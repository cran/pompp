#include <RcppEigen.h>
#include "include/PresenceOnly.hpp"
#include "include/BackgroundVariables.hpp"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "include/safeR.hpp"

// [[Rcpp::plugins(openmp)]]

double PresenceOnly::updateLambdaStar() {
  double a = aL + x.rows() + xprime.rows() + u.rows(),
    b = bL + area;
  lambdaStar = rgamma(a, 1 / b);

  return - lambdaStar * b + (a - 1) * log(lambdaStar);
}

double PresenceOnly::sampleProcesses() {
  // Setting after the resampling
  xObservability.col(xObservability.cols() - 1) =
    bkg->getGPfull(OBSERVABILITY_VARIABLES).head(x.rows());

  // Determining number of points in X' and U
  /*
   * Technically, the correct code would sample from a truncated Poisson.
   * This approximation should only be problematic if the data comes from a
   * homogeneous process, which would defeat the purpose of the analysis anyway
   * so there is no loss.
   */
  long totalPoints = rpois(lambdaStar * area);

  // Sampling from X' and U
  double p, q, uniform;
  long accXp = 0, accU = 0;
  Eigen::VectorXd candidate;
  Eigen::MatrixXd storingCoords(totalPoints, 3); // Put X' on top and U on bottom
  marksPrime = Eigen::VectorXd(totalPoints);
  bkg->startGPs(totalPoints);
  for (int i = 0; i < totalPoints; i++) {
#pragma omp critical
    R_CheckUserInterrupt();
    candidate = bkg->getRandomPoint();
    uniform = log(runif(0, 1));
    q = beta->link(bkg->getVarVec(candidate,
                                  INTENSITY_VARIABLES).transpose())(0);
    if (uniform > q) { // Assign to U
      storingCoords.row(totalPoints - ++accU) = candidate.transpose();
      bkg->acceptNewPoint(INTENSITY_VARIABLES);
    } else {
      p = delta->link(bkg->getVarVec(candidate,
                                     marksPrime(accXp),
                                     marksNugget, marksMu,
                                     OBSERVABILITY_VARIABLES).transpose())(0);
      if (uniform > p + q) { // Assign to X'
        storingCoords.row(accXp++) = candidate.transpose();
        bkg->acceptNewPoint(OBSERVABILITY_VARIABLES);
      } // Else discard candidate and try again.
    }
  }
  bkg->setGPinStone();

  marksPrime.conservativeResize(accXp);

  xxprimeIntensity.resize(x.rows() + accXp, beta->getSize() - 1);
  xxprimeIntensity.topRows(x.rows()) = xIntensity;

  if (accXp) {
    xprime = storingCoords.topRows(accXp);
    xxprimeIntensity.bottomRows(accXp) =
      bkg->getVarMat(xprime, INTENSITY_VARIABLES);
    xprimeObservability.resize(accXp, delta->getSize() - 1);
    xprimeObservability.leftCols(delta->getSize() - 2) =
      bkg->getVarMat(xprime, OBSERVABILITY_VARIABLES);
    xprimeObservability.rightCols(1) =
      bkg->getGP(OBSERVABILITY_VARIABLES);
  } else {
    xprime.resize(0, 3);
    xprimeObservability.resize(0, delta->getSize() - 1);
  }

  if (accU) {
    u = storingCoords.bottomRows(accU);
    uIntensity = bkg->getVarMat(u, INTENSITY_VARIABLES);
  } else {
    u.resize(0, 3);
    uIntensity.resize(0, beta->getSize() - 1);
  }

  return - lgamma(accXp + 1) - lgamma(accU + 1);
}

double PresenceOnly::updateMarksPars(const Eigen::VectorXd& gp) {
  double sqrtMarksNugget = sqrt(marksNugget);

  Eigen::VectorXd logMarks = -gp;
  logMarks.head(x.rows()) += marks.array().log().matrix();
  logMarks.tail(xprime.rows()) += marksPrime.array().log().matrix();

  // Sampling the nugget
  marksNugget = 1 / rgamma(marksNuggetPriora + gp.size() / 2,
                           1 / (marksNuggetPriorb +
                             (logMarks.array() - marksMu).square().sum() / 2));

  // Sampling the mean parameter
  double newVariance = 1 / (1 / marksMuPriors2 + gp.size() / marksNugget);
  marksMu = rnorm(
    newVariance *
      (marksMuPriormu / marksMuPriors2 + logMarks.sum() / marksNugget),
      sqrt(newVariance)
  );

  return 0.;
}

inline double PresenceOnly::applyTransitionKernel() {
  double out, privateOut1, privateOut2;
  out = sampleProcesses();
  out += updateLambdaStar();
#pragma omp parallel
{
#pragma omp sections
{
#pragma omp section
  privateOut1 = beta->sample(xxprimeIntensity, uIntensity);
#pragma omp section
  privateOut2 = delta->sample(xObservability, xprimeObservability);
}
}
  bkg->resampleGPs(marksMu, marksNugget,
                 marks, marksPrime, delta->getNormalMean(),
                 delta->getDataAugmentation(),
                 delta->getBeta()(delta->getSize() - 1));
  out += updateMarksPars(bkg->getGPfull(OBSERVABILITY_VARIABLES));
  return out + privateOut1 + privateOut2;
}

