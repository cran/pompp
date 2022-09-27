#include <RcppEigen.h>
#include "include/PresenceOnly.hpp"
#include "include/BackgroundVariables.hpp"
#include "include/GaussianProcess.hpp"
#include "include/CovarianceFunction.hpp"
#include "include/BinaryRegression.hpp"
#include "include/RegressionPrior.hpp"
#include <progress.hpp>
using namespace Rcpp;

// [[Rcpp::export]]
List cppPOMPP(Eigen::VectorXd beta, Eigen::VectorXd delta,
                   double lambda, Rcpp::String b_updater,
                   Rcpp::String d_updater, Rcpp::String l_updater,
                   Rcpp::List parB, Rcpp::List parD,
                   double lambdaA, double lambdaB,
                   Rcpp::String covsClass, SEXP covariates,
                   double areaD, Rcpp::String xClass,
                   double mu, double nugget,
                   double marksMuMu, double marksMuS2,
                   double marksNuggetA, double marksNuggetB,
                   Eigen::MatrixXd xValues, Eigen::VectorXd xMarks,
                   Eigen::MatrixXd xPositions,
                   Eigen::VectorXi intensityCovs,
                   Eigen::VectorXi observabilityCovs,
                   Eigen::VectorXi xIntensityCovs,
                   Eigen::VectorXi xObservabilityCovs,
                   double maxDist, double sigma2, double phi,
                   int neighborhoodSize,
                   int longCol, int latCol,
                   int burnin, int thin, int iter, int threads, bool verbose) {
  int i, j;

#ifdef _OPENMP
  omp_set_num_threads( threads );
#endif

  // Auxiliary
  Eigen::MatrixXd xInt(xValues.rows(), xIntensityCovs.size());
  for (i = 0; i < xValues.rows(); i++)
    for (j = 0; j < xIntensityCovs.size(); j++)
      xInt(i, j) = xValues(i, xIntensityCovs(j));
  Eigen::MatrixXd xObs(xValues.rows(), xObservabilityCovs.size());
  for (i = 0; i < xValues.rows(); i++)
    for (j = 0; j < xObservabilityCovs.size(); j++)
      xObs(i, j) = xValues(i, xObservabilityCovs(j));

  // Output storage
  int outSize = iter / thin;
  Eigen::MatrixXd outBetas(outSize, beta.size());
  Eigen::MatrixXd outDeltas(outSize, delta.size());
  Eigen::VectorXd outLambdas(outSize);
  Eigen::VectorXd outMus(outSize);
  Eigen::VectorXd outNuggets(outSize);
  Eigen::VectorXd outLogPost(outSize);
  Eigen::VectorXd out_nU(outSize);
  Eigen::VectorXd out_nXp(outSize);
  Eigen::VectorXd outMarksPrimeSum(outSize);
  Eigen::VectorXd outMarksPrimeVariance(outSize);
  Eigen::VectorXd outAllMarksSum(outSize);
  Eigen::VectorXd outAllMarksVariance(outSize);

  // Get prior parameters
  Eigen::VectorXd muB = parB["mean"];
  Eigen::MatrixXd SigmaB = parB["covariance"];
  Eigen::VectorXd muD = parD["mean"];
  Eigen::MatrixXd SigmaD = parD["covariance"];

  PresenceOnly mc(
    xPositions, xInt, xObs,
    new MatrixVariables(
      std::vector<int>(&intensityCovs[0], intensityCovs.data() + intensityCovs.size()),
      std::vector<int>(&observabilityCovs[0], observabilityCovs.data() + observabilityCovs.size()),
      covariates, longCol, latCol,
      new NNGP(xPositions, xPositions.rows(), neighborhoodSize,
               new PowerExponentialCovariance(maxDist, phi, 1, sigma2))
    ), xMarks,
    new LogisticRegression(beta, new NormalPrior(muB, SigmaB)),
    new LogisticRegression(delta, new NormalPrior(muD, SigmaD)),
    lambda, lambdaA, lambdaB,
    mu, nugget,
    areaD, marksMuMu, marksMuS2, marksNuggetA, marksNuggetB
  );

  // Burning in
  if (burnin) {
    if (verbose) Rcpp::Rcout << "Warming up the Markov Chain.\n";
    Progress progr_Burnin(burnin, true);
    for (i = 0; i < burnin; i++)
    {
      progr_Burnin.increment();
      mc.update();
    }
    if (verbose) Rcpp::Rcout << "Warm up complete. ";
  }

  if (verbose) Rcpp::Rcout << "Sampling MCMC.\n";
  double xMarksPrimeSquaredNorm, xMarksSquaredNorm = xMarks.squaredNorm();
  int fullSize;
  Progress progr_Main(outSize, true);
  for (i = 0; i < outSize; i++)
  {
    R_CheckUserInterrupt();
    mc.update(thin);
    outBetas.row(i) = mc.getBeta().transpose();
    outDeltas.row(i) = mc.getDelta().transpose();
    outLambdas[i] = mc.getLambdaStar();
    outMus[i] = mc.getMarksMu();
    outNuggets[i] = mc.getMarksNugget();
    out_nU[i] = mc.getUsize();
    out_nXp[i] = mc.getXpsize();
    fullSize = xMarks.size() + mc.getXpsize();
    outMarksPrimeSum[i] = mc.getMarksPrime().size() ? mc.getMarksPrime().sum() : 0;
    xMarksPrimeSquaredNorm = mc.getMarksPrime().squaredNorm();
    outMarksPrimeVariance[i] = xMarksPrimeSquaredNorm / mc.getXpsize() -
      outMarksPrimeSum[i] * outMarksPrimeSum[i] / (out_nXp[i] * out_nXp[i]);
    outAllMarksSum[i] = outMarksPrimeSum[i] + xMarks.sum();
    outAllMarksVariance[i] = (xMarksPrimeSquaredNorm + xMarksSquaredNorm) / fullSize -
      outAllMarksSum[i] * outAllMarksSum[i] / (fullSize * fullSize);
    outLogPost[i] = mc.getLogPosterior();
    progr_Main.increment();
  }
  if (verbose) Rcpp::Rcout << "MCMC sampling complete.\n";

//  delete mc;

  return Rcpp::List::create(Rcpp::Named("beta") = outBetas,
                            Rcpp::Named("delta") = outDeltas,
                            Rcpp::Named("lambda") = outLambdas,
                            Rcpp::Named("mu") = outMus,
                            Rcpp::Named("nugget") = outNuggets,
                            Rcpp::Named("nU") = out_nU,
                            Rcpp::Named("nXp") = out_nXp,
                            Rcpp::Named("marksPrimeSum") = outMarksPrimeSum,
                            Rcpp::Named("marksPrimeVariance") = outMarksPrimeVariance,
                            Rcpp::Named("allMarksSum") = outAllMarksSum,
                            Rcpp::Named("allMarksVariance") = outAllMarksVariance,
                            Rcpp::Named("logPost") = outLogPost);
}
