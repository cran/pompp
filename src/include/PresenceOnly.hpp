#ifndef __PRESENCE_ONLY_HPP__
#define __PRESENCE_ONLY_HPP__

#include "MarkovChain.hpp"
#include "BinaryRegression.hpp"
#include "BackgroundVariables.hpp"
#include "RegressionPrior.hpp"

#define MAX_ATTEMPTS_GP 1000

class PresenceOnly : public MarkovChain {
  // States
  // Regression members
  BinaryRegression* beta;
  BinaryRegression* delta;

  // Point process members
  const double area;
  Eigen::MatrixXd xprime, u;
  Eigen::MatrixXd xxprimeIntensity, xprimeObservability, uIntensity;
  double sampleProcesses();

  // Data
  const Eigen::MatrixXd x;
  const Eigen::MatrixXd xIntensity;
  Eigen::MatrixXd xObservability;
  BackgroundVariables* bkg;

  // lambda star members
  double lambdaStar;
  double aL, bL;
  double updateLambdaStar();

  // Marks members
  const Eigen::VectorXd marks;
  const double marksMuPriormu, marksMuPriors2,
  marksNuggetPriora, marksNuggetPriorb;
  Eigen::VectorXd marksPrime;
  double marksMu, marksNugget;
  double updateMarksPars(const Eigen::VectorXd& gp);

  // Inherited
  double applyTransitionKernel();

public:
  PresenceOnly(const Eigen::MatrixXd& xPositions,
               const Eigen::MatrixXd& xIntensityCovs,
               const Eigen::MatrixXd& xObservabilityCovs,
               BackgroundVariables* bk,
               const Eigen::VectorXd& observedValues,
               BinaryRegression* b, BinaryRegression* d,
               double lambda, double lambdaA, double lambdaB,
               double mu, double nugget,
               double a, double mmm, double mms2,
               double mna, double mnb) : MarkovChain(), beta(b), delta(d),
               area(a), x(xPositions),
               xIntensity(xIntensityCovs),
               bkg(bk),
               lambdaStar(lambda), aL(lambdaA), bL(lambdaB),
               marks(observedValues), marksMuPriormu(mmm), marksMuPriors2(mms2),
               marksNuggetPriora(mna), marksNuggetPriorb(mnb),
               marksMu(mu), marksNugget(nugget) {
    xObservability = Eigen::MatrixXd(xObservabilityCovs.rows(), xObservabilityCovs.cols() + 1);
    xObservability.leftCols(xObservabilityCovs.cols()) = xObservabilityCovs;
    xxprimeIntensity = xIntensity;
    xprimeObservability = Eigen::MatrixXd(0, 0);
    uIntensity = Eigen::MatrixXd(0, 0);
  }
  ~PresenceOnly() {delete beta; delete delta; delete bkg;}

  // Getters
  Eigen::VectorXd getBeta() {return beta->getBeta();}
  Eigen::VectorXd getDelta() {return delta->getBeta();}
  double getLambdaStar() {return lambdaStar;}
  Eigen::MatrixXd getU() {return u;}
  Eigen::MatrixXd getXprime() {return xprime;}
  int getUsize() {return u.rows();}
  int getXpsize() {return xprime.rows();}
  Eigen::VectorXd getMarksPrime() {return marksPrime;}
  double getMarksMu() {return marksMu;}
  double getMarksNugget() {return marksNugget;}
};

#endif
