#ifndef __pompp_BINARY_REGRESSION_HPP__
#define __pompp_BINARY_REGRESSION_HPP__

#include <RcppEigen.h>
#include "RegressionPrior.hpp"

class BinaryRegression {
protected:
  Eigen::VectorXd betas, normalMean; // normalMean is used in the GP part of the program.
  RegressionPrior* prior;
  const unsigned int n;

public:
  virtual double sample(const Eigen::MatrixXd& onesCovariates,
                        const Eigen::MatrixXd& zerosCovariates) = 0;


  // Link function. Returns in the **LOG** scale. complementaryProb = true
  // calculates the link for 1 - p
  // For the stored beta vector
  Eigen::VectorXd link(const Eigen::MatrixXd& covariates,
                       bool complementaryProb = false) {
    return link(covariates, betas, complementaryProb);
  }

  // For any beta vector
  virtual Eigen::VectorXd link(const Eigen::MatrixXd& covariates,
                               const Eigen::VectorXd& beta,
                               bool complementaryProb = false) = 0;

  // Constructor and destructor
  BinaryRegression(Eigen::VectorXd initialize,
                   RegressionPrior* p) : betas(initialize),
                   normalMean(Eigen::VectorXd(0)), prior(p),
                   n(initialize.size()) {}
  virtual ~BinaryRegression() {delete prior;}

  // Prior setter
  void setPrior(RegressionPrior* p) {prior = p;}

  // Some getters
  Eigen::VectorXd getBeta() {return betas;}
  int getSize() {return n;}
  Eigen::VectorXd getNormalMean() {return normalMean;}
  virtual Eigen::VectorXd getDataAugmentation() = 0; // For data augmentation variables
  // Some setters
  void setNormalMean(Eigen::VectorXd newValue) {normalMean = newValue;}
  void setBeta(Eigen::VectorXd newValue) {betas = newValue;}
};

class LogisticRegression : public BinaryRegression {
  // Data augmentation
  Eigen::VectorXd pg;
public:
  Eigen::VectorXd getDataAugmentation() {return pg;}

  LogisticRegression(Eigen::VectorXd initialize, RegressionPrior* p) :
    BinaryRegression(initialize, p), pg(Eigen::VectorXd(0)) {}

  double sample(const Eigen::MatrixXd& onesCovariates,
                const Eigen::MatrixXd& zerosCovariates);
  Eigen::VectorXd link(const Eigen::MatrixXd& covariates,
                       const Eigen::VectorXd& beta,
                       bool complementaryProb = false);
};

#endif
