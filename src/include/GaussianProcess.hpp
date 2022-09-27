#ifndef __POMPP_GAUSSIANPROCESS_H__
#define __POMPP_GAUSSIANPROCESS_H__

#include <RcppEigen.h>
#include "CovarianceFunction.hpp"
#include "safeR.hpp"

class GaussianProcess {
  // coordinates parameters. Function returns the sampled value
  virtual void sampleNewPoint(Eigen::VectorXd coords,
                              double& mark, double nugget, double mu);

  double updateCovarianceParameters();
  Eigen::MatrixXd recalcPrecision(std::vector<double> newParams); // Used in updateCovarianceParameter()
public:
  // getters
  Eigen::VectorXd getAugmentedValues() {return augmentedValues;}
  Eigen::VectorXd getAugmentedValuesTail() {
    return augmentedValues.tail(augmentedValues.size() - xSize);
  }
  // setters
  void setCovFunction(CovarianceFunction* c) {covFun = c;}

  GaussianProcess(int s) : xSize(s) {}
  GaussianProcess(Eigen::MatrixXd pos, int s,
                  CovarianceFunction* cf) : xSize(s),
  positions(pos.leftCols(2)),
  values(rnorm(xSize)),
  covFun(cf) {augmentedValues = values;}
  virtual ~GaussianProcess() {delete covFun;}

  double getNewPoint(Eigen::VectorXd coords, double& mark,
                     double nugget, double mu)
    {sampleNewPoint(coords, mark, nugget, mu); return propValue;}
  virtual void acceptNewPoint();
  virtual void resampleGP(double marksMu, double marksVariance,
                          const Eigen::VectorXd& xMarks, Eigen::VectorXd& xPrimeMarks,
                          const Eigen::VectorXd& betasPart, const Eigen::VectorXd& pgs,
                          double gamma);

  // Methods to update which points are data augmentation.
  virtual void startUp(int howMany);
  virtual void closeUp();
protected:
  const int xSize; // Used in start up and close up
  int tempAcc; // Used in start up and close up
  int parameterSize, currentIndex;
  Eigen::MatrixXd positions, covariances, augmentedPositions, augmentedCovariances;
  Eigen::VectorXd augmentedValues;
  Eigen::VectorXd values;
  CovarianceFunction* covFun;
  double logDensity;

  // Proposed point
  Eigen::VectorXd propCovariances;
  double propValue;

  double calcDist(Eigen::VectorXd p1, Eigen::VectorXd p2);
};

class NNGP : public GaussianProcess {
  void sampleNewPoint(Eigen::VectorXd coords,
                      double& mark,
                      double nugget, double mu);
  Eigen::MatrixXd recalcPrecision(std::vector<double> newParams); // Used in updateCovarianceParameter()
  void bootUpIminusA();

  // Neighborhood members
  const int neighborhoodSize;
  std::vector<int> neighborhood;
  Eigen::VectorXd distances, D, Arow, propPosition, theseCovariances;
  std::vector<int> getNeighorhood(Eigen::VectorXd coords);
  Eigen::SparseMatrix<double> IminusA, precision;
  std::vector<Eigen::Triplet<double> > trips;
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > sqrtC;
  Eigen::MatrixXi pastCovariancesPositions;
  Eigen::MatrixXd pastCovariances, propPrecision;
  int thisPosition;
  double propD;

public:
  NNGP(Eigen::MatrixXd pos, int s, int M, CovarianceFunction* cf) :
  GaussianProcess(pos, s, cf), neighborhoodSize(M) {}

  void acceptNewPoint();
  // Methods to update which points are data augmentation.
  void startUp(int howMany);
  void closeUp();

  void resampleGP(double marksMu, double marksVariance,
                  const Eigen::VectorXd& xMarks, Eigen::VectorXd& xPrimeMarks,
                  const Eigen::VectorXd& betasPart, const Eigen::VectorXd& pgs,
                  double gamma);
};

#endif
