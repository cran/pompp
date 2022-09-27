#ifndef __BACKGROUND_VARIABLES_HPP__
#define __BACKGROUND_VARIABLES_HPP__

#include "GaussianProcess.hpp"
#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(openmp)]]

const int INTENSITY_VARIABLES = 0;
const int OBSERVABILITY_VARIABLES = 1;

class BackgroundVariables {
  const std::vector<int> intensityCols, observabilityCols;
  GaussianProcess* spatialProcessInt;
  GaussianProcess* spatialProcessObs;
  bool useGPint = false, useGPobs = false;

  Eigen::MatrixXd getVariablesMat(const Eigen::MatrixXd& coordinates,
                               std::vector<int> columns) {
    int n = coordinates.rows();
    Eigen::MatrixXd out(n, columns.size());
#pragma omp parallel for
    for (int i = 0; i < n; i++)
      out.row(i) = getVariablesVec(coordinates.row(i).transpose(), columns).transpose();

    return out;
  }
  // Retrieve variables for any set of columns
  virtual Eigen::VectorXd getVariablesVec(const Eigen::VectorXd& coordinates,
                                          std::vector<int> columns) = 0;
public:
  // Getters
  Eigen::MatrixXd getVarMat(const Eigen::MatrixXd& coordinates, int type) {
    if (type == INTENSITY_VARIABLES) return getVariablesMat(coordinates, intensityCols);
    if (type == OBSERVABILITY_VARIABLES) return getVariablesMat(coordinates, observabilityCols);
    return Eigen::MatrixXd(0, 0);
  }
  Eigen::VectorXd getVarVec(const Eigen::VectorXd& coordinates,
                            double& mark,
                            double nugget, double mu, int type) {
    if (type == INTENSITY_VARIABLES) {
      Eigen::VectorXd out(intensityCols.size() + (useGPint ? 1 : 0));
      out.head(intensityCols.size()) = getVariablesVec(coordinates, intensityCols);
      if (useGPint)
        out(intensityCols.size()) =
          spatialProcessInt->getNewPoint(coordinates.head(2), mark, nugget, mu);
      return out;
    }
    if (type == OBSERVABILITY_VARIABLES) {
      Eigen::VectorXd out(observabilityCols.size() + (useGPobs ? 1 : 0));
      out.head(observabilityCols.size()) = getVariablesVec(coordinates, observabilityCols);
      if (useGPobs)
        out(observabilityCols.size()) =
          spatialProcessObs->getNewPoint(coordinates.head(2), mark, nugget, mu);
      return out;
    }
    return Eigen::VectorXd(0);
  }
  Eigen::VectorXd getVarVec(const Eigen::VectorXd& coordinates, int type) {
    if (type == INTENSITY_VARIABLES) {
      return getVariablesVec(coordinates, intensityCols);
    }
    if (type == OBSERVABILITY_VARIABLES) {
      return getVariablesVec(coordinates, observabilityCols);
    }
    return Eigen::VectorXd(0);
  }
  Eigen::VectorXd getGP(int type) {
    if (type == INTENSITY_VARIABLES) return spatialProcessInt->getAugmentedValuesTail();
    if (type == OBSERVABILITY_VARIABLES) return spatialProcessObs->getAugmentedValuesTail();
    return Eigen::VectorXd(0);
  }
  Eigen::VectorXd getGPfull(int type) {
    if (type == INTENSITY_VARIABLES) return spatialProcessInt->getAugmentedValues();
    if (type == OBSERVABILITY_VARIABLES) return spatialProcessObs->getAugmentedValues();
    return Eigen::VectorXd(0);
  }
  // Setters
  void setGP(GaussianProcess* gp, int type) {
    if (type == INTENSITY_VARIABLES) {
      spatialProcessInt = gp;
      useGPint = true;
      return;
    }
    if (type == OBSERVABILITY_VARIABLES) {
      spatialProcessObs = gp;
      useGPobs = true;
    }
  }

  // Resampler
  void resampleGPs(double marksMu, double marksVariance,
                   const Eigen::VectorXd& xMarks, Eigen::VectorXd& xPrimeMarks,
                   const Eigen::VectorXd& betasPart, const Eigen::VectorXd& pgs,
                   double gamma) {
    if (useGPint) spatialProcessInt->resampleGP(marksMu, marksVariance,
        xMarks, xPrimeMarks, betasPart, pgs, gamma);
    if (useGPobs) spatialProcessObs->resampleGP(marksMu, marksVariance,
        xMarks, xPrimeMarks, betasPart, pgs, gamma);
  }

  // Random point (for point process simulation)
  virtual Eigen::VectorXd getRandomPoint() = 0;
  Eigen::MatrixXd getRandomPoint(int size) {
    Eigen::MatrixXd output(size, 3);
    for (int i = 0; i < size; i++)
      output.row(i) = getRandomPoint().transpose();
    return output;
  }

  // Related to spatial Gaussian Processes
  void startGPs(int howMany) {
    if (useGPint) spatialProcessInt->startUp(howMany);
    if (useGPobs) spatialProcessObs->startUp(howMany);
  }
  void acceptNewPoint(int type) {
    if (type == INTENSITY_VARIABLES && useGPint)
      spatialProcessInt->acceptNewPoint();
    if (type == OBSERVABILITY_VARIABLES && useGPobs)
      spatialProcessObs->acceptNewPoint();
  }
  void setGPinStone() {
    if (useGPint) spatialProcessInt->closeUp();
    if (useGPobs) spatialProcessObs->closeUp();
  }

  // Constructor and destructor
  BackgroundVariables(std::vector<int> intCols, std::vector<int> obsCols,
                      GaussianProcess* gp) :
    intensityCols(intCols), observabilityCols(obsCols),
    spatialProcessObs(gp), useGPobs(true) {}
  virtual ~BackgroundVariables() {delete spatialProcessObs;}

protected:
  double* data;
};

class MatrixVariables : public BackgroundVariables {
  const long rows, cols, longCol, latCol;
  double halfVertSmallestIncrement, halfHorSmallestIncrement;
  // Uses last coordinate (position in covariates matrix) to retrieve covariates.
  Eigen::VectorXd getVariablesVec(const Eigen::VectorXd& coordinates,
                                  std::vector<int> columns) {
    Eigen::VectorXd out(columns.size());
    for (int j = 0; j < columns.size(); j++)
      out(j) = data[columns[j] * rows + (long)coordinates(2)];
    return out;
  }
public:
  MatrixVariables(std::vector<int> intCols, std::vector<int> obsCols,
                  SEXP matrix, int xC, int yC, GaussianProcess* gp) :
                  BackgroundVariables(intCols, obsCols, gp),
                  rows(INTEGER(Rf_getAttrib(matrix, R_DimSymbol))[0]),
                  cols(INTEGER(Rf_getAttrib(matrix, R_DimSymbol))[1]),
                  longCol(xC), latCol(yC) {
    data = REAL(matrix);

    // Used to get a random point by adding a small perturbation around the selected cell center.
    halfVertSmallestIncrement = 0;
    halfHorSmallestIncrement = 0;
    double temp;
    bool testVert = true, testHor = true;
    for (long i = 1; i < rows; i++) {
      temp = fabs(data[latCol * rows] - data[latCol * rows + i]);
      if ((testVert && halfVertSmallestIncrement < temp) ||
          (!testVert && temp && halfVertSmallestIncrement > temp)) {
        testVert = false;
        halfVertSmallestIncrement = temp;
      }
      temp = fabs(data[longCol * rows] - data[longCol * rows + i]);
      if ((testHor && halfHorSmallestIncrement < temp) ||
          (!testHor && temp && halfHorSmallestIncrement > temp)) {
        testHor = false;
        halfHorSmallestIncrement = temp;
      }
    }
    halfVertSmallestIncrement /= 2;
    halfHorSmallestIncrement /= 2;
  }

  // Third coordinate is the random row. Used later to retrieve covariates values easily in getVariablesVec().
  Eigen::VectorXd getRandomPoint() {
    Eigen::VectorXd out(3);
    long selectedRow = long(R::runif(0, 1) * rows);
    out(0) = data[longCol * rows + selectedRow] + R::runif(-1, 1) * halfHorSmallestIncrement;
    out(1) = data[latCol * rows + selectedRow] + R::runif(-1, 1) * halfVertSmallestIncrement;
    out(2) = selectedRow;
    return out;
  }
};

#endif
