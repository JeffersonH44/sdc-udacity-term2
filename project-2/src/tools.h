#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include <cmath>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  static VectorXd calculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);

  /**
   * A helper method to calculate NIS
   */
  static double calculateNIS(Eigen::VectorXd zDiff, Eigen::MatrixXd S);

  /**
   * Angle normalization
   */
  static double angleNormalization(double angle);
};

#endif /* TOOLS_H_ */