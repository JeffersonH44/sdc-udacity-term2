#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

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
  static Eigen::VectorXd calculateRMSE(const std::vector<Eigen::VectorXd> &estimations,
                                       const std::vector<Eigen::VectorXd> &ground_truth);

  /**
  * A helper method to calculate Jacobians.
  */
  static Eigen::MatrixXd calculateJacobian(const Eigen::VectorXd &x_state);

};

#endif /* TOOLS_H_ */
