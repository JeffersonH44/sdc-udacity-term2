#ifndef FusionEKF_H_
#define FusionEKF_H_

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "kalman_filter.h"
#include "tools.h"

class FusionEKF {
public:
  /**
  * Constructor.
  */
  FusionEKF();

  /**
  * Destructor.
  */
  virtual ~FusionEKF();

  /**
  * Run the whole flow of the Kalman Filter from here.
  */
  void ProcessMeasurement(const MeasurementPackage &measurement_pack);

  /**
  * Kalman Filter update and prediction math lives in here.
  */
  KalmanFilter ekf;

private:
  // check whether the tracking toolbox was initiallized or not (first measurement)
  bool isInitialized;

  // previous timestamp
  long long previousTimestamp;

  // tool object used to compute Jacobian and RMSE
  Eigen::MatrixXd RLaser;
  Eigen::MatrixXd RRadar;
  Eigen::MatrixXd HLaser;
  Eigen::MatrixXd Hj;
};

#endif /* FusionEKF_H_ */
