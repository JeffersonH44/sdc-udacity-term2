#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "tools.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <tuple>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool isInitialized;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool useLaser;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool useRadar;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x;

  ///* state covariance matrix
  MatrixXd P;

  ///* predicted sigma points matrix
  MatrixXd XSigPred;

  ///* time when the state is true, in us
  long long time_us_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double stdAcc;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double stdYawdd;

  ///* Laser measurement noise standard deviation position1 in m
  double stdLasPx;

  ///* Laser measurement noise standard deviation position2 in m
  double stdLasPy;

  ///* Radar measurement noise standard deviation radius in m
  double stdRadR;

  ///* Radar measurement noise standard deviation angle in rad
  double stdRadPhi;

  ///* Radar measurement noise standard deviation radius change in m/s
  double stdRadRd ;

  ///* Weights of sigma points
  VectorXd weights;

  ///* State dimension
  int nx;

  ///* Augmented state dimension
  int nAug;

  int nRadar;

  int nLidar;

  ///* Sigma point spreading parameter
  double lambda;

  ///* previous timestamp
  long prevTimestamp;

  double lastRadarNIS;
  double lastLidarNIS;

  MatrixXd RRadar;
  MatrixXd RLidar;


  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void processMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void updateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void updateRadar(MeasurementPackage meas_package);

  MatrixXd augmentedSigmaPoints();
  void sigmaPointPrediction(MatrixXd Xsig_aug, double delta_t);
  void predictMeanAndCovariance();
  std::tuple<VectorXd*, MatrixXd*, MatrixXd*> predictRadarMeasurement();
  void updateRadarState(MatrixXd Zsig, VectorXd z_pred, MatrixXd S, VectorXd z);

  std::tuple<VectorXd*, MatrixXd*, MatrixXd*> predictLidarMeasurement();
  void updateLidarState(MatrixXd Zsig, VectorXd z_pred, MatrixXd S, VectorXd z);
};

#endif /* UKF_H */
