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
   * @param measurementPackage The latest measurement data of either radar or laser
   */
  void processMeasurement(MeasurementPackage measurementPackage);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param deltaT Time between k and k+1 in s
   */
  void prediction(double deltaT);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param measurementPackage The measurement at k+1
   */
  void updateLidar(MeasurementPackage measurementPackage);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param measurementPackage The measurement at k+1
   */
  void updateRadar(MeasurementPackage measurementPackage);

  MatrixXd augmentedSigmaPoints();
  void sigmaPointPrediction(MatrixXd XSigmaAug, double deltaT);
  void predictMeanAndCovariance();
  std::tuple<VectorXd*, MatrixXd*, MatrixXd*> predictMeasurement(bool isRadar);
  void updateState(MatrixXd ZSigmaPoints, VectorXd zPred, MatrixXd S, VectorXd z, bool isRadar);

  const double &getX(int i) const;
private:
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
  // long long time_us_;
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
  double stdRadRd;
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
  ///* NIS for radar and lidar
  double lastRadarNIS;
  double lastLidarNIS;
  ///* Noise matrix for radar and lidar
  MatrixXd RRadar;
  MatrixXd RLidar;
};

#endif /* UKF_H */
