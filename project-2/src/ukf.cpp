#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  useLaser = true;

  // if this is false, radar measurements will be ignored (except during init)
  useRadar = true;

  // predictions variables to make
  nx = 5;
  // augmented by the two random variables
  nAug = 7;

  // initial state vector
  x = VectorXd(nx);

  // initial covariance matrix
  P = MatrixXd(nx, nx);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  stdAcc = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  stdYawdd = 0.7;

  // Laser measurement noise standard deviation position1 in m
  stdLasPx = 0.15;

  // Laser measurement noise standard deviation position2 in m
  stdLasPy = 0.15;

  // Radar measurement noise standard deviation radius in m
  stdRadR = 0.3;

  // Radar measurement noise standard deviation angle in rad
  stdRadPhi = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  stdRadRd = 0.3;

  nRadar = 3;

  nLidar = 2;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  // spreading parameter
  lambda = 3 - nAug;

  // weights for the sigma points
  weights = VectorXd(2 * nAug +1);
  weights.fill(1.0 / (2.0*(lambda + nAug)));
  weights(0) = lambda / (lambda + nAug);

  lastLidarNIS = 0.0;
  lastRadarNIS = 0.0;

  RLidar = MatrixXd(nLidar, nLidar);
  RLidar << stdLasPx * stdLasPx, 0,
             0, stdLasPy * stdLasPy;

  RRadar = MatrixXd(nRadar, nRadar);
  RRadar << stdRadR * stdRadR, 0, 0,
             0, stdRadPhi * stdRadPhi, 0,
             0, 0, stdRadRd * stdRadRd;

  isInitialized = false;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::processMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  if(isInitialized) {
    double dt = (meas_package.timestamp_ - prevTimestamp) / 1000000.0;
    this->prediction(dt);
    if(useLaser && meas_package.sensor_type_ == MeasurementPackage::SensorType::LASER) {
      this->updateLidar(meas_package);
    }
    if(useRadar && meas_package.sensor_type_ == MeasurementPackage::SensorType::RADAR) {
      this->updateRadar(meas_package);
    }
    prevTimestamp = meas_package.timestamp_;
  } else {

    prevTimestamp = meas_package.timestamp_;
    if(meas_package.sensor_type_ == MeasurementPackage::SensorType::LASER) {
      double px = meas_package.raw_measurements_[0];
      double py = meas_package.raw_measurements_[1];
      x << px, py, 0.0, 0.0, 0.0;
    } else {
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double rho_dot = meas_package.raw_measurements_[2];
      x << rho * cos(phi), rho * sin(phi), 0.0, 0.0, 0.0;
    }

    // TODO: do not forget to change those values
    P << 50, 0, 0, 0, 0,
          0, 50, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 0.5, 0,
          0, 0, 0, 0, 0.5;

    isInitialized = true;
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  MatrixXd Xsig_aug = this->augmentedSigmaPoints();
  this->sigmaPointPrediction(Xsig_aug, delta_t);
  this->predictMeanAndCovariance();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::updateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  std::tuple<VectorXd*, MatrixXd*, MatrixXd*> tup = this->predictLidarMeasurement();

  VectorXd *z_pred = std::get<0>(tup);
  MatrixXd *S = std::get<1>(tup);
  MatrixXd *Zsig = std::get<2>(tup);

  /*std::cout << *z_pred << std::endl;
  std::cout << *S << std::endl;
  std::cout << *Zsig << std::endl;*/
  this->updateLidarState(*Zsig, *z_pred, *S, meas_package.raw_measurements_);

  delete(z_pred);
  delete(S);
  delete(Zsig);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::updateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  std::tuple<VectorXd*, MatrixXd*, MatrixXd*> tup = this->predictRadarMeasurement();

  VectorXd *z_pred = std::get<0>(tup);
  MatrixXd *S = std::get<1>(tup);
  MatrixXd *Zsig = std::get<2>(tup);

  /*std::cout << *z_pred << std::endl;
  std::cout << *S << std::endl;
  std::cout << *Zsig << std::endl;*/
  this->updateRadarState(*Zsig, *z_pred, *S, meas_package.raw_measurements_);

  delete(z_pred);
  delete(S);
  delete(Zsig);
}

MatrixXd UKF::augmentedSigmaPoints() {

  //create augmented mean vector
  VectorXd x_aug = VectorXd(nAug);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(nAug, nAug);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(nAug, 2 * nAug + 1);

/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //create augmented mean state
  x_aug.head(5) = x;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P;
  P_aug(5,5) = stdAcc * stdAcc;
  P_aug(6,6) = stdYawdd * stdYawdd;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< nAug; i++)
  {
    Xsig_aug.col(i + 1)       = x_aug + sqrt(lambda + nAug) * L.col(i);
    Xsig_aug.col(i + 1 + nAug) = x_aug - sqrt(lambda + nAug) * L.col(i);
  }

/*******************************************************************************
 * Student part end
 ******************************************************************************/

  //print result
  // std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;

  //write result
  return Xsig_aug;
}

void UKF::sigmaPointPrediction(MatrixXd Xsig_aug, double delta_t) {

  //create matrix with predicted sigma points as columns
  XSigPred = MatrixXd(nx, 2 * nAug + 1);
/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //predict sigma points
  //avoid division by zero
  //write predicted sigma points into right column
  for (int i = 0; i< 2 * nAug + 1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
      px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
      py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
      px_p = p_x + v*delta_t*cos(yaw);
      py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    XSigPred(0,i) = px_p;
    XSigPred(1,i) = py_p;
    XSigPred(2,i) = v_p;
    XSigPred(3,i) = yaw_p;
    XSigPred(4,i) = yawd_p;
  }


/*******************************************************************************
 * Student part end
 ******************************************************************************/
}

void UKF::predictMeanAndCovariance() {
/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //create vector for predicted state
  VectorXd predictedState = VectorXd(nx);

  //create covariance matrix for prediction
  MatrixXd covMatPred = MatrixXd::Zero(nx, nx);
  //predict state mean

  predictedState.fill(0.0);
  for (int i = 0; i < 2 * nAug + 1; i++) {  //iterate over sigma points
    predictedState = predictedState+ weights(i) * XSigPred.col(i);
  }
  //predict state covariance matrix
  for (int i = 0; i < 2 * nAug + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = XSigPred.col(i) - predictedState;
    //angle normalization
    x_diff(3) = Tools::angleNormalization(x_diff(3));

    covMatPred = covMatPred + weights(i) * x_diff * x_diff.transpose() ;
  }

  x = predictedState;
  P = covMatPred;
/*******************************************************************************
 * Student part end
 ******************************************************************************/
}



std::tuple<VectorXd*, MatrixXd*, MatrixXd*> UKF::predictRadarMeasurement() {

  //create matrix for sigma points in measurement space
  MatrixXd *Zsig = new MatrixXd(nRadar, 2 * nAug + 1);

  //mean predicted measurement
  VectorXd *z_pred = new VectorXd(nRadar);

  //measurement covariance matrix S
  MatrixXd *S = new MatrixXd(nRadar, nRadar);

/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //transform sigma points into measurement space
  for(int i = 0; i < 2 * nAug + 1; ++i) {
    VectorXd currentPred = VectorXd(nRadar);
    double px = XSigPred(0, i);
    double py = XSigPred(1, i);
    double v = XSigPred(2, i);
    double yaw = XSigPred(3, i);
    yaw = Tools::angleNormalization(yaw);

    double sqrtVal = sqrt(px*px + py*py);
    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    currentPred << sqrtVal,
            atan2(py, px),
            (px*v1 + py*v2) / sqrtVal;

    (*Zsig).col(i) = currentPred;
  }
  //calculate mean predicted measurement
  (*z_pred).fill(0.0);
  for (int i=0; i < 2*nAug+1; i++) {
    *z_pred = *z_pred + weights(i) * (*Zsig).col(i);
  }
  //calculate measurement covariance matrix S
  (*S).fill(0.0);
  for(int i = 0; i < 2 * nAug + 1; ++i) {
    VectorXd Zdiff = (*Zsig).col(i) - *z_pred;

    Zdiff(1) = Tools::angleNormalization(Zdiff(1));

    *S += weights(i) * Zdiff * Zdiff.transpose();
  }

  *S = *S + RRadar;

/*******************************************************************************
 * Student part end
 ******************************************************************************/

  return std::make_tuple(z_pred, S, Zsig);
}

void UKF::updateRadarState(MatrixXd Zsig, VectorXd z_pred, MatrixXd S, VectorXd z) {

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(nx, nRadar);

/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * nAug + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    z_diff(1) = Tools::angleNormalization(z_diff(1));

    // state difference
    VectorXd x_diff = XSigPred.col(i) - x;
    //angle normalization
    x_diff(3) = Tools::angleNormalization(x_diff(3));

    Tc = Tc + weights(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  z_diff(1) = Tools::angleNormalization(z_diff(1));

  //update state mean and covariance matrix
  x = x + K * z_diff;
  P = P - K*S*K.transpose();

  lastRadarNIS = Tools::calculateNIS(z_diff, S);
/*******************************************************************************
 * Student part end
 ******************************************************************************/
}

std::tuple<VectorXd *, MatrixXd *, MatrixXd *> UKF::predictLidarMeasurement() {
  //create matrix for sigma points in measurement space
  MatrixXd *Zsig = new MatrixXd(nLidar, 2 * nAug + 1);

  //mean predicted measurement
  VectorXd *z_pred = new VectorXd(nLidar);

  //measurement covariance matrix S
  MatrixXd *S = new MatrixXd(nLidar, nLidar);

/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //transform sigma points into measurement space
  for(int i = 0; i < 2 * nAug + 1; ++i) {
    VectorXd currentPred = VectorXd(nLidar);
    double px = XSigPred(0, i);
    double py = XSigPred(1, i);
    double v = XSigPred(2, i);
    double yaw = XSigPred(3, i);
    yaw = Tools::angleNormalization(yaw);

    double sqrtVal = sqrt(px*px + py*py);
    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    currentPred << px, py; // for laser

    (*Zsig).col(i) = currentPred;
  }
  //calculate mean predicted measurement
  (*z_pred).fill(0.0);
  for (int i=0; i < 2*nAug+1; i++) {
    *z_pred = *z_pred + weights(i) * (*Zsig).col(i);
  }
  //calculate measurement covariance matrix S
  (*S).fill(0.0);
  for(int i = 0; i < 2 * nAug + 1; ++i) {
    VectorXd Zdiff = (*Zsig).col(i) - *z_pred;

    // Zdiff(1) = Tools::angleNormalization(Zdiff(1));

    *S += weights(i) * Zdiff * Zdiff.transpose();
  }

  *S = *S + RLidar;

/*******************************************************************************
 * Student part end
 ******************************************************************************/

  return std::make_tuple(z_pred, S, Zsig);
}

void UKF::updateLidarState(MatrixXd Zsig, VectorXd z_pred, MatrixXd S, VectorXd z) {
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(nx, nLidar);

/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * nAug + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization just for radar
    // z_diff(1) = Tools::angleNormalization(z_diff(1));

    // state difference
    VectorXd x_diff = XSigPred.col(i) - x;
    //angle normalization because X contains the 5 variables
    x_diff(3) = Tools::angleNormalization(x_diff(3));

    Tc = Tc + weights(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization in case of radar
  //z_diff(1) = Tools::angleNormalization(z_diff(1));

  //update state mean and covariance matrix
  x = x + K * z_diff;
  P = P - K*S*K.transpose();

  lastLidarNIS = Tools::calculateNIS(z_diff, S);
}


