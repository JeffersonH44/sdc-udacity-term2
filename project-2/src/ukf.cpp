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

  // variables for radar (rho, phi, rho dot)
  nRadar = 3;
  // variables for lidar (px, py)
  nLidar = 2;

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
void UKF::processMeasurement(MeasurementPackage measurementPackage) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  if(isInitialized) {
    bool isRadar = measurementPackage.sensor_type == MeasurementPackage::SensorType::RADAR;
    bool makeUpdate = isRadar ? useRadar : useLaser;

    if(makeUpdate) {
      double dt = (measurementPackage.timestamp - prevTimestamp) / 1000000.0;
      this->prediction(dt);
      this->update(measurementPackage, isRadar);
      prevTimestamp = measurementPackage.timestamp;
    }
  } else {
    prevTimestamp = measurementPackage.timestamp;
    if(measurementPackage.sensor_type == MeasurementPackage::SensorType::LASER) {
      double px = measurementPackage.raw_measurements[0];
      double py = measurementPackage.raw_measurements[1];
      x << px, py, 0.0, 0.0, 0.0;
    } else {
      double rho = measurementPackage.raw_measurements[0];
      double phi = measurementPackage.raw_measurements[1];
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
void UKF::prediction(double deltaT) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  MatrixXd XSigAug = this->augmentedSigmaPoints();
  this->sigmaPointPrediction(XSigAug, deltaT);
  this->predictMeanAndCovariance();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::update(MeasurementPackage measurementPackage, bool isRadar) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  std::tuple<VectorXd*, MatrixXd*, MatrixXd*> tup = this->predictMeasurement(isRadar);

  VectorXd *zPred = std::get<0>(tup);
  MatrixXd *S = std::get<1>(tup);
  MatrixXd *ZSigmaPoints = std::get<2>(tup);

  /*std::cout << *zPred << std::endl;
  std::cout << *S << std::endl;
  std::cout << *ZSigmaPoints << std::endl;*/
  this->updateState(*ZSigmaPoints, *zPred, *S, measurementPackage.raw_measurements, isRadar);

  delete(zPred);
  delete(S);
  delete(ZSigmaPoints);
}

MatrixXd UKF::augmentedSigmaPoints() {

  //create augmented mean vector
  VectorXd xAug = VectorXd(nAug);

  //create augmented state covariance
  MatrixXd PAug = MatrixXd(nAug, nAug);

  //create sigma point matrix
  MatrixXd XSigmaAug = MatrixXd(nAug, 2 * nAug + 1);

/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //create augmented mean state
  xAug.head(5) = x;
  xAug(5) = 0;
  xAug(6) = 0;

  //create augmented covariance matrix
  PAug.fill(0.0);
  PAug.topLeftCorner(5,5) = P;
  PAug(5,5) = stdAcc * stdAcc;
  PAug(6,6) = stdYawdd * stdYawdd;

  //create square root matrix
  MatrixXd L = PAug.llt().matrixL();

  //create augmented sigma points
  XSigmaAug.col(0)  = xAug;
  for (int i = 0; i< nAug; i++)
  {
    XSigmaAug.col(i + 1)       = xAug + sqrt(lambda + nAug) * L.col(i);
    XSigmaAug.col(i + 1 + nAug) = xAug - sqrt(lambda + nAug) * L.col(i);
  }

/*******************************************************************************
 * Student part end
 ******************************************************************************/

  //print result
  // std::cout << "XSigmaAug = " << std::endl << XSigmaAug << std::endl;

  //write result
  return XSigmaAug;
}

void UKF::sigmaPointPrediction(MatrixXd XSigmaAug, double deltaT) {

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
    double px = XSigmaAug(0,i);
    double py = XSigmaAug(1,i);
    double v = XSigmaAug(2,i);
    double yaw = XSigmaAug(3,i);
    double yawd = XSigmaAug(4,i);
    double nuAcc = XSigmaAug(5,i);
    double nuYawdd = XSigmaAug(6,i);

    //predicted state values
    double pxPred, pyPred;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
      pxPred = px + v/yawd * ( sin (yaw + yawd*deltaT) - sin(yaw));
      pyPred = py + v/yawd * ( cos(yaw) - cos(yaw+yawd*deltaT) );
    }
    else {
      pxPred = px + v*deltaT*cos(yaw);
      pyPred = py + v*deltaT*sin(yaw);
    }

    double vPred = v;
    double yawPred = yaw + yawd*deltaT;
    double yawdPred = yawd;

    //add noise
    pxPred = pxPred + 0.5*nuAcc*deltaT*deltaT * cos(yaw);
    pyPred = pyPred + 0.5*nuAcc*deltaT*deltaT * sin(yaw);
    vPred = vPred + nuAcc*deltaT;

    yawPred = yawPred + 0.5*nuYawdd*deltaT*deltaT;
    yawdPred = yawdPred + nuYawdd*deltaT;

    //write predicted sigma point into right column
    XSigPred(0,i) = pxPred;
    XSigPred(1,i) = pyPred;
    XSigPred(2,i) = vPred;
    XSigPred(3,i) = yawPred;
    XSigPred(4,i) = yawdPred;
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
    VectorXd xDiff = XSigPred.col(i) - predictedState;
    //angle normalization
    xDiff(3) = Tools::angleNormalization(xDiff(3));

    covMatPred = covMatPred + weights(i) * xDiff * xDiff.transpose() ;
  }

  x = predictedState;
  P = covMatPred;
/*******************************************************************************
 * Student part end
 ******************************************************************************/
}



std::tuple<VectorXd*, MatrixXd*, MatrixXd*> UKF::predictMeasurement(bool isRadar) {

  int n = isRadar ? nRadar : nLidar;

  //create matrix for sigma points in measurement space
  MatrixXd *ZSigmaPoints = new MatrixXd(n, 2 * nAug + 1);

  //mean predicted measurement
  VectorXd *zPred = new VectorXd(n);

  //measurement covariance matrix S
  MatrixXd *S = new MatrixXd(n, n);

/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //transform sigma points into measurement space
  if(isRadar) {
    for(int i = 0; i < 2 * nAug + 1; ++i) {
      VectorXd currentPred = VectorXd(n);
      double px = XSigPred(0, i);
      double py = XSigPred(1, i);
      double v = XSigPred(2, i);
      double yaw = XSigPred(3, i);
      yaw = Tools::angleNormalization(yaw);

      double sqrtVal = sqrt(px*px + py*py);
      double v1 = cos(yaw) * v;
      double v2 = sin(yaw) * v;

      currentPred << sqrtVal,
              atan2(py, px),
              (px*v1 + py*v2) / sqrtVal;

      (*ZSigmaPoints).col(i) = currentPred;
    }
  } else {
    for(int i = 0; i < 2 * nAug + 1; ++i) {
      VectorXd currentPred = VectorXd(n);
      double px = XSigPred(0, i);
      double py = XSigPred(1, i);

      currentPred << px, py; // for laser

      (*ZSigmaPoints).col(i) = currentPred;
    }
  }

  //calculate mean predicted measurement
  (*zPred).fill(0.0);
  for (int i=0; i < 2*nAug+1; i++) {
    *zPred = *zPred + weights(i) * (*ZSigmaPoints).col(i);
  }
  //calculate measurement covariance matrix S
  (*S).fill(0.0);
  for(int i = 0; i < 2 * nAug + 1; ++i) {
    VectorXd zDiff = (*ZSigmaPoints).col(i) - *zPred;

    if(isRadar) {
      zDiff(1) = Tools::angleNormalization(zDiff(1));
    }

    *S += weights(i) * zDiff * zDiff.transpose();
  }

  *S = *S + (isRadar ? RRadar : RLidar);

/*******************************************************************************
 * Student part end
 ******************************************************************************/

  return std::make_tuple(zPred, S, ZSigmaPoints);
}

void UKF::updateState(MatrixXd ZSigmaPoints, VectorXd zPred, MatrixXd S, VectorXd z, bool isRadar) {

  int n = isRadar ? nRadar : nLidar;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(nx, n);

/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * nAug + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd zDiff = ZSigmaPoints.col(i) - zPred;
    //angle normalization
    if(isRadar) {
      zDiff(1) = Tools::angleNormalization(zDiff(1));
    }

    // state difference
    VectorXd xDiff = XSigPred.col(i) - x;
    //angle normalization
    xDiff(3) = Tools::angleNormalization(xDiff(3));

    Tc = Tc + weights(i) * xDiff * zDiff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd predResidual = z - zPred;

  //angle normalization
  predResidual(1) = Tools::angleNormalization(predResidual(1));

  //update state mean and covariance matrix
  x = x + K * predResidual;
  P = P - K * S * K.transpose();

  double NIS = Tools::calculateNIS(predResidual, S);
  if(isRadar) {
    lastRadarNIS = NIS;
  } else {
    lastLidarNIS = NIS;
  }
/*******************************************************************************
 * Student part end
 ******************************************************************************/
}

const double &UKF::getX(int i) const {
  return x(i);
}


