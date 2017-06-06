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
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // predictions variables to make
  n_x_ = 5;
  // augmented by the two random variables
  n_aug_ = 7;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.7;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  n_radar = 3;

  n_lidar = 2;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  // spreading parameter
  lambda_ = 3 - n_aug_;

  // weights for the sigma points
  weights_ = VectorXd(2 * n_aug_ +1);
  weights_.fill(1.0 / (2.0*(lambda_ + n_aug_)));
  weights_(0) = lambda_ / (lambda_ + n_aug_);

  last_lidar_NIS = 0.0;
  last_radar_NIS = 0.0;

  R_laser = MatrixXd(n_lidar, n_lidar);
  R_laser << std_laspx_ * std_laspx_, 0,
             0, std_laspy_ * std_laspy_;

  R_radar = MatrixXd(n_radar, n_radar);
  R_radar << std_radr_ * std_radr_, 0, 0,
             0, std_radphi_ * std_radphi_, 0,
             0, 0, std_radrd_ * std_radrd_;

  is_initialized_ = false;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  if(is_initialized_) {
    double dt = (meas_package.timestamp_ - prev_timestamp_) / 1000000.0;
    std::cout << "dt: " << dt << std::endl;
    if(meas_package.sensor_type_ == MeasurementPackage::SensorType::LASER) {
      // do nothing right now
    } else {
      this->Prediction(dt);
      this->UpdateRadar(meas_package);
    }
    prev_timestamp_ = meas_package.timestamp_;
  } else {

    prev_timestamp_ = meas_package.timestamp_;
    if(meas_package.sensor_type_ == MeasurementPackage::SensorType::LASER) {
      double px = meas_package.raw_measurements_[0];
      double py = meas_package.raw_measurements_[1];
      x_ << px, py, 0.0, 0.0, 0.0;
    } else {
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double rho_dot = meas_package.raw_measurements_[2];
      x_ << rho * cos(phi), rho * sin(phi), 0.0, 0.0, 0.0;
    }

    // TODO: do not forget to change those values
    P_ << 50, 0, 0, 0, 0,
          0, 50, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 0.5, 0,
          0, 0, 0, 0, 0.5;

    is_initialized_ = true;
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  MatrixXd Xsig_aug = this->AugmentedSigmaPoints();
  this->SigmaPointPrediction(Xsig_aug, delta_t);
  this->PredictMeanAndCovariance();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  std::tuple<VectorXd*, MatrixXd*, MatrixXd*> tup = this->PredictRadarMeasurement();

  VectorXd *z_pred = std::get<0>(tup);
  MatrixXd *S = std::get<1>(tup);
  MatrixXd *Zsig = std::get<2>(tup);

  /*std::cout << *z_pred << std::endl;
  std::cout << *S << std::endl;
  std::cout << *Zsig << std::endl;*/
  this->UpdateRadarState(*Zsig, *z_pred, *S, meas_package.raw_measurements_);

  delete(z_pred);
  delete(S);
  delete(Zsig);
}

MatrixXd UKF::AugmentedSigmaPoints() {

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i + 1)       = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }

/*******************************************************************************
 * Student part end
 ******************************************************************************/

  //print result
  // std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;

  //write result
  return Xsig_aug;
}

void UKF::SigmaPointPrediction(MatrixXd Xsig_aug, double delta_t) {

  //create matrix with predicted sigma points as columns
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //predict sigma points
  //avoid division by zero
  //write predicted sigma points into right column
  for (int i = 0; i< 2 * n_aug_ + 1; i++)
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
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }


/*******************************************************************************
 * Student part end
 ******************************************************************************/
}

void UKF::PredictMeanAndCovariance() {
/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //create vector for predicted state
  VectorXd x = VectorXd(n_x_);

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd::Zero(n_x_, n_x_);
  //predict state mean

  x.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x = x+ weights_(i) * Xsig_pred_.col(i);
  }
  //predict state covariance matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x;
    //angle normalization
    x_diff(3) = Tools::angleNormalization(x_diff(3));

    P = P + weights_(i) * x_diff * x_diff.transpose() ;
  }

  x_ = x;
  P_ = P;
/*******************************************************************************
 * Student part end
 ******************************************************************************/
}



std::tuple<VectorXd*, MatrixXd*, MatrixXd*> UKF::PredictRadarMeasurement() {

  //create matrix for sigma points in measurement space
  MatrixXd *Zsig = new MatrixXd(n_radar, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd *z_pred = new VectorXd(n_radar);

  //measurement covariance matrix S
  MatrixXd *S = new MatrixXd(n_radar, n_radar);

/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //transform sigma points into measurement space
  for(int i = 0; i < 2 * n_aug_ + 1; ++i) {
    VectorXd currentPred = VectorXd(n_radar);
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);
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
  for (int i=0; i < 2*n_aug_+1; i++) {
    *z_pred = *z_pred + weights_(i) * (*Zsig).col(i);
  }
  //calculate measurement covariance matrix S
  (*S).fill(0.0);
  for(int i = 0; i < 2 * n_aug_ + 1; ++i) {
    VectorXd Zdiff = (*Zsig).col(i) - *z_pred;

    Zdiff(1) = Tools::angleNormalization(Zdiff(1));

    *S += weights_(i) * Zdiff * Zdiff.transpose();
  }

  *S = *S + R_radar;

/*******************************************************************************
 * Student part end
 ******************************************************************************/

  return std::make_tuple(z_pred, S, Zsig);
}

void UKF::UpdateRadarState(MatrixXd Zsig, VectorXd z_pred, MatrixXd S, VectorXd z) {

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_radar);

/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    z_diff(1) = Tools::angleNormalization(z_diff(1));

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    x_diff(3) = Tools::angleNormalization(x_diff(3));

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  z_diff(1) = Tools::angleNormalization(z_diff(1));

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  last_radar_NIS = Tools::calculateRadNIS(z_diff, S);
/*******************************************************************************
 * Student part end
 ******************************************************************************/
}
