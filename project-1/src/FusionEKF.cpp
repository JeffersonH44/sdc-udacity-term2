#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  isInitialized = false;

  previousTimestamp = 0;

  // initializing matrices
  RLaser = MatrixXd(2, 2);
  RRadar = MatrixXd(3, 3);
  HLaser = MatrixXd(2, 4);
  Hj = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  RLaser << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  RRadar << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  HLaser << 1, 0, 0, 0,
          0, 1, 0, 0;

  Hj << 1, 1, 0, 0,
          1, 1, 0, 0,
          1, 1, 1, 1;

  // Initialize ekf variables from part 13
  ekf.F = MatrixXd(4, 4);
  ekf.F << 1, 0, 1/*dt*/, 0,
             0, 1, 0, 1/*dt*/,
             0, 0, 1, 0,
             0, 0, 0, 1;

  ekf.P = MatrixXd(4, 4);
  ekf.P << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1000, 0,
             0, 0, 0, 1000;


}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!isInitialized) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf.x = VectorXd(4);
    ekf.x << 1, 1, 1, 1;

    if (measurement_pack.sensorType == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double rho = measurement_pack.rawMeasurements(0);
      double phi = measurement_pack.rawMeasurements(1);
      double rho_dot = measurement_pack.rawMeasurements(2);
      ekf.x(0) = rho * cos(phi);
      ekf.x(1) = rho * sin(phi);
      ekf.x(2) = rho_dot * cos(phi);
      ekf.x(3) = rho_dot * sin(phi);
    }
    else if (measurement_pack.sensorType == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf.x(0) = measurement_pack.rawMeasurements(0);
      ekf.x(1) = measurement_pack.rawMeasurements(1);
    }

    previousTimestamp = measurement_pack.timestamp;
    // done initializing, no need to predict or update
    isInitialized = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
     * using code from lesson 5:13
   */
  // update F matrix
  double dt = (measurement_pack.timestamp - previousTimestamp) / 1000000.0;
  previousTimestamp = measurement_pack.timestamp;
  ekf.F(0, 2) = dt;
  ekf.F(1, 3) = dt;

  // update Q matrix
  double noise_ax = 9, noise_ay = 9;
  double dt_2 = dt * dt;
  double dt_3 = dt_2 * dt;
  double dt_4 = dt_3 * dt;
  ekf.Q = MatrixXd(4, 4);
  ekf.Q <<  (dt_4 / 4) * noise_ax, 0, (dt_3 / 2) * noise_ax, 0,
                0, (dt_4 / 4) * noise_ay, 0, (dt_3 / 2) * noise_ay,
                (dt_3 / 2) * noise_ax, 0, dt_2 * noise_ax, 0,
                0, (dt_3 / 2) * noise_ay, 0, dt_2 * noise_ay;


  ekf.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensorType == MeasurementPackage::RADAR) {
    // Radar updates
    Hj = Tools::calculateJacobian(ekf.x);
    ekf.H = Hj;
    ekf.R = RRadar;
    ekf.UpdateEKF(measurement_pack.rawMeasurements);
  } else {
    // Laser updates
    ekf.H = HLaser;
    ekf.R = RLaser;
    ekf.Update(measurement_pack.rawMeasurements);
  }

  // print the output
  cout << "x_ = " << ekf.x << endl;
  cout << "P_ = " << ekf.P << endl;
}
