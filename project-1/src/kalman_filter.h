#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_
#include "Eigen/Dense"

class KalmanFilter {
public:

  // state vector
  Eigen::VectorXd x;

  // state covariance matrix
  Eigen::MatrixXd P;

  // state transistion matrix
  Eigen::MatrixXd F;

  // process covariance matrix
  Eigen::MatrixXd Q;

  // measurement matrix
  Eigen::MatrixXd H;

  // measurement covariance matrix
  Eigen::MatrixXd R;

  /**
   * Constructor
   */
  KalmanFilter();

  /**
   * Destructor
   */
  virtual ~KalmanFilter();

  /**
   * Init Initializes Kalman filter
   * @param x Initial state
   * @param P Initial state covariance
   * @param F Transition matrix
   * @param H Measurement matrix
   * @param R Measurement covariance matrix
   * @param Q Process covariance matrix
   */
  void Init(Eigen::VectorXd &x, Eigen::MatrixXd &P, Eigen::MatrixXd &F,
      Eigen::MatrixXd &H, Eigen::MatrixXd &R, Eigen::MatrixXd &Q);

  /**
   * Prediction Predicts the state and the state covariance
   * using the process model
   * @param delta_T Time between k and k+1 in s
   */
  void Predict();

  /**
   * Updates the state by using standard Kalman Filter equations
   * @param z The measurement at k+1
   */
  void Update(const Eigen::VectorXd &z);

  /**
   * Updates the state by using Extended Kalman Filter equations
   * @param z The measurement at k+1
   */
  void UpdateEKF(const Eigen::VectorXd &z);

  /**
   * Calculate the h function for the radar measurement.
   * @param z The measurement at k+1
   */
  Eigen::VectorXd transformRadarPred(const Eigen::VectorXd &x);
};

#endif /* KALMAN_FILTER_H_ */
