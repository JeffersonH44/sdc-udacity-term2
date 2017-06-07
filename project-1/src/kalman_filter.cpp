#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x, MatrixXd &P, MatrixXd &F,
                        MatrixXd &H, MatrixXd &R, MatrixXd &Q) {
  this->x = x;
  this->P = P;
  this->F = F;
  this->H = H;
  this->R = R;
  this->Q = Q;
}

void KalmanFilter::Predict() {
  x = F * x;
  MatrixXd Ft = F.transpose();
  P = F * P * Ft + Q;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd zPred = H * x;
  VectorXd y = z - zPred;
  MatrixXd Ht = H.transpose();
  MatrixXd S = H * P * Ht + R;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x = x + (K * y);
  long xSize = x.size();
  MatrixXd I = MatrixXd::Identity(xSize, xSize);
  P = (I - K * H) * P;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  VectorXd y = z - this->transformRadarPred(x);
  while (y(1) > M_PI) {
    y(1) -= 2 * M_PI;
  }
  while (y(1) < -M_PI) {
    y(1) += 2 * M_PI;
  }

  MatrixXd Ht = H.transpose();
  MatrixXd S = H * P * Ht + R;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x = x + (K * y);
  long xSize = x.size();
  MatrixXd I = MatrixXd::Identity(xSize, xSize);
  P = (I - K * H) * P;
}

Eigen::VectorXd KalmanFilter::transformRadarPred(const Eigen::VectorXd &x) {
  VectorXd hPred(3);
  double px = x(0);
  double py = x(1);
  double vx = x(2);
  double vy = x(3);
  hPred(0) = sqrt(px*px + py*py);

  if(fabs(px) > 0.0001){
    hPred(1) = atan2(py, px);
  }

  if(fabs(hPred(0)) < 0.01) {
    hPred(2) = 0.0;
  } else {
    hPred(2) = (px*vx + py*vy) / hPred(0);
  }

  return hPred;
}
