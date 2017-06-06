#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  // ... your code here
  if(estimations.size() == 0 || (estimations.size() != ground_truth.size())) {
    std::cout << "estimations size invalid" << std::endl;
    return rmse;
  }

  //accumulate squared residuals
  for(int i=0; i < estimations.size(); ++i){
    // ... your code here
    VectorXd res = estimations[i] - ground_truth[i];
    res = res.array() * res.array();
    rmse += res;
  }

  //calculate the mean
  // ... your code here
  rmse = rmse.array() / estimations.size();

  //calculate the squared root
  // ... your code here
  rmse = rmse.array().sqrt();

  //return the result
  return rmse;
}

double Tools::calculateRadNIS(Eigen::VectorXd z_diff, Eigen::MatrixXd S) {
  return z_diff.transpose() * S.inverse() * z_diff;
}

double Tools::angleNormalization(double angle) {
  if (angle > M_PI)
    angle = fmod(angle - M_PI, 2*M_PI) - M_PI;
  if (angle < -M_PI)
    angle = fmod(angle + M_PI, 2*M_PI) + M_PI;
  return angle;
}

