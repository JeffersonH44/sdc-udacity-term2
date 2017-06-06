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

double Tools::calculateRadNIS(Eigen::VectorXd z, Eigen::VectorXd z_pred, Eigen::MatrixXd S) {
  VectorXd Zdiff = z - z_pred;

  while(Zdiff(1) > M_PI) Zdiff(1) -= 2 * M_PI;
  while(Zdiff(1) < -M_PI) Zdiff(1) += 2 * M_PI;

  return Zdiff.transpose() * S.inverse() * Zdiff;
}

