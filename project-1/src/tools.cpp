#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::calculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
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

MatrixXd Tools::calculateJacobian(const VectorXd &x_state) {
  MatrixXd Hj(3,4);
  //recover state parameters
  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);

  //TODO: YOUR CODE HERE

  //check division by zero
  if(px < 0.01 && py < 0.01) {
      std::cout << "Division by zero" << std::endl;
  }

  double px_2 = px * px;
  double py_2 = py * py;

  //compute the Jacobian matrix
  double one = px / sqrt(px_2 + py_2);
  double two = py / sqrt(px_2 + py_2);
  double five = -(py / (px_2 + py_2));
  double six = px / (px_2 + py_2);
  double nine = (py*(vx*py - vy*px)) / pow(px_2 + py_2, 1.5);
  double ten = (px*(vy*px - vx*py)) / pow(px_2 + py_2, 1.5);
  double eleven = one;
  double twelve = two;

  Hj << one, two, 0, 0,
          five, six, 0, 0,
          nine, ten, eleven, twelve;

  return Hj;
}
