#include "PID.h"

using namespace std;

/*
* TODO: Complete the PID class.
*/

PID::PID() {}

PID::~PID() {}

void PID::init(double Kp, double Ki, double Kd, double dp, double di, double dd) {
  this->coeff[0] = Kp;
  this->coeff[1] = Ki;
  this->coeff[2] = Kd;
  this->errors[0] = 0.0;
  this->errors[1] = 0.0;
  this->errors[2] = 0.0;
  this->prevCte = 0.0;

  this->currentIndex = 0;
  this->firstStep = true;
  this->d[0] = dp;
  this->d[1] = di;
  this->d[2] = dd;
  this->sumError = 0.0;
  this->bestError = 1e5;
}

void PID::updateError(double cte) {
  this->errors[0] = cte;
  this->errors[1] += cte;
  this->errors[2] = (cte - this->prevCte);

  this->sumError += fabs(cte * cte);

  this->prevCte = cte;
}

double PID::totalError() {
  double totalError = 0.0;
  for(int i = 0; i < 3; ++i) {
    totalError += (this->errors[i] * this->coeff[i]);
  }

  return -totalError;
}

void PID::twiddle() {
  sumError /= 20;
  if(firstStep) {
    coeff[currentIndex] += d[currentIndex];
    if(sumError < bestError) {
      bestError = sumError;
      d[currentIndex] *= 1.1;
      currentIndex = (currentIndex + 1) % 3;
    } else {
      firstStep = !firstStep;
    }
  } else {
    coeff[currentIndex] -= 2*d[currentIndex];
    if(sumError < bestError) {
      bestError = sumError;
      d[currentIndex] *= 1.1;
    } else {
      coeff[currentIndex] += d[currentIndex];
      d[currentIndex] *= 0.9;
    }

    currentIndex = (currentIndex + 1) % 3;
    firstStep = !firstStep;
  }

  cout << "best error: " << bestError << endl;
  cout << "P:" << coeff[0] << " I:" << coeff[1] << " D:" << coeff[2] << endl;
  cout << "current index: " << currentIndex << endl;
  sumError = 0.0;
}
