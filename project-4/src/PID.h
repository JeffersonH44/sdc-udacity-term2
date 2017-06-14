#ifndef PID_H
#define PID_H
#include <iostream>
#include <chrono>
#include <cmath>

using namespace std;

class PID {
public:
  /*
  * Errors
  */
  //P, I, D errors respectively
  double errors[3];

  /*
  * Coefficients
  */
  // P, I, D coefficients
  double coeff[3];

  chrono::high_resolution_clock::time_point startTime;
  chrono::high_resolution_clock::time_point prevTime;

  double prevCte;

  /*
   * params for twiddle
   */
  int currentIndex;
  bool firstStep;
  double d[3];
  double sumError;
  double bestError;
  /*
  * Constructor
  */
  PID();

  /*
  * Destructor.
  */
  virtual ~PID();

  /*
  * Initialize PID.
  */
  void init(double Kp, double Ki, double Kd, double dp, double di, double dd);

  /*
  * Update the PID error variables given cross track error.
  */
  void updateError(double cte);

  /*
  * Calculate the total PID error.
  */
  double totalError();

  void twiddle();
};

#endif /* PID_H */
