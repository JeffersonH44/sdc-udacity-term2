# MPC Controller

## The model

### State vector

For the model we have 6 states for the car:

* x coordinate from map coordinates.
* y coordinate from map coordinates.
* velocity on mph.
* psi in radians for the vehicle orientation.
* cross track error of the vehicle from the desired trajectory.
* error psi degree in radians between trajectory orientation and
the vehicle orientation.

### Actuators

There are two actuators for the model

* the steering angle in radians
* the throttle value that allow to break and accelerate

### Update equations

I used the kinematic equations used in the class that are 
the following:

```
x1 = x0 + v0 * cos(psi0) * dt
y1 = y0 + v0 * sin(psi0) * dt
psi1 = psi0 - v0 * delta0 * dt / Lf
v1 = v0 + a0 * dt;
cte1 = (f0 - y0 + (v0 * sin(epsi0) * dt))
epsi1 = (psi0 - psiDest - (v0 * delta0 * dt / Lf))
```

## Timestep Length and Elapsed Duration (N & dt)

at the beginning I started using the value used on the class
(N=25, dt=0.05), but this was too large and the model was not working
well, so I reduce N=10 and dt=0.2 but the car steering angle doesn't 
stabilize, so to make the model more confident I put N=15, with that
the car start to follow the desired trajectory.

also that dt=0.2 seems to help to deal with latency because I didn't
need add any latency to the model.

## Polynomial Fitting and MPC Preprocessing

The polynomial fitting degree was 3 because it works on most of the road
as said in the class, and for the preprocessing step I transform the 
trajectory and car points from map coordinates to car coordinates:

for the trajectory points I used this:

```
Eigen::VectorXd pointsX(ptsx.size());
Eigen::VectorXd pointsY(ptsy.size());
for(int i = 0; i < ptsx.size(); ++i) {
    double xDiff = ptsx[i] - x;
    double yDiff = ptsy[i] - y;
    
    pointsX[i] = (xDiff * cos(-psi) - yDiff * sin(-psi));
    pointsY[i] = (xDiff * sin(-psi) + yDiff * cos(-psi));
}
``` 

where x and y are the car position in map coordinates and the transformation used was:

```
xDiff = point_trajectory - x_map_car;
yDiff = point_trajectory - y_map_car;
// remember that we go from map orientation to car orientation
// and I assume that map orientation is 0
orientation = mapOrientation - carOrientation

// coordinates from car coordinates
x_car = xDiff * cos(orientation) - yDiff * sin(orientation)
y_car = xDiff * sin(orientation) + yDiff * cos(orientation)
```

then we create a polynomial with the new points from car coordinates to get
the cross track error(cte) and error psi in the following way:

```
auto coeffs = polyfit(pointsX, pointsY, 3);
double cte = polyeval(coeffs, 0);
double epsi = -atan(coeffs[1]);
```

the formula for calculate the cross track error and epsi are the following:
```
cte = y - polynomial(coeffs, 0) // we evaluate x at position 0 because now we are in car coordinates
epsi = -atan(x0*(3*coeffs[3]*x0 + 2*coeffs[2]) + coeffs[1]) // using the derivative inside
// we can simplify the calculation because x = 0
epsi = -atan(coeffs[1])
```

and then we can create the state vector as follows:
```
Eigen::VectorXd state(6);
state << 0, 0, 0, v, cte, epsi;
```

x, y, psi are zero because we are in car coordinates,
the velocity is the same, cte and epsi are the values 
that we calculate before

## Model Predictive Control with Latency

This process was not necessary because the environment doesn't 
change too much between the 100ms and I think that the 'dt' 
parameter of the MPC model for prediction help a lot with that delay,

Also I think that this delay should be bigger because you have 
to consider all the computations that make a self-driving car.

## Parameters for cost function

The way that I calculate the cost function is the following
```
const double paramCte = 3000;
const double paramEpsi = 1000;
const double paramVel = 1.0;
const double paramDelta = 2000.0;
const double paramAcceleration = 1.0;
const double paramDeltaDiff = 5000;
const double paramAccDiff = 10000;

// cost over cte and direction, and velocity should be close to the reference value
for (int t = 0; t < N; ++t) {
  fg[0] += paramCte * CppAD::pow(vars[cte_start + t], 2);
  fg[0] += paramEpsi * CppAD::pow(vars[epsi_start + t], 2);
  fg[0] += paramVel * CppAD::pow(vars[v_start + t] - ref_v, 2);
}

// actuators cost function
for (int t = 0; t < N - 1; ++t) {
  fg[0] += paramDelta * CppAD::pow(vars[delta_start + t], 2);
  fg[0] += paramAcceleration * CppAD::pow(vars[a_start + t], 2);
}

// actuators should act smoothly between timesteps so we take in count the current and the previous timestep
for (int t = 0; t < N - 2; ++t) {
  fg[0] += paramDeltaDiff * CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t], 2);
  fg[0] += paramAccDiff * CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2);
}
```

I used a very high value for 'paramDeltaDiff' and 'paramDeltaDiff'
because a very drastic change of this values make car unstable on 
the road so I penalize these values more than the other ones. I 
don't penalize too much the epsi value to make the car arrive 
smoothly to the desired trajectory.

The velocity wasn't penalized too much otherwise the car
doens't get the 50 mph that is the reference velocity, also 
the acceleration was not necessary because I penalize the 
acceleration between two timesteps.

