[p1]: ./img/p-controller.png "Relation p"
[p2]: ./img/p-controller1.png "Relation p"
[p3]: ./img/p-controller2.png "Relation p"
[d1]: ./img/d-controller.png "Relation d"
[d2]: ./img/d-controller1.png "Relation d"
[i1]: ./img/i-controller.png "Relation i"
[i2]: ./img/i-controller1.png "Relation i"

# Reflection over PID controller

## P - controller

This controller just take the Cross track error(cte) and multiply by a scale 
factor, this makes the steering angle larger or shorter given the cross track 
error:

![alt text][p1]

Higher p gain is better because it will converge to the track line faster, but
if the car is too far away from the track line, it will make run in circles like
the image below.

![alt text][p2]

Also. if there is a dramatic change on the track line it will overshoots the 
trajectory like the image below.

![alt text][p3]

## D - controller

This controller just add a proportional control over the cross track error rate
(derivative term) or how fast the car is moving on the perpendicular direction,
like the image below:

![alt text][d1]

If the d parameter is too low it will underdamp the track line, and if it's too
high, it will take a long time to arrive to the track line, a good gain
 will make the car arrive to the track line faster, like the image below:

![alt text][d2]

## I - controller

This controller takes a sum on line offsets to put the car again in the desired
trajectory, this is used when the car hit some obstacles on the road and makes 
the car goes out of the desired trajectory, a lane offset is like the image below:

![alt text][i2]

A higher gain will make the car unstable around the trajectory because it will
exaggerate the real line offset, a lower gain will underdamp the car over the 
real trajectory.

![alt text][i1]

## Hyperparameters tunning

For hyperparameter tunning I used twiddle using the initial parameters from the
class (p=0.2, d=3.0, i=0.0004), finally I got these parameters that works best
for me (p=0.3, d=2.0, i=0.0004), for twiddle I make 20 iterations from the 
simulator and then I applied it.

I tried manual parameter tunning but is not easy to find a good parameters,
so it's more easy to use twiddle.