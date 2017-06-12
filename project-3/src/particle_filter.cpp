/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[], Map map_landmarks) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	this->num_particles = 1250;
	this->particles.resize(this->num_particles);

  normal_distribution<double> x_gen(x, std[0]);
  normal_distribution<double> y_gen(y, std[1]);
  normal_distribution<double> theta_gen(theta, std[2]);

  for(int i = 0; i < this->num_particles; ++i) {
    Particle p;
    p.id = i;
    p.x = x_gen(eng);
    p.y = y_gen(eng);
    p.theta = theta_gen(eng);
    p.weight = 1;

    p.associations.clear();
    p.sense_x.clear();
    p.sense_y.clear();

    this->particles[i] = p;
  }

  this->is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
  normal_distribution<double> stds[3];
  stds[0] = normal_distribution<double>(0, std_pos[0]);
  stds[1] = normal_distribution<double>(0, std_pos[1]);
  stds[2] = normal_distribution<double>(0, std_pos[2]);

  for(int i = 0; i < this->num_particles; ++i) {
    this->particles[i] = move(this->particles[i], delta_t, stds, velocity, yaw_rate);
  }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

  double stdX = std_landmark[0];
  double stdY = std_landmark[1];
  auto mapLandmarks = map_landmarks.landmark_list;

  for(int i = 0; i < this->num_particles; ++i) {
    double particleX = this->particles[i].x;
    double particleY = this->particles[i].y;
    double particleTheta = this->particles[i].theta;

    vector<LandmarkObs> predictionsPart(observations.size());

    for (int j = 0; j < observations.size(); ++j) {
      auto currentObs = observations[j];
      double x = (cos(particleTheta) * currentObs.getX()) - (sin(particleTheta) * currentObs.getY()) + particleX;
      double y = (sin(particleTheta) * currentObs.getX()) + (cos(particleTheta) * currentObs.getY()) + particleY;

      predictionsPart[j] = LandmarkObs(currentObs.id, x, y);
    }

    double weight = 1.0;

    for (int j = 0; j < predictionsPart.size(); ++j) {
      double predX = predictionsPart[j].getX();
      double predY = predictionsPart[j].getY();
      double closestMarkX = 0.0, closestMarkY = 0.0, minDist = numeric_limits<double>::max();
      for (int k = 0; k < mapLandmarks.size(); ++k) {
        double landX = mapLandmarks[k].x_f;
        double landY = mapLandmarks[k].y_f;

        double currentDist = dist(predX, predY, landX, landY);
        if(currentDist < minDist) {
          minDist = currentDist;
          closestMarkX = landX;
          closestMarkY = landY;
        }
      }

      double currentWeight = multivariateGaussianProb(predX, predY, closestMarkX, closestMarkY, stdX, stdY);
      if(currentWeight > 0.0) {
        weight *= currentWeight;
      }
    }

    this->particles[i].weight = weight;
  }
}

struct myclass {
    bool operator() (Particle i, Particle j) { return i.weight < j.weight; }
} myobj;

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  // using resample wheel there's is no need of normalization

  vector<Particle> resampledParticles(this->particles.size());
  double maxWeight = (*max_element(this->particles.begin(), this->particles.end(), myobj)).weight;

  int index = 0;

  uniform_real_distribution<double> getReal(0, 2.0 * maxWeight);

  double beta = 0.0;
  for(int i = 0; i < this->num_particles; ++i) {
    beta += getReal(eng);
    while(beta > this->particles[index].weight) {
      beta -= this->particles[index].weight;
      index = (index + 1) % this->num_particles;
    }

    resampledParticles[i] = this->particles[index];
  }

  particles = resampledParticles;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
  copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

Particle ParticleFilter::move(Particle p, double delta_t, normal_distribution<double> stds[], double velocity, double yaw_rate) {
  if(yaw_rate < 1e-9) {
    p.x += velocity * cos(p.theta) * delta_t;
    p.y += velocity * sin(p.theta) * delta_t;
  } else {
    p.x += (velocity / yaw_rate) * (sin(p.theta + yaw_rate * delta_t) - sin(p.theta));
    p.y += (velocity / yaw_rate) * (cos(p.theta) - cos(p.theta + yaw_rate * delta_t));
    p.theta += (yaw_rate * delta_t);
  }

  p.x += stds[0](eng);
  p.y += stds[1](eng);
  p.theta += stds[2](eng);

  return p;
}
