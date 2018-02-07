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
#include <limits>
#include <cstring>

#include "particle_filter.h"

using namespace std;

ParticleFilter::ParticleFilter(int nump, double sigp[], double sigl[])
	: num_particles_(nump), is_initialized_(false)
{
	memcpy(sigma_pos_, sigp, sizeof(sigma_pos_));
	memcpy(sigma_landmark_, sigl, sizeof(sigma_landmark_));

	gauss_norm_ = 1 / (2 * M_PI * sigma_landmark_[0] * sigma_landmark_[1]);
}

void ParticleFilter::init(double x, double y, double theta)
{
	// Setting the number of particles. Initialize all particles to first position
	// (based on estimates of x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Adding random Gaussian noise to each particle.

	normal_distribution<double> dist_x(x,     sigma_pos_[0]);
	normal_distribution<double> dist_y(y,     sigma_pos_[1]);
	normal_distribution<double> dist_t(theta, sigma_pos_[2]);

	default_random_engine rndgen;

	int id = 0;
	particles_.resize(num_particles_); // Set of current particles
	for (auto& p : particles_)
	{
		p.id = ++id; // 1, 2, ..., num_particles_
		p.x = dist_x(rndgen);
		p.y = dist_y(rndgen);
		p.theta = dist_t(rndgen);
		p.weight = 1.0f;
		vector<int>().swap(p.associations); // empty by default
		vector<double>().swap(p.sense_x);
		vector<double>().swap(p.sense_y);
	}

	is_initialized_ = true;
}

void ParticleFilter::prediction(double delta_t, double velocity, double yaw_rate)
{
	// Adding measurements to each particle and add random Gaussian noise.

	default_random_engine rndgen;

	for (auto& p : particles_)
	{
		if (fabs(yaw_rate) > EPS)
		{
			p.x += velocity / yaw_rate * (sin(p.theta + yaw_rate * delta_t) - sin(p.theta));
			p.y += velocity / yaw_rate * (cos(p.theta) - cos(p.theta + yaw_rate * delta_t));
			p.theta += yaw_rate * delta_t;
		}
		else // zero yaw_rate
		{
			p.x += velocity * delta_t * cos(p.theta);
			p.y += velocity * delta_t * sin(p.theta);
		}

		normal_distribution<double> dist_x(p.x, sigma_pos_[0]);
    		p.x = dist_x(rndgen);
		normal_distribution<double> dist_y(p.y, sigma_pos_[1]);
    		p.y = dist_y(rndgen);
		normal_distribution<double> dist_t(p.theta, sigma_pos_[2]);
    		p.theta = dist_t(rndgen);

		normalize_angle(p.theta);
	}
}

void ParticleFilter::dataAssociation(const vector<LandmarkObs>& predicted, vector<LandmarkObs>& observations) const
{
	// Finding the predicted measurement that is closest to each observed measurement and assign the 
	// observed measurement to this particular landmark.
	// This method will NOT be called by the grading code. But you will probably find it useful to 
	// implement this method and use it as a helper during the updateWeights phase.

	double dmin;
	double dtmp;
	for (auto& lmo : observations)
	{
		dmin = numeric_limits<double>::max();
		for (const auto& lmp : predicted)
		{
			dtmp = dist(lmo.x, lmo.y, lmp.x, lmp.y);
			if (dtmp < dmin)
			{
				dmin = dtmp;
				lmo.id = lmp.id;
			}
		}
	}
}

void ParticleFilter::updateWeights(double sensor_range, const vector<LandmarkObs>& observations, const Map& map_landmarks)
{
	// Updating the weights of each particle using a mult-variate Gaussian distribution. You can read
	// more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	// according to the MAP'S coordinate system. You will need to transform between the two systems.
	// Keep in mind that this transformation requires both rotation AND translation (but no scaling).

	//Transformed Observation (x_map,y_map) = func(x_particle, y_particle, heading_particle, x_obs, y_obs)

	vector<LandmarkObs> mapObs; // observations in map coordinates
	mapObs.resize(observations.size());

	vector<double>().swap(weights_);

	for (auto& p : particles_)
	{
		for (short n = 0; n < observations.size(); ++n)
		{
			auto& mo = mapObs[n];
			mo.id = observations[n].id;
			mo.x = p.x + cos(p.theta) * observations[n].x - sin(p.theta) * observations[n].y;
			mo.y = p.y + sin(p.theta) * observations[n].x + cos(p.theta) * observations[n].y;
		}

		vector<LandmarkObs> mapLandm; // filtered landmarks (the ones in sensor range)
		for (const auto& ml : map_landmarks.landmark_list)
			if (dist(ml.x_f, ml.y_f, p.x, p.y) <= sensor_range)
				mapLandm.push_back({ ml.id_i, ml.x_f, ml.y_f });


		dataAssociation(mapLandm, mapObs);
		p.weight = 1.0; // resetting weight
		vector<int>().swap(p.associations);
		vector<double>().swap(p.sense_x);
		vector<double>().swap(p.sense_y);

		for (const auto& mo : mapObs)
		{
			double mlx, mly;
			for (const auto& ml : mapLandm) if (ml.id == mo.id) { mlx = ml.x; mly = ml.y; }

			p.weight *= gauss_norm_ * exp(-0.5 * pow(mo.x - mlx, 2.0) / pow(sigma_landmark_[0], 2.0) +
				-0.5 * pow(mo.y - mly, 2.0) / pow(sigma_landmark_[1], 2.0));
			p.associations.push_back(mo.id);
			p.sense_x.push_back(mo.x);
			p.sense_y.push_back(mo.y);
		}
		weights_.push_back(p.weight);
	}
}

void ParticleFilter::resample()
{
	// Resampling particles with replacement with probability proportional to their weight. 

	vector<Particle> particles_new(num_particles_);

	random_device rd;
	mt19937 rndgen(rd());

	discrete_distribution<> dstr(weights_.begin(), weights_.end());

	for (int i = 0; i < num_particles_; ++i)
		particles_new[i] = particles_[dstr(rndgen)];

	particles_.swap(particles_new);
}

//Particle ParticleFilter::SetAssociations(Particle& particle, const vector<int>& associations, 
//                                     const vector<double>& sense_x, const vector<double>& sense_y)
//{
//	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
//	// associations: The landmark id that goes along with each listed association
//	// sense_x: the associations x mapping already converted to world coordinates
//	// sense_y: the associations y mapping already converted to world coordinates
//
//	particle.associations= associations;
//	particle.sense_x = sense_x;
//	particle.sense_y = sense_y;
//}

string ParticleFilter::getAssociations(const Particle& best) const
{
	const vector<int>& v = best.associations;
	stringstream ss;
	copy(v.begin(), v.end(), ostream_iterator<int>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length()-1);  // get rid of the trailing space
	return s;
}

string ParticleFilter::getSenseX(const Particle& best) const
{
	const vector<double>& v = best.sense_x;
	stringstream ss;
	copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length()-1);  // get rid of the trailing space
	return s;
}

string ParticleFilter::getSenseY(const Particle& best) const
{
	const vector<double>& v = best.sense_y;
	stringstream ss;
	copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length()-1);  // get rid of the trailing space
	return s;
}
