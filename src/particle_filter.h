/*
 * particle_filter.h
 *
 * 2D particle filter class.
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#ifndef PARTICLE_FILTER_H_
#define PARTICLE_FILTER_H_

#include <random>
#include "helper_functions.h"

struct Particle
{
	int id;
	double x;
	double y;
	double theta;
	double weight;
	std::vector<int> associations;
	std::vector<double> sense_x;
	std::vector<double> sense_y;
};

class ParticleFilter
{
public:

	// Constructor
	// @param nump Number of particles
	// @param sigp[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
	//   standard deviation of yaw [rad]]
	// @param sigl[] Array of dimension 2 [Landmark measurement uncertainty [x [m], y [m]]]
	ParticleFilter(int nump, double sigp[], double sigl[]);

	// Destructor
	~ParticleFilter() = default;

	/**
	 * init Initializes particle filter by initializing particles to Gaussian
	 *   distribution around first position and all the weights to 1.
	 * @param x Initial x position [m] (simulated estimate from GPS)
	 * @param y Initial y position [m]
	 * @param theta Initial orientation [rad]
	 *   standard deviation of yaw [rad]]
	 */
	void init(double x, double y, double theta);

	/**
	 * prediction Predicts the state for the next time step
	 *   using the process model.
	 * @param delta_t Time between time step t and t+1 in measurements [s]
	 * @param velocity Velocity of car from t to t+1 [m/s]
	 * @param yaw_rate Yaw rate of car from t to t+1 [rad/s]
	 */
	void prediction(double delta_t, double velocity, double yaw_rate);
	
	/**
	 * dataAssociation Finds which observations correspond to which landmarks (likely by using
	 *   a nearest-neighbors data association).
	 * @param predicted Vector of predicted landmark observations
	 * @param observations Vector of landmark observations
	 */
	//void dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations);
	void dataAssociation(const std::vector<LandmarkObs>& predicted, std::vector<LandmarkObs>& observations) const;

	/**
	 * updateWeights Updates the weights for each particle based on the likelihood of the 
	 *   observed measurements. 
	 * @param sensor_range Range [m] of sensor
	 * @param observations Vector of landmark observations
	 * @param map Map class containing map landmarks
	 */
	void updateWeights(double sensor_range, const std::vector<LandmarkObs>& observations, const Map& map_landmarks);
	
	/**
	 * resample Resamples from the updated set of particles to form
	 *   the new set of particles.
	 */
	void resample();

	/*
	 * Set a particles list of associations, along with the associations calculated world x,y coordinates
	 * This can be a very useful debugging tool to make sure transformations are correct and assocations correctly connected
	 */
	//Particle SetAssociations(Particle& particle, const std::vector<int>& associations,
	//	                     const std::vector<double>& sense_x, const std::vector<double>& sense_y);

	std::string getAssociations(const Particle& best) const;
	std::string getSenseX(const Particle& best) const;
	std::string getSenseY(const Particle& best) const;

	/**
	* initialized Returns whether particle filter is initialized yet or not.
	*/
	const bool initialized() const { return is_initialized_; }
	const std::vector<Particle>& getParticles() const { return particles_; }
	
private:

	int num_particles_;           // Number of particles to draw
	bool is_initialized_;         // Flag, if filter is initialized
	std::vector<double> weights_; // Vector of weights of all particles
	double sigma_pos_[3];         // GPS measurement uncertainty [x [m], y [m], theta [rad]]
        double sigma_landmark_[2];    // Landmark measurement uncertainty [x [m], y [m]]
	double gauss_norm_; // precalculating for caching purposes
	std::vector<Particle> particles_; // Set of current particles
};

#endif /* PARTICLE_FILTER_H_ */
