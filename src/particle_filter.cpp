/**
* particle_filter.cpp
*
* Created on: Dec 12, 2016
* Author: Tiffany Huang
*/

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;

std::default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	/**
	* TODO: Set the number of particles. Initialize all particles to
	*   first position (based on estimates of x, y, theta and their uncertainties
	*   from GPS) and all weights to 1.
	* TODO: Add random Gaussian noise to each particle.
	* NOTE: Consult particle_filter.h for more information about this method
	*   (and others in this file).
	* TODO:设置粒子数。将所有粒子初始化到第一个位置（基于对x，y，theta的估计及其来自GPS的不确定性），并将所有权重设置为1。
	* TODO:向每个粒子添加随机高斯噪声。
	* 注意：有关此方法的详细信息，请参阅particle_filter.h（以及此文件中的其他内容）。
	*/


	num_particles = 100;  // TODO: Set the number of particles

	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);
	
	for (int i = 0; i < num_particles; i++){
		Particle particle;
		particle.id = i;
		particle.x = dist_x(gen);
		particle.y = dist_y(gen);
		particle.theta = dist_theta(gen);
		particle.weight = 1.0;

		particles.push_back(particle);
	}

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[],
	double velocity, double yaw_rate) {
	/**
	* TODO: Add measurements to each particle and add random Gaussian noise.
	* TODO: 代入不同的时间步，向每个粒子添加测量值并添加随机高斯噪声。
	* NOTE: When adding noise you may find std::normal_distribution
	*   and std::default_random_engine useful.
	*  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	*  http://www.cplusplus.com/reference/random/default_random_engine/
	*/
	normal_distribution<double> dist_x(0,std_pos[0]);
	normal_distribution<double> dist_y(0,std_pos[1]);
	normal_distribution<double> dist_theta(0,std_pos[2]);

	for (int i = 0; i < num_particles; i++){
		if (fabs(yaw_rate)<0.0001){
			yaw_rate = 0.0001;
		}
		particles[i].x = particles[i].x + velocity*(sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta)) / yaw_rate + dist_x(gen);
		particles[i].y = particles[i].y + velocity*(cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t)) / yaw_rate + dist_y(gen);
		particles[i].theta = particles[i].theta + yaw_rate*delta_t + dist_theta(gen);
	}
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted,
	vector<LandmarkObs>& observations) {
	/**
	* TODO: Find the predicted measurement that is closest to each
	*   observed measurement and assign the observed measurement to this
	*   particular landmark.
	* TODO: 找到与每个观测到的测量值最接近的预测测量值，并将观测到的测量值分配给该特定地标。
	* NOTE: this method will NOT be called by the grading code. But you will
	*   probably find it useful to implement this method and use it as a helper
	*   during the updateWeights phase.
	*/
	for (int j = 0; j < observations.size(); j++){
		double distance;
		double min_distance;
		for (int k = 0; k < predicted.size(); k++){
			distance = dist(observations[j].x, observations[j].y, predicted[k].x, predicted[k].y);
			if ((k == 0) || (distance < min_distance)){
				min_distance = distance;
				observations[j].id = k;
			}
		}
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
	const vector<LandmarkObs> &observations,
	const Map &map_landmarks) {
	/**
	* TODO: Update the weights of each particle using a mult-variate Gaussian
	*   distribution. You can read more about this distribution here:
	*   TODO:使用多变量高斯分布更新每个粒子的权重。你可以在这里阅读更多关于这个分布的信息
	*   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	*   NOTE: The observations are given in the VEHICLE'S coordinate system.
	*   Your particles are located according to the MAP'S coordinate system.
	*   You will need to transform between the two systems. Keep in mind that
	*   this transformation requires both rotation AND translation (but no scaling).
	*   The following is a good resource for the theory:
	*   注：观测值以车辆坐标系给出。粒子是根据地图的坐标系定位的。你需要在两个系统之间转换。请记住，此转换需要旋转和平移（但不需要缩放）。
	*   以下是一个很好的理论资源：
	*   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	*   and the following is a good resource for the actual equation to implement
	*   (look at equation 3.33) 下面是一个很好的资源来实现实际的等式（参见方程式3.33）http://planning.cs.uiuc.edu/node99.html
	*/
	for (int i = 0; i < num_particles; i++){

		vector<LandmarkObs> in_range_landmks;
		for (int j = 0; j < map_landmarks.landmark_list.size(); j++){
			double dx = particles[i].x - map_landmarks.landmark_list[j].x_f;
			double dy = particles[i].y - map_landmarks.landmark_list[j].y_f;
			if ((dx*dx + dy*dy) <= (sensor_range*sensor_range)){
				in_range_landmks.push_back(LandmarkObs{ map_landmarks.landmark_list[j].id_i, map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f });
			}
		}

		vector<LandmarkObs> trans_observations;
		for (int j = 0; j < observations.size(); j++){
			double trans_x = particles[i].x + (cos(particles[i].theta)*observations[j].x - sin(particles[i].theta)*observations[j].y);
			double trans_y = particles[i].y + (sin(particles[i].theta)*observations[j].x + cos(particles[i].theta)*observations[j].y);
			trans_observations.push_back(LandmarkObs{ observations[j].id, trans_x, trans_y });
		}

		dataAssociation(in_range_landmks, trans_observations);

		particles[i].weight = 1.0;

		for (int j = 0; j < trans_observations.size(); j++){
			double single_weight;
			single_weight = multiv_prob(std_landmark[0], std_landmark[1], trans_observations[j].x, trans_observations[j].y,
				in_range_landmks[trans_observations[j].id].x, in_range_landmks[trans_observations[j].id].y);
			particles[i].weight = particles[i].weight*single_weight;
		}
	}


}

void ParticleFilter::resample() {
	/**
	* TODO: Resample particles with replacement with probability proportional
	*   to their weight. 用与粒子重量成比例的概率替换粒子，重新采样。
	* NOTE: You may find std::discrete_distribution helpful here.
	*   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	*/
	double max_weight;

	for (int i = 0; i < num_particles; i++){
		if (i == 0 || particles[i].weight>max_weight){
			max_weight = particles[i].weight;
		}
	}

	std::uniform_real_distribution<float> dist_weight(0.0, max_weight);
	std::uniform_int_distribution<int> dist_index(0, num_particles - 1);

	std::vector<Particle> resample_particles;
	double beta = 0.0;
	int index = dist_index(gen);

	for (int i = 0; i < num_particles; i++){
		beta += dist_weight(gen)*2.0;
		while (beta>particles[index].weight){
			beta = beta - particles[index].weight;
			index = (index + 1) % num_particles;
		}
		resample_particles.push_back(particles[index]);
	}

	particles = resample_particles;
}
	void ParticleFilter::SetAssociations(Particle& particle,
	const vector<int>& associations,
	const vector<double>& sense_x,
	const vector<double>& sense_y) {
	// particle: the particle to which assign each listed association,
	//   and association's (x,y) world coordinates mapping
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates
	// 粒子：指定每个列出的关联的粒子，以及关联的（x，y）世界坐标映射
	// 关联：与每个列出的关联一起使用的地标id
	// sense_x：关联x映射已转换为世界坐标
	// sense_y：关联y映射已转换为世界坐标
	particle.associations= associations;
	particle.sense_x = sense_x;
	particle.sense_y = sense_y;
	}

	string ParticleFilter::getAssociations(Particle best) {
	vector<int> v = best.associations;
	std::stringstream ss;
	copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length()-1);  // get rid of the trailing space
	return s;
	}

	string ParticleFilter::getSenseCoord(Particle best, string coord) {
	vector<double> v;

	if (coord == "X") {
	v = best.sense_x;
	} else {
	v = best.sense_y;
	}

	std::stringstream ss;
	copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length()-1);  // get rid of the trailing space 去掉后面的空格
	return s;
	}