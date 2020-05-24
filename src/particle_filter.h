/**
* particle_filter.h
* 2D particle filter class.
*
* Created on: Dec 12, 2016
* Author: Tiffany Huang
*/

#ifndef PARTICLE_FILTER_H_
#define PARTICLE_FILTER_H_

#include <string>
#include <vector>
#include "helper_functions.h"

struct Particle {
	int id;
	double x;
	double y;
	double theta;
	double weight;
	std::vector<int> associations;
	std::vector<double> sense_x;
	std::vector<double> sense_y;
};


class ParticleFilter {
public:
	// Constructor
	// @param num_particles Number of particles @参数 粒子数
	ParticleFilter() : num_particles(0), is_initialized(false) {}

	// Destructor
	~ParticleFilter() {}

	/**
	* init Initializes particle filter by initializing particles to Gaussian
	*   distribution around first position and all the weights to 1.
	*   init通过将粒子初始化为第一个位置周围的高斯分布并将所有权重设置为1来初始化粒子过滤器。
	* @param x Initial x position [m] (simulated estimate from GPS) x初始x位置[m]（来自GPS的模拟估计）
	* @param y Initial y position [m] y初始y位置[m]
	* @param theta Initial orientation [rad] theta初始方向[rad]
	* @param std[] Array of dimension 3 [standard deviation of x [m],
	*   standard deviation of y [m], standard deviation of yaw [rad]]
	* std[]3维数组[x[m]的标准偏差，y[m]的标准偏差，偏航的标准偏差[rad]]
	*/

	void init(double x, double y, double theta, double std[]);

	/**
	* prediction Predicts the state for the next time step
	*   using the process model. 预测使用流程模型预测下一个时间步骤的状态。
	* @param delta_t Time between time step t and t+1 in measurements [s]
	* @param std_pos[] Array of dimension 3 [standard deviation of x [m],
	*   standard deviation of y [m], standard deviation of yaw [rad]]
	* @param velocity Velocity of car from t to t+1 [m/s]
	* @param yaw_rate Yaw rate of car from t to t+1 [rad/s]
	*/
	void prediction(double delta_t, double std_pos[], double velocity,
		double yaw_rate);

	/**
	* dataAssociation Finds which observations correspond to which landmarks
	*   (likely by using a nearest-neighbors data association).
	* 查找哪些观测值对应于哪些地标（可能是使用最近邻数据关联）。
	* @param predicted Vector of predicted landmark observations
	* @param observations Vector of landmark observations
	*/
	void dataAssociation(std::vector<LandmarkObs> predicted,
		std::vector<LandmarkObs>& observations);

	/**
	* updateWeights Updates the weights for each particle based on the likelihood
	*   of the observed measurements. 根据观察到的测量值的可能性更新每个粒子的权重。
	* @param sensor_range Range [m] of sensor
	* @param std_landmark[] Array of dimension 2
	*   [Landmark measurement uncertainty [x [m], y [m]]]
	* @param observations Vector of landmark observations
	* @param map Map class containing map landmarks
	*/
	void updateWeights(double sensor_range, double std_landmark[],
		const std::vector<LandmarkObs> &observations,
		const Map &map_landmarks);

	/**
	* resample Resamples from the updated set of particles to form
	*   the new set of particles. 从更新的粒子集重新采样以形成新的粒子集。
	*/
	void resample();

	/**
	* Set a particles list of associations, along with the associations'
	*   calculated world x,y coordinates 设置关联的粒子列表，以及关联的计算世界x，y坐标
	* This can be a very useful debugging tool to make sure transformations
	*   are correct and assocations correctly connected
	*这是一个非常有用的调试工具，可以确保转换是正确的，并且关联是正确连接的
	*/
	void SetAssociations(Particle& particle, const std::vector<int>& associations,
		const std::vector<double>& sense_x,
		const std::vector<double>& sense_y);

	/**
	* initialized Returns whether particle filter is initialized yet or not.
	* 返回粒子过滤器是否已初始化。
	*/
	const bool initialized() const {
		return is_initialized;
	}

	/**
	* Used for obtaining debugging information related to particles.
	* 用于获取与粒子相关的调试信息。
	*/
	std::string getAssociations(Particle best);
	std::string getSenseCoord(Particle best, std::string coord);

	// Set of current particles 当前粒子集
	std::vector<Particle> particles;

private:
	// Number of particles to draw 要绘制的粒子数
	int num_particles;

	// Flag, if filter is initialized 筛选器已初始化标志，如果
	bool is_initialized;

	// Vector of weights of all particles 所有粒子的权重向量
	std::vector<double> weights;
};

#endif  // PARTICLE_FILTER_H_