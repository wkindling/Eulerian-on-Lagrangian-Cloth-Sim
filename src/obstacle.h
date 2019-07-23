#pragma once
#ifndef OBSTACLE_H
#define OBSTACLE_H


#include <Eigen/Dense>

#include <vector>
#include <memory>



class PointCloud
{
public:
	PointCloud() { point_num = 0; }
	virtual ~PointCloud() {}

	int point_num;
	Eigen::MatrixXd points_position;
	Eigen::MatrixXd points_normal;

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};


class Box
{
public:
	Box(Eigen::Vector3d _scale, Eigen::Vector3d _position);
	virtual ~Box() {};

private:
	void createBoxBuffer();

public:
	int point_num;
	int edge_num;

	float* buffer;
	int buffer_size;

	Eigen::Vector3d scale;
	Eigen::Vector3d position;

	Eigen::Matrix<double,3, 6> face_norm;
	Eigen::Matrix<double, 2, 12> edge_face;
	Eigen::Matrix<double, 3, 8> vert_edge;
	Eigen::Matrix<double, 1, 12> edge_direction;
	
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class Obstacle
{
public:
	Obstacle();
	virtual ~Obstacle() {};

	double collision_threshold;

	int box_num;
	std::vector<Box*> boxes;

	PointCloud* point_cloud;
};

#endif




