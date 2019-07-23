#pragma once
#ifndef __Collisions__
#define __Collisions__

#include "ArcSim\mesh.hpp"
#include "obstacle.h"
#include <memory>
#include <vector>
#include <Eigen/Dense>

class Edge_
{
public:
	Edge_() {
		vs = Eigen::Vector4i::Zero();
		faces = Eigen::Vector2i::Zero();
		internal = false;
		angle = 0;
		ns[0] = Eigen::Vector3d::Zero();
		ns[1] = Eigen::Vector3d::Zero();
	}
	virtual ~Edge_() {

	}

	Eigen::Vector4i vs;
	Eigen::Vector2i faces;
	bool internal;
	double angle;
	Eigen::Vector3d ns[2];
};

void create_edge(
	std::vector<Edge_* > &es, // output
	const Eigen::MatrixXi &fs, // input
	const Eigen::MatrixXd &vs  // input
);

class Collision
{
public:
	Collision() 
	{
		dist = 0;
		nor1 = Eigen::Vector3d::Zero();
		nor2 = Eigen::Vector3d::Zero();
		pos1 = Eigen::Vector3d::Zero();
		pos2 = Eigen::Vector3d::Zero();
		count1 = 0;
		count2 = 0;
		verts1 = Eigen::Vector3i::Zero();
		verts2 = Eigen::Vector3i::Zero();
		weights1 = Eigen::Vector3d::Zero();
		weights2 = Eigen::Vector3d::Zero();
		tri1 = -1;
		tri2 = -1;
		edge2 = -1;
		edgeDir = Eigen::Vector3d::Zero();
	}
	Collision(double dist, Eigen::Vector3d nor1, Eigen::Vector3d nor2, Eigen::Vector3d pos1, Eigen::Vector3d pos2,
		int count1, int count2, Eigen::Vector3i verts1, Eigen::Vector3i verts2, Eigen::Vector3d weights1, Eigen::Vector3d weights2, int tri1, int tri2) {
		this->dist = dist;
		this->nor1 = nor1;
		this->nor2 = nor2;
		this->pos1 = pos1;
		this->pos2 = pos2;
		this->pos1_ = pos1_;
		this->count1 = count1;
		this->count2 = count2;
		this->verts1 = verts1;
		this->verts2 = verts2;
		this->weights1 = weights1;
		this->weights2 = weights2;
		this->tri1 = tri1;
		this->tri2 = tri2;
	}
	~Collision() {

	}

	double dist; // distance
	Eigen::Vector3d nor1; // normal
	Eigen::Vector3d nor2; // normal on cloth
	Eigen::Vector3d pos1; // position on box
	Eigen::Vector3d pos2; // position on cloth
	Eigen::Vector3d pos1_; // position on box offset to be slightly inside
	int count1; // Collision type for box:   1=vert, 2=edge, 3=face
	int count2; // Collision type for cloth: 1=vert, 2=edge, 3=face
	Eigen::Vector3i verts1; // The collided vertices on box
	Eigen::Vector3i verts2; // The collided vertices on cloth
	Eigen::Vector3d weights1; // The vertex weights on box
	Eigen::Vector3d weights2; // The vertex weights on cloth
	int tri1; // Triangle index for box
	int tri2; // Triangle index for cloth
	std::vector<int> edge1; // all surrounding edge indices for he box
	int edge2; // edge index for cloth
	Eigen::Vector3d edgeDir; // Direction of box edge  :: This a better way
};

void ForBox(
	std::vector<Collision*> &colls,
	double thres,
	const Eigen::Vector3d &scale,
	const Eigen::Matrix4d &E1,
	const Eigen::MatrixXd &v2,
	const Eigen::MatrixXi &f2,
	const Eigen::VectorXi &isEOL2,
	bool EOL);

void ForPoint(
	std::vector<Collision*> &colls,
	double thres,
	const Eigen::MatrixXd &v1,
	const Eigen::MatrixXd &n1,
	const Eigen::MatrixXd &v2,
	const Eigen::MatrixXi &f2,
	bool EOL);

void Find_Collision(const Mesh& mesh, const Obstacle* obs, std::vector<Collision* > &cls);

void Find_Collision_(const Mesh& mesh, const Obstacle* obs, std::vector<Collision*> &cls);

#endif