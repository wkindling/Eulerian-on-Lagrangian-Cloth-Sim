#pragma once
#ifndef PREPROCESSOR_H
#define PREPROCESSOR_H

#include "ArcSim/mesh.hpp"
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "Collisions.h"



class Preprocessor
{
public:
	Preprocessor() {}
	virtual ~Preprocessor() {}

public:
	void preprocess(Mesh& mesh, const Eigen::MatrixXd& boundary, std::vector <Collision* > collisions);
	
	void addEOL(Mesh& mesh, const Eigen::MatrixXd& boundary, std::vector<Collision*> collisions);
	void markPreserve(Mesh& mesh);
	void cleanMesh(Mesh& mesh, const Eigen::MatrixXd& boundary);

private:
	bool onBoundary(const Eigen::MatrixXd& boundary, const Node* node);
	bool onBoundaryA(const Eigen::MatrixXd& boundary, const Eigen::Vector3d& P);
	bool onBoundaryB(const Eigen::MatrixXd& boundary, const Eigen::Vector3d& P, const Vec3& e2);
};





#endif