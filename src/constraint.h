#pragma once
#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include <vector>
#include <memory>
#include <string>

#include "ArcSim/mesh.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/StdVector>

class Obstacle;

class Constraint
{
public:
	Constraint();
	virtual ~Constraint() {}

	void updateConstraint(const Obstacle* obstacle);
	void applyConstraint(const Mesh& mesh, const Obstacle* obstacle, double h ,double t);

public:
	bool occur_collision;

	Eigen::SparseMatrix<double> Aeq;
	Eigen::SparseMatrix<double> Aineq;

	Eigen::VectorXd beq;
	Eigen::VectorXd bineq;

	Eigen::MatrixXd constraint_matrix;

};








#endif