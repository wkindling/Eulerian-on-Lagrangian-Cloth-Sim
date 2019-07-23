#pragma once
#ifndef CLOTH_H
#define CLOTH_H

#include <vector>
#include <memory>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "ArcSim/mesh.hpp"
#include "Mosek/QuadProgMosek.h"

class Force;
class Obstacle;
class Constraint;
class Preprocessor;

struct Material
{
	double density;
	double e;
	double nu;
	double beta;
	double damping1;
	double damping2;
};

struct Remeshing
{
	double refine_angle, refine_compression, refine_velocity;
	double size_min, size_max;
	double aspect_min;
};

class Cloth
{
public:
	Cloth(Eigen::Vector2i res, Eigen::VectorXd& xX00, Eigen::VectorXd& xX01, Eigen::VectorXd& xX10, Eigen::VectorXd& xX11);
	virtual ~Cloth() {};
	
	bool mosekSolve(const Eigen::SparseMatrix<double>& MDK, const Eigen::VectorXd& b,
		const Eigen::SparseMatrix<double>& Aeq, const Eigen::VectorXd& beq,
		const Eigen::SparseMatrix<double>& Aineq, const Eigen::VectorXd& bineq,
		Eigen::VectorXd& v);
	void saveOldMesh();
	void velocityTransfer();
	void step(Obstacle* obstacle, const Eigen::Vector3d& gravity, double h, double t);
	void solve(double h);

	void createClothBuffer();
	

public:
	Mesh mesh;
	Remeshing remeshing;
	Material material;

	Mesh old_mesh;

	Constraint* constraint;
	Force* force;
	Preprocessor* preprocessor;

	Eigen::MatrixXd boundary;
	Eigen::VectorXd v;

	float* buffer;
	float* mesh_buffer;

	int buffer_size;
	int mesh_buffer_size;
};




#endif


