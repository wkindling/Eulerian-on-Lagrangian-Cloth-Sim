#pragma once
#ifndef FORCE_H
#define FORCE_H

#include <vector>
#include <memory>
#include <string>

#include "utility.h"
#include "cloth.h"

#include "ArcSim/mesh.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

class Force
{
public:
	Force() {};
	~Force() {};

	void computeForce(const Mesh& mesh, const Material& mat, const Eigen::Vector3d& gravity, double h);

private:

	void computeBendingForce(const Mesh& mesh, const Material& mat, Eigen::VectorXd& f, std::vector<T>& _MDK, double h);
	void computeBendingForceEOL(const Edge* edge, const Vert* v0, const Vert* v1, const Vert* v2, const Vert* v3,
								Eigen::VectorXd& fb, Eigen::MatrixXd& Kb);

	void computeFaceForce(const Mesh& mesh, Eigen::VectorXd& f, std::vector<T>& _MDK, std::vector<T>& _M, const Eigen::Vector3d& gravity, double h);
	void computeMembraneForceEOL(const Face* face, Eigen::VectorXd& fm, Eigen::MatrixXd& Km);
	void computeInertialForceEOL(const Face* face, Eigen::VectorXd& fi, Eigen::MatrixXd& Mi);


	void computeBending(const Vec3& x0, const Vec3& x1, const Vec3& x2, const Vec3& x3,
		const Vec3& X0, const Vec3& X1, const Vec3& X2, const Vec3& X3,
		double beta, double &W, Eigen::VectorXd &f, Eigen::MatrixXd &K);
	void computeInertial(const Vec3& xa, const Vec3& xb, const Vec3& xc,
		const Vec3& Xa, const Vec3& Xb, const Vec3& Xc, Eigen::Vector3d gravity,
		double rho, double &W, Eigen::VectorXd &f, Eigen::MatrixXd &M);
	void computeMembrane(const Vec3& xa, const Vec3& xb, const Vec3& xc,
		const Vec3& Xa, const Vec3& Xb, const Vec3& Xc,
		double e, double nu, Eigen::Matrix<double, 2, 3>& P, Eigen::Matrix2d& Q,
		double &W, Eigen::VectorXd& f, Eigen::MatrixXd& K);

public:
	Eigen::VectorXd f;
	Eigen::SparseMatrix<double> M;
	Eigen::SparseMatrix<double> K;

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};



#endif