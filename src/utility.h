#pragma once
#ifndef UTILITY_H
#define UTILITY_H

#include "ArcSim/mesh.hpp"
#include "ArcSim/vectors.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>


using namespace Eigen;
using namespace std;

typedef Triplet<double> T;

//Geometry computation functions
MatrixXd computeF(const Face* face);
double vertLineProj(const Vec3& A, const Vec3& B, const Vec3& P);
double vertLineDist(const Vec3& A, const Vec3& B, const Vec3& P);

//Global Matrix Filling Functions
void fillLB	(vector<T>& _MDK, const Matrix3d& KLL, int index);
void fillLLB(vector<T>& _MDK, const Matrix3d& KLL, int index1, int index2);
void fillEB	(vector<T>& _MDK, const Matrix2d& KEE, int index);
void fillEEB(vector<T>& _MDK, const Matrix2d& KEE, int index1, int index2);
void fillELB(vector<T>& _MDK, const MatrixXd& KEL, int index1, int index2);
void fillLEB(vector<T>& _MDK, const MatrixXd& KLE, int index1, int index2);

void fillLM	(vector<T>& _MDK, vector<T>& _M, const Matrix3d& KLL, const Matrix3d& MLL, int index,  double damping, double h);
void fillLLM(vector<T>& _MDK, vector<T>& _M, const Matrix3d& KLL, const Matrix3d& MLL, int index1, int index2,	   double damping, double h);
void fillEM	(vector<T>& _MDK, vector<T>& _M, const Matrix2d& KLL, const Matrix2d& MLL, int index,  double damping, double h);
void fillEEM(vector<T>& _MDK, vector<T>& _M, const Matrix2d& KLL, const Matrix2d& MLL, int index1, int index2,	   double damping, double h);
void fillELM(vector<T>& _MDK, vector<T>& _M, const MatrixXd& KEL, const MatrixXd& MEL, int index1, int index2,     double damping, double h);
void fillLEM(vector<T>& _MDK, vector<T>& _M, const MatrixXd& KLE, const MatrixXd& MLE, int index1, int index2,     double damping, double h);

//Numerical computation
Matrix2d poldec(const Matrix2d& M);
bool PositiveDefinite(MatrixXd G);

//Arcsim Vec & Eigen Vectorxd translation
Vec2 e2v(const Eigen::Vector2d);
Vec3 e2v(const Eigen::Vector3d);
Vector2d v2e(const Vec2);
Vector3d v2e(const Vec3);
Vector2d v322e(const Vec3 v);


















#endif