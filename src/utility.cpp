#include "utility.h"

MatrixXd computeF(const Face* face)
{
	MatrixXd Dx(3, 2);
	Matrix2d DX;
	MatrixXd F;

	Dx << face->v[1]->node->x[0] - face->v[0]->node->x[0], face->v[2]->node->x[0] - face->v[0]->node->x[0],
		  face->v[1]->node->x[1] - face->v[0]->node->x[1], face->v[2]->node->x[1] - face->v[0]->node->x[1],
		  face->v[1]->node->x[2] - face->v[0]->node->x[2], face->v[2]->node->x[2] - face->v[0]->node->x[2];

	DX << face->v[1]->u[0] - face->v[0]->u[0], face->v[2]->u[0] - face->v[0]->u[0],
  		  face->v[1]->u[1] - face->v[0]->u[1], face->v[2]->u[1] - face->v[0]->u[1];

	F = Dx * DX.inverse();
	return F;
}

double vertLineProj(const Vec3& A, const Vec3& B, const Vec3& P)
{
	Vec3 AP = P - A;
	Vec3 AB = B - A;
	
	double proj_coordinate = dot(AP, AB) / dot(AB, AB);

	return proj_coordinate;
}

double vertLineDist(const Vec3& A, const Vec3& B, const Vec3& P)
{
	Vec3 AP = P - A;
	Vec3 AB = B - A;
	
	double proj_coordinate = dot(AP, AB) / dot(AB, AB);

	Vec3 proj = (1.0 - proj_coordinate)*A + proj_coordinate * B;
	Vec3 dist = P - proj;

	return norm(dist);
}

//Global Matrix Filling Functions
void fillLB(std::vector<T>& K, const Eigen::Matrix3d& KLL, int index)
{
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			K.push_back(T(index + i, index + j, KLL(i, j)));
		}
	}
}

void fillLLB(std::vector<T>& K, const Eigen::Matrix3d& KLL, int index1, int index2)
{
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			K.push_back(T(index1 + i, index2 + j, KLL(i, j)));
			K.push_back(T(index2 + j, index1 + i, KLL(i, j)));
		}
	}
}

void fillEB(std::vector<T>& K, const Eigen::Matrix2d& KEE, int index)
{
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			K.push_back(T(index + i, index + j, KEE(i, j)));
		}
	}
}

void fillEEB(std::vector<T>& K, const Eigen::Matrix2d& KEE, int index1, int index2)
{
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			K.push_back(T(index1 + i, index2 + j, KEE(i, j)));
			K.push_back(T(index2 + j, index1 + i, KEE(i, j)));
		}
	}
}

void fillELB(std::vector<T>& K, const Eigen::MatrixXd& KEL, int index1, int index2)
{
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			K.push_back(T(index1 + i, index2 + j, KEL(i, j)));
			K.push_back(T(index2 + j, index1 + i, KEL(i, j)));
		}
	}
}

void fillLEB(std::vector<T>& K, const Eigen::MatrixXd& KLE, int index1, int index2)
{
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			K.push_back(T(index1 + i, index2 + j, KLE(i, j)));
			K.push_back(T(index2 + j, index1 + i, KLE(i, j)));
		}
	}
}

void fillLM(std::vector<T>& K, std::vector<T>&M, const Eigen::Matrix3d& KLL, const Eigen::Matrix3d& MLL, int index, double damping, double h)
{
	Eigen::Matrix3d MDK = MLL + damping * h*h*KLL;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			M.push_back(T(index + i, index + j, MLL(i, j)));
			K.push_back(T(index + i, index + j, MDK(i, j)));
		}
	}

}

void fillLLM(std::vector<T>& K, std::vector<T>& M, const Eigen::Matrix3d& KLL, const Eigen::Matrix3d& MLL, int index1, int index2, double damping, double h)
{
	Eigen::Matrix3d MDK = MLL + damping * h*h*KLL;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			M.push_back(T(index1 + i, index2 + j, MLL(i, j)));
			M.push_back(T(index2 + j, index1 + i, MLL(i, j)));
			K.push_back(T(index1 + i, index2 + j, MDK(i, j)));
			K.push_back(T(index2 + j, index1 + i, MDK(i, j)));
		}
	}
}

void fillEM(std::vector<T>& K, std::vector<T>& M, const Eigen::Matrix2d& KLL, const Eigen::Matrix2d& MLL, int index, double damping, double h)
{
	Eigen::Matrix2d MDK = MLL + damping * h*h*KLL;
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			M.push_back(T(index + i, index + j, MLL(i, j)));
			K.push_back(T(index + i, index + j, MDK(i, j)));
		}
	}
}

void fillEEM(std::vector<T>& K, std::vector<T>& M, const Eigen::Matrix2d& KLL, const Eigen::Matrix2d& MLL, int index1, int index2, double damping, double h)
{
	Eigen::Matrix2d MDK = MLL + damping * h*h*KLL;
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			M.push_back(T(index1 + i, index2 + j, MLL(i, j)));
			M.push_back(T(index2 + j, index1 + i, MLL(i, j)));
			K.push_back(T(index1 + i, index2 + j, MDK(i, j)));
			K.push_back(T(index2 + j, index1 + i, MDK(i, j)));
		}
	}
}

void fillELM(std::vector<T>& K, std::vector<T>& M, const Eigen::MatrixXd& KEL, const Eigen::MatrixXd& MEL, int index1, int index2, double damping, double h)
{
	Eigen::MatrixXd MDK = MEL + damping * h*h*KEL;
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			M.push_back(T(index1 + i, index2 + j, MEL(i, j)));
			M.push_back(T(index2 + j, index1 + i, MEL(i, j)));
			K.push_back(T(index1 + i, index2 + j, MDK(i, j)));
			K.push_back(T(index2 + j, index1 + i, MDK(i, j)));
		}
	}
}

void fillLEM(std::vector<T>& K, std::vector<T>& M, const Eigen::MatrixXd& KLE, const Eigen::MatrixXd& MLE, int index1, int index2, double damping, double h)
{
	Eigen::MatrixXd MDK = MLE + damping * h*h*KLE;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			M.push_back(T(index1 + i, index2 + j, MLE(i, j)));
			M.push_back(T(index2 + j, index1 + i, MLE(i, j)));
			K.push_back(T(index1 + i, index2 + j, MDK(i, j)));
			K.push_back(T(index2 + j, index1 + i, MDK(i, j)));
		}
	}
}

//Numerical computation
Matrix2d poldec(const Matrix2d& M)
{
	double m00 = M(0, 0);
	double m01 = M(0, 1);
	double m10 = M(1, 0);
	double m11 = M(1, 1);
	double detM = m00 * m11 - m01 * m10;
	int sign = 1;
	if (detM < 0)
	{
		sign = -1;
	}
	else if (detM == 0)
	{
		sign = 0;
	}
	Matrix2d TM;
	TM << m11, -m10, -m01, m00;
	Eigen::Matrix2d Q = M + sign * TM;

	double len0 = Q.col(0).norm();
	double len1 = Q.col(1).norm();

	Q.col(0) = Q.col(0) / len0;
	Q.col(1) = Q.col(1) / len1;

	return Q;
}

bool PositiveDefinite(MatrixXd G)
{
	VectorXcd E = G.eigenvalues();

	bool res = true;

	for (int i = 0; i < E.rows(); i++)
	{
		if (E(i).real() <= 0) { res = false; break; }
	}
	return res;
}

//Arcsim Vec & Eigen Vectorxd translation
Vec2 e2v(const Vector2d v)
{
	return Vec2(v(0), v(1));
}

Vec3 e2v(const Vector3d v)
{
	return Vec3(v(0), v(1), v(2));
}

Vector2d v2e(const Vec2 v)
{
	return Vector2d(v[0], v[1]);
}

Vector2d v322e(const Vec3 v)
{
	return Vector2d(v[0], v[1]);
}

Vector3d v2e(const Vec3 v)
{
	return Vector3d(v[0], v[1], v[2]);
}