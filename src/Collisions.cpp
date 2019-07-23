#include "Collisions.h"
#include <Eigen\Dense>
#define _USE_MATH_DEFINES
#include <cmath>
#include <string.h>
#include <iostream>
#include <random>


using namespace std;
using namespace Eigen;

#define EPSILON 1e-6
#define CROSS(dest,v1,v2) \
          dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
          dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
          dest[2]=v1[0]*v2[1]-v1[1]*v2[0];
#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])
#define SUB(dest,v1,v2) \
          dest[0]=v1[0]-v2[0]; \
          dest[1]=v1[1]-v2[1]; \
          dest[2]=v1[2]-v2[2]; 
#define ZERO -EPSILON

int intersect_triangle3_inc(
	const double *orig, const double *dir,
	const double *vert0, const double *vert1, const double *vert2,
	double *t, double *u, double *v)
{
	double edge1[3], edge2[3], tvec[3], pvec[3], qvec[3];
	double det, inv_det;

	/* find vectors for two edges sharing vert0 */
	SUB(edge1, vert1, vert0);
	SUB(edge2, vert2, vert0);

	/* begin calculating determinant - also used to calculate U parameter */
	CROSS(pvec, dir, edge2);

	/* if determinant is near zero, ray lies in plane of triangle */
	det = DOT(edge1, pvec);

	/* calculate distance from vert0 to ray origin */
	SUB(tvec, orig, vert0);
	inv_det = 1.0 / det;

	CROSS(qvec, tvec, edge1);

	if (det > EPSILON)
	{
		*u = DOT(tvec, pvec);
		if (*u < -ZERO || *u > det + ZERO)
			return 0;

		/* calculate V parameter and test bounds */
		*v = DOT(dir, qvec);
		if (*v < -ZERO || *u + *v > det + ZERO)
			return 0;

	}
	else if (det < -EPSILON)
	{
		/* calculate U parameter and test bounds */
		*u = DOT(tvec, pvec);
		if (*u > ZERO || *u < det - ZERO)
			return 0;

		/* calculate V parameter and test bounds */
		*v = DOT(dir, qvec);
		if (*v > ZERO || *u + *v < det - ZERO)
			return 0;
	}
	else return 0;  /* ray is parallell to the plane of the triangle */

	*t = DOT(edge2, qvec) * inv_det;
	(*u) *= inv_det;
	(*v) *= inv_det;

	return 1;
}
// -----------------------------------------------------------

// Box geometry data
double v1_data[] = {
	-1, -1, -1,  1,
	-1, -1,  1,  1,
	-1,  1, -1,  1,
	-1,  1,  1,  1,
	1, -1, -1,  1,
	1, -1,  1,  1,
	1,  1, -1,  1,
	1,  1,  1,  1,
	-1,  0,  0,  1,
	1,  0,  0,  1,
	0, -1,  0,  1,
	0,  1,  0,  1,
	0,  0, -1,  1,
	0,  0,  1,  1,
};

int f1_data[] = {
	0,     8,     2,
	1,     8,     0,
	3,     8,     1,
	2,     8,     3,
	4,    10,     0,
	5,    10,     4,
	1,    10,     5,
	0,    10,     1,
	6,     9,     4,
	7,     9,     6,
	5,     9,     7,
	4,     9,     5,
	2,    11,     6,
	3,    11,     2,
	7,    11,     3,
	6,    11,     7,
	1,    13,     3,
	5,    13,     1,
	7,    13,     5,
	3,    13,     7,
	0,    12,     4,
	2,    12,     0,
	6,    12,     2,
	4,    12,     6,
};


int ev_data[] = {
	0, 1, 8,10,
	2, 0, 8,12,
	1, 3, 8,13,
	3, 2, 8,11,
	0, 4,10,12,
	5, 1,10,13,
	4, 5,10, 9,
	6, 2,11,12,
	4, 6, 9,12,
	3, 7,11,13,
	7, 5, 9,13,
	6, 7, 9,11,
};

// Keeping track of which edges each corner point is made up of
int ve_data[] = {
	0, 1, 4,
	0, 2, 5,
	1, 3, 7,
	2, 3, 9,
	4, 6, 8,
	5, 6, 10,
	7, 8, 11,
	9, 10, 11
};


int ef_data[] = {
	1,7,
	0,21,
	2,16,
	3,13,
	4,20,
	6,17,
	5,11,
	12,22,
	8,23,
	14,19,
	10,18,
	9,15
};

// 1. build a new edge
void create_edge(vector<Edge_*> &es, const MatrixXi &fs, const MatrixXd &vs) {

	struct FaceEdge
	{
		Vector3i verts;
		int face;
		int hash;
	};

	vector<FaceEdge*> tmp;
	for (int k = 0; k < fs.cols(); ++k) {
		for (int i = 0; i < 3; ++i) {
			FaceEdge* faceEdge = new FaceEdge;
			tmp.push_back(faceEdge);
			faceEdge->verts << (i + 0) % 3, (i + 1) % 3, (i + 2) % 3;
			faceEdge->face = k;
			int v0 = faceEdge->verts(0);
			int v1 = faceEdge->verts(1);
			int kmin = min(fs(v0, k), fs(v1, k)) + 1;
			int kmax = max(fs(v0, k), fs(v1, k)) + 1;
			faceEdge->hash = kmin + (3 * fs.cols() + 1)*kmax;
		}
	}

	// Sort by hash number
	struct Compare {
		bool operator() (const FaceEdge* e0, const FaceEdge* e1) const {
			return (e0->hash < e1->hash);
		}
	};
	Compare comp;
	
	sort(tmp.begin(), tmp.end(), comp);

	// Now the twin edges are neighbors in the list
	int k = 0;
	while (k < 3 * fs.cols()) {
		int one_f0 = tmp[k]->face;
		Vector2i e;
		e(0) = fs(tmp[k]->verts(0), one_f0);
		e(1) = fs(tmp[k]->verts(1), one_f0);
		Vector3d a = vs.block<3, 1>(0, fs(0, one_f0));
		Vector3d b = vs.block<3, 1>(0, fs(1, one_f0));
		Vector3d c = vs.block<3, 1>(0, fs(2, one_f0));
		Vector3d n = (b - a).cross(c - a).normalized();
		if (k < 3 * fs.cols() - 1 && tmp[k]->hash == tmp[k + 1]->hash) {
			// k and k+1 are twins
			int one_f1 = tmp[k + 1]->face;
			Vector3d xa1 = vs.block<3, 1>(0, fs(0, one_f1));
			Vector3d xb1 = vs.block<3, 1>(0, fs(1, one_f1));
			Vector3d xc1 = vs.block<3, 1>(0, fs(2, one_f1));
			Vector3d n_ = (xb1 - xa1).cross(xc1 - xa1).normalized();
			double angle = acos(n.dot(n_));
			auto edge = new Edge_;
			es.push_back(edge);
			edge->vs.segment<2>(0) = e; // edge vertex indices
			edge->vs(2) = fs(tmp[k]->verts(2), one_f0); // x2 index in figure
			edge->vs(3) = fs(tmp[k + 1]->verts(2), one_f1); // x3 index in figure
			edge->faces << one_f0, one_f1; // t0 and t1 in figure
			edge->internal = true;
			edge->angle = angle; // ( : concave edges)
			edge->ns[0] = n;
			edge->ns[1] = n_;
			k += 2;
		}
		else {
			// k is an external edge (singleton)
			auto edge = new Edge_;
			es.push_back(edge);
			edge->vs.segment<2>(0) = e; // edge vertex indices
			edge->vs(2) = fs(tmp[k]->verts(2), one_f0); // x2 index in figure
			edge->vs(3) = -1; // x3 index in figure
			edge->faces << one_f0, -1; // t0 and t1 in figure
			edge->internal = false;
			edge->angle = 1e9; // infinite dihedral angle
			edge->ns[0] = n;
			k += 1;
		}
	}
}
// 2. build face normal based on vertices
MatrixXd create_face_n(const MatrixXi &faces, const MatrixXd &verts) {

	MatrixXd ns = MatrixXd::Zero(3, faces.cols());
	for (int k = 0; k < faces.cols(); ++k) {
		const Vector3i &f = faces.col(k);
		const Vector3d &a = verts.block<3, 1>(0, f(0));
		const Vector3d &b = verts.block<3, 1>(0, f(1));
		const Vector3d &c = verts.block<3, 1>(0, f(2));
		ns.col(k) = (b - a).cross(c - a).normalized();
	}
	return ns;
}
// 3. build vertex normals
MatrixXd create_vert_n(const MatrixXi &faces, const MatrixXd &verts) {

	MatrixXd ns = MatrixXd::Zero(3, verts.cols());
	vector<double> angles(verts.cols(), 0); 
	for (int k = 0; k < faces.cols(); ++k) {
		const Vector3i &f = faces.col(k);
		const Vector3d &a = verts.block<3, 1>(0, f(0));
		const Vector3d &b = verts.block<3, 1>(0, f(1));
		const Vector3d &c = verts.block<3, 1>(0, f(2));
		Vector3d ba = b - a;
		Vector3d cb = c - b;
		Vector3d ac = a - c;
		Vector3d N = ba.cross(-ac).normalized();
		ba.normalize();
		cb.normalize();
		ac.normalize();
		double angle1 = acos(ba.dot(-ac));
		double angle2 = acos(cb.dot(-ba));
		double angle3 = acos(ac.dot(-cb));
		ns.col(f(0)) += angle1 * N;
		ns.col(f(1)) += angle2 * N;
		ns.col(f(2)) += angle3 * N;
		angles[f(0)] += angle1;
		angles[f(1)] += angle2;
		angles[f(2)] += angle3;
	}
	for (int k = 0; k < verts.cols(); k++) {
		Vector3d N = ns.col(k) / angles[k];
		ns.col(k) = N.normalized();
	}
	return ns;
}


// 4. AABB 
void AABB(Matrix<double, 6, 1> &aabb, const MatrixXd &V)
{
	aabb(0) = V.row(0).minCoeff();
	aabb(1) = V.row(1).minCoeff();
	aabb(2) = V.row(2).minCoeff();
	aabb(3) = V.row(0).maxCoeff();
	aabb(4) = V.row(1).maxCoeff();
	aabb(5) = V.row(2).maxCoeff();
}

// 5. Check between two AABBs
bool check_AABB(const Matrix<double, 6, 1> &aabb1, const Matrix<double, 6, 1> &aabb2)
{
	double thresh = 1e-3;
	Vector3d threshVec(thresh, thresh, thresh);
	Vector3d min1 = aabb1.segment<3>(0) - threshVec;
	Vector3d max1 = aabb1.segment<3>(3) + threshVec;
	Vector3d min2 = aabb2.segment<3>(0) - threshVec;
	Vector3d max2 = aabb2.segment<3>(3) + threshVec;
	return
		max1(0) >= min2(0) &&
		min1(0) <= max2(0) &&
		max1(1) >= min2(1) &&
		min1(1) <= max2(1) &&
		max1(2) >= min2(2) &&
		min1(2) <= max2(2);
}

// 6. Check rectangle intersection
int intersection(
	const Vector3d &x0, const Vector3d &dx,
	const Vector3d &xa, const Vector3d &xb, const Vector3d &xc,
	double &t)
{
	int intersect = 0;
	double unused1, unused2;

	// Compute xd and xe
	Vector3d xd = xa + 2.0*(xc - xa);
	Vector3d xe = xb + 2.0*(xc - xb);

	// Slow: check all four triangles
	intersect = intersect_triangle3_inc(x0.data(), dx.data(), xa.data(), xb.data(), xc.data(), &t, &unused1, &unused2);
	if (intersect) {
		return 1;
	}
	intersect = intersect_triangle3_inc(x0.data(), dx.data(), xb.data(), xd.data(), xc.data(), &t, &unused1, &unused2);
	if (intersect) {
		return 1;
	}
	intersect = intersect_triangle3_inc(x0.data(), dx.data(), xd.data(), xe.data(), xc.data(), &t, &unused1, &unused2);
	if (intersect) {
		return 1;
	}
	intersect = intersect_triangle3_inc(x0.data(), dx.data(), xe.data(), xa.data(), xc.data(), &t, &unused1, &unused2);
	if (intersect) {
		return 1;
	}

	// No intersection
	t = -1.0;
	return 0;
}


void ForBox(
	vector<Collision*> &colls,
	double thres,
	const Vector3d &scale,
	const Matrix4d &E1,
	const MatrixXd &verts2_,
	const MatrixXi &f2,
	const VectorXi &_,
	bool EOL)
{
	vector<Edge_*> e2;
	create_edge(e2, f2, verts2_);
	// Make a copy first so we can perturb
	MatrixXd v2 = verts2_;

	MatrixXd face1_n;
	MatrixXd verts_n;

	// The first body is always the box
	vector<Edge_*> e1;

	Matrix4d S = Matrix4d::Identity();
	S(0, 0) = 0.5 * scale(0);
	S(1, 1) = 0.5 * scale(1);
	S(2, 2) = 0.5 * scale(2);
	Matrix4d E = E1 * S;
	Map<Matrix<int, 2, 12, 0> > edge_f(ef_data);
	Map<Matrix<int, 4, 12, 0> > edge_v(ev_data);
	Map<Matrix<int, 3, 24, 0> > f1(f1_data);
	Map<Matrix<double, 4, 14, 0> > v1_(v1_data);
	Matrix<double, 4, 14> v1 = E * v1_;

	face1_n = create_face_n(f1, v1);
	verts_n = create_vert_n(f1, v1);
	for (int k = 0; k < 12; ++k) {
		auto edge = new Edge_;
		e1.push_back(edge);
		edge->vs = edge_v.col(k);
		edge->faces = edge_f.col(k);
		edge->internal = false; // for boxes only
		edge->ns[0] = face1_n.col(edge->faces(0));
		edge->ns[1] = face1_n.col(edge->faces(1));
		edge->angle = acos(edge->ns[0].dot(edge->ns[1]));
	}
	// Perturb cloth verts
	std::random_device device;
	std::mt19937 mt;
	std::uniform_real_distribution<> dis(-1.0, 1.0);
	mt.seed(1);
	for (int i2 = 0; i2 < v2.cols(); ++i2) {
		Vector3d r;
		r(0) = dis(mt)*thres*1e-3;
		r(1) = dis(mt)*thres*1e-3;
		r(2) = dis(mt)*thres*1e-3;
		v2.block<3, 1>(0, i2) += r;
	}

	// Precompute the face normals for the cloth
	MatrixXd faceNors2 = create_face_n(f2, v2);

	// Build AABBs
	Matrix<double, 6, 1> aabb_1;
	Matrix<double, 6, 1> aabb_2;
	AABB(aabb_1, v1);
	AABB(aabb_2, v2);

	if (!EOL) {
		// 1: triangle
		// 2: point
		for (int i = 0; i < v2.cols(); i++) {
			Vector3d x2 = v2.block<3, 1>(0, i);
			// AABB test
			Matrix<double, 6, 1> aabbV2;
			aabbV2.segment<3>(0) = x2;
			aabbV2.segment<3>(3) = x2;
			if (!check_AABB(aabbV2, aabb_1)) {
				continue;
			}
			Collision* c = NULL;
			
			bool inside = true;
			for (int j = 0; j < f1.cols(); j++) {
				Vector3i one_f1 = f1.col(j);
				Vector3d x1_a = v1.block<3, 1>(0, one_f1(0));
				Vector3d N = face1_n.col(j);
				Vector3d dx = x2 - x1_a;
				if (N.dot(dx) > 0.0) {
					inside = false;
					break;
				}
			}
			if (!inside)
				continue;
			for (int j = 0; j < f1.cols(); j++) {
				Vector3i one_f1 = f1.col(j);
				Vector3d x1_a = v1.block<3, 1>(0, one_f1(0));
				Vector3d x1_b = v1.block<3, 1>(0, one_f1(1));
				Vector3d x1_c = v1.block<3, 1>(0, one_f1(2));
				// Project x2 onto the triangle 1
				Vector3d N1 = face1_n.col(j);
				Vector3d dx = x2 - x1_a;
				if (N1.dot(dx) > 0.0)
					continue;
				Vector3d x1 = x2 - N1.dot(dx) * N1;
				dx = x2 - x1;
				if (dx.norm() > 5.0*thres) 
					continue;
				// calculate barycentric corrd
				double u, v, w;
				Vector3d t0 = x1_b - x1_a;
				Vector3d t1 = x1_c - x1_a;
				Vector3d t2 = x1 - x1_a;
				double dividend = t0.dot(t0) * t1.dot(t1) - t0.dot(t1) * t0.dot(t1);
				v = (t1.dot(t1) * t2.dot(t0) - t0.dot(t1) * t2.dot(t1)) / dividend;
				w = (t0.dot(t0) * t2.dot(t1) - t0.dot(t1) * t2.dot(t0)) / dividend;
				u = 1 - v - w;
				if (u < 0 || u > 1 || v < 0 || v > 1 || w < 0 || w > 1) continue;

				// Compute cloth normal
				Vector3d N2 = faceNors2.col(i);
				if (N2.dot(N1) < 0.0)
					N2 = -N2;
	
				// Create contact object
				Collision* c_temp = new Collision(dx.norm(), N1, N2, x1, x2, 3, 1, one_f1,
					Vector3i(i, -1, -1), Vector3d(u, v, w), Vector3d(1.0, 0.0, 0.0), j, -1);
				// Is this the closest one so far?
				if (!c || c->dist > dx.norm())
					c = c_temp;
			}
			if (c) colls.push_back(c);
		}
	}


	Map<Matrix<int, 3, 8, ColMajor> > ve_1(ve_data); // ADDED BY NICK
	
	// 1: point
	// 2: triangle
	for (int i = 0; i < 8; i++) { // only up to 8 corner points (ignore face points)
		const Vector3d &x1 = v1.block<3, 1>(0, i);
		// AABB test: check Vertex1 against Body2
		Matrix<double, 6, 1> aabb_1_V;
		aabb_1_V.segment<3>(0) = x1;
		aabb_1_V.segment<3>(3) = x1;
		if (!check_AABB(aabb_1_V, aabb_2))
			continue;
		Collision* c = NULL;
		Vector3d N1 = verts_n.block<3, 1>(0, i); // vertex normal
		for (int j = 0; j < f2.cols(); j++) {
			const Vector3i &one_f2 = f2.col(j);
			const Vector3d &x2_a = v2.block<3, 1>(0, one_f2(0));
			const Vector3d &x2_b = v2.block<3, 1>(0, one_f2(1));
			const Vector3d &x2_c = v2.block<3, 1>(0, one_f2(2));
			Vector3d N2 = faceNors2.col(j);
			// Make sure the triangle normal points outward wrt the box.
			if (N1.dot(N2) < 0.0) N2 = -N2;
			// Is x1 on the correct side?
			Vector3d dx = x1 - x2_a;
			if (dx.dot(N2) < 0.0)
				continue;
			// Project x1 onto the triangle
			Vector3d x2 = x1 - dx.dot(N2) * N2;
			dx = x2 - x1;

			if (dx.norm() > 5 * thres) continue;
			// Compute barycentric coords of x1 wrt tri2
			double u, v, w;
			Vector3d t0 = x2_b - x2_a;
			Vector3d t1 = x2_c - x2_a;
			Vector3d t2 = x1 - x2_a;
			double dividend = t0.dot(t0) * t1.dot(t1) - t0.dot(t1) * t0.dot(t1);
			v = (t1.dot(t1) * t2.dot(t0) - t0.dot(t1) * t2.dot(t1)) / dividend;
			w = (t0.dot(t0) * t2.dot(t1) - t0.dot(t1) * t2.dot(t0)) / dividend;
			u = 1 - v - w;
			if (u < 0 || u > 1 || v < 0 || v > 1 || w < 0 || w > 1) continue;
			// Create contact object
			Collision* c_temp = new Collision(dx.norm(), N1, N2, x1, x2, 1, 3, Vector3i(i, -1, -1), one_f2,
				Vector3d(1.0, 0.0, 0.0), Vector3d(u, v, w), -1, j);
			c_temp->edge1.push_back(ve_1(0, i));
			c_temp->edge1.push_back(ve_1(1, i));
			c_temp->edge1.push_back(ve_1(2, i));
			// Is this the closest one so far?
			if (!c || c->dist > dx.norm()) c = c_temp;
		}
		if (c) colls.push_back(c);
	}

	MatrixXd aabb_1_F(6, f1.cols());

	for (int k = 0; k < f1.cols(); k++) {
		Vector3i f = f1.col(k);
		Vector3d a = v1.block<3, 1>(0, f(0));
		Vector3d b = v1.block<3, 1>(0, f(1));
		Vector3d c = v1.block<3, 1>(0, f(2));
		aabb_1_F.block<3, 1>(0, k) = a;
		aabb_1_F.block<3, 1>(3, k) = a;
		for (int i = 0; i < 3; i++) {
			aabb_1_F(i, k) = min(b(i), aabb_1_F(i, k));
			aabb_1_F(i, k) = min(c(i), aabb_1_F(i, k));
			aabb_1_F(i + 3, k) = max(b(i), aabb_1_F(i + 3, k));
			aabb_1_F(i + 3, k) = max(c(i), aabb_1_F(i + 3, k));
		}
	}

	MatrixXd aabb_2_E(6, e2.size());

	for (int k = 0; k < e2.size(); k++) {
		Vector4i e = e2[k]->vs;
		Vector3d t0 = v2.block<3, 1>(0, e(0));
		Vector3d t1 = v2.block<3, 1>(0, e(1));
		aabb_2_E.block<3, 1>(0, k) = t0;
		aabb_2_E.block<3, 1>(3, k) = t0;
		for (int i = 0; i < 3; i++) {
			aabb_2_E(i, k) = min(t1(i), aabb_2_E(i, k));
			aabb_2_E(i + 3, k) = max(t1(i), aabb_2_E(i + 3, k));
		}
	}

	// 1: Edge
	// 2: Edge
	for (int i = 0; i < e2.size(); i++) {
		auto one_e2 = e2[i];
		const Vector3d &x2_a = v2.block<3, 1>(0, one_e2->vs(0));
		const Vector3d &x2_b = v2.block<3, 1>(0, one_e2->vs(1));
		Vector3d dx2 = x2_b - x2_a;
		double len2 = dx2.norm();
		Vector3d N2 = (one_e2->ns[0] + one_e2->ns[1]).normalized();
		Matrix<double, 6, 1> aabbE2k = aabb_2_E.col(i);
		for (int j = 0; j < e1.size(); j++) {
			Edge_* one_e1 = e1[j];
			if (one_e1->angle < M_PI / 6.0) continue;
			
			// AABB test: 
			if (!check_AABB(aabb_1_F.col(one_e1->faces(0)), aabbE2k) && 
				!check_AABB(aabb_1_F.col(one_e1->faces(1)), aabbE2k)) {
				continue;
			}

			// x1_a and x1_b are the vertices of the edge.
			const Vector3d &x1_a = v1.block<3, 1>(0, one_e1->vs(0));
			const Vector3d &x1_b = v1.block<3, 1>(0, one_e1->vs(1));
			Vector3d dx1 = x1_b - x1_a;
			double len1 = dx1.norm();
			Vector3d tan1 = dx1 / len1;

			double t = acos(tan1.dot(N2));
			if (fabs(t) < M_PI / 90.0 || fabs(M_PI - t) < M_PI / 90.0) {
				continue;
			}
			// Are the two edges parallel?
			t = acos(tan1.dot(dx2) / len2);
			if (fabs(t) < M_PI / 90.0 || fabs(M_PI - t) < M_PI / 90.0) {
				continue;
			}
			Vector3d N = dx1.cross(dx2).normalized();
			// The two triangles of the edge are (a,b,c) and (b,a,d), and n1_c 
			// and n1_d are the two triangle normals
			const Vector3d &x1_c = v1.block<3, 1>(0, one_e1->vs(2));
			const Vector3d &x1d = v1.block<3, 1>(0, one_e1->vs(3));
			const Vector3d &n1_c = one_e1->ns[0];
			const Vector3d &n1_d = one_e1->ns[1];
			// Make the computed normal point along the edge normal.
			Vector3d N1 = (n1_c + n1_d).normalized();
			if (N.dot(N1) < 0.0) N = -N;

			double angleCD = acos(n1_c.dot(n1_d)); 
			double angleCN = acos(n1_c.dot(N)); 
			if (angleCD < 0.0) {
				angleCD = -angleCD;
				angleCN = -angleCN;
			}
			// angleCN must be between 0 and angleCD
			if (angleCN < -M_PI / 90.0 || angleCN - angleCD > M_PI / 90.0) {
				continue;
			}
			// Where does the line intersect the box?
			double u2_c, u2_d, unused1, unused2;
			
			int i2_c = intersection(x2_a, dx2, x1_a, x1_b, x1_c, u2_c);
			int i2_d = intersection(x2_a, dx2, x1_b, x1_a, x1d, u2_d);

			i2_c = i2_c && (0.0 <= u2_c && u2_c <= 1.0);
			i2_d = i2_d && (0.0 <= u2_d && u2_d <= 1.0);

			Vector3d B2_B1 = x2_b - x2_a;
			Vector3d A1_B1 = x1_a - x2_a;
			Vector3d A2_A1 = x1_b - x1_a;
			Vector3d A2A1_B2B1 = A2_A1.cross(B2_B1);
			double u1 = B2_B1.cross(A1_B1).dot(A2A1_B2B1) / A2_A1.cross(B2_B1).dot(A2A1_B2B1);
			double u2 = A2_A1.cross(A1_B1).dot(A2A1_B2B1) / A2_A1.cross(B2_B1).dot(A2A1_B2B1);

			if (u1 < -thres / len1 || u1 > 1.0 + thres / len1 ||
				u2 < -thres / len2 || u2 > 1.0 + thres / len2) {
				continue;
			}
			Vector3d x1, x2, dx;
			if (i2_c && i2_d) {
				// The cloth edge intersects both box triangles
				x1 = (1.0 - u1)*x1_a + u1 * x1_b;
				x2 = (1.0 - u2)*x2_a + u2 * x2_b;
			}
			else if (i2_c && !i2_d) {
				dx = x2_a - x1_a;
				if (dx.dot(n1_c) < 0.0) u2 = max(0.0, min(u2_c, u2));
				else u2 = max(u2_c, min(1.0, u2));
				x2 = (1.0 - u2)*x2_a + u2 * x2_b;
				// Recompute u1 based on the new x2.

				u1 = (x2 - x1_a).dot(x1_b - x1_a) / (x1_b - x1_a).dot(x1_b - x1_a);

				x1 = (1.0 - u1)*x1_a + u1 * x1_b;
			}
			else if (!i2_c && i2_d) {
				
				dx = x2_a - x1_b;
				if (dx.dot(n1_d) < 0.0) u2 = max(0.0, min(u2_d, u2));
				else u2 = max(u2_d, min(1.0, u2));
				x2 = (1.0 - u2)*x2_a + u2 * x2_b;
				// Recompute u1 based on the new x2.

				u1 = (x2 - x1_a).dot(x1_b - x1_a) / (x1_b - x1_a).dot(x1_b - x1_a);

				x1 = (1.0 - u1)*x1_a + u1 * x1_b;
			}
			else continue;
			if (u1 < -thres / len1 || u1 > 1.0 + thres / len1 || 
				u2 < -thres / len2 || u2 > 1.0 + thres / len2) {
				continue;
			}
			// At this point, x2 is the collision point on the cloth and x1 is
			// the collision point on the box.
			dx = x2 - x1;
			if (dx.dot(dx) > 4 * thres * thres) continue;

			Collision* c = new Collision(dx.norm(), N1, N, x1, x2, 2, 2, Vector3i(one_e1->vs(0), one_e1->vs(1), -1),
				Vector3i(one_e2->vs(0), one_e2->vs(1), -1), Vector3d(1.0 - u1, u1, 0.0), Vector3d(1.0 - u2, u2, 0.0), -1, -1);
			
			c->edge1.push_back(j);
			c->edge2 = i;
			c->edgeDir = tan1;
			colls.push_back(c);
		}
	}

	// Corners
	for (int i = 0; i < 8; i++) {
		const Vector3d &x1 = v1.block<3, 1>(0, i);
		// Find the closest collision to the box corner
		int kmin = -1;
		double dmin = INFINITY;
		for (int j = 0; j < colls.size(); j++) {
			if (colls[j]->count1 == 1 && colls[j]->verts1(0) == i) {
				const Vector3d &x2 = colls[j]->pos2;
				Vector3d dx = x2 - x1;
				double d = dx.dot(dx);
				if (d < dmin) {
					kmin = j;
					dmin = d;
				}
			}
		}
		if (kmin != -1) {
			// Make a list of colls to delete
			vector<int> dlist;
			for (int j = 0; j < colls.size(); j++) {
				if (colls[j]->count1 == 1 && colls[j]->verts1(0) == i && j != kmin) {
					dlist.push_back(j);
				}
			}
			// Delete colls in reverse order
			for (int temp : dlist) {
				colls[temp] = colls.back();
				colls.pop_back();
			}
		}
	}

	// Perturb
	for (Collision* one_cl : colls) {
		one_cl->pos1_ = one_cl->pos1 - 0.1 * thres * one_cl->nor1;
	}
}

void ForPoint(
	vector<Collision* > &colls,
	double thres,
	const Eigen::MatrixXd &v1,
	const Eigen::MatrixXd &n1,
	const Eigen::MatrixXd &verts2_,
	const Eigen::MatrixXi &f2,
	bool EOL)
{
	MatrixXd v2 = verts2_;

	// Perturb cloth verts
	std::random_device device;
	std::mt19937 mt;
	std::uniform_real_distribution<> dis(-1.0, 1.0);
	mt.seed(1);
	for (int i = 0; i < v2.cols(); i++) {
		Vector3d r;
		r(0) = dis(mt) * thres * 1e-3;
		r(1) = dis(mt) * thres * 1e-3;
		r(2) = dis(mt) * thres * 1e-3;
		v2.block<3, 1>(0, i) += r;
	}

	// Precompute the face normals for the cloth
	MatrixXd faceNors2 = create_face_n(f2, v2);

	// Build AABBs
	Matrix<double, 6, 1> aabb_2;
	AABB(aabb_2, v2);

	if (EOL) {
		for (int i = 0; i < v2.cols(); i++) {
			Vector3d x2 = v2.block<3, 1>(0, i);
			Collision* c = NULL;
			for (int j = 0; j < v1.cols(); j++) {
				Vector3d x1 = v1.block<3, 1>(0, j);
				Vector3d dx = x2 - x1;

				if (dx.norm() < thres) {
					Vector3d N1 = n1.block<3, 1>(0, j); // vertex normal
					// Create contact object
					Collision* c_temp = new Collision();
					c_temp->dist = dx.norm();
					c_temp->nor1 = N1;
					c_temp->nor2 = N1;
					c_temp->pos1 = x1;
					c_temp->pos2 = x2;
					c_temp->count1 = 3;
					c_temp->count2 = 1;
					c_temp->verts1 << j, -1, -1;
					c_temp->verts2 << i, -1, -1;
					c_temp->weights1 << 1.0, 0.0, 0.0;
					c_temp->weights2 << 1.0, 0.0, 0.0;
					c_temp->tri1 = -1;
					c_temp->tri2 = -1;
					// Is this the closest one so far?
					if (!c || c->dist > dx.norm()) c = c_temp;
				}
			}
			if (c) colls.push_back(c);
		}
	}

	// Vertex1-Triangle2
	for (int i = 0; i < v1.cols(); i++) {
		const Vector3d &x1 = v1.block<3, 1>(0, i);
		// AABB test: check Vertex1 against Body2
		Matrix<double, 6, 1> aabb_1_V;
		aabb_1_V.segment<3>(0) = x1;
		aabb_1_V.segment<3>(3) = x1;
		if (!check_AABB(aabb_1_V, aabb_2)) {
			continue;
		}
		Collision* c = NULL;
		Vector3d N1 = n1.block<3, 1>(0, i); // vertex normal
		for (int j = 0; j < f2.cols(); j++) {
			const Vector3i &one_f2 = f2.col(j);
			const Vector3d &x2_a = v2.block<3, 1>(0, one_f2(0));
			const Vector3d &x2_b = v2.block<3, 1>(0, one_f2(1));
			const Vector3d &x2_c = v2.block<3, 1>(0, one_f2(2));
			Vector3d N2 = faceNors2.col(j);
			// Make sure the triangle normal points outward wrt the box.
			if (N1.dot(N2) < 0.0) {
				N2 = -N2;
			}
			// Is x1 on the correct side?
			Vector3d dx = x1 - x2_a;

			if (dx.dot(N2) < 0.0) continue;
			// Project x1 onto the triangle
			Vector3d x2 = x1 - dx.dot(N2) * N2;
			dx = x2 - x1;
			if (dx.norm() > 5 * thres) continue;
			// Barycentric coord
			double u, v, w;
			Vector3d t0 = x2_b - x2_a;
			Vector3d t1 = x2_c - x2_a;
			Vector3d t2 = x1 - x2_a;
			double dividend = t0.dot(t0) * t1.dot(t1) - t0.dot(t1) * t0.dot(t1);
			v = (t1.dot(t1) * t2.dot(t0) - t0.dot(t1) * t2.dot(t1)) / dividend;
			w = (t0.dot(t0) * t2.dot(t1) - t0.dot(t1) * t2.dot(t0)) / dividend;
			u = 1 - v - w;
			if (u < 0 || u > 1 || v < 0 || v > 1 || w < 0 || w > 1) continue;
			// Create contact object
			Collision* c_temp = new Collision;
			c_temp->dist = dx.norm();
			c_temp->nor1 = N1;
			c_temp->nor2 = N2;
			c_temp->pos1 = x1;
			c_temp->pos2 = x2;
			c_temp->count1 = 1;
			c_temp->count2 = 3;
			c_temp->verts1 << i, -1, -1;
			c_temp->verts2 = one_f2;
			c_temp->weights1 << 1.0, 0.0, 0.0;
			c_temp->weights2 << u, v, w;
			c_temp->tri1 = -1;
			c_temp->tri2 = j;
			// Is this the closest one so far?
			if (!c || c->dist > dx.norm()) c = c_temp;
		}
		if (c) colls.push_back(c);
	}

	// Perturb
	for (Collision* one_cl : colls) {
		one_cl->pos1_ = one_cl->pos1 - 0.1 * thres * one_cl->nor1;
	}
}


void Find_Collision(const Mesh& mesh, const Obstacle* obs, std::vector<Collision*> &colls)
{
	MatrixXd v2(3, mesh.nodes.size());
	MatrixXi f2(3, mesh.faces.size());
	VectorXi EoLs;
	EoLs.resize(mesh.nodes.size());

	for (int i = 0; i < mesh.nodes.size(); i++)
	{
		v2.col(i) = Vector3d(mesh.nodes[i]->x[0], mesh.nodes[i]->x[1], mesh.nodes[i]->x[2]);
		if (mesh.nodes[i]->EoL) EoLs(i) = 1;
	}
	for (int i = 0; i < mesh.faces.size(); i++)
	{
		f2.col(i) = Vector3i(mesh.faces[i]->v[0]->node->index, mesh.faces[i]->v[1]->node->index, mesh.faces[i]->v[2]->node->index);
	}

	ForPoint(colls, obs->collision_threshold, obs->point_cloud->points_position, obs->point_cloud->points_normal, v2, f2, true);

	int c = colls.size();
	for (int b = 0; b < obs->box_num; b++)
	{
		vector<Collision*> colls_list;
		Matrix4d E1;
		E1.setIdentity();
		E1.block<3, 1>(0, 3) = obs->boxes[b]->position;
		ForBox(colls_list, obs->collision_threshold, obs->boxes[b]->scale, E1, v2, f2, EoLs, false);
		colls.insert(colls.end(), colls_list.begin(), colls_list.end());
		

		for (c; c < colls.size(); c++) {
			if (colls[c]->count1 == 1 && colls[c]->count2 == 3) {
				colls[c]->verts1(0) = obs->point_cloud->point_num + (b* obs->boxes[b]->point_num) + (b* obs->boxes[b]->edge_num) + colls[c]->verts1(0);
			}
			for (int e = 0; e < colls[c]->edge1.size(); e++)
				colls[c]->edge1[e] = obs->point_cloud->point_num + (b* obs->boxes[b]->point_num) + (b* obs->boxes[b]->edge_num) + (obs->boxes[b]->point_num + colls[c]->edge1[e]);
		}
	}


}

void Find_Collision_(const Mesh& mesh, const Obstacle* obs, std::vector<Collision*> &colls)
{
	MatrixXd v2(3, mesh.nodes.size());
	MatrixXi f2(3, mesh.faces.size());
	VectorXi EoLs;
	EoLs.resize(mesh.nodes.size());

	for (int i = 0; i < mesh.nodes.size(); i++)
	{
		v2.col(i) = Vector3d(mesh.nodes[i]->x[0], mesh.nodes[i]->x[1], mesh.nodes[i]->x[2]);
		if (mesh.nodes[i]->EoL) EoLs(i) = 1;
	}
	for (int i = 0; i < mesh.faces.size(); i++)
	{
		f2.col(i) = Vector3i(mesh.faces[i]->v[0]->node->index, mesh.faces[i]->v[1]->node->index, mesh.faces[i]->v[2]->node->index);
	}

	ForPoint(colls, obs->collision_threshold, obs->point_cloud->points_position, obs->point_cloud->points_normal, v2, f2, false);

	for (int b = 0; b < obs->box_num; b++)
	{
		vector<Collision*> colls_list;
		Matrix4d E1;
		E1.setIdentity();
		E1.block<3, 1>(0, 3) = obs->boxes[b]->position;
		ForBox(colls_list, obs->collision_threshold, obs->boxes[b]->scale, E1, v2, f2, EoLs, false);
		colls.insert(colls.end(), colls_list.begin(), colls_list.end());
	}
}