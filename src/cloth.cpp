#include "cloth.h"

#include "obstacle.h"
#include "constraint.h"
#include "force.h"
#include "utility.h"
#include "preprocessor.h"

#include "ArcSim/mesh.hpp"
#include "ArcSim/io.hpp"
#include "ArcSim/geometry.hpp"
#include "ArcSim/referenceshape.hpp"

using namespace Eigen;
using namespace std;

Cloth::Cloth(Vector2i res, VectorXd& xX00, VectorXd& xX01, VectorXd& xX10, VectorXd& xX11)
{
	constraint = new Constraint;
	force = new Force;
	preprocessor = new Preprocessor;

	//Build Cloth Mesh
	Vector3d x00, x01, x10, x11;
	Vector3d X00 = Vector3d::Zero(), X01 = Vector3d::Zero(), X10 = Vector3d::Zero(), X11 = Vector3d::Zero();

	x00 = xX00.segment<3>(0);
	x01 = xX01.segment<3>(0);
	x10 = xX10.segment<3>(0);
	x11 = xX11.segment<3>(0);
	X00.segment<2>(0) = xX00.segment<2>(3);
	X01.segment<2>(0) = xX01.segment<2>(3);
	X10.segment<2>(0) = xX10.segment<2>(3);
	X11.segment<2>(0) = xX11.segment<2>(3);

	// Set boundary values
	boundary.resize(3, 4);
	boundary.block<3, 1>(0, 0) = X00;
	boundary.block<3, 1>(0, 1) = X01;
	boundary.block<3, 1>(0, 2) = X11;
	boundary.block<3, 1>(0, 3) = X10;

	for (int i = 0; i < res(0); ++i) 
	{
		double u = i / (res(0) - 1.0);
		Vector3d x0 = (1 - u)*x00 + u * x10;
		Vector3d x1 = (1 - u)*x01 + u * x11;
		Vector3d X0 = (1 - u)*X00 + u * X10;
		Vector3d X1 = (1 - u)*X01 + u * X11;
		for (int j = 0; j < res(1); ++j) 
		{
			double v = j / (res(1) - 1.0);
			Vector3d x = (1 - v)*x0 + v * x1;
			Vector3d X = (1 - v)*X0 + v * X1;
			mesh.add(new Vert(e2v(X), Vec3(0)));
			mesh.add(new Node(e2v(x), e2v(x), Vec3(0), 0, 0, false));
			connect(mesh.verts.back(), mesh.nodes.back());
		}
	}

	for (int i = 0; i < res(0) - 1; ++i) 
	{
		for (int j = 0; j < res(1) - 1; ++j) 
		{
			int k0 = (i * res(1)) + j; // upper right index
			Vector3d x = (v2e(mesh.nodes[k0]->x) + v2e(mesh.nodes[k0 + 1]->x) + v2e(mesh.nodes[k0 + res(1) + 1]->x) + v2e(mesh.nodes[k0 + res(1)]->x)) / 4;
			Vector3d X = (v2e(mesh.verts[k0]->u) + v2e(mesh.verts[k0 + 1]->u) + v2e(mesh.verts[k0 + res(1) + 1]->u) + v2e(mesh.verts[k0 + res(1)]->u)) / 4;
			mesh.add(new Vert(e2v(X), Vec3(0)));
			mesh.add(new Node(e2v(x), e2v(x), Vec3(0), 0, 0, false));
			connect(mesh.verts.back(), mesh.nodes.back());
		}
	}
	int center_index_cnt = 0;
	for (int i = 0; i < res(0) - 1; i++)
	{
		for (int j = 0; j < res(1) - 1; j++)
		{
			int k0 = (i * res(1)) + j; // upper right index
			int kc0 = ((res(0) * res(1)) + center_index_cnt); // center index
			center_index_cnt++;
			vector<Vert*> verts1;
			verts1.push_back(mesh.verts[k0]);
			verts1.push_back(mesh.verts[k0 + 1]);
			verts1.push_back(mesh.verts[kc0]);
			vector<Face*> faces1 = triangulateARC(verts1);
			for (int f = 0; f < faces1.size(); f++)
				mesh.add(faces1[f]);
			vector<Vert*> verts2;
			verts2.push_back(mesh.verts[k0 + 1]);
			verts2.push_back(mesh.verts[k0 + res(1) + 1]);
			verts2.push_back(mesh.verts[kc0]);
			vector<Face*> faces2 = triangulateARC(verts2);
			for (int f = 0; f < faces2.size(); f++)
				mesh.add(faces2[f]);
			vector<Vert*> verts3;
			verts3.push_back(mesh.verts[k0 + res(1) + 1]);
			verts3.push_back(mesh.verts[k0 + res(1)]);
			verts3.push_back(mesh.verts[kc0]);
			vector<Face*> faces3 = triangulateARC(verts3);
			for (int f = 0; f < faces3.size(); f++)
				mesh.add(faces3[f]);
			vector<Vert*> verts4;
			verts4.push_back(mesh.verts[k0 + res(1)]);
			verts4.push_back(mesh.verts[k0]);
			verts4.push_back(mesh.verts[kc0]);
			vector<Face*> faces4 = triangulateARC(verts4);
			for (int f = 0; f < faces4.size(); f++)
				mesh.add(faces4[f]);
		}
	}
	for (int i = 0; i < mesh.faces.size(); i++) 
	{
		mesh.faces[i]->material = &material;
	}

	mark_nodes_to_preserve(mesh);
	compute_ms_data(mesh);

	v.resize(mesh.nodes.size() * 3);

	mesh.parent = this;
	mesh.ref = new ReferenceLinear(mesh);

	material.density = 0.05;
	material.e = 50.0;
	material.nu = 0.01;
	material.beta = 1e-5;
	material.damping1 = 0.0;
	material.damping2 = 1.0;

	remeshing.refine_angle = 0.3;
	remeshing.refine_compression = 0.005;
	remeshing.refine_velocity = 0.5;
	remeshing.size_min = 320e-3;
	remeshing.size_max = 350e-3;
	remeshing.aspect_min = 0.2;

	buffer = NULL;
	createClothBuffer();
}


void Cloth::saveOldMesh()
{
	delete_mesh(old_mesh);
	old_mesh = deep_copy(mesh);
}

/*
VelocityTransfer: for those old EOL points which is caused by the collision in the formaer frame, just put the Eulerian velocity into lag velocity, for other points, just average case by case
1. found in old mesh
 1.1 wasEOl: old EoL point: just transfer Eulerian velocity into lag velocity to let them move away
 1.2 not wasEOL: normal point, don't care
2. not found in old mesh
 2.1 not EoL: new lag point ,just average the traignle's velocity by barycentric coordinates
 2.2 EoL:
    2.2.1 newEolFromSplit: the point is derived from 2 known EoL point, so just average the velocity at the ends
	2.2.2 not newEolFromSplit: the point caused by normal collision, just average the world velocity and store the velocity into Eulerian velocity as much as possible using KKT optimization
*/

void Cloth::velocityTransfer()
{
	v.resize(mesh.nodes.size() * 3 + mesh.EoL_Count * 2);

	v.setZero();

	// Loop through all of our nodes and update their velocities
	for (int n = 0; n < mesh.nodes.size(); n++) 
	{
		Node* node = mesh.nodes[n];
		Vert* vert = node->verts[0];

		bool found = false;
		double how_close = 1e-6;
		int closest = -1;
		for (int j = 0; j < old_mesh.nodes.size(); j++) 
		{
			double dist = unsigned_vv_distance(vert->u, old_mesh.nodes[j]->verts[0]->u);
			if (unsigned_vv_distance(vert->u, old_mesh.nodes[j]->verts[0]->u) < how_close) 
			{
				how_close = dist;
				closest = j;
				found = true;
				break;
			}
		}

		// If the node existed before do one of two things
		if (found) 
		{
			// THIS is the key point bro!
			if (node->EoL_state == Node::WasEOL)
			{
				node->EoL_state = Node::IsLAG;
				MatrixXd F = MatrixXd::Zero(3, 2);
				Vector3d nodev = v2e(node->v);
				Vector2d nodeV = v322e(vert->v);
				for (int j = 0; j < vert->adjf.size(); j++) 
				{
					F += incedent_angle(vert, vert->adjf[j]) * computeF(vert->adjf[j]);
				}
				is_seam_or_boundary(node) ? F *= (1 / M_PI) : F *= (1 / (2 * M_PI));
				Vector3d newvL = nodev - F * nodeV;
				v(3 * n) = newvL(0);
				v(3 * n + 1) = newvL(1);
				v(3 * n + 2) = newvL(2);
			}
			// If it already existed, just pull it's old information, the most noraml point, nothing to care
			else
			{
				v(3 * n) = node->v[0];
				v(3 * n + 1) = node->v[1];
				v(3 * n + 2) = node->v[2];
				if (node->EoL) 
				{
					v(mesh.nodes.size() * 3 + node->EoL_index * 2) = vert->v[0];
					v(mesh.nodes.size() * 3 + node->EoL_index * 2 + 1) = vert->v[1];
				}
			}
		}

		//newly inserted point
		else 
		{
			// If its a LAG point we can just use barycentric averaging
			if (!node->EoL) 
			{
				Face* old_face = get_enclosing_face(old_mesh, Vec2(vert->u[0], vert->u[1]));
				Vec3 bary = get_barycentric_coords(Vec2(vert->u[0], vert->u[1]), old_face);
				Vector3d vwA = v2e(old_face->v[0]->node->v);
				Vector3d vwB = v2e(old_face->v[1]->node->v);
				Vector3d vwC = v2e(old_face->v[2]->node->v);

				if (old_face->v[0]->node->EoL) vwA += -computeF(old_face) * v322e(old_face->v[0]->v);
				if (old_face->v[1]->node->EoL) vwB += -computeF(old_face) * v322e(old_face->v[1]->v);
				if (old_face->v[2]->node->EoL) vwC += -computeF(old_face) * v322e(old_face->v[2]->v);
				Vector3d v_new_world = bary[0] * vwA + bary[1] * vwB + bary[2] * vwC;
				node->v = e2v(v_new_world);
				v(3 * n) = v_new_world(0);
				v(3 * n + 1) = v_new_world(1);
				v(3 * n + 2) = v_new_world(2);
			}
			else 
			{
				if (node->EoL_state == Node::NewEOLFromSplit) // NewEOLFromSplit means this point is derived from 2 known eulerian points, so just average
				{
					Vector3d newvL = Vector3d::Zero();
					Vector3d newvE = Vector3d::Zero();
					for (int e = 0; e < node->adje.size(); e++) 
					{
						if (node->adje[e]->preserve) 
						{
							Node* adjn = other_node(node->adje[e], node);
							newvL += v2e(adjn->v);
							newvE += v2e(adjn->verts[0]->v);
						}
					}
					newvL /= 2.0;
					newvE /= 2.0; //average the two ends of the point

					// We put in the new averaged velocities
					v(mesh.nodes.size() * 3 + node->EoL_index * 2) = newvE(0);
					v(mesh.nodes.size() * 3 + node->EoL_index * 2 + 1) = newvE(1);
					v(3 * n) = newvL(0);
					v(3 * n + 1) = newvL(1);
					v(3 * n + 2) = newvL(2);
				}
				else 
				{
					Face* old_face = get_enclosing_face(old_mesh, Vec2(vert->u[0], vert->u[1]));
					Matrix2d ftf = computeF(old_face).transpose() * computeF(old_face);
					Vector2d dtv = -computeF(old_face).transpose() * v2e(node->v); // the average v has been done in the function split_face()

					MatrixXd KKTl = MatrixXd::Zero(4, 4);
					KKTl.block(0, 0, 2, 2) = ftf;
					KKTl.block<2, 2>(0, 2) = Matrix2d::Identity();
					KKTl.block<2, 2>(2, 0) = Matrix2d::Identity();
					VectorXd KKTr(4);
					KKTr << dtv, VectorXd::Zero(2);
					ConjugateGradient<MatrixXd, Lower | Upper> cg;
					cg.compute(KKTl);
					VectorXd newvE = cg.solve(KKTr);

					// Eulerian component
					vert->v = Vec3(newvE(0), newvE(1), 0.0);
					v(mesh.nodes.size() * 3 + node->EoL_index * 2) = newvE(0);
					v(mesh.nodes.size() * 3 + node->EoL_index * 2 + 1) = newvE(1);

					// Lagrangian component
					Vector3d newvL = v2e(node->v) + -computeF(old_face) * newvE.segment<2>(0);
					node->v = e2v(newvL);
					v(3 * n) = newvL(0);
					v(3 * n + 1) = newvL(1);
					v(3 * n + 2) = newvL(2);
				}
			}
		}
	}
}

void Cloth::solve(double h)
{
	VectorXd b = -(force->M * v + h * force->f);

	MatrixXd M(force-> K);
	MatrixXd I(M.rows(), M.cols());
	double u = 1.0;
	I.setIdentity();
	while (!PositiveDefinite(M))
	{
		u *= 2.0;
		M = M + u * I;
	}

	force->K = M.sparseView();
	bool success = mosekSolve(force->K, b, constraint->Aeq, constraint->beq, constraint->Aineq, constraint->bineq, v);
}

void Cloth::step(Obstacle* obstacle, const Vector3d& gravity, double h, double t)
{
	velocityTransfer();
	force->computeForce(mesh, material, gravity, h);
	constraint->applyConstraint(mesh, obstacle, h, t);
	solve(h);

	for (int n = 0; n < mesh.nodes.size(); n++)
	{
		Node* node = mesh.nodes[n];
		Vert* vert = node->verts[0];

		node->v[0] = v(n * 3);
		node->v[1] = v(n * 3 + 1);
		node->v[2] = v(n * 3 + 2);
		node->x = node->x + h * node->v;
		if (node->EoL)
		{
			vert->v[0] = v(mesh.nodes.size() * 3 + node->EoL_index * 2);
			vert->v[1] = v(mesh.nodes.size() * 3 + node->EoL_index * 2 + 1);
			vert->u[0] = vert->u[0] + h * vert->v[0];
			vert->u[1] = vert->u[1] + h * vert->v[1];;
		}
	}

	createClothBuffer();
}

bool Cloth::mosekSolve(const SparseMatrix<double>& MDK, const VectorXd& b,
	const SparseMatrix<double>& Aeq, const VectorXd& beq,
	const SparseMatrix<double>& Aineq, const VectorXd& bineq,
	VectorXd& v)
{
	QuadProgMosek *program = new QuadProgMosek();
	double inf = numeric_limits<double>::infinity();

	VectorXd xl;
	VectorXd xu;
	xl.setConstant(b.size(), -inf);
	xu.setConstant(b.size(), inf);

	program->setNumberOfVariables(b.size());
	program->setNumberOfEqualities(beq.size());
	program->setNumberOfInequalities(bineq.size());

	program->setObjectiveMatrix(MDK);
	program->setObjectiveVector(b);

	program->setInequalityMatrix(Aineq);
	program->setInequalityVector(bineq);

	program->setEqualityMatrix(Aeq);
	program->setEqualityVector(beq);

	bool success = program->solve();

	v = program->getPrimalSolution();

	return success;
}

void Cloth::createClothBuffer()
{
	delete buffer;
	delete mesh_buffer;
	buffer = (float*)malloc(sizeof(float)*mesh.faces.size() * 3 * 6);
	mesh_buffer = (float*)malloc(sizeof(float)*mesh.faces.size() * 6 * 3);
	

	for (int i = 0; i < mesh.faces.size(); i++)
	{
		buffer[i * 18] = mesh.faces[i]->v[0]->node->x[0];
		buffer[i * 18+1] = mesh.faces[i]->v[0]->node->x[1];
		buffer[i * 18+2] = mesh.faces[i]->v[0]->node->x[2];
		buffer[i * 18+3] = mesh.faces[i]->n[0];
		buffer[i * 18+4] = mesh.faces[i]->n[1];
		buffer[i * 18+5] = mesh.faces[i]->n[2];

		buffer[i * 18 + 6] = mesh.faces[i]->v[1]->node->x[0];
		buffer[i * 18 + 7] = mesh.faces[i]->v[1]->node->x[1];
		buffer[i * 18 + 8] = mesh.faces[i]->v[1]->node->x[2];
		buffer[i * 18 + 9] = mesh.faces[i]->n[0];
		buffer[i * 18 + 10] = mesh.faces[i]->n[1];
		buffer[i * 18 + 11] = mesh.faces[i]->n[2];

		buffer[i * 18 + 12] = mesh.faces[i]->v[2]->node->x[0];
		buffer[i * 18 + 13] = mesh.faces[i]->v[2]->node->x[1];
		buffer[i * 18 + 14] = mesh.faces[i]->v[2]->node->x[2];
		buffer[i * 18 + 15] = mesh.faces[i]->n[0];
		buffer[i * 18 + 16] = mesh.faces[i]->n[1];
		buffer[i * 18 + 17] = mesh.faces[i]->n[2];
	}

	for (int i = 0; i < mesh.faces.size(); i++)
	{
		mesh_buffer[i * 18]   =	 mesh.faces[i]->v[0]->u[0]*0.55 + 0.45;
		mesh_buffer[i * 18 + 1] = -mesh.faces[i]->v[0]->u[1]*0.55 - 0.45;
		mesh_buffer[i * 18 + 2] = 0;
		
		mesh_buffer[i * 18 + 3] = mesh.faces[i]->v[1]->u[0]*0.55 + 0.45;
		mesh_buffer[i * 18 + 4] = -mesh.faces[i]->v[1]->u[1]*0.55 - 0.45;
		mesh_buffer[i * 18 + 5] = 0;

		mesh_buffer[i * 18 + 6] = mesh.faces[i]->v[0]->u[0]*0.55 + 0.45;
		mesh_buffer[i * 18 + 7] = -mesh.faces[i]->v[0]->u[1]*0.55 - 0.45;
		mesh_buffer[i * 18 + 8] = 0;
		
		mesh_buffer[i * 18 + 9]  = mesh.faces[i]->v[2]->u[0]*0.55 + 0.45;
		mesh_buffer[i * 18 + 10] = -mesh.faces[i]->v[2]->u[1]*0.55 - 0.45;
		mesh_buffer[i * 18 + 11] = 0;

		mesh_buffer[i * 18 + 12] = mesh.faces[i]->v[1]->u[0]*0.55 + 0.45;
		mesh_buffer[i * 18 + 13] = -mesh.faces[i]->v[1]->u[1]*0.55 - 0.45;
		mesh_buffer[i * 18 + 14] = 0;
		
		mesh_buffer[i * 18 + 15] = mesh.faces[i]->v[2]->u[0]*0.55 + 0.45;
		mesh_buffer[i * 18 + 16] = -mesh.faces[i]->v[2]->u[1]*0.55 - 0.45;
		mesh_buffer[i * 18 + 17] = 0;
	}


	buffer_size = sizeof(float)*mesh.faces.size() * 6 * 3;
	mesh_buffer_size = sizeof(float)*mesh.faces.size() * 3 * 6;

}