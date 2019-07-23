#include "constraint.h"

#include "obstacle.h"
#include "utility.h"
#include "Collisions.h"

#include "ArcSim/util.hpp"
#include "ArcSim/geometry.hpp"

using namespace Eigen;
using namespace std;

Constraint::Constraint()
{
	occur_collision = false;
}

void Constraint::updateConstraint(const Obstacle* obstacle)
{
	int total_size = obstacle->point_cloud->point_num;

	for (int b = 0; b < obstacle->boxes.size(); b++)
	{
		total_size += (obstacle->boxes[b]->point_num + obstacle->boxes[b]->edge_num);
	}

	constraint_matrix.resize(9, total_size);

	for (int p = 0; p < obstacle->point_cloud->point_num; p++) 
	{
		constraint_matrix.block(0, p, 3, 1) = obstacle->point_cloud->points_normal.col(p);
	}

	for (int b = 0; b < obstacle->boxes.size(); b++) 
	{
		for (int p = 0; p < obstacle->boxes[b]->point_num; p++) 
		{
			int index = obstacle->point_cloud->point_num + (b* obstacle->boxes[b]->point_num) + (b* obstacle->boxes[b]->edge_num) + p;
			Vector3d corner_nor = Vector3d::Zero();
			for (int i = 0; i < 3; i++) 
			{
				for (int j = 0; j < 2; j++) 
				{
					corner_nor += obstacle->boxes[b]->face_norm.col(obstacle->boxes[b]->edge_face(j, obstacle->boxes[b]->vert_edge(i, p)));
				}
			}
			corner_nor /= 6.0;
			constraint_matrix.block(0, index, 3, 1) = corner_nor.normalized();

		}
		for (int e = 0; e < obstacle->boxes[b]->edge_num; e++) 
		{
			int index = obstacle->point_cloud->point_num + (b* obstacle->boxes[b]->point_num) + (b* obstacle->boxes[b]->edge_num) + (obstacle->boxes[b]->point_num + e);
			constraint_matrix.block(0, index, 3, 1) = obstacle->boxes[b]->face_norm.col(obstacle->boxes[b]->edge_face(0, e));
			constraint_matrix.block(3, index, 3, 1) = obstacle->boxes[b]->face_norm.col(obstacle->boxes[b]->edge_face(1, e));
			constraint_matrix.block(6, index, 3, 1) = obstacle->boxes[b]->face_norm.col(obstacle->boxes[b]->edge_direction(e));
		}
	}
}

void Constraint::applyConstraint(const Mesh&mesh, const Obstacle* obstacle, double h, double t)
{
	updateConstraint(obstacle);

	occur_collision = false;

	vector<T> _Aeq;
	vector<T> _Aineq;
	vector< pair<int, double> > _beq;
	vector< pair<int, double> > _bineq;

	int eqsize = 0;
	int ineqsize = 0;

	for (int n = 0; n < mesh.nodes.size(); n++)
	{
		if (mesh.nodes[n]->EoL)
		{
			if (mesh.nodes[n]->cornerID >= 0) // cloth's point has collision with the obstacle
			{
				Vector3d normal = constraint_matrix.block(0, mesh.nodes[n]->cornerID, 3, 1);
				Vector3d ortho1 = Vector3d(0.0, -normal(2), normal(1)); // just to build a frame, maybe can be changed arbitrarily
				Vector3d ortho2 = (ortho1.cross(normal)).normalized();
				ortho2.normalize();

				bool is_flat = true;
				Node* node = mesh.nodes[n];
				Face* face0 = node->verts[0]->adjf[0];
				for (int f = 1; f < node->verts[0]->adjf.size(); f++)
				{
					Face* face1 = node->verts[0]->adjf[f];
					if (get_angle(face0->n, face1->n) > 0.5)
					{
						is_flat = false;
						break;
					}
				}

				// collision  with obstacle's sharp point or the surrounding is is_flat, keep the normal positive and other director equal
				if (is_flat || node->cornerID < obstacle->point_cloud->point_num)
				{

					_Aineq.push_back(T(ineqsize, n * 3, -normal(0)));
					_Aineq.push_back(T(ineqsize, n * 3 + 1, -normal(1)));
					_Aineq.push_back(T(ineqsize, n * 3 + 2, -normal(2)));

					ineqsize++;

					_Aeq.push_back(T(eqsize, n * 3, ortho1(0)));
					_Aeq.push_back(T(eqsize, n * 3 + 1, ortho1(1)));
					_Aeq.push_back(T(eqsize, n * 3 + 2, ortho1(2)));

					eqsize++;

					_Aeq.push_back(T(eqsize, n * 3, ortho2(0)));
					_Aeq.push_back(T(eqsize, n * 3 + 1, ortho2(1)));
					_Aeq.push_back(T(eqsize, n * 3 + 2, ortho2(2)));

					eqsize++;
				}
				//collision with obstacle's box, just ineq constraint can suffice
				else
				{
					_Aineq.push_back(T(ineqsize, n * 3, -normal(0)));
					_Aineq.push_back(T(ineqsize, n * 3 + 1, -normal(1)));
					_Aineq.push_back(T(ineqsize, n * 3 + 2, -normal(2)));

					ineqsize++;

					_Aineq.push_back(T(ineqsize, n * 3, -ortho1(0)));
					_Aineq.push_back(T(ineqsize, n * 3 + 1, -ortho1(1)));
					_Aineq.push_back(T(ineqsize, n * 3 + 2, -ortho1(2)));

					ineqsize++;

					_Aineq.push_back(T(ineqsize, n * 3, -ortho2(0)));
					_Aineq.push_back(T(ineqsize, n * 3 + 1, -ortho2(1)));
					_Aineq.push_back(T(ineqsize, n * 3 + 2, -ortho2(2)));

					ineqsize++;
				}
			}
			// not collision but EOL
			else 
			{
				Node* node = mesh.nodes[n];

				bool is_flat = true;
				Face* face0 = node->verts[0]->adjf[0];
				for (int f = 1; f < node->verts[0]->adjf.size(); f++)
				{
					Face* face1 = node->verts[0]->adjf[f];
					if (get_angle(face0->n, face1->n) > 0.1) 
					{
						is_flat = false;
						break;
					}
				}

				if (is_flat)
				{
					_Aineq.push_back(T(ineqsize, n * 3, -node->n[0]));
					_Aineq.push_back(T(ineqsize, n * 3 + 1, -node->n[1]));
					_Aineq.push_back(T(ineqsize, n * 3 + 2, -node->n[2]));

					ineqsize++;

					Vector3d flatConstraint = v2e(node->n).cross(constraint_matrix.block<3, 1>(6, node->cdEdges[0]));

					_Aeq.push_back(T(eqsize, n * 3, -flatConstraint(0)));
					_Aeq.push_back(T(eqsize, n * 3 + 1, -flatConstraint(1)));
					_Aeq.push_back(T(eqsize, n * 3 + 2, -flatConstraint(2)));

					eqsize++;
				}
				else 
				{
					_Aineq.push_back(T(ineqsize, n * 3, -constraint_matrix(0, node->cdEdges[0])));
					_Aineq.push_back(T(ineqsize, n * 3 + 1, -constraint_matrix(1, node->cdEdges[0])));
					_Aineq.push_back(T(ineqsize, n * 3 + 2, -constraint_matrix(2, node->cdEdges[0])));

					ineqsize++;

					_Aineq.push_back(T(ineqsize, n * 3, -constraint_matrix(3, node->cdEdges[0])));
					_Aineq.push_back(T(ineqsize, n * 3 + 1, -constraint_matrix(4, node->cdEdges[0])));
					_Aineq.push_back(T(ineqsize, n * 3 + 2, -constraint_matrix(5, node->cdEdges[0])));

					ineqsize++;
				}

				// If a boundary, the Eulerian constraint stops it from moving outside
				if (is_seam_or_boundary(node)) 
				{
					Edge* edge=NULL;
					for (int e = 0; e < node->adje.size(); e++) 
					{
						if (is_seam_or_boundary(node->adje[e])) 
						{
							edge = node->adje[e];
							break;
						}
					}
					Node* opp_node = other_node(edge, node);
					// This should be orthogonal to the edge connecting the two nodes
					Vector2d orth_border = Vector2d(node->verts[0]->u[1] - opp_node->verts[0]->u[1], -node->verts[0]->u[0] - opp_node->verts[0]->u[0]).normalized();
					_Aeq.push_back(T(eqsize, mesh.nodes.size() * 3 + node->EoL_index * 2, orth_border(0)));
					_Aeq.push_back(T(eqsize, mesh.nodes.size() * 3 + node->EoL_index * 2 + 1, orth_border(1)));

					eqsize++;
				}
				// If internal, the Eulerian constraint forces tangential motion to realize in the Lagrangian space
				else 
				{
					Vector2d average_tan = Vector2d::Zero();
					for (int e = 0; e < node->adje.size(); e++) 
					{
						if (node->adje[e]->preserve) 
						{
							Edge* edge = node->adje[e];
							if (norm(edge->n[0]->verts[0]->u) > norm(edge->n[1]->verts[0]->u)) 
							{
								average_tan += Vector2d(edge->n[1]->verts[0]->u[0] - edge->n[0]->verts[0]->u[0], edge->n[1]->verts[0]->u[1] - edge->n[0]->verts[0]->u[1]).normalized();
							}
							else 
							{
								average_tan += Vector2d(edge->n[0]->verts[0]->u[0] - edge->n[1]->verts[0]->u[0], edge->n[0]->verts[0]->u[1] - edge->n[1]->verts[0]->u[1]).normalized();
							}
						}
					}
					average_tan.normalize();
					_Aeq.push_back(T(eqsize, mesh.nodes.size() * 3 + node->EoL_index * 2, average_tan(0)));
					_Aeq.push_back(T(eqsize, mesh.nodes.size() * 3 + node->EoL_index * 2 + 1, average_tan(1)));

					eqsize++;
				}
			}
		}
	}

	vector<Collision* > clsLAG;
	Find_Collision_(mesh, obstacle, clsLAG);
	for (int i = 0; i < clsLAG.size(); i++)
	{
		if (clsLAG[i]->count1 == 3 && clsLAG[i]->count2 == 1) 
		{
			if (mesh.nodes[clsLAG[i]->verts2(0)]->EoL) continue;
			_Aineq.push_back(T(ineqsize, clsLAG[i]->verts2(0) * 3, -clsLAG[i]->nor1(0)));
			_Aineq.push_back(T(ineqsize, clsLAG[i]->verts2(0) * 3 + 1, -clsLAG[i]->nor1(1)));
			_Aineq.push_back(T(ineqsize, clsLAG[i]->verts2(0) * 3 + 2, -clsLAG[i]->nor1(2)));

			ineqsize++;
		}
		else if (clsLAG[i]->count1 == 2 && clsLAG[i]->count2 == 2) 
		{
			if (mesh.nodes[clsLAG[i]->verts2(0)]->EoL || mesh.nodes[clsLAG[i]->verts2(1)]->EoL) continue;
			for (int j = 0; j < 2; j++) 
			{
				_Aineq.push_back(T(ineqsize, clsLAG[i]->verts2(j) * 3, -clsLAG[i]->nor2(0) * clsLAG[i]->weights2(j)));
				_Aineq.push_back(T(ineqsize, clsLAG[i]->verts2(j) * 3 + 1, -clsLAG[i]->nor2(1) * clsLAG[i]->weights2(j)));
				_Aineq.push_back(T(ineqsize, clsLAG[i]->verts2(j) * 3 + 2, -clsLAG[i]->nor2(2) * clsLAG[i]->weights2(j)));
			}

			ineqsize++;
		}
		else if (clsLAG[i]->count1 == 1 && clsLAG[i]->count2 == 3) 
		{
			if (mesh.nodes[clsLAG[i]->verts2(0)]->EoL || mesh.nodes[clsLAG[i]->verts2(1)]->EoL || mesh.nodes[clsLAG[i]->verts2(2)]->EoL) continue;
			for (int j = 0; j < 3; j++) 
			{
				_Aineq.push_back(T(ineqsize, clsLAG[i]->verts2(j) * 3, -clsLAG[i]->nor2(0) * clsLAG[i]->weights2(j)));
				_Aineq.push_back(T(ineqsize, clsLAG[i]->verts2(j) * 3 + 1, -clsLAG[i]->nor2(1) * clsLAG[i]->weights2(j)));
				_Aineq.push_back(T(ineqsize, clsLAG[i]->verts2(j) * 3 + 2, -clsLAG[i]->nor2(2) * clsLAG[i]->weights2(j)));
			}

			ineqsize++;
		}
	}

	if (ineqsize > 0) occur_collision = true;
	
	/*
	if (t < 0.345) //Interactive section lol
	{
		
		_Aeq.push_back(T(eqsize, 2 * 3, 1.0));
		_beq.push_back(make_pair(eqsize, mesh.nodes[2]->v[0]));
		eqsize++;

		_Aeq.push_back(T(eqsize, 2 * 3 + 1, 1.0));
		_beq.push_back(make_pair(eqsize, mesh.nodes[2]->v[1]));
		eqsize++;

		_Aeq.push_back(T(eqsize, 2 * 3 + 2, 1.0));
		_beq.push_back(make_pair(eqsize, mesh.nodes[2]->v[2]));
		eqsize++;

		_Aeq.push_back(T(eqsize, 8 * 3, 1.0));
		_beq.push_back(make_pair(eqsize, mesh.nodes[8]->v[0]));
		eqsize++;

		_Aeq.push_back(T(eqsize, 8 * 3 + 1, 1.0));
		_beq.push_back(make_pair(eqsize, mesh.nodes[8]->v[1]));
		eqsize++;

		_Aeq.push_back(T(eqsize, 8 * 3 + 2, 1.0));
		_beq.push_back(make_pair(eqsize, mesh.nodes[8]->v[2]));
		eqsize++;
		
	}
	else
	{
		_Aeq.push_back(T(eqsize, 2 * 3, 1.0));
		_beq.push_back(make_pair(eqsize, mesh.nodes[2]->v[0]+0.2));
		eqsize++;

		_Aeq.push_back(T(eqsize, 2 * 3 + 1, 1.0));
		_beq.push_back(make_pair(eqsize, mesh.nodes[2]->v[1]));
		eqsize++;

		_Aeq.push_back(T(eqsize, 2 * 3 + 2, 1.0));
		_beq.push_back(make_pair(eqsize, mesh.nodes[2]->v[2]));
		eqsize++;
		
		
		_Aeq.push_back(T(eqsize, 8 * 3, 1.0));
		_beq.push_back(make_pair(eqsize, mesh.nodes[8]->v[0]+0.2));
		eqsize++;

		_Aeq.push_back(T(eqsize, 8 * 3 + 1, 1.0));
		_beq.push_back(make_pair(eqsize, mesh.nodes[8]->v[1]));
		eqsize++;

		_Aeq.push_back(T(eqsize, 8 * 3 + 2, 1.0));
		_beq.push_back(make_pair(eqsize, mesh.nodes[8]->v[2]));
		eqsize++;
	}*/
	
	// up node 2  z:0.1
	// forward both 2 & 8 fix -> both 2 & 8 x:0.2
	// fall down nothing

	Aeq.resize(eqsize, mesh.nodes.size() * 3 + mesh.EoL_Count * 2);
	Aineq.resize(ineqsize, mesh.nodes.size() * 3 + mesh.EoL_Count * 2);
	beq.resize(eqsize);
	bineq.resize(ineqsize);

	Aeq.setFromTriplets(_Aeq.begin(), _Aeq.end());
	Aineq.setFromTriplets(_Aineq.begin(), _Aineq.end());

	beq.setZero();
	bineq.setZero();
	for (int i = 0; i < _beq.size(); i++) 
	{
		beq(_beq[i].first) = _beq[i].second;
	}
	for (int i = 0; i < _bineq.size(); i++) 
	{
		bineq(_bineq[i].first) = _bineq[i].second;
	}
}

