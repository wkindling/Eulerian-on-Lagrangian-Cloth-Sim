#include "preprocessor.h"

#include <stdlib.h>

#include "ArcSim/remesh.hpp"
#include "ArcSim/geometry.hpp"
#include "ArcSim/subset.hpp"
#include "ArcSim/dynamicremesh.hpp"

#include "utility.h"

using namespace std;
using namespace Eigen;

bool Preprocessor::onBoundaryA(const MatrixXd &boundary, const Vector3d & P) 
{
	int corner1, corner2;

	for (corner1 = 0; corner1 < boundary.cols(); corner1++)
	{
		corner2 = (corner1 + 1) % boundary.cols();
		Vector3d cor1 = boundary.block<3, 1>(0, corner1);
		Vector3d cor2 = boundary.block<3, 1>(0, corner2);

		if (vertLineDist(e2v(cor1), e2v(cor2), e2v(P)) < 0.025)
		{
			return true;
		}
	}
	return false;
}

bool Preprocessor::onBoundaryB(const MatrixXd& boundary, const Vector3d& P, const Vec3& e2)
{
	int corner1, corner2;
	int corner_count = 0;

	for (corner1 = 0; corner1 < boundary.cols(); corner1++)
	{
		corner2 = (corner1 + 1) % boundary.cols();
		Vector3d cor1 = boundary.block<3, 1>(0, corner1);
		Vector3d cor2 = boundary.block<3, 1>(0, corner2);

		if (vertLineDist(e2v(cor1), e2v(cor2), e2v(P)) < 0.025)
		{
			corner_count++;
			if (corner_count > 1) return true;
			Vector3d boundary_edge = cor1 - cor2;
			double angle = get_angle(e2v(boundary_edge), e2);
			if (angle<M_PI / 5 || angle>(4 * M_PI / 5)) return true;
		}
	}
	return false;
}

bool Preprocessor::onBoundary(const MatrixXd& boundary, const Node* node)
{
	if (!node->EoL) return false;
	int corner1, corner2;
	int corner_count = 0, preserve = 0;
	Vert* vert = node->verts[0];

	for (corner1 = 0; corner1 < boundary.cols(); corner1++)
	{
		corner2 = (corner1 + 1) % boundary.cols();
		Vector3d cor1 = boundary.block<3, 1>(0, corner1);
		Vector3d cor2 = boundary.block<3, 1>(0, corner2);

		if (vertLineDist(e2v(cor1), e2v(cor2), vert->u) < 0.025)
		{
			if (node->cornerID >= 0) return true;
			corner_count++;
			if (corner_count > 1) return true;
			for (int e = 0; e < node->adje.size(); e++)
			{
				Edge* edge = node->adje[e];
				if (edge->preserve)
				{
					preserve++;
					Vector3d boundary_edge = cor1 - cor2;
					Node* node1 = other_node(edge, node);
					double angle = get_angle(e2v(boundary_edge), (vert->u - node1->verts[0]->u));
					if (angle<M_PI / 5 || angle>(4 * M_PI) / 5) return true;
				}
			}
			if (preserve == 0) return true;
		}
	}
	return false;
}

/*
addEOL: detect the collision point, the point can be exactly on the existing vert or on the edge or in the triangle, split the topology and insert new point, mark with flag newEOL
3 types : face - vert , vert - face, edge - edge

Type1. face - vert: just mark the EoL nodes to be  isEOL, in case for EoL transfer

Type2. vert - face: the obstacle's point may puncture the face, just get the barycentric coordinate, split and insert new node mark with newEOL

Type3. edge - edge: the cloth's edge has collision with obstacle's, get the one-dimension coordinate and split the edge
*/


void Preprocessor::addEOL(Mesh& mesh, const MatrixXd& boundary, vector<Collision*> collisions)
{
	for (int n = 0; n < mesh.nodes.size(); n++)
	{
		if (mesh.nodes[n]->EoL)
		{
			mesh.nodes[n]->EoL_state = Node::WasEOL;
		}
	}

	for (int i = 0; i < collisions.size(); i++)
	{		
		/*Obstacle's face -------------  cloth's vert*/
		if (collisions[i]->count1 == 3 && collisions[i]->count2 == 1) 
		{
			Node* node = mesh.nodes[collisions[i]->verts2(0)];
			if (node->EoL)
			{
				node->EoL_state = Node::IsEOL; // mark with isEoL for EoL transfer 
			}
		}
		/*Obstacle's vert -------------- cloth's face*/
		else if (collisions[i]->count1 == 1 && collisions[i]->count2 == 3) 
		{
			if ((mesh.nodes[collisions[i]->verts2(0)]->EoL && mesh.nodes[collisions[i]->verts2(0)]->cornerID == collisions[i]->verts1(0)) ||
				(mesh.nodes[collisions[i]->verts2(1)]->EoL && mesh.nodes[collisions[i]->verts2(1)]->cornerID == collisions[i]->verts1(0)) ||
				(mesh.nodes[collisions[i]->verts2(2)]->EoL && mesh.nodes[collisions[i]->verts2(2)]->cornerID == collisions[i]->verts1(0))) continue; // If the triangle's point touches obstacle's vert

			Vert *v0 = mesh.verts[collisions[i]->verts2(0)],
				 *v1 = mesh.verts[collisions[i]->verts2(1)],
				 *v2 = mesh.verts[collisions[i]->verts2(2)];

			double xX = collisions[i]->weights2(0) * v0->u[0] + collisions[i]->weights2(1) * v1->u[0] + collisions[i]->weights2(2) * v2->u[0];
			double yX = collisions[i]->weights2(0) * v0->u[1] + collisions[i]->weights2(1) * v1->u[1] + collisions[i]->weights2(2) * v2->u[1]; // The coordinate of collision point

			// Boundary
			if (onBoundaryA(boundary, Vector3d(xX, yX, 0.0))) continue; // If the point is on boundary, just skip

			
			Face *f0 = get_enclosing_face(mesh, Vec2(xX, yX));
			Vec3 barycentric = get_barycentric_coords(Vec2(xX, yX), f0);

			bool on_edge = false;
			for (int j = 0; j < 3; j++) 
			{
				if (barycentric[j] < 1e-3) on_edge = true; // The collision point is exactly on the boundary 
			}

			// If this point is on a cloth edge, we should split the edge
			if (on_edge) 
			{
				double least = 1.0;
				int opp_index = -1;
				for (int j = 0; j < 3; j++) 
				{
					if (barycentric[j] < least) 
					{
						least = barycentric[j];
						opp_index = j;
					}
				}
				Edge *e0 = get_opp_edge(f0, f0->v[opp_index]->node);
				double split_weight;
				if (e0->n[0] == f0->v[0]->node) // If v0 is on the edge
				{
					split_weight = 1.0 - barycentric[0];
				}
				else if (e0->n[0] == f0->v[1]->node) // If v1 is on the edge 
				{
					split_weight = 1.0 - barycentric[1];
				}
				else  // If v2 is on the edge
				{
					split_weight = 1.0 - barycentric[2];
				}
				Node *node0 = e0->n[0], *node1 = e0->n[1];
				RemeshOp op = split_edgeForced(e0, split_weight, -1);
				for (size_t v = 0; v < op.added_verts.size(); v++) 
				{
					Vert *new_vert = op.added_verts[v];
					Vert *v0 = adjacent_vert(node0, new_vert),
						 *v1 = adjacent_vert(node1, new_vert);
					new_vert->sizing = 0.5 * (v0->sizing + v1->sizing); // Update the sizing
				}
				op.done();
			}
			else  //The collision point is just in the face, lol, it's easy to deal with
			{
				RemeshOp op = split_face(f0, barycentric); 
				for (size_t v = 0; v < op.added_verts.size(); v++) 
				{
					Vert *new_vert = op.added_verts[v];
					new_vert->sizing = (v0->sizing + v1->sizing + v2->sizing) / 3.0;
				}
				op.done();
			}

			Node *n = mesh.nodes.back(); // Got the newly inserted point
			n->EoL = true;
			n->EoL_state = Node::NewEOL;
			n->preserve = true;
			n->cornerID = collisions[i]->verts1(0);
			n->cdEdges = collisions[i]->edge1;
			n->x = e2v(collisions[i]->pos1_); 
		}
		/*Obstacle's edge ------------ Cloth's edge*/
		if (collisions[i]->count1 == 2 && collisions[i]->count2 == 2)
		{
			if (mesh.nodes[collisions[i]->verts2(0)]->EoL || mesh.nodes[collisions[i]->verts2(1)]->EoL) continue; // If the collision point is exactly on the end

			// The verts won't change since we only ever add verts
			Edge *e0 = get_edge(mesh.nodes[collisions[i]->verts2(0)], mesh.nodes[collisions[i]->verts2(1)]);
			
			double split_weight;

			// We'll need this info for the boundary
			Vert *v0 = mesh.verts[collisions[i]->verts2(0)],
				 *v1 = mesh.verts[collisions[i]->verts2(1)];

			double xX = collisions[i]->weights2(0) * v0->u[0] + collisions[i]->weights2(1) * v1->u[0];
			double yX = collisions[i]->weights2(0) * v0->u[1] + collisions[i]->weights2(1) * v1->u[1]; // Get the position of the collision point

			Face *f0 = get_enclosing_face(mesh, Vec2(xX, yX));

			// Boundary simple
			MatrixXd F = computeF(f0);
			Vector2d e2 = F.transpose() * collisions[i]->edgeDir; // get the edge direction in the material space 
			if (onBoundaryB(boundary, Vector3d(xX, yX, 0.0), Vec3(e2(0), e2(1), 0.0))) continue;

			// Some previous collisions has already split this edge 
			if (e0 == NULL) 
			{
				Vec3 barycentric = get_barycentric_coords(Vec2(xX, yX), f0);
				double least = 1.0;
				int opp_index = -1;
				for (int j = 0; j < 3; j++) 
				{
					if (barycentric[j] < least) 
					{
						least = barycentric[j];
						opp_index = j;
					}
				}
				e0 = get_opp_edge(f0, f0->v[opp_index]->node); //find the new triangle, then the same as before
				if (e0->n[0] == f0->v[0]->node) 
				{
					split_weight = 1.0 - barycentric[0];
				}
				else if (e0->n[0] == f0->v[1]->node) 
				{
					split_weight = 1.0 - barycentric[1];
				}
				else 
				{
					split_weight = 1.0 - barycentric[2];
				}
			}
			else 
			{
				//Just to fit in arcsim, use the weight of n1
				if (e0->n[0]->index == collisions[i]->verts2(1)) 
				{
					split_weight = collisions[i]->weights2(0);
				}
				else {
					split_weight = collisions[i]->weights2(1);
				}
			}

			Node *node0 = e0->n[0], *node1 = e0->n[1];
			RemeshOp op = split_edgeForced(e0, split_weight, -1);
			for (size_t v = 0; v < op.added_verts.size(); v++) 
			{
				Vert *new_vert = op.added_verts[v];
				Vert *v0 = adjacent_vert(node0, new_vert),
					*v1 = adjacent_vert(node1, new_vert);
				new_vert->sizing = 0.5 * (v0->sizing + v1->sizing);
			}
			op.done();

			Node *n = mesh.nodes.back();
			n->EoL = true;
			n->EoL_state = Node::NewEOL;
			n->cdEdges = collisions[i]->edge1;
			n->x = e2v(collisions[i]->pos1_); 
		}
	}
}

void Preprocessor::markPreserve(Mesh& mesh)
{
	for (int i = 0; i < mesh.edges.size(); i++)
	{
		Edge *e = mesh.edges[i];
		e->preserve = false;
		Node *n0 = e->n[0], *n1 = e->n[1];
		if (n0->EoL && n1->EoL)// the edge between 2 EoL points should not be splitted
		{
			if (n0->cornerID >= 0 && n1->cornerID >= 0) continue; //Just for point collision case

			bool match = false;
			if (n0->cornerID >= 0) //corner 2 edge
			{
				for (int j = 0; j < n0->cdEdges.size(); j++)
				{
					if (n1->cdEdges[0] == n0->cdEdges[j])
					{
						match = true;
						break;
					}
				}
			}
			else if (n1->cornerID >= 0)  // corner 2 edge
			{
				for (int j = 0; j < n1->cdEdges.size(); j++)
				{
					if (n0->cdEdges[0] == n1->cdEdges[j]) 
					{
						match = true;
						break;
					}
				}
			}
			else if (n0->cdEdges[0] == n1->cdEdges[0])  //edge 2 edge
			{
				match = true;
			}
			if (match) e->preserve = true;
		}
	}
}

void Preprocessor::cleanMesh(Mesh& mesh, const Eigen::MatrixXd& boundary)
{
	vector<Face*> active_faces = mesh.faces;
	flip_edges(0, active_faces, 0, 0);
	markPreserve(mesh);
	bool allclear = false;

	while (!allclear)
	{
		int iter = 0;
		while (!allclear)
		{
			iter++;
			if (iter == 3) break; // in case be trapped in infinite loop
			allclear = true;
			while (collapse_nonconformal(mesh, allclear));
			markPreserve(mesh); 
			while (collapse_conformal(mesh, allclear));
			markPreserve(mesh);
		}
		allclear = split_illconditioned_faces(mesh);
	}
	markPreserve(mesh);

	for (int n = 0; n < mesh.nodes.size(); n++)
	{
		Node* node = mesh.nodes[n];

		if (onBoundary(boundary, node) || node->EoL_state == Node::WasEOL)
		{
			node->EoL = false;
			node->preserve = false;
			node->cornerID = -1;
			node->cdEdges.clear();
		}
	}

	compute_ws_data(mesh); //external call
}

void Preprocessor::preprocess(Mesh& mesh, const MatrixXd& boundary, vector<Collision*> collisions) // Preserve done!
{
	addEOL(mesh, boundary, collisions);

	cleanMesh(mesh, boundary);
}