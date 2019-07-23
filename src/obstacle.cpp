#include "obstacle.h"

using namespace std;
using namespace Eigen;

Obstacle::Obstacle()
{
	box_num = 0;
	point_cloud = new PointCloud;
	collision_threshold = 1e-2;
}

Box::Box(Vector3d _scale, Vector3d _position)
{
	point_num = 8;
	edge_num = 12;

	scale = _scale;
	position = _position;

	face_norm << 1, -1, 0, 0, 0, 0,
				 0, 0, 1, -1, 0, 0,
				 0, 0, 0, 0, 1, -1;
	edge_face << 1, 1, 1, 1, 3, 3, 0, 2, 0, 2, 0, 0,
				 3, 5, 4, 2, 5, 4, 3, 5, 5, 4, 4, 2;
	vert_edge << 0, 0, 1, 2, 4, 5, 7, 9,
				 1, 2, 3, 3, 6, 6, 8, 10,
				 4, 5, 7, 9, 8, 10, 11, 11;
	edge_direction << 4, 2, 2, 4, 0, 0, 4, 0, 2, 0, 2, 4;

	createBoxBuffer();
}

void Box::createBoxBuffer()
{
	float x1 = position.x() + scale.x() / 2.0 * 0.98;
	float x2 = position.x() - scale.x() / 2.0 * 0.98;
	float y1 = position.y() + scale.y() / 2.0 * 0.98;
	float y2 = position.y() - scale.y() / 2.0 * 0.98;
	float z1 = position.z() + scale.z() / 2.0 * 0.98;
	float z2 = position.z() - scale.z() / 2.0 * 0.98;

	float vertices[] = { x1, y1, z1, 0, 0, -1,
						 x1, y2, z1, 0, 0, -1,
						 x2, y2, z1, 0, 0, -1,

						 x1, y1, z1, 0, 0, -1,	
						 x2, y1, z1, 0, 0, -1,
						 x2, y2, z1, 0, 0, -1,				 

						 x1, y1, z2, 0, 0, 1,
						 x1, y2, z2, 0, 0, 1,
						 x2, y2, z2, 0, 0, 1,
						
						 x1, y1, z2, 0, 0, 1,
						 x2, y1, z2, 0, 0, 1,
						 x2, y2, z2, 0, 0, 1,

						 x1, y1, z1, -1, 0, 0,
						 x1, y1, z2, -1, 0, 0,
						 x1, y2, z2, -1, 0, 0,
						 
						 x1, y1, z1, -1, 0, 0,
						 x1, y2, z1, -1, 0, 0,
					     x1, y2, z2, -1, 0, 0,
			
						 x2, y1, z1, 1, 0, 0,
						 x2, y1, z2, 1, 0, 0,
						 x2, y2, z2, 1, 0, 0,

						 x2, y1, z1, 1, 0, 0,
						 x2, y2, z1, 1, 0, 0,
						 x2, y2, z2, 1, 0, 0,

						 x1, y1, z1, 0, -1, 0,
						 x1, y1, z2, 0, -1, 0,
						 x2, y1, z2, 0, -1, 0,

						 x1, y1, z1, 0, -1, 0,
						 x2, y1, z1, 0, -1, 0,
						 x2, y1, z2, 0, -1, 0,

						 x1, y2, z1, 0, 1, 0,
						 x1, y2, z2, 0, 1, 0,
						 x2, y2, z2, 0, 1, 0,

						 x1, y2, z1, 0, 1, 0,
						 x2, y2, z1, 0, 1, 0,
						 x2, y2, z2, 0, 1, 0};
	buffer = (float*)malloc(sizeof(vertices));
	memcpy(buffer, vertices, sizeof(vertices));
	buffer_size = sizeof(float) * 6 * 6 * 6;
}