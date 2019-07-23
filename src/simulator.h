#pragma once
#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <vector>
#include <memory>
#include <string>

#include <Eigen/Dense>

#include "Collisions.h"



class Cloth;
class Obstacle;
class Constraint;

class Simulator
{
public:
	Simulator();
	virtual ~Simulator(){}

	void step();

public:
	double h;

	Eigen::Vector3d gravity;

	Cloth* cloth;
	Obstacle* obstacle;
	std::vector<Collision*> collisions;

private:
	double t;
};


#endif