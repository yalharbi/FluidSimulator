#pragma once
#include "Grid.h"

class FluidSimulator
{
	Grid* simulationGrid;
	void advect(float timeStep);
	void advectVelocity();
	void advectPressure();
	
	float * approximateVelocity(float * pos);
	float dt;

public:
	FluidSimulator();
	FluidSimulator(int width, int height, int cellSize);
	Vector advectCell(int i, int j);
	void simulateAndDraw();
	void simulateAndDraw(float timeStep);
	float computeTimeStep();
	void draw();
};

