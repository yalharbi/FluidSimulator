#pragma once
#include "Grid.h"

class FluidSimulator
{
	Grid* simulationGrid;
	void advect(float timeStep);
	void advectVelocity();
	void advectPressure();
	void advectCell(int i, int j);
	float * approximateVelocity(float * pos);
	float dt;

public:
	FluidSimulator();
	FluidSimulator(int width, int height, int cellSize);

	void simulateAndDraw();
	void simulateAndDraw(float timeStep);
	float computeTimeStep();
};

