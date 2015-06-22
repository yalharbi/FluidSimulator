#pragma once
#include "Grid.h"

class FluidSimulator
{
	Grid* simulationGrid;
	void advect(float timeStep);


public:
	FluidSimulator();
	FluidSimulator(int width, int height, int cellSize);

	void simulateAndDraw(float timeStep);
};

