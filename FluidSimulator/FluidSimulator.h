#pragma once
#include "Grid.h"
#include <map>
class FluidSimulator
{
	Grid* simulationGrid;
	
	void advectVelocity();
	void advectPressure();
	void addForces(int i, int j);
	float * calculateNegativeDivergence();
	void project();
	float * approximateVelocity(float * pos);
	float dt;
	float * applyPreconditioner(float * Adiag, float * Aplusi, float * Aplusj, float * precon, float * r, std::map<std::pair<int, int>, int> indices);

public:
	FluidSimulator();
	FluidSimulator(int width, int height, int cellSize);
	Cell advectCell(int i, int j);
	void advect(float timeStep);
	void simulateAndDraw();
	void simulateAndDraw(float timeStep);
	float computeTimeStep();
	void draw();
};

