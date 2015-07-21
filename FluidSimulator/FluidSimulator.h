#pragma once
#include "Grid.h"
#include <map>
class FluidSimulator
{
	Grid* simulationGrid;
	double ** pressureCoefficients;
	void advectVelocity();
	void advectPressure();
	void extrapolateVelocity();
	void addForces(int i, int j);
	double * calculateNegativeDivergence(std::map<std::pair<int, int>, int> indices);
	void project();
	double * approximateVelocity(double * pos);
	double dt;
	double * applyPreconditioner(double * Adiag, double * Aplusi, double * Aplusj, double* Aprevi, double * Aprevj, double * precon, double * r, std::map<std::pair<int, int>, int> indices);
	double * FluidSimulator::apply(double * Adiag, double * Aplusi, double * Aplusj, double * Aprevi, double * Aprevj, double * x, int length, std::map<std::pair<int, int>, int> indices);
public:
	FluidSimulator();
	FluidSimulator(int width, int height, int cellSize);
	Cell advectCell(int i, int j);
	void advect(double timeStep);
	void simulateAndDraw();
	void simulateAndDraw(double timeStep);
	double computeTimeStep();
	void draw();
};

