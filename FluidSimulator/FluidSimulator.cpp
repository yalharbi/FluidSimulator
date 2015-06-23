#include "FluidSimulator.h"


FluidSimulator::FluidSimulator()
{
}

FluidSimulator::FluidSimulator(int width, int height, int cellSize)
{
	simulationGrid = new Grid(width, height, cellSize);
}

void FluidSimulator::advect(float timeStep){
	dt = timeStep;
	advectVelocity();
	advectPressure();
}

void FluidSimulator::advectVelocity(){
	for (int i = 0; i < simulationGrid->width; i++){
		for (int j = 0; j < simulationGrid->height; j++){
			advectCell(i, j);
		}
	}
}

void FluidSimulator::advectCell(int i, int j){
	float * position = simulationGrid->getCellPosition(i, j);
	float * apprU = approximateVelocity(position);
}

float * FluidSimulator::approximateVelocity(float *){
	

}

void FluidSimulator::advectPressure(){

}

float FluidSimulator::computeTimeStep(){
	float cellSize = simulationGrid->getCellSize();
	float maxVel = simulationGrid->getMaxVelocity();

	dt = cellSize / maxVel;
	return dt;
}