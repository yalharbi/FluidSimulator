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
	Vector **updatedV = new Vector*[simulationGrid->width];
	for (int i = 0; i < simulationGrid->width; i++)
		updatedV[i] = new Vector[simulationGrid->height];


	for (int i = 0; i < simulationGrid->width; i++){
		for (int j = 0; j < simulationGrid->height; j++){
			updatedV[i][j] = advectCell(i, j);
		}
	}
}

Vector FluidSimulator::advectCell(int i, int j){
	if (i == 0 || j == 0 || i == simulationGrid->width || j == simulationGrid->height)
		return Vector(0,0,0);
	computeTimeStep();
	Vector position = simulationGrid->getCellPosition(i, j);
	Vector midPos = position - simulationGrid->getVelocityVector(i, j)*dt/2;
	Vector tracedParticle = position - simulationGrid->interpolateVelocity(midPos)*dt;

	std::cout << "dt: " << dt << " pos: " << position << " traced: " << tracedParticle << " vel: " << simulationGrid->interpolateVelocity(midPos) << "\n";
	return simulationGrid->interpolateVelocity(tracedParticle);
}


void FluidSimulator::advectPressure(){

}

float FluidSimulator::computeTimeStep(){
	float cellSize = simulationGrid->getCellSize();
	float maxVel = simulationGrid->getMaxVelocity();

	dt = cellSize / maxVel;
	return dt;
}

void FluidSimulator::draw(){
	simulationGrid->setCell(2, 1, FLUID);
	advectCell(2, 1);
	simulationGrid->draw();
}