#include "FluidSimulator.h"


FluidSimulator::FluidSimulator()
{
}

FluidSimulator::FluidSimulator(int width, int height, int cellSize)
{
	simulationGrid = new Grid(width, height, cellSize);
	simulationGrid->setCell(3, 4, FLUID); simulationGrid->setCell(4, 4, FLUID); simulationGrid->setCell(5, 4, FLUID);
}

void FluidSimulator::advect(float timeStep){
	dt = timeStep;
	advectVelocity();
	advectPressure();
}

void FluidSimulator::advectVelocity(){
	Cell ** updatedCells = new Cell*[simulationGrid->width];
	for (int i = 0; i < simulationGrid->width; i++){
		updatedCells[i] = new Cell[simulationGrid->height];
	}



	for (int i = 0; i < simulationGrid->width; i++){
		for (int j = 0; j < simulationGrid->height; j++){
			updatedCells[i][j] = advectCell(i, j);
		}
	}

	simulationGrid->cells = updatedCells;
}

Cell FluidSimulator::advectCell(int i, int j){

	if (i == 0 || j == 0 || i == simulationGrid->width-1 || j == simulationGrid->height-1)
		return simulationGrid->cells[i][j];

	computeTimeStep();
	Vector position = simulationGrid->getCellPosition(i, j);
	Vector midPos = position - simulationGrid->getVelocityVector(i, j)*dt/2;
	Vector tracedParticle = position - simulationGrid->interpolateVelocity(midPos)*dt;

	//std::cout << "dt: " << dt << " pos: " << position << " traced: " << tracedParticle << " vel: " << simulationGrid->interpolateVelocity(midPos) << "\n";
	//simulationGrid->interpolateVelocity(tracedParticle);
	return simulationGrid->getCell(tracedParticle);
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
	simulationGrid->draw();
}