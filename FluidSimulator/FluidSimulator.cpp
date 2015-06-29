#include "FluidSimulator.h"

/* LOG
6/28/2015: According to (Foster and Fedkiw 2001) air cells only have velocity at their boundaries with fluid cells, also (Bridson 2008) p24
		   this means u(i,j) is actually the average of the two boundaries, otherwise fluids won't flow into air cells whos velocity is 0 (because we can't trace them back)
		   TODO: (1) Set all velocities between air cells to 0 (2) Fix getVelocity(i,j) to account for inflow and outflow (two normal velocities)
		   
*/

FluidSimulator::FluidSimulator()
{
}

FluidSimulator::FluidSimulator(int width, int height, int cellSize)
{
	simulationGrid = new Grid(width, height, cellSize);
	simulationGrid->setCell(3, 4, FLUID); simulationGrid->setCell(4, 4, FLUID); simulationGrid->setCell(5, 4, FLUID);
	simulationGrid->setCell(3, 3, FLUID); simulationGrid->setCell(4, 3, FLUID); simulationGrid->setCell(5, 3, FLUID);
	for (int i = 0; i < simulationGrid->width; i++)
		simulationGrid->setCell(i, 0, SOLID);
	for (int i = 0; i < simulationGrid->height; i++){
		simulationGrid->setCell(0, i, SOLID);
		simulationGrid->setCell(simulationGrid->width-1, i, SOLID);
	}
	simulationGrid->showVectorField();
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


	computeTimeStep();
	for (int i = 0; i < simulationGrid->width; i++){
		for (int j = 1; j < simulationGrid->height; j++){
			addForces(i, j);
		}
	}

	for (int i = 0; i < simulationGrid->width; i++){
		for (int j = 0; j < simulationGrid->height; j++){
			updatedCells[i][j] = advectCell(i, j);
		}
	}

	simulationGrid->cells = updatedCells;
	//std::cout << simulationGrid->getVelocityVector(4, 2);
}

void FluidSimulator::addForces(int i, int j){
	if (simulationGrid->getCellType(i, j) == FLUID){
		simulationGrid->cells[i][j].v += dt*(-9.8);
	}

	if (simulationGrid->getCellType(i, j) == FREE && simulationGrid->getCellType(i, j+1) == FLUID){
		simulationGrid->cells[i][j].v += simulationGrid->getVVelocityAt(i, j+1);

	}
}

Cell FluidSimulator::advectCell(int i, int j){

	if (simulationGrid->cells[i][j].cellType == SOLID)
		return simulationGrid->cells[i][j];

	Vector position = simulationGrid->getCellPosition(i, j);
	Vector midPos = position - simulationGrid->getVelocityVector(i, j)*dt/2;
	Vector tracedParticle = position - simulationGrid->interpolateVelocity(midPos)*dt;
	if (i == 4 && j == 4){
		std::cout << "midpos: " << midPos;
		std::cout << "dt: " << dt << " pos: " << position << " traced: " << tracedParticle << " vel: " << simulationGrid->interpolateVelocity(midPos) << "size: " << simulationGrid->getCellSize() << "\n";
		std::cout << "tr: " << simulationGrid->getCell(tracedParticle).v << "\n";
	}
	simulationGrid->interpolateVelocity(tracedParticle);
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