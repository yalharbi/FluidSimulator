#include "FluidSimulator.h"
#define DENSITY 1000
#define MAX_ITERATIONS 1000
#define TOLERANCE 1/1000000

#define square(x) x*x

inline float dotproduct(float * a, float * b, int length){
	float result = 0;
	for (int it = 0; it < length; it++)
		result += a[it] * b[it];
	if (result <0.001)
		return 0.001;
	return result;
}

float * FluidSimulator::apply(float * Adiag, float * Aplusi, float * Aplusj, float * Aprevi, float * Aprevj, float * x, int length, std::map<std::pair<int, int>, int> indices){
	float * ans = new float[length]();

	for (int i = 0; i < simulationGrid->width; i++){
		for (int j = 0; j < simulationGrid->height; j++){
			if (simulationGrid->getCellType(i, j) == FLUID){
				int place = indices[std::make_pair(i, j)];
				int placepj = indices[std::make_pair(i, j - 1)];
				int placepi = indices[std::make_pair(i - 1, j)];
				int placenj = indices[std::make_pair(i, j + 1)];
				int placeni = indices[std::make_pair(i + 1, j)];
				ans[place] = Adiag[place] * x[place];
				ans[place] += Aplusi[place] * x[placeni];
				ans[place] += Aplusj[place] * x[placenj];
				ans[place] += Aprevi[place] * x[placepi];
				ans[place] += Aprevj[place] * x[placepj];
					/*
				for (int k = place - 1; k >= 0; k--){
					int symIndex = indices[std::make_pair(k,j)];
					ans[place] += Aplusi[place-1] * x[place-1];
					ans[place] += Aplusj[place-1] * x[place-1];
					std::cout << "place" << place<< "placeSym" << symIndex<<"in" << "i" << i << "j" << j << "k" << k << "\n";
				}*/
				/*if (simulationGrid->getCellType(i - 1, j) == FLUID)
					ans[place] += Aplusi[placepi] * x[placepi];
				if (simulationGrid->getCellType(i, j - 1) == FLUID)
					ans[place] += Aplusj[placepj] * x[placepj];*/
			}
		}
	}
	

	return ans;
}

inline float * copyV(float * from, int length){
	float * to = new float[length]();
	for (int i = 0; i < length; i++)
		to[i] = from[i];
	return to;
}
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
	/*for (int i = 5; i < 10; i++)
		for (int j = 5; j < 10; j++)
			simulationGrid->setCell(i, j, FLUID);*/
	
	simulationGrid->setCell(3, 5, FLUID); simulationGrid->setCell(4, 5, FLUID); simulationGrid->setCell(5, 5, FLUID);
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
	if (simulationGrid->getCellType(i, j) != SOLID){
		simulationGrid->cells[i][j].v += (-9.8)*dt;
	}

	/*if (simulationGrid->getCellType(i, j) == FREE && simulationGrid->getCellType(i, j-1) == FLUID){
		simulationGrid->cells[i][j].v += (-9.8)*dt;

	}*/
}

Cell FluidSimulator::advectCell(int i, int j){

	if (simulationGrid->cells[i][j].cellType == SOLID)
		return simulationGrid->cells[i][j];


	Vector position = simulationGrid->getCellPosition(i, j);
	Vector midPos = position - simulationGrid->getVelocityVector(i, j)*dt/2;
	Vector tracedParticle = position - simulationGrid->interpolateVelocity(midPos)*dt;
	if ((i == 2 && j == 1) || (i == 4 && j == 2)){
		std::cout << "midpos: " << midPos;
		std::cout << "dt: " << dt << " pos: " << position << " traced: " << tracedParticle << " vel: " << simulationGrid->interpolateVelocity(midPos) << "size: " << simulationGrid->getCellSize() << "\n";
		std::cout << "tr: " << simulationGrid->getCell(tracedParticle).v << "\n";
	}
	simulationGrid->interpolateVelocity(tracedParticle);
	return simulationGrid->getCell(tracedParticle);
}


void FluidSimulator::advectPressure(){

}

float * FluidSimulator::calculateNegativeDivergence(){
	float* rhs = new float[simulationGrid->fluidCellCount];// (float *)malloc(sizeof(float) * simulationGrid->fluidCellCount);
	int place = 0;
	float scale = 1 / sqrt(simulationGrid->getCellSize());
	std::map<std::pair<int, int>, int> indices;

	for (int i = 0; i < simulationGrid->width; i++){
		for (int j = 0; j < simulationGrid->height; j++){
			if (simulationGrid->getCellType(i, j) != FLUID)
				continue;
			Cell cell = simulationGrid->cells[i][j];
			rhs[place++] = (simulationGrid->getHVelocityAt(i + 1, j) - simulationGrid->getHVelocityAt(i, j) +
				simulationGrid->getVVelocityAt(i, j + 1) - simulationGrid->getVVelocityAt(i, j) + 0) *  -scale;
			indices.insert(std::make_pair(std::make_pair(i, j), place - 1));

			std::cout << rhs[place - 1] << "\n";
		}
	}

	for (int i = 0; i < simulationGrid->width; i++){
		for (int j = 0; j < simulationGrid->height; j++){
			if (simulationGrid->getCellType(i, j) != FLUID)
				continue;
			Cell cell = simulationGrid->cells[i][j];
			int place = indices[std::make_pair(i, j)];
			if (simulationGrid->getCellType(i, j + 1) == SOLID){
				rhs[place] = rhs[place] + simulationGrid->getVVelocityAt(i, j+1)*scale;
			}
			if (simulationGrid->getCellType(i, j - 1) == SOLID){
				rhs[place] = rhs[place] - simulationGrid->getVVelocityAt(i, j)*scale;
			}
			if (simulationGrid->getCellType(i+1, j) == SOLID){
				rhs[place] = rhs[place] + simulationGrid->getHVelocityAt(i+1, j)*scale;
			}
			if (simulationGrid->getCellType(i - 1, j) == SOLID){
				rhs[place] = rhs[place] - simulationGrid->getHVelocityAt(i, j)*scale;
			}

			std::cout << "^" << rhs[place] << "\n";
		}
	}
	return rhs;
}

float * FluidSimulator::applyPreconditioner(float * Adiag, float * Aplusi, float * Aplusj, float * precon, float * r, std::map<std::pair<int, int>, int> indices){
	float * q = new float[simulationGrid->fluidCellCount]();

	for (int i = 1; i < simulationGrid->width; i++){ //Lq=r
		for (int j = 1; j < simulationGrid->height; j++){
			int place = indices[std::make_pair(i, j)];
			int placepj = indices[std::make_pair(i, j - 1)];
			int placepi = indices[std::make_pair(i - 1, j)];

			if (simulationGrid->getCellType(i, j) == FLUID){
				float t = r[place] - Aplusi[placepi] * precon[placepi] * q[placepi]
					- Aplusj[placepj] * precon[placepj] * q[placepj];
				q[place] = t*precon[place];
			}
		}
	}

	float * z = new float[simulationGrid->fluidCellCount]();
	for (int i = simulationGrid->width-1; i > 0; i--){ //Ltz = q
		for (int j = simulationGrid->height-1; j > 0; j--){
			int place = indices[std::make_pair(i, j)];
			int placenj = indices[std::make_pair(i, j + 1)];
			int placeni = indices[std::make_pair(i + 1, j)];

			if (simulationGrid->getCellType(i, j) == FLUID){
				float t = q[place] - Aplusi[place] * precon[place] * z[placeni]
					- Aplusj[place] * precon[place] * z[placenj];
				z[place] = t*precon[place];
			}
		}
	}
	return z;
}

void FluidSimulator::project(){
	std::map<std::pair<int, int>, int> indices;
	int index = 0;
	for (int i = 0; i < simulationGrid->width; i++){
		for (int j = 0; j < simulationGrid->height; j++){
			if (simulationGrid->getCellType(i, j) == FLUID)
				indices.insert(std::make_pair(std::make_pair(i, j), index++));
		}
	}


	float * rhs = calculateNegativeDivergence();

	for (int i = 0; i < simulationGrid->fluidCellCount; i++){ // return if pressure is already divergent free
		if (rhs[i] != 0)
			break;
		if (i + 1 == simulationGrid->fluidCellCount){
			std::cout << "divFree\n";
			return;
		}
	}
	float * Adiag = new float[simulationGrid->fluidCellCount]();
	float * Aplusi = new float[simulationGrid->fluidCellCount]();
	float * precon = new float[simulationGrid->fluidCellCount]();
	float * Aplusj = new float[simulationGrid->fluidCellCount]();
	float * Aprevi = new float[simulationGrid->fluidCellCount]();
	float * Aprevj = new float[simulationGrid->fluidCellCount]();

	float scale = dt / (DENSITY*simulationGrid->getCellSize());
	for (int i = 0; i < simulationGrid->width; i++){ //compute A... where A*pressure(unknown) = negative divergence... pressure makes the field divergent free
		for (int j = 0; j < simulationGrid->height; j++){
			int place = indices[std::make_pair(i, j)];
			if (simulationGrid->getCellType(i, j) == FLUID && simulationGrid->getCellType(i + 1, j) == FLUID){
				Adiag[place] += scale;
				Adiag[indices[std::make_pair(i + 1, j)]] += scale;
				Aplusi[place] = -scale;
				Aprevi[indices[std::make_pair(i + 1, j)]] = -scale;
			}
			else if (simulationGrid->getCellType(i, j) == FLUID && simulationGrid->getCellType(i + 1, j) == FREE)
				Adiag[place] += scale;

			if (simulationGrid->getCellType(i, j) == FLUID && simulationGrid->getCellType(i, j+1) == FLUID){
				Adiag[place] += scale;
				Adiag[indices[std::make_pair(i, j+1)]] += scale;
				Aplusj[place] = -scale;
				Aprevj[indices[std::make_pair(i, j+1)]] = -scale;
			}
			else if (simulationGrid->getCellType(i, j) == FLUID && simulationGrid->getCellType(i, j+1) == FREE)
				Adiag[place] += scale;

			
		}
	}
	float tau = 0.97, beta = 0.25;
	for (int i = 1; i < simulationGrid->width; i++){ //compute preconditioner... initial guess for solving the system
		for (int j = 1; j < simulationGrid->height; j++){
			int place = indices[std::make_pair(i, j)];
			int placepj = indices[std::make_pair(i, j-1)];
			int placepi = indices[std::make_pair(i-1, j)];

			if (simulationGrid->getCellType(i, j) == FLUID){
				float e = Adiag[place] - square(Aplusi[placepi] * precon[placepi]) - square(Aplusj[placepj] * precon[placepj])
					- tau*(Aplusi[placepi] * Aplusj[placepi] *square(precon[placepi])
					+ Aplusj[placepj] * Aplusi[placepj]  *square(precon[placepj]));
				if (e < beta * Adiag[place])
					e = Adiag[place];
				if (e == 0)
					e = 0.1;
				precon[place] = 1 / sqrt(e);
			}
		}
	}

	float * p = new float[simulationGrid->fluidCellCount]();
	float * r = copyV(rhs, simulationGrid->fluidCellCount);
	float * z = applyPreconditioner(Adiag, Aplusi, Aplusj, precon, r, indices);
	float * s = copyV(z, simulationGrid->fluidCellCount);

	float dp = dotproduct(z, r, simulationGrid->fluidCellCount);

	bool done = false;

	int iterations = 0;
	while (!done && iterations < MAX_ITERATIONS){
		z = apply(Adiag, Aplusi, Aplusj, Aprevi, Aprevj, s, simulationGrid->fluidCellCount, indices);
		float alph = dp / dotproduct(z, s, simulationGrid->fluidCellCount);

		float maxR = 0;
	//	std::cout << "s0" << s[0];
		for (int i = 0; i < simulationGrid->fluidCellCount; i++){
			p[i] = p[i] + alph*s[i];
			r[i] = r[i] - alph*z[i];
			if (r[i] > maxR)
				maxR = r[i];
		}
		if (abs(maxR) <= TOLERANCE){
			done = true;
			continue;
		}
		z = applyPreconditioner(Adiag, Aplusi, Aplusj, precon, r, indices);
		float dpNew = dotproduct(z, r, simulationGrid->fluidCellCount);
		float beeta = dpNew / dp;
	//	std::cout << "beeta" << beeta << "dp" << dp << "dpnew" << dpNew;
		for (int i = 0; i < simulationGrid->fluidCellCount; i++)
			s[i] = z[i] + beeta*s[i];
		dp = dpNew;
		iterations++;
		if (iterations + 1 == MAX_ITERATIONS || done) std::cout << "MX" << maxR << "#";
	}
	std::cout << "uterations: " << iterations << "\n";

	for (int i = 0; i < simulationGrid->fluidCellCount; i++)
		std::cout << "p: " << p[i] << "\n";

	for (int i = 0; i < simulationGrid->width; i++){ 
		for (int j = 0; j < simulationGrid->height; j++){
			float sc =  dt / (DENSITY * sqrt(simulationGrid->getCellSize()));
			int place = indices[std::make_pair(i, j)];
			int placenj = indices[std::make_pair(i, j + 1)];
			int placeni = indices[std::make_pair(i + 1, j)];

			if (simulationGrid->getCellType(i, j) == FLUID){
				simulationGrid->cells[i][j].u -= sc*p[place];
				simulationGrid->cells[i + 1][j].u += sc*p[place];

				simulationGrid->cells[i][j].v -= sc*p[place];
				simulationGrid->cells[i][j+1].v += sc*p[place];
			}
		}
	}
	for (int i = 0; i < simulationGrid->width; i++){
		for (int j = 0; j < simulationGrid->height; j++){

			if (simulationGrid->getCellType(i, j) == SOLID){
				simulationGrid->cells[i][j].u = 0;
				if (i + 1 <simulationGrid->width)
					simulationGrid->cells[i + 1][j].u = 0;

				simulationGrid->cells[i][j].v = 0;
				if (j + 1 <simulationGrid->height)
					simulationGrid->cells[i][j + 1].v = 0;
			}

		}
	}
	
}

float FluidSimulator::computeTimeStep(){
	float cellSize = simulationGrid->getCellSize();
	float maxVel = simulationGrid->getMaxVelocity();

	dt = cellSize / maxVel;
	dt = 0.5;
	return dt;
}

void FluidSimulator::draw(){
	simulationGrid->draw();
}

void FluidSimulator::simulateAndDraw(){
	project();
	advect(9.8);
	draw();
}