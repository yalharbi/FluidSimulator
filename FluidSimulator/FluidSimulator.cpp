#include "FluidSimulator.h"
#define DENSITY 1000
#define MAX_ITERATIONS 100
#define TOLERANCE 1/1000000

#define square(x) x*x
#define getIndex(x,y) indices[std::make_pair(x,y)]

inline double dotproduct(double * a, double * b, int length){
	double result = 0;
	for (int it = 0; it < length; it++){
		result += a[it] * b[it];
	}

	return result;
}

double * FluidSimulator::apply(double * Adiag, double * Aplusi, double * Aplusj, double * Aprevi, double * Aprevj, double * x, int length, std::map<std::pair<int, int>, int> indices){
	double * ans = new double[length]();

	for (int i = 0; i < simulationGrid->width; i++){
		for (int j = 0; j < simulationGrid->height; j++){
			if (simulationGrid->getCellType(i, j) == FLUID){
				int place = indices[std::make_pair(i, j)];

				ans[place] = Adiag[place] * x[place];
				if (indices.count(std::make_pair(i+1, j)) > 0)
					ans[place] += Aplusi[place] * x[getIndex(i + 1, j)];
				if (indices.count(std::make_pair(i, j + 1)) > 0)
					ans[place] += Aplusj[place] * x[getIndex(i, j + 1)];
				if (indices.count(std::make_pair(i-1, j)) > 0)
					ans[place] += Aprevi[place] * x[getIndex(i - 1, j)];
				if (indices.count(std::make_pair(i, j - 1)) > 0)
					ans[place] += Aprevj[place] * x[getIndex(i, j - 1)];

				/*if (simulationGrid->getCellType(i - 1, j) == FLUID)
				ans[place] += Aplusi[placepi] * x[placepi];
				if (simulationGrid->getCellType(i, j - 1) == FLUID)
				ans[place] += Aplusj[placepj] * x[placepj];*/
			}
		}
	}


	return ans;
}

inline double * copyV(double * from, int length){
	double * to = new double[length]();
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
	int mid = simulationGrid->width / 2;
	for (int i = mid-5; i < mid+5; i++)
		for (int j = 11; j <25; j++)
			simulationGrid->setCell(i, j, FLUID);
	/*simulationGrid->setCell(3, 5, FLUID); simulationGrid->setCell(4, 5, FLUID); simulationGrid->setCell(2, 5, FLUID);
	simulationGrid->setCell(3, 4, FLUID); simulationGrid->setCell(4, 4, FLUID); simulationGrid->setCell(2, 4, FLUID);
	simulationGrid->setCell(3, 3, FLUID); simulationGrid->setCell(4, 3, FLUID); simulationGrid->setCell(2, 3, FLUID);*/
	for (int i = 0; i < simulationGrid->width; i++)
		simulationGrid->setCell(i, 0, SOLID);
	for (int i = 0; i < simulationGrid->height; i++){
		simulationGrid->setCell(0, i, SOLID);
		simulationGrid->setCell(simulationGrid->width - 1, i, SOLID);
	}
	//simulationGrid->showVectorField();
}

void FluidSimulator::advect(double timeStep){
	dt = timeStep;
	extrapolateVelocity();
	advectVelocity();
	advectPressure();
	
}

void FluidSimulator::extrapolateVelocity(){

}

void FluidSimulator::advectVelocity(){


	Cell ** updatedCells = new Cell*[simulationGrid->width];
	for (int i = 0; i < simulationGrid->width; i++){
		updatedCells[i] = new Cell[simulationGrid->height];
	}


	computeTimeStep();

	for (int i = 0; i < simulationGrid->width; i++){
		for (int j = 0; j < simulationGrid->height; j++){
			updatedCells[i][j] = advectCell(i, j);
		}
	}
	simulationGrid->cells = updatedCells;
	for (int i = 0; i < simulationGrid->width; i++){
		for (int j = 0; j < simulationGrid->height; j++){
			addForces(i, j);
		}
	}


	//std::cout << simulationGrid->getVelocityVector(4, 2);
}

void FluidSimulator::addForces(int i, int j){
	if (i <= 0 || j <= 0 || i >= simulationGrid->width - 1 || j >= simulationGrid->width - 1)
		return;
	if (simulationGrid->getCellType(i, j) == SOLID) return;
	if (simulationGrid->getCellType(i, j) != FLUID  && simulationGrid->getCellType(i, j + 1) != FLUID && simulationGrid->getCellType(i, j - 1) != FLUID && simulationGrid->getCellType(i-1, j) != FLUID && simulationGrid->getCellType(i+1, j) != FLUID){
		simulationGrid->cells[i][j].v = 0;
		simulationGrid->cells[i][j].u = 0;
		simulationGrid->cells[i][j].cellType = FREE;
	}
	if (simulationGrid->getCellType(i, j) != SOLID){
		simulationGrid->cells[i][j].v += (-9.8)*dt;
	}
	
	



}

Cell FluidSimulator::advectCell(int i, int j){

	if (simulationGrid->cells[i][j].cellType == SOLID)
		return simulationGrid->cells[i][j];

	Vector vel = Vector((simulationGrid->cells[i][j].u + simulationGrid->cells[i + 1][j].u) / 2.0, (simulationGrid->cells[i][j].v + simulationGrid->cells[i][j + 1].v) / 2.0, 0);
	Vector position = simulationGrid->getCellPosition(i, j);
	Vector midPos = position - vel*dt / 2;
	Vector tracedParticle = position - simulationGrid->interpolateVelocity(midPos)*dt;

	
	simulationGrid->interpolateVelocity(tracedParticle);
	return simulationGrid->getCell(tracedParticle);
}


void FluidSimulator::advectPressure(){

}

double * FluidSimulator::calculateNegativeDivergence(std::map<std::pair<int, int>, int> indices){
	double* rhs = new double[simulationGrid->getFluidCellCount()];// (double *)malloc(sizeof(double) * simulationGrid->getFluidCellCount());
	int place = 0;
	double scale = 1 / simulationGrid->getCellSize();


	for (int i = 0; i < simulationGrid->width; i++){
		for (int j = 0; j < simulationGrid->height; j++){
			if (simulationGrid->getCellType(i, j) != FLUID)
				continue;
			place = indices[std::make_pair(i, j)];
			Cell cell = simulationGrid->cells[i][j];
			rhs[place] = (simulationGrid->cells[i + 1][j].u - simulationGrid->cells[i][j].u +
				simulationGrid->cells[i][j + 1].v - simulationGrid->cells[i][j].v + 0) *  -scale;


			//std::cout << rhs[place - 1] << "\n";
		}
	}

	for (int i = 0; i < simulationGrid->width; i++){
		for (int j = 0; j < simulationGrid->height; j++){
			if (simulationGrid->getCellType(i, j) != FLUID)
				continue;
			Cell cell = simulationGrid->cells[i][j];
			int place = indices[std::make_pair(i, j)];
			if (simulationGrid->getCellType(i, j + 1) == SOLID){
				rhs[place] = rhs[place] + simulationGrid->cells[i][j+1].v*scale;
			}
			if (simulationGrid->getCellType(i, j - 1) == SOLID){
				rhs[place] = rhs[place] - simulationGrid->cells[i][j].v*scale;
			}
			if (simulationGrid->getCellType(i + 1, j) == SOLID){
				rhs[place] = rhs[place] +simulationGrid->cells[i + 1][j].u*scale;
			}
			if (simulationGrid->getCellType(i - 1, j) == SOLID){
				rhs[place] = rhs[place] - simulationGrid->cells[i][j].u*scale;
			}

			//std::cout << "^" << rhs[place] << "\n";
		}
	}
	return rhs;
}

double * FluidSimulator::applyPreconditioner(double * Adiag, double * Aplusi, double * Aplusj, double* Aprevi, double * Aprevj, double * precon, double * r, std::map<std::pair<int, int>, int> indices){
	double * q = new double[simulationGrid->getFluidCellCount()]();

	for (int i = 1; i < simulationGrid->width; i++){ //Lq=r
		for (int j = 1; j < simulationGrid->height; j++){
			int place = indices[std::make_pair(i, j)];
			
			if (simulationGrid->getCellType(i, j) == FLUID){
				//std::cout << "p" << place << "placepj" << placepj << "placepi" << placepi << "\n";
				double t = r[place];
				if (indices.count(std::make_pair(i-1, j)) > 0)
					t -= Aprevi[place] * precon[getIndex(i - 1, j)] * q[getIndex(i - 1, j)];
				if (indices.count(std::make_pair(i, j - 1)) > 0)
					t -= Aprevj[place] * precon[getIndex(i, j - 1)] * q[getIndex(i, j - 1)];
				q[place] = t*precon[place];
			}
		}
	}

	double * z = new double[simulationGrid->getFluidCellCount()]();
	for (int i = simulationGrid->width - 1; i >= 0; i--){ //Ltz = q
		for (int j = simulationGrid->height - 1; j >= 0; j--){
			int place = indices[std::make_pair(i, j)];

			if (simulationGrid->getCellType(i, j) == FLUID){
				double t = q[place];
				if (indices.count(std::make_pair(i+1, j)) > 0)
					t -= Aplusi[place] * precon[place] * z[getIndex(i + 1, j)];
				if (indices.count(std::make_pair(i, j + 1)) > 0)
					t -= Aplusj[place] * precon[place] * z[getIndex(i, j + 1)];
				z[place] = t*precon[place];
			}
		}
	}
	return z;
}

void FluidSimulator::project(){
	
	int fluidCellsCount = simulationGrid->getFluidCellCount();
	std::cout <<"N" << fluidCellsCount << "\n";

	std::map<std::pair<int, int>, int> indices;

	int index = 0;
	for (int i = 0; i < simulationGrid->width; i++){
		for (int j = 0; j < simulationGrid->height; j++){
			if (simulationGrid->getCellType(i, j) == FLUID){
				indices.insert(std::make_pair(std::make_pair(i, j), index++));

			}
		}
	}


	pressureCoefficients = new double*[fluidCellsCount];
	double * rhs = calculateNegativeDivergence(indices);

	for (int i = 0; i < fluidCellsCount; i++){ // return if pressure is already divergent free
		if (rhs[i] != 0)
			break;
		if (i + 1 == fluidCellsCount){
			std::cout << "divFree\n";
			return;
		}
	}
	double * Adiag = new double[fluidCellsCount]();
	double * Aplusi = new double[fluidCellsCount]();
	double * precon = new double[fluidCellsCount]();
	double * Aplusj = new double[fluidCellsCount]();
	double * Aprevi = new double[fluidCellsCount]();
	double * Aprevj = new double[fluidCellsCount]();

	double scale = dt / (DENSITY*simulationGrid->getCellSize()*simulationGrid->getCellSize());
	for (int i = 0; i < simulationGrid->width; i++){ //compute A... where A*pressure(unknown) = negative divergence... pressure makes the field divergent free
		for (int j = 0; j < simulationGrid->height; j++){
			if (simulationGrid->getCellType(i, j) != FLUID)		continue;

			int place = indices[std::make_pair(i, j)];
			if (simulationGrid->getCellType(i, j) == FLUID && simulationGrid->getCellType(i + 1, j) == FLUID){
				Adiag[place] += scale;
				Aplusi[place] = -scale;
			}
			else if (simulationGrid->getCellType(i, j) == FLUID && simulationGrid->getCellType(i + 1, j) == FREE){
				Adiag[place] += scale;
			}

			if (simulationGrid->getCellType(i, j) == FLUID && simulationGrid->getCellType(i, j + 1) == FLUID){
				Adiag[place] += scale;
				Aplusj[place] = -scale;
			}
			else if (simulationGrid->getCellType(i, j) == FLUID && simulationGrid->getCellType(i, j + 1) == FREE)
				Adiag[place] += scale;

			if (simulationGrid->getCellType(i, j) == FLUID && simulationGrid->getCellType(i - 1, j) == FLUID){
				Adiag[place] += scale;
				Aprevi[place] = -scale;
			}
			else if (simulationGrid->getCellType(i, j) == FLUID && simulationGrid->getCellType(i - 1, j) == FREE){
				Adiag[place] += scale;

			}
			if (simulationGrid->getCellType(i, j) == FLUID && simulationGrid->getCellType(i, j-1) == FLUID){
				Adiag[place] += scale;
				Aprevj[place] = -scale;
			}
			else if (simulationGrid->getCellType(i, j) == FLUID && simulationGrid->getCellType(i, j-1) == FREE){
				Adiag[place] += scale;
			}

		/*	if (simulationGrid->getCellType(i, j) == FREE && simulationGrid->getCellType(i, j + 1) == FLUID)
				Adiag[indices[std::make_pair(i, j + 1)]] += scale;
			if (simulationGrid->getCellType(i, j) == FREE && simulationGrid->getCellType(i + 1, j) == FLUID)
				Adiag[indices[std::make_pair(i + 1, j)]] += scale;*/

		}
	}

	/*for (int i = 15; i < simulationGrid->width; i=i+9){
		for (int j = 0; j < simulationGrid->height; j++){
			if (simulationGrid->getCellType(i, j) != FLUID) continue;
			int place = indices[std::make_pair(i, j)];


			std::cout << "***place:" << i << " " << j << " " << place << " Adiag: " << Adiag[place] / scale << "\n";
			std::cout << " Aprevi: " << Aprevi[place] / scale << " " << getIndex(i-1, j) << " " << indices.count(std::make_pair(i-1, j)) << "\n";
			std::cout << " Aprevj: " << Aprevj[place] / scale << " " << getIndex(i, j-1) << " " << indices.count(std::make_pair(i, j - 1)) << "\n";
			std::cout << " Aplusi: " << Aplusi[place] / scale << " " << getIndex(i + 1, j) << " " << indices.count(std::make_pair(i + 1, j)) << "\n";
			std::cout << " Aplusj: " << Aplusj[place] / scale << " " << getIndex(i, j+1) << " " << indices.count(std::make_pair(i, j + 1)) << "\n";
		}
	}*/
	double tau = 0.97, beta = 0.25;
	for (int i = 0; i < simulationGrid->width; i++){ //compute preconditioner... initial guess for solving the system
		for (int j = 0; j < simulationGrid->height; j++){
			
			

			if (simulationGrid->getCellType(i, j) == FLUID){
				int place = indices[std::make_pair(i, j)];

				double e = Adiag[place];
				if (indices.count(std::make_pair(i - 1, j)) > 0)
					e -= square(Aprevi[place] * precon[getIndex(i-1,j)]);
				if (indices.count(std::make_pair(i, j - 1)) > 0)
					e -= square(Aprevj[place] * precon[getIndex(i, j - 1)]);
				if (indices.count(std::make_pair(i - 1, j)) > 0)
					e -= tau*(Aplusi[getIndex(i - 1, j)] * Aplusj[getIndex(i - 1, j)] * square(precon[getIndex(i - 1, j)]));
				if (indices.count(std::make_pair(i, j - 1)) > 0)
					e -= tau*(Aplusj[getIndex(i, j - 1)] * Aplusi[getIndex(i, j - 1)] * square(precon[getIndex(i, j - 1)]));

				if (e < beta * Adiag[place])
					e = Adiag[place];

				precon[place] = 1 / sqrt(e);
			}
		}
	}

	double * p = new double[fluidCellsCount]();
	double * r = copyV(rhs, fluidCellsCount);
	double * z = applyPreconditioner(Adiag, Aplusi, Aplusj, Aprevi, Aprevj, precon, r, indices);
	double * s = copyV(z, fluidCellsCount);


	double dp = dotproduct(z, r, fluidCellsCount);
	bool done = false;


	double maxR;
	int iterations = 0;
	while (!done && iterations < MAX_ITERATIONS){
		z = apply(Adiag, Aplusi, Aplusj, Aprevi, Aprevj, s, fluidCellsCount, indices);
		double den = dotproduct(z, s, fluidCellsCount);

		double alph = dp / den;

		maxR = 0;

		for (int i = 0; i < fluidCellsCount; i++){
			p[i] = p[i] + alph*s[i];
			r[i] = r[i] - alph*z[i];
			if (r[i] > maxR)
				maxR = r[i];
		}
		if (abs(maxR) <= TOLERANCE){
			done = true;
			//std::cout << "MX" << maxR << "#";
			break;
		}
		z = applyPreconditioner(Adiag, Aplusi, Aplusj, Aprevi, Aprevj, precon, r, indices);
		double dpNew = dotproduct(z, r, fluidCellsCount);
		
		double beeta = dpNew / dp;

		//std::cout << "beeta" << beeta << "dp" << dp << "dpnew" << dpNew;
		for (int i = 0; i < fluidCellsCount; i++)
			s[i] = z[i] + beeta*s[i];
		dp = dpNew;
		iterations++;
		if (iterations + 1 == MAX_ITERATIONS || done) std::cout << "MX" << maxR << "#";
	}
	double *  res = apply(Adiag, Aplusi, Aplusj, Aprevi, Aprevj, p, fluidCellsCount, indices);
	for (int i = 0; i < fluidCellsCount; i++){
		//std::cout << "MX" << maxR << "#";
		std::cout << "p" << p[i] << "\n";
	}
	//std::cout << "uterations: " << iterations << "\n";

	/*for (int i = 0; i < fluidCellsCount; i++)
	std::cout << "p: " << p[i] << "\n";*/
	
	for (int i = 0; i < simulationGrid->width; i++){
		for (int j = 0; j < simulationGrid->height; j++){
			double sc = dt / (DENSITY * simulationGrid->getCellSize());
			

			if (simulationGrid->getCellType(i, j) == FLUID){
				int place = indices[std::make_pair(i, j)];
				simulationGrid->cells[i][j].u -= sc*p[place];
				simulationGrid->cells[i + 1][j].u += sc*p[place];

				simulationGrid->cells[i][j].v -= sc*p[place];
				simulationGrid->cells[i][j + 1].v += sc*p[place];
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

	delete[] rhs;
	delete[] Adiag;
	delete[] Aplusi;
	delete[] precon;
	delete[] Aplusj;
	delete[] Aprevi;
	delete[] Aprevj;

}

double FluidSimulator::computeTimeStep(){
	double cellSize = 5*sqrt(simulationGrid->getCellSize());
	double maxVel = simulationGrid->getMaxVelocity();

	dt = cellSize / maxVel;
	dt = 0.5;
	//std::cout << dt << "\n";
	return dt;
}

void FluidSimulator::draw(){
	simulationGrid->draw();
}

void FluidSimulator::simulateAndDraw(){
	advect(0.01);
	project();
	

	draw();
}