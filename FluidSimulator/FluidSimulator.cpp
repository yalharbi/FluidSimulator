/* LOG
6/28/2015: According to (Foster and Fedkiw 2001) air cells only have velocity at their boundaries with fluid cells, also (Bridson 2008) p24
this means u(i,j) is actually the average of the two boundaries, otherwise fluids won't flow into air cells whos velocity is 0 (because we can't trace them back)
TODO: (1) Set all velocities between air cells to 0 (2) Fix getVelocity(i,j) to account for inflow and outflow (two normal velocities)

7/21/2015: Still debugging... My next step should be to implement u(i+1/2)(j) as an average, just like in the book.. currently they do not account for
neighboring cells
*/

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

double * FluidSimulator::apply(double * Adiag, double * Aplusi, double * Aplusj, double * Aprevi, double * Aprevj, double * x, int length){
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


FluidSimulator::FluidSimulator()
{
}

FluidSimulator::FluidSimulator(int width, int height, int cellSize)
{
	simulationGrid = new Grid(width, height, cellSize);
	int mid = simulationGrid->width / 2;
	for (int i = mid-5; i <mid+5; i++)
		for (int j = 10; j <25; j++)
			simulationGrid->setCell(i, j, FLUID);
	/*simulationGrid->setCell(3, 5, FLUID); simulationGrid->setCell(4, 5, FLUID); simulationGrid->setCell(2, 5, FLUID);
	simulationGrid->setCell(3, 4, FLUID); simulationGrid->setCell(4, 4, FLUID); simulationGrid->setCell(2, 4, FLUID);
	simulationGrid->setCell(3, 3, FLUID); simulationGrid->setCell(4, 3, FLUID); simulationGrid->setCell(2, 3, FLUID);*/
	for (int i = 0; i < simulationGrid->width; i++)
		simulationGrid->setCell(i, 0, SOLID);
	for (int i = 0; i < simulationGrid->height; i++){
		simulationGrid->setCell(0, i, SOLID);
		simulationGrid->setCell(simulationGrid->width - 1, i, SOLID);
		simulationGrid->setCell(simulationGrid->width - 2, i, SOLID);
	}
	simulationGrid->showVectorField();
}

void FluidSimulator::advect(double timeStep){
	dt = timeStep;
	extrapolateVelocity();
	advectVelocity();
	advectPressure();
	
}

void FluidSimulator::extrapolateVelocity(){
	for (int i = 0; i < simulationGrid->width; i++){
		for (int j = 0; j < simulationGrid->height; j++){
			if (simulationGrid->getCellType(i, j) != FLUID && simulationGrid->getCellType(i, j)!= SOLID && (simulationGrid->getCellType(i, j - 1) == FLUID ||simulationGrid->getCellType(i, j + 1) == FLUID || simulationGrid->getCellType(i - 1, j) == FLUID || simulationGrid->getCellType(i + 1, j) == FLUID)){
				std::cout << i << " " << j;
				if (simulationGrid->getCellType(i, j - 1) == FLUID){
				simulationGrid->cells[i][j+1].v = simulationGrid->cells[i][j - 1].v;

				}
				if (simulationGrid->getCellType(i, j + 1) == FLUID){
				simulationGrid->cells[i][j].v = simulationGrid->cells[i][j + 1].v;

				}
				if (simulationGrid->getCellType(i+1, j) == FLUID){
				simulationGrid->cells[i][j].u = simulationGrid->cells[i+1][j].u;

				}
				if (simulationGrid->getCellType(i-1, j) == FLUID){
					if (simulationGrid->getCellType(i + 1, j) != SOLID && simulationGrid->getCellType(i + 1, j) != FLUID)
						simulationGrid->cells[i+1][j].u = simulationGrid->cells[i-1][j].u;

				}
			}
		}
	}
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

	if (simulationGrid->getCellType(i, j) != FLUID && simulationGrid->getCellType(i, j - 1) != FLUID && simulationGrid->getCellType(i, j + 1) != FLUID && simulationGrid->getCellType(i - 1, j) != FLUID && simulationGrid->getCellType(i+1, j) != FLUID){
		simulationGrid->cells[i][j].v = 0;
		simulationGrid->cells[i][j].u = 0;
		simulationGrid->cells[i][j].cellType = FREE;
	}
	if (simulationGrid->getCellType(i, j) == FLUID ){
		simulationGrid->cells[i][j].v += (-9.8)*dt;
	}
	else {
		/*if (simulationGrid->getCellType(i, j - 1) == FLUID){
			simulationGrid->cells[i][j].v = simulationGrid->cells[i][j - 1].v;
		}
		if (simulationGrid->getCellType(i, j + 1) == FLUID){
			simulationGrid->cells[i][j].v = simulationGrid->cells[i][j + 1].v;
		}
		/*if (simulationGrid->getCellType(i+1, j) == FLUID){
			simulationGrid->cells[i][j].u = simulationGrid->cells[i+1][j].u;
		}
		if (simulationGrid->getCellType(i-1, j) == FLUID){
			simulationGrid->cells[i][j].u = simulationGrid->cells[i-1][j].u;
		}*/
	}
	



}

Cell FluidSimulator::advectCell(int i, int j){

	if (simulationGrid->cells[i][j].cellType == SOLID)
		return simulationGrid->cells[i][j];

	Vector vel(simulationGrid->getHVelocityAt(i, j), simulationGrid->getVVelocityAt(i, j), 0);
	if (vel[0] == 0 && vel[1]==0) return simulationGrid->cells[i][j];
	Vector position = simulationGrid->getCellPosition(i, j);
	Vector midPos = position - vel*dt / 2;
	Vector tracedParticle = position - simulationGrid->interpolateVelocity(midPos)*dt;

	/*std::cout << "\n" << "i" << i << " j" << j << "\n";
		std::cout << "vel" << vel << "\n";
		std::cout << "pos" << position << "\n";
		std::cout << "midpos" << midPos << "\n";
		std::cout << "tracedP" << tracedParticle << "\n";
	

	std::cout<< "intrp" << simulationGrid->interpolateVelocity(tracedParticle) << "\n";*/
	if (simulationGrid->getCell(tracedParticle).cellType == SOLID)
		return simulationGrid->cells[i][j];
	Cell cell = simulationGrid->getCell(tracedParticle);
	cell.u = simulationGrid->interpolateVelocity(tracedParticle)[0];
	cell.v = simulationGrid->interpolateVelocity(tracedParticle)[1];
	return cell;
}


void FluidSimulator::advectPressure(){

}

double * FluidSimulator::calculateNegativeDivergence(){
	double* rhs = new double[simulationGrid->getFluidCellCount()];// (double *)malloc(sizeof(double) * simulationGrid->getFluidCellCount());
	int place = 0;
	double scale = 1 / simulationGrid->getCellSize();


	for (int i = 0; i < simulationGrid->width; i++){
		for (int j = 0; j < simulationGrid->height; j++){
			if (simulationGrid->getCellType(i, j) != FLUID)
				continue;
			place = indices[std::make_pair(i, j)];
			Cell cell = simulationGrid->cells[i][j];
			rhs[place] = (simulationGrid->cells[i+1][j].u - simulationGrid->cells[i][j].u +
				simulationGrid->cells[i][j + 1].v - simulationGrid->cells[i][j].v + 0) *  -scale;


			//std::cout << rhs[place] << "\n";
		}
	}

	for (int i = 0; i < simulationGrid->width; i++){
		for (int j = 0; j < simulationGrid->height; j++){
			if (simulationGrid->getCellType(i, j) != FLUID)
				continue;
			Cell cell = simulationGrid->cells[i][j];
			int place = indices[std::make_pair(i, j)];
			if (simulationGrid->getCellType(i, j + 1) == SOLID){
				rhs[place] = rhs[place] + simulationGrid->getAveragedVVelocity(i, j + 1)[1] * scale;
			}
			if (simulationGrid->getCellType(i, j - 1) == SOLID){
				rhs[place] = rhs[place] - simulationGrid->getAveragedVVelocity(i, j)[1] * scale;
			}
			if (simulationGrid->getCellType(i + 1, j) == SOLID){
				rhs[place] = rhs[place] + simulationGrid->getAveragedHVelocity(i+1, j)[0] * scale;
			}
			if (simulationGrid->getCellType(i - 1, j) == SOLID){
				rhs[place] = rhs[place] - simulationGrid->getAveragedHVelocity(i, j)[0] * scale;
			}
			//if (rhs[place] < 1 && rhs[place] > -1) rhs[place] = 0;
		/*	std::cout << i << " " << j;
			std::cout << "^" << rhs[place] << "\n";*/
		}
	}
	return rhs;
}

double * FluidSimulator::applyPreconditioner(double * Adiag, double * Aplusi, double * Aplusj, double* Aprevi, double * Aprevj, double * precon, double * r){
	double * q = new double[simulationGrid->getFluidCellCount()]();

	for (int i = 1; i < simulationGrid->width; i++){ //Lq=r
		for (int j = 1; j < simulationGrid->height; j++){
		
			
			if (simulationGrid->getCellType(i, j) == FLUID){
				int place = indices[std::make_pair(i, j)];
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
			

			if (simulationGrid->getCellType(i, j) == FLUID){
				int place = indices[std::make_pair(i, j)];
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
	
	
	/*for (int j = simulationGrid->height - 3; j >=1; j--){
		for (int i = simulationGrid->width - 3; i >= 1; i--){
			if (simulationGrid->getCellType(i, j) == FLUID){

				std::cout << "i" << i << "j" << j << ": " << simulationGrid->cells[i-1][j].u << " " << simulationGrid->cells[i][j-1].v << " " << simulationGrid->cells[i][j].u << " " << simulationGrid->cells[i][j].v << " " << simulationGrid->cells[i + 1][j].u << " " << simulationGrid->cells[i][j + 1].v << "\n";
				
			}
		}
	}*/

	int fluidCellsCount = simulationGrid->getFluidCellCount();
	std::cout <<"N" << fluidCellsCount << "\n";

	

	int index = 0;
	for (int i = 0; i < simulationGrid->width; i++){
		for (int j = 0; j < simulationGrid->height; j++){
			if (simulationGrid->getCellType(i, j) == FLUID){
				indices.insert(std::make_pair(std::make_pair(i, j), index++));
			}
		}
	}
	std::cout << indices.size();


	pressureCoefficients = new double*[fluidCellsCount];
	double * rhs = calculateNegativeDivergence();

	for (int i = 0; i < fluidCellsCount; i++){ // return if pressure is already divergent free
		if (rhs[i] != 0)
			break;
		if (i + 1 == fluidCellsCount){
			indices.clear();
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

				precon[place] = 1 / sqrt(e+0.000001);
			}
		}
	}

	double * p = new double[fluidCellsCount]();
	double * r = copyV(rhs, fluidCellsCount);
	double * z = applyPreconditioner(Adiag, Aplusi, Aplusj, Aprevi, Aprevj, precon, r);
	double * s = copyV(z, fluidCellsCount);


	double dp = dotproduct(z, r, fluidCellsCount);
	bool done = false;


	double maxR;
	int iterations = 0;
	while (!done && iterations < MAX_ITERATIONS){
		z = apply(Adiag, Aplusi, Aplusj, Aprevi, Aprevj, s, fluidCellsCount);
		double den = dotproduct(z, s, fluidCellsCount);
		if (den == 0)
			break;
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
		z = applyPreconditioner(Adiag, Aplusi, Aplusj, Aprevi, Aprevj, precon, r);
		double dpNew = dotproduct(z, r, fluidCellsCount);
		
		double beeta = dpNew / dp;

		//std::cout << "beeta" << beeta << "dp" << dp << "dpnew" << dpNew;
		for (int i = 0; i < fluidCellsCount; i++)
			s[i] = z[i] + beeta*s[i];
		dp = dpNew;
		iterations++;
		//if (iterations + 1 == MAX_ITERATIONS || done) std::cout << "MX" << maxR << "#";
	}
	for (int i = 0; i < fluidCellsCount; i++)
		p[i] *= dt / (DENSITY * simulationGrid->getCellSize());
	double *  res = apply(Adiag, Aplusi, Aplusj, Aprevi, Aprevj, p, fluidCellsCount);
	for (int i = 0; i < simulationGrid->width; i++){ 
		for (int j = 0; j < simulationGrid->height; j++){
			if (simulationGrid->getCellType(i, j) != FLUID) continue;
			int place = indices[std::make_pair(i, j)];
			if (simulationGrid->getCellType(i, j) == FLUID && ( simulationGrid->getCellType(i - 1, j)!= FLUID || simulationGrid->getCellType(i + 1, j) != FLUID)){
				std::cout << dt / (DENSITY * simulationGrid->getCellSize()) << "p " << p[place] << " ";
				std::cout << i << " " << j << "\n";
				//std::cout << "^" << rhs[place] << "\n";
			}
		}
	}
	//std::cout << "uterations: " << iterations << "\n";

	/*for (int i = 0; i < fluidCellsCount; i++)
	std::cout << "p: " << p[i] << "\n";*/
	
	for (int i = 0; i < simulationGrid->width; i++){
		for (int j = 0; j < simulationGrid->height; j++){
			double sc = 1;
			

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

	/*for (int i = 0; i < simulationGrid->width; i++){ //compute preconditioner... initial guess for solving the system
		for (int j = 0; j < simulationGrid->height; j++){
			if (simulationGrid->getCellType(i, j) != FLUID) continue;
			int place = indices[std::make_pair(i, j)];
			if (simulationGrid->getCellType(i, j) == FLUID && (simulationGrid->getCellType(i - 1, j) != FLUID || simulationGrid->getCellType(i + 1, j) != FLUID)){
				
				//std::cout << "u " << simulationGrid->cells[i][j].u << " v" << simulationGrid->cells[i][j].v ;
				std::cout << i << " " << j << "\n";
				if (simulationGrid->getCellType(i - 1, j) == FLUID)
					std::cout << "* " << (simulationGrid->cells[i + 1][j].u + simulationGrid->cells[i][j].u) / 2 << " " << (simulationGrid->cells[i][j].v + simulationGrid->cells[i][j + 1].v) / 2 << "\n";
				if (simulationGrid->getCellType(i + 1, j) == FLUID)
					std::cout << "* " << (simulationGrid->cells[i + 1][j].u + simulationGrid->cells[i][j].u) / 2 << " " << (simulationGrid->cells[i][j].v + simulationGrid->cells[i][j+1].v) / 2 <<"\n";
				//std::cout << "^" << rhs[place] << "\n";
			}
		}
	}*/
	indices.clear();
	delete[] rhs;
	delete[] Adiag;
	delete[] Aplusi;
	delete[] precon;
	delete[] Aplusj;
	delete[] Aprevi;
	delete[] Aprevj;
	std::cout << "aft" << indices.size();
}

double FluidSimulator::computeTimeStep(){
	double cellSize = simulationGrid->getCellSize();
	double maxVel = simulationGrid->getMaxVelocity();

	dt = cellSize / maxVel;
	dt = 0.5;
	//if (dt < 0.15) dt = 0.15;
	std::cout << dt << "\n";
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