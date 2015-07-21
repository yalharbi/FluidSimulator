#include "Grid.h"


Grid::Grid()
{
}


Grid::Grid(int w, int h, double cSize){
	fluidCellCount = 0;
	width = w / cSize;
	height = h / cSize;
	cellSize = cSize;
	drawVectorField = false;
	offset = cellSize / 10;
	cells = new Cell*[width];
	for (int i = 0; i < width; i++){
		cells[i] = new Cell[height];
	}


	for (int j = 0; j < height; j++){
		cells[0][j].u = 0;
		cells[width - 1][j].u = 0;
		//cells[width - 1][j].v = 0;
	}

}

void Grid::draw(){


	for (int i = 0; i < width; i++){
		for (int j = 0; j < height; j++){
			Cell cell = cells[i][j];
			glBegin(GL_QUADS);
			if (cell.cellType == FREE)
				glColor3f(0.5, 0.5, 0.5);
			else if (cell.cellType == SOLID)
				glColor3f(1, 1, 0);
			else
				glColor3f(0, 0, 1);

			glVertex3f(i*cellSize + offset, j*cellSize + offset, 0);
			glVertex3f(i*cellSize - offset + cellSize, j*cellSize + offset, 0);
			glVertex3f(i*cellSize - offset + cellSize, j*cellSize - offset + cellSize, 0);
			glVertex3f(i*cellSize + offset, j*cellSize - offset + cellSize, 0);
			glEnd();

			if (drawVectorField){
				glBegin(GL_LINES);
				glColor3f(1, 1, 1);
				Vector vv = getVelocityVector(i, j);
				double length = vv.length();
				vv = vv.normalize();
				glVertex3f(i*cellSize + cellSize / 2, j*cellSize + cellSize / 2, 0);
				glVertex3f(i*cellSize + cellSize / 2 + vv[0] * length, j*cellSize + cellSize / 2 + vv[1] * length, 0);
				glEnd();

				glPointSize(2);
				glBegin(GL_POINTS);
				glColor3f(0.75, 0.75, 0);
				glVertex3f(i*cellSize + cellSize / 2 + vv[0] * length, j*cellSize + cellSize / 2 + vv[1] * length, 0);
				glEnd();

			}
		}
	}
	float xyz[3];
	Vector pos = getCellPosition(4, 2);
	xyz[0] = pos[0] + cellSize / 2; xyz[1] = pos[1] + cellSize / 2; xyz[2] = pos[2];
	//std::cout << pos[0] << "-" << pos[1] << "\n";
	glPointSize(5);
	glBegin(GL_POINTS);
	glColor3f(1, 1, 1);
	glVertex3fv(xyz);
	glEnd();
}

void Grid::setCell(int i, int j, CellType cType){
	cells[i][j].setType(cType);
	if (cType == FLUID)
		fluidCellCount++;
}

Cell Grid::getCell(Vector pos){
	int i = pos[0] / cellSize;
	int j = pos[1] / cellSize;
	int w = 0;// pos[2] / cellSize;

	return cells[i][j];
}

CellType Grid::getCellType(int i, int j){
	return cells[i][j].cellType;
}

double Grid::getHVelocityAt(int i, int j){
	if (i == width - 1)
		return cells[i][j].u/2;
	double centerHVelocity = (cells[i][j].u + cells[i + 1][j].u) / 2;
	return centerHVelocity;
}

double Grid::getVVelocityAt(int i, int j){
	if (j == height - 1)
		return cells[i][j].v/2;
	double centerVVelocity = (cells[i][j].v + cells[i][j + 1].v) / 2;
	return centerVVelocity;
}

Vector Grid::getVelocityVector(int i, int j){
	return Vector(getHVelocityAt(i, j), getVVelocityAt(i, j), 0);
}

void Grid::showVectorField(){
	drawVectorField = true;
}

void Grid::hideVectorField(){
	drawVectorField = false;
}

double Grid::getCellSize(){
	return cellSize;
}

double Grid::getMaxVelocity(){
	// TODO
	double max = 0;
	for (int i = 0; i < width - 1; i++){
		for (int j = 0; j < height - 1; j++){
			//if (getCellType(i, j) != FLUID) continue;
			if (abs(cells[i][j].u) > max)
				max = abs(cells[i][j].u);
			if (abs(cells[i][j].v) > max)
				max = abs(cells[i][j].v);
		}
	}

	return max + sqrt(5 * sqrt(cellSize)*9.8);
}

Vector Grid::interpolateVelocity(Vector position){


	int x0 = position[0] / cellSize;
	int x1 = x0 + 1;
	int y0 = position[1] / cellSize;
	int y1 = y0 + 1;

	if (x1 >= width - 1 || x0 <= 0 || y1 >= height - 1 || y0 <= 0)
		return Vector(0, 0, 0);

	double alpha = (position[0] - x0*cellSize) / cellSize;
	double approximateU = (1 - alpha)*cells[x0][y0].u + (alpha)*cells[x1][y0].u;

	double beta = (position[1] - y0*cellSize) / cellSize;
	double approximateV = (1 - beta)*cells[x0][y0].v + (beta)*cells[x0][y1].v;

	Vector x0y0(cells[x0][y0].u, cells[x0][y0].v, 0);
	Vector x1y0(cells[x1][y0].u, cells[x1][y0].v, 0);
	Vector x0y1(cells[x0][y1].u, cells[x0][y1].v, 0);
	Vector x1y1(cells[x1][y1].u, cells[x1][y1].v, 0);

	double den = (y1*cellSize - y0*cellSize)*(x1*cellSize - x0*cellSize);
	Vector bilinear = x0y0*(x1*cellSize - position[0])*(y1*cellSize - position[1])
		+ x1y0*(position[0] - x0*cellSize)*(y1*cellSize - position[1])
		+ x0y1*(x1*cellSize - position[0])*(position[1] - y0*cellSize)
		+ x1y1*(position[0] - x0*cellSize)*(position[1] - y0*cellSize);
	bilinear = bilinear / den;
	//return getVelocityVector(x0, y0);
	return Vector(approximateU, approximateV, 0);
	return bilinear;
}

Vector Grid::getCellPosition(int i, int j){
	double xyz[3];

	xyz[0] = i*cellSize + cellSize / 2;
	xyz[1] = j*cellSize + cellSize / 2;
	xyz[2] = 0; //3D

	return Vector(xyz[0], xyz[1], xyz[2]);
}

int Grid::getFluidCellCount(){
	fluidCellCount = 0;
	for (int i = 0; i < width; i++){
		for (int j = 0; j < height; j++){
			if (getCellType(i, j) != FLUID) continue;
			fluidCellCount++;
		}
	}
	return fluidCellCount;
}