#include "Grid.h"


Grid::Grid()
{
}


Grid::Grid(int w, int h, float cSize){
	width = w/cSize;
	height = h/cSize;
	cellSize = cSize;
	drawVectorField = false;
	offset = cellSize / 10;
	cells = new Cell*[width];
	for (int i = 0; i < width; i++){
		cells[i] = new Cell[height];
	}

}

void Grid::draw(){
	

	for (int i = 0; i < width; i++){
		for (int j = 0; j < height; j++){
			Cell cell = cells[i][j];
			glBegin(GL_QUADS);
			if (cell.cellType == FREE)
				glColor3f(0.5,0.5,0.5);
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
				if (cell.u>cellSize/2) cell.u = cellSize / 2;
				if (cell.v>cellSize/2) cell.v = cellSize / 2;
				glVertex3f(i*cellSize + cellSize / 2 - cell.u, j*cellSize + cellSize / 2 - cell.v, 0);
				glVertex3f(i*cellSize + cellSize / 2 + cell.u/2,  j*cellSize  + cellSize / 2 + cell.v/2, 0);
				glEnd();

				glPointSize(2);
				glBegin(GL_POINTS);
				glColor3f(1, 1, 0);
				glVertex3f(i*cellSize + cellSize / 2 + cell.u / 2, j*cellSize + cellSize / 2 + cell.v / 2, 0);
				glEnd();

			}
		}
	}

	
}

void Grid::setCell(int i, int j, CellType cType){
	cells[i][j].cellType = cType;
}

void Grid::showVectorField(){
	drawVectorField = true;
}

void Grid::hideVectorField(){
	drawVectorField = false;
}