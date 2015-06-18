#include "Grid.h"
#define OFFSET 10

Grid::Grid()
{
}


Grid::Grid(int w, int h, float cSize){
	width = w/cSize;
	height = h/cSize;
	cellSize = cSize;
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

			glVertex3f(i*cellSize + OFFSET, j*cellSize + OFFSET, 0);
			glVertex3f(i*cellSize - OFFSET + cellSize, j*cellSize + OFFSET, 0);
			glVertex3f(i*cellSize - OFFSET + cellSize, j*cellSize - OFFSET + cellSize, 0);
			glVertex3f(i*cellSize + OFFSET, j*cellSize - OFFSET + cellSize, 0);
			glEnd();
		}
	}

}

void Grid::setCell(int i, int j, CellType cType){
	cells[i][j].cellType = cType;
}