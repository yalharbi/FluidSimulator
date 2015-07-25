#pragma once
#include "Cell.h"
#include <Windows.h>
#include <gl\GL.h>
#include <glut.h>
#include <iostream>
#include "Vector.h"
#include <math.h>

class Grid
{
	double cellSize, offset;
	Vector * markers;
	bool drawVectorField;

public:
	int fluidCellCount;
	Cell ** cells;
	int width, height;
	Grid();
	Grid(int w, int h, double cellSize);
	void draw();
	void setCell(int i, int j, CellType cellType);
	Cell getCell(Vector position);
	CellType getCellType(int i, int j);
	void showVectorField();
	void hideVectorField();
	double getCellSize();
	double getMaxVelocity();
	Vector getCellPosition(int i, int j);
	double getHVelocityAt(int i, int j);
	double getVVelocityAt(int i, int j);
	Vector getVelocityVector(int i, int j);
	Vector interpolateVelocity(Vector position);
	int getFluidCellCount();
	Vector Grid::getAveragedVVelocity(int i, int j);
	Vector Grid::getAveragedHVelocity(int i, int j);
	Vector Grid::getAveragedVelocity(int i, int j);
};

