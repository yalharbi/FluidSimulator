#pragma once
#include "Cell.h"
#include <Windows.h>
#include <gl\GL.h>
#include <glut.h>
#include <iostream>
#include "Vector.h"

class Grid
{
	float cellSize, offset;
	
	bool drawVectorField;

public:
	Cell ** cells;
	int width, height;
	Grid();
	Grid(int w, int h, float cellSize);
	void draw();
	void setCell(int i, int j, CellType cellType);
	Cell getCell(Vector position);
	void showVectorField();
	void hideVectorField();
	float getCellSize();
	float getMaxVelocity();
	Vector getCellPosition(int i, int j);
	float getHVelocityAt(int i, int j);
	float getVVelocityAt(int i, int j);
	Vector getVelocityVector(int i, int j);
	Vector interpolateVelocity(Vector position);
};

