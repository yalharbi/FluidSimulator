#pragma once
#include "Cell.h"
#include <Windows.h>
#include <gl\GL.h>
#include <glut.h>
#include <iostream>


class Grid
{
	float cellSize, offset;
	Cell ** cells;
	bool drawVectorField;

public:
	int width, height;
	Grid();
	Grid(int w, int h, float cellSize);
	void draw();
	void setCell(int i, int j, CellType cellType);
	void showVectorField();
	void hideVectorField();
	float getCellSize();
	float getMaxVelocity();
	float * getCellPosition(int i, int j);
};

