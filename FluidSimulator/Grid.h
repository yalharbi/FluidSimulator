#pragma once
#include "Cell.h"
#include <Windows.h>
#include <gl\GL.h>
#include <glut.h>
#include <iostream>


class Grid
{
	int width, height;
	float cellSize;
	Cell ** cells;

public:
	Grid();
	Grid(int w, int h, float cellSize);
	void draw();
	void setCell(int i, int j, CellType cellType);
};
