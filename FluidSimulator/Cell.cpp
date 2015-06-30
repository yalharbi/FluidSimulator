#include "Cell.h"

Cell::Cell(){
	cellType = FREE;
	u = 0;
	v = 0;
	p = 0;
	density = 0;
}

Cell::Cell(CellType cType, float u0, float v0, float p0){
	cellType = cType;
	u = u0;
	v = v0;
	p = p0;
}

void Cell::setVelocity(float uUpdated, float vUpdated){
	u = uUpdated;
	v = vUpdated;
}

void Cell::setPressure(float pUpdated){
	p = pUpdated;
}

void Cell::setType(CellType cType){
	cellType = cType;
	if (cType == FLUID)
		density = 1000;
}
