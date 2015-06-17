enum CellType {FREE, FLUID, SOLID};

class Cell{
	CellType cellType;
	float u, v, p;

public:
	Cell();
	Cell(CellType type, float u0, float v0, float p0);

	void setVelocity(float u, float v);
	void setPressure(float p);
};