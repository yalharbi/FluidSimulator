

enum CellType {FREE, FLUID, SOLID};

class Cell{
	


public:
	float u, v, p, density;
	CellType cellType;
	Cell();
	Cell(CellType type, float u0, float v0, float p0);

	void setVelocity(float u, float v);
	void setPressure(float p);
	void setType(CellType cellType);
};