#pragma once
#include <stdio.h>
#include <tchar.h>
#include <string>
#include <iostream>
#include <ostream>
#include <sstream>
#include <math.h>

class Vector
{
public:
	double xyz[3];
	Vector();
	Vector(double x, double y, double z);
	~Vector(void);
	double& operator[](int i);
	Vector normalize();
	double length();
	friend std::ostream& operator<<(std::ostream& os, Vector a){
		std::string result = "(" + std::to_string((long double)a[0]) + ", " + std::to_string((long double)a[1]) + ", " + std::to_string((long double)a[2]) + ")";
		os<<result;
		return os;
	}
	double operator*(Vector LHS);
	Vector operator^(Vector LHS);
	Vector operator*(double scalar);
	Vector operator/(double scalar);
	Vector operator+(Vector LHS);
	Vector operator-(Vector LHS);
	Vector rotateP(Vector p1, Vector p2, double angle);
	Vector rotateV(Vector aD, double thetad);
	Vector intersectRayWithTriangle(Vector ray, Vector pts[3]);
	Vector intersectRayWithPlane(Vector ray, Vector pt, Vector normal);
};

