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
	float xyz[3];
	Vector();
	Vector(float x, float y, float z);
	~Vector(void);
	float& operator[](int i);
	Vector normalize();
	float length();
	friend std::ostream& operator<<(std::ostream& os, Vector a){
		std::string result = "(" + std::to_string((long double)a[0]) + ", " + std::to_string((long double)a[1]) + ", " + std::to_string((long double)a[2]) + ")";
		os<<result;
		return os;
	}
	float operator*(Vector LHS);
	Vector operator^(Vector LHS);
	Vector operator*(float scalar);
	Vector operator/(float scalar);
	Vector operator+(Vector LHS);
	Vector operator-(Vector LHS);
	Vector rotateP(Vector p1, Vector p2, float angle);
	Vector rotateV(Vector aD, float thetad);
	Vector intersectRayWithTriangle(Vector ray, Vector pts[3]);
	Vector intersectRayWithPlane(Vector ray, Vector pt, Vector normal);
};

