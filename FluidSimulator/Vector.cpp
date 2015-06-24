
#include "Vector.h"

float TORAD(float x){
	return x * 3.14159265 / (float)180 ;
}

Vector::Vector(){
	xyz[0] = 0;
	xyz[1] = 0;
	xyz[2] = 0;
}
Vector::Vector(float x, float y, float z)
{
	xyz[0] = x;
	xyz[1] = y;
	xyz[2] = z;
}

float& Vector::operator[](int i){
	return xyz[i];
}


float Vector::length(){
	return sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2]);
}

Vector Vector::normalize(){
	float length = this->length();
	return Vector(xyz[0]/length, xyz[1]/length, xyz[2]/length);
}

float Vector::operator*(Vector LHS){
	return xyz[0]*LHS[0]+xyz[1]*LHS[1]+xyz[2]*LHS[2];
}
Vector Vector::operator^(Vector LHS){
	float i = xyz[1]*LHS[2]-xyz[2]*LHS[1];
	float j = -(xyz[0]*LHS[2]-xyz[2]*LHS[0]);
	float k = xyz[0]*LHS[1]-xyz[1]*LHS[0];
	return Vector(i, j, k);
}
Vector Vector::operator*(float scalar){
	return Vector(xyz[0]*scalar, xyz[1]*scalar, xyz[2]*scalar);
}
Vector Vector::operator/(float scalar){
	if(scalar==0) return *this;
	return (*this)*(1/scalar);
}

Vector Vector::operator+(Vector LHS){
	Vector result(xyz[0]+LHS[0], xyz[1]+LHS[1], xyz[2]+LHS[2]);
	return result;
}

Vector Vector::operator-(Vector LHS){
	Vector result(xyz[0]-LHS[0], xyz[1]-LHS[1], xyz[2]-LHS[2]);
	return result;
}


Vector Vector::intersectRayWithPlane(Vector ray, Vector pt, Vector normal){
	float t = ((pt - ray)*normal) / (ray*normal);
	//std::cout << t << "\n";
	return ray + ray*t;
}

Vector::~Vector(void)
{
}
