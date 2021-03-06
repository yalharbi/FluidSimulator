#pragma once
#include <Windows.h>
#include <gl\GL.h>
#include <glut.h>
#include "Grid.h"
#include "FluidSimulator.h"

#define WIDTH 400
#define HEIGHT 400

FluidSimulator *simulator;
Grid* grid;
int time = 0;;
void display(){
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	simulator->draw();


	if (time==0){
		simulator->simulateAndDraw();
		//time = 1;
	}




	glutSwapBuffers();
	glFlush();
}

void reshape(int width, int height) {
	simulator = new FluidSimulator(width, height,10);
	/*grid = new Grid(width, height, 100);
	grid->showVectorField();
	/*for (int i = 0; i < width; i++){
		for (int j = 0; j < height; j++){
			if (i == 0 || i == width - 1 || j == 0 || j == height - 1)
				grid->setCell(i, j, SOLID);
		}
	}
	grid->setCell(2, 1, FLUID);*/
	glViewport(0, 0, width, height);
	gluOrtho2D(0,width,0, height);
}

void idle(){
	glutPostRedisplay();
}

void key(unsigned char key, int x, int y){
	time--;
	//glutPostRedisplay();
}

int main(int argc, char** argv){
	//initialization
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowPosition(300, 300);
	glutInitWindowSize(WIDTH, HEIGHT);
	glutCreateWindow("Grid-Based Fluid Simulator - Yazeed Alharbi - Professor Tricoche - Summer 2015");
	glutReshapeFunc(reshape);
	glutDisplayFunc(display);
	glutIdleFunc(idle);
	glutKeyboardFunc(key);

	glutMainLoop();

}