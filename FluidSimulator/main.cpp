#include <Windows.h>
#include <gl\GL.h>
#include <GL\glut.h>
#define WIDTH 480
#define HEIGHT 480

void display(){
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glPointSize(5);
	glBegin(GL_POINTS);
	glColor3f(0, 1, 0);
	glVertex3f(0,0,0);
	glEnd();


	glutSwapBuffers();
	glFlush();
}

void reshape(int width, int height) {
	glViewport(0, 0, width, height);
	gluOrtho2D(-2,2,-2,2);
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
	

	glutMainLoop();

}