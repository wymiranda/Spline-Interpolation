#include <Eigen/Dense>
#include <GLFW\glfw3.h>
#include <iostream>
#include <vector>
#include "includes/splines.h"

double WIDTH = 800, HEIGHT = 600;
double XMIN = -3.5, XMAX = 3.5, YMIN = -1.5, YMAX = 1.5;
Eigen::VectorXd x(0), y(0);

//Coordenadas da tela para o plano de desenho(cartesiano)
void coordinateConversion(double& x, double& y) {
	x = x * ((XMAX - XMIN) / WIDTH) + XMIN;
	y = -y * ((YMAX - YMIN) / HEIGHT) - YMIN;
}

void setCoordenates(double xCoord, double yCoord) {
	x.conservativeResize(x.size() + 1);
	y.conservativeResize(y.size() + 1);
	coordinateConversion(xCoord, yCoord);
	x(x.size() - 1) = xCoord;
	y(y.size() - 1) = yCoord;
}

void deleteLastPoint() {
	x.conservativeResize(x.size() - 1);
	y.conservativeResize(y.size() - 1);
}

static void mouseClick(GLFWwindow* window, int button, int action, int mods)
{
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
	{
		double xmouse, ymouse;
		glfwGetCursorPos(window, &xmouse, &ymouse);
		setCoordenates(xmouse, ymouse);
	}
	if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS)
	{
		deleteLastPoint();
	}
}

int main() {

	glfwInit();
	GLFWwindow* window = glfwCreateWindow(WIDTH, HEIGHT, "Splines", NULL, NULL);
	glfwMakeContextCurrent(window);
	glClearColor(0.8f, 0.8f, 1.0f, 1.0f);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(XMIN, XMAX, YMIN, YMAX, 1, -1);

	glfwSetMouseButtonCallback(window, mouseClick);

	while (!glfwWindowShouldClose(window))
	{


		if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
			glfwSetWindowShouldClose(window, GLFW_TRUE);

		glClear(GL_COLOR_BUFFER_BIT);

		spline::drawAxis(XMIN, XMAX, YMIN, YMAX);
		//spline::quadratic::quadraticSpline(x, y);
		//spline::cubic::cubicSpline(x, y);
		//spline::parametric1::parametricSpline(x, y);
		spline::parametric::parametricSpline(x, y);
		glfwPollEvents();
		spline::drawPoints(x, y);

		glfwSwapBuffers(window);
	}
	glfwTerminate();
	return 0;
}

