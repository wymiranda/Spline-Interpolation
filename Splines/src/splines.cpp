#include <iostream>
#include <GLFW\glfw3.h>
#include <Eigen\Dense>
#include "includes/splines.h"

namespace spline
{
	void drawPoints(Eigen::MatrixXd points) {
		glColor3f(0.0f, 0.0f, 0.0f);
		glPointSize(5.0f);

		glBegin(GL_POINTS);
		for (int i = 0; i < points.rows(); i++) glVertex3d(points(i, 0), points(i, 1), 0);
		glEnd();
	}

	void drawPoints(Eigen::VectorXd x, Eigen::VectorXd y) {
		if (x.size() != y.size())
			std::cout << "Os vetores X e Y não devem ter tamanhos diferentes!" << std::endl;

		glColor3f(0.f, 0.f, 0.f);
		glPointSize(5.f);

		glBegin(GL_POINTS);
		for (int i = 0; i < x.size(); i++) glVertex3d(x(i), y(i), 0);
		glEnd();
	}

	void drawAxis(double xmin, double xmax, double ymin, double ymax) {
		glColor3f(0.f, 0.f, 1.f);
		glLineWidth(0.5f);

		glBegin(GL_LINES);
		glVertex3d(xmin, 0, 0);
		glVertex3d(xmax, 0, 0);
		glVertex3d(0, ymin, 0);
		glVertex3d(0, ymax, 0);
		glEnd();
	}

	void drawFunction(function f, double xa, double xb) {
		glColor3f(1.f, 1.f, 0.f);
		glLineWidth(5.0f);

		int points = 100;
		double increment = (xa - xb) / ((double)points - 1);
		double num = xa;

		glBegin(GL_LINE_STRIP);
		for (int i = 0; i < points; i++) {
			glVertex3d(num, f(num), 0);
			num += increment;
		}
		glEnd();
	}

	namespace quadratic {

		uint64_t N;							//Número de equações
		Eigen::VectorXd xCoord, yCoord;		// Coordenadas x e y dos pontos de entrada
		Eigen::MatrixXd A;					// AX=B
		Eigen::VectorXd X, B;				// AX=B

		Eigen::MatrixXd coefficients;		// Coeficientes das parábolas encontradas

		void setPointCoord(Eigen::VectorXd x, Eigen::VectorXd y) {
			N = x.size() - 1;
			xCoord = x;
			yCoord = y;
		}

		void setCoefficients(Eigen::VectorXd coeff) {
			coefficients = Eigen::MatrixXd::Zero(N, 3);

			for (uint64_t i = 0; i < N; i++) {
				coefficients.row(i) = coeff.segment(i * 3, 3);
			}
		}

		void setMatrixA(Eigen::VectorXd x, Eigen::MatrixXd& out) {
			out = Eigen::MatrixXd::Zero(3 * N, 3 * N);

			for (uint64_t n = 0; n < N; n++) {

				out.block<1, 3>(2 * n, 3 * n) << x(n) * x(n), x(n), 1;
				out.block<1, 3>(2 * n + 1, 3 * n) << x(n + 1) * x(n + 1), x(n + 1), 1;
			}

			for (uint64_t n = 0; n < N - 1; n++)
				out.block<1, 6>(2 * N + n, 3 * n) << 2 * x(n + 1), 1, 0, -2 * x(n + 1), -1, 0;

			out(3 * N - 1, 0) = 1;

		}

		void setVectorB(Eigen::VectorXd y, Eigen::VectorXd& out) {
			out = Eigen::VectorXd::Zero(3 * N);

			for (uint64_t i = 0; i < N; i++)
				out.segment(2 * i, 2) << y(i), y(i + 1);
		}

		void solveEquations(Eigen::MatrixXd A, Eigen::VectorXd& X, Eigen::VectorXd B) {
			if (A.rows() != B.size())
				std::cout << "A matriz de A e o vetor B não devem ter quantidade de linhas diferentes";

			X = A.colPivHouseholderQr().solve(B);
		}

		void quadraticSpline(Eigen::VectorXd x, Eigen::VectorXd y) {

			if (x.size() < 3) return;

			if (x.size() != y.size())
				std::cout << "O tamanho de X e Y devem ser iguais" << std::endl;


			setPointCoord(x, y);
			setMatrixA(xCoord, A);
			setVectorB(yCoord, B);
			solveEquations(A, X, B);
			setCoefficients(X);
			drawSpline(xCoord, coefficients);

		}

		double quadFunction(double x, double a, double b, double c) {
			return x * x * a + x * b + c;
		}

		void drawSpline(Eigen::VectorXd x, Eigen::MatrixXd coeff) {

			glColor3f(0.f, 1.f, 0.f);
			glLineWidth(2.0f);

			int sections = 100;
			double increment;
			double xa, xb;

			for (uint64_t i = 0; i < N; i++) {

				xa = x(i);
				xb = x(i + 1);
				increment = (double)(xb - xa) / sections;

				double num = xa;

				glBegin(GL_LINE_STRIP);
				for (int j = 0; j <= sections; j++) {
					glVertex3f(num, (float)quadFunction(num, coeff(i, 0), coeff(i, 1), coeff(i, 2)), 0.f);
					num += increment;
				}
				glEnd();
			}
		}
	}


	namespace cubic {

		uint64_t N;							//Número de equações
		Eigen::VectorXd xCoord, yCoord;		// Coordenadas x e y dos pontos de entrada
		Eigen::MatrixXd A;					// AX=B
		Eigen::VectorXd X, B;				// AX=B
		Eigen::VectorXd h;					// h é o vetor que contém "X(i) - X(i-i)"

		Eigen::MatrixXd coefficients;		// Coeficientes dos polinômios

		void setPointCoord(Eigen::VectorXd x, Eigen::VectorXd y) {
			N = x.size() - 1;
			xCoord = x;
			yCoord = y;
		}

		void setDifferenceH(Eigen::VectorXd& x, Eigen::VectorXd& out) {
			out = Eigen::VectorXd::Zero(N);

			for (uint64_t i = 0; i < N; i++)
				out(i) = (x(i + 1) - x(i));
		}

		void setMatrixA(Eigen::VectorXd h, Eigen::MatrixXd& out) {
			out = Eigen::MatrixXd::Zero(N - 1, N - 1);

			out.block<1, 2>(0, 0) << 2 * (h(0) + h(1)), h(1);

			for (uint64_t i = 1; i < N - 2; i++)
				out.block<1, 3>(i, i - 1) << h(i), 2 * (h(i) + h(i + 1)), h(i + 1);

			out.block<1, 2>(N - 2, N - 3) << h(N - 2), 2 * (h(N - 2) + h(N - 1));
		}

		void setVectorB(Eigen::VectorXd h, Eigen::VectorXd y, Eigen::VectorXd& out) {
			out = Eigen::VectorXd::Zero(N - 1);

			for (uint64_t i = 0; i < N - 1; i++)
				out(i) = (y(i + 2) - y(i + 1)) / h(i + 1) - (y(i + 1) - y(i)) / h(i);

			out *= 6;
		}

		void setCoefficients(Eigen::VectorXd h, Eigen::VectorXd y, Eigen::VectorXd X, Eigen::MatrixXd& coeff) {
			coeff = Eigen::MatrixXd::Zero(N, 4);

			for (uint64_t i = 0; i < N; i++) {
				coeff(i, 0) = (X(i + 1) - X(i)) / (6 * h(i));
				coeff(i, 1) = X(i + 1) / 2;
				coeff(i, 2) = (y(i + 1) - y(i)) / (h(i)) + h(i) * (2 * X(i + 1) + X(i)) / 6;
				coeff(i, 3) = y(i + 1);
			}

		}

		void solveEquations(Eigen::MatrixXd A, Eigen::VectorXd& X, Eigen::VectorXd B) {
			X = Eigen::VectorXd::Zero(N + 1);

			if (A.rows() != B.size())
				std::cout << "A matriz de A e o vetor B não devem ter quantidade de linhas diferentes";

			X.segment(1, N - 1) = A.colPivHouseholderQr().solve(B);
		}

		void cubicSpline(Eigen::VectorXd x, Eigen::VectorXd y) {

			if (x.size() < 4) return;

			if (xCoord.size() != yCoord.size())
				std::cout << "O tamanho dos vetores X e Y devem ser iguais" << std::endl;

			setPointCoord(x, y);
			setDifferenceH(xCoord, h);
			setMatrixA(h, A);
			setVectorB(h, yCoord, B);
			solveEquations(A, X, B);
			setCoefficients(h, y, X, coefficients);
			drawSpline(xCoord, coefficients);
		}

		double cubicFunction(double x, double xk, double a, double b, double c, double d) {
			return pow(x - xk, 3) * a + pow(x - xk, 2) * b + (x - xk) * c + d;
		}

		void drawSpline(Eigen::VectorXd x, Eigen::MatrixXd coeff) {

			glColor3f(1.f, 0.f, 0.f);
			glLineWidth(2.0f);

			int sections = 40;
			double increment;
			double xa, xb;

			for (uint64_t i = 0; i < x.size() - 1; i++) {

				xa = x(i);
				xb = x(i + 1);
				increment = (xb - xa) / sections;

				double num = xa;
				for (int j = 0; j < sections; j++) {

					glBegin(GL_LINES);
					glVertex3f(num, cubicFunction(num, xb, coeff(i, 0), coeff(i, 1), coeff(i, 2), coeff(i, 3)), 0.f);
					glVertex3f(num + increment, cubicFunction(num + increment, xb, coeff(i, 0), coeff(i, 1), coeff(i, 2), coeff(i, 3)), 0.f);
					glEnd();

					num += increment;
				}
			}
		}
	}
	/*
	namespace parametric1 {
		//A(0) = 0; B(0) = 0
		uint64_t N;							//Número de equações
		Eigen::VectorXd x, y;				// Coordenadas x e y dos pontos de entrada
		Eigen::VectorXd h;					// h é o vetor que contém "t(i) - t(i-1)"

		Eigen::MatrixXd xCoefficients, yCoefficients;		// Coeficientes dos polinômios


		void setPointCoord(Eigen::VectorXd x_, Eigen::VectorXd y_) {
			N = x_.size() - 1;
			x = x_;
			y = y_;
		}

		void setDifferenceH(double interval_t) {
			//Considerando igualmente espaçado, tal que:
			//t(i) = i/N onde (i = 0, 1, 2, 3, ..., N)
			// Deste modo 0 < t < 1 sempre será verdadeiro, e "t(i) - t(i-1)" = 1

			h = Eigen::VectorXd::Constant(N, 1, interval_t);
		}

		void cubicFunction(double t, double tk, int line, double &x, double &y) {
			x =  pow(t - tk, 3) * xCoefficients(line, 0) + pow(t - tk, 2) * xCoefficients(line, 1) + (t - tk) * xCoefficients(line, 2) + xCoefficients(line, 3);
			y = pow(t - tk, 3) * yCoefficients(line, 0) + pow(t - tk, 2) * yCoefficients(line, 1) + (t - tk) * yCoefficients(line, 2) + yCoefficients(line, 3);
		}

		void calculeCoefficients() {

			xCoefficients = Eigen::MatrixXd::Zero(N, 4);
			yCoefficients = Eigen::MatrixXd::Zero(N, 4);

			xCoefficients(0, 2) = (x(1) - x(0)) / h(0);
			yCoefficients(0, 2) = (y(1) - y(0)) / h(0);
			xCoefficients(0, 3) = x(0);
			yCoefficients(0, 3) = y(0);

			for (int i = 1; i < N; i++) {
				xCoefficients(i, 3) = x(i);
				yCoefficients(i, 3) = y(i);

				xCoefficients(i, 1) = (3 * xCoefficients(i - 1, 0) * h(i - 1) + xCoefficients(i - 1, 1));
				yCoefficients(i, 1) = (3 * yCoefficients(i - 1, 0) * h(i - 1) + yCoefficients(i - 1, 1));

				xCoefficients(i, 2) = 3 * xCoefficients(i - 1, 0) * pow(h(i - 1), 2) + 2 * xCoefficients(i - 1, 1) * h(i - 1) + xCoefficients(i - 1, 2);
				yCoefficients(i, 2) = 3 * yCoefficients(i - 1, 0) * pow(h(i - 1), 2) + 2 * yCoefficients(i - 1, 1) * h(i - 1) + yCoefficients(i - 1, 2);

				xCoefficients(i, 0) = (x(i + 1) - xCoefficients(i, 1) * pow(h(i), 2) - xCoefficients(i, 2) * h(i) - xCoefficients(i, 3)) / pow(h(i), 3);
				yCoefficients(i, 0) = (y(i + 1) - yCoefficients(i, 1) * pow(h(i), 2) - yCoefficients(i, 2) * h(i) - yCoefficients(i, 3)) / pow(h(i), 3);
			}

		}

		void drawSplines() {

			glColor3f(1.f, 1.f, 0.f);
			glLineWidth(2.0f);

			int sections = 20;

			for (uint64_t i = 0; i < N; i++) {

				double tk = (double)i/N;
				double t = tk;
				double coord_x, coord_y;

				glBegin(GL_LINE_STRIP);
				for (int j = 0; j <= sections; j++) {
					double increment = h(i)/ sections;
					cubicFunction(t, tk, i, coord_x, coord_y);
					glVertex3f(coord_x, coord_y, 0.0f);
					t += increment;
				}
				glEnd();

			}
		}

		void parametricSpline(Eigen::VectorXd x_, Eigen::VectorXd y_) {

			if (x_.size() < 2) return;
			if (x_.size() != y_.size()) {
				std::cout << "O tamanho dos vetores X e Y devem ser iguais" << std::endl;
				return;
			}

			setPointCoord(x_, y_);
			setDifferenceH(1/N);
			calculeCoefficients();
			drawSplines();

		}
	}
	*/

	namespace parametric {
		//A(0) = 0; B(0) = 0
		uint64_t N;							//Número de equações
		Eigen::VectorXd x, y;				// Coordenadas x e y dos pontos de entrada
		Eigen::MatrixXd A;					// AX=B
		Eigen::VectorXd X_x, X_y, B_x, B_y;				// AX=B
		Eigen::VectorXd h;					// h é o vetor que contém "t(i) - t(i-1)"

		Eigen::MatrixXd xCoefficients, yCoefficients;		// Coeficientes dos polinômios


		void setPointCoord(Eigen::VectorXd x_, Eigen::VectorXd y_) {
			N = x_.size() - 1;
			x = x_;
			y = y_;
		}

		void setDifferenceH(double interval_t) {
			//Considerando igualmente espaçado, tal que:
			//t(i) = i/N onde (i = 0, 1, 2, 3, ..., N)
			// Deste modo 0 < t < 1 sempre será verdadeiro, e "t(i) - t(i-1)" = 1

			h = Eigen::VectorXd::Constant(N, 1, interval_t);
		}

		void setMatrixA() {
			A = Eigen::MatrixXd::Zero(N - 1, N - 1);

			for (int i = 0; i < N - 1; i++) {
				if (i > 0) A(i, i - 1) = h(i);
				A(i, i) = 2 * (h(i + 1) + h(i));
				if (i < N - 2) A(i, i + 1) = h(i + 1);
			}
		}

		void setVectorB() {
			B_x = Eigen::VectorXd::Zero(N - 1);
			B_y = Eigen::VectorXd::Zero(N - 1);

			for (uint64_t i = 1; i < N; i++) {
				B_x(i - 1) = (x(i + 1) - x(i)) / h(i) - (x(i) - x(i - 1)) / h(i - 1);
				B_y(i - 1) = (y(i + 1) - y(i)) / h(i) - (y(i) - y(i - 1)) / h(i - 1);
			}


			B_x *= 6;
			B_y *= 6;

		}

		void solveEquations() {
			X_x = Eigen::VectorXd::Zero(N + 1);
			X_y = Eigen::VectorXd::Zero(N + 1);

			if (A.rows() != B_x.size())
				std::cout << "A matriz de A e o vetor B não devem ter quantidade de linhas diferentes";

			X_x.segment(1, N - 1) = A.colPivHouseholderQr().solve(B_x);
			X_y.segment(1, N - 1) = A.colPivHouseholderQr().solve(B_y);
		}

		void cubicFunction(double t, double tk, int line, double& x, double& y) {
			x = pow(t - tk, 3) * xCoefficients(line, 0) + pow(t - tk, 2) * xCoefficients(line, 1) + (t - tk) * xCoefficients(line, 2) + xCoefficients(line, 3);
			y = pow(t - tk, 3) * yCoefficients(line, 0) + pow(t - tk, 2) * yCoefficients(line, 1) + (t - tk) * yCoefficients(line, 2) + yCoefficients(line, 3);
		}

		void calculeCoefficients() {

			xCoefficients = Eigen::MatrixXd::Zero(N, 4);
			yCoefficients = Eigen::MatrixXd::Zero(N, 4);



			for (int i = 1; i <= N; i++) {
				xCoefficients(i - 1, 0) = (X_x(i) - X_x(i - 1)) / (6 * h(i - 1));
				yCoefficients(i - 1, 0) = (X_y(i) - X_y(i - 1)) / (6 * h(i - 1));

				xCoefficients(i - 1, 1) = (X_x(i)) / 2;
				yCoefficients(i - 1, 1) = (X_y(i)) / 2;

				xCoefficients(i - 1, 2) = (x(i) - x(i - 1)) / h(i - 1) + h(i - 1) * (2 * X_x(i) + X_x(i - 1)) / 6;
				yCoefficients(i - 1, 2) = (y(i) - y(i - 1)) / h(i - 1) + h(i - 1) * (2 * X_y(i) + X_y(i - 1)) / 6;

				xCoefficients(i - 1, 3) = x(i);
				yCoefficients(i - 1, 3) = y(i);
			}

		}

		void drawSplines() {

			glColor3f(1.f, 0.f, 1.f);
			glLineWidth(2.0f);

			int sections = 20;

			for (uint64_t i = 0; i < N; i++) {

				double tk = (double)(i + 1) / N;
				double t = (double)i / N;
				double coord_x, coord_y;
				double increment = h(i) / sections;

				glBegin(GL_LINE_STRIP);
				for (int j = 0; j <= sections; j++) {
					cubicFunction(t, tk, i, coord_x, coord_y);
					glVertex3f(coord_x, coord_y, 0.0f);
					t += increment;
				}
				glEnd();

			}
		}

		void parametricSpline(Eigen::VectorXd x_, Eigen::VectorXd y_) {

			if (x_.size() < 3) return;
			if (x_.size() != y_.size()) {
				std::cout << "O tamanho dos vetores X e Y devem ser iguais" << std::endl;
				return;
			}

			setPointCoord(x_, y_);
			setDifferenceH((double)1 / N);
			setMatrixA();
			setVectorB();
			solveEquations();
			calculeCoefficients();
			drawSplines();
		}
	}
}
