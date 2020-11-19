#pragma once

namespace spline {


	/// <summary>
		/// Ponteiro para uma função.
		/// </summary>
	typedef double (*function)(double);


	/// <summary>
		/// Desenha os pontos na tela criada pelo OpenGL.
		/// </summary>
		/// <param name="points"> vetor com todos os pontos que devem ser desenhados</param>
	void drawPoints(Eigen::MatrixXd points);

	/// <summary>
		/// Desenha os pontos na tela criada pelo OpenGL
		///	Note: x e y devem ter o mesmo tamanho!
		/// </summary>
		/// <param name="x">coordenada x de todos os pontos</param>
		/// <param name="y">coordenada y de todos os pontos</param>
	void drawPoints(Eigen::VectorXd x, Eigen::VectorXd y);

	/// <summary>
		/// Desenha os eixos na tela criada pelo OpenGL.
		///	Note: Os limites da tela são definidos considerandos coordenadas cartesianas.
		/// </summary>
		/// <param name="xmin">o limite esquerdo da tela</param>
		/// <param name="xmax">o limite direito da tela</param>
		/// <param name="ymin">o limite inferior da tela</param>
		/// <param name="ymax">o limite superior da tela</param>
	void drawAxis(double xmin, double xmax, double ymin, double ymax);

	/// <summary>
		/// Desenha uma função qualquer na tela criada pelo OpenGL.
		/// </summary>
		/// <param name="f">ponteiro para a função a ser desenhada</param>
		/// <param name="xa">o limite onde a função começa a ser desenhada</param>
		/// <param name="xb">o limite onde a função termina de ser  desenhada</param>
	void drawFunction(function f, double xa, double xb);


	namespace quadratic {

		/// <summary>
		/// Armazena as coordenadas x e y dos pontos recebidos.
		///Note: Tem de ter o mesmo número de coordenadas.
		/// </summary>
		/// <param name="x">coordenada x dos pontos</param>
		/// <param name="y">coordenada y dos pontos</param>
		void setPointCoord(Eigen::VectorXd x, Eigen::VectorXd y);

		/// <summary>
		/// Os coeficientes das parábolas são recebidos em um vetor,
		/// esta função os organiza em uma matriz nx3.
		/// </summary>
		/// <param name="coeff">vetor com os coeficientes a, b e c das parábolas  </param>
		void setCoefficients(Eigen::VectorXd coeff);

		/// <summary>
		/// Organiza a matriz A para poder ser posteriormente utilizado na na resolução do sistema linear.
		/// </summary>
		/// <param name="x">coordenada x de todos os pontos</param>
		/// <param name="out">matriz A, que será utilizada para a resolução do sistema linear</param>
		void setMatrixA(Eigen::VectorXd x, Eigen::MatrixXd& out);

		/// <summary>
		/// Organiza o vetor B para poder ser posteriormente utilizado na na resolução do sistema linear.
		/// </summary>
		/// <param name="x">coordenada y de todos os pontos</param>
		/// <param name="out">vetor B, que será utilizada para a resolução do sistema linear</param>
		void setVectorB(Eigen::VectorXd y, Eigen::VectorXd& out);


		/// <summary>
		/// Resolução do sistema linear Ax = B.
		/// </summary>
		/// <param name="A">matriz A</param>
		/// <param name="X">Vetor que vai receber a resolução do sistema</param>
		/// <param name="B">vetor B</param>
		void solveEquations(Eigen::MatrixXd A, Eigen::VectorXd& X, Eigen::VectorXd B);

		/// <summary>
		/// Calcula a resolução das equações para a implementação de splines quadráticas,
		/// tem como resultado o desenho das equações em uma tela do OpenGL.
		/// </summary>
		/// <param name="x">coordenada x dos pontos</param>
		/// <param name="y">coordenada y dos pontos</param>
		void quadraticSpline(Eigen::VectorXd x, Eigen::VectorXd y);

		/// <summary>
		/// Equação matemática para calcular o resultado de uma entrada
		/// em um polinômio do segundo grau.
		/// </summary>
		/// <param name="x">incógnita x</param>
		/// <param name="a">coeficiente a</param>
		/// <param name="b">coeficiente b</param>
		/// <param name="c">coeficiente c</param>
		/// <returns>valor de y ou f(x)</returns>
		double quadFunction(double x, double a, double b, double c);

		/// <summary>
		/// Desenha as parábolas que formam a spline quadrática em uma janela do OpenGL.
		/// </summary>
		/// <param name="x">coordenada x dos pontos</param>
		/// <param name="coeff">coeficientes de todas as equações</param>
		void drawSpline(Eigen::VectorXd x, Eigen::MatrixXd coeff);
	}


	namespace cubic {

		/// <summary>
		/// Armazena as coordenadas x e y dos pontos recebidos.
		///Note: Tem de ter o mesmo número de coordenadas.
		/// </summary>
		/// <param name="x">coordenada x dos pontos</param>
		/// <param name="y">coordenada y dos pontos</param>
		void setPointCoord(Eigen::VectorXd x, Eigen::VectorXd y);

		/// <summary>
		/// Calcula a diferença entre a coordenada x de um ponto e a coordenada anterior.
		/// As diferenças calculadas são armazenadas em um vetor.
		/// "x(i) - x(i-1)"
		/// </summary>
		/// <param name="x">coordenada x dos pontos</param>
		/// <param name="out"> Vetor de saída, que irá receber o resultado das diferenças</param>
		void setDifferenceH(Eigen::VectorXd& x, Eigen::VectorXd& out);

		/// <summary>
		/// Organiza a matriz A para poder ser posteriormente utilizado na na resolução do sistema linear.
		/// </summary>
		/// <param name="h">vetor com as diferenças entre as coordenadas x</param>
		/// <param name="out">matriz A, que será utilizada para a resolução do sistema linear</param>
		void setMatrixA(Eigen::VectorXd h, Eigen::MatrixXd& out);

		/// <summary>
		/// Organiza o vetor B para poder ser posteriormente utilizado na na resolução do sistema linear.
		/// </summary>
		/// <param name="h">vetor com as diferenças entre as coordenadas x</param>
		/// <param name="y">coordenada y de todos os pontos</param>
		/// <param name="out">vetor B, que será utilizada para a resolução do sistema linear</param>
		void setVectorB(Eigen::VectorXd h, Eigen::VectorXd y, Eigen::VectorXd& out);

		/// <summary>
		/// Os coeficientes dos polinômios de grau 3 são recebidos em um vetor.
		/// Esta função os organiza em uma matriz nx4.
		/// </summary>
		/// <param name="h">vetor com as diferenças entre as coordenadas x</param>
		/// <param name="y">coordenada y de todos os pontos</param>
		/// <param name="X"></param>
		/// <param name="coeff">vetor com os coeficientes a, b, c e d das dos polinômios  </param>
		void setCoefficients(Eigen::VectorXd h, Eigen::VectorXd y, Eigen::VectorXd X, Eigen::MatrixXd& coeff);

		/// <summary>
		/// Calcula a resolução das equações para a implementação de splines cúbicas,
		/// tem como resultado o desenho das equações em uma tela do OpenGL.
		/// </summary>
		/// <param name="x">coordenada x dos pontos</param>
		/// <param name="y">coordenada y dos pontos</param>
		void cubicSpline(Eigen::VectorXd x, Eigen::VectorXd y);

		/// <summary>
		/// Equação matemática para calcular o resultado de uma entrada
		/// em um polinômio do terceiro grau.
		/// </summary>
		/// <param name="x">incógnita x</param>
		/// <param name="xk">uma constante</param>
		/// <param name="a">coeficiente a</param>
		/// <param name="b">coeficiente b</param>
		/// <param name="c">coeficiente c</param>
		/// <param name="d">coeficiente d</param>
		/// <returns>valor de y ou f(x)</returns>
		double cubicFunction(double x, double xk, double a, double b, double c, double d);

		/// <summary>
		/// Desenha os polinômios que formam a spline cóbica em uma janela do OpenGL.
		/// </summary>
		/// <param name="x">coordenada x dos pontos</param>
		/// <param name="coeff">coeficientes de todas as equações</param>
		void drawSpline(Eigen::VectorXd x, Eigen::MatrixXd coeff);
	}

	/*
	namespace parametric1 {
		void setPointCoord(Eigen::VectorXd x, Eigen::VectorXd y);
		void setDifferenceH(double interval_t = 1);
		void cubicFunction(double t, double tk, int line, double& x, double& y);
		void calculeCoefficients();
		void drawSplines();
		void parametricSpline(Eigen::VectorXd x_, Eigen::VectorXd y_);
	}
	*/

	namespace parametric {
		void setPointCoord(Eigen::VectorXd x, Eigen::VectorXd y);
		void setDifferenceH(double interval_t);
		void setMatrixA();
		void setVectorB();
		void cubicFunction(double t, double tk, int line, double& x, double& y);
		void solveEquations();
		void calculeCoefficients();
		void drawSplines();
		void parametricSpline(Eigen::VectorXd x_, Eigen::VectorXd y_);
	}
}
