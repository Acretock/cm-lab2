#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <cstdlib>   
#include <ctime> 
#include <map>
#include <omp.h>

const double PI = 3.141592653589793238;

using namespace std;

double function(double x) {
	return x * sin(x) * cos(2.0 * x);
}
const double exactValue = -1.206375395410046;
const double epsilon = 1e-7;
bool compare_float(double a, double b) {
	return fabs(a - b) < epsilon;
}
double RectIntegral(double a, double b, int n) {
	double result = 0, h = (b - a) / n;

	for (int i = 0; i < n; i++) {
		result += function(a + h * (i + 0.5));
	}
	result *= h;
	return result;
}
double TrapezoidalIntegral(double a, double b, int n) {
	const double width = (b - a) / n;

	double trapezoidal_integral = 0;
	for (int step = 0; step < n; step++) {
		const double x1 = a + step * width;
		const double x2 = a + (step + 1) * width;

		trapezoidal_integral += 0.5 * (x2 - x1) * (function(x1) + function(x2));

	}

	return trapezoidal_integral;
}
double TrapezoidalIntegral5(int n, std::vector<double>& xValues) {


	double trapezoidal_integral = 0;
	for (int step = 0; step < n - 1; step++) {
		double x1 = xValues[step];
		double x2 = xValues[step + 1];

		trapezoidal_integral += 0.5 * (x2 - x1) * (function(x1) + function(x2));

	}

	return trapezoidal_integral;
}

double SimpsonIntegral(double a, double b, int n) {
	double h = (b - a) / n;

	// Internal sample points, there should be n - 1 of them
	double sum_odds = 0.0;
	for (int i = 1; i < n; i += 2)
	{
		sum_odds += function(a + i * h);
	}
	double sum_evens = 0.0;
	for (int i = 2; i < n; i += 2)
	{
		sum_evens += function(a + i * h);
	}

	return (function(a) + function(b) + 2 * sum_evens + 4 * sum_odds) * h / 3;
}

double romberg(double a, double b) {
	const int N = 5;
	double h[N + 1], r[N + 1][N + 1];
	for (int i = 1; i < N + 1; ++i) {
		h[i] = (b - a) / pow(2, i - 1);
	}
	r[1][1] = h[1] / 2 * (function(a) + function(b));
	for (int i = 2; i < N + 1; ++i) {
		double coeff = 0;
		for (int k = 1; k <= pow(2, i - 2); ++k) {
			coeff += function(a + (2 * k - 1) * h[i]);
		}
		r[i][1] = 0.5 * (r[i - 1][1] + h[i - 1] * coeff);
	}

	for (int i = 2; i < N + 1; ++i) {
		for (int j = 2; j <= i; ++j) {
			r[i][j] = r[i][j - 1] + (r[i][j - 1] - r[i - 1][j - 1]) / (pow(4, j - 1) - 1);
		}
	}
	return r[N][N];
}

double chebyshevNode(double i, double n0) {
	auto a = 2.0;
	auto b = 0.0;
	return ((a + b) / 2.0 + (a - b) / 2.0 * cos((PI * (2.0 * i + 1.0) / (2.0 * n0))));
}

void task_1() {
	std::ofstream out_1("Rect.txt");
	std::ofstream out_2("Trap.txt");
	double i = 1;
	double a = RectIntegral(0, 2, i);
	double b = TrapezoidalIntegral(0, 2, i);
	while (!compare_float(a, exactValue)) {
		i++;
		a = RectIntegral(0, 2, i);
		out_1 << std::fixed << std::setprecision(12) << abs(a - exactValue) << std::endl;

	}
	int k = 1;
	std::cout << std::fixed << std::setprecision(12) << "Exact value:\t\t\t" << exactValue << std::endl;
	std::cout << std::fixed << std::setprecision(12) << "Rectangles method:\t\t" << a << std::endl << "Rectangles method error margin:\t" << abs(a - exactValue) << std::endl;
	std::cout << std::setprecision(0) << "Rectangles method iterations:\t " << i << std::endl;
	while (!compare_float(b, exactValue)) {
		k++;
		b = TrapezoidalIntegral(0, 2, k);
		out_2 << std::fixed << std::setprecision(12) << abs(b - exactValue) << std::endl;
	}

	std::cout << std::fixed << std::setprecision(12) << "Trapezoid metod:\t\t" << b << std::endl << "Trapezoid metod error magin:\t" << abs(b - exactValue) << std::endl;

	std::cout << "Trapezoid metod iterations:\t" << k << std::endl;
	out_1.close();
	out_2.close();
}

void task_2() {
	std::ofstream out_1("Simpson.txt");
	std::ofstream out_2("Delta.txt");
	double i = 1;
	double a = SimpsonIntegral(0, 2, i);
	while (!compare_float(a, exactValue)) {
		i++;
		a = SimpsonIntegral(0, 2, i);
		out_1 << std::fixed << std::setprecision(12) << a << std::endl;
		out_2 << std::fixed << std::setprecision(12) << abs(a - exactValue) << std::endl;
	}

	std::cout << std::fixed << std::setprecision(12) << "Exact value:\t" << exactValue << std::endl;
	std::cout << std::fixed << std::setprecision(12) << "Simpson's method :\t" << a << std::endl << "Error margin :\t" << abs(a - exactValue) << std::endl;
	std::cout << "Simpson's method iteration number: " << i << std::endl;

	out_1.close();
	out_2.close();
}

void task_3() {

	std::cout << "                 Aitken process               " << std::endl;
	double q = 100;
	double n = 10;
	double F1 = TrapezoidalIntegral(0, 2, n);
	double F2 = TrapezoidalIntegral(0, 2, q * n);
	double F3 = TrapezoidalIntegral(0, 2, q * q * n);

	 //Calculation of the order of the main term of the error 
	
	double p = 1.0 / log(q) * log((F3 - F2) / (F2 - F1));
	double I = F1 + (F1 - F2) * (F1 - F2) / (2 * F2 - F1 - F3);
	std::cout << std::fixed << std::setprecision(12) << "Exact value: \t\t\t" << exactValue << std::endl;
	std::cout << std::fixed << std::setprecision(12) << "Clarified value : \t\t" << I << std::endl;
	std::cout << std::fixed << std::setprecision(12) << "Effective order of precision : \t" << p << std::endl;
	std::cout << std::fixed << std::setprecision(12) << "Error margin: \t\t\t" << abs(I - exactValue) << std::endl << std::endl << std::endl;
	std::cout << "                  Runge  Method" << std::endl;

	double n1 = 1000;
	double n2 = 2 * n1;
	double n3 = 4 * n2;
	double Sm1 = TrapezoidalIntegral(0, 2, n1);
	//std::cout << std::fixed << std::setprecision(12) << Sm1 << "\t" << abs(exactValue - Sm1);


	double Sm2 = TrapezoidalIntegral(0, 2, n2);
	double Sm3 = TrapezoidalIntegral(0, 2, n3);
	I = Sm2 + 1 / 3 * (Sm2 - Sm1);


	p = 1.0 / log(2) * log((Sm3 - Sm2) / (Sm2 - Sm1));

std::cout << std::fixed << std::setprecision(12) << "Exact value: \t\t\t" << exactValue << std::endl;
std::cout << std::fixed << std::setprecision(12) << "Clarified value: \t\t" << I << std::endl;
std::cout << std::fixed << std::setprecision(12) << "Effective order of precision \t" << p << std::endl;
std::cout << std::fixed << std::setprecision(12) << "Error margin \t\t\t" << abs(I - exactValue) << std::endl << std::endl << std::endl;


	std::cout << "                 Romberg method             " << std::endl;
	double a = romberg(0, 2);
	std::cout << std::fixed << std::setprecision(12) << "Exact value: \t\t" << exactValue << std::endl;
	std::cout << std::fixed << std::setprecision(12) << "Romberg method value: \t" << a << std::endl;
	std::cout << std::fixed << std::setprecision(12) << "Error margin: \t\t" << abs(a - exactValue) << std::endl << std::endl << std::endl;


	std::cin.clear();
}

void task_5() {
	double answer = 0;
	double delta = 0;
	double iterations = 0;
	for (int N = 10; N < 1000; N++) {
		std::vector<double> xValues;
		std::vector<double> yValues;
		double temp;

		for (int i = 0; i < N; i++)
		{
			xValues.push_back(chebyshevNode(i, N));
			yValues.push_back(function(xValues[i]));
		}
		// chebyshevNode has reversed order
		std::reverse(xValues.begin(), xValues.end());
		std::reverse(yValues.begin(), yValues.end());

		temp = TrapezoidalIntegral5(N, xValues);
		delta = abs(temp - answer);
		if (delta < epsilon) {
			answer = temp;
			iterations = N;
			break;
		}
		answer = temp;
	}
	std::cout << std::fixed << std::setprecision(12) << "Exact value: \t\t" << exactValue << std::endl;
	std::cout << std::fixed << std::setprecision(12) << "Resulted value: \t" << answer << std::endl;
	std::cout << std::fixed << std::setprecision(12) << "Error margin: \t\t" << abs(exactValue - answer) << std::endl;
	std::cout << "Nedded nodes : \t\t" << iterations << std::endl;

	std::cin.clear();
}

void task_6() {

	std::cout << "             Monte-Carlo Metod               " << std::endl;
	std::ofstream out("task_6_out.txt");
	std::ofstream out2("task_6.1_out.txt");

	srand(time(NULL));
	double Integral = 0;
	double x1 = 0.0, x2 = 2.0;

	for (int N = 10; N < 2000; N++) {

		Integral = 0;
		for (int k = 1; k <= 100; k++) {

			double Sum = 0;
			for (int i = 0; i < N; i++) {
				Sum += function(static_cast<double>(rand()) / function(static_cast<double>(RAND_MAX) * (x2 - x1) + x1) * (fabs(x2 - x1)));
			}

			double meanSum = Sum / N;
			Integral += meanSum;
		}
		Integral /= 100.0;

		double error = (Integral - exactValue);
		out << std::fixed << std::setprecision(12) << Integral << std::endl;
		out2 << std::fixed << std::setprecision(12) << error << std::endl;
	}


	std::cout << std::fixed << std::setprecision(12) << "Exact value: \t\t\t" << exactValue << std::endl;
	std::cout << std::fixed << std::setprecision(12) << "Resulted value: \t\t" << Integral << std::endl;
	std::cout << std::fixed << std::setprecision(12) << "Error margin: \t\t\t" << abs(exactValue - Integral) << std::endl;
	std::cin.clear();
	out.close();
	out2.close();
}

int main() {
	int c = -1;
	while (c != 0) {
		cout << "enter number for task, 0 for exit" << endl;
		cin >> c;
		switch (c) {
		case 1:
			task_1();
			break;
		case 2:
			task_2();
			break;
		case 3:
			task_3();
			break;
		case 5:
			task_5();
			break;
		case 6:
			task_6();
			break;
		default:
			break;
		}
	}
}