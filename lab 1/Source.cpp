#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>

long double a = -1.4;
long double b = 0.2;
long double Eps = 0.00000001;
long double ans = -0.111832559158963; //Посчитано в matlab
long double ans2 = -0.66666666666666667; //Посчитано в desmos
long double sup_dF = 0.40552;
long double sup_dF2 = 0.728;
long double startEps = 1;
int cnt = 12;
int sampling = 40;
long double x0 = (a + b) / 2;




long double F(long double x) {
	return  x * std::expl(x) + 0.1;
}
long double F_(long double x) {
	return  - 0.1 / std::expl(x);
}

long double F2(long double x) {
	return  3 * std::pow(x, 3) + 5 * std::pow(x, 2) + 5 * x + 2;
}
long double F_2(long double x) {
	return (- 3 * std::pow(x, 3) - 5 * std::pow(x, 2) - 2)/ 5;
}


long double Error(long double point, long double ans) {
	return  ans - point;
}

void sPrint(std::ofstream& stream, long double num) {
	if (stream.is_open()) { stream << num << std::endl; }
}

long double X_kPlusOne(long double x0, long double x1, long double x2) {
	return ((x1*x1 - x2*x0)/(2.*x1 - x2 - x0));
}

void BisectionMethod(std::ofstream& result, std::ofstream& error,long double input_x0,long double a,long double b, long double Eps, long double ans, long double func(long double)) {
	long double c = input_x0;
	sPrint(result, c);
	sPrint(error, Error(c, ans));
	while (abs(b - a) >= (2. * Eps)) {
		if (func(a) * func(c) < 0) {
			b = c;
		}
		else {
			a = c;
		}
		c = (a + b) / 2;
		sPrint(result, c);
		sPrint(error, Error(c, ans));
	}
}
int BisectionMethod(long double input_x0,long double a, long double b, long double Eps, long double func(long double)) {
	long double c = input_x0;
	int cnt = 1;
	while (abs(b - a) > (2. * Eps)) {
		if (func(a) * func(c) < 0) {
			b = c;
		}
		else {
			a = c;
		}
		c = (a + b) / 2;
		++cnt;
	}
	return cnt;
}

std::pair<int, long double> BisectionMethod(long double input_x0, long double a, long double b, long double Eps, long double ans, long double func(long double)) {
	long double c = input_x0;
	int cnt = 1;
	while (abs(b - a) > (2. * Eps)) {
		if (func(a) * func(c) < 0) {
			b = c;
		}
		else {
			a = c;
		}
		c = (a + b) / 2;
		++cnt;
	}
	return { cnt, Error(c, ans) };
}

void SimpleIterationMethod(std::ofstream& result, std::ofstream& error, long double input_x0, long double Eps, long double ans, long double sup_dF ,long double func(long double), long double func_(long double)) {
	long double x0 = input_x0;
	sPrint(result, x0);
	sPrint(error, Error(x0, ans));
	long double x1 = func_(x0);
	sPrint(result, x1);
	sPrint(error, Error(x1, ans));
	long double x2 = func_(x1);
	sPrint(result, x2);
	sPrint(error, Error(x2, ans));
	while (1) {
		long double x2_ = X_kPlusOne(x0, x1, x2);
		long double x3 = func_(x2_);
		sPrint(result, x3);
		sPrint(error, Error(x3, ans));
		if (abs(x3 - x2_) > (sup_dF / (1. - sup_dF)) * Eps) {
			x0 = x2_;
			x1 = x3;
			x2 = func_(x1);
		}
		else {
			break;
		}
	}
}

int SimpleIterationMethod(long double input_x0, long double Eps, long double sup_dF, long double func(long double), long double func_(long double)) {
	int cnt = 0;
	long double x0 = input_x0;
	++cnt;
	long double x1 = func_(x0);
	++cnt;
	long double x2 = func_(x1);
	++cnt;
	while (1) {
		long double x2_ = X_kPlusOne(x0, x1, x2);
		long double x3 = func_(x2_);
		++cnt;
		if (abs(x3 - x2_) > (sup_dF / (1. - sup_dF)) * Eps) {
			x0 = x2_;
			x1 = x3;
			x2 = func_(x1);
		}
		else { 
			return cnt;
		}
	}
}

std::pair<int, long double> SimpleIterationMethod(long double input_x0, long double Eps, long double ans, long double sup_dF, long double func(long double), long double func_(long double)) {
	int cnt = 0;
	long double x0 = input_x0;
	++cnt;
	long double x1 = func_(x0);
	++cnt;
	long double x2 = func_(x1);
	++cnt;
	while (1) {
		long double x2_ = X_kPlusOne(x0, x1, x2);
		long double x3 = func_(x2_);
		++cnt;
		if (abs(x3 - x2_) > (sup_dF / (1. - sup_dF)) * Eps) {
			x0 = x2_;
			x1 = x3;
			x2 = func_(x1);
		}
		else {
			return { cnt, Error(x3, ans) };
		}
	}
}


void ComputationalCostBisectionMethod(std::ofstream& cost, std::ofstream& achievability, long double input_x0, long double a, long double b, long double startEps, long double ans ,int cnt, long double func(long double)) {
	std::pair<int, long double>  currentPair;
	currentPair = BisectionMethod(input_x0, a, b, startEps, ans, func);
	sPrint(cost, currentPair.first);
	sPrint(achievability, currentPair.second);
	for (int i = 1; i != cnt; ++i) {
		currentPair  = BisectionMethod(input_x0, a, b, startEps/(std::pow(10,i)), ans, func);
		sPrint(cost, currentPair.first);
		sPrint(achievability, currentPair.second);
	}
}

void ComputationalCostSimpleIterationMethod(std::ofstream& cost, std::ofstream& achievability , long double x0, long double Eps, long double ans, long double sup_dF, long double func(long double), long double func_(long double)) {
	std::pair<int, long double> currentPair;
	currentPair = SimpleIterationMethod(x0, startEps, ans, sup_dF, func, func_);
	sPrint(cost, currentPair.first);
	sPrint(achievability, currentPair.second);
	for (int i = 1; i != cnt; ++i) {
		currentPair = SimpleIterationMethod(x0, startEps / (std::pow(10, i)), ans, sup_dF, func, func_);
		sPrint(cost, currentPair.first);
		sPrint(achievability, currentPair.second);
	}
}


void NumberIterationsOfChoiceApproximationBisectionMethod(std::ofstream& IterNumApproximate, int sampling, long double ans,long double Eps, long double func(long double)) {
	for (int i = 0; i != (sampling + 1); ++i) {
		long double x0 = a + ((b - a) / sampling) * i;
		if (x0 > ans) {
			break;
		}
		sPrint(IterNumApproximate, BisectionMethod((x0+b)/2, x0, b, Eps, func));
	}
}


void NumberIterationsOfChoiceApproximationSimpleIter(std::ofstream& IterNumApproximate, int sampling, long double Eps, long double sup_dF, long double func(long double), long double func_(long double)) {
	for (int i = 0; i != (sampling + 1); ++i) {
		long double x0 = a + ((b - a) / sampling) * i;
		sPrint(IterNumApproximate, SimpleIterationMethod(x0, Eps, sup_dF, func, func_));
	}
}


int main() {
	std::ofstream result;
	std::ofstream error;
	std::ofstream cost;
	std::ofstream achievability;
	std::ofstream choiceapr;
	std::ofstream result2;
	std::ofstream error2;
	std::ofstream cost2;
	std::ofstream achievability2;
	std::ofstream choiceapr2;
	std::ofstream result3;
	std::ofstream error3;
	std::ofstream cost3;
	std::ofstream achievability3;
	std::ofstream choiceapr3;
	std::ofstream result4;
	std::ofstream error4;
	std::ofstream cost4;
	std::ofstream achievability4;
	std::ofstream choiceapr4;

	result.open("Results.txt");
	error.open("ErrorRate.txt");
	cost.open("Cost.txt");
	achievability.open("achievability.txt");
	choiceapr.open("choiceapr.txt");
	result2.open("Results2.txt");
	error2.open("ErrorRate2.txt");
	cost2.open("Cost2.txt");
	achievability2.open("achievability2.txt");
	choiceapr2.open("choiceapr2.txt");
	result3.open("Results3.txt");
	error3.open("ErrorRate3.txt");
	cost3.open("Cost3.txt");
	achievability3.open("achievability3.txt");
	choiceapr3.open("choiceapr3.txt");
	result4.open("Results4.txt");
	error4.open("ErrorRate4.txt");
	cost4.open("Cost4.txt");
	achievability4.open("achievability4.txt");
	choiceapr4.open("choiceapr4.txt");

	BisectionMethod(result, error,x0, a, b, Eps, ans, F);
	ComputationalCostBisectionMethod(cost, achievability,x0, a, b, startEps, ans, cnt, F);
	NumberIterationsOfChoiceApproximationBisectionMethod(choiceapr, sampling, ans, Eps, F);
	SimpleIterationMethod(result2, error2, x0, Eps, ans, sup_dF, F, F_);
	ComputationalCostSimpleIterationMethod(cost2, achievability2, x0, startEps, ans, sup_dF, F, F_);
	NumberIterationsOfChoiceApproximationSimpleIter(choiceapr2, sampling, Eps, sup_dF, F, F_);

	BisectionMethod(result3, error3, x0, a, b, Eps, ans2, F2);
	ComputationalCostBisectionMethod(cost3, achievability3, x0, a, b, startEps, ans2, cnt, F2);
	NumberIterationsOfChoiceApproximationBisectionMethod(choiceapr3, sampling, ans2, Eps, F2);
	SimpleIterationMethod(result4, error4, x0, Eps, ans2, sup_dF2, F2, F_2);
	ComputationalCostSimpleIterationMethod(cost4, achievability4, x0, startEps, ans2, sup_dF2, F2, F_2);
	NumberIterationsOfChoiceApproximationSimpleIter(choiceapr4, sampling, Eps, sup_dF, F2, F_2);


	result.close();
	error.close();
	cost.close();
	achievability.close();
	choiceapr.close();
	result2.close();
	error2.close();
	cost2.close();
	achievability2.close();
	choiceapr2.close();
	result3.close();
	error3.close();
	cost3.close();
	achievability3.close();
	choiceapr3.close();
	result4.close();
	error4.close();
	cost4.close();
	achievability4.close();
	choiceapr4.close();
}





/*
#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>

double a = -1.4;
double b = 0.2;
double Eps = 0.001;
long double ans = -0.111832559158963; //Посчитано в matlab
long double ans2 = -0.66666666666666667; //Посчитано в desmos
long double sup_dF = 0.45;
long double sup_dF2 = 0.6;
long double startEps = 1;
int cnt = 16;
long double x0 = (a + b) / 2;

long double F(long double x) {
	return  x * std::expl(x) + 0.1;
}
long double F_(long double x) {
	return  - 0.1 / std::expl(x);
}

long double F2(long double x) {
	return  3 * std::pow(x, 3) + 5 * std::pow(x, 2) + 5 * x + 2;
}
long double F_2(long double x) {
	return (- 3 * std::pow(x, 3) - 5 * std::pow(x, 2) - 2)/ 5;
}


long double Error(long double point, long double ans) {
	return  ans - point;
}

void sPrint(std::ofstream& stream, long double num) {
	if (stream.is_open()) { stream << num << std::endl; }
}

void BisectionMethod(std::ofstream& result, std::ofstream& error,long double a,long double b, long double Eps, long double ans, long double func(long double)) {
	long double c = (a + b) / 2;
	sPrint(result, c);
	sPrint(error, Error(c, ans));
	while (abs(b - a) >= (2. * Eps)) {
		if (func(a) * func(c) < 0) {
			b = c;
		}
		else {
			a = c;
		}
		c = (a + b) / 2;
		sPrint(result, c);
		sPrint(error, Error(c, ans));
	}
}
int BisectionMethod(long double a, long double b, long double Eps, long double func(long double)) {
	long double c = (a + b) / 2;
	int cnt = 1;
	while (abs(b - a) > (2. * Eps)) {
		if (func(a) * func(c) < 0) {
			b = c;
		}
		else {
			a = c;
		}
		c = (a + b) / 2;
		++cnt;
	}
	return cnt;
}
std::pair<int, long double> BisectionMethod(long double a, long double b, long double Eps, long double ans, long double func(long double)) {
	long double c = (a + b) / 2;
	int cnt = 1;
	while (abs(b - a) > (2. * Eps)) {
		if (func(a) * func(c) < 0) {
			b = c;
		}
		else {
			a = c;
		}
		c = (a + b) / 2;
		++cnt;
	}
	return { cnt, Error(c, ans) };
}

void SimpleIterationMethod(std::ofstream& result, std::ofstream& error, long double x0, long double Eps, long double ans, long double sup_dF ,long double func(long double), long double func_(long double)) {
	long double currentX = x0;
	sPrint(result, currentX);
	sPrint(error, Error(currentX, ans));
	long double newX = func_(currentX);
	sPrint(result, newX);
	sPrint(error, Error(newX, ans));
	while((sup_dF/(1. - sup_dF))*abs(newX-currentX) > Eps){
		currentX = newX;
		newX = func_(currentX);
		sPrint(result, newX);
		sPrint(error, Error(newX, ans));
	}
}
int SimpleIterationMethod(long double x0, long double Eps, long double sup_dF, long double func(long double), long double func_(long double)) {
	int cnt = 0;
	long double currentX = x0;
	++cnt;
	long double newX = func_(currentX);
	++cnt;
	while ((sup_dF / (1. - sup_dF)) * abs(newX - currentX) > Eps) {
		long double delta = func(newX + Eps) * func(newX - Eps);
		currentX = newX;
		newX = func_(currentX);
		++cnt;
	}
	return cnt;
}
std::pair<int, long double> SimpleIterationMethod(long double x0, long double Eps, long double ans, long double sup_dF, long double func(long double), long double func_(long double)) {
	int cnt = 0;
	long double currentX = x0;
	++cnt;
	long double newX = func_(currentX);
	++cnt;
	while ((sup_dF / (1. - sup_dF)) * abs(newX - currentX) > Eps) {
		long double delta = func(newX + Eps) * func(newX - Eps);
		currentX = newX;
		newX = func_(currentX);
		++cnt;
	}
	return { cnt, Error(newX, ans) };
}


void ComputationalCostBisectionMethod(std::ofstream& cost, std::ofstream& achievability, long double a, long double b, long double startEps, long double ans ,int cnt, long double func(long double)) {
	std::pair<int, long double>  currentPair;
	currentPair = BisectionMethod(a, b, startEps, ans, func);
	sPrint(cost, currentPair.first);
	sPrint(achievability, currentPair.second);
	for (int i = 1; i != cnt; ++i) {
		currentPair  = BisectionMethod(a, b, startEps/(std::pow(10,i)), ans, func);
		sPrint(cost, currentPair.first);
		sPrint(achievability, currentPair.second);
	}
}

void ComputationalCostSimpleIterationMethod(std::ofstream& cost, std::ofstream& achievability , long double x0, long double Eps, long double ans, long double sup_dF, long double func(long double), long double func_(long double)) {
	std::pair<int, long double> currentPair;
	currentPair = SimpleIterationMethod(x0, startEps, ans, sup_dF, func, func_);
	sPrint(cost, currentPair.first);
	sPrint(achievability, currentPair.second);
	for (int i = 1; i != cnt; ++i) {
		currentPair = SimpleIterationMethod(x0, startEps / (std::pow(10, i)), ans, sup_dF, func, func_);
		sPrint(cost, currentPair.first);
		sPrint(achievability, currentPair.second);
	}
}


int main() {
	std::ofstream result;
	std::ofstream error;
	std::ofstream cost;
	std::ofstream achievability;
	std::ofstream result2;
	std::ofstream error2;
	std::ofstream cost2;
	std::ofstream achievability2;
	std::ofstream result3;
	std::ofstream error3;
	std::ofstream cost3;
	std::ofstream achievability3;
	std::ofstream result4;
	std::ofstream error4;
	std::ofstream cost4;
	std::ofstream achievability4;

	result.open("Results.txt");
	error.open("ErrorRate.txt");
	cost.open("Cost.txt");
	achievability.open("achievability.txt");
	result2.open("Results2.txt");
	error2.open("ErrorRate2.txt");
	cost2.open("Cost2.txt");
	achievability2.open("achievability2.txt");
	result3.open("Results3.txt");
	error3.open("ErrorRate3.txt");
	cost3.open("Cost3.txt");
	achievability3.open("achievability3.txt");
	result4.open("Results4.txt");
	error4.open("ErrorRate4.txt");
	cost4.open("Cost4.txt");
	achievability4.open("achievability4.txt");

	BisectionMethod(result, error, a, b, Eps, ans, F);
	ComputationalCostBisectionMethod(cost, achievability, a, b, startEps, ans, cnt, F);
	SimpleIterationMethod(result2, error2, x0, Eps, ans, sup_dF, F, F_);
	ComputationalCostSimpleIterationMethod(cost2, achievability2, x0, startEps, ans, sup_dF, F, F_);

	BisectionMethod(result3, error3, a, b, Eps, ans2, F2);
	ComputationalCostBisectionMethod(cost3, achievability3, a, b, startEps, ans2, cnt, F2);
	SimpleIterationMethod(result4, error4, x0, Eps, ans2, sup_dF2, F2, F_2);
	ComputationalCostSimpleIterationMethod(cost4, achievability4, x0, startEps, ans2, sup_dF2, F2, F_2);


	result.close();
	error.close();
	cost.close();
	achievability.close();
	result2.close();
	error2.close();
	cost2.close();
	achievability2.close();
	result3.close();
	error3.close();
	cost3.close();
	achievability3.close();
	result4.close();
	error4.close();
	cost4.close();
	achievability4.close();
}

*/