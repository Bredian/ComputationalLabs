//Soft works for 2nd order Runge–Kutta methods
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#define FTYPE long double

//Settings
//y' = DYX function
#define FUNC (2 * y / (x * ((2 * pow(x, 2) * y * log(y)) - 1)))
//y(X0) = Y0
#define X0_VALUE 1.0
#define Y0_VALUE 1
// BORDER_LEFT < X < BORDER_RIGHT
#define BORDER_LEFT 1.0
#define BORDER_RIGHT 1.2
#define EPSILON 1E-4


FTYPE func(FTYPE x, FTYPE y) {
	return FUNC;
}
//Get f_1 and f_2
FTYPE* evaluateFunc(int funcNum, FTYPE xn, FTYPE yn, FTYPE MethodMatrix[3][3], FTYPE h) {
	FTYPE* f = (FTYPE*)malloc(sizeof(FTYPE) * 2);
	for (int i = 0; i < 2; i++) {
		f[i] = 0;
	}
	f[0] = func(xn, yn);
	for (int i = 1; i < funcNum - 1; i++) {
		FTYPE x = xn + MethodMatrix[i][0] * h;
		FTYPE y = yn;
		for (int j = 1; j < i + 2; j++) {
			y += MethodMatrix[i][j] * f[j - 1] * h;
		}
		f[i] = func(x, y);
	}
	return f;
}
//Get y_n+1 from y_n
FTYPE evaluateY(FTYPE yn, FTYPE methodMatrix[3][3], FTYPE* functions, int length, FTYPE h) {
	FTYPE ynp1 = yn;
	for (int i = 1; i < length; i++) {
		FTYPE a = methodMatrix[length - 1][i];
		ynp1 += functions[i - 1] * methodMatrix[length - 1][i] * h;
	}
	return ynp1;
}
//Reassign X0, Y0, borders and epsilon
FTYPE X0 = X0_VALUE;
FTYPE Y0 = Y0_VALUE;

FTYPE borders[2] = { BORDER_LEFT, BORDER_RIGHT };

FTYPE epsilon = EPSILON;

//Get h value
FTYPE evalGridStepH(int steps) {
	return (borders[1] - borders[0]) / steps;
}
//check if our error is less than epsilon
bool checkIfErrorAcceptable(FTYPE delta, int method) {
	FTYPE val = (fabsl(delta)) / (powl(2.0, method) - 1);
	return val <= epsilon ? false : true;
}


int main(int argc, const char* argv[]) {
	//function on h grid
	FTYPE* f = NULL;
	//function on 2h grid
	FTYPE* fProxy = NULL;
	//Butcher tableau
	FTYPE methodMatrix[3][3] = {
		{0, 0, 0},
		{1, 1, 0},
		{0.5, 0.5, 0}
	};
	//points in uniform grid
	FTYPE calcPoints[11];
	//y in uniform grid with h step
	FTYPE valuesInPoints[11];
	//y in uniform grid with 2h step
	FTYPE valuesInPointsProxy[11];
	//y(2h) - y(h)
	FTYPE deltas[11];
	//initial setup
	for (int i = 0; i < 11; i++) {
		calcPoints[i] = borders[0] + ((borders[1] - borders[0]) / 10) * i;
	}
	valuesInPoints[0] = Y0;
	valuesInPointsProxy[0] = Y0;
	int stepM = 3;
	//counting appropriate number of steps
	int steps = 5;
	FTYPE maxDelta;
	FTYPE h = 0.0;
	FTYPE x = X0;
	FTYPE xProxy = X0;
	FTYPE value = Y0;
	FTYPE valueProxy = Y0;
	//finding approproate step
	do {
		x = X0;
		xProxy = X0;
		value = Y0;
		valueProxy = Y0;
		for (int i = 0; i <= steps; i++) {
			h = evalGridStepH(steps);
			x = h * i + X0;
			xProxy = 2 * h * i + X0;
			for (int j = 0; j < 11; j++) {
				if (x == calcPoints[j]) {
					valuesInPoints[j] = value;
				}
				if (xProxy == calcPoints[j]) {
					valuesInPointsProxy[j] = valueProxy;
					break;
				}
			}
			//count values for output, will be changed if step is too small
			f = evaluateFunc(stepM, x, value, methodMatrix, h);
			fProxy = evaluateFunc(stepM, xProxy, valueProxy, methodMatrix, 2 * h);
			value = evaluateY(value, methodMatrix, f, stepM, h);
			valueProxy = evaluateY(valueProxy, methodMatrix, fProxy, stepM, 2 * h);
		}
		free(f);
		free(fProxy);
		for (int i = 0; i < 11; i++) {
			deltas[i] = fabsl(valuesInPointsProxy[i] - valuesInPoints[i]);
		}
		maxDelta = deltas[0];
		for (int i = 1; i < 11; i++) {
			if (maxDelta < deltas[i]) {
				maxDelta = deltas[i];
			}
		}
		steps = steps * 2;
	} while (checkIfErrorAcceptable(maxDelta, stepM));
	FTYPE step = evalGridStepH(steps/2);
	FILE* results = fopen( "results.txt", "w+");
	FTYPE finalValue = Y0;
	fprintf(results, "x = %Lf\t y = %Lf\n", X0, finalValue);
	for (int i = 1; i <= steps/2; i++) {
		f = evaluateFunc(stepM, step * i + X0, Y0, methodMatrix, step);
		finalValue = evaluateY(finalValue, methodMatrix, f, stepM, step);
		fprintf(results, "x = %Lf\t y = %Lf\n", X0 + i * step, finalValue);
		free(f);
	}
	fclose(results);
	printf("Calculated with %d steps\nAccuracy: %Lf\nh=%Lf\nDiff=%Le\n\n", steps / 2, epsilon, h, maxDelta);
	for (int i = 0; i < 11; i++) {
		printf("%d) x=%Lf\n\ty(h)=%Lf\n\ty(2h)=%Lf\n\t--------------\n\ty(2h)-y(h)=%Lf\n\n", i + 1, calcPoints[i], valuesInPoints[i], valuesInPointsProxy[i], valuesInPointsProxy[i] - valuesInPoints[i]);
	}
	system("pause");
	return 0;
}
