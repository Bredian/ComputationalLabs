#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define ftype long double

ftype getMax(ftype* array, int size) {
	if (size < 1) {
		return -1;
	}
	else {
		ftype max = array[0];
		for (int i = 1; i < size; i++) {
			if (array[i] > max) {
				max = array[i];
			}
		}
		return max;
	}
}

ftype sign(ftype a) {
	if (a > 0) {
		return 1;
	}
	else if (a == 0) {
		return 0;
	}
	else {
		return -1;
	}
}
ftype a0, a1, a2, a3, a4, a5, a6, n;
ftype evaluate(ftype X) {
	ftype result = a0 * powl(X, 2.0 * n) + a1 * powl(X, n + 2.0) + a2 * powl(X, n + 1.0) + a3 * powl(X, n) + a4 * powl(X, 2) + a5 * X + a6;
	return result;
}

typedef struct {
	ftype left;
	ftype right;
} interval;

void printInterval(interval a) {
	printf("\t%Lf<=Z<=%Lf\n", a.left, a.right);
}

int main()
{
	printf("Computational Mathematics Lab #1\nVariant II.2\n\n");
	
	//Initial values
	ftype gamma0 = 5.0 / 3.0;
	ftype rho0 = 1.694 * powl(10.0, -4.0);
	ftype U0 = 0.0;
	ftype P0 = 1.013 * powl(10.0, 6.0);

	ftype gamma3 = 7.0 / 5.0;
	ftype C3 = 3.6537 * powl(10.0, 4.0);
	ftype U3 = 1.229 * powl(10, 4.0);
	ftype P3 = 1.6768 * powl(10.0, 6.0);

	printf("Initial values:\n\tgamma_0=%Lf\n\trho_0=%Lf\n\tU_0=%Lf\n\tP_0=%Lf\n\n\tgamma_3=%Lf\n\tC_3=%Lf\n\tU_3=%Lf\n\tP_3=%Lf\n", gamma0, rho0, U0, P0, gamma3, C3, U3, P3);
	//Data to get coefficients
	n = 2.0 * gamma3 / (gamma3 - 1.0);

	ftype alpha0 = (gamma0 - 1.0) / (gamma0 + 1.0);
	ftype rho3 = (gamma3 * P3) / (powl(C3, 2.0));
	ftype X = P3/P0;
	ftype mu = (U3-U0) * sqrtl(rho0 * (gamma0 - 1.0)/(2.0 * P0));
	ftype nu = (2.0 / (gamma3 - 1.0)) * sqrtl((gamma3 * (gamma0 - 1.0) / 2.0) * (P3 / P0) * (rho0 / rho3));

	//Get coefficients
	a0 = powl(X, 2.0);
	a1 = -1.0 * alpha0 * powl(nu, 2.0);
	a2 = 2.0 * alpha0 * nu * (mu + nu) * X;
	a3 = -1.0 * (2.0 + powl(mu + nu, 2.0)) * X;
	a4 = -1.0 * powl(nu, 2.0);
	a5 = 2.0 * nu * (mu + nu);
	a6 = -1.0 * powl((mu + nu), 2.0) + 1;

	printf("\nFirst step results(coefficients):\n\ta0=%Lf\n\ta1=%Lf\n\ta2=%Lf\n\ta3=%Lf\n\ta4=%Lf\n\ta5=%Lf\n\ta6=%Lf\n", a0, a1, a2, a3, a4, a5, a6);

	//Roots localization in ring
	ftype arrayOneToN[6] = { a1, a2, a3, a4, a5, a6 };
	ftype arrayZeroToNMinOne[6] = { a0, a1, a2, a3, a4, a5 };

	ftype A = getMax(arrayOneToN, 6);
	ftype B = getMax(arrayZeroToNMinOne, 6);

	ftype leftEdge = fabsl(a6) / (fabsl(a6) + B);
	ftype rightEdge = 1 + (A / fabsl(a0));

	printf("\nSecond step result(localization):\nRing:\n\t%Lf<=|Z|<=%Lf\n", leftEdge, rightEdge);
	//Root separating
	ftype N = 1024;
	ftype step = (rightEdge - leftEdge) / N;
	interval intervals[6];
	int intervalsAmount = 0;
	int lastSign = sign(evaluate(leftEdge));
		
	N = N * 2;
	step = (rightEdge - leftEdge) / N;
	for (ftype X = leftEdge; X <= rightEdge; X += step) {
		if (sign(evaluate(X)) != lastSign) {
			lastSign = sign(evaluate(X));
			intervalsAmount++;
			intervals[intervalsAmount - 1].left = X - step;
			intervals[intervalsAmount - 1].right = X;
		}
	}
	
	printf("Intervals:\n");
	for (int i = 0; i < intervalsAmount; i++) {
		printInterval(intervals[i]);
	}
	//Getting roots with required precision
	printf("Third step result(clarification by dichotomy method):\nInput epsilon: ");
	ftype epsilon = 0.0;
	ftype* roots = (ftype*) malloc(intervalsAmount * sizeof(ftype));
	scanf_s("%Lf", &epsilon);
	for (int i = 0; i < intervalsAmount; i++) {
		ftype l = intervals[i].left;
		ftype r = intervals[i].right;
		while ((r - l) > epsilon) {
			ftype middle = (l + r) / 2.0;

			ftype leftValue = evaluate(l);
			ftype middleValue = evaluate(middle);

			ftype testVal = leftValue * middleValue;

			if (testVal < 0) {
				r = middle;
			}
			else {
				l = middle;
			}
		}
		ftype root = (l + r) / 2.0;
		roots[i] = root;
		printf("Root number %d Z = %Lf\n", i, root);
	}
	printf("\nFourth step result(D0):\n");
	for (int i = 0; i < intervalsAmount; i++) {
		if (roots[i] > 0) {
			printf("Root Z[%d]=%Lf\n", i, roots[i]);
			printf("---------------------------\n");
			
			ftype P1 = powl(roots[i], n) * P3;

			ftype rho1 = rho0 * ((gamma0 -1.0) + (gamma0 + 1.0) * (P1 / P0)) / ((gamma0 + 1.0) + (gamma0 - 1.0) * (P1 / P0));

			ftype D0 = U0 - sqrtl((rho1 * (P1 - P0)) / (rho0 * (rho1 - rho0)));
			
			printf("\tD0=%Lf\n", D0);
		}
	}
	system("pause");
	return 0;
}
