#include <iostream>
#include <cmath>
#include <vector>
#include "fluxes.h"
using namespace std;

int main() {
	const double gamma = 1.4;
	structFlux structRoe;

	vector<double> UL, UR, n;
	UL = {0.9, 0.5, 0.1, 2};
	UR = {1, 0.6, 0.15, 2.3};
	n = {0.5, sqrt(3)/2};
	structFlux result = HLLE(UL, UR, gamma, n);
	vector<double> F = result.F;
	for (double i : F) {
		cout << i << "\n";
	}
	double smag = result.smag;
	cout << smag << "\n";
	return 0;
}