#include <iostream>
#include <cmath>
#include <vector>
#include "fluxes.h"
using namespace std;

int main() {
	const double gamma = 1.4;
	structFlux structRoe;

	vector<double> UL, UR(4, 2.0), n(2, 1.0);
	UL = {1, 2, 6, 7};
	structFlux result = rusanov(UL, UR, gamma, n);
	vector<double> F = result.F;
	for (double i : F) {
		cout << i << "\n";
	}
	cout << F.size() << "\n";
	return 0;
}