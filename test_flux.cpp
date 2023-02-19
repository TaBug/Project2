#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <vector>
#include "fluxes.h"
using namespace std;

void unitTest() {
	double rho0, a0, M, gamma, alpha, p0;

	rho0 = 1; a0 = 1; M = 0.1; gamma = 1.4; alpha = 8 * M_PI / 180;
	p0 = pow(a0, 2) * rho0 / gamma;

	// calculate initial state
	double rho, p, nvMag, u, v, rhoE;
	rho = rho0 * pow((1 + (gamma - 1) / 2 * pow(M, 2)), -1 / (gamma - 1));
	p = p0 / pow((1 + (gamma - 1) / 2 * pow(M, 2)), gamma / (gamma - 1));
	nvMag = M * a0;
	u = nvMag * cos(alpha);
	v = nvMag * sin(alpha);
	rhoE = p / (gamma - 1) + .5 * rho * (pow(u, 2) + pow(v, 2));
	vector<double> U1 = { rho, rho * u, rho * v, rhoE };

	vector<vector<double>> output(3, vector<double>(2));
	// left BC (Roe Flux)
	vector<double> nleft = { -1, 0 };
	structFlux roeFlux = roe(U1, U1, gamma, nleft);
	output[0] = roeFlux.F;

	// bot BC (wall flux)
	vector<double> nbot = { 0, -1 };
	structFlux wallflux = wallFlux(U1, nbot, gamma);
	output[1] = wallflux.F;

	// diagonal BC (Roe Flux)
	vector<double> ndiag = { 1 / sqrt(2), 1 / sqrt(2) };
	structFlux diagFlux = roe(U1, U1, gamma, ndiag);
	output[2] = diagFlux.F;
	for (int i = 0; i < output.size(); i++) {
		cout << "[" << " ";
		for (int j = 0; j < output[i].size(); j++) {
			cout << output[i][j] << " ";
		}
		cout << "]" << "\n";
	}
}

int main() {
	unitTest();
	/*
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
	*/
	return 0;
}