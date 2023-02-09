#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include "tools.h"
using namespace std;

struct struct {
	vector<double> F;
	double smag;
};

structFlux rusanov(vector<double>& UL, vector<double>& UR, double gamma, vector<double>& n, structFlux output) {
	// process left state
	double rL, uL, vL, unL, qL, pL, rHL, HL, cL;

	rL = UL[0];
	uL = UL[1] / rL;
	vL = UL[2] / rL;
	unL = uL * n[0] + vL * n[1];
	qL = sqrt(pow(UL[1], 2) + pow(UL[2], 2)) / rL;
	pL = (gamma - 1) * (UL[3] - 0.5 * rL * pow(qL, 2));
	if (pL < 0 || rL < 0) {
		cout << "Non-physical state!" << "\n";
		cout << UL[0] << "\n";
		cout << UL[1] << "\n";
	}
	if (pL < 0) {
		pL = -pL;
	}
	if (rL < 0) {
		rL = -rL;
	}
	rHL = UL[3] + pL;
	HL = rHL / rL;
	cL = sqrt(gamma * pL / rL);

	// left flux 
	vector<double> FL(4, 0.0);
	FL[0] = rL * unL;
	FL[1] = UL[1] * unL + pL * n[0];
	FL[2] = UL[2] * unL + pL * n[1];
	FL[3] = rHL * unL;

	// process right state
	double rR, uR, vR, unR, qR, pR, rHR, HR, cR;

	rR = UR[0];
	uR = UR[1] / rR;
	vR = UR[2] / rR;
	unR = uR * n[0] + vR * n[1];
	qR = sqrt(pow(UR[1], 2) + pow(UR[2], 2)) / rR;
	pR = (gamma - 1) * (UR[3] - 0.5 * rR * pow(qR, 2));
	if (pR < 0 || rR < 0) {
		cout << "Non-physical state!" << "\n";
	}
	if (pR < 0) {
		pR = -pR;
	}
	if (rR < 0) {
		rR = -rR;
	}
	rHR = UR[3] + pR;
	HR = rHR / rR;
	cR = sqrt(gamma * pR / rR);

	// right flux
	vector<double> FR(4, 0.0);
	FR[0] = rR * unR;
	FR[1] = UR[1] * unR + pR * n[0];
	FR[2] = UR[2] * unR + pR * n[1];
	FR[3] = rHR * unR;

	// difference in states
	vector<double> du;
	du = subtractVectors(UR, UL);

	// max wave speed
	double smag;
	double sLmax = uL + cL;
	double sRmax = uR + cR;
	if (sLmax > sRmax) {
		smag = sLmax;
	}
	else {
		smag = sRmax;
	}

	// flux assembly
	vector<double> F(4, 0.0);
	F[0] = .5 * (FL[0] + FR[0]) - 0.5 * smag * (UR[0] - UL[0]);
	F[1] = .5 * (FL[1] + FR[1]) - 0.5 * smag * (UR[1] - UL[1]);
	F[2] = .5 * (FL[2] + FR[2]) - 0.5 * smag * (UR[2] - UL[2]);
	F[3] = .5 * (FL[3] + FR[3]) - 0.5 * smag * (UR[3] - UL[3]);

	output.F = F;
	output.smag = smag;
	return output;
}