#include <iostream>
#include <vector>
#include <cmath>
#include "fluxes.h"
#include "getRes.h"

using namespace std;

void FVM_1st(vector<vector<double>>& bounds, vector<vector<double>>& nodes, vector<vector<double>>& interiorFaces, vector<vector<double>>& u,
    vector<vector<double>>& B2E, vector<vector<double>>& Bn, vector<vector<double>>& In, int& nelem) {
    double CFL;
    int opt;
    int time;

    // User input for which Flux Function
    cout << "Choose Flux Function:\n" << "1 - Roe\n" << "2 - Rusanov\n" << "3 - HLLE\n" << "Enter Option: ";
    cin >> opt;

    // User input for CFL value
    cout << "Enter CFL: ";
    cin >> CFL;

    // User input for maximum time iterations
    cout << "Enter Max. Time Iterations: ";
    cin >> time;

    // Loop through for maximum number of time iterations (100,000)
    timeStep: for (int t = 1; t < time; t++) {
        // initialize L1 Residual and Residual to 0 every time iteration
        vector<vector<double>> residual(nelem, vector<double>(4));
        double resL1 = 0;

        // calculate residual of each element
        residual = getRes(nelem, opt, u, interiorFaces, In, bounds, Bn, B2E, nodes);

        // calculate L1 Residual and stop calculation if < 10e-5
        for (int i = 0; i < nelem; i++) {
            for (int j = 0; j < 4; j++) {
                resL1 += abs(residual[i][j]);
            }
        }
        //cout << resL1 << "\n";
        if (resL1 < pow(10, -5)) {
            break;
        }

        // update state using Forward Euler (first order accurate)
        for (int i = 0; i < nelem; i++) {
            double dt = (2 * CFL) / residual[i][4]; // calculate local time step
            for (int j = 0; j < 4; j++) {
                u[i][j] = u[i][j] - (dt * residual[i][j]); // update state
            }
        }
        //cout << t << "\n";
    }

return;
}