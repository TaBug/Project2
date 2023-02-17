#include <iostream>
#include <vector>
#include <cmath>
#include <solver.h>

using namespace std;

void rk2(vector<vector<double>> &u, vector<vector<double>> &residual,int &nelem){
    vector<vector<double>> f0(nelem,vector<double>(4));
    vector<vector<double>> f1(nelem,vector<double>(4));

    for{
        // initialize L1 Residual and Residual to 0 every time iteration
        vector<vector<double>> residual(nelem,vector<double>(4));
        vector<vector<double>> residual2(nelem,vector<double>(4));
        double resL1 = 0;

        // calculate residual of each element
        residual = secondOrderFV(opt,u,Area,nodes,elem,Minf,alphaDeg,Bn,In,elemBounds,bounds,interiorFaces,globalEdge,I2E,B2E);
        
        // calculate L1 Residual and stop calculation if < 10e-5
        for (int i = 0; i < nelem; i++){
            for (int j = 0; j < 4; j++){
                resL1 += abs(residual[i][j]);
            }
        }
        if (resL1 < pow(10,-5)){
                break;
        }
        
        // first step of RK2
        for (int i = 0; i < nelem; i++){
            double dt = (2*Area[i]*CFL)/residual[i][4]; // calculate local time step
            
            // find f0 and use it to find uf0 which is used for next step of RK2
            for (int j = 0; j < 4; j++){
                f0[i][j] = -residual[i][j]/Area[i];

                uf0[i][j] = u[i][j] + dt*f0[j];
            }
        }

        // calculate residuals for the state uf0
        residual2 = secondOrderFV();

        // second step of RK2
        for (int i = 0; i < nelem; i++){
            double dt = (2*Area[i]*CFL)/residual[i][4]; // calculate local time step
            
            // find f1 and update state using f0 and f1
            for (int j = 0; j < 4; j++){
                f1[i][j] = -residual2[i][j]/Area[i];

                u[i][j] = u[i][j] + .5*dt*(f0[i][j] + f1[i][j]);
            }
        }

    }

}