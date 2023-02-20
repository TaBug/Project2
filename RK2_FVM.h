//
//  RK2_FVM.h
//  p2-AEROSP-623
//
//  Created by Jake Yeaman on 2/17/23.
//

#ifndef RK2_FVM_h
#define RK2_FVM_h

#include <iostream>
#include <vector>
#include <cmath>
#include "solver.h"

using namespace std;

void rk2(int opt, vector<vector<double>> &u, vector<double> const &area ,vector<vector<double>> const &nodes, vector<vector<double>> const &elem, double Minf, double alphaDeg, vector<vector<double>> const &Bn, vector<vector<double>> const &In, vector<vector<int>> const &elemBounds, vector<vector<double>> const &bounds, vector<vector<double>> const &interiorFaces, vector<vector<int>> const &globalEdge, vector<vector<double>> const &I2E, vector<vector<double>> const &B2E, string limiterType, double convergedVal, vector<double>const &Area, int nelem, double CFL, vector<double> &L1_norms, int &niter){
    vector<vector<double>> f0(nelem,vector<double>(4));
    vector<vector<double>> f1(nelem,vector<double>(4));
    
    double residSum = DBL_MAX;
    niter = 0;

    while(residSum > convergedVal){
        niter++;
        // initialize L1 Residual and Residual to 0 every time iteration
        vector<vector<double>> residual(nelem,vector<double>(4));
        vector<vector<double>> residual2(nelem,vector<double>(4));
        double resL1 = 0;

        // calculate residual of each element
        residual = secondOrderFV(opt, u, Area, nodes, elem, Minf, alphaDeg, Bn, In, elemBounds, bounds, interiorFaces, globalEdge, I2E, B2E, limiterType);
        
        // calculate L1 Residual and stop calculation if < 10e-5
        for (int i = 0; i < nelem; i++){
            for (int j = 0; j < 4; j++){
                resL1 += abs(residual[i][j]);
            }
        }
        
        for (int i = 0; i < nelem; i++){
            bool hasNan = false;
            for (int j = 0; j < 4; j++){
                if(isnan(abs(residual[i][j]))){
                    hasNan = true;
                }
            }
            if(hasNan == true){
                cout << i + 1 << "\n";
            }
        }
        
        if (resL1 < pow(10,-5)){
                break;
        }
        
        vector<vector<double>> f0(nelem,vector<double>(4,0));
        vector<vector<double>> uf0(nelem,vector<double>(4,0));
        // first step of RK2
        for (int i = 0; i < nelem; i++){

            double dt = (2*Area[i]*CFL)/residual[i][4]; // calculate local time step
            
            // find f0 and use it to find uf0 which is used for next step of RK2
            for (int j = 0; j < 4; j++){
                f0[i][j] = -residual[i][j]/Area[i];
                uf0[i][j] = u[i][j] + dt*f0[i][j];
            }
        }

        // calculate residuals for the state uf0
        residual2 = secondOrderFV(opt, uf0, Area, nodes, elem, Minf, alphaDeg, Bn, In, elemBounds, bounds, interiorFaces, globalEdge, I2E, B2E, limiterType);
        
        double sum = 0;
        for(int i = 0; i < residual2.size(); i++){
            for(int j = 0; j < residual2[0].size() - 1;j++){
                sum += abs(residual2[i][j]);
            }
        }
        
        residSum = sum;
        L1_norms.emplace_back(residSum);
        
//        cout << residSum << "\n";
        
        if(niter % 2500 == 0){
            cout << "\n\nIteration " << niter << " residual is " << residSum << "\n\n\n";
        }
        
//        if(niter % 25000 == 0){ // asking for termination every 25000 iterations
//            cout << "Should terminate (1 == YES & 0 == NO)? ";
//            int termFlag = - 1;
//            cin >> termFlag;
//            if(termFlag == 1){
//                break;
//            }
//            else{
//                cout << "\n";
//            }
//            
//        }

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
    
    // int stopHere = -1;

}

#endif /* RK2_FVM_h */


