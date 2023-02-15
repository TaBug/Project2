#include <iostream>
#include <vector>
#include <cmath>
#include "fluxes.h"

using namespace std;
vector<vector<double>> getRes(int &nelem,int &opt,vector<vector<double>> &u, vector<vector<double>> &interiorFaces,vector<vector<double>> &In,vector<vector<double>> &bounds,vector<vector<double>> &Bn,vector<vector<double>> &B2E,vector<vector<double>> &nodes){
    vector<vector<double>> residual(nelem,vector<double>(5));
    vector<double> F(4);
    double s;
    vector<double> n(2);
    vector<double> F(4);

    // set flux function based on user choice
    /*structFlux (*flux)(vector<double>& UL, vector<double>& UR, double gamma, vector<double>& n);
    structFlux output;
    if (opt == 1){
        flux = roe;
    }
    else if (opt == 2){
        flux = rusanov;
    }
    else if (opt == 3){
        flux = hlle;
    }*/

    // loop over interior edges
    for (int i = 0; i < interiorFaces.size(); i++){
            int elemL = interiorFaces[i][2] - 1; // set element left of edge
            int elemR = interiorFaces[i][3] - 1; // set element right of edge

            // initialize left/right states
            vector<double> uL(4);
            vector<double> uR(4);
            
            // get normal vector for edge
            n[0] = In[i][0];
            n[1] = In[i][1];

            for (int j = 0; j < 4; j++){
                uL[j] = u[elemL][j]; // set left state
                uR[j] = u[elemR][j]; // set right state
            }
            
            // calculate length of edge
            int node1 = interiorFaces[i][0] - 1;
            int node2 = interiorFaces[i][1] - 1;\
            double length = sqrt(pow(nodes[node1][0]+nodes[node2][0],2) + pow(nodes[node1][1]+nodes[node2][1],2));

            // call chosen flux function to compute flux and wave speed using
            // left state, right state, gamma (1.4), and normal vector
            if (opt == 1){
                output = roe(uL,uR,1.4,n);
            }
            else if (opt == 2){
                output = rusanov(uL,uR,1.4,n);
            }
            else if (opt == 3){
                output = HLLE(uL,uR,1.4,n);
            }
            F = output.F;
            s = output.smag;

            // add F*length to residual of left element, subtract it from right element
            for (int j = 0; j < 4; j++){
                residual[elemL][j] += F[j]*length;
                residual[elemR][j] -= F[j]*length;
            }
            // add to the wave speed of left/right elements
            residual[elemL][4] += s*length;
            residual[elemR][4] += s*length;
    }
    // Loop over boundary edges
    for (int i = 0; i < B2E.size(); i++){
        int elem = B2E[i][0] - 1; // set element connected to boundary

        vector<double> uTemp(4); // initialize temporary holder for state vector


        // get normal vector for edge
        n[0] = Bn[i][0];
        n[1] = Bn[i][1];

        // set temporary holder to element state
        for (int j = 0; j < 4; j++){
            uTemp[j] = u[elem][j];
        }
        
        // calculate length of edge
        int node1 = bounds[i][0] - 1;
        int node2 = bounds[i][1] - 1;
        double length = sqrt(pow(nodes[node1][0]+nodes[node2][0],2) + pow(nodes[node1][1]+nodes[node2][1],2));

        // call flux function depending on boundary type, computes Flux and wave speed
        if (B2E[i][2] == 4 && (nodes[node1][0] == -100 && nodes[node2][0] == -100)){
            output = inflowFlux(uTemp,n);
        }
        else if (B2E[i][2] == 4 && (nodes[node1][0] != -100 && nodes[node2][0] != -100)){
            output = outFlowFlux(uTemp,n);
        }
        else {
            output = wallFlux(uTemp,n);
        }
        F = output.F;
        s = output.smag;

        // add F*length to residual
        for (int j = 0; j < 4; j++){
            residual[elem][j] += F[j]*length;
        }
        // add wave speed
        residual[elem][4] += s*length;
    }
    return residual;
}
