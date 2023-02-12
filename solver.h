//
//  solver.h
//  p2-AEROSP-623
//
//  Created by Jake Yeaman on 2/6/23.
//

#ifndef solver_h
#define solver_h

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_set>
#include <numeric>
#include <cmath>
#include <cfloat>
#include <algorithm>
#include "elem2Edge.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace std;
using namespace Eigen;

//double computeL1ResidualNorm(){
//
//}

vector<double> computeFreestreamState(double Minf, double alphaDeg){
    double alphaRad = alphaDeg * M_PI / 180.0;
    double gamma = 1.4;
    vector<double> uInf = {1.0, Minf*cos(alphaRad), Minf*sin(alphaRad), (1/(gamma*gamma-gamma))+0.5*(Minf*Minf)};
    return uInf;
}

/// Limiter Functions
///

//vector<double> computeNeighborIndices(){
//
//}

vector<Vector3d> computeP(vector<vector<double>> const &nodes, vector<vector<double>> const &elem, vector<double> const &U_cell, int const iCell){
    double xCentroid = (1/3)*((nodes[ elem[iCell][0]-1 ][0]) + (nodes[ elem[iCell][1]-1 ][0]) + (nodes[ elem[iCell][2]-1 ][0]));
    double yCentroid = (1/3)*((nodes[ elem[iCell][0]-1 ][1]) + (nodes[ elem[iCell][1]-1 ][1]) + (nodes[ elem[iCell][2]-1 ][1]));
    vector<Vector3d> P;
    P.reserve(U_cell.size());
    
    for(int i = 0; i < U_cell.size(); i++){
        P.push_back(Vector3d(xCentroid, yCentroid, U_cell[i]));
    }
    
    return P;
}

vector<double> computeBoundaryState(vector<vector<double>> const &nodes, vector<vector<double>> const &elem, vector<double> const &U_cell, int const iCell, const int iFace, bool isWall, double Minf, double alphaDeg, vector<vector<double>> const &Bn, int iBound){
    
    vector<double> U_ghost;
    U_ghost.reserve(4);
    
    if(isWall == true){ // Wall Boundary
        double rho = U_cell[0];
        Vector2d vCell = {U_cell[1]/rho,U_cell[2]/rho};
        Vector2d n = {Bn[iBound][0],Bn[iBound][1]};
        Vector2d vb = vCell - (vCell.dot(n))*n;
        U_ghost = {rho, rho*vb[0], rho*vb[1], U_cell[3]};
    }
    else{ // Freestream
        U_ghost = computeFreestreamState(Minf, alphaDeg);
    }
    return U_ghost;
}

vector<Vector3d> computeP_boundary(vector<vector<double>> const &nodes, vector<vector<double>> const &elem, vector<double> const &U_cell, int const iCell, const int iFace, bool isWall, double Minf, double alphaDeg, vector<vector<double>> const &Bn, int iBound){
    
    // Finding the boundary edge midpoint
    vector<size_t> nodeIndices = {0,1,2}; // representing localface 1, 2, and 3.
    nodeIndices.erase(nodeIndices.begin() + (iFace - 1)); // removing local face index so our two remaining indices correspond to the index of edge nodes
    vector<double> edgeMidpoint = {(nodes[elem[iCell][nodeIndices[0]]-1][0] + nodes[elem[iCell][nodeIndices[1]]-1][0])/2, (nodes[elem[iCell][nodeIndices[0]]-1][1] + nodes[elem[iCell][nodeIndices[1]]-1][1])/2}; // cell0 and cellk edge interface midpoint
    
    // Finding current, non-ghost cell center
    double xCentroid = (1/3)*((nodes[ elem[iCell][0]-1 ][0]) + (nodes[ elem[iCell][1]-1 ][0]) + (nodes[ elem[iCell][2]-1 ][0]));
    double yCentroid = (1/3)*((nodes[ elem[iCell][0]-1 ][1]) + (nodes[ elem[iCell][1]-1 ][1]) + (nodes[ elem[iCell][2]-1 ][1]));
    
    // Finding ghost cell center
    double dx = edgeMidpoint[0] - xCentroid;
    double dy = edgeMidpoint[1] - yCentroid;
    double xCentroidGhost = edgeMidpoint[0] + dx;
    double yCentroidGhost = edgeMidpoint[1] + dy;
    
//    vector<double> U_ghost;
//    U_ghost.reserve(4);
//
//    if(isWall == true){ // Wall Boundary
//        double rho = U_cell[0];
//        Vector2d vCell = {U_cell[1]/rho,U_cell[2]/rho};
//        Vector2d n = {Bn[iBound][0],Bn[iBound][1]};
//        Vector2d vb = vCell - (vCell.dot(n))*n;
//        U_ghost = {rho, rho*vb[0], rho*vb[1], U_cell[3]};
//    }
//    else{ // Freestream
//        U_ghost = computeFreestreamState(Minf, alphaDeg);
//    }
    
    vector<double> U_ghost = computeBoundaryState(nodes, elem, U_cell, iCell, iFace, isWall, Minf, alphaDeg, Bn, iBound);
    
    // Computing P
    vector<Vector3d> P;
    P.reserve(4);
    
    for(int iU = 0; iU < 4; iU++){
        P.emplace_back(Vector3d(xCentroidGhost,yCentroidGhost,U_ghost[iU]));
    }
    
    return P; 
    
}

// iNeighbor: a vector of neightboring cell indeces (size of two on boundaries and size of three on interior)
vector<Vector2d> computeL(vector<vector<double>> const &nodes, vector<vector<double>> const &elem, vector<vector<double>> const &U, int iCell, vector<int> const iNeighbor, vector<int> const iFaces, double Minf, double alphaDeg, vector<vector<double>> const &Bn, vector<vector<double>> const &bounds, vector<vector<double>> const &interiorFaces, vector<vector<int>> const &elemBounds){
    // iFaces: vector consisting of the local face index of the interface edge between cell0 and cellk, same ordering as iNeighbor
    
    // TODO: acount for boundary edge states: farfield = free-stream and wall = density and energy equivalent to the state inside the domain
    // Assuming that for boundary cases, L_123 are corrected to the only possible L vector value. Therefore, size should always be zero
    
    vector<Vector3d> Pi;
    vector<Vector3d> Pj;
    vector<Vector3d> Pk;
    vector<Vector3d> zerosVec(4,Vector3d::Zero());

    
    // LET iNeighbor BE NEGATIVE IFF ON A WALL BOUNDARY
    // TODO: IF iNeighbor is NEGATIVE, corresponding iFace should be the local face of iCell. Ensure this is passed in correctly to the code
    // If not on a wall boundary, compute P with the neighboring cell coordinates and neighboring cell states
    if(iNeighbor[0] < 0){
        
        int iGlobal = elemBounds[iCell][iFaces[0]-1];
        iGlobal2Local currFace = iG2L(iGlobal, bounds, interiorFaces);
        bool isWall = false;
        if(bounds[currFace.index][2] == 4){ // TODO: double check to make sure wall boundary groups are equal to 4
            isWall = true;
        }
        Pi = computeP_boundary(nodes, elem, U[iCell], iCell, iFaces[0], isWall, Minf, alphaDeg, Bn, currFace.index);
        
    }
    
    else{
        Pi = computeP(nodes, elem, U[iNeighbor[0]], iNeighbor[0]);
    }
    
    if(iNeighbor[1] < 0){
        int iGlobal = elemBounds[iCell][iFaces[1]-1];
        iGlobal2Local currFace = iG2L(iGlobal, bounds, interiorFaces);
        bool isWall = false;
        if(bounds[currFace.index][2] == 4){ // TODO: double check to make sure wall boundary groups are equal to 4
            isWall = true;
        }
        Pj = computeP_boundary(nodes, elem, U[iCell], iCell, iFaces[1], isWall, Minf, alphaDeg, Bn, currFace.index);
    }
    
    else{
        Pj = computeP(nodes, elem, U[iNeighbor[1]], iNeighbor[1]);
    }
    
    if(iNeighbor[2] < 0){
        int iGlobal = elemBounds[iCell][iFaces[2]-1];
        iGlobal2Local currFace = iG2L(iGlobal, bounds, interiorFaces);
        bool isWall = false;
        if(bounds[currFace.index][2] == 4){ // TODO: double check to make sure wall boundary groups are equal to 4
            isWall = true;
        }
        Pk = computeP_boundary(nodes, elem, U[iCell], iCell, iFaces[2], isWall, Minf, alphaDeg, Bn, currFace.index);
        
    }
    
    else{
        Pk = computeP(nodes, elem, U[iNeighbor[2]], iNeighbor[2]);
    }
    
    
    vector<Vector3d> n;
    n.reserve(Pi.size());
    
    for(int i = 0; i < Pi.size(); i++){
        n.emplace_back(Vector3d((Pj[i] - Pi[i]).cross(Pk[i] - Pi[i])));
    }
    
    vector<Vector2d> Lijk;
    Lijk.reserve(Pi.size());
    
    for(int i = 0; i < Pi.size(); i++){
        Lijk.emplace_back(Vector2d(-n[i][0]/n[i][2], -n[i][1]/n[i][2]));
    }

    
    return Lijk;
    
}

vector<Vector2d> computeL_LCD(vector<Vector2d> L, vector<vector<double>> const &U, vector<vector<double>> const &nodes, vector<vector<double>> const &elem, int iCell, vector<int> const iNeighbor, vector<int> const iFaces, double Minf, double alphaDeg, vector<vector<double>> const &Bn, vector<vector<int>> const &elemBounds, vector<vector<double>> const &bounds, vector<vector<double>> const &interiorFaces){
    // TODO: IF iNeighbor is NEGATIVE, corresponding iFace should be the local face of iCell. Ensure this is passed in correctly to the code
    // iNeighbor: vector consisting of the index of the neighboring elements
    // iFaces: vector consisting of the local face index of the interface edge between cell0 and cellk, same ordering as iNeighbor
    
    vector<Vector2d> Llcd;
    Llcd.reserve(U[0].size());
    
    vector<double> alpha(U[0].size(),DBL_MAX); // storing an independent alpha value for each state variable
    vector<double> cellCentroid = {(nodes[elem[iCell][0]-1][0] + nodes[elem[iCell][1]-1][0] + nodes[elem[iCell][2]-1][0])/3,(nodes[elem[iCell][0]-1][1] + nodes[elem[iCell][1]-1][1] + nodes[elem[iCell][2]-1][1])/3};
    
    for(int i = 0; i < 3; i++){
        
        vector<size_t> nodeIndices = {0,1,2}; // representing localface 1, 2, and 3.
        nodeIndices.erase(nodeIndices.begin() + (iFaces[i] - 1)); // removing local face index so our two remaining indices correspond to the index of edge nodes

        int ik = iNeighbor[i]; // cellk index
        vector<double> edgeMidpoint;
        
        if(ik < 0){ // if on a boudnary
            edgeMidpoint = {(nodes[elem[iCell][nodeIndices[0]]-1][0] + nodes[elem[iCell][nodeIndices[1]]-1][0])/2, (nodes[elem[iCell][nodeIndices[0]]-1][1] + nodes[elem[iCell][nodeIndices[1]]-1][1])/2}; // cell0 and boundary edge interface midpoint
        }
        else{ // If not on a boundary
            edgeMidpoint = {(nodes[elem[ik][nodeIndices[0]]-1][0] + nodes[elem[ik][nodeIndices[1]]-1][0])/2, (nodes[elem[ik][nodeIndices[0]]-1][1] + nodes[elem[ik][nodeIndices[1]]-1][1])/2}; // cell0 and cellk edge interface midpoint
        }
        
        Vector2d r0k(edgeMidpoint[0] - cellCentroid[0],edgeMidpoint[1] - cellCentroid[1]);
        
        for(int j = 0; j < alpha.size(); j++){ // populating alpha for each state
         
            double alpha_k;
            double stateDifference;
            double rDotL;
            
            if(ik < 0){ // if on a boundary
                int iGlobal = elemBounds[iCell][iFaces[i]-1];
                iGlobal2Local currFace = iG2L(iGlobal, bounds, interiorFaces);
                bool isWall = false;
                if(bounds[currFace.index][2] == 4){ // TODO: double check to make sure wall boundary groups are equal to 4
                    isWall = true;
                }
                vector<double> Ub = computeBoundaryState(nodes, elem, U[iCell], iCell, iFaces[i], isWall, Minf, alphaDeg, Bn, currFace.index);
                stateDifference = Ub[j] - U[iCell][j];
                rDotL = r0k.dot(L[j]);
            }
            else{ // if not on a boundary
                stateDifference = U[ik][j] - U[iCell][j];
                rDotL = r0k.dot(L[j]);
            }
            
            if(rDotL > max(stateDifference,0.0)){
                alpha_k = max(stateDifference,0.0)/rDotL;
            }
            else if(r0k.dot(L[j]) < min(stateDifference,0.0)){
                alpha_k = min(stateDifference,0.0)/rDotL;
            }
            else{
                alpha_k = 1;
            }
            
            // Updating alpha if alphak < alphaCurrent
            
            if(alpha_k < alpha[j]){
                alpha[j] = alpha_k;
            }
        }
        
    }
    
    for(int i = 0; i < alpha.size(); i++){
        Llcd.emplace_back(alpha[i]*L[i]);
    }
    
    return Llcd;
    
}

/// BARTH JESPERSON

vector<Vector2d> compute_rN(vector<vector<double>> const &nodes, vector<vector<double>> const &elem, int iCell){
    
    vector<double> cellCentroid = {(nodes[elem[iCell][0]-1][0] + nodes[elem[iCell][1]-1][0] + nodes[elem[iCell][2]-1][0])/3,(nodes[elem[iCell][0]-1][1] + nodes[elem[iCell][1]-1][1] + nodes[elem[iCell][2]-1][1])/3};
    
    vector<Vector2d> rN(3,Vector2d(2,0.0));
    
    // Iterate through each node
    for(int i = 0; i < 3; i++){
        double xN = nodes[elem[iCell][i]-1][0];
        double yN = nodes[elem[iCell][i]-1][1];
        rN[i][0] = xN - cellCentroid[0]; // xComponent of rN
        rN[i][1] = yN - cellCentroid[1]; // yComponent of rN
    }
    
    return rN;
    
}

vector<Vector2d> barthJespersen(vector<vector<double>> const &nodes, vector<vector<double>> const &elem, vector<double> const &area, vector<vector<double>> const &U, int iCell, vector<int> const iNeighbor){
    // vector<vector<double>> uALL(4, vector<double>(4,0.0)); // EACH ROW REPRESENTS A SPECIFIC STATE VARIABLE (E.G DENSTIY)
    vector<vector<double>> uALL(4,vector<double>(1 + iNeighbor.size(),0.0));
    
    // Populating uALL for iCell
    uALL[0][0] = U[iCell][0];
    uALL[1][0] = U[iCell][1];
    uALL[2][0] = U[iCell][2];
    uALL[3][0] = U[iCell][3];
    
    // Populating uALL for neighbors
    for(int i = 0; i < iNeighbor.size(); i++){
        uALL[0][i+1] = U[iNeighbor[i]][0];
        uALL[1][i+1] = U[iNeighbor[i]][1];
        uALL[2][i+1] = U[iNeighbor[i]][2];
        uALL[3][i+1] = U[iNeighbor[i]][3];
    }
    
    vector<double> uMAX(4,0.0);
    vector<double> uMIN(4,0.0);
    
    // Iterate through each state variable
    for(int i = 0; i < 4; i++){
        uMAX[i] = *(max_element(uALL[i].begin(), uALL[i].end()));
        uMIN[i] = *(min_element(uALL[i].begin(), uALL[i].end()));
    }
    
    // Find cell0 state
    vector<double> u0 = U[iCell];
    
    // Obtain all rNs
    vector<Vector2d> rN = compute_rN(nodes, elem, iCell);
    
    // Optain all L
    vector<Vector2d> L = computeL(nodes, elem, U, iCell, iNeighbor);
    
    // Obtain node states
    // Each row is a node and each col is a state variable
    vector<vector<double>> uiN;
    uiN.reserve(3);
    
    // Iterate through each node
    for(int iN = 0; iN < 3; iN++){
        vector<double> currState;
        currState.reserve(4);
        // Iterate through each state var.
        for(int iU = 0; iU < 4; iU++){
            currState.emplace_back(u0[iU] + rN[iN].dot(L[iU]));
        }
        uiN.emplace_back(currState);
    }
    
    // Computing alphaNs
    vector<double> alphas(4, DBL_MAX); // Final, minimul alpha values
    
    for(int iN = 0; iN < 3; iN++){ // Iterate through each node
        
        for(int iU = 0; iU < 4; iU++){ // Iterate through each state variable
            
            if((uiN[iN][iU]-u0[iU])>0){
                
                double alpha_iN = min(1.0,(uMAX[iU]-u0[iU])/(uiN[iN][iU]-u0[iU]));
                
                if(alpha_iN < alphas[iU]){
                    alphas[iU] = alpha_iN;
                }
                
            }
            else if ((uiN[iN][iU]-u0[iU])<0){
                
                double alpha_iN = max(1.0,(uMAX[iU]-u0[iU])/(uiN[iN][iU]-u0[iU]));
                
                if(alpha_iN < alphas[iU]){
                    alphas[iU] = alpha_iN;
                }
                
            }
            else{
                
                double alpha_iN = 1.0;
                
                if(alpha_iN < alphas[iU]){
                    alphas[iU] = alpha_iN;
                }
            }
            
        } // end loop for each state variable
    } // end loop for each node
    
    // Scale the limeters
    for(int iU = 0; iU < 4; iU++){
        
        L[iU] *= alphas[iU];
        
    }
    
    return L;
    
}




///


// 2nd Order Finite Volume Driver
void secondOrderFV( vector<vector<double>> &U, int nelem, int nfaces){
    
    double residual;;
    // double residual = computeL1ResidualNorm();
    
    while(abs(residual) > pow(10,-5)){
        
        vector<double> grad_u(nelem,0);
        for(int i = 0; i < nfaces; i++){
            
            
            
        }
        
        
    }
    
}


#endif /* solver_h */
