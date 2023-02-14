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
#include "fluxes.h"

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

vector<Vector2d> barthJespersen(vector<vector<double>> const &nodes, vector<vector<double>> const &elem, vector<double> const &area, vector<vector<double>> const &U, int iCell, vector<int> const iNeighbor, double Minf, double alphaDeg, vector<vector<double>> const &Bn, vector<vector<int>> const &elemBounds, vector<vector<double>> const &bounds, vector<vector<double>> const &interiorFaces, vector<int> const iFaces){
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
    vector<Vector2d> L = computeL(nodes, elem, U, iCell, iNeighbor, iFaces, Minf, alphaDeg, Bn, bounds, interiorFaces, elemBounds);
    
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


vector<int> findAdjElem(vector<vector<double>> const &elem, vector<vector<double>> const &I2E, vector<vector<double>> const &B2E, vector<vector<int>> const &elemBounds, vector<vector<double>> const &bounds, vector<vector<double>> const &interiorFaces, vector<vector<int>> const &globalEdge, int iFaceGlobal){
    
    // COMPUTES THE L and R elems for an input face (if no adj element exists due to being on a boundary, negative value is used)
    iGlobal2Local iLocal = iG2L(iFaceGlobal, bounds, interiorFaces);
    vector<int> elemLR = {-1, -1};
    
    if(iLocal.isBound == false){ // interior edge
        elemLR[0] = int(I2E[iLocal.index][0] - 1); // Left Element
        elemLR[1] = int(I2E[iLocal.index][1] - 1); // Right Element
    }
    
    else{ // boundary edge
        elemLR[0] = int(B2E[iLocal.index][0] - 1); // Left Element
        // No right element exists
    }
    
    return elemLR; // RETURNING L and R ELEMENT TO EDGE IN BASE-0 INDEXING
    
}


vector<int> findNeighbors(vector<vector<double>> const &elem, vector<vector<double>> const &I2E, vector<vector<double>> const &B2E, vector<vector<int>> const &elemBounds, vector<vector<double>> const &bounds, vector<vector<double>> const &interiorFaces, vector<vector<int>> const &globalEdge, int iElem){
    // iElem is 0-based
    unordered_set<int> neighbors;
    
    // finding all edges surrounding an element
    vector<int> edges = {elemBounds[iElem][0], elemBounds[iElem][1], elemBounds[iElem][1]};
    
    // populating the neighbor set
    for(int iEdge = 0; iEdge < 3; iEdge++){
        neighbors.insert(edges[iEdge]);
    }
    
    // Populating vector containing all neighbor indices
    vector<int> iNeighbors;
    iNeighbors.insert(iNeighbors.end(), neighbors.begin(), neighbors.end());
    
    // If we are on a boundary, some edges do not have an adjacent element, therefore indicate this with negative values
    if(iNeighbors.size() == 1){
        iNeighbors.push_back(-1);
        iNeighbors.push_back(-1);
    }
    else if (iNeighbors.size() == 2){
        iNeighbors.push_back(-1);
    }
    
    return iNeighbors;
    
}

double computeEdgeLength(int iEdgeGlobal, vector<vector<int>> const &globalEdge, vector<vector<double>> const &nodes){
    double n1x = nodes[globalEdge[iEdgeGlobal][0] - 1][0];
    double n1y = nodes[globalEdge[iEdgeGlobal][0] - 1][1];
    double n2x = nodes[globalEdge[iEdgeGlobal][1] - 1][0];
    double n2y = nodes[globalEdge[iEdgeGlobal][1] - 1][1];
    
    double deltaL = sqrt(pow(n2x-n1x, 2) + pow(n2y-n1y, 2));
    return deltaL;
}


// 2nd Order Finite Volume Driver
void secondOrderFV(int opt, vector<vector<double>> &U, vector<double> const &area ,vector<vector<double>> const &nodes, vector<vector<double>> const &elem, double Minf, double alphaDeg, vector<vector<double>> const &Bn, vector<vector<double>> const &In, vector<vector<int>> const &elemBounds, vector<vector<double>> const &bounds, vector<vector<double>> const &interiorFaces, vector<vector<int>> const &globalEdge, vector<vector<double>> const &I2E, vector<vector<double>> const &B2E){
    
    int nelem = int(elem.size());
    int nfaces = int(globalEdge.size());
    double residual;;
    // double residual = computeL1ResidualNorm();
    
    vector<double> F(4);
    double s;
    
    // set the flux to use
//    structFlux (*flux)(vector<double>& UL, vector<double>& UR, double gamma, vector<double>& n);
    structFlux output;
    
//    if(opt == 1){
//        flux = roe;
//    }
//
//    else if (opt == 2){
//        flux = rusanov;
//    }
//
//    else if (opt == 3){
//        flux = hlle;
//    }
    
    while(abs(residual) > pow(10,-5)){ // TODO: Change
        
        Vector2d zeros = {0.0,0.0};
        vector<vector<Vector2d>> grad_u(nelem,vector<Vector2d>(4,zeros));
        
        for(int iElem = 0; iElem < nelem; iElem++){
            
            vector<int> iNeighbor;
            // Finding iNeighbor (all neighboring elements) and iFaces (local indices of edges adjacent to neighbor)
            vector<int> globalFaceIndices = {elemBounds[iElem][0], elemBounds[iElem][1], elemBounds[iElem][2]}; // global face indices for edges with local face 1, 2, and 3, respectively
            
            for(int iFace = 0; iFace < 3; iFace++){
                vector<int> adjElems = findAdjElem(elem, I2E, B2E, elemBounds, bounds, interiorFaces, globalEdge, globalFaceIndices[iFace]);
                auto it = find(adjElems.begin(), adjElems.end(), iElem);
                adjElems.erase(it);
                if(!adjElems.empty()){
                    iNeighbor.push_back(adjElems[0]);
                }
                else{
                    iNeighbor.push_back(-1);
                }
            }
            
            vector<int> iFaces = {1,2,3};
            
            
            for(int iF = 0; iF < 3; iF++){ // Iterating through each element's three faces
                
                // Generating appropriate indices
                int iFaceGlobal = elemBounds[iElem][iF];
                iGlobal2Local iLocal = iG2L(iFaceGlobal, bounds, interiorFaces);
                int iFaceLocal = iLocal.index;
                
                // Finding the Left and Right elements
                vector<int> elemLR = findAdjElem(elem, I2E, B2E, elemBounds, bounds, interiorFaces, globalEdge, iFaceGlobal);
                
                // Initializing the normal vector pointing from L to R
                Vector2d n;
                
                if(iLocal.isBound == false){ // interior edge
                    n = {In[iFaceLocal][0], In[iFaceLocal][1]};
                }
                else{ // boundary edge
                    n = {Bn[iFaceLocal][0], Bn[iFaceLocal][1]};
                }
                
                // Gradient calculation
                
                
                /// BARTH JESPERSEN
                
                // vector<Vector2d> L = barthJespersen(nodes, elem, area, U, iElem, iNeighbor, Minf, alphaDeg, Bn, elemBounds, bounds, interiorFaces, iFaces)
                
                /// MP LIMITER (LCD)
                vector<Vector2d> L_input = computeL(nodes, elem, U, iElem, iNeighbor, iFaces, Minf, alphaDeg, Bn, bounds, interiorFaces, elemBounds);
                
                vector<Vector2d> L = computeL_LCD(L_input, U, nodes, elem, iElem, iNeighbor, iFaces, Minf, alphaDeg, Bn, elemBounds, bounds, interiorFaces);
                
                /// NO LIMITER
                
                vector<Vector2d> L_TEMP = L_input;
                /// END GRADIENT COMPUTATION
                
                // Compute face length
                double delta_l = computeEdgeLength(iFaceGlobal, globalEdge, nodes);
                
                // Find u_hat (the average of the L and R cell averages)
                    // TODO: do this
                vector<double> UL_limiting = U[elemLR[0]];
                vector<double> UR_limiting = U[elemLR[1]];
                vector<double> Uhat;
                Uhat.reserve(UL_limiting.size());
                
                for(int iU = 0; iU < UL_limiting.size(); iU++){
                    Uhat.emplace_back(0.5*(UL_limiting[iU]+UR_limiting[iU]));
                }
                
                
                // Add u_hat*n*delta_l to gradient value in vector for L elem and subtract this for the right element (if existing)
                    // TODO: do this
                
                for(int iU = 0; iU < UL_limiting.size(); iU++){
                    Uhat.emplace_back(0.5*(UL_limiting[iU]+UR_limiting[iU]));
                }
                
                // Adding u_hat*n*delta_l to elemL and subtracting from elemR
                for(int iU = 0; iU < UL_limiting.size(); iU++){
                    int iElemL = elemLR[0];
                    Vector2d gradU_add = {n[0]*Uhat[iU]*delta_l,n[1]*Uhat[iU]*delta_l};
                    grad_u[iElemL][iU] += gradU_add;
                    
                    if(iLocal.isBound == false){ // if right element exists
                        int iElemR = elemLR[1];
                        grad_u[iElemR][iU] -= gradU_add;
                    }
    
                }
                
                
            } // end for iFace
            
        } // end for iElem
        
        // Dividing each grad_u by its element area
        for(int iElem = 0; iElem < nelem; iElem++){
            for(int iU = 0; iU < 4; iU++){
                grad_u[iElem][iU] /= area[iElem];
            }
        }
        
        // Initialize Residual and Wave Speed on Each Cell to be Zero
        vector<vector<double>> residuals(nelem, vector<double>(5,0)); // residual for each state and then the wave speed
        
        
        for(int iElem = 0; iElem < nelem; iElem++){
            for(int iF = 0; iF < 3; iF++){
                
                // Generating appropriate indices
                int iFaceGlobal = elemBounds[iElem][iF];
                iGlobal2Local iLocal = iG2L(iFaceGlobal, bounds, interiorFaces);
                int iFaceLocal = iLocal.index;
                
                // Finding the Left and Right elements
                vector<int> elemLR = findAdjElem(elem, I2E, B2E, elemBounds, bounds, interiorFaces, globalEdge, iFaceGlobal);
                
                // Initializing the normal vector pointing from L to R
                vector<double> n;
                
                // Computing normal vectors of vector<double> type
                if(iLocal.isBound == false){ // interior edge
                    n = {In[iFaceLocal][0], In[iFaceLocal][1]};
                }
                
                else{ // boundary edge
                    n = {Bn[iFaceLocal][0], Bn[iFaceLocal][1]};
                }
                
                // Compute face length
                double length = computeEdgeLength(iFaceGlobal, globalEdge, nodes);
                
                // Computing face midpoint vector and cell centroid vector
                // midpoint
                double n1x = nodes[globalEdge[iFaceGlobal][0] - 1][0];
                double n1y = nodes[globalEdge[iFaceGlobal][0] - 1][1];
                double n2x = nodes[globalEdge[iFaceGlobal][1] - 1][0];
                double n2y = nodes[globalEdge[iFaceGlobal][1] - 1][1];
                
                Vector2d midpoint = {0.5*(n1x+n2x),0.5*(n1y+n2y)};
                
                // centroid for left element
                int iElemL = elemLR[0];
                double xCentroidL = (1/3)*((nodes[ elem[iElemL][0]-1 ][0]) + (nodes[ elem[iElemL][1]-1 ][0]) + (nodes[ elem[iElemL][2]-1 ][0]));
                double yCentroidL = (1/3)*((nodes[ elem[iElemL][0]-1 ][1]) + (nodes[ elem[iElemL][1]-1 ][1]) + (nodes[ elem[iElemL][2]-1 ][1]));
                
                Vector2d centroidL = {xCentroidL, yCentroidL};
                
                // centroid for right element
                Vector2d centroidR;
                
                if(elemLR.size()>1){ // not on a boundary
                    int iElemR = elemLR[1];
                    double xCentroidR = (1/3)*((nodes[ elem[iElemR][0]-1 ][0]) + (nodes[ elem[iElemR][1]-1 ][0]) + (nodes[ elem[iElemR][2]-1 ][0]));
                    double yCentroidR = (1/3)*((nodes[ elem[iElemR][0]-1 ][1]) + (nodes[ elem[iElemR][1]-1 ][1]) + (nodes[ elem[iElemR][2]-1 ][1]));
                    centroidR = {xCentroidR, yCentroidR};
                }
                else{
                    double dx = midpoint[0] - xCentroidL;
                    double dy = midpoint[1] - yCentroidL;
                    double xCentroidGhost = midpoint[0] + dx;
                    double yCentroidGhost = midpoint[1] + dy;
                    centroidR = {xCentroidGhost, yCentroidGhost};
                }
                
                // Compute UL at the face midpoint
                
                vector<double> UL = U[iElemL];
                
                for(int iU = 0; iU < 4; iU++){
                    double UL_add = grad_u[iElemL][iU].dot(midpoint-centroidL);
                    UL[iU]+=UL_add;
                }
                
                // Compute UR at the face midpoint
                vector<double> UR;
                if(elemLR.size()>1){ // not on a boundary
                    int iElemR = elemLR[1];
                    UR = U[iElemR];
                    for(int iU = 0; iU < 4; iU++){
                        double UR_add = grad_u[iElemR][iU].dot(midpoint-centroidR);
                        UR[iU]+=UR_add;
                    }
                }
                else{
                    int dummy = 0;
                    bool isWall = false;
                    if(abs(n1x) < 100 && abs(n2x) < 100 && abs(n1y) < 100 && abs(n2y) < 100){
                        isWall = true;
                    }
                    UR = computeBoundaryState(nodes, elem, U[iElem], iElem, dummy, isWall, Minf, alphaDeg, Bn, iLocal.index);
                    for(int iU = 0; iU < 4; iU++){
                        double UR_add = grad_u[iElemL][iU].dot(midpoint-centroidR); // same grad_u as left element ????
                        UR[iU]+=UR_add;
                    }
                }
                
                // Call appropriate flux function
                
                if(opt == 1){
                    output = roe(UL, UR, 1.4, n);
                }
            
                else if (opt == 2){
                    output = rusanov(UL, UR, 1.4, n);
                }
            
                else if (opt == 3){
                    output = HLLE(UL, UR, 1.4, n);
                }
                
                vector<double> F = output.F;
                double s = output.smag;
                
                // Increment and decrement the residuals
                for(int j = 0; j < 4; j++){
                    residuals[iElemL][j] += F[j]*length;
                }
                residuals[iElemL][4] += s*length;
                
                if(elemLR.size() > 1){
                    int iElemR = elemLR[1];
                    for(int j = 0; j < 4; j++){
                        residuals[iElemR][j] -= F[j]*length;
                    }
                    residuals[iElemR][4] -= s*length; // I dont know what it means to add wave speed tallies on L and R cells
                }
            } // end for iElem
        } // end for iFace
        
    }
    
}


#endif /* solver_h */
