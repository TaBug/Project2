//
//  elem2Edge.h
//  p2-AEROSP-623
//
//  Created by Jake Yeaman on 2/10/23.
//

#ifndef elem2Edge_h
#define elem2Edge_h

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

using namespace std;

struct iGlobal2Local{
    bool isBound;
    int index;
};

vector<vector<int>> genElemBounds(vector<vector<double>> const &elem, vector<vector<double>> const &I2E, vector<vector<double>> const &B2E){
   
    int niedge = int(I2E.size());
    int nbedge = int(B2E.size());
    int ngedge = niedge + nbedge; // Number of global edges (total number of edges)
    int nelem = int(elem.size());
    
    vector<vector<int>> elemBounds(nelem,vector<int>(3,0));
    
    // Iterate through each edge and add the information into elemBounds. This will double count and overwrite elements but this is fine
    
    // Start with interior edges
    for(int iie = 0; iie < niedge; iie++){
        
        int iE1 = int(I2E[iie][0]);
        int iFace1 = int(I2E[iie][1]);
        elemBounds[iE1-1][iFace1-1] = iie;
        
        int iE2 = int(I2E[iie][2]);
        int iFace2 = int(I2E[iie][3]);
        elemBounds[iE2-1][iFace2-1] = iie;
        
    }
    
    // Now do boundary edges
    for(int ibe = 0; ibe < nbedge; ibe++){
        
        int iE = int(B2E[ibe][0]);
        int iFace = int(B2E[ibe][1]);
        elemBounds[iE-1][iFace-1] = ibe + niedge;
        
    }
    
    return elemBounds;
}

vector<vector<int>> genGlobalEdge(vector<vector<double>> const &bounds, vector<vector<double>> const &interiorFaces){
    
    int niedge = int(interiorFaces.size());
    int nbedge = int(bounds.size());
    int ngedge = niedge + nbedge; // Number of global edges (total number of edges)
    
    vector<vector<int>> globalEdge(ngedge,vector<int>(2,0));
    
    // Loop through interior edges
    for(int iie = 0; iie < niedge; iie++){
        
        globalEdge[iie][0] = interiorFaces[iie][0]; // node 1 for edge #iie
        globalEdge[iie][1] = interiorFaces[iie][1]; // node 2 for edge #iie
        
    }
    
    // Loop through boundary edges
    for(int ibe = 0; ibe < nbedge; ibe++){
        
        globalEdge[ibe + niedge][0] = bounds[ibe][0];
        globalEdge[ibe + niedge][1] = bounds[ibe][1];
        
    }
    
    return globalEdge;
    
}

iGlobal2Local iG2L(int i, vector<vector<double>> const &bounds, vector<vector<double>> const &interiorFaces){
    
    int niedge = int(interiorFaces.size());
    int nbedge = int(bounds.size());
    int ngedge = niedge + nbedge; // Number of global edges (total number of edges)
    
    if(i>=niedge){ // boundary index
        iGlobal2Local iLocal = {true,i-niedge};
        return iLocal;
    }
    else{ // interior index
        iGlobal2Local iLocal = {false,i};
        return iLocal;
    }
    
}

#endif /* elem2Edge_h */
