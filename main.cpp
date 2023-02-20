//
//  main.cpp
//  p2-AEROSP-623
//
//  Created by Jake Yeaman on 2/6/23.
//

#include "solver.h"
#include "processMesh.h"
#include "FE_FVM.h"
#include "RK2_FVM.h"
#include "elem2Edge.h"



using namespace std;

int main(){
    
    // TODO: FILL
    int opt = 2;
    double CFL = 1;
    double Minf = 0.25;
    double alphaDeg = 8;
 
    cout << "Computational Fluid Dynamics II Driver...\n";
    int nnode, nbedge, nelem;
    string fileName = "/Users/jakeyeaman/Desktop/Comp Fluid Dyn II/Project-1/CoarseMesh";
    getMatSizes(fileName, nnode,nbedge,nelem);
        
    int niedge = 0.5*(3*nelem - nbedge);
    
    vector<vector<double>> nodes(nnode, vector<double>(2));
    vector<vector<double>> elem(nelem, vector<double>(3)); // [node1, node2, node3]
    vector<vector<double>> bounds(nbedge, vector<double>(3)); // [node1, node2, boundaryGroup]
    
    readGmshFile(fileName,nnode,nbedge,nelem,nodes,elem,bounds);
    generateGri("/Users/jakeyeaman/Desktop/Comp Fluid Dyn II/Project-2/CoarseMesh.gri", nnode, nbedge, nelem, nodes, elem, bounds);
    vector<vector<double>> interiorFaces = genInteriorFaceVec(niedge, nelem, nbedge, bounds, elem);
        
    niedge = int(interiorFaces.size());

    // GENERATING TASK 2 DATA STRUCTURES
    vector<vector<double>> I2E = genI2E(niedge, nelem, interiorFaces, elem);
    vector<vector<double>> B2E = genB2E(nbedge, nelem, bounds, elem);
    vector<vector<double>> In = genIn(niedge, interiorFaces, elem, nodes, I2E);
    vector<vector<double>> Bn = genBn(nbedge, bounds, nodes, elem, B2E);
    vector<double> Area = genArea(nelem, elem, nodes);

    vector<vector<double>> u(nelem, vector<double> (4));
    for(int i=0;i<nelem;++i)  {
        u[i][0]=1;
        u[i][1]=.25*cos(8*3.14159/180);
        u[i][2]=.25*sin(8*3.14159/180);
        u[i][3]=1/(.4*1.4)+.25*.25/2;
    }

    FVM_1st(bounds, nodes, interiorFaces, u, B2E, Bn, In, nelem);
    
    vector<vector<double>> U_firstOrder = u;
    
    vector<vector<int>> elemBounds = genElemBounds(elem, I2E, B2E);
    vector<vector<int>> globalEdges = genGlobalEdge(bounds, interiorFaces);
    
    vector<double> L1_norms;
    L1_norms.reserve(10000);
    
    rk2(opt, U_firstOrder, Area, nodes, elem, Minf, alphaDeg, Bn, In, elemBounds, bounds, interiorFaces, globalEdges, I2E, B2E, "BJ", pow(10,-5), Area, elem.size(), CFL,L1_norms);
    
    ofstream outputFileStates;
    outputFileStates.open("/Users/jakeyeaman/Desktop/Comp Fluid Dyn II/Project-2/outputStates.txt");

    for(int iElem = 0; iElem < elem.size(); iElem++){
        for(int iU = 0; iU < 4; iU++){
            outputFileStates << U_firstOrder[iElem][iU] << " ";
        }
        outputFileStates << "\n";
    }
    outputFileStates.close();
    
    ofstream outputFileNorms;
    outputFileNorms.open("/Users/jakeyeaman/Desktop/Comp Fluid Dyn II/Project-2/outputL1Norms.txt");
    for(int iNorm = 0; iNorm < L1_norms.size(); iNorm++){
        outputFileNorms << L1_norms[iNorm] << "\n";
    }
    outputFileNorms.close();
    
    return 0;

}
