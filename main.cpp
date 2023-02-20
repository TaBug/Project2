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
#include "read_gri.h"



using namespace std;

int main(){
    
    // TODO: FILL OUT (also make sure Minf if picked correctly in FE_FVM.h as well)
    int opt = 2;
    string limiterType = "BJ";
    double CFL = 0.5;
    double Minf = 0.25;
    double alphaDeg = 8;
    double convergenceThreshold = pow(10,-5);
    // TODO: end
 
    cout << "Computational Fluid Dynamics II Driver...\n";
//    int nnode, nbedge, nelem;
    string fileName = "/Users/jakeyeaman/Desktop/Comp Fluid Dyn II/Project-2/c1.gri";
//    getMatSizes(fileName, nnode,nbedge,nelem);
//
//    int niedge = 0.5*(3*nelem - nbedge);
//
//    vector<vector<double>> nodes(nnode, vector<double>(2));
//    vector<vector<double>> elem(nelem, vector<double>(3)); // [node1, node2, node3]
//    vector<vector<double>> bounds(nbedge, vector<double>(3)); // [node1, node2, boundaryGroup]
    
    int nnode = getnnode(fileName);
    int nbedge = getnbedge(fileName);
    int nelem = getnelem(fileName);
    int niedge = 0.5 * (3 * nelem - nbedge); // number of interior edges
    
    vector<vector<double>> nodes=getnodes(fileName);
    vector<vector<double>> elem=getelement(fileName);
    vector<vector<double>> bounds=getboundaries(fileName);
    
    // readGmshFile(fileName,nnode,nbedge,nelem,nodes,elem,bounds);
    // generateGri("/Users/jakeyeaman/Desktop/Comp Fluid Dyn II/Project-2/CoarseMesh.gri", nnode, nbedge, nelem, nodes, elem, bounds);
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
        u[i][1]=Minf*cos(8*3.14159/180);
        u[i][2]=Minf*sin(8*3.14159/180);
        u[i][3]=1/(.4*1.4)+Minf*Minf/2;
    }

    FVM_1st(bounds, nodes, interiorFaces, u, B2E, Bn, In, nelem,Minf);
    
    vector<vector<double>> U_firstOrder = u;
    
    vector<vector<int>> elemBounds = genElemBounds(elem, I2E, B2E);
    vector<vector<int>> globalEdges = genGlobalEdge(bounds, interiorFaces);
    
    vector<double> L1_norms;
    L1_norms.reserve(10000);
    
    int niter = 0;
    rk2(opt, U_firstOrder, Area, nodes, elem, Minf, alphaDeg, Bn, In, elemBounds, bounds, interiorFaces, globalEdges, I2E, B2E, limiterType, convergenceThreshold, Area, elem.size(), CFL,L1_norms, niter);
    
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
    
    cout << "Total Number of Iterations: " << niter << "\n";
    
    return 0;

}
