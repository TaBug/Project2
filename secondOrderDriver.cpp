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
    double gamma = 1.4;
    double CFL = 1;
    double Minf = 0.25;
    double alphaDeg = 8;
    vector<double> uinf = computeFreestreamState(Minf, alphaDeg);
 
    cout << "Computational Fluid Dynamics II Driver...\n";
    int nnode, nbedge, nelem;
    string fileName = "C:/Users/david/OneDrive - Umich/UM/Courses/AE 623/Projects/Project 2/Project2/CoarseMesh.msh";
    getMatSizes(fileName, nnode, nbedge, nelem); // finds the number of nodes, boundary edges, and elements
        
    int niedge = 0.5 * (3 * nelem - nbedge); // number of interior edges
    
    vector<vector<double>> nodes(nnode, vector<double>(2)); // node coordinates matrix (#nodes x 2)
    vector<vector<double>> elem(nelem, vector<double>(3)); // [node1, node2, node3]
    vector<vector<double>> bounds(nbedge, vector<double>(3)); // [node1, node2, boundaryGroup]
    readGmshFile(fileName,nnode,nbedge,nelem,nodes,elem,bounds); // construct the matrices
    generateGri("C:/Users/david/OneDrive - Umich/UM/Courses/AE 623/Projects/Project 2/Project2/CoarseMesh.txt", nnode, nbedge, nelem, nodes, elem, bounds); // convert mesh to .gri file
    vector<vector<double>> interiorFaces = genInteriorFaceVec(niedge, nelem, nbedge, bounds, elem);

    // GENERATING TASK 2 DATA STRUCTURES
    vector<vector<double>> I2E = genI2E(niedge, nelem, interiorFaces, elem);
    vector<vector<double>> B2E = genB2E(nbedge, nelem, bounds, elem);
    vector<vector<double>> In = genIn(niedge, interiorFaces, elem, nodes, I2E);
    vector<vector<double>> Bn = genBn(nbedge, bounds, nodes, elem, B2E);
    vector<double> Area = genArea(nelem, elem, nodes);

    // initialize the state to free-stream condition
    vector<vector<double>> u(nelem, uinf);

    // 1st-order finite volume
    // update the state to the converged state using 1st-order FVM
    int opt;
    FVM_1st(bounds, nodes, interiorFaces, u, B2E, Bn, In, nelem, opt, uinf);
    
    vector<vector<double>> U_firstOrder = u;
    vector<vector<int>> elemBounds = genElemBounds(elem, I2E, B2E);
    vector<vector<int>> globalEdges = genGlobalEdge(bounds, interiorFaces);
    
    rk2(opt, U_firstOrder, Area, nodes, elem, Minf, alphaDeg, Bn, In, elemBounds, bounds, interiorFaces, globalEdges, I2E, B2E, "BJ", pow(10,-2), Area, elem.size(), CFL);
    
//    ofstream outputFile;
//    outputFile.open("/Users/jakeyeaman/Desktop/Comp Fluid Dyn II/Project-2/output.txt");
//
//    for(int iElem = 0; iElem < u.size(); iElem++){
//        for(int iU = 0; iU < 1; iU++){
//            outputFile << u[iElem][iU] << " ";
//        }
//        outputFile << "\n";
//    }
//    outputFile.close();
    return 0;
}
