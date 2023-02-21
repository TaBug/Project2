#include <iostream>
#include <vector>
#include <cmath>
#include "FE_FVM.h"
#include <iostream>
#include "processMesh (2).h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_set>
#include <numeric>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include "adaptation.h"
#include "solver.h"
#include "RK2_FVM.h"
#include "spline.h"
#include "elem2Edge.h"
using namespace std;

int main () {
    
        double gamma = 1.4;
    int opt = 2;
    double CFL = .1;
    double Minf = 0.5;
    double alphaDeg = 8;
 cout << "Computational Fluid Dynamics II Driver...\n";
 int nnode, nbedge, nelem;
 string fileName = "/mnt/c/Users/sasha/Documents/MATLAB/aero623/.vscode/Project2-main/CoarseMesh (1)";
getMatSizes(fileName, nnode,nbedge,nelem);
    
    int niedge = 0.5*(3*nelem - nbedge);
    
    vector<vector<double>> nodes(nnode, vector<double>(2));
    vector<vector<double>> elem(nelem, vector<double>(3)); // [node1, node2, node3]
    vector<vector<double>> bounds(nbedge, vector<double>(3)); // [node1, node2, boundaryGroup]
    
    readGmshFile(fileName,nnode,nbedge,nelem,nodes,elem,bounds);
generateGri("/mnt/c/Users/sasha/Documents/MATLAB/aero623/.vscode/Project2-main/CoarseMesh.gri", nnode, nbedge, nelem, nodes, elem, bounds);
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
    u[i][1]=.5*cos(8*3.14159/180);
    u[i][2]=.5*sin(8*3.14159/180);
    u[i][3]=1/(.4*1.4)+(.5*.5)/2;
}

FVM_1st(bounds,nodes,interiorFaces,u,B2E,Bn,In,nelem,Area);

cout<<"first order"<<endl;

    
    
   vector<vector<int>> elemBounds = genElemBounds(elem, I2E, B2E);
    vector<vector<int>> globalEdges = genGlobalEdge(bounds, interiorFaces);
    
    vector<double> L1_norms;
    L1_norms.reserve(10000);
    int niter=0;
    rk2(opt, u, Area, nodes, elem, Minf, alphaDeg, Bn, In, elemBounds, bounds, interiorFaces, globalEdges, I2E, B2E, "BJ", pow(10,-5), Area, elem.size(), CFL,L1_norms,niter);
   ofstream outputFileStates;
    outputFileStates.open("/mnt/c/Users/sasha/Documents/MATLAB/aero623/.vscode/Project2-main/outputs.txt");
    for(int iElem = 0; iElem < elem.size(); iElem++){
        for(int iU = 0; iU < 4; iU++){
            outputFileStates << u[iElem][iU] << " ";
            
        }
        outputFileStates << "\n";
    }
  
    outputFileStates.close();
    outputFileStates.open("/mnt/c/Users/sasha/Documents/MATLAB/aero623/.vscode/Project2-main/elem.txt");
    for(int iElem = 0; iElem < elem.size(); iElem++){
        for(int iU = 0; iU < 3; iU++){
            outputFileStates << elem[iElem][iU] << " ";
        }
        outputFileStates << "\n";
    }
    outputFileStates.close();
     
    outputFileStates.open("/mnt/c/Users/sasha/Documents/MATLAB/aero623/.vscode/Project2-main/nodes.txt");
    for(int iElem = 0; iElem < nodes.size(); iElem++){
        for(int iU = 0; iU < 2; iU++){
            outputFileStates << nodes[iElem][iU] << " ";
        }
        outputFileStates << "\n";
    }
    outputFileStates.close();

  
   vector<vector<double>> testSpline = spline(nodes,bounds,B2E);
   
  adapt(testSpline,elemBounds,Bn,u,niedge, nnode, nbedge, nelem, I2E, B2E, elem, nodes, bounds,interiorFaces);
  
    interiorFaces = genInteriorFaceVec(niedge, nelem, nbedge, bounds, elem);
     
    I2E = genI2E(niedge, nelem, interiorFaces, elem);
    B2E = genB2E(nbedge, nelem, bounds, elem);
    In = genIn(niedge, interiorFaces, elem, nodes, I2E);
    Bn = genBn(nbedge, bounds, nodes, elem, B2E);
    Area = genArea(nelem, elem, nodes);

    FVM_1st(bounds,nodes,interiorFaces,u,B2E,Bn,In,nelem,Area);

    elemBounds = genElemBounds(elem, I2E, B2E);
    globalEdges = genGlobalEdge(bounds, interiorFaces);
    L1_norms.reserve(10000);
    niter=0;
    rk2(opt, u, Area, nodes, elem, Minf, alphaDeg, Bn, In, elemBounds, bounds, interiorFaces, globalEdges, I2E, B2E, "BJ", pow(10,-5), Area, elem.size(), CFL,L1_norms,niter);
generateGri("/mnt/c/Users/sasha/Documents/MATLAB/aero623/.vscode/Project2-main/adapt1.gri", nnode, nbedge, nelem, nodes, elem, bounds);
    outputFileStates.open("/mnt/c/Users/sasha/Documents/MATLAB/aero623/.vscode/Project2-main/outputs1.txt");
    for(int iElem = 0; iElem < elem.size(); iElem++){
        for(int iU = 0; iU < 4; iU++){
            outputFileStates << u[iElem][iU] << " ";
            
        }
        outputFileStates << "\n";
    }
    outputFileStates.close();

    outputFileStates.open("/mnt/c/Users/sasha/Documents/MATLAB/aero623/.vscode/Project2-main/elem1.txt");
    for(int iElem = 0; iElem < elem.size(); iElem++){
        for(int iU = 0; iU < 3; iU++){
            outputFileStates << elem[iElem][iU] << " ";
        }
        outputFileStates << "\n";
    }
    outputFileStates.close();
     
    outputFileStates.open("/mnt/c/Users/sasha/Documents/MATLAB/aero623/.vscode/Project2-main/nodes1.txt");
    for(int iElem = 0; iElem < nodes.size(); iElem++){
        for(int iU = 0; iU < 2; iU++){
            outputFileStates << nodes[iElem][iU] << " ";
        }
        outputFileStates << "\n";
    }
    outputFileStates.close();


    testSpline = spline(nodes,bounds,B2E);
  adapt(testSpline,elemBounds,Bn,u,niedge, nnode, nbedge, nelem, I2E, B2E, elem, nodes, bounds,interiorFaces);
    interiorFaces = genInteriorFaceVec(niedge, nelem, nbedge, bounds, elem);
    
    I2E = genI2E(niedge, nelem, interiorFaces, elem);
    B2E = genB2E(nbedge, nelem, bounds, elem);
    In = genIn(niedge, interiorFaces, elem, nodes, I2E);
    Bn = genBn(nbedge, bounds, nodes, elem, B2E);
    Area = genArea(nelem, elem, nodes);

    FVM_1st(bounds,nodes,interiorFaces,u,B2E,Bn,In,nelem,Area);
    
    elemBounds = genElemBounds(elem, I2E, B2E);
    globalEdges = genGlobalEdge(bounds, interiorFaces);
    L1_norms.reserve(10000);
    niter=0;
    rk2(opt, u, Area, nodes, elem, Minf, alphaDeg, Bn, In, elemBounds, bounds, interiorFaces, globalEdges, I2E, B2E, "BJ", pow(10,-5), Area, elem.size(), CFL,L1_norms,niter);
generateGri("/mnt/c/Users/sasha/Documents/MATLAB/aero623/.vscode/Project2-main/adapt2.gri", nnode, nbedge, nelem, nodes, elem, bounds);
    outputFileStates.open("/mnt/c/Users/sasha/Documents/MATLAB/aero623/.vscode/Project2-main/outputs2.txt");
    for(int iElem = 0; iElem < elem.size(); iElem++){
        for(int iU = 0; iU < 4; iU++){
            outputFileStates << u[iElem][iU] << " ";
            
        }
        outputFileStates << "\n";
    }
    outputFileStates.close();

    outputFileStates.open("/mnt/c/Users/sasha/Documents/MATLAB/aero623/.vscode/Project2-main/elem2.txt");
    for(int iElem = 0; iElem < elem.size(); iElem++){
        for(int iU = 0; iU < 3; iU++){
            outputFileStates << elem[iElem][iU] << " ";
        }
        outputFileStates << "\n";
    }
    outputFileStates.close();
     
    outputFileStates.open("/mnt/c/Users/sasha/Documents/MATLAB/aero623/.vscode/Project2-main/nodes2.txt");
    for(int iElem = 0; iElem < nodes.size(); iElem++){
        for(int iU = 0; iU < 2; iU++){
            outputFileStates << nodes[iElem][iU] << " ";
        }
        outputFileStates << "\n";
    }
    outputFileStates.close();

testSpline = spline(nodes,bounds,B2E);
  adapt(testSpline,elemBounds,Bn,u,niedge, nnode, nbedge, nelem, I2E, B2E, elem, nodes, bounds,interiorFaces);
    interiorFaces = genInteriorFaceVec(niedge, nelem, nbedge, bounds, elem);
    
    I2E = genI2E(niedge, nelem, interiorFaces, elem);
    B2E = genB2E(nbedge, nelem, bounds, elem);
    In = genIn(niedge, interiorFaces, elem, nodes, I2E);
    Bn = genBn(nbedge, bounds, nodes, elem, B2E);
    Area = genArea(nelem, elem, nodes);

    FVM_1st(bounds,nodes,interiorFaces,u,B2E,Bn,In,nelem,Area);
    
    elemBounds = genElemBounds(elem, I2E, B2E);
    globalEdges = genGlobalEdge(bounds, interiorFaces);
    L1_norms.reserve(10000);
    niter=0;
    rk2(opt, u, Area, nodes, elem, Minf, alphaDeg, Bn, In, elemBounds, bounds, interiorFaces, globalEdges, I2E, B2E, "BJ", pow(10,-5), Area, elem.size(), CFL,L1_norms,niter);
generateGri("/mnt/c/Users/sasha/Documents/MATLAB/aero623/.vscode/Project2-main/adapt3.gri", nnode, nbedge, nelem, nodes, elem, bounds);
    outputFileStates.open("/mnt/c/Users/sasha/Documents/MATLAB/aero623/.vscode/Project2-main/outputs3.txt");
    for(int iElem = 0; iElem < elem.size(); iElem++){
        for(int iU = 0; iU < 4; iU++){
            outputFileStates << u[iElem][iU] << " ";
            
        }
        outputFileStates << "\n";
    }
    outputFileStates.close();

    outputFileStates.open("/mnt/c/Users/sasha/Documents/MATLAB/aero623/.vscode/Project2-main/elem3.txt");
    for(int iElem = 0; iElem < elem.size(); iElem++){
        for(int iU = 0; iU < 3; iU++){
            outputFileStates << elem[iElem][iU] << " ";
        }
        outputFileStates << "\n";
    }
    outputFileStates.close();
     
    outputFileStates.open("/mnt/c/Users/sasha/Documents/MATLAB/aero623/.vscode/Project2-main/nodes3.txt");
    for(int iElem = 0; iElem < nodes.size(); iElem++){
        for(int iU = 0; iU < 2; iU++){
            outputFileStates << nodes[iElem][iU] << " ";
        }
        outputFileStates << "\n";
    }
    outputFileStates.close();


testSpline = spline(nodes,bounds,B2E);
  adapt(testSpline,elemBounds,Bn,u,niedge, nnode, nbedge, nelem, I2E, B2E, elem, nodes, bounds,interiorFaces);
    interiorFaces = genInteriorFaceVec(niedge, nelem, nbedge, bounds, elem);
    
    I2E = genI2E(niedge, nelem, interiorFaces, elem);
    B2E = genB2E(nbedge, nelem, bounds, elem);
    In = genIn(niedge, interiorFaces, elem, nodes, I2E);
    Bn = genBn(nbedge, bounds, nodes, elem, B2E);
    Area = genArea(nelem, elem, nodes);

    FVM_1st(bounds,nodes,interiorFaces,u,B2E,Bn,In,nelem,Area);
    
    elemBounds = genElemBounds(elem, I2E, B2E);
    globalEdges = genGlobalEdge(bounds, interiorFaces);
    L1_norms.reserve(10000);
    niter=0;
    rk2(opt, u, Area, nodes, elem, Minf, alphaDeg, Bn, In, elemBounds, bounds, interiorFaces, globalEdges, I2E, B2E, "BJ", pow(10,-5), Area, elem.size(), CFL,L1_norms,niter);
generateGri("/mnt/c/Users/sasha/Documents/MATLAB/aero623/.vscode/Project2-main/adapt4.gri", nnode, nbedge, nelem, nodes, elem, bounds);
    outputFileStates.open("/mnt/c/Users/sasha/Documents/MATLAB/aero623/.vscode/Project2-main/outputs4.txt");
    for(int iElem = 0; iElem < elem.size(); iElem++){
        for(int iU = 0; iU < 4; iU++){
            outputFileStates << u[iElem][iU] << " ";
            
        }
        outputFileStates << "\n";
    }
    outputFileStates.close();

    outputFileStates.open("/mnt/c/Users/sasha/Documents/MATLAB/aero623/.vscode/Project2-main/elem4.txt");
    for(int iElem = 0; iElem < elem.size(); iElem++){
        for(int iU = 0; iU < 3; iU++){
            outputFileStates << elem[iElem][iU] << " ";
        }
        outputFileStates << "\n";
    }
    outputFileStates.close();
     
    outputFileStates.open("/mnt/c/Users/sasha/Documents/MATLAB/aero623/.vscode/Project2-main/nodes4.txt");
    for(int iElem = 0; iElem < nodes.size(); iElem++){
        for(int iU = 0; iU < 2; iU++){
            outputFileStates << nodes[iElem][iU] << " ";
        }
        outputFileStates << "\n";
    }
    outputFileStates.close();

testSpline = spline(nodes,bounds,B2E);
  adapt(testSpline,elemBounds,Bn,u,niedge, nnode, nbedge, nelem, I2E, B2E, elem, nodes, bounds,interiorFaces);
    interiorFaces = genInteriorFaceVec(niedge, nelem, nbedge, bounds, elem);
    
    I2E = genI2E(niedge, nelem, interiorFaces, elem);
    B2E = genB2E(nbedge, nelem, bounds, elem);
    In = genIn(niedge, interiorFaces, elem, nodes, I2E);
    Bn = genBn(nbedge, bounds, nodes, elem, B2E);
    Area = genArea(nelem, elem, nodes);

    FVM_1st(bounds,nodes,interiorFaces,u,B2E,Bn,In,nelem,Area);

    elemBounds = genElemBounds(elem, I2E, B2E);
    globalEdges = genGlobalEdge(bounds, interiorFaces);
    L1_norms.reserve(10000);
    niter=0;
    rk2(opt, u, Area, nodes, elem, Minf, alphaDeg, Bn, In, elemBounds, bounds, interiorFaces, globalEdges, I2E, B2E, "BJ", pow(10,-5), Area, elem.size(), CFL,L1_norms,niter);
generateGri("/mnt/c/Users/sasha/Documents/MATLAB/aero623/.vscode/Project2-main/adapt5.gri", nnode, nbedge, nelem, nodes, elem, bounds);
    outputFileStates.open("/mnt/c/Users/sasha/Documents/MATLAB/aero623/.vscode/Project2-main/outputs5.txt");
    for(int iElem = 0; iElem < elem.size(); iElem++){
        for(int iU = 0; iU < 4; iU++){
            outputFileStates << u[iElem][iU] << " ";
            
        }
        outputFileStates << "\n";
    }
    outputFileStates.close();

    outputFileStates.open("/mnt/c/Users/sasha/Documents/MATLAB/aero623/.vscode/Project2-main/elem5.txt");
    for(int iElem = 0; iElem < elem.size(); iElem++){
        for(int iU = 0; iU < 3; iU++){
            outputFileStates << elem[iElem][iU] << " ";
        }
        outputFileStates << "\n";
    }
    outputFileStates.close();
     
    outputFileStates.open("/mnt/c/Users/sasha/Documents/MATLAB/aero623/.vscode/Project2-main/nodes5.txt");
    for(int iElem = 0; iElem < nodes.size(); iElem++){
        for(int iU = 0; iU < 2; iU++){
            outputFileStates << nodes[iElem][iU] << " ";
        }
        outputFileStates << "\n";
    }
    outputFileStates.close();





    ofstream outputFileNorms;
    outputFileNorms.open("/mnt/c/Users/sasha/Documents/MATLAB/aero623/.vscode/Project2-main/outputL1Norms.txt");
    for(int iNorm = 0; iNorm < L1_norms.size(); iNorm++){
        outputFileNorms << L1_norms[iNorm] << "\n";
    }
    outputFileNorms.close();
    
    return 0;
   
}