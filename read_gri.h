#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

vector<vector<double>> getnodes(string fileName){

    // Open the input file
    ifstream inputFile(fileName);

    // Read the first line
    int numnodes;
    inputFile >> numnodes; // nodes
        int numelem;
    inputFile >> numelem; // elemnts
        int b4;
    inputFile >> b4; // farfield boundary
        int b1;
    inputFile >> b1; // slat boundary
        int b2;
    inputFile >> b2; // main boundary
        int b3;
    inputFile >> b3; // flap boundary    

    // Create a vector of vectors to store nodes
    vector<vector<double>> nodes(numnodes, vector<double>(2));

    // Read in nodes from the file
    for (int i = 0; i < numnodes; i++) {
        for (int j = 0; j < 2; j++) {
            inputFile >> nodes[i][j];
        }
    }

    // Close the input file
    inputFile.close();

return nodes;

}



vector<vector<double>> getelement(string fileName){

    // Open the input file
    ifstream inputFile(fileName);

    // Read the first line
    int numnodes;
    inputFile >> numnodes; // nodes
        int numelem;
    inputFile >> numelem; // elemnts
        int b4;
    inputFile >> b4; // farfield boundary
        int b1;
    inputFile >> b1; // slat boundary
        int b2;
    inputFile >> b2; // main boundary
        int b3;
    inputFile >> b3; // flap boundary    

    int skip=numnodes*2+b1*2+b2*2+b3*2+b4*2;

    
    // skip the values that I dont want
     vector<double> skipvec(skip);

    for (int i = 0; i < skip; i++) {
     inputFile >> skipvec[i];
    
    }   

    // Create a vector of vectors to store elem
    vector<vector<double>> elem(numelem, vector<double>(3));

    // Read in elem from the file
    for (int i = 0; i < numelem; i++) {
        for (int j = 0; j < 3; j++) {
            inputFile >> elem[i][j];
        }
    }

    // Close the input file
    inputFile.close();

return elem;

}

vector<vector<double>> getboundaries(string fileName){

    // Open the input file
    ifstream inputFile(fileName);

    // Read the first line
        int numnodes;
    inputFile >> numnodes; // nodes
        int numelem;
    inputFile >> numelem; // elemnts
        int b4;
    inputFile >> b4; // farfield boundary
        int b1;
    inputFile >> b1; // slat boundary
        int b2;
    inputFile >> b2; // main boundary
        int b3;
    inputFile >> b3; // flap boundary    

    int skip=numnodes*2;

   
    // skip the values that I dont want
     vector<double> skipvec(skip);

    for (int i = 0; i < skip; i++) {
     inputFile >> skipvec[i];
    
    }   

    // Create a vector of vectors to store bounds
    vector<vector<double>> bounds((b1+b2+b3+b4), vector<double>(3));

    // farfield
    for (int i = 0; i < b4; i++) {
        for (int j = 0; j < 2; j++) {
            inputFile >> bounds[i][j];
        }
        bounds[i][2]=4;
    }
        // slat
    for (int i = 0; i < b1; i++) {
        for (int j = 0; j < 2; j++) {
            inputFile >> bounds[i+b4][j];
        }
        bounds[i+b4][2]=1;
    }
        // main
    for (int i = 0; i < b2; i++) {
        for (int j = 0; j < 2; j++) {
            inputFile >> bounds[i+b4+b1][j];
        }
        bounds[i+b4+b1][2]=2;
    }
            // flap
    for (int i = 0; i < b3; i++) {
        for (int j = 0; j < 2; j++) {
            inputFile >> bounds[i+b4+b1+b2][j];
        }
        bounds[i+b4+b1+b2][2]=3;
    }
    // Close the input file
    inputFile.close();

return bounds;

}

int getnnode(string fileName){

    // Open the input file
    ifstream inputFile(fileName);

    // Read the first line
    int nnode;
    inputFile >> nnode; // nodes
        int numelem;
    inputFile >> numelem; // elemnts
        int b4;
    inputFile >> b4; // farfield boundary
        int b1;
    inputFile >> b1; // slat boundary
        int b2;
    inputFile >> b2; // main boundary
        int b3;
    inputFile >> b3; // flap boundary    

return nnode;

}

int getnelem(string fileName){

    // Open the input file
    ifstream inputFile(fileName);

    // Read the first line
    int nnode;
    inputFile >> nnode; // nodes
        int nelem;
    inputFile >> nelem; // elemnts
        int b4;
    inputFile >> b4; // farfield boundary
        int b1;
    inputFile >> b1; // slat boundary
        int b2;
    inputFile >> b2; // main boundary
        int b3;
    inputFile >> b3; // flap boundary    

return nelem;

}

int getnbedge(string fileName){

    // Open the input file
    ifstream inputFile(fileName);

    // Read the first line
    int nnode;
    inputFile >> nnode; // nodes
        int numelem;
    inputFile >> numelem; // elemnts
        int b4;
    inputFile >> b4; // farfield boundary
        int b1;
    inputFile >> b1; // slat boundary
        int b2;
    inputFile >> b2; // main boundary
        int b3;
    inputFile >> b3; // flap boundary    

    int nbedge = b1+b2+b3+b4;

return nbedge;

}