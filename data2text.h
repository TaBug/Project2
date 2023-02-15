//
//  data2text.h
//  p2-AEROSP-623
//
//  Created by Jake Yeaman on 2/14/23.
//

#ifndef data2text_h
#define data2text_h

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>


void state2text(vector<vector<double>> &U, string filename){
    
    ofstream outputFile;
    outputFile.open(filename);
    
    for(int iElem = 0; iElem < U.size(); iElem++){
        for(int iU = 0; iU < 4; iU++){
            outputFile << U[iElem][iU] << " ";
        }
        outputFile << "\n";
    }
    
    return;
    
}

void L1Residual2text(vector<double> &L1Residual, string filename){
    ofstream outputFile;
    outputFile.open(filename);
    for(int iR = 0; iR < L1Residual.size(); iR++){
        outputFile << L1Residual[iR] << "\n";
    }
    return;
}

#endif /* data2text_h */
