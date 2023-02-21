#include <iostream>
#include <vector>
#include <cmath>
#include <bits/stdc++.h>
#include "elem2Edge.h"
using namespace std;

bool compareInterval(vector<double> i1, vector<double> i2)
{
    return (i1[2] > i2[2]);
}

  void adapt(vector<vector<double>> &splineCoord,vector<vector<int>> &elemBounds,vector<vector<double>> &Bn,vector<vector<double>> &u,int &niedge, int &nnode, int &nbedge, int &nelem, vector<vector<double>> &I2E, vector<vector<double>> &B2E, vector<vector<double>> &elem, vector<vector<double>> &nodes, vector<vector<double>> &bounds, vector<vector<double>> &interiorFaces) {
    vector<double> mach(nelem);
    vector<vector<int>> elembounds=genElemBounds(elem,I2E,B2E);
    for (int i=0;i<nelem;++i) {
        double r=u[i][0];
        double U=u[i][1]/r;
        double v=u[i][2]/r;
        double p=(1.4-1)*(u[i][3]-.5*r*(pow(U,2)+pow(v,2)));
       double c=sqrt(1.4*p/r);
        mach[i]=sqrt(pow(U,2)+pow(v,2))/c;

    }
    vector<vector<double>> interior_error(niedge,vector<double>(3));
    for (int i=0;i<niedge;++i) {
        interior_error[i][0]=1;
        interior_error[i][1]=i+1;
        int node1 = interiorFaces[i][0] - 1;
        int node2 = interiorFaces[i][1] - 1;
        double length = sqrt(pow(nodes[node1][0]-nodes[node2][0],2) + pow(nodes[node1][1]-nodes[node2][1],2));
        int elemL= interiorFaces[i][2] - 1;
        int elemR=interiorFaces[i][3] - 1;
        interior_error[i][2]=abs(mach[elemL]-mach[elemR])*sqrt(length);
    }
 
    vector<vector<double>> boundary_error(nbedge,vector<double>(3));
    for (int i=0;i<nbedge;++i) {
        boundary_error[i][0]=0;
        boundary_error[i][1]=i+1;
        int node1 = bounds[i][0] - 1;
        int node2 = bounds[i][1] - 1;
        double length = sqrt(pow(nodes[node1][0]-nodes[node2][0],2) + pow(nodes[node1][1]-nodes[node2][1],2));
        int elem= B2E[i][0] - 1;
        double r=u[elem][0];
       double U=u[elem][1]/r;
       double v=u[elem][2]/r;
       double p=(1.4-1)*(u[i][3]-.5*r*(pow(U,2)+pow(v,2)));
      double  c=sqrt(1.4*p/r);
      if (B2E[i][2]==4){
         boundary_error[i][2]=0;
      }
      else {
        boundary_error[i][2]=abs(U*Bn[i][0]+v*Bn[i][1])*sqrt(length)/c;
      }
    }
	
	 
	 interior_error.insert( interior_error.end(), boundary_error.begin(), boundary_error.end());
     sort(interior_error.begin(), interior_error.end(), compareInterval);
	 
	 
     int num_edges=niedge+nbedge;
     int num_flag=num_edges*.03;
     vector<vector<double>> flag_elem(nelem,vector<double>(4));
	 
     for (int i=0;i<num_flag;++i) {
        if (interior_error[i][0]) {
            int edge=interior_error[i][1]-1;
            int elemL=I2E[edge][0]-1;
            int faceL=I2E[edge][1];
            int elemR=I2E[edge][2]-1;
            int faceR=I2E[edge][3];
            flag_elem[elemL][faceL]=-1;
            flag_elem[elemR][faceR]=-1;

        }
        else {
            int edge=interior_error[i][1]-1;
            int elemL=B2E[edge][0]-1;
            int faceL=B2E[edge][1];
            flag_elem[elemL][faceL]=edge+1;
        }
     }
	 

     vector<int> globalFlag(num_edges);
     for (int i=0;i<nelem;++i){
        if(flag_elem[i][1]!=0 || flag_elem[i][2]!=0 || flag_elem[i][3]!=0) {
            int edge1 = elemBounds[i][0];
            int edge2 = elemBounds[i][1];
            int edge3 = elemBounds[i][2];

            globalFlag[edge1] = 1;
            globalFlag[edge2] = 1;
            globalFlag[edge3] = 1;
        }
     }
	
    // matrix mapping new node coordinates that are on the Spline, from spline output
	vector<vector<double>> splineNode(nbedge, vector<double>(2));
    // process Spline output into 2D vector that is usable for this function
	double size1 = splineCoord[0].size();
	double size2 = splineCoord[2].size();
	double size3 = splineCoord[4].size();

	// populates first size1 rows with x,y of Group 1
	for (int i = 0; i < size1; i++){
		splineNode[i][0] = splineCoord[0][i];
		splineNode[i][1] = splineCoord[1][i];
	}
	// populates first size1 rows with x,y of Group 2
	for (int i = 0; i < size2; i++){
		splineNode[size1+i][0] = splineCoord[2][i];
		splineNode[size1+i][1] = splineCoord[3][i];
	}
	// populates first size1 rows with x,y of Group 3
	for (int i = 0; i < size3; i++){
		splineNode[size1+size2+i][0] = splineCoord[4][i];
		splineNode[size1+size2+i][1] = splineCoord[5][i];
	}
    vector<vector<double>> midPoints(nelem,vector<double>(3));
    vector<vector<double>> boundMod(nbedge,vector<double>(2));
    int tally = 0;

	
    for (int i=0; i<num_edges;i++){
        if (globalFlag[i] == 1){
            if (i<niedge){
                
                int elemL=I2E[i][0]-1;
                int faceL=I2E[i][1];
                int elemR=I2E[i][2]-1;
                int faceR=I2E[i][3];

                flag_elem[elemL][faceL]=-1;
                flag_elem[elemR][faceR]=-1;

                flag_elem[elemL][0] += 1;
                flag_elem[elemR][0] += 1;

                // create new nodes at midpoints of interior faces, map to elements
                double x1 = nodes[interiorFaces[i][0] - 1][0];
                double x2 = nodes[interiorFaces[i][1] - 1][0];
                double y1 = nodes[interiorFaces[i][0] - 1][1];
                double y2 = nodes[interiorFaces[i][1] - 1][1];

                nodes.push_back(vector<double>{(x1 + x2)/2,(y1 + y2)/2}); // creates new node and adds to nodes matrix
                int newNode = nodes.size(); // gets the node number of newest node

                midPoints[elemL][faceL-1] = newNode; // map midpoint of interior face i to adjacent left element and corresponding local face
                midPoints[elemR][faceR-1] = newNode; // map midpoint of interior face i to adjacent right element and corresponding local face
            }
            else {
                int edge=i-niedge;
                int elemL=B2E[edge][0]-1;
                int faceL=B2E[edge][1];
                flag_elem[elemL][faceL]=edge+1;

                flag_elem[elemL][0] += 1;
                if (B2E[edge][2] != 4){
                    nodes.push_back(vector<double>{splineNode[edge][0],splineNode[edge][1]}); // create new node using coordinates from splines
                }
                else{
                    double x1 = nodes[bounds[i-niedge][0] - 1][0];
                    double x2 = nodes[bounds[i-niedge][1] - 1][0];
                    double y1 = nodes[bounds[i-niedge][0] - 1][1];
                    double y2 = nodes[bounds[i-niedge][1] - 1][1];

                    nodes.push_back(vector<double>{(x1 + x2)/2,(y1 + y2)/2}); // creates new node and adds to nodes matrix
                }
                int newNode = nodes.size(); // gets node number of newest node

                boundMod[edge][0] = 1;
                boundMod[edge][1] = newNode;

			    midPoints[elemL][faceL-1] = newNode; // map midpoint of boundary face i to adjacent element and corresponding local face
            }
            
        }
    }
	
	
	for(int j = 0; j<boundMod.size();j++){
		
                if(boundMod[j][0]){ // && bounds[i][2]!=4
				
                    double hold = bounds[j+tally][1];
                    bounds[j+tally][1]=boundMod[j][1];
                    bounds.insert(bounds.begin()+j+tally+1,vector<double>{boundMod[j][1],hold,bounds[j+tally][2]});
                    tally++;
                }
				
	}
	
    // Loop over elements, determine if element has flagged edges, how many flagged edges, and which edges are flagged,
	// perform local refinement correspondingly
	for (int i = 0; i < nelem; i++){
		// store current nodes for element i in placeholder variables
		double place1 = elem[i][0];
		double place2 = elem[i][1];
		double place3 = elem[i][2];
        double node1;
        double node2;
        double node3;
		vector<double> newNodeCoord(2);

		// FOr 0 flagged edges
		if (flag_elem[i][0] == 0){
			elem[i][0] = place1;
			elem[i][1] = place2;
			elem[i][2] = place3;
		}
		// For 1 flagged edge: determine which edge is flagged, nodes 1 + 2 are on edge, crossNode is across from edge
		else if (flag_elem[i][0] == 1){
			double crossNode;
			double newNode;
			//I2E_size++; // nothing to do with I2E, just adds 1 to niedge since 1 new interior edge is created
	
			// following conditional statements determine naming of nodes based on which local face is flagged
			if (flag_elem[i][1] != 0){
				node1 = place2;
				node2 = place3;
				crossNode = place1;
				newNode = midPoints[i][0];
			}
			else if (flag_elem[i][2] != 0){
				node1 = place1;
				node2 = place3;
				crossNode = place2;
				newNode = midPoints[i][1]; 
			}
			else if (flag_elem[i][3] != 0){
				node1 = place2;
				node2 = place1;
				crossNode = place3;
				newNode = midPoints[i][2]; 
			}
			
			// create 1st new element
			elem.push_back(vector<double>{newNode,crossNode,node2});

			// replace current element i with other new element
			elem[i][0] = newNode;
			elem[i][1] = node1;
			elem[i][2] = crossNode;

            u.push_back(vector<double>{u[i][0],u[i][1],u[i][2],u[i][3]});
		}
		// For 2 flagged edges: determine which edges are flagged, nodes 1 + 2 are on edge, crossNode is across from edge
		else if (flag_elem[i][0] == 2){
			vector<double> edge1(2);
			vector<double> edge2(2);
			vector<double> edge3(2);
			vector<double> edge4(2);

			double mag1;
			double mag2;
			double mag3;
			double mag4;

			double nodeShare;
			double newNode1;
			double newNode2;
			
			// determine which 2 local faces are flagged, determine naming of nodes
			if ((flag_elem[i][1] != 0) & (flag_elem[i][2] != 0)){
				node1 = place1;
				node2 = place2;
				nodeShare = place3;
				newNode1 = midPoints[i][0];
				newNode2 = midPoints[i][1];
			}
			else if ((flag_elem[i][1] != 0) & (flag_elem[i][3] != 0)){
				node1 = place3;
				node2 = place1;
				nodeShare = place2;
				newNode1 = midPoints[i][2];
				newNode2 = midPoints[i][0];
			}
			else if ((flag_elem[i][2] != 0) & (flag_elem[i][3] != 0)){
				node1 = place2;
				node2 = place3;
				nodeShare = place1;
				newNode1 = midPoints[i][1];
				newNode2 = midPoints[i][2];
			}

			// compute edge vectors
			for (int j = 0; j < 2; j++){
				edge1[j] = nodes[nodeShare-1][j] - nodes[node1-1][j];
				edge3[j] = nodes[node2-1][j] - nodes[node1-1][j];

				edge2[j] = nodes[nodeShare-1][j] - nodes[node2-1][j];
				edge4[j] = nodes[node1-1][j] - nodes[node2-1][j];
			}
			// magnitude/length of each edge
			mag1 = sqrt(pow(edge1[0],2)+pow(edge1[1],2));
			mag2 = sqrt(pow(edge2[0],2)+pow(edge2[1],2));
			mag3 = sqrt(pow(edge3[0],2)+pow(edge3[1],2));
			mag4 = sqrt(pow(edge4[0],2)+pow(edge4[1],2));

			// cos/dot product theory to determine angles between edges
			double theta1 = acos(((edge1[0]*edge3[0]) + (edge1[1]*edge3[1]))/(mag1*mag3));
			double theta2 = acos(((edge2[0]*edge4[0]) + (edge2[1]*edge4[1]))/(mag2*mag4));
			
			// determine how element will be refined based on the angles between edges
			if (theta1 >= theta2){
				// create 1st new element
				elem.push_back(vector<double>{node1,newNode1,newNode2});

				// create 2nd new element
				elem.push_back(vector<double>{node1,node2,newNode1});

				// replace current element i with 3rd new element
				elem[i][0] = newNode1;
				elem[i][1] = nodeShare;
				elem[i][2] = newNode2;
			}
			else{
				// create 1st new element
				elem.push_back(vector<double>{node2,newNode1,newNode2});

				// create 2nd new element
				elem.push_back(vector<double>{node1,node2,newNode2});

				// replace current element i with 3rd new element
				elem[i][0] = newNode1;
				elem[i][1] = nodeShare;
				elem[i][2] = newNode2;
			}

            u.push_back(vector<double>{u[i][0],u[i][1],u[i][2],u[i][3]});
            u.push_back(vector<double>{u[i][0],u[i][1],u[i][2],u[i][3]});
		}
		// For 3 flagged edges
		else if (flag_elem[i][0] == 3){
			// store current nodes for element i in placeholder variables
			double place1 = elem[i][0];
			double place2 = elem[i][1];
			double place3 = elem[i][2];
			vector<double> newNodeCoord(2);

			double newNode1 = midPoints[i][0];
			double newNode2 = midPoints[i][1];;
			double newNode3 = midPoints[i][2];;

			node1 = place1;
			node2 = place2;
			node3 = place3;

			// create 1st new element
			elem.push_back(vector<double>{node1,newNode3,newNode2});

			// create 2nd new element
			elem.push_back(vector<double>{node2,newNode1,newNode3});

			// create 3rd new element
			elem.push_back(vector<double>{node3,newNode2,newNode1});

			// replace current element i with other new element
			elem[i][0] = newNode1;
			elem[i][1] = newNode2;
			elem[i][2] = newNode3;

            // set states
            u.push_back(vector<double>{u[i][0],u[i][1],u[i][2],u[i][3]});
            u.push_back(vector<double>{u[i][0],u[i][1],u[i][2],u[i][3]});
            u.push_back(vector<double>{u[i][0],u[i][1],u[i][2],u[i][3]});
		}
	}
	// find new number of nodes, elements, and boundary edges
	nnode = nodes.size();
	nelem = elem.size();
	nbedge = bounds.size();
	niedge = .5*(3*nelem - nbedge);
	
}