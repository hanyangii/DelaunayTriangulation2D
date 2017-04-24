#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <list>
#include <algorithm>

#include "Delaunay.h"

using namespace std;


list<Triangle> MergeTriangles(list<Triangle> LList, list<Triangle> RList){
	list<Triangle> merged_list;
	
	return merged_list;
}

list<Triangle> MergeSets(list<Point> PointList){
	
	if(PointList.size() > 6){ // divide Points to two subsets

	}
	else{ // Make LR edge for Delaunay Triangle
		Point[3] Lset;
		Point[3] Rset;
		for(int i=0; i<3;i++){
			Lset[0] = PointList.pop_front();
		}
		Traingle LTriangle = MakeTriangle(Lset[0],Lset[1],Lset[2]);
		int i=0;
		while(!PointList.empty()){
			Rset[i++]=PointList.pop_front();
		}


		Triangle RTriangle;
		if(i>2){//Rset has three points -> make triangle
			RTriangle = MakeTriangle(Rset[0],Rset[1],Rset[2]);
			Rminy = MinYpoint(RTriangle);
			Lminy = MinYpoint(LTriangle);
			Edge base_LR = linkEdge(Rminy, Lminy);
			
			//Select potential LR
			Edge potential_LR = PotentialLR(base_LR, RTriangle.
		}
		else{


		}

	}// endif
}

list<Triangle> Delaunay(list<Point> PointList){
	PointList = PointList.sort(compare_y());
	PointList = PointList.sort(compare_x());
	int idx=0;
	for(int i =PointList.begin(); i!= PointList.end(); ++i){
			i->num = idx++;
	}

	return MergeSets(PointList);
}


list<Point> InputPoints(){
	//Input points
	ifstream PointFile("Points.txt");
	list<Point> PointList;
	
	string line;
	while(!PointFile.eof()){
		getline(PointFile, line);
		stringstream linestream(line);
		string word;
		int idx = 0;
		Point new_pt;
		while(getline(linestream, word, ' ')){
			if(idx%2==0) new_pt.x= strtof((word).c_str(),0);
			else new_pt.y= strtof((word).c_str(),0);
			if(++idx ==2) break;
		} 
		PointList.push_back(new_pt);
	}
	PointFile.close();
	
	return PointList;
	
}

int main(){
	
	list<Point> PointList=InputPoints();
	list<Traingle> result_mesh = Delaunay(PointList);

	return 0;
}
