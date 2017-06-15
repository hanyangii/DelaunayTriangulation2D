#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "Header.h"
#include "Source.cpp"

int main(){
	//InputPoints
	vector<Point> PointList=InputPoints();

	//Make Matrix for Edges
	int** Edges  = EdgeMatrix(PointList);

	//Link Edges 
	Mesh result_mesh = Delaunay(PointList);

	result_mesh.PrintAllTrinagle();
	
	return 0;
}