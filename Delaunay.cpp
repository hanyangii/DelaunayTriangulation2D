#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <list>

using namespace std;

struct Point{
	float x=0.0;
	float y=0.0;
};

struct Edge{
	// start.x < end.x
	Point start;
	Point end;
};

struct Triangle{
	//clockwise from a point with the least x value
	Point points[3];
	list<Edge> EdgeList;	
};

list<Triangle> TriangleList;

list<Triangle> MergeTriangles(list<Triangle> LList, list<Triangle> RList){
	list<Triangle> merged_list;
	
	return merged_list;
}

list<Triangle> Delaunay(list<Point> PointList){
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
