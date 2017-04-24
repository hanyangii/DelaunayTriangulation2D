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
		int num =0;
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

struct compare_x{
	  bool operator(Point* a, Point* b){
		     return a -> x < b->x;
		}
}
struct compare_y{
	  bool operator(Point* a, Point* b){
		   return a->y < b->y;
		}
}

Edge linkEdge(Point start, Point end){
	Edge e;
	e.start = start;
	e.end = end;
	return e;
}

Triangle MakeTriangle(Point p1, Point p2, Point p3){
	Triangle triangle;
	
	triangle.points[0]=p1;
	triangle.points[1]=p2;
	triangle.points[2]=p3;
	
	for(int i=0; i<3; i++){
		Edge e;
		e.start = triangle.points[i];
		e.start = triangle.points[(i+1)%3];
		triangle.EdgeList.push_back(e)
	}

	return triangle;
}

Point MinYpoint(Triangle tri){
	float a=tri.points[0].y, b=tri.points[1].y, c=tri.points[2].y;
	if(a>b && a>c) return a;
	else if(b>a&&b>c) return b;
	else return c;
}
