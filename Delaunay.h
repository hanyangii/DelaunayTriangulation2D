#pragma once
#define DELAUNAY_HEADER

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <list>
#include <algorithm>
#include <time.h>
#include <cmath>
#include <assert.h>

#define PI 3.14159265
#define INF 1000

using namespace std;


bool display = false;

class Point
{
	public:
		double x;
		double y;
		double angle;// angle between this point and min y point 
		int num;

		bool operator==(const Point &rhs) const{
			return x==rhs.x && y==rhs.y;
		}

		void addNum(int _num){
			num=_num;
		}
};


bool cmp(const pair<double, double> &u, const pair<double, double>&v){
	if(u.first<v.first) return true;
	else if(u.first==v.first) return u.second<v.second;
	else return false;
}


class Edge{
	// start point is located at left of end point in every edge
	public:
		Point* start;
		Point* end;
		
		friend bool operator== (const Edge &rhs, const Edge &lhs);

		void AddPoint(Point* _start, Point* _end){
			start=_start;
			end=_end;
		}

};


bool operator==(const Edge &rhs, const Edge&lhs){
	return rhs.start == lhs.start && rhs.end == lhs.end;
}

class Mesh
{

	public:
		vector<Point*> PointList;

		size_t point_num(){
			return PointList.size();
		}

		Mesh(vector<Point*> _PointList){
			for(vector<Point*>:: iterator it =_PointList.begin();it!=_PointList.end();it++){
				Point* ppointer  = *it;
				PointList.push_back(ppointer);
			}
		}

	
		
};




bool compare_x(const Point& a, const Point& b){
    return a.x < b.x;
}

bool compare_y(const Point& a, const Point&b){
    return a.y < b.y;
}

bool compare_point_y(Point* a, Point* b){
    return a->y < b->y;
}

bool compare_angle(Point* a, Point* b){
			return a->angle < b->angle;
}



Point* MinYpoint(Mesh mesh){
	vector <Point*> Plist(mesh.PointList.begin(), mesh.PointList.end());
	sort(Plist.begin(), Plist.end(),compare_point_y);
	return Plist.front();
}

int arg_min(const double *arr, size_t length) {
    // returns the minimum value of array
    size_t i;
	int arg=0;	
    double minimum = arr[0];
    for (i = 1; i < length; ++i) {
        if (minimum > arr[i]) {
            minimum = arr[i];
			arg =i;
        }
    }
    return arg;
}

Point SortingAngle(Point minp, vector<Point> llist , bool isLmin){
	double *pp = new double [ llist.size() ]();
	int t=0;
	for(vector<Point> :: iterator i = llist.begin(); i!=llist.end();i++){
		pp[t++] = (i->y - minp.y)/(i->x - minp.x);
		if(!isLmin){
			pp[t-1]=-pp[t-1];
		}
		if(display){
			cout << i->num << "±â¿ï±â : "<<pp[t-1]<<endl; 
		}
	}

	int arg = arg_min(pp, llist.size());
	Point n_min_point = llist.at(arg);
	
	return n_min_point;
}


bool IsIntersection(double* circle, Point base_p, Point n){
	double vec2_x = circle[0]-base_p.x;
	double vec2_y = circle[1]-base_p.y;
	double vec1_x = n.x-base_p.x;
	double vec1_y = n.y-base_p.y;
  
	double det = (vec1_x*vec2_y)-(vec1_y*vec2_x);
	double dott = (vec1_x*vec2_x)+(vec1_y*vec2_y);

	if(dott>0 && det >0) return true;
	else return false;
}



bool isInter(Point a1, Point a2, Point b1, Point b2){
	if(display) cout<<a1.x<<" "<<a1.y<<"-"<<a2.x<<" "<<a2.y<<endl;
	if(display) cout<<b1.x<<" "<<b1.y<<"-"<<b2.x<<" "<<b2.y<<endl;
	if(a1.y<b1.y && a2.y>b2.y&&a2.x>b1.x&&a2.y>b1.y) return true;
	else return false;
}
//#endif