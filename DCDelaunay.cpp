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

#define INF 100000000000000
#define PI 3.14159265

using namespace std;

ofstream coutFile;
bool display = false;

class Point
{
	public:
		double x;
		double y;
		double angle;// angle between this point and min y point 
		int idx;

		bool operator==(const Point &rhs) const{
			return x==rhs.x && y==rhs.y;
		}
};

typedef vector<const Point*> PointList;

class Edge{
	// start point is located at left of end point in every edge
	public:
		const Point* start;
		const Point* end;
		
		Edge() : start(NULL), end(NULL) {}
		Edge(const Point* s, const Point* e) : start(s), end(e) {}
		void set(const Point* s, const Point* e) { start = s, end = e; }

};

bool operator==(const Edge &rhs, const Edge&lhs){
	return rhs.start == lhs.start && rhs.end == lhs.end;
}

class Mesh{
	public:
		vector<const Point*> PointList;

		size_t point_num(){
			return PointList.size();
		}

		void InputList(vector<const Point*> _PointList){
			for(vector<const Point*>:: iterator it =_PointList.begin();it!=_PointList.end();it++){
				const Point* ppointer  = *it;
				PointList.push_back(ppointer);
			}
		}	
};

bool compare_x(const Point& a, const Point& b){
    return a.x < b.x;
}


bool compare_point_y(const Point* a, const Point* b){
    return a->y < b->y;
}

void AddEdge(int n1, int n2, vector<int>* adj){
	assert(n1>=0&&n2>=0);
	if(display){
		coutFile << "add edge : " << n1 << " " << n2 <<"\n" << endl;
	}
	adj[n1].push_back(n2);
	adj[n2].push_back(n1);
}

void DeleteEdge(int n1, int n2, vector<int>* adj){
	if(true) coutFile <<"\n****** Delete : "<<n1<<" "<<n2<<"******\n"<<endl;
	for(int i = 0; i<adj[n1].size();i++){
		if(adj[n1].at(i)==n2){
			adj[n1].erase(adj[n1].begin()+i);
			break;
		}
	}

	for(int i = 0; i<adj[n2].size();i++){
		if(adj[n2].at(i)==n1){
			adj[n2].erase(adj[n2].begin()+i);
			break;
		}
	}
}

double det(const Point* central, const Point* start, const Point* end){
	double start_x = start->x - central->x, start_y = start->y - central->y;
	double end_x = end->x - central->x, end_y = end->y - central->y;

	return ((start_x * end_y) - (start_y*end_x));
}

double dot(const Point* central, const Point* start, const Point* end){
	double start_x = start->x - central->x, start_y = start->y - central->y;
	double end_x = end->x - central->x, end_y = end->y - central->y;

	return ((start_x * end_x) + (start_y*end_y));
}

double tan(const Point* central, const Point* start, const Point* end){
	double start_x = start->x - central->x, start_y = start->y - central->y;
	double end_x = end->x - central->x, end_y = end->y - central->y;
	double dott= ((start_x * end_x) + (start_y*end_y));
	double dist = sqrt((start_x*start_x+start_y*start_y)*(end_x*end_x+end_y*end_y));
	return -dott/dist;
	//return det(central, start, end)/dot(central, start, end);
}

//Find baseLR from two sets
Edge FindInitialBaseLR(const PointList& left, const PointList& right){
	vector<const Point*> Lpoints(left.begin(), left.end());
	vector<const Point*> Rpoints(right.begin(), right.end());
	
	sort(Lpoints.begin(), Lpoints.end(), compare_point_y);
	sort(Rpoints.begin(), Rpoints.end(), compare_point_y);

	const Point* Lmin = Lpoints.front();
	const Point* Rmin = Rpoints.front();

	for(int i =1;i<Rpoints.size();i++){
			if(det(Lmin, Rmin, Rpoints[i])<0) Rmin=Rpoints[i];
			if(Rpoints[i]->y>Rmin->y && Rpoints[i]->y > Lmin->y ) break;
	}
	for(int i=1;i<Lpoints.size();i++){
			if(det(Rmin, Lmin, Lpoints[i])>0) Lmin=Lpoints[i];
			if(Lpoints[i]->y>Rmin->y && Lpoints[i]->y > Lmin->y ) break;
	}

	return Edge(Lmin, Rmin);

}


void SortCandidates(const PointList pts, vector<const Point*>* newList, Edge baseLR,bool is_left_mesh){
	const Point* Lpoint = baseLR.start, *Rpoint = baseLR.end;
	vector<pair<double, const Point*> > data;
	///////////////////////HERE///////////////////
	if (is_left_mesh) {
		for (int i = 0; i < pts.size(); i++) {
			const Point* p = pts[i];
			if (det(Lpoint, Rpoint, p)>0) data.push_back(make_pair(tan(Lpoint, Rpoint, p), p));
		}
	} 
	else {  // !is_left_mesh
		for (int i = 0; i < pts.size(); i++) {
			const Point* p = pts[i];
			if (det(Rpoint, p, Lpoint)>0) data.push_back(make_pair(tan(Rpoint, Lpoint, p), p));
		}
	}
	/*
	for(int i=0;i<OriginalList.size();i++){
		if(isLMesh&&det(Lpoint, Rpoint, OriginalList[i])>0){
			OriginalList[i]->angle=tan(Lpoint, Rpoint, OriginalList[i]);
			newList.push_back(OriginalList[i]);
		}
		else if(!isLMesh&&det(Rpoint, OriginalList[i],Lpoint)>0){
			OriginalList[i]->angle=tan(Rpoint, Lpoint, OriginalList[i]);
			newList.push_back(OriginalList[i]);
		}
	}
	*/
	sort(data.begin(), data.end());
	newList->resize(data.size());
	for(int i = 0;i<data.size();++i) newList ->at(i) = data[i].second;
	if(display){
		for(int i=0;i<newList->size();i++){
			coutFile<<newList->at(i)->idx<<" : "<<newList->at(i)->x<<" "<<newList->at(i)->y<<" tan : "<<newList->at(i)->angle<<endl;
		}
		coutFile<<"\n"<<endl;
	}

}

//if circumcircle of present point and baseLR points contain next point, return ture
bool isContain(Edge baseLR, const Point* now, const Point* next){
	const Point *c1=baseLR.start, *c2=baseLR.end, *c3=now;
	const double dA = c1->x * c1->x + c1->y * c1->y;
    const double dB = c2->x * c2->x + c2->y * c2->y;
    const double dC = c3->x * c3->x + c3->y * c3->y;
  
    const double aux1 = (dA*(c3->y - c2->y) + dB*(c1->y - c3->y) + dC*(c2->y - c1->y));
    const double aux2 = -(dA*(c3->x - c2->x) + dB*(c1->x - c3->x) + dC*(c2->x - c1->x));
    const double div = (2*(c1->x*(c3->y - c2->y) + c2->x*(c1->y-c3->y) + c3->x*(c2->y - c1->y)));
 
    if(div == 0){ 
        return false;
    }

	//Circumcircle
	const double center_x = aux1/div;
	const double center_y = aux2/div;
	const double double_radius = (center_x - c1->x)*(center_x - c1->x) + (center_y - c1->y)*(center_y - c1->y);
	
	//if(display) cout<<"CircumCircle : x="<<a[0]<<" y="<<a[1]<<" r="<<a[2]<<endl;

	const double x_dist = next->x - center_x;
	const double y_dist = next->y - center_y;
	const double dist = x_dist*x_dist + y_dist*y_dist;

	if(display) coutFile<<"Tri : "<<c1->idx<<" "<<c2->idx<<" "<<c3->idx<<"  Point : "<<next->idx<<" :::: "<<double_radius<<" "<<dist<<endl;
	if(double_radius > dist){
		if(display) coutFile<<"True"<<endl;
		return true;
	}
	else{
		if(display) coutFile<<"False"<<endl;
		return false;
	}
}

//Make Next LR edge (reculsively)
void MakeLRedge(const PointList& left, const PointList& right,
                const Edge& baseLR,
                Mesh* mesh, vector<int>* adj){
	if(display) coutFile << "++++++++++++++++++++ newLR : " << baseLR.start->idx << " " << baseLR.end->idx <<"++++++++++++++++++++++++++++\n" << endl;
	
	//Sort points by angle
	PointList Lcandidates, Rcandidates;
	SortCandidates(left, &Lcandidates, baseLR , true);
	SortCandidates(right, &Rcandidates, baseLR , false);
	
	//Select L candidate
	bool isLcandid = false;
	const Point* Lcandid;
	for(int i=0;i<Lcandidates.size();i++){
		if(i==Lcandidates.size()-1){
			Lcandid = Lcandidates[i];
			isLcandid=true;
			if(display) coutFile<<"Lcandid : "<<Lcandid->idx<<endl;
			break;
		}
		else{
			bool del=false;
			for(int j=i+1;j<Lcandidates.size();j++){
				if(isContain(baseLR, Lcandidates[i], Lcandidates[j])){
					DeleteEdge(baseLR.start->idx, Lcandidates[i]->idx, adj);
					del=true;
					break;
				}
			}

			if(!del){
				isLcandid = true;
				Lcandid = Lcandidates[i];
				if(display) coutFile<<"Lcandid : "<<Lcandid->idx<<endl;
				break;
			}
		}
	}

	//Select R candidate
	bool isRcandid = false;
	const Point* Rcandid;
	for(int i=0;i<Rcandidates.size();i++){
		if(i==Rcandidates.size()-1){
			Rcandid = Rcandidates[i];
			isRcandid=true;
			if(display) coutFile<<"Rcandid : "<<Rcandid->idx<<endl;
			break;
		}
		else{
			bool del=false;
			for(int j=i+1;j<Rcandidates.size();j++){
				if(isContain(baseLR, Rcandidates[i], Rcandidates[j])){
					DeleteEdge(baseLR.end->idx, Rcandidates[i]->idx, adj);
					del=true;
					break;
				}
			}

			if(!del){
				isRcandid = true;
				Rcandid = Rcandidates[i];
				if(display) coutFile<<"Rcandid : "<<Rcandid->idx<<endl;
				break;
			}
		}
	}
	

	//Choose one in two candidates
	bool isRcontainL =false, isLcontainR = false;
	Edge newLR;
	if(isLcandid&&isRcandid){
		isRcontainL = isContain(baseLR, Rcandid, Lcandid);
		isLcontainR = isContain(baseLR, Lcandid, Rcandid);

		assert(isRcontainL || isLcontainR == true);
		if(isRcontainL){
			newLR.set(Lcandid,baseLR.end);
			AddEdge(newLR.start->idx, newLR.end->idx, adj);
			
		}
		else{
			newLR.set(baseLR.start, Rcandid);
			AddEdge(newLR.start->idx, newLR.end->idx, adj);
		}
	}
	else if(isLcandid){
		newLR.set(Lcandid,baseLR.end);
		AddEdge(newLR.start->idx, newLR.end->idx, adj);
	}
	else if(isRcandid){
		newLR.set(baseLR.start, Rcandid);
		AddEdge(newLR.start->idx, newLR.end->idx, adj);
	}
	else{
		//No candidates
		return ;
	}

	 MakeLRedge(left, right, newLR, mesh, adj);
}

//Divide and Conquer
void MergeSets(PointList left, PointList right,
               Mesh* mesh, vector<int>* adj){
	if(display){
			coutFile << "Lset : " <<endl;
			for(PointList::iterator i =left.begin(); i!=left.end();i++){
				coutFile << (*i)->idx << " " << endl;
			}
			coutFile << "Rset : " <<endl;
			for(PointList::iterator i =right.begin(); i!=right.end();i++){
				coutFile << (*i)->idx << " " << endl;
			}
			coutFile << "\n"<<endl;
		}


	// Divide each mesh to subsets OR Make Edge and triangle
	if(left.size()>3){
		const size_t half_size = ceil(left.size() / 2.0);
		PointList new_LList(left.begin(), left.begin()+half_size);
		PointList new_RList(left.begin()+half_size, left.end());
		MergeSets(new_LList, new_RList, mesh, adj);
	}
	else{
		//if number of points <= 3
		for (int i = 0; i < left.size(); ++i) {
			for (int j = i + 1; j < left.size(); ++j) {
				AddEdge(left[i]->idx, left[j]->idx, adj);
			}
		}
	}

	if(right.size()>3){
		const size_t half_size = ceil(right.size() / 2.0);
		PointList new_LList(right.begin(), right.begin()+half_size);
		PointList new_RList(right.begin()+half_size, right.end());
		MergeSets(new_LList, new_RList, mesh, adj);
	}
	else{
		//if number of points <= 3
		for (int i = 0; i < right.size(); ++i) {
			for (int j = i + 1; j < right.size(); ++j) {
				AddEdge(right[i]->idx, right[j]->idx, adj);
			}
		}
	}

	//Conquer two sets
	/*ector <Point*> new_PointList(LL.begin(), Lmesh.PointList.end());
	new_PointList.insert(new_PointList.end(), Rmesh.PointList.begin(), Rmesh.PointList.end());
	Mesh new_mesh(new_PointList);
	*/

	//Select the lowest line for baseLR
	Edge baseLR = FindInitialBaseLR(left, right);
	AddEdge(baseLR.start->idx, baseLR.end->idx, adj);

	//Link LR edges
	MakeLRedge(left, right, baseLR, mesh, adj);
	
	if(display) coutFile << "\n===============================================Two sets merged====================================\n" << endl;
}



void Delaunay(const vector<Point>& pts, Mesh* mesh, vector<int>* adj){
		
	bool erase = false ;
	assert(pts.size()>=2);

	//Devide Points to two sets and merge
	
	if(pts.size()>3){
		if(display) cout << "number of points : " << pts.size()<<"\n"<<endl;
		size_t const half_size = ceil(pts.size() / 2.0);
		vector<const Point*> LList;
		vector<const Point*> RList;
		for(int i = 0;i<half_size;i++) LList.push_back(&pts[i]);
		for(int i=half_size;i<pts.size();i++) RList.push_back(&pts[i]);
		MergeSets(LList, RList, mesh, adj);
	}
	else if(pts.size()==3){
		vector<const Point*> Plist;
		for(int i=0;i<pts.size();i++) Plist.push_back(&pts[i]);
		Mesh mesh;
		mesh.InputList(Plist);
		//vector<Point>::iterator it = pts.begin();
		for(int i=0;i<3;i++){
			AddEdge(pts[i].idx, pts[(i+1)%3].idx, adj);
			
		}
	}
	else{
		vector<const Point*> Plist;
		for(int i=0;i<pts.size();i++) Plist.push_back(&pts[i]);
		Mesh mesh;
		mesh.InputList(Plist);
		AddEdge(pts.begin()->idx, pts.end()->idx, adj);
	}
	
}

// Save data from input file to Point List
vector<Point> InputPoints(){
	//Input points
	ifstream PointFile("Points.txt");
	vector<Point> PointList;

	
	string line;
	
	int tt=0;
	while(!PointFile.eof()){
		getline(PointFile, line);
		stringstream linestream(line);
		string word;
		int idx = 0;
		Point new_pt;
		
		while(getline(linestream, word, '\t')){
			double inputnum = atof((word).c_str());
			if(idx%2==0){ 
				new_pt.x= inputnum;//*((double)rand()/(double)RAND_MAX);
			}
			else{
				new_pt.y= inputnum;//*((double)rand()/(double)RAND_MAX);
			}
			if(++idx ==2) break;
		} 
		PointList.push_back(new_pt);

		tt++;
	}
	PointFile.close();


	Point last_pt;
	last_pt.x = INF;
	last_pt.y = INF;
	double epsilon = 0.001;
	
	//Sort Points
	sort(PointList.begin(), PointList.end(), compare_x);
	
	//If there are two same points, move one point little bit. 
	for(int i=0;i<PointList.size();i++){
		Point current  = PointList[i];
		if(last_pt.x == current.x || last_pt.y==current.y){
			PointList[i].x+=((double)rand()/(double)RAND_MAX*epsilon );
			PointList[i].y+=((double)rand()/(double)RAND_MAX*epsilon );
		}
		last_pt=current;
	}
	
	return PointList;
	
}


int main(){
	//InputPoints	
	coutFile.open("Cout.txt");
	srand((unsigned int)time(NULL));
	vector<Point> pts;
	pts=InputPoints();

	//Make Matrix for Edges
	int point_num = pts.size();
	vector<int>* adj;
	adj = new vector<int> [ point_num ];

	//give indicies
	ofstream newFile;
	newFile.open("newPoints.txt");
	for(int i = 0; i<pts.size(); i++){
		newFile << pts[i].x << "\t" << pts[i].y << endl;
		pts[i].idx = i;
	}
	newFile.close();

	//Link Edges 
	clock_t start, finish;
	Mesh result_mesh;
	start=clock();
	Delaunay(pts, &result_mesh, adj);
	finish=clock();

	double duration = (double)(finish - start)/CLOCKS_PER_SEC;
	cout << "Duration : "<<duration<<endl;

	ofstream linkFile;
	linkFile.open("Link.txt");
	for(int i=0; i<point_num;i++){
		for(int j=0; j < adj[i].size();j++){
			linkFile << i << " " << adj[i].at(j) << endl;
		}
	}
	
	linkFile.close();
	coutFile.close();
	return 0;
}