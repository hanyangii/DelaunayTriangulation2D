#include "Delaunay.h"

#define INF 100000000000000

ofstream coutFile;
vector<int>* adj;
vector<Point> FullPointList;
int cnt=0;
clock_t start, finish;
	

void AddEdge(int n1, int n2){
	assert(n1>=0&&n2>=0);
	if(display){
		coutFile << "add edge : " << n1 << " " << n2 <<"\n" << endl;
	}
	adj[n1].push_back(n2);
	adj[n2].push_back(n1);
}

void DeleteEdge(int n1, int n2){
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

double det(Point* central, Point* start, Point* end){
	double start_x = start->x - central->x, start_y = start->y - central->y;
	double end_x = end->x - central->x, end_y = end->y - central->y;

	return ((start_x * end_y) - (start_y*end_x));
}

double dot(Point* central, Point* start, Point* end){
	double start_x = start->x - central->x, start_y = start->y - central->y;
	double end_x = end->x - central->x, end_y = end->y - central->y;

	return ((start_x * end_x) + (start_y*end_y));
}

double tan(Point* central, Point* start, Point* end){
	double start_x = start->x - central->x, start_y = start->y - central->y;
	double end_x = end->x - central->x, end_y = end->y - central->y;
	double dott= ((start_x * end_x) + (start_y*end_y));
	double dist = sqrt((start_x*start_x+start_y*start_y)*(end_x*end_x+end_y*end_y));
	return -dott/dist;
	//return det(central, start, end)/dot(central, start, end);
}

//Find baseLR from two sets
Edge initialBaseLR(Mesh Lmesh, Mesh Rmesh){
	vector<Point*> Lpoints(Lmesh.PointList.begin(), Lmesh.PointList.end());
	vector<Point*> Rpoints(Rmesh.PointList.begin(), Rmesh.PointList.end());

	if(Lpoints[0]->num==9 && Rpoints[0]->num==17){
		cout<<"now!!"<<endl;
	}

	sort(Lpoints.begin(), Lpoints.end(), compare_point_y);
	sort(Rpoints.begin(), Rpoints.end(), compare_point_y);

	Point* Lmin = Lpoints.front();
	Point* Rmin = Rpoints.front();

	//if(Lmin->y < Rmin->y){
		for(int i =1;i<Rpoints.size();i++){
			if(det(Lmin, Rmin, Rpoints[i])<0) Rmin=Rpoints[i];
			if(Rpoints[i]->y>Rmin->y && Rpoints[i]->y > Lmin->y ) break;
		}
	//}
	//else{
		for(int i=1;i<Lpoints.size();i++){
			if(det(Rmin, Lmin, Lpoints[i])>0) Lmin=Lpoints[i];
			if(Lpoints[i]->y>Rmin->y && Lpoints[i]->y > Lmin->y ) break;
		}
	//}

	Edge baseLR;
	baseLR.AddPoint(Lmin, Rmin);
	return baseLR;

}

vector<Point*> SortCandidates(vector<Point*> OriginalList, Edge baseLR, bool isLMesh){
	Point* Lpoint = baseLR.start, *Rpoint = baseLR.end;
	vector<Point*> newList;
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

	sort(newList.begin(), newList.end(), compare_angle);
	if(display){
		for(int i=0;i<newList.size();i++){
			coutFile<<newList[i]->num<<" : "<<newList[i]->x<<" "<<newList[i]->y<<" tan : "<<newList[i]->angle<<endl;
		}
		coutFile<<"\n"<<endl;
	}
	return newList;
}

//if circumcircle of present point and baseLR points contain next point, return ture
bool isContain(Edge baseLR, Point* now, Point* next){
	Point *c1=baseLR.start, *c2=baseLR.end, *c3=now;
	double dA = c1->x * c1->x + c1->y * c1->y;
    double dB = c2->x * c2->x + c2->y * c2->y;
    double dC = c3->x * c3->x + c3->y * c3->y;
  
    double aux1 = (dA*(c3->y - c2->y) + dB*(c1->y - c3->y) + dC*(c2->y - c1->y));
    double aux2 = -(dA*(c3->x - c2->x) + dB*(c1->x - c3->x) + dC*(c2->x - c1->x));
    double div = (2*(c1->x*(c3->y - c2->y) + c2->x*(c1->y-c3->y) + c3->x*(c2->y - c1->y)));
 
    if(div == 0){ 
        return false;
    }

	//Circumcircle
	double center_x = aux1/div;
	double center_y = aux2/div;
	double double_radius = (center_x - c1->x)*(center_x - c1->x) + (center_y - c1->y)*(center_y - c1->y);
	
	//if(display) cout<<"CircumCircle : x="<<a[0]<<" y="<<a[1]<<" r="<<a[2]<<endl;

	double x_dist = next->x - center_x;
	double y_dist = next->y - center_y;
	double dist = x_dist*x_dist + y_dist*y_dist;

	if(display) coutFile<<"Tri : "<<c1->num<<" "<<c2->num<<" "<<c3->num<<"  Point : "<<next->num<<" :::: "<<double_radius<<" "<<dist<<endl;
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
Mesh MakeLRedge(Mesh new_mesh, Mesh Lmesh, Mesh Rmesh, Edge baseLR){
	if(display) coutFile << "++++++++++++++++++++ newLR : " << baseLR.start->num << " " << baseLR.end->num <<"++++++++++++++++++++++++++++\n" << endl;
	
	//Sort points by angle
	vector<Point*> Lcandidates = SortCandidates(Lmesh.PointList, baseLR , true);
	vector<Point*> Rcandidates = SortCandidates(Rmesh.PointList, baseLR , false);
	
	//Select L candidate
	bool isLcandid = false;
	Point* Lcandid;
	for(int i=0;i<Lcandidates.size();i++){
		if(i==Lcandidates.size()-1){
			Lcandid = Lcandidates[i];
			isLcandid=true;
			if(display) coutFile<<"Lcandid : "<<Lcandid->num<<endl;
			break;
		}
		else{
			bool del=false;
			for(int j=i+1;j<Lcandidates.size();j++){
				if(isContain(baseLR, Lcandidates[i], Lcandidates[j])){
					DeleteEdge(baseLR.start->num, Lcandidates[i]->num);
					del=true;
					break;
				}
			}

			if(!del){
				isLcandid = true;
				Lcandid = Lcandidates[i];
				if(display) coutFile<<"Lcandid : "<<Lcandid->num<<endl;
				break;
			}
		}
	}

	//Select R candidate
	bool isRcandid = false;
	Point* Rcandid;
	for(int i=0;i<Rcandidates.size();i++){
		if(i==Rcandidates.size()-1){
			Rcandid = Rcandidates[i];
			isRcandid=true;
			if(display) coutFile<<"Rcandid : "<<Rcandid->num<<endl;
			break;
		}
		else{
			bool del=false;
			for(int j=i+1;j<Rcandidates.size();j++){
				if(isContain(baseLR, Rcandidates[i], Rcandidates[j])){
					DeleteEdge(baseLR.end->num, Rcandidates[i]->num);
					del=true;
					break;
				}
			}

			if(!del){
				isRcandid = true;
				Rcandid = Rcandidates[i];
				if(display) coutFile<<"Rcandid : "<<Rcandid->num<<endl;
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
			newLR.AddPoint(Lcandid,baseLR.end);
			AddEdge(newLR.start->num, newLR.end->num);
			
		}
		else{
			newLR.AddPoint(baseLR.start, Rcandid);
			AddEdge(newLR.start->num, newLR.end->num);
		}
	}
	else if(isLcandid){
		newLR.AddPoint(Lcandid,baseLR.end);
		AddEdge(newLR.start->num, newLR.end->num);
	}
	else if(isRcandid){
		newLR.AddPoint(baseLR.start, Rcandid);
		AddEdge(newLR.start->num, newLR.end->num);
	}
	else{
		//No candidates
		return new_mesh;
	}

	return MakeLRedge(new_mesh, Lmesh, Rmesh, newLR);
}

//Divide and Conquer
Mesh MergeSets(Mesh Lmesh, Mesh Rmesh){
	if(display){
			coutFile << "Lset : " <<endl;
			for(vector<Point*>::iterator i =Lmesh.PointList.begin(); i!=Lmesh.PointList.end();i++){
				coutFile << (*i)->num << " " << endl;
			}
			coutFile << "Rset : " <<endl;
			for(vector<Point*>::iterator i =Rmesh.PointList.begin(); i!=Rmesh.PointList.end();i++){
				coutFile << (*i)->num << " " << endl;
			}
			coutFile << "\n"<<endl;
		}


	// Divide each mesh to subsets OR Make Edge and triangle
	if(Lmesh.point_num()>3){
		size_t list_size;
		if(Lmesh.point_num()%2 == 1) list_size = Lmesh.point_num()/2+1;
		else list_size = Lmesh.point_num()/2;
		
		size_t const half_size = list_size;

		vector <Point*> new_LList(Lmesh.PointList.begin(), Lmesh.PointList.begin()+half_size);
		vector <Point*> new_RList(Lmesh.PointList.begin()+half_size, Lmesh.PointList.end());
		Mesh new_Lmesh(new_LList);
		Mesh new_Rmesh(new_RList);
		Lmesh = MergeSets(new_Lmesh, new_Rmesh);
	}
	else{
		//if number of points <= 3
		for(vector<Point*>::iterator i=Lmesh.PointList.begin(); i!=Lmesh.PointList.end(); i++){
			for(vector<Point*>::iterator j = i+1; j!= Lmesh.PointList.end();j++){
				AddEdge((*i)->num, (*j)->num);
			}
		}
	}

	if(Rmesh.point_num()>3){
		size_t list_size;
		if(Rmesh.point_num()%2 == 1) list_size = Rmesh.point_num()/2+1;
		else list_size = Rmesh.point_num()/2;
		
		size_t const half_size = list_size;
		vector <Point*> new_LList(Rmesh.PointList.begin(), Rmesh.PointList.begin()+half_size);
		vector <Point*> new_RList(Rmesh.PointList.begin()+half_size, Rmesh.PointList.end());
		Mesh new_Lmesh(new_LList);
		Mesh new_Rmesh(new_RList);
		Rmesh = MergeSets(new_Lmesh, new_Rmesh);
	}
	else{
		//if number of points <=3
		for(vector<Point*>::iterator i=Rmesh.PointList.begin(); i!=Rmesh.PointList.end(); i++){
			for(vector<Point*>::iterator j = i+1; j!= Rmesh.PointList.end();j++){
				AddEdge((*i)->num, (*j)->num);
			}
		}
	}

	//Conquer two sets
	vector <Point*> new_PointList(Lmesh.PointList.begin(), Lmesh.PointList.end());
	new_PointList.insert(new_PointList.end(), Rmesh.PointList.begin(), Rmesh.PointList.end());
	Mesh new_mesh(new_PointList);

	//Select the lowest line for baseLR
	Edge baseLR = initialBaseLR(Lmesh, Rmesh);
	AddEdge(baseLR.start->num, baseLR.end->num);

	//Link LR edges
	cnt =0;
	new_mesh = MakeLRedge(new_mesh, Lmesh, Rmesh, baseLR);
	
	if(display) coutFile << "\n===============================================Two sets merged====================================\n" << endl;
	return new_mesh;
}



Mesh Delaunay(){
	
	ofstream newFile;
	newFile.open("newPoints.txt");
	
	//Sort points from left to right & give indicies
	int num=0;
	for(vector<Point>::iterator i=FullPointList.begin();i!=FullPointList.end();i++){
		newFile<< i->x << "\t" << i->y << endl;
		i->addNum(num++);	
	}
	newFile.close();
	
	bool erase = false ;
	assert(FullPointList.size()>=2);

	//Devide Points to two sets and merge
	start=clock();
	if(FullPointList.size()>3){
		if(display) cout << "number of points : " << FullPointList.size()<<"\n"<<endl;
		int list_size =0;
		if(FullPointList.size()%2 == 1) list_size = FullPointList.size()/2+1;
		else list_size = FullPointList.size()/2;
		
		size_t const half_size = list_size;
		vector<Point*> LList;
		vector<Point*> RList;
		for(int i = 0;i<half_size;i++) LList.push_back(&FullPointList[i]);
		for(int i=half_size;i<FullPointList.size();i++) RList.push_back(&FullPointList[i]);

		Mesh LMesh(LList);
		Mesh RMesh(RList);
		return MergeSets(LMesh, RMesh);
	}
	else if(FullPointList.size()==3){
		vector<Point*> Plist;
		for(int i=0;i<FullPointList.size();i++) Plist.push_back(&FullPointList[i]);
		Mesh mesh(Plist);
		//vector<Point>::iterator it = FullPointList.begin();
		for(int i=0;i<3;i++){
			AddEdge(FullPointList[i].num, FullPointList[(i+1)%3].num);
			
		}
		return mesh;
	}
	else{
		vector<Point*> Plist;
		for(int i=0;i<FullPointList.size();i++) Plist.push_back(&FullPointList[i]);
		Mesh mesh(Plist);
		AddEdge(FullPointList.begin()->num, FullPointList.end()->num);
		return mesh;
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
	
	//If there are two same points, move one point little bit. 
	
	sort(PointList.begin(), PointList.end(), compare_x);
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
	FullPointList=InputPoints();

	//Make Matrix for Edges
	int point_num = FullPointList.size();
	adj = new vector<int> [ point_num ];

	
	//Link Edges 
	Mesh result_mesh = Delaunay();
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