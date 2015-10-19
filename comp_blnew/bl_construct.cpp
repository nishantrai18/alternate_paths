#include <iostream>
#include <algorithm>
#include <cstdio>
#include <queue>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <map>
#include <stack>
#include <cmath>

using namespace std;

#define INF 1000000000
#define N 310000

struct edge{
    int end_node;
    float weight;
};

typedef struct edge edge;

typedef pair< int , int> epair;

typedef pair< double, double > ponts;


bool operator <( edge a, edge b ) 
{
	return a.weight > b.weight;
}

struct sort_pred {
    bool operator()(const std::pair<int,int> &left, const std::pair<int,int> &right) {
        return left.second < right.second;
    }
};

////////////////

map<epair, float> lengths, olen, revlengths, revolen;
map<long long int, int> vertices;
map<int, pair<double,double> > coords;

vector <vector <edge> > graph,revgraph,prevtree,foratree;
vector <int> prev,fora,origpath;

float landmarks_for[40][N+10],landmarks_rev[40][N+10],forcost[N+10],revcost[N+10],forpi[N+10],revpi[N+10];
float dista[N+10],distb[N+10],opt,param_share,param_stretch,param_localt;
int visa[N+10],visb[N+10],stat[N+10],sz,size,cnt=0;

double altp=0,totp=0,giters=0,initp=0,altp_dist[1100],origopt,avglt=0,avgshr=0,avgp=0;
int num_landmarks=32;

double shr_pt[1100],st_pt[1100];

int sz_pt=0;

//////////////////

int graph_print(vector <vector <edge> > graph, int size)
{
	int i,j,t;
	for(i=0;i<=size;i++)
    {
    	cout<<i<<":";
    	for(j=0;j<graph[i].size();j++)
    	{
    		cout<<graph[i][j].end_node<<" ";
    	}
    	if(graph[i].size())
    		cout<<"\n";
    }
    return 0;
}

void marks_plt(vector <int> marks, FILE* f)
{
	int i,j;
	fprintf(f,"[");
	for(i=0;i<marks.size()-1;i++)
	{
		fprintf(f,"[%lf,%lf],",coords[marks[i]].first,coords[marks[i]].second);
	}
	fprintf(f,"[%lf,%lf]",coords[marks[i]].first,coords[marks[i]].second);
	fprintf(f,"]\n");
}

double norm(double n)
{
	return n/1000000;
}

double distance(double a,double b,double c,double d)
{
	return sqrt(((a-c)*(a-c))+((b-d)*(b-d)));
}

double mod(double x)
{
	if(x<0)
		return -x;
	return x;
}

string tstring(int n)
{
	string str="";
	if(n==0)
		return "0";
	int i=0;
	while(n>0)
	{
		str+=string(1,(char)('0'+n%10));
		n/=10;
	}
	return str;
}

void dijkstra( vector < vector < edge > > graph, float dist[], int size, int node );
float biastar_forlength(float forcost[], float revcost[], float dista[], float distb[], int visa[], int visb[], int na, int nb, int* opv );
int penalise(vector <int> optpath, float alpha);
int findpath(vector <int> &path, int na, int nc, int nb);
double find_dist(int na, int nb, int* opv);
double comp_share(vector <int> optpath, vector <int> path, double* st);
void plotty(int na, int nb, vector <double> &alphas);
float bidijkstra(int na, int nb, int *opv);
double compute_mid(double l, double r, double* shr, double* st, int na, int nb);
ponts values(double p, int na, int nb);
void trivial_plotty(int na, int nb);	
void revert();

int main()
{
	srand (static_cast <unsigned> (time(0)));

	int m,i,j,x,y,iters;

	char c;

	cin>>c;
	
    cout<<"Input Size, edges : ";
    cin>>sz>>m;
    
    size=sz;
    float share,stretch,localt,threshold;
	
	share=0.8,stretch=0.25,localt=0.15;
	threshold=0.05;

	int sum=0,tmp=iters,v=1;

    graph.resize(sz+2);
    revgraph.resize(sz+2);
    prev.resize(sz+2);
    fora.resize(sz+2);

    for(i=0;i<m;i++)
    {
    	cin>>c;
    	int x,y,mx,my;
    	float a;
        scanf("%d %d %f",&x,&y,&a);
        if(vertices.find(x)==vertices.end())
        {
        	vertices.insert(make_pair(x,v));
        	v++;
        }
        if(vertices.find(y)==vertices.end())
        {
        	vertices.insert(make_pair(y,v));
        	v++;
        }
        edge t;
        t.end_node=vertices[y];
        t.weight=a;
        graph[vertices[x]].push_back(t);
    }

    cout<<"INPUT DONE, Actual Size is "<<v<<"\n";
    size=sz=v-1;
    //cin>>v;

    //undir_random_graph(graph,size,threshold);   

    for(i=1;i<=size;i++)
    {
    	for(j=0;j<graph[i].size();j++)
    	{
    		int x=i,y=graph[i][j].end_node;
    		float a=graph[i][j].weight;
    		epair p1 (x,y);
	        if(lengths.find(p1)==lengths.end())
	        {
	        	lengths.insert(make_pair(p1,a));
	        	olen.insert(make_pair(p1,a));;
	        }
	        else
	        {
				lengths[p1]=min(lengths[p1],a);
				olen[p1]=min(olen[p1],a);
			}
	        epair p2 (y,x);
	        if(revlengths.find(p2)==revlengths.end())
	        {
	        	revlengths.insert(make_pair(p2,a));
	        	revolen.insert(make_pair(p2,a));;
	        }
	        else
	        {
				revlengths[p2]=min(revlengths[p2],a);
				revolen[p2]=min(revolen[p2],a);
			}
	        /*epair p2 (y,x);
	        lengths.insert(make_pair(p2,a));
    		*/
    	}
    }

    for(i=0;i<=size;i++)
    {
        epair p1 (0,i);
        lengths.insert(make_pair(p1,0));
        olen.insert(make_pair(p1,0));
        epair p2 (i,0);
        lengths.insert(make_pair(p2,0));
        olen.insert(make_pair(p2,0));
    	for(j=0;j<graph[i].size();j++)
    	{
    		edge t;
    		t.end_node=i;
    		t.weight=graph[i][j].weight;
    		revgraph[graph[i][j].end_node].push_back(t);
    	}
    }
	
    cout<<"REVERSING DONE, COOL\n";

	int jump=(size/num_landmarks);
	
	for(i=1;i<=num_landmarks+1;i++)
	{
		for(x=0;x<=size;x++)
		{
			landmarks_for[i][x]=0;
		}	
	}

	for(x=0;x<=size;x++)
	{
		stat[x]=1;
	}

	cout<<"REACHED THE PREPROCESSING FOR ALT\n";

	for(i=1;i<=num_landmarks;i++)
	{
		float mx=0;
		j=1;
		for(x=1;x<=size;x++)
		{
			if(stat[x]&&(landmarks_for[num_landmarks+1][x]>mx))
			{
				j=x;
				mx=landmarks_for[num_landmarks+1][x];
			}
		}	
		stat[j]=0;

		cout<<"J IS "<<j<<", ";

		for(x=0;x<=size;x++)
		{
			landmarks_for[i][x]=landmarks_rev[i][x]=INF;
		}

		dijkstra(graph,landmarks_for[i],size,j);
		dijkstra(revgraph,landmarks_rev[i],size,j);

		for(x=0;x<=size;x++)
		{
			landmarks_for[num_landmarks+1][x]+=sqrt(landmarks_for[i][x]);
		}
		cout<<i<<" is done\n";
	}

	FILE* f;
	f=fopen("newyork.co.gr","r");
	int tn;
	double cx,cy;
	fscanf(f,"%d",&tn);
	for(i=0;i<tn;i++)
	{
		fscanf(f,"%c",&c);
		fscanf(f,"%c",&c);
		fscanf(f,"%d %lf %lf",&j,&cx,&cy);
		//cout<< j<< "\t";
		j=vertices[j];
		//cout<< j<< endl;
		cx=norm(cx);
		cy=norm(cy);
		coords.insert(make_pair(j,make_pair(cy,cx)));
	}

	fclose(f);
	
	int chck=5,tcase[6][2]={{246102,113929},{63180,195353},{220182,88009},{86316,218489},{158970,26797},{19295,151468}};

	for(int gt=0;gt<=chck;gt++)
		cout<<"add_tree("<<"paths_"+tstring(gt)<<", \""<<"H_"+tstring(gt)+"\") ;\n";
	for(int gt=0;gt<=chck;gt++)
		cout<<"<script src=\""<<"pth_files/paths_"+tstring(gt)+".js\"></script>\n";

	while(chck>=0)
	{
		int S=tcase[chck][0],T=tcase[chck][1],opv;

		//int S=rand()%size;
		//int T=(S+(size/2))%size,opv;

		//fscanf(f,"%d %d",&S,&T);

		FILE* bl=fopen(("bl_files/blackline"+tstring(chck)+".bl").c_str(),"w");
		FILE* tbl=fopen(("bl_files/trivial_blackline"+tstring(chck)+".bl").c_str(),"w");
		FILE* pth=fopen(("pth_files/paths_"+tstring(chck)+".js").c_str(),"w");

		for ( i = 0; i <= size; ++i )
		{
			dista[i] = INF;
			distb[i] = INF;
			visa[i]=visb[i]=0;
		}

		origopt=find_dist( S, T, &opv);
	    origpath.clear();
		findpath(origpath,S,opv,T);

		cout<<S<<" "<<T<<endl;

		vector <double> alphas;

		plotty(S,T,alphas);

		priority_queue < ponts > q;	

		for(i=0;i<sz_pt;i++)
			q.push((ponts) {shr_pt[i],st_pt[i]});

		fprintf(pth,"var %s = [\n",("paths_"+tstring(chck)).c_str());
		for(i=0;i<alphas.size();i++)
		{
			penalise(origpath,alphas[i]);
			int opv;
			find_dist(S, T, &opv);
			vector <int> path;
			findpath(path, S, opv, T);
			marks_plt(path,pth);
			if(i!=(alphas.size()-1))
				fprintf(pth,",");
		}
		fprintf(pth,"];\n");

		while(!q.empty())
		{
			ponts t=q.top();
			q.pop();
			fprintf(bl,"%lf %lf\n",t.first,t.second);
		}

		sz_pt=0;

		trivial_plotty(S,T);

		priority_queue < ponts > qu;	

		for(i=0;i<sz_pt;i++)
		{
			qu.push((ponts) {shr_pt[i],st_pt[i]});
		}

		while(!qu.empty())
		{
			ponts t=qu.top();
			qu.pop();
			fprintf(tbl,"%lf %lf\n",t.first,t.second);
		}

		fclose(bl);
		fclose(tbl);
		fclose(pth);

		chck--;
		revert();
	}

	return 0;
}

void trivial_plotty(int na, int nb)
{
	int i,j,t;
	double a=0.02,b,shr,st;	
	for(i=0;(float)(i*a)<1;i++)
	{
		ponts c=values((i*a),na,nb);
		shr=c.first;
		st=c.second;
		shr_pt[sz_pt]=shr;
		st_pt[sz_pt]=st;
		sz_pt++;
		if(i%5==0)
			cout<<i<<" "<<shr<<" and "<<st<<endl;
		if(shr<0.00001)
			break;
	}
}

void plotty(int na, int nb, vector <double> &alphas)
{
	sz_pt=0;
	int i,j;
	double l=0,r=sqrt((double)INF*10),m;
	
	ponts a=values(l,na,nb);
	shr_pt[sz_pt]=a.first;
	st_pt[sz_pt]=a.second;
	sz_pt++;

	a=values(r,na,nb);
	shr_pt[sz_pt]=a.first;
	st_pt[sz_pt]=a.second;
	sz_pt++;

	double shr,st;
	stack <ponts> stk;
	stk.push(make_pair(l,r));
	
	alphas.push_back(l);
	alphas.push_back(r);

	i=0;
	while(!stk.empty())
	{
		i++;
		ponts c,d;
		c=stk.top();
		//cout<<c.first<<" with the value "<<c.second<<endl;
		stk.pop();
		m=compute_mid(c.first,c.second,&shr,&st,na,nb);
		if(m>(-1))
		{
			alphas.push_back(m);
			stk.push(make_pair(m,c.second));
			//cout<<m<<" and other "<<c.second<<endl;
			stk.push(make_pair(c.first,m));
			//cout<<c.first<<" and other "<<m<<endl;
			shr_pt[sz_pt]=shr;
			st_pt[sz_pt]=st;
			sz_pt++;
		}	
	}

	sort(alphas.begin(), alphas.end());
	vector <double> tmp;
	for(i=0;i<alphas.size();i++)
	{
		tmp.push_back(alphas[i]);
		while((i<alphas.size())&&(alphas[i]==alphas[i+1]))
			i++;	
	}
	alphas.clear();
	alphas=tmp;
}

double compute_mid(double l, double r, double* shr, double* st, int na, int nb)
{
	//cout<<"L AND  R ARE "<<l<<" "<<r<<endl; 

	ponts a=values(l,na,nb);
	ponts b=values(r,na,nb);
	double s=mod((b.second-a.second)/(b.first-a.first));
	ponts c=values(s,na,nb);
	*shr=c.first;
	*st=c.second;
	
	//cout<<b.first<< " "<< b.second<< " "<< a.first<< " "<< a.second<<endl;

	if(s>sqrt(INF))
		s=sqrt(INF);

	//cout<<*shr<<" hi "<<*st<< " si "<<" "<<s <<endl;

	
	double x=distance(*shr,*st,a.first,a.second);
	double y=distance(*shr,*st,b.first,b.second);
	//cout<<x <<" hodor "<<y<<endl;

	if((x<0.001)||(y<0.001))
		return -3;

	//cout<<"THE SLOPE IS WONDERFULLY AND GLORIOUSLY "<<s<<endl;

	return s;
}

ponts values(double p, int na, int nb)
{
	//cout<<"FIRST STEP\n";

	//cout<<"FOR THE VALUE "<<p<<endl;

	penalise(origpath,p);
	
	int S,T,i,opv;

	//cout<<"SECOND STEP\n";

	for ( i = 0; i <= size; ++i )
	{
		dista[i] = INF;
		distb[i] = INF;
		visa[i]=visb[i]=0;
	}

	//cout<<"THIRD STEP\n";

	double st=find_dist( na,nb, &opv);
	
	vector <int> path;

	findpath(path,na,opv,nb);

	double shr=comp_share(origpath,path,&st);

	st=(st/origopt);

	shr=(shr/origopt);

	//cout<<shr << " and "<< st << endl;

	return make_pair(shr,st);
}

void dijkstra( vector < vector < edge > > graph, float dist[], int size, int node ) 
{
	//cout<<" REAK\n";
	dist[node] = 0;
	priority_queue < edge > q;
	q.push( ( edge ) { node, 0 } );

	while ( !q.empty() ) {
		edge p = q.top();
		q.pop();
		//cout<<p.end_node<<" : "<<p.weight<<" : "<<dist[p.end_node]<<" AHEM\n";

		for ( int i = 0; i<graph[p.end_node].size(); ++i )
		{
			int u = p.end_node;
			int v = graph[p.end_node][i].end_node;
			
			if ((dist[u]+graph[p.end_node][i].weight)<dist[v])
			{
				dist[ v ] = dist[ u ] + graph[ p.end_node ][ i ].weight;
				q.push( (edge) {v,dist[v]} );
				//q.push( graph[ p.end_node ][ i ] );
			}
		}
	}
}

float biastar_forlength(float forcost[], float revcost[], float dista[], float distb[], int visa[], int visb[], int na, int nb, int* opv )
{
	//cout<<"BREAK\n";
	float ans=INF,a,b,ra=0,rb=0,corr;

	prev[na] = fora[nb] = dista[na] = distb[nb] = 0;
	
	if(na==nb)
		return 0;
	
	//cout<<na<<" and "<<nb<<"\n";
	int iters=0;

	priority_queue < edge > qa;
	priority_queue < edge > qb;
	
	qa.push( (edge) { na, dista[na]} );
	qb.push( (edge) { nb, distb[nb]} );

	int common=0;

	while (((ra+rb)<(ans+revcost[nb]+forcost[nb]))&&(!qa.empty())&&(!qb.empty())&&(!common)) 
	{
		edge p = qa.top();
		int u = p.end_node;
		//cout<<p.end_node<<" : "<<p.weight<<" A\n";

		if(visa[u])
		{
			qa.pop();
			continue;
		}

		ra=p.weight;

		visa[u] = 1;
		if(visb[u])
			{
			common=1;
			a=dista[u]+distb[u];
			if(a<ans)
			{
				ans=a;	
				*opv=u;
			}	
		}	
		
		qa.pop();
		for ( int i = 0; i<graph[u].size(); ++i )
		{
			int v = graph[u][i].end_node;
			corr=forcost[v]-forcost[u];
			double tk=lengths[make_pair(u,v)];
			if((dista[u]+tk+corr)<dista[v])
			{
				dista[ v ] = dista[ u ] + tk + corr ;
				prev[v]=u;
				if(visb[v]==0)
					qa.push((edge) {v,dista[v]});
			}
			if(visb[v]==1)
			{
				a=dista[u]+distb[v] + tk + corr ;
				if(a<ans)
				{
					ans=a;	
					*opv=u;
				}
			}	
		}

		iters++;

		p = qb.top();
		u = p.end_node;
		//cout<<p.end_node<<" : "<<p.weight<<" B\n";

		if(visb[u])
		{
			qb.pop();
			continue;
		}

		rb=p.weight;		

		visb[u] = 1;
		if(visa[u])
			{
			common=1;
			a=dista[u]+distb[u];
			if(a<ans)
			{
				ans=a;	
				*opv=u;
			}	
		}	

		qb.pop();
		for ( int i = 0; i<revgraph[u].size(); ++i )
		{
			int v = revgraph[u][i].end_node;
			corr=revcost[v]-revcost[u];
			double tk=revlengths[make_pair(u,v)];
			if ((distb[u]+tk+corr)<distb[v])
			{
				distb[ v ] = distb[ u ] + tk + corr;
				fora[v]=u;
				if(visa[v]==0)
					qb.push((edge) {v,distb[v]});
			}
			if(visa[v]==1)
			{
				a=distb[u]+dista[v]+tk+corr;
				if(a<ans)
				{
					ans=a;	
					*opv=u;
				}	
			}
		}

	}
	return ans+forcost[na]-forcost[nb];
}

int penalise(vector <int> optpath, float alpha)
{
	int i,j,t,index_nc=0;
	for(i=0;i<optpath.size()-1;i++)
	{
		lengths[make_pair(optpath[i],optpath[i+1])]=(1+alpha)*olen[make_pair(optpath[i],optpath[i+1])];
		revlengths[make_pair(optpath[i+1],optpath[i])]=(1+alpha)*revolen[make_pair(optpath[i+1],optpath[i])];
	}		
	return 0;
}

int findpath(vector <int> &path, int na, int nc, int nb)
{
	int i,j=nc,t,index_nc=0;
	while(j!=0)
	{
		path.push_back(j);
		j=prev[j];
		index_nc++;
	}
	reverse(path.begin(), path.end());
	j=fora[nc];
	while(j!=0)
	{
		path.push_back(j);
		j=fora[j];
	}
	return index_nc-1;
}

double find_dist(int na, int nb, int* opv)
{	
	int cla = na, clb = nb, i;
	
	float ans;

	float* tmp_dista=(float*)malloc((N+10)*sizeof(float));
	float* tmp_distb=(float*)malloc((N+10)*sizeof(float));
	int* tmp_visa=(int*)malloc((N+10)*sizeof(int));
	int* tmp_visb=(int*)malloc((N+10)*sizeof(int));

	for (int i = 0; i <= size; ++i )
	{
		tmp_dista[i] = INF;
		tmp_distb[i] = INF;
		tmp_visa[i]=tmp_visb[i]=0;
		forpi[i]=revpi[i]=0;
	}

	int x,j;

	for(i=0;i<20;i++)
		stat[i]=1;
	
	for(i=1;i<=3;i++)
	{
		float mx=0,a,b,cst=0;
		j=1;
		for(x=1;x<=num_landmarks;x++)
		{
			a=max(landmarks_rev[x][na]-landmarks_rev[x][nb],landmarks_for[x][nb]-landmarks_for[x][na]);
			b=max(landmarks_for[x][nb]-landmarks_for[x][na],landmarks_rev[x][na]-landmarks_rev[x][nb]);
			cst=a+b;	
			if((cst>mx)&&stat[x])
			{
				j=i;
				mx=cst;
			}
		}	
		stat[j]=0;
		for(x=0;x<=size;x++)
		{
			forpi[x]=max(forpi[x],landmarks_for[j][nb]-landmarks_for[j][x]);
			forpi[x]=max(forpi[x],landmarks_rev[j][x]-landmarks_rev[j][nb]);
			revpi[x]=max(revpi[x],landmarks_for[j][x]-landmarks_for[j][na]);
			revpi[x]=max(revpi[x],landmarks_rev[j][na]-landmarks_rev[j][x]);
		}
	}

	for(x=0;x<=size;x++)
	{
		forcost[x]=((forpi[x]-revpi[x]+revpi[nb])/2);
		revcost[x]=((revpi[x]-forpi[x]+forpi[na])/2);
	}		
	
	ans=biastar_forlength( forcost, revcost, tmp_dista, tmp_distb, tmp_visa, tmp_visb, cla, clb, opv);
	free(tmp_dista);
	free(tmp_distb);
	free(tmp_visb);
	free(tmp_visa);

	return ans;
} 

double comp_share(vector <int> optpath, vector <int> path, double* st)
{
	map<epair,int> e1, e2;
	double sum=0;
	int i;
	*st=0;
	for(i=0;i<optpath.size()-1;i++)
	{
		epair tmp=make_pair(optpath[i],optpath[i+1]);
		e2.insert(make_pair(tmp,1));
	}
	for(i=0;i<path.size()-1;i++)
	{
		epair tmp=make_pair(path[i],path[i+1]);
		*st+=olen[make_pair(path[i],path[i+1])];
		e1.insert(make_pair(tmp,1));
	}
	for(map<epair,int>::iterator itr=e1.begin(); itr!=e1.end();itr++)
	{
		if(e2.find(itr->first)!=e2.end())
		{
			sum+=olen[itr->first];
		}
	}
		return sum;
}

float bidijkstra(int na, int nb, int *opv)
{
	//cout<<"BREAK\n";
	float ans=INF,a,b,ra=0,rb=0,ma=0,mb=0;

	prev[na] = fora[nb] = dista[na] = distb[nb] = 0;
	
	if(na==nb)
		return 0;

	priority_queue < edge > qa, qb;
	
	qa.push( (edge) { na, 0 } );
	qb.push( (edge) { nb, 0 } );

	int iters=0;

	//Possible stopping conditions:
	//((!qa.empty())&&(!qb.empty())&&(((ma<(limit*ans))&&(mb<(limit*ans)))||(((ra+rb)<ans))))
	//((!qa.empty())&&(!qb.empty())&&(((ra+rb)<ans)))	
	//((!qa.empty())&&(!qb.empty()))

	while ((!qa.empty())&&(!qb.empty())&&((ra+rb)<ans))
	{
		edge p = qa.top();
		int u = p.end_node;
		//cout<<p.end_node<<" : "<<p.weight<<" A\n";

		if(visa[u])
		{	
			qa.pop();
			continue;
		}

		ra=dista[u];

		visa[u] = 1;
		if(visb[u])
			{
			a=dista[u]+distb[u];
			if(a<ans)
			{
				ans=a;	
				*opv=u;
			}
		}	
		
		for ( int i = 0; i<graph[u].size(); ++i )
		{
			int v = graph[u][i].end_node;
			graph[u][i].weight=lengths[make_pair(u,v)];
			if((dista[u]+graph[u][i].weight)<dista[v])
			{
				dista[ v ] = dista[ u ] + graph[u][ i ].weight;
				prev[v]=u;
				a=0;
				if(visb[v])
					{
					a=dista[v]+distb[v];
				}	
				if(a<ans)
					qa.push((edge) {v,dista[v]});
			}
			if(visb[v]==1)
			{
				a=dista[u]+distb[v]+graph[u][i].weight;
				if(a<ans)
				{
					ans=a;	
					*opv=u;
				}
			}
			ma=max(ma,dista[v]);	
		}

		p = qb.top();
		u = p.end_node;
		//cout<<p.end_node<<" : "<<p.weight<<" B\n";
		
		if(visb[u])
		{
			qb.pop();
			continue;
		}

		rb=distb[u];		
		visb[u] = 1;
		if(visa[u])
			{
			a=dista[u]+distb[u];
			if(a<ans)
				ans=a;	
		}	

		for ( int i = 0; i<revgraph[u].size(); ++i )
		{
			int v = revgraph[u][i].end_node;
			revgraph[u][i].weight=revlengths[make_pair(u,v)];
			if ((distb[u]+revgraph[u][i].weight)<distb[v])
			{
				distb[ v ] = distb[ u ] + revgraph[u][ i ].weight;
				fora[v]=u;
				a=0;
				if(visa[v])
					{
					a=dista[v]+distb[v];
				}	
				if(a<ans)
					qb.push((edge) {v,distb[v]});
			}
			if(visa[v]==1)
			{
				a=distb[u]+dista[v]+revgraph[u][i].weight;
				if(a<ans)
				{
					ans=a;	
					*opv=u;
				}
			}
			mb=max(mb,distb[v]);
		}

		iters++;
	}		

	//cout<<"THE TOTAL ITERATIONS WERE: "<<iters<<"\n";
	return ans;
} 

void revert()
{
	for( map<epair, float>::iterator it=lengths.begin();it!=lengths.end(); it++)
	{
		lengths[it->first]=olen[it->first];
	}
	for( map<epair, float>::iterator it=revlengths.begin();it!=revlengths.end(); it++)
	{
		revlengths[it->first]=revolen[it->first];
	}
}