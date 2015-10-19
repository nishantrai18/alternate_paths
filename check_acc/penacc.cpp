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
#include "random_graph.h"

using namespace std;

#define INF 100000000
#define N 265100

typedef pair< int , int> epair;

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

map<epair, float> lengths, revlengths, olen, revolen;
map<long long int, int> vertices;
map<int, pair<double,double> > coords;

vector <vector <edge> > graph,revgraph,prevtree,foratree;
vector <int> prev,fora,origpath;

float landmarks_for[40][N+10],landmarks_rev[40][N+10],forcost[N+10],revcost[N+10],forpi[N+10],revpi[N+10];
float dista[N+10],distb[N+10],opt,param_share,param_stretch,param_localt;
int visa[N+10],visb[N+10],stat[N+10],sz,size,cnt=0;

double altp=0,gshr=0,totp=0,giters=0,initp=0,gstretch=0,glocalt=0,altp_dist[1100],origopt,avglt=0,avgshr=0,avgp=0,avgstr=0;
		
int num_landmarks=32;

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

double norm(double n)
{
	return n/1000000;
}

void dijkstra( vector < vector < edge > > graph, float dist[], int size, int node );
int check_accuracy(float share,float stretch,float localt);
float biastar_forlength(float forcost[], float revcost[], float dista[], float distb[], int visa[], int visb[], int na, int nb, int *opv );
int penalise(vector <int> path, vector <double> &tdist, int na, int nb, float alpha);
int findpath(vector <int> &path, int na, int nc, int nb);
float chshare_node(vector <int> patha, vector <int> pathb);
int chlocopt(vector <int> path, int index_nc, int na, int nc, int nb, float param_t);
double alt_path(int na, int nb, float alpha);
double findlocalt(vector <int> path, vector <double> tdist, int na, int nb);
double find_dist(int na, int nb, int *opv );
int proper_ttest(vector <int> path, vector <double> tdist, int na, int nb, float param_t);
double comp_share(vector <int> optpath, vector <int> path);
float biastar_forlength_orig(float forcost[], float revcost[], float dista[], float distb[], int visa[], int visb[], int na, int nb, int* opv );
double find_dist_orig(int na, int nb, int *opv );
void revert();

int main()
{
	srand (static_cast <unsigned> (time(0)));

	int m,i,j,x,y,iters;

	char c;
	
	/*for(i=0;i<4;i++)
	{
		c='0';
		while(c!='\n')
			cin>>c;
	}
	*/
	cin>>c;
	
    cout<<"Input Size, edges : ";
    cin>>sz>>m;
    
    size=sz;
    float share,stretch,localt,threshold;
	
	share=0.8,stretch=0.25,localt=0.25;
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
        
        /*t.end_node=vertices[x];
        t.weight=a;
        graph[vertices[y]].push_back(t);
        */

        if(i%500000==0)
        	cout<<i<<endl;

        //epair p1 (x,y);
        //lengths.insert(make_pair(p1,a));
        /*t.end_node=x;
        t.weight=a;
        graph[y].push_back(t);
        epair p2 (y,x);
        lengths.insert(make_pair(p2,a));
    	*/
    }

    cout<<"INPUT DONE, Actual Size is "<<v<<"\n";
    size=sz=v-1;

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
				revlengths[p2]=min(lengths[p2],a);
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

	char* files[]={"NYin"};

	int cot=0;
	float alpha=0.25;

	for ( i = 0; i <= size; ++i )
	{
		dista[i] = INF;
		distb[i] = INF;
		visa[i]=visb[i]=0;
	}

	while(cot>=0)
	{
		f=fopen(files[cot],"r");

		int S,T,chck=1000;

		int topa=chck;
		altp=0,totp=0,giters=0,initp=0,gstretch=0,gshr=0;

		for(i=0;i<1100;i++)
		{
			altp_dist[i]=0;
		}

		clock_t begin = clock();
		
		int fail=0;

		vector <pair < int, double > > distas;

		double maxim=0;

		int csk=0;

		while(chck>0)
		{	
			
			fscanf(f,"%d %d",&S,&T);
			avglt=0,avgshr=0,avgp=0,avgstr=0;
			int na=S,nb=T;

			int opv;
		    
			for ( i = 0; i <= size; ++i )
			{
				dista[i] = INF;
				distb[i] = INF;
				fora[i]=prev[i]=visa[i]=visb[i]=0;
			}

		    origpath.clear();
		    origopt=find_dist_orig( na, nb, &opv);
			
			findpath(origpath,na,opv,nb);

			int cd=0;	
			for(i=0;i<9;i++)
			{
				//cout<<i<<endl;
				double a=alt_path(na,nb,alpha);
				if(a>0.99)
					cd++;
				if(cd>=2)
					break;
			}

			int a=avgp;
			distas.push_back(make_pair(a,origopt));

			if(origopt>maxim)
				maxim=origopt;

			if(a==0)
				fail++;

			if(chck%10==0)
			{
				cout<<chck<<", "<<flush;
				cout<<"AVERAGE LOCAL OPTIMALITY, STRETCH, SHARE, PATHS ARE "<<(avglt/avgp)<<" "<<(avgstr/avgp)<<" "<<(avgshr/avgp)<<" "<<avgp<< endl;
			}

			if(avgp)
			{
				glocalt+=(avglt/avgp);
				gshr+=(avgshr/avgp);
				altp+=avgp;
				gstretch+=(avgstr/avgp);
				csk++;
			}

			revert();

			chck--;
		}

		cout<<endl;

		clock_t end = clock();

		sort(distas.begin(), distas.end(), sort_pred() );
		
		double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

		cout<<"THE TIME ELAPSED IS : "<<elapsed_secs<<" AND THE AVERAGE IS : "<<elapsed_secs/topa<<"\n";
		cout<<"THE AVERAGE ALT PATHS FOUND IS : "<<altp/topa<<"\n";
		cout<<"THE AVERAGE STRETCH IS : "<<gstretch/csk<<"\n";
		cout<<"THE AVERAGE LOCAL OPTIMALITY IS : "<<glocalt/altp<<"\n";
		cout<<"THE SUCCESS PERCENTAGE IS : "<<1-((fail*(1.0))/topa)<<"\n";		
	
		cot--;	
		fclose(f);
	}

    return 0;
}

int check_accuracy(float share,float stretch,float localt)
{
	int S,T,chck=100,i,j;

	int topa=chck;
	altp=0,totp=0,giters=0,initp=0;

	for(i=0;i<1100;i++)
	{
		altp_dist[i]=0;
	}

	clock_t begin = clock();
	
	int fail=0;

	vector <pair < int, double > > distas;

	double maxim=0;

	while(chck>0)
	{	
		S=rand()%size+1;
		T=rand()%size+1;
		
		//cout<<S<<" and "<<T<<"\n";
		if(S==T)
			continue;
		
		int a=alt_path(S,T,0.05);

		distas.push_back(make_pair(a,opt));

		if(opt>maxim)
			maxim=opt;

		if(a==0)
			fail++;

		//if(chck%10==0)
		//	cout<<chck<<" "<<flush;
		chck--;
	}

	cout<<"\n";

	clock_t end = clock();

	sort(distas.begin(), distas.end(), sort_pred() );
	
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

	cout<<"THE TIME ELAPSED IS : "<<elapsed_secs<<" AND THE AVERAGE IS : "<<elapsed_secs/topa<<"\n";
	cout<<"THE AVERAGE ALT PATHS FOUND IS : "<<altp/topa<<"\n";
	cout<<"THE AVERAGE NODES CHECKED IS : "<<totp/topa<<"\n";
	cout<<"THE AVERAGE INITIAL NODES ARE : "<<initp/topa<<"\n";
	cout<<"THE AVERAGE ITERATIONS ARE : "<<giters/topa<<"\n";
	cout<<"THE SUCCESS PERCENTAGE IS : "<<1-((fail*(1.0))/topa)<<"\n";
	cout<<"THE DISTRIBUTION IS :\n";

	for(i=0;i<200;i++)
	{
		cout<<i<<":"<<altp_dist[i]<<"  ";
	}
	cout<<"\n";

	/*cout<<"THE OTHER DISTRIBUTION IS :\n";

	for(i=0;i<distas.size();i++)
	{
		cout<<distas[i].first<<":"<<distas[i].second<<"  ";
	}
	cout<<"\n";
	*/
	/*
	cout<<"YET ANOTHER DISTRIBUTION IS :\n";

	j=0;
	for(i=1;i<=10;i++)
	{
		double tot=0,ct=0;
		while(((j<i*distas.size())&&(distas[j].second<(i*maxim/10))))
		{
			tot=tot+distas[j].first;
			ct=ct+1;
			j++;
		}
		cout<<(int) (i*maxim/10)<<":"<<(tot/ct)<<"  ";
	}
	cout<<"\n";
	*/
	
	cout<<"YET ANOTHER DISTRIBUTION IS :\n";

	j=0;
	for(i=1;i<=10;i++)
	{
		double tot=0,ct=0,dt=0,fl=0;
		while(j<(i*distas.size()/10))
		{
			dt=dt+distas[j].second;
			tot=tot+distas[j].first;
			if(distas[j].first==0)
				fl=fl+1;
			ct=ct+1;
			j++;
		}
		cout<<(int) (dt/ct)<<" :  "<<(tot/ct)<<"  "<<1-(fl/ct)<<endl;
	}
	cout<<"\n";
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
	
	*opv=na;
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

	//cout<<"THE NUMBER FOR A IS : "<<iters<<"\n";
	return ans+forcost[na]-forcost[nb];
}

float biastar_forlength_orig(float forcost[], float revcost[], float dista[], float distb[], int visa[], int visb[], int na, int nb, int* opv )
{
	//cout<<"BREAK\n";
	float ans=INF,a,b,ra=0,rb=0,corr;

	prev[na] = fora[nb] = dista[na] = distb[nb] = 0;
	
	if(na==nb)
		return 0;
	
	*opv=na;

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
			double tk=olen[make_pair(u,v)];
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
			double tk=revolen[make_pair(u,v)];
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

	//cout<<"THE NUMBER FOR A IS : "<<iters<<"\n";
	return ans+forcost[na]-forcost[nb];
}

double alt_path(int na, int nb, float alpha)
{
	int i,j,t,x,opv=0;
	
	float mina=INF,corr=0;

	map<int,int> status;
	vector <int> candidates, optpath, init_cands;

	prevtree.clear();
	foratree.clear();
	fora.clear();
	prev.clear();

	for ( i = 0; i <= size; ++i )
	{
		dista[i] = INF;
		distb[i] = INF;
		fora[i]=prev[i]=visa[i]=visb[i]=0;
	}

	opt=find_dist( na, nb, &opv);

	findpath(optpath,na,opv,nb);
	
	vector <double> tdist;
	penalise(optpath,tdist,na,nb,alpha);

	double minlt=findlocalt(optpath,tdist,na,nb);

	double shr=comp_share(origpath,optpath);

	if((shr/origopt)<0.8)
	{
		avgshr+=(shr/origopt);
		avgp+=1;
		avglt+=minlt;
		avgstr+=(tdist[tdist.size()-1]/origopt);
	}

	optpath.clear();
	tdist.clear();

	return (shr/origopt);
}

int penalise(vector <int> optpath, vector <double> &tdist, int na, int nb, float alpha)
{
	int i,j,t,index_nc=0;
	double sum=0;
	for(i=0;i<optpath.size()-1;i++)
	{
		lengths[make_pair(optpath[i],optpath[i+1])]=(1+alpha)*olen[make_pair(optpath[i],optpath[i+1])];
		revlengths[make_pair(optpath[i+1],optpath[i])]=(1+alpha)*revolen[make_pair(optpath[i+1],optpath[i])];
		tdist.push_back(sum);
		sum+=olen[make_pair(optpath[i],optpath[i+1])];
	}
	tdist.push_back(sum);
		
	return 0;
}

int findpath(vector <int> &path, int na, int nc, int nb)
{
	int i,j=nc,t,index_nc=0,k=10000;
	while((j!=0)&&(k>0))
	{
		path.push_back(j);
		j=prev[j];
		index_nc++;
		k--;
	}
	reverse(path.begin(), path.end());
	j=fora[nc];
	k=10000;
	while((j!=0)&&k)
	{
		path.push_back(j);
		j=fora[j];
		k--;
	}
	if(k==0)
		cout<<"GOLO\n";
	return index_nc-1;
}

float chshare_node(vector <int> patha, vector <int> pathb)
{
	int i=0,j=0,t=0,la=patha.size(),lb=pathb.size();
	sort(patha.begin(), patha.end());
	sort(pathb.begin(), pathb.end());
	
	float shr=0;
	
	while((i<la)&&(j<lb))
	{
		if(patha[i]>patha[j])
			j++;
		else if(patha[i]<patha[j])
			i++;
		else
		{	
			shr=shr+1;	
			i++;
		}	
	}

	return shr;
}

int proper_ttest(vector <int> path, vector <double> tdist, int na, int nb, float param_t)
{	
	int i, j=0 ;
	float da, db,ta,tb;

	j=0;
	for (i=0;i<path.size();i++)
	{
		if(j<i)
			j=i;
		while((j<path.size())&&((tdist[j]-tdist[i])<param_t*opt))
		{
			//cout<<i<<" and "<<j<<" , "<<path[i]<<" and "<<path[j]<<endl;
			j++;
		}
		//cout<<i<<" and "<<j<<" , "<<path[i]<<" and "<<path[j]<<endl;
		if(j>=path.size())
			break;
		if((tdist[j]-tdist[i])>param_t*opt)
			j--;
		double tmp=(tdist[j]-tdist[i]), opta;
		if(tmp==0)
			continue;
		//cout<<"THE QUERY IS FOR "<<path[i]<<" and "<<path[j]<<" with distance "<<(tdist[j]-tdist[i])<<endl;
		int opv;
		opta=find_dist_orig(path[i],path[j],&opv);
		
		/*if(opta!=opta2)
			cout<<opta<<" and "<<opta2<<endl;
		*/
		//cout<<tmp<<" AND THE THE DISTANCE "<<opta<<endl;
		if((tmp-opta)*(tmp-opta)>0.0001)
			return 0;
	}
	return 1;
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
		//cout<<x<<" : "<<forcost[x]<<" , "<<revcost[x]<<"\n";
	}		
	
	ans=biastar_forlength( forcost, revcost, tmp_dista, tmp_distb, tmp_visa, tmp_visb, cla, clb, opv);

	//ans=INF;
	//cout<<"chlocopt"<<" "<<na<<" "<<nc<<" "<<nb<<" "<<ans<<"\n";
	//cout<<da+db<<","<<ans<<" :: "<<"nc is "<<nc<<","<<cla<<" and "<<clb<<"\n";

	free(tmp_dista);
	free(tmp_distb);
	free(tmp_visb);
	free(tmp_visa);

	return ans;
} 

double find_dist_orig(int na, int nb, int* opv)
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
		//cout<<x<<" : "<<forcost[x]<<" , "<<revcost[x]<<"\n";
	}		
	
	ans=biastar_forlength_orig( forcost, revcost, tmp_dista, tmp_distb, tmp_visa, tmp_visb, cla, clb, opv);

	//ans=INF;
	//cout<<"chlocopt"<<" "<<na<<" "<<nc<<" "<<nb<<" "<<ans<<"\n";
	//cout<<da+db<<","<<ans<<" :: "<<"nc is "<<nc<<","<<cla<<" and "<<clb<<"\n";

	free(tmp_dista);
	free(tmp_distb);
	free(tmp_visb);
	free(tmp_visa);

	return ans;
} 

double findlocalt(vector <int> path, vector <double> tdist, int na, int nb)
{
	return 0;

	double l=0,r=1,m;
	int i,j;

	//cout<<l<<" and "<<r<<endl;

	while((r-l)>0.01)
	{
		//cout<<l<<" and "<<r<<endl;
		m=((l+r)/2);
		i=proper_ttest(path,tdist,na,nb,m);
		if(i==0)
		{
			r=m;
		}
		else
		{
			l=m;
		}
	}
	return ((l+r)/2);
}

double comp_share(vector <int> optpath, vector <int> path)
{
	map<epair,int> e1, e2;
	double sum=0;
	int i;
	for(i=0;i<optpath.size()-1;i++)
	{
		epair tmp=make_pair(optpath[i],optpath[i+1]);
		e2.insert(make_pair(tmp,1));
	}
	for(i=0;i<path.size()-1;i++)
	{
		epair tmp=make_pair(path[i],path[i+1]);
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