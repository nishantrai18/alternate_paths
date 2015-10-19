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

#define INF 1000000000
#define N 2100000

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

map<epair, float> lengths, olen;
map<long long int, int> vertices;
map<int, pair<double,double> > coords;

vector <vector <edge> > graph,revgraph,prevtree,foratree;
vector <int> prev,fora;

float landmarks_for[40][N+10],landmarks_rev[40][N+10],forcost[N+10],revcost[N+10],forpi[N+10],revpi[N+10];
float dista[N+10],distb[N+10],opt,param_share,param_stretch,param_localt;
int visa[N+10],visb[N+10],stat[N+10],sz,size,cnt=0;

double altp=0,totp=0,giters=0,initp=0,altp_dist[1100];
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
int check_accuracy(float share,float stretch,float localt, FILE* f);
float biastar_forlength(float forcost[], float revcost[], float dista[], float distb[], int visa[], int visb[], int na, int nb );
float biastar(float forcost[], float revcost[], int na, int nb, int *opv );
int penalise(vector <int> path, int na, int nb, float alpha);
int findpath(vector <int> &path, int na, int nc, int nb);
float chshare_node(vector <int> patha, vector <int> pathb);
int chlocopt(vector <int> path, int index_nc, int na, int nc, int nb, float param_t);
int alt_path(int na, int nb, FILE* f, float alpha);
void pplt(vector <int> path, FILE* f);

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
    /*cout<<"Enter parameters : share, stretch, localt: ";
    cin>>share>>stretch>>localt;

    cout<<"Enter threshold for random graph: ";
    cin>>threshold;
	*/

	/*for(i=0;i<2;i++)
	{
		c='0';
		while(c!='\n')
			cin>>c;
	}
	*/

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
	        /*epair p2 (y,x);
	        lengths.insert(make_pair(p2,a));
    		*/
    	}
    }


	/*for(i=1;i<=size;i++)
    {
    	for(j=0;j<graph[i].size();j++)
    	{
    		if(graph[i][j].end_node>=i)
    			break;
    		printf("%d %d %d\n",i,graph[i][j].end_node, (int) (graph[i][j].weight));
    	}
    }
    	
	
	for(i=1;i<=size;i++)
    {
    	for(j=0;j<graph[i].size();j++)
    	{
    		printf("%d ",graph[i][j].end_node);
    	}
    	cout<<"\n";
    }
	*/

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

    //check_accuracy(share,stretch,localt);

	f=fopen("tcpen0.js","w");
	
	float alpha=0.05;

	for(i=0;i<20;i++)
	{
		cout<<i<<endl;
		alt_path( 70994, 14, f, alpha);
		if(i<19)
			fprintf(f,",");
		fprintf(f,"\n");
		alpha+=0.05;
	}
	fclose(f);

    return 0;
}

void pplt(vector <int> path, FILE* f)
{
	int i,j;
	fprintf(f,"[");
	for(i=0;i<path.size()-1;i++)
	{
		fprintf(f,"[%lf,%lf],",coords[path[i]].first,coords[path[i]].second);
	}
	fprintf(f,"[%lf,%lf]]",coords[path[i]].first,coords[path[i]].second);
}

int check_accuracy(float share,float stretch,float localt, FILE* f)
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
		
		int a=alt_path(S,T,f,0.05);

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

float biastar_forlength(float forcost[], float revcost[], float dista[], float distb[], int visa[], int visb[], int na, int nb )
{
	//cout<<"BREAK\n";
	float ans=INF,a,b,ra=0,rb=0,corr;

	dista[na] = 0;
	distb[nb] = 0;
	
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
				ans=a;	
		}	
		
		qa.pop();
		for ( int i = 0; i<graph[u].size(); ++i )
		{
			int v = graph[u][i].end_node;
			corr=forcost[v]-forcost[u];
			graph[u][i].weight=lengths[make_pair(u,v)];
			if((dista[u]+graph[u][i].weight+corr)<dista[v])
			{
				dista[ v ] = dista[ u ] + graph[u][ i ].weight + corr ;
				if(visb[v]==0)
					qa.push((edge) {v,dista[v]});
			}
			if(visb[v]==1)
			{
				a=dista[u]+distb[v]+graph[u][i].weight + corr ;
				if(a<ans)
					ans=a;	
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
				ans=a;	
		}	

		qb.pop();
		for ( int i = 0; i<revgraph[u].size(); ++i )
		{
			int v = revgraph[u][i].end_node;
			corr=revcost[v]-revcost[u];
			graph[u][i].weight=lengths[make_pair(u,v)];
			if ((distb[u]+revgraph[u][i].weight+corr)<distb[v])
			{
				distb[ v ] = distb[ u ] + revgraph[u][ i ].weight + corr;
				if(visa[v]==0)
					qb.push((edge) {v,distb[v]});
			}
			if(visa[v]==1)
			{
				a=distb[u]+dista[v]+revgraph[u][i].weight+corr;
				if(a<ans)
					ans=a;	
			}
		}

		/*cout<<ra+rb<<":"<<ans+revcost[nb]<<"\n"	;
		
		for(int i=1;i<=size;i++)
		{
			printf("%d:%d, ",i,(int)dista[i]);
		}
		printf("\tThis was dist a\n");
		for(int i=1;i<=size;i++)
		{
			printf("%d:%d, ",i,(int)distb[i]);
		}
		printf("\tThis was dist b\n");


		for(int i=1;i<=size;i++)
		{
			if(visa[i])
			printf("%d, ",i);
		}
		printf("\tThis was visa\n");
		for(int i=1;i<=size;i++)
		{
			if(visb[i])
			printf("%d, ",i);
		}
		printf("\tThis was visb\n");
		*/
	}

	//cout<<"THE NUMBER FOR A IS : "<<iters<<"\n";
	return ans+forcost[na]-forcost[nb];
} 

float biastar(float forcost[], float revcost[], int na, int nb, int *opv )
{
	//cout<<"BREAK\n";
	float ans=INF,a,b,ra=0,rb=0,corr=0,ma=0,mb=0;

	float* tmp_fora=(float*)malloc((N+10)*sizeof(float));
	float* tmp_prev=(float*)malloc((N+10)*sizeof(float));

	prev[na] = fora[nb] = dista[na] = distb[nb] = 0;
		
	if(na==nb)
		return 0;
	
	for(int i=0;i<=sz;i++)
	{
		tmp_prev[i]=tmp_fora[i]=0;
	}

	//cout<<na<<" and "<<nb<<"\n";
	int iters=0;

	priority_queue < edge > qa, qb;
	
	qa.push( (edge) { na, dista[na]} );
	qb.push( (edge) { nb, distb[nb]} );

	cout<<na<<nb<<endl;

	int common=0;

	while (((ra+rb)<(ans+forcost[na]-forcost[nb]))&&(!qa.empty())&&(!qb.empty())&&(!common)) 
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
			graph[u][i].weight=lengths[make_pair(u,v)];
			if((dista[u]+graph[u][i].weight+corr)<dista[v])
			{
				dista[ v ] = dista[ u ] + graph[u][ i ].weight + corr ;
				prev[v]=u;
				tmp_prev[v]=graph[u][i].weight;
				
				/*a=0;
				if(visb[v])
					{
					a=dista[v]+distb[v];
				}	
				if(a<limit*(ans+forcost[na]-forcost[nb]))
					qa.push((edge) {v,dista[v]});
				*/
				if(visb[v]==0)
					qa.push((edge) {v,dista[v]});
				
			}
			if(visb[v]==1)
			{
				a=dista[u]+distb[v]+graph[u][i].weight + corr ;
				if(a<ans)
				{
					ans=a;
					*opv=u;	
				}
			}	
			ma=max(ma,dista[v]);
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
			graph[u][i].weight=lengths[make_pair(u,v)];
			if ((distb[u]+revgraph[u][i].weight+corr)<distb[v])
			{
				distb[ v ] = distb[ u ] + revgraph[u][ i ].weight + corr;
				fora[v]=u;
				tmp_fora[v] = revgraph[u][i].weight;
				
				if(visa[v]==0)
					qb.push((edge) {v,distb[v]});
				
			}
			if(visa[v]==1)
			{
				a=distb[u]+dista[v]+revgraph[u][i].weight+corr;
				if(a<ans)
				{
					ans=a;
					*opv=u;	
				}
			}
			mb=max(mb,distb[v]);
		}
	}

	giters+=iters;

	free(tmp_prev);
	free(tmp_fora);

	//cout<<"THE NUMBER FOR A IS : "<<iters<<"\n";
	return ans+forcost[na]-forcost[nb];
} 

int alt_path(int na, int nb, FILE* f, float alpha)
{
	int i,j,t,x,opv=0;
	
	float mina=INF,corr=0;

	map<int,int> status;
	vector <int> candidates, optpath, init_cands;

	cout<<"COOL\n";
	cout<<na<<" "<<nb<<endl;

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
	
    cout<<"COOL1\n";

    //corr=forcost[na]-forcost[nb];

	opt=biastar(forcost, revcost, na, nb, &opv);
	
	/*for ( i = 0; i <= size; ++i )
	{
		if(visa[i]*visb[i]==1)
		{
			//cout<<i<<" ";
			if((dista[i]+distb[i])<mina)
			{
				mina=dista[i]+distb[i];
				opv=i;
			}
		}
	}*/	

	cout<<"COOL\n";
	findpath(optpath,na,opv,nb);
		
	cout<<opv<<" "<<na<<" "<<nb<<endl;

	/*cout<<"The shortest path is : "<<"\n";
	for(int y=0;y<optpath.size();y++)
	{
		cout<<optpath[y]<<" ";
	}	
	cout<<"COOL\n";
	*/
	penalise(optpath,na,nb,alpha);

	cout<<"The shortest between "<<na<<" and "<<nb<<" is "<< opt<<" or "<<mina<<"\n";
	cout<<"It goes through "<<opv<<"\n";

	cout<<"The shortest path is : "<<"\n";
	for(int y=0;y<optpath.size();y++)
	{
		cout<<optpath[y]<<" ";
	}	
	cout<<"\n";

	pplt(optpath,f);

	/*	for(int i=1;i<=size;i++)
		{
			printf("%d:%d, ",i,(int)dista[i]);
		}
		printf("\tThis was dist a\n");
		for(int i=1;i<=size;i++)
		{
			printf("%d:%d, ",i,(int)distb[i]);
		}
		printf("\tThis was dist b\n");


    cout<<"COOL3\n";
	*/	

	return cnt;
}

int penalise(vector <int> optpath, int na, int nb, float alpha)
{
	int i,j,t,index_nc=0;
	
	for(i=0;i<optpath.size()-1;i++)
	{
		lengths[make_pair(optpath[i],optpath[i+1])]=(1+alpha)*olen[make_pair(optpath[i],optpath[i+1])];
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

int chlocopt(vector <int> path, int index_nc, int na, int nc, int nb, float param_t)
{
	int cla = na, clb = nb, i;
	float da = dista[nc], db = distb[nc],ta,tb;

	//cout<<"The parameters is "<<param_t<<" "<<opt<<" "<<param_t*opt<<"\n";

	int l=0,r=index_nc,m;
	while(l<=r)
	{
		m=((l+r)/2);
		ta=(dista[nc]-dista[path[m]]);
		if(ta>=(param_t*opt))
		{
			l=m+1;
		}
		else
		{
			r=m-1;
		}		
	}
	if(l>index_nc)
		l=index_nc;

	cla = path[l];
	da = (dista[nc]-dista[path[l]]);
	
	//cout<<index_nc-l<<" ";

	
	r=path.size()-1,l=index_nc,m;
	while(l<=r)
	{
		m=((l+r)/2);
		tb=(distb[nc]-distb[path[m]]);
		if(tb>=(param_t*opt))
		{
			r=m-1;
		}
		else
		{
			l=m+1;
		}		
	}
	if(l>=path.size())
		l=path.size()-1;

	clb = path[l];
	db = (distb[nc]-distb[path[l]]);
	
	//cout<<l-index_nc<<" GOGO " <<path.size()<<"\n";

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
	
	ans=biastar_forlength( forcost, revcost, tmp_dista, tmp_distb, tmp_visa, tmp_visb, cla, clb);
		
	/*for (int i = 0; i <= size; ++i )
	{
		tmp_dista[i] = INF;
		tmp_distb[i] = INF;
		tmp_visa[i]=tmp_visb[i]=0;
	}

	ans=bidijkstra_forlength( tmp_dista, tmp_distb, tmp_visa, tmp_visb, cla, clb, da+db );
	*/

	//ans=INF;
	//cout<<"chlocopt"<<" "<<na<<" "<<nc<<" "<<nb<<" "<<ans<<"\n";
	//cout<<da+db<<","<<ans<<" :: "<<"nc is "<<nc<<","<<cla<<" and "<<clb<<"\n";

	free(tmp_dista);
	free(tmp_distb);
	free(tmp_visb);
	free(tmp_visa);

	if(((da+db)-ans)*((da+db)-ans)<0.0001)
		return 1;
	else
		return 0;
}