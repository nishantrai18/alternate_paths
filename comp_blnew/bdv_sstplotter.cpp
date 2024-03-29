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

#define INF 100000000
#define N 21000

struct edge{
    int end_node;
    float weight;
};

typedef struct edge edge;


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

map<epair, float> lengths;
map<long long int, int> vertices;

vector <vector <edge> > graph,revgraph,prevtree,foratree;
vector <int> prev,fora;

float landmarks_for[120][N+10],landmarks_rev[120][N+10],forcost[N+10],revcost[N+10],forpi[N+10],revpi[N+10];
float dista[N+10],distb[N+10],opt,param_share,param_stretch,param_localt;
int visa[N+10],visb[N+10],stat[N+10],sz,size,cnt=0;

double altp=0,totp=0,giters=0,initp=0,altp_dist[1100];
int num_landmarks=32;
int ttest[3][3],total_calls=0;

double thest[1100], thes[1100];
int thesz=1;

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

void dijkstra( vector < vector < edge > > graph, float dist[], int size, int node );
float biastar_forlength(float forcost[], float revcost[], float dista[], float distb[], int visa[], int visb[], int na, int nb, float val );
float bidijkstra_forlength(float dista[], float distb[], int visa[], int visb[], int na, int nb, float val  );
float bidijkstra(float limit, int na, int nb, int* opv );
int alt_path(int na, int nb, float share, float stretch, float localt);
int check_validity(vector <int> pathab, vector <int> pathacb, int index_nc, int na, int nc, int nb, float share_amt_prev[], float share_amt_fora[], 
				   float share, float stretch, float localt);
int findpath(vector <int> &path, int na, int nc, int nb);
float chshare(vector <int> patha, vector <int> pathb);
int chlocopt(vector <int> pathacb, int index_nc, int na, int nc, int nb, float param_t);
int chstretch(vector <int> pathacb, int index_nc, int na, int nc, int nb, float stretch);
int findshare(vector <int> optpath, float share_amt_prev[], float share_amt_fora[], int na, int nb);
void dfs(int fg, int n, vector <int> optpath, vector <vector <edge> > arr, float share_amt[]);
void stack_dfs(int fg, int n, vector <int> optpath, vector <vector <edge> > arr, float share_amt[]);
void rem_plateau(vector <int> &partner, map<int,int> &status,int n);
int proper_ttest(vector <int> path, int index_nc, int na, int nc, int nb, float param_t);
double find_dist(int na, int nb);

int main(int argc, char **argv)
{
	srand (static_cast <unsigned> (time(0)));

	int m,i,j,x,y,iters;

	char c;
	cin>>c;

    cout<<"Input Size, edges : ";
    cin>>sz>>m;
    
    size=sz;
    float share,stretch,localt,threshold;
	
	share=0.8,stretch=0.4,localt=0.15;
	threshold=0.05;
    /*cout<<"Enter parameters : share, stretch, localt: ";
    cin>>share>>stretch>>localt;

    cout<<"Enter threshold for random graph: ";
    cin>>threshold;
	*/

	cout<<"The parameters are : share: "<<share<<", stretch: "<<stretch<<", localt: "<<localt<<endl;

    int sum=0,tmp=iters,v=1;

    graph.resize(sz+2);
    revgraph.resize(sz+2);
    prev.resize(sz+2);
    fora.resize(sz+2);

    for(i=0;i<m;i++)
    {
    	cin>>c;
		long long int x,y,mx,my;
    	float a;
        scanf("%lld %lld %f",&x,&y,&a);
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
    size=sz=v;
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
	        	lengths.insert(make_pair(p1,a));
	        else
	        	lengths[p1]=min(lengths[p1],a);
	        /*epair p2 (y,x);
	        if(lengths.find(p2)==lengths.end())
	        	lengths.insert(make_pair(p2,a));
	        else
	        	lengths[p2]=min(lengths[p2],a);
    		*/
    	}
    }

    for(i=0;i<=size;i++)
    {
        epair p1 (0,i);
        lengths.insert(make_pair(p1,0));
        epair p2 (i,0);
        lengths.insert(make_pair(p2,0));
        
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
	
	FILE* fp=fopen(argv[1],"r");							//argument 1 is the file to be read from

	FILE *wr=fopen(argv[2],"w");							//argument 2 is the file to be written into

	int S,T;

	int cnts = 100;
	fprintf(wr,"%d\n",cnts);
		
	for(i=0;i<cnts;i++)
	{
		fscanf(fp,"%d %d",&S,&T);

		thes[0]=INF;
		thesz=1;

		//cout<<S<<" "<<T<<endl;

		//check_accuracy(share,stretch,localt);
		alt_path( S, T, share, stretch, localt);

		int fg=1;
		if(thesz<=1)
			fg=0;

		fprintf(wr,"%d %d %d\n",S,T,fg);
		fprintf(wr,"%d\n",thesz-1);
		for(int i=1;i<thesz;i++)
		{
			fprintf(wr,"%lf %lf\n",thes[i],thest[i]);
		}	
	}
	
	fclose(fp);
	fclose(wr);
	
    return 0;
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

float biastar_forlength(float forcost[], float revcost[], float dista[], float distb[], int visa[], int visb[], int na, int nb, float val )
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

	while (((ra+rb)<(ans+revcost[nb]+forcost[nb]))&&(!qa.empty())&&(!qb.empty())&&(!common)&&((ans+forcost[na]-forcost[nb])>=val)) 
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

float bidijkstra_forlength(float dista[], float distb[], int visa[], int visb[], int na, int nb, float val )
{
	float ans=INF,a,b,ra=0,rb=0,ma=0,mb=0;

	dista[na] = distb[nb] = 0;
	
	if(na==nb)
		return 0;
	
	priority_queue < edge > qa, qb;
	
	qa.push( (edge) { na, 0 } );
	qb.push( (edge) { nb, 0 } );

	int iters=0;

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
			if(a<=ans)
				ans=a;	
			else 
			{
				qa.pop();
				continue;
			}
		}	
		
		qa.pop();
		for ( int i = 0; i<graph[u].size(); ++i )
		{
			int v = graph[u][i].end_node;
			if((dista[u]+graph[u][i].weight)<dista[v])
			{
				dista[ v ] = dista[ u ] + graph[u][ i ].weight;
				a=0;
				if(visb[v])
					{
					a=dista[v]+distb[v];
				}	
				if(a<ans)
					qa.push((edge) {v,distb[v]});
			}
			if(visb[v]==1)
			{
				a=dista[u]+distb[v]+graph[u][i].weight;
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
		
		rb=distb[u];

		visb[u] = 1;
		if(visa[u])
			{
			a=dista[u]+distb[u];
			if(a<=ans)
				ans=a;	
			else 
			{
				qb.pop();
				continue;
			}
		}	

		qb.pop();
		for ( int i = 0; i<revgraph[u].size(); ++i )
		{
			int v = revgraph[u][i].end_node;
			if ((distb[u]+revgraph[u][i].weight)<distb[v])
			{
				distb[ v ] = distb[ u ] + revgraph[u][ i ].weight;
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
					ans=a;	
			}
		}
	}

	cout<<"THE NUMBER FOR B IS : "<<iters<<"\n";
	return ans;
} 

float bidijkstra( float limit, int na, int nb, int *opv )
{
	//cout<<"BREAK\n";
	float ans=INF,a,b,ra=0,rb=0,ma=0,mb=0;

	float* tmp_fora=(float*)malloc((N+10)*sizeof(float));
	float* tmp_prev=(float*)malloc((N+10)*sizeof(float));

	prev[na] = fora[nb] = dista[na] = distb[nb] = 0;
	
	if(na==nb)
		return 0;
	
	for(int i=0;i<=sz;i++)
	{
		tmp_prev[i]=tmp_fora[i]=0;
	}

	priority_queue < edge > qa, qb;
	
	qa.push( (edge) { na, 0 } );
	qb.push( (edge) { nb, 0 } );

	int iters=0;

	//Possible stopping conditions:
	//((!qa.empty())&&(!qb.empty())&&(((ma<(limit*ans))&&(mb<(limit*ans)))||(((ra+rb)<ans))))
	//((!qa.empty())&&(!qb.empty())&&(((ra+rb)<ans)))	
	//((!qa.empty())&&(!qb.empty()))

	while ((!qa.empty())&&(!qb.empty())&&(((ra+rb)<limit*ans)))
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
			if((dista[u]+graph[u][i].weight)<dista[v])
			{
				dista[ v ] = dista[ u ] + graph[u][ i ].weight;
				prev[v]=u;
				tmp_prev[v]=graph[u][i].weight;
				a=0;
				if(visb[v])
					{
					a=dista[v]+distb[v];
				}	
				if(a<limit*ans)
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
			{
				ans=a;	
				*opv=u;
			}	
		}	

		for ( int i = 0; i<revgraph[u].size(); ++i )
		{
			int v = revgraph[u][i].end_node;
			if ((distb[u]+revgraph[u][i].weight)<distb[v])
			{
				distb[ v ] = distb[ u ] + revgraph[u][ i ].weight;
				fora[v]=u;
				tmp_fora[v] = revgraph[u][i].weight;
				a=0;
				if(visa[v])
					{
					a=dista[v]+distb[v];
				}	
				if(a<limit*ans)
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
		
		//cout<<ma<<":"<<mb<<" "<<limit*ans<<"\n"	;
		
		//cout<<ra+rb<<":"<<ans<<"\n"	;
		/*
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

		/*for(int i=1;i<=size;i++)
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
	giters+=iters;

	prevtree.clear();
	foratree.clear();
	prevtree.resize(sz+2);
	foratree.resize(sz+2);

	for(int i=1;i<=sz;i++)
	{
		edge t1 = {i,tmp_prev[i]};
		prevtree[prev[i]].push_back(t1);
		edge t2 = {i,tmp_fora[i]};
		foratree[fora[i]].push_back(t2);
	}

	//cout<<"THE TOTAL ITERATIONS WERE: "<<iters<<"\n";

	free(tmp_prev);
	free(tmp_fora);

	return ans;
} 

int alt_path(int na, int nb, float share, float stretch, float localt)
{
	int i,j,t,opv=0;
	
	float mina=INF;

	map<int,int> status;
	vector <int> candidates, optpath;
	float* share_amt_prev=(float*)malloc((N+10)*sizeof(float));
	float* share_amt_fora=(float*)malloc((N+10)*sizeof(float));
	float* tshare_amt_prev=(float*)malloc((N+10)*sizeof(float));
	float* tshare_amt_fora=(float*)malloc((N+10)*sizeof(float));

    //cout<<"COOL\n";

	for ( i = 0; i <= size; ++i )
	{
		dista[i] = INF;
		distb[i] = INF;
		share_amt_prev[i]=share_amt_fora[i]=0;
		tshare_amt_prev[i]=tshare_amt_fora[i]=0;
		fora[i]=prev[i]=visa[i]=visb[i]=0;
	}

    //cout<<"COOL1\n";

	opt=bidijkstra((1+stretch), na, nb, &opv);

    //cout<<"COOL2\n";

    int tots=0;
	priority_queue < edge > sortcand;

	for ( i = 0; i <= size; ++i )
	{
		if(visa[i]*visb[i]==1)
		{
			sortcand.push((edge) {i,(dista[i]+distb[i])});
			tots++;
			status.insert(make_pair(i,1));
			/*if((dista[i]+distb[i])<mina)
			{
				mina=dista[i]+distb[i];
				opv=i;
			}*/
		}
	}	

	mina=dista[opv]+distb[opv];

	findpath(optpath,na,opv,nb);
	
	findshare(optpath,share_amt_prev,share_amt_fora,na,nb);
	
	int cands=0,cnt=0;

	cout<<"THE TOTAL INITIAL CANDIDATES ARE "<<tots<<"\n";
	
	initp+=tots;
	i=0;
	while(!sortcand.empty())
	{
		j=sortcand.top().end_node;
		sortcand.pop();
		if((status[j]==1)&&((dista[j]+distb[j])<=((1+stretch)*opt)))
		{
			i++;
			vector <int> pathacb;
			int index_nc=findpath(pathacb,na,j,nb);
			vector <int> partner;
			partner.clear();
			rem_plateau(partner,status,j);

			t=0;
			for(int tr=0;(!t)&&(tr<partner.size());tr++)
			{
				cands++;
				j=partner[tr];
				t=check_validity(optpath,pathacb,index_nc,na,j,nb,share_amt_prev,share_amt_fora,share,stretch,localt);
				param_stretch=((dista[j]+distb[j])/opt);
				if(t==1)
				{
					altp+=1;
					cnt++;					
					thest[thesz]=param_stretch;
					thes[thesz]=param_share;
					thesz++;
					//cout<<"The params are : "<<param_localt<<" "<<param_stretch<<" "<<param_share<<"\n";
				}
				//cout<<"The params are : "<<param_localt<<" "<<param_stretch<<" "<<param_share<<"\n";
			}			
		}
		if(cands%10000==0)
		{
			cout<<cands<<" AND "<<i<<"\n";
			cout<<"The params are : "<<param_localt<<" "<<param_stretch<<" "<<param_share<<"\n";
		}

	}
	//cout<<"THE TOTAL CANDIDATES WERE "<<cands<<"\n";
	altp_dist[cnt]+=1;
	
	totp+=cands;

	free(share_amt_prev);
	free(share_amt_fora);
	free(tshare_amt_prev);
	free(tshare_amt_fora);

	//cout<<"GOOD"<<endl;
	return cnt;
}

void rem_plateau(vector <int> &partner, map<int,int> &status,int n)
{
	int i,j,t,tmp=n;
	//cout<<"INITIAL N IS "<<n<<"\n";
	partner.push_back(n);
	while(n)
	{
		t=fora[n];
		int fg=0;
		for(i=0;i<prevtree[n].size();i++)
		{
			if((t==prevtree[n][i].end_node)&&(status.find(t)!=status.end()))
        	{
				//cout<<"PARTNERS IN PLATEAU IS "<<t<<"\n";
				status[t]=0;
				fg=1,n=t;
				partner.push_back(n);
				break;
			}	
		}
		if(fg==0)
			n=0;
	}
	n=tmp;
	while(n)
	{
		t=prev[n];
		int fg=0;
		for(i=0;i<foratree[n].size();i++)
		{
			if((t==foratree[n][i].end_node)&&(status.find(t)!=status.end()))
			{
				//cout<<"PARTNERS IN PLATEAU IS "<<t<<"\n";
				status[t]=0;
				fg=1,n=t;
				partner.push_back(n);
				break;
			}	
		}
		if(fg==0)
			n=0;
	}
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

int check_validity(vector <int> pathab, vector <int> pathacb, int index_nc, int na, int nc, int nb, float share_amt_prev[], float share_amt_fora[], 
				   float share, float stretch, float localt)
{
	int i,j,t;
	float shr,str,loc;
	
	shr=share_amt_fora[nc]+share_amt_prev[nc];

	param_share=(shr/opt);
	
	if((shr/opt)>share)
		shr=0;
	else
		shr=1;

	if(shr==0)
		return 0;

	str=chstretch(pathacb,index_nc,na,nc,nb,stretch);
	//str=1;
	//cout<<"STRETCH: "<<str<<"\n";
	int heur=chlocopt(pathacb,index_nc,na,nc,nb,localt);
	
	loc=heur;

	int prop = 1;

	//proper_ttest( pathacb, index_nc, na, nc, nb, localt);

	ttest[prop][heur]++;
	total_calls++;
	

	//loc=1; 
	//cout<<"LOCAL: "<<loc<<"\n";

	//cout<<"SHARE : "<<shr<<"\n";

	return shr*str*loc;

}

int findshare(vector <int> optpath, float share_amt_prev[], float share_amt_fora[], int na, int nb)
{
	int i,j,t;
	vector <int> sorted_path;
	sorted_path=optpath;
	sort(sorted_path.begin(), sorted_path.end());

	/*
	cout<<"The sorted path is : "<<"\n";
	for(int y=0;y<sorted_path.size();y++)
	{
		cout<<sorted_path[y]<<" ";
	}	
	cout<<"\n";
	*/

	share_amt_prev[0]=0;
	share_amt_fora[0]=0;
	
	cnt=0;
	stack_dfs(1,na,sorted_path,prevtree,share_amt_prev);

	//stack_dfs(1,na,sorted_path,prevtree,share_amt_prev);
	cnt=0;
	stack_dfs(0,nb,sorted_path,foratree,share_amt_fora);
	//stack_dfs(0,nb,sorted_path,foratree,share_amt_fora);
	//cout<<"\n\n";
	return 0;
}

void dfs(int fg, int n, vector <int> optpath, vector <vector <edge> > arr, float share_amt[])
{
	int j,t=1;
	//cout<<n<<" ";
	if(fg)
		share_amt[n]=share_amt[prev[n]];
	else
		share_amt[n]=share_amt[fora[n]];
	if(binary_search(optpath.begin(), optpath.end(), n))
	{
		epair p (0,n);
		if(fg)
			p.first=prev[n];
		else
		{
			p.first=n;
			p.second=fora[n];
		}	
		share_amt[n]+=lengths[p];
	}
	for(j=0;j<arr[n].size();j++)
	{
	    dfs(fg,arr[n][j].end_node,optpath,arr,share_amt);
	}
}

void stack_dfs(int fg, int n, vector <int> optpath, vector <vector <edge> > arr, float share_amt[])
{
	int j,t=1;
	stack <int> stk;
	stk.push(n);
	while(!stk.empty())
	{
		n=stk.top();
		//cout<<n<<" ";
		stk.pop();
		if(fg)
			share_amt[n]=share_amt[prev[n]];
		else
			share_amt[n]=share_amt[fora[n]];
		if(binary_search(optpath.begin(), optpath.end(), n))
		{
			epair p (0,n);
			if(fg)
				p.first=prev[n];
			else
				p.first=fora[n];
			share_amt[n]+=lengths[p];
		}
		for(j=0;j<arr[n].size();j++)
		{
		    stk.push(arr[n][j].end_node);
		}
	}
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
			a=max(landmarks_rev[x][cla]-landmarks_rev[x][clb],landmarks_for[x][clb]-landmarks_for[x][cla]);
			b=max(landmarks_for[x][clb]-landmarks_for[x][cla],landmarks_rev[x][cla]-landmarks_rev[x][clb]);
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
			forpi[x]=max(forpi[x],landmarks_for[j][clb]-landmarks_for[j][x]);
			forpi[x]=max(forpi[x],landmarks_rev[j][x]-landmarks_rev[j][clb]);
			revpi[x]=max(revpi[x],landmarks_for[j][x]-landmarks_for[j][cla]);
			revpi[x]=max(revpi[x],landmarks_rev[j][cla]-landmarks_rev[j][x]);
		}
	}

	for(x=0;x<=size;x++)
	{
		forcost[x]=((forpi[x]-revpi[x]+revpi[clb])/2);
		revcost[x]=((revpi[x]-forpi[x]+forpi[cla])/2);
		//cout<<x<<" : "<<forcost[x]<<" , "<<revcost[x]<<"\n";
	}		
	
	ans=biastar_forlength( forcost, revcost, tmp_dista, tmp_distb, tmp_visa, tmp_visb, cla, clb, da+db );
		
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

/*If a via path Pv has stretch (1 + E) and passes the T-test for T = B · dist(s, t) (with 0 < β < 1), then Pv is a B/B-E -UBS path.*/

int chstretch(vector <int> path, int index_nc, int na, int nc, int nb, float stretch)
{
	//cout<<"chstretch : "<<na<<" "<<nc<<" "<<nb<<"\n";
	int i,j,t;
	float a,b;
	a=((dista[nc]+distb[nc])/opt);
	param_localt=a;
	//cout<<(dista[nc]+distb[nc])<<" "<<opt<<" "<<a<<"\n";
	b=(a-1)+((a-1)/stretch);
	//cout<<"Value is : "<<b<<"\n";
	param_stretch=b;
	return 1;
	//return chlocopt(path,index_nc,na,nc,nb,b);
}

int proper_ttest(vector <int> path, int index_nc, int na, int nc, int nb, float param_t)
{	
	int cla = na, clb = nb, i, j ,rl;
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
	rl=l;
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

	//cout<<rl<<" "<<index_nc<<" "<<na<<" "<<nc<<" "<<nb<<endl;

	j=index_nc+1;
	for (i=rl;i<index_nc;i++)
	{
		while((j<=l)&&((distb[nc]-distb[path[j]]+dista[nc]-dista[path[i]])<param_t))
		{
			j++;
		}
		if(j>l)
			break;
		if((distb[nc]-distb[path[j]]+dista[nc]-dista[path[i]])>param_t)
			j--;
		double tmp=(distb[nc]-distb[path[j]]+dista[nc]-dista[path[i]]), opta;
		opta=find_dist(path[i],path[j]);
		if((tmp-opta)*(tmp-opta)>0.0001)
			return 0;
	}
	return 1;
}

double find_dist(int na, int nb)
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
	
	ans=biastar_forlength( forcost, revcost, tmp_dista, tmp_distb, tmp_visa, tmp_visb, cla, clb, INF );

	//ans=INF;
	//cout<<"chlocopt"<<" "<<na<<" "<<nc<<" "<<nb<<" "<<ans<<"\n";
	//cout<<da+db<<","<<ans<<" :: "<<"nc is "<<nc<<","<<cla<<" and "<<clb<<"\n";

	free(tmp_dista);
	free(tmp_distb);
	free(tmp_visb);
	free(tmp_visa);

	return ans;
} 