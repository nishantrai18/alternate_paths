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
#define N 41000

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

double det(double a,double b,double c,double d)
{
    return ((a*d*(1.0))-(b*c*(1.0)));
}

double area(vector <ponts> arr)
{
    double t=0;
    int i,j,m=arr.size();
    for(j=0;j<m;j++)
    {   
        t+=det(arr[j].first,arr[j].second,arr[(j+1)%m].first,arr[(j+1)%m].second);
    }
    t=mod(t)/2;
    return t;
}

double dista(double x1, double y1, double x2, double y2, double px, double py);
ponts nearest(ponts a, ponts b, ponts arr[], int m, vector <int> &ids);
double dist_nearest(ponts mar[], ponts arr[], int m, vector <int> ids, double &mina); 

int main(int argc, char** argv)
{
    int a,n,m,i,j,x,y,k,iters;

    double mat[4][2];

    for(i=0;i<3;i++)
    {
        for(j=0;j<2;j++)
        {
            mat[i][j]=0;
        }
    }

    //scanf("%d",&k);
    k=100;

    if(*argc <= 3)
    {
        cout<<"NOT ENOUGH ARGUMENTS\n";
        return 0;
    }

    FILE* f=fopen(argv[1],"r");             //Black line
    FILE* fb=fopen(argv[2],"r");            //BDV 
    FILE* fa=fopen(argv[3],"r");            //A star

    fscanf(f,"%d",&k);
    fscanf(fa,"%d",&k);
    fscanf(fb,"%d",&k);
        
    while(k>0)
    {
        int S,T;
        fscanf(f,"%d%d",&S,&T);
        fscanf(fa,"%d%d",&S,&T);
        fscanf(fb,"%d%d",&S,&T);
        
        fscanf(f,"%d",&a);
        fscanf(fa,"%d",&n);
        fscanf(fb,"%d",&m);
        
        ponts bp[1100],rp[1100];
        vector <ponts> nrst;
        vector <int> ids;
        
        vector <ponts> bl, bdv, as;
        bdv.resize(m);
        as.resize(n);
        bl.resize(a);

        for(i=0;i<a;i++)
            fscanf(f,"%lf %lf",&bl[i].first,&bl[i].second);

        for(i=0;i<n;i++)
            fscanf(fa,"%lf %lf",&as[i].first,&as[i].second);

        for(i=0;i<m;i++)
            fscanf(fb,"%lf %lf",&bdv[i].first,&bdv[i].second);

        /*COMPARE THE POINTS AS YOU LIKE IT*/

        ids.clear(),nrst.clear();
        for(i=0;i<m;i++)
        {
            fscanf(fb,"%lf %lf",&rp[i].first,&rp[i].second);
        }
        for(i=a-1;i>0;i--)
        {
            ponts tmp=nearest(bp[i],bp[i-1],rp,m,ids);
            nrst.push_back(tmp);
        }

        mina=100000;

        mat[0][1]+=area(nrst);
        mat[1][1]+=dist_nearest(rp,bp,a,ids,mina);
        mat[2][1]+=mina;
        
        /*cout<<"AR: "<<area(nrst)<<endl;
        cout<<"COST B: "<<dist_nearest(rp,bp,a,ids,mina)<<endl;
        cout<<"MINA: "<<mina<<endl;
        */

        k--;    
    }

    for(i=0;i<3;i++)
    {
        for(j=0;j<2;j++)
        {
            cout<<mat[i][j]<<" ";
        }
        cout<<endl;
    }

    fclose(f);
    fclose(fa);
    fclose(fb);   

    return 0;
}

double dista(double x1, double y1, double x2, double y2, double px, double py)
{
    double dx = x2 - x1;
    double dy = y2 - y1;
    if ((dx == 0) && (dy == 0))
    {
        dx = px - x1;
        dy = py - y1;
        return sqrt(dx * dx + dy * dy);
    }

    double t = ((px - x1) * dx + (py - y1) * dy) / (dx * dx + dy * dy);

    if (t < 0)
    {
        dx = px - x1;
        dy = py - y1;
    }
    else if (t > 1)
    {
        dx = px - x2;
        dy = py - y2;
    }
    else
    {
        dx = px - (x1 + t * dx);
        dy = py - (y1 + t * dy);
    }
    return sqrt(dx * dx + dy * dy);
}

ponts nearest(ponts a, ponts b, ponts arr[], int m, vector <int> &ids)
{
    int i,j,t,o=0;
    ponts ans;
    double min=100000,x,y;
    for(i=0;i<m;i++)
    {
        x=dista(a.first,a.second,b.first,b.second,arr[i].first,arr[i].second);
        if(x<min)
        {
            min=x;
            ans=arr[i];
            o=i;
        }
    }
    ids.push_back(o);
    return ans;
}

double dist_nearest(ponts mar[], ponts arr[], int m, vector <int> ids, double &mina)
{
    int i,j,t;
    ponts ans;
    double mins=100000,an=0,x,y;
    int stat[1100];
    for(i=0;i<m;i++)
    {
        stat[i]=0;
    }
    for(j=0;j<ids.size();j++)   
    {
        mins=100000;
        if(stat[ids[j]]==1)
            continue;
        stat[ids[j]]=1;
        ponts a=mar[ids[j]];
        for(i=0;i<m-1;i++)
        {
            x=dista(arr[i].first,arr[i].second,arr[i+1].first,arr[i+1].second,a.first,a.second);
            if(x<mins)
                mins=x;
        }
        an+=mins;
        mina=min(mins,mina);
        //cout<<"IT IS "<<ids[j]<<" "<<min<<endl;
    }
    return an;
}