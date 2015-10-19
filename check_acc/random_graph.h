#ifndef RANDOM_GRAPH
#define RANDOM_GRAPH

#include<iostream>
#include<algorithm>
#include<cstdio>
#include<vector>
#include<cstdlib>
#include<ctime>

struct edge{
    int end_node;
    float weight;
};

typedef struct edge edge;

int random_graph(std::vector <std::vector <edge> > &graph, int size, float threshold)
{
    int i,j;
    graph.clear();
    graph.resize(size+2);
    for(i=1;i<=size;i++)
    {
        for(j=1;j<=size;j++)
        {
            float prob=(static_cast <float> (rand()) / static_cast <float> (RAND_MAX));
            //std::cout<<prob<<" ";
            if(prob<threshold)
            {
                edge t;
                t.end_node=j;
                t.weight=threshold;
                graph[i].push_back(t);
            }
        }
        //std::cout<<"\n";
    }
    return 0;
}

#endif