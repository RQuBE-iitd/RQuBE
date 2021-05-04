#include "include/ARRIVAL/parseRegEx/regexToAutomata.cc"
#include "fsa/fsa.cc"

#ifndef random_H
#define random_H
#include "include/ARRIVAL/random/random.cc"
#endif

#ifndef Graph_H
#define Graph_H
#include "include/ARRIVAL/Graph/Graph.cc"
#endif

#ifdef __MACH__
#include <mach/mach.h>
#endif

#include<bits/stdc++.h>
vector<unordered_map<int,unordered_set<int>>> transition_table;

double supp=0.3;
int topCentNodes=50;
int sample=2000;
double r=0.8;
int dir=0;
int relDest=10;
double decimals=1000;
int p_count=0;
double eta=0.001;
int d=2; 
int tailSet=-1;

template<class A, class B, class C>
class triple{
    public:
        A first;
        B second;
        C third;
        triple(A a,B b,C c){
            first=a;
            second=b;
            third=c;
        }
};

inline int findLabel(Graph *g,int v){
    if(g->nodes[v]->labels.size()==1) return -1;
	return g->nodes[v]->labels[1];
}

void horizonFWD(Graph *g,int src,int horizon,unordered_map<int,set<vector<int>>>& fwd){
    queue<vector<int>> q;
    q.push({src});
    while(!q.empty()){
        auto t = q.front();
        q.pop();
        fwd[t.back()].insert(t);
        if(t.size()>=horizon)
            continue;

        Node *currNode = g->nodes[t.back()];
        for(int i=0; i<currNode->allFwdEdges.size();i++){
            if(find(t.begin(),t.end(),currNode->allFwdEdges[i])!=t.end())  // SIMPLE PATH CHECK
                continue;
            t.push_back(currNode->allFwdEdges[i]);
            q.push(t);
            t.pop_back(); 
        }
    }
}


void horizonBWD(Graph *g,int dest,int horizon,unordered_map<int,set<vector<int>>>& fwd,map<vector<int>,int>& pathStrings){
    set<vector<int>> distPaths;
    queue<pair<vector<int>,unordered_set<int>>> q;
    q.push({{dest},{dest}});

    while(!q.empty()){
        auto t = q.front();
        q.pop();
        if(fwd.find(t.first.back())!=fwd.end()){
            int temp=t.first.back();
            t.first.pop_back();
            for(auto i:fwd[temp]){
                bool flag=0;
                for(auto j:i){
                    if(j!=temp && t.second.find(j)!=t.second.end()){
                        flag=1;
                        break;
                    }
                }
                if(flag)
                    continue;
                copy(t.first.rbegin(),t.first.rend(),back_inserter(i));
                auto it = distPaths.insert(i);
                if(it.second){
                    vector<int> path;
                    for(int j=1;j<i.size();j++){
                        path.push_back(g->nodes[i[j-1]]->fwdEdge[i[j]]);
                        // cout<<path.back()<<" ";
                    }
                //     for(auto j:i)
                //         path.push_back(findLabel(g,j));
                    pathStrings[path]++;
                }
            }
            t.first.push_back(temp);
        }

        if(t.first.size()>=horizon)
            continue;

        Node *currNode = g->nodes[t.first.back()];
        for(int i=0; i<currNode->allBwdEdges.size();i++){
            if(find(t.first.begin(),t.first.end(),currNode->allBwdEdges[i])!=t.first.end())  // SIMPLE PATH CHECK
                continue;
            t.first.push_back(currNode->allBwdEdges[i]);
            t.second.insert(currNode->allBwdEdges[i]);
            q.push(t);
            t.first.pop_back();
            t.second.erase(currNode->allBwdEdges[i]);
        }
    }
}

void initFSA( vector<vector<string>>& input, vector<vector<string>> &train,vector<vector<int>>& validation)
{
    if(validation.empty()){
        transition_table.clear();
        construct(input,1,transition_table);
        tailSet=1;
        return;
    }

    int tail;
    bool flag=1;
	for(tail=maxStringLength;tail>1 && !validation.empty();tail--){
        construct(train,tail,transition_table);
        double acc=0;
        for(auto j:validation)
            if(isAccepted(transition_table,j,0,0))
                acc++;
        if((acc/validation.size())>=0.9){
            flag=0;
            break;
        }
        transition_table.clear();
    }

    transition_table.clear();
    cout<<"k "<<tail<<endl;
    tailSet=tail;
    construct(input,tail,transition_table);
}

void nrwr(Graph* g,unordered_set<int>& vis_nodes,int src,int horizon,double tol=1e-5,int maxiter=1e6){
    vector<double> ranks(g->nodes.size(),0);
    vector<unordered_set<int>> rank_states(g->nodes.size());
    set<pair<int,int>> top50;

    int tempLabel = findLabel(g,src);
    vector<int> startStates;
    int l=0;
    if(present(transition_table[0],tempLabel))
        tr(transition_table[0][tempLabel],j){
            startStates.push_back(*j);
        }

    int counter = 1000;
    auto start = std::chrono::system_clock::now();
    int node=src,iter=0;
    int state=startStates[rand()%startStates.size()];
    double prevnorm=0;
    srand(time(NULL));
    int c=0;
    while(true){
        l++;
        ranks[node]++;
        rank_states[node].insert(state);
        auto pos = top50.find({ranks[node]-1,node}); 
        if(pos!=top50.end()){
            top50.erase(pos);
            top50.insert({ranks[node],node});
            c++;
        }
        else if(top50.size()<=topCentNodes)
            top50.insert({ranks[node],node});
        else if(top50.begin()->first<ranks[node]){
            top50.erase(top50.begin());
            top50.insert({ranks[node],node});
            c=0;
        }
        iter++;
        if(iter%1000==0 || c > counter ){
            vector<double> temp_ranks(g->nodes.size(),0);
            transform(ranks.begin(), ranks.end(), temp_ranks.begin(),[&iter](double i) { return i / iter; });
            double curnorm = (double)sqrt(inner_product(temp_ranks.begin(),temp_ranks.end(),temp_ranks.begin(),0.000000));
            if(abs(curnorm-prevnorm)<tol || iter>maxiter || c>counter){
                break;
            }
            prevnorm = curnorm;
        }
        int prob = rand()%100;    // restart prob
        if(prob>90 || l >= 2*horizon){
            state=startStates[rand()%startStates.size()];
            node=src;
            l=0;
            continue;
        }
        Node *currNode = g->nodes[node];
        vector<pair<int,int>> temp;            //  node, state
        for(int i=0; i < currNode->allFwdEdges.size(); i++){
            int d = currNode->allFwdEdges[i];
            int label = findLabel(g,d);
            if(present(transition_table[state],label))
                tr(transition_table[state][label],j){
                    temp.push_back({d,*j});
                }
        }
        if(temp.size()==0){
            state=startStates[rand()%startStates.size()];
            node=src;
            l=0;
            continue;
        }
        int randomEle = rand()%temp.size();
        node = temp[randomEle].first;
        state = temp[randomEle].second;
    }

    vector<pair<int,double>> topNodes;
    iter-=ranks[src];
    double sum=0;
    for(auto i:top50){
        if(i.second==src)
            continue;
        topNodes.push_back({i.second,i.first/(double)iter});
        sum+=(i.first/(double)iter);
    }
    for(int i=topNodes.size()-2;i>=0;i--)
         topNodes[i].second+=topNodes[i+1].second;
    double dec = topNodes.back().second,mult=1;
    while(dec<0.1)
        dec*=10,mult*=10;

    // Neighbourhood Visit
    int searchDepth=d;
    int count=0;
    for(int ii=topNodes.size()-1;ii>=0;ii--){
        count++;
        int node = topNodes[ii].first;
    	vis_nodes.insert(node);
        queue<triple<int, int, int>> q;  // triple = node,state,length
	    for(int state: rank_states[node])
	    	q.push(triple<int,int,int>(node,state,0));

        while(q.size()!=0){
	        triple<int,int,int> temp = q.front();
	        q.pop();
	        if(temp.third >= searchDepth)   
	            continue;
	        Node *currNode = g->nodes[temp.first];
	        for(int i=0; i< currNode->allFwdEdges.size(); i++){
	            int d = currNode->allFwdEdges[i];
	            int label = findLabel(g,d);
	            if(present(transition_table[temp.second],label)){
	                tr(transition_table[temp.second][label],j){
	                    q.push(triple<int,int,int>(d,*j,temp.third+1));
	                    if(finalStates.find(*j)!= finalStates.end())
	                    	vis_nodes.insert(d);
	                }
	            }
	        }
	    }
        if(ii<(topNodes.size()-1) && (double)(abs(topNodes[ii+1].second-topNodes[ii].second)*mult)<eta) {
            break;
        }
    }
}

void RandomPaths(Graph *g,int dest,int horizon,unordered_map<int,set<vector<int>>>& fwd,map<vector<int>,int>& pathStrings,  set<vector<int>> &distPaths, int maxIter){

    vector<int> t;
    unordered_set<int> s;
    t.push_back(dest);
    s.insert(dest);
    int c=0;
    for(int ii=0; ii< maxIter; ii++){
        if(fwd.count(t.back())){
            int temp=t.back();
            t.pop_back();
            for(auto i:fwd[temp]){
                bool flag=0;
                for(auto j:i){
                    if(j!=temp && s.find(j)!=s.end()){
                        flag=1;
                        break;
                    }
                }
                if(flag)
                    continue;
                copy(t.rbegin(),t.rend(),back_inserter(i));
                auto it = distPaths.insert(i);
                if(it.second){
                    vector<int> path;
                    for(int j=1;j<i.size();j++){
                        path.push_back(g->nodes[i[j-1]]->fwdEdge[i[j]]);
                    }
                    pathStrings[path]++;
                }
            }
            t.push_back(temp);
        }

        if(t.size()>=horizon){
 	    t.clear();
            t.push_back(dest);
            s.clear();
            s.insert(dest);
	}

        Node *currNode = g->nodes[t.back()];
        int sz = currNode->allBwdEdges.size();

        if(sz==0){
            t.clear();
            t.push_back(dest);
            s.clear();
            s.insert(dest);
        }
        else{
            int nextNode = currNode->allBwdEdges[rand()%sz];
            if(find(s.begin(),s.end(),nextNode) != s.end()){  // SIMPLE PATH CHECK
                t.clear();
                t.push_back(dest);
                s.clear();
                s.insert(dest);
            }
            else {
                t.push_back(nextNode);
                s.insert(nextNode);
            }
        }
    }
}


int findAnswerSet(Graph *g,int src,unordered_set<int>& vis_nodes,vector<pair<int,double>> &answerSet,unordered_map<int,set<vector<int>>>& fwd,int relDest,int horizon,int support){
    int t=0,c=0,br=0;
    int minIter=500;
    int currSample = sample;
    currSample = max(currSample, 2* support);
    minIter = max(minIter, support);
    unordered_map<int, set<vector<int>>> paths;
    unordered_map<int, int> acc_count;
    unordered_set<int> prev,curr;
    int iter=0;
    while(true){
        iter++;
        vector<pair<int,double>> recommendedRankings;
        for(auto i:vis_nodes){
            map<vector<int>, int> pathStrings;
            RandomPaths(g,i,horizon,fwd,pathStrings,paths[i],currSample);
            int acceptedCount=0;
            int totalCount=0;
            tr(pathStrings, j){
                if(isAccepted(transition_table,j->first,0,0))
                    acceptedCount += j->second;
            }
            
            acc_count[i] += acceptedCount;
            totalCount = paths[i].size();
            if(acc_count[i] > support ){
                int t = acc_count[i]/(totalCount*1.0) * decimals;
                recommendedRankings.push_back({i,((double)t)/decimals});
            }
        }

        sort(recommendedRankings.begin(),recommendedRankings.end(),[](pair<int,float> a,pair<int,float> b){
            return a.second>b.second;
        });
        if(recommendedRankings.size()>relDest && recommendedRankings[relDest-1].second == recommendedRankings.back().second){
            answerSet = recommendedRankings;
            break;
        }

        vis_nodes = unordered_set<int>();
        br=0;

        int n = recommendedRankings.size();
        if(n <= relDest){
            answerSet = recommendedRankings;
            break;
        }
     
        if((int)(r*n) < relDest){
            for(int i=0;i<n;i++){
                if(i<relDest || recommendedRankings[relDest-1].second == recommendedRankings[i].second)
                    answerSet.pb(recommendedRankings[i]);
                else break;
            }
            return iter;
        }
        else if(recommendedRankings[(int)(r*n)].second == recommendedRankings.back().second){
            for(int i=0;i<n;i++){
                if(i<relDest || recommendedRankings[i].second > recommendedRankings[(int)(r*n)].second)
                    vis_nodes.insert(recommendedRankings[i].first);
                else break;
            }
        }
        else {
            for(int i=0;i<n;i++){
                if(i>(int)(r*n)){
                    if(recommendedRankings[(int)(r*n)].second == recommendedRankings[i].second)
                        vis_nodes.insert(recommendedRankings[i].first);
                    else
                        break;
                }
                else
                    vis_nodes.insert(recommendedRankings[i].first);
            }
        }

        for(int i=0;i<n;i++){
            if(i<relDest || recommendedRankings[relDest-1].second == recommendedRankings[i].second)
                curr.insert(recommendedRankings[i].first);
            else break;
        }

        if(prev.size()==curr.size()){
            bool flag =1;
            for(auto i:curr){
                if(prev.find(i)==prev.end()){
                    flag =0;
                    break;
                }
            }
            if(flag){
                for(int i=0;i<n;i++){
                    if(i<relDest || recommendedRankings[relDest-1].second == recommendedRankings[i].second)
                    answerSet.pb(recommendedRankings[i]);
                }
                break;
            }
        }
        currSample = max((int)(0.8*currSample) , minIter);
        prev=curr;
        curr= unordered_set<int>();
        // cout<<vis_nodes.size()<<endl;
    }
    return iter;

}

int main(int argc, char *argv[])
{
    char *edgeFile = argv[1];    // Edge File
    char *labelFile = argv[2];   // Frequent label file
    int src = atoi(argv[3]);     // Source Node
    int dest = atoi(argv[4]);    // Destination Node
    int horizon = atoi(argv[5]); // Path Length

    supp = atof(argv[6]);        // Support threshold
    sample = atoi(argv[7]);      // z, samples to take in each iteration
    r = atof(argv[8]);           // Retain percentage
    dir =atoi(argv[9]);         // Directed/undirected
    relDest = atoi(argv[10]);    // Top k ranks required
    eta = atof(argv[11]);        // Eta threshold
    d = stoi(argv[12]);          // Number of hops to be taken in NRWR

    string attrFile = "include/files/att.txt";
    float walkLength;
    float numWalks;
    Graph *newG = new Graph(edgeFile, labelFile, &attrFile[0], dir);
    Random *rands = new Random(newG->numEdges,1);

    ofstream out1("out/topk.csv");
    ofstream out2("out/nfa.csv");

    alpha.erase(all(alpha));
    maxStringLength = INT_MIN;
    minStringLength = INT_MAX;
    transition_table.clear();

    horizon++;
    map<vector<int>,int> pathStrings;
    unordered_map<int,set<vector<int>>> fwd;
    horizonFWD(newG,src,(horizon/2),fwd);
    horizonBWD(newG,dest,(horizon/2)+(horizon%2),fwd,pathStrings);

    vector<vector<string>> input,train;
    vector<vector<int>> validation;
    int totalPaths=0,n=pathStrings.size();
    for(auto i:pathStrings){
        vector<string> temp;
        bool flag=1;
        totalPaths += i.second;
        for(auto j:i.first){
            temp.push_back(to_string(j));
            if(alpha.find(to_string(j))==alpha.end())
                flag=0;
            alpha.insert(to_string(j));
        }
        int len=temp.size();
        maxStringLength = max(len,maxStringLength);
        minStringLength = min(len,minStringLength);
        input.push_back(temp);
        if(flag && validation.size()<(n*0.2))               // 80-20 split
            validation.push_back(i.first);
        else train.push_back(temp);
    }
    cout<<"input "<<input.size()<<endl;
    if(input.empty()){
        cout<<"Empty set encountered at "<<src<<" "<<dest<<endl;
        return 0;
    }

    initFSA(input,train,validation);                               // Initialise nfa with 90% validation set acceptance
    int support = supp*totalPaths;
    cout<<"support "<<support<<" "<<totalPaths<<endl;

    // NFA Guided Random Walks with Restarts   
    unordered_set<int> vis_nodes;
    nrwr(newG,vis_nodes,src,horizon);
    cout<<"VISITED NODES: "<<vis_nodes.size()<<endl;

    // Creating the answer set
    vector<pair<int, double>> answerSet;
    findAnswerSet(newG,src,vis_nodes,answerSet,fwd,relDest,horizon,support);
    cout<<"answer found "<<answerSet.size()<<endl;

    vector<pair<string,int>> temp;
    out1<<"SRC-"<<src<<",";
    out1<<endl;

    out1<<"DEST-"<<dest<<",";
    out1<<endl;

    for(auto a : answerSet){
        if(a.first==src || a.first==dest)
            continue;
        out1<<a.first<<","<<a.second<<endl;
    }

    // WRITING NFA
    out2<<"SRC-"<<endl;
    out2<<"DEST-"<<", "<<", "<<"Tail Learned,"<<tailSet<<endl;
    out2<<",";
    for(auto l:alpha)
        out2<<l<<",";
    out2<<endl;
    for(int state=0;state<transition_table.size();state++){
        out2<<state<<",";
        for(auto iter:alpha){
            if(transition_table[state].find(stoi(iter))!=transition_table[state].end())
                for(auto nextState:transition_table[state][stoi(iter)])
                    out2<<nextState<<" ";
            out2<<",";
        }
        if(finalStates.find(state)!=finalStates.end())
            out2<<"FINAL";
        out2<<endl;
    }
    out1.close();
    out2.close();
    return 0;
}
