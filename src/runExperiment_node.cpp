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
int dir=1;
int relDest=10;
double decimals=1000;
int p_count=0;
double eta =0.001;
int d =2; 

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
                    for(auto j:i)
                        path.push_back(findLabel(g,j));
                    pathStrings[path]++;
                }
            }
            t.first.push_back(temp);
        }

        if(t.first.size()>=horizon)
            continue;

        Node *currNode = g->nodes[t.first.back()];
        for(int i=0; i<currNode->allBwdEdges.size();i++){
            if(find(t.second.begin(),t.second.end(),currNode->allBwdEdges[i])!=t.second.end())  // SIMPLE PATH CHECK
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
        if(prob>90 || l >= horizon){
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
        if(ii<(topNodes.size()-1) && (double)(abs(topNodes[ii+1].second-topNodes[ii].second)*mult)< eta) {
            cout<<"top cen "<<ii<<endl;
            break;
        }
    }
}

void RandomPaths(Graph *g,int dest,int horizon,unordered_map<int,set<vector<int>>>& fwd,map<vector<int>,int>& pathStrings,set<vector<int>> &distPaths, int maxIter){
    vector<int> t;
    unordered_set<int> s;
    t.push_back(dest);
    s.insert(dest);
    int c=0;
    for(int ii=0; ii<maxIter; ii++){
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
                    for(auto j:i)
                        path.push_back(findLabel(g,j));
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
    currSample = max(currSample, support);
    minIter = max(minIter, support);
    unordered_map<int, set<vector<int>>> paths;
    unordered_map<int, int> acc_count;
    unordered_set<int> prev,curr;
    int iter =0;
    while(true){
        iter++;
        vector<pair<int,double>> recommendedRankings;
        for(auto i:vis_nodes){
            map<vector<int>, int> pathStrings;
            RandomPaths(g,i,(horizon/2)+(horizon%2),fwd,pathStrings,paths[i],currSample);
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
     
        // if(recommendedRankings[(int)(r*n)].second == recommendedRankings.back().second) 
        //     c++;
        // else
        //     c=0;
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
    }
    return iter;

}

vector<double> metrics(Graph *g,int src,int dest,vector<pair<int, double>> &answerSet,vector<triple<int,double,int>>& node_accep,unordered_map<int,set<vector<int>>>& fwd,int horizon,int support, int relDest){
    double diff=0,count=0,acc=0,act=0,time=0;
    unordered_map<int,double> m;

    for(auto i:node_accep){
        m[i.first]=i.second;
    }

    int destlabel = findLabel(g,dest);
    int com=0,lowerc=0;

    for(int i=0; i< (int)answerSet.size(); i++){
        if(m.find(answerSet[i].first)==m.end()){
            map<vector<int>,int> pathStrings;
		    horizonBWD(g,answerSet[i].first,(horizon/2)+(horizon%2),fwd,pathStrings);
		    int acceptedCount=0;
		    int totalCount=0;
		    for(auto j :pathStrings){
		        totalCount+= j.second;
		        if(isAccepted(transition_table,j.first,0, 0))
		            acceptedCount+= j.second;
		    }
		    double tempAcc = acceptedCount/(double)totalCount;
            tempAcc = (double)((int)(tempAcc *decimals))/decimals;
		    if(acceptedCount > support){
		    	node_accep.push_back(triple<int,double,int>(answerSet[i].first,tempAcc,acceptedCount));
                m[answerSet[i].first] = tempAcc;
            }
        }
        else com++;
    }


    sort(node_accep.begin(),node_accep.end(),[](triple<int,double,int> a,triple<int,double,int> b){
        return a.second>b.second;
    });

    unordered_set<int> GT;
    vector<double> scores;
    for(int i=0;i<relDest;i++){
        while(i<node_accep.size()-1 && node_accep[i].second == node_accep[i+1].second){
            // cout<<node_accep[i].second<<endl;
            GT.insert(node_accep[i].first);
            scores.push_back(node_accep[i++].second);
        }
        GT.insert(node_accep[i].first);
        scores.push_back(node_accep[i].second);
    }
    cout<<GT.size()<<" GT SIZE "<<answerSet.size()<<endl;

    double tp=0,fp=0;
    double normDiff=0;
    int size=answerSet.size();
    int it=min(relDest,size);
   
    double x_mean =0;
    double y_mean =0;
    double n =0;
    double x,y, x_2=0, y_2=0, xy=0;
    vector<double> Xs, Ys;
    unordered_set<double> confX, confY;
    vector<double> rankX, rankY;

    for(int i=0;i<min((int)answerSet.size(),(int)scores.size());i++)
        normDiff+=abs(m[answerSet[i].first]-scores[i]);
    normDiff = 1-(normDiff/min((int)answerSet.size(),(int)scores.size()));

    for(auto i:answerSet){
        if(GT.find(i.first) != GT.end()){
            tp++;
            Xs.pb(i.second);   //already ranked by our algo
            Ys.pb(m[i.first]);
            confX.insert(i.second);
            confY.insert(m[i.first]);
        }
    }
    if(tp ==0) return {0, 0, 0,(double)GT.size()};
    double pearson =0;
    if(confX.size() == 1 && confY.size() ==1){
         pearson = 1; 
         cout<<"pearson "<<*confX.begin()<<" "<<*confY.begin()<<endl;
    } 
    else{
        n = Xs.size();
        double curr =0;
        unordered_map<double , double> temp;
        for(int i =0; i< n; i++){
            if(temp.count(Xs[i])==0){
                if(i > 0){
                    double prev = Xs[i-1];
                    temp[prev] =  ((curr*(curr+1))/2 - (temp[prev]*(temp[prev]-1))/2)/ (curr - temp[prev]+1);
                }
                temp[Xs[i]]= ++curr;
            }
            else  curr++;
            
        }
        temp[Xs[n-1]] =  ((curr*(curr+1))/2 - (temp[Xs[n-1]]*(temp[Xs[n-1]]-1))/2)/ (curr - temp[Xs[n-1]]+1);
        for(int i =0; i< n; i++){
            rankX.pb(temp[Xs[i]]);
            // cout<<Xs[i]<<" "<<temp[Xs[i]]<<endl;
        }
        
        vector<double> sortedY = Ys;
        sort(sortedY.begin(), sortedY.end(), [](double a , double b){ return a > b;});
        curr =0;
        temp.clear();
        for(int i =0; i< n; i++){
            if(temp.count(sortedY[i])==0){
                if(i > 0){
                    double prev = sortedY[i-1];
                    temp[prev] =  ((curr*(curr+1))/2 - (temp[prev]*(temp[prev]-1))/2)/ (curr - temp[prev]+1);
                }
                temp[sortedY[i]]= ++curr;
            }
            else curr++;
        }
        temp[sortedY[n-1]] =  ((curr*(curr+1))/2 - (temp[sortedY[n-1]]*(temp[sortedY[n-1]]-1))/2) / (curr - temp[sortedY[n-1]]+1);
        for(int i =0; i< n; i++){
            rankY.pb(temp[Ys[i]]);
            // cout<<Ys[i]<<" "<<temp[Ys[i]]<<endl;
        }

        double num=0, den=0;
        den = n*(n*n-1);
        for(int i =0; i< n; i++){
            num += (rankX[i]-rankY[i])* (rankX[i]-rankY[i]);
        }
        pearson = 1 - (6*num/den);
        cout<<"pearson "<<pearson<<endl;
    }
    

    
    double precision = min(1.0,tp/min(relDest,(int)answerSet.size()));
    cout<<"precision  "<<precision<<endl;
    return { precision,normDiff,pearson,(double)GT.size()};
}

int main(int argc, char *argv[])
{
    char *edgeFile = argv[1];
    char *labelFile = argv[2];
    char *attrFile = argv[3];
    string res = argv[4];           // Ground truth directory
    int st = atoi(argv[5]);         // Start pair number from ground turth (created using bbfs)
    int en = atoi(argv[6]);         // End pair number from ground turth (created using bbfs)

    // Algorithm parameters
    supp = atof(argv[7]);
    sample = atoi(argv[8]);
    r = atof(argv[9]);
    dir =atoi(argv[10]);
    relDest = atoi(argv[11]);
    eta = atof(argv[12]);
    d = stoi(argv[13]);
    int horizon = atoi(argv[14]);
    string exp = argv[15];           // Experiment performing
    
    int nn = en-st;
    float walkLength;
    float numWalks;
    Graph *newG = new Graph(edgeFile, labelFile, attrFile, dir);
    Random *rands = new Random(newG->numEdges,1);
    double countZero=0,acc=0,act=0,time1=0,time2=0,time3=0,precision=0,avgSup=0,totalVis=0,avgAnsSet=0,avgGt=0,pearson=0;
    double csim=0;
    double iter =0;
    while(st<en){
        alpha.erase(all(alpha));
        maxStringLength = INT_MIN;
        minStringLength = INT_MAX;

        ifstream fn(res+"/"+to_string(st)+".csv");
        string line,source,destination,node,accep,accepCount,K;
        int maxAccCount=0;
        vector<triple<int,double,int>> node_accep,temp;
        getline(fn, line);
        stringstream s(line);
        getline(s, source, ',');
        getline(s, destination, ',');
        getline(s,K,',');
        getline(s,K,',');
        while(getline(fn,line))
	        if(line !=""){
	        	stringstream s(line);
		        getline(s, node, ',');
		        getline(s, accep, ',');
		        getline(s,accepCount, ',');
		        temp.push_back(triple<int,double,int>(stoi(node),(double)((int)(stod(accep)*decimals))/decimals,stoi(accepCount)));
		        maxAccCount= max(maxAccCount,stoi(accepCount));
	        }
        fn.close();
        
        transition_table.clear();
        int src = stoi(source),dest=stoi(destination);
        // ------------------------Time in NFA Creation--------------------------------------------- 
        // create the unordered_set of alphabets - alpha
        // create the vector<vector<string>> labelSequences - input
        // also set minStringLength and maxStringLength 
        
        auto start = chrono::system_clock::now();
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

        initFSA(input,train,validation);                                // Initialise nfa with 90% validation set acceptance
	    // transition_table.clear();
        // construct(input,stoi(K),transition_table);
        auto end = std::chrono::system_clock::now();
        chrono::duration<double> elapsed_seconds = end-start;
        time1+=elapsed_seconds.count();
        // NFA CREATED

        int support = supp*totalPaths;
        cout<<"support "<<support<<" "<<totalPaths<<endl;
        tr(temp,i){
        	if(i->third > support)
        	   node_accep.push_back(*i);
            // else node_accep.push_back(triple<int,double,int>(i->first,0,i->third));
        }
        cout<<"nodes "<<node_accep.size()<<endl;
        sort(node_accep.begin(), node_accep.end(), [](triple<int, double, int> &a, triple<int, double, int> &b)
            {return a.second > b.second ||
            	(a.second==b.second && 
            	(a.third > b.third || (a.third == b.third && a.first > b.first)));
            });

        // Ground truth initialized

        // NFA Guided Random Walks with Restarts
        vector<double> ranks(newG->nodes.size(),0);
        vector<unordered_set<int>> rank_states(newG->nodes.size());
       
        unordered_set<int> vis_nodes;
        start = chrono::system_clock::now();
        nrwr(newG,vis_nodes,src,horizon);
        end = std::chrono::system_clock::now();
        totalVis += vis_nodes.size();
        elapsed_seconds = end-start;
        time2+=elapsed_seconds.count();
        cout<<"VISITED NODES: "<<vis_nodes.size()<<endl;

        // Creating the answer set
        vector<pair<int, double>> answerSet;
        start = chrono::system_clock::now();
        iter += findAnswerSet(newG,src,vis_nodes,answerSet,fwd,relDest,horizon,support);
        end = std::chrono::system_clock::now();
        avgAnsSet += answerSet.size();
        elapsed_seconds = end-start;
        time3+=elapsed_seconds.count();
        cout<<"answer found "<<answerSet.size()<<endl;
        /////////////////////////////////// metrics ////////////////////////////////////

        if(answerSet.size()!=0){
	        auto t = metrics(newG,src,dest,answerSet,node_accep,fwd,horizon,support, relDest);
            precision+=t[0];
       		csim +=t[1];
        	pearson +=max(0.0, t[2]);
            avgGt+=t[3];
            
        }
        else {
            countZero++;
        }

        ///////////////////////////DONE/////////////////////////

        cout<<st<<" :DONE\n";
        st++;
    }

    string name="";
    if(exp=="supp"){
        if(supp==0.1)
            name="1";
        else if(supp==0.3)
            name="3";
        else if(supp==0.5)
            name="5";
        else if(supp==0.7)
            name="7";
        else name="9";
        name+=".supp";
    }
    else if(exp=="retain"){
        if(r==0.6)
            name="6";
        else if(r==0.4)
            name="4";
        else if(r==0.2)
            name="2";
        name+=".retain";
    }
    else if(exp=="k"){
        if(relDest==20)
            name="20";
        else if(relDest==50)
            name="50";
        else if(relDest==100)
            name="100";
        name+=".k";
    }
    else if(exp=="z"){
	if(sample==500)
	    name="0.5";
        if(sample==1000)
            name="1";
        else if(sample==2000)
            name="2";
        else if(sample==3000)
            name="3";
        else if(sample==4000)
            name="4";
        else if(sample==5000)
            name="5";
        name+=".iters";
    }
    else if(exp=="eta"){
        if(eta==0.1)
            name="1";
        else if(eta==0.01)
            name="2";
        else if(eta==0.001)
            name="3";
        else if(eta==0.0001)
            name="4";
        else if(eta==0.00001)
            name="5";
        name+=".eta";
    }
    else if(exp=="d"){
        if(d==0)
            name="0";
        else if(d==1)
            name="1";
        else if(d==2)
            name="2";
        else if(d==3)
            name="3";
        name+=".d";
    }
    cout<<(res+"/"+name)<<endl;
    ofstream myOutFile(res+"/"+name);
    myOutFile<<"Pearson's coefficient "<<pearson/nn<<endl; 
    myOutFile<<"precision "<<precision/nn<<endl;
    myOutFile<<"csim "<<csim/nn<<endl;
    myOutFile<<"total time "<<(time1+time2+time3)/nn<<endl;
    myOutFile<<"time in nfa creation "<<time1/nn<<endl;
    myOutFile<<"time in finding all vis nodes "<<time2/nn<<endl;
    myOutFile<<"Avg answer Set size "<<avgAnsSet/nn<<endl;
    myOutFile<<"Avg gt size "<<avgGt/nn<<endl;  
    myOutFile<<"Avg iter "<<iter/nn<<endl;  
    myOutFile<<"zeros "<<countZero<<endl;
    myOutFile.close();

    cout<<"Pearson's coefficient "<<pearson/nn<<endl; 
    cout<<"precision "<<precision/nn<<endl;
    cout<<"csim "<<csim/nn<<endl;
    cout<<"total time "<<(time1+time2+time3)/nn<<endl;
    cout<<"time in nfa creation "<<time1/nn<<endl;
    cout<<"time in finding all vis nodes "<<time2/nn<<endl;
    cout<<"total time "<<(time1+time2+time3)/nn<<endl;
    cout<<"count (pearson) "<<p_count<<endl;
    cout<<"Avg answer Set size "<<avgAnsSet/nn<<endl;
    cout<<"Avg gt size "<<avgGt/nn<<endl;
    cout<<"Avg iter "<<iter/nn<<endl; 
    cout<<"zeros "<<countZero<<endl; 
    return 0;
}