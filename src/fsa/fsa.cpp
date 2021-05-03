#include<bits/stdc++.h>
using namespace std;

#define sz(a) int((a).size())
#define pb push_back
#define all(c) (c).begin(),(c).end()
#define tr(c,i) for(typeof((c).begin()) i = (c).begin(); i != (c).end(); i++)
#define present(c,x) ((c).find(x) != (c).end())
#define cpresent(c,x) (find(all(c),x) != (c).end())

long int alpha_max;
long int k=1,maxStringLength=-1,minStringLength=INT_MAX;
unordered_map<string,long int> prefix, sets;  //prefix = prefix to state number, sets = set to state number
long int state_index=1;
set<long int> finalStates;
set<string> totalFinalStrings;
set<string> alpha;

struct trieNode {
    unordered_map<string,trieNode *> children;
    vector<pair<long int,string>> states;
    bool EOW;
};

trieNode* newNode(){
    trieNode* node=new trieNode;
    // for(long int i=0;i<alpha_max;i++)
    //     node->children.pb(NULL);
    node->EOW=0;
    return node;
}


void insert(trieNode *root,vector<string> &s){
    for(long int i=0;i<s.size();i++){
        if(root->children.find(s[i])==root->children.end()){
            trieNode *child=newNode();
            if(i==s.size()-1)
                child->EOW=1;
            root->children[s[i]]=child;
        }
        else if(i==s.size()-1)
            root->children[s[i]]->EOW=1;
        root=root->children[s[i]];
    }
}

void print(trieNode *root){
    if(root==NULL)
        return;
    tr(alpha,i){
        if(!present(root->children,*i))
            continue;
        cout<<*i<<" "<<root->children[*i]->EOW<<" ";
        for(long int j=0;j<root->children[*i]->states.size();j++){
            cout<<"("<<root->children[*i]->states[j].first<<","<<root->children[*i]->states[j].second<<") ";
        }
        cout<<endl;
        print(root->children[*i]);
    }
}

void calPrefix(trieNode *root,long int k){
    if(root==NULL)
        return;
    tr(alpha,i){
        if(!present(root->children,*i))
            continue;
        calPrefix(root->children[*i],k);
        for(long int j=0;j<root->children[*i]->states.size();j++){
            if(root->children[*i]->states[j].first+1<=k)
                if(root->children[*i]->states[j].second=="E")
                    root->states.pb(make_pair(root->children[*i]->states[j].first+1,*i));
                else
                    root->states.pb(make_pair(root->children[*i]->states[j].first+1,*i+"."+root->children[*i]->states[j].second)); 
        }
    }
    if(root->EOW==1)
        root->states.pb(make_pair(0,"E"));
    if(root->states.size()==0)
        root->states.pb(make_pair(2*k,"T"));
}

void findMap(trieNode *root, string pre){
    if(root==NULL)
       return;

    tr(root->children,i){
         string state_string ="";
         for(long int j =0; j< i->second->states.size(); j++){
            state_string += i->second->states[j].second;
            
            if(j<i->second->states.size()-1)
                state_string += ".";
         }
         if(state_string!="" && !present(sets,state_string)){
            sets.insert({state_string,state_index++});
        }
        if(state_string!="")
            prefix.insert(make_pair(pre +"." +i->first, sets[state_string]));
        findMap(i->second, pre+ "."+i->first);
    }


}

void initialFinalStates(vector<vector<string>> &input){
    string str="";
    tr(input,i){
        string temp="";
        tr(*i, j)
            temp = temp + "."+*j;
        if(i->size()<=k)
            str=str+"."+temp;
        finalStates.insert(prefix[temp]);
    }
    if(str!="")
        str=str.substr(1,str.length()-1);
    // sets.insert({str,0});
    prefix.insert({"",0});
}

void transitions(vector<unordered_map<string, set<long int>>> &transition_table){
    tr(prefix, i){
        tr(alpha, j){
            if(present(prefix,i->first+ "."+ *j))
                 transition_table[i->second][*j].insert(prefix[i->first+ "."+ *j]);
        } 
    }
}

long int totalStrings(vector<unordered_map<string,set<long int>>> &transition_table,int curLength,int curState,string curString, int maxStrings){
    if(totalFinalStrings.size() > maxStrings)
        return -1;
    if(present(finalStates,curState)){
        totalFinalStrings.insert(curString);
    }

    tr(alpha,i){
        if(present(transition_table[curState],*i)){
            for(auto j = transition_table[curState][*i].begin();j!=transition_table[curState][*i].end();j++){
                if(curLength<maxStringLength){
                    long int temp = totalStrings(transition_table,curLength+1,*j,curString+"."+*i, maxStrings);    
                    if(temp ==-1)
                        return -1;
                }
            }
        }
    }
    return totalFinalStrings.size();
}


bool isAccepted(vector<unordered_map<string,set<long int>>> &transition_table, vector<string> &path, int pos, int curState){    
    if(pos == path.size())
        if present(finalStates,curState) 
            return true;
        else
            return false;
    string s= path[pos];
    bool a  = false;
    if(present(transition_table[curState], s)){
        for(auto j = transition_table[curState][s].begin();j!=transition_table[curState][s].end();j++){
            a = a || isAccepted(transition_table, path,pos+1, *j);
            if(a)
                return a; 
        }
    }
    return false;
}


int main(int argc, char *argv[]){
    string line;
    int numDes = atoi(argv[1]);
    string outfile = "results.csv";
    ifstream myfile("output/0.txt");
    ofstream myOutFile(outfile);
    vector<vector<string>> input;

   while (getline(myfile, line)){
        if(line=="")
            continue;
        vector<string> temp;
        char *cstr = &line[0u];
        char *t = strtok(cstr, " ");
        while(t){
            temp.push_back(string(t));
            alpha.insert(string(t));
            t = strtok(NULL, " ");
        }
        long int len = temp.size();
        maxStringLength = max(len,maxStringLength);
        minStringLength = min(len,minStringLength);
        input.pb(temp);
   
    }

    myfile.close();
    alpha_max = alpha.size(); 

    int k = min(minStringLength+2,  maxStringLength);
    
    for(;k>=1;k--){
        state_index=1;
        trieNode *root=newNode();

        tr(input,i)
            insert(root,*i);

        calPrefix(root,k);
        // print(root);

        findMap(root,"");

        initialFinalStates(input);

        vector<unordered_map<string, set<long int>> > transition_table(state_index , unordered_map<string,set<long int>>(alpha_max)); 
        transitions(transition_table);

        long int stringCount = totalStrings(transition_table,0,0,"", 150 * input.size());
        
        if(stringCount ==-1)
            break;
        myOutFile<<k<<","<<totalFinalStrings.size()/(float)input.size();
        
        for (int j =1; j <= numDes ; j++){
            ifstream myfile("output/"+ to_string(j)+ ".txt");
            int acceptedCount=0;
            int totalCount =0;
            while (getline(myfile, line)){
                if(line=="")
                    continue;
                totalCount++;
                vector<string> temp;
                char *cstr = &line[0u];
                char *t = strtok(cstr, " ");
                while(t){
                    temp.push_back(string(t));
                    t = strtok(NULL, " ");
                }
                acceptedCount += isAccepted(transition_table,temp,0, 0);
            }
            if(totalCount!=0)
                myOutFile<<","<<acceptedCount/(float)totalCount;
            else
                myOutFile<<",0";
            myfile.close();
        }
        myOutFile<<endl;
        // cout<<endl;
        // tr(transition_table,i){
        //     tr(*i,j){
        //         cout<<j->first<<" "<<"(";
        //         tr(j->second, l)
        //             cout<<*l<<",";
        //         cout<<") ";
        //     }
        //     cout<<endl;
        // }

        // tr(prefix,i)
        //     cout<<i->first<<" "<<i->second<<endl;
        // cout<<endl;

        // tr(sets,i)
        //     cout<<i->first<<" "<<i->second<<endl;
        // cout<<endl;

        // tr(finalStates,i)
        //     cout<<*i<<" ";
        // cout<<endl;

        // cout<<totalFinalStrings.size()<<endl;
    //     // tr(input,i)
    //     //     cout<<*i<<" ";
    //     // cout<<endl;

        cout<<k<<" "<<totalFinalStrings.size()<<" "<<input.size()<<endl;   
        cout<<maxStringLength<<" "<<minStringLength<<" "<<alpha_max<<endl;

    //     myfile<<k<<","<<totalFinalStrings.size()<<","<<input.size()<<endl;
        totalFinalStrings.erase(all(totalFinalStrings));
        prefix.erase(all(prefix));
        sets.erase(all(sets));
        finalStates.erase(all(finalStates));
    }
    myOutFile.close();
    return 0;
}