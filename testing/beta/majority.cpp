#include "majority.h"
extern int numTrees,numTaxas,rooted;
extern tree *T;

tree majorityConsensus(){
  tree ret;
  int sum[numTaxas+5];
  memset(sum,0,sizeof(sum));
  for(int i=0;i<numTrees;++i){
    if(!i){
      ret = T[i];
      for(int j=0;j<ret.cnt;++j) ret.goodLabel[j] = 1;
    }
    else ret = majorityMerge(T[i],ret);
	}
  ret.precompute();
	for(int i=0;i<numTrees;++i){
		for(int j=1;j<=numTaxas;++j) ret.good[j] = (ret.H[j] != -1);
		days(ret,T[i]);
		for(int j=1;j<=numTaxas;++j) sum[j] += ret.good[j];
	}
	for(int i=1;i<=numTaxas;++i){
		if(sum[i] > numTrees/2) ret.good[i] = 1;
		else ret.good[i] = 0;
	}
	//Remove non-majority nodes from T1
	ret = contract(ret);
  return ret;
}

tree majContract(tree X){
	int label[X.cnt+5], anc[X.cnt+5],tmp = 0, ai = 0;
	memset(label,-1,sizeof(label));
	memset(anc,-1,sizeof(anc));
	pair<int,int> S[X.cnt+5];
  tree ret;
  for(int i=0;i<X.cnt;++i) tmp += (X.goodLabel[i] > 0);
  ret = tree(X.N,tmp,rooted);
  tmp = 0;
  for(int i=0;i<X.cnt;++i){
    if(X.goodLabel[i] > 0) label[i] = tmp++, ret.goodLabel[label[i]] = X.goodLabel[i];
    if(X.leaf[i] > 0) ret.leaf[label[i]] = X.leaf[i], ret.idx[X.leaf[i]] = label[i];
  }
  ret.root = label[X.root];
  S[ai++] = pair<int,int>(X.root,-1);
  while(ai){
    pair<int,int> t1 = S[--ai];
    if(t1.second == -1){
      if(X.goodLabel[t1.first] > 0 && anc[t1.first] >= 0)
        ret.G[label[anc[t1.first]]].push_back(label[t1.first]);
      if(X.goodLabel[t1.first] > 0) anc[t1.first] = t1.first;
    }
    ++t1.second;
    if(t1.second < X.G[t1.first].size()){
      S[ai++] = t1;
      S[ai++] = pair<int,int>(X.G[t1.first][t1.second],-1);
      anc[X.G[t1.first][t1.second]] = anc[t1.first];
    }
  }
  for(int i=0;i<ret.cnt;++i)
		for(vector<int>::iterator it=ret.G[i].begin();it!=ret.G[i].end();++it) ret.parent[*it] = i;
	return ret;
}

/*
Given two trees A and B (B is R(i-1) )
-Merge A to B and update the counters.
-Find clusters in B incompatible with A and reduce counters
-Contract B and output it 
NOTE: goodLabel is the special counter I am using
*/

tree majorityMerge(tree A,tree B){
	  pair<int,int> S[max(A.cnt,B.cnt)+5];
  int ai = 0, cnt = 0, POS[numTaxas+5], DEPTH[B.cnt+5], pid[B.cnt+5], newNodes = 0;
  bool can[B.cnt+5]; memset(can,0,sizeof(can));
  memset(A.size,0,(A.cnt+5));
  //Relabelling leaves in A so that each subtree in A can be represented as a consecutive range of leaves
  A.precompute();
  
  for(int i=1;i<=numTaxas;++i) B.label[i] = A.label[i];
  
  vector <int> L[numTaxas+5],R[numTaxas+5],BEFORE[B.cnt+5],AFTER[B.cnt+5];
  S[ai++] = pair<int,int>(B.root,-1);
  //Relabelling leaves of B using leaf labels in A and change B into an ordered tree.
  while(ai){
    pair<int,int> t1 = S[--ai];
    if(t1.second == -1){
      if(B.leaf[t1.first] > 0) B.minL[t1.first] = B.label[B.leaf[t1.first]];
      else B.minL[t1.first] = numTaxas+1;
    }
    ++t1.second;
    if(t1.second < B.G[t1.first].size()){
      S[ai++] = t1;
      S[ai++] = pair<int,int>(B.G[t1.first][t1.second],-1);  
    }
    else{
      if(B.parent[t1.first] >= 0) L[B.minL[t1.first]].push_back(t1.first);
      B.G[t1.first].clear();
      if(B.parent[t1.first] >= 0){
        if(B.minL[t1.first] < B.minL[B.parent[t1.first]]) B.minL[B.parent[t1.first]] = B.minL[t1.first];  
      }
    }
  }
  for(int i=1;i<=numTaxas;L[i].clear(),++i)
    for(vector<int>::iterator it=L[i].begin();it!=L[i].end();++it)
      B.G[B.parent[*it]].push_back(*it);
  //Construct the left and right paths for the leaves in B
  cnt = 0, POS[numTaxas+1] = numTaxas+1, POS[numTaxas+2] = -1, DEPTH[B.root] = 0; 
  S[ai++] = pair<int,int>(B.root,-1);
  while(ai){
    pair<int,int> t1 = S[--ai];
    if(t1.second == -1){
      if(B.leaf[t1.first] > 0)
        B.minL[t1.first] = B.maxL[t1.first] = B.leaf[t1.first], POS[B.leaf[t1.first]] = ++cnt;
      else B.minL[t1.first] = numTaxas+1, B.maxL[t1.first] = numTaxas+2;
    }
    ++t1.second;
    if(t1.second < B.G[t1.first].size()){
      DEPTH[B.G[t1.first][t1.second]] = DEPTH[t1.first] + 1;
      S[ai++] = t1;
      S[ai++] = pair<int,int>(B.G[t1.first][t1.second],-1);
    }
    else{
      L[B.minL[t1.first]].push_back(t1.first);
      R[B.maxL[t1.first]].push_back(t1.first);
      int cur = 0;
      for(vector<int>::iterator it=B.G[t1.first].begin();it!=B.G[t1.first].end();++it){
        BEFORE[t1.first].push_back(0);
        AFTER[t1.first].push_back(0);
      }
      if(B.parent[t1.first] >= 0){
        int p = B.parent[t1.first];
        if(POS[B.minL[t1.first]] < POS[B.minL[p]]) B.minL[p] = B.minL[t1.first];
        if(POS[B.maxL[t1.first]] > POS[B.maxL[p]]) B.maxL[p] = B.maxL[t1.first];
      }
      cur = 0;
      for(vector<int>::iterator it=B.G[t1.first].begin();it!=B.G[t1.first].end();++it)
        pid[*it] = cur++;
    }
  }
  //Querying subtrees with A to see if they are compatible with B
  S[ai++] = pair<int,int>(A.root,-1);
  while(ai){
    pair<int,int> t1 = S[--ai];
    if(t1.second == -1){
      if(A.leaf[t1.first]) A.minL[t1.first] = A.maxL[t1.first] = A.leaf[t1.first];
      else A.minL[t1.first] = numTaxas+1, A.maxL[t1.first] = numTaxas+2;
    }
    ++t1.second;
    if(t1.second < A.G[t1.first].size()){
      S[ai++] = t1;
      S[ai++] = pair<int,int>(A.G[t1.first][t1.second],-1);
    }
    else{
      int a = A.minL[t1.first], b = A.maxL[t1.first];
    	if(POS[b] - POS[a] + 1 != A.size[t1.first]) goto end;
      if(DEPTH[L[a].back()] >= DEPTH[R[b].back()]){
        if(DEPTH[L[a].back()] - DEPTH[B.idx[b]] > 0) goto end;
        a = L[a].back();
        b = R[b][DEPTH[B.idx[b]]-DEPTH[a]];
      }
      else{
        if(DEPTH[R[b].back()] - DEPTH[B.idx[a]] > 0) goto end;
        b = R[b].back();
        a = L[a][DEPTH[B.idx[a]]-DEPTH[b]];
      }
   		if(a == b) can[a] = 1, ++B.goodLabel[a];
      else if(B.parent[a] == B.parent[b]){
        int t1 = pid[a], t2 = pid[b];
        ++newNodes;
        ++BEFORE[B.parent[a]][t1];
        ++AFTER[B.parent[a]][t2];
      }
      end:;
      if(A.parent[t1.first] >= 0){
        int p = A.parent[t1.first];
        if(POS[A.minL[t1.first]] < POS[A.minL[p]]) A.minL[p] = A.minL[t1.first];
        if(POS[A.maxL[t1.first]] > POS[A.maxL[p]]) A.maxL[p] = A.maxL[t1.first];
      }
    }
  }
  for(int i=0;i<B.cnt;++i) B.goodLabel[i] -= (can[i] == 0);
  tree ret;
  if(newNodes > 0){
    ret = tree(numTaxas,B.cnt+newNodes,rooted);
    ret.root = B.root;
    ret.goodLabel[ret.root] = B.goodLabel[B.root];
    ret.parent[ret.root] = -1;
    int cnt = B.cnt,cur;
    for(int a=0;a<B.cnt;++a){
      ret.leaf[a] = B.leaf[a];
      cur = a;
      for(int i=0;i<(int)B.G[a].size();++i){
        while(BEFORE[a][i]){
          ret.parent[cnt] = cur;
          ret.G[cur].push_back(cnt);
          ret.goodLabel[cnt] = 1;
          cur = cnt++;
          --BEFORE[a][i];
        }
        ret.G[cur].push_back(B.G[a][i]);
        ret.goodLabel[B.G[a][i]] = B.goodLabel[B.G[a][i]];
        ret.parent[B.G[a][i]] = cur;
        while(AFTER[a][i]){	
          cur = ret.parent[cur];
          --AFTER[a][i];
        }
      }
    }
    for(int i=1;i<=numTaxas;++i) ret.idx[i] = B.idx[i];
    assert(cnt == ret.cnt);
	}
	else ret = B;
  ret = majContract(ret);
  return ret;
}
