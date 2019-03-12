#include <iostream>
#include <cstdio>
#include <string>
#include <cstdlib>
#include <vector>
#include <fstream>
#include "wrapper.h"
#include "tree.h"
using namespace std;

extern tree *T;
extern int numTrees,numTaxas,rooted;

void nexReader(string s){
  //freopen(s.c_str(),"r",stdin); //Read in data
  ifstream fin(s.c_str());
	string S;
  vector <string> trees;
  getline(fin,S);
  getline(fin,S);
	while(S != "translate"){
		if(S == "translate") break;
		getline(fin,S);
	}
  printf("Inputting Taxas ...\n");
  while(1){
    getline(fin,S);
    ++numTaxas;
    if(S[S.length()-1] == ';') break;
  }
  printf("Inputting Trees ...\n");
  while(1){
    getline(fin,S);
    if(S == "END;") break;
    for(unsigned int i=0;i<S.length();++i)
      if(S[i] == '('){
        trees.push_back(S.substr(i,S.length()-i));
        break;
      }
  }
  numTrees = trees.size();
  T = new tree[numTrees];
  for(int i=0;i<numTrees;++i){
		T[i] = tree(numTaxas,trees[i],rooted);
	}
	return;
}
