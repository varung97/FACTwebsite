#include <iostream>
#include <cstdio>
#include <vector>
#include <string>
#include <cstdlib>
#include <ctime>
#include <stack>
#include "tree.h"
#include "wrapper.h"
#include "strict.h"
#include "majority.h"
#include "loose.h"
#include "greedy.h"  
#include "majorityplus.h"
#define NUM_ALGO 5
using namespace std;

int numTrees,numTaxas,rooted, chosenAlgo, startTime, endTime;
tree *T;
tree ans;
string input;
string algoName[NUM_ALGO] = {
	"Strict Consensus", 
	"Majority Consensus", 
	"Loose Consensus", 
	"Greedy Consensus",
	"Majority Plus Consensus"
};

void endProg(string s){
	cout << s << endl;
	exit(0);
}

int main(){
	srand(time(NULL));
	printf("Key in filename: ");
	cin >> input;
	printf("Enter 1 for Rooted Trees and 0 for Unrooted Trees: ");
	scanf("%d",&rooted);
  if(rooted != 1 && rooted != 0) endProg("Invalid Input. Ending Program ... \n");
	nexReader(input);

  printf("Data set has %d trees, and each tree has %d taxas.\n",numTrees,numTaxas);
  
  printf("Enter the algorithm you want to run. Please choose the following options:\n\n");
  for(int i=0;i<NUM_ALGO;++i) printf("Press %d to choose %s\n",i,algoName[i].c_str());
  printf("\nEnter your Choice: ");
  scanf("%d",&chosenAlgo);
  if(!(chosenAlgo >= 0 && chosenAlgo < NUM_ALGO)) endProg("Invalid Option Chosen. Ending Program ... \n");
  printf("You have chosen %s\n\n",algoName[chosenAlgo].c_str());
  startTime = clock();
  switch(chosenAlgo){
		case 0:
			ans = strictConsensus();
			break;
		case 1:
			ans = majorityConsensus();
			break;
		case 2:
			ans = looseConsensus();
			break;
		case 3:
			ans = greedyConsensus();
			break;
		case 4:
			ans = majorityPlusConsensus();
		default:
			break;
	}
	ans.printNex();
	printf("%s tree has %d nodes\n",algoName[chosenAlgo].c_str(),ans.cnt);
	endTime = clock();
	printf("\nProcessing done. Time taken = %.3lf\n",(double)(endTime - startTime)/(double)(CLOCKS_PER_SEC));
	return 0;
}
