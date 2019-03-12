#ifndef MAJORITY_H
#define MAJORITY_H
#include "tree.h"
#include "strict.h"
#include "loose.h"
#include <vector>
#include <list>
#include <cstring>
using namespace std;

tree majorityConsensus();
tree majContract(tree X);
tree majorityMerge(tree A,tree B);  //A is R(i-1), B is the T(i)

#endif
