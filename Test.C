#include <vector>
#include "GFW.h"
using std::vector;
int main() {
  GFW *fGFW = new GFW();
  fGFW->AddRegion("pos",3,3,-1,1,1,1);
  fGFW->CreateRegions();
  vector<int> oneVec = vector<int>{ 1, 2, 4 };
  for(auto i: oneVec) printf("%i, ",i);
  printf("\n");
  return 0;
}
