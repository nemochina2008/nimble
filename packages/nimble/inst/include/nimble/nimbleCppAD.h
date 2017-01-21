/* Definitions only to be included when a nimbleFunction needs CppAD */

#include<vector>

class nimbleCppADinfoClass {
 public:
  std::vector<CppAD::AD<double> > ADindependentVars;
};
