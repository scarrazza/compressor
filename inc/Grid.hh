// Compressor - January 2015
// Author:  Stefano Carrazza
// Contact: stefano.carrazza@mi.infn.it
#pragma once

#include <vector>
using std::vector;

class Grid
{
public:
  Grid();
  ~Grid();
  int size() const { return _x.size(); }
  double at(int i) const { return _x[i]; }
private:
  vector<double> _x;
};
