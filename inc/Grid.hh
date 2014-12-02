#pragma once

#include <vector>
using std::vector;

class Grid
{
public:
  Grid();
  int size() const { return _x.size(); }
  double at(int i) const { return _x[i]; }
private:
  vector<double> _x;
};
