#pragma once

#include <vector>
using std::vector;

class Grid
{
public:
  Grid();
  ~Grid();
  int size() const { return (int) _x.size(); }
  double at(int i) const { return _x[i]; }
private:
  vector<double> _x;
};
