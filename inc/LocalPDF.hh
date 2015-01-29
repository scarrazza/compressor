// Compressor - January 2015
// Author:  Stefano Carrazza
// Contact: stefano.carrazza@mi.infn.it
#pragma once

#include <vector>
#include "LHAPDF/LHAPDF.h"
using std::vector;
class Grid;

class LocalPDF
{
public:
  LocalPDF(vector<LHAPDF::PDF*> pdf, int nf, Grid* const& x, int Q);
  ~LocalPDF();
  double xfxQ(int const& r, int const& ifl, int const& ix);
  double size() { return _rep; }
private:
  int _rep;
  int _nf;
  int _nx;
  double ***_grid;
};
