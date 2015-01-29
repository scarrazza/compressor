// Compressor - January 2015
// Author:  Stefano Carrazza
// Contact: stefano.carrazza@mi.infn.it

#include "LocalPDF.hh"
#include "Grid.hh"

LocalPDF::LocalPDF(vector<LHAPDF::PDF*> pdf, int nf, Grid* const& x, int Q):
  _rep(pdf.size()),
  _nf(nf),
  _nx(x->size())
{
  _grid = new double**[pdf.size()];
  for (int r = 0; r < (int) pdf.size(); r++)
    {
      _grid[r] = new double*[2*nf+1];
      for (int f = -nf; f <= nf; f++)
	{
	  _grid[r][f+nf] = new double[x->size()];
	  for (int i = 0; i < x->size(); i++)
	    _grid[r][f+nf][i] = pdf[r]->xfxQ(f, x->at(i), Q);
	}
    }
}

LocalPDF::~LocalPDF()
{
  for (int r = 0; r < _rep; r++)
    {
      for (int f = 0; f < 2*_nf+1; f++)
	if (_grid[r][f]) delete[] _grid[r][f];
      if (_grid[r]) delete[] _grid[r];
    }
  delete[] _grid;
}

double LocalPDF::xfxQ(int const& r, int const& ifl, int const& ix)
{
  return _grid[r][ifl+_nf][ix];
}
