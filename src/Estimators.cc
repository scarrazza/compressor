// Compressor - January 2015
// Author:  Stefano Carrazza
// Contact: stefano.carrazza@mi.infn.it

#include <iostream>
#include "Estimators.hh"
#include "Grid.hh"
#include "TMatrixD.h"
#include "TDecompSVD.h"
#include "TMatrixDBase.h"
#include "TMatrixDSymEigen.h"
#include "TDecompLU.h"
using namespace std;

double CentralValue::Evaluate(LocalPDF* const& pdf, const int &fl, const vector<int> &index,
                              vector<double> const& w, const int &x) const
{
  double res = 0, wtot = 0;
  const int n = index.size();
  for (int i = 0; i < n; i++)
    {
      res += w[i]*pdf->xfxQ(index[i],fl,x);
      wtot+= w[i];
    }
  return res / wtot;
}

double StdDeviation::Evaluate(LocalPDF* const& pdf, const int &fl, const vector<int> &index,
                              vector<double> const& w, const int &x) const
{
  double sum = 0, sq_sum = 0, wtot = 0, w2tot = 0;

  const int n = index.size();
  for (int i = 0; i < n; i++)
    {
      const double v = pdf->xfxQ(index[i],fl,x);
      sum += w[i]*v;
      sq_sum += w[i]*v*v;
      wtot += w[i];
      w2tot += w[i]*w[i];
    }

  return sqrt( (sq_sum * wtot - sum * sum) / (wtot*wtot - w2tot) );
}

double Skewness::Evaluate(LocalPDF* const& pdf, const int &fl, const vector<int> &index,
                          vector<double> const& w, const int &x) const
{

  CentralValue ecv;
  StdDeviation est;

  double cv = ecv.Evaluate(pdf,fl,index,w,x);
  double st = est.Evaluate(pdf,fl,index,w,x);

  double sum = 0, wtot = 0;
  const int n = index.size();
  for (int i = 0; i < n; i++)
    {
      const double v = pdf->xfxQ(index[i],fl,x) - cv;
      sum += w[i]*v*v*v;
      wtot += w[i];
    }

  return (sum / wtot) / (st*st*st);
}

double Kurtosis::Evaluate(LocalPDF* const &pdf, const int &fl, const vector<int> &index,
                          vector<double> const& w, const int &x) const
{
  CentralValue ecv;
  StdDeviation est;

  double cv = ecv.Evaluate(pdf,fl,index,w,x);
  double st = est.Evaluate(pdf,fl,index,w,x);

  double sum = 0, wtot = 0;
  const int n = index.size();
  for (int i = 0; i < n; i++)
    {
      const double v = pdf->xfxQ(index[i],fl,x) - cv;
      sum += w[i]*v*v*v*v;
      wtot += w[i];
    }

  return (sum / wtot) / (st*st*st*st);
}

double moment5th::Evaluate(LocalPDF* const& pdf, const int &fl, const vector<int> &index,
                           vector<double> const& w, const int &x) const
{
  CentralValue ecv;
  StdDeviation est;

  double cv = ecv.Evaluate(pdf,fl,index,w,x);
  double st = est.Evaluate(pdf,fl,index,w,x);

  double sum = 0, wtot = 0;
  const int n = index.size();
  for (int i = 0; i < n; i++)
    {
      const double v = pdf->xfxQ(index[i],fl,x) - cv;
      sum += w[i]*v*v*v*v*v;
      wtot+= w[i];
    }

  return (sum / wtot) / (st*st*st*st*st);
}

double moment6th::Evaluate(LocalPDF* const &pdf, const int &fl, const vector<int> &index,
                           vector<double> const& w, const int &x) const
{
  CentralValue ecv;
  StdDeviation est;

  double cv = ecv.Evaluate(pdf,fl,index,w,x);
  double st = est.Evaluate(pdf,fl,index,w,x);

  double sum = 0, wtot = 0;
  const int n = index.size();
  for (int i = 0; i < n; i++)
    {
      const double v = pdf->xfxQ(index[i],fl,x) - cv;
      sum += w[i]*v*v*v*v*v*v;
    }

  return (sum / wtot) / (st*st*st*st*st*st);
}

vector<double> Kolmogorov::Evaluate(LocalPDF* const& pdf, int const& fl, vector<int> const& index,
                                    vector<double> const& w, int const& x) const
{
  vector<double> res(_regions,0);

  CentralValue ecv;
  StdDeviation est;

  double cv = ecv.Evaluate(pdf,fl,index,w,x);
  double st = est.Evaluate(pdf,fl,index,w,x);

  double wtot = 0;
  const int n = index.size();
  for (int i = 0; i < n; i++)
    {
      const double v = pdf->xfxQ(index[i],fl,x);
      if (v <= cv -2*st)
        res[0] += w[i];
      else if (v <= cv - st)
        res[1] += w[i];
      else if (v <= cv)
        res[2] += w[i];
      else if (v <= cv + st)
        res[3] += w[i];
      else if (v <= cv + 2*st)
        res[4] += w[i];
      else if (v > cv + 2*st)
        res[5] += w[i];
      else
        cout << "ERROR" << endl;
      wtot += w[i];
    }

  for (int l = 0; l < _regions; l++) res[l] /= wtot;

  return res;
}

TMatrixD Correlation::Evaluate(LocalPDF* const& pdf, const vector<int> &ids, vector<int> const& index,
                               vector<double> const& w, Grid * const &x) const
{
  const int n  = index.size();
  const int nx = _size / (int) ids.size();
  /*
  const int xx[3] = { 
    (int) (1/5.*x->size()), 
    (int) (1/2.*x->size()), 
    (int) (4/5.*x->size())
  };
  */

  const int xx[5] = { 
    (int) (1/5.*x->size()), 
    (int) (2/5.*x->size()), 
    (int) (3/5.*x->size()), 
    (int) (4/5.*x->size()), 
    (int) (5/5.*x->size())-1, 
  };

  TMatrixD m(_size,_size);

  for (size_t fl1 = 0; fl1 < ids.size(); fl1++)
    for (int ix1 = 0; ix1 < nx; ix1++)
      {
        const int i = nx*fl1+ix1;
        for (size_t fl2 = 0; fl2 < ids.size(); fl2++)
          for (int ix2 = 0; ix2 < nx; ix2++)
            {
              const int j = nx*fl2+ix2;
              double ab = 0, a = 0, b = 0, wtot = 0, w2tot = 0;
              double sq_a = 0, sq_b = 0;
              for (int r = 0; r < n; r++)
                {
                  const double v1 = pdf->xfxQ(index[r], ids[fl1], xx[ix1]);
                  const double v2 = pdf->xfxQ(index[r], ids[fl2], xx[ix2]);
                  ab +=  w[r]*v1*v2;
                  a += w[r]*v1;
                  b += w[r]*v2;
                  sq_a += w[r]*v1*v1;
                  sq_b += w[r]*v2*v2;
                  wtot += w[r];
		  w2tot += w[r]*w[r];
                }
	    
              double sig1 = sqrt( (sq_a * wtot - a * a) / ( wtot*wtot - w2tot) );
              double sig2 = sqrt( (sq_b * wtot - b * b) / ( wtot*wtot - w2tot) );

              m(i,j) = wtot*wtot/(wtot*wtot-w2tot)*(ab/wtot - a/wtot*b/wtot)/(sig1*sig2);
            }
      }

  /*
  TMatrixDSym mtm(TMatrixDSym::kAtA,m);
  TMatrixDSymEigen eigen(mtm);
  const TVectorD eigenVal = eigen.GetEigenValues();

  for (int i = 0; i < eigenVal.GetNoElements(); i++)
    res[i] = eigenVal(i);
  */

  /*
  TMatrixD mp(m);
  mp.Invert();
  TMatrixD r = m*mp;

  res.resize(_size,0);
  for (int i = 0; i < _size; i++) res[0] += r(i,i);
  */

  return m;
}
