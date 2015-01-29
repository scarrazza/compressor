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

double CentralValue::Evaluate(LocalPDF* const& pdf, const int &fl,
                              const vector<int> &index, const int &x) const
{
  double res = 0;
  const int n = index.size();
  for (int i = 0; i < n; i++)
    res += pdf->xfxQ(index[i],fl,x);
  return res / n;
}

double StdDeviation::Evaluate(LocalPDF* const& pdf, const int &fl,
                              const vector<int> &index, const int &x) const
{
  double sum = 0, sq_sum = 0;

  const int n = index.size();
  for (int i = 0; i < n; i++)
    {
      const double v = pdf->xfxQ(index[i],fl,x);
      sum += v;
      sq_sum += v*v;
    }

  return sqrt(sq_sum / (n-1.0) - n/(n-1.0) * sum/n * sum/n);
}

double Skewness::Evaluate(LocalPDF* const& pdf, const int &fl,
                              const vector<int> &index, const int &x) const
{

  CentralValue ecv;
  StdDeviation est;

  double cv = ecv.Evaluate(pdf,fl,index,x);
  double st = est.Evaluate(pdf,fl,index,x);

  double sum = 0;
  const int n = index.size();
  for (int i = 0; i < n; i++)
    {
      const double v = pdf->xfxQ(index[i],fl,x) - cv;
      sum += v*v*v;
    }

  return (sum / n) / (st*st*st);
}

double Kurtosis::Evaluate(LocalPDF* const &pdf, const int &fl,
                              const vector<int> &index, const int &x) const
{
  CentralValue ecv;
  StdDeviation est;

  double cv = ecv.Evaluate(pdf,fl,index,x);
  double st = est.Evaluate(pdf,fl,index,x);

  double sum = 0;
  const int n = index.size();
  for (int i = 0; i < n; i++)
    {
      const double v = pdf->xfxQ(index[i],fl,x) - cv;
      sum += v*v*v*v;
    }

  return (sum / n) / (st*st*st*st);
}

double moment5th::Evaluate(LocalPDF* const& pdf, const int &fl,
                              const vector<int> &index, const int &x) const
{
  CentralValue ecv;
  StdDeviation est;

  double cv = ecv.Evaluate(pdf,fl,index,x);
  double st = est.Evaluate(pdf,fl,index,x);

  double sum = 0;
  const int n = index.size();
  for (int i = 0; i < n; i++)
    {
      const double v = pdf->xfxQ(index[i],fl,x) - cv;
      sum += v*v*v*v*v;
    }

  return (sum / n) / (st*st*st*st*st);
}

double moment6th::Evaluate(LocalPDF* const &pdf, const int &fl,
                              const vector<int> &index, const int &x) const
{
  CentralValue ecv;
  StdDeviation est;

  double cv = ecv.Evaluate(pdf,fl,index,x);
  double st = est.Evaluate(pdf,fl,index,x);

  double sum = 0;
  const int n = index.size();
  for (int i = 0; i < n; i++)
    {
      const double v = pdf->xfxQ(index[i],fl,x) - cv;
      sum += v*v*v*v*v*v;
    }

  return (sum / n) / (st*st*st*st*st*st);
}

vector<double> Kolmogorov::Evaluate(LocalPDF* const& pdf, int const& fl,
                        vector<int> const& index, int const& x) const
{
  vector<double> res(_regions,0);

  CentralValue ecv;
  StdDeviation est;

  double cv = ecv.Evaluate(pdf,fl,index,x);
  double st = est.Evaluate(pdf,fl,index,x);

  const int n = index.size();
  for (int i = 0; i < n; i++)
    {
      const double v = pdf->xfxQ(index[i],fl,x);
      if (v <= cv -2*st)
        res[0] += 1;
      else if (v <= cv - st)
        res[1] += 1;
      else if (v <= cv)
        res[2] += 1;
      else if (v <= cv + st)
        res[3] += 1;
      else if (v <= cv + 2*st)
        res[4] += 1;
      else if (v > cv + 2*st)
        res[5] += 1;
      else
        cout << "ERROR" << endl;
    }

  for (int l = 0; l < _regions; l++) res[l] /= (double) n;

  return res;
}

TMatrixD Correlation::Evaluate(LocalPDF* const& pdf, const vector<int> &ids, vector<int> const& index, Grid * const &x) const
{
  const int n  = index.size();
  const int nx = _size / (int) ids.size();
  const int xx[3] = { 
    (int) (1/5.*x->size()), 
    (int) (1/2.*x->size()), 
    (int) (4/5.*x->size())
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
              double ab = 0, a = 0, b = 0;
              double sq_a = 0, sq_b = 0;
              for (int r = 0; r < n; r++)
                {
                  const double v1 = pdf->xfxQ(index[r], ids[fl1], xx[ix1]);
                  const double v2 = pdf->xfxQ(index[r], ids[fl2], xx[ix2]);
                  ab +=  v1*v2;
                  a += v1;
                  b += v2;
                  sq_a += v1*v1;
                  sq_b += v2*v2;
                }
              ab /= n;
              a  /= n;
              b  /= n;

              double sig1 = sqrt(sq_a/(n-1.0) - n/(n-1.0)*a*a);
              double sig2 = sqrt(sq_b/(n-1.0) - n/(n-1.0)*b*b);

              m(i,j) = n/(n-1.0)*(ab - a*b)/(sig1*sig2);
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
