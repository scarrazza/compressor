
#include <iostream>
#include "Estimators.hh"
#include "Grid.hh"
#include "TMatrixD.h"
#include "TDecompSVD.h"
#include "TMatrixDBase.h"
#include "TMatrixDSymEigen.h"
using namespace std;

double CentralValue::Evaluate(const vector<LHAPDF::PDF *> &pdf, const int &fl,
                              const vector<int> &index, const double &x, const double &Q) const
{
  double res = 0;
  const int n = index.size();
  for (int i = 0; i < n; i++)
    res += pdf[index[i]]->xfxQ(fl,x,Q);
  return res / n;
}

double StdDeviation::Evaluate(const vector<LHAPDF::PDF *> &pdf, const int &fl,
                              const vector<int> &index, const double &x, const double &Q) const
{
  double sum = 0, sq_sum = 0;

  const int n = index.size();
  for (int i = 0; i < n; i++)
    {
      const double v = pdf[index[i]]->xfxQ(fl,x,Q);
      sum += v;
      sq_sum += v*v;
    }

  return sqrt(sq_sum / n - sum/n * sum/n);
}

double Skewness::Evaluate(const vector<LHAPDF::PDF *> &pdf, const int &fl,
                              const vector<int> &index, const double &x, const double &Q) const
{

  CentralValue ecv;
  StdDeviation est;

  double cv = ecv.Evaluate(pdf,fl,index,x,Q);
  double st = est.Evaluate(pdf,fl,index,x,Q);

  double sum = 0;
  const int n = index.size();
  for (int i = 0; i < n; i++)
    {
      const double v = pdf[index[i]]->xfxQ(fl,x,Q) - cv;
      sum += v*v*v;
    }

  return (sum / n) / (st*st*st);
}

double Kurtosis::Evaluate(const vector<LHAPDF::PDF *> &pdf, const int &fl,
                              const vector<int> &index, const double &x, const double &Q) const
{
  CentralValue ecv;
  StdDeviation est;

  double cv = ecv.Evaluate(pdf,fl,index,x,Q);
  double st = est.Evaluate(pdf,fl,index,x,Q);

  double sum = 0;
  const int n = index.size();
  for (int i = 0; i < n; i++)
    {
      const double v = pdf[index[i]]->xfxQ(fl,x,Q) - cv;
      sum += v*v*v*v;
    }

  return (sum / n) / (st*st*st*st);
}

double moment5th::Evaluate(const vector<LHAPDF::PDF *> &pdf, const int &fl,
                              const vector<int> &index, const double &x, const double &Q) const
{
  CentralValue ecv;
  StdDeviation est;

  double cv = ecv.Evaluate(pdf,fl,index,x,Q);
  double st = est.Evaluate(pdf,fl,index,x,Q);

  double sum = 0;
  const int n = index.size();
  for (int i = 0; i < n; i++)
    {
      const double v = pdf[index[i]]->xfxQ(fl,x,Q) - cv;
      sum += v*v*v*v*v;
    }

  return (sum / n) / (st*st*st*st*st);
}

double moment6th::Evaluate(const vector<LHAPDF::PDF *> &pdf, const int &fl,
                              const vector<int> &index, const double &x, const double &Q) const
{
  CentralValue ecv;
  StdDeviation est;

  double cv = ecv.Evaluate(pdf,fl,index,x,Q);
  double st = est.Evaluate(pdf,fl,index,x,Q);

  double sum = 0;
  const int n = index.size();
  for (int i = 0; i < n; i++)
    {
      const double v = pdf[index[i]]->xfxQ(fl,x,Q) - cv;
      sum += v*v*v*v*v*v;
    }

  return (sum / n) / (st*st*st*st*st*st);
}

vector<double> Kolmogorov::Evaluate(vector<LHAPDF::PDF*> const& pdf, int const& fl,
                        vector<int> const& index,double const& x, double const& Q) const
{
  vector<double> res(_regions,0);

  CentralValue ecv;
  StdDeviation est;

  double cv = ecv.Evaluate(pdf,fl,index,x,Q);
  double st = est.Evaluate(pdf,fl,index,x,Q);

  const int n = index.size();
  for (int i = 0; i < n; i++)
    {
      const double v = pdf[index[i]]->xfxQ(fl,x,Q);
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

vector<double> EigCorrelation::Evaluate(vector<LHAPDF::PDF*> const& pdf, const vector<int> &ids,
                        vector<int> const& index, Grid * const &x, double const& Q) const
{
  vector<double> res(_size,0);

  const int n  = index.size();
  const int nx = _size / (int) ids.size();
  const int xx[3] = { 10, (int) (x->size()/2.0), 60};

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
                  const double v1 = pdf[index[r]]->xfxQ(ids[fl1], x->at(xx[ix1]), Q);
                  const double v2 = pdf[index[r]]->xfxQ(ids[fl2], x->at(xx[ix2]), Q);
                  ab +=  v1*v2;
                  a += v1;
                  b += v2;
                  sq_a += v1*v1;
                  sq_b += v2*v2;
                }
              ab /= n;
              a  /= n;
              b  /= n;

              double sig1 = sqrt(sq_a/n - a*a);
              double sig2 = sqrt(sq_b/n - b*b);

              m(i,j) = (ab - a*b)/(sig1*sig2);
            }
      }

  TMatrixDSym mtm(TMatrixDSym::kAtA,m);
  TMatrixDSymEigen eigen(mtm);
  const TVectorD eigenVal = eigen.GetEigenValues();

  for (int i = 0; i < eigenVal.GetNoElements(); i++)
    res[i] = eigenVal(i);

  return res;
}
