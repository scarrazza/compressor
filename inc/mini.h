#ifndef MINI_H
#define MINI_H

#include <vector>
#include <string>
#include "LHAPDF/LHAPDF.h"
using std::vector;
using std::string;

class RandomGenerator;

class Mini
{

public:
  Mini(int rep, vector<LHAPDF::PDF*> pdf);
  ~Mini();
  double iterate(int* index);
  double ComputeAVG(double,double,int,int*);
  void   Save(int* index, string);

private:
  int fRep;
  int fNMut;
  RandomGenerator *fRg;
  vector<double> fX;
  vector<LHAPDF::PDF*> fPDF;
  vector<double*> fCV;
  vector<int*> fMut;
  vector<int> fPids;
};

#endif // MINI_H
