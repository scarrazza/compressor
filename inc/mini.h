#ifndef MINI_H
#define MINI_H

#include <vector>
#include <string>
#include "LHAPDF/LHAPDF.h"
using std::vector;
using std::string;
#define Q 1.0

class RandomGenerator;

class Mini
{

public:
  Mini(int rep, vector<LHAPDF::PDF*> pdf, int seed);
  ~Mini();
  double iterate(int* index);
  void ComputeEstimators(int,double,double,int,int*,
			 double*,double*,double*,double*,double*);
  double GetCV(int i,int j) { return fCV[i][j]; }
  double GetSD(int i,int j) { return fSD[i][j]; }
  double GetSK(int i,int j) { return fSK[i][j]; }
  double GetKU(int i,int j) { return fKU[i][j]; }
  double GetKO(int i,int j, int l) { return fKO[i][j][l]; }
  vector<double> GetX() { return fX; }

private:
  int fRep;
  int fNMut;
  RandomGenerator *fRg;
  vector<double> fX;
  vector<LHAPDF::PDF*> fPDF;
  vector<double*> fCV;
  vector<double*> fSD;
  vector<double*> fKU;
  vector<double*> fSK;
  vector<double**> fKO;
  vector<int*> fMut;
  vector<int> fPids;
};

#endif // MINI_H
