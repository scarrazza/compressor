#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <cstdlib>
#include <sys/stat.h>
#include "LHAPDF/LHAPDF.h"
#include "utils.hh"
#include "Minimizer.hh"
#include "Grid.hh"
#include "Estimators.hh"
using namespace std;

void splash();

int main(int argc, char** argv)
{
  int    rep;
  bool   compress = true;
  string priorname;
  double Q = 1.0;
  if (argc > 2) { rep = atoi(argv[1]); priorname.assign(argv[2]); }
  else { 
    cout << "\n usage: ./compressor [REP] [PDF prior name] [energy Q=1] [compress=true]\n" << endl; 
    exit(-1);
  }
  if (argc >= 4) { Q = atof(argv[3]); }
  if (argc == 5) { compress = atoi(argv[4]); }

  splash();

  cout << "* Prior: " << priorname << endl;
  cout << "* Desired compression: " << rep << endl;
  cout << "* Input energy Q = " << Q << endl;
  cout << "* Creating folder " << priorname.c_str() << endl;
  mkdir(priorname.c_str(),0777);
  fstream f;

  // allocate LHAPDF set
  const LHAPDF::PDFSet set(priorname.c_str());
  vector<LHAPDF::PDF*> pdf = set.mkPDFs();

  // allocate common rg for random testing before fit
  RandomGenerator *rg = new RandomGenerator(0,0);
  Grid *x = new Grid();

  cout << "* X grid size: " << x->size() << " points, x=["
       << x->at(0) << ", " << x->at(x->size()-1) << "]" << endl;
  Minimizer min(pdf,x,Q);

  vector<EstimatorsM*> estM = min.GetMomentEstimators();
  vector<EstimatorsS*> estS = min.GetStatEstimators();
  vector<EstimatorsC*> estC = min.GetCorrEstimators();

  // Computing error function for random set
  const int trials = 1000;
  cout << "* Random trials: " << trials << endl;
  vector<int> index;
  double*   estCval = new double[estC[0]->getSize()];
  double**  estMval = new double*[min.GetIDS().size()];
  double*** estSval = new double**[min.GetIDS().size()];
  for (size_t fl = 0; fl < min.GetIDS().size(); fl++) {
      estMval[fl] = new double[x->size()];
      estSval[fl] = new double*[x->size()];
      for (int ix = 0; ix < x->size(); ix++)
        estSval[fl][ix] = new double[estS[0]->getRegions()];
    }

  vector< vector<double> > erfMs;
  vector< vector<double> > erfSs;
  vector< vector<double> > erfCs;
  vector<double> t(trials,0.0);
  for (size_t es = 0; es < estM.size(); es++) erfMs.push_back(t);
  for (size_t es = 0; es < estS.size(); es++) erfSs.push_back(t);
  for (size_t es = 0; es < estC.size(); es++) erfCs.push_back(t);

  for (int t = 0; t < trials; t++)
    {
      cout << "* Processing trial " << t+1 << "/" << trials << "\r";
      cout.flush();

      index.resize(rep,0);
      // generate random replicas, no duplicates
      randomize(pdf.size()-1,rg,index);

      // computing estimators
      for (size_t es = 0; es < estM.size(); es++)
        {
          for (size_t fl = 0; fl < min.GetIDS().size(); fl++)
            for (int ix = 0; ix < x->size(); ix++)
              estMval[fl][ix] = estM[es]->Evaluate(pdf,min.GetIDS()[fl],index,x->at(ix),Q);
          erfMs[es][t] += ERF(min.GetIDS().size(), x->size(), estMval, min.GetPriorMomentEstValues()[es]);
        }

      // computing estimators
      for (size_t es = 0; es < estS.size(); es++)
        {
          for (size_t fl = 0; fl < min.GetIDS().size(); fl++)
            for (int ix = 0; ix < x->size(); ix++)
              {
                vector<double> res = estS[es]->Evaluate(pdf,min.GetIDS()[fl],index,x->at(ix),Q);
                for (int l = 0; l < estS[es]->getRegions(); l++) estSval[fl][ix][l] = res[l];
              }
          erfSs[es][t] += ERFS(min.GetIDS().size(),x->size(),estS[es]->getRegions(),
                              estSval,min.GetPriorStatEstValues()[es]);
        }
      
      for (size_t es = 0; es < estC.size(); es++)
        {
	  vector<double> res = estC[es]->Evaluate(pdf,min.GetIDS(),index,x,Q);
	  for (int l = 0; l < estC[es]->getSize(); l++) estCval[l] = res[l];      
      
          erfCs[es][t] += ERFC(estC[es]->getSize(), estCval, min.GetPriorCorrEstValues()[es]);
        }
      
    }
  cout << endl;

  // printing results
  vector<double> N;
  double cv = 0, md = 0, up50 = 0, dn50 = 0;
  double up68 = 0, dn68 = 0, up90 = 0, dn90 = 0;
  cout.precision(4);
  f.precision(4);
  f << scientific;

  stringstream file0("");
  file0 << priorname.c_str() << "/erf_random.dat";
  f.open(file0.str().c_str(), ios::app | ios::out);
  f << rep << "\t";

  for (size_t es = 0; es < erfMs.size(); es++)
    {
      cout << "\n* Estimator: " << estM[es]->getName() << endl;
      ComputeCV(erfMs[es],cv, md, dn50, up50, dn68, up68, dn90, up90);
      cout << scientific << "*   mean,median = " << cv << " " << md << endl;
      cout << scientific << "* -50%cl,+50%cl = " << dn50 << " " << up50 << endl;
      cout << scientific << "* -68%cl,+68%cl = " << dn68 << " " << up68 << endl;
      cout << scientific << "* -90%cl,+90%cl = " << dn90 << " " << up90 << endl;
      f << cv << "\t" << md << "\t" << dn50 << "\t" << up50 << "\t"
        << dn68 << "\t" << up68 << "\t" << dn90 << "\t" << up90 << "\t";

      if (es < 4) N.push_back(up68);
    }

  for (size_t es = 0; es < erfSs.size(); es++)
    {
      cout << "\n* Estimator: " << estS[es]->getName() << endl;
      ComputeCV(erfSs[es],cv, md, dn50, up50, dn68, up68, dn90, up90);
      cout << scientific << "*   mean,median = " << cv << " " << md << endl;
      cout << scientific << "* -50%cl,+50%cl = " << dn50 << " " << up50 << endl;
      cout << scientific << "* -68%cl,+68%cl = " << dn68 << " " << up68 << endl;
      cout << scientific << "* -90%cl,+90%cl = " << dn90 << " " << up90 << endl;
      f << cv << "\t" << md << "\t" << dn50 << "\t" << up50 << "\t"
        << dn68 << "\t" << up68 << "\t" << dn90 << "\t" << up90 << "\t";

      N.push_back(up68);
    }

  for (size_t es = 0; es < erfCs.size(); es++)
    {
      cout << "\n* Estimator: " << estC[es]->getName() << endl;
      ComputeCV(erfCs[es],cv, md, dn50, up50, dn68, up68, dn90, up90);
      cout << scientific << "*   mean,median = " << cv << " " << md << endl;
      cout << scientific << "* -50%cl,+50%cl = " << dn50 << " " << up50 << endl;
      cout << scientific << "* -68%cl,+68%cl = " << dn68 << " " << up68 << endl;
      cout << scientific << "* -90%cl,+90%cl = " << dn90 << " " << up90 << endl;
      f << cv << "\t" << md << "\t" << dn50 << "\t" << up50 << "\t"
        << dn68 << "\t" << up68 << "\t" << dn90 << "\t" << up90 << "\t";

      N.push_back(up68);
    }

  f << endl;
  f.close();

  if (compress)
    {
      cout << "\n* Error Function Normalizations" << endl;
      cout << "*   CV: " << N[0] << endl;
      cout << "*   SD: " << N[1] << endl;
      cout << "*   SK: " << N[2] << endl;
      cout << "*   KU: " << N[3] << endl;
      cout << "*   KO: " << N[4] << endl;
      cout << "*   CORR: "<< N[5] << endl;

      rg->SetSeed(0);
      min.setupminimizer(rep,N,rg);

      cout << "\n* Compressing:" << endl;
      const int Nite = 15000;
      double e;
      for (int i = 0; i < Nite; i++)
        {
          e = min.iterate();
          cout << "* Iteration " << i+1 << "/" << Nite << "\t ERF = " << e << "\r";
          cout.flush();
        }
      cout << endl;

      index = min.getIndex();

      for (size_t es = 0; es < estM.size(); es++) erfMs[es][0] = 0.0;
      for (size_t es = 0; es < estS.size(); es++) erfSs[es][0] = 0.0;
      for (size_t es = 0; es < estS.size(); es++) erfCs[es][0] = 0.0;

      // computing estimators
      for (size_t es = 0; es < estM.size(); es++)
        {
          for (size_t fl = 0; fl < min.GetIDS().size(); fl++)
            for (int ix = 0; ix < x->size(); ix++)
              estMval[fl][ix] = estM[es]->Evaluate(pdf,min.GetIDS()[fl],index,x->at(ix),Q);
          erfMs[es][0] += ERF(min.GetIDS().size(), x->size(), estMval, min.GetPriorMomentEstValues()[es]);
        }

      // computing estimators
      for (size_t es = 0; es < estS.size(); es++)
        {
          for (size_t fl = 0; fl < min.GetIDS().size(); fl++)
            for (int ix = 0; ix < x->size(); ix++)
              {
                vector<double> res = estS[es]->Evaluate(pdf,min.GetIDS()[fl],index,x->at(ix),Q);
                for (int l = 0; l < estS[es]->getRegions(); l++) estSval[fl][ix][l] = res[l];
              }
          erfSs[es][0] += ERFS(min.GetIDS().size(),x->size(),estS[es]->getRegions(),
                              estSval,min.GetPriorStatEstValues()[es]);
        }
      
      // computing c estimators
      for (size_t es = 0; es < estC.size(); es++)
        {
          vector<double> res = estC[es]->Evaluate(pdf,min.GetIDS(),index,x,Q);	  
	  for (int l = 0; l < estC[es]->getSize(); l++) estCval[l] = res[l];      

          erfCs[es][0] += ERFC(estC[es]->getSize(), estCval, min.GetPriorCorrEstValues()[es]);
        }      

      stringstream file1("");
      file1 << priorname.c_str() << "/erf_compression.dat";
      f.open(file1.str().c_str(), ios::out | ios::app);
      f << rep << "\t";
      for (size_t es = 0; es < erfMs.size(); es++)
        {
          cout << "\n* Estimator: " << estM[es]->getName() << endl;
          cout << scientific << "*   mean = " << erfMs[es][0] << endl;
          f << erfMs[es][0] << "\t";
        }

      for (size_t es = 0; es < erfSs.size(); es++)
        {
          cout << "\n* Estimator: " << estS[es]->getName() << endl;
          cout << scientific << "*   mean = " << erfSs[es][0] << endl;
          f << erfSs[es][0] << "\t";
        }

      for (size_t es = 0; es < erfCs.size(); es++)
        {
          cout << "\n* Estimator: " << estC[es]->getName() << endl;
          cout << scientific << "*   mean = " << erfCs[es][0] << endl;
          f << erfCs[es][0] << "\t";
        }
      f << endl;
      f.close();

      stringstream file2("");
      file2 << priorname.c_str() << "/replica_compression_" << rep << ".dat";
      f.open(file2.str().c_str(), ios::out);
      for (int i = 0; i < rep; i++) f << index[i] << endl;
      f.close();

      cout << "\n* Please run now: ./buildgrid.py " << rep << " " << priorname << endl;
      cout << "* Use validation.C for ERFs plots" << endl;
    }

  delete rg;

  delete[] estCval;

  for (int i = 0; i < (int) min.GetIDS().size(); i++)
    if (estMval[i]) delete[] estMval[i];
  delete[] estMval;

  for (int i = 0; i < (int) min.GetIDS().size(); i++)
    {
      for (int j = 0; j < (int)x->size(); j++)
	if (estSval[i][j]) delete[] estSval[i][j];
      if (estSval[i]) delete[] estSval[i];
    }
  delete[] estSval;

  delete x;

  for (size_t i = 0; i < pdf.size(); i++)
    if (pdf[i]) delete pdf[i];
  pdf.clear();
  

  return 0;
}

void splash()
{
  cout << "   ___                                                    \n" <<
    "  / __\\___  _ __ ___  _ __  _ __ ___  ___ ___  ___  _ __  \n" <<
    " / /  / _ \\| '_ ` _ \\| '_ \\| '__/ _ \\/ __/ __|/ _ \\| '__| \n" <<
    "/ /__| (_) | | | | | | |_) | | |  __/\\__ \\__ \\ (_) | |    \n" <<
    "\\____/\\___/|_| |_| |_| .__/|_|  \\___||___/___/\\___/|_|    \n" <<
    "                     |_|                                  \n" << endl;
  cout << " __v" << VERSION <<"__ Authors: S. Carrazza, J. I. Latorre\n" << endl;
}
