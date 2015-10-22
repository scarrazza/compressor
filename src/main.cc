// Compressor - January 2015
// Author:  Stefano Carrazza
// Contact: stefano.carrazza@mi.infn.it

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
#include "TMatrixD.h"
using namespace std;

void splash();

int main(int argc, char** argv)
{
  int    rep;
  bool   nocompress = false;
  unsigned long int seed = 0;
  string priorname;
  double Q = 1.0;  
  if (argc > 3) { rep = atoi(argv[1]); priorname.assign(argv[2]); Q = atof(argv[3]); }
  else { 
    cout << "\n usage: ./compressor [rep] [pdf-name] [Q] [seed] [no-grid]" << endl;
    cout << "   optional: [seed]      seed value" << endl;
    cout << "             [no-grid]   1 do not export grid\n" << endl;
    exit(-1);
  }
  if (argc >= 5) { seed = atoi(argv[4]); }
  if (argc >= 6) { nocompress = atoi(argv[5]); }

  splash();
  
  cout << "-------------------------------------------" << endl;
  cout << "- Setup summary                           -" << endl;
  cout << "-------------------------------------------" << endl;
  cout << "- Input PDF set      : " << priorname << endl;
  cout << "- Desired compression: " << rep << endl;
  cout << "- Input energy Qin   : " << Q << " GeV" << endl;
  cout << "- Creating out.folder: " << priorname.c_str() << "/" << endl;
  mkdir(priorname.c_str(),0777);  

  // allocate common rg for random testing before fit
  cout << "- Seed               : " << seed << endl;
  RandomGenerator *rg = new RandomGenerator(0,seed);
  Grid *x = new Grid();

  cout << "- X grid size        : " << x->size() << " points, x=["
       << x->at(0) << ", " << x->at(x->size()-1) << "]" << endl;
  cout << "-------------------------------------------" << endl;

  // allocate LHAPDF set
  cout << "\n-------------------------------------------" << endl;
  cout << "- Loading grid with LHAPDF6               -" << endl;
  cout << "-------------------------------------------" << endl;
  const LHAPDF::PDFSet set(priorname.c_str());
  vector<LHAPDF::PDF*> lhapdf = set.mkPDFs();

  cout << "\n-------------------------------------------" << endl;
  cout << "- Preloading grid in memory               -" << endl;
  cout << "-------------------------------------------" << endl;  
  int nf = 3;  
  LocalPDF *pdf = new LocalPDF(lhapdf,nf,x,Q);  
  Minimizer min(pdf,x,nf);

  vector<EstimatorsM*> estM = min.GetMomentEstimators();
  vector<EstimatorsS*> estS = min.GetStatEstimators();
  vector<EstimatorsC*> estC = min.GetCorrEstimators();
  TMatrixD invPrior(min.GetPriorInvMatrix());

  // Computing error function for random set
  const int trials = 1000;
  cout << "\n* Random trials: " << trials << endl;
  vector<int> index;
  vector<double> w;
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
      w.resize(rep,1);
      // generate random replicas, no duplicates
      randomize(pdf->size()-1,rg,index);

      // computing estimators
      for (size_t es = 0; es < estM.size(); es++)
        {
          for (size_t fl = 0; fl < min.GetIDS().size(); fl++)
            for (int ix = 0; ix < x->size(); ix++)
              estMval[fl][ix] = estM[es]->Evaluate(pdf,min.GetIDS()[fl],index,w,ix);
          erfMs[es][t] += ERF(min.GetIDS().size(), x->size(), estMval, min.GetPriorMomentEstValues()[es]);
        }

      // computing estimators
      for (size_t es = 0; es < estS.size(); es++)
        {
          for (size_t fl = 0; fl < min.GetIDS().size(); fl++)
            for (int ix = 0; ix < x->size(); ix++)
              {
                vector<double> res = estS[es]->Evaluate(pdf,min.GetIDS()[fl],index,w,ix);
                for (int l = 0; l < estS[es]->getRegions(); l++) estSval[fl][ix][l] = res[l];
              }
          erfSs[es][t] += ERFS(min.GetIDS().size(),x->size(),estS[es]->getRegions(),
                              estSval,min.GetPriorStatEstValues()[es]);
        }

      for (size_t es = 0; es < estC.size(); es++)
        {
          TMatrixD m = estC[es]->Evaluate(pdf,min.GetIDS(),index,w,x);
          TMatrixD r = m*invPrior;

	  estCval[0] = 0;
          for (int l = 0; l < estC[es]->getSize(); l++) estCval[0] += r(l,l);

          erfCs[es][t] += ERFC(1, estCval, min.GetPriorCorrEstValues()[es]);
        }

    }
  cout << endl;

  // printing results
  vector<double> N;
  double cv = 0, md = 0, up50 = 0, dn50 = 0;
  double up68 = 0, dn68 = 0, up90 = 0, dn90 = 0;
  cout.precision(4);
  fstream f, f2;
  f.precision(4);
  f << scientific;

  stringstream file0("");
  file0 << priorname.c_str() << "/erf_random.dat";
  f.open(file0.str().c_str(), ios::app | ios::out);
  f << rep << "\t";

  for (size_t es = 0; es < erfMs.size(); es++)
    {
      cout << "\n- Estimator: " << estM[es]->getName() << endl;
      ComputeCV(erfMs[es],cv, md, dn50, up50, dn68, up68, dn90, up90);
      cout << scientific << "-   mean,median = " << cv << " " << md << endl;
      cout << scientific << "- [-50%cl,+50%cl] = " << dn50 << " " << up50 << endl;
      cout << scientific << "- [-68%cl,+68%cl] = " << dn68 << " " << up68 << endl;
      cout << scientific << "- [-90%cl,+90%cl] = " << dn90 << " " << up90 << endl;
      f << cv << "\t" << md << "\t" << dn50 << "\t" << up50 << "\t"
        << dn68 << "\t" << up68 << "\t" << dn90 << "\t" << up90 << "\t";

      //if (es < 4) N.push_back(up68);
      N.push_back(up68);
    }

  for (size_t es = 0; es < erfSs.size(); es++)
    {
      cout << "\n- Estimator: " << estS[es]->getName() << endl;
      ComputeCV(erfSs[es],cv, md, dn50, up50, dn68, up68, dn90, up90);
      cout << scientific << "-   mean,median = " << cv << " " << md << endl;
      cout << scientific << "- [-50%cl,+50%cl] = " << dn50 << " " << up50 << endl;
      cout << scientific << "- [-68%cl,+68%cl] = " << dn68 << " " << up68 << endl;
      cout << scientific << "- [-90%cl,+90%cl] = " << dn90 << " " << up90 << endl;
      f << cv << "\t" << md << "\t" << dn50 << "\t" << up50 << "\t"
        << dn68 << "\t" << up68 << "\t" << dn90 << "\t" << up90 << "\t";

      N.push_back(up68);
    }

  for (size_t es = 0; es < erfCs.size(); es++)
    {
      cout << "\n- Estimator: " << estC[es]->getName() << endl;
      ComputeCV(erfCs[es],cv, md, dn50, up50, dn68, up68, dn90, up90);
      cout << scientific << "-   mean,median = " << cv << " " << md << endl;
      cout << scientific << "- [-50%cl,+50%cl] = " << dn50 << " " << up50 << endl;
      cout << scientific << "- [-68%cl,+68%cl] = " << dn68 << " " << up68 << endl;
      cout << scientific << "- [-90%cl,+90%cl] = " << dn90 << " " << up90 << endl;
      f << cv << "\t" << md << "\t" << dn50 << "\t" << up50 << "\t"
        << dn68 << "\t" << up68 << "\t" << dn90 << "\t" << up90 << "\t";

      N.push_back(up68);
    }

  f << endl;
  f.close();

  if (!nocompress)
    {
      cout << "\n- Error Function Normalizations" << endl;
      cout << "-   CV : " << N[0] << endl;
      cout << "-   SD : " << N[1] << endl;
      cout << "-   SK : " << N[2] << endl;
      cout << "-   KU : " << N[3] << endl;
      cout << "-   KO : " << N[4] << endl;
      cout << "-   EIG: " << N[5] << endl;

      rg->SetSeed(0);
      min.setupminimizer(rep,N,rg);

      cout << "\n- Compressing:" << endl;
      const int Nite = 30000;
      double e;
      for (int i = 0; i < Nite; i++)
        {
          e = min.iterate();
          e = min.iterate_w();
          cout << "* Iteration " << i+1 << "/" << Nite << "\t ERF = " << e << "\r";
          cout.flush();
        }
      cout << endl;      

      index = min.getIndex();
      w = min.getW();

      for (size_t es = 0; es < estM.size(); es++) erfMs[es][0] = 0.0;
      for (size_t es = 0; es < estS.size(); es++) erfSs[es][0] = 0.0;
      for (size_t es = 0; es < estS.size(); es++) erfCs[es][0] = 0.0;

      // computing estimators
      for (size_t es = 0; es < estM.size(); es++)
        {
          for (size_t fl = 0; fl < min.GetIDS().size(); fl++)
            for (int ix = 0; ix < x->size(); ix++)
              estMval[fl][ix] = estM[es]->Evaluate(pdf,min.GetIDS()[fl],index,w,ix);
          erfMs[es][0] += ERF(min.GetIDS().size(), x->size(), estMval, min.GetPriorMomentEstValues()[es]);
        }

      // computing estimators
      for (size_t es = 0; es < estS.size(); es++)
        {
          for (size_t fl = 0; fl < min.GetIDS().size(); fl++)
            for (int ix = 0; ix < x->size(); ix++)
              {
                vector<double> res = estS[es]->Evaluate(pdf,min.GetIDS()[fl],index,w,ix);
                for (int l = 0; l < estS[es]->getRegions(); l++) estSval[fl][ix][l] = res[l];
              }
          erfSs[es][0] += ERFS(min.GetIDS().size(),x->size(),estS[es]->getRegions(),
                              estSval,min.GetPriorStatEstValues()[es]);
        }

      // computing c estimators
      for (size_t es = 0; es < estC.size(); es++)
        {
          TMatrixD m = estC[es]->Evaluate(pdf,min.GetIDS(),index,w,x);
          TMatrixD r = m*invPrior;

          estCval[0] = 0;
          for (int l = 0; l < estC[es]->getSize(); l++) estCval[0] += r(l,l);

          erfCs[es][0] += ERFC(1, estCval, min.GetPriorCorrEstValues()[es]);
        }


      stringstream file1("");
      file1 << priorname.c_str() << "/erf_compression.dat";
      f.open(file1.str().c_str(), ios::out | ios::app);
      f << rep << "\t";
      for (size_t es = 0; es < erfMs.size(); es++)
        {
          cout << "\n- Estimator: " << estM[es]->getName() << endl;
          cout << scientific << "-   mean = " << erfMs[es][0] << endl;
          f << erfMs[es][0] << "\t";
        }

      for (size_t es = 0; es < erfSs.size(); es++)
        {
          cout << "\n- Estimator: " << estS[es]->getName() << endl;
          cout << scientific << "-   mean = " << erfSs[es][0] << endl;
          f << erfSs[es][0] << "\t";
        }

      for (size_t es = 0; es < erfCs.size(); es++)
        {
          cout << "\n- Estimator: " << estC[es]->getName() << endl;
          cout << scientific << "-   mean = " << erfCs[es][0] << endl;
          f << erfCs[es][0] << "\t";
        }
      f << endl;
      f.close();

      stringstream file2(""), file3("");
      file2 << priorname.c_str() << "/replica_compression_" << rep << ".dat";
      file3 << priorname.c_str() << "/replica_weights_" << rep << ".dat";
      f.open(file2.str().c_str(), ios::out);
      f2.open(file3.str().c_str(), ios::out);
      for (int i = 0; i < rep; i++)
        {
          f << index[i] << endl;
          f2 << w[i] << endl;
        }
      f.close();
      f2.close();

      cout << "\n- In order to create the compressed grid in the LHAPDF6 format" << endl;
      cout << "  Please run now: ./compressor_buildgrid " << rep << " " << priorname << endl;
      cout << "- Use compressor_validate.C in order to generate the ERFs plots\n" << endl;
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
  for (size_t i = 0; i < lhapdf.size(); i++)
    if (lhapdf[i]) delete lhapdf[i];
  lhapdf.clear();

  return 0;
}

void splash()
{
  cout << "    ___                                                    \n" <<
    "   / __\\___  _ __ ___  _ __  _ __ ___  ___ ___  ___  _ __  \n" <<
    "  / /  / _ \\| '_ ` _ \\| '_ \\| '__/ _ \\/ __/ __|/ _ \\| '__| \n" <<
    " / /__| (_) | | | | | | |_) | | |  __/\\__ \\__ \\ (_) | |    \n" <<
    " \\____/\\___/|_| |_| |_| .__/|_|  \\___||___/___/\\___/|_|    \n" <<
    "                      |_|                                  \n" << endl;
  cout << "  __v" << VERSION <<"__ Author: S. Carrazza\n" << endl;
}
