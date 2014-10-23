#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "mini.h"
#include "randomgenerator.h"
#include "LHAPDF/LHAPDF.h"
using namespace std;

const double Ncv[] = { 
  2.497476e+01+2.027556e+01,
  1.616557e+01+1.377080e+01,
  1.173817e+01+1.008291e+01,
  9.756332e+00+8.513520e+00,
  7.968408e+00+6.916763e+00,
  6.725073e+00+5.932574e+00,
  5.624087e+00+4.909153e+00,
  5.517243e+00+4.868520e+00,
  4.440808e+00+3.905348e+00,
  4.095582e+00+3.591852e+00
};

const double Nsd[] = {
  2.206626e+02+1.443792e+02,
  1.878154e+02+1.357397e+02,
  1.593149e+02+1.354671e+02,
  1.500699e+02+1.316032e+02,
  1.400564e+02+1.231716e+02,
  1.246468e+02+1.093228e+02,
  1.147630e+02+1.013076e+02,
  1.132564e+02+1.007084e+02,
  9.756910e+01+8.582596e+01,
  9.428089e+01+8.356309e+01
};

const double Nsk[] = {
  1.902686e+03+2.588718e+02,
  1.603385e+03+2.911056e+02,
  1.394172e+03+2.962958e+02,
  1.270620e+03+2.899914e+02,
  1.164156e+03+2.930363e+02,
  1.055383e+03+2.935948e+02,
  9.762340e+02+2.952764e+02,
  9.210893e+02+3.199038e+02,
  8.476124e+02+3.134301e+02,
  8.190550e+02+3.147638e+02
};

const double Nku[] = {
  2.998481e+05+4.841705e+03,
  2.828328e+05+1.210993e+04,
  2.668503e+05+2.028426e+04,
  2.553422e+05+2.476674e+04,
  2.443311e+05+3.038502e+04,
  2.306691e+05+3.821358e+04,
  2.203701e+05+4.453629e+04,
  2.122242e+05+4.988864e+04,
  2.032739e+05+5.451479e+04,
  1.968498e+05+5.814638e+04
};

const double Nko[] = {
  5.158859e+01+3.530100e+00,
  2.725654e+01+2.090950e+00,
  1.900916e+01+1.612000e+00,
  1.485242e+01+1.204225e+00,
  1.231963e+01+1.091860e+00,
  1.056806e+01+9.173944e-01,
  9.293225e+00+7.847041e-01,
  8.393598e+00+7.721688e-01,
  7.593241e+00+6.450926e-01,
  6.978772e+00+6.071000e-01
};

// compute
void Mini::ComputeEstimators(int n, double x, double q, int f, int* index, 
			     double &cv, double &std, double &ske, 
			     double &kur, double *res)
{
  double sum = 0, sq_sum = 0;
  double *pdf = new double[n];
  for (size_t i = 0; i < n; i++)
    {
      pdf[i] = fPDF[index[i]]->xfxQ(f, x, q);
      sum += pdf[i];
      sq_sum += pdf[i]*pdf[i]; 
    }

  cv = sum / (double) n;
  std= sqrt(sq_sum / (double) n - cv*cv);

  double sq_sum2 = 0, cub_sum = 0;
  for (size_t i = 0; i < n; i++)
    {
      const double v = (pdf[i] - cv);
      cub_sum += v*v*v;
      sq_sum2 += v*v*v*v;
    }

  ske = (cub_sum / n) / (std * std * std);
  kur = (sq_sum2 / n) / (std * std * std * std) - 3.0;

  // Kolmogoroz
  for (int l = 0; l < 6; l++) res[l] = 0;
 
  for (size_t i = 0; i < n; i++)
    {
      const double v = pdf[i];
      if (v <= cv -2*std)
	res[0] += 1;
      else if (v <= cv - std)
	res[1] += 1;
      else if (v <= cv)
	res[2] += 1;
      else if (v <= cv + std)
	res[3] += 1;
      else if (v <= cv + 2*std)
	res[4] += 1;
      else if (v > cv + 2*std)
	res[5] += 1;
      else
	cout << "ERROR" << endl;
    }

  for (int l = 0; l < 6; l++) res[l] /= (double) n;

  delete[] pdf;
}

Mini::Mini(int rep, vector<LHAPDF::PDF*> pdf, int seed):
  fRep(rep),
  fNMut(5),
  fRg(NULL),
  fPDF(pdf)
{
  fRg = new RandomGenerator(0,seed);

  for (int i = 0; i < (int) fNMut; i++)
    fMut.push_back(new int[fRep]);
  
  /*
  fX.push_back(1.0000000000000001E-009);
  fX.push_back(1.4508287784959398E-009);
  fX.push_back(2.1049041445120207E-009);
  fX.push_back(3.0538555088334157E-009);
  fX.push_back(4.4306214575838816E-009);
  fX.push_back(6.4280731172843209E-009);
  fX.push_back(9.3260334688321995E-009);
  fX.push_back(1.3530477745798068E-008);
  fX.push_back(1.9630406500402714E-008);
  fX.push_back(2.8480358684358022E-008);
  fX.push_back(4.1320124001153370E-008);
  fX.push_back(5.9948425031894094E-008);
  fX.push_back(8.6974900261778356E-008);  
  fX.push_back(1.2618568830660210E-007);
  fX.push_back(1.8307382802953678E-007);
  fX.push_back(2.6560877829466870E-007);
  fX.push_back(3.8535285937105315E-007);
  fX.push_back(5.5908101825122239E-007);
  fX.push_back(8.1113083078968731E-007);  
  fX.push_back(1.1768119524349981E-006);
  fX.push_back(1.7073526474706905E-006);
  fX.push_back(2.4770763559917115E-006);
  fX.push_back(3.5938136638046262E-006);
  fX.push_back(5.2140082879996849E-006);
  fX.push_back(7.5646332755462914E-006);
  */
  fX.push_back(1.0974987654930569E-005);
  fX.push_back(1.5922827933410941E-005);
  fX.push_back(2.3101297000831580E-005);
  fX.push_back(3.3516026509388410E-005);
  fX.push_back(4.8626015800653536E-005);
  fX.push_back(7.0548023107186455E-005);
  
  fX.push_back(1.0235310218990269E-004);
  fX.push_back(1.4849682622544667E-004);
  fX.push_back(2.1544346900318823E-004);
  fX.push_back(3.1257158496882353E-004);
  fX.push_back(4.5348785081285824E-004);
  fX.push_back(6.5793322465756835E-004);
  fX.push_back(9.5454845666183481E-004);
  fX.push_back(1.3848863713938717E-003);
  fX.push_back(2.0092330025650459E-003);
  fX.push_back(2.9150530628251760E-003);
  fX.push_back(4.2292428743894986E-003);
  fX.push_back(6.1359072734131761E-003);
  fX.push_back(8.9021508544503934E-003);
  fX.push_back(1.2915496650148829E-002);
  fX.push_back(1.8738174228603830E-002);
  fX.push_back(2.7185882427329403E-002);
  fX.push_back(3.9442060594376556E-002);
  fX.push_back(5.7223676593502207E-002);
  fX.push_back(8.3021756813197525E-002);
  fX.push_back(0.10000000000000001);
  fX.push_back(0.11836734693877551);
  fX.push_back(0.13673469387755102);
  fX.push_back(0.15510204081632653);
  fX.push_back(0.17346938775510204);
  fX.push_back(0.19183673469387758);
  fX.push_back(0.21020408163265308);
  fX.push_back(0.22857142857142856);
  fX.push_back(0.24693877551020407);
  fX.push_back(0.26530612244897961);
  fX.push_back(0.28367346938775512);
  fX.push_back(0.30204081632653063);
  fX.push_back(0.32040816326530613);
  fX.push_back(0.33877551020408170);
  fX.push_back(0.35714285714285710);
  fX.push_back(0.37551020408163271);
  fX.push_back(0.39387755102040811);
  fX.push_back(0.41224489795918373);
  fX.push_back(0.43061224489795924);
  fX.push_back(0.44897959183673475);
  fX.push_back(0.46734693877551026);
  fX.push_back(0.48571428571428565);
  fX.push_back(0.50408163265306127);
  fX.push_back(0.52244897959183678);
  fX.push_back(0.54081632653061229);
  fX.push_back(0.55918367346938780);
  fX.push_back(0.57755102040816331);
  fX.push_back(0.59591836734693870);
  fX.push_back(0.61428571428571421);
  fX.push_back(0.63265306122448983);
  fX.push_back(0.65102040816326534);
  
  fX.push_back(0.66938775510204085);
  fX.push_back(0.68775510204081625);
  fX.push_back(0.70612244897959175);
  fX.push_back(0.72448979591836737);
  fX.push_back(0.74285714285714288);
  fX.push_back(0.76122448979591839);
  fX.push_back(0.77959183673469379);
  fX.push_back(0.79795918367346941);
  fX.push_back(0.81632653061224492);
  fX.push_back(0.83469387755102042);
  fX.push_back(0.85306122448979593);
  fX.push_back(0.87142857142857133);
  fX.push_back(0.88979591836734695);
  fX.push_back(0.90816326530612246);
  /*
  fX.push_back(0.92653061224489797);
  fX.push_back(0.94489795918367347);
  fX.push_back(0.96326530612244898);
  fX.push_back(0.98163265306122449);
  fX.push_back(1.0000000000000000);
  */

  fPids = pdf[0]->flavors();

  int *index = new int[pdf.size()-1];
  for (int i = 0; i < (int) pdf.size()-1; i++) index[i] = i+1;

  for (int i = 0; i < (int) fPids.size(); i++)
    {
      fCV.push_back(new double[fX.size()]);
      fSD.push_back(new double[fX.size()]);
      fSK.push_back(new double[fX.size()]);
      fKU.push_back(new double[fX.size()]);
      fKO.push_back(new double*[fX.size()]);

      for (int j = 0; j < (int) fX.size(); j++)
	{
	  double cv = 0, std = 0, kur = 0, sk = 0;
	  double *res = new double[6];
	  ComputeEstimators((int) pdf.size() - 1,fX[j], Q, fPids[i], index, 
			    cv, std, sk, kur, res);
	  fCV[i][j] = cv;
	  fSD[i][j] = std;
	  fSK[i][j] = sk;
	  fKU[i][j] = kur;
	  
	  fKO[i][j] = new double[6];
	  for (int l = 0; l < 6; l++) fKO[i][j][l] = res[l];

	  delete[] res;
	}
    }

  delete[] index;

}

Mini::~Mini()
{
  delete fRg;
  for (int i = 0; i < (int) fCV.size(); i++)
    if (fCV[i]) delete[] fCV[i];
  fCV.clear();

  for (int i = 0; i < (int) fSD.size(); i++)
    if (fSD[i]) delete[] fSD[i];
  fSD.clear();

  for (int i = 0; i < (int) fKU.size(); i++)
    if (fKU[i]) delete[] fKU[i];
  fKU.clear();

  for (int i = 0; i < (int) fSK.size(); i++)
    if (fSK[i]) delete[] fSK[i];
  fSK.clear();

  for (int i = 0; i < (int) fKO.size(); i++)
    {
      for (int j = 0; j < (int) fX.size(); j++)
	if (fKO[i][j]) delete[] fKO[i][j];
      if (fKO[i]) delete[] fKO[i];
    }  
  fKO.clear();

  for (int i = 0; i < (int) fMut.size(); i++)
    if (fMut[i]) delete[] fMut[i];
  fMut.clear();
}

double Mini::iterate(int* index)
{
  // ERF
  double berf = 0;
  double *res = new double[6];
  for (int f = 0; f < (int) fPids.size(); f++)
    for (int i = 0; i < (int) fX.size(); i++)
      {
	double cv = 0, std = 0, sk = 0, kur = 0;
	ComputeEstimators(fRep, fX[i], Q, fPids[f], index, 
			  cv, std, sk, kur,res);
	berf += pow(fCV[f][i] - cv,  2.0) / Ncv[fRep / 10 - 1];
	berf += pow(fSD[f][i] - std, 2.0) / Nsd[fRep / 10 - 1];
	berf += pow(fSK[f][i] - sk,  2.0) / Nsk[fRep / 10 - 1];
	berf += pow(fKU[f][i] - kur, 2.0) / Nku[fRep / 10 - 1];
	
	for (int l = 0; l < 6; l++) {
	  const double v = fKO[f][i][l] - res[l];
	  berf += v*v / Nko[fRep / 10 - 1];
	}
      }

  // set mut
  for (int i = 0; i < fNMut; i++)
    for (int j = 0; j < (int) fRep; j++)
        fMut[i][j] = index[j];

  // GA mutation
  for (int i = 0; i < fNMut; i++)
    {
      const double g = fRg->GetRandomUniform();
      
      int nmut = 4; // 10%
      if (g <= 0.3) nmut = 1; // 30%
      else if (g > 0.3 && g <= 0.6) nmut = 2; // 30%
      else if (g > 0.6 && g <= 0.7) nmut = 3; // 10%
     
      for (int t = 0; t < nmut; t++)
        {
          const int pos = fRg->GetRandomUniform(fPDF.size()-1)+1;
          fMut[i][fRg->GetRandomUniform(fRep)] = pos;
        }
    }

  // Compute ERF
  double *erf = new double[fNMut];
  for (int i = 0; i < fNMut; i++)
    {
      erf[i] = 0.0;
      for (int f = 0; f < (int) fPids.size(); f++)
        for (int j = 0; j < (int) fX.size(); j++)
	  {
	    double cv = 0, std = 0, kur = 0, sk = 0;
	    ComputeEstimators(fRep, fX[j], Q, fPids[f], fMut[i],
			      cv,std,sk,kur,res);
	    erf[i] += pow(fCV[f][j] - cv,  2.0) / Ncv[fRep / 10 - 1];
	    erf[i] += pow(fSD[f][j] - std, 2.0) / Nsd[fRep / 10 - 1];
	    erf[i] += pow(fSK[f][j] - sk,  2.0) / Nsk[fRep / 10 - 1];
	    erf[i] += pow(fKU[f][j] - kur, 2.0) / Nku[fRep / 10 - 1];
	    
	    for (int l = 0; l < 6; l++) {
	      const double v = fKO[f][j][l] - res[l];
	      erf[i] += v*v / Nko[fRep / 10 - 1];
	    }
	  }
    }

  // Selection
  int id = 0;
  double bestchi2 = erf[0];
  for (int i = 0; i < fNMut; i++)
    if (erf[i] < bestchi2)
      {
        bestchi2 = erf[i];
        id = i;
      }

  if (bestchi2 < berf)
    for (int i = 0; i < (int) fRep; i++) index[i] = fMut[id][i];   
  else bestchi2 = berf;

  delete[] erf;
  delete[] res;
  return bestchi2;
}
