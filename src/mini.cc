#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "mini.h"
#include "randomgenerator.h"
#include "LHAPDF/LHAPDF.h"
#define Q 1.0
using namespace std;

// compute
double Mini::ComputeAVG(double x, double q, int f, int* index)
{
  double sum = 0;
  for (size_t i = 0; i < fRep; i++)
    sum += fPDF[index[i]]->xfxQ(f, x, q);

  return sum / (double) fRep;
}

Mini::Mini(int rep, vector<LHAPDF::PDF*> pdf):
  fRep(rep),
  fNMut(5),
  fRg(NULL),
  fPDF(pdf)
{
  fRg = new RandomGenerator(0,0);

  for (int i = 0; i < (int) fNMut; i++)
    fMut.push_back(new int[fRep]);
  
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
  fX.push_back(0.92653061224489797);
  fX.push_back(0.94489795918367347);
  fX.push_back(0.96326530612244898);
  fX.push_back(0.98163265306122449);
  fX.push_back(1.0000000000000000);

  fPids = pdf[0]->flavors();

  for (int i = 0; i < (int) fPids.size(); i++)
    {
      fCV.push_back(new double[fX.size()]);
      for (int j = 0; j < (int) fX.size(); j++)
        fCV[i][j] = pdf[0]->xfxQ(fPids[i], fX[j], Q);
    }

}

Mini::~Mini()
{
  delete fRg;
  for (int i = 0; i < (int) fCV.size(); i++)
    if (fCV[i]) delete[] fCV[i];
  fCV.clear();

  for (int i = 0; i < (int) fMut.size(); i++)
    if (fMut[i]) delete[] fMut[i];
  fMut.clear();
}

double Mini::iterate(int* index)
{
  // ERF
  double berf = 0;
  for (int f = 0; f < (int) fPids.size(); f++)
    for (int i = 0; i < (int) fX.size(); i++)
      berf += pow(fCV[f][i] - ComputeAVG(fX[i], Q, fPids[f], index), 2.0);

  //
  for (int i = 0; i < fNMut; i++)
    for (int j = 0; j < (int) fRep; j++)
        fMut[i][j] = index[j];

  // GA mutation
  for (int i = 0; i < fNMut; i++)
    {
      const double g = fRg->GetRandomUniform();

      int nmut = 3;
      if (g < 0.7) nmut = 1;
      else if (g < 0.9 && g >= 0.7) nmut = 2;

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
          erf[i] += pow(fCV[f][j] - ComputeAVG(fX[j], Q, fPids[f], fMut[i]), 2.0);
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
    {
      for (int i = 0; i < (int) fRep; i++)
        index[i] = fMut[id][i];
    }
  else bestchi2 = berf;

  delete[] erf;

  return bestchi2;
}

void Mini::Save(int *index, string dir)
{
  fstream f;

  for (int i = 0; i < (int) fPids.size(); i++)
    {
      stringstream file("");
      file << dir << "/scatter_" << fPids[i] << ".dat";
      f.open(file.str().c_str(), ios::out);
      for (int ix = 0; ix < (int) fX.size(); ix++)
        f << scientific << fCV[i][ix] << "\t" << ComputeAVG(fX[ix],Q,fPids[i],index) << endl;
      f.close();
    }
}
