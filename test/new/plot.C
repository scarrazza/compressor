void plot()
{
  vector<string> filename;
  filename.push_back("cv_erf.dat");  
  filename.push_back("sd_erf.dat");  
  filename.push_back("sk_erf.dat");  
  filename.push_back("ku_erf.dat");  
  filename.push_back("ko_erf.dat");
  
  TMultiGraph** mg = new TMultiGraph*[filename.size()];

  TLegend *l = new TLegend(0.5,0.67,0.88,0.88);
  l->SetBorderSize(0);
  l->SetFillColor(0);
  l->SetFillStyle(0);

  int rep;
  double cv, md, dn50, up50, dn68, up68, dn90, up90;
  for (int i = 0; i < (int) filename.size(); i++)
    {
      stringstream file("");
      file << "set_r_erf/" << filename[i];

      stringstream file2("");
      file2 << "set_c_erf/" << filename[i];

      TGraphErrors *comp = new TGraphErrors(file2.str().c_str(),"%lg %lg %lg");
      comp->SetMarkerStyle(20);
      comp->SetMarkerColor(kRed);
      comp->SetLineColor(kRed);
      TGraphErrors *g1 = new TGraphErrors(22);
      TGraphErrors *g2 = new TGraphErrors(22);
      TGraphErrors *g3 = new TGraphErrors(22);
      TGraphErrors *g4 = new TGraphErrors(22);
      TGraphErrors *g5 = new TGraphErrors(22);

      fstream f;
      f.open(file.str().c_str(),ios::in);
      int index = 0;
      for(;;)
	{
	  f >> rep >> cv >> md >> dn50 >> up50 >> dn68 >> up68 >> dn90 >> up90;
	  if (f.eof()) break;
	  
	  g1->SetPoint(index, rep, (up68+dn68)/2);
	  g1->SetPointError(index, 0, up68-(up68+dn68)/2);
	  g2->SetPoint(index, rep, md);
	  g3->SetPoint(index, rep, (up50+dn50)/2);
	  g3->SetPointError(index, 0, up50-(up50+dn50)/2);
	  g4->SetPoint(index, rep, (up90+dn90)/2);
	  g4->SetPointError(index, 0, up90-(up90+dn90)/2);
	  g5->SetPoint(index, rep, cv);
	  
	  index++;
	}
      
      g1->SetMarkerColor(kBlue);
      g1->SetLineColor(kGreen+2);
      g2->SetMarkerStyle(20);
      g2->SetMarkerColor(kBlue);
      g3->SetLineColor(kBlue);
      g4->SetLineColor(kOrange);
      g5->SetMarkerStyle(24);
      g5->SetMarkerColor(kBlue);
      
      mg[i] = new TMultiGraph();

      mg[i]->Add(g4,"p0");
      mg[i]->Add(g1,"p0");
      mg[i]->Add(g2,"p0");
      mg[i]->Add(g3,"p0");
      mg[i]->Add(g5,"p0");
      mg[i]->Add(comp,"p");

      if (i == 0)
	{
	  l->AddEntry(comp,"Compressed (20)","p");
	  l->AddEntry(g5,"Random Mean (1k)","p");
	  l->AddEntry(g2,"Random Median (1k)","p");
	  l->AddEntry(g3,"Random 50% c.l. (1k)","l");
	  l->AddEntry(g1,"Random 68% c.l. (1k)","l");
	  l->AddEntry(g4,"Random 90% c.l. (1k)","l");	  
	}
    }
 
  TCanvas *c = new TCanvas("c","",1200,700);

  c->Divide(3,2);

  vector<string> label;
  label.push_back("Error Function: Central Value");
  label.push_back("Error Function: Standard Deviation");
  label.push_back("Error Function: Skewness");
  label.push_back("Error Function: Kurtosis");
  label.push_back("Error Function: Kolmogorov");

  double Ncv[] = {
    4.095582e+00-3.591852e+00,
    9.428089e+01-8.356309e+01,
    8.190550e+02-3.147638e+02,
    1.968498e+05-5.814638e+04,
    6.978772e+00-6.071000e-01
  };
  
  for (int i = 0; i < (int) filename.size(); i++)
    {
      c->cd(i+1);
      c->cd(i+1)->SetLogx();     
      c->cd(i+1)->SetTickx();
      c->cd(i+1)->SetTicky();
      
      mg[i]->SetTitle(label[i].c_str());
      mg[i]->Draw("A");
      mg[i]->GetXaxis()->SetTitle("Number of Replicas");
      mg[i]->GetXaxis()->CenterTitle(true);      
      l->Draw("same");
               
      TGraph *g = new TGraph(2);
      g->SetPoint(0, 10, Ncv[i]);
      g->SetPoint(1, 1000, Ncv[i]);
      g->SetLineColor(kBlue);
      g->SetLineStyle(2);
      g->Draw("same,L");
      if (i == 0) l->AddEntry(g,"Lower 68% band","l");
      /*
      TGraph *gd = new TGraph(2);
      gd->SetPoint(0, 10, cv2[i]);
      gd->SetPoint(1, 1000, cv2[i]);
      gd->SetLineColor(kBlue);
      gd->SetLineStyle(2);
      //gd->Draw("same,L");
      */
    }

}
