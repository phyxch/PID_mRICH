#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <cmath>
#include <cassert>
#include <string>
#include <iostream>

void calMeanRadiaus()
{
  string date = "BeamTest";
  string inputfile = Form("/work/eic/xusun/output/database/BeamTest/database_%s_1.root",date.c_str());
  cout << "read in file: " << inputfile.c_str() << endl;
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  assert(File_InPut);

  string key_proton = "h_mPhotonDist_proton_vx_0_vy_0_mom_0_theta_0_phi_0";
  TH2D *h_PDF_proton = (TH2D*)File_InPut->Get(key_proton.c_str());
  assert(h_PDF_proton);

  TCanvas *c_Identified = new TCanvas("c_Identified","c_Identified",10,10,800,800);
  c_Identified->cd();
  c_Identified->cd()->SetLeftMargin(0.15);
  c_Identified->cd()->SetBottomMargin(0.15);
  c_Identified->cd()->SetGrid(0,0);
  c_Identified->cd()->SetTicks(1,1);

  h_PDF_proton->SetTitle("120 GeV/c proton");
  h_PDF_proton->SetStats(0);

  h_PDF_proton->GetXaxis()->SetTitle("out_x (mm)");
  h_PDF_proton->GetXaxis()->CenterTitle();
  h_PDF_proton->GetXaxis()->SetTitleOffset(1.25);
  h_PDF_proton->GetXaxis()->SetNdivisions(505);

  h_PDF_proton->GetYaxis()->SetTitle("out_y (mm)");
  h_PDF_proton->GetYaxis()->CenterTitle();
  h_PDF_proton->GetYaxis()->SetTitleOffset(1.25);
  h_PDF_proton->GetYaxis()->SetNdivisions(505);

  h_PDF_proton->Draw("colz");

  const int mNumOfPixels = 37; // 16*2 3mm-pixels + 2*2 2mm-glasswindow + 1 1mm-gap
  const double mPixels[mNumOfPixels+1] = {-52.5,-50.5,-47.5,-44.5,-41.5,-38.5,-35.5,-32.5,-29.5,-26.5,-23.5,-20.5,-17.5,-14.5,-11.5,-8.5,-5.5,-2.5,-0.5,0.5,2.5,5.5,8.5,11.5,14.5,17.5,20.5,23.5,26.5,29.5,32.5,35.5,38.5,41.5,44.5,47.5,50.5,52.5};
  TH1D *h_radiaus = new TH1D("h_radiaus","h_radiaus",mNumOfPixels,mPixels);

  for(int i_binX = 0; i_binX < mNumOfPixels; ++i_binX)
  {
    for(int i_binY = 0; i_binY < mNumOfPixels; ++i_binY)
    {
      int globalBin = h_PDF_proton->GetBin(i_binX+1,i_binY+1);
      float weight = h_PDF_proton->GetBinContent(globalBin);
      float out_x = h_PDF_proton->GetXaxis()->GetBinCenter(i_binX+1);
      float out_y = h_PDF_proton->GetYaxis()->GetBinCenter(i_binY+1);
      float radiaus = TMath::Sqrt(out_x*out_x+out_y*out_y);
      cout << "i_binX = " << i_binX << ", out_x = " << out_x << endl;
      cout << "i_binY = " << i_binY << ", out_y = " << out_y << endl;
      cout << "globalBin = " << globalBin << ", radiaus =" << radiaus << ", weight = " << weight << endl;
      cout << endl;
      h_radiaus->Fill(radiaus,weight);
    }
  }

  TCanvas *c_radiaus = new TCanvas("c_radiaus","c_radiaus",10,10,800,800);
  c_radiaus->cd();
  c_radiaus->cd()->SetLeftMargin(0.15);
  c_radiaus->cd()->SetBottomMargin(0.15);
  c_radiaus->cd()->SetGrid(0,0);
  c_radiaus->cd()->SetTicks(1,1);

  h_radiaus->SetTitle("radiaus resolution");
  h_radiaus->SetStats(0);

  h_radiaus->GetXaxis()->SetTitle("radiaus (mm)");
  h_radiaus->GetXaxis()->CenterTitle();
  h_radiaus->GetXaxis()->SetTitleOffset(1.25);
  h_radiaus->GetXaxis()->SetNdivisions(505);
  h_radiaus->GetXaxis()->SetRangeUser(2.5,50.5);

  h_radiaus->GetYaxis()->SetTitle("counts");
  h_radiaus->GetYaxis()->CenterTitle();
  h_radiaus->GetYaxis()->SetTitleOffset(1.25);
  h_radiaus->GetYaxis()->SetNdivisions(505);

  h_radiaus->Draw("hE");

  TF1 *f_gaus = new TF1("f_gaus","gaus",0.0,100.0);
  f_gaus->SetParameter(0,1.0);
  f_gaus->SetParameter(1,40.0);
  f_gaus->SetParameter(2,10.0);
  f_gaus->SetRange(2.5,50.5);

  h_radiaus->Fit("f_gaus","MNR");
  f_gaus->SetLineColor(2);
  f_gaus->SetLineWidth(2);
  f_gaus->SetLineStyle(2);
  f_gaus->Draw("l same");
  float mean = f_gaus->GetParameter(1);
  float sigma = f_gaus->GetParameter(2);

  string meanR = Form("R = %2.2f",mean);
  string deltaR = Form("#DeltaR = %2.2f",sigma);

  TLegend *leg = new TLegend(0.2,0.6,0.5,0.8);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(h_radiaus,"simulation","P");
  leg->AddEntry(f_gaus,meanR.c_str(),"l");
  leg->AddEntry(f_gaus,deltaR.c_str(),"l");
  leg->Draw("same");

}
