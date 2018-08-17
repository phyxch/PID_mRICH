#include "TString.h"
#include "TFile.h"
#include "TEllipse.h"
#include "TCanvas.h"

void plotQA_HoughTransform()
{
  TFile *File_InPut = TFile::Open("../out.root");
  TH3D *h_mCherenkovRing = (TH3D*)File_InPut->Get("h_mCherenkovRing")->Clone();

  TCanvas *c_CherenkovRing = new TCanvas("c_HoughTransform","c_HoughTransform",10,10,1600,800);
  c_CherenkovRing->Divide(2,1);
  for(int i_pad = 0; i_pad < 2; ++i_pad)
  {
    c_CherenkovRing->cd(i_pad+1)->SetLeftMargin(0.15);
    c_CherenkovRing->cd(i_pad+1)->SetBottomMargin(0.15);
    c_CherenkovRing->cd(i_pad+1)->SetTicks(1,1);
    c_CherenkovRing->cd(i_pad+1)->SetGrid(0,0);
  }
  c_CherenkovRing->cd(1);
  h_mCherenkovRing->Project3D("xy")->Draw("colz");
  c_CherenkovRing->cd(2);
  TH1D *h_radiaus = (TH1D*)h_mCherenkovRing->Project3D("z")->Clone();
  h_radiaus->Draw();
  h_radiaus->Fit("gaus");
  c_CherenkovRing->SaveAs("../figures/c_CherenkovRing.eps");

  TH2D *h_mNumOfCherenkovPhotons = (TH2D*)File_InPut->Get("h_mNumOfCherenkovPhotons")->Clone();
  TCanvas *c_NumOfPhotons = new TCanvas("c_NumOfPhotons","c_NumOfPhotons",10,10,1600,800);
  c_NumOfPhotons->Divide(2,1);
  for(int i_pad = 0; i_pad < 2; ++i_pad)
  {
    c_NumOfPhotons->cd(i_pad+1)->SetLeftMargin(0.15);
    c_NumOfPhotons->cd(i_pad+1)->SetBottomMargin(0.15);
    c_NumOfPhotons->cd(i_pad+1)->SetTicks(1,1);
    c_NumOfPhotons->cd(i_pad+1)->SetGrid(0,0);
  }
  c_NumOfPhotons->cd(1);
  h_mNumOfCherenkovPhotons->ProjectionX()->Draw();
  c_NumOfPhotons->cd(2);
  h_mNumOfCherenkovPhotons->ProjectionY()->Draw();
  c_NumOfPhotons->SaveAs("../figures/c_NumOfPhotons.eps");
}
