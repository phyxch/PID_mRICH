#include "TString.h"
#include "TFile.h"
#include "TEllipse.h"
#include "TCanvas.h"

void plotQA_ringFinder()
{
  TFile *File_InPut = TFile::Open("../out.root");
  TH2D *h_mPhotonDist = (TH2D*)File_InPut->Get("h_mPhotonDist")->Clone();
  TH3D *h_mCherenkovRing = (TH3D*)File_InPut->Get("h_mCherenkovRing")->Clone();

  int hBin_x = -1;
  int hBin_y = -1;
  int hBin_r = -1;

  int globalBin = h_mHoughTransform->GetMaximumBin(hBin_x,hBin_y,hBin_r);
  int maxVote = h_mHoughTransform->GetBinContent(globalBin);
  double x_Cherenkov = h_mCherenkovRing->GetXaxis()->GetBinCenter(hBin_x);
  double y_Cherenkov = h_mCherenkovRing->GetYaxis()->GetBinCenter(hBin_y);
  double r_Cherenkov = h_mCherenkovRing->GetZaxis()->GetBinCenter(hBin_r);

  TEllipse *circle = new TEllipse(x_Cherenkov,y_Cherenkov,r_Cherenkov,r_Cherenkov);
  TCanvas *c_ringFinder = new TCanvas("c_ringFinder","c_ringFinder",10,10,800,800);
  c_ringFinder->SetLeftMargin(0.15);
  c_ringFinder->SetBottomMargin(0.15);
  c_ringFinder->SetTicks(1,1);
  c_ringFinder->SetGrid(0,0);
  h_mPhotonDist->Draw("colz");

  circle->SetFillColor(0);
  circle->SetFillStyle(0);
  circle->SetLineWidth(4);
  circle->SetLineColor(4);
  circle->Draw("same");

  c_ringFinder->SaveAs("../figures/c_ringFinder.eps");
}
