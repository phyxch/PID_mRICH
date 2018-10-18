#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
using namespace std;

void plotQA_Timing(const int PMTid = 182, const int SiPMid = 649)
{
  string inputPMT = Form("/home/xusun/Data/mRICH/BeamTest/tdc/richTDC_run%d/sspRich.root",PMTid);
  cout << "read in file: " << inputPMT.c_str() << endl;
  TFile *File_InPutPMT = TFile::Open(inputPMT.c_str());
  assert(File_InPutPMT);

  TCanvas *c_PMT = new TCanvas("c_PMT","c_PMT",10,10,800,800);
  c_PMT->cd();
  c_PMT->cd()->SetLeftMargin(0.15);
  c_PMT->cd()->SetBottomMargin(0.15);
  c_PMT->cd()->SetGrid(0,0);
  c_PMT->cd()->SetTicks(1,1);

  TH1F *h_PMT = new TH1F("h_PMT","h_PMT",200,1900.5,2100.5);
  TTree * tree_mRICH_PMT = (TTree*)File_InPutPMT->Get("data");
  tree_mRICH_PMT->Draw("time>>h_PMT","pol==1");

  string inputSiPM = Form("/home/xusun/Data/mRICH/BeamTest/tdc/sipmTDC_run%d/sspRich.root",SiPMid);
  cout << "read in file: " << inputSiPM.c_str() << endl;
  TFile *File_InPutSiPM = TFile::Open(inputSiPM.c_str());
  assert(File_InPutSiPM);

  TH1F *h_SiPM = new TH1F("h_SiPM","h_SiPM",200,1900.5,2100.5);
  h_SiPM->SetLineColor(2);
  TTree * tree_mRICH_SiPM = (TTree*)File_InPutSiPM->Get("data");
  tree_mRICH_SiPM->SetAlias("shift",Form("1.0*%f",1490.0));
  tree_mRICH_SiPM->SetAlias("TimeShift","time+shift");
  tree_mRICH_SiPM->Draw("TimeShift>>h_SiPM","pol==1","same");
  // tree_mRICH_SiPM->Draw("time");

  TLegend *leg = new TLegend(0.2,0.5,0.5,0.65);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(h_PMT,"PMT time","L");
  leg->AddEntry(h_SiPM,"SiPM time + 1490","L");
  leg->Draw("same");
}
