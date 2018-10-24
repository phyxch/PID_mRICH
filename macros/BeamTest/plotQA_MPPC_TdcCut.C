#include "string"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "../draw.h"

void plotQA_MPPC_TdcCut()
{
  string inputfile = "/home/xusun/Data/mRICH/BeamTest/QA/MPPC/sipmTimeCuts.root";
  TFile *File_InPut = TFile::Open(inputfile.c_str());

  TH2F *h_mTimeCuts = (TH2F*)File_InPut->Get("h_mTimeCuts")->Clone();
  TH1F *h_mTdcStart = (TH1F*)h_mTimeCuts->ProjectionY("h_mTdcStart",1,1)->Clone();
  TH1F *h_mTdcStop  = (TH1F*)h_mTimeCuts->ProjectionY("h_mTdcStop",2,2)->Clone();

  TCanvas *c_TimeCuts = new TCanvas("c_TimeCuts","c_TimeCuts",10,10,1200,400);
  c_TimeCuts->Divide(3,1);
  for(int i_pad = 0; i_pad < 3; ++i_pad)
  {
    c_TimeCuts->cd(i_pad+1);
    c_TimeCuts->cd(i_pad+1)->SetLeftMargin(0.15);
    c_TimeCuts->cd(i_pad+1)->SetBottomMargin(0.15);
    c_TimeCuts->cd(i_pad+1)->SetTicks(1,1);
    c_TimeCuts->cd(i_pad+1)->SetGrid(0,0);
  }

  c_TimeCuts->cd(1);
  h_mTimeCuts->SetTitle("PMT TDC Cuts");
  h_mTimeCuts->SetStats(0);
  h_mTimeCuts->Draw("colz");
  // PlotLine(-0.5,2.5,15,15,2,2,2);
  // PlotLine(-0.5,2.5,317,317,2,2,2);
  // plotTopLegend((char*)"proton",0.0,140,0.06,1,0.0,42,0,1);

  c_TimeCuts->cd(2);
  h_mTdcStart->SetTitle("proton tdc Start");
  h_mTdcStart->SetStats(0);
  h_mTdcStart->SetMarkerStyle(24);
  h_mTdcStart->SetMarkerSize(1.5);
  h_mTdcStart->SetMarkerColor(2);
  h_mTdcStart->GetYaxis()->SetTitle("tdc");
  h_mTdcStart->GetYaxis()->SetTitleSize(0.06);
  h_mTdcStart->GetYaxis()->CenterTitle();
  h_mTdcStart->GetYaxis()->SetRangeUser(490,590);
  h_mTdcStart->DrawCopy("p");

  c_TimeCuts->cd(3);
  h_mTdcStop->SetTitle("proton tdc Stop");
  h_mTdcStop->SetStats(0);
  h_mTdcStop->SetMarkerStyle(24);
  h_mTdcStop->SetMarkerSize(1.5);
  h_mTdcStop->SetMarkerColor(1);
  h_mTdcStop->GetYaxis()->SetTitle("tdc");
  h_mTdcStop->GetYaxis()->SetTitleSize(0.06);
  h_mTdcStop->GetYaxis()->CenterTitle();
  h_mTdcStop->GetYaxis()->SetRangeUser(490,590);
  h_mTdcStop->DrawCopy("p");

  c_TimeCuts->SaveAs("../../figures/BeamTest_QA/c_TdcCuts_MPPC.eps");
}
