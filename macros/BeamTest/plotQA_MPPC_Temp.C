#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
using namespace std;

void plotQA_MPPC_Temp()
{
  string input_30 = "/home/xusun/Data/mRICH/BeamTest/QA/MPPC/sipmTDC_run649.root";
  TFile *File_30 = TFile::Open(input_30.c_str());
  assert(File_30);
  TH1F *h_tdc_30 = (TH1F*)File_30->Get("h_mTDC_pixelX_23_pixelY_7");

  string input_20 = "/home/xusun/Data/mRICH/BeamTest/QA/MPPC/sipmTDC_run639.root";
  TFile *File_20 = TFile::Open(input_20.c_str());
  assert(File_20);
  TH1F *h_tdc_20 = (TH1F*)File_20->Get("h_mTDC_pixelX_23_pixelY_7");

  string input_10 = "/home/xusun/Data/mRICH/BeamTest/QA/MPPC/sipmTDC_run674.root";
  TFile *File_10 = TFile::Open(input_10.c_str());
  assert(File_10);
  TH1F *h_tdc_10 = (TH1F*)File_10->Get("h_mTDC_pixelX_23_pixelY_7");

  string input_0 = "/home/xusun/Data/mRICH/BeamTest/QA/MPPC/sipmTDC_run686.root";
  TFile *File_0 = TFile::Open(input_0.c_str());
  assert(File_0);
  TH1F *h_tdc_0 = (TH1F*)File_0->Get("h_mTDC_pixelX_23_pixelY_7");

  string input_room = "/home/xusun/Data/mRICH/BeamTest/QA/MPPC/sipmTDC_run697.root";
  TFile *File_room  = TFile::Open(input_room.c_str());
  assert(File_room);
  TH1F *h_tdc_room = (TH1F*)File_room->Get("h_mTDC_pixelX_23_pixelY_7");


  TCanvas *c_MPPC_temp = new TCanvas("c_MPPC_temp","c_MPPC_temp",10,10,800,800);
  c_MPPC_temp->SetLeftMargin(0.15);
  c_MPPC_temp->SetBottomMargin(0.15);
  c_MPPC_temp->SetRightMargin(0.15);
  c_MPPC_temp->SetTicks(1,1);
  c_MPPC_temp->SetGrid(0,0);

  h_tdc_30->SetTitle("BV 65V & Threshold 50");
  h_tdc_30->SetStats(0);
  h_tdc_30->GetXaxis()->SetTitle("tdc");
  h_tdc_30->GetXaxis()->CenterTitle();
  h_tdc_30->GetXaxis()->SetRangeUser(500,650);
  h_tdc_30->GetYaxis()->SetTitle("counts");
  h_tdc_30->GetYaxis()->CenterTitle();
  h_tdc_30->SetLineColor(1);
  h_tdc_30->Draw();

  h_tdc_room->SetLineColor(kOrange+6);
  h_tdc_room->SetLineWidth(3);
  h_tdc_room->Draw("same");

  h_tdc_20->SetLineColor(2);
  h_tdc_20->Draw("same");

  h_tdc_10->SetLineColor(4);
  h_tdc_10->Draw("same");

  h_tdc_0->SetLineColor(6);
  h_tdc_0->Draw("same");

  TLegend *leg = new TLegend(0.6,0.5,0.9,0.65);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(h_tdc_30," -30","L");
  leg->AddEntry(h_tdc_20," -20","L");
  leg->AddEntry(h_tdc_10," -10","L");
  leg->AddEntry(h_tdc_0," 0","L");
  leg->AddEntry(h_tdc_room," room temp","L");
  leg->Draw("same");

  string c_mppc_temp = "../../figures/BeamTest_QA/c_TDC_MPPC_temp.eps";
  c_MPPC_temp->SaveAs(c_mppc_temp.c_str());


  TCanvas *c_MPPC_temp_div = new TCanvas("c_MPPC_temp_div","c_MPPC_temp_div",10,10,1600,400);
  c_MPPC_temp_div->Divide(4,1,0,0);
  for(int i_pad = 0; i_pad < 4; ++i_pad)
  {
    c_MPPC_temp_div->cd(i_pad+1)->SetTopMargin(0.1);
    c_MPPC_temp_div->cd(i_pad+1)->SetBottomMargin(0.15);
    c_MPPC_temp_div->cd(i_pad+1)->SetTicks(1,1);
    c_MPPC_temp_div->cd(i_pad+1)->SetGrid(0,0);
  }

  c_MPPC_temp_div->cd(1);
  h_tdc_30->SetTitle("BV 65V & Threshold 50");
  h_tdc_30->SetStats(0);
  h_tdc_30->GetXaxis()->SetTitle("tdc");
  h_tdc_30->GetXaxis()->CenterTitle();
  h_tdc_30->GetXaxis()->SetRangeUser(500,650);
  h_tdc_30->GetYaxis()->SetRangeUser(0,300);
  h_tdc_30->SetLineColor(1);
  h_tdc_30->Draw();

  h_tdc_room->SetLineColor(kOrange+6);
  h_tdc_room->SetLineWidth(3);
  h_tdc_room->Draw("same");

  TLegend *leg_30 = new TLegend(0.6,0.5,0.85,0.7);
  leg_30->SetFillColor(0);
  leg_30->SetBorderSize(0);
  leg_30->AddEntry(h_tdc_30,"  -30","L");
  leg_30->AddEntry(h_tdc_room,"room temp","L");
  leg_30->Draw("same");

  c_MPPC_temp_div->cd(2);
  h_tdc_20->SetTitle("BV 65V & Threshold 50");
  h_tdc_20->SetStats(0);
  h_tdc_20->GetXaxis()->SetTitle("tdc");
  h_tdc_20->GetXaxis()->CenterTitle();
  h_tdc_20->GetXaxis()->SetRangeUser(500,650);
  h_tdc_20->GetYaxis()->SetRangeUser(0,300);
  h_tdc_20->SetLineColor(2);
  h_tdc_20->Draw();

  TLegend *leg_20 = new TLegend(0.6,0.5,0.85,0.6);
  leg_20->SetFillColor(0);
  leg_20->SetBorderSize(0);
  leg_20->AddEntry(h_tdc_20,"  -20","L");
  leg_20->Draw("same");

  c_MPPC_temp_div->cd(3);
  h_tdc_10->SetTitle("BV 65V & Threshold 50");
  h_tdc_10->SetStats(0);
  h_tdc_10->GetXaxis()->SetTitle("tdc");
  h_tdc_10->GetXaxis()->CenterTitle();
  h_tdc_10->GetXaxis()->SetRangeUser(500,650);
  h_tdc_10->GetYaxis()->SetRangeUser(0,300);
  h_tdc_10->SetLineColor(4);
  h_tdc_10->Draw();

  TLegend *leg_10 = new TLegend(0.6,0.5,0.85,0.6);
  leg_10->SetFillColor(0);
  leg_10->SetBorderSize(0);
  leg_10->AddEntry(h_tdc_10,"  -10","L");
  leg_10->Draw("same");

  c_MPPC_temp_div->cd(4);
  h_tdc_0->SetTitle("BV 65V & Threshold 50");
  h_tdc_0->SetStats(0);
  h_tdc_0->GetXaxis()->SetTitle("tdc");
  h_tdc_0->GetXaxis()->CenterTitle();
  h_tdc_0->GetXaxis()->SetRangeUser(500,650);
  h_tdc_0->GetYaxis()->SetRangeUser(0,300);
  h_tdc_0->SetLineColor(6);
  h_tdc_0->Draw();

  TLegend *leg_0 = new TLegend(0.6,0.5,0.85,0.6);
  leg_0->SetFillColor(0);
  leg_0->SetBorderSize(0);
  leg_0->AddEntry(h_tdc_0,"  0","L");
  leg_0->Draw("same");

  string c_mppc_temp_div = "../../figures/BeamTest_QA/c_TDC_MPPC_temp_div.eps";
  c_MPPC_temp_div->SaveAs(c_mppc_temp_div.c_str());
}
