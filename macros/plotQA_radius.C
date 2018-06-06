#include "TF1.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TMath.h"
#include "TLegend.h"

double radius(double *x_val, double *par)
{
  double p = x_val[0];
  double m0 = par[0];
  double n = par[1];
  double f = par[2];

  double numerator = (n*n-1.0)*p*p-m0*m0;
  double denominator = (2.0-n*n)*p*p+m0*m0;

  double r = 0.0;
  if( numerator > 0) r = f*TMath::Sqrt(numerator/denominator);

  return r;
}

void plotQA_radius()
{
  double f = 6.0*25.4; // mm
  // double f = 1.0; // mm
  double n = 1.03; // refractive index
  double mass[4] = {0.138,0.494,0.938,0.00051}; // GeV
  TString ParType[4] = {"pion","Kaon","proton","electron"};
  int color[4] = {2,4,1,6};
  int style[4] = {2,4,1,6};

  TH1F *h_frame = new TH1F("h_frame","h_frame",170,-1.2,15.8);
  for(int i_bin = 0; i_bin < 160; ++i_bin)
  {
    h_frame->SetBinContent(i_bin+1,-100.0);
    h_frame->SetBinError(i_bin+1,1.0);
  }
  h_frame->SetTitle("");
  h_frame->SetStats(0);
  h_frame->GetXaxis()->SetNdivisions(505,'N');
  h_frame->GetXaxis()->SetLabelSize(0.03);
  h_frame->GetXaxis()->SetTitle("momentum (GeV/c)");
  h_frame->GetXaxis()->SetTitleSize(0.05);
  h_frame->GetXaxis()->SetTitleOffset(1.2);
  h_frame->GetXaxis()->CenterTitle();

  h_frame->GetYaxis()->SetRangeUser(0.0,44.0);
  h_frame->GetYaxis()->SetNdivisions(505,'N');
  h_frame->GetYaxis()->SetTitle("radius (mm)");
  h_frame->GetYaxis()->SetTitleSize(0.05);
  h_frame->GetYaxis()->SetLabelSize(0.03);
  h_frame->GetYaxis()->CenterTitle();

  TCanvas *c_play = new TCanvas("c_play","c_play",10,10,800,800);
  c_play->SetLeftMargin(0.15);
  c_play->SetBottomMargin(0.15);
  c_play->SetGrid(0,0);
  c_play->SetTicks(1,1);
  h_frame->DrawCopy("pE");

  TF1 *f_radius[4];
  for(int i_particle = 0; i_particle < 4; ++i_particle)
  {
    TString FuncName = Form("f_%s",ParType[i_particle].Data());
    double p_min = TMath::Sqrt(mass[i_particle]*mass[i_particle]/(n*n-1.0));
    // cout << "mass = " << mass[i_particle] << ", p_min = " << p_min<< endl;
    f_radius[i_particle] = new TF1(FuncName.Data(),radius,0,15.0,3);
    f_radius[i_particle]->SetNpx(10000);
    f_radius[i_particle]->FixParameter(0,mass[i_particle]);
    f_radius[i_particle]->FixParameter(1,n);
    f_radius[i_particle]->FixParameter(2,f);
    f_radius[i_particle]->SetLineColor(color[i_particle]);
    f_radius[i_particle]->SetLineStyle(style[i_particle]);
    f_radius[i_particle]->SetLineWidth(4);
    f_radius[i_particle]->SetRange(p_min*0.8,15.0);
    f_radius[i_particle]->Draw("l same");
  }

  TLegend *leg = new TLegend(0.5,0.4,0.85,0.6);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  for(int i_particle = 0; i_particle < 4; ++i_particle)
  {
    leg->AddEntry(f_radius[i_particle],ParType[i_particle].Data(),"L");
  }
  leg->Draw("same");

  c_play->SaveAs("../figures/c_radius.eps");
}
