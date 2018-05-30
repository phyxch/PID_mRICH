#include "TString.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGaxis.h"
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TColor.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TExec.h"
#include "TPolyMarker.h"
#include "TVirtualPad.h"
#include "TPolyLine.h"
#include "TVector3.h"
#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"
#include "TVirtualFitter.h"
#include "Math/MinimizerOptions.h"
#include "TLorentzVector.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"
#include "TGeoMatrix.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include <vector>
#include "TF1.h"
#include "TProfile.h"

//----------------------------------------------------------------------------------------
TLatex* plotTopLegend(char* label,Float_t x=-1,Float_t y=-1,Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1)
{
    // coordinates in NDC!
    // plots the string label in position x and y in NDC coordinates
    // size is the text size
    // color is the text color

    if((x<0||y<0) && NDC == 1)
    {   // defaults
      x=gPad->GetLeftMargin()*1.15;
      y=(1-gPad->GetTopMargin())*1.04;
    }
    TLatex* text=new TLatex(x,y,label);
    text->SetTextFont(font);
    text->SetTextSize(size);
    if(NDC == 1) text->SetNDC();
    text->SetTextColor(color);
    text->SetTextAngle(angle);
    text->SetTextAlign(align);
    text->Draw();
    return text;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void PlotLine(Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)
{
    TLine* Zero_line = new TLine();
    Zero_line -> SetX1(x1_val);
    Zero_line -> SetX2(x2_val);
    Zero_line -> SetY1(y1_val);
    Zero_line -> SetY2(y2_val);
    Zero_line -> SetLineWidth(LineWidth);
    Zero_line -> SetLineStyle(LineStyle);
    Zero_line -> SetLineColor(Line_Col);
    Zero_line -> Draw();
    //delete Zero_line;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void PlotHistLine2(TH1F* Histo, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle, Float_t x_start, Float_t x_stop, Float_t y_min, Float_t y_max)
{
    TPolyLine* Hist_line = new TPolyLine();
    Hist_line -> SetLineWidth(LineWidth);
    Hist_line -> SetLineStyle(LineStyle);
    Hist_line -> SetLineColor(Line_Col);
    for(Int_t i = 1; i < Histo->GetNbinsX(); i++)
    {
        Double_t x = Histo->GetBinCenter(i);
        Double_t y = Histo->GetBinContent(i);
        if(y >= y_min && y < y_max && y != 0 && x >= x_start && x <= x_stop)
        {
            //cout << "x = " << x << ", y = " << y << endl;
            Hist_line->SetNextPoint(x,y);
        }
    }
    Hist_line -> Draw();
    //delete Hist_line;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void PlotHistErrorBand(TH1F* Histo, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle, Int_t FillStyle, Int_t FillColor, Float_t x_start, Float_t x_stop)
{
    // Modified March 14th
    Int_t N_points = 0;
    for(Int_t i = 1; i < Histo->GetNbinsX(); i++)
    {
        Double_t x = Histo->GetBinCenter(i);
        Double_t y = Histo->GetBinContent(i);
        if(y != 0 && x >= x_start && x <= x_stop)
        {
            N_points++;
        }
    }
    Int_t N_total_points = N_points*2+1;
    //cout << "N_total_points = " << N_total_points << endl;
    TGraph* Hist_line = new TGraph(N_total_points);
    Hist_line -> SetLineWidth(LineWidth);
    Hist_line -> SetLineStyle(LineStyle);
    Hist_line -> SetLineColor(Line_Col);
    Hist_line -> SetFillStyle(FillStyle);
    Hist_line -> SetFillColor(FillColor);
    Int_t N_point = 0;
    for(Int_t i = 1; i < Histo->GetNbinsX(); i++)
    {
        Double_t y     = Histo->GetBinContent(i);
        Double_t y_err = Histo->GetBinError(i);
        if(y != 0)
        {
            Double_t x = Histo->GetBinCenter(i);
            if(x >= x_start && x <= x_stop)
            {
                //cout << "N_point = " << N_point << ", x = " << x << ", y = " << y << ", y_err = " << y_err << endl;
                Hist_line->SetPoint(N_point,x,y-y_err);
                Hist_line->SetPoint(N_total_points-2-N_point,x,y+y_err);
                if(N_point == 0) Hist_line->SetPoint(N_total_points-1,x,y-y_err);
                N_point++;
            }
        }
    }
    Hist_line -> Draw("f");
    //delete Hist_line;
}
//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
Int_t Hist_interpolate_and_error(TH1F* hist, Double_t x, Double_t &Int_val, Double_t &Int_err)
{
    // Linear interpolation of a histogram
    // Be careful, extrapolation is not done! -> closest data point is returned
    // Error calculation included
    // x-bin error is not taken into account

    Int_t return_val        = 0;
    Double_t bin_entries[2] = {0,0};
    Double_t bin_error[2]   = {0,0};
    Double_t bin_x_val[2]   = {0,0};
    Int_t binx_high;
    Int_t flag_max          = 0;

    if(hist->GetEntries() == 0) // no data poins to interpolate
    {
        return_val = -1;
        Int_val = 0;
        Int_err = 0;
        //cout << "No entries in histogram" << endl;
        return 0;
    }
    else
    {
        Int_t err_counter = 0;
        for(Int_t binx = 1; binx < hist->GetNbinsX(); binx++)
        {
            Double_t bin_error_val   = hist->GetBinError(binx);
            Double_t bin_x_pos       = hist->GetBinCenter(binx);
            if(bin_error_val != 0)
            {
                err_counter++;
                bin_entries[1] = hist->GetBinContent(binx);
                bin_error[1]   = hist->GetBinError(binx);
                bin_x_val[1]   = hist->GetBinCenter(binx);
                if(bin_x_pos >= x)
                {
                    binx_high = binx;
                    flag_max  = 1;
                    break;
                }
                else flag_max = 0;
            }
        }
        if(err_counter == 1 || flag_max == 0) // There is no lower/uppper data point -> extrapolation, return closest values
        {
            return_val = 0;
            Int_val = bin_entries[1];
            Int_err = bin_error[1];
            return return_val;
        }
        for(Int_t binx_low = binx_high; binx_low > 0; binx_low--)
        {
            bin_entries[0] = hist->GetBinContent(binx_low);
            bin_error[0]   = hist->GetBinError(binx_low);
            bin_x_val[0]   = hist->GetBinCenter(binx_low);
            if(bin_x_val[0] < x && bin_error[0] != 0)
            {
                break;
            }
        }

        if(bin_error[0] != 0 && bin_error[1] != 0)
        {
            return_val = 1;
            if(bin_x_val[0] != bin_x_val[1])
            {
                Double_t slope = (bin_entries[1]-bin_entries[0])/(bin_x_val[1]-bin_x_val[0]);
                Double_t t_val = bin_entries[1]-slope*bin_x_val[1];
                Int_val        = slope*x+t_val;

                Double_t x1     = bin_x_val[0];
                Double_t x2     = bin_x_val[1];
                Double_t y1     = bin_entries[0];
                Double_t y2     = bin_entries[1];
                Double_t x1_err = 0.0;
                Double_t x2_err = 0.0;
                Double_t y1_err = bin_error[0];
                Double_t y2_err = bin_error[1];

                Double_t termA  = ((y2-y1)/((x2-x1)*(x2-x1)))*(x-x2);
                Double_t termB  = -((y2-y1)/((x2-x1)*(x2-x1)))*(x-x2) - (y2-y1)/(x2-x1);
                Double_t termC  = -(x-x2)/(x2-x1);
                Double_t termD  = (x-x2)/(x2-x1) + 1.0;
                termA *= x1_err;
                termB *= x2_err;
                termC *= y1_err;
                termD *= y2_err;
                Int_err = TMath::Sqrt(termA*termA+termB*termB+termC*termC+termD*termD);
                return return_val;
            }
            else
            {
                Int_val = bin_entries[0];
                Int_err = bin_error[0];
                return return_val;
            }
        }
    }

    return return_val;

}
//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
Int_t TGraphAsymmErrors_interpolate_and_error(TGraphAsymmErrors* hist, Double_t x, Double_t &Int_val, Double_t &Int_err_low,Double_t &Int_err_high)
{
    // V2: 03.12.2012 -> bug fixed to calculate error bars
    // Linear interpolation of a TGraphAsymmErrors
    // Be careful, extrapolation is not done! -> closest data point is returned
    // Error calculation included
    // x-bin error is not taken into account

    Int_t return_val        = 0;
    Double_t bin_entries[2] = {0,0};
    Double_t bin_error[4]   = {0,0,0,0};
    Double_t bin_x_val[2]   = {0,0};
    Int_t binx_high;
    Int_t flag_max          = 0;

    if(hist->GetN() == 0) // no data poins to interpolate
    {
        return_val   = -1;
        Int_val      = 0;
        Int_err_low  = 0;
        Int_err_high = 0;
        //cout << "No entries in TGraphAsymmErrors" << endl;
        return 0;
    }
    else
    {
        Int_t err_counter = 0;
        for(Int_t epoint = 0; epoint < hist->GetN(); epoint++)
        {
            hist->GetPoint(epoint,bin_x_val[1],bin_entries[1]);
            bin_error[2] = hist->GetErrorYlow(epoint);
            bin_error[3] = hist->GetErrorYhigh(epoint);

            err_counter++;
            if(bin_x_val[1] >= x)
            {
                binx_high = epoint;
                flag_max  = 1;
                break;
            }
            else flag_max = 0;
        }
        if(err_counter == 1 || flag_max == 0) // There is no lower/uppper data point -> extrapolation, return closest values
        {
            return_val = 0;
            Int_val      = bin_entries[1];
            Int_err_low  = bin_error[2];
            Int_err_high = bin_error[3];
            return return_val;
        }
        for(Int_t epoint = binx_high; epoint >= 0; epoint--)
        {
            hist->GetPoint(epoint,bin_x_val[0],bin_entries[0]);
            bin_error[0] = hist->GetErrorYlow(epoint);
            bin_error[1] = hist->GetErrorYhigh(epoint);

            if(bin_x_val[0] < x && bin_error[0] != 0)
            {
                break;
            }
        }

        //cout << "bin0 = " << bin_error[0] << ", bin2 = " << bin_error[2] << endl;

        if(bin_error[0] != 0 && bin_error[2] != 0)
        {
            return_val = 1;
            if(bin_x_val[0] != bin_x_val[1])
            {
                Double_t slope = (bin_entries[1]-bin_entries[0])/(bin_x_val[1]-bin_x_val[0]);
                Double_t t_val = bin_entries[1]-slope*bin_x_val[1];
                Int_val        = slope*x+t_val;

                Double_t x1     = bin_x_val[0];
                Double_t x2     = bin_x_val[1];
                Double_t y1     = bin_entries[0];
                Double_t y2     = bin_entries[1];
                Double_t x1_err = 0.0;
                Double_t x2_err = 0.0;

                //cout << "x1 = " << x1 << ", x2 = " << x2 << ", y1 = " << y1 << ", y2 = " << y2 << endl;

                Double_t y1_err = bin_error[0];
                Double_t y2_err = bin_error[2];

                Double_t termA  = ((y2-y1)/((x2-x1)*(x2-x1)))*(x-x2);
                Double_t termB  = -((y2-y1)/((x2-x1)*(x2-x1)))*(x-x2) - (y2-y1)/(x2-x1);
                Double_t termC  = -(x-x2)/(x2-x1);
                Double_t termD  = (x-x2)/(x2-x1) + 1.0;

                termA *= x1_err;
                termB *= x2_err;
                termC *= y1_err;
                termD *= y2_err;
                Int_err_low = TMath::Sqrt(termA*termA+termB*termB+termC*termC+termD*termD);

                y1_err = bin_error[1];
                y2_err = bin_error[3];

                termA  = ((y2-y1)/((x2-x1)*(x2-x1)))*(x-x2);
                termB  = -((y2-y1)/((x2-x1)*(x2-x1)))*(x-x2) - (y2-y1)/(x2-x1);
                termC  = -(x-x2)/(x2-x1);
                termD  = (x-x2)/(x2-x1) + 1.0;

                termA *= x1_err;
                termB *= x2_err;
                termC *= y1_err;
                termD *= y2_err;
                Int_err_high = TMath::Sqrt(termA*termA+termB*termB+termC*termC+termD*termD);
                return return_val;
            }
            else
            {
                Int_val      = bin_entries[0];
                Int_err_low  = bin_error[0];
                Int_err_high = bin_error[1];
                return return_val;
            }
        }
    }

    return return_val;

}
//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
void Draw_TGAE_new_Symbol(TGraphAsymmErrors* tgae, Int_t style, Int_t color, Float_t size)
{
    TString HistName;
    Float_t size_A = 1.35*size;
    Float_t size_B = size;
    Float_t size_C = size;
    Int_t alt_marker = style;
    Int_t style_in = style;
    if(style == 24)
    {
        alt_marker = 20;
        size_A = 1.35*size;
    }
    if(style == 25)
    {
        alt_marker = 21;
        size_A = 1.35*size;
    }
    if(style == 26)
    {
        alt_marker = 22;
        size_A = 1.5*size;
    }
    if(style == 23)
    {
        alt_marker = 23;
        size_A = 1.35*size;
    }
    if(style == 30 || style == 29)
    {
        alt_marker = 29;
        size_A = 1.55*size;
    }
    if(style == 260)
    {
        alt_marker = 26;
        size_A = 1.15*size;
        style = 26;
    }
    if(style == 300)
    {
        alt_marker = 30;
        size_A = 1.3*size;
        style = 30;
    }

    // black and filled outer marker
    HistName = "tgae_dummy_A";
    TGraphAsymmErrors* ge_clone_A = (TGraphAsymmErrors*)tgae->Clone(HistName.Data());
    ge_clone_A->SetMarkerSize(size_A);
    ge_clone_A->SetMarkerStyle(alt_marker);
    ge_clone_A->SetMarkerColor(1);
    ge_clone_A->SetLineColor(10);
    if(style_in == 260 || style_in == 300) ge_clone_A->SetLineColor(1);
    ge_clone_A->Draw("same PZ0");

    if(!(style_in == 260 || style_in == 300))
    {
        // white and filled inner marker
        HistName = "tgae_dummy_B";
        TGraphAsymmErrors* ge_clone_B = (TGraphAsymmErrors*)tgae->Clone(HistName.Data());
        ge_clone_B->SetMarkerSize(size_B);
        ge_clone_B->SetMarkerStyle(alt_marker);
        ge_clone_B->SetMarkerColor(10);
        ge_clone_B->SetLineColor(10);
        ge_clone_B->Draw("same PZ0");
    }

    // color inner marker
    HistName = "tgae_dummy_C";
    TGraphAsymmErrors* ge_clone_C = (TGraphAsymmErrors*)tgae->Clone(HistName.Data());
    ge_clone_C->SetMarkerSize(size_C);
    ge_clone_C->SetMarkerStyle(style);
    ge_clone_C->SetMarkerColor(color);
    ge_clone_C->SetLineColor(1);
    ge_clone_C->Draw("same PE1");
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Draw_TGAE_Point_new_Symbol(Double_t x_val, Double_t y_val, Double_t x_min_err, Double_t x_max_err,
                                Double_t y_min_err, Double_t y_max_err,
                                Int_t style, Int_t color, Float_t size)
{
    TGraphAsymmErrors* tgae = new TGraphAsymmErrors();
    tgae->SetPoint(0,x_val,y_val);
    tgae->SetPointError(0,x_min_err,x_max_err,y_min_err,y_max_err);

    TString HistName;
    Float_t size_A = 1.35*size;
    Float_t size_B = size;
    Float_t size_C = size;
    Int_t alt_marker = style;
    Int_t style_in = style;
    if(style == 24)
    {
        alt_marker = 20;
        size_A = 1.35*size;
    }
    if(style == 25)
    {
        alt_marker = 21;
        size_A = 1.35*size;
    }
    if(style == 26)
    {
        alt_marker = 22;
        size_A = 1.5*size;
    }
    if(style == 23)
    {
        alt_marker = 23;
        size_A = 1.35*size;
    }
    if(style == 30 || style == 29)
    {
        alt_marker = 29;
        size_A = 1.55*size;
    }
    if(style == 260)
    {
        alt_marker = 26;
        size_A = 1.15*size;
        style = 26;
    }
    if(style == 300)
    {
        alt_marker = 30;
        size_A = 1.3*size;
        style = 30;
    }

    // black and filled outer marker
    HistName = "tgae_dummy_A";
    TGraphAsymmErrors* ge_clone_A = (TGraphAsymmErrors*)tgae->Clone(HistName.Data());
    ge_clone_A->SetMarkerSize(size_A);
    ge_clone_A->SetMarkerStyle(alt_marker);
    ge_clone_A->SetMarkerColor(1);
    ge_clone_A->SetLineColor(10);
    if(style_in == 260 || style_in == 300) ge_clone_A->SetLineColor(1);
    ge_clone_A->Draw("same PZ0");

    if(!(style_in == 260 || style_in == 300))
    {
        // white and filled inner marker
        HistName = "tgae_dummy_B";
        TGraphAsymmErrors* ge_clone_B = (TGraphAsymmErrors*)tgae->Clone(HistName.Data());
        ge_clone_B->SetMarkerSize(size_B);
        ge_clone_B->SetMarkerStyle(alt_marker);
        ge_clone_B->SetMarkerColor(10);
        ge_clone_B->SetLineColor(10);
        ge_clone_B->Draw("same PZ0");
    }

    // color inner marker
    HistName = "tgae_dummy_C";
    TGraphAsymmErrors* ge_clone_C = (TGraphAsymmErrors*)tgae->Clone(HistName.Data());
    ge_clone_C->SetMarkerSize(size_C);
    ge_clone_C->SetMarkerStyle(style);
    ge_clone_C->SetMarkerColor(color);
    ge_clone_C->SetLineColor(1);
    ge_clone_C->Draw("same PZ0");
}
//----------------------------------------------------------------------------------------
