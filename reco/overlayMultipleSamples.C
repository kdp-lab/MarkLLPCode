#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"
#include <TError.h>
#include "TSystem.h"

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

void SetPlotStyle();
void mySmallText(Double_t x, Double_t y, Color_t color, char* text);






void openRootFiles() {
    SetPlotStyle();
    // Create a TChain for each file
    TChain* StauMCPsChain_1000_005 = new TChain("StauMCPs");
    StauMCPsChain_1000_005->Add("1000_0.05_reco_TIMING32ns64ns.root");
    TChain* StauMCPsChain_1000_01 = new TChain("StauMCPs");
    StauMCPsChain_1000_01->Add("1000_0.1_reco_TIMING32ns64ns.root");
    TChain* StauMCPsChain_1000_1 = new TChain("StauMCPs");
    StauMCPsChain_1000_1->Add("1000_1_reco_TIMING32ns64ns.root");
    TChain* StauMCPsChain_4500_01 = new TChain("StauMCPs");
    StauMCPsChain_4500_01->Add("4500_0.1_reco_TIMING32ns64ns.root");
    TChain* StauMCPsChain_4500_1 = new TChain("StauMCPs");
    StauMCPsChain_4500_1->Add("4500_1_reco_TIMING32ns64ns.root");
    TChain* StauMCPsChain_4500_10 = new TChain("StauMCPs");
    StauMCPsChain_4500_10->Add("4500_10_reco_TIMING32ns64ns.root");
    std::vector<float >* prod_stau_endpoint_r_1000_005;
    TBranch* b_prod_stau_endpoint_r_1000_005;
    prod_stau_endpoint_r_1000_005 = 0;
    std::vector<float >* prod_stau_endpoint_r_1000_01;
    TBranch* b_prod_stau_endpoint_r_1000_01;
    prod_stau_endpoint_r_1000_01 = 0;
    std::vector<float >* prod_stau_endpoint_r_1000_1;
    TBranch* b_prod_stau_endpoint_r_1000_1;
    prod_stau_endpoint_r_1000_1 = 0;
    std::vector<float >* prod_stau_endpoint_r_4500_01;
    TBranch* b_prod_stau_endpoint_r_4500_01;
    prod_stau_endpoint_r_4500_01 = 0;
    std::vector<float >* prod_stau_endpoint_r_4500_1;
    TBranch* b_prod_stau_endpoint_r_4500_1;
    prod_stau_endpoint_r_4500_1 = 0;
    std::vector<float >* prod_stau_endpoint_r_4500_10;
    TBranch* b_prod_stau_endpoint_r_4500_10;
    prod_stau_endpoint_r_4500_10 = 0;

    StauMCPsChain_1000_005->SetBranchAddress("prod_stau_endpoint_r", &prod_stau_endpoint_r_1000_005, &b_prod_stau_endpoint_r_1000_005);
    StauMCPsChain_1000_01->SetBranchAddress("prod_stau_endpoint_r", &prod_stau_endpoint_r_1000_01, &b_prod_stau_endpoint_r_1000_01);
    StauMCPsChain_1000_1->SetBranchAddress("prod_stau_endpoint_r", &prod_stau_endpoint_r_1000_1, &b_prod_stau_endpoint_r_1000_1);
    StauMCPsChain_4500_01->SetBranchAddress("prod_stau_endpoint_r", &prod_stau_endpoint_r_4500_01, &b_prod_stau_endpoint_r_4500_01);
    StauMCPsChain_4500_1->SetBranchAddress("prod_stau_endpoint_r", &prod_stau_endpoint_r_4500_1, &b_prod_stau_endpoint_r_4500_1);
    StauMCPsChain_4500_10->SetBranchAddress("prod_stau_endpoint_r", &prod_stau_endpoint_r_4500_10, &b_prod_stau_endpoint_r_4500_10);

    // Create a canvas
    TCanvas* c1 = new TCanvas("c1", "Overlayed Histograms", 800, 600);
    c1->SetLogy();

    // Create histograms
    TH1F* h1000_005 = new TH1F("h1000_005", "Stau Decay Endpoint Radius", 100, 0, 2000); // Adjust bins and range as needed
    TH1F* h1000_01 = new TH1F("h1000_01", "Stau Decay Endpoint Radius", 100, 0, 2000);
    TH1F* h1000_1 = new TH1F("h1000_1", "Stau Decay Endpoint Radius", 100, 0, 2000);
    TH1F* h4500_01 = new TH1F("h4500_01", "Stau Decay Endpoint Radius", 100, 0, 2000);
    TH1F* h4500_1 = new TH1F("h4500_1", "Stau Decay Endpoint Radius", 100, 0, 2000);
    TH1F* h4500_10 = new TH1F("h4500_10", "Stau Decay Endpoint Radius", 100, 0, 2000);

    for (int iEvt = 0 ; iEvt < StauMCPsChain_1000_005->GetEntries(); ++iEvt){
        StauMCPsChain_1000_005->GetEntry(iEvt);
        StauMCPsChain_1000_01->GetEntry(iEvt);
        StauMCPsChain_1000_1->GetEntry(iEvt);
        StauMCPsChain_4500_01->GetEntry(iEvt);
        StauMCPsChain_4500_1->GetEntry(iEvt);
        StauMCPsChain_4500_10->GetEntry(iEvt);

        for (unsigned int iStau = 0; iStau < prod_stau_endpoint_r_1000_005->size(); ++iStau){
            h1000_005->Fill(prod_stau_endpoint_r_1000_005->at(iStau));
        }
        for (unsigned int iStau = 0; iStau < prod_stau_endpoint_r_1000_01->size(); ++iStau){
            h1000_01->Fill(prod_stau_endpoint_r_1000_01->at(iStau));
        }
        for (unsigned int iStau = 0; iStau < prod_stau_endpoint_r_1000_1->size(); ++iStau){
            h1000_1->Fill(prod_stau_endpoint_r_1000_1->at(iStau));
        }
        for (unsigned int iStau = 0; iStau < prod_stau_endpoint_r_4500_01->size(); ++iStau){
            h4500_01->Fill(prod_stau_endpoint_r_4500_01->at(iStau));
        }
        for (unsigned int iStau = 0; iStau < prod_stau_endpoint_r_4500_1->size(); ++iStau){
            h4500_1->Fill(prod_stau_endpoint_r_4500_1->at(iStau));
        }
        for (unsigned int iStau = 0; iStau < prod_stau_endpoint_r_4500_10->size(); ++iStau){
            h4500_10->Fill(prod_stau_endpoint_r_4500_10->at(iStau));
        }
    }

    // Set histogram styles
    h1000_005->SetLineColor(kRed);
    h1000_01->SetLineColor(kBlue);
    h1000_1->SetLineColor(kGreen);
    h4500_01->SetLineColor(kMagenta);
    h4500_1->SetLineColor(kCyan);
    h4500_10->SetLineColor(kOrange);

    h1000_005->SetLineWidth(2);
    h1000_01->SetLineWidth(2);
    h1000_1->SetLineWidth(2);
    h4500_01->SetLineWidth(2);
    h4500_1->SetLineWidth(2);
    h4500_10->SetLineWidth(2);

    // Remove stat box
    h1000_005->SetStats(0);
    h1000_005->SetMaximum(1800); // Adjust this value as needed
    h1000_01->SetStats(0);
    h1000_01->SetMaximum(1800); // Adjust this value as needed
    h1000_1->SetStats(0);
    h1000_1->SetMaximum(1800); // Adjust this value as needed
    h4500_01->SetStats(0);
    h4500_01->SetMaximum(1800); // Adjust this value as needed
    h4500_1->SetStats(0);
    h4500_1->SetMaximum(1800); // Adjust this value as needed
    h4500_10->SetStats(0);
    h4500_10->SetMaximum(1800); // Adjust this value as needed

    // Set axis labels
    h1000_005->GetXaxis()->SetTitle("Decay Length (R_{xy}) [mm]");
    h1000_005->GetYaxis()->SetTitle("Num. Staus");

    // Create a legend
    TLegend* legend = new TLegend(0.5, 0.6, 0.95, 0.95);  // Adjusted position
    legend->AddEntry(h1000_005, "m_{#tilde{#tau}}: 1 TeV, c#tau: 15 mm", "l");
    legend->AddEntry(h1000_01, "m_{#tilde{#tau}}: 1 TeV, c#tau: 30 mm", "l");
    legend->AddEntry(h1000_1, "m_{#tilde{#tau}}: 1 TeV, c#tau: 300 mm", "l");
    legend->AddEntry(h4500_01, "m_{#tilde{#tau}}: 4.5 TeV, c#tau: 30 mm", "l");
    legend->AddEntry(h4500_1, "m_{#tilde{#tau}}: 4.5 TeV, c#tau: 300 mm", "l");
    legend->AddEntry(h4500_10, "m_{#tilde{#tau}}: 4.5 TeV, c#tau: 3000 mm", "l");
    // Set text size
    legend->SetTextSize(0.0435);  // Adjust the value to your preference

    // Draw histograms
    h1000_005->Draw();
    h1000_01->Draw("SAME");
    h1000_1->Draw("SAME");
    h4500_01->Draw("SAME");
    h4500_1->Draw("SAME");
    h4500_10->Draw("SAME");

    // Draw legend
    legend->Draw();

    // Update canvas
    c1->Update();
    c1->SaveAs("overlayed_histograms_log.pdf"); // Save the canvas as an image file
}


















// Declare helper functions 
void SetPlotStyle() {
  // from ATLAS plot style macro

  // use plain black on white colors
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetHistLineColor(1);

  gStyle->SetPalette(1);

  // set the paper & margin sizes
  gStyle->SetPaperSize(20, 26);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadLeftMargin(0.16);

  // set title offsets (for axis label)
  gStyle->SetTitleXOffset(1.4);
  gStyle->SetTitleYOffset(1.4);

  // use large fonts
  gStyle->SetTextFont(42);
  gStyle->SetTextSize(0.05);
  gStyle->SetLabelFont(42, "x");
  gStyle->SetTitleFont(42, "x");
  gStyle->SetLabelFont(42, "y");
  gStyle->SetTitleFont(42, "y");
  gStyle->SetLabelFont(42, "z");
  gStyle->SetTitleFont(42, "z");
  gStyle->SetLabelSize(0.05, "x");
  gStyle->SetTitleSize(0.05, "x");
  gStyle->SetLabelSize(0.05, "y");
  gStyle->SetTitleSize(0.05, "y");
  gStyle->SetLabelSize(0.05, "z");
  gStyle->SetTitleSize(0.05, "z");

  // use bold lines and markers
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(2.);
  gStyle->SetHistLineWidth(3.);
  gStyle->SetLineStyleString(2, "[12 12]");

  // get rid of error bar caps
  gStyle->SetEndErrorSize(0.);

  // do not display any of the standard histogram decorations
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
}

void mySmallText(Double_t x, Double_t y, Color_t color, char* text) {
  Double_t tsize = 0.044;
  TLatex l;
  l.SetTextSize(tsize);
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x, y, text);
}
