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
#include <cmath>

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;


//TString massPoint = "1TeV";
//TString massPoint = "2.5TeV";
//TString massPoint = "4TeV";
//TString massPoint = "4.5TeV";
//TString massPoint = "misc"; // For multiple mass points, custom defined filenames

std::map<TString, TH1F*> displaced_track_pt_map;
std::map<TString, TH1F*> displaced_track_d0_map;
std::map<TString, TH1F*> displaced_track_z0_map;
std::map<TString, TH1F*> displaced_track_eta_map;
std::map<TString, TH1F*> displaced_resz0_map;
std::map<TString, TH1F*> displaced_resd0_map;
std::map<TString, TH1F*> displaced_chi2_reduced_map;
std::map<TString, TH1F*> displaced_nhits_map;
std::map<TString, TH1F*> displaced_matched_d0_map;
std::map<TString, TH1F*> displaced_matched_rxy_map;
std::map<TString, TH1F*> displaced_matched_eta_map;
std::map<TString, TH1F*> displaced_tp_d0_map;
std::map<TString, TH1F*> displaced_tp_rxy_map;
std::map<TString, TH1F*> displaced_tp_eta_map;
std::map<TString, TH1F*> displaced_mcp_d0_map;
std::map<TString, TH1F*> displaced_mcp_rxy_map;
std::map<TString, TH1F*> displaced_mcp_eta_map;
std::map<TString, TH1F*> displaced_mcp_rxy_map_clone;
std::map<TString, TH1F*> displaced_tp_pt_map;
std::map<TString, TH1F*> displaced_matched_pt_map;
std::map<TString, TH1F*> displaced_mcp_pt_map;
std::map<TString, TH1F*> displaced_eff_acc_d0_map;
std::map<TString, TH1F*> displaced_eff_d0_map;
std::map<TString, TH1F*> displaced_eff_acc_rxy_map;
std::map<TString, TH1F*> displaced_eff_rxy_map;
std::map<TString, TH1F*> displaced_eff_acc_pt_map;
std::map<TString, TH1F*> displaced_eff_pt_map;
std::map<TString, TH1F*> displaced_eff_acc_z0_map; // FIXME make eff vs z0, eta plots
std::map<TString, TH1F*> displaced_eff_z0_map;
std::map<TString, TH1F*> displaced_eff_acc_eta_map;
std::map<TString, TH1F*> displaced_eff_eta_map;

std::map<TString, TH1F*> stau_mcp_pt_map;
std::map<TString, TH1F*> stau_mcp_eta_map;
std::map<TString, TH1F*> stau_mcp_phi_map;
std::map<TString, TH1F*> stau_mcp_d0_map;
std::map<TString, TH1F*> stau_mcp_z0_map;
std::map<TString, TH1F*> stau_chi2_reduced_map;

std::map<TString, TH1F*> stau_tp_pt_map;
std::map<TString, TH1F*> stau_tp_eta_map;
std::map<TString, TH1F*> stau_tp_phi_map;
std::map<TString, TH1F*> stau_tp_d0_map;
std::map<TString, TH1F*> stau_tp_z0_map;
std::map<TString, TH1F*> stau_matched_pt_map;
std::map<TString, TH1F*> stau_velores_map;
std::map<TString, TH1F*> stau_velores_nhits_map;
std::map<TString, TH1F*> stau_matched_pt_all_map;
std::map<TString, TH1F*> stau_track_pt_map;
std::map<TString, TH1F*> stau_nhits_map;

std::map<TString, TH1F*> stau_matched_phi_map;
std::map<TString, TH1F*> stau_matched_eta_map;
std::map<TString, TH1F*> stau_track_eta_map;

std::map<TString, TH1F*> stau_eff_acc_pt_map;
std::map<TString, TH1F*> stau_eff_acc_phi_map;
std::map<TString, TH1F*> stau_eff_pt_map;
std::map<TString, TH1F*> stau_eff_acc_z0_map; // FIXME make eff vs z0, eta plots
std::map<TString, TH1F*> dstau_eff_z0_map;
std::map<TString, TH1F*> stau_eff_acc_eta_map;
std::map<TString, TH1F*> stau_eff_eta_map;

std::map<TString, TH1F*> stau_rxy_decay_map;
std::map<TString, TH1F*> stau_velo_residuals_map;

std::map<TString, TH1F*> fake_track_pt_map;
std::map<TString, TH1F*> fake_track_eta_map;
std::map<TString, TH1F*> fake_track_d0_map;
std::map<TString, TH1F*> fake_track_z0_map;
std::map<TString, TH1F*> fake_track_chi2_reduced_map;
std::map<TString, TH1F*> fake_track_nhits_map;

// For overlaying displaced tracks (wtih same binning) with fake track distributions
std::map<TString, TH1F*> displaced_overlay_track_nhits_map;
std::map<TString, TH1F*> displaced_overlay_track_pt_map;
std::map<TString, TH1F*> displaced_overlay_track_eta_map;
std::map<TString, TH1F*> displaced_overlay_track_d0_map;
std::map<TString, TH1F*> displaced_overlay_track_z0_map;
std::map<TString, TH1F*> displaced_overlay_track_chi2_reduced_map;

  std::map<TString, TString> legend_map = {
    {"1000_0.05", "1 TeV, 15 mm"},
    {"1000_0.1", "1 TeV, 30 mm"},
    {"1000_1", "1 TeV, 300 mm"},
    {"1000_10", "1 TeV, 3000 mm"},
    {"2500_0.1", "2.5 TeV, 30 mm"},
    {"2500_1", "2.5 TeV, 300 mm"},
    {"2500_10", "2.5 TeV, 3000 mm"},
    {"4000_0.1", "4 TeV, 30 mm"},
    {"4000_1", "4 TeV, 300 mm"},
    {"4000_10", "4 TeV, 3000 mm"},
    {"4500_0.1", "4.5 TeV, 30 mm"},
    {"4500_1", "4.5 TeV, 300 mm"},
    {"4500_10", "4.5 TeV, 3000 mm"},
    {"1000_0.05_osgcomparison", "1 TeV, 15 mm, v2.9"},
    {"1000_0.1_osgcomparison", "1 TeV, 30 mm, v2.9"},
    {"1000_1_osgcomparison", "1 TeV, 300 mm, v2.9"},
    {"1000_10_osgcomparison", "1 TeV, 3000 mm, v2.9"},
    {"2500_0.1_osgcomparison", "2.5 TeV, 30 mm, v2.9"},
    {"2500_1_osgcomparison", "2.5 TeV, 300 mm, v2.9"},
    {"2500_10_osgcomparison", "2.5 TeV, 3000 mm, v2.9"},
    {"4000_0.1_osgcomparison", "4 TeV, 30 mm, v2.9"},
    {"4000_1_osgcomparison", "4 TeV, 300 mm, v2.9"},
    {"4000_10_osgcomparison", "4 TeV, 3000 mm, v2.9"},
    {"4500_0.1_osgcomparison", "4.5 TeV, 30 mm, v2.9"},
    {"4500_1_osgcomparison", "4.5 TeV, 300 mm, v2.9"},
    {"4500_10_osgcomparison", "4.5 TeV, 3000 mm, v2.9"},
    {"1000_0.05_bib", "1 TeV, 15 mm w/ BIB"},
    {"1000_0.1_bib", "1 TeV, 30 mm w/ BIB"},
    {"1000_1_bib", "1 TeV, 300 mm w/ BIB"},
    {"1000_10_bib", "1 TeV, 3000 mm w/ BIB"},
    {"2500_0.1_bib", "2.5 TeV, 30 mm w/ BIB"},
    {"2500_1_bib", "2.5 TeV, 300 mm w/ BIB"},
    {"2500_10_bib", "2.5 TeV, 3000 mm w/ BIB"},
    {"4000_0.1_bib", "4 TeV, 30 mm w/ BIB"},
    {"4000_1_bib", "4 TeV, 300 mm w/ BIB"},
    {"4000_10_bib", "4 TeV, 3000 mm w/ BIB"},
    {"4500_0.1_bib", "4.5 TeV, 30 mm w/ BIB"},
    {"4500_1_bib", "4.5 TeV, 300 mm w/ BIB"},
    {"4500_10_bib", "4.5 TeV, 3000 mm w/ BIB"},
    {"2500_10_loose", "2.5 TeV, 3000 mm, Loose Timing"},
    {"4000_10_loose", "4 TeV, 3000 mm, Loose Timing"},
    {"2500_10_medium", "2.5 TeV, 3000 mm, Medium Timing"},
    {"4000_10_medium", "4 TeV, 3000 mm, Medium Timing"},
    {"2500_10_tight", "2.5 TeV, 3000 mm, Nominal Timing"},
    {"4000_10_tight", "4 TeV, 3000 mm, Nominal Timing"},
    {"4000_10_loose_velores", "4 TeV, 3000 mm, Loose Timing"},
    {"4000_10_medium_velores", "4 TeV, 3000 mm, Medium Timing"},
    {"2500_10_medium_velores", "2.5 TeV, 3000 mm, Medium Timing"},
    {"4000_10_tight_velores", "4 TeV, 3000 mm, Nominal Timing"}
  };

// Populate the maps with histograms
void initialize_histograms(std::vector<TString> fileNames) {
  for (TString fileName : fileNames){


    stau_mcp_pt_map[fileName] = new TH1F("stau_mcp_pt", ";Monte Carlo Stau p_{T} [GeV]; Tracking particles / 100.0 GeV", 11, 1900.0, 3000.0);
    stau_mcp_eta_map[fileName] = new TH1F("stau_mcp_eta", ";Monte Carlo Stau Eta; Monte Carlo Stau / 0.1", 10, -1.25, 1.25);
    stau_mcp_phi_map[fileName] = new TH1F("stau_mcp_phi", ";Monte Carlo Stau Eta; Monte Carlo Stau / 0.1 radians", 16, -3.2, -3.2);
    stau_mcp_d0_map[fileName] = new TH1F("stau_mcp_d0", ";Monte Carlo Stau d_{0}; Monte Carlo Stau / 10 mm", 50, -250.0, 250.0);
    stau_mcp_z0_map[fileName] = new TH1F("stau_mcp_z0", ";Monte Carlo Stau z_{0}; Monte Carlo Stau / 10 mm", 50, -250.0, 250.0);
    stau_chi2_reduced_map[fileName] = new TH1F("stau_chi2_reduced", "; #chi^{2} / ndf; Tracks / 0.5", 10, 0, 5);

    stau_tp_pt_map[fileName] = new TH1F("stau_tp_pt", ";Tracking Particle Stau p_{T} [GeV]; Tracking particles / 250.0 GeV", 11, 1900.0, 3000.0);
    stau_tp_eta_map[fileName] = new TH1F("stau_tp_eta", ";Tracking Particle Stau Eta; Tracking particle / 0.2", 10, -1.25, 1.25);
    stau_tp_phi_map[fileName] = new TH1F("stau_tp_phi", ";Tracking Particle Stau #phi; Tracking particle / 0.2 radians", 16, -3.2, -3.2);
    stau_tp_d0_map[fileName] = new TH1F("stau_tp_d0", ";Tracking Particle Stau d_{0}; Tracking Particle Stau / 10 mm", 50, -250.0, 250.0);
    stau_tp_z0_map[fileName] = new TH1F("stau_tp_z0", ";Tracking Particle Stau z_{0}; Tracking Particle Stau / 10 mm", 50, -250.0, 250.0);
    stau_matched_pt_map[fileName] = new TH1F("stau_matched_pt", ";Tracking particle p_{T} [GeV]; Tracking particles / 250.0 GeV", 11, 1900.0, 3000.0);
    stau_velores_map[fileName] = new TH1F("stau_velores", "; NHits; |#Delta Stau Velo.| / True Stau Velo. ", 15, 6, 21);
    stau_matched_pt_all_map[fileName] = new TH1F("stau_matched_pt_all", ";Tracking particle p_{T} [GeV]; Tracking particles (all) / 100.0 GeV", 22, 1900.0, 3000.0);

    stau_nhits_map[fileName] = new TH1F("stau_nhits_map_" + fileName, ";Number of Tracker Hits; Tracks / Hit", 15, 6, 21);
    stau_track_pt_map[fileName] = new TH1F("stau_track_pt", ";Stau Track p_{T} [GeV]; Tracks / 100.0 GeV", 50, 0, 5000.0);

    stau_matched_phi_map[fileName] = new TH1F("stau_matched_phi", ";Stau #phi; Stau / 0.2 radians", 16, -3.2, -3.2);
    stau_matched_eta_map[fileName] = new TH1F("stau_matched_eta", ";Stau #eta; Tracking particle / 0.2", 10, -1.25, 1.25);
    stau_track_eta_map[fileName] = new TH1F("stau_track_eta", ";Stau Track #eta; Tracks / 0.1", 50, -2.5, 2.5);
    stau_rxy_decay_map[fileName] = new TH1F("stau_rxy_decay_map_" + fileName, ";Monte Carlo Stau Decay Point r_{xy}; Monte Carlo Stau / 100 mm", 20, 0, 4000.0);
    stau_velo_residuals_map[fileName] = new TH1F("stau_velo_residuals_map_" + fileName, "; Stau True Velocity - Measured Velocity [mm/ns]; Tracks / 0.5 mm/ns", 20, -5, 5);


    displaced_track_pt_map[fileName] = new TH1F("displaced_track_pt_" + fileName, ";Stau Decay Product Track p_{T} [GeV]; Tracks / 500.0 GeV", 10, 0, 5000.0);
    displaced_track_d0_map[fileName] = new TH1F("displaced_track_d0_" + fileName, ";Stau Decay Product Track d_{0} [mm]; Tracks / 10 mm", 30, -150.0, 150.0);
    displaced_track_z0_map[fileName] = new TH1F("displaced_track_z0_" + fileName, ";Stau Decay Product Track z_{0} [mm]; Tracks / 10 mm", 30, -150.0, 150.0);
    displaced_track_eta_map[fileName] = new TH1F("displaced_track_eta_" + fileName, ";Stau Decay Product Track #eta; Tracks / 0.2", 30, -1.5, 1.5);
    displaced_resz0_map[fileName] = new TH1F("displaced_resz0_" + fileName, ";|Track z_{0} - Matched z_{0}|; Tracks / 0.2 mm", 25, 0, 5);
    displaced_resd0_map[fileName] = new TH1F("displaced_resd0_" + fileName, ";|Track d_{0} - Matched d_{0}|; Tracks / 0.2 mm", 20, 0, 2);
    displaced_chi2_reduced_map[fileName] = new TH1F("displaced_chi2_reduced_" + fileName, "; #chi^{2} / ndf; Tracks / 0.5", 10, 0, 5);
    displaced_nhits_map[fileName] = new TH1F("displaced_nhits_" + fileName, ";Number of Tracker Hits; Tracks / Hit", 15, 0, 15);
    displaced_matched_d0_map[fileName] = new TH1F("displaced_matched_d0_" + fileName, "; Tracking Particle d_{0} [mm]; Tracking particles / 6 mm", 10, 0, 150.0);
    displaced_matched_rxy_map[fileName] = new TH1F("displaced_matched_rxy_" + fileName, "; Tracking Particle Vertex r_{xy} [mm]; Tracking particles / 20 mm", 15, 0, 600.0);
    displaced_matched_eta_map[fileName] = new TH1F("displaced_matched_eta_" + fileName, "; Tracking Particle #eta; Tracking particles / 0.2", 24, -2.4, 2.4);
    displaced_tp_d0_map[fileName] = new TH1F("displaced_tp_d0_" + fileName, "; |Displaced Tracking Tracking Particle d_{0}| [mm]; Tracking Particle Stau / 6 mm", 10, 0, 150.0);
    displaced_tp_rxy_map[fileName] = new TH1F("displaced_tp_rxy_" + fileName, "; |Displaced Tracking Tracking Particle d_{0}| [mm]; Tracking Particle Stau / 20 mm", 15, 0, 600.0);
    displaced_tp_eta_map[fileName] = new TH1F("displaced_tp_eta_" + fileName, ";  Tracking Particle #eta; Tracking particles / 0.2", 24, -2.4, 2.4);
    displaced_mcp_d0_map[fileName] = new TH1F("displaced_mcp_d0_" + fileName, ";Monte Carlo Stau Decay Product d_{0}; Monte Carlo Stau / 20 mm", 10, 0, 150.0);
    displaced_mcp_rxy_map[fileName] = new TH1F("displaced_mcp_rxy_" + fileName, ";Monte Carlo Stau Decay Product d_{0}; Monte Carlo Stau / 20 mm", 15, 0, 600.0);
    displaced_mcp_eta_map[fileName] = new TH1F("displaced_mcp_eta_" + fileName, "; Tracking Particle #eta; Tracking particles / 0.2", 24, -2.4, 2.4);
    displaced_mcp_rxy_map_clone[fileName] = new TH1F("displaced_mcp_rxy_" + fileName + "_clone", ";Monte Carlo Stau Decay Product d_{0}; Monte Carlo Stau / 20 mm", 15, 0, 600.0);
    displaced_tp_pt_map[fileName] = new TH1F("displaced_tp_pt_" + fileName, ";Displaced Tracking Particle p_{T} [GeV]; Displaced tracking particles / 500.0 GeV", 10, 0, 5000.0);
    displaced_matched_pt_map[fileName] = new TH1F("displaced_matched_pt_" + fileName, ";Tracking Particle p_{T} [GeV]; Tracking particles / 500.0 GeV", 10, 0, 5000.0);
    displaced_mcp_pt_map[fileName] = new TH1F("displaced_mcp_pt_" + fileName, ";Monte Carlo Stau Decay Product p_{T} [GeV]; Tracking particles / 500.0 GeV", 10, 0, 5000.0);
    fake_track_pt_map[fileName] = new TH1F("fake_track_pt_" + fileName, ";Fake Track p_{T} [GeV]; Tracks / 10.0 GeV", 50, 0, 500.0);
    fake_track_eta_map[fileName] = new TH1F("fake_track_eta_" + fileName, ";Fake Track #eta; Tracks / 0.2", 30, -1.5, 1.5);
    fake_track_d0_map[fileName] = new TH1F("fake_track_d0_" + fileName, ";Fake Track d_{0}; Tracking particles / 25 mm", 20, -250.0, 250.0);
    fake_track_z0_map[fileName] = new TH1F("fake_track_z0_" + fileName, ";Fake Track z_{0}; Tracking particles / 25 mm", 20, -250.0, 250.0);
    fake_track_chi2_reduced_map[fileName] = new TH1F("displaced_overlay_track_chi2_reduced_" + fileName, "; #chi^{2} / ndf; Tracks / 0.5", 20, 0, 10);
    fake_track_nhits_map[fileName] = new TH1F("fake_track_nhits_map_" + fileName, ";Number of Tracker Hits; Tracks / Hit", 15, 0, 15);
    displaced_overlay_track_nhits_map[fileName] = new TH1F("displaced_overlay_track_nhits_" + fileName, ";Number of Tracker Hits; Tracks / Hit", 15, 0, 15);
    displaced_overlay_track_pt_map[fileName] = new TH1F("displaced_overlay_track_pt_" + fileName, ";Fake Track p_{T} [GeV]; Tracks / 10.0 GeV", 50, 0, 500.0);
    displaced_overlay_track_eta_map[fileName] = new TH1F("displaced_overlay_track_eta_" + fileName, ";Fake Track #eta; Tracks / 0.2", 30, -1.5, 1.5);
    displaced_overlay_track_d0_map[fileName] = new TH1F("displaced_overlay_track_d0_" + fileName, ";Fake Track d_{0}; Tracking particles / 25 mm", 20, -250.0, 250.0);
    displaced_overlay_track_z0_map[fileName] = new TH1F("displaced_overlay_track_z0_" + fileName, ";Fake Track z_{0}; Tracking particles / 25 mm", 20, -250.0, 250.0);
    displaced_overlay_track_chi2_reduced_map[fileName] = new TH1F("fake_track_chi2_reduced_" + fileName, "; #chi^{2} / ndf; Tracks / 0.5", 20, 0, 10);
  }
}

void SetPlotStyle();
void mySmallText(Double_t x, Double_t y, Color_t color, char* text);
bool approximatelyEqual(float a, float b, float epsilon); // helper functions to identify duplicate stau tracks
bool isApproximatelyEqualToAny(const std::vector<float>& vec, float value, float epsilon);
void rootAnalyzer(const TString fileName, std::string timing);

// Function to compute the mean of a vector
float computeMean(const std::vector<float>& vec) {
    if (vec.empty()) return 0.0f;
    float sum = 0.0f;
    for (float value : vec) {
        sum += value;
    }
    return sum / vec.size();
}

// Function to compute the standard deviation of a vector
float computeStandardDeviation(const std::vector<float>& vec, float mean) {
    if (vec.size() < 2) return 0.0f;
    float sum = 0.0f;
    for (float value : vec) {
        float diff = value - mean;
        sum += diff * diff;
    }
    return std::sqrt(sum / (vec.size() - 1));  // Using (n-1) for sample standard deviation
}

// Function to compute the standard error of the mean
float computeStandardError(const std::vector<float>& vec, float stdDev) {
    if (vec.size() == 0) return 0.0f;
    return stdDev / std::sqrt(vec.size());
}

// Example Usage: callRootAnalyzer() (change filenames in this file)
// Will output plots as pdfs to a subdirectory named according to the filename

void rootAnalyzer(
const TString fileName,
std::string timing
){  

  double stau_mass = (fileName[0] == '1') ? 1000.0 :
                       (fileName[0] == '2') ? 2500.0 :
                       (fileName[0] == '4') ? 4000.0 : 0.0;
    
    // Example Usage: rootAnalyzer("1000_0.05_reco") 
    // Will output plots as pdfs to a subdirectory named according to the filename
    // ----------------------------------------------------------------------------------------------------------------
    // Set plot style and other configurations
    SetPlotStyle();



    // ----------------------------------------------------------------------------------------------------------------
    // Read NTuples from one file
    TChain* MCPs = new TChain("MCPs");
    std::string rootFileDir = "recoRootFiles"; 
    if (timing == "tight"){
      rootFileDir = "tightTimingRecoRootFiles/"; 
    }
    if (timing == "medium"){
      rootFileDir = "mediumTimingRecoRootFiles/";
    }
    if (timing == "loose"){
      rootFileDir = "looseTimingRecoRootFiles/";
    }
    if (timing == "mediumFix"){
      rootFileDir = "timingFixMediumRootFiles/";
    }
    if (timing == "ITMedium"){
      rootFileDir = "timingMediumITSeedingRootFiles/";
    }
    if (timing == "ITMediumNoBIB"){
      rootFileDir = "mediumTimingNoBIBITSeedingRecoRootFiles/";
    }
    if (timing == "ITOTMediumNoBIB"){
      rootFileDir = "mediumTimingNoBIBITOTSeedingRecoRootFiles/";
    }
    if (timing == "mediumNoBIB"){
      rootFileDir = "mediumTimingNoBIBRecoRootFilesStripsOn/";
    }
    if (timing == "mediumNoBIBPropBack"){
      rootFileDir = "mediumTimingNoBIBPropBackStripsOn/";
    }
    if (timing == "mediumNoBIBStripsOff"){
      rootFileDir = "mediumTimingNoBIBStripsOff/";
    }
    if (timing == "mediumPrompt"){
      rootFileDir = "MediumTimingPrompt10pBIB/";
    }
    if (timing == "tightPrompt"){
      rootFileDir = "TightTimingPrompt10pBIB/";
    }
    if (timing == "loosePrompt"){
      rootFileDir = "LooseTimingPrompt10pBIB/";
    }
    if (timing == "mediumDisplaced"){
      rootFileDir = "MediumTimingTwoPasses0pBIB/";
    }
    if (timing == "openHouse"){
      rootFileDir = "openHouseRootFiles/";
    }
    if (timing == "medium10pBIBLRTTest"){
      rootFileDir = "DisplacedMedium10pBIBNoCut/";
    }

    MCPs->Add(rootFileDir + fileName + "_reco.root");
    TChain* StauMCPs = new TChain("StauMCPs");
    StauMCPs->Add(rootFileDir + fileName + "_reco.root");
    TChain* StauDecayProductsMCPs = new TChain("StauDecayProductsMCPs");
    StauDecayProductsMCPs->Add(rootFileDir + fileName + "_reco.root");
    TChain* AllTracks = new TChain("AllTracks");
    AllTracks->Add(rootFileDir + fileName + "_reco.root");
    TChain* AllStauTracks = new TChain("AllStauTracks");
    AllStauTracks->Add(rootFileDir + fileName + "_reco.root");
    TChain* AllDecayProductTracks = new TChain("AllDecayProductTracks");
    AllDecayProductTracks->Add(rootFileDir + fileName + "_reco.root");
    TChain* AllFakeTracks = new TChain("AllFakeTracks");
    AllFakeTracks->Add(rootFileDir + fileName + "_reco.root");
    TChain* AllHits = new TChain("AllHits");
    AllHits->Add(rootFileDir + fileName + "_reco.root");
    TChain* AllHitsActual = new TChain("AllHitsActual");
    AllHitsActual->Add(rootFileDir + fileName + "_reco.root");

    if (MCPs->GetEntries() == 0 || StauMCPs->GetEntries() == 0 || StauDecayProductsMCPs->GetEntries() == 0 || AllTracks->GetEntries() == 0 || AllStauTracks->GetEntries() == 0 || AllDecayProductTracks->GetEntries() == 0 || AllFakeTracks->GetEntries() == 0 || AllHits->GetEntries() == 0 || AllHitsActual->GetEntries() == 0) {
        cout << "A file / branch doesn't exist or is empty, returning..." << "\n";
        return;
    }

    // ----------------------------------------------------------------------------------------------------------------
    // Define leaves and branches
    
    std::vector<float>* mcp_pt;
    std::vector<float>* mcp_phi;
    std::vector<float>* mcp_eta;
    std::vector<int>* pdgid;
    std::vector<int>* status;
    std::vector<float>* prod_vertex_r;
    std::vector<float>* prod_vertex_z;
    std::vector<float>* prod_endpoint_r;
    std::vector<float>* prod_endpoint_z;
    std::vector<float>* prod_traveldist;
    std::vector<float>* prod_time;
    std::vector<int>* id;

    std::vector<float>* mcp_stau_pt;
    std::vector<float>* mcp_stau_phi;
    std::vector<float>* mcp_stau_eta;
    std::vector<float>* mcp_stau_d0;
    std::vector<float>* mcp_stau_z0;
    std::vector<float>* mcp_daughter_r_vertex;
    std::vector<float>* mcp_daughter_r_endpoint;
    std::vector<float>* mcp_daughter_z_vertex;
    std::vector<float>* mcp_daughter_z_endpoint;
    std::vector<float>* prod_stau_vertex_r;
    std::vector<float>* prod_stau_vertex_z;
    std::vector<float>* prod_stau_endpoint_r;
    std::vector<float>* prod_stau_endpoint_z;
    std::vector<float>* prod_traveldist_stau;
    std::vector<int>* mcp_stau_track_bool;
    std::vector<int>* mcp_stau_track_reconstructable_bool;

    std::vector<float>* mcp_daughter_pt;
    std::vector<float>* mcp_daughter_phi;
    std::vector<float>* mcp_daughter_eta;
    std::vector<float>* mcp_daughter_d0;
    std::vector<float>* mcp_daughter_z0;
    std::vector<int>* mcp_daughter_track_bool;
    std::vector<int>* mcp_daughter_track_reconstructable_bool;

    std::vector<int>* nhits;
    std::vector<int>* pixel_nhits;
    std::vector<int>* inner_nhits;
    std::vector<int>* outer_nhits;
    std::vector<float>* track_pt;
    std::vector<float>* track_eta;
    std::vector<float>* track_theta;
    std::vector<int>* ndf;
    std::vector<float>* chi2;
    std::vector<float>* chi2_red;

    std::vector<float>* LC_stau_pt_match;
    std::vector<float>* LC_stau_reco_velo;
    std::vector<float>* LC_stau_true_velo;
    std::vector<float>* LC_stau_velo_res;
    std::vector<float>* LC_stau_track_pt;
    std::vector<float>* LC_stau_track_eta;
    std::vector<float>* LC_stau_eta_match;
    std::vector<float>* LC_stau_track_theta;
    std::vector<float>* LC_stau_phi_match;
    std::vector<int>* LC_stau_ndf;
    std::vector<float>* LC_stau_chi2;
    std::vector<float>* LC_stau_d0;
    std::vector<float>* LC_stau_z0;
    std::vector<int>* LC_stau_nhits;
    std::vector<int>* LC_stau_pixel_nhits;
    std::vector<int>* LC_stau_inner_nhits;
    std::vector<int>* LC_stau_outer_nhits;
    std::vector<float>* LC_stau_pt_res;
    std::vector<float>* LC_stau_dr;
    std::vector<float>* LC_stau_hit_r;
    std::vector<float>* LC_stau_hit_x;
    std::vector<float>* LC_stau_hit_y;
    std::vector<float>* LC_stau_hit_z;

    std::vector<float>* LC_daughter_pt_match;
    std::vector<float>* LC_daughter_track_pt;
    std::vector<float>* LC_daughter_track_eta;
    std::vector<float>* LC_daughter_eta_match;
    std::vector<float>* LC_daughter_track_theta;
    std::vector<float>* LC_daughter_phi_match;
    std::vector<int>* LC_daughter_ndf;
    std::vector<float>* LC_daughter_chi2;
    std::vector<float>* LC_daughter_d0;
    std::vector<float>* LC_daughter_z0;
    std::vector<float>* LC_daughter_d0_match;
    std::vector<float>* LC_daughter_z0_match;
    std::vector<float>* LC_daughter_r_vertex_match;
    std::vector<float>* LC_daughter_r_endpoint_match;
    std::vector<float>* LC_daughter_z_vertex_match;
    std::vector<float>* LC_daughter_z_endpoint_match;
    std::vector<int>* LC_daughter_nhits;
    std::vector<int>* LC_daughter_pixel_nhits;
    std::vector<int>* LC_daughter_inner_nhits;
    std::vector<int>* LC_daughter_outer_nhits;
    std::vector<float>* LC_daughter_pt_res;
    std::vector<float>* LC_daughter_dr;

    std::vector<float>* fake_theta;
    std::vector<float>* fake_eta;
    std::vector<float>* fake_pt;
    std::vector<float>* fake_phi;
    std::vector<float>* fake_d0;
    std::vector<float>* fake_z0;
    std::vector<int>* fake_ndf;
    std::vector<float>* fake_chi2;
    std::vector<float>* fake_chi2_reduced;
    std::vector<int>* fake_nhits;
    std::vector<int>* fake_pixel_nhits;
    std::vector<int>* fake_inner_nhits;
    std::vector<int>* fake_outer_nhits;

    std::vector<float>* x;
    std::vector<float>* y;
    std::vector<float>* z;
    std::vector<int>* hit_pdgid;
    std::vector<float>* time;
    std::vector<float>* corrected_time;
    std::vector<int>* hit_layer;
    std::vector<int>* hit_detector;
    std::vector<int>* hit_side;


    std::vector<float>* hitR;
    std::vector<float>* hitZ;

    TBranch* b_mcp_pt;
    TBranch* b_mcp_phi;
    TBranch* b_mcp_eta;
    TBranch* b_pdgid;
    TBranch* b_status;
    TBranch* b_prod_vertex_r;
    TBranch* b_prod_vertex_z;
    TBranch* b_prod_endpoint_r;
    TBranch* b_prod_endpoint_z;
    TBranch* b_prod_traveldist;
    TBranch* b_prod_time;
    TBranch* b_id;

    TBranch* b_mcp_stau_pt;
    TBranch* b_mcp_stau_phi;
    TBranch* b_mcp_stau_eta;
    TBranch* b_mcp_stau_d0;
    TBranch* b_mcp_stau_z0;
    TBranch* b_prod_stau_vertex_r;
    TBranch* b_prod_stau_vertex_z;
    TBranch* b_prod_stau_endpoint_r;
    TBranch* b_prod_stau_endpoint_z;
    TBranch* b_prod_traveldist_stau;
    TBranch* b_mcp_stau_track_bool;
    TBranch* b_mcp_stau_track_reconstructable_bool;

    TBranch* b_mcp_daughter_pt;
    TBranch* b_mcp_daughter_phi;
    TBranch* b_mcp_daughter_eta;
    TBranch* b_mcp_daughter_d0;
    TBranch* b_mcp_daughter_z0;
    TBranch* b_mcp_daughter_r_vertex;
    TBranch* b_mcp_daughter_r_endpoint;
    TBranch* b_mcp_daughter_z_vertex;
    TBranch* b_mcp_daughter_z_endpoint;
    TBranch* b_mcp_daughter_track_bool;
    TBranch* b_mcp_daughter_track_reconstructable_bool;

    TBranch* b_nhits;
    TBranch* b_pixel_nhits;
    TBranch* b_inner_nhits;
    TBranch* b_outer_nhits;
    TBranch* b_track_pt;
    TBranch* b_track_eta;
    TBranch* b_track_theta;
    TBranch* b_ndf;
    TBranch* b_chi2;
    TBranch* b_chi2_red;

    TBranch* b_LC_stau_pt_match;
    TBranch* b_LC_stau_reco_velo;
    TBranch* b_LC_stau_true_velo;
    TBranch* b_LC_stau_velo_res;
    TBranch* b_LC_stau_track_pt;
    TBranch* b_LC_stau_track_eta;
    TBranch* b_LC_stau_eta_match;
    TBranch* b_LC_stau_track_theta;
    TBranch* b_LC_stau_phi_match;
    TBranch* b_LC_stau_ndf;
    TBranch* b_LC_stau_chi2;
    TBranch* b_LC_stau_d0;
    TBranch* b_LC_stau_z0;
    TBranch* b_LC_stau_nhits;
    TBranch* b_LC_stau_pixel_nhits;
    TBranch* b_LC_stau_inner_nhits;
    TBranch* b_LC_stau_outer_nhits;
    TBranch* b_LC_stau_pt_res;
    TBranch* b_LC_stau_dr;
    TBranch* b_LC_stau_hit_r;
    TBranch* b_LC_stau_hit_x;
    TBranch* b_LC_stau_hit_y;
    TBranch* b_LC_stau_hit_z;

    TBranch* b_LC_daughter_pt_match;
    TBranch* b_LC_daughter_track_pt;
    TBranch* b_LC_daughter_track_eta;
    TBranch* b_LC_daughter_eta_match;
    TBranch* b_LC_daughter_track_theta;
    TBranch* b_LC_daughter_phi_match;
    TBranch* b_LC_daughter_ndf;
    TBranch* b_LC_daughter_chi2;
    TBranch* b_LC_daughter_d0;
    TBranch* b_LC_daughter_z0;
    TBranch* b_LC_daughter_d0_match;
    TBranch* b_LC_daughter_z0_match;
    TBranch* b_LC_daughter_r_vertex_match;
    TBranch* b_LC_daughter_r_endpoint_match;
    TBranch* b_LC_daughter_z_vertex_match;
    TBranch* b_LC_daughter_z_endpoint_match;
    TBranch* b_LC_daughter_nhits;
    TBranch* b_LC_daughter_pixel_nhits;
    TBranch* b_LC_daughter_inner_nhits;
    TBranch* b_LC_daughter_outer_nhits;
    TBranch* b_LC_daughter_pt_res;
    TBranch* b_LC_daughter_dr;

    TBranch* b_fake_theta;
    TBranch* b_fake_eta;
    TBranch* b_fake_pt;
    TBranch* b_fake_phi;
    TBranch* b_fake_d0;
    TBranch* b_fake_z0;
    TBranch* b_fake_ndf;
    TBranch* b_fake_chi2;
    TBranch* b_fake_chi2_reduced;
    TBranch* b_fake_nhits;
    TBranch* b_fake_pixel_nhits;
    TBranch* b_fake_inner_nhits;
    TBranch* b_fake_outer_nhits;

    TBranch* b_x;
    TBranch* b_y;
    TBranch* b_z;
    TBranch* b_hit_pdgid;
    TBranch* b_time;
    TBranch* b_corrected_time;
    TBranch* b_hit_layer;
    TBranch* b_hit_detector;
    TBranch* b_hit_side;

    TBranch* b_hitR;
    TBranch* b_hitZ;

    mcp_pt = 0;
    mcp_phi = 0;
    mcp_eta = 0;
    pdgid = 0;
    status = 0;
    prod_vertex_r = 0;
    prod_vertex_z = 0;
    prod_endpoint_r = 0;
    prod_endpoint_z = 0;
    prod_traveldist = 0;
    prod_time = 0;
    id = 0;

    mcp_stau_pt = 0;
    mcp_stau_phi = 0;
    mcp_stau_eta = 0;
    mcp_stau_d0 = 0;
    mcp_stau_z0 = 0;
    prod_stau_vertex_r = 0;
    prod_stau_vertex_z = 0;
    prod_stau_endpoint_r = 0;
    prod_stau_endpoint_z = 0;
    prod_traveldist_stau = 0;
    mcp_stau_track_bool = 0;
    mcp_stau_track_reconstructable_bool = 0;

    mcp_daughter_pt = 0;
    mcp_daughter_phi = 0;
    mcp_daughter_eta = 0;
    mcp_daughter_d0 = 0;
    mcp_daughter_z0 = 0;
    mcp_daughter_r_vertex = 0;
    mcp_daughter_r_endpoint = 0;
    mcp_daughter_z_vertex = 0;
    mcp_daughter_z_endpoint = 0;

    mcp_daughter_track_bool = 0;
    mcp_daughter_track_reconstructable_bool = 0;

    nhits = 0;
    pixel_nhits = 0;
    inner_nhits = 0;
    outer_nhits = 0;
    track_pt = 0;
    track_eta = 0;
    track_theta = 0;
    ndf = 0;
    chi2 = 0;
    chi2_red = 0;

    LC_stau_pt_match = 0;
    LC_stau_reco_velo = 0;
    LC_stau_true_velo = 0;
    LC_stau_velo_res = 0;
    LC_stau_track_pt = 0;
    LC_stau_track_eta = 0;
    LC_stau_eta_match = 0;
    LC_stau_track_theta = 0;
    LC_stau_phi_match = 0;
    LC_stau_ndf = 0;
    LC_stau_chi2 = 0;
    LC_stau_d0 = 0;
    LC_stau_z0 = 0;
    LC_stau_nhits = 0;
    LC_stau_pixel_nhits = 0;
    LC_stau_inner_nhits = 0;
    LC_stau_outer_nhits = 0;
    LC_stau_pt_res = 0;
    LC_stau_dr = 0;
    LC_stau_hit_r = 0;
    LC_stau_hit_x = 0;
    LC_stau_hit_y = 0;
    LC_stau_hit_z = 0;

    LC_daughter_pt_match = 0;
    LC_daughter_track_pt = 0;
    LC_daughter_track_eta = 0;
    LC_daughter_eta_match = 0;
    LC_daughter_track_theta = 0;
    LC_daughter_phi_match = 0;
    LC_daughter_ndf = 0;
    LC_daughter_chi2 = 0;
    LC_daughter_d0 = 0;
    LC_daughter_z0 = 0;
    LC_daughter_d0_match = 0;
    LC_daughter_z0_match = 0;
    LC_daughter_r_vertex_match = 0;
    LC_daughter_r_endpoint_match = 0;
    LC_daughter_z_vertex_match = 0;
    LC_daughter_z_endpoint_match = 0;
    LC_daughter_nhits = 0;
    LC_daughter_pixel_nhits = 0;
    LC_daughter_inner_nhits = 0;
    LC_daughter_outer_nhits = 0;
    LC_daughter_pt_res = 0;
    LC_daughter_dr = 0;

    fake_theta = 0;
    fake_eta = 0;
    fake_pt = 0;
    fake_phi = 0;
    fake_d0 = 0;
    fake_z0 = 0;
    fake_ndf = 0;
    fake_chi2 = 0;
    fake_chi2_reduced = 0;
    fake_nhits = 0;
    fake_pixel_nhits = 0;
    fake_inner_nhits = 0;
    fake_outer_nhits = 0;

    x = 0;
    y = 0;
    z = 0;
    hit_pdgid = 0;
    time = 0;
    corrected_time = 0;
    hit_layer = 0;
    hit_detector = 0;
    hit_side = 0;

    hitR = 0;
    hitZ = 0;
    
    // ----------------------------------------------------------------------------------------------------------------
    // Load the files and set branch addresses
    MCPs->SetBranchAddress("mcp_pt", &mcp_pt, &b_mcp_pt);
    MCPs->SetBranchAddress("mcp_phi", &mcp_phi, &b_mcp_phi);
    MCPs->SetBranchAddress("mcp_eta", &mcp_eta, &b_mcp_eta);
    MCPs->SetBranchAddress("pdgid", &pdgid, &b_pdgid);
    MCPs->SetBranchAddress("status", &status, &b_status);
    MCPs->SetBranchAddress("prod_vertex_r", &prod_vertex_r, &b_prod_vertex_r);
    MCPs->SetBranchAddress("prod_vertex_z", &prod_vertex_z, &b_prod_vertex_z);
    MCPs->SetBranchAddress("prod_endpoint_r", &prod_endpoint_r, &b_prod_endpoint_r);
    MCPs->SetBranchAddress("prod_endpoint_z", &prod_endpoint_z, &b_prod_endpoint_z);
    MCPs->SetBranchAddress("prod_traveldist", &prod_traveldist, &b_prod_traveldist);
    MCPs->SetBranchAddress("prod_time", &prod_time, &b_prod_time);
    MCPs->SetBranchAddress("id", &id, &b_id);

    StauMCPs->SetBranchAddress("mcp_stau_pt", &mcp_stau_pt, &b_mcp_stau_pt);
    StauMCPs->SetBranchAddress("mcp_stau_phi", &mcp_stau_phi, &b_mcp_stau_phi);
    StauMCPs->SetBranchAddress("mcp_stau_eta", &mcp_stau_eta, &b_mcp_stau_eta);
    StauMCPs->SetBranchAddress("mcp_stau_d0", &mcp_stau_d0, &b_mcp_stau_d0);
    StauMCPs->SetBranchAddress("mcp_stau_z0", &mcp_stau_z0, &b_mcp_stau_z0);
    StauMCPs->SetBranchAddress("prod_stau_vertex_r", &prod_stau_vertex_r, &b_prod_stau_vertex_r);
    StauMCPs->SetBranchAddress("prod_stau_vertex_z", &prod_stau_vertex_z, &b_prod_stau_vertex_z);
    StauMCPs->SetBranchAddress("prod_stau_endpoint_r", &prod_stau_endpoint_r, &b_prod_stau_endpoint_r);
    StauMCPs->SetBranchAddress("prod_stau_endpoint_z", &prod_stau_endpoint_z, &b_prod_stau_endpoint_z);
    StauMCPs->SetBranchAddress("prod_traveldist_stau", &prod_traveldist_stau, &b_prod_traveldist_stau);
    StauMCPs->SetBranchAddress("mcp_stau_track_bool", &mcp_stau_track_bool, &b_mcp_stau_track_bool);
    StauMCPs->SetBranchAddress("mcp_stau_track_reconstructable_bool", &mcp_stau_track_reconstructable_bool, &b_mcp_stau_track_reconstructable_bool);

    StauDecayProductsMCPs->SetBranchAddress("mcp_daughter_pt", &mcp_daughter_pt, &b_mcp_daughter_pt);
    StauDecayProductsMCPs->SetBranchAddress("mcp_daughter_phi", &mcp_daughter_phi, &b_mcp_daughter_phi);
    StauDecayProductsMCPs->SetBranchAddress("mcp_daughter_eta", &mcp_daughter_eta, &b_mcp_daughter_eta);
    StauDecayProductsMCPs->SetBranchAddress("mcp_daughter_d0", &mcp_daughter_d0, &b_mcp_daughter_d0);
    StauDecayProductsMCPs->SetBranchAddress("mcp_daughter_z0", &mcp_daughter_z0, &b_mcp_daughter_z0);
    StauDecayProductsMCPs->SetBranchAddress("mcp_daughter_r_vertex", &mcp_daughter_r_vertex, &b_mcp_daughter_r_vertex);
    StauDecayProductsMCPs->SetBranchAddress("mcp_daughter_r_endpoint", &mcp_daughter_r_endpoint, &b_mcp_daughter_r_endpoint);
    StauDecayProductsMCPs->SetBranchAddress("mcp_daughter_z_vertex", &mcp_daughter_z_vertex, &b_mcp_daughter_z_vertex);
    StauDecayProductsMCPs->SetBranchAddress("mcp_daughter_z_endpoint", &mcp_daughter_z_endpoint, &b_mcp_daughter_z_endpoint);
    StauDecayProductsMCPs->SetBranchAddress("mcp_daughter_track_bool", &mcp_daughter_track_bool, &b_mcp_daughter_track_bool);
    StauDecayProductsMCPs->SetBranchAddress("mcp_daughter_track_reconstructable_bool", &mcp_daughter_track_reconstructable_bool, &b_mcp_daughter_track_reconstructable_bool);

    AllTracks->SetBranchAddress("nhits", &nhits, &b_nhits);
    AllTracks->SetBranchAddress("pixel_nhits", &pixel_nhits, &b_pixel_nhits);
    AllTracks->SetBranchAddress("inner_nhits", &inner_nhits, &b_inner_nhits);
    AllTracks->SetBranchAddress("outer_nhits", &outer_nhits, &b_outer_nhits);
    AllTracks->SetBranchAddress("track_pt", &track_pt, &b_track_pt);
    AllTracks->SetBranchAddress("track_eta", &track_eta, &b_track_eta);
    AllTracks->SetBranchAddress("track_theta", &track_theta, &b_track_theta);
    AllTracks->SetBranchAddress("ndf", &ndf, &b_ndf);
    AllTracks->SetBranchAddress("chi2", &chi2, &b_chi2);
    AllTracks->SetBranchAddress("chi2_red", &chi2_red, &b_chi2_red);

    AllStauTracks->SetBranchAddress("LC_stau_pt_match", &LC_stau_pt_match, &b_LC_stau_pt_match);
    AllStauTracks->SetBranchAddress("LC_stau_reco_velo", &LC_stau_reco_velo, &b_LC_stau_reco_velo);
    AllStauTracks->SetBranchAddress("LC_stau_true_velo", &LC_stau_true_velo, &b_LC_stau_true_velo);
    AllStauTracks->SetBranchAddress("LC_stau_velo_res", &LC_stau_velo_res, &b_LC_stau_velo_res);
    AllStauTracks->SetBranchAddress("LC_stau_track_pt", &LC_stau_track_pt, &b_LC_stau_track_pt);
    AllStauTracks->SetBranchAddress("LC_stau_track_eta", &LC_stau_track_eta, &b_LC_stau_track_eta);
    AllStauTracks->SetBranchAddress("LC_stau_eta_match", &LC_stau_eta_match, &b_LC_stau_eta_match);
    AllStauTracks->SetBranchAddress("LC_stau_track_theta", &LC_stau_track_theta, &b_LC_stau_track_theta);
    AllStauTracks->SetBranchAddress("LC_stau_phi_match", &LC_stau_phi_match, &b_LC_stau_phi_match);
    AllStauTracks->SetBranchAddress("LC_stau_ndf", &LC_stau_ndf, &b_LC_stau_ndf);
    AllStauTracks->SetBranchAddress("LC_stau_chi2", &LC_stau_chi2, &b_LC_stau_chi2);
    AllStauTracks->SetBranchAddress("LC_stau_d0", &LC_stau_d0, &b_LC_stau_d0);
    AllStauTracks->SetBranchAddress("LC_stau_z0", &LC_stau_z0, &b_LC_stau_z0);
    AllStauTracks->SetBranchAddress("LC_stau_nhits", &LC_stau_nhits, &b_LC_stau_nhits);
    AllStauTracks->SetBranchAddress("LC_stau_pixel_nhits", &LC_stau_pixel_nhits, &b_LC_stau_pixel_nhits);
    AllStauTracks->SetBranchAddress("LC_stau_inner_nhits", &LC_stau_inner_nhits, &b_LC_stau_inner_nhits);
    AllStauTracks->SetBranchAddress("LC_stau_outer_nhits", &LC_stau_outer_nhits, &b_LC_stau_outer_nhits);
    AllStauTracks->SetBranchAddress("LC_stau_pt_res", &LC_stau_pt_res, &b_LC_stau_pt_res);
    AllStauTracks->SetBranchAddress("LC_stau_dr", &LC_stau_dr, &b_LC_stau_dr);
    AllStauTracks->SetBranchAddress("LC_stau_hit_r", &LC_stau_hit_r, &b_LC_stau_hit_r);
    AllStauTracks->SetBranchAddress("LC_stau_hit_x", &LC_stau_hit_x, &b_LC_stau_hit_x);
    AllStauTracks->SetBranchAddress("LC_stau_hit_y", &LC_stau_hit_y, &b_LC_stau_hit_y);
    AllStauTracks->SetBranchAddress("LC_stau_hit_z", &LC_stau_hit_z, &b_LC_stau_hit_z);

    AllDecayProductTracks->SetBranchAddress("LC_daughter_pt_match", &LC_daughter_pt_match, &b_LC_daughter_pt_match);
    AllDecayProductTracks->SetBranchAddress("LC_daughter_track_pt", &LC_daughter_track_pt, &b_LC_daughter_track_pt);
    AllDecayProductTracks->SetBranchAddress("LC_daughter_track_eta", &LC_daughter_track_eta, &b_LC_daughter_track_eta);
    AllDecayProductTracks->SetBranchAddress("LC_daughter_eta_match", &LC_daughter_eta_match, &b_LC_daughter_eta_match);
    AllDecayProductTracks->SetBranchAddress("LC_daughter_track_theta", &LC_daughter_track_theta, &b_LC_daughter_track_theta);
    AllDecayProductTracks->SetBranchAddress("LC_daughter_phi_match", &LC_daughter_phi_match, &b_LC_daughter_phi_match);
    AllDecayProductTracks->SetBranchAddress("LC_daughter_ndf", &LC_daughter_ndf, &b_LC_daughter_ndf);
    AllDecayProductTracks->SetBranchAddress("LC_daughter_chi2", &LC_daughter_chi2, &b_LC_daughter_chi2);
    AllDecayProductTracks->SetBranchAddress("LC_daughter_d0", &LC_daughter_d0, &b_LC_daughter_d0);
    AllDecayProductTracks->SetBranchAddress("LC_daughter_z0", &LC_daughter_z0, &b_LC_daughter_z0);
    AllDecayProductTracks->SetBranchAddress("LC_daughter_d0_match", &LC_daughter_d0_match, &b_LC_daughter_d0_match);
    AllDecayProductTracks->SetBranchAddress("LC_daughter_z0_match", &LC_daughter_z0_match, &b_LC_daughter_z0_match);
    AllDecayProductTracks->SetBranchAddress("LC_daughter_r_vertex_match", &LC_daughter_r_vertex_match, &b_LC_daughter_r_vertex_match);
    AllDecayProductTracks->SetBranchAddress("LC_daughter_r_endpoint_match", &LC_daughter_r_endpoint_match, &b_LC_daughter_r_endpoint_match);
    AllDecayProductTracks->SetBranchAddress("LC_daughter_z_vertex_match", &LC_daughter_z_vertex_match, &b_LC_daughter_z_vertex_match);
    AllDecayProductTracks->SetBranchAddress("LC_daughter_z_endpoint_match", &LC_daughter_z_endpoint_match, &b_LC_daughter_z_endpoint_match);
    AllDecayProductTracks->SetBranchAddress("LC_daughter_nhits", &LC_daughter_nhits, &b_LC_daughter_nhits);
    AllDecayProductTracks->SetBranchAddress("LC_daughter_pixel_nhits", &LC_daughter_pixel_nhits, &b_LC_daughter_pixel_nhits);
    AllDecayProductTracks->SetBranchAddress("LC_daughter_inner_nhits", &LC_daughter_inner_nhits, &b_LC_daughter_inner_nhits);
    AllDecayProductTracks->SetBranchAddress("LC_daughter_outer_nhits", &LC_daughter_outer_nhits, &b_LC_daughter_outer_nhits);
    AllDecayProductTracks->SetBranchAddress("LC_daughter_pt_res", &LC_daughter_pt_res, &b_LC_daughter_pt_res);
    AllDecayProductTracks->SetBranchAddress("LC_daughter_dr", &LC_daughter_dr, &b_LC_daughter_dr);

    AllFakeTracks->SetBranchAddress("fake_theta", &fake_theta, &b_fake_theta);
    AllFakeTracks->SetBranchAddress("fake_eta", &fake_eta, &b_fake_eta);
    AllFakeTracks->SetBranchAddress("fake_pt", &fake_pt, &b_fake_pt);
    AllFakeTracks->SetBranchAddress("fake_phi", &fake_phi, &b_fake_phi);
    AllFakeTracks->SetBranchAddress("fake_d0", &fake_d0, &b_fake_d0);
    AllFakeTracks->SetBranchAddress("fake_z0", &fake_z0, &b_fake_z0);
    AllFakeTracks->SetBranchAddress("fake_ndf", &fake_ndf, &b_fake_ndf);
    AllFakeTracks->SetBranchAddress("fake_chi2", &fake_chi2, &b_fake_chi2);
    AllFakeTracks->SetBranchAddress("fake_chi2_reduced", &fake_chi2_reduced, &b_fake_chi2_reduced);
    AllFakeTracks->SetBranchAddress("fake_nhits", &fake_nhits, &b_fake_nhits);
    AllFakeTracks->SetBranchAddress("fake_pixel_nhits", &fake_pixel_nhits, &b_fake_pixel_nhits);
    AllFakeTracks->SetBranchAddress("fake_inner_nhits", &fake_inner_nhits, &b_fake_inner_nhits);
    AllFakeTracks->SetBranchAddress("fake_outer_nhits", &fake_outer_nhits, &b_fake_outer_nhits);

    AllHits->SetBranchAddress("x", &x, &b_x);
    AllHits->SetBranchAddress("y", &y, &b_y);
    AllHits->SetBranchAddress("z", &z, &b_z);
    AllHits->SetBranchAddress("hit_pdgid", &hit_pdgid, &b_hit_pdgid);
    AllHits->SetBranchAddress("time", &time, &b_time);
    AllHits->SetBranchAddress("corrected_time", &corrected_time, &b_corrected_time);
    AllHits->SetBranchAddress("hit_layer", &hit_layer, &b_hit_layer);
    AllHits->SetBranchAddress("hit_detector", &hit_detector, &b_hit_detector);
    AllHits->SetBranchAddress("hit_side", &hit_side, &b_hit_side);

    AllHitsActual->SetBranchAddress("r", &hitR, &b_hitR);
    AllHitsActual->SetBranchAddress("z", &hitZ, &b_hitZ);



    // ----------------------------------------------------------------------------------------------------------------
    // Declare histograms to be written to and plotted

    TH2F* hit_rz = new TH2F("hit_rz", ";|z| [mm] ; r [mm]", 11500, 0, 2350, 8000, 0, 1650);

    TH1F* stau_mcp_pt = new TH1F("stau_mcp_pt", ";Monte Carlo Stau p_{T} [GeV]; Tracking particles / 100.0 GeV", 11, 1900.0, 3000.0);
    TH1F* stau_mcp_p = new TH1F("stau_mcp_p", ";Monte Carlo Stau p [GeV]; Tracking particles / 100.0 GeV", 20, 0, 5000.0);
    TH1F* stau_mcp_E = new TH1F("stau_mcp_E", ";Monte Carlo Stau Energy [GeV]; Tracking particles / 100.0 GeV", 24, 4000.0, 5200.0);
    TH2F* stau_mcp_pt_eta = new TH2F("stau_mcp_pt_eta", ";#eta ; p_{T} [GeV]", 15, -1.5, 1.5, 20, 0, 5000.0);
    TH1F* stau_mcp_eta = new TH1F("stau_mcp_eta", ";Monte Carlo Stau Eta; Monte Carlo Stau / 0.1", 15, -1.25, 1.25);
    TH1F* stau_mcp_phi = new TH1F("stau_mcp_phi", ";Monte Carlo Stau Eta; Monte Carlo Stau / 0.1 radians", 16, -3.2, -3.2);
    TH1F* stau_mcp_d0 = new TH1F("stau_mcp_d0", ";Monte Carlo Stau d_{0}; Monte Carlo Stau / 10 mm", 50, -250.0, 250.0);
    TH1F* stau_mcp_z0 = new TH1F("stau_mcp_z0", ";Monte Carlo Stau z_{0}; Monte Carlo Stau / 10 mm", 50, -250.0, 250.0);

    TH1F* stau_tp_pt = new TH1F("stau_tp_pt", ";Tracking Particle Stau p_{T} [GeV]; Tracking particles / 250.0 GeV", 11, 1900.0, 3000.0);
    TH1F* stau_tp_eta = new TH1F("stau_tp_eta", ";Tracking Particle Stau Eta; Tracking particle / 0.2", 15, -1.25, 1.25);
    TH1F* stau_tp_phi = new TH1F("stau_tp_phi", ";Tracking Particle Stau #phi; Tracking particle / 0.2 radians", 16, -3.2, -3.2);
    TH1F* stau_tp_d0 = new TH1F("stau_tp_d0", ";Tracking Particle Stau d_{0}; Tracking Particle Stau / 10 mm", 50, -250.0, 250.0);
    TH1F* stau_tp_z0 = new TH1F("stau_tp_z0", ";Tracking Particle Stau z_{0}; Tracking Particle Stau / 10 mm", 50, -250.0, 250.0);

    TH1F* displaced_mcp_pt = new TH1F("displaced_mcp_pt", ";Monte Carlo Stau Decay Product p_{T} [GeV]; Tracking particles / 500.0 GeV", 10, 0, 5000.0);
    TH1F* displaced_mcp_eta = new TH1F("displaced_mcp_eta", ";Monte Carlo Stau Decay Product #eta; Monte Carlo Stau / 0.2", 30, -1.5, 1.5);
    TH1F* displaced_mcp_phi = new TH1F("displaced_mcp_phi", ";Monte Carlo Stau Decay Product #eta; Monte Carlo Stau / 0.2 radians", 16, -3.2, -3.2);
    TH1F* displaced_mcp_d0 = new TH1F("displaced_mcp_d0", ";Monte Carlo Stau Decay Product |d_{0}|; Monte Carlo Stau / 10 mm", 10, 0, 100.0);
    TH1F* displaced_mcp_rxy = new TH1F("displaced_mcp_rxy", ";Monte Carlo Stau Decay Product r_{xy}; Monte Carlo Stau / 10 mm", 25, 0, 500.0);
    TH1F* displaced_mcp_z0 = new TH1F("displaced_mcp_z0", ";Monte Carlo Stau Decay Product z_{0}; Monte Carlo Stau / 25 mm", 20, -250.0, 250.0);

    TH1F* displaced_tp_pt = new TH1F("displaced_tp_pt", ";Displaced Tracking particle p_{T} [GeV]; Displaced tracking particles / 500.0 GeV", 10, 0, 5000.0);
    TH1F* displaced_tp_eta = new TH1F("displaced_tp_eta", ";Displaced Tracking Particle #eta; Displaced tracking particle / 0.2", 15, -1.5, 1.5);
    TH1F* displaced_tp_phi = new TH1F("displaced_tp_phi", ";Displaced Tracking Particle Phi; Displaced tracking particle / 0.2 radians", 16, -3.2, -3.2);
    TH1F* displaced_tp_d0 = new TH1F("displaced_tp_d0", ";Displaced Tracking Tracking Particle d_{0}; Tracking Particle Stau / 25 mm", 20, -250.0, 250.0);
    TH1F* displaced_tp_z0 = new TH1F("displaced_tp_z0", ";Displaced Tracking Tracking Particle z_{0}; Tracking Particle Stau / 25 mm", 20, -250.0, 250.0);

    TH1F* stau_matched_pt = new TH1F("stau_matched_pt", ";Tracking particle p_{T} [GeV]; Tracking particles / 250.0 GeV", 11, 1900.0, 3000.0);
    TH1F* stau_matched_pt_all = new TH1F("stau_matched_pt_all", ";Tracking particle p_{T} [GeV]; Tracking particles (all) / 100.0 GeV", 11, 1900.0, 3000.0);
    TH1F* displaced_matched_pt = new TH1F("displaced_matched_pt", ";Tracking particle p_{T} [GeV]; Tracking particles / 500.0 GeV", 10, 0, 5000.0);
    TH1F* stau_track_pt = new TH1F("stau_track_pt", ";Stau Track p_{T} [GeV]; Tracks / 100.0 GeV", 50, 0, 5000.0);
    TH1F* displaced_track_pt = new TH1F("displaced_track_pt", ";Displaced Stau Decay Product Track p_{T} [GeV]; Tracks / 100.0 GeV", 50, 0, 5000.0);

    TH1F* stau_matched_eta = new TH1F("stau_matched_eta", ";Stau #eta; Tracking particle / 0.2", 15, -1.25, 1.25);
    TH1F* stau_matched_eta_all = new TH1F("stau_matched_eta_all", ";Stau #eta; Tracking particle (all) / 0.1", 50, -2.5, 2.5);
    TH1F* displaced_matched_eta = new TH1F("displaced_matched_eta", ";Tracking Particle #eta; Tracking particle / 0.2", 15, -1.5, 1.5);
    TH1F* stau_track_eta = new TH1F("stau_track_eta", ";Stau Track #eta; Tracks / 0.1", 50, -2.5, 2.5);
    TH1F* displaced_track_eta = new TH1F("displaced_track_eta", ";Displaced Stau Decay Product Track #eta; Tracks / 0.1", 50, -2.5, 2.5);

    TH1F* displaced_matched_phi = new TH1F("displaced_matched_phi", ";Tracking Particle #phi; Tracking particle / 0.4 radians", 16, -3.2, 3.2);

    TH1F* displaced_track_d0 = new TH1F("displaced_track_d0", ";Displaced Stau Decay Product Track d_{0}; Tracking particles / 25 mm", 100, -250.0, 250.0);
    TH1F* displaced_track_z0 = new TH1F("displaced_track_z0", ";Displaced Stau Decay Product Track z_{0}; Tracking particles / 25 mm", 100, -250.0, 250.0);
    TH1F* displaced_matched_d0 = new TH1F("displaced_matched_d0", ";Tracking particle d_{0}; Tracking particles / 25 mm", 20, -250.0, 250.0);
    TH1F* displaced_matched_z0 = new TH1F("displaced_matched_z0", ";Tracking particle z_{0}; Tracking particles / 25 mm", 20, -250.0, 250.0);

    TH1F* displaced_resPt = new TH1F("displaced_resPt", ";Displaced Track p_{T} - Matched p_{T}; Tracks / 100.0 GeV", 100, -5000, 5000);
    TH1F* displaced_resEta = new TH1F("displaced_resEta", ";Displaced Track #eta - Matched #eta; Tracks / 0.01", 20, -0.1, 0.1);
    TH1F* displaced_resz0 = new TH1F("displaced_resz0", ";|Track z_{0} - Matched z_{0}|; Tracks / 0.2 mm", 25, 0, 5);
    TH1F* displaced_resd0 = new TH1F("displaced_resd0", ";|Track d_{0} - Matched d_{0}|; Tracks / 0.1 mm", 20, 0, 2);

    TH1F* displaced_PtRel = new TH1F("displaced_PtRel", ";|((Displaced Track p_{T} - Matched p_{T}) / Matched p_{T})|; Tracks / 0.1", 20, 0, 2);
    TH1F* PtRelvsPt = new TH1F("PtRelvsPt", ";Tracking particle p_{T} [GeV];p_{T} resolution ", 10, 0, 5000.0);

    TH1F* nhits_vs_d0 = new TH1F("nhits_vs_d0", ";d_{0} [mm] ; NHits", 10, 0, 100.0);
    TH1F* nhits_vs_rxy = new TH1F("nhits_vs_rxy", ";r_{xy} [mm]; NHits ", 25, 0, 500.0);

    TH1F* stau_chi2_reduced = new TH1F("stau_chi2_reduced", "; #chi^{2} / ndf; Tracks / 0.5", 20, 0, 10);
    TH1F* displaced_chi2_reduced = new TH1F("displaced_chi2_reduced", "; #chi^{2} / ndf; Tracks / 0.5", 20, 0, 10);

    TH1F* stau_resPt = new TH1F("stau_resPt", ";stau Track p_{T} - Matched p_{T}; Tracks / 100.0 GeV", 100, -5000, 5000);
    TH1F* stau_resEta = new TH1F("stau_resEta", ";stau Track #eta - Matched #eta; Tracks / 0.01", 20, -0.1, 0.1);
    TH1F* stau_true_velocity = new TH1F("stau_true_velocity", "; Stau True Velocity [mm/ns]; Matched Staus / 2 mm/ns", 17, 150, 184);
    TH1F* stau_velo_residuals = new TH1F("stau_velo_residuals", "; Stau True Velocity - Measured Velocity [mm/ns]; Tracks / 0.5 mm/ns", 20, -5, 5);
    TH1F* stau_resz0 = new TH1F("stau_resz0", ";stau Track z_{0} - Matched z_{0}; Tracks / 1 mm", 20, -10, 10);
    TH1F* stau_resd0 = new TH1F("stau_resd0", ";stau Track d_{0} - Matched d_{0}; Tracks / 1 mm", 20, -10, 10);

    TH1F* stau_PtRel = new TH1F("stau_PtRel", ";|((stau Track p_{T} - Matched p_{T}) / Matched p_{T})|; Tracks / 0.5", 10, 0, 5);

    TH1F* fake_track_pt = new TH1F("fake_track_pt", ";Fake Track p_{T} [GeV]; Tracks / 10.0 GeV", 50, 0, 500.0);
    TH1F* fake_track_eta = new TH1F("fake_track_eta", ";Fake Track #eta; Tracks / 0.1", 50, -2.5, 2.5);
    TH1F* fake_track_d0 = new TH1F("fake_track_d0", ";Fake Track d_{0}; Tracking particles / 25 mm", 20, -250.0, 250.0);
    TH1F* fake_track_z0 = new TH1F("fake_track_z0", ";Fake Track z_{0}; Tracking particles / 25 mm", 20, -250.0, 250.0);
    TH1F* fake_track_chi2_reduced = new TH1F("fake_track_chi2_reduced", "; #chi^{2} / ndf; Tracks / 0.5", 20, 0, 10);

    TH1F* displaced_resEta_z0 = new TH1F("displaced_resEta_z0", ";Matched Displaced track #eta ; z_{0} res [mm] ", 20, -10, 10);
    TH1F* displaced_resPt_z0 = new TH1F("displaced_resPt_z0", ";Matched Displaced track p_{T}; z_{0} res [mm] ", 20, -10, 10);
    TH1F* displaced_resz0_z0 = new TH1F("displaced_resz0_z0", "; z{0} [mm]; z_{0} res [mm] ", 20, -10, 10);

    TH1F* displaced_resEta_PtRel = new TH1F("displaced_resEta_PtRel", ";stau Track d_{0} - Matched d_{0}; p_{T} res / p_{T} ", 20, -10, 10);
    TH1F* displaced_resPt_PtRel = new TH1F("displaced_resPt_PtRel", ";Matched Displaced track p_{T}; p_{T} res / p_{T}", 20, -10, 10);
    TH1F* displaced_resz0_PtRel = new TH1F("displaced_resz0_PtRel", ";z_{0} [mm]; p_{T} res / p_{T}", 20, -10, 10);

    TH1F* displaced_resEta_Eta = new TH1F("displaced_resEta_Eta", ";Matched Displaced track #eta; #eta res", 20, -10, 10);
    TH1F* displaced_resPt_Eta = new TH1F("displaced_resPt_Eta", ";Matched Displaced track p_{T}; #eta res", 20, -10, 10);
    TH1F* displaced_resz0_Eta = new TH1F("displaced_resz0_Eta", ";z_{0} [mm]; eta res", 20, -10, 10);

    TH1F* displaced_nhits = new TH1F("displaced_nhits", ";Number of Tracker Hits; Tracks / 1 Hit", 15, 0, 15);
    
    // resolution (z0, pt) per different variables
    // truth level and track level distributions (pt, eta, phi, z0, d0, etc.)
    
    // ----------------------------------------------------------------------------------------------------------------
    // Calculate tracking efficiency, fake/genuine track. rates, etc. and fill histograms

    char ctxt[500];
    TCanvas c;
    TCanvas cLog1;
    cLog1.SetLogy();
    TString dir = rootFileDir + "Plots_" + fileName + "";
    gSystem->mkdir(dir);
    TString DIR = dir + "/";

    // Counters
    int numRecoableStaus = 0;
    int numMatchedStauTracks = 0;
    int nFakeTrack = 0;
    int nEvents = 0; 
    int numRecoableDisplacedParticles = 0;
    int numTotalDisplacedParticles = 0;
    int numMatchedDisplacedTracks = 0;
    int numGoodQualityMatchedDisplacedTracks = 0;
    int displacedChi2RedCut = 5.0;
    int nPtBins = 10;
    int nVeloBins = 15;
    int nNHitsBins_rxy = 25;
    int nNHitsBins_d0 = 10;
    int nStauPrompt = 0;
    std::vector<std::vector<float> > rxyByBin(nNHitsBins_rxy + 1);
    std::vector<std::vector<float> > d0ByBin(nNHitsBins_d0 + 1);
    std::vector<std::vector<float>> ptResByBin(nPtBins);
    std::vector<std::vector<float>> veloResByBin(nVeloBins);

    // Event loop
    for (int iEvt = 0; iEvt < StauMCPs->GetEntries(); ++iEvt) {
      ++nEvents;
      // Load the entry data into the branch
      StauMCPs->GetEntry(iEvt);
      AllStauTracks->GetEntry(iEvt);
      StauDecayProductsMCPs->GetEntry(iEvt);
      AllDecayProductTracks->GetEntry(iEvt);
      AllFakeTracks->GetEntry(iEvt);
      AllHitsActual->GetEntry(iEvt);
      
      std::cout << "hitR size: " << hitR->size() << "\n";
      for (unsigned int iHit = 0; iHit < hitR->size(); ++iHit){
        hit_rz->Fill(hitZ->at(iHit), hitR->at(iHit));
      }

      std::vector<bool > trueReconstructable; // for each stau determines 
      std::vector<bool > trueHasTrack;// FIXME no longer needed with dividing histograms

      
      
      /// PROCESS STAUS FIRST
      for (unsigned int iStau = 0; iStau < mcp_stau_pt->size(); ++iStau){ // Loop through staus in event, will either be 2, 4, or 6
        trueReconstructable.push_back(mcp_stau_track_reconstructable_bool->at(iStau)); // only first two indices will be used for efficiency
        trueHasTrack.push_back(mcp_stau_track_bool->at(iStau));
        if (mcp_stau_pt->size() == 4){ // if any of the staus are either reconstructable or have a track, ensure this is counted in efficiency
          if (mcp_stau_track_reconstructable_bool->at(0) || mcp_stau_track_reconstructable_bool->at(2)){
            trueReconstructable[0] = true;
          }
          if (mcp_stau_track_bool->at(0) || mcp_stau_track_bool->at(2)){
            trueHasTrack[0] = true;
          }
          if (mcp_stau_track_reconstructable_bool->at(1) || mcp_stau_track_reconstructable_bool->at(3)){
            trueReconstructable[1] = true;
          }
          if (mcp_stau_track_bool->at(1) || mcp_stau_track_bool->at(3)){
            trueHasTrack[1] = true;
          }
        }
        else if (mcp_stau_pt->size() == 6){
          if (mcp_stau_track_reconstructable_bool->at(0) || mcp_stau_track_reconstructable_bool->at(2) || mcp_stau_track_reconstructable_bool->at(4)){
            trueReconstructable[0] = true;
          }
          if (mcp_stau_track_bool->at(0) || mcp_stau_track_bool->at(2) || mcp_stau_track_bool->at(4)){
            trueHasTrack[0] = true;
          }
          if (mcp_stau_track_reconstructable_bool->at(1) || mcp_stau_track_reconstructable_bool->at(3) || mcp_stau_track_reconstructable_bool->at(5)){
            trueReconstructable[1] = true;
          }
          if (mcp_stau_track_bool->at(1) || mcp_stau_track_bool->at(3) || mcp_stau_track_bool->at(5)){
            trueHasTrack[1] = true;
          }
        }
        if (prod_stau_vertex_r->at(iStau) < 1.0){
          nStauPrompt++;
          stau_mcp_pt->Fill(mcp_stau_pt->at(iStau));
          stau_mcp_p->Fill(mcp_stau_pt->at(iStau)*std::cosh(mcp_stau_eta->at(iStau)));
          stau_mcp_E->Fill(sqrt(mcp_stau_pt->at(iStau)*std::cosh(mcp_stau_eta->at(iStau)) * mcp_stau_pt->at(iStau)*std::cosh(mcp_stau_eta->at(iStau)) + stau_mass*stau_mass)); // change for diff masses FIXME read in by file name
          stau_mcp_pt_eta->Fill(mcp_stau_eta->at(iStau), mcp_stau_pt->at(iStau));
          stau_mcp_pt_map[fileName]->Fill(mcp_stau_pt->at(iStau));
          stau_mcp_eta->Fill(mcp_stau_eta->at(iStau));
          stau_mcp_eta_map[fileName]->Fill(mcp_stau_eta->at(iStau));
          stau_mcp_phi->Fill(mcp_stau_phi->at(iStau));
          stau_mcp_phi_map[fileName]->Fill(mcp_stau_phi->at(iStau));
          stau_mcp_d0->Fill(mcp_stau_d0->at(iStau));
          stau_mcp_z0->Fill(mcp_stau_z0->at(iStau));
          stau_rxy_decay_map[fileName]->Fill(prod_stau_endpoint_r->at(iStau));
        }
        
        /*if (mcp_stau_pt->size() > 1 && iStau > 1){
          stau_mcp_p->Fill(mcp_stau_pt->at(1)*std::cosh(mcp_stau_eta->at(1)));
          stau_mcp_E->Fill(sqrt(mcp_stau_pt->at(1)*std::cosh(mcp_stau_eta->at(1)) * mcp_stau_pt->at(1)*std::cosh(mcp_stau_eta->at(1)) + 1000.0*1000.0));
          stau_mcp_pt->Fill(mcp_stau_pt->at(1));
          stau_mcp_pt_eta->Fill(mcp_stau_eta->at(1), mcp_stau_pt->at(1));
          stau_mcp_pt_map[fileName]->Fill(mcp_stau_pt->at(1));
          stau_mcp_eta->Fill(mcp_stau_eta->at(1));
          stau_mcp_eta_map[fileName]->Fill(mcp_stau_eta->at(1));
          stau_mcp_phi->Fill(mcp_stau_phi->at(1));
          stau_mcp_phi_map[fileName]->Fill(mcp_stau_phi->at(1));
          stau_mcp_d0->Fill(mcp_stau_d0->at(1));
          stau_mcp_z0->Fill(mcp_stau_z0->at(1));
        }*/
        

        // FIXME fill stau mcp histograms 
      }

      if (trueReconstructable[0]){ // fill reconstructable staus
        numRecoableStaus++;
        stau_tp_pt->Fill(mcp_stau_pt->at(0));
        stau_tp_pt_map[fileName]->Fill(mcp_stau_pt->at(0));
        stau_tp_eta->Fill(mcp_stau_eta->at(0));
        stau_tp_eta_map[fileName]->Fill(mcp_stau_eta->at(0));
        stau_tp_phi->Fill(mcp_stau_phi->at(0));
        stau_tp_phi_map[fileName]->Fill(mcp_stau_phi->at(0));
        stau_tp_d0->Fill(mcp_stau_d0->at(0));
        stau_tp_z0->Fill(mcp_stau_z0->at(0));
      } 
      if (trueReconstructable[1] && trueReconstructable.size() > 1){
        numRecoableStaus++;
        stau_tp_pt->Fill(mcp_stau_pt->at(1));
        stau_tp_pt_map[fileName]->Fill(mcp_stau_pt->at(1));
        stau_tp_eta->Fill(mcp_stau_eta->at(1));
        stau_tp_eta_map[fileName]->Fill(mcp_stau_eta->at(1));
        stau_tp_phi->Fill(mcp_stau_phi->at(1));
        stau_tp_phi_map[fileName]->Fill(mcp_stau_phi->at(1));
        stau_tp_d0->Fill(mcp_stau_d0->at(1));
        stau_tp_z0->Fill(mcp_stau_z0->at(1));
      } 

      std::vector<float > previousEtas;
      for (unsigned int iStauTrack = 0; iStauTrack < LC_stau_pt_match->size(); ++iStauTrack){ 
        // Use std::find to check if the value exists in the vector
        if (isApproximatelyEqualToAny(previousEtas, LC_stau_eta_match->at(iStauTrack), 0.001)) { // FIXME do this with parent / daughter mcp relations later!! 
          //std::cout << "Value " << LC_stau_eta_match->at(iStauTrack) << " exists in the vector." << std::endl;
        }
        else{
          numMatchedStauTracks++;
          stau_matched_pt_map[fileName]->Fill(LC_stau_pt_match->at(iStauTrack));
          stau_matched_pt->Fill(LC_stau_pt_match->at(iStauTrack));
          stau_matched_eta->Fill(LC_stau_eta_match->at(iStauTrack));
          stau_matched_eta_map[fileName]->Fill(LC_stau_eta_match->at(iStauTrack));
          stau_matched_phi_map[fileName]->Fill(LC_stau_phi_match->at(iStauTrack));
        }
        stau_nhits_map[fileName]->Fill(LC_stau_nhits->at(iStauTrack));
        stau_chi2_reduced_map[fileName]->Fill(LC_stau_chi2->at(iStauTrack) / LC_stau_ndf->at(iStauTrack));
        
        stau_matched_pt_all->Fill(LC_stau_pt_match->at(iStauTrack));
        previousEtas.push_back(LC_stau_eta_match->at(iStauTrack));
        stau_track_pt->Fill(LC_stau_track_pt->at(iStauTrack));
        
        int veloBin = stau_nhits_map[fileName]->FindBin(LC_stau_nhits->at(iStauTrack));
        //std::cout << "velobin: " << veloBin << "\n";
        //std::cout << "LC_stau_nhits->at(iStauTrack): " << LC_stau_nhits->at(iStauTrack) << "\n";
        //stau_true_velocity->Fill(LC_stau_true_velo->at(iStauTrack));
        //stau_velo_residuals->Fill(LC_stau_true_velo->at(iStauTrack) - LC_stau_reco_velo->at(iStauTrack));
        //stau_velo_residuals_map[fileName]->Fill(LC_stau_true_velo->at(iStauTrack) - LC_stau_reco_velo->at(iStauTrack));
        //if (veloBin > 0) veloResByBin[veloBin - 1].push_back(LC_stau_velo_res->at(iStauTrack) / LC_stau_true_velo->at(iStauTrack)); // COMMENT OUT WHEN RUNNING OVER OLDER SAMPLES
        
        
        stau_track_pt_map[fileName]->Fill(LC_stau_track_pt->at(iStauTrack));
        
        stau_matched_eta_all->Fill(LC_stau_eta_match->at(iStauTrack));
        stau_track_eta->Fill(LC_stau_track_eta->at(iStauTrack)); // fixme fill other track properties

        stau_resPt->Fill(LC_stau_track_pt->at(iStauTrack) - LC_stau_pt_match->at(iStauTrack));
        stau_resEta->Fill(LC_stau_track_eta->at(iStauTrack) - LC_stau_eta_match->at(iStauTrack));

        stau_PtRel->Fill(abs(LC_stau_track_pt->at(iStauTrack) - LC_stau_pt_match->at(iStauTrack)) / LC_stau_pt_match->at(iStauTrack));

        stau_chi2_reduced->Fill(LC_stau_chi2->at(iStauTrack) / LC_stau_ndf->at(iStauTrack)); 

      }

      /// PROCESS STAU DECAY PRODUCTS

      for (unsigned int iDP = 0; iDP < mcp_daughter_pt->size(); ++iDP){ // Loop through staus in event, will either be 2, 4, or 6
        numTotalDisplacedParticles++;
        displaced_mcp_pt->Fill(mcp_daughter_pt->at(iDP));
        displaced_mcp_eta->Fill(mcp_daughter_eta->at(iDP));
        displaced_mcp_eta_map[fileName]->Fill(mcp_daughter_eta->at(iDP));
        displaced_mcp_phi->Fill(mcp_daughter_phi->at(iDP));
        displaced_mcp_d0->Fill(mcp_daughter_d0->at(iDP));
        displaced_mcp_z0->Fill(mcp_daughter_z0->at(iDP));

        if(mcp_daughter_track_reconstructable_bool->at(iDP)){
          numRecoableDisplacedParticles++;
          displaced_tp_pt->Fill(mcp_daughter_pt->at(iDP));
          displaced_tp_eta->Fill(mcp_daughter_eta->at(iDP));
          
          displaced_tp_phi->Fill(mcp_daughter_phi->at(iDP));
          displaced_tp_d0->Fill(mcp_daughter_d0->at(iDP));
          displaced_tp_z0->Fill(mcp_daughter_z0->at(iDP));
          displaced_tp_eta_map[fileName]->Fill(mcp_daughter_eta->at(iDP));
          displaced_tp_pt_map[fileName]->Fill(mcp_daughter_pt->at(iDP));
          displaced_tp_rxy_map[fileName]->Fill(mcp_daughter_r_vertex->at(iDP));
          displaced_tp_d0_map[fileName]->Fill(abs(mcp_daughter_d0->at(iDP)));
        }

        displaced_mcp_pt_map[fileName]->Fill(mcp_daughter_pt->at(iDP));
        displaced_mcp_d0_map[fileName]->Fill(abs(mcp_daughter_d0->at(iDP)));
        displaced_mcp_rxy_map[fileName]->Fill(abs(mcp_daughter_r_vertex->at(iDP)));

      }
      
      for (unsigned int iDPTrack = 0; iDPTrack < LC_daughter_pt_match->size(); ++iDPTrack){ 
        numMatchedDisplacedTracks++;
        displaced_matched_pt->Fill(LC_daughter_pt_match->at(iDPTrack));
        displaced_matched_eta->Fill(LC_daughter_eta_match->at(iDPTrack));
        displaced_matched_d0->Fill(LC_daughter_d0_match->at(iDPTrack));
        displaced_matched_z0->Fill(LC_daughter_z0_match->at(iDPTrack));
        displaced_matched_phi->Fill(LC_daughter_phi_match->at(iDPTrack));
        
        displaced_track_pt->Fill(LC_daughter_track_pt->at(iDPTrack));
        displaced_track_eta->Fill(LC_daughter_track_eta->at(iDPTrack)); // fixme fill other track properties
        displaced_track_d0->Fill(LC_daughter_d0->at(iDPTrack));
        displaced_track_z0->Fill(LC_daughter_z0->at(iDPTrack));

        displaced_resPt->Fill(LC_daughter_track_pt->at(iDPTrack) - LC_daughter_pt_match->at(iDPTrack));
        displaced_resEta->Fill(LC_daughter_track_eta->at(iDPTrack) - LC_daughter_eta_match->at(iDPTrack));
        displaced_resd0->Fill(abs(LC_daughter_d0->at(iDPTrack) - LC_daughter_d0_match->at(iDPTrack)));
        displaced_resz0->Fill(abs(LC_daughter_z0->at(iDPTrack) - LC_daughter_z0_match->at(iDPTrack)));

        displaced_PtRel->Fill(abs(LC_daughter_track_pt->at(iDPTrack) - LC_daughter_pt_match->at(iDPTrack)) / LC_daughter_pt_match->at(iDPTrack));
        int ptBin = displaced_matched_pt->FindBin(LC_daughter_pt_match->at(iDPTrack));
        int rxyBin = displaced_mcp_rxy->FindBin(LC_daughter_r_vertex_match->at(iDPTrack));
        int d0Bin = displaced_mcp_d0->FindBin(abs(LC_daughter_d0->at(iDPTrack)));
        //std::cout << "rxyBin: " << rxyBin << " and d0Bin: " << d0Bin << "\n";
        rxyByBin[rxyBin - 1].push_back(LC_daughter_nhits->at(iDPTrack));
        d0ByBin[d0Bin - 1].push_back(LC_daughter_nhits->at(iDPTrack));
        ptResByBin[ptBin - 1].push_back(abs(LC_daughter_track_pt->at(iDPTrack) - LC_daughter_pt_match->at(iDPTrack)) / LC_daughter_pt_match->at(iDPTrack));
        
        displaced_chi2_reduced->Fill(LC_daughter_chi2->at(iDPTrack) / LC_daughter_ndf->at(iDPTrack)); 
        if (LC_daughter_chi2->at(iDPTrack) / LC_daughter_ndf->at(iDPTrack) > displacedChi2RedCut){
          numGoodQualityMatchedDisplacedTracks++;
        }

        displaced_nhits->Fill(LC_daughter_nhits->at(iDPTrack));

        displaced_track_pt_map[fileName]->Fill(LC_daughter_track_pt->at(iDPTrack));
        displaced_track_d0_map[fileName]->Fill(LC_daughter_d0->at(iDPTrack));
        displaced_track_z0_map[fileName]->Fill(LC_daughter_z0->at(iDPTrack));
        displaced_track_eta_map[fileName]->Fill(LC_daughter_track_eta->at(iDPTrack));
        displaced_chi2_reduced_map[fileName]->Fill(LC_daughter_chi2->at(iDPTrack) / LC_daughter_ndf->at(iDPTrack)); 
        displaced_nhits_map[fileName]->Fill(LC_daughter_nhits->at(iDPTrack));

        displaced_overlay_track_nhits_map[fileName]->Fill(LC_daughter_nhits->at(iDPTrack)); // FIXME add vertex, IT, OT hits too!
        displaced_overlay_track_pt_map[fileName]->Fill(LC_daughter_track_pt->at(iDPTrack));
        //std::cout << "LC_daughter_track_eta->at(iDPTrack): " << LC_daughter_track_eta->at(iDPTrack) << "\n";
        displaced_overlay_track_eta_map[fileName]->Fill(LC_daughter_track_eta->at(iDPTrack));
        displaced_overlay_track_d0_map[fileName]->Fill(LC_daughter_d0->at(iDPTrack));
        displaced_overlay_track_z0_map[fileName]->Fill(LC_daughter_z0->at(iDPTrack));
        displaced_overlay_track_chi2_reduced_map[fileName]->Fill(LC_daughter_chi2->at(iDPTrack) / LC_daughter_ndf->at(iDPTrack));

        displaced_matched_eta_map[fileName]->Fill(LC_daughter_eta_match->at(iDPTrack));
        displaced_matched_d0_map[fileName]->Fill(abs(LC_daughter_d0_match->at(iDPTrack)));
        displaced_matched_rxy_map[fileName]->Fill(abs(LC_daughter_r_vertex_match->at(iDPTrack)));
        displaced_matched_pt_map[fileName]->Fill(LC_daughter_pt_match->at(iDPTrack));
        displaced_resd0_map[fileName]->Fill(abs(LC_daughter_d0->at(iDPTrack) - LC_daughter_d0_match->at(iDPTrack)));
        displaced_resz0_map[fileName]->Fill(abs(LC_daughter_z0->at(iDPTrack) - LC_daughter_z0_match->at(iDPTrack)));
      }
      //std::cout << "displaced_matched_pt_map[filename]->GetEntries(): " << displaced_matched_pt_map[fileName]->GetEntries() << " for: " << fileName << "\n";
      //std::cout << "fake_pt->size(): " << fake_pt->size() << "\n";
      //std::cout << "fake_pt->size(): " << fake_eta->size() << "\n";
      //std::cout << "fake_pt->size(): " << fake_d0->size() << "\n";
      //std::cout << "fake_pt->size(): " << fake_z0->size() << "\n";
      for (unsigned int iFakeTrk = 0; iFakeTrk < fake_pt->size(); ++iFakeTrk){ // Fill fake track properties
        ++nFakeTrack;
        fake_track_pt->Fill(fake_pt->at(iFakeTrk));
        fake_track_eta->Fill(fake_eta->at(iFakeTrk));
        fake_track_d0->Fill(fake_d0->at(iFakeTrk));
        fake_track_z0->Fill(fake_z0->at(iFakeTrk));
        fake_track_chi2_reduced->Fill(fake_chi2_reduced->at(iFakeTrk));
        fake_track_nhits_map[fileName]->Fill(fake_nhits->at(iFakeTrk));
        fake_track_pt_map[fileName]->Fill(fake_pt->at(iFakeTrk));
        fake_track_eta_map[fileName]->Fill(fake_eta->at(iFakeTrk));
        fake_track_d0_map[fileName]->Fill(fake_d0->at(iFakeTrk));
        fake_track_z0_map[fileName]->Fill(fake_z0->at(iFakeTrk));
        fake_track_chi2_reduced_map[fileName]->Fill(fake_chi2_reduced->at(iFakeTrk));
      }
    

  } // OUTSIDE EVENT LOOP

  for (int bin = 0; bin < nPtBins; ++bin){
    PtRelvsPt->SetBinContent(bin + 1, computeMean(ptResByBin[bin]));
    PtRelvsPt->SetBinError(bin + 1, computeStandardError(ptResByBin[bin], computeStandardDeviation(ptResByBin[bin], computeMean(ptResByBin[bin]))));
  }

  for (int bin = 0; bin < nNHitsBins_rxy; ++bin){
    nhits_vs_rxy->SetBinContent(bin + 1, computeMean(rxyByBin[bin]));
    nhits_vs_rxy->SetBinError(bin + 1, computeStandardError(rxyByBin[bin], computeStandardDeviation(rxyByBin[bin], computeMean(rxyByBin[bin]))));
  }

  for (int bin = 0; bin < nNHitsBins_d0; ++bin){
    nhits_vs_d0->SetBinContent(bin + 1, computeMean(d0ByBin[bin]));
    nhits_vs_d0->SetBinError(bin + 1, computeStandardError(d0ByBin[bin], computeStandardDeviation(d0ByBin[bin], computeMean(d0ByBin[bin]))));
  }
  
  for (int bin = 0; bin < nVeloBins; ++bin){
    //std::cout << "test 1" << "\n";
    stau_velores_map[fileName]->SetBinContent(bin + 1, computeMean(veloResByBin[bin]));
    stau_velores_map[fileName]->SetBinError(bin + 1, computeStandardError(veloResByBin[bin], computeStandardDeviation(veloResByBin[bin], computeMean(veloResByBin[bin]))));
  }

  /// Draw and save histograms
  c.cd();
  stau_mcp_pt->Draw();
  c.SaveAs(DIR + "stau_mcp_pt.pdf");
  stau_mcp_p->Draw();
  c.SaveAs(DIR + "stau_mcp_p.pdf");
  cLog1.cd();
  stau_mcp_E->Draw();
  cLog1.SaveAs(DIR + "stau_mcp_E.pdf");
  stau_true_velocity->Scale(1.0 / stau_true_velocity->Integral());
  stau_true_velocity->Draw("HIST"); 
  cLog1.SaveAs(DIR + "stau_true_velo.pdf");
  c.cd();
  stau_mcp_pt_eta->Draw();
  c.SaveAs(DIR + "stau_mcp_pt_eta.pdf");
  stau_mcp_eta->Draw();
  c.SaveAs(DIR + "stau_mcp_eta.pdf");
  stau_mcp_phi->Draw();
  c.SaveAs(DIR + "stau_mcp_phi.pdf");
  stau_mcp_d0->Draw();
  c.SaveAs(DIR + "stau_mcp_d0.pdf");
  stau_mcp_z0->Draw();
  c.SaveAs(DIR + "stau_mcp_z0.pdf");
  stau_matched_pt->Draw();
  c.SaveAs(DIR + "stau_matched_pt.pdf");
  stau_matched_eta->Draw();
  c.SaveAs(DIR + "stau_matched_eta.pdf");
  stau_track_pt->Draw();
  c.SaveAs(DIR + "stau_track_pt.pdf");

  TCanvas* c2 = new TCanvas("c", "r-z View", 1400, 700);  // Wider canvas for proportional layout
  c2->SetRightMargin(0.15);  // Leave space for η labels
  c2->SetTopMargin(0.15);
  c2->cd();
  hit_rz->Draw();

  std::vector<double> eta_lines = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6};
  double z_max = 2300;
  double r_max = 1600;
  double label_offset = 75;

  TLatex latex;
  latex.SetTextSize(0.03);
  latex.SetTextFont(42);

  for (double eta : eta_lines) {
      double theta = 2.0 * atan(exp(-eta));

      // Draw the η line
      const int N = 500;
      double z_vals[N], r_vals[N];
      for (int i = 0; i < N; ++i) {
          double z = z_max * i / (N - 1);
          z_vals[i] = z;
          r_vals[i] = z * tan(theta);
      }

      TGraph* gr = new TGraph(N, z_vals, r_vals);
      gr->SetLineColor(kGray + 1);
      gr->SetLineStyle(2);
      gr->Draw("L SAME");

      // Draw label: position depends on η
      if (eta <= 1.0) {
          double z = r_max / tan(theta);
          latex.DrawLatex(z, r_max + label_offset, Form("%.1f", eta));
      } else {
          double r = z_max * tan(theta);
          latex.DrawLatex(z_max + label_offset, r, Form("%.1f", eta));
      }
  }

  c2->SaveAs("hitrz_eta_overlay.pdf");
  
  c.cd();
  stau_velo_residuals->Draw();
  c.SaveAs(DIR + "stau_track_velo_residuals.pdf");
  stau_track_eta->Draw();
  c.SaveAs(DIR + "stau_track_eta.pdf");

  stau_chi2_reduced->Draw();
  c.SaveAs(DIR + "stau_chi2_reduced.pdf");

  displaced_matched_pt->Draw();
  c.SaveAs(DIR + "displaced_matched_pt.pdf");
  displaced_matched_eta->Draw();
  c.SaveAs(DIR + "displaced_matched_eta.pdf");
  displaced_track_pt->Draw();
  c.SaveAs(DIR + "displaced_track_pt.pdf");
  /*
  double d0_min = 0, d0_max = 90;  // Adjust based on data
  double r_min = 0, r_max = 110;  // Adjust based on data
  const int nBins = 10;  

  // Histogram to accumulate NHits
  TH1D h_mean_rxy("h_mean", "NHits vs. r_xy", nBins, r_min, r_max);
  TH1D h_mean_d0("h_mean", "NHits vs. d_0", nBins, d0_min, d0_max);

  // Fill histogram
  for (size_t i = 0; i < r_xy.size(); i++) {
    h_mean_rxy.Fill(r_xy[i], nHits[i]);
    h_mean_d0.Fill(d_0[i], nHits[i]);
  }

  // Store binned values
  std::vector<double> bin_centers_rxy, mean_values_rxy, std_devs_rxy;
  std::vector<double> bin_centers_d0, mean_values_d0, std_devs_d0;
    
  for (int i = 1; i <= nBins; i++) {
    double count = h_mean_rxy.GetBinContent(i);
    if (count > 0) { // Avoid empty bins
          bin_centers_rxy.push_back(h_mean_rxy.GetBinCenter(i));
          mean_values_rxy.push_back(h_mean_rxy.GetBinContent(i));  // Mean NHits
          std_devs_rxy.push_back(h_mean_rxy.GetBinError(i));  // Standard deviation
      }
  }

  for (int i = 1; i <= nBins; i++) {
    double count = h_mean_d0.GetBinContent(i);
    if (count > 0) { // Avoid empty bins
        bin_centers_d0.push_back(h_mean_d0.GetBinCenter(i));
        mean_values_d0.push_back(h_mean_d0.GetBinContent(i));  // Mean NHits
        std_devs_d0.push_back(h_mean_d0.GetBinError(i));  // Standard deviation
    }
}


  int n_rxy = bin_centers_rxy.size();
  TGraphErrors nhits_rxy(n_rxy, bin_centers_rxy.data(), mean_values_rxy.data(), nullptr, std_devs_rxy.data());
  

  nhits_rxy.Draw("AP");
  c.SaveAs(DIR + "nhits_vs_rxy.pdf");


  int n_d0 = bin_centers_d0.size();
  TGraphErrors nhits_d0(n_d0, bin_centers_d0.data(), mean_values_d0.data(), nullptr, std_devs_d0.data());
  nhits_d0.Draw("AP");
  c.SaveAs(DIR + "nhits_vs_d0.pdf");
  */
  displaced_track_eta->Draw();
  c.SaveAs(DIR + "displaced_track_eta.pdf");
  displaced_matched_d0->Draw();
  c.SaveAs(DIR + "displaced_matched_d0.pdf");
  displaced_matched_z0->Draw();
  c.SaveAs(DIR + "displaced_matched_z0.pdf");
  displaced_matched_phi->Draw();
  c.SaveAs(DIR + "displaced_matched_phi.pdf");
  displaced_track_d0->Draw();
  c.SaveAs(DIR + "displaced_track_d0.pdf");
  displaced_track_z0->Draw();
  c.SaveAs(DIR + "displaced_track_z0.pdf");
  displaced_chi2_reduced->Draw();
  c.SaveAs(DIR + "displaced_chi2_reduced.pdf");
  displaced_nhits->Draw();
  c.SaveAs(DIR + "displaced_nhits.pdf");

  displaced_mcp_pt->Draw();
  c.SaveAs(DIR + "displaced_mcp_pt.pdf");
  displaced_mcp_eta->Draw();
  c.SaveAs(DIR + "displaced_mcp_eta.pdf");
  displaced_mcp_phi->Draw();
  c.SaveAs(DIR + "displaced_mcp_phi.pdf");
  displaced_mcp_d0->Draw();
  c.SaveAs(DIR + "displaced_mcp_d0.pdf");
  displaced_mcp_z0->Draw();
  c.SaveAs(DIR + "displaced_mcp_z0.pdf");

  displaced_tp_pt->Draw();
  c.SaveAs(DIR + "displaced_tp_pt.pdf");
  displaced_tp_eta->Draw();
  c.SaveAs(DIR + "displaced_tp_eta.pdf");
  displaced_tp_phi->Draw();
  c.SaveAs(DIR + "displaced_tp_phi.pdf");
  displaced_tp_d0->Draw();
  c.SaveAs(DIR + "displaced_tp_d0.pdf");
  displaced_tp_z0->Draw();
  c.SaveAs(DIR + "displaced_tp_z0.pdf");

  displaced_resPt->Draw();
  c.SaveAs(DIR + "displaced_resPt.pdf");
  displaced_resEta->Draw();
  c.SaveAs(DIR + "displaced_resEta.pdf");
  //displaced_resd0->Scale(1.0 / displaced_resd0->GetEntries());
  displaced_resd0->Draw();
  c.SaveAs(DIR + "displaced_resd0.pdf");
  //displaced_resz0->Scale(1.0 / displaced_resz0->GetEntries());
  displaced_resz0->Draw();
  c.SaveAs(DIR + "displaced_resz0.pdf");

  displaced_PtRel->Draw();
  c.SaveAs(DIR + "displaced_PtRel.pdf");

  PtRelvsPt->Draw();
  c.SaveAs(DIR + "PtRelvsPt.pdf");

  nhits_vs_d0->Draw();
  c.SaveAs(DIR + "nhits_vs_d0.pdf");

  nhits_vs_rxy->Draw();
  c.SaveAs(DIR + "nhits_vs_rxy.pdf");

  fake_track_pt->Draw();
  c.SaveAs(DIR + "fake_track_pt.pdf");
  fake_track_eta->Draw();
  c.SaveAs(DIR + "fake_track_eta.pdf");
  fake_track_d0->Draw();
  c.SaveAs(DIR + "fake_track_d0.pdf");
  fake_track_z0->Draw();
  c.SaveAs(DIR + "fake_track_z0.pdf");
  fake_track_chi2_reduced->Draw();
  c.SaveAs(DIR + "fake_track_chi2_reduced.pdf");


  // Calculate efficiency by dividing histograms
  // staus first
  stau_matched_pt->Sumw2();
  stau_tp_pt->Sumw2();
  TH1F* stau_eff_acc_pt = (TH1F*)stau_matched_pt->Clone();
  stau_eff_acc_pt->SetName("eff_pt");
  stau_eff_acc_pt->GetYaxis()->SetTitle("Efficiency");
  stau_eff_acc_pt->Divide(stau_matched_pt, stau_tp_pt, 1.0, 1.0, "B");
  
  stau_eff_acc_pt->SetAxisRange(0, 1.1, "Y");

  stau_eff_acc_pt->Draw();
  stau_eff_acc_pt->Write();
  c.SaveAs(DIR + "stau_eff_acc_pt.pdf");

  stau_matched_eta->Sumw2();
  stau_tp_eta->Sumw2();
  TH1F* stau_eff_acc_eta = (TH1F*)stau_matched_eta->Clone();
  stau_eff_acc_eta->SetName("eff_eta");
  stau_eff_acc_eta->GetYaxis()->SetTitle("Efficiency");
  stau_eff_acc_eta->Divide(stau_matched_eta, stau_tp_eta, 1.0, 1.0, "B");
  
  stau_eff_acc_eta->SetAxisRange(0, 1.1, "Y");

  stau_eff_acc_eta->Draw();
  stau_eff_acc_eta->Write();
  c.SaveAs(DIR + "stau_eff_acc_eta.pdf");

  // next displaced stau decay products
  // with acceptance
  displaced_matched_pt->Sumw2();
  displaced_tp_pt->Sumw2();
  TH1F* displaced_eff_acc_pt = (TH1F*)displaced_matched_pt->Clone();
  displaced_eff_acc_pt->SetName("eff_pt");
  displaced_eff_acc_pt->GetYaxis()->SetTitle("Efficiency");
  displaced_eff_acc_pt->Divide(displaced_matched_pt, displaced_tp_pt, 1.0, 1.0, "B");
  
  displaced_eff_acc_pt->SetAxisRange(0, 1.1, "Y");

  displaced_eff_acc_pt->Draw();
  displaced_eff_acc_pt->Write();
  c.SaveAs(DIR + "displaced_eff_acc_pt.pdf");
  // w/o acceptance
  displaced_matched_pt->Sumw2();
  displaced_mcp_pt->Sumw2();
  TH1F* displaced_eff_pt = (TH1F*)displaced_matched_pt->Clone();
  displaced_eff_pt->SetName("eff_pt");
  displaced_eff_pt->GetYaxis()->SetTitle("Efficiency");
  displaced_eff_pt->Divide(displaced_matched_pt, displaced_mcp_pt, 1.0, 1.0, "B");
  
  displaced_eff_pt->SetAxisRange(0, 1.1, "Y");

  displaced_eff_pt->Draw();
  displaced_eff_pt->Write();
  c.SaveAs(DIR + "displaced_eff_pt.pdf");

  displaced_matched_eta->Sumw2();
  displaced_tp_eta->Sumw2();
  TH1F* displaced_eff_acc_eta = (TH1F*)displaced_matched_eta->Clone();
  displaced_eff_acc_eta->SetName("eff_eta");
  displaced_eff_acc_eta->GetYaxis()->SetTitle("Efficiency");
  displaced_eff_acc_eta->Divide(displaced_matched_eta, displaced_tp_eta, 1.0, 1.0, "B");
  
  displaced_eff_acc_eta->SetAxisRange(0, 1.1, "Y");

  displaced_eff_acc_eta->Draw();
  displaced_eff_acc_eta->Write();
  c.SaveAs(DIR + "displaced_eff_acc_eta.pdf");

  displaced_matched_eta->Sumw2();
  displaced_mcp_eta->Sumw2();
  TH1F* displaced_eff_eta = (TH1F*)displaced_matched_eta->Clone();
  displaced_eff_eta->SetName("eff_eta");
  displaced_eff_eta->GetYaxis()->SetTitle("Efficiency");
  displaced_eff_eta->Divide(displaced_matched_eta, displaced_mcp_eta, 1.0, 1.0, "B");
  
  displaced_eff_eta->SetAxisRange(0, 1.1, "Y");

  displaced_eff_eta->Draw();
  displaced_eff_eta->Write();
  c.SaveAs(DIR + "displaced_eff_eta.pdf");

  

  displaced_matched_d0->Sumw2();
  displaced_tp_d0->Sumw2();
  TH1F* displaced_eff_acc_d0 = (TH1F*)displaced_matched_d0->Clone();
  displaced_eff_acc_d0->SetName("eff_eta");
  displaced_eff_acc_d0->GetYaxis()->SetTitle("Efficiency");
  displaced_eff_acc_d0->Divide(displaced_matched_d0, displaced_tp_d0, 1.0, 1.0, "B");
  
  displaced_eff_acc_d0->SetAxisRange(0, 1.1, "Y");

  displaced_eff_acc_d0->Draw();
  c.SaveAs(DIR + "displaced_eff_acc_d0.pdf");

  displaced_matched_d0->Sumw2();
  displaced_mcp_d0->Sumw2();
  TH1F* displaced_eff_d0 = (TH1F*)displaced_matched_d0->Clone();
  displaced_eff_d0->SetName("eff_eta");
  displaced_eff_d0->GetYaxis()->SetTitle("Efficiency");
  displaced_eff_d0->Divide(displaced_matched_d0, displaced_mcp_d0, 1.0, 1.0, "B");
  
  displaced_eff_d0->SetAxisRange(0, 1.1, "Y");

  displaced_eff_d0->Draw();
  c.SaveAs(DIR + "displaced_eff_d0.pdf");

  displaced_matched_z0->Sumw2();
  displaced_tp_z0->Sumw2();
  TH1F* displaced_eff_acc_z0 = (TH1F*)displaced_matched_z0->Clone();
  displaced_eff_acc_z0->SetName("eff_eta");
  displaced_eff_acc_z0->GetYaxis()->SetTitle("Efficiency");
  displaced_eff_acc_z0->Divide(displaced_matched_z0, displaced_tp_z0, 1.0, 1.0, "B");
  
  displaced_eff_acc_z0->SetAxisRange(0, 1.1, "Y");

  displaced_eff_acc_z0->Draw();
  c.SaveAs(DIR + "displaced_eff_acc_z0.pdf");

  displaced_matched_z0->Sumw2();
  displaced_mcp_z0->Sumw2();
  TH1F* displaced_eff_z0 = (TH1F*)displaced_matched_z0->Clone();
  displaced_eff_z0->SetName("eff_eta");
  displaced_eff_z0->GetYaxis()->SetTitle("Efficiency");
  displaced_eff_z0->Divide(displaced_matched_z0, displaced_mcp_z0, 1.0, 1.0, "B");
  
  displaced_eff_z0->SetAxisRange(0, 1.1, "Y");

  displaced_eff_z0->Draw();
  c.SaveAs(DIR + "displaced_eff_z0.pdf");

  displaced_matched_phi->Sumw2();
  displaced_tp_phi->Sumw2();
  TH1F* displaced_eff_acc_phi = (TH1F*)displaced_matched_phi->Clone();
  displaced_eff_acc_phi->SetName("eff_eta");
  displaced_eff_acc_phi->GetYaxis()->SetTitle("Efficiency");
  displaced_eff_acc_phi->Divide(displaced_matched_phi, displaced_tp_phi, 1.0, 1.0, "B");
  
  displaced_eff_acc_phi->SetAxisRange(0, 1.1, "Y");

  displaced_eff_acc_phi->Draw();
  c.SaveAs(DIR + "displaced_eff_acc_phi.pdf");

  displaced_matched_phi->Sumw2();
  displaced_mcp_phi->Sumw2();
  TH1F* displaced_eff_phi = (TH1F*)displaced_matched_phi->Clone();
  displaced_eff_phi->SetName("eff_eta");
  displaced_eff_phi->GetYaxis()->SetTitle("Efficiency");
  displaced_eff_phi->Divide(displaced_matched_phi, displaced_mcp_phi, 1.0, 1.0, "B");
  
  displaced_eff_phi->SetAxisRange(0, 1.1, "Y");

  displaced_eff_phi->Draw();
  c.SaveAs(DIR + "displaced_eff_phi.pdf");

  
  // Create a canvas to draw the histograms
  TCanvas *c_over = new TCanvas("c", "Overlay Histograms", 800, 600);

  // Draw the first histogram
  displaced_eff_acc_pt->SetLineColor(kRed); // Set color for better distinction
  displaced_eff_acc_pt->Draw();

  // Draw the second histogram on the same canvas
  displaced_eff_pt->SetLineColor(kBlue); // Set a different color
  displaced_eff_pt->Draw("SAME"); // "SAME" option overlays it on the existing canvas

  // Add a legend to distinguish between the histograms
  TLegend *legend = new TLegend(0.7, 0.8, 0.9, 0.9); // Adjust the coordinates as needed
  legend->AddEntry(displaced_eff_acc_pt, "Eff. w/ Acceptance", "l");
  legend->AddEntry(displaced_eff_pt, "Total Eff. ", "l");
  legend->Draw();

  // Save the canvas as a PDF
  c_over->SaveAs(DIR + "displaced_eff_total_acc_pt.pdf");

  // Draw the first histogram
  displaced_eff_acc_phi->SetLineColor(kRed); // Set color for better distinction
  displaced_eff_acc_phi->Draw();

  // Draw the second histogram on the same canvas
  displaced_eff_phi->SetLineColor(kBlue); // Set a different color
  displaced_eff_phi->Draw("SAME"); // "SAME" option overlays it on the existing canvas

  // Add a legend to distinguish between the histograms
  legend->Draw();

  // Save the canvas as a PDF
  c_over->SaveAs(DIR + "displaced_eff_total_acc_phi.pdf");

  // Draw the first histogram
  displaced_eff_acc_eta->SetLineColor(kRed); // Set color for better distinction
  displaced_eff_acc_eta->Draw();

  // Draw the second histogram on the same canvas
  displaced_eff_eta->SetLineColor(kBlue); // Set a different color
  displaced_eff_eta->Draw("SAME"); // "SAME" option overlays it on the existing canvas

  legend->Draw();

  // Save the canvas as a PDF
  c_over->SaveAs(DIR + "displaced_eff_total_acc_eta.pdf");

  // Draw the first histogram
  displaced_eff_acc_d0->SetLineColor(kRed); // Set color for better distinction
  displaced_eff_acc_d0->Draw();

  // Draw the second histogram on the same canvas
  displaced_eff_d0->SetLineColor(kBlue); // Set a different color
  displaced_eff_d0->Draw("SAME"); // "SAME" option overlays it on the existing canvas

  // Add a legend to distinguish between the histograms
  legend->Draw();

  // Save the canvas as a PDF
  c_over->SaveAs(DIR + "displaced_eff_total_acc_d0.pdf");


  // Draw the first histogram
  displaced_eff_acc_z0->SetLineColor(kRed); // Set color for better distinction
  displaced_eff_acc_z0->Draw();

  // Draw the second histogram on the same canvas
  displaced_eff_z0->SetLineColor(kBlue); // Set a different color
  displaced_eff_z0->Draw("SAME"); // "SAME" option overlays it on the existing canvas

  legend->Draw();

  // Save the canvas as a PDF
  c_over->SaveAs(DIR + "displaced_eff_total_acc_z0.pdf");

  std::cout << "---------------------------------------" << "\n";
  std::cout << "--------   SAMPLE SUMMARY     ---------" << "\n";
  std::cout << "---------------------------------------" << "\n";

  std::cout << "numRecoableStaus: " << numRecoableStaus << "\n";
  std::cout << "numMatchedStauTracks: " << numMatchedStauTracks << "\n";
  std::cout << "total stau tracking efficiency: " << float(numMatchedStauTracks) / float(numRecoableStaus) << "\n";
  std::cout << "fake rate: " << float(nFakeTrack) / float(nEvents) << "\n";
  std::cout << "numEvents: " << nEvents << "\n";
  std::cout << "nStauPrompt: " << nStauPrompt << "\n";
  std::cout << "numRecoableDisplacedParticles: " << numRecoableDisplacedParticles << "\n";
  std::cout << "numTotalDisplacedParticles: " << numTotalDisplacedParticles << "\n";
  std::cout << "numMatchedDisplacedTracks: " << numMatchedDisplacedTracks << "\n";
  std::cout << "num good quality (Chi2 red < 5.0) displaced tracks: " << numGoodQualityMatchedDisplacedTracks << "\n";
  std::cout << "ratio of good quality to all matched displaced tracks: " << float(numGoodQualityMatchedDisplacedTracks) / float(numMatchedDisplacedTracks) << "\n";
  std::cout << "total displaced tracking efficiency: " << float(numMatchedDisplacedTracks) / float(numRecoableDisplacedParticles) << "\n";
  std::cout << "total displaced tracking AxE: " << float(numMatchedDisplacedTracks) / float(numTotalDisplacedParticles) << "\n";
  }


  void overlayMultipleSamples(std::vector<TString> fileNames, std::string timing){

    std::string rootFileDir = "recoRootFiles"; 
    if (timing == "tight"){
      rootFileDir = "tightTimingRecoRootFiles/"; 
    }
    if (timing == "medium"){
      rootFileDir = "mediumTimingRecoRootFiles/";
    }
    if (timing == "loose"){
      rootFileDir = "looseTimingRecoRootFiles/";
    }
    if (timing == "mediumFix"){
      rootFileDir = "timingFixMediumRootFiles/";
    }
    if (timing == "ITMedium"){
      rootFileDir = "timingMediumITSeedingRootFiles/";
    }
    if (timing == "ITMediumNoBIB"){
      rootFileDir = "mediumTimingNoBIBITSeedingRecoRootFiles/";
    }
    if (timing == "ITOTMediumNoBIB"){
      rootFileDir = "mediumTimingNoBIBITOTSeedingRecoRootFiles/";
    }
    if (timing == "mediumNoBIB"){
      rootFileDir = "mediumTimingNoBIBRecoRootFilesStripsOn/";
    }
    if (timing == "mediumNoBIBPropBack"){
      rootFileDir = "mediumTimingNoBIBPropBackStripsOn/";
    }
    if (timing == "mediumNoBIBStripsOff"){
      rootFileDir = "mediumTimingNoBIBStripsOff/";
    }
    if (timing == "mediumPrompt"){
      rootFileDir = "MediumTimingPrompt10pBIB/";
    }
    if (timing == "tightPrompt"){
      rootFileDir = "TightTimingPrompt10pBIB/";
    }
    if (timing == "loosePrompt"){
      rootFileDir = "LooseTimingPrompt10pBIB/";
    }
    if (timing == "mediumDisplaced"){
      rootFileDir = "MediumTimingTwoPasses0pBIB/";
    }
    if (timing == "openHouse"){
      rootFileDir = "openHouseRootFiles/";
    }
    if (timing == "medium10pBIBLRTTest"){
      rootFileDir = "DisplacedMedium10pBIBNoCut/";
    }

    

  //////////////////////////////////////////
  /// overlay multiple sample histograms ///
  //////////////////////////////////////////

  int colorCounter = 1;

  TCanvas* cTrackD0 = new TCanvas("cTrackD0", "Overlayed Histograms", 800, 600);
  cTrackD0->SetLogy();
  TCanvas* cTrackPt = new TCanvas("cTrackPt", "Overlayed Histograms", 800, 600);
  TCanvas* cNHits = new TCanvas("cNHits", "Overlayed Histograms", 800, 600);
  TCanvas* cChi2R = new TCanvas("cChi2R", "Overlayed Histograms", 800, 600);
  

  TCanvas* cFakeTrackD0 = new TCanvas("cFakeTrackD0", "Overlayed Histograms", 800, 600);
  cFakeTrackD0->SetLogy();
  TCanvas* cFakeTrackZ0 = new TCanvas("cFakeTrackZ0", "Overlayed Histograms", 800, 600);
  cFakeTrackZ0->SetLogy();
  TCanvas* cFakeTrackPt = new TCanvas("cFakeTrackPt", "Overlayed Histograms", 800, 600);
  cFakeTrackPt->SetLogy();
  TCanvas* cFakeNHits = new TCanvas("cFakeNHits", "Overlayed Histograms", 800, 600);
  cFakeNHits->SetLogy();
  TCanvas* cFakeChi2R = new TCanvas("cFakeChi2R", "Overlayed Histograms", 800, 600);
  cFakeChi2R->SetLogy();
  TCanvas* cFakeEta = new TCanvas("cFakeEta", "Overlayed Histograms", 800, 600);
  cFakeEta->SetLogy();

  TCanvas* cResZ0 = new TCanvas("cResZ0", "Overlayed Histograms", 800, 600);
  TCanvas* cResD0 = new TCanvas("cResD0", "Overlayed Histograms", 800, 600);

  TCanvas* cEffAccD0 = new TCanvas("cEffAccD0", "Overlayed Histograms", 800, 600);
  TCanvas* cEffD0 = new TCanvas("cEffD0", "Overlayed Histograms", 800, 600);
  TCanvas* cEffAccRXY = new TCanvas("cEffAccRXY", "Overlayed Histograms", 800, 600);
  TCanvas* cEffRXY = new TCanvas("cEffRXY", "Overlayed Histograms", 800, 600);
  TCanvas* cEffAccPt = new TCanvas("cEffAccPt", "Overlayed Histograms", 800, 600);
  TCanvas* cEffAccEta = new TCanvas("cEffAccEta", "Overlayed Histograms", 800, 600);
  TCanvas* cEffPt = new TCanvas("cEffPt", "Overlayed Histograms", 800, 600);
  TCanvas* cEffEta = new TCanvas("cEffEta", "Overlayed Histograms", 800, 600);
  TCanvas* cStEffEta = new TCanvas("cStEffEta", "Overlayed Histograms", 800, 600);
  TCanvas* cStEffAccEta = new TCanvas("cStEffAccEta", "Overlayed Histograms", 800, 600);
  TCanvas* cStEffPt = new TCanvas("cStEffPt", "Overlayed Histograms", 800, 600);
  TCanvas* cStEffAccPt = new TCanvas("cStEffAccPt", "Overlayed Histograms", 800, 600);
  TCanvas* cStEffAccPhi = new TCanvas("cStEffAccPhi", "Overlayed Histograms", 800, 600);
  TCanvas* cStTrackPt = new TCanvas("cStTrackPt", "Overlayed Histograms", 800, 600);
  TCanvas* cStVeloRes = new TCanvas("cStVeloRes", "Overlayed Histograms", 800, 600);
  TCanvas* cStNHits = new TCanvas("cStNHits", "Overlayed Histograms", 800, 600);
  TCanvas* cStChi2 = new TCanvas("cStChi2", "Overlayed Histograms", 800, 600);
  TCanvas* cStVeloResid = new TCanvas("cStVeloResid", "Overlayed Histograms", 800, 600);
  

  TLegend* legendFake0 = new TLegend(0.6, 0.8, 0.95, 0.95);
  legendFake0->SetTextSize(0.0215); // Adjust text size
  TLegend* legendFake1 = new TLegend(0.7, 0.8, 0.95, 0.95);
  legendFake1->SetTextSize(0.0215); // Adjust text size
  TLegend* legendFake2 = new TLegend(0.7, 0.8, 0.95, 0.95);
  legendFake2->SetTextSize(0.0215); // Adjust text size
  TLegend* legendFake3 = new TLegend(0.7, 0.8, 0.95, 0.95);
  legendFake3->SetTextSize(0.0215); // Adjust text size
  TLegend* legendFake4 = new TLegend(0.6, 0.8, 0.95, 0.95);
  legendFake4->SetTextSize(0.0215); // Adjust text size
  TLegend* legendFake5 = new TLegend(0.6, 0.8, 0.95, 0.95);
  legendFake5->SetTextSize(0.0215); // Adjust text size

  TLegend* legendLog = new TLegend(0.7, 0.7, 0.95, 0.95);
  legendLog->SetTextSize(0.03); // Adjust text size
  TLegend* legendStau = new TLegend(0.24, 0.95, 0.925, 0.8);
  legendStau->SetTextSize(0.025); // Adjust text size
  legendStau->SetNColumns(2); // Set legend to two columns
  TLegend* legend1 = new TLegend(0.6, 0.6, 0.95, 0.95);
  legend1->SetTextSize(0.05); // Adjust text size
  TLegend* legend2 = new TLegend(0.6, 0.6, 0.95, 0.95);
  legend2->SetTextSize(0.05); // Adjust text size
  TLegend* legend3 = new TLegend(0.6, 0.8, 0.95, 0.95);
  legend3->SetTextSize(0.0215); // Adjust text size
  TCanvas* cRxy = new TCanvas("cRxy", "Overlayed Histograms", 800, 600);
  TCanvas* cStRxy = new TCanvas("cStRxy", "Overlayed Histograms", 800, 600);
  TLegend* legendRxy = new TLegend(0.65, 0.75, 0.95, 0.95);
  legendRxy->SetTextSize(0.04); // Adjust text size
  cRxy->SetLogy();
  cStRxy->SetLogy();
  for (TString fileName : fileNames){
    std::cout << "fileName: " << fileName << "\n";

    TString dir = rootFileDir + "Overlay_" + fileName + "";
    gSystem->mkdir(dir);
    TString DIR = dir + "/";
    TString OverlayDIR = rootFileDir + "Overlay_" + fileName + "/";
    if (colorCounter == fileNames.size()){
      gSystem->mkdir(OverlayDIR);
    }

    // Mcp info
    cRxy->cd();
    displaced_mcp_rxy_map[fileName]->SetLineColor(colorCounter);
    displaced_mcp_rxy_map[fileName]->SetMaximum(500);
    displaced_mcp_rxy_map[fileName]->GetYaxis()->SetTitle("Num. Charged Decay Products");
    displaced_mcp_rxy_map[fileName]->GetXaxis()->SetTitle("Prod. R_{xy} [mm]");
    displaced_mcp_rxy_map[fileName]->GetXaxis()->SetTitleSize(0.05);
    displaced_mcp_rxy_map[fileName]->GetYaxis()->SetTitleSize(0.05);
    if (colorCounter == 1) displaced_mcp_rxy_map[fileName]->Draw("HIST");
    else displaced_mcp_rxy_map[fileName]->Draw("HIST SAME");
    if (colorCounter == fileNames.size()){
      legendRxy->Draw();
      cRxy->SaveAs(OverlayDIR + "displaced_mcp_prod_rxy.pdf");
    }

    cStRxy->cd();
    stau_rxy_decay_map[fileName]->SetLineColor(colorCounter);
    stau_rxy_decay_map[fileName]->SetMaximum(5000);
    stau_rxy_decay_map[fileName]->GetYaxis()->SetTitle("Num. Staus");
    stau_rxy_decay_map[fileName]->GetXaxis()->SetTitle("Decay R_{xy} [mm]");
    stau_rxy_decay_map[fileName]->GetXaxis()->SetTitleSize(0.05);
    stau_rxy_decay_map[fileName]->GetYaxis()->SetTitleSize(0.05);
    legendRxy->AddEntry(stau_rxy_decay_map[fileName], legend_map[fileName], "l");
    if (colorCounter == 1) stau_rxy_decay_map[fileName]->Draw("HIST");
    else stau_rxy_decay_map[fileName]->Draw("HIST SAME");
    if (colorCounter == fileNames.size()){
      legendRxy->Draw();
      cStRxy->SaveAs(OverlayDIR + "stau_mcp_decay_rxy.pdf");
    }

    // Overlay fake (meant for BIB) and displaced track propertes
    if (fileNames.size() == 1){ // Only do if working with one sample
      cFakeTrackPt->cd();
      fake_track_pt_map[fileName]->SetLineColor(colorCounter);
      displaced_overlay_track_pt_map[fileName]->SetLineColor(colorCounter + 1);
      fake_track_pt_map[fileName]->Scale(1.0 / fake_track_pt_map[fileName]->Integral());
      displaced_overlay_track_pt_map[fileName]->Scale(1.0 / displaced_overlay_track_pt_map[fileName]->Integral());
      fake_track_pt_map[fileName]->GetXaxis()->SetTitle("Track p_{T} [GeV]");
      fake_track_pt_map[fileName]->Draw("HIST");
      displaced_overlay_track_pt_map[fileName]->Draw("HIST SAME"); 
      legendFake0->AddEntry(displaced_overlay_track_pt_map[fileName], ("Displaced Trks, w/ " + std::to_string(static_cast<int>(displaced_overlay_track_pt_map[fileName]->GetEntries())) + " tracks").c_str(), "l");
      legendFake0->AddEntry(fake_track_pt_map[fileName], ("Fake Tracks, w/ " + std::to_string(static_cast<int>(fake_track_pt_map[fileName]->GetEntries())) + " tracks").c_str(), "l");
      legendFake0->Draw();
      cFakeTrackPt->SaveAs(OverlayDIR + "fake_displaced_pt.pdf");

      cFakeEta->cd();
      std::cout << "displaced_overlay_track_eta_map[fileName] entries: " << displaced_overlay_track_eta_map[fileName]->GetEntries() << "\n";
      fake_track_eta_map[fileName]->SetLineColor(colorCounter);
      displaced_overlay_track_eta_map[fileName]->SetLineColor(colorCounter + 1);
      fake_track_eta_map[fileName]->Scale(1.0 / fake_track_eta_map[fileName]->Integral());
      displaced_overlay_track_eta_map[fileName]->Scale(1.0 / displaced_overlay_track_eta_map[fileName]->Integral());
      fake_track_eta_map[fileName]->GetXaxis()->SetTitle("Track #eta");
      fake_track_eta_map[fileName]->Draw("HIST");
      displaced_overlay_track_eta_map[fileName]->Draw("HIST SAME");
      legendFake1->AddEntry(displaced_overlay_track_eta_map[fileName], "Displaced Trks", "l");
      legendFake1->AddEntry(fake_track_eta_map[fileName], "Fake Tracks", "l");
      legendFake1->Draw();
      cFakeEta->SaveAs(OverlayDIR + "fake_displaced_eta.pdf");

      cFakeTrackD0->cd();
      fake_track_d0_map[fileName]->SetLineColor(colorCounter);
      displaced_overlay_track_d0_map[fileName]->SetLineColor(colorCounter + 1);
      fake_track_d0_map[fileName]->Scale(1.0 / fake_track_d0_map[fileName]->Integral());
      displaced_overlay_track_d0_map[fileName]->Scale(1.0 / displaced_overlay_track_d0_map[fileName]->Integral());
      fake_track_d0_map[fileName]->GetXaxis()->SetTitle("Track d_{0} [mm]");
      fake_track_d0_map[fileName]->Draw("HIST");
      displaced_overlay_track_d0_map[fileName]->Draw("HIST SAME");
      legendFake2->AddEntry(displaced_overlay_track_d0_map[fileName], "Displaced Trks", "l");
      legendFake2->AddEntry(fake_track_d0_map[fileName], "Fake Tracks", "l");
      legendFake2->Draw();
      cFakeTrackD0->SaveAs(OverlayDIR + "fake_displaced_d0.pdf");

      cFakeTrackZ0->cd();
      fake_track_z0_map[fileName]->SetLineColor(colorCounter);
      displaced_overlay_track_z0_map[fileName]->SetLineColor(colorCounter + 1);
      fake_track_z0_map[fileName]->Scale(1.0 / fake_track_z0_map[fileName]->Integral());
      displaced_overlay_track_z0_map[fileName]->Scale(1.0 / displaced_overlay_track_z0_map[fileName]->Integral());
      fake_track_z0_map[fileName]->GetXaxis()->SetTitle("Track # [mm]");
      fake_track_z0_map[fileName]->Draw("HIST");
      displaced_overlay_track_z0_map[fileName]->Draw("HIST SAME");
      legendFake3->AddEntry(displaced_overlay_track_z0_map[fileName], "Displaced Trks", "l");
      legendFake3->AddEntry(fake_track_z0_map[fileName], "Fake Tracks", "l");
      legendFake3->Draw();
      cFakeTrackZ0->SaveAs(OverlayDIR + "fake_displaced_z0.pdf");

      cFakeNHits->cd();
      fake_track_nhits_map[fileName]->SetLineColor(colorCounter);
      displaced_overlay_track_nhits_map[fileName]->SetLineColor(colorCounter + 1);
      fake_track_nhits_map[fileName]->Scale(1.0 / fake_track_nhits_map[fileName]->Integral());
      displaced_overlay_track_nhits_map[fileName]->Scale(1.0 / displaced_overlay_track_nhits_map[fileName]->Integral());
      fake_track_nhits_map[fileName]->GetXaxis()->SetTitle("Number of Tracker Hits");
      
      fake_track_nhits_map[fileName]->Draw("HIST");
      displaced_overlay_track_nhits_map[fileName]->Draw("HIST SAME");
      legendFake4->AddEntry(displaced_overlay_track_nhits_map[fileName], ("Displaced Trks, w/ " + std::to_string(static_cast<int>(displaced_overlay_track_pt_map[fileName]->GetEntries())) + " tracks").c_str(), "l");
      legendFake4->AddEntry(fake_track_nhits_map[fileName], ("Fake Tracks, w/ " + std::to_string(static_cast<int>(fake_track_pt_map[fileName]->GetEntries())) + " tracks").c_str(), "l");
      legendFake4->Draw();
      cFakeNHits->SaveAs(OverlayDIR + "fake_displaced_nhits.pdf"); // FIXME add n vxd, inner, outer hits

      cFakeChi2R->cd();
      fake_track_chi2_reduced_map[fileName]->SetLineColor(colorCounter);
      fake_track_chi2_reduced_map[fileName]->SetMaximum(1.0);
      fake_track_chi2_reduced_map[fileName]->Scale(1.0 / fake_track_chi2_reduced_map[fileName]->Integral());
      displaced_overlay_track_chi2_reduced_map[fileName]->Scale(1.0 / displaced_overlay_track_chi2_reduced_map[fileName]->Integral());
      
      displaced_overlay_track_chi2_reduced_map[fileName]->SetLineColor(colorCounter + 1);
      fake_track_chi2_reduced_map[fileName]->GetXaxis()->SetTitle("#chi^{2} / ndf");
      
      fake_track_chi2_reduced_map[fileName]->Draw("HIST");
      displaced_overlay_track_chi2_reduced_map[fileName]->Draw("HIST SAME");
      legendFake5->AddEntry(displaced_overlay_track_chi2_reduced_map[fileName], ("Displaced Trks, w/ " + std::to_string(static_cast<int>(displaced_overlay_track_pt_map[fileName]->GetEntries())) + " tracks").c_str(), "l");
      legendFake5->AddEntry(fake_track_chi2_reduced_map[fileName], ("Fake Tracks, w/ " + std::to_string(static_cast<int>(fake_track_pt_map[fileName]->GetEntries())) + " tracks").c_str(), "l");
      legendFake5->Draw();
      cFakeChi2R->SaveAs(OverlayDIR + "fake_displaced_Chi2R.pdf");
    }
    
    cTrackD0->cd();

    // Draw histograms on the same canvas
    displaced_track_d0_map[fileName]->SetLineColor(colorCounter);
    displaced_track_d0_map[fileName]->GetXaxis()->SetTitleSize(0.05);
    displaced_track_d0_map[fileName]->GetYaxis()->SetTitleSize(0.05);
    displaced_track_d0_map[fileName]->SetLineWidth(2); // Set line width to 3
    displaced_track_d0_map[fileName]->SetMarkerSize(1.5); // Set marker size to 1.5
    displaced_track_d0_map[fileName]->SetMaximum(240);

    if (colorCounter == 1) displaced_track_d0_map[fileName]->Draw("HIST");
    else displaced_track_d0_map[fileName]->Draw("HIST SAME");

    legendLog->AddEntry(displaced_track_d0_map[fileName], legend_map[fileName], "l");
    legendLog->Draw();
    if (colorCounter == fileNames.size()) cTrackD0->SaveAs(OverlayDIR + "displaced_track_d0_overlay.pdf");
    
    cTrackPt->cd();
    // Draw histograms on the same canvas
    displaced_track_pt_map[fileName]->SetLineColor(colorCounter);
    displaced_track_pt_map[fileName]->GetXaxis()->SetTitleSize(0.05);
    displaced_track_pt_map[fileName]->GetYaxis()->SetTitleSize(0.05);
    displaced_track_pt_map[fileName]->SetMaximum(60);
    displaced_track_pt_map[fileName]->SetLineWidth(2); // Set line width to 3
    displaced_track_pt_map[fileName]->SetMarkerSize(1.5); // Set marker size to 1.5

    // Draw the first histogram and the rest on the same canvas
    if (colorCounter == 1) displaced_track_pt_map[fileName]->Draw("HIST");
    else displaced_track_pt_map[fileName]->Draw("HIST SAME");

    legend1->AddEntry(displaced_track_pt_map[fileName], legend_map[fileName], "l");

    // Save canvas
    if (colorCounter == fileNames.size()){
      legend1->Draw();
      cTrackPt->SaveAs(OverlayDIR + "displaced_track_pt_overlay.pdf");
    } 

    cNHits->cd();
    // Draw histograms on the same canvas
    displaced_nhits_map[fileName]->SetLineColor(colorCounter);
    displaced_nhits_map[fileName]->GetXaxis()->SetTitleSize(0.05);
    displaced_nhits_map[fileName]->GetYaxis()->SetTitleSize(0.05);
    displaced_nhits_map[fileName]->SetMaximum(300);
    displaced_nhits_map[fileName]->SetLineWidth(2); // Set line width to 3
    displaced_nhits_map[fileName]->SetMarkerSize(1.5); // Set marker size to 1.5

    // Draw the first histogram and the rest on the same canvas
    if (colorCounter == 1) displaced_nhits_map[fileName]->Draw("HIST");
    else displaced_nhits_map[fileName]->Draw("HIST SAME");

    
    if (colorCounter == fileNames.size()){
      legend1->Draw();
      cNHits->SaveAs(OverlayDIR + "displaced_track_nhits_overlay.pdf");
    } 

    cChi2R->cd();
    // Draw histograms on the same canvas
    displaced_chi2_reduced_map[fileName]->SetLineColor(colorCounter);
    displaced_chi2_reduced_map[fileName]->GetXaxis()->SetTitleSize(0.05);
    displaced_chi2_reduced_map[fileName]->GetYaxis()->SetTitleSize(0.05);
    displaced_chi2_reduced_map[fileName]->SetMaximum(160);
    displaced_chi2_reduced_map[fileName]->SetLineWidth(2); // Set line width to 3
    displaced_chi2_reduced_map[fileName]->SetMarkerSize(1.5); // Set marker size to 1.5

    // Draw the first histogram and the rest on the same canvas
    if (colorCounter == 1) displaced_chi2_reduced_map[fileName]->Draw("HIST");
    else displaced_chi2_reduced_map[fileName]->Draw("HIST SAME");

    if (colorCounter == fileNames.size()){
      legend1->Draw();
      cChi2R->SaveAs(OverlayDIR + "displaced_track_chi2red_overlay.pdf");
    } 


    cResZ0->cd();
    displaced_resz0_map[fileName]->SetLineColor(colorCounter);
    displaced_resz0_map[fileName]->GetXaxis()->SetTitleSize(0.05);
    displaced_resz0_map[fileName]->GetYaxis()->SetTitleSize(0.05);
    displaced_resz0_map[fileName]->SetLineWidth(2); // Set line width to 3
    displaced_resz0_map[fileName]->SetMarkerSize(1.5); // Set marker size to 1.5

    // Draw the first histogram and the rest on the same canvas
    if (colorCounter == 1) displaced_resz0_map[fileName]->Draw("HIST");
    else displaced_resz0_map[fileName]->Draw("HIST SAME");
    
    if (colorCounter == fileNames.size()){
      legend1->Draw();
      cResZ0->SaveAs(OverlayDIR + "displaced_resz0_overlay.pdf");
    } 

    cResD0->cd();
    displaced_resd0_map[fileName]->SetLineColor(colorCounter);
    displaced_resd0_map[fileName]->GetXaxis()->SetTitleSize(0.05);
    displaced_resd0_map[fileName]->GetYaxis()->SetTitleSize(0.05);
    displaced_resd0_map[fileName]->SetLineWidth(2); // Set line width to 3
    displaced_resd0_map[fileName]->SetMarkerSize(1.5); // Set marker size to 1.5
    displaced_resd0_map[fileName]->SetMaximum(400);

    if (colorCounter == 1) displaced_resd0_map[fileName]->Draw("HIST");
    else displaced_resd0_map[fileName]->Draw("SAME");

    if (colorCounter == fileNames.size()){
      legend1->Draw();
      cResD0->SaveAs(OverlayDIR + "displaced_resd0_overlay.pdf");
    } 


    // efficiency vs. d0
    cEffAccD0->cd();
    displaced_matched_d0_map[fileName]->Sumw2();
    displaced_tp_d0_map[fileName]->Sumw2();
    std::cout << "entries for : " << fileName << " matched d0: " << displaced_matched_d0_map[fileName]->GetEntries() << "\n";
    std::cout << "entries for : " << fileName << " tp d0: " << displaced_tp_d0_map[fileName]->GetEntries() << "\n";
    displaced_eff_acc_d0_map[fileName] = (TH1F*)displaced_matched_d0_map[fileName]->Clone();
    displaced_eff_acc_d0_map[fileName]->SetName("eff_eta");
    displaced_eff_acc_d0_map[fileName]->GetYaxis()->SetTitle("Efficiency");
    displaced_eff_acc_d0_map[fileName]->Divide(displaced_matched_d0_map[fileName], displaced_tp_d0_map[fileName], 1.0, 1.0, "B");
    displaced_eff_acc_d0_map[fileName]->SetAxisRange(0, 1.1, "Y");
    displaced_eff_acc_d0_map[fileName]->SetLineColor(colorCounter);
    displaced_eff_acc_d0_map[fileName]->GetXaxis()->SetTitleSize(0.05);
    displaced_eff_acc_d0_map[fileName]->GetYaxis()->SetTitleSize(0.05);
    displaced_eff_acc_d0_map[fileName]->SetLineWidth(2); // Set line width to 3
    displaced_eff_acc_d0_map[fileName]->SetMarkerSize(1.5); // Set marker size to 1.5
    if (colorCounter == 1) displaced_eff_acc_d0_map[fileName]->Draw();
    else displaced_eff_acc_d0_map[fileName]->Draw("SAME");
    
    legend2->AddEntry(displaced_eff_acc_d0_map[fileName], legend_map[fileName], "l");
    
    if (colorCounter == fileNames.size()){
      legend2->Draw();
      cEffAccD0->SaveAs(OverlayDIR + "displaced_eff_acc_d0_overlay.pdf");
    } 

    cEffD0->cd();
    displaced_matched_d0_map[fileName]->Sumw2();
    displaced_mcp_d0_map[fileName]->Sumw2();
    std::cout << "entries for : " << fileName << " matched d0: " << displaced_matched_d0_map[fileName]->GetEntries() << "\n";
    std::cout << "entries for : " << fileName << " mcp d0: " << displaced_mcp_d0_map[fileName]->GetEntries() << "\n";
    displaced_eff_d0_map[fileName] = (TH1F*)displaced_matched_d0_map[fileName]->Clone();
    displaced_eff_d0_map[fileName]->SetName("eff_d0");
    displaced_eff_d0_map[fileName]->GetYaxis()->SetTitle("Efficiency");
    displaced_eff_d0_map[fileName]->Divide(displaced_matched_d0_map[fileName], displaced_mcp_d0_map[fileName], 1.0, 1.0, "B");
    displaced_eff_d0_map[fileName]->SetAxisRange(0, 1.1, "Y");
    displaced_eff_d0_map[fileName]->SetLineColor(colorCounter);
    displaced_eff_d0_map[fileName]->GetXaxis()->SetTitleSize(0.05);
    displaced_eff_d0_map[fileName]->GetYaxis()->SetTitleSize(0.05);
    displaced_eff_d0_map[fileName]->SetLineWidth(2); // Set line width to 3
    displaced_eff_d0_map[fileName]->SetMarkerSize(1.5); // Set marker size to 1.5
    if (colorCounter == 1) displaced_eff_d0_map[fileName]->Draw();
    else displaced_eff_d0_map[fileName]->Draw("SAME");
    
    if (colorCounter == fileNames.size()){
      legend2->Draw();
      cEffD0->SaveAs(OverlayDIR + "displaced_eff_d0_overlay.pdf");
    } 

    cEffAccRXY->cd();
    // efficiency vs. rxy
    displaced_matched_rxy_map[fileName]->Sumw2();
    displaced_tp_rxy_map[fileName]->Sumw2();
    displaced_eff_acc_rxy_map[fileName] = (TH1F*)displaced_matched_rxy_map[fileName]->Clone();
    displaced_eff_acc_rxy_map[fileName]->SetName("eff_eta");
    displaced_eff_acc_rxy_map[fileName]->GetYaxis()->SetTitle("Efficiency");
    displaced_eff_acc_rxy_map[fileName]->Divide(displaced_matched_rxy_map[fileName], displaced_tp_rxy_map[fileName], 1.0, 1.0, "B");
    displaced_eff_acc_rxy_map[fileName]->SetAxisRange(0, 1.1, "Y");
    displaced_eff_acc_rxy_map[fileName]->SetLineColor(colorCounter);
    displaced_eff_acc_rxy_map[fileName]->GetXaxis()->SetTitleSize(0.05);
    displaced_eff_acc_rxy_map[fileName]->GetYaxis()->SetTitleSize(0.05);
    displaced_eff_acc_rxy_map[fileName]->SetLineWidth(2); // Set line width to 3
    displaced_eff_acc_rxy_map[fileName]->SetMarkerSize(1.5); // Set marker size to 1.5
    if (colorCounter == 1) displaced_eff_acc_rxy_map[fileName]->Draw();
    else displaced_eff_acc_rxy_map[fileName]->Draw("SAME");
    
   
    if (colorCounter == fileNames.size()){
      legend2->Draw();
      cEffAccRXY->SaveAs(OverlayDIR + "displaced_eff_acc_rxy_overlay.pdf");
    } 

    cEffRXY->cd();
    displaced_matched_rxy_map[fileName]->Sumw2();
    displaced_mcp_rxy_map[fileName]->Sumw2();
    displaced_eff_rxy_map[fileName] = (TH1F*)displaced_matched_rxy_map[fileName]->Clone();
    displaced_eff_rxy_map[fileName]->SetName("eff_eta");
    displaced_eff_rxy_map[fileName]->GetYaxis()->SetTitle("Efficiency");
    displaced_eff_rxy_map[fileName]->Divide(displaced_matched_rxy_map[fileName], displaced_mcp_rxy_map[fileName], 1.0, 1.0, "B");
    displaced_eff_rxy_map[fileName]->SetAxisRange(0, 1.1, "Y");
    displaced_eff_rxy_map[fileName]->SetLineColor(colorCounter);
    displaced_eff_rxy_map[fileName]->GetXaxis()->SetTitleSize(0.05);
    displaced_eff_rxy_map[fileName]->GetYaxis()->SetTitleSize(0.05);
    displaced_eff_rxy_map[fileName]->SetLineWidth(2); // Set line width to 3
    displaced_eff_rxy_map[fileName]->SetMarkerSize(1.5); // Set marker size to 1.5
    if (colorCounter == 1) displaced_eff_rxy_map[fileName]->Draw();
    else displaced_eff_rxy_map[fileName]->Draw("SAME");
    if (colorCounter == fileNames.size()){
      legend2->Draw();
      cEffRXY->SaveAs(OverlayDIR + "displaced_eff_rxy_overlay.pdf");
    } 
    cEffAccPt->cd();
    displaced_matched_pt_map[fileName]->Sumw2();
    displaced_tp_pt_map[fileName]->Sumw2();
    displaced_eff_acc_pt_map[fileName] = (TH1F*)displaced_matched_pt_map[fileName]->Clone();
    displaced_eff_acc_pt_map[fileName]->SetName("eff_eta");
    displaced_eff_acc_pt_map[fileName]->GetYaxis()->SetTitle("Efficiency");
    displaced_eff_acc_pt_map[fileName]->Divide(displaced_matched_pt_map[fileName], displaced_tp_pt_map[fileName], 1.0, 1.0, "B");
    displaced_eff_acc_pt_map[fileName]->SetAxisRange(0, 1.1, "Y");
    displaced_eff_acc_pt_map[fileName]->SetLineColor(colorCounter);
    displaced_eff_acc_pt_map[fileName]->GetXaxis()->SetTitleSize(0.05);
    displaced_eff_acc_pt_map[fileName]->GetYaxis()->SetTitleSize(0.05);
    displaced_eff_acc_pt_map[fileName]->SetLineWidth(2); // Set line width to 3
    displaced_eff_acc_pt_map[fileName]->SetMarkerSize(1.5); // Set marker size to 1.5
    if (colorCounter == 1) displaced_eff_acc_pt_map[fileName]->Draw();
    else displaced_eff_acc_pt_map[fileName]->Draw("SAME");
    legend3->AddEntry(displaced_eff_acc_pt_map[fileName], legend_map[fileName], "l");
    if (colorCounter == fileNames.size()){
      legend3->Draw();
      cEffAccPt->SaveAs(OverlayDIR + "displaced_eff_acc_pt_overlay.pdf");
    } 

    cEffPt->cd();

    displaced_matched_pt_map[fileName]->Sumw2();
    displaced_mcp_pt_map[fileName]->Sumw2();
    displaced_eff_pt_map[fileName] = (TH1F*)displaced_matched_pt_map[fileName]->Clone();
    displaced_eff_pt_map[fileName]->SetName("eff_eta");
    displaced_eff_pt_map[fileName]->GetYaxis()->SetTitle("Efficiency");
    displaced_eff_pt_map[fileName]->Divide(displaced_matched_pt_map[fileName], displaced_mcp_pt_map[fileName], 1.0, 1.0, "B");
    displaced_eff_pt_map[fileName]->SetAxisRange(0, 1.1, "Y");
    displaced_eff_pt_map[fileName]->SetLineColor(colorCounter);
    displaced_eff_pt_map[fileName]->GetXaxis()->SetTitleSize(0.05);
    displaced_eff_pt_map[fileName]->GetYaxis()->SetTitleSize(0.05);
    displaced_eff_pt_map[fileName]->SetLineWidth(2); // Set line width to 3
    displaced_eff_pt_map[fileName]->SetMarkerSize(1.5); // Set marker size to 1.5
    if (colorCounter == 1) displaced_eff_pt_map[fileName]->Draw();
    else displaced_eff_pt_map[fileName]->Draw("SAME");
    if (colorCounter == fileNames.size()){
      legend3->Draw();
      cEffPt->SaveAs(OverlayDIR + "displaced_eff_pt_overlay.pdf");
    } 

    cEffAccEta->cd();
    displaced_matched_eta_map[fileName]->Sumw2();
    displaced_tp_eta_map[fileName]->Sumw2();
    displaced_eff_acc_eta_map[fileName] = (TH1F*)displaced_matched_eta_map[fileName]->Clone();
    displaced_eff_acc_eta_map[fileName]->SetName("eff_eta");
    displaced_eff_acc_eta_map[fileName]->GetYaxis()->SetTitle("Efficiency");
    displaced_eff_acc_eta_map[fileName]->Divide(displaced_matched_eta_map[fileName], displaced_tp_eta_map[fileName], 1.0, 1.0, "B");
    displaced_eff_acc_eta_map[fileName]->SetAxisRange(0, 1.1, "Y");
    displaced_eff_acc_eta_map[fileName]->SetLineColor(colorCounter);
    displaced_eff_acc_eta_map[fileName]->GetXaxis()->SetTitleSize(0.05);
    displaced_eff_acc_eta_map[fileName]->GetYaxis()->SetTitleSize(0.05);
    displaced_eff_acc_eta_map[fileName]->SetLineWidth(2); // Set line width to 3
    displaced_eff_acc_eta_map[fileName]->SetMarkerSize(1.5); // Set marker size to 1.5
    if (colorCounter == 1) displaced_eff_acc_eta_map[fileName]->Draw();
    else displaced_eff_acc_eta_map[fileName]->Draw("SAME");
    if (colorCounter == fileNames.size()){
      legend3->Draw();
      cEffAccEta->SaveAs(OverlayDIR + "displaced_eff_acc_eta_overlay.pdf");
    } 

    cEffEta->cd();

    displaced_matched_eta_map[fileName]->Sumw2();
    displaced_mcp_eta_map[fileName]->Sumw2();
    displaced_eff_eta_map[fileName] = (TH1F*)displaced_matched_eta_map[fileName]->Clone();
    displaced_eff_eta_map[fileName]->SetName("eff_eta");
    displaced_eff_eta_map[fileName]->GetYaxis()->SetTitle("Efficiency");
    displaced_eff_eta_map[fileName]->Divide(displaced_matched_eta_map[fileName], displaced_mcp_eta_map[fileName], 1.0, 1.0, "B");
    displaced_eff_eta_map[fileName]->SetAxisRange(0, 1.1, "Y");
    displaced_eff_eta_map[fileName]->SetLineColor(colorCounter);
    displaced_eff_eta_map[fileName]->GetXaxis()->SetTitleSize(0.05);
    displaced_eff_eta_map[fileName]->GetYaxis()->SetTitleSize(0.05);
    displaced_eff_eta_map[fileName]->SetLineWidth(2); // Set line width to 3
    displaced_eff_eta_map[fileName]->SetMarkerSize(1.5); // Set marker size to 1.5
    if (colorCounter == 1) displaced_eff_eta_map[fileName]->Draw();
    else displaced_eff_eta_map[fileName]->Draw("SAME");
    if (colorCounter == fileNames.size()){
      legend3->Draw();
      cEffEta->SaveAs(OverlayDIR + "displaced_eff_eta_overlay.pdf");
    } 

    int lineStyle = 1; // Default to solid
    int color = kBlack; // Default color

    std::string strFileName = fileName.Data();

    // Determine color
    if (strFileName.find("4000") != std::string::npos) {
      color = kBlue;
    } else if (strFileName.find("2500") != std::string::npos) {
      color = kRed;
    }

    // Determine line style
    if (strFileName.find("tight") != std::string::npos) {
      lineStyle = 1; // Solid
    } else if (strFileName.find("medium") != std::string::npos) {
        lineStyle = 2; // Dashed
    } else if (strFileName.find("loose") != std::string::npos) {
        lineStyle = 3; // Dotted
    }

    

    cStEffEta->cd();
    stau_matched_eta_map[fileName]->Sumw2();
    stau_mcp_eta_map[fileName]->Sumw2();
    stau_eff_eta_map[fileName] = (TH1F*)stau_matched_eta_map[fileName]->Clone();
    stau_eff_eta_map[fileName]->SetName("st_eff_eta");
    stau_eff_eta_map[fileName]->GetYaxis()->SetTitle("Efficiency");
    stau_eff_eta_map[fileName]->Divide(stau_matched_eta_map[fileName], stau_mcp_eta_map[fileName], 1.0, 1.0, "B");
    stau_eff_eta_map[fileName]->SetAxisRange(0, 1.1, "Y");
    stau_eff_eta_map[fileName]->SetLineColor(color);
    stau_eff_eta_map[fileName]->SetLineStyle(lineStyle);
    stau_eff_eta_map[fileName]->GetXaxis()->SetTitleSize(0.05);
    stau_eff_eta_map[fileName]->GetYaxis()->SetTitleSize(0.05);
    stau_eff_eta_map[fileName]->SetLineWidth(2); // Set line width to 3
    stau_eff_eta_map[fileName]->SetMarkerSize(1.5); // Set marker size to 1.5
    if (colorCounter == 1) stau_eff_eta_map[fileName]->Draw();
    else stau_eff_eta_map[fileName]->Draw("SAME");
    legendStau->AddEntry(stau_eff_eta_map[fileName], legend_map[fileName], "l");
    if (colorCounter == fileNames.size()){
      legendStau->Draw();
      cStEffEta->SaveAs(OverlayDIR + "stau_eff_eta_overlay.pdf");
    } 

    cStEffAccEta->cd();
    stau_matched_eta_map[fileName]->Sumw2();
    stau_tp_eta_map[fileName]->Sumw2();
    stau_eff_acc_eta_map[fileName] = (TH1F*)stau_matched_eta_map[fileName]->Clone();
    stau_eff_acc_eta_map[fileName]->SetName("st_eff_eta");
    stau_eff_acc_eta_map[fileName]->GetYaxis()->SetTitle("Efficiency");
    stau_eff_acc_eta_map[fileName]->Divide(stau_matched_eta_map[fileName], stau_tp_eta_map[fileName], 1.0, 1.0, "B");
    stau_eff_acc_eta_map[fileName]->SetAxisRange(0, 1.3, "Y");

    stau_eff_acc_eta_map[fileName]->SetLineColor(color);
    stau_eff_acc_eta_map[fileName]->SetLineStyle(lineStyle);
    stau_eff_acc_eta_map[fileName]->GetXaxis()->SetTitleSize(0.05);
    stau_eff_acc_eta_map[fileName]->GetYaxis()->SetTitleSize(0.05);
    stau_eff_acc_eta_map[fileName]->SetLineWidth(2); // Set line width to 3
    stau_eff_acc_eta_map[fileName]->SetMarkerSize(1.5); // Set marker size to 1.5
    if (colorCounter == 1) stau_eff_acc_eta_map[fileName]->Draw();
    else stau_eff_acc_eta_map[fileName]->Draw("SAME");
    if (colorCounter == fileNames.size()){
      legendStau->Draw();
      cStEffAccEta->SaveAs(OverlayDIR + "stau_eff_acc_eta_overlay.pdf");
    } 

    cStEffPt->cd();

    stau_matched_pt_map[fileName]->Sumw2();
    stau_mcp_pt_map[fileName]->Sumw2();
    stau_eff_pt_map[fileName] = (TH1F*)stau_matched_pt_map[fileName]->Clone();
    stau_eff_pt_map[fileName]->SetName("st_eff_pt");
    stau_eff_pt_map[fileName]->GetYaxis()->SetTitle("Efficiency");
    stau_eff_pt_map[fileName]->Divide(stau_matched_pt_map[fileName], stau_mcp_pt_map[fileName], 1.0, 1.0, "B");
    stau_eff_pt_map[fileName]->SetAxisRange(0, 1.1, "Y");
    stau_eff_pt_map[fileName]->SetLineColor(color);
    stau_eff_pt_map[fileName]->GetXaxis()->SetTitleSize(0.05);
    stau_eff_pt_map[fileName]->GetYaxis()->SetTitleSize(0.05);
    stau_eff_pt_map[fileName]->SetLineWidth(2); // Set line width to 3
    stau_eff_pt_map[fileName]->SetMarkerSize(1.5); // Set marker size to 1.5
    if (colorCounter == 1) stau_eff_pt_map[fileName]->Draw();
    else stau_eff_pt_map[fileName]->Draw("SAME");
    if (colorCounter == fileNames.size()){
      legendStau->Draw();
      cStEffPt->SaveAs(OverlayDIR + "stau_eff_pt_overlay.pdf");
    } 
    
    cStEffAccPt->cd();

    stau_matched_pt_map[fileName]->Sumw2();
    stau_tp_pt_map[fileName]->Sumw2();
    stau_eff_acc_pt_map[fileName] = (TH1F*)stau_matched_pt_map[fileName]->Clone();
    stau_eff_acc_pt_map[fileName]->SetName("eff_acc_pt");
    stau_eff_acc_pt_map[fileName]->GetYaxis()->SetTitle("Efficiency");
    stau_eff_acc_pt_map[fileName]->Divide(stau_matched_pt_map[fileName], stau_tp_pt_map[fileName], 1.0, 1.0, "B");
    stau_eff_acc_pt_map[fileName]->SetAxisRange(0, 1.1, "Y");
    stau_eff_acc_pt_map[fileName]->SetLineColor(color);
    stau_eff_acc_pt_map[fileName]->SetLineStyle(lineStyle);
    stau_eff_acc_pt_map[fileName]->GetXaxis()->SetTitleSize(0.05);
    stau_eff_acc_pt_map[fileName]->GetYaxis()->SetTitleSize(0.05);
    stau_eff_acc_pt_map[fileName]->SetLineWidth(2); // Set line width to 3
    stau_eff_acc_pt_map[fileName]->SetMarkerSize(1.5); // Set marker size to 1.5
    if (colorCounter == 1) stau_eff_acc_pt_map[fileName]->Draw();
    else stau_eff_acc_pt_map[fileName]->Draw("SAME");
    if (colorCounter == fileNames.size()){
      legendStau->Draw();
      cStEffAccPt->SaveAs(OverlayDIR + "stau_eff_acc_pt_overlay.pdf");
    } 

    cStEffAccPhi->cd();

    stau_matched_phi_map[fileName]->Sumw2();
    stau_tp_phi_map[fileName]->Sumw2();
    stau_eff_acc_phi_map[fileName] = (TH1F*)stau_matched_phi_map[fileName]->Clone();
    stau_eff_acc_phi_map[fileName]->SetName("eff_acc_phi");
    stau_eff_acc_phi_map[fileName]->GetYaxis()->SetTitle("Efficiency");
    stau_eff_acc_phi_map[fileName]->Divide(stau_matched_phi_map[fileName], stau_tp_phi_map[fileName], 1.0, 1.0, "B");
    stau_eff_acc_phi_map[fileName]->SetAxisRange(0, 1.1, "Y");
    stau_eff_acc_phi_map[fileName]->SetLineColor(color);
    stau_eff_acc_phi_map[fileName]->SetLineStyle(lineStyle);
    stau_eff_acc_phi_map[fileName]->GetXaxis()->SetTitleSize(0.05);
    stau_eff_acc_phi_map[fileName]->GetYaxis()->SetTitleSize(0.05);
    stau_eff_acc_phi_map[fileName]->SetLineWidth(2); // Set line width to 3
    stau_eff_acc_phi_map[fileName]->SetMarkerSize(1.5); // Set marker size to 1.5
    if (colorCounter == 1) stau_eff_acc_phi_map[fileName]->Draw();
    else stau_eff_acc_phi_map[fileName]->Draw("SAME");
    if (colorCounter == fileNames.size()){
      legendStau->Draw();
      cStEffAccPhi->SaveAs(OverlayDIR + "stau_eff_acc_phi_overlay.pdf");
    } 


    cStTrackPt->cd();
    // Draw histograms on the same canvas
    stau_track_pt_map[fileName]->SetLineColor(color);
    stau_track_pt_map[fileName]->SetLineStyle(lineStyle);
    stau_track_pt_map[fileName]->GetXaxis()->SetTitleSize(0.05);
    stau_track_pt_map[fileName]->GetYaxis()->SetTitleSize(0.05);
    stau_track_pt_map[fileName]->SetMaximum(625);
    stau_track_pt_map[fileName]->SetLineWidth(2); // Set line width to 3
    stau_track_pt_map[fileName]->SetMarkerSize(1.5); // Set marker size to 1.5

    // Draw the first histogram and the rest on the same canvas
    if (colorCounter == 1) stau_track_pt_map[fileName]->Draw("HIST");
    else stau_track_pt_map[fileName]->Draw("HIST SAME");

    // Save canvas
    if (colorCounter == fileNames.size()){
      legendStau->Draw();
      cStTrackPt->SaveAs(OverlayDIR + "stau_track_pt_overlay.pdf");
    } 

    cStVeloRes->cd();
    // Draw histograms on the same canvas
    stau_velores_map[fileName]->SetLineColor(color);
    stau_velores_map[fileName]->SetLineStyle(lineStyle);
    stau_velores_map[fileName]->GetXaxis()->SetTitleSize(0.05);
    stau_velores_map[fileName]->GetYaxis()->SetTitleSize(0.05);
    stau_velores_map[fileName]->SetMaximum(0.05);
    stau_velores_map[fileName]->SetLineWidth(2); // Set line width to 3
    stau_velores_map[fileName]->SetMarkerSize(1.5); // Set marker size to 1.5

    // Draw the first histogram and the rest on the same canvas
    if (colorCounter == 1) stau_velores_map[fileName]->Draw();
    else stau_velores_map[fileName]->Draw("SAME");

    // Save canvas
    if (colorCounter == fileNames.size()){
      legendStau->Draw();
      cStVeloRes->SaveAs(OverlayDIR + "stau_velores_nhits_overlay.pdf");
    } 

    

    cStNHits->cd();
    // Draw histograms on the same canvas
    stau_nhits_map[fileName]->SetLineColor(color);
    stau_nhits_map[fileName]->SetLineStyle(lineStyle);
    stau_nhits_map[fileName]->GetXaxis()->SetTitleSize(0.05);
    stau_nhits_map[fileName]->GetYaxis()->SetTitleSize(0.05);
    stau_nhits_map[fileName]->SetMaximum(1200);
    stau_nhits_map[fileName]->SetLineWidth(2); // Set line width to 3
    stau_nhits_map[fileName]->SetMarkerSize(1.5); // Set marker size to 1.5

    // Draw the first histogram and the rest on the same canvas
    if (colorCounter == 1) stau_nhits_map[fileName]->Draw("HIST");
    else stau_nhits_map[fileName]->Draw("HIST SAME");
    
    if (colorCounter == fileNames.size()){
      legendStau->Draw();
      cStNHits->SaveAs(OverlayDIR + "stau_track_nhits_overlay.pdf");
    } 

    cStVeloResid->cd();
    // Draw histograms on the same canvas
    stau_velo_residuals_map[fileName]->SetLineColor(color);
    stau_velo_residuals_map[fileName]->SetLineStyle(lineStyle);
    stau_velo_residuals_map[fileName]->GetXaxis()->SetTitleSize(0.05);
    stau_velo_residuals_map[fileName]->GetYaxis()->SetTitleSize(0.05);
    stau_velo_residuals_map[fileName]->SetMaximum(200);
    stau_velo_residuals_map[fileName]->SetLineWidth(2); // Set line width to 3
    stau_velo_residuals_map[fileName]->SetMarkerSize(1.5); // Set marker size to 1.5

    // Draw the first histogram and the rest on the same canvas
    if (colorCounter == 1) stau_velo_residuals_map[fileName]->Draw("HIST");
    else stau_velo_residuals_map[fileName]->Draw("HIST SAME");
    
    if (colorCounter == fileNames.size()){
      legendStau->Draw();
      cStVeloResid->SaveAs(OverlayDIR + "stau_track_velo_residuals_overlay.pdf");
    } 

    colorCounter++; 


    }
    
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
  gStyle->SetTextSize(0.07);
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
  gStyle->SetMarkerSize(1.2);
  gStyle->SetHistLineWidth(2.);
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
  gStyle->SetTickLength(0.02, "X");  // default is usually ~0.03
  gStyle->SetTickLength(0.02, "Y");
}

void mySmallText(Double_t x, Double_t y, Color_t color, char* text) {
  Double_t tsize = 0.044;
  TLatex l;
  l.SetTextSize(tsize);
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x, y, text);
}

bool approximatelyEqual(float a, float b, float epsilon) {
    return std::abs(a - b) <= epsilon;
}

bool isApproximatelyEqualToAny(const std::vector<float>& vec, float value, float epsilon) {
    for (float element : vec) {
        if (approximatelyEqual(element, value, epsilon)) {
            return true;
        }
    }
    return false;
}

void callRootAnalyzer(){
  std::vector<TString> fileNames;
  //if (massPoint == "1TeV") fileNames = {"1000_0.1", "1000_1", "1000_10"};
  //if (massPoint == "2.5TeV") fileNames = {"2500_0.1", "2500_1", "2500_10"};
  //if (massPoint == "4TeV") fileNames = {"4000_0.1", "4000_1", "4000_10"};
  //if (massPoint == "4.5TeV") fileNames = {"4500_0.1", "4500_1", "4500_10"};
  //fileNames = {"2500_1_bib"};
  //std::string timing = "ITMediumNoBIB"; 
  //std::string timing = "ITOTMediumNoBIB"; 
  //std::string timing = "mediumNoBIB"; 
  //std::string timing = "mediumNoBIBStripsOff"; 
  //std::string timing = "mediumNoBIBPropBack"; 
  //std::string timing = "mediumDisplaced"; 
  //std::string timing = "mediumPrompt"; 
  std::string timing = "medium10pBIBLRTTest";
  //fileNames = {"1000_0.1"};
  //fileNames = {"1000_1", "1000_1_osgcomparison"};
  //fileNames = {"1000_1_osgcomparison", "2500_1_osgcomparison", "4000_1_osgcomparison"};
  //fileNames = {"4000_10", "4000_10_osgcomparison"};
  //fileNames = {"1000_1", "2500_1", "4000_1"};
  //fileNames = {"2500_10_tight", "4000_10_tight"};
  //fileNames = {"4000_10_tight", "4000_10_medium", "4000_10_loose"};
  fileNames = {"4000_10"};
  //fileNames = {"2500_1", "4000_1", "2500_10_medium", "4000_10_medium"};
  //fileNames = {"4000_10_tight", "4000_10_medium", "4000_10_loose", "2500_10_tight", "2500_10_medium", "2500_10_loose"};
  initialize_histograms(fileNames);
  for (int i = 0; i < fileNames.size(); i++){
    std::cout << "fileNames[i]: " << fileNames[i];
    rootAnalyzer(fileNames[i], timing);
  }
  if (fileNames.size() > 1){
    overlayMultipleSamples(fileNames, timing);
  }
}
