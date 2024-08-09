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

// Declare global histograms (to overlay multiple samples) FIXME find better way of doing this!! 
TH1F* displaced_track_pt_1_01_ = new TH1F("displaced_track_pt_1_01", ";Stau Decay Product Track p_{T} [GeV]; Tracks / 250.0 GeV", 20, 0, 5000.0);
TH1F* displaced_track_pt_1_1_ = new TH1F("displaced_track_pt_1_1", ";Stau Decay Product Track p_{T} [GeV]; Tracks / 250.0 GeV", 20, 0, 5000.0);
TH1F* displaced_track_pt_45_01_ = new TH1F("displaced_track_pt_45_01", ";Stau Decay Product Track p_{T} [GeV]; Tracks / 250.0 GeV", 20, 0, 5000.0);
TH1F* displaced_track_pt_45_1_ = new TH1F("displaced_track_pt_45_1", ";Stau Decay Product Track p_{T} [GeV]; Tracks / 250.0 GeV", 20, 0, 5000.0);

TH1F* displaced_track_d0_1_01_ = new TH1F("displaced_track_d0_1_01", ";Stau Decay Product Track d_{0} [mm]; Tracking particles / 10 mm", 30, -150.0, 150.0);
TH1F* displaced_track_d0_1_1_ = new TH1F("displaced_track_d0_1_1", ";Stau Decay Product Track d_{0} [mm]; Tracking particles / 10 mm", 30, -150.0, 150.0);
TH1F* displaced_track_d0_45_01_ = new TH1F("displaced_track_d0_45_01", ";Stau Decay Product Track d_{0} [mm]; Tracking particles / 10 mm", 30, -150.0, 150.0);
TH1F* displaced_track_d0_45_1_ = new TH1F("displaced_track_d0_45_1", ";Stau Decay Product Track d_{0} [mm]; Tracking particles / 10 mm", 30, -150.0, 150.0);

TH1F* displaced_chi2_reduced_1_01_ = new TH1F("displaced_chi2_reduced_1_01", "; #chi^{2} / ndf; Tracks / 0.5", 10, 0, 5);
TH1F* displaced_chi2_reduced_1_1_ = new TH1F("displaced_chi2_reduced_1_1", "; #chi^{2} / ndf; Tracks / 0.5", 10, 0, 5);
TH1F* displaced_chi2_reduced_45_01_ = new TH1F("displaced_chi2_reduced_45_01", "; #chi^{2} / ndf; Tracks / 0.5", 10, 0, 5);
TH1F* displaced_chi2_reduced_45_1_ = new TH1F("displaced_chi2_reduced_45_1", "; #chi^{2} / ndf; Tracks / 0.5", 10, 0, 5);

TH1F* displaced_nhits_1_01_ = new TH1F("displaced_nhits_1_01", ";Number of Tracker Hits; Tracks / Hit", 15, 0, 15);
TH1F* displaced_nhits_1_1_ = new TH1F("displaced_nhits_1_1", ";Number of Tracker Hits; Tracks / Hit", 15, 0, 15);
TH1F* displaced_nhits_45_01_ = new TH1F("displaced_nhits_45_01", ";Number of Tracker Hits; Tracks / Hit", 15, 0, 15);
TH1F* displaced_nhits_45_1_ = new TH1F("displaced_nhits_45_1", ";Number of Tracker Hits; Tracks / Hit", 15, 0, 15);

// for efficiency & AxE need matched, tp, mcp FIXME RENAME D0 NOW RXY
TH1F* displaced_matched_d0_1_01_ = new TH1F("displaced_matched_d0_1_01", "; Tracking Particle Vertex r_{xy} [mm]; Tracking particles / 50 mm", 10, 0, 500.0);
TH1F* displaced_matched_d0_1_1_ = new TH1F("displaced_matched_d0_1_1", "; Tracking Particle Vertex r_{xy} [mm]; Tracking particles / 50 mm", 10, 0, 500.0);
TH1F* displaced_matched_d0_45_01_ = new TH1F("displaced_matched_d0_45_01", "; Tracking Particle Vertex r_{xy} [mm]; Tracking particles / 50 mm", 10, 0, 500.0);
TH1F* displaced_matched_d0_45_1_ = new TH1F("displaced_matched_d0_45_1", "; Tracking Particle Vertex r_{xy} [mm]; Tracking particles / 50 mm", 10, 0, 500.0);

TH1F* displaced_tp_d0_1_01_ = new TH1F("displaced_tp_d0_1_01", "; |Displaced Tracking Tracking Particle d_{0}| [mm]; Tracking Particle Stau / 50 mm", 10, 0, 500.0);
TH1F* displaced_tp_d0_1_1_ = new TH1F("displaced_tp_d0_1_1", "; |Displaced Tracking Tracking Particle d_{0}| [mm]; Tracking Particle Stau / 50 mm", 10, 0, 500.0);
TH1F* displaced_tp_d0_45_01_ = new TH1F("displaced_tp_d0_45_01", "; |Displaced Tracking Tracking Particle d_{0}| [mm]; Tracking Particle Stau / 50 mm", 10, 0, 500.0);
TH1F* displaced_tp_d0_45_1_ = new TH1F("displaced_tp_d0_45_1", "; |Displaced Tracking Tracking Particle d_{0}| [mm]; Tracking Particle Stau / 50 mm", 10, 0, 500.0);

TH1F* displaced_mcp_d0_1_01_ = new TH1F("displaced_mcp_d0_1_01", ";Monte Carlo Stau Decay Product d_{0}; Monte Carlo Stau / 50 mm", 10, 0, 500.0);
TH1F* displaced_mcp_d0_1_1_ = new TH1F("displaced_mcp_d0_1_1", ";Monte Carlo Stau Decay Product d_{0}; Monte Carlo Stau / 50 mm", 10, 0, 500.0);
TH1F* displaced_mcp_d0_45_01_ = new TH1F("displaced_mcp_d0_45_01", ";Monte Carlo Stau Decay Product d_{0}; Monte Carlo Stau / 50 mm", 10, 0, 500.0);
TH1F* displaced_mcp_d0_45_1_ = new TH1F("displaced_mcp_d0_45_1", ";Monte Carlo Stau Decay Product d_{0}; Monte Carlo Stau / 50 mm", 10, 0, 500.0);

TH1F* displaced_tp_pt_1_01_ = new TH1F("displaced_tp_pt_1_01", ";Displaced Tracking Particle p_{T} [GeV]; Displaced tracking particles / 500.0 GeV", 10, 0, 5000.0);
TH1F* displaced_tp_pt_1_1_ = new TH1F("displaced_tp_pt_1_1", ";Displaced Tracking Particle p_{T} [GeV]; Displaced tracking particles / 500.0 GeV", 10, 0, 5000.0);
TH1F* displaced_tp_pt_45_01_ = new TH1F("displaced_tp_pt_45_01", ";Displaced Tracking Particle p_{T} [GeV]; Displaced tracking particles / 500.0 GeV", 10, 0, 5000.0);
TH1F* displaced_tp_pt_45_1_ = new TH1F("displaced_tp_pt_45_1", ";Displaced Tracking Particle p_{T} [GeV]; Displaced tracking particles / 500.0 GeV", 10, 0, 5000.0);

TH1F* displaced_matched_pt_1_01_ = new TH1F("displaced_matched_pt_1_01", ";Tracking Particle p_{T} [GeV]; Tracking particles / 500.0 GeV", 10, 0, 5000.0);
TH1F* displaced_matched_pt_1_1_ = new TH1F("displaced_matched_pt_1_1", ";Tracking Particle p_{T} [GeV]; Tracking particles / 500.0 GeV", 10, 0, 5000.0);
TH1F* displaced_matched_pt_45_01_ = new TH1F("displaced_matched_pt_45_01", ";Tracking Particle p_{T} [GeV]; Tracking particles / 500.0 GeV", 10, 0, 5000.0);
TH1F* displaced_matched_pt_45_1_ = new TH1F("displaced_matched_pt_45_1", ";Tracking Particle p_{T} [GeV]; Tracking particles / 500.0 GeV", 10, 0, 5000.0);

TH1F* displaced_mcp_pt_1_01_ = new TH1F("displaced_mcp_pt_1_01", ";Monte Carlo Stau Decay Product p_{T} [GeV]; Tracking particles / 500.0 GeV", 10, 0, 5000.0);
TH1F* displaced_mcp_pt_1_1_ = new TH1F("displaced_mcp_pt_1_1", ";Monte Carlo Stau Decay Product p_{T} [GeV]; Tracking particles / 500.0 GeV", 10, 0, 5000.0);
TH1F* displaced_mcp_pt_45_01_ = new TH1F("displaced_mcp_pt_45_01", ";Monte Carlo Stau Decay Product p_{T} [GeV]; Tracking particles / 500.0 GeV", 10, 0, 5000.0);
TH1F* displaced_mcp_pt_45_1_ = new TH1F("displaced_mcp_pt_45_1", ";Monte Carlo Stau Decay Product p_{T} [GeV]; Tracking particles / 500.0 GeV", 10, 0, 5000.0);

void SetPlotStyle();
void mySmallText(Double_t x, Double_t y, Color_t color, char* text);
bool approximatelyEqual(float a, float b, float epsilon); // helper functions to identify duplicate stau tracks
bool isApproximatelyEqualToAny(const std::vector<float>& vec, float value, float epsilon);
void rootAnalyzer(const TString fileName);

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

void callRootAnalyzer(){
  //std::vector<TString> fileNames = {"1000_0.05_reco", "1000_0.1_reco", "1000_1_reco", "4500_0.1_reco", "4500_1_reco", "4500_10_reco"};
  std::vector<TString> fileNames = {"1000_0.1_reco", "1000_1_reco", "4500_0.1_reco", "4500_1_reco"};
  for (int i = 0; i < fileNames.size(); i++){
    rootAnalyzer(fileNames[i]);
  }
}

// Example Usage: callRootAnalyzer() (change filenames in this file)
// Will output plots as pdfs to a subdirectory named according to the filename

void rootAnalyzer(
const TString fileName
){  
    
    // Example Usage: rootAnalyzer("1000_0.05_reco") 
    // Will output plots as pdfs to a subdirectory named according to the filename
    // ----------------------------------------------------------------------------------------------------------------
    // Set plot style and other configurations
    SetPlotStyle();



    // ----------------------------------------------------------------------------------------------------------------
    // Read NTuples from one file
    TChain* MCPs = new TChain("MCPs");
    MCPs->Add(fileName + "_TIMING32ns64ns.root");
    TChain* StauMCPs = new TChain("StauMCPs");
    StauMCPs->Add(fileName + "_TIMING32ns64ns.root");
    TChain* StauDecayProductsMCPs = new TChain("StauDecayProductsMCPs");
    StauDecayProductsMCPs->Add(fileName + "_TIMING32ns64ns.root");
    TChain* AllTracks = new TChain("AllTracks");
    AllTracks->Add(fileName + "_TIMING32ns64ns.root");
    TChain* AllStauTracks = new TChain("AllStauTracks");
    AllStauTracks->Add(fileName + "_TIMING32ns64ns.root");
    TChain* AllDecayProductTracks = new TChain("AllDecayProductTracks");
    AllDecayProductTracks->Add(fileName + "_TIMING32ns64ns.root");
    TChain* AllFakeTracks = new TChain("AllFakeTracks");
    AllFakeTracks->Add(fileName + "_TIMING32ns64ns.root");
    TChain* AllHits = new TChain("AllHits");
    AllHits->Add(fileName + "_TIMING32ns64ns.root");

    if (MCPs->GetEntries() == 0 || StauMCPs->GetEntries() == 0 || StauDecayProductsMCPs->GetEntries() == 0 || AllTracks->GetEntries() == 0 || AllStauTracks->GetEntries() == 0 || AllDecayProductTracks->GetEntries() == 0 || AllFakeTracks->GetEntries() == 0 || AllHits->GetEntries() == 0) {
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

    // ----------------------------------------------------------------------------------------------------------------
    // Declare histograms to be written to and plotted

    TH1F* stau_mcp_pt = new TH1F("stau_mcp_pt", ";Monte Carlo Stau p_{T} [GeV]; Tracking particles / 250.0 GeV", 10, 0, 2500.0);
    TH1F* stau_mcp_eta = new TH1F("stau_mcp_eta", ";Monte Carlo Stau Eta; Monte Carlo Stau / 0.1", 25, -2.5, 2.5);
    TH1F* stau_mcp_phi = new TH1F("stau_mcp_phi", ";Monte Carlo Stau Eta; Monte Carlo Stau / 0.1 radians", 16, -3.2, -3.2);
    TH1F* stau_mcp_d0 = new TH1F("stau_mcp_d0", ";Monte Carlo Stau d_{0}; Monte Carlo Stau / 10 mm", 50, -250.0, 250.0);
    TH1F* stau_mcp_z0 = new TH1F("stau_mcp_z0", ";Monte Carlo Stau z_{0}; Monte Carlo Stau / 10 mm", 50, -250.0, 250.0);

    TH1F* stau_tp_pt = new TH1F("stau_tp_pt", ";Tracking Particle Stau p_{T} [GeV]; Tracking particles / 250.0 GeV", 10, 0, 2500.0);
    TH1F* stau_tp_eta = new TH1F("stau_tp_eta", ";Tracking Particle Stau Eta; Tracking particle / 0.2", 25, -2.5, 2.5);
    TH1F* stau_tp_phi = new TH1F("stau_tp_phi", ";Tracking Particle Stau Phi; Tracking particle / 0.2 radians", 16, -3.2, -3.2);
    TH1F* stau_tp_d0 = new TH1F("stau_tp_d0", ";Tracking Particle Stau d_{0}; Tracking Particle Stau / 10 mm", 50, -250.0, 250.0);
    TH1F* stau_tp_z0 = new TH1F("stau_tp_z0", ";Tracking Particle Stau z_{0}; Tracking Particle Stau / 10 mm", 50, -250.0, 250.0);

    TH1F* displaced_mcp_pt = new TH1F("displaced_mcp_pt", ";Monte Carlo Stau Decay Product p_{T} [GeV]; Tracking particles / 500.0 GeV", 10, 0, 5000.0);
    TH1F* displaced_mcp_eta = new TH1F("displaced_mcp_eta", ";Monte Carlo Stau Decay Product #eta; Monte Carlo Stau / 0.2", 25, -2.5, 2.5);
    TH1F* displaced_mcp_phi = new TH1F("displaced_mcp_phi", ";Monte Carlo Stau Decay Product #eta; Monte Carlo Stau / 0.2 radians", 16, -3.2, -3.2);
    TH1F* displaced_mcp_d0 = new TH1F("displaced_mcp_d0", ";Monte Carlo Stau Decay Product d_{0}; Monte Carlo Stau / 25 mm", 20, -250.0, 250.0);
    TH1F* displaced_mcp_z0 = new TH1F("displaced_mcp_z0", ";Monte Carlo Stau Decay Product z_{0}; Monte Carlo Stau / 25 mm", 20, -250.0, 250.0);

    TH1F* displaced_tp_pt = new TH1F("displaced_tp_pt", ";Displaced Tracking particle p_{T} [GeV]; Displaced tracking particles / 500.0 GeV", 10, 0, 5000.0);
    TH1F* displaced_tp_eta = new TH1F("displaced_tp_eta", ";Displaced Tracking Particle #eta; Displaced tracking particle / 0.2", 25, -2.5, 2.5);
    TH1F* displaced_tp_phi = new TH1F("displaced_tp_phi", ";Displaced Tracking Particle Phi; Displaced tracking particle / 0.2 radians", 16, -3.2, -3.2);
    TH1F* displaced_tp_d0 = new TH1F("displaced_tp_d0", ";Displaced Tracking Tracking Particle d_{0}; Tracking Particle Stau / 25 mm", 20, -250.0, 250.0);
    TH1F* displaced_tp_z0 = new TH1F("displaced_tp_z0", ";Displaced Tracking Tracking Particle z_{0}; Tracking Particle Stau / 25 mm", 20, -250.0, 250.0);

    TH1F* stau_matched_pt = new TH1F("stau_matched_pt", ";Tracking particle p_{T} [GeV]; Tracking particles / 250.0 GeV", 10, 0, 2500.0);
    TH1F* stau_matched_pt_all = new TH1F("stau_matched_pt_all", ";Tracking particle p_{T} [GeV]; Tracking particles (all) / 100.0 GeV", 30, 2000.0, 5000.0);
    TH1F* displaced_matched_pt = new TH1F("displaced_matched_pt", ";Tracking particle p_{T} [GeV]; Tracking particles / 500.0 GeV", 10, 0, 5000.0);
    TH1F* stau_track_pt = new TH1F("stau_track_pt", ";Stau Track p_{T} [GeV]; Tracks / 10.0 GeV", 50, 0, 5000.0);
    TH1F* displaced_track_pt = new TH1F("displaced_track_pt", ";Displaced Stau Decay Product Track p_{T} [GeV]; Tracks / 100.0 GeV", 50, 0, 5000.0);

    TH1F* stau_matched_eta = new TH1F("stau_matched_eta", ";Stau #eta; Tracking particle / 0.2", 25, -2.5, 2.5);
    TH1F* stau_matched_eta_all = new TH1F("stau_matched_eta_all", ";Stau #eta; Tracking particle (all) / 0.1", 50, -2.5, 2.5);
    TH1F* displaced_matched_eta = new TH1F("displaced_matched_eta", ";Tracking Particle #eta; Tracking particle / 0.2", 25, -2.5, 2.5);
    TH1F* stau_track_eta = new TH1F("stau_track_eta", ";Stau Track #eta; Tracks / 0.1", 50, -2.5, 2.5);
    TH1F* displaced_track_eta = new TH1F("displaced_track_eta", ";Displaced Stau Decay Product Track #eta; Tracks / 0.1", 50, -2.5, 2.5);

    TH1F* displaced_matched_phi = new TH1F("displaced_matched_phi", ";Tracking Particle #phi; Tracking particle / 0.4 radians", 16, -3.2, 3.2);

    TH1F* displaced_track_d0 = new TH1F("displaced_track_d0", ";Displaced Stau Decay Product Track d_{0}; Tracking particles / 25 mm", 100, -250.0, 250.0);
    TH1F* displaced_track_z0 = new TH1F("displaced_track_z0", ";Displaced Stau Decay Product Track z_{0}; Tracking particles / 25 mm", 100, -250.0, 250.0);
    TH1F* displaced_matched_d0 = new TH1F("displaced_matched_d0", ";Tracking particle d_{0}; Tracking particles / 25 mm", 20, -250.0, 250.0);
    TH1F* displaced_matched_z0 = new TH1F("displaced_matched_z0", ";Tracking particle z_{0}; Tracking particles / 25 mm", 20, -250.0, 250.0);

    TH1F* displaced_resPt = new TH1F("displaced_resPt", ";Displaced Track p_{T} - Matched p_{T}; Tracks / 100.0 GeV", 100, -5000, 5000);
    TH1F* displaced_resEta = new TH1F("displaced_resEta", ";Displaced Track #eta - Matched #eta; Tracks / 0.01", 20, -0.1, 0.1);
    TH1F* displaced_resz0 = new TH1F("displaced_resz0", ";Displaced Track z_{0} - Matched z_{0}; Tracks / 1 mm", 20, -10, 10);
    TH1F* displaced_resd0 = new TH1F("displaced_resd0", ";Displaced Track d_{0} - Matched d_{0}; Tracks / 1 mm", 20, -10, 10);

    TH1F* displaced_PtRel = new TH1F("displaced_PtRel", ";|((Displaced Track p_{T} - Matched p_{T}) / Matched p_{T})|; Tracks / 0.1", 20, 0, 2);
    TH1F* PtRelvsPt = new TH1F("PtRelvsPt", ";Tracking particle p_{T} [GeV];p_{T} resolution ", 10, 0, 5000.0);

    TH1F* stau_chi2_reduced = new TH1F("stau_chi2_reduced", "; #chi^{2} / ndf; Tracks / 0.5", 20, 0, 10);
    TH1F* displaced_chi2_reduced = new TH1F("displaced_chi2_reduced", "; #chi^{2} / ndf; Tracks / 0.5", 20, 0, 10);

    TH1F* stau_resPt = new TH1F("stau_resPt", ";stau Track p_{T} - Matched p_{T}; Tracks / 100.0 GeV", 100, -5000, 5000);
    TH1F* stau_resEta = new TH1F("stau_resEta", ";stau Track #eta - Matched #eta; Tracks / 0.01", 20, -0.1, 0.1);
    TH1F* stau_resz0 = new TH1F("stau_resz0", ";stau Track z_{0} - Matched z_{0}; Tracks / 1 mm", 20, -10, 10);
    TH1F* stau_resd0 = new TH1F("stau_resd0", ";stau Track d_{0} - Matched d_{0}; Tracks / 1 mm", 20, -10, 10);

    TH1F* stau_PtRel = new TH1F("stau_PtRel", ";|((stau Track p_{T} - Matched p_{T}) / Matched p_{T})|; Tracks / 0.5", 10, 0, 5);

    TH1F* fake_track_pt = new TH1F("fake_track_pt", ";Fake Track p_{T} [GeV]; Tracks / 100.0 GeV", 50, 0, 5000.0);
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


    //hisRZTTStubs_ = dir.make<TH2F>("RZ TTStubs", ";Stub z; Stub r", 400, -300, 300, 400, 0., 120); taken from iranalyzer

    
    
    // ----------------------------------------------------------------------------------------------------------------
    // Calculate tracking efficiency, fake/genuine track. rates, etc. and fill histograms

    char ctxt[500];
    TCanvas c;
    TString dir = "Plots_" + fileName + "_TIMING32ns64ns";
    gSystem->mkdir(dir);
    TString DIR = dir + "/";

    // Counters
    int numRecoableStaus = 0;
    int numMatchedStauTracks = 0;
    int numRecoableDisplacedParticles = 0;
    int numMatchedDisplacedTracks = 0;
    int nPtBins = 10;
    std::vector<std::vector<float>> ptResByBin(nPtBins);

    // Event loop
    for (int iEvt = 0; iEvt < StauMCPs->GetEntries(); ++iEvt) {
      // Load the entry data into the branch
      StauMCPs->GetEntry(iEvt);
      AllStauTracks->GetEntry(iEvt);
      StauDecayProductsMCPs->GetEntry(iEvt);
      AllDecayProductTracks->GetEntry(iEvt);
      AllFakeTracks->GetEntry(iEvt);

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
        stau_mcp_pt->Fill(mcp_stau_pt->at(0));
        stau_mcp_eta->Fill(mcp_stau_eta->at(0));
        stau_mcp_phi->Fill(mcp_stau_phi->at(0));
        stau_mcp_d0->Fill(mcp_stau_d0->at(0));
        stau_mcp_z0->Fill(mcp_stau_z0->at(0));
        if (mcp_stau_pt->size() > 1){
          stau_mcp_pt->Fill(mcp_stau_pt->at(1));
          stau_mcp_eta->Fill(mcp_stau_eta->at(1));
          stau_mcp_phi->Fill(mcp_stau_phi->at(1));
          stau_mcp_d0->Fill(mcp_stau_d0->at(1));
          stau_mcp_z0->Fill(mcp_stau_z0->at(1));
        }
        

        // FIXME fill stau mcp histograms 
      }

      if (trueReconstructable[0]){ // fill reconstructable staus
        numRecoableStaus++;
        stau_tp_pt->Fill(mcp_stau_pt->at(0));
        stau_tp_eta->Fill(mcp_stau_eta->at(0));
        stau_tp_phi->Fill(mcp_stau_phi->at(0));
        stau_tp_d0->Fill(mcp_stau_d0->at(0));
        stau_tp_z0->Fill(mcp_stau_z0->at(0));
      } 
      if (trueReconstructable[1] && trueReconstructable.size() > 1){
        numRecoableStaus++;
        stau_tp_pt->Fill(mcp_stau_pt->at(1));
        stau_tp_eta->Fill(mcp_stau_eta->at(1));
        stau_tp_phi->Fill(mcp_stau_phi->at(1));
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
          stau_matched_pt->Fill(LC_stau_pt_match->at(iStauTrack));
          stau_matched_eta->Fill(LC_stau_eta_match->at(iStauTrack));
        }
        
        stau_matched_pt_all->Fill(LC_stau_pt_match->at(iStauTrack));
        previousEtas.push_back(LC_stau_eta_match->at(iStauTrack));
        stau_track_pt->Fill(LC_stau_track_pt->at(iStauTrack));
        
        stau_matched_eta_all->Fill(LC_stau_eta_match->at(iStauTrack));
        stau_track_eta->Fill(LC_stau_track_eta->at(iStauTrack)); // fixme fill other track properties

        stau_resPt->Fill(LC_stau_track_pt->at(iStauTrack) - LC_stau_pt_match->at(iStauTrack));
        stau_resEta->Fill(LC_stau_track_eta->at(iStauTrack) - LC_stau_eta_match->at(iStauTrack));

        stau_PtRel->Fill(abs(LC_stau_track_pt->at(iStauTrack) - LC_stau_pt_match->at(iStauTrack)) / LC_stau_pt_match->at(iStauTrack));

        stau_chi2_reduced->Fill(LC_stau_chi2->at(iStauTrack) / LC_stau_ndf->at(iStauTrack)); 

      }

      /// PROCESS STAU DECAY PRODUCTS

      for (unsigned int iDP = 0; iDP < mcp_daughter_pt->size(); ++iDP){ // Loop through staus in event, will either be 2, 4, or 6
        displaced_mcp_pt->Fill(mcp_daughter_pt->at(iDP));
        displaced_mcp_eta->Fill(mcp_daughter_eta->at(iDP));
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
        }
        if(fileName == "1000_0.1_reco"){ // fixme use to fill specific histograms
          if(mcp_daughter_track_reconstructable_bool->at(iDP)){
            displaced_tp_pt_1_01_->Fill(mcp_daughter_pt->at(iDP));
            displaced_tp_d0_1_01_->Fill(mcp_daughter_r_vertex->at(iDP));
          }
          displaced_mcp_pt_1_01_->Fill(mcp_daughter_pt->at(iDP));
          displaced_mcp_d0_1_01_->Fill(mcp_daughter_r_vertex->at(iDP));
        }
        
        else if(fileName == "1000_1_reco"){
          if(mcp_daughter_track_reconstructable_bool->at(iDP)){
            displaced_tp_pt_1_1_->Fill(mcp_daughter_pt->at(iDP));
            displaced_tp_d0_1_1_->Fill(mcp_daughter_r_vertex->at(iDP));
          }
          displaced_mcp_pt_1_1_->Fill(mcp_daughter_pt->at(iDP));
          displaced_mcp_d0_1_1_->Fill(mcp_daughter_r_vertex->at(iDP));
        }
        
        else if(fileName == "4500_0.1_reco"){
          if(mcp_daughter_track_reconstructable_bool->at(iDP)){
            displaced_tp_pt_45_01_->Fill(mcp_daughter_pt->at(iDP));
            displaced_tp_d0_45_01_->Fill(mcp_daughter_r_vertex->at(iDP));
          }
          displaced_mcp_pt_45_01_->Fill(mcp_daughter_pt->at(iDP));
          displaced_mcp_d0_45_01_->Fill(mcp_daughter_r_vertex->at(iDP));
          
        }
        else if(fileName == "4500_1_reco"){
          if(mcp_daughter_track_reconstructable_bool->at(iDP)){
            displaced_tp_pt_45_1_->Fill(mcp_daughter_pt->at(iDP));
            displaced_tp_d0_45_1_->Fill(mcp_daughter_r_vertex->at(iDP));
          }
          displaced_mcp_pt_45_1_->Fill(mcp_daughter_pt->at(iDP));
          displaced_mcp_d0_45_1_->Fill(mcp_daughter_r_vertex->at(iDP));
        }
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
        displaced_resd0->Fill(LC_daughter_d0->at(iDPTrack) - LC_daughter_d0_match->at(iDPTrack));
        displaced_resz0->Fill(LC_daughter_z0->at(iDPTrack) - LC_daughter_z0_match->at(iDPTrack));

        displaced_PtRel->Fill(abs(LC_daughter_track_pt->at(iDPTrack) - LC_daughter_pt_match->at(iDPTrack)) / LC_daughter_pt_match->at(iDPTrack));
        int ptBin = displaced_matched_pt->FindBin(LC_daughter_pt_match->at(iDPTrack));
        ptResByBin[ptBin - 1].push_back(abs(LC_daughter_track_pt->at(iDPTrack) - LC_daughter_pt_match->at(iDPTrack)) / LC_daughter_pt_match->at(iDPTrack));
        
        displaced_chi2_reduced->Fill(LC_daughter_chi2->at(iDPTrack) / LC_daughter_ndf->at(iDPTrack)); 

        displaced_nhits->Fill(LC_daughter_nhits->at(iDPTrack));

        if(fileName == "1000_0.1_reco"){ 
          displaced_track_pt_1_01_->Fill(LC_daughter_track_pt->at(iDPTrack));
          displaced_track_d0_1_01_->Fill(LC_daughter_d0->at(iDPTrack));
          displaced_chi2_reduced_1_01_->Fill(LC_daughter_chi2->at(iDPTrack) / LC_daughter_ndf->at(iDPTrack)); 
          displaced_nhits_1_01_->Fill(LC_daughter_nhits->at(iDPTrack));
          displaced_matched_d0_1_01_->Fill(LC_daughter_r_vertex_match->at(iDPTrack));
          displaced_matched_pt_1_01_->Fill(LC_daughter_pt_match->at(iDPTrack));
        }
        else if(fileName == "1000_1_reco"){
          displaced_track_pt_1_1_->Fill(LC_daughter_track_pt->at(iDPTrack));
          displaced_track_d0_1_1_->Fill(LC_daughter_d0->at(iDPTrack));
          displaced_chi2_reduced_1_1_->Fill(LC_daughter_chi2->at(iDPTrack) / LC_daughter_ndf->at(iDPTrack)); 
          displaced_nhits_1_1_->Fill(LC_daughter_nhits->at(iDPTrack));
          displaced_matched_d0_1_1_->Fill(LC_daughter_r_vertex_match->at(iDPTrack));
          displaced_matched_pt_1_1_->Fill(LC_daughter_pt_match->at(iDPTrack));
        }
        else if(fileName == "4500_0.1_reco"){
          displaced_track_pt_45_01_->Fill(LC_daughter_track_pt->at(iDPTrack));
          displaced_track_d0_45_01_->Fill(LC_daughter_d0->at(iDPTrack));
          displaced_chi2_reduced_45_01_->Fill(LC_daughter_chi2->at(iDPTrack) / LC_daughter_ndf->at(iDPTrack)); 
          displaced_nhits_45_01_->Fill(LC_daughter_nhits->at(iDPTrack));
          displaced_matched_d0_45_01_->Fill(LC_daughter_r_vertex_match->at(iDPTrack));
          displaced_matched_pt_45_01_->Fill(LC_daughter_pt_match->at(iDPTrack));
        }
        else if(fileName == "4500_1_reco"){
          displaced_track_pt_45_1_->Fill(LC_daughter_track_pt->at(iDPTrack));
          displaced_track_d0_45_1_->Fill(LC_daughter_d0->at(iDPTrack));
          displaced_chi2_reduced_45_1_->Fill(LC_daughter_chi2->at(iDPTrack) / LC_daughter_ndf->at(iDPTrack)); 
          displaced_nhits_45_1_->Fill(LC_daughter_nhits->at(iDPTrack));
          displaced_matched_d0_45_1_->Fill(LC_daughter_r_vertex_match->at(iDPTrack));
          displaced_matched_pt_45_1_->Fill(LC_daughter_pt_match->at(iDPTrack));
        }
      }

      for (unsigned int iFakeTrk = 0; iFakeTrk < fake_pt->size(); ++iFakeTrk){ // Fill fake track properties
        fake_track_pt->Fill(fake_pt->at(iFakeTrk));
        fake_track_eta->Fill(fake_eta->at(iFakeTrk));
        fake_track_d0->Fill(fake_d0->at(iFakeTrk));
        fake_track_z0->Fill(fake_z0->at(iFakeTrk));
        fake_track_chi2_reduced->Fill(fake_chi2_reduced->at(iFakeTrk));
      }
    

  } // OUTSIDE EVENT LOOP

  for (int bin = 0; bin < nPtBins; ++bin){
    PtRelvsPt->SetBinContent(bin, computeMean(ptResByBin[bin]));
    PtRelvsPt->SetBinError(bin, computeStandardError(ptResByBin[bin], computeStandardDeviation(ptResByBin[bin], computeMean(ptResByBin[bin]))));
  }

  /// Draw and save histograms
  stau_mcp_pt->Draw();
  c.SaveAs(DIR + "stau_mcp_pt.pdf");
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
  displaced_resd0->Draw();
  c.SaveAs(DIR + "displaced_resd0.pdf");
  displaced_resz0->Draw();
  c.SaveAs(DIR + "displaced_resz0.pdf");

  displaced_PtRel->Draw();
  c.SaveAs(DIR + "displaced_PtRel.pdf");

  PtRelvsPt->Draw();
  c.SaveAs(DIR + "PtRelvsPt.pdf");

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
  legend->AddEntry(displaced_eff_acc_phi, "Eff. w/ Acceptance", "l");
  legend->AddEntry(displaced_eff_phi, "Total Eff. ", "l");
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


  // overlay multiple sample histograms :

  TCanvas* cLog = new TCanvas("cLog", "Overlayed Histograms", 800, 600);
  cLog->SetLogy();

  TLegend* legendLog = new TLegend(0.7, 0.7, 0.95, 0.95);
  legendLog->SetTextSize(0.03); // Adjust text size
  legendLog->AddEntry(displaced_track_pt_1_01_, "1 TeV, 30 mm", "l");
  legendLog->AddEntry(displaced_track_pt_1_1_, "1 TeV, 300 mm", "l");
  legendLog->AddEntry(displaced_track_pt_45_01_, "4.5 TeV, 30 mm", "l");
  legendLog->AddEntry(displaced_track_pt_45_1_, "4.5 TeV, 300 mm", "l");
  legendLog->Draw();

  // Draw histograms on the same canvas
  displaced_track_d0_1_01_->SetLineColor(kRed);
  displaced_track_d0_1_01_->GetXaxis()->SetTitleSize(0.05);
  displaced_track_d0_1_01_->GetYaxis()->SetTitleSize(0.05);
  displaced_track_d0_1_01_->SetLineWidth(2); // Set line width to 3
  displaced_track_d0_1_01_->SetMarkerSize(1.5); // Set marker size to 1.5
  displaced_track_d0_1_01_->SetMaximum(240);
  displaced_track_d0_1_1_->SetLineColor(kBlue);
  displaced_track_d0_1_1_->GetXaxis()->SetTitleSize(0.05);
  displaced_track_d0_1_1_->GetYaxis()->SetTitleSize(0.05);
  displaced_track_d0_1_1_->SetLineWidth(2); // Set line width to 3
  displaced_track_d0_1_1_->SetMarkerSize(1.5); // Set marker size to 1.5
  displaced_track_d0_45_01_->SetLineColor(kGreen);
  displaced_track_d0_45_01_->GetXaxis()->SetTitleSize(0.05);
  displaced_track_d0_45_01_->GetYaxis()->SetTitleSize(0.05);
  displaced_track_d0_45_01_->SetLineWidth(2); // Set line width to 3
  displaced_track_d0_45_01_->SetMarkerSize(1.5); // Set marker size to 1.5
  displaced_track_d0_45_1_->SetLineColor(kMagenta);
  displaced_track_d0_45_1_->GetXaxis()->SetTitleSize(0.05);
  displaced_track_d0_45_1_->GetYaxis()->SetTitleSize(0.05);
  displaced_track_d0_45_1_->SetLineWidth(2); // Set line width to 3
  displaced_track_d0_45_1_->SetMarkerSize(1.5); // Set marker size to 1.5

  // Draw the first histogram and the rest on the same canvas
  displaced_track_d0_1_01_->Draw("HIST");
  displaced_track_d0_1_1_->Draw("HIST SAME");
  displaced_track_d0_45_01_->Draw("HIST SAME");
  displaced_track_d0_45_1_->Draw("HIST SAME");

  // Create and draw legend
  legendLog->Draw();

  // Save canvas
  cLog->SaveAs("displaced_track_d0_overlay.pdf");

  // Create canvas
  TCanvas* c1 = new TCanvas("c1", "Overlayed Histograms", 800, 600);

  // Draw histograms on the same canvas
  displaced_track_pt_1_01_->SetLineColor(kRed);
  displaced_track_pt_1_01_->GetXaxis()->SetTitleSize(0.05);
  displaced_track_pt_1_01_->GetYaxis()->SetTitleSize(0.05);
  displaced_track_pt_1_01_->SetLineWidth(2); // Set line width to 3
  displaced_track_pt_1_01_->SetMarkerSize(1.5); // Set marker size to 1.5
  displaced_track_pt_1_1_->SetLineColor(kBlue);
  displaced_track_pt_1_1_->GetXaxis()->SetTitleSize(0.05);
  displaced_track_pt_1_1_->GetYaxis()->SetTitleSize(0.05);
  displaced_track_pt_1_1_->SetLineWidth(2); // Set line width to 3
  displaced_track_pt_1_1_->SetMarkerSize(1.5); // Set marker size to 1.5
  displaced_track_pt_45_01_->SetLineColor(kGreen);
  displaced_track_pt_45_01_->GetXaxis()->SetTitleSize(0.05);
  displaced_track_pt_45_01_->GetYaxis()->SetTitleSize(0.05);
  displaced_track_pt_45_01_->SetLineWidth(2); // Set line width to 3
  displaced_track_pt_45_01_->SetMarkerSize(1.5); // Set marker size to 1.5
  displaced_track_pt_45_1_->SetLineColor(kMagenta);
  displaced_track_pt_45_1_->GetXaxis()->SetTitleSize(0.05);
  displaced_track_pt_45_1_->GetYaxis()->SetTitleSize(0.05);
  displaced_track_pt_45_1_->SetLineWidth(2); // Set line width to 3
  displaced_track_pt_45_1_->SetMarkerSize(1.5); // Set marker size to 1.5

  // Draw the first histogram and the rest on the same canvas
  displaced_track_pt_1_01_->Draw("HIST");
  displaced_track_pt_1_1_->Draw("HIST SAME");
  displaced_track_pt_45_01_->Draw("HIST SAME");
  displaced_track_pt_45_1_->Draw("HIST SAME");

  // Create and draw legend
  TLegend* legend1 = new TLegend(0.7, 0.7, 0.95, 0.95);
  legend1->SetTextSize(0.03); // Adjust text size
  legend1->AddEntry(displaced_track_pt_1_01_, "1 TeV, 30 mm", "l");
  legend1->AddEntry(displaced_track_pt_1_1_, "1 TeV, 300 mm", "l");
  legend1->AddEntry(displaced_track_pt_45_01_, "4.5 TeV, 30 mm", "l");
  legend1->AddEntry(displaced_track_pt_45_1_, "4.5 TeV, 300 mm", "l");
  legend1->Draw();

  // Save canvas
  c1->SaveAs("displaced_track_pt_overlay.pdf");

    // Draw histograms on the same canvas
  displaced_nhits_1_01_->SetLineColor(kRed);
  displaced_nhits_1_01_->GetXaxis()->SetTitleSize(0.05);
  displaced_nhits_1_01_->GetYaxis()->SetTitleSize(0.05);
  displaced_nhits_1_01_->SetLineWidth(2); // Set line width to 3
  displaced_nhits_1_01_->SetMarkerSize(1.5); // Set marker size to 1.5
  displaced_nhits_1_1_->SetLineColor(kBlue);
  displaced_nhits_1_1_->GetXaxis()->SetTitleSize(0.05);
  displaced_nhits_1_1_->GetYaxis()->SetTitleSize(0.05);
  displaced_nhits_1_1_->SetLineWidth(2); // Set line width to 3
  displaced_nhits_1_1_->SetMarkerSize(1.5); // Set marker size to 1.5
  displaced_nhits_45_01_->SetLineColor(kGreen);
  displaced_nhits_45_01_->GetXaxis()->SetTitleSize(0.05);
  displaced_nhits_45_01_->GetYaxis()->SetTitleSize(0.05);
  displaced_nhits_45_01_->SetLineWidth(2); // Set line width to 3
  displaced_nhits_45_01_->SetMarkerSize(1.5); // Set marker size to 1.5
  displaced_nhits_45_1_->SetLineColor(kMagenta);
  displaced_nhits_45_1_->GetXaxis()->SetTitleSize(0.05);
  displaced_nhits_45_1_->GetYaxis()->SetTitleSize(0.05);
  displaced_nhits_45_1_->SetLineWidth(2); // Set line width to 3
  displaced_nhits_45_1_->SetMarkerSize(1.5); // Set marker size to 1.5

  // Draw the first histogram and the rest on the same canvas
  displaced_nhits_1_01_->Draw("HIST");
  displaced_nhits_1_1_->Draw("HIST SAME");
  displaced_nhits_45_01_->Draw("HIST SAME");
  displaced_nhits_45_1_->Draw("HIST SAME");

  // Create and draw legend
  legend1->Draw();

  // Save canvas
  c1->SaveAs("displaced_track_nhits_overlay.pdf");

  // Draw histograms on the same canvas
  displaced_chi2_reduced_1_01_->SetLineColor(kRed);
  displaced_chi2_reduced_1_01_->GetXaxis()->SetTitleSize(0.05);
  displaced_chi2_reduced_1_01_->GetYaxis()->SetTitleSize(0.05);
  displaced_chi2_reduced_1_01_->SetLineWidth(2); // Set line width to 3
  displaced_chi2_reduced_1_01_->SetMarkerSize(1.5); // Set marker size to 1.5
  displaced_chi2_reduced_1_1_->SetLineColor(kBlue);
  displaced_chi2_reduced_1_1_->GetXaxis()->SetTitleSize(0.05);
  displaced_chi2_reduced_1_1_->GetYaxis()->SetTitleSize(0.05);
  displaced_chi2_reduced_1_1_->SetLineWidth(2); // Set line width to 3
  displaced_chi2_reduced_1_1_->SetMarkerSize(1.5); // Set marker size to 1.5
  displaced_chi2_reduced_45_01_->SetLineColor(kGreen);
  displaced_chi2_reduced_45_01_->GetXaxis()->SetTitleSize(0.05);
  displaced_chi2_reduced_45_01_->GetYaxis()->SetTitleSize(0.05);
  displaced_chi2_reduced_45_01_->SetLineWidth(2); // Set line width to 3
  displaced_chi2_reduced_45_01_->SetMarkerSize(1.5); // Set marker size to 1.5
  displaced_chi2_reduced_45_1_->SetLineColor(kMagenta);
  displaced_chi2_reduced_45_1_->GetXaxis()->SetTitleSize(0.05);
  displaced_chi2_reduced_45_1_->GetYaxis()->SetTitleSize(0.05);
  displaced_chi2_reduced_45_1_->SetLineWidth(2); // Set line width to 3
  displaced_chi2_reduced_45_1_->SetMarkerSize(1.5); // Set marker size to 1.5

  // Draw the first histogram and the rest on the same canvas
  displaced_chi2_reduced_1_01_->Draw("HIST");
  displaced_chi2_reduced_1_1_->Draw("HIST SAME");
  displaced_chi2_reduced_45_01_->Draw("HIST SAME");
  displaced_chi2_reduced_45_1_->Draw("HIST SAME");

  // Create and draw legend
  legend1->Draw();

  // Save canvas
  c1->SaveAs("displaced_track_chi2red_overlay.pdf");

  // efficiency vs. d0
  displaced_matched_d0_1_01_->Sumw2();
  displaced_tp_d0_1_01_->Sumw2();
  TH1F* displaced_eff_acc_d0_1_01 = (TH1F*)displaced_matched_d0_1_01_->Clone();
  displaced_eff_acc_d0_1_01->SetName("eff_eta");
  displaced_eff_acc_d0_1_01->GetYaxis()->SetTitle("Efficiency");
  displaced_eff_acc_d0_1_01->Divide(displaced_matched_d0_1_01_, displaced_tp_d0_1_01_, 1.0, 1.0, "B");
  displaced_eff_acc_d0_1_01->SetAxisRange(0, 1.1, "Y");
  displaced_eff_acc_d0_1_01->SetLineColor(kRed);
  displaced_eff_acc_d0_1_01->GetXaxis()->SetTitleSize(0.05);
  displaced_eff_acc_d0_1_01->GetYaxis()->SetTitleSize(0.05);
  displaced_eff_acc_d0_1_01->SetLineWidth(2); // Set line width to 3
  displaced_eff_acc_d0_1_01->SetMarkerSize(1.5); // Set marker size to 1.5
  displaced_eff_acc_d0_1_01->Draw();

  displaced_matched_d0_1_1_->Sumw2();
  displaced_tp_d0_1_1_->Sumw2();
  TH1F* displaced_eff_acc_d0_1_1 = (TH1F*)displaced_matched_d0_1_1_->Clone();
  displaced_eff_acc_d0_1_1->SetName("eff_eta");
  displaced_eff_acc_d0_1_1->GetYaxis()->SetTitle("Efficiency");
  displaced_eff_acc_d0_1_1->Divide(displaced_matched_d0_1_1_, displaced_tp_d0_1_1_, 1.0, 1.0, "B");
  displaced_eff_acc_d0_1_1->SetAxisRange(0, 1.1, "Y");
  displaced_eff_acc_d0_1_1->SetLineColor(kBlue);
  displaced_eff_acc_d0_1_1->GetXaxis()->SetTitleSize(0.05);
  displaced_eff_acc_d0_1_1->GetYaxis()->SetTitleSize(0.05);
  displaced_eff_acc_d0_1_1->SetLineWidth(2); // Set line width to 3
  displaced_eff_acc_d0_1_1->SetMarkerSize(1.5); // Set marker size to 1.5
  displaced_eff_acc_d0_1_1->Draw("SAME");

  displaced_matched_d0_45_01_->Sumw2();
  displaced_tp_d0_45_01_->Sumw2();
  TH1F* displaced_eff_acc_d0_45_01 = (TH1F*)displaced_matched_d0_45_01_->Clone();
  displaced_eff_acc_d0_45_01->SetName("eff_eta");
  displaced_eff_acc_d0_45_01->GetYaxis()->SetTitle("Efficiency");
  displaced_eff_acc_d0_45_01->Divide(displaced_matched_d0_45_01_, displaced_tp_d0_45_01_, 1.0, 1.0, "B");
  displaced_eff_acc_d0_45_01->SetAxisRange(0, 1.1, "Y");
  displaced_eff_acc_d0_45_01->SetLineColor(kGreen);
  displaced_eff_acc_d0_45_01->GetXaxis()->SetTitleSize(0.05);
  displaced_eff_acc_d0_45_01->GetYaxis()->SetTitleSize(0.05);
  displaced_eff_acc_d0_45_01->SetLineWidth(2); // Set line width to 3
  displaced_eff_acc_d0_45_01->SetMarkerSize(1.5); // Set marker size to 1.5
  //displaced_eff_acc_d0_45_01->Draw("SAME");

  displaced_matched_d0_45_1_->Sumw2();
  displaced_tp_d0_45_1_->Sumw2();
  TH1F* displaced_eff_acc_d0_45_1 = (TH1F*)displaced_matched_d0_45_1_->Clone();
  displaced_eff_acc_d0_45_1->SetName("eff_eta");
  displaced_eff_acc_d0_45_1->GetYaxis()->SetTitle("Efficiency");
  displaced_eff_acc_d0_45_1->Divide(displaced_matched_d0_45_1_, displaced_tp_d0_45_1_, 1.0, 1.0, "B");
  displaced_eff_acc_d0_45_1->SetAxisRange(0, 1.1, "Y");
  displaced_eff_acc_d0_45_1->SetLineColor(kMagenta);
  displaced_eff_acc_d0_45_1->GetXaxis()->SetTitleSize(0.05);
  displaced_eff_acc_d0_45_1->GetYaxis()->SetTitleSize(0.05);
  displaced_eff_acc_d0_45_1->SetLineWidth(2); // Set line width to 3
  displaced_eff_acc_d0_45_1->SetMarkerSize(1.5); // Set marker size to 1.5
  //displaced_eff_acc_d0_45_1->Draw("SAME");
  
  TLegend* legend2 = new TLegend(0.75, 0.75, 0.95, 0.95);
  legend2->SetTextSize(0.0275); // Adjust text size
  legend2->AddEntry(displaced_eff_acc_d0_1_01, "1 TeV, 30 mm", "l");
  legend2->AddEntry(displaced_eff_acc_d0_1_1, "1 TeV, 300 mm", "l");
  legend2->Draw();

  c1->SaveAs("displaced_eff_acc_r_vertex_overlay.pdf");

  displaced_matched_d0_1_01_->Sumw2();
  displaced_mcp_d0_1_01_->Sumw2();
  TH1F* displaced_eff_d0_1_01 = (TH1F*)displaced_matched_d0_1_01_->Clone();
  displaced_eff_d0_1_01->SetName("eff_eta");
  displaced_eff_d0_1_01->GetYaxis()->SetTitle("Efficiency");
  displaced_eff_d0_1_01->Divide(displaced_matched_d0_1_01_, displaced_mcp_d0_1_01_, 1.0, 1.0, "B");
  displaced_eff_d0_1_01->SetAxisRange(0, 1.1, "Y");
  displaced_eff_d0_1_01->SetLineColor(kRed);
  displaced_eff_d0_1_01->GetXaxis()->SetTitleSize(0.05);
  displaced_eff_d0_1_01->GetYaxis()->SetTitleSize(0.05);
  displaced_eff_d0_1_01->SetLineWidth(2); // Set line width to 3
  displaced_eff_d0_1_01->SetMarkerSize(1.5); // Set marker size to 1.5
  displaced_eff_d0_1_01->Draw();

  displaced_matched_d0_1_1_->Sumw2();
  displaced_mcp_d0_1_1_->Sumw2();
  TH1F* displaced_eff_d0_1_1 = (TH1F*)displaced_matched_d0_1_1_->Clone();
  displaced_eff_d0_1_1->SetName("eff_eta");
  displaced_eff_d0_1_1->GetYaxis()->SetTitle("Efficiency");
  displaced_eff_d0_1_1->Divide(displaced_matched_d0_1_1_, displaced_mcp_d0_1_1_, 1.0, 1.0, "B");
  displaced_eff_d0_1_1->SetAxisRange(0, 1.1, "Y");
  displaced_eff_d0_1_1->SetLineColor(kBlue);
  displaced_eff_d0_1_1->GetXaxis()->SetTitleSize(0.05);
  displaced_eff_d0_1_1->GetYaxis()->SetTitleSize(0.05);
  displaced_eff_d0_1_1->SetLineWidth(2); // Set line width to 3
  displaced_eff_d0_1_1->SetMarkerSize(1.5); // Set marker size to 1.5
  displaced_eff_d0_1_1->Draw("SAME");

  displaced_matched_d0_45_01_->Sumw2();
  displaced_mcp_d0_45_01_->Sumw2();
  TH1F* displaced_eff_d0_45_01 = (TH1F*)displaced_matched_d0_45_01_->Clone();
  displaced_eff_d0_45_01->SetName("eff_eta");
  displaced_eff_d0_45_01->GetYaxis()->SetTitle("Efficiency");
  displaced_eff_d0_45_01->Divide(displaced_matched_d0_45_01_, displaced_mcp_d0_45_01_, 1.0, 1.0, "B");
  displaced_eff_d0_45_01->SetAxisRange(0, 1.1, "Y");
  displaced_eff_d0_45_01->SetLineColor(kGreen);
  displaced_eff_d0_45_01->GetXaxis()->SetTitleSize(0.05);
  displaced_eff_d0_45_01->GetYaxis()->SetTitleSize(0.05);
  displaced_eff_d0_45_01->SetLineWidth(2); // Set line width to 3
  displaced_eff_d0_45_01->SetMarkerSize(1.5); // Set marker size to 1.5
  //displaced_eff_d0_45_01->Draw("SAME");

  displaced_matched_d0_45_1_->Sumw2();
  displaced_mcp_d0_45_1_->Sumw2();
  TH1F* displaced_eff_d0_45_1 = (TH1F*)displaced_matched_d0_45_1_->Clone();
  displaced_eff_d0_45_1->SetName("eff_eta");
  displaced_eff_d0_45_1->GetYaxis()->SetTitle("Efficiency");
  displaced_eff_d0_45_1->Divide(displaced_matched_d0_45_1_, displaced_mcp_d0_45_1_, 1.0, 1.0, "B");
  displaced_eff_d0_45_1->SetAxisRange(0, 1.1, "Y");
  displaced_eff_d0_45_1->SetLineColor(kMagenta);
  displaced_eff_d0_45_1->GetXaxis()->SetTitleSize(0.05);
  displaced_eff_d0_45_1->GetYaxis()->SetTitleSize(0.05);
  displaced_eff_d0_45_1->SetLineWidth(2); // Set line width to 3
  displaced_eff_d0_45_1->SetMarkerSize(1.5); // Set marker size to 1.5
  //displaced_eff_d0_45_1->Draw("SAME");
  
  legend2->Draw();

  c1->SaveAs("displaced_eff_r_vertex_overlay.pdf");

  displaced_matched_pt_1_01_->Sumw2();
  displaced_tp_pt_1_01_->Sumw2();
  TH1F* displaced_eff_acc_pt_1_01 = (TH1F*)displaced_matched_pt_1_01_->Clone();
  displaced_eff_acc_pt_1_01->SetName("eff_eta");
  displaced_eff_acc_pt_1_01->GetYaxis()->SetTitle("Efficiency");
  displaced_eff_acc_pt_1_01->Divide(displaced_matched_pt_1_01_, displaced_tp_pt_1_01_, 1.0, 1.0, "B");
  displaced_eff_acc_pt_1_01->SetAxisRange(0, 1.1, "Y");
  displaced_eff_acc_pt_1_01->SetLineColor(kRed);
  displaced_eff_acc_pt_1_01->GetXaxis()->SetTitleSize(0.05);
  displaced_eff_acc_pt_1_01->GetYaxis()->SetTitleSize(0.05);
  displaced_eff_acc_pt_1_01->SetLineWidth(2); // Set line width to 3
  displaced_eff_acc_pt_1_01->SetMarkerSize(1.5); // Set marker size to 1.5
  displaced_eff_acc_pt_1_01->Draw();

  displaced_matched_pt_1_1_->Sumw2();
  displaced_tp_pt_1_1_->Sumw2();
  TH1F* displaced_eff_acc_pt_1_1 = (TH1F*)displaced_matched_pt_1_1_->Clone();
  displaced_eff_acc_pt_1_1->SetName("eff_eta");
  displaced_eff_acc_pt_1_1->GetYaxis()->SetTitle("Efficiency");
  displaced_eff_acc_pt_1_1->Divide(displaced_matched_pt_1_1_, displaced_tp_pt_1_1_, 1.0, 1.0, "B");
  displaced_eff_acc_pt_1_1->SetAxisRange(0, 1.1, "Y");
  displaced_eff_acc_pt_1_1->SetLineColor(kBlue);
  displaced_eff_acc_pt_1_1->GetXaxis()->SetTitleSize(0.05);
  displaced_eff_acc_pt_1_1->GetYaxis()->SetTitleSize(0.05);
  displaced_eff_acc_pt_1_1->SetLineWidth(2); // Set line width to 3
  displaced_eff_acc_pt_1_1->SetMarkerSize(1.5); // Set marker size to 1.5
  displaced_eff_acc_pt_1_1->Draw("SAME");

  displaced_matched_pt_45_01_->Sumw2();
  displaced_tp_pt_45_01_->Sumw2();
  TH1F* displaced_eff_acc_pt_45_01 = (TH1F*)displaced_matched_pt_45_01_->Clone();
  displaced_eff_acc_pt_45_01->SetName("eff_eta");
  displaced_eff_acc_pt_45_01->GetYaxis()->SetTitle("Efficiency");
  displaced_eff_acc_pt_45_01->Divide(displaced_matched_pt_45_01_, displaced_tp_pt_45_01_, 1.0, 1.0, "B");
  displaced_eff_acc_pt_45_01->SetAxisRange(0, 1.1, "Y");
  displaced_eff_acc_pt_45_01->SetLineColor(kGreen);
  displaced_eff_acc_pt_45_01->GetXaxis()->SetTitleSize(0.05);
  displaced_eff_acc_pt_45_01->GetYaxis()->SetTitleSize(0.05);
  displaced_eff_acc_pt_45_01->SetLineWidth(2); // Set line width to 3
  displaced_eff_acc_pt_45_01->SetMarkerSize(1.5); // Set marker size to 1.5
  displaced_eff_acc_pt_45_01->Draw("SAME");

  displaced_matched_pt_45_1_->Sumw2();
  displaced_tp_pt_45_1_->Sumw2();
  TH1F* displaced_eff_acc_pt_45_1 = (TH1F*)displaced_matched_pt_45_1_->Clone();
  displaced_eff_acc_pt_45_1->SetName("eff_eta");
  displaced_eff_acc_pt_45_1->GetYaxis()->SetTitle("Efficiency");
  displaced_eff_acc_pt_45_1->Divide(displaced_matched_pt_45_1_, displaced_tp_pt_45_1_, 1.0, 1.0, "B");
  displaced_eff_acc_pt_45_1->SetAxisRange(0, 1.1, "Y");
  displaced_eff_acc_pt_45_1->SetLineColor(kMagenta);
  displaced_eff_acc_pt_45_1->GetXaxis()->SetTitleSize(0.05);
  displaced_eff_acc_pt_45_1->GetYaxis()->SetTitleSize(0.05);
  displaced_eff_acc_pt_45_1->SetLineWidth(2); // Set line width to 3
  displaced_eff_acc_pt_45_1->SetMarkerSize(1.5); // Set marker size to 1.5
  displaced_eff_acc_pt_45_1->Draw("SAME");
  
  TLegend* legend3 = new TLegend(0.8, 0.8, 0.95, 0.95);
  legend3->SetTextSize(0.0215); // Adjust text size
  legend3->AddEntry(displaced_eff_acc_pt_1_01, "1 TeV, 30 mm", "l");
  legend3->AddEntry(displaced_eff_acc_pt_1_1, "1 TeV, 300 mm", "l");
  legend3->AddEntry(displaced_eff_acc_pt_45_01, "4.5 TeV, 30 mm", "l");
  legend3->AddEntry(displaced_eff_acc_pt_45_1, "4.5 TeV, 300 mm", "l");
  legend3->Draw();

  c1->SaveAs("displaced_eff_acc_pt_overlay.pdf");

  displaced_matched_pt_1_01_->Sumw2();
  displaced_mcp_pt_1_01_->Sumw2();
  TH1F* displaced_eff_pt_1_01 = (TH1F*)displaced_matched_pt_1_01_->Clone();
  displaced_eff_pt_1_01->SetName("eff_eta");
  displaced_eff_pt_1_01->GetYaxis()->SetTitle("Efficiency");
  displaced_eff_pt_1_01->Divide(displaced_matched_pt_1_01_, displaced_mcp_pt_1_01_, 1.0, 1.0, "B");
  displaced_eff_pt_1_01->SetAxisRange(0, 1.1, "Y");
  displaced_eff_pt_1_01->SetLineColor(kRed);
  displaced_eff_pt_1_01->GetXaxis()->SetTitleSize(0.05);
  displaced_eff_pt_1_01->GetYaxis()->SetTitleSize(0.05);
  displaced_eff_pt_1_01->SetLineWidth(2); // Set line width to 3
  displaced_eff_pt_1_01->SetMarkerSize(1.5); // Set marker size to 1.5
  displaced_eff_pt_1_01->Draw();

  displaced_matched_pt_1_1_->Sumw2();
  displaced_mcp_pt_1_1_->Sumw2();
  TH1F* displaced_eff_pt_1_1 = (TH1F*)displaced_matched_pt_1_1_->Clone();
  displaced_eff_pt_1_1->SetName("eff_eta");
  displaced_eff_pt_1_1->GetYaxis()->SetTitle("Efficiency");
  displaced_eff_pt_1_1->Divide(displaced_matched_pt_1_1_, displaced_mcp_pt_1_1_, 1.0, 1.0, "B");
  displaced_eff_pt_1_1->SetAxisRange(0, 1.1, "Y");
  displaced_eff_pt_1_1->SetLineColor(kBlue);
  displaced_eff_pt_1_1->GetXaxis()->SetTitleSize(0.05);
  displaced_eff_pt_1_1->GetYaxis()->SetTitleSize(0.05);
  displaced_eff_pt_1_1->SetLineWidth(2); // Set line width to 3
  displaced_eff_pt_1_1->SetMarkerSize(1.5); // Set marker size to 1.5
  displaced_eff_pt_1_1->Draw("SAME");

  displaced_matched_pt_45_01_->Sumw2();
  displaced_mcp_pt_45_01_->Sumw2();
  TH1F* displaced_eff_pt_45_01 = (TH1F*)displaced_matched_pt_45_01_->Clone();
  displaced_eff_pt_45_01->SetName("eff_eta");
  displaced_eff_pt_45_01->GetYaxis()->SetTitle("Efficiency");
  displaced_eff_pt_45_01->Divide(displaced_matched_pt_45_01_, displaced_mcp_pt_45_01_, 1.0, 1.0, "B");
  displaced_eff_pt_45_01->SetAxisRange(0, 1.1, "Y");
  displaced_eff_pt_45_01->SetLineColor(kGreen);
  displaced_eff_pt_45_01->GetXaxis()->SetTitleSize(0.05);
  displaced_eff_pt_45_01->GetYaxis()->SetTitleSize(0.05);
  displaced_eff_pt_45_01->SetLineWidth(2); // Set line width to 3
  displaced_eff_pt_45_01->SetMarkerSize(1.5); // Set marker size to 1.5
  displaced_eff_pt_45_01->Draw("SAME");

  displaced_matched_pt_45_1_->Sumw2();
  displaced_mcp_pt_45_1_->Sumw2();
  TH1F* displaced_eff_pt_45_1 = (TH1F*)displaced_matched_pt_45_1_->Clone();
  displaced_eff_pt_45_1->SetName("eff_eta");
  displaced_eff_pt_45_1->GetYaxis()->SetTitle("Efficiency");
  displaced_eff_pt_45_1->Divide(displaced_matched_pt_45_1_, displaced_mcp_pt_45_1_, 1.0, 1.0, "B");
  displaced_eff_pt_45_1->SetAxisRange(0, 1.1, "Y");
  displaced_eff_pt_45_1->SetLineColor(kMagenta);
  displaced_eff_pt_45_1->GetXaxis()->SetTitleSize(0.05);
  displaced_eff_pt_45_1->GetYaxis()->SetTitleSize(0.05);
  displaced_eff_pt_45_1->SetLineWidth(2); // Set line width to 3
  displaced_eff_pt_45_1->SetMarkerSize(1.5); // Set marker size to 1.5
  displaced_eff_pt_45_1->Draw("SAME");
  
  legend3->Draw();

  c1->SaveAs("displaced_eff_pt_overlay.pdf");



  std::cout << "---------------------------------------" << "\n";
  std::cout << "--------   SAMPLE SUMMARY     ---------" << "\n";
  std::cout << "---------------------------------------" << "\n";

  std::cout << "numRecoableStaus: " << numRecoableStaus << "\n";
  std::cout << "numMatchedStauTracks: " << numMatchedStauTracks << "\n";
  std::cout << "total stau tracking efficiency: " << float(numMatchedStauTracks) / float(numRecoableStaus) << "\n";

  std::cout << "numRecoableDisplacedParticles: " << numRecoableDisplacedParticles << "\n";
  std::cout << "numMatchedDisplacedTracks: " << numMatchedDisplacedTracks << "\n";
  std::cout << "total displaced tracking efficiency: " << float(numMatchedDisplacedTracks) / float(numRecoableDisplacedParticles) << "\n";
  
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
