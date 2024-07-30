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


void rootAnalyzer(
const TString fileName
){  
    // Example Usage: rootAnalyzer("1000_10_reco_rt")
    // Will output plots as pdfs to a subdirectory named according to the filename

    // ----------------------------------------------------------------------------------------------------------------
    // Set plot style and other configurations
    SetPlotStyle();



    // ----------------------------------------------------------------------------------------------------------------
    // Read NTuples from one file
    TChain* MCPs = new TChain("MCPs");
    MCPs->Add(fileName + ".root");
    TChain* StauMCPs = new TChain("StauMCPs");
    StauMCPs->Add(fileName + ".root");
    TChain* StauDecayProductsMCPs = new TChain("StauDecayProductsMCPs");
    StauDecayProductsMCPs->Add(fileName + ".root");
    TChain* AllTracks = new TChain("AllTracks");
    AllTracks->Add(fileName + ".root");
    TChain* AllStauTracks = new TChain("AllStauTracks");
    AllStauTracks->Add(fileName + ".root");
    TChain* AllDecayProductTracks = new TChain("AllDecayProductTracks");
    AllDecayProductTracks->Add(fileName + ".root");
    TChain* AllFakeTracks = new TChain("AllFakeTracks");
    AllFakeTracks->Add(fileName + ".root");
    TChain* AllHits = new TChain("AllHits");
    AllHits->Add(fileName + ".root");

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
    // Loop over the entries
    for (int i = 0; i < MCPs->GetEntries(); ++i) {
      // Load the entry data into the branch
      MCPs->GetEntry(i);
      // Now mcp_pt should be populated with data
      std::cout << "Entry " << i << " mcp_pt->size(): " << mcp_pt->size() << "\n";
      for (int it = 0; it < (int)mcp_pt->size(); ++it) {
          std::cout << "mcp_pt: " << mcp_pt->at(it) << " at: " << it << "\n";
      }
    }
    

    // ----------------------------------------------------------------------------------------------------------------
    // Declare histograms to be written to and plotted

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