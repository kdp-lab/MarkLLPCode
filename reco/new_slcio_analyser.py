import pyLCIO
import ROOT
import glob
import json
from math import *
import numpy as np
import uproot

# ############## SETUP #############################
# Prevent ROOT from drawing while you're running -- good for slow remote servers
# Instead, save files and view them with an sftp client like Fetch (feel free to ask me for my UTK license)
#ROOT.gROOT.SetBatch()

# Set up some options, constants
max_events = 100 # Set to -1 to run over all events
Bfield = 5 #T, 3.57 for legacy

def check_hard_radiation(mcp, fractional_threshold):
    had_hard_rad = False
    daughters = mcp.getDaughters() 
    for d in daughters:
        if d.getPDG() == 22 or d.getPDG() == 23 or d.getPDG() == 24:        
            if d.getEnergy() > fractional_threshold*mcp.getEnergy():
                had_hard_rad = True   
    return had_hard_rad

def acceptanceCutsOld(mcp):
    r_vertex = sqrt(mcp.getVertex()[0] ** 2 + mcp.getVertex()[1] ** 2)
    r_endpoint = sqrt(mcp.getEndpoint()[0] ** 2 + mcp.getEndpoint()[1] ** 2)
    print("r_vertex, r_endpoint:", r_vertex, r_endpoint)
    z_vertex = mcp.getVertex()[2]
    z_endpoint = mcp.getEndpoint()[2]
    trvl_dist = sqrt(r_endpoint ** 2 + z_vertex ** 2) - sqrt(r_vertex ** 2 + z_vertex ** 2)
    if r_vertex < 552.0 and r_endpoint > 105.0 and abs(z_vertex) < 1000.0: # produced within itk, decay outside of vxd
        reconstructable = True
        
        print("stau is possibly reconstructable, with z_vertex, endpoint of:", z_vertex, z_endpoint)
    else:
        reconstructable = False
        print("stau is not reconstructable (type 1), with z_vertex, endpoint of:", z_vertex, z_endpoint)
    if trvl_dist < 20.0: 
        print("TOO SHORT LIVED TO CROSS TWO LAYERS")
        reconstructable = False
    if r_vertex > 51.0 and travel_dist < 350.0: 
            print("stau is not reconstructable (type 2), with z_vertex, endpoint of:", z_vertex, z_endpoint)
            #reconstructable = False
    return reconstructable
"""
def acceptanceCutsPrompt(mcp):
    r_vertex = sqrt(mcp.getVertex()[0] ** 2 + mcp.getVertex()[1] ** 2)
    r_endpoint = sqrt(mcp.getEndpoint()[0] ** 2 + mcp.getEndpoint()[1] ** 2)
    print("r_vertex, r_endpoint:", r_vertex, r_endpoint)
    z_vertex = mcp.getVertex()[2]
    z_endpoint = mcp.getEndpoint()[2]
    #trvl_dist = sqrt(r_endpoint ** 2 + z_vertex ** 2) - sqrt(r_vertex ** 2 + z_vertex ** 2)
    if r_endpoint > 105.0: # produced within itk, decay outside of vxd
        reconstructable = True
        print("stau is reconstructable, with z_vertex, endpoint of:", z_vertex, z_endpoint)
    else:
        reconstructable = False
        print("stau is not reconstructable (type 1), with z_vertex, endpoint of:", z_vertex, z_endpoint)
    return reconstructable

def acceptanceCutsDisplaced(mcp):
    r_vertex = sqrt(mcp.getVertex()[0] ** 2 + mcp.getVertex()[1] ** 2)
    r_endpoint = sqrt(mcp.getEndpoint()[0] ** 2 + mcp.getEndpoint()[1] ** 2)
    print("r_vertex, r_endpoint:", r_vertex, r_endpoint)
    z_vertex = mcp.getVertex()[2]
    z_endpoint = mcp.getEndpoint()[2]
    #trvl_dist = sqrt(r_endpoint ** 2 + z_vertex ** 2) - sqrt(r_vertex ** 2 + z_vertex ** 2)
    if r_vertex < 552.0 and r_endpoint > 1490.0: # produced within itk, decay outside ot
        reconstructable = True
        print("displaced product is reconstructable, with z_vertex, endpoint of:", z_vertex, z_endpoint)
    else:
        reconstructable = False
        print("stau is not reconstructable (type 1), with z_vertex, endpoint of:", z_vertex, z_endpoint)
    return reconstructable
"""
# Gather input files
# Note: these are using the path convention from the singularity command in the MuCol tutorial (see README)
base_path = "/home/larsonma/MarkLLPCode/reco/"
sampleNames = ["1000_0.05", "1000_0.1", "1000_1", "1000_10"]
#sampleNames = ["1000_0.05"]
# Loop over sample names to construct file paths
for sample in sampleNames:
    # Use glob to get the file paths
    file_pattern = f"{base_path}{sample}_reco.slcio"
    output_root = f"{base_path}{sample}_reco.root"
    print("file_pattern: ", file_pattern)
    fnames = glob.glob(file_pattern)
    print("Found %i files."%len(fnames))


    import ROOT
    # Create a new ROOT file
    file = ROOT.TFile(output_root, "RECREATE")

    # Create a new tree
    mcp_tree = ROOT.TTree("MCPs", "All MCPs")
    mcp_stau_tree = ROOT.TTree("StauMCPs", "All Stau MCPs")
    mcp_stau_decay_products_tree = ROOT.TTree("StauDecayProducts", "All Stau Charged Decay Product MCPs")
    track_tree = ROOT.TTree("MCPs", "All MCPs")
    stau_track_tree = ROOT.TTree("StauTracks", "All Stau Tracks")
    daughter_track_tree = ROOT.TTree("DaughterTracks", "All Charged Stau Decay Product Tracks")
    fake_track_tree = ROOT.TTree("FakeTracks", "All fake tracks, w/o mcp or multiple non-duplicate MCPs")
    hits_tree = ROOT.TTree("Hits", "Info about all tracker hits")

    # Create empty lists for each variable, and empty root variables to fill root file, also create branches
    mcp_pt = [] #mcp = MCParticle (truth)
    mcp_pt_rt = ROOT.std.vector('float')()
    mcp_tree.Branch('mcp_pt', mcp_pt_rt)

    mcp_phi = []
    mcp_phi_rt = ROOT.std.vector('float')()
    mcp_tree.Branch('mcp_phi', mcp_phi_rt)

    mcp_eta = []
    mcp_eta_rt = ROOT.std.vector('float')()
    mcp_tree.Branch('mcp_eta', mcp_eta_rt)

    pdgid = []
    pdgid_rt = ROOT.std.vector('int')()
    mcp_tree.Branch('pdgid', pdgid_rt)

    status = []
    status_rt = ROOT.std.vector('int')()
    mcp_tree.Branch('status', status_rt)

    prod_vertex_r = []
    prod_vertex_r_rt = ROOT.std.vector('float')()
    mcp_tree.Branch('prod_vertex_r', prod_vertex_r_rt)

    prod_vertex_z = []
    prod_vertex_z_rt = ROOT.std.vector('float')()
    mcp_tree.Branch('prod_vertex_z', prod_vertex_z_rt)

    prod_endpoint_r = []
    prod_endpoint_r_rt = ROOT.std.vector('float')()
    mcp_tree.Branch('prod_endpoint_r', prod_endpoint_r_rt)

    prod_endpoint_z = []
    prod_endpoint_z_rt = ROOT.std.vector('float')()
    mcp_tree.Branch('prod_endpoint_z', prod_endpoint_z_rt)

    prod_traveldist = []
    prod_traveldist_rt = ROOT.std.vector('float')()
    mcp_tree.Branch('prod_traveldist', prod_traveldist_rt)

    prod_time = []
    prod_time_rt = ROOT.std.vector('float')()
    mcp_tree.Branch('prod_time', prod_time_rt)

    id = []
    id_rt = ROOT.std.vector('int')()
    mcp_tree.Branch('id', id_rt)


    mcp_stau_pt = [] #mcp_stau = MCParticle stau
    mcp_stau_pt_rt = ROOT.std.vector('float')()
    mcp_stau_tree.Branch('mcp_stau_pt', mcp_stau_pt_rt)

    mcp_stau_phi = []
    mcp_stau_phi_rt = ROOT.std.vector('float')()
    mcp_stau_tree.Branch('mcp_stau_phi', mcp_stau_phi_rt)

    mcp_stau_eta = []
    mcp_stau_eta_rt = ROOT.std.vector('float')()
    mcp_stau_tree.Branch('mcp_stau_eta', mcp_stau_eta_rt)
    
    mcp_stau_d0 = []
    mcp_stau_d0_rt = ROOT.std.vector('float')()
    mcp_stau_tree.Branch('mcp_stau_d0', mcp_stau_d0_rt)

    mcp_stau_z0 = []
    mcp_stau_z0_rt = ROOT.std.vector('float')()
    mcp_stau_tree.Branch('mcp_stau_z0', mcp_stau_z0_rt)

    mcp_stau_track_bool = []
    mcp_stau_track_bool_rt = ROOT.std.vector('int')()
    mcp_stau_tree.Branch('mcp_stau_track_bool', mcp_stau_track_bool_rt)

    mcp_stau_track_reconstructable_bool = []
    mcp_stau_track_reconstructable_bool_rt = ROOT.std.vector('int')()
    mcp_stau_tree.Branch('mcp_stau_track_reconstructable_bool', mcp_stau_track_reconstructable_bool_rt)

    mcp_two_stau_track_bool = [] ### Dont need in root file


    mcp_daughter_pt = [] #daughter_mcp = charged decay products of tau from stau
    mcp_daughter_pt_rt = ROOT.std.vector('float')()
    mcp_stau_decay_products_tree.Branch('mcp_daughter_pt_rt', mcp_daughter_pt_rt)

    mcp_daughter_phi = []
    mcp_daughter_phi_rt = ROOT.std.vector('float')()
    mcp_stau_decay_products_tree.Branch('mcp_daughter_phi', mcp_daughter_phi_rt)

    mcp_daughter_eta = []
    mcp_daughter_eta_rt = ROOT.std.vector('float')()
    mcp_stau_decay_products_tree.Branch('mcp_daughter_eta', mcp_daughter_eta_rt)

    mcp_daughter_d0 = []
    mcp_daughter_d0_rt = ROOT.std.vector('float')()
    mcp_stau_decay_products_tree.Branch('mcp_daughter_d0', mcp_daughter_d0_rt)

    mcp_daughter_z0 = []
    mcp_daughter_z0_rt = ROOT.std.vector('float')()
    mcp_stau_decay_products_tree.Branch('mcp_daughter_z0', mcp_daughter_z0_rt)

    mcp_daughter_track_bool = []
    mcp_daughter_track_bool_rt = ROOT.std.vector('int')()
    mcp_stau_decay_products_tree.Branch('mcp_daughter_track_bool', mcp_daughter_track_bool_rt)

    mcp_daughter_track_reconstructable_bool = []
    mcp_daughter_track_reconstructable_bool_rt = ROOT.std.vector('int')()
    mcp_stau_decay_products_tree.Branch('mcp_daughter_track_reconstructable_bool', mcp_daughter_track_reconstructable_bool_rt)

    # TRACK - (FOR ALL TRACKS in TRACKCOLLECTION) ### FIXME add resolutions for real tracks (associated to 1 mcp) and matching (but already done for decay product and stau tracks)
    nhits = []
    nhits_rt = ROOT.std.vector('int')()
    track_tree.Branch('nhits', nhits_rt)

    pixel_nhits = []
    pixel_nhits_rt = ROOT.std.vector('int')()
    track_tree.Branch('pixel_nhits', pixel_nhits_rt)

    inner_nhits = []
    inner_nhits_rt = ROOT.std.vector('int')()
    track_tree.Branch('inner_nhits', inner_nhits_rt)

    outer_nhits = [] 
    outer_nhits = ROOT.std.vector('int')()
    track_tree.Branch('outer_nhits', outer_nhits_rt)

    track_pt = [] #This is track pt
    track_pt_rt = ROOT.std.vector('float')()
    track_tree.Branch('track_pt', track_pt_rt)

    track_eta = [] #This is track eta
    track_eta_rt = ROOT.std.vector('float')()
    track_tree.Branch('track_eta', track_eta_rt)

    track_theta = []
    track_theta_rt = ROOT.std.vector('float')()
    track_tree.Branch('track_theta', track_theta_rt)

    ndf = []
    ndf_rt = ROOT.std.vector('int')()
    track_tree.Branch('ndf', ndf_rt)

    chi2 = []
    chi2_rt = ROOT.std.vector('float')()
    track_tree.Branch('chi2', chi2_rt)

    chi2_red = []
    chi2_red_rt = ROOT.std.vector('float')()
    track_tree.Branch('chi2_red', chi2_red_rt)


    ### FIRST FOR STAUS
    LC_stau_pt_match = []
    LC_stau_pt_match_rt = ROOT.std.vector('float')()
    stau_track_tree.Branch('LC_stau_pt_match', LC_stau_pt_match_rt)

    LC_stau_track_pt = []
    LC_stau_track_pt_rt = ROOT.std.vector('float')()
    stau_track_tree.Branch('LC_stau_track_pt', LC_stau_track_pt_rt)

    LC_stau_track_eta = []
    LC_stau_track_eta_rt = ROOT.std.vector('float')()
    stau_track_tree.Branch('LC_stau_track_eta', LC_stau_track_eta_rt)

    LC_stau_eta_match = []
    LC_stau_eta_match_rt = ROOT.std.vector('float')()
    stau_track_tree.Branch('LC_stau_eta_match', LC_stau_eta_match_rt)

    LC_stau_track_theta = []
    LC_stau_track_theta_rt = ROOT.std.vector('float')()
    stau_track_tree.Branch('LC_stau_track_theta', LC_stau_track_theta_rt)

    LC_stau_phi_match = []
    LC_stau_phi_match_rt = ROOT.std.vector('float')()
    stau_track_tree.Branch('LC_stau_phi_match', LC_stau_phi_match_rt)
    
    LC_stau_ndf = []
    LC_stau_ndf_rt = ROOT.std.vector('int')()
    stau_track_tree.Branch('LC_stau_ndf', LC_stau_ndf_rt)

    LC_stau_chi2 = []
    LC_stau_chi2_rt = ROOT.std.vector('float')()
    stau_track_tree.Branch('LC_stau_chi2', LC_stau_chi2_rt)
    
    LC_stau_d0 = []
    LC_stau_d0_rt = ROOT.std.vector('float')()
    stau_track_tree.Branch('LC_stau_d0', LC_stau_d0_rt)

    LC_stau_z0 = []
    LC_stau_z0_rt = ROOT.std.vector('float')()
    stau_track_tree.Branch('LC_stau_z0', LC_stau_z0_rt)

    LC_stau_nhits = []
    LC_stau_nhits_rt = ROOT.std.vector('int')()
    stau_track_tree.Branch('LC_stau_nhits', LC_stau_nhits_rt)
    
    LC_stau_pixel_nhits = []
    LC_stau_pixel_nhits_rt = ROOT.std.vector('int')()
    stau_track_tree.Branch('LC_stau_pixel_nhits', LC_stau_pixel_nhits_rt)

    LC_stau_inner_nhits = []
    LC_stau_inner_nhits_rt = ROOT.std.vector('int')()
    stau_track_tree.Branch('LC_stau_inner_nhits', LC_stau_inner_nhits_rt)

    LC_stau_outer_nhits = []
    LC_stau_outer_nhits_rt = ROOT.std.vector('int')()
    stau_track_tree.Branch('LC_stau_outer_nhits', LC_stau_outer_nhits_rt)

    LC_stau_pt_res = [] 
    LC_stau_pt_res_rt = ROOT.std.vector('float')()
    stau_track_tree.Branch('LC_stau_pt_res', LC_stau_pt_res_rt)

    LC_stau_dr = []
    LC_stau_dr_rt = ROOT.std.vector('float')()
    stau_track_tree.Branch('LC_stau_dr', LC_stau_dr_rt)
    
    LC_stau_hit_r = []
    LC_stau_hit_r_rt = ROOT.std.vector('float')()
    stau_track_tree.Branch('LC_stau_hit_r', LC_stau_hit_r_rt)

    LC_stau_hit_x = []
    LC_stau_hit_x_rt = ROOT.std.vector('float')()
    stau_track_tree.Branch('LC_stau_hit_x', LC_stau_hit_x_rt)

    LC_stau_hit_y = []
    LC_stau_hit_y_rt = ROOT.std.vector('float')()
    stau_track_tree.Branch('LC_stau_hit_y', LC_stau_hit_y_rt)

    LC_stau_hit_z = []
    LC_stau_hit_z_rt = ROOT.std.vector('float')()
    stau_track_tree.Branch('LC_stau_hit_z', LC_stau_hit_z_rt)

    ### NEXT FOR DAUGHTERS OF TAUS (FROM STAUS) (antcipate tracks to be displaced)
    LC_daughter_pt_match = []
    LC_daughter_pt_match_rt = ROOT.std.vector('float')()
    daughter_track_tree.Branch('LC_daughter_pt_match', LC_daughter_pt_match_rt)
    
    LC_daughter_track_pt = []
    LC_daughter_track_pt_rt = ROOT.std.vector('float')()
    daughter_track_tree.Branch('LC_daughter_track_pt', LC_daughter_track_pt_rt)
    
    LC_daughter_track_eta = []
    LC_daughter_track_eta_rt = ROOT.std.vector('float')()
    daughter_track_tree.Branch('LC_daughter_track_eta', LC_daughter_track_eta_rt)

    LC_daughter_eta_match = []
    LC_daughter_eta_match_rt = ROOT.std.vector('float')()
    daughter_track_tree.Branch('LC_daughter_eta_match', LC_daughter_eta_match_rt)

    LC_daughter_track_theta = []
    LC_daughter_track_theta_rt = ROOT.std.vector('float')()
    daughter_track_tree.Branch('LC_daughter_track_theta', LC_daughter_track_theta_rt)

    LC_daughter_phi_match = []
    LC_daughter_phi_match_rt = ROOT.std.vector('float')()
    daughter_track_tree.Branch('LC_daughter_phi_match', LC_daughter_phi_match_rt)

    LC_daughter_ndf = []
    LC_daughter_ndf_rt = ROOT.std.vector('int')()
    daughter_track_tree.Branch('LC_daughter_ndf', LC_daughter_ndf_rt)

    LC_daughter_chi2 = []
    LC_daughter_chi2_rt = ROOT.std.vector('float')()
    daughter_track_tree.Branch('LC_daughter_chi2', LC_daughter_chi2_rt)

    LC_daughter_d0 = []
    LC_daughter_d0_rt = ROOT.std.vector('float')()
    daughter_track_tree.Branch('LC_daughter_d0', LC_daughter_d0_rt)

    LC_daughter_z0 = []
    LC_daughter_z0_rt = ROOT.std.vector('float')()
    daughter_track_tree.Branch('LC_daughter_z0', LC_daughter_z0_rt)

    LC_daughter_nhits = []
    LC_daughter_nhits_rt = ROOT.std.vector('int')()
    daughter_track_tree.Branch('LC_daughter_nhits', LC_daughter_nhits_rt)

    LC_daughter_pixel_nhits = []
    LC_daughter_pixel_nhits_rt = ROOT.std.vector('int')()
    daughter_track_tree.Branch('LC_daughter_pixel_nhits',LC_daughter_pixel_nhits_rt)

    LC_daughter_inner_nhits = []
    LC_daughter_inner_nhits_rt = ROOT.std.vector('int')()
    daughter_track_tree.Branch('LC_daughter_inner_nhits', LC_daughter_inner_nhits_rt)

    LC_daughter_outer_nhits = []
    LC_daughter_outer_nhits_rt = ROOT.std.vector('int')()
    daughter_track_tree.Branch('LC_daughter_outer_nhits', LC_daughter_outer_nhits_rt)

    LC_daughter_pt_res = [] 
    LC_daughter_pt_res_rt = ROOT.std.vector('float')()
    daughter_track_tree.Branch('LC_daughter_pt_res', LC_daughter_pt_res_rt)

    LC_daughter_dr = []
    LC_daughter_dr_rt = ROOT.std.vector('float')()
    daughter_track_tree.Branch('LC_daughter_dr', LC_daughter_dr_rt)

    ### NEXT FOR ALL TRACKS, indiscriminate of type of particle ### FIXME perform association from track to mcp looping only through track collection (so as to find faketracks too)
    """LC_pt_match = []
    LC_track_pt = []
    LC_track_eta = []
    LC_eta_match = []
    LC_track_theta = []
    LC_phi_match = []
    LC_ndf = []
    LC_chi2 = []
    LC_d0 = []
    LC_z0 = []
    LC_nhits = []
    LC_pixel_nhits = []
    LC_inner_nhits = []
    LC_outer_nhits = []
    LC_pt_res = [] 
    LC_dr = []"""

    # Fake Tracks FIXME how to add faketracks?, loop through remaining tracks without any MC association? 
    fake_theta = []
    fake_theta_rt = ROOT.std.vector('float')()
    fake_track_tree.Branch('fake_theta', fake_theta_rt)
    
    fake_eta = []
    fake_eta_rt = ROOT.std.vector('float')()
    fake_track_tree.Branch('fake_eta', fake_eta_rt)

    fake_pt = []
    fake_pt_rt = ROOT.std.vector('float')()
    fake_track_tree.Branch('fake_pt', fake_pt_rt)

    fake_phi = []
    fake_phi_rt = ROOT.std.vector('float')()
    fake_track_tree.Branch('fake_phi', fake_phi_rt)

    fake_d0 = []
    fake_d0_rt = ROOT.std.vector('float')()
    fake_track_tree.Branch('fake_d0', fake_d0_rt)

    fake_z0 = []
    fake_z0_rt = ROOT.std.vector('float')()
    fake_track_tree.Branch('fake_z0', fake_z0_rt)

    fake_ndf = []
    fake_ndf_rt = ROOT.std.vector('int')()
    fake_track_tree.Branch('fake_ndf', fake_ndf_rt)

    fake_chi2 = []
    fake_chi2_rt = ROOT.std.vector('float')()
    fake_track_tree.Branch('fake_chi2', fake_chi2_rt)

    fake_chi2_reduced = []
    fake_chi2_reduced_rt = ROOT.std.vector('float')()
    fake_track_tree.Branch('fake_chi2_reduced', fake_chi2_reduced_rt)

    fake_nhits = []
    fake_nhits_rt = ROOT.std.vector('int')()
    fake_track_tree.Branch('fake_nhits', fake_nhits_rt)

    fake_pixel_nhits = []
    fake_pixel_nhits_rt = ROOT.std.vector('int')()
    fake_track_tree.Branch('fake_pixel_nhits', fake_pixel_nhits_rt)

    fake_inner_nhits = []
    fake_inner_nhits_rt = ROOT.std.vector('int')()
    fake_track_tree.Branch('fake_inner_nhits', fake_inner_nhits_rt)

    fake_outer_nhits = []
    fake_outer_nhits_rt = ROOT.std.vector('int')()
    fake_track_tree.Branch('fake_outer_nhits', fake_outer_nhits_rt)

    # HITS
    x = []
    x_rt = ROOT.std.vector('float')()
    hits_tree.Branch('x', x_rt)
    
    y = []
    y_rt = ROOT.std.vector('float')()
    hits_tree.Branch('y', y_rt)
    
    z = []
    z_rt = ROOT.std.vector('float')()
    hits_tree.Branch('z', z_rt)

    hit_pdgid = []
    hit_pdgid_rt = ROOT.std.vector('int')() ### FIXME does not work atm
    hits_tree.Branch('hit_pdgid', hit_pdgid_rt)

    time = []
    time_rt = ROOT.std.vector('float')()
    hits_tree.Branch('time', time_rt)

    corrected_time = []
    corrected_time_rt = ROOT.std.vector('float')()
    hits_tree.Branch('corrected_time', corrected_time_rt)

    hit_layer = []
    hit_layer_rt = ROOT.std.vector('int')()
    hits_tree.Branch('hit_layer', hit_layer_rt)

    hit_detector = []
    hit_detector_rt = ROOT.std.vector('int')()
    hits_tree.Branch('hit_detector', hit_detector_rt)

    hit_side = []
    hit_side_rt = ROOT.std.vector('int')()
    hits_tree.Branch('hit_side', hit_side_rt)

    sim_VB_x, sim_VB_y, sim_VB_z, sim_VB_time, sim_VB_pdg, sim_VB_mcpid, sim_VB_layer = [], [], [], [], [], [], []
    sim_VE_x, sim_VE_y, sim_VE_z, sim_VE_time, sim_VE_pdg, sim_VE_mcpid, sim_VE_layer = [], [], [], [], [], [], []
    sim_IB_x, sim_IB_y, sim_IB_z, sim_IB_time, sim_IB_pdg, sim_IB_mcpid, sim_IB_layer = [], [], [], [], [], [], []
    sim_IE_x, sim_IE_y, sim_IE_z, sim_IE_time, sim_IE_pdg, sim_IE_mcpid, sim_IE_layer = [], [], [], [], [], [], []
    sim_OB_x, sim_OB_y, sim_OB_z, sim_OB_time, sim_OB_pdg, sim_OB_mcpid, sim_OB_layer = [], [], [], [], [], [], []
    sim_OE_x, sim_OE_y, sim_OE_z, sim_OE_time, sim_OE_pdg, sim_OE_mcpid, sim_OE_layer = [], [], [], [], [], [], []

    reco_VB_x, reco_VB_y, reco_VB_z, reco_VB_time, reco_VB_layer = [], [], [], [], []
    reco_VE_x, reco_VE_y, reco_VE_z, reco_VE_time, reco_VE_layer = [], [], [], [], []
    reco_IB_x, reco_IB_y, reco_IB_z, reco_IB_time, reco_IB_layer = [], [], [], [], []
    reco_IE_x, reco_IE_y, reco_IE_z, reco_IE_time, reco_IE_layer = [], [], [], [], []
    reco_OB_x, reco_OB_y, reco_OB_z, reco_OB_time, reco_OB_layer = [], [], [], [], []
    reco_OE_x, reco_OE_y, reco_OE_z, reco_OE_time, reco_OE_layer = [], [], [], [], []

    speedoflight = 299792458/1000000  # mm/ns

    # Make counter variables
    n_mcp_stau = 0
    n_charged_mcp_daughter = 0 ### daughters of staus 
    n_recoable_daughter = 0 ### recoable stau decay products
    no_inner_hits = 0
    i = 0
    num_matched_tracks = 0 # total matched tracks
    num_matched_stau_tracks = 0
    num_reconstructable_stau_tracks = 0

    num_matched_daughter_tracks = 0
    num_dupes = 0 ### FIXME how to get # dupe tracks, perhaps if one mcp has two tracks associated to it? 
    num_fake_tracks = 0 ### FIXME how to get # fake, loop through all track collections and see if a track doesn't have mcp associated?? 
    hard_rad_discard = 0
    total_n_pfo_mu = 0
    reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
    reader.setReadCollectionNames(["MCParticle", "PandoraPFOs", "SiTracks", "SiTracks_Refitted", "MCParticle_SiTracks", "MCParticle_SiTracks_Refitted", "IBTrackerHits", "IETrackerHits", "OBTrackerHits", "OETrackerHits", "VBTrackerHits", "VETrackerHits"])
    # ############## LOOP OVER EVENTS AND FILL HISTOGRAMS  #############################
    # Loop over events
    for f in fnames:
        print("---------------------")
        print("processing: " + f)
        print("---------------------")
        if max_events > 0 and i >= max_events: break
        reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
        reader.open(f)
        for ievt,event in enumerate(reader): 
            if max_events > 0 and i >= max_events: break
            #if i%10 == 0: 
            print("Processing event %i."%i)

            # Print all the collection names in the event
            collection_names = event.getCollectionNames()
            if (ievt == 0):
                print("COLLECTION NAMES: ")
                for name in collection_names:
                    print(name)
            

            #has_pfo_stau = False
            
            # Get the collections we care about
            relationCollection = event.getCollection('MCParticle_SiTracks_Refitted')
            lcRelation = pyLCIO.UTIL.LCRelationNavigator(relationCollection)

            # Get the collections we care about
            mcpCollection = event.getCollection("MCParticle")
            trackCollection = event.getCollection("SiTracks_Refitted") ### NOTE use SiTracks_Refitted

            hit_collections = []
            IBTrackerHits = event.getCollection('ITBarrelHits')
            hit_collections.append(IBTrackerHits)
            IETrackerHits = event.getCollection('ITEndcapHits')
            hit_collections.append(IETrackerHits)
            OBTrackerHits = event.getCollection('OTBarrelHits')
            hit_collections.append(OBTrackerHits)
            OETrackerHits = event.getCollection('OTEndcapHits')
            hit_collections.append(OETrackerHits)
            VBTrackerHits = event.getCollection('VXDBarrelHits')
            hit_collections.append(VBTrackerHits)
            VETrackerHits = event.getCollection('VXDEndcapHits')
            hit_collections.append(VETrackerHits)

            #hitRelation = pyLCIO.UTIL.LCRelationNavigator(hit_collections)
            
            # Relations
            VBrelationCollection = event.getCollection('VXDBarrelHitsRelations')
            VBrelation = pyLCIO.UTIL.LCRelationNavigator(VBrelationCollection)
            
            VErelationCollection = event.getCollection('VXDEndcapHitsRelations')
            VErelation = pyLCIO.UTIL.LCRelationNavigator(VErelationCollection)
            
            IBrelationCollection = event.getCollection('ITBarrelHitsRelations')
            IBrelation = pyLCIO.UTIL.LCRelationNavigator(IBrelationCollection)
            
            IErelationCollection = event.getCollection('ITEndcapHitsRelations')
            IErelation = pyLCIO.UTIL.LCRelationNavigator(IErelationCollection)
            
            OBrelationCollection = event.getCollection('OTBarrelHitsRelations')
            OBrelation = pyLCIO.UTIL.LCRelationNavigator(OBrelationCollection)
            
            OErelationCollection = event.getCollection('OTEndcapHitsRelations')
            OErelation = pyLCIO.UTIL.LCRelationNavigator(OErelationCollection)

            # MCPs
            imcp_pt = []
            imcp_eta = []
            imcp_phi = []
            ipdgid = []
            istatus = []
            iprod_vertex_r = []
            iprod_vertex_z = []
            iprod_endpoint_r = []
            iprod_endpoint_z = []
            iprod_traveldist = []
            iprod_time = []
            iid = []
            trackAssociatedStau = []

            # Tracks
            itrack_pt = []
            itrack_eta = []
            itrack_theta = []
            indf = []
            ichi2 = []
            inhits = []

            # HITS
            ix = []
            iy = []
            iz = []
            ihit_pdg = []
            itime = []
            icorrected_time = []
            ihit_layer = []
            ihit_detector = []
            ihit_side = []

            # Sim Hits
            isim_VB_x, isim_VB_y, isim_VB_z, isim_VB_time, isim_VB_pdg, isim_VB_mcpid, isim_VB_layer = [], [], [], [], [], [], []
            isim_VE_x, isim_VE_y, isim_VE_z, isim_VE_time, isim_VE_pdg, isim_VE_mcpid, isim_VE_layer = [], [], [], [], [], [], []
            isim_IB_x, isim_IB_y, isim_IB_z, isim_IB_time, isim_IB_pdg, isim_IB_mcpid, isim_IB_layer = [], [], [], [], [], [], []
            isim_IE_x, isim_IE_y, isim_IE_z, isim_IE_time, isim_IE_pdg, isim_IE_mcpid, isim_IE_layer = [], [], [], [], [], [], []
            isim_OB_x, isim_OB_y, isim_OB_z, isim_OB_time, isim_OB_pdg, isim_OB_mcpid, isim_OB_layer = [], [], [], [], [], [], []
            isim_OE_x, isim_OE_y, isim_OE_z, isim_OE_time, isim_OE_pdg, isim_OE_mcpid, isim_OE_layer = [], [], [], [], [], [], []

            # Reco. Hits

            ireco_VB_x, ireco_VB_y, ireco_VB_z, ireco_VB_time, ireco_VB_layer = [], [], [], [], []
            ireco_VE_x, ireco_VE_y, ireco_VE_z, ireco_VE_time, ireco_VE_layer = [], [], [], [], []
            ireco_IB_x, ireco_IB_y, ireco_IB_z, ireco_IB_time, ireco_IB_layer = [], [], [], [], []
            ireco_IE_x, ireco_IE_y, ireco_IE_z, ireco_IE_time, ireco_IE_layer = [], [], [], [], []
            ireco_OB_x, ireco_OB_y, ireco_OB_z, ireco_OB_time, ireco_OB_layer = [], [], [], [], []
            ireco_OE_x, ireco_OE_y, ireco_OE_z, ireco_OE_time, ireco_OE_layer = [], [], [], [], []

            for coll_name, relation, (ix, iy, iz, itime, ipdg, imcpid, ilayer), (irel_x, irel_y, irel_z, irel_time, irel_layer) in [
                ("VertexBarrelCollection", VBrelation, (isim_VB_x, isim_VB_y, isim_VB_z, isim_VB_time, isim_VB_pdg, isim_VB_mcpid, isim_VB_layer), (ireco_VB_x, ireco_VB_y, ireco_VB_z, ireco_VB_time, ireco_VB_layer)),
                ("VertexEndcapCollection", VErelation, (isim_VE_x, isim_VE_y, isim_VE_z, isim_VE_time, isim_VE_pdg, isim_VE_mcpid, isim_VE_layer), (ireco_VE_x, ireco_VE_y, ireco_VE_z, ireco_VE_time, ireco_VE_layer)),
                ("InnerTrackerBarrelCollection", IBrelation, (isim_IB_x, isim_IB_y, isim_IB_z, isim_IB_time, isim_IB_pdg, isim_IB_mcpid, isim_IB_layer), (ireco_IB_x, ireco_IB_y, ireco_IB_z, ireco_IB_time, ireco_IB_layer)),
                ("InnerTrackerEndcapCollection", IErelation, (isim_IE_x, isim_IE_y, isim_IE_z, isim_IE_time, isim_IE_pdg, isim_IE_mcpid, isim_IE_layer), (ireco_IE_x, ireco_IE_y, ireco_IE_z, ireco_IE_time, ireco_IE_layer)),
                ("OuterTrackerBarrelCollection", OBrelation, (isim_OB_x, isim_OB_y, isim_OB_z, isim_OB_time, isim_OB_pdg, isim_OB_mcpid, isim_OB_layer), (ireco_OB_x, ireco_OB_y, ireco_OB_z, ireco_OB_time, ireco_OB_layer)),
                ("OuterTrackerEndcapCollection", OErelation, (isim_OE_x, isim_OE_y, isim_OE_z, isim_OE_time, isim_OE_pdg, isim_OE_mcpid, isim_OE_layer), (ireco_OE_x, ireco_OE_y, ireco_OE_z, ireco_OE_time, ireco_OE_layer))
            ]:
                try:
                    collection = event.getCollection(coll_name)
                    for hit in collection:
                        # print("hit",collection.getTypeName(), collection.getParameters())
                        position = hit.getPosition()
                        ix.append(position[0])
                        iy.append(position[1])
                        iz.append(position[2])
                        itime.append(hit.getTime())

                        # Retrieve the MCParticle and its PDG code
                        mcp = hit.getMCParticle()
                        hit_pdg = mcp.getPDG() if mcp else None
                        mcpid = mcp.id() if mcp else None
                        ipdg.append(hit_pdg)
                        imcpid.append(mcpid)
                        
                        # Get layer
                        encoding = collection.getParameters().getStringVal(pyLCIO.EVENT.LCIO.CellIDEncoding)
                        decoder = pyLCIO.UTIL.BitField64(encoding)
                        cellID = int(hit.getCellID0())
                        decoder.setValue(cellID)
                        detector = decoder["system"].value()
                        layer = decoder['layer'].value()
                        side = decoder["side"].value()
                        
                        ilayer.append(layer)
                        
                        recohit = relation.getRelatedFromObjects(hit)
                        if len(recohit) > 0:
                            recohit = recohit[0]
                            rel_position = recohit.getPosition()
                            irel_x.append(rel_position[0])
                            irel_y.append(rel_position[1])
                            irel_z.append(rel_position[2])
                            irel_time.append(recohit.getTime())
                            irel_layer.append(layer)
                            
                except Exception as e:
                    print(f"Error accessing {coll_name}: {e}")   # Storing this to use for matching in the next loop

            ### perform truth association 

            imcp_stau_pt = []
            imcp_stau_eta = []
            imcp_stau_phi = []
            imcp_stau_d0 = []
            imcp_stau_z0 = []
            imcp_stau_track_bool = []
            imcp_stau_track_reconstructable_bool = []
            imcp_two_stau_track_bool = []

            imcp_daughter_pt = []
            imcp_daughter_eta = []
            imcp_daughter_phi = []
            imcp_daughter_d0 = []
            imcp_daughter_z0 = []
            imcp_daughter_track_bool = []
            imcp_daughter_track_reconstructable_bool = []

            # Loop over the truth objects and fill histograms
            print("len si tracks refitted: ", len(trackCollection))
            print("len mcp collec: ", len(mcpCollection))
            for mcp in mcpCollection:

                mcp_p = mcp.getMomentum()
                mcp_tlv = ROOT.TLorentzVector()
                mcp_tlv.SetPxPyPzE(mcp_p[0], mcp_p[1], mcp_p[2], mcp.getEnergy())
                pdg = mcp.getPDG()

                imcp_pt.append(mcp_tlv.Perp())
                mcp_pt_rt.clear()
                mcp_pt_rt.push_back(mcp_tlv.Perp())
                imcp_eta.append(mcp_tlv.Eta())
                mcp_eta_rt.clear()
                mcp_eta_rt.push_back(mcp_tlv.Eta())
                imcp_phi.append(mcp_tlv.Phi())
                mcp_phi_rt.clear()
                mcp_phi_rt.push_back(mcp_tlv.Phi())
                ipdgid.append(pdg)
                pdgid_rt.clear()
                pdgid_rt.push_back(pdg)
                istatus.append(mcp.getGeneratorStatus())
                status_rt.clear()
                status_rt.push_back(mcp.getGeneratorStatus())
                iprod_vertex_r.append(sqrt(mcp.getVertex()[0] ** 2 + mcp.getVertex()[1]**2))
                prod_vertex_r_rt.clear()
                prod_vertex_r_rt.push_back(sqrt(mcp.getVertex()[0] ** 2 + mcp.getVertex()[1]**2))
                iprod_vertex_z.append(mcp.getVertex()[2])
                prod_vertex_z_rt.clear()
                prod_vertex_z_rt.push_back(mcp.getVertex()[2])
                iprod_endpoint_r.append(sqrt(mcp.getEndpoint()[0] ** 2 + mcp.getEndpoint()[1]**2))
                prod_endpoint_r_rt.clear()
                prod_endpoint_r_rt.push_back(sqrt(mcp.getEndpoint()[0] ** 2 + mcp.getEndpoint()[1]**2))
                iprod_endpoint_z.append(mcp.getEndpoint()[2])
                prod_endpoint_z_rt.clear()
                prod_endpoint_z_rt.push_back(mcp.getEndpoint()[2])
                travel_dist = sqrt(mcp.getEndpoint()[0]**2 + mcp.getEndpoint()[1]**2 + mcp.getEndpoint()[2]**2) - sqrt(mcp.getVertex()[0]**2 + mcp.getVertex()[1]**2 + mcp.getVertex()[2]**2)
                iprod_traveldist.append(travel_dist)
                prod_traveldist_rt.clear()
                prod_traveldist_rt.push_back(travel_dist)
                print("--------------------------------")
                print("vertex x, y, z: ", mcp.getVertex()[0], mcp.getVertex()[1], mcp.getVertex()[2])
                print("endpoint x, y, z: ", mcp.getEndpoint()[0], mcp.getEndpoint()[1], mcp.getEndpoint()[2])
                
                print("travel_dist: ", travel_dist)
                iprod_time.append(mcp.getTime())
                prod_time_rt.clear()
                prod_time_rt.push_back(mcp.getTime())
                iid.append(mcp.id())
                id_rt.clear()
                id_rt.push_back(mcp.id())
                mcp_tree.Fill()
                print("PID, Status, pT, eta, phi: ", pdg, mcp.getGeneratorStatus(), mcp_tlv.Perp(), mcp_tlv.Eta(), mcp_tlv.Phi())
                print("--------------------------------")
                #print("mcp pdgid: ", mcp.getPDG())
                if abs(mcp.getPDG())==1000015 or abs(mcp.getPDG())==2000015: #### STAUS
                    #print("length of related to for mcp: ", len(lcRelation.getRelatedToObjects(mcp)))
                    #print("length of related from for mcp: ", len(lcRelation.getRelatedToObjects(mcp)))
                    if (travel_dist == 0):
                        #print("travel 0")
                        continue ### skip staus with decay length 0 

                    imcp_stau_pt.append(mcp_tlv.Perp()) 
                    mcp_stau_pt_rt.clear()
                    mcp_stau_pt_rt.push_back(mcp_tlv.Perp())
                    imcp_stau_eta.append(mcp_tlv.Eta())
                    mcp_stau_eta_rt.clear()
                    mcp_stau_eta_rt.push_back(mcp_tlv.Perp())
                    imcp_stau_phi.append(mcp_tlv.Phi())
                    mcp_stau_phi_rt.clear()
                    mcp_stau_phi_rt.push_back(mcp_tlv.Perp())
                    n_mcp_stau += 1

                    # Get the vertex position (true vertex from mcp)
                    # Assuming mcp.getVertex() returns a list or tuple with three elements
                    stau_vx, stau_vy, stau_vz = mcp.getVertex()[0], mcp.getVertex()[1], mcp.getVertex()[2]
                    
                    # Get the momentum
                    stau_px, stau_py, stau_pz = mcp.getMomentum()[0], mcp.getMomentum()[1], mcp.getMomentum()[2]
                    
                    # Calculate transverse impact parameter (d0)
                    stau_pt = sqrt(stau_px**2 + stau_py**2)  # Transverse momentum
                    stau_d0 = (stau_vx * stau_py - stau_vy * stau_px) / stau_pt
                    
                    # Calculate longitudinal impact parameter (z0)
                    stau_z0 = stau_vz - (stau_vz * stau_pt / sqrt(stau_px**2 + stau_py**2 + stau_pz**2))

                    imcp_stau_d0.append(stau_d0)
                    mcp_stau_d0_rt.clear()
                    mcp_stau_d0_rt.push_back(mcp_tlv.Perp())
                    imcp_stau_z0.append(stau_z0)
                    mcp_stau_z0_rt.clear()
                    mcp_stau_z0_rt.push_back(mcp_tlv.Perp())
                    mcp_stau_tree.Fill() ### FIXME may get appended twice when filling later

                    #stau_reconstructable = False
                    r_vertex = sqrt(stau_vx ** 2 + stau_vy ** 2)
                    r_endpoint = sqrt(mcp.getEndpoint()[0] ** 2 + mcp.getEndpoint()[1] ** 2)
                    print("r_vertex, r_endpoint:", r_vertex, r_endpoint)
                    z_vertex = mcp.getVertex()[2]
                    z_endpoint = mcp.getEndpoint()[2]
                    stau_reconstructable = acceptanceCutsOld(mcp)
                    if stau_reconstructable:
                        num_reconstructable_stau_tracks += 1
                    print("Truth pt, eta, phi:", mcp_tlv.Perp(), mcp_tlv.Eta(), mcp_tlv.Phi())
                    
                    stauHasTrack = False ### Will remain false if doesn't pass low pt cut, but each stau definitely should.. 
                    if (mcp_tlv.Perp() > 0.5): # Remove ultra-low pt tracks    
                        tracks = lcRelation.getRelatedToObjects(mcp)
                        if (len(tracks) > 0):
                            stauHasTrack = True ### if there is a track related to this stau 
                            if (len(tracks) > 1):
                                num_dupes += 1 ### FIXME only stau duplicate tracks right now
                        two_stau_track_bool = False
                        for track in tracks:
                            ii_staux = []
                            ii_stauy = []
                            ii_stauz = []
                            ii_staur = []
                            LC_pixel_nhit = 0
                            LC_inner_nhit = 0
                            LC_outer_nhit = 0
                            lastLayer = -1
                            nLayersCrossed = 0
                            for hit in track.getTrackerHits():
                                # now decode hits
                                encoding = hit_collections[0].getParameters().getStringVal(pyLCIO.EVENT.LCIO.CellIDEncoding)
                                decoder = pyLCIO.UTIL.BitField64(encoding)
                                cellID = int(hit.getCellID0())
                                decoder.setValue(cellID)
                                detector = decoder["system"].value()
                                layer = decoder['layer'].value()
                                if (lastLayer != layer):
                                    print("increment nlayerscrossed")
                                    nLayersCrossed += 1 ### NOTE counting each of vertex doublet layers as individual layer
                                print("hit layer: ", layer)
                                position = hit.getPosition()
                                pos_x = position[0]
                                pos_y = position[1]
                                pos_z = position[2]
                                print("pos_x, pos_y, pos_z:", pos_x, pos_y, pos_z)
                                pos_r = sqrt(pos_x ** 2 + pos_y ** 2)
                                if detector == 1 or detector == 2:
                                    LC_pixel_nhit += 1
                                if detector == 3 or detector == 4:
                                    LC_inner_nhit += 1
                                if detector == 5 or detector == 6:
                                    LC_outer_nhit += 1
                                lastLayer = layer
                                
                                ii_staur.append(pos_r)
                                LC_stau_hit_r.clear()
                                LC_stau_hit_r.push_back(pos_r)
                                ii_staux.append(pos_x)
                                LC_stau_hit_x.clear()
                                LC_stau_hit_x.push_back(pos_x)
                                ii_stauy.append(pos_y) 
                                LC_stau_hit_y.clear()
                                LC_stau_hit_y.push_back(pos_y)
                                ii_stauz.append(pos_z)
                                LC_stau_hit_z.clear()
                                LC_stau_hit_z.push_back(pos_z)
                                stau_track_tree.Fill() ### FIXME may get appended twice when filling later
                            
                            theta = np.pi/2- np.arctan(track.getTanLambda())
                            phi = track.getPhi()
                            eta = -np.log(np.tan(theta/2))
                            pt  = 0.2998 * Bfield / fabs(track.getOmega() * 1000.)
                            track_tlv = ROOT.TLorentzVector()
                            track_tlv.SetPtEtaPhiE(pt, eta, phi, 0)
                            dr = mcp_tlv.DeltaR(track_tlv)
                            nhitz = track.getTrackerHits().size()
                            print("nhits for stau track: ", nhitz)
                            ptres = abs(mcp_tlv.Perp() - pt) / (mcp_tlv.Perp())
                            d0 = track.getD0()
                            z0 = track.getZ0()

                            stau_track_mcps = lcRelation.getRelatedFromObjects(track)
                            print("num mcps associated to stau track:", len(stau_track_mcps))
                            print("stau track pdgs:", [mcp.getPDG() for mcp in stau_track_mcps])
                            if (len(stau_track_mcps) > 1 & len(set(stau_track_mcps)) == 1):
                                two_stau_track_bool = True
                                print("two staus associated to this track")
                            
                            if (nLayersCrossed < 4):
                                print("stau track doesn't pass n layers / n hits cut")
                            else:
                                print("found matched track")
                                num_matched_stau_tracks += 1
                                imcp_two_stau_track_bool.append(two_stau_track_bool)
                                LC_stau_pt_match.append(mcp_tlv.Perp())
                                LC_stau_pt_match_rt.clear()
                                LC_stau_pt_match_rt.push_back(mcp_tlv.Perp())
                                LC_stau_track_pt.append([pt])
                                LC_stau_track_pt.clear()
                                LC_stau_track_pt.push_back(pt)
                                LC_stau_track_eta.append([eta])
                                LC_stau_track_eta.clear()
                                LC_stau_track_eta.push_back(eta)
                                LC_stau_eta_match.append(mcp_tlv.Eta())
                                LC_stau_eta_match.clear()
                                LC_stau_eta_match.push_back(mcp_tlv.Eta())
                                LC_stau_track_theta.append([theta])
                                LC_stau_track_theta.clear()
                                LC_stau_track_theta.push_back(theta)
                                LC_stau_phi_match.append([phi])
                                LC_stau_phi_match.clear()
                                LC_stau_phi_match.push_back(phi)
                                LC_stau_ndf.append([track.getNdf()])
                                LC_stau_ndf.clear()
                                LC_stau_ndf.push_back(track.getNdf())
                                LC_stau_chi2.append([track.getChi2()])
                                LC_stau_chi2.clear()
                                LC_stau_chi2.push_back(track.getChi2())
                                LC_stau_d0.append([d0])
                                LC_stau_d0.clear()
                                LC_stau_d0.push_back(d0)
                                LC_stau_z0.append([z0])
                                LC_stau_z0.clear()
                                LC_stau_z0.push_back(z0)
                                LC_stau_nhits.append([nhitz])
                                LC_stau_nhits.clear()
                                LC_stau_nhits.push_back(nhitz)
                                LC_stau_pt_res.append([ptres])
                                LC_stau_pt_res.clear()
                                LC_stau_pt_res.push_back(ptres)
                                LC_stau_dr.append([dr])
                                LC_stau_dr.clear()
                                LC_stau_dr.push_back(dr)
                                LC_stau_pixel_nhits.append([LC_pixel_nhit])
                                LC_stau_pixel_nhits.clear()
                                LC_stau_pixel_nhits.push_back(LC_pixel_nhit)
                                LC_stau_inner_nhits.append([LC_inner_nhit])
                                LC_stau_inner_nhits.clear()
                                LC_stau_inner_nhits.push_back(LC_inner_nhit)
                                LC_stau_outer_nhits.append([LC_outer_nhit])
                                LC_stau_outer_nhits.clear()
                                LC_stau_outer_nhits.push_back(LC_outer_nhit)
                                LC_stau_hit_r.append(ii_staur)
                                LC_stau_hit_x.append(ii_staux)
                                LC_stau_hit_y.append(ii_stauy)
                                LC_stau_hit_z.append(ii_stauz)
                                stau_track_tree.Fill()

                    imcp_stau_track_bool.append(stauHasTrack)
                    mcp_stau_track_bool_rt.clear()
                    mcp_stau_track_bool_rt.push_back(stauHasTrack)
                    imcp_stau_track_reconstructable_bool.append(stau_reconstructable)
                    mcp_stau_track_reconstructable_bool_rt.clear()
                    mcp_stau_track_reconstructable_bool_rt.push_back(stau_reconstructable)
                    mcp_stau_tree.Fill()
                    print("truth stau d0, z0, has matched track: ", stau_d0, stau_z0, stauHasTrack) ### end mcp loop 


                    ### START LOOP THROUGH DAUGHTER PARTICLES
                    mcpDaughters = mcp.getDaughters()
                    print("num stau daughters:", len(mcpDaughters))
                    for mcpDaughter in mcpDaughters: 
                        mcp_daughter_p = mcpDaughter.getMomentum()
                        mcp_daughter_tlv = ROOT.TLorentzVector()
                        mcp_daughter_tlv.SetPxPyPzE(mcp_daughter_p[0], mcp_daughter_p[1], mcp_daughter_p[2], mcpDaughter.getEnergy())
                        daughterHasTrack = False
                        daughter_reconstructable = False
                        if abs(mcpDaughter.getPDG())==1000015 or abs(mcpDaughter.getPDG())==2000015: ### IGNORE STAUS
                            continue
                        print("pdgid of daughter:", mcpDaughter.getPDG())
                        if (mcpDaughter.getCharge() != 0):
                            if mcpDaughter.getGeneratorStatus() == 0:
                                continue
                            print("pdgid of charged daughter:", mcpDaughter.getPDG())
                            print("gen status daughter 1:", mcpDaughter.getGeneratorStatus())
                            n_charged_mcp_daughter += 1 
                            imcp_daughter_pt.append(mcp_daughter_tlv.Perp()) 
                            print("daughter_pt:", imcp_daughter_pt)
                            imcp_daughter_eta.append(mcp_daughter_tlv.Eta())
                            imcp_daughter_phi.append(mcp_daughter_tlv.Phi())
                            #n_charged_mcp_daughter += 1
                            # Get the vertex position
                            daughter_vx, daughter_vy, daughter_vz = mcpDaughter.getVertex()[0], mcpDaughter.getVertex()[1], mcpDaughter.getVertex()[2]

                            print("daughter vertex position: ", daughter_vx, daughter_vy, daughter_vz)
                            
                            # Get the momentum
                            daughter_px, daughter_py, daughter_pz = mcpDaughter.getMomentum()[0], mcpDaughter.getMomentum()[1], mcpDaughter.getMomentum()[2]
                            
                            # Calculate transverse impact parameter (d0)
                            daughter_pt = sqrt(daughter_px**2 + daughter_py**2)  # Transverse momentum
                            print("daughter_pt:", daughter_pt)
                            daughter_d0 = (daughter_vx * daughter_py - daughter_vy * daughter_px) / daughter_pt
                            
                            # Calculate longitudinal impact parameter (z0)
                            daughter_z0 = daughter_vz - (daughter_vz * daughter_pt / sqrt(daughter_px**2 + daughter_py**2 + daughter_pz**2))

                            imcp_daughter_d0.append(daughter_d0)
                            imcp_daughter_z0.append(daughter_z0)

                            print("daughter d0, z0: ", daughter_d0, daughter_z0)

                            if acceptanceCutsOld(mcpDaughter):
                                daughter_reconstructable = True
                                print("recoable daughter 1")
                                n_recoable_daughter += 1
                                
                                ### Perform truth association to tracks (Note: likely won't work without additional seeding layers / loosening nhits)
                                tracks = lcRelation.getRelatedToObjects(mcpDaughter)
                                if (len(tracks) > 0):
                                    daughterHasTrack = True ### if there is a track related to this stau 
                                    if (len(tracks) > 1):
                                        num_dupes += 1 ### FIXME add for all tracks too
                                for track in tracks:
                                    theta = np.pi/2- np.arctan(track.getTanLambda())
                                    phi = track.getPhi()
                                    eta = -np.log(np.tan(theta/2))
                                    pt  = 0.2998 * Bfield / fabs(track.getOmega() * 1000.)
                                    track_tlv = ROOT.TLorentzVector()
                                    track_tlv.SetPtEtaPhiE(pt, eta, phi, 0)
                                    dr = mcp_daughter_tlv.DeltaR(track_tlv)
                                    nhitz = track.getTrackerHits().size()
                                    print("daughter hits: ", nhitz)
                                    ptres = abs(mcp_daughter_tlv.Perp() - pt) / (mcp_daughter_tlv.Perp())
                                    d0 = track.getD0()
                                    z0 = track.getZ0()
                                    LC_daughter_pt_match.append(mcp_daughter_tlv.Perp())
                                    LC_daughter_track_pt.append([pt])
                                    LC_daughter_track_eta.append([eta])
                                    LC_daughter_eta_match.append(mcp_daughter_tlv.Eta())
                                    LC_daughter_track_theta.append([theta])
                                    LC_daughter_phi_match.append([phi])
                                    LC_daughter_ndf.append([track.getNdf()])
                                    LC_daughter_chi2.append([track.getChi2()])
                                    LC_daughter_d0.append([d0])
                                    LC_daughter_z0.append([z0])
                                    LC_daughter_nhits.append([nhitz])
                                    LC_daughter_pt_res.append([ptres])
                                    LC_daughter_dr.append([dr])

                                    LC_pixel_nhit = 0
                                    LC_inner_nhit = 0
                                    LC_outer_nhit = 0
                                    for hit in track.getTrackerHits():
                                    # now decode hits
                                        encoding = hit_collections[0].getParameters().getStringVal(pyLCIO.EVENT.LCIO.CellIDEncoding)
                                        decoder = pyLCIO.UTIL.BitField64(encoding)
                                        cellID = int(hit.getCellID0())
                                        decoder.setValue(cellID)
                                        detector = decoder["system"].value()
                                        if detector == 1 or detector == 2:
                                            LC_pixel_nhit += 1
                                        if detector == 3 or detector == 4:
                                            LC_inner_nhit += 1
                                        if detector == 5 or detector == 6:
                                            LC_outer_nhit += 1
                                    LC_daughter_pixel_nhits.append([LC_pixel_nhit])
                                    LC_daughter_inner_nhits.append([LC_inner_nhit])
                                    LC_daughter_outer_nhits.append([LC_outer_nhit])
                                    print("found matched first daughter track")
                                    print("--------------------------------")
                                    num_matched_daughter_tracks += 1
                        imcp_daughter_track_bool.append(daughterHasTrack)
                        imcp_daughter_track_reconstructable_bool.append(daughter_reconstructable)
                        
                        daughterHasTrack = False
                        daughter_reconstructable = False
                        mcpDaughterDaughters = mcpDaughter.getDaughters()
                        for mcpDaughterDaughter in mcpDaughterDaughters:
                            mcp_daughterDaughter_p = mcpDaughter.getMomentum()
                            mcp_daughterDaughter_tlv = ROOT.TLorentzVector()
                            mcp_daughterDaughter_tlv.SetPxPyPzE(mcp_daughterDaughter_p[0], mcp_daughterDaughter_p[1], mcp_daughterDaughter_p[2], mcpDaughterDaughter.getEnergy())
                            if (mcpDaughterDaughter.getCharge() != 0):
                                if mcpDaughter.getGeneratorStatus() == 23: 
                                    continue
                                print("pdgid of daughter's charged daughter:", mcpDaughterDaughter.getPDG())
                                print("gen status daughter 2:", mcpDaughter.getGeneratorStatus())
                                n_charged_mcp_daughter += 1 
                                imcp_daughter_pt.append(mcp_daughterDaughter_tlv.Perp()) 
                                print("daughter_pt:", imcp_daughter_pt)
                                imcp_daughter_eta.append(mcp_daughterDaughter_tlv.Eta())
                                imcp_daughter_phi.append(mcp_daughterDaughter_tlv.Phi())

                                #n_charged_mcp_daughter += 1
                                # Get the vertex position
                                daughter_vx, daughter_vy, daughter_vz = mcpDaughterDaughter.getVertex()[0], mcpDaughterDaughter.getVertex()[1], mcpDaughterDaughter.getVertex()[2]

                                print("daughter vertex position: ", daughter_vx, daughter_vy, daughter_vz)
                                
                                # Get the momentum
                                daughter_px, daughter_py, daughter_pz = mcpDaughterDaughter.getMomentum()[0], mcpDaughterDaughter.getMomentum()[1], mcpDaughterDaughter.getMomentum()[2]
                                
                                # Calculate transverse impact parameter (d0)
                                daughter_pt = sqrt(daughter_px**2 + daughter_py**2)  # Transverse momentum
                                daughter_d0 = (daughter_vx * daughter_py - daughter_vy * daughter_px) / daughter_pt
                                
                                # Calculate longitudinal impact parameter (z0)
                                daughter_z0 = daughter_vz - (daughter_vz * daughter_pt / sqrt(daughter_px**2 + daughter_py**2 + daughter_pz**2))

                                imcp_daughter_d0.append(daughter_d0)
                                imcp_daughter_z0.append(daughter_z0)

                                print("daughter d0, z0: ", daughter_d0, daughter_z0)
                                if acceptanceCutsOld(mcpDaughterDaughter):
                                    daughter_reconstructable = True
                                    n_recoable_daughter += 1
                                    print("recoable daughter 2")
                                    print("is stau displaced product")

                                    ### Perform truth association to tracks (Note: likely won't work without additional seeding layers / loosening nhits)

                                    tracks = lcRelation.getRelatedToObjects(mcpDaughterDaughter)
                                    if (len(tracks) > 0):
                                        daughterHasTrack = True ### if there is a track related to this stau 
                                        if (len(tracks) > 1):
                                            num_dupes += 1 ### FIXME add for all tracks too
                                    for track in tracks: 
                                        theta = np.pi/2- np.arctan(track.getTanLambda())
                                        phi = track.getPhi()
                                        eta = -np.log(np.tan(theta/2))
                                        pt  = 0.2998 * Bfield / fabs(track.getOmega() * 1000.)
                                        track_tlv = ROOT.TLorentzVector()
                                        track_tlv.SetPtEtaPhiE(pt, eta, phi, 0)
                                        dr = mcp_daughterDaughter_tlv.DeltaR(track_tlv)
                                        nhitz = track.getTrackerHits().size()
                                        print("nhitz: ", nhitz)
                                        ptres = abs(mcp_daughterDaughter_tlv.Perp() - pt) / (mcp_daughterDaughter_tlv.Perp())
                                        d0 = track.getD0()
                                        z0 = track.getZ0()
                                        LC_daughter_pt_match.append(mcp_daughterDaughter_tlv.Perp())
                                        LC_daughter_track_pt.append([pt])
                                        LC_daughter_track_eta.append([eta])
                                        LC_daughter_eta_match.append(mcp_daughterDaughter_tlv.Eta())
                                        LC_daughter_track_theta.append([theta])
                                        LC_daughter_phi_match.append([phi])
                                        LC_daughter_ndf.append([track.getNdf()])
                                        LC_daughter_chi2.append([track.getChi2()])
                                        LC_daughter_d0.append([d0])
                                        LC_daughter_z0.append([z0])
                                        LC_daughter_nhits.append([nhitz])
                                        LC_daughter_pt_res.append([ptres])
                                        LC_daughter_dr.append([dr])

                                        LC_pixel_nhit = 0
                                        LC_inner_nhit = 0
                                        LC_outer_nhit = 0
                                        for hit in track.getTrackerHits():
                                        # now decode hits
                                            encoding = hit_collections[0].getParameters().getStringVal(pyLCIO.EVENT.LCIO.CellIDEncoding)
                                            decoder = pyLCIO.UTIL.BitField64(encoding)
                                            cellID = int(hit.getCellID0())
                                            decoder.setValue(cellID)
                                            detector = decoder["system"].value()
                                            if detector == 1 or detector == 2:
                                                LC_pixel_nhit += 1
                                            if detector == 3 or detector == 4:
                                                LC_inner_nhit += 1
                                            if detector == 5 or detector == 6:
                                                LC_outer_nhit += 1
                                        LC_daughter_pixel_nhits.append([LC_pixel_nhit])
                                        LC_daughter_inner_nhits.append([LC_inner_nhit])
                                        LC_daughter_outer_nhits.append([LC_outer_nhit])
                                        print("found matched second daughter track")
                                        print("--------------------------------")
                                        num_matched_daughter_tracks += 1
                        imcp_daughter_track_bool.append(daughterHasTrack)
                        imcp_daughter_track_reconstructable_bool.append(daughter_reconstructable)        
                ### FIXME add another recursion to account for pi0 decay products? kaon decay products? 

            if n_charged_mcp_daughter > 0:
                mcp_daughter_pt.append(imcp_daughter_pt)
                print("len(imcp_daughter_pt): ", len(imcp_daughter_pt))
                mcp_daughter_eta.append(imcp_daughter_eta)
                mcp_daughter_phi.append(imcp_daughter_phi) 
                mcp_daughter_d0.append(imcp_daughter_d0)
                mcp_daughter_z0.append(imcp_daughter_z0) 
                mcp_daughter_track_bool.append(imcp_daughter_track_bool)
                mcp_daughter_track_reconstructable_bool.append(imcp_daughter_track_reconstructable_bool)
            
            ### Fill mcp stau values after looping through MCPs            
            if n_mcp_stau > 0:
                mcp_stau_pt.append(imcp_stau_pt)
                print("len(imcp_stau_pt): ", len(imcp_stau_pt))
                mcp_stau_eta.append(imcp_stau_eta)
                mcp_stau_phi.append(imcp_stau_phi)          
                mcp_stau_d0.append(imcp_stau_d0)
                mcp_stau_z0.append(imcp_stau_z0)  
                mcp_stau_track_bool.append(imcp_stau_track_bool)
                mcp_two_stau_track_bool.append(imcp_two_stau_track_bool)
                mcp_stau_track_reconstructable_bool.append(imcp_stau_track_reconstructable_bool)

            # Loop over the track objects
            for track in trackCollection:
                isFakeTrack = False
                # Hit per track
                iix = []
                iiy = []
                iiz = []
                iihit_pdg = []
                iitime = []
                iicorrected_time = []
                iihit_layer = []
                iihit_detector = []
                iihit_side = []
                
                pixel_nhit = 0
                inner_nhit = 0
                outer_nhit = 0
                lastLayer = -1
                nLayersCrossed = 0 
                for hit in track.getTrackerHits():
                    try:
                        mcp = hit.getMCParticle()
                        hit_pdg = mcp.getPDG()
                    except:
                        hit_pdg = 0
                    # now decode hits
                    encoding = hit_collections[0].getParameters().getStringVal(pyLCIO.EVENT.LCIO.CellIDEncoding)
                    decoder = pyLCIO.UTIL.BitField64(encoding)
                    cellID = int(hit.getCellID0())
                    decoder.setValue(cellID)
                    detector = decoder["system"].value()
                    layer = decoder['layer'].value()
                    side = decoder["side"].value()
                    if (lastLayer != layer):
                        nLayersCrossed += 1 ### NOTE counting each of vertex doublet layers as individual layer
                    if detector == 1 or detector == 2:
                        pixel_nhit += 1
                    if detector == 3 or detector == 4:
                        inner_nhit += 1
                    if detector == 5 or detector == 6:
                        outer_nhit += 1
                    position = hit.getPosition()
                    pos_x = position[0]
                    pos_y = position[1]
                    pos_z = position[2]

                    d = sqrt(position[0]*position[0] + position[1]
                    * position[1] + position[2]*position[2])
                    tof = d/speedoflight

                    resolution = 0.03
                    if detector > 2:
                        resolution = 0.06

                    corrected_t = hit.getTime()*(1.+ROOT.TRandom3(ievt).Gaus(0., resolution)) - tof
                    
                    iix.append(pos_x)
                    iiy.append(pos_y)
                    iiz.append(pos_z)
                    iihit_pdg.append(hit_pdg)
                    iitime.append(hit.getTime())
                    iicorrected_time.append(corrected_t)
                    iihit_detector.append(detector)
                    iihit_layer.append(layer)
                    iihit_side.append(side)
                    lastLayer = layer
                if(nLayersCrossed < 4):
                    print("track did not pass n layers / n hits cut")
                    continue
                track_mcps = lcRelation.getRelatedFromObjects(track)
                if len(track_mcps) < 1:
                    isFakeTrack = True
                if len(track_mcps) > 1: 
                    pdg_ids = [mcp.getPDG() for mcp in track_mcps]
                    print("fake pdgids?:", pdg_ids)
                    # Check if all PDG IDs are the same
                    print(len(set(pdg_ids)))
                    if len(set(pdg_ids)) > 1:
                        num_fake_tracks += 1
                        isFakeTrack = True
                        print("fake track with pdgids:", pdg_ids)
                    
                if isFakeTrack: ### FAKE TRACK (NO MCPS RELATED)
                    theta = np.pi/2- np.arctan(track.getTanLambda())
                    phi = track.getPhi()
                    eta = -np.log(np.tan(theta/2))
                    pt  = 0.2998 * Bfield / fabs(track.getOmega() * 1000.)
                    track_tlv = ROOT.TLorentzVector()
                    track_tlv.SetPtEtaPhiE(pt, eta, phi, 0)
                    nhitz = track.getTrackerHits().size()
                    d0 = track.getD0()
                    z0 = track.getZ0()
                    fake_pt.append([pt])
                    fake_eta.append([eta])
                    fake_theta.append([theta])
                    fake_ndf.append([track.getNdf()])
                    fake_chi2.append([track.getChi2()])
                    fake_d0.append([d0])
                    fake_z0.append([z0])
                    print("fake pt:", pt)
                    print("fake nhitz:", nhitz)
                    print("fake_chi2:", track.getChi2())
                    print("fake reduced chi2", track.getChi2() / float(track.getNdf()))
                    fake_nhits.append([nhitz])

                    LC_pixel_nhit = 0
                    LC_inner_nhit = 0
                    LC_outer_nhit = 0
                    for hit in track.getTrackerHits():
                    # now decode hits
                        encoding = hit_collections[0].getParameters().getStringVal(pyLCIO.EVENT.LCIO.CellIDEncoding)
                        decoder = pyLCIO.UTIL.BitField64(encoding)
                        cellID = int(hit.getCellID0())
                        decoder.setValue(cellID)
                        detector = decoder["system"].value()
                        if detector == 1 or detector == 2:
                            LC_pixel_nhit += 1
                        if detector == 3 or detector == 4:
                            LC_inner_nhit += 1
                        if detector == 5 or detector == 6:
                            LC_outer_nhit += 1
                    fake_pixel_nhits.append([LC_pixel_nhit])
                    fake_inner_nhits.append([LC_inner_nhit])
                    fake_outer_nhits.append([LC_outer_nhit])
                    num_fake_tracks += 1
                ### FIXME find way to get mcp from reco. tracks (not just staus?), use 'from' relation? 
                print("len(mcps_track):", len(track_mcps))
                for track_mcp in track_mcps:
                    print("track mcps: ", track_mcp.getPDG())
                
                theta = np.pi/2- np.arctan(track.getTanLambda())
                phi = track.getPhi()
                eta = -np.log(np.tan(theta/2))
                pt  = 0.2998 * Bfield / fabs(track.getOmega() * 1000.)
                track_tlv = ROOT.TLorentzVector()
                track_tlv.SetPtEtaPhiE(pt, eta, phi, 0)
                nhitz = track.getTrackerHits().size()

                d0 = track.getD0()
                z0 = track.getZ0()

                inhits.append(nhitz)
                itrack_pt.append(pt)
                itrack_eta.append(eta)
                itrack_theta.append(theta)
                indf.append(track.getNdf())
                ichi2.append(track.getChi2())
                print("track reduced chi2", track.getChi2() / float(track.getNdf()))
                # print("Reco pt, eta, phi, nhits, dr:", pt, eta, phi, nhitz, dr)

                pixel_nhits.append([pixel_nhit])
                inner_nhits.append([inner_nhit])
                outer_nhits.append([outer_nhit])
                ix.append(iix)
                iy.append(iiy)
                iz.append(iiz)
                ihit_pdg.append(iihit_pdg)
                itime.append(iitime)
                icorrected_time.append(iicorrected_time)
                ihit_detector.append(iihit_detector)
                ihit_layer.append(iihit_layer)
                ihit_side.append(iihit_side)
                            
            # print("End of event \n")
            # This is here to check that we never reconstruct multiple muons
            # If we did, we'd have to match the correct muon to the MCP object to do eff/res plots
            # But since we don't, we can skip that step
            i+=1
            nhits.append(inhits)
            track_pt.append(itrack_pt)
            track_eta.append(itrack_eta)
            track_theta.append(itrack_theta)
            ndf.append(indf)
            chi2.append(ichi2)
            x.append(ix)
            y.append(iy)
            z.append(iz)
            hit_pdgid.append(ihit_pdg)
            time.append(itime)
            corrected_time.append(icorrected_time)
            hit_layer.append(ihit_layer)
            hit_detector.append(ihit_detector)
            hit_side.append(ihit_side)
            mcp_pt.append(imcp_pt)
            mcp_eta.append(imcp_eta)
            mcp_phi.append(imcp_phi)
            pdgid.append(ipdgid)
            status.append(istatus)
            prod_time.append(iprod_time)
            id.append(iid)

            """reco_VB_x.append(ireco_VB_x); reco_VB_y.append(ireco_VB_y); reco_VB_z.append(ireco_VB_z); reco_VB_time.append(ireco_VB_time); reco_VB_layer.append(ireco_VB_layer)
            reco_VE_x.append(ireco_VE_x); reco_VE_y.append(ireco_VE_y); reco_VE_z.append(ireco_VE_z); reco_VE_time.append(ireco_VE_time); reco_VE_layer.append(ireco_VE_layer)
            reco_IB_x.append(ireco_IB_x); reco_IB_y.append(ireco_IB_y); reco_IB_z.append(ireco_IB_z); reco_IB_time.append(ireco_IB_time); reco_IB_layer.append(ireco_IB_layer)
            reco_IE_x.append(ireco_IE_x); reco_IE_y.append(ireco_IE_y); reco_IE_z.append(ireco_IE_z); reco_IE_time.append(ireco_IE_time); reco_IE_layer.append(ireco_IE_layer)
            reco_OB_x.append(ireco_OB_x); reco_OB_y.append(ireco_OB_y); reco_OB_z.append(ireco_OB_z); reco_OB_time.append(ireco_OB_time); reco_OB_layer.append(ireco_OB_layer)
            reco_OE_x.append(ireco_OE_x); reco_OE_y.append(ireco_OE_y); reco_OE_z.append(ireco_OE_z); reco_OE_time.append(ireco_OE_time); reco_OE_layer.append(ireco_OE_layer)"""
        
        reader.close()

        # ############## MANIPULATE, PRETTIFY, AND SAVE HISTOGRAMS #############################
        print("\nSummary statistics:")
        print("for: ", sample)
        print("Ran over %i events."%i)
        print("Found:")
        #print("\t%i MCPs"%len(np.ravel(mcp_pt)))
        print("\t%i Stau MCPs"%n_mcp_stau)
        print('\t%i reconstructable stau tracks'%(num_reconstructable_stau_tracks))
        print('\t%i matched stau tracks'%(num_matched_stau_tracks))
        stau_eff = float(num_matched_stau_tracks) / float(n_mcp_stau)
        stau_eff_no_acc = float(num_matched_stau_tracks) / (float(num_reconstructable_stau_tracks))
        print("Approx. Stau Total Tracking Eff: ", stau_eff)
        print("Stau Total Tracking Eff w/ Acceptance: ", stau_eff_no_acc) 
        print("\t num charged stau decay products per event:", n_charged_mcp_daughter / max_events)
        print("\t n recoable stau decay products:", n_recoable_daughter / max_events)
        print('\t%i matched daughter tracks per event'%(num_matched_daughter_tracks))
        daughter_eff = float(num_matched_daughter_tracks) / float(n_charged_mcp_daughter) ### FIXME divide by ALL (not just charged) decay products for total eff w/o acc?
        daughter_eff_acc = float(num_matched_daughter_tracks) / float(n_recoable_daughter)
        print("Approx. Stau Total Tracking Eff: ", daughter_eff)
        print("Stau Total Tracking Eff w/ Acceptance: ", daughter_eff_acc) 
        print('\t%i fake tracks per event'%num_fake_tracks)
        #print('\t%i GeV'%np.max(mcp_stau_pt))

        # Make a list of all the data you want to save

        data_list = {
            #"sim_VB_x": sim_VB_x, "sim_VB_y": sim_VB_y, "sim_VB_z": sim_VB_z, "sim_VB_time": sim_VB_time, "sim_VB_pdg": sim_VB_pdg, "sim_VB_mcpid": sim_VB_mcpid, "sim_VB_layer": sim_VB_layer,
            #"sim_VE_x": sim_VE_x, "sim_VE_y": sim_VE_y, "sim_VE_z": sim_VE_z, "sim_VE_time": sim_VE_time, "sim_VE_pdg": sim_VE_pdg, "sim_VE_mcpid": sim_VE_mcpid, "sim_VE_layer": sim_VE_layer,
            #"sim_IB_x": sim_IB_x, "sim_IB_y": sim_IB_y, "sim_IB_z": sim_IB_z, "sim_IB_time": sim_IB_time, "sim_IB_pdg": sim_IB_pdg, "sim_IB_mcpid": sim_IB_mcpid, "sim_IB_layer": sim_IB_layer,
            #"sim_IE_x": sim_IE_x, "sim_IE_y": sim_IE_y, "sim_IE_z": sim_IE_z, "sim_IE_time": sim_IE_time, "sim_IE_pdg": sim_IE_pdg, "sim_IE_mcpid": sim_IE_mcpid, "sim_IE_layer": sim_IE_layer,
            #"sim_OB_x": sim_OB_x, "sim_OB_y": sim_OB_y, "sim_OB_z": sim_OB_z, "sim_OB_time": sim_OB_time, "sim_OB_pdg": sim_OB_pdg, "sim_OB_mcpid": sim_OB_mcpid, "sim_OB_layer": sim_OB_layer,
            #"sim_OE_x": sim_OE_x, "sim_OE_y": sim_OE_y, "sim_OE_z": sim_OE_z, "sim_OE_time": sim_OE_time, "sim_OE_pdg": sim_OE_pdg, "sim_OE_mcpid": sim_OE_mcpid, "sim_OE_layer": sim_OE_layer,
            #"reco_VB_x": reco_VB_x, "reco_VB_y": reco_VB_y, "reco_VB_z": reco_VB_z, "reco_VB_time": reco_VB_time, "reco_VB_layer": reco_VB_layer,
            #"reco_VE_x": reco_VE_x, "reco_VE_y": reco_VE_y, "reco_VE_z": reco_VE_z, "reco_VE_time": reco_VE_time, "reco_VE_layer": reco_VE_layer,
            #"reco_IB_x": reco_IB_x, "reco_IB_y": reco_IB_y, "reco_IB_z": reco_IB_z, "reco_IB_time": reco_IB_time, "reco_IB_layer": reco_IB_layer,
            #"reco_IE_x": reco_IE_x, "reco_IE_y": reco_IE_y, "reco_IE_z": reco_IE_z, "reco_IE_time": reco_IE_time, "reco_IE_layer": reco_IE_layer,
            #"reco_OB_x": reco_OB_x, "reco_OB_y": reco_OB_y, "reco_OB_z": reco_OB_z, "reco_OB_time": reco_OB_time, "reco_OB_layer": reco_OB_layer,
            #"reco_OE_x": reco_OE_x, "reco_OE_y": reco_OE_y, "reco_OE_z": reco_OE_z, "reco_OE_time": reco_OE_time, "reco_OE_layer": reco_OE_layer
        }

        ### FOR ALL MCPS
        data_list["mcp_pt"] = mcp_pt
        data_list["mcp_eta"] = mcp_eta
        data_list["mcp_phi"] = mcp_phi
        data_list["pdgid"] = pdgid
        data_list["status"] = status
        data_list["prod_vertex"] = prod_vertex
        data_list["prod_endpoint"] = prod_endpoint
        data_list["prod_time"] = prod_time
        data_list["id"] = id

        ### FOR STAU MCPS
        data_list["mcp_stau_pt"] = mcp_stau_pt
        data_list["mcp_stau_eta"] = mcp_stau_eta
        data_list["mcp_stau_phi"] = mcp_stau_phi
        data_list["mcp_stau_d0"] = mcp_stau_d0
        data_list["mcp_stau_z0"] = mcp_stau_z0
        data_list["mcp_stau_track_bool"] = mcp_stau_track_bool
        data_list["mcp_stau_track_reconstructable_bool"] = mcp_stau_track_reconstructable_bool

        ### FOR STAU CHARGED DECAY PRODUCTS
        data_list["mcp_daughter_pt"] = mcp_daughter_pt
        data_list["mcp_daughter_eta"] = mcp_daughter_eta
        data_list["mcp_daughter_phi"] = mcp_daughter_phi
        data_list["mcp_daughter_d0"] = mcp_daughter_d0
        data_list["mcp_daughter_z0"] = mcp_daughter_z0
        data_list["mcp_daughter_track_bool"] = mcp_daughter_track_bool
        data_list["mcp_daughter_track_reconstructable_bool"] = mcp_daughter_track_reconstructable_bool

        ### FOR STAU MATCHED TRACKS
        data_list["LC_stau_pt_match"] = LC_stau_pt_match
        data_list["LC_stau_track_pt"] = LC_stau_track_pt
        data_list["LC_stau_track_eta"] = LC_stau_track_eta
        data_list["LC_stau_eta_match"] = LC_stau_eta_match
        data_list["LC_stau_track_theta"] = LC_stau_track_theta
        data_list["LC_stau_phi_match"] = LC_stau_phi_match
        data_list["LC_stau_ndf"] = LC_stau_ndf
        data_list["LC_stau_chi2"] = LC_stau_chi2
        data_list["LC_stau_d0"] = LC_stau_d0
        data_list["LC_stau_z0"] = LC_stau_z0
        data_list["LC_stau_nhits"] = LC_stau_nhits
        data_list["LC_stau_pixel_nhits"] = LC_stau_pixel_nhits
        data_list["LC_stau_inner_nhits"] = LC_stau_inner_nhits
        data_list["LC_stau_outer_nhits"] = LC_stau_outer_nhits
        data_list["LC_stau_pt_res"] = LC_stau_pt_res
        data_list["LC_stau_dr"] = LC_stau_dr
        data_list["LC_stau_hit_r"] = LC_stau_hit_r
        data_list["LC_stau_hit_x"] = LC_stau_hit_x
        data_list["LC_stau_hit_y"] = LC_stau_hit_y
        data_list["LC_stau_hit_z"] = LC_stau_hit_z
        data_list["mcp_two_stau_track_bool"] = mcp_two_stau_track_bool

        ### FOR STAU CHARGED DECAY PRODUCTS MATCHED TRACKS
        data_list["LC_daughter_pt_match"] = LC_daughter_pt_match
        data_list["LC_daughter_track_pt"] = LC_daughter_track_pt
        data_list["LC_daughter_track_eta"] = LC_daughter_track_eta
        data_list["LC_daughter_eta_match"] = LC_daughter_eta_match
        data_list["LC_daughter_track_theta"] = LC_daughter_track_theta
        data_list["LC_daughter_phi_match"] = LC_daughter_phi_match
        data_list["LC_daughter_ndf"] = LC_daughter_ndf
        data_list["LC_daughter_chi2"] = LC_daughter_chi2
        data_list["LC_daughter_d0"] = LC_daughter_d0
        data_list["LC_daughter_z0"] = LC_daughter_z0
        data_list["LC_daughter_nhits"] = LC_daughter_nhits
        data_list["LC_daughter_pixel_nhits"] = LC_daughter_pixel_nhits
        data_list["LC_daughter_inner_nhits"] = LC_daughter_inner_nhits
        data_list["LC_daughter_outer_nhits"] = LC_daughter_outer_nhits
        data_list["LC_daughter_pt_res"] = LC_daughter_pt_res
        data_list["LC_daughter_dr"] = LC_daughter_dr
        
        ### FOR ALL TRACKS
        data_list["nhits"] = nhits
        data_list["pixel_nhits"] = pixel_nhits
        data_list["inner_nhits"] = inner_nhits
        data_list["outer_nhits"] = outer_nhits
        data_list["track_pt"] = track_pt
        data_list["track_eta"] = track_eta
        data_list["track_theta"] = track_theta
        data_list["ndf"] = ndf
        data_list["chi2"] = chi2

        ### For fake tracks
        data_list["fake_pt"] = fake_theta
        data_list["fake_theta"] = fake_theta
        data_list["fake_eta"] = fake_theta
        data_list["fake_phi"] = fake_theta
        data_list["fake_d0"] = fake_theta
        data_list["fake_z0"] = fake_theta
        data_list["fake_ndf"] = fake_theta
        data_list["fake_chi2"] = fake_theta
        data_list["fake_pixel_nhits"] = fake_theta
        data_list["fake_inner_nhits"] = fake_theta
        data_list["fake_outer_nhits"] = fake_theta

        ### FOR HITS
        data_list["x"] = x
        data_list["y"] = y
        data_list["z"] = z
        #data_list["hit_pdgid"] = hit_pdgid ### FIXME all values are 0 
        data_list["time"] = time
        data_list["corrected_time"] = corrected_time
        data_list["hit_layer"] = hit_layer
        data_list["hit_detector"] = hit_detector
        data_list["hit_side"] = hit_side

        # After the loop is finished, save the data_list to a .json file
        output_pattern = f"{base_path}{sample}_reco.json"
        with open(output_pattern, 'w') as fp:
            json.dump(data_list, fp)

        # Write the tree to the file
        file.Write()
        file.Close()