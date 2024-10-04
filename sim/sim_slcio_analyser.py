import pyLCIO
import ROOT
import glob
import json
from math import *
import numpy as np


# ############## SETUP #############################
# Prevent ROOT from drawing while you're running -- good for slow remote servers
# Instead, save files and view them with an sftp client like Fetch (feel free to ask me for my UTK license)
ROOT.gROOT.SetBatch()

# Set up some options
max_events = -1 # Set to -1 to run over all events

# Gather input files
fnames = glob.glob("/home/larsonma/MarkLLPCode/sim/2500_0.1_sim.slcio") 
print("Found %i files."%len(fnames))

# Create empty lists for each variable
mcp_pt = [] #mcp = monte-carlo particle (truth)
mcp_phi = []
mcp_eta = []
pdgid = []
status = []
prod_vertex = []
prod_time = []
id = []

VB_x, VB_y, VB_z, VB_time, VB_pdg, VB_mcpid = [], [], [], [], [], []
VE_x, VE_y, VE_z, VE_time, VE_pdg, VE_mcpid = [], [], [], [], [], []
IB_x, IB_y, IB_z, IB_time, IB_pdg, IB_mcpid = [], [], [], [], [], []
IE_x, IE_y, IE_z, IE_time, IE_pdg, IE_mcpid = [], [], [], [], [], []
OB_x, OB_y, OB_z, OB_time, OB_pdg, OB_mcpid = [], [], [], [], [], []
OE_x, OE_y, OE_z, OE_time, OE_pdg, OE_mcpid = [], [], [], [], [], []

i = 0
# reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
# reader.setReadCollectionNames(["MCParticle", "PandoraPFOs", "SiTracks", "SiTracks_Refitted", "MCParticle_SiTracks", "MCParticle_SiTracks_Refitted", "IBTrackerHits", "IETrackerHits", "OBTrackerHits", "OETrackerHits", "VBTrackerHits", "VETrackerHits"])
# ############## LOOP OVER EVENTS AND FILL HISTOGRAMS  #############################
# Loop over events
for f in fnames:
    if max_events > 0 and i >= max_events: break
    reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(f)
    for event in reader: 
        if max_events > 0 and i >= max_events: break
        if i%10 == 0: print("Processing event %i."%i)

        mcpCollection = event.getCollection("MCParticle")

        # Make counter variables
        n_mcp_mu = 0

        # MCPs
        imcp_pt = []
        imcp_eta = []
        imcp_phi = []
        ipdgid = []
        istatus = []
        iprod_vertex = []
        iprod_time = []
        iid = []

         # Sim Hits
        iVB_x, iVB_y, iVB_z, iVB_time, iVB_pdg, iVB_mcpid = [], [], [], [], [], []
        iVE_x, iVE_y, iVE_z, iVE_time, iVE_pdg, iVE_mcpid = [], [], [], [], [], []
        iIB_x, iIB_y, iIB_z, iIB_time, iIB_pdg, iIB_mcpid = [], [], [], [], [], []
        iIE_x, iIE_y, iIE_z, iIE_time, iIE_pdg, iIE_mcpid = [], [], [], [], [], []
        iOB_x, iOB_y, iOB_z, iOB_time, iOB_pdg, iOB_mcpid = [], [], [], [], [], []
        iOE_x, iOE_y, iOE_z, iOE_time, iOE_pdg, iOE_mcpid = [], [], [], [], [], []

        # Loop over the truth objects and fill histograms
        for mcp in mcpCollection:
            mcp_p = mcp.getMomentum()
            mcp_tlv = ROOT.TLorentzVector()
            mcp_tlv.SetPxPyPzE(mcp_p[0], mcp_p[1], mcp_p[2], mcp.getEnergy())
            pdg = mcp.getPDG()

            imcp_pt.append(mcp_tlv.Perp())
            imcp_eta.append(mcp_tlv.Eta())
            imcp_phi.append(mcp_tlv.Phi())
            ipdgid.append(pdg)
            istatus.append(mcp.getGeneratorStatus())
            iprod_vertex.append([mcp.getVertex()[i] for i in range(3)])
            iprod_time.append(mcp.getTime())
            iid.append(mcp.id())
            # print("PID, Status, pT, eta, phi, prod_vertex, prod_time, id: ", pdg, mcp.getGeneratorStatus(), mcp_tlv.Perp(), mcp_tlv.Eta(), mcp_tlv.Phi(), [mcp.getVertex()[i] for i in range(3)], mcp.getTime(), mcp.id())

        for coll_name, (ix, iy, iz, itime, ipdg, imcpid) in [
            ("VertexBarrelCollection", (iVB_x, iVB_y, iVB_z, iVB_time, iVB_pdg, iVB_mcpid)),
            ("VertexEndcapCollection", (iVE_x, iVE_y, iVE_z, iVE_time, iVE_pdg, iVE_mcpid)),
            ("InnerTrackerBarrelCollection", (iIB_x, iIB_y, iIB_z, iIB_time, iIB_pdg, iIB_mcpid)),
            ("InnerTrackerEndcapCollection", (iIE_x, iIE_y, iIE_z, iIE_time, iIE_pdg, iIE_mcpid)),
            ("OuterTrackerBarrelCollection", (iOB_x, iOB_y, iOB_z, iOB_time, iOB_pdg, iOB_mcpid)),
            ("OuterTrackerEndcapCollection", (iOE_x, iOE_y, iOE_z, iOE_time, iOE_pdg, iOE_mcpid))
        ]:
            try:
                collection = event.getCollection(coll_name)
                for hit in collection:
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

            except Exception as e:
                print(f"Error accessing {coll_name}: {e}")                
        # print("End of event \n")
        i += 1
        mcp_pt.append(imcp_pt)
        mcp_eta.append(imcp_eta)
        mcp_phi.append(imcp_phi)
        pdgid.append(ipdgid)
        status.append(istatus)
        prod_vertex.append(iprod_vertex)
        prod_time.append(iprod_time)
        id.append(iid)
        VB_x.append(iVB_x); VB_y.append(iVB_y); VB_z.append(iVB_z); VB_time.append(iVB_time); VB_pdg.append(iVB_pdg); VB_mcpid.append(iVB_mcpid)
        VE_x.append(iVE_x); VE_y.append(iVE_y); VE_z.append(iVE_z); VE_time.append(iVE_time); VE_pdg.append(iVE_pdg); VE_mcpid.append(iVE_mcpid)
        IB_x.append(iIB_x); IB_y.append(iIB_y); IB_z.append(iIB_z); IB_time.append(iIB_time); IB_pdg.append(iIB_pdg); IB_mcpid.append(iIB_mcpid)
        IE_x.append(iIE_x); IE_y.append(iIE_y); IE_z.append(iIE_z); IE_time.append(iIE_time); IE_pdg.append(iIE_pdg); IE_mcpid.append(iIE_mcpid)
        OB_x.append(iOB_x); OB_y.append(iOB_y); OB_z.append(iOB_z); OB_time.append(iOB_time); OB_pdg.append(iOB_pdg); OB_mcpid.append(iOB_mcpid)
        OE_x.append(iOE_x); OE_y.append(iOE_y); OE_z.append(iOE_z); OE_time.append(iOE_time); OE_pdg.append(iOE_pdg); OE_mcpid.append(iOE_mcpid)


       
    reader.close()

# ############## MANIPULATE, PRETTIFY, AND SAVE HISTOGRAMS #############################
print("\nSummary statistics:")
print("Ran over %i events."%i)
print("Found:")
print("\tSanity check # of staus:", len(mcp_pt[pdgid == 1000015]))

# Make a list of all the data you want to save
data_list = {
    "VB_x": VB_x, "VB_y": VB_y, "VB_z": VB_z, "VB_time": VB_time, "VB_pdg": VB_pdg, "VB_mcpid": VB_mcpid,
    "VE_x": VE_x, "VE_y": VE_y, "VE_z": VE_z, "VE_time": VE_time, "VE_pdg": VE_pdg, "VE_mcpid": VE_mcpid,
    "IB_x": IB_x, "IB_y": IB_y, "IB_z": IB_z, "IB_time": IB_time, "IB_pdg": IB_pdg, "IB_mcpid": IB_mcpid, 
    "IE_x": IE_x, "IE_y": IE_y, "IE_z": IE_z, "IE_time": IE_time, "IE_pdg": IE_pdg, "IE_mcpid": IE_mcpid,
    "OB_x": OB_x, "OB_y": OB_y, "OB_z": OB_z, "OB_time": OB_time, "OB_pdg": OB_pdg, "OB_mcpid": OB_mcpid,
    "OE_x": OE_x, "OE_y": OE_y, "OE_z": OE_z, "OE_time": OE_time, "OE_pdg": OE_pdg, "OE_mcpid": OE_mcpid
}
data_list["mcp_pt"] = mcp_pt
data_list["mcp_eta"] = mcp_eta
data_list["mcp_phi"] = mcp_phi
data_list["pdgid"] = pdgid
data_list["status"] = status
data_list["prod_vertex"] = prod_vertex
data_list["prod_time"] = prod_time
data_list["id"] = id

# After the loop is finished, save the data_list to a .json file
output_json = "/home/larsonma/MarkLLPCode/sim/2500_0.1_sim.json"
with open(output_json, 'w') as fp:
    json.dump(data_list, fp)