import sys, os, time
start = time.time()
import argparse
import ROOT 
from ROOT import * #RDataFrame

# Enable multi-threading
ROOT.ROOT.EnableImplicitMT()

# Batch mode
ROOT.gROOT.SetBatch(True)

ROOT.gInterpreter.Declare("""
    #include "Utilities.h"
""")

#MC
#path = "/lustre/cms/store/user/fsimone/DsToPhiPi_ToMuMu_MuFilter_TuneCP5_13TeV-pythia8-evtgen/SkimPhiPi_Summer20UL18_DsPhiPi_ModFilter_Mini_v0/210131_142127/0000/Tree_PhiPi_1.root"
#tree_name = "Tree3Mu/ntuple"

#Data
path = "/lustre/cms/store/user/fsimone/DoubleMuonLowMass/SkimPhiPi_UL2018_Run2018C_ModFilter_Mini_v1/210131_134455/0000/" 
tree_name = "Tree3Mu/ntuple"

#generating the list of all .root files in given directory and subdirectories
chain = []
for r, d, f in os.walk(path): # r=root, d=directories, f = files
    for file in f:
        if '.root' in file:
            chain.append(os.path.join(r, file))

df = ROOT.RDataFrame(tree_name, chain)

#entries = df.Count()
#print("total ",entries.GetValue())
#print(df.GetColumnNames())
#print(df.GetColumnType("Trigger_hltname"))

# DsPhiPi ANALYSIS CUTFLOW

# 1 -> Event fires L1 and HLT

df = df.Define("HLT_mask", "Trigger_hltdecision==1 && contains(Trigger_hltname,\"HLT_DoubleMu3_Trk_Tau3mu_v\")").Filter("ROOT::VecOps::Sum(HLT_mask) >0")
df = df.Define("L1_mask", "Trigger_l1decision==1 && ( contains(Trigger_l1name,\"L1_DoubleMu\") || contains(Trigger_l1name,\"L1_TripleMu\"))").Filter("ROOT::VecOps::Sum(L1_mask) >0")

# Selections on triplets
# 2 -> 2mu+track candidate mass in (1.62-2.02)GeV
# 3 -> at least 2 track associated with PV
# 4 -> Significance of BS-SV distance in the transverse plane > 2
triplet_selection = "Triplet2_Mass>1.62 && Triplet2_Mass<2.02 && \
                     RefittedPV2_NTracks > 1 && \
                     FlightDistBS_SV_Significance > 2 "

# Events with at least one good candidate
df = df.Define("triplet_mask1", triplet_selection).Filter("ROOT::VecOps::Sum(triplet_mask1) >0") 

# 5 -> Muons and track within CMS acceptance
acceptance_selection = "((abs(Mu01_Eta)<1.2 && Mu01_Pt>3.5) || (abs(Mu01_Eta)>=1.2 && abs(Mu01_Eta)<2.4 && Mu01_Pt>2.0)) &&\
                        ((abs(Mu02_Eta)<1.2 && Mu02_Pt>3.5) || (abs(Mu02_Eta)>=1.2 && abs(Mu02_Eta)<2.4 && Mu02_Pt>2.0)) &&\
                        Tr_Pt>1.2"
# Events with at least one good candidate
df = df.Define("triplet_mask2", acceptance_selection).Filter("ROOT::VecOps::Sum(triplet_mask2)>0")

# Compute dR and dZ between muons/tracks
df = df.Define("dR12", "deltaR_vec(Mu01_Eta, Mu02_Eta, Mu01_Phi, Mu02_Phi)")
df = df.Define("dR13", "deltaR_vec(Mu01_Eta, Tr_Eta, Mu01_Phi, Tr_Phi)")
df = df.Define("dR23", "deltaR_vec(Mu02_Eta, Tr_Eta, Mu02_Phi, Tr_Phi)")

# 6 -> min and max deltaR requirement
dR_selection = "dR12>DELTAR_MIN && dR13>DELTAR_MIN && dR23>DELTAR_MIN &&\
                dR12<DELTAR_MAX && dR13<DELTAR_MAX && dR23<DELTAR_MAX"
df = df.Define("triplet_mask3", dR_selection).Filter("ROOT::VecOps::Sum(triplet_mask3)>0")

# Find index in "Muon_" and "Track_" branches
df = df.Define("Mu01_index", "match(MuonPt, Mu01_Pt)")
df = df.Define("Mu02_index", "match(MuonPt, Mu02_Pt)")
df = df.Define("Tr_index",   "match(MuonPt, Tr_Pt)")

# 7 -> Apply Muon ID Global and Particle Flow
df = df.Define("Mu01_ID", "muon_id(Mu01_index, Muon_isGlobal && Muon_isPF)")
df = df.Define("Mu02_ID", "muon_id(Mu02_index, Muon_isGlobal && Muon_isPF)")

# 8 -> IP(track, BS) z direction < 20 cm and xy direction < 0.3 cm
df = df.Define("Tr_IPcut", "muon_id(Tr_index, (Track_dz<20 && Track_dxy<0.3) )")
df = df.Define("triplet_mask4", "Mu01_ID && Mu02_ID && Tr_IPcut").Filter("ROOT::VecOps::Sum(triplet_mask4)>0")

# 9 -> dimuon mass compatible with phi(1020) RefTrack1_Eta
df = df.Define("Dimu_mass", "dimu_mass(RefTrack1_Pt, RefTrack1_Eta, RefTrack1_Phi, RefTrack2_Pt, RefTrack2_Eta, RefTrack2_Phi)")
df = df.Define("triplet_mask5", "Dimu_mass>1.0 && Dimu_mass<1.04").Filter("ROOT::VecOps::Sum(triplet_mask5)>0")

# 10 -> Trigger Matching (to-do)

# Keep best candidate based on vertex chi2
df = df.Define("BestTriplet_index", "bestcandidate(TripletVtx2_Chi2)")
df = df.Define("BestTriplet_mass", "flattening(Triplet2_Mass, BestTriplet_index)")


# Create a histogram from `x` and draw it
## h = df.Histo1D(("h_mass", "h_mass", 80, 1.65, 2.05), "BestTriplet_mass")
## c =  ROOT.TCanvas("canvas", "canvas", 1000,1000)
## h.Draw()
## c.SaveAs("myhisto.png")


#fitting invariant mass :)

x = RooRealVar("BestTriplet_mass", "2mu+1trk inv. mass (GeV)", 1.65, 2.05)
x.setBins(80)

# We first declare the RooDataHistHelper
rdhMaker = ROOT.RooDataHistHelper("dataset", "Title of dataset", ROOT.RooArgSet(x))
 
# Then, we move it into an RDataFrame action:
data_result = df.Book(ROOT.std.move(rdhMaker), ("BestTriplet_mass"))
data = roo_data_hist_result.GetValue()
entries = data.numEntries()

#data = RooDataHist("data", "h_mass", RooArgSet(x), RooFit.Import(h, ROOT.kFALSE))

x.setRange("R1", 1.70, 1.80)
x.setRange("R2", 1.89, 1.925)
x.setRange("R3", 1.99, 2.02)

meanCB = RooRealVar("mean", "meanCB", 1.97, 1.94, 2.1)
sigmaCB1 = RooRealVar("#sigma_{CB}", "sigmaCB1", 0.02, 0.001, 0.1)    
alpha1 = RooRealVar("#alpha1", "alpha1", 1.0, 0.5, 10.0)
nSigma1 = RooRealVar("n1", "n1", 1.0, 0.1, 25.0)
sig_right = RooCBShape("sig_right", "sig_right", x, meanCB, sigmaCB1, alpha1, nSigma1)

meanCB2 = RooRealVar("mean2", "meanCB2", 1.87, 1.82, 1.89)
sigmaCB2 = RooRealVar("#sigma2_{CB}", "sigmaCB2", 0.05, 0.001, 0.05)
alpha2 = RooRealVar("#alpha2", "alpha2", 1.0, 0.5, 10.0)
nSigma2 = RooRealVar("n2", "n2", 1.0, 0.1, 25.0)
sig_left = RooCBShape("sig_left", "sig_left", x, meanCB2, sigmaCB2, alpha2, nSigma2)

gamma = RooRealVar("#Gamma", "Gamma", -1, -2.0, -1e-2)
exp_bkg = RooExponential("exp_bkg", "exp_bkg", x, gamma)
exp_bkg.fitTo(data, RooFit.Range("R1,R2,R3"))

nSig_right = RooRealVar("nSig_R", "Number of signal candidates", entries*0.05, 1.0, 1e+6)
nSig_left = RooRealVar("nSig_L", "Number of signal 2 candidates", entries*0.02, 1.0, 1e+6)
nBkg = RooRealVar("nBkg", "Bkg component", entries*0.8, 1.0, 1e+6)

totalPDF = RooAddPdf("totalPDF", "totalPDF", RooArgList(sig_right, sig_left, exp_bkg), RooArgList(nSig_right, nSig_left, nBkg))

r = totalPDF.fitTo(data, RooFit.Extended(ROOT.kTRUE), RooFit.Save(ROOT.kTRUE))  

xframe = x.frame()
xframe.SetTitle("")
xframe.SetXTitle("2mu +1trk inv. mass (GeV)")
#totalPDF.paramOn(xframe, RooFit.Parameters(RooArgSet(meanCB, meanCB2, sigmaCB1, sigmaCB2, gamma, nSig_right, nSig_left, nBkg)), RooFit.Layout(0.2, 0.2, 0.6))
data.plotOn(xframe)
totalPDF.plotOn(xframe, RooFit.Components(RooArgSet(sig_right, sig_left)), RooFit.LineColor(ROOT.kRed), RooFit.LineStyle(ROOT.kDashed))
totalPDF.plotOn(xframe, RooFit.Components(RooArgSet(exp_bkg)), RooFit.LineColor(ROOT.kGreen), RooFit.LineStyle(ROOT.kDashed))
totalPDF.plotOn(xframe)

c1 = ROOT.TCanvas("c1", "c1", 900, 900)
xframe.Draw()
c1.SaveAs("myfit.png")

