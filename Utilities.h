/// Header file with functions needed to execute the Python version of
///// postselection step of the analysis. The header is declared to the
///// ROOT C++ interpreter prior to the start of the analysis via the
///// `ROOT.gInterpreter.Declare()` function.
/////

#include <iostream>
#include <vector>
#include <algorithm>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TH2F.h>
#include <string>

using namespace ROOT::VecOps;
//using RNode = ROOT::RDF::RNode;
using rvec_f = const RVec<float> &;
using rvec_d = const RVec<double> &;
using rvec_i = const RVec<int> &;
using rvec_b = const RVec<bool> &;

//values for cuts and constant

const float MU_MASS=0.10566;
const float PT_CUT_MU_ENDCAP=2;
const float PT_CUT_MU_BARREL=3.5;
const float ETA_CUT_MU=2.4;
const float DELTAR_MAX=0.8;
const float DELTAR_MIN=0.0001;
const float DELTAZ_MAX=0.8;

RVec<bool> contains(const RVec<std::string> vec, const std::string &y) {
   RVec<bool> out;
   //check if a vector element contains a substring and returns vector of booleans
   for (auto str : vec){
     out.push_back(str.find(y) != std::string::npos);
   }
   return out;
}

double deltaR(float eta1, float eta2, float phi1, float phi2){
    auto dp = std::abs(phi1 - phi2);
    auto deta = std::abs(eta1 - eta2);
    if (dp > float(M_PI))
        dp -= float(2 * M_PI);
    Double_t n = TMath::Sqrt(dp*dp + deta*deta);
    return n;
}

RVec<double> deltaR_vec(rvec_d eta1, rvec_d eta2, rvec_d phi1, rvec_d phi2){
    RVec<double> out;
    for(unsigned i=0; i<eta1.size(); ++i)
       out.push_back(deltaR(eta1.at(i), eta2.at(i), phi1.at(i), phi2.at(i)));
    return out;
}

Bool_t isPairDeltaRGood(rvec_f MuonEta, rvec_f MuonPhi, vector<int> index){
    // The function returns 'true' if all of the 3 possible pairs of muons have dR<DELTAR_CUT
    Double_t dR12 = deltaR(MuonEta[index[0]], MuonEta[index[1]], MuonPhi[index[0]], MuonPhi[index[1]]);
    Double_t dR13 = deltaR(MuonEta[index[0]], MuonEta[index[2]], MuonPhi[index[0]], MuonPhi[index[2]]);
    Double_t dR23 = deltaR(MuonEta[index[1]], MuonEta[index[2]], MuonPhi[index[1]], MuonPhi[index[2]]);
    
    if (dR12<DELTAR_MAX && dR13<DELTAR_MAX && dR23<DELTAR_MAX) return true;
    else return false;
}

Bool_t isPairDeltaZGood(double vz1, double vz2, double vz3, double DeltaZmax){
    double dZ12 = TMath::Abs(vz2 - vz1);
    double dZ13 = TMath::Abs(vz3 - vz1);
    double dZ23 = TMath::Abs(vz3 - vz2);
    
    if (dZ12<DeltaZmax && dZ13<DeltaZmax && dZ23<DeltaZmax) return true;
    else return false;
}

double Cos3D_(double TripletVtx_x, double TripletVtx_y, double TripletVtx_z, double RefittedPV_x, double RefittedPV_y, double RefittedPV_z, double Triplet_Pt, double Triplet_Eta,double Triplet_Phi){
    // Computes the angle between the momentum vector of the candidatet (b) and the vector from the primary vertex (a)
    double a_x = TripletVtx_x - RefittedPV_x;
    double a_y = TripletVtx_y - RefittedPV_y;
    double a_z = TripletVtx_z - RefittedPV_z;
    TVector3 b;
    b.SetPtEtaPhi(Triplet_Pt, Triplet_Eta, Triplet_Phi);
    double b_x = b.Px();
    double b_y = b.Py();
    double b_z = b.Pz();
    double a_mod = sqrt(a_x*a_x + a_y*a_y + a_z*a_z);
    double b_mod = abs(b.Mag());
    double cos_ang = ((a_x*b_x)+(a_y*b_y)+(a_z*b_z))/(a_mod*b_mod);
    return cos_ang;
}

double Cos2D_(double TripletVtx_x, double TripletVtx_y, double RefittedPV_x, double RefittedPV_y, double Triplet_Pt, double Triplet_Eta,double Triplet_Phi){
    // Computes the angle between the momentum vector of the 4mu quadruplet (b) and the vector from the primary vertex (a)
    double a_x = TripletVtx_x - RefittedPV_x;
    double a_y = TripletVtx_y - RefittedPV_y;
    TVector3 b;
    b.SetPtEtaPhi(Triplet_Pt, Triplet_Eta, Triplet_Phi);
    double b_x = b.Px();
    double b_y = b.Py();
    double a_mod = sqrt(a_x*a_x + a_y*a_y);
    double b_mod = sqrt(b_x*b_x + b_y*b_y);
    double cos_ang = ((a_x*b_x)+(a_y*b_y))/(a_mod*b_mod);
    return cos_ang;
}

vector<int> get_3index(ROOT::VecOps::RVec<float> MuonPt, double pt1, double pt2, double pt3){
    vector<int> index;
    int i=0;
    auto i1 = std::find(MuonPt.begin(), MuonPt.end(), pt1);
    auto i2 = std::find(MuonPt.begin(), MuonPt.end(), pt2);
    auto i3 = std::find(MuonPt.begin(), MuonPt.end(), pt3);
    if (i1 != MuonPt.end() && i2 != MuonPt.end() && i3 != MuonPt.end()) {
        index.push_back(std::distance(MuonPt.begin(), i1));
        index.push_back(std::distance(MuonPt.begin(), i2));
        index.push_back(std::distance(MuonPt.begin(), i3));
    }
    else{
        index.push_back(-1);
        index.push_back(-1);
        index.push_back(-1);
    }
    return index;
}

vector<int> get_2index(ROOT::VecOps::RVec<float> MuonPt, double pt1, double pt2){
    vector<int> index;
    int i=0;
    auto i1 = std::find(MuonPt.begin(), MuonPt.end(), pt1);
    auto i2 = std::find(MuonPt.begin(), MuonPt.end(), pt2);
    if (i1 != MuonPt.end() && i2 != MuonPt.end()) {
        index.push_back(std::distance(MuonPt.begin(), i1));
        index.push_back(std::distance(MuonPt.begin(), i2));
    }
    else{
        index.push_back(-1);
        index.push_back(-1);
    }
    return index;
}

RVec<int> match(ROOT::VecOps::RVec<double> branch1, ROOT::VecOps::RVec<double> branch2){
//returns vector of indeces such that branch2[index]=branch1
    RVec<int> index;
    for(unsigned i = 0; i<branch1.size(); i++){
      auto idx = std::find(branch2.begin(), branch2.end(), branch1.at(i));
      if( idx != branch2.end()) index.push_back(std::distance(branch2.begin(), idx));
      else index.push_back(-99);
    }
    return index;
}

RVec<bool> muon_id(RVec<int> index, RVec<bool> branch2){
//returns booleans from branch2, aligned to branch1 based on index
    RVec<bool> out;
    for(unsigned i = 0; i<index.size(); i++){
       out.push_back(branch2.at(i));
    }
    return out;
}

RVec<double> dimu_mass(RVec<double> pt1, RVec<double> eta1, RVec<double> phi1, RVec<double> pt2, RVec<double> eta2, RVec<double> phi2){
   TLorentzVector mu1, mu2, dimu;
   RVec<double> out;
   for(unsigned i = 0; i<pt1.size(); i++){
      mu1.SetPtEtaPhiM(pt1.at(i), eta1.at(i), phi1.at(i), MU_MASS);
      mu2.SetPtEtaPhiM(pt2.at(i), eta2.at(i), phi2.at(i), MU_MASS);
      dimu = mu1 + mu2;
      out.push_back(dimu.M());
   }
   return out;
}

int bestcandidate(RVec<double> TripletChi2){
//returns index pointing to best candidate (one per event) based on Vtx fit Chi2
  auto ptr = std::min_element(TripletChi2.begin(), TripletChi2.end());
  return std::distance(TripletChi2.begin(), ptr);
}

double flattening(ROOT::VecOps::RVec<double> var, int index){
    double value = -99;
    try {
        value = var.at(index);
    } catch (const std::out_of_range& e) {
        std::cout << "Not valid index " << std::endl;
        return -99;
    }
    return value;
}

