/// Header file with functions needed to execute the Python version of
/// postselection step of the analysis. The header is declared to the
/// ROOT C++ interpreter prior to the start of the analysis via the
/// `ROOT.gInterpreter.Declare()` function.
///

#ifndef POST_H
#define POST_H

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TFile.h"
#include "TH2D.h"
#include "TLatex.h"
#include "Math/Vector4D.h"
#include "TStyle.h"
#include <map>

#include "TDavixFile.h"

using namespace ROOT::VecOps;
using RNode = ROOT::RDF::RNode;
using rvec_f = const RVec<float> &;
using rvec_i = const RVec<int> &;
using rvec_b = const RVec<bool> &;

//values for cuts and constant 

const size_t ONLYELE=1;
const size_t ONLYMU=0;

const float PT_CUT_MU=  30;
const float ETA_CUT_MU= 2.4;
const float ISO_CUT_MU= 0.15;

const float PT_CUT_ELE_UL2016=  30;
const float PT_CUT_ELE_UL2017=  38;
const float PT_CUT_ELE_UL2018=  35;

const float ETA_CUT_ELE= 2.4;
const float ISO_CUT_ELE= 0.08;

const float REL_ISO_CUT_LEP_VETO_ELE=   0.2;
const float PT_CUT_LEP_VETO_ELE=        15;
const float ETA_CUT_LEP_VETO_ELE=       2.4;

const float REL_ISO_CUT_LEP_VETO_MU=    0.4;
const float PT_CUT_LEP_VETO_MU=         10;
const float ETA_CUT_LEP_VETO_MU=        2.4;

const float DR_OVERLAP_CONE_TAU=        0.5;
const float DR_OVERLAP_CONE_OTHER=      0.4;

const float PT_CUT_JET= 30;
const float ETA_CUT_JET=5;

const float DELTAETA_JJ_CUT=2.5;

const float BTAG_PT_CUT =   30;
const float BTAG_ETA_CUT=   5;
string BTAG_ALGO   =   "DeepFlv";
string BTAG_WP     =   "M";
const float BTAG_WP_VALUE = 0.3033;
const float BTAG_WP_VALUE_LOOSE = 0.0521;

const size_t ID_TAU_RECO_DEEPTAU_VSJET_LOOSE_ELE_UL2016APV = 8;
const size_t ID_TAU_RECO_DEEPTAU_VSJET_VETO_ELE_UL2016APV = 8;
const size_t ID_TAU_RECO_DEEPTAU_VSJET_LOOSE_MU_UL2016APV = 4;
const size_t ID_TAU_RECO_DEEPTAU_VSJET_VETO_MU_UL2016APV = 4;
const size_t ID_TAU_RECO_DEEPTAU_VSJET_LOOSE_ELE_UL2016 = 4 ;
const size_t ID_TAU_RECO_DEEPTAU_VSJET_VETO_ELE_UL2016 = 4;
const size_t ID_TAU_RECO_DEEPTAU_VSJET_LOOSE_MU_UL2016 = 8;
const size_t ID_TAU_RECO_DEEPTAU_VSJET_VETO_MU_UL2016 = 8;
const size_t ID_TAU_RECO_DEEPTAU_VSJET_LOOSE_ELE_UL2017 = 16;
const size_t ID_TAU_RECO_DEEPTAU_VSJET_VETO_ELE_UL2017 = 16;
const size_t ID_TAU_RECO_DEEPTAU_VSJET_LOOSE_MU_UL2017 = 8;
const size_t ID_TAU_RECO_DEEPTAU_VSJET_VETO_MU_UL2017 = 8;
const size_t ID_TAU_RECO_DEEPTAU_VSJET_LOOSE_ELE_UL2018 = 16;
const size_t ID_TAU_RECO_DEEPTAU_VSJET_VETO_ELE_UL2018 = 16;
const size_t ID_TAU_RECO_DEEPTAU_VSJET_LOOSE_MU_UL2018 = 16;
const size_t ID_TAU_RECO_DEEPTAU_VSJET_VETO_MU_UL2018 = 16;

const size_t ID_TAU_RECO_DEEPTAU_VSJET=  64;
const size_t ID_TAU_RECO_DEEPTAU_VSELE=  4; 
const size_t ID_TAU_RECO_DEEPTAU_VSMU=   8; 

const size_t ID_TAU_RECO_MVA=            8; 
const size_t ID_TAU_ANTIMU=              1; 
const size_t ID_TAU_ANTIELE=             2; 
const float PT_CUT_TAU=30;
const float ETA_CUT_TAU=2.3;
const float M_JJ_CUT=   500;
const float MET_CUT=    40;

const string vsJetwp = "VTight";
const string vsElewp = "VLoose";
const string vsMuwp = "Tight";

const string remote_storage = "https://vbs-pg-support.web.cern.ch/nanoAOD-tools/python/postprocessing/";



TFile *TauID_SF_pt_DeepTau2017v2p1VSjet_UL2017 = TFile::Open(TString(remote_storage) + TString("data/tauSF/TauID_SF_pt_") + TString("DeepTau2017v2p1VSjet") + TString("_") + TString("UL2017") + TString(".root"));
TString path_down = TString(TString(vsJetwp) + TString("_down"));
TString path_cent = TString(TString(vsJetwp) + TString("_cent"));
TString path_up = TString(TString(vsJetwp) + TString("_up"));
TF1 * TauID_SF_pt_DeepTau2017v2p1VSjet_UL2017_h_down = (TF1*)TauID_SF_pt_DeepTau2017v2p1VSjet_UL2017->Get(path_down);
TF1 * TauID_SF_pt_DeepTau2017v2p1VSjet_UL2017_h_cent =  (TF1*)TauID_SF_pt_DeepTau2017v2p1VSjet_UL2017->Get(path_cent);
TF1 * TauID_SF_pt_DeepTau2017v2p1VSjet_UL2017_h_up =  (TF1*)TauID_SF_pt_DeepTau2017v2p1VSjet_UL2017->Get(path_up);

TFile *TauID_SF_eta_DeepTau2017v2p1VSe_UL2017 = TFile::Open(TString(remote_storage) + TString("data/tauSF/TauID_SF_eta_") + TString("DeepTau2017v2p1VSe") + TString("_") + TString("UL2017") + TString(".root"));
TString histoname_ele = TString(vsElewp);
TH1F * TauID_SF_eta_DeepTau2017v2p1VSe_UL2017_hist = (TH1F *) TauID_SF_eta_DeepTau2017v2p1VSe_UL2017->Get(histoname_ele);

TFile *TauID_SF_eta_DeepTau2017v2p1VSmu_UL2017 = TFile::Open(TString(remote_storage) + TString("data/tauSF/TauID_SF_eta_") + TString("DeepTau2017v2p1VSmu") + TString("_") + TString("2017ReReco") + TString(".root"));
TString histoname_mu = TString(vsMuwp);
TH1F * TauID_SF_eta_DeepTau2017v2p1VSmu_UL2017_hist = (TH1F *) TauID_SF_eta_DeepTau2017v2p1VSmu_UL2017->Get(histoname_mu);

TFile *TauES_dm_DeepTau2017v2p1VSjet_UL2017 = TFile::Open(TString(remote_storage) + TString("data/tauSF/TauES_dm_") + TString("DeepTau2017v2p1VSjet") + TString("_") + TString("UL2017") +  TString(".root"));        
TH1F * TauES_dm_DeepTau2017v2p1VSjet_UL2017_hist_low = (TH1F *) TauES_dm_DeepTau2017v2p1VSjet_UL2017->Get("tes");

TFile *TauES_dm_DeepTau2017v2p1VSjet_UL2017_ptgt100 = TFile::Open(TString(remote_storage) + TString("data/tauSF/TauES_dm_") + TString("DeepTau2017v2p1VSjet") + TString("_") + TString("2017ReReco") + TString("_ptgt100.root"));
TH1F * TauES_dm_DeepTau2017v2p1VSjet_UL2017_ptgt100_hist_high = (TH1F *) TauES_dm_DeepTau2017v2p1VSjet_UL2017_ptgt100->Get("tes");

TFile *TauFES_eta_dm_DeepTau2017v2p1VSe_UL2017 = TFile::Open(TString(remote_storage) + TString("data/tauSF/TauFES_eta-dm_") + TString("DeepTau2017v2p1VSe") + TString("_") + TString("2017")  + TString("ReReco") + TString(".root"));
TGraphAsymmErrors * TauFES_eta_dm_DeepTau2017v2p1VSe_UL2017_graph = (TGraphAsymmErrors *) TauFES_eta_dm_DeepTau2017v2p1VSe_UL2017->Get("fes");

float deltaPhi (float phi1, float phi2){
    float dphi = (phi1-phi2);
    while(dphi >  M_PI) dphi -= 2*M_PI;
    while(dphi < -M_PI) dphi += 2*M_PI;
    return dphi;
}

float deltaR(float eta1, float phi1, float eta2, float phi2){
    return hypot(eta1 - eta2, deltaPhi(phi1, phi2)); 
}

RVec<float> getFlattenedMatrixColumn(rvec_f flattened_matrix, int nColumns, int column_index){
    RVec<float> result;
    for (int i = 0; i < flattened_matrix.size()/nColumns; i++) result.emplace_back(flattened_matrix[column_index + i*nColumns]);
    return result;
}


bool LepVetoEle(rvec_i Electron_idx, rvec_f Electron_pt, rvec_f Electron_eta, rvec_b Iso_WPL, rvec_f jetRelIso, rvec_f Muon_pt, rvec_f Muon_eta, rvec_f Iso04_all, rvec_b Muon_looseId )
{
    bool IsEleVetoPassed = true;
    bool IsMuVetoPassed = true;
    for (size_t i = 0; i < Electron_pt.size(); i++) {
        if(i != Electron_idx[0] && Iso_WPL[i] && Electron_pt[i] > PT_CUT_LEP_VETO_ELE && abs(Electron_eta[i]) < ETA_CUT_LEP_VETO_ELE && !(abs(Electron_eta[i])>1.4442 && abs(Electron_eta[i])<1.566) && jetRelIso[i] < REL_ISO_CUT_LEP_VETO_ELE) IsEleVetoPassed = false;
    }
    if(IsEleVetoPassed == true){
        for (size_t i = 0; i < Muon_pt.size(); i++) {
            if(Muon_looseId[i] && Muon_pt[i] > PT_CUT_LEP_VETO_MU && abs(Muon_eta[i]) < ETA_CUT_LEP_VETO_MU && Iso04_all[i] < REL_ISO_CUT_LEP_VETO_MU) IsMuVetoPassed = false;
        }
    }
    return IsEleVetoPassed && IsMuVetoPassed;
}



bool LepVetoMu(rvec_i Muon_idx, rvec_f Electron_pt, rvec_f Electron_eta, rvec_b Electron_mvaFall17V2Iso_WPL, rvec_f Electron_jetRelIso, rvec_f Muon_pt, rvec_f Muon_eta, rvec_f Muon_pfRelIso04_all, rvec_b Muon_looseId)
{
    bool IsEleVetoPassed = true;
    bool IsMuVetoPassed = true;
    for (size_t i = 0; i < Electron_pt.size(); i++) {
        if(Electron_mvaFall17V2Iso_WPL[i] && Electron_pt[i] > PT_CUT_LEP_VETO_ELE && abs(Electron_eta[i]) < ETA_CUT_LEP_VETO_ELE && !(abs(Electron_eta[i])>1.4442 && abs(Electron_eta[i])<1.566) && Electron_jetRelIso[i] < REL_ISO_CUT_LEP_VETO_ELE) IsEleVetoPassed = false;
    }
    if(IsEleVetoPassed == true){
        for (size_t i = 0; i < Muon_pt.size(); i++) {
            if(i != Muon_idx[0] && Muon_looseId[i] && Muon_pt[i] > PT_CUT_LEP_VETO_MU && abs(Muon_eta[i]) < ETA_CUT_LEP_VETO_MU && Muon_pfRelIso04_all[i] < REL_ISO_CUT_LEP_VETO_MU) IsMuVetoPassed = false;
        }
    }
    return IsEleVetoPassed && IsMuVetoPassed;
}

RVec<size_t> GoodJets(rvec_i jetId, rvec_f eta, rvec_f pt, rvec_i puId){
   RVec<int> idx;
   for (size_t i = 0; i < pt.size(); i++) {
      if (jetId[i] >= 2 && abs(eta[i]) < 5. && pt[i] > PT_CUT_JET && (pt[i] > 50. || (pt[i] <= 50. && puId[i] >= 7))) idx.emplace_back(i);
   }
   return idx;
}

bool atleast2GoodJets(rvec_i GoodJets_idx){
    if (GoodJets_idx.size() >= 2) return true;
    else return false;
}


RVec<size_t> SelectVBSJets_invmass(rvec_f pt, rvec_f eta, rvec_f phi, rvec_f mass, rvec_i GoodJets_idx)
{
    
    RVec<size_t> idx;
    // Find first lepton pair with invariant mass closest to Z mass
    auto idx_cmb = Combinations(GoodJets_idx, 2);
    //auto best_mass = -1.;
    float best_mass = -1.;
    size_t best_i1 = 0; size_t best_i2 = 0;
    for (size_t i = 0; i < idx_cmb[0].size(); i++) {
        const auto i1 = idx_cmb[0][i];
        const auto i2 = idx_cmb[1][i];
        if (abs(eta[GoodJets_idx[i1]] - eta[GoodJets_idx[i2]]) >= DELTAETA_JJ_CUT) {
            ROOT::Math::PtEtaPhiMVector p1(pt[GoodJets_idx[i1]], eta[GoodJets_idx[i1]], phi[GoodJets_idx[i1]], mass[GoodJets_idx[i1]]);
            ROOT::Math::PtEtaPhiMVector p2(pt[GoodJets_idx[i2]], eta[GoodJets_idx[i2]], phi[GoodJets_idx[i2]], mass[GoodJets_idx[i2]]);
            float this_mass = (p1 + p2).M();
            if (this_mass > best_mass) {
                best_mass = this_mass;
                best_i1 = GoodJets_idx[i1];
                best_i2 = GoodJets_idx[i2];
            }
        }
    } 
    idx.emplace_back(best_i1);
    idx.emplace_back(best_i2);
    return idx;
}

float GetInvMass(rvec_f pt, rvec_f eta, rvec_f phi, rvec_f mass, rvec_i VBSJets_idx)
{
    ROOT::Math::PtEtaPhiMVector p1(pt[VBSJets_idx[0]], eta[VBSJets_idx[0]], phi[VBSJets_idx[0]], mass[VBSJets_idx[0]]);
    ROOT::Math::PtEtaPhiMVector p2(pt[VBSJets_idx[1]], eta[VBSJets_idx[1]], phi[VBSJets_idx[1]], mass[VBSJets_idx[1]]);
    return (p1 + p2).M();
}

float GetInvMassNoIndex(float pt1, float eta1, float phi1, float mass1, float pt2, float eta2, float phi2, float mass2)
{
    ROOT::Math::PtEtaPhiMVector p1(pt1, eta1, phi1, mass1);
    ROOT::Math::PtEtaPhiMVector p2(pt2, eta2, phi2, mass2);
    return (p1 + p2).M();
}

float GetInvMassNoIndex3(float pt1, float eta1, float phi1, float mass1, float pt2, float eta2, float phi2, float mass2, float pt3, float eta3, float phi3, float mass3)
{
    ROOT::Math::PtEtaPhiMVector p1(pt1, eta1, phi1, mass1);
    ROOT::Math::PtEtaPhiMVector p2(pt2, eta2, phi2, mass2);
    ROOT::Math::PtEtaPhiMVector p3(pt3, eta3, phi3, mass3);
    return (p1 + p2 + p3).M();
}

float GetInvMassNoIndex4(float pt1, float eta1, float phi1, float mass1, float pt2, float eta2, float phi2, float mass2, float pt3, float eta3, float phi3, float mass3, float pt4, float eta4, float phi4, float mass4)
{
    ROOT::Math::PtEtaPhiMVector p1(pt1, eta1, phi1, mass1);
    ROOT::Math::PtEtaPhiMVector p2(pt2, eta2, phi2, mass2);
    ROOT::Math::PtEtaPhiMVector p3(pt3, eta3, phi3, mass3);
    ROOT::Math::PtEtaPhiMVector p4(pt4, eta4, phi4, mass4);
    return (p1 + p2 + p3 + p4).M();
}

float GetLeading(rvec_f Jet_pt, rvec_i VBSJet_idx){
    return Jet_pt[VBSJet_idx[0]];
}

float GetSubLeading(rvec_f Jet_pt, rvec_i VBSJet_idx){
    return Jet_pt[VBSJet_idx[1]];
}

float GetLepton(rvec_f Electron_pt, rvec_i Electron_idx, rvec_f Muon_pt, rvec_i Muon_idx, int GoodLeptonFamily){
    if (GoodLeptonFamily == 0) return Electron_pt[Electron_idx[0]];
    else return Muon_pt[Muon_idx[0]];
}

int GetLeptonTightFlag(rvec_i Electron_idx, rvec_i Muon_idx, int GoodLeptonFamily){
    if (GoodLeptonFamily == 0) return Electron_idx[1];
    else return Muon_idx[1];
}

float GetTau(rvec_f pt, rvec_i idx){
    return pt[idx[0]];
}

RVec<int> SelectElectron(rvec_f lepton_pt, rvec_f lepton_eta, rvec_f lepton_phi, rvec_f lepton_jetRelIso, rvec_b lepton_mvaFall17V2Iso_WPL, rvec_f lepton_mvaFall17V2Iso_WP90, rvec_f jet_eta, rvec_f jet_phi, rvec_i VBSJets_idx, string Year){
    
    float PT_CUT_ELE;
    
    if(Year == "UL2016" || Year == "UL2016APV") PT_CUT_ELE = PT_CUT_ELE_UL2016;
    else if (Year == "UL2017") PT_CUT_ELE = PT_CUT_ELE_UL2017;
    else PT_CUT_ELE = PT_CUT_ELE_UL2018;
    
    
    float jet1_idx = VBSJets_idx[0];
    float jet2_idx = VBSJets_idx[1];
    //const auto jet1eta = lepton_pt[0]
    float jet1eta = jet_eta[jet1_idx];
    float jet2eta = jet_eta[jet2_idx];
    float jet1phi = jet_phi[jet1_idx];
    float jet2phi = jet_phi[jet2_idx];
    
    const float isocone = DR_OVERLAP_CONE_OTHER;

    RVec<size_t> Tleptons_idx;
    RVec<size_t> LnTleptons_idx;
    
    bool IsLooseID, IsTightID, IsTightIso, IsLooseIso, IsInEtaRegion, IsInPtRegion;
    
    for (size_t i = 0; i < lepton_pt.size(); i++) {
        //setting loose and tight, eta, and pt criteria for leptons depending on lepton flavour
        IsLooseID = lepton_mvaFall17V2Iso_WPL[i];
        IsTightID = lepton_mvaFall17V2Iso_WP90[i];
        IsTightIso = lepton_jetRelIso[i]<ISO_CUT_ELE && lepton_jetRelIso[i]>=0.;
        IsLooseIso = lepton_jetRelIso[i]<1. && lepton_jetRelIso[i]>=0.;
        IsInEtaRegion = abs(lepton_eta[i])<ETA_CUT_ELE && !(abs(lepton_eta[i])>1.4442 && abs(lepton_eta[i])<1.566);
        IsInPtRegion = lepton_pt[i] > PT_CUT_ELE;

        //find tight and loose-not-tight leptons filtering with jet-lep isolation criteria
        if (IsInEtaRegion && IsInPtRegion){
            if(IsLooseID && IsLooseIso){
                if(deltaR(lepton_eta[i], lepton_phi[i], jet1eta,  jet1phi) > isocone && deltaR(lepton_eta[i], lepton_phi[i], jet2eta,  jet2phi) > isocone){
                    if(IsTightID && IsTightIso) Tleptons_idx.emplace_back(i);
                    else LnTleptons_idx.emplace_back(i);
                }
            }
        }
    }
 
    RVec<int> idx(2);
    //select leading tight/loose-not-tight lepton
    if (Tleptons_idx.size() > 0){
        idx[0] = Tleptons_idx[0];
        idx[1] = 1;
    }
    else if (LnTleptons_idx.size() > 0){
        idx[0] = LnTleptons_idx[0];
        idx[1] = 0;
    }
    else{
        idx[0] = -1;
        idx[1] = -1;  
    }
    return idx;
}


RVec<int> SelectMuon(rvec_f lepton_pt, rvec_f lepton_eta, rvec_f lepton_phi, rvec_b lepton_tightId, rvec_b lepton_looseId, rvec_f Iso04_all, rvec_f jet_eta, rvec_f jet_phi, rvec_i VBSJets_idx){
    
    float jet1_idx = VBSJets_idx[0];
    float jet2_idx = VBSJets_idx[1];
    float jet1eta = jet_eta[jet1_idx];
    float jet2eta = jet_eta[jet2_idx];
    float jet1phi = jet_phi[jet1_idx];
    float jet2phi = jet_phi[jet2_idx];
    
    const float isocone = DR_OVERLAP_CONE_OTHER;

    RVec<size_t> Tleptons_idx;
    RVec<size_t> LnTleptons_idx;
    
    bool IsLooseID, IsTightID, IsTightIso, IsLooseIso, IsInEtaRegion, IsInPtRegion;
    
    for (size_t i = 0; i < lepton_pt.size(); i++) {
        //setting loose and tight, eta, and pt criteria for leptons depending on lepton flavour
        IsTightID = lepton_tightId[i];
        IsLooseID = lepton_looseId[i];
        IsTightIso = Iso04_all[i]<ISO_CUT_MU && Iso04_all[i]>=0.;
        IsLooseIso = Iso04_all[i]<1. && Iso04_all[i]>=0.;
        IsInEtaRegion = abs(lepton_eta[i]) < ETA_CUT_MU;
        IsInPtRegion = lepton_pt[i] > PT_CUT_MU;
        

        //find tight and loose-not-tight leptons filtering with jet-lep isolation criteria
        if (IsInEtaRegion && IsInPtRegion){
            if(IsLooseID && IsLooseIso){
                if(deltaR(lepton_eta[i], lepton_phi[i], jet1eta,  jet1phi) > isocone && deltaR(lepton_eta[i], lepton_phi[i], jet2eta,  jet2phi) > isocone){
                    if(IsTightID && IsTightIso) Tleptons_idx.emplace_back(i);
                    else LnTleptons_idx.emplace_back(i);
                }
            }
        }
    }
    RVec<int> idx(2);
    //select leading tight/loose-not-tight lepton
    if (Tleptons_idx.size() > 0){
        idx[0] = Tleptons_idx[0];
        idx[1] = 1;
    }
    else if (LnTleptons_idx.size() > 0){
        idx[0] = LnTleptons_idx[0];
        idx[1] = 0;
    }
    else{
        idx[0] = -1;
        idx[1] = -1;  
    }
    
    return idx;
}


int DetermineGoodLepton(bool HLT_IsoMu27, bool HLT_Mu50, bool HLT_Ele35_WPTight_Gsf, bool HLT_Ele32_WPTight_Gsf_L1DoubleEG, bool HLT_Photon200, bool HLT_PFHT250, bool HLT_PFHT350, rvec_i Electron_idx, rvec_f Electron_pt, rvec_f Electron_eta, rvec_b Electron_mvaFall17V2Iso_WPL, rvec_f Electron_jetRelIso, rvec_i Muon_idx, rvec_f Muon_pt, rvec_f Muon_eta, rvec_f Muon_pfRelIso04_all, rvec_b Muon_looseId){
    bool passMu = false;
    bool passEle = false;
    bool passHT = false;
    int GoodLeptonFamily;
    
    if(HLT_IsoMu27) passMu = true;
    if(HLT_Ele35_WPTight_Gsf) passEle = true;
    
    bool ele_lepton_veto = false;
    bool mu_lepton_veto = false;
    
    if(Electron_idx[1] != -1) ele_lepton_veto = LepVetoEle(Electron_idx, Electron_pt, Electron_eta, Electron_mvaFall17V2Iso_WPL, Electron_jetRelIso, Muon_pt, Muon_eta, Muon_pfRelIso04_all, Muon_looseId);
    if(Muon_idx[1] != -1) mu_lepton_veto = LepVetoMu(Muon_idx, Electron_pt, Electron_eta, Electron_mvaFall17V2Iso_WPL, Electron_jetRelIso, Muon_pt, Muon_eta, Muon_pfRelIso04_all, Muon_looseId);

    bool SingleEle=false;
    bool SingleMu=false;
    if(passEle && !passMu){
        if(Electron_idx[1] != -1 && ele_lepton_veto){
            GoodLeptonFamily = 0;
            //lepton_TightRegion[0] = copy.deepcopy(ele_TightRegion)
            SingleEle = true;
            SingleMu = false;
        }
        else GoodLeptonFamily = -1;
    }
    else if(!passEle && passMu){
        if(Muon_idx[1] != -1 && mu_lepton_veto){
            GoodLeptonFamily = 1;
            //lepton_TightRegion[0] = copy.deepcopy(mu_TightRegion)
            SingleEle = false;
            SingleMu = true;
        }
        else GoodLeptonFamily = -1;
    }
    else if(passEle && passMu){ 
        if(Muon_idx[1] == -1 && Electron_idx[1] != -1 && ele_lepton_veto){
            GoodLeptonFamily = 0;
            //lepton_TightRegion[0] = copy.deepcopy(ele_TightRegion)
            SingleEle = true;
            SingleMu = false;
        }
        else if(Electron_idx[1] == -1 && Muon_idx[1] != -1 && mu_lepton_veto){
            GoodLeptonFamily = 1;
            //lepton_TightRegion[0] = copy.deepcopy(mu_TightRegion)
            SingleEle = false;
            SingleMu = true;
        }
                
        else if(Muon_idx[1] != -1 && Electron_idx[1] != -1){
            if(ele_lepton_veto && !mu_lepton_veto){
                GoodLeptonFamily = 0;
                //lepton_TightRegion[0] = copy.deepcopy(ele_TightRegion)
                SingleEle = false;
                SingleMu = true;
            }
            else if(!ele_lepton_veto && mu_lepton_veto){           
                GoodLeptonFamily = 1;
                //lepton_TightRegion[0] = copy.deepcopy(mu_TightRegion)
                SingleEle = false;
                SingleMu = true;
            }

            else if(ele_lepton_veto && mu_lepton_veto){
                if(Electron_pt[Electron_idx[0]] > Muon_pt[Muon_idx[0]]){
                    GoodLeptonFamily = 0;
                    //lepton_TightRegion[0] = copy.deepcopy(ele_TightRegion)
                    SingleEle = true;
                    SingleMu = false;
                }
                else{
                    GoodLeptonFamily = 1;
                    //lepton_TightRegion[0] = copy.deepcopy(mu_TightRegion)
                    SingleEle = false;
                    SingleMu = true;
                }
            }
            else GoodLeptonFamily = -1;
        }
        else GoodLeptonFamily = -1;
    }
    if(!(SingleEle || SingleMu)) GoodLeptonFamily = -1;
    
    return GoodLeptonFamily;
}


RVec<int> SelectAndVetoTaus(rvec_f Tau_pt, rvec_f Tau_eta, rvec_f Tau_phi, RVec<UChar_t> &Tau_idDeepTau2017v2p1VSjet, RVec<UChar_t> &Tau_idDeepTau2017v2p1VSe, RVec<UChar_t> &Tau_idDeepTau2017v2p1VSmu, int GoodLeptonFamily, rvec_i Electron_idx, rvec_f Electron_eta, rvec_f Electron_phi, rvec_i Muon_idx, rvec_f Muon_eta, rvec_f Muon_phi, rvec_f Jet_eta, rvec_f Jet_phi, rvec_i VBSJet_idx, string Year)
{
    size_t ID_TAU_RECO_DEEPTAU_VSJET_VETO_ELE,  ID_TAU_RECO_DEEPTAU_VSJET_VETO_MU, ID_TAU_RECO_DEEPTAU_VSJET_LOOSE_ELE, ID_TAU_RECO_DEEPTAU_VSJET_LOOSE_MU;
    
    if(Year == "UL2016" || Year == "UL2016APV") ID_TAU_RECO_DEEPTAU_VSJET_VETO_ELE = ID_TAU_RECO_DEEPTAU_VSJET_VETO_ELE_UL2016;
    else if (Year == "UL2017") ID_TAU_RECO_DEEPTAU_VSJET_VETO_ELE = ID_TAU_RECO_DEEPTAU_VSJET_VETO_ELE_UL2017;
    else ID_TAU_RECO_DEEPTAU_VSJET_VETO_ELE = ID_TAU_RECO_DEEPTAU_VSJET_VETO_ELE_UL2018;
    
    if(Year == "UL2016" || Year == "UL2016APV") ID_TAU_RECO_DEEPTAU_VSJET_VETO_MU = ID_TAU_RECO_DEEPTAU_VSJET_VETO_MU_UL2016;
    else if (Year == "UL2017") ID_TAU_RECO_DEEPTAU_VSJET_VETO_MU = ID_TAU_RECO_DEEPTAU_VSJET_VETO_MU_UL2017;
    else ID_TAU_RECO_DEEPTAU_VSJET_VETO_MU = ID_TAU_RECO_DEEPTAU_VSJET_VETO_MU_UL2018;
    
    if(Year == "UL2016" || Year == "UL2016APV") ID_TAU_RECO_DEEPTAU_VSJET_LOOSE_ELE = ID_TAU_RECO_DEEPTAU_VSJET_LOOSE_ELE_UL2016;
    else if (Year == "UL2017") ID_TAU_RECO_DEEPTAU_VSJET_LOOSE_ELE = ID_TAU_RECO_DEEPTAU_VSJET_LOOSE_ELE_UL2017;
    else ID_TAU_RECO_DEEPTAU_VSJET_LOOSE_ELE = ID_TAU_RECO_DEEPTAU_VSJET_LOOSE_ELE_UL2018;
    
    if(Year == "UL2016" || Year == "UL2016APV") ID_TAU_RECO_DEEPTAU_VSJET_LOOSE_MU = ID_TAU_RECO_DEEPTAU_VSJET_LOOSE_MU_UL2016;
    else if (Year == "UL2017") ID_TAU_RECO_DEEPTAU_VSJET_LOOSE_MU = ID_TAU_RECO_DEEPTAU_VSJET_LOOSE_MU_UL2017;
    else ID_TAU_RECO_DEEPTAU_VSJET_LOOSE_MU = ID_TAU_RECO_DEEPTAU_VSJET_LOOSE_MU_UL2018;
    
    //setting jet-related quantities if isolation from them is needed
    
    float jet1eta = Jet_eta[VBSJet_idx[0]];
    float jet2eta = Jet_eta[VBSJet_idx[1]];
    float jet1phi = Jet_phi[VBSJet_idx[0]];
    float jet2phi = Jet_phi[VBSJet_idx[1]];
    float isocone = DR_OVERLAP_CONE_OTHER;

    size_t nTau=0;
    RVec<int> idx(2);
    
    if (Tau_eta.size()==0){
        idx[0] = -1;
        idx[1] = -1;
        return idx;
    } 
    float cutloose_vsjet;
    bool isAtLeastLoose = false;
    for (size_t i = 0; i < Tau_eta.size(); i++) {
        
        if (GoodLeptonFamily == 0) {
            cutloose_vsjet = ID_TAU_RECO_DEEPTAU_VSJET_VETO_ELE;
            if ((Tau_idDeepTau2017v2p1VSjet[i]>=cutloose_vsjet && Tau_idDeepTau2017v2p1VSe[i]>=ID_TAU_RECO_DEEPTAU_VSELE && Tau_idDeepTau2017v2p1VSmu[i]>=ID_TAU_RECO_DEEPTAU_VSMU) && deltaR(Tau_eta[i], Tau_phi[i], Electron_eta[Electron_idx[0]], Electron_phi[Electron_idx[0]])>DR_OVERLAP_CONE_TAU && deltaR(Tau_eta[i], Tau_phi[i], jet1eta, jet1phi)>isocone && deltaR(Tau_eta[i], Tau_phi[i], jet2eta, jet2phi)>isocone && Tau_pt[i] >=PT_CUT_TAU && abs(Tau_eta[i])<=ETA_CUT_TAU){
                nTau++;
                
                if(Tau_idDeepTau2017v2p1VSjet[i]>=ID_TAU_RECO_DEEPTAU_VSJET){
                    idx[0] = i;
                    //idx[1] = 0;
                    idx[1] = 1;
                    isAtLeastLoose = true;
                } 
                else{
                    if(Tau_idDeepTau2017v2p1VSjet[i]>=ID_TAU_RECO_DEEPTAU_VSJET_LOOSE_ELE){
                        idx[0] = i;
                        //idx[1] = 1;
                        idx[1] = 0;
                        isAtLeastLoose = true;
                    }
                }
            }   
        }

        else if (GoodLeptonFamily == 1) {
            cutloose_vsjet = ID_TAU_RECO_DEEPTAU_VSJET_VETO_MU;
            if ((Tau_idDeepTau2017v2p1VSjet[i]>=cutloose_vsjet && Tau_idDeepTau2017v2p1VSe[i]>=ID_TAU_RECO_DEEPTAU_VSELE && Tau_idDeepTau2017v2p1VSmu[i]>=ID_TAU_RECO_DEEPTAU_VSMU) && deltaR(Tau_eta[i], Tau_phi[i], Muon_eta[Muon_idx[0]], Muon_phi[Muon_idx[0]])>DR_OVERLAP_CONE_TAU && deltaR(Tau_eta[i], Tau_phi[i], jet1eta, jet1phi)>isocone && deltaR(Tau_eta[i], Tau_phi[i], jet2eta, jet2phi)>isocone && Tau_pt[i] >=PT_CUT_TAU && abs(Tau_eta[i])<=ETA_CUT_TAU){
                nTau++;
                if(Tau_idDeepTau2017v2p1VSjet[i]>=ID_TAU_RECO_DEEPTAU_VSJET){
                    idx[0] = i;
                    idx[1] = 1;
                    isAtLeastLoose = true;
                } 
                else{
                    if(Tau_idDeepTau2017v2p1VSjet[i]>=ID_TAU_RECO_DEEPTAU_VSJET_LOOSE_MU){
                        idx[0] = i;
                        idx[1] = 0;
                        isAtLeastLoose = true;
                    }
                }
            }
        }
    }
    if(nTau!=1 || isAtLeastLoose == false) idx[1] = -1;                                                                                               
    return idx;
}

bool SameCharge(int GoodLeptonFamily, rvec_i Electron_idx, rvec_i Electron_charge, rvec_i Muon_idx, rvec_i Muon_charge, rvec_i Tau_idx, rvec_i Tau_charge){
    if(GoodLeptonFamily == 0){
        if(Electron_charge[Electron_idx[0]] == Tau_charge[Tau_idx[0]]) return true;
        else return false;
    }
    else if(GoodLeptonFamily == 1){
        if(Muon_charge[Muon_idx[0]] == Tau_charge[Tau_idx[0]]) return true;
        else return false;
    }
    return false;
}

int CountBJets(rvec_f Jet_pt, rvec_f Jet_eta, rvec_f Jet_btagDeepFlavB){
    int nb=0;
    for (int i = 0; i < Jet_pt.size(); i++) {
        if (Jet_btagDeepFlavB[i]>=BTAG_WP_VALUE && Jet_pt[i]>BTAG_PT_CUT && abs(Jet_eta[i])<BTAG_ETA_CUT) nb++;
    }
    return nb;
}

bool BVeto(rvec_f Jet_pt, rvec_f Jet_eta, rvec_f Jet_btagDeepFlavB, rvec_i GoodJet_idx)
{
    bool veto = false;
    for (size_t i = 0; i < GoodJet_idx.size(); i++) {
        //if (Jet_btagDeepFlavB[i]>=WP_btagger[BTAG_ALGO][BTAG_WP])*(Jet_pt[i]>BTAG_PT_CUT)*(abs(Jet_eta[i])<BTAG_ETA_CUT) return true;
        //if ((Jet_btagDeepFlavB[i]>=BTAG_WP_VALUE) && (Jet_pt[i]>BTAG_PT_CUT) && (abs(Jet_eta[i])<BTAG_ETA_CUT)) return false;
        if ((Jet_btagDeepFlavB[GoodJet_idx[i]]>=BTAG_WP_VALUE) && (Jet_pt[GoodJet_idx[i]]>BTAG_PT_CUT) && (abs(Jet_eta[GoodJet_idx[i]])<BTAG_ETA_CUT)) return false;
    }
    return true;
}



bool BVeto_loose(rvec_f Jet_pt, rvec_f Jet_eta, rvec_f Jet_btagDeepFlavB, rvec_i GoodJet_idx)
{
    bool veto = false;
    for (size_t i = 0; i < GoodJet_idx.size(); i++) {
        if ((Jet_btagDeepFlavB[GoodJet_idx[i]]>=BTAG_WP_VALUE_LOOSE) && (Jet_pt[GoodJet_idx[i]]>BTAG_PT_CUT) && (abs(Jet_eta[GoodJet_idx[i]])<BTAG_ETA_CUT)) return false;
    }
    return true;
}


float GetLog2(float x){
    if(x > 0.) return log2(x);
    else return 0;
}

RVec<float> getTES(rvec_f Tau_pt, rvec_i Tau_decayMode, const RVec<UChar_t> &Tau_genPartFlav, bool IsMC, string year){
    string id = "DeepTau2017v2p1VSjet";
    
    float pt_low  = 34;
    float pt_high = 170;
    
    RVec<float> result_all;
    
    for (int i = 0; i < Tau_pt.size(); i++){
        RVec<float> result(3);
                
        int SelectedTau_decayMode = Tau_decayMode[i];
        float SelectedTau_pt = Tau_pt[i];
        //auto SelectedTau_genPartFlav = Tau_genPartFlav[i];
        
        if(IsMC == false){
            result[0] = 1.;
            result[1] = 1.;
            result[2] = 1.;
        }

        else if((SelectedTau_decayMode == 0 || SelectedTau_decayMode == 1 || SelectedTau_decayMode == 10 || SelectedTau_decayMode == 11) && Tau_genPartFlav[i] == 5){ 
            int bin = TauES_dm_DeepTau2017v2p1VSjet_UL2017_hist_low->GetXaxis()->FindBin(SelectedTau_decayMode);
            float tes = TauES_dm_DeepTau2017v2p1VSjet_UL2017_hist_low->GetBinContent(bin);
            float err;
            if (SelectedTau_pt >= pt_high){
                int bin_high = TauES_dm_DeepTau2017v2p1VSjet_UL2017_ptgt100_hist_high->GetXaxis()->FindBin(SelectedTau_decayMode);
                err = TauES_dm_DeepTau2017v2p1VSjet_UL2017_ptgt100_hist_high->GetBinError(bin_high);
            }
            else if (SelectedTau_pt > pt_low){
                int bin_high = TauES_dm_DeepTau2017v2p1VSjet_UL2017_ptgt100_hist_high->GetXaxis()->FindBin(SelectedTau_decayMode);
                float err_high = TauES_dm_DeepTau2017v2p1VSjet_UL2017_ptgt100_hist_high->GetBinError(bin_high);
                float err_low  = TauES_dm_DeepTau2017v2p1VSjet_UL2017_hist_low->GetBinError(bin);
                err      = err_low + (err_high-err_low)/(pt_high-pt_low)*(SelectedTau_pt-pt_low);
            }
            else err = TauES_dm_DeepTau2017v2p1VSjet_UL2017_hist_low->GetBinError(bin);

            result[0] = tes-err;
            result[1] = tes;
            result[2] = tes+err;
        }
        else{
            result[0] = 1.;
            result[1] = 1.;
            result[2] = 1.;
        }
        result_all.emplace_back(result[0]);
        result_all.emplace_back(result[1]);
        result_all.emplace_back(result[2]);
    }
    return result_all;
}



RVec<float> getFES(rvec_f Tau_eta, rvec_i Tau_decayMode, rvec_i Tau_genPartFlav, bool IsMC, string year){
    RVec<float> result_all;
    string id = "DeepTau2017v2p1VSe";
    for (int i = 0; i < Tau_eta.size(); i++){
        RVec<float> result(3);
                
        int SelectedTau_decayMode = Tau_decayMode[i];
        float SelectedTau_eta = Tau_eta[i];
        int SelectedTau_genPartFlav = Tau_genPartFlav[i];
        
        if(IsMC == false){
            result[0] = 1.;
            result[1] = 1.;
            result[2] = 1.;
        }
        
        //else if((SelectedTau_decayMode == 0 || SelectedTau_decayMode == 1) && (SelectedTau_genPartFlav == '0x01' || SelectedTau_genPartFlav == '0x03')){ 
        else if((SelectedTau_decayMode == 0 || SelectedTau_decayMode == 1) && (SelectedTau_genPartFlav == 1 || SelectedTau_genPartFlav == 3)){ 

            int endcap_index_shift = 0;
            if(abs(SelectedTau_eta) >= 1.5) endcap_index_shift = 2;

            float y = TauFES_eta_dm_DeepTau2017v2p1VSe_UL2017_graph->GetY()[SelectedTau_decayMode + endcap_index_shift];
            float yup  = TauFES_eta_dm_DeepTau2017v2p1VSe_UL2017_graph->GetErrorYhigh(SelectedTau_decayMode + endcap_index_shift);
            float ylow = TauFES_eta_dm_DeepTau2017v2p1VSe_UL2017_graph->GetErrorYlow(SelectedTau_decayMode + endcap_index_shift);

            result[0] = y-ylow;
            result[1] = y;
            result[2] = y + yup;
        }
        else{
            result[0] = 1.;
            result[1] = 1.;
            result[2] = 1.;
        }
        result_all.emplace_back(result[0]);
        result_all.emplace_back(result[1]);
        result_all.emplace_back(result[2]);
    }
    return result_all;
}

float M1T(float Lepton_pt, float Lepton_eta, float Lepton_phi, float Lepton_mass, float SelectedTau_pt, float SelectedTau_eta, float SelectedTau_phi, float SelectedTau_mass, float MET_pt, float MET_phi){
    ROOT::Math::PtEtaPhiMVector lep_p4(Lepton_pt, Lepton_eta, Lepton_phi, Lepton_mass);
    ROOT::Math::PtEtaPhiMVector tau_p4(SelectedTau_pt, SelectedTau_eta, SelectedTau_phi, SelectedTau_mass); //*taucorr
    auto leptau_p4 = lep_p4 + tau_p4;
    auto leptau_pt2 = leptau_p4.Perp2();
    auto leptau_px = leptau_p4.Px();
    auto leptau_py = leptau_p4.Py();
    auto leptau_mass2 = leptau_p4.M2();
    auto leptau_et = sqrt(leptau_mass2 + leptau_pt2);
    auto MET_px = MET_pt*cos(MET_phi);
    auto MET_py = MET_pt*sin(MET_phi);
    auto sys_et2 = pow(leptau_et + MET_pt,2.);
    auto sys_pt2 = pow(leptau_px + MET_px, 2.) + pow(leptau_py + MET_py, 2.);
    auto M1T2 = sys_et2 - sys_pt2;
    auto sign_M1T2 = M1T2/abs(M1T2);

    return sign_M1T2*sqrt(sign_M1T2*M1T2);
}

float Mo1(float Lepton_pt, float Lepton_eta, float Lepton_phi, float Lepton_mass, float SelectedTau_pt, float SelectedTau_eta, float SelectedTau_phi, float SelectedTau_mass, float MET_pt, float MET_phi){
    ROOT::Math::PtEtaPhiMVector lep_p4(Lepton_pt, Lepton_eta, Lepton_phi, Lepton_mass);
    ROOT::Math::PtEtaPhiMVector tau_p4(SelectedTau_pt, SelectedTau_eta, SelectedTau_phi, SelectedTau_mass); //*taucorr
    auto leptau_p4 = lep_p4 + tau_p4;
    auto lep_pt = Lepton_pt;
    auto tau_pt = SelectedTau_pt; //*taucorr
    auto leptau_px = leptau_p4.Px();
    auto leptau_py = leptau_p4.Py();
    auto MET_px = MET_pt*cos(MET_phi);
    auto MET_py = MET_pt*sin(MET_phi);
    auto sys_eo2 = pow(lep_pt + tau_pt + MET_pt, 2.);
    auto sys_pt2 = pow(leptau_px + MET_px, 2.) + pow(leptau_py + MET_py, 2.);
    auto Mo12 = sys_eo2 - sys_pt2;
    auto sign_Mo12 = Mo12/abs(Mo12);

    return sign_Mo12*sqrt(sign_Mo12*Mo12);
}
    

unordered_map<int,float> xsecs({
{208,35.85},
{25,51.1},
{26,191.3},
{18,3.697},
{19,0.5297},
{20,0.2529},
{21,0.4062},
{22,0.216},
{23,0.0758},
{28,12.178},
{50,0.2086},
{51,0.1651},
{52,0.05565},
{53,0.01398},
{54,0.2147},
{15,72.1},
{64,27.59},
{67,7181.0},
{28,12.178},
{29,0.014193333333333334},
{30,0.007096666666666667},
{31,0.007096666666666667},
{32,0.007096666666666667},
{33,0.014193333333333334},
{34,0.007096666666666667},
{35,0.007096666666666667},
{36,0.007096666666666667},
{37,0.014193333333333334},
{38,35.85},
{40,1.0315},
{42,0.0118},
{43,2.7757},
{44,0.0896},
{45,0.237},
{46,0.212},
{47,0.952},
{1,0.564},
{2,13.74},
{11,0.003194},
{4,0.003194},
{5,0.003194},
{6,0.003194},
{8,0.003194},
{3,0.001586},
{7,0.001586},
{10,0.001586},
{71,0.002014},
{72,0.01036},
{73,0.01595},
{363,1.},
{364,1.},
{365,1.},
{366,1.},
{367,1.},
{341,1.},
{342,1.},
{343,1.},
{344,1.},
{345,1.},
{110,51.1},
{111,191.3},
{103,3.697},
{104,0.5297},
{105,0.2529},
{106,0.4062},
{107,0.216},
{108,0.0758},
{113,12.178},
{134,0.2086},
{135,0.1651},
{136,0.05565},
{137,0.01398},
{138,0.2147},
{100,72.1},
{148,27.59},
{151,7181.0},
{113,12.178},
{114,0.014193333333333334},
{115,0.007096666666666667},
{116,0.007096666666666667},
{117,0.007096666666666667},
{118,0.014193333333333334},
{119,0.007096666666666667},
{120,0.007096666666666667},
{121,0.007096666666666667},
{122,0.014193333333333334},
{123,35.85},
{125,1.0315},
{127,0.0118},
{128,2.7757},
{129,0.0896},
{130,0.237},
{131,0.212},
{132,0.952},
{86,0.9738},
{87,13.74},
{96,0.003194},
{89,0.003194},
{90,0.003194},
{91,0.003194},
{93,0.003194},
{88,0.001586},
{92,0.001586},
{95,0.001586},
{155,0.002014},
{156,0.01036},
{157,0.01595},
{369,1.},
{370,1.},
{371,1.},
{347,1.},
{348,1.},
{349,1.},
{194,51.1},
{195,191.3},
{187,3.697},
{188,0.5297},
{189,0.2529},
{190,0.4062},
{191,0.216},
{192,0.0758},
{197,12.178},
{219,0.2086},
{220,0.1651},
{221,0.05565},
{222,0.01398},
{223,0.2147},
{184,72.1},
{233,27.59},
{236,7181.0},
{197,12.178},
{198,0.014193333333333334},
{199,0.007096666666666667},
{200,0.007096666666666667},
{201,0.007096666666666667},
{202,0.014193333333333334},
{203,0.007096666666666667},
{204,0.007096666666666667},
{205,0.007096666666666667},
{206,0.014193333333333334},
{207,35.85},
{209,1.0315},
{211,0.0118},
{212,2.7757},
{213,0.0896},
{214,0.237},
{215,0.212},
{216,0.952},
{170,0.9738},
{171,13.74},
{180,0.003194},
{173,0.003194},
{174,0.003194},
{175,0.003194},
{177,0.003194},
{172,0.001586},
{176,0.001586},
{179,0.001586},
{240,0.002014},
{241,0.01036},
{242,0.01595},
{373,1.},
{374,1.},
{375,1.},
{376,1.},
{377,1.},
{351,1.},
{352,1.},
{353,1.},
{354,1.},
{355,1.},
{279,51.1},
{280,191.3},
{272,3.697},
{273,0.5297},
{274,0.2529},
{275,0.4062},
{276,0.216},
{277,0.0758},
{282,12.178},
{304,0.2086},
{305,0.1651},
{306,0.05565},
{307,0.01398},
{308,0.2147},
{269,72.1},
{318,27.59},
{321,7181.0},
{282,12.178},
{283,0.014193333333333334},
{284,0.007096666666666667},
{285,0.007096666666666667},
{286,0.007096666666666667},
{287,0.014193333333333334},
{288,0.007096666666666667},
{289,0.007096666666666667},
{290,0.007096666666666667},
{291,0.014193333333333334},
{292,35.85},
{294,1.0315},
{296,0.0118},
{297,2.7757},
{298,0.0896},
{299,0.237},
{300,0.212},
{301,0.952},
{255,0.9738},
{256,13.74},
{265,0.003194},
{258,0.003194},
{259,0.003194},
{260,0.003194},
{262,0.003194},
{257,0.001586},
{261,0.001586},
{264,0.001586},
{325,0.002014},
{326,0.01036},
{327,0.01595},
{379,1.},
{380,1.},
{381,1.},
{382,1.},
{357,1.},
{358,1.},
{359,1.},
{360,1.},
});


unordered_map<int,float> Nevents({
{208,5674000},
{25,27805647},
{26,53848477},
{18,1511805},
{19,6277000},
{20,5792000},
{21,308442},
{22,1264826},
{23,3723000},
{28,3018000},
{50,5190000},
{51,5072000},
{52,5394000},
{53,5302000},
{54,800000},
{15,37505000},
{64,7934000},
{67,90947213},
{28,3018000},
{29,1994000},
{30,1991000},
{31,1906000},
{32,1997000},
{33,1953000},
{34,1914000},
{35,1958000},
{36,1994000},
{37,1956000},
{38,2300000},
{40,4666982},
{42,1000000},
{43,6134000},
{44,2204000},
{45,1452000},
{46,1977996},
{47,440780},
{1,16862000},
{2,19622315},
{11,500000},
{4,500000},
{5,500000},
{6,500000},
{8,500000},
{3,972000},
{7,927966},
{10,497032},
{71,1921000},
{72,1983000},
{73,1994000},
{363,0},
{364,0},
{365,0},
{366,0},
{367,0},
{341,0},
{342,0},
{343,0},
{344,0},
{345,0},
{110,31562465},
{111,55939475},
{103,1416230},
{104,5401000},
{105,6017000},
{106,308983},
{107,1264826},
{108,3967000},
{113,2900000},
{134,4159000},
{135,4595000},
{136,4554000},
{137,4534000},
{138,698000},
{100,43546000},
{148,7584000},
{151,71839442},
{113,2900000},
{114,1959000},
{115,1927000},
{116,1982000},
{117,1991000},
{118,1984000},
{119,1859000},
{120,1996000},
{121,1973000},
{122,1991000},
{123,2491000},
{125,4882983},
{127,1000000},
{128,6443000},
{129,2893000},
{130,1500000},
{131,2240994},
{132,330462},
{86,15928000},
{87,15890000},
{96,500000},
{89,499000},
{90,500000},
{91,493000},
{93,500000},
{88,992608},
{92,997445},
{95,499183},
{155,1926000},
{156,1956000},
{157,1953000},
{369,0},
{370,0},
{371,0},
{347,0},
{348,0},
{349,0},
{194,29890946},
{195,60212926},
{187,3534208},
{188,13822000},
{189,14036000},
{190,655018},
{191,2891483},
{192,9530000},
{197,7098000},
{219,9854000},
{220,178000},
{221,9898000},
{222,9524000},
{223,1736000},
{184,106724000},
{233,7889000},
{236,195529774},
{197,7098000},
{198,3780000},
{199,3876000},
{200,3998000},
{201,3962000},
{202,3978000},
{203,3985000},
{204,3840000},
{205,3984000},
{206,1910000},
{207,5649000},
{209,6828983},
{211,998000},
{212,12974000},
{213,6440000},
{214,2811630},
{215,5070989},
{216,869559},
{170,40839000},
{171,41708429},
{180,497000},
{173,500000},
{174,498000},
{175,500000},
{177,500000},
{172,997000},
{176,975090},
{179,499000},
{240,1978000},
{241,1952000},
{242,1902000},
{373,0},
{374,0},
{375,0},
{376,0},
{377,0},
{351,0},
{352,0},
{353,0},
{354,0},
{355,0},
{279,29919798},
{280,61613294},
{272,4437068},
{273,19816000},
{274,19608000},
{275,970179},
{276,10450000},
{277,14000000},
{282,9994000},
{304,9894000},
{305,9961999},
{306,9994000},
{307,9889000},
{308,2500000},
{269,145020000},
{318,7940000},
{321,195510810},

{282,9994000},
{283,4904000},
{284,4928000},
{285,4784000},
{286,4998000},
{287,4986000},
{288,4864000},
{289,4948000},
{290,4958000},
{291,4678000},
{292,7956000},
{294,9865972},
{296,940000},
{297,12966000},
{298,9492608},
{299,2987000},
{300,7328993},
{301,1208288},
{255,56886000},
{256,64107525},
{265,500000},
{258,500000},
{259,500000},
{260,500000},
{262,496000},
{257,974000},
{261,994626},
{264,493998},
{325,1927000},
{326,1939000},
{327,1930000},
{379,0},
{380,0},
{381,0},
{382,0},
{357,0},
{358,0},
{359,0},
{360,0},
});


unordered_set<int> data_flags({339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, });
unordered_set<int> dataEle_flags({361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, });
unordered_set<int> dataMu_flags({339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360,});

bool isMC(int SampleFlag){
    bool is_in = data_flags.find(SampleFlag) != data_flags.end();
    if (is_in == true) return false;
    else return true;
}

float getLumi(string Year, bool IsMC){
    if (IsMC == false) return 1.;
    else if(Year == "UL2016APV") return 36.33*0.5373;
    else if(Year == "UL2016") return 36.33*4627;
    //else if (Year == 2017) return 41.48;
    else if (Year == "UL2017") return 41.48;
    else return 59.83;
}

float getXSec(int Sample, bool IsMC){
    if (IsMC == false) return 1.;
    else{
        return  xsecs[Sample];
    }
}

float getNevents(int Sample, bool IsMC){
    if (IsMC == false) return 1.;
    else{
        return  Nevents[Sample];
    }
}


unordered_map<int,bool> IsPdfHessian({
{ 194 , true},
{ 195 , true},
{ 187 , true},
{ 188 , true},
{ 189 , true},
{ 190 , true},
{ 191 , true},
{ 192 , true},
{ 219 , true},
{ 220 , true},
{ 221 , true},
{ 222 , true},
{ 223 , true},
{ 184 , true},
{ 233 , false},
{ 236 , true},
{ 197 , false},
{ 198 , false},
{ 199 , false},
{ 200 , false},
{ 201 , false},
{ 202 , false},
{ 203 , false},
{ 204 , false},
{ 205 , false},
{ 206 , false},
{ 207 , true},
{ 208 , true},
{ 209 , false},
{ 211 , true},
{ 212 , true},
{ 213 , true},
{ 214 , true},
{ 215 , true},
{ 216 , true},
{ 170 , false},
{ 171 , false},
{ 180 , false},
{ 173 , false},
{ 174 , false},
{ 175 , false},
{ 177 , false},
{ 172 , false},
{ 176 , false},
{ 179 , false},
{ 240 , true},
{ 241 , true},
{ 242 , true}, 
});

float Get_isolation(rvec_f Electron_jetRelIso, rvec_i Electron_idx, rvec_f Muon_pfRelIso04_all, rvec_i Muon_idx, int GoodLeptonFamily){
    if(GoodLeptonFamily == 0) return Electron_jetRelIso[Electron_idx[0]];
    else return Muon_pfRelIso04_all[Muon_idx[0]];
}

#endif
