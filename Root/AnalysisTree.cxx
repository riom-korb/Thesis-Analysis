//////////////////////////////////////////////////////////
// AnalysisTree Tree Class - Based on Root MakeClass
// Used to create smearing functions
// Authors: Lorraine Courneyea (lorraine.courneyea@cern.ch)
//          Brock Moir        (brock.moir@cern.ch)
//////////////////////////////////////////////////////////

#define AnalysisTree_cxx
#include "Smearing/AnalysisTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

// =================================================================

AnalysisTree::AnalysisTree(TTree *tree)
{
   Init(tree);
}

// =================================================================

AnalysisTree::~AnalysisTree()
{
  return;
}

// =================================================================

Int_t AnalysisTree::GetEntry(Long64_t entry)
{
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

// =================================================================

void AnalysisTree::Init(TTree *tree)
{
   // Set object pointers
   ph_eta = 0;
   ph_phi = 0;
   ph_pt = 0;
   jet_eta = 0;
   jet_phi = 0;
   jet_pt = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   fChain->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
   fChain->SetBranchAddress("actualIntPerXing", &actualIntPerXing, &b_actualIntPerXing);
   fChain->SetBranchAddress("averageIntPerXing", &averageIntPerXing, &b_averageIntPerXing);
   fChain->SetBranchAddress("isSimulation", &isSimulation, &b_isSimulation);
   fChain->SetBranchAddress("lbn", &lbn, &b_lbn);
   
   fChain->SetBranchAddress("MET_RefFinal_et", &MET_RefFinal_et, &b_MET_RefFinal_et);
   fChain->SetBranchAddress("MET_RefFinal_etx", &MET_RefFinal_etx, &b_MET_RefFinal_etx);
   fChain->SetBranchAddress("MET_RefFinal_ety", &MET_RefFinal_ety, &b_MET_RefFinal_ety);
   fChain->SetBranchAddress("MET_RefFinal_phi", &MET_RefFinal_phi, &b_MET_RefFinal_phi);
   fChain->SetBranchAddress("ph_eta", &ph_eta, &b_ph_eta);
   fChain->SetBranchAddress("ph_phi", &ph_phi, &b_ph_phi);
   fChain->SetBranchAddress("ph_pt", &ph_pt, &b_ph_pt);
   fChain->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("Res", &Res, &b_Res);
   fChain->SetBranchAddress("TrigPrescale", &TrigPrescale, &b_TrigPrescale);
   fChain->SetBranchAddress("TrigThreshold", &TrigThreshold, &b_TrigThresh);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("pileup", &pileup, &b_pileup);
   fChain->SetBranchAddress("MCeventweight", &MCeventweight, &b_MCeventweight);

}

// =================================================================

