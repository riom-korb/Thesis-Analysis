//////////////////////////////////////////////////////////
// Response Tree Class - Based on Root MakeClass
// Used to create smearing functions
// Authors: Lorraine Courneyea (lorraine.courneyea@cern.ch)
//          Brock Moir         (brock.moir@cern.ch) 
//////////////////////////////////////////////////////////

#ifndef AnalysisTree_h
#define AnalysisTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>


class AnalysisTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UInt_t          EventNumber;
   UInt_t          RunNumber;
   Float_t         actualIntPerXing;
   Float_t         averageIntPerXing;
   Bool_t          isSimulation;
   UInt_t          lbn;
   UInt_t          TrigPrescale;
   Float_t         Res;
   Float_t         TrigThreshold;
   float           weight;
   float           lumi;
   float           pileup;
   float           MCeventweight;
   Float_t         MET_RefFinal_et;
   Float_t         MET_RefFinal_etx;
   Float_t         MET_RefFinal_ety;
   Float_t         MET_RefFinal_phi;
   std::vector<float>   *ph_eta;
   std::vector<float>   *ph_phi;
   std::vector<float>   *ph_pt;
   std::vector<float>   *jet_eta;
   std::vector<float>   *jet_phi;
   std::vector<float>   *jet_pt;

   // List of branches
   TBranch        *b_EventNumber;   //!
   TBranch        *b_RunNumber;   //!
   TBranch        *b_actualIntPerXing;   //!
   TBranch        *b_averageIntPerXing;   //!
   TBranch        *b_isSimulation;   //!
   TBranch        *b_lbn;   //!
   TBranch        *b_MET_RefFinal_et;   //!
   TBranch        *b_MET_RefFinal_etx;   //!
   TBranch        *b_MET_RefFinal_ety;   //!
   TBranch        *b_MET_RefFinal_phi;   //!
   TBranch        *b_ph_eta;   //!
   TBranch        *b_ph_phi;   //!
   TBranch        *b_ph_pt;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_Res;   //!
   TBranch        *b_TrigPrescale;   //!
   TBranch        *b_TrigThresh;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_lumi;
   TBranch        *b_pileup;
   TBranch        *b_MCeventweight;

   AnalysisTree(TTree *tree=0);
   virtual ~AnalysisTree();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual void     Init(TTree *tree);
};

#endif
