/****************************************************
*****************************************************
SkimSlim
Author: Brock Moir

This class iterates over the full 2011 ATLAS data and MC sets
and identifies direct photon events to be stored locally.  

*****************************************************
*****************************************************/

// RootCore includes
#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <EventLoop/OutputStream.h>

// My Packages
#include <Smearing/SkimSlim.h>
#include <Smearing/BinningTool.h>

// Root includes
#include <TLorentzVector.h>
#include <TFile.h>
#include <TTree.h>

// stl includes
#include <iostream>
#include <math.h>
#include <sstream>

#define jets jet_AntiKt4LCTopo
#define MET MET_RefFinal

/*preprocessor magic for debugging*/
#define XXX std::cout<<" I am here: "<<__FILE__<<":"<<__LINE__<<std::endl;

// this is needed to distribute the algorithm to the workers
ClassImp(SkimSlim);

SkimSlim :: SkimSlim ()
{
  // put any variable initialization code here
  // e.g. initialize all pointers to 0
  // note that advanced variable intialization (e.g. histogram
  // creation) should be done in initialize

  // object related variables
  m_jetName = "jet_AntiKt4LCTopo_";
  m_METName = "MET_RefFinal_";
  final_photons = new D3PDReader::PhotonD3PDObject();
  final_jets = new D3PDReader::JetD3PDObject();
  final_electrons = new D3PDReader::ElectronD3PDObject();
  final_muons = new D3PDReader::MuonD3PDObject();
  final_taus = new D3PDReader::TauD3PDObject();
  elptcut = 10000;
  phptcut = 10000;
  jetptcut = 10000;
  muptcut = 10000;
  tauptcut = 10000;
  
  // output tree variables
  out_photons = new D3PDReader::PhotonD3PDObject();
  out_jets = new D3PDReader::JetD3PDObject();
  
  // tools
  m_tdt = new D3PD::TrigDecisionToolD3PD();
  m_grl =0;
  m_JetCleaningTool = 0;
  m_JES = 0;
  m_pileupTool = 0;

  // configuration
  m_isData = true; 
  m_outputStreamName="output";
  m_deltaPhiB2B = 3.0;
  m_doTrigger = true;
  //m_directph = true;

};

EL::StatusCode SkimSlim :: setupJob (EL::Job& job)
{

  // manipulate the job before it gets submitted
  // e.g. add output datasets

  EL::OutputStream stream (m_outputStreamName.Data());
  job.outputAdd	(stream);
  job.useD3PDReader();

  return EL::StatusCode::SUCCESS;
};



EL::StatusCode SkimSlim :: changeInput (bool firstFile)
{
  // do everything you need to do when we change
  // to a new input file, e.g. reset branch addresses

  //setup trig decision tool
  TTree * tree = wk()->tree();
  if(!tree){
    return EL::StatusCode::FAILURE;
  }
  TFile * file = tree->GetCurrentFile();
  if(!file){
    return EL::StatusCode::FAILURE;
  }

  TString trigMetaTree(tree->GetName());
  trigMetaTree+="Meta/TrigConfTree";
  TTree* confTree = dynamic_cast< TTree* >( file->Get(trigMetaTree ) );
  if( ! confTree ) {
    std::cout << "Couldn't retrieve configuration metadata tree!" << std::endl;
    return EL::StatusCode::FAILURE;
  }
  if( ! m_tdt->SetEventTree( tree ) ) {
    std::cout << "Problems with setting the event tree to the TDT" << std::endl;
    return EL::StatusCode::FAILURE;
  }
  if( ! m_tdt->SetConfigTree( confTree ) ) {
    std::cout << "Problems with setting the config tree to the TDT" << std::endl;
    return EL::StatusCode::FAILURE;
  }

  return EL::StatusCode::SUCCESS;
};



EL::StatusCode SkimSlim :: initialize ()
{

  // do everything you need to do to set up initially
  // on the worker node, e.g. create histograms and
  // output trees.

  event = wk()->d3pdreader();
  event->jets.SetPrefix( m_jetName );
  event->MET.SetPrefix( m_METName );

  // Setup a TTree in a new output file
  outputFile= wk()->getOutputFile (m_outputStreamName.Data());
  outputTree= new TTree("AnalysisTree","Jet Response");
  outputTree->SetDirectory(outputFile);
  outputTree->SetAutoSave(100000000);
  outputTree->SetAutoFlush( 30000000 );
  TTree::SetBranchStyle(1);
 
  //Objects from the input Tree to copy to output
  event->eventinfo.SetActive(kTRUE, "^thr.*|^lbn$|^isSimulation$|^a.*PerXing$|^RunNumber$|^EventNumber$");
  event->MET.SetActive(kTRUE, "^thr.*|^et$|^etx$|^ety$|^phi$");
  event->WriteTo( outputTree );

  // Setup selected objects to copy to output
  if(m_directph){
    out_photons->SetPrefix( "ph_" );
    out_photons->SetActive( kTRUE, "^thr.*|^pt$|^eta$|^phi$" );
    out_photons->WriteTo( outputTree );
  }

  out_jets->SetPrefix( "jet_" );
  out_jets->SetActive( kTRUE, "^thr.*|^pt$|^eta$|^phi$" );
  out_jets->WriteTo( outputTree );

  // Add individual branches to the output tree
  outputTree->Branch("TrigPrescale",     &trigPre,        "TrigPrescale/i");
  outputTree->Branch("TrigThreshold",    &trigThresh,     "TrigThresh/f");
  outputTree->Branch("Res",              &res,            "Res/f");
  outputTree->Branch("weight",           &event_weight,   "weight/f");
  outputTree->Branch("lumi",             &MClumi,         "lumi/f");
  outputTree->Branch("pileup",           &MCpileup,       "pileup/f");
  outputTree->Branch("MCeventweight",    &MCeventweight,  "MCeventweight/f");
  event->WriteTo( outputTree );

  this->BookCutflow();

  // EnergyRescaler
  energyRescaler = new eg2011::EnergyRescaler();
  energyRescaler->useDefaultCalibConstants("2011");  // this is the default...

  // Pilup Reweighting
  m_pileupTool->SetUnrepresentedDataAction(2); 
  m_pileupTool->Initialize();


  if (!m_JetCleaningTool){
    throw std::string ("No JetCleaningTool configured");
  }
  m_JetCleaningTool->initialize();

  return EL::StatusCode::SUCCESS;

};



EL::StatusCode SkimSlim :: execute ()
{
  // initialize event variables
  final_electrons->Clear();
  final_photons->Clear();
  final_jets->Clear();
  final_muons->Clear();
  final_taus->Clear();

  out_jets->Clear();
  out_photons->Clear();

  trigPre=0;
  trigThresh=0;
  res=-666;
  MClumi = 1;
  MCpileup = 1;
  MCeventweight = 1;
  event_weight=1;

  this->FillCutflow_event("All");

  /** check GRL */
  if (!(event->eventinfo.isSimulation() ||
      m_grl.HasRunLumiBlock(event->eventinfo.RunNumber(),
                            event->eventinfo.lbn()))) {
    return EL::StatusCode::SUCCESS;
  }
  this->FillCutflow_event("GRL");
    
  /** check Trigger */
  bool passTrig=false;
  m_tdt->GetEntry(wk()->treeEntry() );
  for(unsigned int i=0; i<m_investigateTrig.size(); i++){
    std::string chain(m_investigateTrig[i]);
    //    m_chainGroup = new D3PD::ChainGroup( m_tdt->GetChainGroup( chain ) );
    passTrig=m_tdt->IsPassed(chain);
    if(passTrig){
      trigPre=m_tdt->GetChainGroup( chain ).GetPrescale();
      trigThresh=m_trigThreshold[i];
      break;
    }
  }
  if(m_doTrigger&&!passTrig) return EL::StatusCode::SUCCESS;
  this->FillCutflow_event("Trig");
  
  if (m_directph) phptcut = (trigThresh + 10.) * 1000.;
  else jetptcut = (trigThresh + 10.) * 1000.;

  // LAR remove events affected by LAr bursts and data corruption
  if (event->eventinfo.larError() >= 2) return EL::StatusCode::SUCCESS;
  this->FillCutflow_event("LArError");
  
  // clean objects and overlap removal
  this->getElectrons();
  this->getPhotons();
  this->getJets();
  this->getMuons();
  this->getTaus();

  this->FillCutflow_directph("All");
  // I want clean events with just one jet and one photon
  if(! final_electrons->n()==0) return EL::StatusCode::SUCCESS;
  this->FillCutflow_directph("Electrons");
  if(! final_muons->n()==0)     return EL::StatusCode::SUCCESS;
  this->FillCutflow_directph("Muons");
  if(! final_taus->n()==0)      return EL::StatusCode::SUCCESS;
  this->FillCutflow_directph("Taus");
  if(! final_photons->n()==1)   return EL::StatusCode::SUCCESS;
  this->FillCutflow_directph("Photons");
  if(! final_jets->n()==1)      return EL::StatusCode::SUCCESS;
  this->FillCutflow_directph("Jets");

  //check for photons in the EM Calo Barrel
  if(abs((*final_photons)[0].eta())>2.5) return EL::StatusCode::SUCCESS; 
  this->FillCutflow_directph("#eta_{#gamma}");

  //but not in the crack region
  if(abs((*final_photons)[0].eta())>1.37 
     && abs((*final_photons)[0].eta())<1.52) return EL::StatusCode::SUCCESS;
  this->FillCutflow_directph("Crack");

  if(this->deltaR(0.0, 0.0,(*final_jets)[0].phi(),(*final_photons)[0].phi())
     < m_deltaPhiB2B) return EL::StatusCode::SUCCESS;  
  this->FillCutflow_directph("#Delta#phi");
  
  res = ((*final_jets)[0].pt() - (*final_photons)[0].pt())
        /(*final_photons)[0].pt();

  out_jets->Add((*final_jets)[0]);
  out_photons->Add((*final_photons)[0]);

  if(m_isData) event_weight = trigPre;
  else {
    double cross_sec = wk()->metaData()->getDouble("cross_section");
    double num_evnt  = wk()->metaData()->getDouble("num_events");
    MClumi = num_evnt/cross_sec;
    MCeventweight  = event->mcevt[0].weight()[0];
    MCpileup = m_pileupTool->GetCombinedWeight(
			       event->eventinfo.RunNumber(),
			       event->eventinfo.mc_channel_number(),
			       event->eventinfo.averageIntPerXing());
    event_weight     = MCpileup * MCeventweight * trigPre / MClumi;
  }
  
  // write out the event
  event->ReadAllActive();
  outputTree->Fill();

  return EL::StatusCode::SUCCESS;
};



EL::StatusCode SkimSlim :: finalize ()
{
  // do everything you need to do before the job on the
  // worker node ends.  so far I have no example of what
  // that would be, but the hook is here just in case.
  return EL::StatusCode::SUCCESS;
};

void SkimSlim::getElectrons()
{

  for ( int i = 0; i < event->el.n(); i++ ) {

    if (event->el[i].pt() < elptcut) continue;      
    final_electrons->Add(event->el[i]);
  }     
}

void SkimSlim::getPhotons()
{
  // cleaning and overlap removal
  int k=0;
  for ( int i=0; i < event->ph.n(); i++ ) {
    float delR = 0;
    bool overlap = false;
    // energy rescaling, this might be the reccommended way
    // https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/EnergyScaleResolutionRecommendations
    // https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/EnergyRescaler
    double E = event->ph[i].E();
    double pt= 0;
    if(m_isData) E = 1000.0 * energyRescaler->applyEnergyCorrectionGeV(
	  				  event->ph[i].eta(),
	  				  event->ph[i].phi(),  
					  event->ph[i].E()/1000.0,  
					  event->ph[i].Et()/1000.0,  
					  0, "UNCONVERTED_PHOTON");
    /*else E = 1000.0 * energyRescaler->getSmearingCorrection(
	  				  event->ph[i].eta(),
					  event->ph[i].E()/1000.0, 
					  false, "2011");
    */
    pt=E/cosh(event->ph[i].eta());

    // photon pt cut
    if (pt < phptcut) continue;

    // photon id, using tight
    // https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/PhotonReconstruction#PID_variables
    // https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/IsEMIdentification#2_1_2012_data_analyses_with_rele
    if(event->ph[i].isEMTight() == 0) continue;

    // overlap removal
    for(int j=0; j < (int)final_electrons->n() && !overlap; j++){ 
      delR = this->deltaR((*final_electrons)[j].eta(), 
			  event->ph[i].eta(), 
			  (*final_electrons)[j].phi(),
			  event->ph[i].phi()); 
      if(delR < 0.2) overlap = true;
    }
    if(overlap) continue;
 
    // after overlap removal fill the photons, and apply the scaling
    (*final_photons) += event->ph[i];  
    (*final_photons)[k].E()=E;
    (*final_photons)[k].pt()=pt;
    k++;
  } // end of cleaning and overlap removal
}


void SkimSlim::getJets()
{
  double m_mu = event->eventinfo.averageIntPerXing();
  int npv = 0;
  
  // this is for the jet rescaler
  for (int i = 0; i < (int)event->vx.n(); ++i) {
    // count the number of vertices with 2 or more tracks
    if ( event->vx[i].nTracks() >= 2 ) npv++;
  }
  
  // loop over jets for cleaning, scaling and overlap removal
  int k=0;
  for (int i = 0; i < (int)event->jets.n(); ++i) {
    //bool overlap = false;
    
    // apply rescaling to the jets 
    double Eraw    = event->jets[i].constscale_E();
    double eta_det = event->jets[i].constscale_eta();
    double eta, phi, m;
    eta     = event->jets[i].EtaOrigin();
    phi     = event->jets[i].PhiOrigin();
    m       = event->jets[i].MOrigin();
    
    TLorentzVector jet = m_JES->ApplyOffsetEtaJES(Eraw,eta_det,
                                                  eta,phi,m,m_mu,npv);
    
    if (jet.Pt() < jetptcut) continue;
    if (!(m_JetCleaningTool->accept(jet.Eta(),
				  event->jets[i].NegativeE(),
                                  event->jets[i].hecf(),
                                  event->jets[i].HECQuality(),
                                  event->jets[i].emfrac(),
                                  event->jets[i].sumPtTrk(),
                                  event->jets[i].LArQuality(),
                                  event->jets[i].Timing(),
                                  event->jets[i].fracSamplingMax(),
				  event->jets[i].AverageLArQF())) ) continue;

    
    // overlap removal
   bool overlap = false; 
   for(int j=0; j<final_electrons->n() && !overlap; j++){ 
      if (this->deltaR((*final_electrons)[j].eta(), jet.Eta(), 
		       (*final_electrons)[j].phi(), jet.Phi()) 
	  < 0.3) overlap = true;
    }
    
   
    for(int j=0; j<final_photons->n() && !overlap; j++){ 
      if (this->deltaR((*final_photons)[j].eta(), jet.Eta(), 
		       (*final_photons)[j].phi(), jet.Phi()) 
	  < 0.3) overlap = true;
    }
    if (overlap) continue;
    

    // fill the jets and scale energy
    (*final_jets) += event->jets[i];
    (*final_jets)[k].E()=jet.E();
    (*final_jets)[k].eta()=jet.Eta();
    (*final_jets)[k].pt()=jet.Pt();
    (*final_jets)[k].phi()=jet.Phi();
    k++;
    
  }
}

void SkimSlim::getMuons(){
  for(int i=0; i < (int)event->mu_staco.n(); i++){
    // overlap removal & pt cut
    if(event->mu_staco[i].pt()<muptcut) continue;
    bool overlap=false;
    for(int j=0; j < (int)final_electrons->n() && !overlap; j++){ 
      if (this->deltaR((*final_electrons)[j].eta(), event->mu_staco[i].eta(), 
		       (*final_electrons)[j].phi(), event->mu_staco[i].phi()) 
	  < 0.2) overlap = true;
    }
    for(int j=0; j < (int)final_photons->n() && !overlap; j++){ 
      if (this->deltaR((*final_photons)[j].eta(), event->mu_staco[i].eta(), 
		       (*final_photons)[j].phi(), event->mu_staco[i].phi()) 
	  < 0.2) overlap = true;
    }
    for(int j=0; j < (int)final_jets->n() && !overlap; j++){ 
      if (this->deltaR((*final_jets)[j].eta(), event->mu_staco[i].eta(), 
		       (*final_jets)[j].phi(), event->mu_staco[i].phi()) 
	  < 0.3) overlap = true;
    }
    if(overlap) continue;
    (*final_muons) += event->mu_staco[i];
  }
}

void SkimSlim::getTaus(){
  for(int i=0; i < (int)event->tau.n(); i++){
    if(event->tau[i].pt()<tauptcut) continue;
    // overlap removal
    bool overlap=false;
    for(int j=0; j < (int)final_electrons->n() && !overlap; j++){ 
      if (this->deltaR((*final_electrons)[j].eta(), event->tau[i].eta(), 
		       (*final_electrons)[j].phi(), event->tau[i].phi()) 
	  < 0.2) overlap = true;
    }
    for(int j=0; j < (int)final_photons->n() && !overlap; j++){ 
      if (this->deltaR((*final_photons)[j].eta(), event->tau[i].eta(), 
		       (*final_photons)[j].phi(), event->tau[i].phi()) 
	  < 0.2) overlap = true;
    }
    for(int j=0; j < (int)final_jets->n() && !overlap; j++){ 
      if (this->deltaR((*final_jets)[j].eta(), event->tau[i].eta(), 
		       (*final_jets)[j].phi(), event->tau[i].phi()) 
	  < 0.3) overlap = true;
    }
    for(int j=0; j < (int)final_muons->n() && !overlap; j++){ 
      if (this->deltaR((*final_muons)[j].eta(), event->tau[i].eta(), 
		       (*final_muons)[j].phi(), event->tau[i].phi()) 
	  < 0.2) overlap = true;
    }
    if(overlap) continue;
    (*final_taus) += event->tau[i];
  }

}

float SkimSlim::deltaR (float eta1, float eta2, float phi1, float phi2)
{
  float deta = eta2 - eta1;
  float dphi = phi2 - phi1;
  float delR = sqrt(deta*deta + dphi*dphi);
  return delR; 
}

void SkimSlim::BookCutflow()
{
  cutflow_event_names.push_back("All");
  cutflow_event_names.push_back("GRL");
  cutflow_event_names.push_back("Trig");
  cutflow_event_names.push_back("LArError");

  cutflow_event    = new TH1F("cutflow", "cutflow", 
			      (int)cutflow_event_names.size(), 
			      -0.5, (int)cutflow_event_names.size()-0.5);
  for (int i = 0; i < (int)cutflow_event_names.size(); ++i)
    {
      cutflow_event->GetXaxis()->SetBinLabel(i+1, cutflow_event_names[i]);
    }
  wk()->addOutput(cutflow_event);

  cutflow_directph_names.push_back("All");
  cutflow_directph_names.push_back("Electrons");
  cutflow_directph_names.push_back("Muons");
  cutflow_directph_names.push_back("Taus");
  cutflow_directph_names.push_back("Photons");
  cutflow_directph_names.push_back("Jets");
  if(m_directph){
    cutflow_directph_names.push_back("#eta_{#gamma}");
    cutflow_directph_names.push_back("Crack"); 
  }
  cutflow_directph_names.push_back("#Delta#phi");
 
  cutflow_directph = new TH1F("cutflow", "cutflow", 
			      (int)cutflow_directph_names.size(), 
			      -0.5, (int)cutflow_directph_names.size()-0.5);
  for (int i = 0; i < (int)cutflow_directph_names.size(); ++i)
    {
     cutflow_directph->GetXaxis()->SetBinLabel(i+1, cutflow_directph_names[i]);
    }

  wk()->addOutput(cutflow_directph);

}

void SkimSlim::FillCutflow_event(TString setName)
{
  cutflow_event->Fill(cutflow_event->GetXaxis()->GetBinCenter(cutflow_event->GetXaxis()->FindBin(setName.Data())));
}

void SkimSlim::FillCutflow_directph(TString setName)
{
  cutflow_directph->Fill(cutflow_directph->GetXaxis()->GetBinCenter(cutflow_directph->GetXaxis()->FindBin(setName.Data())));
}
