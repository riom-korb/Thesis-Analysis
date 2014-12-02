/****************************************************
*****************************************************
Author: Brock Moir
This class iterates over the full 2011 ATLAS data set or a similar MC sample
and locally stores physics events variables from identified events.

Skimming refers to grabbing the events of interest while 
slimming refers to storing only the variables of interest
*****************************************************
*****************************************************/

// stl includes
#include <iostream>
#include <math.h>
#include <sstream>

// header file
#include <Smearing/SkimSlim.h>

// Library files for 
// the ROOTCore analysis framework
// and the ROOT analysis framework
#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <EventLoop/OutputStream.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TTree.h>

// Mostly how I did my debugging
#define XXX std::cout<<" I am here: "<<__FILE__<<":"<<__LINE__<<std::endl;

// Needed to distribute the algorithm to the workers
ClassImp(SkimSlim);

// Initialize everything in the constructor.
SkimSlim::SkimSlim()
{
	// The definition and method of clustering
	// of jet and MET physics objects that I want.
	// Different analysis might use different definitions
	m_jetName = "jet_AntiKt4LCTopo_";
	m_METName = "MET_RefFinal_";

	// Standard physics object classes defined in RootCore
	final_photons = new D3PDReader::PhotonD3PDObject();
	final_jets = new D3PDReader::JetD3PDObject();
	final_electrons = new D3PDReader::ElectronD3PDObject();
	final_muons = new D3PDReader::MuonD3PDObject();
	final_taus = new D3PDReader::TauD3PDObject();

	// Default minimum cuts on each physics objects' transverse momentum 
	elptcut = 10000;  // electrons
	phptcut = 10000;  // photons
	jetptcut = 10000; // hadronic jets
	muptcut = 10000;  // muons
	tauptcut = 10000; // taus

	//  Selected physics objects from my analysis
	out_photons = new D3PDReader::PhotonD3PDObject();
	out_jets = new D3PDReader::JetD3PDObject();

	// Tools defined in RootCore
	m_tdt = new D3PD::TrigDecisionToolD3PD();
	m_grl = 0;
	m_JetCleaningTool = 0;
	m_JES = 0;
	m_pileupTool = 0;

	// Default configuration variables,
	// can be changed in steering macro
	m_isData = true;
	m_outputStreamName = "output";
	m_deltaPhiB2B = 3.0;
	m_doTrigger = true;
};

// Most of the methods in my SkimSlim class are based 
// on a template from the RootCore tutorial 
EL::StatusCode SkimSlim::setupJob(EL::Job& job)
{
	// manipulate the job before it gets submitted
	// e.g. add output datasets
	EL::OutputStream stream(m_outputStreamName.Data());
	job.outputAdd(stream);
	job.useD3PDReader();
	return EL::StatusCode::SUCCESS;
};

// There is a lot of configuration stuff that needs to be done
// everytime a new data file is read
EL::StatusCode SkimSlim::changeInput(bool firstFile)
{
	// do everything you need to do when we change
	// to a new input file, e.g. reset branch addresses
	//setup trig decision tool
	TTree * tree = wk()->tree();
	if (!tree){
		return EL::StatusCode::FAILURE;
	}
	TFile * file = tree->GetCurrentFile();
	if (!file){
		return EL::StatusCode::FAILURE;
	}
	TString trigMetaTree(tree->GetName());
	trigMetaTree += "Meta/TrigConfTree";
	TTree* confTree = dynamic_cast< TTree* >(file->Get(trigMetaTree));
	if (!confTree) {
		std::cout << "Couldn't retrieve configuration metadata tree!" << std::endl;
		return EL::StatusCode::FAILURE;
	}
	if (!m_tdt->SetEventTree(tree)) {
		std::cout << "Problems with setting the event tree to the TDT" << std::endl;
		return EL::StatusCode::FAILURE;
	}
	if (!m_tdt->SetConfigTree(confTree)) {
		std::cout << "Problems with setting the config tree to the TDT" << std::endl;
		return EL::StatusCode::FAILURE;
	}
	return EL::StatusCode::SUCCESS;
};

// This analysis will take place in parrallel on 
// hundreds of computers world-wide.  These are all the 
// initializations that need to be done on each worker node
// such that the resultant analysis files can be properly merged 
// and stored locally
EL::StatusCode SkimSlim::initialize()
{

	// Each worker needs an output file
	outputFile = wk()->getOutputFile(m_outputStreamName.Data());
 
	// The TTree is a ROOT defined object that 
	// allows for efficient storage of physics events
	outputTree = new TTree("AnalysisTree", "Jet Response");
	outputTree->SetDirectory(outputFile);
	outputTree->SetAutoSave(100000000);
	outputTree->SetAutoFlush(30000000);
	TTree::SetBranchStyle(1);

	// A physics event is a snap shot of a single bunch crossing
	// in the proton-proton collider known as the Large Hadron Collider.
	// Objects in the event correspond to real physics objects 
	// such as high energy electrons and photons 
	event = wk()->d3pdreader();
	event->jets.SetPrefix(m_jetName);
	event->MET.SetPrefix(m_METName);

	// Physics objects from the input Tree to copy to the output Tree
	// Defines which variables are to be saved during slimming
	event->eventinfo.SetActive(kTRUE, "^thr.*|^lbn$|^isSimulation$|^a.*PerXing$|^RunNumber$|^EventNumber$");
	event->MET.SetActive(kTRUE, "^thr.*|^et$|^etx$|^ety$|^phi$");
	event->WriteTo(outputTree);
    out_photons->SetPrefix("ph_");
	out_photons->SetActive(kTRUE, "^thr.*|^pt$|^eta$|^phi$");
	out_photons->WriteTo(outputTree);
	out_jets->SetPrefix("jet_");
	out_jets->SetActive(kTRUE, "^thr.*|^pt$|^eta$|^phi$");
	out_jets->WriteTo(outputTree);

	// Add branches for event variables of interest
	outputTree->Branch("TrigPrescale", &trigPre, "TrigPrescale/i");
	outputTree->Branch("TrigThreshold", &trigThresh, "TrigThresh/f");
	outputTree->Branch("Res", &res, "Res/f");
	outputTree->Branch("weight", &event_weight, "weight/f");
	outputTree->Branch("lumi", &MClumi, "lumi/f");
	outputTree->Branch("pileup", &MCpileup, "pileup/f");
	outputTree->Branch("MCeventweight", &MCeventweight, "MCeventweight/f");
	event->WriteTo(outputTree);

	// Physics analysis tools defined in RootCore
	energyRescaler = new eg2011::EnergyRescaler();
	energyRescaler->useDefaultCalibConstants("2011"); // this is the default...

	m_pileupTool->SetUnrepresentedDataAction(2);
	m_pileupTool->Initialize();

	if (!m_JetCleaningTool){
		throw std::string("No JetCleaningTool configured");
	}
	m_JetCleaningTool->initialize();

	return EL::StatusCode::SUCCESS;
};

// The meat of my analysis, 
// here is where the method of Skimming is defined
// Execute is called for each event,
// optimization is important for turn around time on 
// petabytes of data
EL::StatusCode SkimSlim::execute()
{

	// Variables that need to be initiallized for each event 
	final_electrons->Clear();
	final_photons->Clear();
	final_jets->Clear();
	final_muons->Clear();
	final_taus->Clear();
	out_jets->Clear();
	out_photons->Clear();
	trigPre = 0;
	trigThresh = 0;
	MClumi = 1;
	MCpileup = 1;
	MCeventweight = 1;
	event_weight = 1;

    // A good runs list is a list of periods 
	// when the detectors were operating in a sub-optimal state
	if (!(event->eventinfo.isSimulation() ||
		m_grl.HasRunLumiBlock(event->eventinfo.RunNumber(),
		event->eventinfo.lbn()))) {
		return EL::StatusCode::SUCCESS;
	}
 
	// The trigger system classifies the physics event type 
	// Different analysis are interested in different triggers
	bool passTrig = false;
	m_tdt->GetEntry(wk()->treeEntry());
	for (unsigned int i = 0; i<m_investigateTrig.size(); i++){
		std::string chain(m_investigateTrig[i]);
		// m_chainGroup = new D3PD::ChainGroup( m_tdt->GetChainGroup( chain ) );
		passTrig = m_tdt->IsPassed(chain);
		if (passTrig){
			trigPre = m_tdt->GetChainGroup(chain).GetPrescale();
			trigThresh = m_trigThreshold[i];
			break;
		}
	}
	if (m_doTrigger&&!passTrig) return EL::StatusCode::SUCCESS;
	phptcut = (trigThresh + 10.) * 1000.;
	
	// Clean physics objects and remove overlapping duplicates
    // Each of these methods performs a slightly different task
	this->getElectrons();
	this->getPhotons();
	this->getJets();

	// I want events with only one jet and one photon
	if (!final_electrons->n() == 0) return EL::StatusCode::SUCCESS;
	if (!final_photons->n() == 1) return EL::StatusCode::SUCCESS;
	if (!final_jets->n() == 1) return EL::StatusCode::SUCCESS;

	// Only consider high precision photons
	if (abs((*final_photons)[0].eta())>2.5) return EL::StatusCode::SUCCESS;
    
	// This region is uninstrumented
	if (abs((*final_photons)[0].eta())>1.37
		&& abs((*final_photons)[0].eta())<1.52) return EL::StatusCode::SUCCESS;

	// My analysis requires the jet and photon objects to be back to back
	// in the transverse plane
	if (this->deltaR(0.0, 0.0, (*final_jets)[0].phi(), (*final_photons)[0].phi())
		< m_deltaPhiB2B) return EL::StatusCode::SUCCESS;

	// Only write the leading jet and photon
	out_jets->Add((*final_jets)[0]);
	out_photons->Add((*final_photons)[0]);

	// There is a difference if the input file is a data file or a simulation file
	if (m_isData) event_weight = trigPre;
	else {
		double cross_sec = wk()->metaData()->getDouble("cross_section");
		double num_evnt = wk()->metaData()->getDouble("num_events");
		MClumi = num_evnt / cross_sec;
		MCeventweight = event->mcevt[0].weight()[0];
		MCpileup = m_pileupTool->GetCombinedWeight(
			event->eventinfo.RunNumber(),
			event->eventinfo.mc_channel_number(),
			event->eventinfo.averageIntPerXing());
		event_weight = MCpileup * MCeventweight * trigPre / MClumi;
	}

	// Write out the event
	event->ReadAllActive();
	outputTree->Fill();
	return EL::StatusCode::SUCCESS;
};

// Template method
EL::StatusCode SkimSlim::finalize()
{
	// do everything you need to do before the job on the
	// worker node ends. so far I have no example of what
	// that would be, but the hook is here just in case.
	return EL::StatusCode::SUCCESS;
};

// The electrons should be the purest physics object
// All electrons must pass a minimum transverse momentum
void SkimSlim::getElectrons()
{
	for (int i = 0; i < event->el.n(); i++) {  // loop over all electrons in event
		if (event->el[i].pt() < elptcut) continue;
		final_electrons->Add(event->el[i]);
	}
}

// Photons can be mistaken for electrons,
// it is called overlap when there are two different physics objects
// assigned to the same energy deposit in the detector
// Similar to the electrons, the photons must pass a PT cut
void SkimSlim::getPhotons()
{
	// cleaning and overlap removal
	int k = 0;
	for (int i = 0; i < event->ph.n(); i++) {  // loop over all photons in event
		float delR = 0;
		bool overlap = false;

		// energy rescaling, this might be the reccommended way
		// https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/EnergyScaleResolutionRecommendations
		// https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/EnergyRescaler
		double E = event->ph[i].E();
		double pt = 0;
		if (m_isData) E = 1000.0 * energyRescaler->applyEnergyCorrectionGeV(
			event->ph[i].eta(),
			event->ph[i].phi(),
			event->ph[i].E() / 1000.0,
			event->ph[i].Et() / 1000.0,
			0, "UNCONVERTED_PHOTON");

		pt = E / cosh(event->ph[i].eta());
		if (pt < phptcut) continue;  // photon pt cut 

		// photon id, using tight
		// https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/PhotonReconstruction#PID_variables
		// https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/IsEMIdentification#2_1_2012_data_analyses_with_rele
		if (event->ph[i].isEMTight() == 0) continue;

		// overlap removal
		for (int j = 0; j < (int)final_electrons->n() && !overlap; j++){
			delR = this->deltaR((*final_electrons)[j].eta(),
				event->ph[i].eta(),
				(*final_electrons)[j].phi(),
				event->ph[i].phi());
			if (delR < 0.2) overlap = true;
		}
		if (overlap) continue;

		// after overlap removal fill the photons, and apply the scaling
		(*final_photons) += event->ph[i];
		(*final_photons)[k].E() = E;
		(*final_photons)[k].pt() = pt;
		k++;
	} // end of loop over event photons
}

// Jets are the most complicated objects.  They are the fragments of hadrons
// that are torn into being as quarks and gluons collide at very high energies
void SkimSlim::getJets()
{
	double m_mu = event->eventinfo.averageIntPerXing();
	int npv = 0;

	// this is for the jet rescaler
	for (int i = 0; i < (int)event->vx.n(); ++i) {
		// count the number of vertices with 2 or more tracks
		if (event->vx[i].nTracks() >= 2) npv++;
	}

	// loop over jets for cleaning, scaling and overlap removal
	int k = 0;
	for (int i = 0; i < (int)event->jets.n(); ++i) {
		double Eraw = event->jets[i].constscale_E();
		double eta_det = event->jets[i].constscale_eta();
		double eta, phi, m;

		eta = event->jets[i].EtaOrigin();
		phi = event->jets[i].PhiOrigin();
		m = event->jets[i].MOrigin();

		TLorentzVector jet = m_JES->ApplyOffsetEtaJES(Eraw, eta_det,
			eta, phi, m, m_mu, npv);

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
			event->jets[i].AverageLArQF()))) continue;

		// overlap removal
		bool overlap = false;
		for (int j = 0; j<final_electrons->n() && !overlap; j++){
			if (this->deltaR((*final_electrons)[j].eta(), jet.Eta(),
				(*final_electrons)[j].phi(), jet.Phi())
				< 0.3) overlap = true;
		}
		for (int j = 0; j<final_photons->n() && !overlap; j++){
			if (this->deltaR((*final_photons)[j].eta(), jet.Eta(),
				(*final_photons)[j].phi(), jet.Phi())
				< 0.3) overlap = true;
		}
		if (overlap) continue;

		// fill the jets and scale energy
		(*final_jets) += event->jets[i];
		(*final_jets)[k].E() = jet.E();
		(*final_jets)[k].eta() = jet.Eta();
		(*final_jets)[k].pt() = jet.Pt();
		(*final_jets)[k].phi() = jet.Phi();
		k++;
	}
}

// deltaR is essentially the angular distance between two objects
float SkimSlim::deltaR(float eta1, float eta2, float phi1, float phi2)
{
	float deta = eta2 - eta1;
	float dphi = phi2 - phi1;
	float delR = sqrt(deta*deta + dphi*dphi);
	return delR;
}

