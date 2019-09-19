// -*- C++ -*-
//
// Package:    L1Trigger/phase2L1TauAnalyzer
// Class:      phase2L1TauAnalyzer
// 
/**\class phase2L1TauAnalyzer phase2L1TauAnalyzer.cc L1Trigger/phase2L1TauAnalyzer/plugins/phase2L1TauAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Isobel Ojalvo
//         Created:  Fri, 26 May 2017 15:44:30 GMT
//
//


// system include files
#include <memory>

// Math Include
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include <TLorentzVector.h>
#include <memory>
#include <math.h>
#include <vector>
#include <list>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

//#include "L1Trigger/phase2Demonstrator/interface/triggerGeometryTools.hh"
#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"

//Vertex and gen particle
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"


#include "L1Trigger/phase2L1TauAnalyzer/plugins/helpers.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/L1Trigger/interface/L1PFTau.h"

#include "DataFormats/Candidate/interface/Candidate.h"

//#include "DataFormats/L1Trigger/interface/L1PFObject.h"

#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/L1TrackTrigger/interface/L1TkPrimaryVertex.h"

#include "L1Trigger/Phase2L1Taus/interface/L1PFTauProducer.hh"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

using namespace l1t;

class phase2L1TauAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit phase2L1TauAnalyzer(const edm::ParameterSet&);
  ~phase2L1TauAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  // ----------member data ---------------------------
  typedef std::vector<reco::GenParticle> GenParticleCollectionType;

  struct genVisTau{
    reco::Candidate::LorentzVector p4;
    int decayMode;
  };

  const CaloSubdetectorGeometry * ebGeometry;
  edm::ESHandle<CaloGeometry> caloGeometry_;

  //edm::EDGetTokenT< L1CaloClusterCollection > L1ClustersToken_;
  //  edm::EDGetTokenT<l1t::PFCandidateCollection> L1PFToken_;
  edm::EDGetTokenT< vector<l1t::PFCandidate> > L1PFToken_;
  edm::EDGetTokenT<L1TkPrimaryVertexCollection>    pvToken_;
  edm::EDGetTokenT< L1PFTauCollection > L1PFTauToken_;
  edm::EDGetTokenT<std::vector<reco::GenParticle> > genToken_;
  edm::EDGetTokenT< vector<pat::Tau>  > MiniTausToken_;
  edm::EDGetTokenT< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > ttTrackToken_;
  edm::EDGetTokenT<std::vector<pat::PackedCandidate> > PackedCands_;
  edm::EDGetTokenT<EcalEBTrigPrimDigiCollection> ecalTPGBToken_;

  edm::InputTag genSrc_;
  edm::InputTag L1TrackInputTag;

  TTree* pi0Tree;
  TTree* piPlusTree;

  TTree* oneProngTree;
  TTree* oneProngPi0Tree;
  TTree* threeProngTree;
  
  TTree* efficiencyTree;

  const reco::Vertex *pv;

  double genPt, genEta, genPhi;
  int decayMode, run, lumi, event;
  double l1Pt, l1Eta, l1Phi, l1DM;
  double zVTX, l1TauZ, l1PVDZ;
  double l1TauChargedIso,l1TauNeutralIso,l1TauRawIso;
  double gen1ProngPt, gen1ProngEta, gen1ProngPhi;
  double gen3ProngPt, gen3ProngEta, gen3ProngPhi;
  double gen1ProngPi0Pt, gen1ProngPi0Eta, gen1ProngPi0Phi;
  double genPiZeroPt, genPiZeroEta, genPiZeroPhi;

  double track1_pt, track1_eta, track1_phi;
  double track2_pt, track2_eta, track2_phi;
  double track3_pt, track3_eta, track3_phi;
  double track2_dR, track2_dz;
  double track3_dR, track3_dz;

  double l1Track_pt, l1Track_eta, l1Track_phi;
  double l1Track_dR, l1Track_dz, l1Track_dRmin;
  double l1TrackRecoTauPtDiff;

  double pfCand_pt, pfCand_eta, pfCand_phi;
  double pfCand_dR, pfCand_dRmin, pfCand_dz;

  double l1StripPt, l1StripEta, l1StripPhi;
  double l1StripDR;

  double other_track_pt, other_track_eta, other_track_phi, other_track_dR, other_track_dz;

  double recoPt, recoEta, recoPhi; 
  float recoChargedIso, recoNeutralIso, recoRawIso; 
  int l1IsoVLoose, l1IsoLoose, l1IsoMedium, l1IsoTight;
  int l1RelIsoVLoose, l1RelIsoLoose, l1RelIsoMedium, l1RelIsoTight;
  int recoDecayMode; 

  int gen1ProngDecayMode;
  int gen1ProngPi0DecayMode;
  int gen3ProngDecayMode;

  double genPi0Pt, genPi0Eta,  genPi0Phi;
  double objectPi0Pt, objectPi0Eta, objectPi0Phi;
  double objectPi0DeltaEta, objectPi0DeltaPhi;

  double genPiPlusPt, genPiPlusEta,  genPiPlusPhi;
  double objectPiPlusPt, objectPiPlusEta, objectPiPlusPhi;
  double objectPiPlusDeltaEta, objectPiPlusDeltaPhi;
  int objectPiPlusIsChargedHadron;

  TH1F* nEvents;
  TH1F* track_pt;
  TH1F* track_pt_eta2p1;
  TH1F* track_pt_eta2p4;

  TH1F* l1Tau_pt;     
  TH1F* l1Tau_pt_eta2p4;     
  TH1F* l1Tau_pt_eta2p1;

  TH1F* l1SingleProngTau_pt;     
  TH1F* l1SingleProngTau_pt_eta2p4;     
  TH1F* l1SingleProngTau_pt_eta2p1;

  TH1F* l1SingleProngPi0Tau_pt;     
  TH1F* l1SingleProngPi0Tau_pt_eta2p4;     
  TH1F* l1SingleProngPi0Tau_pt_eta2p1;

  TH1F* l1SingleProngTauIso_pt;     
  TH1F* l1SingleProngTauIso_pt_eta2p4;     
  TH1F* l1SingleProngTauIso_pt_eta2p1;

  TH1F* l1ThreeProngTau_pt;     
  TH1F* l1ThreeProngTau_pt_eta2p4;     
  TH1F* l1ThreeProngTau_pt_eta2p1;

  TH1F* l1ThreeProngTauIso_pt;     
  TH1F* l1ThreeProngTauIso_pt_eta2p4;     
  TH1F* l1ThreeProngTauIso_pt_eta2p1;

  //TH2F* l1OtherTracks;
  //TH2F* l1EcalCrystals;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
phase2L1TauAnalyzer::phase2L1TauAnalyzer(const edm::ParameterSet& cfg):
  //L1ClustersToken_( consumes< L1CaloClusterCollection >(cfg.getParameter<edm::InputTag>("L1Clusters"))),
  L1PFToken_(       consumes< vector<l1t::PFCandidate> >(cfg.getParameter<edm::InputTag>("L1PFObjects"))),
  pvToken_(         consumes<L1TkPrimaryVertexCollection> (cfg.getParameter<edm::InputTag>("L1VertexInputTag"))),
  L1PFTauToken_(    consumes< L1PFTauCollection    >(cfg.getParameter<edm::InputTag>("l1TauObjects"))),
  MiniTausToken_(   consumes< vector<pat::Tau>     >(cfg.getParameter<edm::InputTag>("miniTaus"))),
  PackedCands_(     consumes< std::vector<pat::PackedCandidate> >(cfg.getParameter<edm::InputTag>("packedCandidates"))),
  ecalTPGBToken_(   consumes<EcalEBTrigPrimDigiCollection>(cfg.getParameter<edm::InputTag>("ecalTPGsBarrel"))),
  genSrc_ ((        cfg.getParameter<edm::InputTag>( "genParticles")))
{
  //now do what ever initialization is needed
  usesResource("TFileService");
  genToken_ =     consumes<std::vector<reco::GenParticle> >(genSrc_);
  
  L1TrackInputTag = cfg.getParameter<edm::InputTag>("L1TrackInputTag");
  ttTrackToken_ = consumes< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > >(L1TrackInputTag);  

  edm::Service<TFileService> fs;
  efficiencyTree = fs->make<TTree>("efficiencyTree", "Crystal cluster individual crystal pt values");
  efficiencyTree->Branch("run",    &run,     "run/I");
  efficiencyTree->Branch("lumi",   &lumi,    "lumi/I");
  efficiencyTree->Branch("event",  &event,   "event/I");
  
  efficiencyTree->Branch("genPt",  &genPt,   "genPt/D");
  efficiencyTree->Branch("genEta", &genEta,   "genEta/D");
  efficiencyTree->Branch("genPhi", &genPhi,   "genPhi/D");
  efficiencyTree->Branch("genDM",  &decayMode,   "genDM/I");

  efficiencyTree->Branch("recoPt",  &recoPt,   "recoPt/D");
  efficiencyTree->Branch("recoEta", &recoEta,   "recoEta/D");
  efficiencyTree->Branch("recoPhi", &recoPhi,   "recoPhi/D");
  efficiencyTree->Branch("recoChargedIso", &recoChargedIso,   "recoChargedIso/D");
  efficiencyTree->Branch("recoNeutralIso", &recoNeutralIso,   "recoNeutralIso/D");
  efficiencyTree->Branch("recoRawIso", &recoRawIso,   "recoDBIso/D");
  efficiencyTree->Branch("recoDM",  &recoDecayMode,   "recoDM/I");

  efficiencyTree->Branch("l1Track_pt",  &l1Track_pt, "l1Track_pt/D");
  efficiencyTree->Branch("l1Track_eta", &l1Track_eta, "l1Track_eta/D");
  efficiencyTree->Branch("l1Track_phi", &l1Track_phi, "l1Track_phi/D");
  efficiencyTree->Branch("l1Track_dR", &l1Track_dR, "l1Track_dR/D");
  efficiencyTree->Branch("l1Track_dz", &l1Track_dz, "l1Track_dz/D");
  efficiencyTree->Branch("l1TrackRecoTauPtDiff", &l1TrackRecoTauPtDiff, "l1TrackRecoTauPtDiff/D");

  efficiencyTree->Branch("pfCand_pt",  &pfCand_pt, "pfCand_pt/D");
  efficiencyTree->Branch("pfCand_eta", &pfCand_eta, "pfCand_eta/D");
  efficiencyTree->Branch("pfCand_phi", &pfCand_phi, "pfCand_phi/D");
  efficiencyTree->Branch("pfCand_dz",  &pfCand_dz, "pfCand_dz/D");
  efficiencyTree->Branch("pfCand_dR",  &pfCand_dR, "pfCand_dR/D");

  efficiencyTree->Branch("l1Pt",  &l1Pt,   "l1Pt/D");
  efficiencyTree->Branch("l1Eta", &l1Eta,   "l1Eta/D");
  efficiencyTree->Branch("l1Phi", &l1Phi,   "l1Phi/D");
  efficiencyTree->Branch("l1TauZ", &l1TauZ, "l1TauZ/D");
  efficiencyTree->Branch("zVTX",   &zVTX,   "zVTX/D");
  efficiencyTree->Branch("l1PVDZ", &l1PVDZ, "l1PVDZ/D");

  efficiencyTree->Branch("l1StripPt",  &l1StripPt,  "l1StripPt/D");
  efficiencyTree->Branch("l1StripEta", &l1StripEta, "l1StripEta/D");
  efficiencyTree->Branch("l1StripPhi", &l1StripPhi, "l1StripPhi/D");
  efficiencyTree->Branch("l1StripDR",  &l1StripDR,  "l1StripDR/D");

  efficiencyTree->Branch("l1RawIso",     &l1TauRawIso,     "l1RawIso/D");
  efficiencyTree->Branch("l1ChargedIso", &l1TauChargedIso, "l1ChargedIso/D");
  efficiencyTree->Branch("l1NeutralIso", &l1TauNeutralIso, "l1NeutralIso/D");

  efficiencyTree->Branch("l1IsoVLoose", &l1IsoVLoose, "l1IsoVLoose/I");
  efficiencyTree->Branch("l1IsoLoose",  &l1IsoLoose,  "l1IsoLoose/I");
  efficiencyTree->Branch("l1IsoMedium", &l1IsoMedium, "l1IsoMedium/I");
  efficiencyTree->Branch("l1IsoTight",  &l1IsoTight,  "l1IsoTight/I");

  efficiencyTree->Branch("l1RelIsoVLoose", &l1RelIsoVLoose, "l1RelIsoVLoose/I");
  efficiencyTree->Branch("l1RelIsoLoose",  &l1RelIsoLoose,  "l1RelIsoLoose/I");
  efficiencyTree->Branch("l1RelIsoMedium", &l1RelIsoMedium, "l1RelIsoMedium/I");
  efficiencyTree->Branch("l1RelIsoTight",  &l1RelIsoTight,  "l1RelIsoTight/I");

  efficiencyTree->Branch("l1DM", &l1DM,   "l1DM/D");


  pi0Tree = fs->make<TTree>("pi0Tree", "Crystal cluster individual crystal pt values");
  pi0Tree->Branch("run",    &run,     "run/I");
  pi0Tree->Branch("lumi",   &lumi,    "lumi/I");
  pi0Tree->Branch("event",  &event,   "event/I");
  
  pi0Tree->Branch("genPt",     &genPi0Pt,           "genPt/D");
  pi0Tree->Branch("genEta",    &genPi0Eta,          "genEta/D");
  pi0Tree->Branch("genPhi",    &genPi0Phi,          "genPhi/D");

  pi0Tree->Branch("objectPt",  &objectPi0Pt,        "objectPt/D");
  pi0Tree->Branch("objectEta", &objectPi0Eta,       "objectEta/D");
  pi0Tree->Branch("objectPhi", &objectPi0Phi,       "objectPhi/D");

  pi0Tree->Branch("deltaEta",  &objectPi0DeltaEta,   "deltaEta/D");
  pi0Tree->Branch("deltaPhi",  &objectPi0DeltaPhi,   "deltaPhi/D");

  piPlusTree = fs->make<TTree>("piPlusTree", "Crystal cluster individual crystal pt values");
  piPlusTree->Branch("run",    &run,     "run/I");
  piPlusTree->Branch("lumi",   &lumi,    "lumi/I");
  piPlusTree->Branch("event",  &event,   "event/I");
  
  piPlusTree->Branch("genPt",     &genPiPlusPt,           "genPt/D");
  piPlusTree->Branch("genEta",    &genPiPlusEta,          "genEta/D");
  piPlusTree->Branch("genPhi",    &genPiPlusPhi,          "genPhi/D");

  piPlusTree->Branch("objectPt",  &objectPiPlusPt,        "objectPt/D");
  piPlusTree->Branch("objectEta", &objectPiPlusEta,       "objectEta/D");
  piPlusTree->Branch("objectPhi", &objectPiPlusPhi,       "objectPhi/D");

  piPlusTree->Branch("objectPiPlusIsChargedHadron", &objectPiPlusIsChargedHadron,"objectPiPlusIsChargedHadron/I");

  piPlusTree->Branch("deltaEta",  &objectPi0DeltaEta,   "deltaEta/D");
  piPlusTree->Branch("deltaPhi",  &objectPi0DeltaPhi,   "deltaPhi/D");


  oneProngTree = fs->make<TTree>("oneProngTree", "Crystal cluster individual crystal pt values");
  oneProngTree->Branch("run",    &run,     "run/I");
  oneProngTree->Branch("lumi",   &lumi,    "lumi/I");
  oneProngTree->Branch("event",  &event,   "event/I");
  
  oneProngTree->Branch("genPt",   &gen1ProngPt,   "genPt/D");
  oneProngTree->Branch("genEta",  &gen1ProngEta,   "genEta/D");
  oneProngTree->Branch("genPhi",  &gen1ProngPhi,   "genPhi/D");
  oneProngTree->Branch("genDM",   &gen1ProngDecayMode,   "genDM/I");


  /*
  oneProngTree->Branch("trackPt",  &oneProngTau.trackPt,   "trackPt/D");
  oneProngTree->Branch("trackEta", &oneProngTau.trackEta,   "trackEta/D");
  oneProngTree->Branch("trackPhi", &oneProngTau.trackPhi,   "trackPhi/D");

  oneProngTree->Branch("objectPt",  &oneProngTau.objectPt,   "objectPt/D");
  oneProngTree->Branch("objectEta", &oneProngTau.objectEta,   "objectEta/D");
  oneProngTree->Branch("objectPhi", &oneProngTau.objectPhi,   "objectPhi/D");

  oneProngTree->Branch("HCALEnergy",  &oneProngTau.HCALEnergy,   "HCALEnergy/D");
  oneProngTree->Branch("ECALEnergy",  &oneProngTau.ECALEnergy,   "ECALEnergy/D");
  
  oneProngTree->Branch("isoRaw",  &oneProngTau.iso,   "isoRaw/D");
  oneProngTree->Branch("decayMode",  &oneProngTau.tauDecayMode,   "decayMode/I");
  */

  oneProngPi0Tree = fs->make<TTree>("oneProngPi0Tree", "Crystal cluster individual crystal pt values");
  oneProngPi0Tree->Branch("run",    &run,     "run/I");
  oneProngPi0Tree->Branch("lumi",   &lumi,    "lumi/I");
  oneProngPi0Tree->Branch("event",  &event,   "event/I");
  
  oneProngPi0Tree->Branch("genPt",   &gen1ProngPi0Pt,   "genPt/D");
  oneProngPi0Tree->Branch("genEta",  &gen1ProngPi0Eta,   "genEta/D");
  oneProngPi0Tree->Branch("genPhi",  &gen1ProngPi0Phi,   "genPhi/D");
  oneProngPi0Tree->Branch("genDM",   &gen1ProngPi0DecayMode,   "genDM/I");

  oneProngPi0Tree->Branch("genPiPlusPt",   &genPiPlusPt,   "genPiPlusPt/D");
  oneProngPi0Tree->Branch("genPiPlusEta",  &genPiPlusEta,   "genPiPlusEta/D");
  oneProngPi0Tree->Branch("genPiPlusPhi",  &genPiPlusPhi,   "genPiPlusPhi/D");

  oneProngPi0Tree->Branch("genPiZeroPt",   &genPiZeroPt,   "genPiZeroPt/D");
  oneProngPi0Tree->Branch("genPiZeroEta",  &genPiZeroEta,   "genPiZeroEta/D");
  oneProngPi0Tree->Branch("genPiZeroPhi",  &genPiZeroPhi,   "genPiZeroPhi/D");


  /*  oneProngPi0Tree->Branch("trackPt",  &oneProngPi0Tau.trackPt,   "trackPt/D");
  oneProngPi0Tree->Branch("trackEta", &oneProngPi0Tau.trackEta,   "trackEta/D");
  oneProngPi0Tree->Branch("trackPhi", &oneProngPi0Tau.trackPhi,   "trackPhi/D");
  
  oneProngPi0Tree->Branch("objectPt",  &oneProngPi0Tau.objectPt,   "objectPt/D");
  oneProngPi0Tree->Branch("objectEta", &oneProngPi0Tau.objectEta,   "objectEta/D");
  oneProngPi0Tree->Branch("objectPhi", &oneProngPi0Tau.objectPhi,   "objectPhi/D");

  oneProngPi0Tree->Branch("stripPt",  &oneProngPi0Tau.stripPt,   "stripPt/D");
  oneProngPi0Tree->Branch("stripEta", &oneProngPi0Tau.stripEta,   "stripEta/D");
  oneProngPi0Tree->Branch("stripPhi", &oneProngPi0Tau.stripPhi,   "stripPhi/D");

  oneProngPi0Tree->Branch("HCALEunergy",  &oneProngPi0Tau.HCALEnergy,   "HCALEnergy/D");
  oneProngPi0Tree->Branch("ECALEnergy",  &oneProngPi0Tau.ECALEnergy,   "ECALEnergy/D");
  
  oneProngPi0Tree->Branch("isoRaw",     &oneProngPi0Tau.iso,   "isoRaw/D");
  oneProngPi0Tree->Branch("decayMode",  &oneProngPi0Tau.tauDecayMode,   "decayMode/I");
  */

  threeProngTree = fs->make<TTree>("threeProngTree", "Crystal cluster individual crystal pt values");
  threeProngTree->Branch("run",    &run,     "run/I");
  threeProngTree->Branch("lumi",   &lumi,    "lumi/I");
  threeProngTree->Branch("event",  &event,   "event/I");
  
  threeProngTree->Branch("genPt",   &gen3ProngPt,   "genPt/D");
  threeProngTree->Branch("genEta",  &gen3ProngEta,   "genEta/D");
  threeProngTree->Branch("genPhi",  &gen3ProngPhi,   "genPhi/D");
  threeProngTree->Branch("genDM",   &gen3ProngDecayMode,   "genDM/I");


  threeProngTree->Branch("track1_pt",  &track1_pt,   "track1_pt/D");
  threeProngTree->Branch("track2_pt",  &track2_pt,   "track2_pt/D");
  threeProngTree->Branch("track3_pt",  &track3_pt,   "track3_pt/D");
  threeProngTree->Branch("other_track_pt",  &other_track_pt,   "other_track_pt/D");

  threeProngTree->Branch("track1_eta",  &track1_eta,   "track1_eta/D");
  threeProngTree->Branch("track2_eta",  &track2_eta,   "track2_eta/D");
  threeProngTree->Branch("track3_eta",  &track3_eta,   "track3_eta/D");  
  threeProngTree->Branch("other_track_eta",  &other_track_eta,   "other_track_eta/D");  
  
  threeProngTree->Branch("track1_phi",  &track1_phi,   "track1_phi/D");  
  threeProngTree->Branch("track2_phi",  &track2_phi,   "track2_phi/D");  
  threeProngTree->Branch("track3_phi",  &track3_phi,   "track3_phi/D");  
  threeProngTree->Branch("other_track_phi",  &other_track_phi,   "other_track_phi/D");  

  threeProngTree->Branch("track2_dR",  &track2_dR,   "track2_dR/D");  
  threeProngTree->Branch("track3_dR",  &track3_dR,   "track3_dR/D");  
  threeProngTree->Branch("other_track_dR",  &other_track_dR,   "other_track_dR/D");  

  threeProngTree->Branch("track2_dz",  &track2_dz,   "track2_dz/D");  
  threeProngTree->Branch("track3_dz",  &track3_dz,   "track3_dz/D");  
  threeProngTree->Branch("other_track_dz",  &other_track_dz,   "other_track_dz/D");  


  /*
  threeProngTree->Branch("objectPt",  &threeProngTau.objectPt,   "objectPt/D");
  threeProngTree->Branch("objectEta", &threeProngTau.objectEta,   "objectEta/D");
  threeProngTree->Branch("objectPhi", &threeProngTau.objectPhi,   "objectPhi/D");


  threeProngTree->Branch("HCALEnergy",  &threeProngTau.HCALEnergy,   "HCALEnergy/D");
  threeProngTree->Branch("ECALEnergy",  &threeProngTau.ECALEnergy,   "ECALEnergy/D");

  threeProngTree->Branch("isoRaw",     &threeProngTau.iso,   "isoRaw/D");
  threeProngTree->Branch("decayMode",  &threeProngTau.tauDecayMode,   "decayMode/I");
  */

  nEvents              = fs->make<TH1F>( "nEvents"  , "nEvents", 2,  0., 1. );
  track_pt             = fs->make<TH1F>( "track_pt"  , "p_{t}", 200,  0., 200. );
  track_pt_eta2p4      = fs->make<TH1F>( "track_pt_eta2p4"  , "p_{t}", 200,  0., 200. );
  track_pt_eta2p1      = fs->make<TH1F>( "track_pt_eta2p1"  , "p_{t}", 200,  0., 200. );

  l1Tau_pt             = fs->make<TH1F>( "l1Tau_ptl"  , "p_{t}", 200,  0., 200. );
  l1Tau_pt_eta2p4      = fs->make<TH1F>( "l1Tau_pt_eta2p4"  , "p_{t}", 200,  0., 200. );
  l1Tau_pt_eta2p1      = fs->make<TH1F>( "l1Tau_pt_eta2p1"  , "p_{t}", 200,  0., 200. );

  l1SingleProngTau_pt            = fs->make<TH1F>( "l1SingleProngTau"  , "p_{t}", 200,  0., 200. );
  l1SingleProngTau_pt_eta2p4     = fs->make<TH1F>( "l1SingleProngTau_eta2p4"  , "p_{t}", 200,  0., 200. );
  l1SingleProngTau_pt_eta2p1     = fs->make<TH1F>( "l1SingleProngTau_eta2p1"  , "p_{t}", 200,  0., 200. );

  l1SingleProngPi0Tau_pt            = fs->make<TH1F>( "l1SingleProngTau"  , "p_{t}", 200,  0., 200. );
  l1SingleProngPi0Tau_pt_eta2p4     = fs->make<TH1F>( "l1SingleProngTau_eta2p4"  , "p_{t}", 200,  0., 200. );
  l1SingleProngPi0Tau_pt_eta2p1     = fs->make<TH1F>( "l1SingleProngTau_eta2p1"  , "p_{t}", 200,  0., 200. );

  l1SingleProngTauIso_pt         = fs->make<TH1F>( "l1SingleProngTauIso"  , "p_{t}", 200,  0., 200. );
  l1SingleProngTauIso_pt_eta2p4  = fs->make<TH1F>( "l1SingleProngTauIso_eta2p4"  , "p_{t}", 200,  0., 200. );
  l1SingleProngTauIso_pt_eta2p1  = fs->make<TH1F>( "l1SingleProngTauIso_eta2p1"  , "p_{t}", 200,  0., 200. );

  l1ThreeProngTau_pt             = fs->make<TH1F>( "l1ThreeProngTau_pt"  , "p_{t}", 200,  0., 200. );
  l1ThreeProngTau_pt_eta2p4      = fs->make<TH1F>( "l1ThreeProngTau_pt_eta2p4"  , "p_{t}", 200,  0., 200. );
  l1ThreeProngTau_pt_eta2p1      = fs->make<TH1F>( "l1ThreeProngTau_pt_eta2p1"  , "p_{t}", 200,  0., 200. );

  l1ThreeProngTauIso_pt          = fs->make<TH1F>( "l1ThreeProngTauIso_pt"  , "p_{t}", 200,  0., 200. );
  l1ThreeProngTauIso_pt_eta2p4   = fs->make<TH1F>( "l1ThreeProngTauIso_pt_eta2p4"  , "p_{t}", 200,  0., 200. );
  l1ThreeProngTauIso_pt_eta2p1   = fs->make<TH1F>( "l1ThreeProngTauIso_pt_eta2p1"  , "p_{t}", 200,  0., 200. );

  //l1OtherTracks    = fs->make<TH2F>( "l1OtherTracks"  , "p_{t}", 40, 0., 0.3, 40, 0., 100. );
  //l1EcalCrystals   = fs->make<TH2F>( "ecal_crystals"  , "p_{t}", 100, -4., 4, 100, -4., 4. );

}


phase2L1TauAnalyzer::~phase2L1TauAnalyzer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)


}


//
// member functions
//

// ------------ method called for each event  ------------
void
phase2L1TauAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  
  run   = iEvent.id().run();
  lumi  = iEvent.id().luminosityBlock();
  event = iEvent.id().event();
  
  //edm::Handle< std::vector<L1CaloCluster> > l1Clusters;
  //iEvent.getByToken( L1ClustersToken_, l1Clusters);
  
  //  edm::Handle< std::vector<L1PFObject> > l1PFChargedCandidates;
  //  iEvent.getByToken( L1PFToken_, l1PFChargedCandidates);
  
  // PF Candidates
  edm::Handle<l1t::PFCandidateCollection> l1PFCandidates;
  iEvent.getByToken(L1PFToken_, l1PFCandidates);
  l1t::PFCandidateCollection pfChargedHadrons;
  l1t::PFCandidateCollection l1PFCandidates_sort;

  for(auto l1PFCand : *l1PFCandidates)
    l1PFCandidates_sort.push_back(l1PFCand);

  std::sort(l1PFCandidates_sort.begin(), l1PFCandidates_sort.end(), [](l1t::PFCandidate i,l1t::PFCandidate j){return(i.pt() > j.pt());});   

  // get all PF Charged Hadrons
  for(auto l1PFCand : l1PFCandidates_sort)
    {
      // this saves some low pt 3 prong taus, but should be optimized with respect to 1 prong pi0
      // Commenting out this requirement for now:
      //      if(l1PFCand.id() == l1t::PFCandidate::ChargedHadron || (l1PFCand.pt()<5 && l1PFCand.id() == l1t::PFCandidate::Electron)){
      pfChargedHadrons.push_back(l1PFCand);
      continue;
      
    }
  // end of PF Candidates

  edm::Handle<L1TkPrimaryVertexCollection> L1VertexHandle;
  iEvent.getByToken(pvToken_, L1VertexHandle);

  edm::Handle< std::vector<L1PFTau> > l1PFTaus;
  iEvent.getByToken( L1PFTauToken_, l1PFTaus);

  Handle<std::vector<pat::PackedCandidate> > packedcands;
  iEvent.getByToken(PackedCands_, packedcands);
  
  edm::Handle< std::vector<pat::Tau> > miniTaus;
  if(!iEvent.getByToken( MiniTausToken_, miniTaus))
    std::cout<<"No miniAOD particles found"<<std::endl;
  
  //control plot for ecal crystals
  edm::Handle<EcalEBTrigPrimDigiCollection> ecaltpgCollection;
  iEvent.getByToken( ecalTPGBToken_, ecaltpgCollection);

  edm::ESHandle<CaloGeometry> caloGeometryHandle;
  iSetup.get<CaloGeometryRecord>().get(caloGeometryHandle);
  const CaloGeometry* caloGeometry_ = caloGeometryHandle.product();  
  ebGeometry = caloGeometry_->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);

  /*
  for(auto& tpg : *ecaltpgCollection.product())
    {

      if(tpg.encodedEt() > 0) 
      {

        GlobalVector position;
	  auto cell = ebGeometry->getGeometry(tpg.id());

	    float et = tpg.encodedEt()/8.;

	      if(et<0.001) continue;//
	        //float energy = et / sin(position.theta());
		  float eta = cell->getPosition().eta();
		    float phi = cell->getPosition().phi();
		      //l1EcalCrystals->Fill(eta,phi,et);
		      }
    }
  */
  // L1 tracks  
  std::vector<TTTrack< Ref_Phase2TrackerDigi_ > > l1Tracks;
  edm::Handle< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > l1trackHandle;
  iEvent.getByToken(ttTrackToken_, l1trackHandle);
  l1Tracks.clear();

  //Find and sort the tracks
  for(size_t track_index=0; track_index<l1trackHandle->size(); ++track_index)
    {

      edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > ptr(l1trackHandle, track_index);
      //double pt  = ptr->getMomentum().perp();
      //double eta = ptr->getMomentum().eta();

      //only using tracks with eta less than 1.5 and pt greater than 2.5 GeV
      //if(abs(eta)<1.5 && pt > 2.5)
      l1Tracks.push_back(l1trackHandle->at(track_index));       

    }

  std::sort(l1Tracks.begin(), l1Tracks.end(), [](TTTrack< Ref_Phase2TrackerDigi_ > i,TTTrack< Ref_Phase2TrackerDigi_ > j){return(i.getMomentum().perp() > j.getMomentum().perp());});   

  if(l1Tracks.size()>0)
    track_pt->Fill(l1Tracks.at(0).getMomentum().perp());

  for(auto l1Track: l1Tracks){
    if(abs(l1Track.getMomentum().eta())<2.1){
      track_pt_eta2p1->Fill(l1Track.getMomentum().perp());
      break;
    }
  }

  for(auto l1Track: l1Tracks){
    if(abs(l1Track.getMomentum().eta())<2.4){
      track_pt_eta2p4->Fill(l1Track.getMomentum().perp());
      break;
    }
  }
  

  std::vector<L1PFTau> l1PFTaus_sorted;

  for(unsigned int i = 0; i < l1PFTaus->size(); i++){
    l1PFTaus_sorted.push_back(l1PFTaus->at(i));
  }
  
  std::sort(l1PFTaus_sorted.begin(), l1PFTaus_sorted.end(), [](L1PFTau i,L1PFTau j){return(i.p4().pt() > j.p4().pt());});   


  //filling general rate tree
  if(l1PFTaus_sorted.size()>0)
    l1Tau_pt->Fill(l1PFTaus_sorted.at(0).pt());
  
  bool f_2p1=false;
  bool f_2p4=false;

  bool f0=false;
  bool f0_2p1=false;
  bool f0_2p4=false;

  bool f1=false;
  bool f1_2p1=false;
  bool f1_2p4=false;

  bool f10 = false;
  bool f10_2p1=false;
  bool f10_2p4=false;

  for(unsigned int i = 0; i < l1PFTaus_sorted.size(); i++){
    
    if(abs(l1PFTaus->at(i).eta())<2.1 && !f_2p1){
      l1Tau_pt_eta2p1->Fill(l1PFTaus->at(i).pt());
      f_2p1=true;
    }

    if(abs(l1PFTaus->at(i).eta())<2.4 && !f_2p4){
      l1Tau_pt_eta2p4->Fill(l1PFTaus->at(i).pt());
      f_2p4=true;
    }

    if(l1PFTaus->at(i).tauType()==0){
      
      if(!f0){
	l1SingleProngTau_pt->Fill(l1PFTaus->at(i).pt());
	f0 = true;
      }

      if(abs(l1PFTaus->at(i).eta())<2.1 && !f0_2p1){
        l1SingleProngTau_pt_eta2p1->Fill(l1PFTaus->at(i).pt());
	f0_2p1=true;
      }
      if(abs(l1PFTaus->at(i).eta())<2.4 && !f0_2p4){
	l1SingleProngTau_pt_eta2p4->Fill(l1PFTaus->at(i).pt());
	f0_2p4=true;
      }
      
    }//close 1prong

    if(l1PFTaus->at(i).tauType()==1){
      
      if(!f1){
	l1SingleProngPi0Tau_pt->Fill(l1PFTaus->at(i).pt());
	f1=true;
      }

      if(abs(l1PFTaus->at(i).eta())<2.1 && !f1_2p1){
        l1SingleProngPi0Tau_pt_eta2p1->Fill(l1PFTaus->at(i).pt());
	f1_2p1=true;
      }
      if(abs(l1PFTaus->at(i).eta())<2.4 && !f1_2p4){
	l1SingleProngPi0Tau_pt_eta2p4->Fill(l1PFTaus->at(i).pt());
	f1_2p4=true;
      }
      
    }//close 1prongPi0

    if(l1PFTaus->at(i).tauType()==10){
      if(!f10){
	l1ThreeProngTau_pt->Fill(l1PFTaus->at(i).pt());
	f10=true;
      }
      if(abs(l1PFTaus->at(i).eta())<2.1 && !f10_2p1){
        l1ThreeProngTau_pt_eta2p1->Fill(l1PFTaus->at(i).pt());
	f10_2p1=true;
      }
      if(abs(l1PFTaus->at(i).eta())<2.4 && !f10_2p4){
	l1ThreeProngTau_pt_eta2p4->Fill(l1PFTaus->at(i).pt());
	f10_2p4=true;
      }
      
    }//close 3prong


  }


  // Get genParticles
  edm::Handle<GenParticleCollectionType> genParticleHandle;
  if(!iEvent.getByToken(genToken_,genParticleHandle))
    std::cout<<"No gen Particles Found "<<std::endl;
  else
    std::cout<<"Gen Particles size "<<genParticleHandle->size()<<std::endl;

  std::vector<reco::GenParticle> genTaus;
  std::vector<reco::GenParticle> genPiZeros;
  std::vector<reco::GenParticle> genPiPluss;
  std::vector<reco::GenParticle> genParticles;

  //reco::GenParticle* genTau;
  for(unsigned int i = 0; i< genParticleHandle->size(); i++){
    edm::Ptr<reco::GenParticle> ptr(genParticleHandle, i);
    genParticles.push_back(*ptr);

    if(abs(ptr->pdgId())==111 && abs(ptr->eta()<1.74)){
      genPiZeros.push_back(*ptr);
      //std::cout<<"Found PiZero PDGID 111 pt: "<<ptr->pt()<<" eta: "<<ptr->eta()<<" phi: "<<ptr->phi()<<std::endl;
    }
    if(abs(ptr->pdgId())==211 && abs(ptr->eta()<1.74)){
      genPiPluss.push_back(*ptr);
      //std::cout<<"Found PiPlus PDGID 111 pt: "<<ptr->pt()<<" eta: "<<ptr->eta()<<" phi: "<<ptr->phi()<<std::endl;
    }
    if(abs(ptr->pdgId())==15){
      genTaus.push_back(*ptr);
    }
  }

  std::vector<genVisTau> GenOneProngTaus;
  std::vector<genVisTau> GenOneProngPi0Taus;
  std::vector<genVisTau> GenThreeProngTaus;
  //Find and Sort the 1 Prong, 1 Prong + pi0 and 3 Prong Taus

  //std::cout<<"starting gen taus"<<std::endl;

  for(auto genTau: genTaus){

    reco::Candidate::LorentzVector visGenTau= getVisMomentum(&genTau, &genParticles);
    genVisTau Temp;
    genPt = visGenTau.pt();
    genEta = visGenTau.eta();
    genPhi = visGenTau.phi();
    decayMode = GetDecayMode(&genTau);
    Temp.p4 = visGenTau;
    Temp.decayMode = decayMode;

    //std::cout<<"Tau Decay Mode "<<decayMode<<std::endl;
    //std::cout<<"tau vis pt: "<<genPt<<" genEta: "<<genEta<<" genPhi: "<<genPhi<<std::endl;

    if(decayMode >21 ){
      std::cout<<"found 3 prong tau: "<<decayMode<<std::endl;
      GenThreeProngTaus.push_back(Temp);
    }

    if(decayMode == 10 ){
      GenOneProngTaus.push_back(Temp);
    }

    if(decayMode > 10 && decayMode < 20 ){
      GenOneProngPi0Taus.push_back(Temp);
    }
  }

  for(unsigned int i = 0; i < miniTaus->size(); i++){
    recoEta        = -99;
    recoPhi        = -99;
    recoPt         = -99;
    recoChargedIso = -99;
    recoNeutralIso = -99;
    recoRawIso     = -99;
    recoDecayMode  = -99;
     
    if(miniTaus->at(i).tauID("decayModeFinding")>0){
      recoEta        = miniTaus->at(i).p4().Eta();
      recoPhi        = miniTaus->at(i).p4().Phi();
      recoPt         = miniTaus->at(i).p4().Pt();
      recoChargedIso = miniTaus->at(i).tauID("chargedIsoPtSum");
      recoNeutralIso = miniTaus->at(i).tauID("neutralIsoPtSum");
      recoRawIso     = miniTaus->at(i).tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
      recoDecayMode  = miniTaus->at(i).decayMode();
      std::cout<<"=== Found recoTau: "<<recoPt<<" Eta: "<<recoEta<<" Phi: "<< recoPhi <<std::endl;
    }
    else{
      continue;
    }
        
    // Get the L1 track that is closest to the recoTau.
    l1Track_pt = 0;
    l1Track_eta = -10;
    l1Track_phi = -10;
    l1Track_dz = 0;
    l1Track_dR = 99.99;
    l1Track_dRmin = 99.99;
    l1TrackRecoTauPtDiff = 99.99;

    for(auto l1Track: l1Tracks){
      // Calculate dR with respect to the recoTau.
      l1Track_dR = reco::deltaR(l1Track.getMomentum().eta(), l1Track.getMomentum().phi(),
				recoEta, recoPhi);

      if ((l1Track_dR < l1Track_dRmin)
	  && (l1Track_dR < 0.3))
	{
	  l1Track_pt = l1Track.getMomentum().perp();
	  l1Track_eta = l1Track.getMomentum().eta();
	  l1Track_phi = l1Track.getMomentum().phi();
	  l1Track_dz = l1Track.getPOCA().z() - L1VertexHandle->at(0).getZvertex();
	  l1Track_dRmin = l1Track_dR;
	  l1TrackRecoTauPtDiff = recoPt - l1Track_pt;
	}
    }

    printf("l1Track_dR is %f for the closest l1 Track to the recoTau.\n", l1Track_dR);
    printf("l1TrackRecoTauPtDiff is %f for the closest l1 Track (l1Track_pt = %f) to the recoTau with recoPt %f\n",
	   l1TrackRecoTauPtDiff, l1Track_pt, recoPt);

    // Get the PF Cand that is closest to the recoTau.
    pfCand_pt = 0;
    pfCand_eta = -10;
    pfCand_phi = -10;
    pfCand_dz = 0;
    pfCand_dR = 99.99;
    pfCand_dRmin = 99.99;
    
    for(auto l1PFCand : pfChargedHadrons)
      {
	// Calculate dR with respect to the recoTau.                                                                                             
	pfCand_dR = reco::deltaR(l1PFCand.p4().eta(), l1PFCand.p4().phi(),
				 recoEta, recoPhi);
	if ((pfCand_dR < pfCand_dRmin)
	    && (pfCand_dR < 0.3))
	  {
	    pfCand_pt = l1PFCand.pt();
	    pfCand_eta = l1PFCand.p4().eta();
	    pfCand_phi = l1PFCand.p4().phi();
	    pfCand_dz = l1PFCand.vz() - L1VertexHandle->at(0).getZvertex();
	    pfCand_dRmin = pfCand_dR;

	  }
      }
    
    printf("pfCand_dR is %f for the closest PF Cand to recoTau.\n", pfCand_dR);
    printf("pfCand_pt is %f, recoPt is %f (difference is %f).\n", pfCand_pt, recoPt, pfCand_pt - recoPt);
     
    l1Pt = 0;
    l1Eta = -10;
    l1Phi = -10;
    l1DM = -10;
    l1TauRawIso = 100;
    l1TauNeutralIso = 100;
    l1TauChargedIso = 100;
    l1IsoVLoose = -10;
    l1IsoLoose = -10;
    l1IsoMedium = -10;
    l1IsoTight = -10;

    zVTX = 0;
    l1TauZ = 0;
    l1PVDZ = 0;  // delta Z = tau's z minus zVTX
     
    l1RelIsoVLoose = -10;
    l1RelIsoLoose = -10;
    l1RelIsoMedium = -10;
    l1RelIsoTight = -10;

    std::cout << "l1PFTaus size: " << l1PFTaus->size() << std::endl;

    // Set primary vertex z position
    if (L1VertexHandle->size() > 0)
      {
	zVTX = L1VertexHandle->at(0).getZvertex();
      }

    for(unsigned int i = 0; i < l1PFTaus->size(); i++){
      /*       std::cout << "l1PFTaus->at(i).p4().Eta() = " << l1PFTaus->at(i).p4().Eta() << "   "
	       << "l1PFTaus->at(i).p4().Phi() = " << l1PFTaus->at(i).p4().Phi() << "   " 
	       << "recoEta = " << recoEta << "   " 
	       << "recoPhi = " << recoPhi << "   "
	       << "l1PFTaus->at(i).p4().pt() = " << l1PFTaus->at(i).p4().pt() << "   " 
	       << "l1TauPt  " << l1TauPt << std::endl;*/
       
      if(( reco::deltaR(l1PFTaus->at(i).eta(), l1PFTaus->at(i).phi(), 
			recoEta, recoPhi) < 0.5 )
	 && (l1PFTaus->at(i).pt() > l1Pt))
	{
	  l1Eta = l1PFTaus->at(i).eta();
	  l1Phi = l1PFTaus->at(i).phi();
	  l1Pt  = l1PFTaus->at(i).pt();

	  l1TauZ = l1PFTaus->at(i).p4().z();   // adding in z position of tau

	  l1StripPt  = l1PFTaus->at(i).strip_p4().pt();
	  l1StripEta = l1PFTaus->at(i).strip_p4().eta();
	  l1StripPhi = l1PFTaus->at(i).strip_p4().phi();
	  l1StripDR  = reco::deltaR(l1PFTaus->at(i).eta(),
				    l1PFTaus->at(i).phi(),
				    l1PFTaus->at(i).strip_p4().eta(),
				    l1PFTaus->at(i).strip_p4().phi());
	  if (L1VertexHandle->size() > 0)
	    {
	      l1PVDZ = l1TauZ - zVTX;
	      std::cout << "l1TauPVDZ (without reco/l1 matching): " << l1PVDZ << std::endl;
	    }


	  l1TauChargedIso = l1PFTaus->at(i).chargedIso();
	  l1TauNeutralIso = l1PFTaus->at(i).neutralIso();
	  l1TauRawIso = l1PFTaus->at(i).rawIso();
	  l1DM        = l1PFTaus->at(i).tauType();
	  l1IsoVLoose =  l1PFTaus->at(i).passVLooseIso();
	  l1IsoLoose  =  l1PFTaus->at(i).passLooseIso();
	  l1IsoMedium =  l1PFTaus->at(i).passMediumIso();
	  l1IsoTight  =  l1PFTaus->at(i).passTightIso();

	  l1RelIsoVLoose =  l1PFTaus->at(i).passVLooseRelIso();
	  l1RelIsoLoose  =  l1PFTaus->at(i).passLooseRelIso();
	  l1RelIsoMedium =  l1PFTaus->at(i).passMediumRelIso();
	  l1RelIsoTight  =  l1PFTaus->at(i).passTightRelIso();

	     
	  std::cout<<" Match found l1Pt: "<< l1Pt <<
	    " Eta: "<< l1Eta  <<
	    " Phi: "<< l1Phi  << 
	    " Pass tight iso: " << l1PFTaus->at(i).passTightIso() <<
	    " Pass VLoose Iso: " << l1PFTaus->at(i).passVLooseIso() <<
	    std::endl;
	}
       
    } // end of loop over L1 taus
   
     
    genPt = 0;
    genEta = -100;
    genPhi = -100;
     
    for(auto genTau: genTaus){
       
      reco::Candidate::LorentzVector visGenTau= getVisMomentum(&genTau, &genParticles);
      genVisTau Temp;
       
      if( reco::deltaR(recoEta, 
		       recoPhi, 
		       visGenTau.eta(), visGenTau.phi()) < 0.5){
	genPt = visGenTau.pt();
	genEta = visGenTau.eta();
	genPhi = visGenTau.phi();
	decayMode = GetDecayMode(&genTau);
	std::cout<<"    tau vis pt: "<<genPt<<" genEta: "<<genEta<<" genPhi: "<<genPhi<<std::endl;
      }
    }
       
    efficiencyTree->Fill();

  } // end of loop over reco (miniTaus)
   
  for(auto genTau: GenThreeProngTaus){
    genPt  = genTau.p4.pt();
    genEta = genTau.p4.eta();
    genPhi = genTau.p4.phi();
     
    for(unsigned int i = 0; i < miniTaus->size(); i++){
      if( reco::deltaR(miniTaus->at(i).p4().Eta(), 
		       miniTaus->at(i).p4().Phi(), 
		       genEta, genPhi) < 0.5)
	if(miniTaus->at(i).tauID("decayModeFinding")>0 && miniTaus->at(i).decayMode()==10){
	  recoEta        = miniTaus->at(i).p4().Eta();
	  recoPhi        = miniTaus->at(i).p4().Phi();
	  recoPt         = miniTaus->at(i).p4().Pt();
	  recoDecayMode  = miniTaus->at(i).decayMode();
	  //std::cout<<"accessing the signal tracks size"<<std::endl;
	  //std::cout<<"size: "<<miniTaus->at(i).signalChargedHadrCands().size()<<std::endl;
	  const reco::CandidatePtrVector chargedHadrons = miniTaus->at(i).signalChargedHadrCands();
	  //std::cout<<"size CH: "<< (*chargedHadrons)->size()<<std::endl;
	  //for (size_t j = 0; i < .size(); ++i) {
	  //if(miniTaus->at(i).signalChargedHadrCands().size()>2){
	  int j = 0;

	  pat::PackedCandidate const* pfCand_1 = 0;
	  pat::PackedCandidate const* pfCand_2 = 0;
	  pat::PackedCandidate const* pfCand_3 = 0;

	  math::PtEtaPhiMLorentzVector tempP4_track1;
	  math::PtEtaPhiMLorentzVector tempP4_track2;
	  math::PtEtaPhiMLorentzVector tempP4_track3;
	     
	  for(reco::CandidatePtrVector::const_iterator iter = chargedHadrons.begin(); iter != chargedHadrons.end(); iter++){
	    //pat::PackedCandidate const* packedCand = dynamic_cast<pat::PackedCandidate const*>(iter->get());
	       
	    if(j==0){
	      pfCand_1 = dynamic_cast<pat::PackedCandidate const*>(iter->get());
	      std::cout<<"track 1 pt: "<<(pfCand_1)->pt()<<" eta: "<<(pfCand_1)->eta()<<" phi: "<<(pfCand_1)->phi()<<std::endl;
	      tempP4_track1.SetPt((pfCand_1)->pt()); 
	      tempP4_track1.SetEta((pfCand_1)->eta()); 
	      tempP4_track1.SetPhi((pfCand_1)->phi());

	      track1_pt      = (pfCand_1)->pt();
	      track1_eta     = (pfCand_1)->eta();
	      track1_phi     = (pfCand_1)->phi();
	    }

	    if(j==1){
	      pfCand_2 = dynamic_cast<pat::PackedCandidate const*>(iter->get());
	      //std::cout<<"track 2 pt: "<<track2_pt<<" eta: "<<track2_eta<<" phi: "<<track2_phi<<std::endl;
	      std::cout<<"track 2 pt: "<<(pfCand_2)->pt()<<" eta: "<< (pfCand_2)->eta() <<" phi: "<< (pfCand_2)->phi() <<std::endl;
	      tempP4_track2.SetPt((pfCand_2)->pt()); 
	      tempP4_track2.SetEta((pfCand_2)->eta()); 
	      tempP4_track2.SetPhi((pfCand_2)->phi());
	      //tempP4_track2.SetPtEtaPhiM((pfCand_2)->pt(), (pfCand_2)->eta(), (pfCand_2)->phi(), 1.335);

	      track2_pt      = (pfCand_2)->pt();
	      track2_eta     = (pfCand_2)->eta();
	      track2_phi     = (pfCand_2)->phi();

	      track2_dR      = ROOT::Math::VectorUtil::DeltaR(tempP4_track1,tempP4_track2);
	      track2_dz      = abs((pfCand_1)->pz() - (pfCand_2)->pz());

	    }

	    if(j==2){
	      pfCand_3 = dynamic_cast<pat::PackedCandidate const*>(iter->get());
	      std::cout<<"track 3 pt: "<<(pfCand_3)->pt()<<" eta: "<< (pfCand_3)->eta() <<" phi: "<< (pfCand_3)->phi() <<std::endl;
	      tempP4_track3.SetPt((pfCand_3)->pt()); 
	      tempP4_track3.SetEta((pfCand_3)->eta()); 
	      tempP4_track3.SetPhi((pfCand_3)->phi());

	      //tempP4_track3.SetPtEtaPhiM((pfCand_3)->pt(), (pfCand_3)->eta(), (pfCand_3)->phi(), 1.335);
	      track3_pt      = (pfCand_3)->pt();
	      track3_eta     = (pfCand_3)->eta();
	      track3_phi     = (pfCand_3)->phi();

	      track3_dR      = ROOT::Math::VectorUtil::DeltaR(tempP4_track1,tempP4_track3);
	      track3_dz      = abs((pfCand_1)->pz() - (pfCand_3)->pz());
	      
	    }
	    j++;
	  }
	     

	  float highestPt_other_track = 0;
	  for(unsigned int i=0;i<packedcands->size();i++){
	    bool isASignalCand = false;
	    const pat::PackedCandidate & c = (*packedcands)[i];

	    //check if it is nearby
	    if( reco::deltaR(track1_eta, track1_phi, c.eta(), c.phi()) < 0.3){
	      for(reco::CandidatePtrVector::const_iterator iter = chargedHadrons.begin(); iter != chargedHadrons.end(); iter++){
		pat::PackedCandidate const* signalcand = dynamic_cast<pat::PackedCandidate const*>(iter->get());
		//std::cout<<"c.sourceCandidatePtr(0).key() "<<c.sourceCandidatePtr(0).key()<<" signalcand->sourceCandidatePtr(0).key() "<<signalcand->sourceCandidatePtr(0).key()<<std::endl;
		if(c.pt() == signalcand->pt()){
		  //std::cout<<"MATCHED c.sourceCandidatePtr(0).key() "<<c.sourceCandidatePtr(0).key()<<" signalcand->sourceCandidatePtr(0).key() "<<signalcand->sourceCandidatePtr(0).key()<<std::endl;
		  isASignalCand = true;
		}
	      }

	      //if(!isASignalCand)
	      //l1OtherTracks->Fill(c.pt(),ROOT::Math::VectorUtil::DeltaR(tempP4_track1,c.p4()));

	      if(!isASignalCand && c.pt() > highestPt_other_track){
		other_track_pt = c.pt();
		other_track_eta = c.eta();
		other_track_phi = c.phi();
		other_track_dR  = ROOT::Math::VectorUtil::DeltaR(tempP4_track1,c.p4());
		other_track_dz  = abs((pfCand_1)->pz() - c.pz());
	      }
	    }
	  }
	  std::cout<<"Filling the tree"<<std::endl;
	  threeProngTree->Fill();
	}
    }
     
  }
   
  std::cout<<"Finished Analyzing the Taus"<<std::endl;
}


// ------------ method called once each job just before starting event loop  ------------
void 
phase2L1TauAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
phase2L1TauAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
phase2L1TauAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(phase2L1TauAnalyzer);
