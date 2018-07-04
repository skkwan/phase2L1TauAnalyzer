// -*- C++ -*-
//
// Package:    L1Trigger/phase2L1TauAnalyzerRates
// Class:      phase2L1TauAnalyzerRates
// 
/**\class phase2L1TauAnalyzerRates phase2L1TauAnalyzerRates.cc L1Trigger/phase2L1TauAnalyzerRates/plugins/phase2L1TauAnalyzerRates.cc

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

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"


#include "L1Trigger/phase2Demonstrator/interface/triggerGeometryTools.hh"
#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"

#include "L1Trigger/phase2L1TauAnalyzer/plugins/helpers.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/L1Trigger/interface/L1PFTau.h"

#include "DataFormats/Phase2L1CaloTrig/interface/L1CaloCluster.h"
#include "DataFormats/L1Trigger/interface/L1PFObject.h"

#include "DataFormats/PatCandidates/interface/Tau.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

using namespace l1t;

class phase2L1TauAnalyzerRates : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit phase2L1TauAnalyzerRates(const edm::ParameterSet&);
      ~phase2L1TauAnalyzerRates();

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

  edm::EDGetTokenT< L1CaloClusterCollection > L1ClustersToken_;
  edm::EDGetTokenT< L1PFObjectCollection > L1PFToken_;
  edm::EDGetTokenT< L1PFTauCollection > L1PFTauToken_;
  edm::EDGetTokenT<std::vector<reco::GenParticle> > genToken_;
  edm::EDGetTokenT< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > ttTrackToken_;
  edm::EDGetTokenT<EcalEBTrigPrimDigiCollection> ecalTPGBToken_;
  edm::InputTag genSrc_;
  edm::InputTag L1TrackInputTag;

  const CaloSubdetectorGeometry * ebGeometry;
  edm::ESHandle<CaloGeometry> caloGeometry_;

  TTree* pi0Tree;
  TTree* piPlusTree;

  TTree* oneProngTree;
  TTree* oneProngPi0Tree;
  TTree* threeProngTree;
  
  TTree* efficiencyTree;

  double genPt, genEta, genPhi;
  int decayMode, run, lumi, event;
  double l1TauPt, l1TauEta, l1TauPhi, l1TauDecayMode;
  double l1TauChargedIso,l1TauNeutralIso,l1TauRawIso;
  double gen1ProngPt, gen1ProngEta, gen1ProngPhi;
  double gen3ProngPt, gen3ProngEta, gen3ProngPhi;
  double gen1ProngPi0Pt, gen1ProngPi0Eta, gen1ProngPi0Phi;
  double genPiZeroPt, genPiZeroEta, genPiZeroPhi;

  double recoPt, recoEta, recoPhi; 
  float recoChargedIso, recoNeutralIso, recoRawIso; 
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

  TH1F* l1TauIso_pt;     
  TH1F* l1TauIso_pt_eta2p4;     
  TH1F* l1TauIso_pt_eta2p1;

  TH1F* l1SingleProngTau_pt;     
  TH1F* l1SingleProngTau_pt_eta2p4;     
  TH1F* l1SingleProngTau_pt_eta2p1;

  TH1F* l1SingleProngTauIso_pt;     
  TH1F* l1SingleProngTauIso_pt_eta2p4;     
  TH1F* l1SingleProngTauIso_pt_eta2p1;

  TH1F* l1SingleProngPi0Tau_pt;     
  TH1F* l1SingleProngPi0Tau_pt_eta2p4;     
  TH1F* l1SingleProngPi0Tau_pt_eta2p1;

  TH1F* l1SingleProngPi0TauIso_pt;     
  TH1F* l1SingleProngPi0TauIso_pt_eta2p4;     
  TH1F* l1SingleProngPi0TauIso_pt_eta2p1;

  TH1F* l1ThreeProngTau_pt;     
  TH1F* l1ThreeProngTau_pt_eta2p4;     
  TH1F* l1ThreeProngTau_pt_eta2p1;

  TH1F* l1ThreeProngTauIso_pt;     
  TH1F* l1ThreeProngTauIso_pt_eta2p4;     
  TH1F* l1ThreeProngTauIso_pt_eta2p1;

  TH1F* l1DiTau_pt;     
  TH1F* l1DiTau_pt_eta2p4;     
  TH1F* l1DiTau_pt_eta2p1;

  TH1F* l1DiTauIso_pt;     
  TH1F* l1DiTauIso_pt_eta2p4;     
  TH1F* l1DiTauIso_pt_eta2p1;

  TH2F* l1EcalCrystals;
  TH2F* l1EcalCrystals_2;
  TH2F* l1EcalCrystals_5;
  TH2F* l1EcalCrystals_10;
  TH2F* l1EcalCrystals_15;
  TH2F* l1EcalCrystals_20;
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
phase2L1TauAnalyzerRates::phase2L1TauAnalyzerRates(const edm::ParameterSet& cfg):
  L1ClustersToken_( consumes< L1CaloClusterCollection >(cfg.getParameter<edm::InputTag>("L1Clusters"))),
  L1PFToken_(       consumes< L1PFObjectCollection >(cfg.getParameter<edm::InputTag>("l1PFObjects"))),
  L1PFTauToken_(    consumes< L1PFTauCollection    >(cfg.getParameter<edm::InputTag>("l1TauObjects"))),
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


  /*
  efficiencyTree->Branch("trackPt",  &trackPt,   "trackPt/D");
  efficiencyTree->Branch("trackEta", &trackEta,   "trackEta/D");
  efficiencyTree->Branch("trackPhi", &trackPhi,   "trackPhi/D");
  */
  efficiencyTree->Branch("l1TauPt",  &l1TauPt,   "l1TauPt/D");
  efficiencyTree->Branch("l1TauEta", &l1TauEta,   "l1TauEta/D");
  efficiencyTree->Branch("l1TauPhi", &l1TauPhi,   "l1TauPhi/D");

  efficiencyTree->Branch("l1RawIso",     &l1TauRawIso,     "l1RawIso/D");
  efficiencyTree->Branch("l1ChargedIso", &l1TauChargedIso, "l1ChargedIso/D");
  efficiencyTree->Branch("l1NeutralIso", &l1TauNeutralIso, "l1NeutralIso/D");

  efficiencyTree->Branch("l1TauDecayMode", &l1TauDecayMode,   "l1TauDecayMode/D");

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

  oneProngPi0Tree->Branch("HCALEnergy",  &oneProngPi0Tau.HCALEnergy,   "HCALEnergy/D");
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

  /*
  threeProngTree->Branch("objectPt",  &threeProngTau.objectPt,   "objectPt/D");
  threeProngTree->Branch("objectEta", &threeProngTau.objectEta,   "objectEta/D");
  threeProngTree->Branch("objectPhi", &threeProngTau.objectPhi,   "objectPhi/D");


  threeProngTree->Branch("HCALEnergy",  &threeProngTau.HCALEnergy,   "HCALEnergy/D");
  threeProngTree->Branch("ECALEnergy",  &threeProngTau.ECALEnergy,   "ECALEnergy/D");

  threeProngTree->Branch("isoRaw",     &threeProngTau.iso,   "isoRaw/D");
  threeProngTree->Branch("decayMode",  &threeProngTau.tauDecayMode,   "decayMode/I");
  */

  nEvents              = fs->make<TH1F>( "nEvents"  , "nEvents", 2,  0., 2. );
  track_pt             = fs->make<TH1F>( "track_pt"  , "p_{t}", 300,  0., 300. );
  track_pt_eta2p4      = fs->make<TH1F>( "track_pt_eta2p4"  , "p_{t}", 300,  0., 300. );
  track_pt_eta2p1      = fs->make<TH1F>( "track_pt_eta2p1"  , "p_{t}", 300,  0., 300. );

  l1Tau_pt             = fs->make<TH1F>( "l1Tau_pt"  , "p_{t}", 300,  0., 300. );
  l1Tau_pt_eta2p4      = fs->make<TH1F>( "l1Tau_pt_eta2p4"  , "p_{t}", 300,  0., 300. );
  l1Tau_pt_eta2p1      = fs->make<TH1F>( "l1Tau_pt_eta2p1"  , "p_{t}", 300,  0., 300. );

  l1TauIso_pt            = fs->make<TH1F>( "l1TauIso_pt"  , "p_{t}", 300,  0., 300. );
  l1TauIso_pt_eta2p4     = fs->make<TH1F>( "l1TauIso_pt_eta2p4"  , "p_{t}", 300,  0., 300. );
  l1TauIso_pt_eta2p1     = fs->make<TH1F>( "l1TauIso_pt_eta2p1"  , "p_{t}", 300,  0., 300. );

  l1DiTau_pt             = fs->make<TH1F>( "l1DiTau_pt"  , "p_{t}", 300,  0., 300. );
  l1DiTau_pt_eta2p4      = fs->make<TH1F>( "l1DiTau_pt_eta2p4"  , "p_{t}", 300,  0., 300. );
  l1DiTau_pt_eta2p1      = fs->make<TH1F>( "l1DiTau_pt_eta2p1"  , "p_{t}", 300,  0., 300. );

  l1DiTauIso_pt          = fs->make<TH1F>( "l1DiTauIso_pt"  , "p_{t}", 300,  0., 300. );
  l1DiTauIso_pt_eta2p4   = fs->make<TH1F>( "l1DiTauIso_pt_eta2p4"  , "p_{t}", 300,  0., 300. );
  l1DiTauIso_pt_eta2p1   = fs->make<TH1F>( "l1DiTauIso_pt_eta2p1"  , "p_{t}", 300,  0., 300. );

  l1SingleProngTau_pt            = fs->make<TH1F>( "l1SingleProngTau"  , "p_{t}", 300,  0., 300. );
  l1SingleProngTau_pt_eta2p4     = fs->make<TH1F>( "l1SingleProngTau_eta2p4"  , "p_{t}", 300,  0., 300. );
  l1SingleProngTau_pt_eta2p1     = fs->make<TH1F>( "l1SingleProngTau_eta2p1"  , "p_{t}", 300,  0., 300. );

  l1SingleProngTauIso_pt         = fs->make<TH1F>( "l1SingleProngTauIso"  , "p_{t}", 300,  0., 300. );
  l1SingleProngTauIso_pt_eta2p4  = fs->make<TH1F>( "l1SingleProngTauIso_eta2p4"  , "p_{t}", 300,  0., 300. );
  l1SingleProngTauIso_pt_eta2p1  = fs->make<TH1F>( "l1SingleProngTauIso_eta2p1"  , "p_{t}", 300,  0., 300. );

  l1SingleProngPi0Tau_pt            = fs->make<TH1F>( "l1SingleProngPi0Tau"  , "p_{t}", 300,  0., 300. );
  l1SingleProngPi0Tau_pt_eta2p4     = fs->make<TH1F>( "l1SingleProngPi0Tau_eta2p4"  , "p_{t}", 300,  0., 300. );
  l1SingleProngPi0Tau_pt_eta2p1     = fs->make<TH1F>( "l1SingleProngPi0Tau_eta2p1"  , "p_{t}", 300,  0., 300. );

  l1SingleProngPi0TauIso_pt         = fs->make<TH1F>( "l1SingleProngPi0TauIso"  , "p_{t}", 300,  0., 300. );
  l1SingleProngPi0TauIso_pt_eta2p4  = fs->make<TH1F>( "l1SingleProngPi0TauIso_eta2p4"  , "p_{t}", 300,  0., 300. );
  l1SingleProngPi0TauIso_pt_eta2p1  = fs->make<TH1F>( "l1SingleProngPi0TauIso_eta2p1"  , "p_{t}", 300,  0., 300. );

  l1ThreeProngTau_pt             = fs->make<TH1F>( "l1ThreeProngTau_pt"  , "p_{t}", 300,  0., 300. );
  l1ThreeProngTau_pt_eta2p4      = fs->make<TH1F>( "l1ThreeProngTau_pt_eta2p4"  , "p_{t}", 300,  0., 300. );
  l1ThreeProngTau_pt_eta2p1      = fs->make<TH1F>( "l1ThreeProngTau_pt_eta2p1"  , "p_{t}", 300,  0., 300. );

  l1ThreeProngTauIso_pt          = fs->make<TH1F>( "l1ThreeProngTauIso_pt"  , "p_{t}", 300,  0., 300. );
  l1ThreeProngTauIso_pt_eta2p4   = fs->make<TH1F>( "l1ThreeProngTauIso_pt_eta2p4"  , "p_{t}", 300,  0., 300. );
  l1ThreeProngTauIso_pt_eta2p1   = fs->make<TH1F>( "l1ThreeProngTauIso_pt_eta2p1"  , "p_{t}", 300,  0., 300. );

  l1EcalCrystals   = fs->make<TH2F>( "ecal_crystals"  , "p_{t}", 100, -4., 4, 100, -4., 4. );
  l1EcalCrystals_2   = fs->make<TH2F>( "ecal_crystals_2"  , "p_{t}", 100, -4., 4, 100, -4., 4. );
  l1EcalCrystals_5   = fs->make<TH2F>( "ecal_crystals_5"  , "p_{t}", 100, -4., 4, 100, -4., 4. );
  l1EcalCrystals_10  = fs->make<TH2F>( "ecal_crystals_10"  , "p_{t}", 100, -4., 4, 100, -4., 4. );
  l1EcalCrystals_15  = fs->make<TH2F>( "ecal_crystals_15"  , "p_{t}", 100, -4., 4, 100, -4., 4. );
  l1EcalCrystals_20  = fs->make<TH2F>( "ecal_crystals_20"  , "p_{t}", 100, -4., 4, 100, -4., 4. );

}


phase2L1TauAnalyzerRates::~phase2L1TauAnalyzerRates()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)


}


//
// member functions
//

// ------------ method called for each event  ------------
void
phase2L1TauAnalyzerRates::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   
   nEvents->Fill(1);

   run   = iEvent.id().run();
   lumi  = iEvent.id().luminosityBlock();
   event = iEvent.id().event();
   
   edm::Handle< std::vector<L1CaloCluster> > l1Clusters;
   iEvent.getByToken( L1ClustersToken_, l1Clusters);

  edm::Handle< std::vector<L1PFObject> > l1PFChargedCandidates;
  iEvent.getByToken( L1PFToken_, l1PFChargedCandidates);

   edm::Handle< std::vector<L1PFTau> > l1PFTaus;
  iEvent.getByToken( L1PFTauToken_, l1PFTaus);

  //control plot for ecal crystals
  edm::Handle<EcalEBTrigPrimDigiCollection> ecaltpgCollection;
  iEvent.getByToken( ecalTPGBToken_, ecaltpgCollection);

  edm::ESHandle<CaloGeometry> caloGeometryHandle;
  iSetup.get<CaloGeometryRecord>().get(caloGeometryHandle);
  const CaloGeometry* caloGeometry_ = caloGeometryHandle.product();  
  ebGeometry = caloGeometry_->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);

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
	  l1EcalCrystals->Fill(eta,phi,et);
	  if(et>2)
	    l1EcalCrystals_2->Fill(eta,phi,et);
	  if(et>5)
	    l1EcalCrystals_5->Fill(eta,phi,et);
	  if(et>10)
	    l1EcalCrystals_10->Fill(eta,phi,et);
	  if(et>15)
	    l1EcalCrystals_15->Fill(eta,phi,et);
	  if(et>20)
	    l1EcalCrystals_20->Fill(eta,phi,et);

	}
    }


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

  std::vector<L1PFTau> l1PFTaus_sorted_iso;

  for(unsigned int i = 0; i < l1PFTaus->size(); i++){
    if(l1PFTaus->at(i).chargedIso()<100)
      if(l1PFTaus->at(i).chargedIso()/l1PFTaus->at(i).pt()<0.2)
	l1PFTaus_sorted_iso.push_back(l1PFTaus->at(i));
  }
  
  std::sort(l1PFTaus_sorted_iso.begin(), l1PFTaus_sorted_iso.end(), [](L1PFTau i,L1PFTau j){return(i.p4().pt() > j.p4().pt());});   


  //filling general rate tree
  if(l1PFTaus_sorted.size()>0)
    l1Tau_pt->Fill(l1PFTaus_sorted.at(0).pt());

  if(l1PFTaus_sorted.size()>1)
    l1DiTau_pt->Fill(l1PFTaus_sorted.at(1).pt());
  
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
    
    if(abs(l1PFTaus_sorted.at(i).eta())<2.1 && !f_2p1){
      l1Tau_pt_eta2p1->Fill(l1PFTaus_sorted.at(i).pt());
      f_2p1=true;
    }

    if(abs(l1PFTaus_sorted.at(i).eta())<2.4 && !f_2p4){
      l1Tau_pt_eta2p4->Fill(l1PFTaus_sorted.at(i).pt());
      f_2p4=true;
    }

    if(l1PFTaus_sorted.at(i).tauType()==0){
      
      if(!f0){
	l1SingleProngTau_pt->Fill(l1PFTaus_sorted.at(i).pt());
	f0 = true;
      }

      if(abs(l1PFTaus_sorted.at(i).eta())<2.1 && !f0_2p1){
        l1SingleProngTau_pt_eta2p1->Fill(l1PFTaus_sorted.at(i).pt());
	f0_2p1=true;
      }
      if(abs(l1PFTaus_sorted.at(i).eta())<2.4 && !f0_2p4){
	l1SingleProngTau_pt_eta2p4->Fill(l1PFTaus_sorted.at(i).pt());
	f0_2p4=true;
      }
      
    }//close 1prong

    if(l1PFTaus_sorted.at(i).tauType()==1){
      
      if(!f1){
	l1SingleProngPi0Tau_pt->Fill(l1PFTaus_sorted.at(i).pt());
	f1=true;
      }

      if(abs(l1PFTaus_sorted.at(i).eta())<2.1 && !f1_2p1){
        l1SingleProngPi0Tau_pt_eta2p1->Fill(l1PFTaus_sorted.at(i).pt());
	f1_2p1=true;
      }
      if(abs(l1PFTaus_sorted.at(i).eta())<2.4 && !f1_2p4){
	l1SingleProngPi0Tau_pt_eta2p4->Fill(l1PFTaus_sorted.at(i).pt());
	f1_2p4=true;
      }
      
    }//close 1prongPi0

    if(l1PFTaus_sorted.at(i).tauType()==10){
      if(!f10){
	l1ThreeProngTau_pt->Fill(l1PFTaus_sorted.at(i).pt());
	f10=true;
      }
      if(abs(l1PFTaus_sorted.at(i).eta())<2.1 && !f10_2p1){
        l1ThreeProngTau_pt_eta2p1->Fill(l1PFTaus_sorted.at(i).pt());
	f10_2p1=true;
      }
      if(abs(l1PFTaus_sorted.at(i).eta())<2.4 && !f10_2p4){
	l1ThreeProngTau_pt_eta2p4->Fill(l1PFTaus_sorted.at(i).pt());
	f10_2p4=true;
      }
      
    }//close 3prong


  }


  //filling iso rate tree
  if(l1PFTaus_sorted_iso.size()>0)
    l1TauIso_pt->Fill(l1PFTaus_sorted_iso.at(0).pt());

  if(l1PFTaus_sorted_iso.size()>1)
    l1DiTauIso_pt->Fill(l1PFTaus_sorted_iso.at(1).pt());
  
  f_2p1=false;
  f_2p4=false;

  f0=false;
  f0_2p1=false;
  f0_2p4=false;

  f1=false;
  f1_2p1=false;
  f1_2p4=false;

  f10 = false;
  f10_2p1=false;
  f10_2p4=false;

  for(unsigned int i = 0; i < l1PFTaus_sorted_iso.size(); i++){
    
    if(abs(l1PFTaus_sorted_iso.at(i).eta())<2.1 && !f_2p1){
      l1TauIso_pt_eta2p1->Fill(l1PFTaus_sorted_iso.at(i).pt());
      f_2p1=true;
    }

    if(abs(l1PFTaus_sorted_iso.at(i).eta())<2.4 && !f_2p4){
      l1TauIso_pt_eta2p4->Fill(l1PFTaus_sorted_iso.at(i).pt());
      f_2p4=true;
    }

    if(l1PFTaus_sorted_iso.at(i).tauType()==0){
      
      if(!f0){
	l1SingleProngTauIso_pt->Fill(l1PFTaus_sorted_iso.at(i).pt());
	f0 = true;
      }

      if(abs(l1PFTaus_sorted_iso.at(i).eta())<2.1 && !f0_2p1){
        l1SingleProngTauIso_pt_eta2p1->Fill(l1PFTaus_sorted_iso.at(i).pt());
	f0_2p1=true;
      }
      if(abs(l1PFTaus_sorted_iso.at(i).eta())<2.4 && !f0_2p4){
	l1SingleProngTauIso_pt_eta2p4->Fill(l1PFTaus_sorted_iso.at(i).pt());
	f0_2p4=true;
      }
      
    }//close 1prong

    if(l1PFTaus_sorted_iso.at(i).tauType()==1){
      
      if(!f1){
	l1SingleProngPi0TauIso_pt->Fill(l1PFTaus_sorted_iso.at(i).pt());
	f1=true;
      }

      if(abs(l1PFTaus_sorted_iso.at(i).eta())<2.1 && !f1_2p1){
        l1SingleProngPi0TauIso_pt_eta2p1->Fill(l1PFTaus_sorted_iso.at(i).pt());
	f1_2p1=true;
      }
      if(abs(l1PFTaus_sorted_iso.at(i).eta())<2.4 && !f1_2p4){
	l1SingleProngPi0TauIso_pt_eta2p4->Fill(l1PFTaus_sorted_iso.at(i).pt());
	f1_2p4=true;
      }
      
    }//close 1prongPi0

    if(l1PFTaus_sorted_iso.at(i).tauType()==10){
      if(!f10){
	l1ThreeProngTauIso_pt->Fill(l1PFTaus_sorted_iso.at(i).pt());
	f10=true;
      }
      if(abs(l1PFTaus_sorted_iso.at(i).eta())<2.1 && !f10_2p1){
        l1ThreeProngTauIso_pt_eta2p1->Fill(l1PFTaus_sorted_iso.at(i).pt());
	f10_2p1=true;
      }
      if(abs(l1PFTaus_sorted_iso.at(i).eta())<2.4 && !f10_2p4){
	l1ThreeProngTauIso_pt_eta2p4->Fill(l1PFTaus_sorted_iso.at(i).pt());
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



   for(auto genPiZero: genPiZeros){

     genPi0Pt = genPiZero.pt();
     genPi0Eta = genPiZero.eta();
     genPi0Phi = genPiZero.phi();
     
     for(unsigned int i = 0; i < l1Clusters->size(); i++){
       if( reco::deltaR(l1Clusters->at(i).p4().Eta(), l1Clusters->at(i).p4().Phi(),genPiZero.p4().eta(),genPiZero.p4().phi()) < 0.1){
	 if(l1Clusters->at(i).isPi0()){
	   objectPi0Eta += (objectPi0Eta*objectPi0Pt + l1Clusters->at(i).p4().Eta()*l1Clusters->at(i).p4().Pt() )/(objectPi0Pt + l1Clusters->at(i).p4().Pt());
	   objectPi0Phi += (objectPi0Phi*objectPi0Pt + l1Clusters->at(i).p4().Phi()*l1Clusters->at(i).p4().Pt() )/(objectPi0Pt + l1Clusters->at(i).p4().Pt());
	   objectPi0Pt += l1Clusters->at(i).p4().Pt();
	   //std::cout<<l1Clusters->at(i)<<std::endl;
	 }
       }
     }
     pi0Tree->Fill();
   }

   for(auto genPiPlus: genPiPluss){

     genPiPlusPt = genPiPlus.pt();
     genPiPlusEta = genPiPlus.eta();
     genPiPlusPhi = genPiPlus.phi();
     float maxPt = 0;

     for(unsigned int iCand = 0; iCand < l1PFChargedCandidates->size(); iCand++){
       if( reco::deltaR(l1PFChargedCandidates->at(iCand).p4().Eta(), l1PFChargedCandidates->at(iCand).p4().Phi(),genPiPlus.p4().eta(),genPiPlus.p4().phi()) < 0.02 && l1PFChargedCandidates->at(iCand).p4().Pt() > maxPt){
	 objectPiPlusEta = l1PFChargedCandidates->at(iCand).p4().Eta();
	 objectPiPlusPhi = l1PFChargedCandidates->at(iCand).p4().Phi();
	 objectPiPlusPt  = l1PFChargedCandidates->at(iCand).p4().Pt();
	 maxPt = l1PFChargedCandidates->at(iCand).p4().Pt();
	 //std::cout<<l1PFChargedCandidates->at(iCand)<<std::endl;	   
	 objectPiPlusIsChargedHadron = l1PFChargedCandidates->at(iCand).isChargedHadron();
       }
     }
       /*for(unsigned int i = 0; i < l1Clusters->size(); i++){
       if( reco::deltaR(l1Clusters->at(i).p4().Eta(), l1Clusters->at(i).p4().Phi(),genPiPlus.p4().eta(),genPiPlus.p4().phi()) < 0.15){
	 if(!l1Clusters->at(i).isPi0()){
	   objectPiPlusEta += (objectPiPlusEta*objectPiPlusPt + l1Clusters->at(i).p4().Eta()*l1Clusters->at(i).p4().Pt() )/(objectPiPlusPt + l1Clusters->at(i).p4().Pt());
	   objectPiPlusPhi += (objectPiPlusPhi*objectPiPlusPt + l1Clusters->at(i).p4().Phi()*l1Clusters->at(i).p4().Pt() )/(objectPiPlusPt + l1Clusters->at(i).p4().Pt());
	   objectPiPlusPt += l1Clusters->at(i).p4().Pt();
	   //std::cout<<l1Clusters->at(i)<<std::endl;
	 }
	 }*/
     //std::cout<<"filling piplus tree"<<std::endl;
     piPlusTree->Fill();
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
     std::cout<<"Tau Decay Mode "<<decayMode<<std::endl;
     std::cout<<"tau vis pt: "<<genPt<<" genEta: "<<genEta<<" genPhi: "<<genPhi<<std::endl;

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

     l1TauPt = 0;
     l1TauEta = -99;
     l1TauPhi = -99;
     l1TauDecayMode = -99;
     l1TauRawIso = 100;
     l1TauNeutralIso = 100;
     l1TauChargedIso = 100;

     recoEta        = -99;
     recoPhi        = -99;
     recoPt         = -99;
     recoChargedIso = -99;
     recoNeutralIso = -99;
     recoRawIso     = -99;
     recoDecayMode  = -99;

     for(unsigned int i = 0; i < l1PFTaus->size(); i++){
       if( reco::deltaR(l1PFTaus->at(i).p4().Eta(), 
			l1PFTaus->at(i).p4().Phi(), 
			genEta, genPhi) < 0.5 &&
	   l1PFTaus->at(i).p4().Pt()>l1TauPt){

	 l1TauEta = l1PFTaus->at(i).p4().Eta();
	 l1TauPhi = l1PFTaus->at(i).p4().Phi();
	 l1TauPt = l1PFTaus->at(i).p4().Pt();
	 l1TauChargedIso = l1PFTaus->at(i).chargedIso();
	 l1TauNeutralIso = l1PFTaus->at(i).neutralIso();
	 l1TauRawIso = l1PFTaus->at(i).rawIso();
	 l1TauDecayMode = l1PFTaus->at(i).tauType();
	 std::cout<<"Match found l1Pt: "<<l1TauPt<<" Eta: "<<l1TauEta<<" Phi: "<<l1TauPhi<<std::endl;
       }
     }
     efficiencyTree->Fill();
   }

   std::cout<<"Finished Analyzing the Taus"<<std::endl;
}


// ------------ method called once each job just before starting event loop  ------------
void 
phase2L1TauAnalyzerRates::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
phase2L1TauAnalyzerRates::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
phase2L1TauAnalyzerRates::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(phase2L1TauAnalyzerRates);
