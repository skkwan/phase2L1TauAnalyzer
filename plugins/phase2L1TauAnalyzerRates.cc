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


#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"

#include "L1Trigger/phase2L1TauAnalyzer/plugins/helpers.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/L1Trigger/interface/L1PFTau.h"

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
  TH1F* track_eta;
  TH1F* track_pt_eta2p4;
  TH1F* track_pt_eta2p1;
  TH1F* track_pt_eta1p4;

  TH1F* l1Tau_pt;     
  TH1F* l1Tau_pt_eta2p4;     
  TH1F* l1Tau_pt_eta2p1;
  TH1F* l1Tau_pt_eta1p4;

  TH1F* l1TauIsoVLoose_pt;     
  TH1F* l1TauIsoVLoose_pt_eta2p4;     
  TH1F* l1TauIsoVLoose_pt_eta2p1;
  TH1F* l1TauIsoVLoose_pt_eta1p4;

  TH1F* l1TauIsoLoose_pt;     
  TH1F* l1TauIsoLoose_pt_eta2p4;     
  TH1F* l1TauIsoLoose_pt_eta2p1;
  TH1F* l1TauIsoLoose_pt_eta1p4;

  TH1F* l1TauIsoMedium_pt;     
  TH1F* l1TauIsoMedium_pt_eta2p4;     
  TH1F* l1TauIsoMedium_pt_eta2p1;
  TH1F* l1TauIsoMedium_pt_eta1p4;

  TH1F* l1TauIsoTight_pt;     
  TH1F* l1TauIsoTight_pt_eta2p4;     
  TH1F* l1TauIsoTight_pt_eta2p1;
  TH1F* l1TauIsoTight_pt_eta1p4;

  TH1F* l1TauIsoTightV0_pt;     
  TH1F* l1TauIsoTightV0_pt_eta2p4;     
  TH1F* l1TauIsoTightV0_pt_eta2p1;
  TH1F* l1TauIsoTightV0_pt_eta1p4;
  /// rel iso
  TH1F* l1TauIsoRelVLoose_pt;     
  TH1F* l1TauIsoRelVLoose_pt_eta2p4;     
  TH1F* l1TauIsoRelVLoose_pt_eta2p1;
  TH1F* l1TauIsoRelVLoose_pt_eta1p4;

  TH1F* l1TauIsoRelLoose_pt;     
  TH1F* l1TauIsoRelLoose_pt_eta2p4;     
  TH1F* l1TauIsoRelLoose_pt_eta2p1;
  TH1F* l1TauIsoRelLoose_pt_eta1p4;

  TH1F* l1TauIsoRelMedium_pt;     
  TH1F* l1TauIsoRelMedium_pt_eta2p4;     
  TH1F* l1TauIsoRelMedium_pt_eta2p1;
  TH1F* l1TauIsoRelMedium_pt_eta1p4;

  TH1F* l1TauIsoRelTight_pt;     
  TH1F* l1TauIsoRelTight_pt_eta2p4;     
  TH1F* l1TauIsoRelTight_pt_eta2p1;
  TH1F* l1TauIsoRelTight_pt_eta1p4;

  TH1F* l1TauIsoRelTightV0_pt;     
  TH1F* l1TauIsoRelTightV0_pt_eta2p4;     
  TH1F* l1TauIsoRelTightV0_pt_eta2p1;
  TH1F* l1TauIsoRelTightV0_pt_eta1p4;
  ////

  TH1F* l1Tau_eta;     
  TH1F* l1Tau_eta_eta2p4;     
  TH1F* l1Tau_eta_eta2p1;
  TH1F* l1Tau_eta_eta1p4;

  /// iso
  TH1F* l1TauIsoVLoose_eta;     
  TH1F* l1TauIsoVLoose_eta_eta2p4;     
  TH1F* l1TauIsoVLoose_eta_eta2p1;
  TH1F* l1TauIsoVLoose_eta_eta1p4;

  TH1F* l1TauIsoLoose_eta;     
  TH1F* l1TauIsoLoose_eta_eta2p4;     
  TH1F* l1TauIsoLoose_eta_eta2p1;
  TH1F* l1TauIsoLoose_eta_eta1p4;

  TH1F* l1TauIsoMedium_eta;     
  TH1F* l1TauIsoMedium_eta_eta2p4;     
  TH1F* l1TauIsoMedium_eta_eta2p1;
  TH1F* l1TauIsoMedium_eta_eta1p4;

  TH1F* l1TauIsoTight_eta;     
  TH1F* l1TauIsoTight_eta_eta2p4;     
  TH1F* l1TauIsoTight_eta_eta2p1;
  TH1F* l1TauIsoTight_eta_eta1p4;

  TH1F* l1TauIsoTightV0_eta;     
  TH1F* l1TauIsoTightV0_eta_eta2p4;     
  TH1F* l1TauIsoTightV0_eta_eta2p1;
  TH1F* l1TauIsoTightV0_eta_eta1p4;

  /// rel iso
  TH1F* l1TauIsoRelVLoose_eta;     
  TH1F* l1TauIsoRelVLoose_eta_eta2p4;     
  TH1F* l1TauIsoRelVLoose_eta_eta2p1;
  TH1F* l1TauIsoRelVLoose_eta_eta1p4;

  TH1F* l1TauIsoRelLoose_eta;     
  TH1F* l1TauIsoRelLoose_eta_eta2p4;     
  TH1F* l1TauIsoRelLoose_eta_eta2p1;
  TH1F* l1TauIsoRelLoose_eta_eta1p4;

  TH1F* l1TauIsoRelMedium_eta;     
  TH1F* l1TauIsoRelMedium_eta_eta2p4;     
  TH1F* l1TauIsoRelMedium_eta_eta2p1;
  TH1F* l1TauIsoRelMedium_eta_eta1p4;

  TH1F* l1TauIsoRelTight_eta;     
  TH1F* l1TauIsoRelTight_eta_eta2p4;     
  TH1F* l1TauIsoRelTight_eta_eta2p1;
  TH1F* l1TauIsoRelTight_eta_eta1p4;

  TH1F* l1TauIsoRelTightV0_eta;     
  TH1F* l1TauIsoRelTightV0_eta_eta2p4;     
  TH1F* l1TauIsoRelTightV0_eta_eta2p1;
  TH1F* l1TauIsoRelTightV0_eta_eta1p4;

  TH1F* l1SingleProngTau_pt;     
  TH1F* l1SingleProngTau_pt_eta2p4;     
  TH1F* l1SingleProngTau_pt_eta2p1;

  TH1F* l1SingleProngTauIso_pt;     
  TH1F* l1SingleProngTauIso_pt_eta2p4;     
  TH1F* l1SingleProngTauIso_pt_eta2p1;

  TH1F* l1SingleProngTauIsoTight_pt;     
  TH1F* l1SingleProngTauIsoTight_pt_eta2p4;     
  TH1F* l1SingleProngTauIsoTight_pt_eta2p1;

  TH1F* l1SingleProngPi0Tau_pt;     
  TH1F* l1SingleProngPi0Tau_pt_eta2p4;     
  TH1F* l1SingleProngPi0Tau_pt_eta2p1;

  TH1F* l1SingleProngPi0TauIso_pt;     
  TH1F* l1SingleProngPi0TauIso_pt_eta2p4;     
  TH1F* l1SingleProngPi0TauIso_pt_eta2p1;

  TH1F* l1SingleProngPi0TauIsoTight_pt;     
  TH1F* l1SingleProngPi0TauIsoTight_pt_eta2p4;     
  TH1F* l1SingleProngPi0TauIsoTight_pt_eta2p1;

  TH1F* l1ThreeProngTau_pt;     
  TH1F* l1ThreeProngTau_pt_eta2p4;     
  TH1F* l1ThreeProngTau_pt_eta2p1;

  TH1F* l1ThreeProngTauIso_pt;     
  TH1F* l1ThreeProngTauIso_pt_eta2p4;     
  TH1F* l1ThreeProngTauIso_pt_eta2p1;

  TH1F* l1ThreeProngTauIsoTight_pt;     
  TH1F* l1ThreeProngTauIsoTight_pt_eta2p4;     
  TH1F* l1ThreeProngTauIsoTight_pt_eta2p1;

  //iso
  TH1F* l1DiTau_pt;     
  TH1F* l1DiTau_pt_eta2p4;     
  TH1F* l1DiTau_pt_eta2p1;
  TH1F* l1DiTau_pt_eta1p4;

  TH1F* l1DiTauIsoVLoose_pt;     
  TH1F* l1DiTauIsoVLoose_pt_eta2p4;     
  TH1F* l1DiTauIsoVLoose_pt_eta2p1;
  TH1F* l1DiTauIsoVLoose_pt_eta1p4;

  TH1F* l1DiTauIsoLoose_pt;     
  TH1F* l1DiTauIsoLoose_pt_eta2p4;     
  TH1F* l1DiTauIsoLoose_pt_eta2p1;
  TH1F* l1DiTauIsoLoose_pt_eta1p4;

  TH1F* l1DiTauIsoMedium_pt;     
  TH1F* l1DiTauIsoMedium_pt_eta2p4;     
  TH1F* l1DiTauIsoMedium_pt_eta2p1;
  TH1F* l1DiTauIsoMedium_pt_eta1p4;

  TH1F* l1DiTauIsoTight_pt;     
  TH1F* l1DiTauIsoTight_pt_eta2p4;     
  TH1F* l1DiTauIsoTight_pt_eta2p1;
  TH1F* l1DiTauIsoTight_pt_eta1p4;

  TH1F* l1DiTauIsoTightV0_pt;     
  TH1F* l1DiTauIsoTightV0_pt_eta2p4;     
  TH1F* l1DiTauIsoTightV0_pt_eta2p1;
  TH1F* l1DiTauIsoTightV0_pt_eta1p4;

  //rel iso
  TH1F* l1DiTauIsoRelVLoose_pt;     
  TH1F* l1DiTauIsoRelVLoose_pt_eta2p4;     
  TH1F* l1DiTauIsoRelVLoose_pt_eta2p1;
  TH1F* l1DiTauIsoRelVLoose_pt_eta1p4;

  TH1F* l1DiTauIsoRelLoose_pt;     
  TH1F* l1DiTauIsoRelLoose_pt_eta2p4;     
  TH1F* l1DiTauIsoRelLoose_pt_eta2p1;
  TH1F* l1DiTauIsoRelLoose_pt_eta1p4;

  TH1F* l1DiTauIsoRelMedium_pt;     
  TH1F* l1DiTauIsoRelMedium_pt_eta2p4;     
  TH1F* l1DiTauIsoRelMedium_pt_eta2p1;
  TH1F* l1DiTauIsoRelMedium_pt_eta1p4;

  TH1F* l1DiTauIsoRelTight_pt;     
  TH1F* l1DiTauIsoRelTight_pt_eta2p4;     
  TH1F* l1DiTauIsoRelTight_pt_eta2p1;
  TH1F* l1DiTauIsoRelTight_pt_eta1p4;

  TH1F* l1DiTauIsoRelTightV0_pt;     
  TH1F* l1DiTauIsoRelTightV0_pt_eta2p4;     
  TH1F* l1DiTauIsoRelTightV0_pt_eta2p1;
  TH1F* l1DiTauIsoRelTightV0_pt_eta1p4;


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


  threeProngTree = fs->make<TTree>("threeProngTree", "Crystal cluster individual crystal pt values");
  threeProngTree->Branch("run",    &run,     "run/I");
  threeProngTree->Branch("lumi",   &lumi,    "lumi/I");
  threeProngTree->Branch("event",  &event,   "event/I");
  
  threeProngTree->Branch("genPt",   &gen3ProngPt,   "genPt/D");
  threeProngTree->Branch("genEta",  &gen3ProngEta,   "genEta/D");
  threeProngTree->Branch("genPhi",  &gen3ProngPhi,   "genPhi/D");
  threeProngTree->Branch("genDM",   &gen3ProngDecayMode,   "genDM/I");


  nEvents              = fs->make<TH1F>( "nEvents"  , "nEvents", 2,  0., 2. );
  track_pt             = fs->make<TH1F>( "track_pt"  , "p_{t}", 300,  0., 300. );
  track_eta             = fs->make<TH1F>( "track_eta"  , "#eta", 100,  -3.5, 3.5 );
  track_pt_eta2p4      = fs->make<TH1F>( "track_pt_eta2p4"  , "p_{t}", 300,  0., 300. );
  track_pt_eta2p1      = fs->make<TH1F>( "track_pt_eta2p1"  , "p_{t}", 300,  0., 300. );
  track_pt_eta1p4      = fs->make<TH1F>( "track_pt_eta1p4"  , "p_{t}", 300,  0., 300. );

  l1Tau_pt             = fs->make<TH1F>( "l1Tau_pt"  , "p_{t}", 300,  0., 300. );
  l1Tau_pt_eta2p4      = fs->make<TH1F>( "l1Tau_pt_eta2p4"  , "p_{t}", 300,  0., 300. );
  l1Tau_pt_eta2p1      = fs->make<TH1F>( "l1Tau_pt_eta2p1"  , "p_{t}", 300,  0., 300. );
  l1Tau_pt_eta1p4      = fs->make<TH1F>( "l1Tau_pt_eta1p4"  , "p_{t}", 300,  0., 300. );

  /// Iso Taus
  l1TauIsoVLoose_pt            = fs->make<TH1F>( "l1TauIsoVLoose_pt"  , "p_{t}", 300,  0., 300. );
  l1TauIsoVLoose_pt_eta2p4     = fs->make<TH1F>( "l1TauIsoVLoose_pt_eta2p4"  , "p_{t}", 300,  0., 300. );
  l1TauIsoVLoose_pt_eta2p1     = fs->make<TH1F>( "l1TauIsoVLoose_pt_eta2p1"  , "p_{t}", 300,  0., 300. );
  l1TauIsoVLoose_pt_eta1p4     = fs->make<TH1F>( "l1TauIsoVLoose_pt_eta1p4"  , "p_{t}", 300,  0., 300. );

  l1TauIsoLoose_pt            = fs->make<TH1F>( "l1TauIsoLoose_pt"  , "p_{t}", 300,  0., 300. );
  l1TauIsoLoose_pt_eta2p4     = fs->make<TH1F>( "l1TauIsoLoose_pt_eta2p4"  , "p_{t}", 300,  0., 300. );
  l1TauIsoLoose_pt_eta2p1     = fs->make<TH1F>( "l1TauIsoLoose_pt_eta2p1"  , "p_{t}", 300,  0., 300. );
  l1TauIsoLoose_pt_eta1p4     = fs->make<TH1F>( "l1TauIsoLoose_pt_eta1p4"  , "p_{t}", 300,  0., 300. );

  l1TauIsoMedium_pt            = fs->make<TH1F>( "l1TauIsoMedium_pt"  , "p_{t}", 300,  0., 300. );
  l1TauIsoMedium_pt_eta2p4     = fs->make<TH1F>( "l1TauIsoMedium_pt_eta2p4"  , "p_{t}", 300,  0., 300. );
  l1TauIsoMedium_pt_eta2p1     = fs->make<TH1F>( "l1TauIsoMedium_pt_eta2p1"  , "p_{t}", 300,  0., 300. );
  l1TauIsoMedium_pt_eta1p4     = fs->make<TH1F>( "l1TauIsoMedium_pt_eta1p4"  , "p_{t}", 300,  0., 300. );

  l1TauIsoTight_pt            = fs->make<TH1F>( "l1TauIsoTight_pt"  , "#pt", 300,  0., 300. );
  l1TauIsoTight_pt_eta2p4     = fs->make<TH1F>( "l1TauIsoTight_pt_eta2p4"  , "#pt", 300,  0., 300. );  
  l1TauIsoTight_pt_eta2p1     = fs->make<TH1F>( "l1TauIsoTight_pt_eta2p1"  , "#pt", 300,  0., 300. );  
  l1TauIsoTight_pt_eta1p4     = fs->make<TH1F>( "l1TauIsoTight_pt_eta1p4"  , "#pt", 300,  0., 300. );  

  l1TauIsoTightV0_pt            = fs->make<TH1F>( "l1TauIsoTightV0_pt"  , "#pt", 300,  0., 300. );
  l1TauIsoTightV0_pt_eta2p4     = fs->make<TH1F>( "l1TauIsoTightV0_pt_eta2p4"  , "#pt", 300,  0., 300. );  
  l1TauIsoTightV0_pt_eta2p1     = fs->make<TH1F>( "l1TauIsoTightV0_pt_eta2p1"  , "#pt", 300,  0., 300. );  
  l1TauIsoTightV0_pt_eta1p4     = fs->make<TH1F>( "l1TauIsoTightV0_pt_eta1p4"  , "#pt", 300,  0., 300. );  

  ///Rel Iso Taus
  l1TauIsoRelVLoose_pt            = fs->make<TH1F>( "l1TauIsoRelVLoose_pt"  , "p_{t}", 300,  0., 300. );
  l1TauIsoRelVLoose_pt_eta2p4     = fs->make<TH1F>( "l1TauIsoRelVLoose_pt_eta2p4"  , "p_{t}", 300,  0., 300. );
  l1TauIsoRelVLoose_pt_eta2p1     = fs->make<TH1F>( "l1TauIsoRelVLoose_pt_eta2p1"  , "p_{t}", 300,  0., 300. );
  l1TauIsoRelVLoose_pt_eta1p4     = fs->make<TH1F>( "l1TauIsoRelVLoose_pt_eta1p4"  , "p_{t}", 300,  0., 300. );

  l1TauIsoRelLoose_pt            = fs->make<TH1F>( "l1TauIsoRelLoose_pt"  , "p_{t}", 300,  0., 300. );
  l1TauIsoRelLoose_pt_eta2p4     = fs->make<TH1F>( "l1TauIsoRelLoose_pt_eta2p4"  , "p_{t}", 300,  0., 300. );
  l1TauIsoRelLoose_pt_eta2p1     = fs->make<TH1F>( "l1TauIsoRelLoose_pt_eta2p1"  , "p_{t}", 300,  0., 300. );
  l1TauIsoRelLoose_pt_eta1p4     = fs->make<TH1F>( "l1TauIsoRelLoose_pt_eta1p4"  , "p_{t}", 300,  0., 300. );

  l1TauIsoRelMedium_pt            = fs->make<TH1F>( "l1TauIsoRelMedium_pt"  , "p_{t}", 300,  0., 300. );
  l1TauIsoRelMedium_pt_eta2p4     = fs->make<TH1F>( "l1TauIsoRelMedium_pt_eta2p4"  , "p_{t}", 300,  0., 300. );
  l1TauIsoRelMedium_pt_eta2p1     = fs->make<TH1F>( "l1TauIsoRelMedium_pt_eta2p1"  , "p_{t}", 300,  0., 300. );
  l1TauIsoRelMedium_pt_eta1p4     = fs->make<TH1F>( "l1TauIsoRelMedium_pt_eta1p4"  , "p_{t}", 300,  0., 300. );

  l1TauIsoRelTight_pt            = fs->make<TH1F>( "l1TauIsoRelTight_pt"  , "#pt", 300,  0., 300. );
  l1TauIsoRelTight_pt_eta2p4     = fs->make<TH1F>( "l1TauIsoRelTight_pt_eta2p4"  , "#pt", 300,  0., 300. );  
  l1TauIsoRelTight_pt_eta2p1     = fs->make<TH1F>( "l1TauIsoRelTight_pt_eta2p1"  , "#pt", 300,  0., 300. );  
  l1TauIsoRelTight_pt_eta1p4     = fs->make<TH1F>( "l1TauIsoRelTight_pt_eta1p4"  , "#pt", 300,  0., 300. );  

  l1TauIsoRelTightV0_pt            = fs->make<TH1F>( "l1TauIsoRelTightV0_pt"  , "#pt", 300,  0., 300. );
  l1TauIsoRelTightV0_pt_eta2p4     = fs->make<TH1F>( "l1TauIsoRelTightV0_pt_eta2p4"  , "#pt", 300,  0., 300. );  
  l1TauIsoRelTightV0_pt_eta2p1     = fs->make<TH1F>( "l1TauIsoRelTightV0_pt_eta2p1"  , "#pt", 300,  0., 300. );  
  l1TauIsoRelTightV0_pt_eta1p4     = fs->make<TH1F>( "l1TauIsoRelTightV0_pt_eta1p4"  , "#pt", 300,  0., 300. );  

  l1Tau_eta             = fs->make<TH1F>( "l1Tau_eta"  , "#eta", 100,  -3.5, 3.5 );
  l1Tau_eta_eta2p4      = fs->make<TH1F>( "l1Tau_eta_eta2p4"  , "#eta", 100,  -3.5, 3.5 );
  l1Tau_eta_eta2p1      = fs->make<TH1F>( "l1Tau_eta_eta2p1"  , "#eta", 100,  -3.5, 3.5 );
  l1Tau_eta_eta1p4      = fs->make<TH1F>( "l1Tau_eta_eta1p4"  , "#eta", 100,  -3.5, 3.5 );

  // Iso Taus eta
  l1TauIsoVLoose_eta            = fs->make<TH1F>( "l1TauIsoVLoose_eta"  , "#eta", 100,  -3.5, 3.5 );
  l1TauIsoVLoose_eta_eta2p4     = fs->make<TH1F>( "l1TauIsoVLoose_eta_eta2p4"  , "#eta", 100,  -3.5, 3.5 );
  l1TauIsoVLoose_eta_eta2p1     = fs->make<TH1F>( "l1TauIsoVLoose_eta_eta2p1"  , "#eta", 100,  -3.5, 3.5 );
  l1TauIsoVLoose_eta_eta1p4     = fs->make<TH1F>( "l1TauIsoVLoose_eta_eta1p4"  , "#eta", 100,  -3.5, 3.5 );

  l1TauIsoLoose_eta            = fs->make<TH1F>( "l1TauIsoLoose_eta"  , "#eta", 100,  -3.5, 3.5 );
  l1TauIsoLoose_eta_eta2p4     = fs->make<TH1F>( "l1TauIsoLoose_eta_eta2p4"  , "#eta", 100,  -3.5, 3.5 );
  l1TauIsoLoose_eta_eta2p1     = fs->make<TH1F>( "l1TauIsoLoose_eta_eta2p1"  , "#eta", 100,  -3.5, 3.5 );
  l1TauIsoLoose_eta_eta1p4     = fs->make<TH1F>( "l1TauIsoLoose_eta_eta1p4"  , "#eta", 100,  -3.5, 3.5 );

  l1TauIsoMedium_eta            = fs->make<TH1F>( "l1TauIsoMedium_eta"  , "#eta", 100,  -3.5, 3.5 );
  l1TauIsoMedium_eta_eta2p4     = fs->make<TH1F>( "l1TauIsoMedium_eta_eta2p4"  , "#eta", 100,  -3.5, 3.5 );
  l1TauIsoMedium_eta_eta2p1     = fs->make<TH1F>( "l1TauIsoMedium_eta_eta2p1"  , "#eta", 100,  -3.5, 3.5 );
  l1TauIsoMedium_eta_eta1p4     = fs->make<TH1F>( "l1TauIsoMedium_eta_eta1p4"  , "#eta", 100,  -3.5, 3.5 );

  l1TauIsoTight_eta            = fs->make<TH1F>( "l1TauIsoTight_eta"  , "#eta", 100,  -3.5, 3.5 );
  l1TauIsoTight_eta_eta2p4     = fs->make<TH1F>( "l1TauIsoTight_eta_eta2p4"  , "#eta", 100,  -3.5, 3.5 );
  l1TauIsoTight_eta_eta2p1     = fs->make<TH1F>( "l1TauIsoTight_eta_eta2p1"  , "#eta", 100,  -3.5, 3.5 );
  l1TauIsoTight_eta_eta1p4     = fs->make<TH1F>( "l1TauIsoTight_eta_eta1p4"  , "#eta", 100,  -3.5, 3.5 );

  l1TauIsoTightV0_eta            = fs->make<TH1F>( "l1TauIsoTightV0_eta"  , "#eta", 100,  -3.5, 3.5 );
  l1TauIsoTightV0_eta_eta2p4     = fs->make<TH1F>( "l1TauIsoTightV0_eta_eta2p4"  , "#eta", 100,  -3.5, 3.5 );
  l1TauIsoTightV0_eta_eta2p1     = fs->make<TH1F>( "l1TauIsoTightV0_eta_eta2p1"  , "#eta", 100,  -3.5, 3.5 );
  l1TauIsoTightV0_eta_eta1p4     = fs->make<TH1F>( "l1TauIsoTightV0_eta_eta1p4"  , "#eta", 100,  -3.5, 3.5 );

  // Rel Iso Taus eta
  l1TauIsoRelVLoose_eta            = fs->make<TH1F>( "l1TauIsoRelVLoose_eta"  , "#eta", 100,  -3.5, 3.5 );
  l1TauIsoRelVLoose_eta_eta2p4     = fs->make<TH1F>( "l1TauIsoRelVLoose_eta_eta2p4"  , "#eta", 100,  -3.5, 3.5 );
  l1TauIsoRelVLoose_eta_eta2p1     = fs->make<TH1F>( "l1TauIsoRelVLoose_eta_eta2p1"  , "#eta", 100,  -3.5, 3.5 );
  l1TauIsoRelVLoose_eta_eta1p4     = fs->make<TH1F>( "l1TauIsoRelVLoose_eta_eta1p4"  , "#eta", 100,  -3.5, 3.5 );

  l1TauIsoRelLoose_eta            = fs->make<TH1F>( "l1TauIsoRelLoose_eta"  , "#eta", 100,  -3.5, 3.5 );
  l1TauIsoRelLoose_eta_eta2p4     = fs->make<TH1F>( "l1TauIsoRelLoose_eta_eta2p4"  , "#eta", 100,  -3.5, 3.5 );
  l1TauIsoRelLoose_eta_eta2p1     = fs->make<TH1F>( "l1TauIsoRelLoose_eta_eta2p1"  , "#eta", 100,  -3.5, 3.5 );
  l1TauIsoRelLoose_eta_eta1p4     = fs->make<TH1F>( "l1TauIsoRelLoose_eta_eta1p4"  , "#eta", 100,  -3.5, 3.5 );

  l1TauIsoRelMedium_eta            = fs->make<TH1F>( "l1TauIsoRelMedium_eta"  , "#eta", 100,  -3.5, 3.5 );
  l1TauIsoRelMedium_eta_eta2p4     = fs->make<TH1F>( "l1TauIsoRelMedium_eta_eta2p4"  , "#eta", 100,  -3.5, 3.5 );
  l1TauIsoRelMedium_eta_eta2p1     = fs->make<TH1F>( "l1TauIsoRelMedium_eta_eta2p1"  , "#eta", 100,  -3.5, 3.5 );
  l1TauIsoRelMedium_eta_eta1p4     = fs->make<TH1F>( "l1TauIsoRelMedium_eta_eta1p4"  , "#eta", 100,  -3.5, 3.5 );

  l1TauIsoRelTight_eta            = fs->make<TH1F>( "l1TauIsoRelTight_eta"  , "#eta", 100,  -3.5, 3.5 );
  l1TauIsoRelTight_eta_eta2p4     = fs->make<TH1F>( "l1TauIsoRelTight_eta_eta2p4"  , "#eta", 100,  -3.5, 3.5 );
  l1TauIsoRelTight_eta_eta2p1     = fs->make<TH1F>( "l1TauIsoRelTight_eta_eta2p1"  , "#eta", 100,  -3.5, 3.5 );
  l1TauIsoRelTight_eta_eta1p4     = fs->make<TH1F>( "l1TauIsoRelTight_eta_eta1p4"  , "#eta", 100,  -3.5, 3.5 );

  l1TauIsoRelTightV0_eta            = fs->make<TH1F>( "l1TauIsoRelTightV0_eta"  , "#eta", 100,  -3.5, 3.5 );
  l1TauIsoRelTightV0_eta_eta2p4     = fs->make<TH1F>( "l1TauIsoRelTightV0_eta_eta2p4"  , "#eta", 100,  -3.5, 3.5 );
  l1TauIsoRelTightV0_eta_eta2p1     = fs->make<TH1F>( "l1TauIsoRelTightV0_eta_eta2p1"  , "#eta", 100,  -3.5, 3.5 );
  l1TauIsoRelTightV0_eta_eta1p4     = fs->make<TH1F>( "l1TauIsoRelTightV0_eta_eta1p4"  , "#eta", 100,  -3.5, 3.5 );

  l1DiTau_pt             = fs->make<TH1F>( "l1DiTau_pt"  , "p_{t}", 300,  0., 300. );
  l1DiTau_pt_eta2p4      = fs->make<TH1F>( "l1DiTau_pt_eta2p4"  , "p_{t}", 300,  0., 300. );
  l1DiTau_pt_eta2p1      = fs->make<TH1F>( "l1DiTau_pt_eta2p1"  , "p_{t}", 300,  0., 300. );
  l1DiTau_pt_eta1p4      = fs->make<TH1F>( "l1DiTau_pt_eta1p4"  , "p_{t}", 300,  0., 300. );

  // Di Iso taus
  l1DiTauIsoVLoose_pt          = fs->make<TH1F>( "l1DiTauIsoVLoose_pt"  , "p_{t}", 300,  0., 300. );
  l1DiTauIsoVLoose_pt_eta2p4   = fs->make<TH1F>( "l1DiTauIsoVLoose_pt_eta2p4"  , "p_{t}", 300,  0., 300. );
  l1DiTauIsoVLoose_pt_eta2p1   = fs->make<TH1F>( "l1DiTauIsoVLoose_pt_eta2p1"  , "p_{t}", 300,  0., 300. );
  l1DiTauIsoVLoose_pt_eta1p4   = fs->make<TH1F>( "l1DiTauIsoVLoose_pt_eta1p4"  , "p_{t}", 300,  0., 300. );

  l1DiTauIsoLoose_pt          = fs->make<TH1F>( "l1DiTauIsoLoose_pt"  , "p_{t}", 300,  0., 300. );
  l1DiTauIsoLoose_pt_eta2p4   = fs->make<TH1F>( "l1DiTauIsoLoose_pt_eta2p4"  , "p_{t}", 300,  0., 300. );
  l1DiTauIsoLoose_pt_eta2p1   = fs->make<TH1F>( "l1DiTauIsoLoose_pt_eta2p1"  , "p_{t}", 300,  0., 300. );
  l1DiTauIsoLoose_pt_eta1p4   = fs->make<TH1F>( "l1DiTauIsoLoose_pt_eta1p4"  , "p_{t}", 300,  0., 300. );

  l1DiTauIsoMedium_pt          = fs->make<TH1F>( "l1DiTauIsoMedium_pt"  , "p_{t}", 300,  0., 300. );
  l1DiTauIsoMedium_pt_eta2p4   = fs->make<TH1F>( "l1DiTauIsoMedium_pt_eta2p4"  , "p_{t}", 300,  0., 300. );
  l1DiTauIsoMedium_pt_eta2p1   = fs->make<TH1F>( "l1DiTauIsoMedium_pt_eta2p1"  , "p_{t}", 300,  0., 300. );
  l1DiTauIsoMedium_pt_eta1p4   = fs->make<TH1F>( "l1DiTauIsoMedium_pt_eta1p4"  , "p_{t}", 300,  0., 300. );

  l1DiTauIsoTight_pt          = fs->make<TH1F>( "l1DiTauIsoTight_pt"  , "p_{t}", 300,  0., 300. );
  l1DiTauIsoTight_pt_eta2p4   = fs->make<TH1F>( "l1DiTauIsoTight_pt_eta2p4"  , "p_{t}", 300,  0., 300. );
  l1DiTauIsoTight_pt_eta2p1   = fs->make<TH1F>( "l1DiTauIsoTight_pt_eta2p1"  , "p_{t}", 300,  0., 300. );
  l1DiTauIsoTight_pt_eta1p4   = fs->make<TH1F>( "l1DiTauIsoTight_pt_eta1p4"  , "p_{t}", 300,  0., 300. );

  l1DiTauIsoTightV0_pt          = fs->make<TH1F>( "l1DiTauIsoTightV0_pt"  , "p_{t}", 300,  0., 300. );
  l1DiTauIsoTightV0_pt_eta2p4   = fs->make<TH1F>( "l1DiTauIsoTightV0_pt_eta2p4"  , "p_{t}", 300,  0., 300. );
  l1DiTauIsoTightV0_pt_eta2p1   = fs->make<TH1F>( "l1DiTauIsoTightV0_pt_eta2p1"  , "p_{t}", 300,  0., 300. );
  l1DiTauIsoTightV0_pt_eta1p4   = fs->make<TH1F>( "l1DiTauIsoTightV0_pt_eta1p4"  , "p_{t}", 300,  0., 300. );

  // Di Rel Iso taus
  l1DiTauIsoRelVLoose_pt          = fs->make<TH1F>( "l1DiTauIsoRelVLoose_pt"  , "p_{t}", 300,  0., 300. );
  l1DiTauIsoRelVLoose_pt_eta2p4   = fs->make<TH1F>( "l1DiTauIsoRelVLoose_pt_eta2p4"  , "p_{t}", 300,  0., 300. );
  l1DiTauIsoRelVLoose_pt_eta2p1   = fs->make<TH1F>( "l1DiTauIsoRelVLoose_pt_eta2p1"  , "p_{t}", 300,  0., 300. );
  l1DiTauIsoRelVLoose_pt_eta1p4   = fs->make<TH1F>( "l1DiTauIsoRelVLoose_pt_eta1p4"  , "p_{t}", 300,  0., 300. );

  l1DiTauIsoRelLoose_pt          = fs->make<TH1F>( "l1DiTauIsoRelLoose_pt"  , "p_{t}", 300,  0., 300. );
  l1DiTauIsoRelLoose_pt_eta2p4   = fs->make<TH1F>( "l1DiTauIsoRelLoose_pt_eta2p4"  , "p_{t}", 300,  0., 300. );
  l1DiTauIsoRelLoose_pt_eta2p1   = fs->make<TH1F>( "l1DiTauIsoRelLoose_pt_eta2p1"  , "p_{t}", 300,  0., 300. );
  l1DiTauIsoRelLoose_pt_eta1p4   = fs->make<TH1F>( "l1DiTauIsoRelLoose_pt_eta1p4"  , "p_{t}", 300,  0., 300. );

  l1DiTauIsoRelMedium_pt          = fs->make<TH1F>( "l1DiTauIsoRelMedium_pt"  , "p_{t}", 300,  0., 300. );
  l1DiTauIsoRelMedium_pt_eta2p4   = fs->make<TH1F>( "l1DiTauIsoRelMedium_pt_eta2p4"  , "p_{t}", 300,  0., 300. );
  l1DiTauIsoRelMedium_pt_eta2p1   = fs->make<TH1F>( "l1DiTauIsoRelMedium_pt_eta2p1"  , "p_{t}", 300,  0., 300. );
  l1DiTauIsoRelMedium_pt_eta1p4   = fs->make<TH1F>( "l1DiTauIsoRelMedium_pt_eta1p4"  , "p_{t}", 300,  0., 300. );

  l1DiTauIsoRelTight_pt          = fs->make<TH1F>( "l1DiTauIsoRelTight_pt"  , "p_{t}", 300,  0., 300. );
  l1DiTauIsoRelTight_pt_eta2p4   = fs->make<TH1F>( "l1DiTauIsoRelTight_pt_eta2p4"  , "p_{t}", 300,  0., 300. );
  l1DiTauIsoRelTight_pt_eta2p1   = fs->make<TH1F>( "l1DiTauIsoRelTight_pt_eta2p1"  , "p_{t}", 300,  0., 300. );
  l1DiTauIsoRelTight_pt_eta1p4   = fs->make<TH1F>( "l1DiTauIsoRelTight_pt_eta1p4"  , "p_{t}", 300,  0., 300. );

  l1DiTauIsoRelTightV0_pt          = fs->make<TH1F>( "l1DiTauIsoRelTightV0_pt"  , "p_{t}", 300,  0., 300. );
  l1DiTauIsoRelTightV0_pt_eta2p4   = fs->make<TH1F>( "l1DiTauIsoRelTightV0_pt_eta2p4"  , "p_{t}", 300,  0., 300. );
  l1DiTauIsoRelTightV0_pt_eta2p1   = fs->make<TH1F>( "l1DiTauIsoRelTightV0_pt_eta2p1"  , "p_{t}", 300,  0., 300. );
  l1DiTauIsoRelTightV0_pt_eta1p4   = fs->make<TH1F>( "l1DiTauIsoRelTightV0_pt_eta1p4"  , "p_{t}", 300,  0., 300. );

  l1SingleProngTau_pt            = fs->make<TH1F>( "l1SingleProngTau"  , "p_{t}", 300,  0., 300. );
  l1SingleProngTau_pt_eta2p4     = fs->make<TH1F>( "l1SingleProngTau_eta2p4"  , "p_{t}", 300,  0., 300. );
  l1SingleProngTau_pt_eta2p1     = fs->make<TH1F>( "l1SingleProngTau_eta2p1"  , "p_{t}", 300,  0., 300. );

  l1SingleProngTauIso_pt         = fs->make<TH1F>( "l1SingleProngTauIso"  , "p_{t}", 300,  0., 300. );
  l1SingleProngTauIso_pt_eta2p4  = fs->make<TH1F>( "l1SingleProngTauIso_eta2p4"  , "p_{t}", 300,  0., 300. );
  l1SingleProngTauIso_pt_eta2p1  = fs->make<TH1F>( "l1SingleProngTauIso_eta2p1"  , "p_{t}", 300,  0., 300. );

  l1SingleProngTauIsoTight_pt         = fs->make<TH1F>( "l1SingleProngTauIsoTight"  , "p_{t}", 300,  0., 300. );
  l1SingleProngTauIsoTight_pt_eta2p4  = fs->make<TH1F>( "l1SingleProngTauIsoTight_eta2p4"  , "p_{t}", 300,  0., 300. );
  l1SingleProngTauIsoTight_pt_eta2p1  = fs->make<TH1F>( "l1SingleProngTauIsoTight_eta2p1"  , "p_{t}", 300,  0., 300. );

  l1SingleProngPi0Tau_pt            = fs->make<TH1F>( "l1SingleProngPi0Tau"  , "p_{t}", 300,  0., 300. );
  l1SingleProngPi0Tau_pt_eta2p4     = fs->make<TH1F>( "l1SingleProngPi0Tau_eta2p4"  , "p_{t}", 300,  0., 300. );
  l1SingleProngPi0Tau_pt_eta2p1     = fs->make<TH1F>( "l1SingleProngPi0Tau_eta2p1"  , "p_{t}", 300,  0., 300. );

  l1SingleProngPi0TauIso_pt         = fs->make<TH1F>( "l1SingleProngPi0TauIso"  , "p_{t}", 300,  0., 300. );
  l1SingleProngPi0TauIso_pt_eta2p4  = fs->make<TH1F>( "l1SingleProngPi0TauIso_eta2p4"  , "p_{t}", 300,  0., 300. );
  l1SingleProngPi0TauIso_pt_eta2p1  = fs->make<TH1F>( "l1SingleProngPi0TauIso_eta2p1"  , "p_{t}", 300,  0., 300. );

  l1SingleProngPi0TauIsoTight_pt         = fs->make<TH1F>( "l1SingleProngPi0TauIsoTight"  , "p_{t}", 300,  0., 300. );
  l1SingleProngPi0TauIsoTight_pt_eta2p4  = fs->make<TH1F>( "l1SingleProngPi0TauIsoTight_eta2p4"  , "p_{t}", 300,  0., 300. );
  l1SingleProngPi0TauIsoTight_pt_eta2p1  = fs->make<TH1F>( "l1SingleProngPi0TauIsoTight_eta2p1"  , "p_{t}", 300,  0., 300. );

  l1ThreeProngTau_pt             = fs->make<TH1F>( "l1ThreeProngTau_pt"  , "p_{t}", 300,  0., 300. );
  l1ThreeProngTau_pt_eta2p4      = fs->make<TH1F>( "l1ThreeProngTau_pt_eta2p4"  , "p_{t}", 300,  0., 300. );
  l1ThreeProngTau_pt_eta2p1      = fs->make<TH1F>( "l1ThreeProngTau_pt_eta2p1"  , "p_{t}", 300,  0., 300. );

  l1ThreeProngTauIso_pt          = fs->make<TH1F>( "l1ThreeProngTauIso_pt"  , "p_{t}", 300,  0., 300. );
  l1ThreeProngTauIso_pt_eta2p4   = fs->make<TH1F>( "l1ThreeProngTauIso_pt_eta2p4"  , "p_{t}", 300,  0., 300. );
  l1ThreeProngTauIso_pt_eta2p1   = fs->make<TH1F>( "l1ThreeProngTauIso_pt_eta2p1"  , "p_{t}", 300,  0., 300. );

  l1ThreeProngTauIsoTight_pt          = fs->make<TH1F>( "l1ThreeProngTauIsoTight_pt"  , "p_{t}", 300,  0., 300. );
  l1ThreeProngTauIsoTight_pt_eta2p4   = fs->make<TH1F>( "l1ThreeProngTauIsoTight_pt_eta2p4"  , "p_{t}", 300,  0., 300. );
  l1ThreeProngTauIsoTight_pt_eta2p1   = fs->make<TH1F>( "l1ThreeProngTauIsoTight_pt_eta2p1"  , "p_{t}", 300,  0., 300. );

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
   
   edm::Handle< std::vector<L1PFTau> > l1PFTaus;
   iEvent.getByToken( L1PFTauToken_, l1PFTaus);


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
  
  if(l1Tracks.size()>0){
    track_pt->Fill(l1Tracks.at(0).getMomentum().perp());
    track_eta->Fill(l1Tracks.at(0).getMomentum().eta());
  }
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

  for(auto l1Track: l1Tracks){
    if(fabs(l1Track.getMomentum().eta())<1.4){
      track_pt_eta1p4->Fill(l1Track.getMomentum().perp());
      break;
    }
  }
  

  std::vector<L1PFTau> l1PFTaus_sorted;
  std::vector<L1PFTau> l1PFTaus_sorted_eta2p4;
  std::vector<L1PFTau> l1PFTaus_sorted_eta2p1;
  std::vector<L1PFTau> l1PFTaus_sorted_eta1p4;

  l1PFTaus_sorted.clear(); l1PFTaus_sorted_eta2p4.clear(); l1PFTaus_sorted_eta2p1.clear();  l1PFTaus_sorted_eta1p4.clear();

  for(unsigned int i = 0; i < l1PFTaus->size(); i++){
    l1PFTaus_sorted.push_back(l1PFTaus->at(i));

    if(fabs(l1PFTaus->at(i).eta())<2.4){
      l1PFTaus_sorted_eta2p4.push_back(l1PFTaus->at(i));}

    if(fabs(l1PFTaus->at(i).eta())<2.1){
      l1PFTaus_sorted_eta2p1.push_back(l1PFTaus->at(i));}

    if(fabs(l1PFTaus->at(i).eta())<1.4){
      l1PFTaus_sorted_eta1p4.push_back(l1PFTaus->at(i));}
  }

  std::sort(l1PFTaus_sorted.begin(),        l1PFTaus_sorted.end(),        [](L1PFTau i,L1PFTau j){return(i.p4().pt() > j.p4().pt());});   
  std::sort(l1PFTaus_sorted_eta2p4.begin(), l1PFTaus_sorted_eta2p4.end(), [](L1PFTau i,L1PFTau j){return(i.p4().pt() > j.p4().pt());});   
  std::sort(l1PFTaus_sorted_eta2p1.begin(), l1PFTaus_sorted_eta2p1.end(), [](L1PFTau i,L1PFTau j){return(i.p4().pt() > j.p4().pt());});   
  std::sort(l1PFTaus_sorted_eta1p4.begin(), l1PFTaus_sorted_eta1p4.end(), [](L1PFTau i,L1PFTau j){return(i.p4().pt() > j.p4().pt());});   

  std::vector<L1PFTau> l1PFTaus_sorted_iso_VLoose;
  std::vector<L1PFTau> l1PFTaus_sorted_iso_VLoose_eta2p4;
  std::vector<L1PFTau> l1PFTaus_sorted_iso_VLoose_eta2p1;
  std::vector<L1PFTau> l1PFTaus_sorted_iso_VLoose_eta1p4;

  std::vector<L1PFTau> l1PFTaus_sorted_iso_Loose;
  std::vector<L1PFTau> l1PFTaus_sorted_iso_Loose_eta2p4;
  std::vector<L1PFTau> l1PFTaus_sorted_iso_Loose_eta2p1;
  std::vector<L1PFTau> l1PFTaus_sorted_iso_Loose_eta1p4;

  std::vector<L1PFTau> l1PFTaus_sorted_iso_Medium;
  std::vector<L1PFTau> l1PFTaus_sorted_iso_Medium_eta2p4;
  std::vector<L1PFTau> l1PFTaus_sorted_iso_Medium_eta2p1;
  std::vector<L1PFTau> l1PFTaus_sorted_iso_Medium_eta1p4;

  std::vector<L1PFTau> l1PFTaus_sorted_iso_Tight;
  std::vector<L1PFTau> l1PFTaus_sorted_iso_Tight_eta2p4;
  std::vector<L1PFTau> l1PFTaus_sorted_iso_Tight_eta2p1;
  std::vector<L1PFTau> l1PFTaus_sorted_iso_Tight_eta1p4;

  // rel iso taus
  std::vector<L1PFTau> l1PFTaus_sorted_iso_rel_VLoose;
  std::vector<L1PFTau> l1PFTaus_sorted_iso_rel_VLoose_eta2p4;
  std::vector<L1PFTau> l1PFTaus_sorted_iso_rel_VLoose_eta2p1;
  std::vector<L1PFTau> l1PFTaus_sorted_iso_rel_VLoose_eta1p4;

  std::vector<L1PFTau> l1PFTaus_sorted_iso_rel_Loose;
  std::vector<L1PFTau> l1PFTaus_sorted_iso_rel_Loose_eta2p4;
  std::vector<L1PFTau> l1PFTaus_sorted_iso_rel_Loose_eta2p1;
  std::vector<L1PFTau> l1PFTaus_sorted_iso_rel_Loose_eta1p4;

  std::vector<L1PFTau> l1PFTaus_sorted_iso_rel_Medium;
  std::vector<L1PFTau> l1PFTaus_sorted_iso_rel_Medium_eta2p4;
  std::vector<L1PFTau> l1PFTaus_sorted_iso_rel_Medium_eta2p1;
  std::vector<L1PFTau> l1PFTaus_sorted_iso_rel_Medium_eta1p4;

  std::vector<L1PFTau> l1PFTaus_sorted_iso_rel_Tight;
  std::vector<L1PFTau> l1PFTaus_sorted_iso_rel_Tight_eta2p4;
  std::vector<L1PFTau> l1PFTaus_sorted_iso_rel_Tight_eta2p1;
  std::vector<L1PFTau> l1PFTaus_sorted_iso_rel_Tight_eta1p4;
  
  // here
  for(unsigned int i = 0; i < l1PFTaus->size(); i++){
    if(l1PFTaus->at(i).passVLooseIso()){
      l1PFTaus_sorted_iso_VLoose.push_back(l1PFTaus->at(i));
      
      if(fabs(l1PFTaus->at(i).eta())<2.4){
	l1PFTaus_sorted_iso_VLoose_eta2p4.push_back(l1PFTaus->at(i));
      }
      if(fabs(l1PFTaus->at(i).eta())<2.1){
	l1PFTaus_sorted_iso_VLoose_eta2p1.push_back(l1PFTaus->at(i));
      }
      if(fabs(l1PFTaus->at(i).eta())<1.4){
	l1PFTaus_sorted_iso_VLoose_eta1p4.push_back(l1PFTaus->at(i));
      }
    }
    //loose
    if(l1PFTaus->at(i).passLooseIso()){
      l1PFTaus_sorted_iso_Loose.push_back(l1PFTaus->at(i));
      
      if(fabs(l1PFTaus->at(i).eta())<2.4){
	l1PFTaus_sorted_iso_Loose_eta2p4.push_back(l1PFTaus->at(i));
      }
      if(fabs(l1PFTaus->at(i).eta())<2.1){
	l1PFTaus_sorted_iso_Loose_eta2p1.push_back(l1PFTaus->at(i));
      }
      if(fabs(l1PFTaus->at(i).eta())<1.4){
	l1PFTaus_sorted_iso_Loose_eta1p4.push_back(l1PFTaus->at(i));
      }
    }
    //med
    if(l1PFTaus->at(i).passMediumIso()){
      l1PFTaus_sorted_iso_Medium.push_back(l1PFTaus->at(i));
      
      if(fabs(l1PFTaus->at(i).eta())<2.4){
	l1PFTaus_sorted_iso_Medium_eta2p4.push_back(l1PFTaus->at(i));
      }
      if(fabs(l1PFTaus->at(i).eta())<2.1){
	l1PFTaus_sorted_iso_Medium_eta2p1.push_back(l1PFTaus->at(i));
      }
      if(fabs(l1PFTaus->at(i).eta())<1.4){
	l1PFTaus_sorted_iso_Medium_eta1p4.push_back(l1PFTaus->at(i));
      }
    }
    //tight
    if(l1PFTaus->at(i).passTightIso()){
      l1PFTaus_sorted_iso_Tight.push_back(l1PFTaus->at(i));
      
      if(fabs(l1PFTaus->at(i).eta())<2.4){
	l1PFTaus_sorted_iso_Tight_eta2p4.push_back(l1PFTaus->at(i));
      }
      if(fabs(l1PFTaus->at(i).eta())<2.1){
	l1PFTaus_sorted_iso_Tight_eta2p1.push_back(l1PFTaus->at(i));
      }
      if(fabs(l1PFTaus->at(i).eta())<1.4){
	l1PFTaus_sorted_iso_Tight_eta1p4.push_back(l1PFTaus->at(i));
      }
    }

    //rel iso
    if(l1PFTaus->at(i).passVLooseRelIso()){
      l1PFTaus_sorted_iso_rel_VLoose.push_back(l1PFTaus->at(i));
      
      if(fabs(l1PFTaus->at(i).eta())<2.4){
	l1PFTaus_sorted_iso_rel_VLoose_eta2p4.push_back(l1PFTaus->at(i));
      }
      if(fabs(l1PFTaus->at(i).eta())<2.1){
	l1PFTaus_sorted_iso_rel_VLoose_eta2p1.push_back(l1PFTaus->at(i));
      }
      if(fabs(l1PFTaus->at(i).eta())<1.4){
	l1PFTaus_sorted_iso_rel_VLoose_eta1p4.push_back(l1PFTaus->at(i));
      }
    }
    //loose
    if(l1PFTaus->at(i).passLooseRelIso()){
      l1PFTaus_sorted_iso_rel_Loose.push_back(l1PFTaus->at(i));
      
      if(fabs(l1PFTaus->at(i).eta())<2.4){
	l1PFTaus_sorted_iso_rel_Loose_eta2p4.push_back(l1PFTaus->at(i));
      }
      if(fabs(l1PFTaus->at(i).eta())<2.1){
	l1PFTaus_sorted_iso_rel_Loose_eta2p1.push_back(l1PFTaus->at(i));
      }
      if(fabs(l1PFTaus->at(i).eta())<1.4){
	l1PFTaus_sorted_iso_rel_Loose_eta1p4.push_back(l1PFTaus->at(i));
      }
    }
    //med
    if(l1PFTaus->at(i).passMediumRelIso()){
      l1PFTaus_sorted_iso_rel_Medium.push_back(l1PFTaus->at(i));
      
      if(fabs(l1PFTaus->at(i).eta())<2.4){
	l1PFTaus_sorted_iso_rel_Medium_eta2p4.push_back(l1PFTaus->at(i));
      }
      if(fabs(l1PFTaus->at(i).eta())<2.1){
	l1PFTaus_sorted_iso_rel_Medium_eta2p1.push_back(l1PFTaus->at(i));
      }
      if(fabs(l1PFTaus->at(i).eta())<1.4){
	l1PFTaus_sorted_iso_rel_Medium_eta1p4.push_back(l1PFTaus->at(i));
      }
    }
    //tight
    if(l1PFTaus->at(i).passTightRelIso()){
      l1PFTaus_sorted_iso_rel_Tight.push_back(l1PFTaus->at(i));
      
      if(fabs(l1PFTaus->at(i).eta())<2.4){
	l1PFTaus_sorted_iso_rel_Tight_eta2p4.push_back(l1PFTaus->at(i));
      }
      if(fabs(l1PFTaus->at(i).eta())<2.1){
	l1PFTaus_sorted_iso_rel_Tight_eta2p1.push_back(l1PFTaus->at(i));
      }
      if(fabs(l1PFTaus->at(i).eta())<1.4){
	l1PFTaus_sorted_iso_rel_Tight_eta1p4.push_back(l1PFTaus->at(i));
      }
    }

  }
  
  ///check if sorting is really needed... taus arrive pre-sorted  
  //std::sort(l1PFTaus_sorted_iso.begin(),        l1PFTaus_sorted_iso.end(),        [](L1PFTau i,L1PFTau j){return(i.p4().pt() > j.p4().pt());});   
  //std::sort(l1PFTaus_sorted_iso_eta2p4.begin(), l1PFTaus_sorted_iso_eta2p4.end(), [](L1PFTau i,L1PFTau j){return(i.p4().pt() > j.p4().pt());});   
  //std::sort(l1PFTaus_sorted_iso_eta2p1.begin(), l1PFTaus_sorted_iso_eta2p1.end(), [](L1PFTau i,L1PFTau j){return(i.p4().pt() > j.p4().pt());});   
  //std::sort(l1PFTaus_sorted_iso_eta1p4.begin(), l1PFTaus_sorted_iso_eta1p4.end(), [](L1PFTau i,L1PFTau j){return(i.p4().pt() > j.p4().pt());});   


  // 
  std::vector<L1PFTau> l1PFTaus_sorted_iso_2;
  std::vector<L1PFTau> l1PFTaus_sorted_iso_2_eta2p4;
  std::vector<L1PFTau> l1PFTaus_sorted_iso_2_eta2p1;
  std::vector<L1PFTau> l1PFTaus_sorted_iso_2_eta1p4;

  for(unsigned int i = 0; i < l1PFTaus->size(); i++){
    if(l1PFTaus->at(i).chargedIso()/l1PFTaus->at(i).p4().pt()<0.1){
      l1PFTaus_sorted_iso_2.push_back(l1PFTaus->at(i));
      
      if(fabs(l1PFTaus->at(i).eta())<2.4){
	l1PFTaus_sorted_iso_2_eta2p4.push_back(l1PFTaus->at(i));
      }
      if(fabs(l1PFTaus->at(i).eta())<2.1){
	l1PFTaus_sorted_iso_2_eta2p1.push_back(l1PFTaus->at(i));
      }
      if(fabs(l1PFTaus->at(i).eta())<1.4){
	l1PFTaus_sorted_iso_2_eta1p4.push_back(l1PFTaus->at(i));
      }
    }
  }
  
  std::sort(l1PFTaus_sorted_iso_2.begin(),        l1PFTaus_sorted_iso_2.end(),        [](L1PFTau i,L1PFTau j){return(i.p4().pt() > j.p4().pt());});   
  std::sort(l1PFTaus_sorted_iso_2_eta2p4.begin(), l1PFTaus_sorted_iso_2_eta2p4.end(), [](L1PFTau i,L1PFTau j){return(i.p4().pt() > j.p4().pt());});   
  std::sort(l1PFTaus_sorted_iso_2_eta2p1.begin(), l1PFTaus_sorted_iso_2_eta2p1.end(), [](L1PFTau i,L1PFTau j){return(i.p4().pt() > j.p4().pt());});   
  std::sort(l1PFTaus_sorted_iso_2_eta1p4.begin(), l1PFTaus_sorted_iso_2_eta1p4.end(), [](L1PFTau i,L1PFTau j){return(i.p4().pt() > j.p4().pt());});   


  //filling general rate tree
  if(l1PFTaus_sorted.size()>0){
    l1Tau_pt->Fill(l1PFTaus_sorted.at(0).pt());
    l1Tau_eta->Fill(l1PFTaus_sorted.at(0).eta());
  }

  if(l1PFTaus_sorted.size()>1)
    l1DiTau_pt->Fill(l1PFTaus_sorted.at(1).pt());

  //eta restricted 2.4
  if(l1PFTaus_sorted_eta2p4.size()>0){
    l1Tau_pt_eta2p4->Fill(  l1PFTaus_sorted_eta2p4.at(0).pt());
    l1Tau_eta_eta2p4->Fill(  l1PFTaus_sorted_eta2p4.at(0).eta());
  }

  if(l1PFTaus_sorted_eta2p4.size()>1)
    l1DiTau_pt_eta2p4->Fill(l1PFTaus_sorted_eta2p4.at(1).pt());

  //eta restricted 2.1
  if(l1PFTaus_sorted_eta2p1.size()>0){
    l1Tau_pt_eta2p1->Fill(  l1PFTaus_sorted_eta2p1.at(0).pt());
    l1Tau_eta_eta2p1->Fill(  l1PFTaus_sorted_eta2p1.at(0).eta());
  }

  if(l1PFTaus_sorted_eta2p1.size()>1)
    l1DiTau_pt_eta2p1->Fill(l1PFTaus_sorted_eta2p1.at(1).pt());

  //eta restricted 1.4
  if(l1PFTaus_sorted_eta1p4.size()>0){
    l1Tau_pt_eta1p4->Fill(  l1PFTaus_sorted_eta1p4.at(0).pt());
    l1Tau_eta_eta1p4->Fill(  l1PFTaus_sorted_eta1p4.at(0).eta());
  }

  if(l1PFTaus_sorted_eta1p4.size()>1)
    l1DiTau_pt_eta1p4->Fill(l1PFTaus_sorted_eta1p4.at(1).pt());
  
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
    /*  
    if(abs(l1PFTaus_sorted.at(i).eta())<2.1 && !f_2p1){
      l1Tau_pt_eta2p1->Fill(l1PFTaus_sorted.at(i).pt());
      f_2p1=true;
    }

    if(abs(l1PFTaus_sorted.at(i).eta())<2.4 && !f_2p4){
      l1Tau_pt_eta2p4->Fill(l1PFTaus_sorted.at(i).pt());
      f_2p4=true;
     }
    */
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
  // VLoose
  if(l1PFTaus_sorted_iso_VLoose.size()>0){
    l1TauIsoVLoose_pt->Fill(l1PFTaus_sorted_iso_VLoose.at(0).pt());
    l1TauIsoVLoose_eta->Fill(l1PFTaus_sorted_iso_VLoose.at(0).eta());
  }

  if(l1PFTaus_sorted_iso_VLoose.size()>1)
    l1DiTauIsoVLoose_pt->Fill(l1PFTaus_sorted_iso_VLoose.at(1).pt());
  
  //eta restricted 2.4
  if(l1PFTaus_sorted_iso_VLoose_eta2p4.size()>0){
    l1TauIsoVLoose_pt_eta2p4->Fill(  l1PFTaus_sorted_iso_VLoose_eta2p4.at(0).pt());
    l1TauIsoVLoose_eta_eta2p4->Fill(  l1PFTaus_sorted_iso_VLoose_eta2p4.at(0).eta());
  }

  if(l1PFTaus_sorted_iso_VLoose_eta2p4.size()>1)
    l1DiTauIsoVLoose_pt_eta2p4->Fill(l1PFTaus_sorted_iso_VLoose_eta2p4.at(1).pt());

  //eta restricted 2.1
  if(l1PFTaus_sorted_iso_VLoose_eta2p1.size()>0){
    l1TauIsoVLoose_pt_eta2p1->Fill(  l1PFTaus_sorted_iso_VLoose_eta2p1.at(0).pt());
    l1TauIsoVLoose_eta_eta2p1->Fill(  l1PFTaus_sorted_iso_VLoose_eta2p1.at(0).eta());
  }

  if(l1PFTaus_sorted_iso_VLoose_eta2p1.size()>1)
    l1DiTauIsoVLoose_pt_eta2p1->Fill(l1PFTaus_sorted_iso_VLoose_eta2p1.at(1).pt());

  //eta restricted 1.4
  if(l1PFTaus_sorted_iso_VLoose_eta1p4.size()>0){
    l1TauIsoVLoose_pt_eta1p4->Fill(  l1PFTaus_sorted_iso_VLoose_eta1p4.at(0).pt());
    l1TauIsoVLoose_eta_eta1p4->Fill(  l1PFTaus_sorted_iso_VLoose_eta1p4.at(0).eta());
  }

  if(l1PFTaus_sorted_iso_VLoose_eta1p4.size()>1)
    l1DiTauIsoVLoose_pt_eta1p4->Fill(l1PFTaus_sorted_iso_VLoose_eta1p4.at(1).pt());


  // Loose finish mee
  if(l1PFTaus_sorted_iso_Loose.size()>0){
    l1TauIsoLoose_pt->Fill(l1PFTaus_sorted_iso_Loose.at(0).pt());
    l1TauIsoLoose_eta->Fill(l1PFTaus_sorted_iso_Loose.at(0).eta());
  }

  if(l1PFTaus_sorted_iso_Loose.size()>1)
    l1DiTauIsoLoose_pt->Fill(l1PFTaus_sorted_iso_Loose.at(1).pt());
  
  //eta restricted 2.4
  if(l1PFTaus_sorted_iso_Loose_eta2p4.size()>0){
    l1TauIsoLoose_pt_eta2p4->Fill(  l1PFTaus_sorted_iso_Loose_eta2p4.at(0).pt());
    l1TauIsoLoose_eta_eta2p4->Fill(  l1PFTaus_sorted_iso_Loose_eta2p4.at(0).eta());
  }

  if(l1PFTaus_sorted_iso_Loose_eta2p4.size()>1)
    l1DiTauIsoLoose_pt_eta2p4->Fill(l1PFTaus_sorted_iso_Loose_eta2p4.at(1).pt());

  //eta restricted 2.1
  if(l1PFTaus_sorted_iso_Loose_eta2p1.size()>0){
    l1TauIsoLoose_pt_eta2p1->Fill(  l1PFTaus_sorted_iso_Loose_eta2p1.at(0).pt());
    l1TauIsoLoose_eta_eta2p1->Fill(  l1PFTaus_sorted_iso_Loose_eta2p1.at(0).eta());
  }

  if(l1PFTaus_sorted_iso_Loose_eta2p1.size()>1)
    l1DiTauIsoLoose_pt_eta2p1->Fill(l1PFTaus_sorted_iso_Loose_eta2p1.at(1).pt());

  //eta restricted 1.4
  if(l1PFTaus_sorted_iso_Loose_eta1p4.size()>0){
    l1TauIsoLoose_pt_eta1p4->Fill(  l1PFTaus_sorted_iso_Loose_eta1p4.at(0).pt());
    l1TauIsoLoose_eta_eta1p4->Fill(  l1PFTaus_sorted_iso_Loose_eta1p4.at(0).eta());
  }

  if(l1PFTaus_sorted_iso_Loose_eta1p4.size()>1)
    l1DiTauIsoLoose_pt_eta1p4->Fill(l1PFTaus_sorted_iso_Loose_eta1p4.at(1).pt());


  // Medium finish mee
  if(l1PFTaus_sorted_iso_Medium.size()>0){
    l1TauIsoMedium_pt->Fill(l1PFTaus_sorted_iso_Medium.at(0).pt());
    l1TauIsoMedium_eta->Fill(l1PFTaus_sorted_iso_Medium.at(0).eta());
  }

  if(l1PFTaus_sorted_iso_Medium.size()>1)
    l1DiTauIsoMedium_pt->Fill(l1PFTaus_sorted_iso_Medium.at(1).pt());
  
  //eta restricted 2.4
  if(l1PFTaus_sorted_iso_Medium_eta2p4.size()>0){
    l1TauIsoMedium_pt_eta2p4->Fill(  l1PFTaus_sorted_iso_Medium_eta2p4.at(0).pt());
    l1TauIsoMedium_eta_eta2p4->Fill(  l1PFTaus_sorted_iso_Medium_eta2p4.at(0).eta());
  }

  if(l1PFTaus_sorted_iso_Medium_eta2p4.size()>1)
    l1DiTauIsoMedium_pt_eta2p4->Fill(l1PFTaus_sorted_iso_Medium_eta2p4.at(1).pt());

  //eta restricted 2.1
  if(l1PFTaus_sorted_iso_Medium_eta2p1.size()>0){
    l1TauIsoMedium_pt_eta2p1->Fill(  l1PFTaus_sorted_iso_Medium_eta2p1.at(0).pt());
    l1TauIsoMedium_eta_eta2p1->Fill(  l1PFTaus_sorted_iso_Medium_eta2p1.at(0).eta());
  }

  if(l1PFTaus_sorted_iso_Medium_eta2p1.size()>1)
    l1DiTauIsoMedium_pt_eta2p1->Fill(l1PFTaus_sorted_iso_Medium_eta2p1.at(1).pt());

  //eta restricted 1.4
  if(l1PFTaus_sorted_iso_Medium_eta1p4.size()>0){
    l1TauIsoMedium_pt_eta1p4->Fill(  l1PFTaus_sorted_iso_Medium_eta1p4.at(0).pt());
    l1TauIsoMedium_eta_eta1p4->Fill(  l1PFTaus_sorted_iso_Medium_eta1p4.at(0).eta());
  }

  if(l1PFTaus_sorted_iso_Medium_eta1p4.size()>1)
    l1DiTauIsoMedium_pt_eta1p4->Fill(l1PFTaus_sorted_iso_Medium_eta1p4.at(1).pt());


// Tight finish mee
  if(l1PFTaus_sorted_iso_Tight.size()>0){
    l1TauIsoTight_pt->Fill(l1PFTaus_sorted_iso_Tight.at(0).pt());
    l1TauIsoTight_eta->Fill(l1PFTaus_sorted_iso_Tight.at(0).eta());
  }

  if(l1PFTaus_sorted_iso_Tight.size()>1)
    l1DiTauIsoTight_pt->Fill(l1PFTaus_sorted_iso_Tight.at(1).pt());
  
  //eta restricted 2.4
  if(l1PFTaus_sorted_iso_Tight_eta2p4.size()>0){
    l1TauIsoTight_pt_eta2p4->Fill(  l1PFTaus_sorted_iso_Tight_eta2p4.at(0).pt());
    l1TauIsoTight_eta_eta2p4->Fill(  l1PFTaus_sorted_iso_Tight_eta2p4.at(0).eta());
  }

  if(l1PFTaus_sorted_iso_Tight_eta2p4.size()>1)
    l1DiTauIsoTight_pt_eta2p4->Fill(l1PFTaus_sorted_iso_Tight_eta2p4.at(1).pt());

  //eta restricted 2.1
  if(l1PFTaus_sorted_iso_Tight_eta2p1.size()>0){
    l1TauIsoTight_pt_eta2p1->Fill(  l1PFTaus_sorted_iso_Tight_eta2p1.at(0).pt());
    l1TauIsoTight_eta_eta2p1->Fill(  l1PFTaus_sorted_iso_Tight_eta2p1.at(0).eta());
  }

  if(l1PFTaus_sorted_iso_Tight_eta2p1.size()>1)
    l1DiTauIsoTight_pt_eta2p1->Fill(l1PFTaus_sorted_iso_Tight_eta2p1.at(1).pt());

  //eta restricted 1.4
  if(l1PFTaus_sorted_iso_Tight_eta1p4.size()>0){
    l1TauIsoTight_pt_eta1p4->Fill(  l1PFTaus_sorted_iso_Tight_eta1p4.at(0).pt());
    l1TauIsoTight_eta_eta1p4->Fill(  l1PFTaus_sorted_iso_Tight_eta1p4.at(0).eta());
  }

  if(l1PFTaus_sorted_iso_Tight_eta1p4.size()>1)
    l1DiTauIsoTight_pt_eta1p4->Fill(l1PFTaus_sorted_iso_Tight_eta1p4.at(1).pt());

  //////////////////REL ISO

  // VLoose
  if(l1PFTaus_sorted_iso_rel_VLoose.size()>0){
    l1TauIsoRelVLoose_pt->Fill(l1PFTaus_sorted_iso_rel_VLoose.at(0).pt());
    l1TauIsoRelVLoose_eta->Fill(l1PFTaus_sorted_iso_rel_VLoose.at(0).eta());
  }

  if(l1PFTaus_sorted_iso_rel_VLoose.size()>1)
    l1DiTauIsoRelVLoose_pt->Fill(l1PFTaus_sorted_iso_rel_VLoose.at(1).pt());
  
  //eta restricted 2.4
  if(l1PFTaus_sorted_iso_rel_VLoose_eta2p4.size()>0){
    l1TauIsoRelVLoose_pt_eta2p4->Fill(  l1PFTaus_sorted_iso_rel_VLoose_eta2p4.at(0).pt());
    l1TauIsoRelVLoose_eta_eta2p4->Fill(  l1PFTaus_sorted_iso_rel_VLoose_eta2p4.at(0).eta());
  }

  if(l1PFTaus_sorted_iso_rel_VLoose_eta2p4.size()>1)
    l1DiTauIsoRelVLoose_pt_eta2p4->Fill(l1PFTaus_sorted_iso_rel_VLoose_eta2p4.at(1).pt());

  //eta restricted 2.1
  if(l1PFTaus_sorted_iso_rel_VLoose_eta2p1.size()>0){
    l1TauIsoRelVLoose_pt_eta2p1->Fill(  l1PFTaus_sorted_iso_rel_VLoose_eta2p1.at(0).pt());
    l1TauIsoRelVLoose_eta_eta2p1->Fill(  l1PFTaus_sorted_iso_rel_VLoose_eta2p1.at(0).eta());
  }

  if(l1PFTaus_sorted_iso_rel_VLoose_eta2p1.size()>1)
    l1DiTauIsoRelVLoose_pt_eta2p1->Fill(l1PFTaus_sorted_iso_rel_VLoose_eta2p1.at(1).pt());

  //eta restricted 1.4
  if(l1PFTaus_sorted_iso_rel_VLoose_eta1p4.size()>0){
    l1TauIsoRelVLoose_pt_eta1p4->Fill(  l1PFTaus_sorted_iso_rel_VLoose_eta1p4.at(0).pt());
    l1TauIsoRelVLoose_eta_eta1p4->Fill(  l1PFTaus_sorted_iso_rel_VLoose_eta1p4.at(0).eta());
  }

  if(l1PFTaus_sorted_iso_rel_VLoose_eta1p4.size()>1)
    l1DiTauIsoRelVLoose_pt_eta1p4->Fill(l1PFTaus_sorted_iso_rel_VLoose_eta1p4.at(1).pt());


  // Loose finish mee
  if(l1PFTaus_sorted_iso_rel_Loose.size()>0){
    l1TauIsoRelLoose_pt->Fill(l1PFTaus_sorted_iso_rel_Loose.at(0).pt());
    l1TauIsoRelLoose_eta->Fill(l1PFTaus_sorted_iso_rel_Loose.at(0).eta());
  }

  if(l1PFTaus_sorted_iso_rel_Loose.size()>1)
    l1DiTauIsoRelLoose_pt->Fill(l1PFTaus_sorted_iso_rel_Loose.at(1).pt());
  
  //eta restricted 2.4
  if(l1PFTaus_sorted_iso_rel_Loose_eta2p4.size()>0){
    l1TauIsoRelLoose_pt_eta2p4->Fill(  l1PFTaus_sorted_iso_rel_Loose_eta2p4.at(0).pt());
    l1TauIsoRelLoose_eta_eta2p4->Fill(  l1PFTaus_sorted_iso_rel_Loose_eta2p4.at(0).eta());
  }

  if(l1PFTaus_sorted_iso_rel_Loose_eta2p4.size()>1)
    l1DiTauIsoRelLoose_pt_eta2p4->Fill(l1PFTaus_sorted_iso_rel_Loose_eta2p4.at(1).pt());

  //eta restricted 2.1
  if(l1PFTaus_sorted_iso_rel_Loose_eta2p1.size()>0){
    l1TauIsoRelLoose_pt_eta2p1->Fill(  l1PFTaus_sorted_iso_rel_Loose_eta2p1.at(0).pt());
    l1TauIsoRelLoose_eta_eta2p1->Fill(  l1PFTaus_sorted_iso_rel_Loose_eta2p1.at(0).eta());
  }

  if(l1PFTaus_sorted_iso_rel_Loose_eta2p1.size()>1)
    l1DiTauIsoRelLoose_pt_eta2p1->Fill(l1PFTaus_sorted_iso_rel_Loose_eta2p1.at(1).pt());

  //eta restricted 1.4
  if(l1PFTaus_sorted_iso_rel_Loose_eta1p4.size()>0){
    l1TauIsoRelLoose_pt_eta1p4->Fill(  l1PFTaus_sorted_iso_rel_Loose_eta1p4.at(0).pt());
    l1TauIsoRelLoose_eta_eta1p4->Fill(  l1PFTaus_sorted_iso_rel_Loose_eta1p4.at(0).eta());
  }

  if(l1PFTaus_sorted_iso_rel_Loose_eta1p4.size()>1)
    l1DiTauIsoRelLoose_pt_eta1p4->Fill(l1PFTaus_sorted_iso_rel_Loose_eta1p4.at(1).pt());


  // Medium finish mee
  if(l1PFTaus_sorted_iso_rel_Medium.size()>0){
    l1TauIsoRelMedium_pt->Fill(l1PFTaus_sorted_iso_rel_Medium.at(0).pt());
    l1TauIsoRelMedium_eta->Fill(l1PFTaus_sorted_iso_rel_Medium.at(0).eta());
  }

  if(l1PFTaus_sorted_iso_rel_Medium.size()>1)
    l1DiTauIsoRelMedium_pt->Fill(l1PFTaus_sorted_iso_rel_Medium.at(1).pt());
  
  //eta restricted 2.4
  if(l1PFTaus_sorted_iso_rel_Medium_eta2p4.size()>0){
    l1TauIsoRelMedium_pt_eta2p4->Fill(  l1PFTaus_sorted_iso_rel_Medium_eta2p4.at(0).pt());
    l1TauIsoRelMedium_eta_eta2p4->Fill(  l1PFTaus_sorted_iso_rel_Medium_eta2p4.at(0).eta());
  }

  if(l1PFTaus_sorted_iso_rel_Medium_eta2p4.size()>1)
    l1DiTauIsoRelMedium_pt_eta2p4->Fill(l1PFTaus_sorted_iso_rel_Medium_eta2p4.at(1).pt());

  //eta restricted 2.1
  if(l1PFTaus_sorted_iso_rel_Medium_eta2p1.size()>0){
    l1TauIsoRelMedium_pt_eta2p1->Fill(  l1PFTaus_sorted_iso_rel_Medium_eta2p1.at(0).pt());
    l1TauIsoRelMedium_eta_eta2p1->Fill(  l1PFTaus_sorted_iso_rel_Medium_eta2p1.at(0).eta());
  }

  if(l1PFTaus_sorted_iso_rel_Medium_eta2p1.size()>1)
    l1DiTauIsoRelMedium_pt_eta2p1->Fill(l1PFTaus_sorted_iso_rel_Medium_eta2p1.at(1).pt());

  //eta restricted 1.4
  if(l1PFTaus_sorted_iso_rel_Medium_eta1p4.size()>0){
    l1TauIsoRelMedium_pt_eta1p4->Fill(  l1PFTaus_sorted_iso_rel_Medium_eta1p4.at(0).pt());
    l1TauIsoRelMedium_eta_eta1p4->Fill(  l1PFTaus_sorted_iso_rel_Medium_eta1p4.at(0).eta());
  }

  if(l1PFTaus_sorted_iso_rel_Medium_eta1p4.size()>1)
    l1DiTauIsoRelMedium_pt_eta1p4->Fill(l1PFTaus_sorted_iso_rel_Medium_eta1p4.at(1).pt());


// Tight finish mee
  if(l1PFTaus_sorted_iso_rel_Tight.size()>0){
    l1TauIsoRelTight_pt->Fill(l1PFTaus_sorted_iso_rel_Tight.at(0).pt());
    l1TauIsoRelTight_eta->Fill(l1PFTaus_sorted_iso_rel_Tight.at(0).eta());
  }

  if(l1PFTaus_sorted_iso_rel_Tight.size()>1)
    l1DiTauIsoRelTight_pt->Fill(l1PFTaus_sorted_iso_rel_Tight.at(1).pt());
  
  //eta restricted 2.4
  if(l1PFTaus_sorted_iso_rel_Tight_eta2p4.size()>0){
    l1TauIsoRelTight_pt_eta2p4->Fill(  l1PFTaus_sorted_iso_rel_Tight_eta2p4.at(0).pt());
    l1TauIsoRelTight_eta_eta2p4->Fill(  l1PFTaus_sorted_iso_rel_Tight_eta2p4.at(0).eta());
  }

  if(l1PFTaus_sorted_iso_rel_Tight_eta2p4.size()>1)
    l1DiTauIsoRelTight_pt_eta2p4->Fill(l1PFTaus_sorted_iso_rel_Tight_eta2p4.at(1).pt());

  //eta restricted 2.1
  if(l1PFTaus_sorted_iso_rel_Tight_eta2p1.size()>0){
    l1TauIsoRelTight_pt_eta2p1->Fill(  l1PFTaus_sorted_iso_rel_Tight_eta2p1.at(0).pt());
    l1TauIsoRelTight_eta_eta2p1->Fill(  l1PFTaus_sorted_iso_rel_Tight_eta2p1.at(0).eta());
  }

  if(l1PFTaus_sorted_iso_rel_Tight_eta2p1.size()>1)
    l1DiTauIsoRelTight_pt_eta2p1->Fill(l1PFTaus_sorted_iso_rel_Tight_eta2p1.at(1).pt());

  //eta restricted 1.4
  if(l1PFTaus_sorted_iso_rel_Tight_eta1p4.size()>0){
    l1TauIsoRelTight_pt_eta1p4->Fill(  l1PFTaus_sorted_iso_rel_Tight_eta1p4.at(0).pt());
    l1TauIsoRelTight_eta_eta1p4->Fill(  l1PFTaus_sorted_iso_rel_Tight_eta1p4.at(0).eta());
  }

  if(l1PFTaus_sorted_iso_rel_Tight_eta1p4.size()>1)
    l1DiTauIsoRelTight_pt_eta1p4->Fill(l1PFTaus_sorted_iso_rel_Tight_eta1p4.at(1).pt());


  //////////////////single prong
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
  //l1PFTaus_sorted_iso_VLoose
  for(unsigned int i = 0; i < l1PFTaus_sorted_iso_VLoose.size(); i++){
    /*
    if(abs(l1PFTaus_sorted_iso_VLoose.at(i).eta())<2.1 && !f_2p1){
      l1TauIso_pt_eta2p1->Fill(l1PFTaus_sorted_iso_VLoose.at(i).pt());
      f_2p1=true;
    }

    if(abs(l1PFTaus_sorted_iso_VLoose.at(i).eta())<2.4 && !f_2p4){
      l1TauIso_pt_eta2p4->Fill(l1PFTaus_sorted_iso_VLoose.at(i).pt());
      f_2p4=true;
    }
    */
    if(l1PFTaus_sorted_iso_VLoose.at(i).tauType()==0){
      
      if(!f0){
	l1SingleProngTauIso_pt->Fill(l1PFTaus_sorted_iso_VLoose.at(i).pt());
	f0 = true;
      }

      if(abs(l1PFTaus_sorted_iso_VLoose.at(i).eta())<2.1 && !f0_2p1){
        l1SingleProngTauIso_pt_eta2p1->Fill(l1PFTaus_sorted_iso_VLoose.at(i).pt());
	f0_2p1=true;
      }
      if(abs(l1PFTaus_sorted_iso_VLoose.at(i).eta())<2.4 && !f0_2p4){
	l1SingleProngTauIso_pt_eta2p4->Fill(l1PFTaus_sorted_iso_VLoose.at(i).pt());
	f0_2p4=true;
      }
      
    }//close 1prong

    if(l1PFTaus_sorted_iso_VLoose.at(i).tauType()==1){
      
      if(!f1){
	l1SingleProngPi0TauIso_pt->Fill(l1PFTaus_sorted_iso_VLoose.at(i).pt());
	f1=true;
      }

      if(abs(l1PFTaus_sorted_iso_VLoose.at(i).eta())<2.1 && !f1_2p1){
        l1SingleProngPi0TauIso_pt_eta2p1->Fill(l1PFTaus_sorted_iso_VLoose.at(i).pt());
	f1_2p1=true;
      }
      if(abs(l1PFTaus_sorted_iso_VLoose.at(i).eta())<2.4 && !f1_2p4){
	l1SingleProngPi0TauIso_pt_eta2p4->Fill(l1PFTaus_sorted_iso_VLoose.at(i).pt());
	f1_2p4=true;
      }
      
    }//close 1prongPi0

    if(l1PFTaus_sorted_iso_VLoose.at(i).tauType()==10){
      if(!f10){
	l1ThreeProngTauIso_pt->Fill(l1PFTaus_sorted_iso_VLoose.at(i).pt());
	f10=true;
      }
      if(abs(l1PFTaus_sorted_iso_VLoose.at(i).eta())<2.1 && !f10_2p1){
        l1ThreeProngTauIso_pt_eta2p1->Fill(l1PFTaus_sorted_iso_VLoose.at(i).pt());
	f10_2p1=true;
      }
      if(abs(l1PFTaus_sorted_iso_VLoose.at(i).eta())<2.4 && !f10_2p4){
	l1ThreeProngTauIso_pt_eta2p4->Fill(l1PFTaus_sorted_iso_VLoose.at(i).pt());
	f10_2p4=true;
      }
      
    }//close 3prong

  }


  if(l1PFTaus_sorted_iso_2.size()>0){
    l1TauIsoTightV0_pt->Fill( l1PFTaus_sorted_iso_2.at(0).pt() );
    l1TauIsoTightV0_eta->Fill(l1PFTaus_sorted_iso_2.at(0).eta());
  }

  if(l1PFTaus_sorted_iso_2.size()>1){
    l1DiTauIsoTightV0_pt->Fill(l1PFTaus_sorted_iso_2.at(1).pt());
  }

  //eta restricted 2.4
  if(l1PFTaus_sorted_iso_2_eta2p4.size()>0){
    l1TauIsoTightV0_pt_eta2p4->Fill(   l1PFTaus_sorted_iso_2_eta2p4.at(0).pt());
    l1TauIsoTightV0_eta_eta2p4->Fill(  l1PFTaus_sorted_iso_2_eta2p4.at(0).eta());
  }

  if(l1PFTaus_sorted_iso_2_eta2p4.size()>1)
    l1DiTauIsoTightV0_pt_eta2p4->Fill(l1PFTaus_sorted_iso_2_eta2p4.at(1).pt());

  //eta restricted 2.1
  if(l1PFTaus_sorted_iso_2_eta2p1.size()>0){
    l1TauIsoTightV0_pt_eta2p1->Fill(   l1PFTaus_sorted_iso_2_eta2p1.at(0).pt());
    l1TauIsoTightV0_eta_eta2p1->Fill(  l1PFTaus_sorted_iso_2_eta2p1.at(0).eta());
  }

  if(l1PFTaus_sorted_iso_2_eta2p1.size()>1)
    l1DiTauIsoTightV0_pt_eta2p1->Fill(l1PFTaus_sorted_iso_2_eta2p1.at(1).pt());

  //eta restricted 1p4
  if(l1PFTaus_sorted_iso_2_eta1p4.size()>0){
    l1TauIsoTightV0_pt_eta1p4->Fill(   l1PFTaus_sorted_iso_2_eta1p4.at(0).pt());
    l1TauIsoTightV0_eta_eta1p4->Fill(  l1PFTaus_sorted_iso_2_eta1p4.at(0).eta());
  }

  if(l1PFTaus_sorted_iso_2_eta1p4.size()>1)
    l1DiTauIsoTightV0_pt_eta1p4->Fill(l1PFTaus_sorted_iso_2_eta1p4.at(1).pt());

  //std::cout<<"here 5"<<std::endl;  
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

  for(unsigned int i = 0; i < l1PFTaus_sorted_iso_2.size(); i++){
    /*
    if(abs(l1PFTaus_sorted_iso_2.at(i).eta())<2.1 && !f_2p1){
      l1TauIsoTight_pt_eta2p1->Fill(l1PFTaus_sorted_iso_2.at(i).pt());
      f_2p1=true;
    }

    if(abs(l1PFTaus_sorted_iso_2.at(i).eta())<2.4 && !f_2p4){
      l1TauIsoTight_pt_eta2p4->Fill(l1PFTaus_sorted_iso_2.at(i).pt());
      f_2p4=true;
    }
    */
    if(l1PFTaus_sorted_iso_2.at(i).tauType()==0){
      
      if(!f0){
	l1SingleProngTauIsoTight_pt->Fill(l1PFTaus_sorted_iso_2.at(i).pt());
	f0 = true;
      }

      if(abs(l1PFTaus_sorted_iso_2.at(i).eta())<2.1 && !f0_2p1){
        l1SingleProngTauIsoTight_pt_eta2p1->Fill(l1PFTaus_sorted_iso_2.at(i).pt());
	f0_2p1=true;
      }
      if(abs(l1PFTaus_sorted_iso_2.at(i).eta())<2.4 && !f0_2p4){
	l1SingleProngTauIsoTight_pt_eta2p4->Fill(l1PFTaus_sorted_iso_2.at(i).pt());
	f0_2p4=true;
      }
      
    }//close 1prong

    if(l1PFTaus_sorted_iso_2.at(i).tauType()==1){
      
      if(!f1){
	l1SingleProngPi0TauIsoTight_pt->Fill(l1PFTaus_sorted_iso_2.at(i).pt());
	f1=true;
      }

      if(abs(l1PFTaus_sorted_iso_2.at(i).eta())<2.1 && !f1_2p1){
        l1SingleProngPi0TauIsoTight_pt_eta2p1->Fill(l1PFTaus_sorted_iso_2.at(i).pt());
	f1_2p1=true;
      }
      if(abs(l1PFTaus_sorted_iso_2.at(i).eta())<2.4 && !f1_2p4){
	l1SingleProngPi0TauIsoTight_pt_eta2p4->Fill(l1PFTaus_sorted_iso_2.at(i).pt());
	f1_2p4=true;
      }
      
    }//close 1prongPi0

    if(l1PFTaus_sorted_iso_2.at(i).tauType()==10){
      if(!f10){
	l1ThreeProngTauIsoTight_pt->Fill(l1PFTaus_sorted_iso_2.at(i).pt());
	f10=true;
      }
      if(abs(l1PFTaus_sorted_iso_2.at(i).eta())<2.1 && !f10_2p1){
        l1ThreeProngTauIsoTight_pt_eta2p1->Fill(l1PFTaus_sorted_iso_2.at(i).pt());
	f10_2p1=true;
      }
      if(abs(l1PFTaus_sorted_iso_2.at(i).eta())<2.4 && !f10_2p4){
	l1ThreeProngTauIsoTight_pt_eta2p4->Fill(l1PFTaus_sorted_iso_2.at(i).pt());
	f10_2p4=true;
      }
      
    }//close 3prong


  }

  //std::cout<<"here 6"<<std::endl;

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
