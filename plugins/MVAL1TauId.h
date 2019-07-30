#ifndef L1Trigger_phase2L1TauAnalyzer_plugins_MVAL1TauId_h
#define L1Trigger_phase2L1TauAnalyzer_plugins_MVAL1TauId_h

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/L1Trigger/interface/L1PFTau.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/L1TrackTrigger/interface/L1TkPrimaryVertex.h"


class MVAL1TauId {

 public:
  MVAL1TauId(const std::string &tmvaWeight = "",
             const std::string &tmvaMethod = "",
             Float_t impactParTkThreshod_ = 1.,
             const std::vector<std::string> &tmvaVariables = std::vector<std::string>());

  MVAL1TauId(const edm::ParameterSet &ps);
  ~MVAL1TauId();

  const std::string method() const { return tmvaMethod_; }

  typedef std::map<std::string, std::pair<float *, float> > variables_list_t;

 protected:
  void setup();
  void runMva();
  void bookReader();
  void resetVariables();
  void initVariables();

  variables_list_t variables_;

  TMVA::Reader *reader_;
  std::string tmvaWeights_, tmvaMethod_;
  std::vector<std::string> tmvaVariables_;
  std::vector<std::string> tmvaSpectators_;
  std::map<std::string, std::string> tmvaNames_;

  Int_t version_;
  Float_t impactParTkThreshod_;
  bool cutBased_;



};

#endif

