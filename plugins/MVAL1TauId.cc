#include "MVAL1TauId.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

MVAL1TauId::MVAL1TauId(const edm::ParameterSet &ps) {

  impactParTkThreshod_ = 1.;
  tmvaWeights_ = edm::FileInPath(ps.getParameter<std::string>("tmvaWeights")).fullPath();
  tmvaMethod_ = ps.getParameter<std::string>("tmvaMethod");
  tmvaVariables_ = ps.getParameter<std::vector<std::string> >("tmvaVariables");
  tmvaSpectators_ = ps.getParameter<std::vector<std::string> >("tmvaSpectators");
  reader_ = nullptr;

  setup();

}

void MVAL1TauId::setup() {
  initVariables();

  tmvaVariables_.clear();
  tmvaVariables_.push_back("l1Pt");
  
  tmvaNames_["l1Pt"] = "l1Pt";

  
}


MVAL1TauId::~MVAL1TauId() {
  if (reader_)
    delete reader_;
}

void MVAJetPuId::bookReader() {
  reader_ = new TMVA::Reader("!Color:Silent");
  assert(!tmvaMethod_.empty() && !tmvaWeights_.empty());
  for (std::vector<std::string>::iterator it = tmvaVariables_.begin(); it != tmvaVariables_.end(); ++it) {
    if (tmvaNames_[*it].empty()) {
      tmvaNames_[*it] = *it;
    }
    reader_->AddVariable(*it, variables_[tmvaNames_[*it]].first);
  }
  for (std::vector<std::string>::iterator it = tmvaSpectators_.begin(); it != tmvaSpectators_.end(); ++it) {
    if (tmvaNames_[*it].empty()) {
      tmvaNames_[*it] = *it;
    }
    reader_->AddSpectator(*it, variables_[tmvaNames_[*it]].first);
  }
  reco::details::loadTMVAWeights(reader_, tmvaMethod_, tmvaWeights_);
}

void MVAJetPuId::runMva() {
  if (!reader_) {
    bookReader();
  }
  if (fabs(internalId_.jetEta_) < 5.0)
    internalId_.mva_ = reader_->EvaluateMVA(tmvaMethod_.c_str());
  if (fabs(internalId_.jetEta_) >= 5.0)
    internalId_.mva_ = -2.;
  internalId_.idFlag_ = computeIDflag(internalId_.mva_, internalId_.jetPt_, internalId_.jetEta_);
}
