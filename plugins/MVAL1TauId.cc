#include "MVAL1TauId.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"


#include "CommonTools/MVAUtils/interface/TMVAZipReader.h"

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

void MVAL1TauId::bookReader() {
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

void MVAJetPuId::set(const PileupJetIdentifier &id) { internalId_ = id; }

void MVAL1TauId::runMva() {
  if (!reader_) {
    bookReader();
  }
  
  internalId_.mva = reader_->EvaluateMVA(tmvaMethod_.c_str());
  
}

#define INIT_VARIABLE(NAME, TMVANAME, VAL) \
  internalId_.NAME##_ = VAL;               \
  variables_[#NAME] = std::make_pair(&internalId_.NAME##_, VAL);


void MVAL1TauId::initVariables() {
  INIT_VARIABLE(mva, "", -100.);

  INIT_VARIABLE(l1Pt, "l1Pt", 0);
  INIT_VARIABLE(genPt, "genPt", 0);

}

#undef INIT_VARIABLE
