#include "AliMuonCompactQuickAccEff.h"

#include <cassert>
#include <iostream>
#include "AliAnalysisRunList.h"
#include "AliMuonCompactEvent.h"
#include "AliMuonCompactManuStatus.h"
#include "AliMuonCompactMapping.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TMath.h"
#include "TParameter.h"
#include "TTree.h"

namespace
{
void GetDropPercentage(int npairs, int referenceNofJpsi, double& drop, double& dropError)
{
  drop = 100.0 * (1.0 - npairs * 1.0 / referenceNofJpsi);
  Double_t relativeError = TMath::Sqrt(1.0 / npairs + 1.0 / referenceNofJpsi);
  dropError = drop * relativeError;
}

std::vector<AliMuonCompactEvent> GetEvents(TTree* tree)
{
  /// Read events from the tree
  std::vector<AliMuonCompactEvent> events;
  AliMuonCompactEvent* compactEvent{ nullptr };
  tree->SetBranchAddress("event", &compactEvent);

  for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);
    events.push_back(*compactEvent);
  }
  return events;
}

std::vector<AliMuonCompactEvent> GetEvents(const char* treeFile)
{
  std::unique_ptr<TFile> f{ TFile::Open(treeFile) };
  if (!f->IsOpen())
    return {};

  TTree* tree = static_cast<TTree*>(f->Get("compactevents"));
  if (!tree)
    return {};

  return GetEvents(tree);
}

} // namespace

/// \ingroup compact
AliMuonCompactQuickAccEff::AliMuonCompactQuickAccEff(const char* eventFile, ULong64_t maxevents,
                                                     bool rejectMonoCathodeClusters)
  : fMaxEvents{ maxevents }, fRejectMonoCathodeClusters{ rejectMonoCathodeClusters }, mEvents{ GetEvents(eventFile) }
{
}

Bool_t AliMuonCompactQuickAccEff::ValidateCluster(const AliMuonCompactCluster& cl,
                                                  const AliMuonCompactManuStatus& manuStatus, UInt_t causeMask)
{
  UInt_t bendingMask = 0;
  UInt_t nonBendingMask = 0;

  if (cl.BendingManuIndex() >= 0) {
    assert(cl.BendingManuIndex() <= (int)manuStatus.size());
    bendingMask = manuStatus[cl.BendingManuIndex()];
  }
  if (cl.NonBendingManuIndex() >= 0) {
    assert(cl.NonBendingManuIndex() <= (int)manuStatus.size());
    nonBendingMask = manuStatus[cl.NonBendingManuIndex()];
  }

  Bool_t station12 = (cl.BendingManuIndex() >= 0 && cl.BendingManuIndex() < 7152) ||
                     (cl.NonBendingManuIndex() >= 0 && cl.NonBendingManuIndex() < 7152);

  Bool_t bendingIsOK = ((bendingMask & causeMask) == 0);
  Bool_t nonBendingIsOK = ((nonBendingMask & causeMask) == 0);

  if (fRejectMonoCathodeClusters) {
    // note that in most of the cases removing the mono-cathode clusters
    // lead to a worst reproduction of the full simulation, simply
    // because mono-cathode are not explicitely killed by the
    // regular reconstruction.
    // we keep the "option" here however for the record.
    if (station12) {
      // for station12 it's ok to have a monocathode cluster
      return (bendingIsOK || nonBendingIsOK);
    } else {
      // for station345 it is *not* ok to have only non bending
      return (bendingIsOK && nonBendingIsOK);
    }
  }

  return (bendingIsOK || nonBendingIsOK);
}

Bool_t AliMuonCompactQuickAccEff::ValidateTrack(const AliMuonCompactTrack& track,
                                                const AliMuonCompactManuStatus& manuStatus, UInt_t causeMask)
{
  /// We remove from the track all the clusters
  /// located on a bad manu.
  /// Then we consider the track survived if we get
  /// at least one cluster per station

  if (manuStatus.AllClean() || causeMask == 0)
    return kTRUE;

  Int_t currentCh;
  Int_t currentSt;
  Int_t previousCh = -1;
  Int_t nChHitInSt4 = 0;
  Int_t nChHitInSt5 = 0;
  UInt_t presentStationMask = 0;
  const UInt_t requestedStationMask = 0x1F;
  const Bool_t request2ChInSameSt45 = kTRUE;

  /* CompactMapping* cm = GetCompactMapping(); */

  for (std::vector<AliMuonCompactCluster>::size_type i = 0; i < track.mClusters.size(); ++i) {
    const AliMuonCompactCluster& cl = track.mClusters[i];

    if (!ValidateCluster(cl, manuStatus, causeMask)) {
      continue;
    }

    currentCh = cl.DetElemId() / 100 - 1;
    currentSt = currentCh / 2;

    // build present station mask
    presentStationMask |= (1 << currentSt);

    // count the number of chambers hit in station 4 that contain cluster(s)
    if (currentSt == 3 && currentCh != previousCh) {
      ++nChHitInSt4;
      previousCh = currentCh;
    }

    // count the number of chambers hit in station 5 that contain cluster(s)
    if (currentSt == 4 && currentCh != previousCh) {
      ++nChHitInSt5;
      previousCh = currentCh;
    }
  }

  // at least one cluster per requested station
  if ((requestedStationMask & presentStationMask) != requestedStationMask) {
    return kFALSE;
  }

  if (request2ChInSameSt45) {
    // 2 chambers hit in the same station (4 or 5)
    return (nChHitInSt4 == 2 || nChHitInSt5 == 2);
  } else {
    // or 2 chambers hit in station 4 & 5 together
    return (nChHitInSt4 + nChHitInSt5 >= 2);
  }

  return kTRUE;
}

Int_t AliMuonCompactQuickAccEff::ComputeNofPairs(const AliMuonCompactManuStatus& manustatus, UInt_t causeMask)
{
  Int_t npairs = 0;

  const double m2 = 0.1056584 * 0.1056584;
  Int_t nTracks = 0;
  Int_t nValidatedTracks = 0;

  uint64_t maxevents = fMaxEvents;

  if (!maxevents) {
    maxevents = mEvents.size();
  }

  for (auto i = 0; i < maxevents; ++i) {
    const AliMuonCompactEvent& e = mEvents[i];

    for (std::vector<AliMuonCompactTrack>::size_type j = 0; j < e.mTracks.size(); ++j) {
      const AliMuonCompactTrack& t1 = e.mTracks[j];

      ++nTracks;
      if (!ValidateTrack(t1, manustatus, causeMask))
        continue;
      ++nValidatedTracks;

      for (std::vector<AliMuonCompactTrack>::size_type k = j + 1; k < e.mTracks.size(); ++k) {
        const AliMuonCompactTrack& t2 = e.mTracks[k];

        if (!ValidateTrack(t2, manustatus, causeMask))
          continue;

        double p1square = t1.mPx * t1.mPx + t1.mPy * t1.mPy + t1.mPz * t1.mPz;

        double p2square = t2.mPx * t2.mPx + t2.mPy * t2.mPy + t2.mPz * t2.mPz;

        double e = sqrt(m2 + p1square + p2square + 2.0 * sqrt(p1square) * sqrt(p2square));
        double pz = t1.mPz + t2.mPz;

        double y = 0.5 * log((e + pz) / (e - pz));

        if (y >= -4 && y <= -2.5) {
          ++npairs;
        }
      }
    }
  }

  std::cout << Form("nTracks %d nValidated %d npairs %d", nTracks, nValidatedTracks, npairs) << std::endl;

  return npairs;
}

void AliMuonCompactQuickAccEff::ComputeEvolution(std::vector<int>& vrunlist,
                                                 const std::map<int, AliMuonCompactManuStatus>& manuStatusForRuns,
                                                 const char* outputfile)
{
  std::cout << "ComputeEvolution(const std::vector<AliMuonCompactEvent>& events,...)" << std::endl;
  RefValues ref = GetRef();

  std::vector<TGraphErrors*> gdrop;
  // one graph for each "bad" cause (but on
  // manu level only)
  // - ped is a bit ill-defined for manu, we consider a manu
  // "bad for ped" if 70% of its channels are bad
  // - manu occupancy
  // - hv
  // - lv
  // - missing (i.e. buspatch removed from configuration)

  std::vector<UInt_t> causes;

  causes.push_back(AliMuonCompactManuStatus::MANUOUTOFCONFIGMASK);
  causes.push_back(AliMuonCompactManuStatus::MANUOUTOFCONFIGMASK | AliMuonCompactManuStatus::MANUBADHVMASK);
  causes.push_back(AliMuonCompactManuStatus::MANUOUTOFCONFIGMASK | AliMuonCompactManuStatus::MANUBADPEDMASK);
  causes.push_back(AliMuonCompactManuStatus::MANUOUTOFCONFIGMASK | AliMuonCompactManuStatus::MANUBADOCCMASK);
  causes.push_back(AliMuonCompactManuStatus::MANUOUTOFCONFIGMASK | AliMuonCompactManuStatus::MANUBADPEDMASK |
                   AliMuonCompactManuStatus::MANUBADOCCMASK | AliMuonCompactManuStatus::MANUBADHVMASK);
  causes.push_back(AliMuonCompactManuStatus::MANUOUTOFCONFIGMASK | AliMuonCompactManuStatus::MANUBADPEDMASK |
                   AliMuonCompactManuStatus::MANUBADOCCMASK | AliMuonCompactManuStatus::MANUBADHVMASK |
                   AliMuonCompactManuStatus::MANUREJECTMASK);
  causes.push_back(AliMuonCompactManuStatus::MANUOUTOFCONFIGMASK | AliMuonCompactManuStatus::MANUBADPEDMASK |
                   AliMuonCompactManuStatus::MANUBADOCCMASK | AliMuonCompactManuStatus::MANUBADHVMASK |
                   AliMuonCompactManuStatus::MANUBADLVMASK);
  for (std::vector<UInt_t>::size_type i = 0; i < causes.size(); ++i) {
    TGraphErrors* g = new TGraphErrors(vrunlist.size());
    gdrop.push_back(g);
    g->SetName(Form("acceffdrop%s", AliMuonCompactManuStatus::CauseAsString(causes[i]).c_str()));
    g->SetMarkerStyle(20);
    g->SetMarkerSize(1.5);
  }

  for (std::vector<int>::size_type i = 0; i < vrunlist.size(); ++i) {
    Int_t runNumber = vrunlist[i];

    std::cout << Form("---- RUN %6d", runNumber) << std::endl;

    auto it = manuStatusForRuns.find(runNumber);

    const auto& manuStatus = it->second;

    for (auto icause = 0; icause < causes.size(); ++icause) {
      auto nbad = manuStatus.CountBad(causes[icause]);

      std::cout << Form("RUN %6d %30s rejected manus = %6d => ", runNumber,
                        AliMuonCompactManuStatus::CauseAsString(causes[icause]).c_str(), nbad);
      Int_t npairs = ComputeNofPairs(manuStatus, causes[icause]);
      double drop, dropError;
      GetDropPercentage(npairs, ref.nJpsi, drop, dropError);
      std::cout << Form("RUN %6d %30s AccxEff drop %7.2f %% +- %5.2f %%", runNumber, " ", drop, dropError) << std::endl;
      gdrop[icause]->SetPoint(i, runNumber, drop);
      gdrop[icause]->SetPointError(i, 0.0, dropError);
    }
  }

  TFile* fout = TFile::Open(outputfile, "recreate");
  for (std::vector<UInt_t>::size_type icause = 0; icause < causes.size(); ++icause) {
    gdrop[icause]->Write();
  }

  std::cout << "RefNofJpsi      = " << rv.nJpsi << std::endl;
  std::cout << "RefNofInputJpsi = " << rv.nInputJpsi << std::endl;
  std::cout << "NofEvents       = " << rv.nEvents << std::endl;
  std::cout << "RefAccEff       = " << rv.AccEff << " +- " << rv.AccEffError << std::endl;

  TParameter<Double_t>("RefNofJpsi", rv.nofJpsi).Write();
  TParameter<Double_t>("RefNofInputJpsi", rv.nInputJpsi).Write();
  TParameter<Double_t>("NofEvents", rv.nEvents).Write();
  TParameter<Double_t>("RefAccEff", rv.AccEff).Write();
  TParameter<Double_t>("RefAccEffError", rv.AccEffError).Write();

  delete fout;
}

void AliMuonCompactQuickAccEff::ComputeEvolutionFromManuStatus(const char* treeFile, const char* runlist,
                                                               const char* outputfile, const char* manustatusfile,
                                                               const char* ocdbPath, Int_t runNumber)
{
  AliMuonCompactMapping::GetCompactMapping(ocdbPath, runNumber);

  std::map<int, AliMuonCompactManuStatus> manuStatusForRuns;
  AliMuonCompactManuStatus::ReadManuStatus(manustatusfile, manuStatusForRuns);

  AliAnalysisRunList rl(runlist);
  std::vector<int> vrunlist = rl.AsVector();
  ComputeEvolution(vrunlist, manuStatusForRuns, outputfile);
}

void AliMuonCompactQuickAccEff::ComputeLoss(const AliMuonCompactManuStatus& ref, const AliMuonCompactManuStatus& ms,
                                            UInt_t refCauseMask)
{
  int npairsRef = ComputeNofPairs(ref, 0);
  int npairs = ComputeNofPairs(ms, refCauseMask);

  double drop, dropError;
  GetDropPercentage(npairs, npairsRef, drop, dropError);
  std::cout << Form("drop = %7.2f %% +- %5.2f %%\n", drop, dropError);
}

void AliMuonCompactQuickAccEff::PrintEvents(int nmax)
{
  int nevents = nmax ? nmax : mEvents.size();

  for (auto i = 0; i < nevents; i++) {
    std::cout << mEvents[i] << std::endl;
  }
}

RefValues AliMuonCompactQuickAccEff::AliMuonCompactQuickAccEff::GetRef()
{
  RefValues rv;

  rv.nofJpsi = ComputeNofPairs(AliMuonCompactManuStatus{}, 0);

  // we count the number of input Jpsi which are in the correct rapidity range
  rv.nInputJpsi =
    std::count_if(mEvents.begin(), end, [](const AliMuonCompactEvent& e) { return e.mY >= -4 && e.mY <= -2.5; });

  rv.nEvents = end - mEvents.begin();
  rv.AccEff = rv.nJpsi / (1.0 * rv.nInputJpsi);
  rv.AccEffErr = TMath::Sqrt(1.0 / rv.nofJpsi + 1.0 / rv.nInputJpsi) * rv.AccEff;

  return rv;
}
