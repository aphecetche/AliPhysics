#ifndef ALIMUONCOMPACTQUICKACCEFF_H
#define ALIMUONCOMPACTQUICKACCEFF_H

#include <map>
#include <vector>
#include "Rtypes.h"

struct AliMuonCompactCluster;
struct AliMuonCompactEvent;
struct AliMuonCompactTrack;
class AliMuonCompactManuStatus;
class TH1;
class TTree;

/**

  @ingroup pwg_muondep_compact

  @class AliMuonCompactQuickAccEff

  @brief Quick Acc x Eff processor

  This class is meant to get a quick computation of
  the evolution of the Acc x Eff for some runs.

*/

class AliMuonCompactQuickAccEff
{
 public:
  AliMuonCompactQuickAccEff(const char* compactEventFile, ULong64_t maxevents = 0,
                            bool rejectMonoCathodeClusters = false);

  void ComputeEvolution(std::vector<int>& vrunlist, const std::map<int, AliMuonCompactManuStatus>& manuStatusForRuns,
                        const char* outputfile);

  Bool_t ValidateCluster(const AliMuonCompactCluster& cl, const AliMuonCompactManuStatus& manuStatus, UInt_t causeMask);

  Bool_t ValidateTrack(const AliMuonCompactTrack& track, const AliMuonCompactManuStatus& manuStatus, UInt_t causeMask);

  Int_t ComputeNofPairs(const AliMuonCompactManuStatus& manustatus, UInt_t causeMask);

  void ComputeEvolutionFromManuStatus(const char* treeFile, const char* runList, const char* outputfile,
                                      const char* manustatusfile, const char* ocdbPath = "raw://", Int_t runNumber = 0);

  void ComputeLoss(const AliMuonCompactManuStatus& ref, const AliMuonCompactManuStatus& ms, UInt_t refCauseMask = 63);

  void PrintEvents(int nmax = 0);

  struct RefValues {
    int nJpsi;
    int nInputJpsi;
    int nEvents;
    double AccEff;
    double AccEffErr;
  };

  RefValues GetRef();

 private:
  ULong64_t fMaxEvents;
  bool fRejectMonoCathodeClusters;
  std::vector<AliMuonCompactEvent> mEvents;
};

#endif
