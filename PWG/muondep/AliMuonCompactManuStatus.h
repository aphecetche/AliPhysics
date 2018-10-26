#ifndef ALIMUONCOMPACTMANUSTATUS_H
#define ALIMUONCOMPACTMANUSTATUS_H

#include <map>
#include <array>
#include <vector>
#include "Rtypes.h"

/**

@ingroup pwg_muondep_compact

@class AliMuonCompactManuStatus

@brief Utility class to compute status of MCH manus

*/

class AliMuonCompactManuStatus
{
public:

    static const UInt_t MANUBADPEDMASK;
    static const UInt_t MANUBADHVMASK;
    static const UInt_t MANUBADLVMASK;
    static const UInt_t MANUBADOCCMASK;
    static const UInt_t MANUOUTOFCONFIGMASK;
    static const UInt_t MANUREJECTMASK;

    AliMuonCompactManuStatus() {}

    AliMuonCompactManuStatus(int runNumber, const char* ocdbPath="raw://");

    AliMuonCompactManuStatus(std::array<UInt_t, 16828> manuStatus) : mManuStatus(manuStatus) {}

    AliMuonCompactManuStatus(const AliMuonCompactManuStatus&) = default;
    AliMuonCompactManuStatus& operator=(const AliMuonCompactManuStatus&) = default;

    static std::string CauseAsString(UInt_t cause);

    bool AllClean() const;

    UInt_t operator[](int i) const { return mManuStatus[i]; }

    void Print(bool all = false);

    AliMuonCompactManuStatus RemoveBusPatches(std::vector<int> busPatchIds) const;

    static void ReadManuStatus(const char* inputfile, std::map<int,AliMuonCompactManuStatus>& manuStatusForRuns);
    static void WriteToBinaryFile(const char* runlist, const char* outputfile, const char* ocdbPath="raw://");

    int CountBad(UInt_t cause) const;

private:
    std::array<UInt_t,16828> mManuStatus; 
};

#endif

