#include "AliMuonCompactManuStatus.h"

#include "AliAnalysisTriggerScalers.h"
#include "AliCDBManager.h"
#include "AliMUONCDB.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONPadStatusMaker.h"
#include "AliMUONRecoParam.h"
#include "AliMUONRejectList.h"
#include "AliMpCDB.h"
#include "AliMpConstants.h"
#include "AliMpDDLStore.h"
#include "AliMpDetElement.h"
#include "AliMpManuIterator.h"
#include "AliMuonCompactMapping.h"
#include <cassert>
#include <algorithm>

/// \ingroup compact
const UInt_t AliMuonCompactManuStatus::MANUBADPEDMASK = ( 1 << 0 );
const UInt_t AliMuonCompactManuStatus::MANUBADHVMASK = ( 1 << 1 );
const UInt_t AliMuonCompactManuStatus::MANUBADLVMASK =  (1 << 2 );
const UInt_t AliMuonCompactManuStatus::MANUBADOCCMASK = ( 1 << 3 );
const UInt_t AliMuonCompactManuStatus::MANUOUTOFCONFIGMASK = ( 1 << 4 );
const UInt_t AliMuonCompactManuStatus::MANUREJECTMASK = ( 1 << 5 );

namespace {

std::array<UInt_t,16828> BuildFromOCDB(Int_t runNumber, const char* ocdbPath)
{
    std::array<UInt_t,16828> vManuStatus;

    AliCDBManager* man = AliCDBManager::Instance();
    man->SetDefaultStorage(ocdbPath);
    man->SetRun(runNumber);

    if (!AliMpDDLStore::Instance())
    {
        AliMpCDB::LoadAll();
    }
   
    AliMuonCompactMapping* cm = AliMuonCompactMapping::GetCompactMapping(ocdbPath,runNumber);

    AliMUONCalibrationData cd(runNumber,true);

    AliMUONRejectList* rl = cd.RejectList();
    assert(rl->IsBinary());

    AliMUONPadStatusMaker statusMaker(cd);

    AliMUONRecoParam* recoParam = AliMUONCDB::LoadRecoParam();

    statusMaker.SetLimits(*recoParam);

    AliMpManuIterator it;
    Int_t detElemId, manuId;

    Int_t pedCheck = (
            AliMUONPadStatusMaker::kPedMeanZero |
            AliMUONPadStatusMaker::kPedMeanTooLow |
            AliMUONPadStatusMaker::kPedMeanTooHigh |
            AliMUONPadStatusMaker::kPedSigmaTooLow |
            AliMUONPadStatusMaker::kPedSigmaTooHigh ); /* 62 */
    
    Int_t hvCheck = (
            AliMUONPadStatusMaker::kHVError |
            AliMUONPadStatusMaker::kHVTooLow |
            AliMUONPadStatusMaker::kHVTooHigh |
            AliMUONPadStatusMaker::kHVChannelOFF |
            AliMUONPadStatusMaker::kHVSwitchOFF ); /* 31 */
    
    Int_t occCheck = (
            AliMUONPadStatusMaker::kManuOccupancyTooHigh 
            ); /* 4 */
    
    Int_t lvCheck = ( AliMUONPadStatusMaker::kLVTooLow ); /* 8 */

    Int_t ntotal(0);

    while ( it.Next(detElemId,manuId) )
    {
        AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);
        Int_t busPatchId = AliMpDDLStore::Instance()->GetBusPatchId(detElemId,manuId);
        
        UInt_t manuStatus = 0;

        Int_t manubadped=0;
        Int_t manubadocc=0;
        Int_t manubadhv=0;
        Int_t manubadlv=0;
        Int_t manumissing=0;
        Int_t manureject=0;

        for ( Int_t manuChannel = 0; manuChannel < AliMpConstants::ManuNofChannels(); ++manuChannel )
        {
            if ( de->IsConnectedChannel(manuId,manuChannel) )
            {
                ++ntotal;

                UInt_t status = statusMaker.PadStatus(detElemId, manuId, manuChannel);

                //if (!status) continue;

                if ( status & AliMUONPadStatusMaker::BuildStatus(pedCheck,0,0,0) )
                {
                    ++manubadped;
                }

                if ( status & AliMUONPadStatusMaker::BuildStatus(0,hvCheck,0,0) )
                {
                    ++manubadhv;
                }

                if ( status & AliMUONPadStatusMaker::BuildStatus(0,0,lvCheck,0) )
                {
                    ++manubadlv;
                }

                if ( status & AliMUONPadStatusMaker::BuildStatus(0,0,0,occCheck) )
                {
                    ++manubadocc;
                }

                if ( status & AliMUONPadStatusMaker::BuildStatus(128 /*AliMUONPadStatusMaker::kMissing*/,0,0,0) )
                {
                    ++manumissing;
                }

                Float_t proba = TMath::Max(rl->DetectionElementProbability(detElemId),rl->BusPatchProbability(busPatchId));
                proba = TMath::Max(proba,rl->ManuProbability(detElemId,manuId));
                proba = TMath::Max(proba,rl->ChannelProbability(detElemId,manuId,manuChannel));

                if ( proba > 0 )
                {
                    ++manureject;
                }
            }

            if ( manubadped>=0.9*de->NofChannelsInManu(manuId) )
            {
                manuStatus |= AliMuonCompactManuStatus::MANUBADPEDMASK;
            }

            if ( manubadhv )
            {
                manuStatus |= AliMuonCompactManuStatus::MANUBADHVMASK;
            }

            if ( manubadlv )
            {
                manuStatus |= AliMuonCompactManuStatus::MANUBADLVMASK;
            }
            
            if ( manubadocc ) 
            {
                manuStatus |= AliMuonCompactManuStatus::MANUBADOCCMASK;
            }
            
            if ( manumissing) 
            {
                manuStatus |= AliMuonCompactManuStatus::MANUOUTOFCONFIGMASK;
            }
            
            if ( manureject >= 0.9*de->NofChannelsInManu(manuId) )
            {
                manuStatus |= AliMuonCompactManuStatus::MANUREJECTMASK;

            }
            
            Int_t manuAbsIndex = cm->FindManuAbsIndex(detElemId,manuId);
            vManuStatus[manuAbsIndex] = manuStatus;
        }
    }

    assert(ntotal==1064008);

    man->ClearCache();

    return vManuStatus;
}
}

AliMuonCompactManuStatus::AliMuonCompactManuStatus(int runNumber, const char* ocdbPath) : mManuStatus(BuildFromOCDB(runNumber,ocdbPath))
{
}

std::string AliMuonCompactManuStatus::CauseAsString(UInt_t cause)
{
    std::string rv = "";

    if ( cause & MANUBADPEDMASK ) rv += "_ped";
    if ( cause & MANUBADHVMASK ) rv += "_hv";
    if ( cause & MANUBADLVMASK ) rv += "_lv";
    if ( cause & MANUBADOCCMASK ) rv += "_occ";
    if ( cause & MANUOUTOFCONFIGMASK ) rv += "_config";
    if ( cause & MANUREJECTMASK ) rv += "_reject";

    return rv;
}

void AliMuonCompactManuStatus::Print(bool all)
{
    AliMuonCompactMapping* cm = AliMuonCompactMapping::GetCompactMapping();

    for ( auto i = 0; i < mManuStatus.size(); ++i )
    {
        Int_t absManuId = cm->AbsManuId(i);
        Int_t detElemId = cm->GetDetElemIdFromAbsManuId(absManuId);
        Int_t manuId = cm->GetManuIdFromAbsManuId(absManuId);
        Int_t busPatchId = AliMpDDLStore::Instance()->GetBusPatchId(detElemId,manuId);
        if ( mManuStatus[i] || all )
        {
            std::cout << Form("status[%04d]=%6x (DE %04d BP %04d MANU %04d) %s",i,mManuStatus[i],detElemId,busPatchId,manuId,CauseAsString(mManuStatus[i]).c_str()) << std::endl;
        }
    }
}

void AliMuonCompactManuStatus::WriteToBinaryFile(const char* runlist, const char* outputfile, const char* ocdbPath)
{
    AliAnalysisTriggerScalers ts(runlist,ocdbPath);

    std::vector<int> vrunlist = ts.GetRunList(); // FIXME: should need to bring in all the AliAnalysisTriggerScalers class just to read the runlist... 

    std::ofstream out(outputfile,std::ios::binary);

    std::vector<int>::size_type nruns = vrunlist.size();

    out.write((char*)&nruns,sizeof(int));
    out.write((char*)&vrunlist[0],nruns*sizeof(int));
    
    for ( std::vector<int>::size_type i = 0; i < vrunlist.size(); ++i )
    {
        Int_t runNumber = vrunlist[i];
        auto manuStatus = BuildFromOCDB(runNumber,ocdbPath);
        out.write((char*)&manuStatus[0],
                manuStatus.size()*sizeof(int));
        assert(manuStatus.size()==16828);
        std::cout << Form("RUN %6d", runNumber) << std::endl;
    }
    out.close();
}

void AliMuonCompactManuStatus::ReadManuStatus(const char* inputfile,
    std::map<int,AliMuonCompactManuStatus>& manuStatusForRuns)
{
    std::ifstream in(inputfile,std::ios::binary);

    int nruns;

    in.read((char*)&nruns,sizeof(int));

    std::cout << "nruns=" << nruns << std::endl;

    std::vector<int> vrunlist;

    vrunlist.resize(nruns,0);

    std::array<UInt_t,16828> manuStatus;

    in.read((char*)&vrunlist[0],sizeof(int)*nruns);

    for ( auto i = 0; i < vrunlist.size(); ++i ) 
    {
        Int_t runNumber = vrunlist[i];
        std::cout << runNumber << " ";
        in.read((char*)&manuStatus[0],sizeof(UInt_t)*manuStatus.size());
        manuStatusForRuns[runNumber] = AliMuonCompactManuStatus(manuStatus);
    }
    std::cout << std::endl;
}

// Return a new manustatus object with the manus for the given buspatch removed from the configuration
AliMuonCompactManuStatus AliMuonCompactManuStatus::RemoveBusPatches(std::vector<int> busPatchIds) const
{
  AliMuonCompactMapping* cm = AliMuonCompactMapping::GetCompactMapping();

  auto n{*this};

  for (auto i = 0; i < mManuStatus.size(); ++i) {
    Int_t absManuId = cm->AbsManuId(i);
    Int_t detElemId = cm->GetDetElemIdFromAbsManuId(absManuId);
    Int_t manuId = cm->GetManuIdFromAbsManuId(absManuId);
    Int_t bp = AliMpDDLStore::Instance()->GetBusPatchId(detElemId, manuId);
    if (std::find(std::begin(busPatchIds),std::end(busPatchIds),bp) != busPatchIds.end()) {
      n.mManuStatus[i] |= MANUOUTOFCONFIGMASK;
    }
  }

  return n;
}

/// True if all manu status are at zero
bool AliMuonCompactManuStatus::AllClean() const {
  return std::count(std::begin(mManuStatus),std::end(mManuStatus),0)==mManuStatus.size();
}

int AliMuonCompactManuStatus::CountBad(UInt_t cause) const
{
  return std::count_if(mManuStatus.begin(), mManuStatus.end(), [&](int n) { return (n & cause); });
}