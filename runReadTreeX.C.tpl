#include "/afs/cern.ch/user/c/chiw/cernbox/JpsiJpsiUps/preCut/ParticleCand/src/ParticleCand.C"
#include "ReadTree.C"

// #define RUN_MINOR
#define RUN_ALL
// #define RUN_JOB

void runReadTreeX(){
    TChain *chain = new TChain("mkcands/X_data","");
    // chain->Add("/home/storage0/users/chiwang/storage2/CMS-Analysis/JJP/Analysis/MultiCandChecker/TripleMeson-ntuple.root");
    // chain->Add("/home/storage0/users/chiwang/storage2/CMS-Analysis/JJP/TripleMeson-ntuple.root");
    #ifdef RUN_JOB
    chain->Add("JOB_DATA");
    #else

    #ifdef RUN_MINOR
    //chain->Add("/home/storage0/users/chiwang/storage2/CMS-Analysis/JJP/Analysis/MultiCandChecker/JJU-ntuple.root");
    chain->Add("~/storage2/CMS-Analysis/JpsiJpsiUps/rootNtuple/P_Run2022/P0_Run2022Cv1_0000.root");
    //chain->Add("/afs/cern.ch/user/c/chiw/cernbox/JpsiJpsiUps/rootNtuple/P_Run2024/P0_Run2024Cv1_0000.root");
    #endif

    

    #ifdef RUN_ALL
    chain->Add("~/storage2/CMS-Analysis/JpsiJpsiUps/rootNtuple/P_Run2022/*.root");
    chain->Add("~/storage2/CMS-Analysis/JpsiJpsiUps/rootNtuple/P_Run2023/*.root");
    chain->Add("~/storage2/CMS-Analysis/JpsiJpsiUps/rootNtuple/P_Run2024/*.root");
    #endif
    #endif

    ReadTree myReadTree(chain);
    myReadTree.Loop();
}
