#include"MyClass.C"
#include<TChain.h>


void LoopPlots()
{
  TChain chain("Events");
  //chain.Add("MLNanoAODv9/Haa4e_mass20MeV_ctau010um_MLNanoAODv9_0.root");
  chain.Add("MLNanoAODv9/Haa4e_mass20MeV_ctau010um_MLNanoAODv9_10*.root");
  chain.Add("MLNanoAODv9/Haa4e_mass20MeV_ctau010um_MLNanoAODv9_20*.root");
  //chain.Add("MLNanoAODv9/Haa4e_mass20MeV_ctau010um_MLNanoAODv9_30*.root");
  //chain.Add("MLNanoAODv9/Haa4e_mass20MeV_ctau010um_MLNanoAODv9_40*.root");

  MyClass O(&chain);
  O.Loop(); 
  return;
}
