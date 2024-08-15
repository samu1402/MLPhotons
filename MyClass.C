#define MyClass_cxx
#include "MyClass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveStats.h>
#include <TMath.h>
#include <TGraph.h>
#include <THStack.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include "Math/Vector4D.h"
#include "TVector2.h"
#include <TEfficiency.h>

#include <iostream>

#define ETA_MAX 1.5 //Max eta for barrel
#define dEta 0.0174 * 1.5 // Crystal size, CFR sec. II article
#define dPhi 0.0174 * 1.5
#define dR2 dEta * dEta + dPhi * dPhi //match radius squared
#define M_MIN 5 // MeV, mass cut for a

using LorentzVector = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>;
using LorentzVectorE = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<double>>;

void HistoSave(TH1F *histo, TString filename, TString XLabel, TString YLabel = "");
void Histo2Save(TH2F *histo, TString filename, TString XLabel, TString YLabel);
void CompareHisto(TH1F *histo1, TH1F *histo2, Float_t x1 = 0.8, Float_t y1 = 0.8, Float_t x2 = 0.9, Float_t y2 = 0.9, TString lab1 = "", TString lab2 = "", TString Xlab = "", TString title = "", TString path = "");
void CompareHisto(TH1F *histo1, TH1F *histo2, TH1F *histo3, Float_t x1 = 0.8, Float_t y1 = 0.8, Float_t x2 = 0.9, Float_t y2 = 0.9, TString lab1 = "", TString lab2 = "", TString lab3 = "", TString Xlab = "", TString title = "", TString path = "");
void CompareHisto(TH1F *histo1, TH1F *histo2, TH1F *histo3, TH1F *histo4, Float_t x1, Float_t y1, Float_t x2, Float_t y2, TString lab1, TString lab2, TString lab3, TString lab4, TString Xlab, TString title, TString path);
void SaveHistoDiff(TH1F *hist1, TH1F *hist2, TString title, TString stats, TString path);
int isDaughter(Int_t *partPdgIds, Int_t *partMothIds, Int_t id);
void CompareHistoZoom(TH1F *histo1, TH1F *histo2, Float_t x1, Float_t y1, Float_t x2, Float_t y2, TString lab1, TString lab2, TString Xlab, TString title, TString path, Float_t xz_1, Float_t yz_1, Float_t xz_2, Float_t yz_2, Float_t xmin, Float_t xmax);
void EffSave(TEfficiency *histo, TString filename);


void MyClass::Loop()
{
  //   In a ROOT session, you can do:
  //      root> .L MyClass.C
  //      root> MyClass t
  //      root> t.GetEntry(12); // Fill t data members with entry number 12
  //      root> t.Show();       // Show values of entry 12
  //      root> t.Show(16);     // Read and show values of entry 16
  //      root> t.Loop();       // Loop on all entries
  //

  //     This is the loop skeleton where:
  //    jentry is the global entry number in the chain
  //    ientry is the entry number in the current Tree
  //  Note that the argument to GetEntry must be:
  //    jentry for TChain::GetEntry
  //    ientry for TTree::GetEntry and TBranch::GetEntry
  //
  //       To read only selected branches, Insert statements like:
  // METHOD1:
  //    fChain->SetBranchStatus("*",0);  // disable all branches
  //    fChain->SetBranchStatus("branchname",1);  // activate branchname
  // METHOD2: replace line
  //    fChain->GetEntry(jentry);       //read all branches
  // by  b_branchname->GetEntry(ientry); //read only this branch
  if (fChain == 0)
    return;
  Long64_t nentries = fChain->GetEntriesFast(); 
  Long64_t nbytes = 0, nb = 0;

  TDatabasePDG *pdg = TDatabasePDG::Instance();
  TParticlePDG *electron = pdg->GetParticle(11);
  Double_t e_mass = electron->Mass(); //Electron mass from pdg in GeV

  const float eta_min = -1.5; //Plot limits for different quantities
  const float eta_max = 1.5;

  const float pt_min = 0;   // GeV
  const float pt_max = 150; // GeV

  const float mee_min = 0;   // MeV Inv mass of ee
  const float mee_max = 120; // MeV

  const float m4e_min = 115; // GeV Inv mass of 4e
  const float m4e_max = 135; // GeV

  const int eta_bin = 100;  //Plot binning for different quantities
  const int phi_bin = 100;
  const int pt_bin = 100;
  const int mee_bin = 100;
  const int m4e_bin = 100;

  //const Double_t dRs[] = {0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035};
  //const Int_t n_dRs = 7;

  //const Int_t n_EleNeigh = 4;
  //const Int_t n_PhoNeigh = 4;

  // Gen histograms: kinematics
  TH1F *H_phi_GenEle = new TH1F("#phi e", "#phi e", phi_bin, -TMath::Pi(), TMath::Pi());
  TH1F *H_phi_GenPho = new TH1F("#phi #gamma", "#phi #gamma", phi_bin, -TMath::Pi(), TMath::Pi());

  TH1F *H_eta_GenEle = new TH1F("#eta e", "#eta e", eta_bin, eta_min, eta_max);
  TH1F *H_eta_GenPho = new TH1F("#eta #gamma", "#eta #gamma", eta_bin, eta_min, eta_max);

  TH1F *H_pt_GenEle = new TH1F("pt e", "pt# e", pt_bin, pt_min, pt_max);
  TH1F *H_pt_GenPho = new TH1F("pt #gamma", "pt #gamma", pt_bin, pt_min, pt_max);

  TH1F *H_s_ee = new TH1F("#sqrt(s) ee", "#sqrt(s) ee", mee_bin, mee_min, mee_max);
  TH1F *H_s_4e = new TH1F("#sqrt(s) 4e", "#sqrt(s) 4e", m4e_bin, m4e_min, m4e_max);

  TH2F *H_aa_dphideta = new TH2F("", "#Delta#phi Vs #Delta#eta aa", 100, -TMath::Pi(), TMath::Pi(), 100, -3, 3);
  TH2F *H_eem_dphideta = new TH2F("", "#Delta#phi Vs #Delta#eta ee low m", 100, -0.02, 0.02, 100, -0.02, 0.02);
  TH2F *H_eeM_dphideta = new TH2F("", "#Delta#phi Vs #Delta#eta ee high m", 100, -0.02, 0.02, 100, -0.02, 0.02);

  TH1F *H_EleMass_matched = new TH1F("Ele mass", "Ele mass", 100, -0.1, 0.1);

  // Reco histograms: kinematics
  TH1F *H_phi_Ele = new TH1F("reco #phi e", "reco #phi e", phi_bin, -TMath::Pi(), TMath::Pi());
  TH1F *H_phi_Pho = new TH1F("reco #phi #gamma", "reco #phi #gamma", phi_bin, -TMath::Pi(), TMath::Pi());

  TH1F *H_eta_Ele = new TH1F("reco #eta e", "reco #eta e", eta_bin, eta_min, eta_max);
  TH1F *H_eta_Pho = new TH1F("reco #eta #gamma", "reco #eta #gamma", eta_bin, eta_min, eta_max);

  TH1F *H_pt_Ele = new TH1F("reco pt e", "reco pt e", pt_bin, pt_min, pt_max);
  TH1F *H_pt_Pho = new TH1F("reco pt #gamma", "reco pt #gamma", pt_bin, pt_min, pt_max);

  TH1F *H_s_eeReco = new TH1F("#sqrt(s) ee reco", "#sqrt(s) ee reco", mee_bin, mee_min, mee_max);
  TH1F *H_s_4eReco = new TH1F("sqrt{s} 4e", "sqrt{s} 4e", m4e_bin, m4e_min, m4e_max);
  TH1F *H_s_eeReco_M = new TH1F("sqrt{s} ee M", "sqrt{s} ee M", mee_bin, mee_min, mee_max);
  TH1F *H_s_eeReco_E = new TH1F("sqrt{s} ee E", "sqrt{s} ee E", mee_bin, mee_min, mee_max);
  TH1F *H_s_4eReco_E = new TH1F("sqrt{s} 4e E", "sqrt{s} 4e E", m4e_bin, m4e_min, m4e_max);

  TH1F *H_s_ggReco = new TH1F("sqrt{s} gg reco", "#sqrt(s) gg reco", mee_bin, mee_min, mee_max);
  TH1F *H_s_4gReco = new TH1F("sqrt{s} 4g", "sqrt{s} 4g", m4e_bin, m4e_min, m4e_max);

  // GenReco histogram: kinematics
  TH2F *H_phiEle_GenReco = new TH2F("phi Scatter Ele", "phi Scatter Ele", phi_bin, -TMath::Pi(), TMath::Pi(), phi_bin, -TMath::Pi(), TMath::Pi());
  TH2F *H_phiPho_GenReco = new TH2F("phi Scatter Pho", "phi Scatter Pho", phi_bin, -TMath::Pi(), TMath::Pi(), phi_bin, -TMath::Pi(), TMath::Pi());

  TH2F *H_etaEle_GenReco = new TH2F("eta Scatter Ele", "eta Scatter Ele", eta_bin, eta_min, eta_max, eta_bin, eta_min, eta_max);
  TH2F *H_etaPho_GenReco = new TH2F("eta Scatter Pho", "eta Scatter Pho", eta_bin, eta_min, eta_max, eta_bin, eta_min, eta_max);

  TH2F *H_ptEle_GenReco = new TH2F("pt Scatter Ele", "pt Scatter Ele", pt_bin, pt_min, pt_max, pt_bin, pt_min, pt_max);
  TH2F *H_ptPho_GenReco = new TH2F("pt Scatter Pho", "pt Scatter Pho", pt_bin, pt_min, pt_max, pt_bin, pt_min, pt_max);

  // Gen histogram: misidentified diPho/diEle pt
  TH1F *H_pt_EleSum = new TH1F("Gen pt e summed", "Gen pt e sum", pt_bin, pt_min, pt_max);
  TH1F *H_pt_PhoSum = new TH1F("Gen pt #gamma summed", "Gen pt #gamma sum", pt_bin, pt_min, pt_max);

  // Gen histogram: # of Reco matches
  TH1F *H_GenMatchEle = new TH1F("Ele neighbors", "", 5, -0.5, 4.5);
  TH1F *H_GenMatchPho = new TH1F("Pho neighbors", "", 5, -0.5, 4.5);

  // GenReco histogram: dR of matches
  TH1F *H_dRMatchEle = new TH1F("Ele dR match", "Ele dR match", 50, 0, sqrt(dR2));
  TH1F *H_dRMatchPho = new TH1F("Pho dR match", "Pho dR match", 50, 0, sqrt(dR2));

  // ML histograms: kinematics
  TH1F *H_phi_GenMLPho = new TH1F("ML #phi #gamma", "ML #phi #gamma", phi_bin, -TMath::Pi(), TMath::Pi());
  TH1F *H_eta_GenMLPho = new TH1F("ML #eta #gamma", "ML #eta #gamma", eta_bin, eta_min, eta_max);
  TH1F *H_pt_GenMLPho = new TH1F("ML pt #gamma", "ML pt #gamma", pt_bin, pt_min, pt_max);
  TH1F *H_phi_MLPho = new TH1F("ML reco #phi #gamma", "ML reco #phi #gamma", phi_bin, -TMath::Pi(), TMath::Pi());
  TH1F *H_eta_MLPho = new TH1F("ML reco #eta #gamma", "ML reco #eta #gamma", eta_bin, eta_min, eta_max);
  TH1F *H_pt_MLPho = new TH1F("ML reco pt #gamma", "ML reco pt #gamma", pt_bin, pt_min, pt_max);
  TH1F *H_pt_MLPhoSum = new TH1F("ML Gen pt #gamma summed", "ML Gen pt #gamma sum", pt_bin, pt_min, pt_max);

  // ML histograms: GenReco kinematics, matches
  TH1F *H_GenMatchMLPho = new TH1F("MLPho neighbors", "", 5, -0.5, 4.5);
  TH2F *H_phiMLPho_GenReco = new TH2F("#phi Scatter MLPho", "#phi Scatter MLPho", phi_bin, -TMath::Pi(), TMath::Pi(), phi_bin, -TMath::Pi(), TMath::Pi());
  TH2F *H_etaMLPho_GenReco = new TH2F("#eta Scatter MLPho", "#eta Scatter MLPho", eta_bin, eta_min, eta_max, eta_bin, eta_min, eta_max);
  TH2F *H_ptMLPho_GenReco = new TH2F("pt Scatter MLPho", "pt Scatter MLPho", pt_bin, pt_min, pt_max, pt_bin, pt_min, pt_max);
  TH1F *H_dRMatchMLPho = new TH1F("MLPho dR match", "MLPho dR match", 50, 0, sqrt(dR2));

  // GenReco histograms: # of matches per reco part
  TH1F *H_EleMatches = new TH1F("Ele Matches", "Ele Matches", 4, -0.5, 3.5);
  TH1F *H_PhoMatches = new TH1F("Pho Matches", "Pho Matches", 4, -0.5, 3.5);

  // Gen histograms: dR of diEle
  TH1F *H_GenEledR = new TH1F("Gen Ele dR", "Gen Ele dR", 100, 0, sqrt(dR2));

  // Histogram for ML Vs nonML Pho
  TH2F *H_etaPho_RecoVsML = new TH2F("#eta pho ML Vs non ML", "#eta pho ML Vs non ML", eta_bin, eta_min, eta_max, eta_bin, eta_min, eta_max);

  // GenReco histogram: invariant mass
  TH2F *H_s_eeGenReco = new TH2F("#sqrt(s) ee Gen Vs Reco", "#sqrt(s) ee Gen Vs Reco", mee_bin, mee_min, mee_max, mee_bin, mee_min, mee_max);

  // Debug histograms
  TH2F *H = new TH2F("", "", mee_bin, mee_min, mee_max, mee_bin, mee_min, mee_max);
  TH1F *H_phi5_1 = new TH1F("", "", phi_bin, -TMath::Pi(), TMath::Pi());
  TH1F *H_eta5_1 = new TH1F("", "", eta_bin, eta_min, eta_max);
  TH1F *H_pt5_1 = new TH1F("", "", pt_bin, pt_min, pt_max);
  TH1F *H_phi5_2 = new TH1F("", "", phi_bin, -TMath::Pi(), TMath::Pi());
  TH1F *H_eta5_2 = new TH1F("", "", eta_bin, eta_min, eta_max);
  TH1F *H_pt5_2 = new TH1F("", "", pt_bin, pt_min, pt_max);
  TH1F *H_phi5_3 = new TH1F("", "", phi_bin, -TMath::Pi(), TMath::Pi());
  TH1F *H_eta5_3 = new TH1F("", "", eta_bin, eta_min, eta_max);
  TH1F *H_pt5_3 = new TH1F("", "", pt_bin, pt_min, pt_max);
  TH1F *H_phi5_4 = new TH1F("", "", phi_bin, -TMath::Pi(), TMath::Pi());
  TH1F *H_eta5_4 = new TH1F("", "", eta_bin, eta_min, eta_max);
  TH1F *H_pt5_4 = new TH1F("", "", pt_bin, pt_min, pt_max);

  // ML debug histograms
  TH1F *H_MLPho_dphi = new TH1F("(#phi Pho - MLPho) / #phi Pho", "(#phi Pho - MLPho) / #phi Pho", 50, -1, 1);
  TH1F *H_MLPho_deta = new TH1F("(#eta Pho - MLPho) / #eta Ph", "(#eta Pho - MLPho) / #eta Pho", 50, -1, 1);
  TH1F *H_MLPho_dpt = new TH1F("(p_{t} Pho - MLPho) / p_{t} Pho", "(p_{t} Pho - MLPho) / p_{t} Pho", 50, -1, 1);

  // Histograms: counting neighbours
  //TH1F **H_EleNeigh = new TH1F *[n_dRs];
  //TH1F **H_PhoNeigh = new TH1F *[n_dRs];

  // Histogram for efficiency
  TEfficiency *H_EleEff_pt = new TEfficiency("Ele reco efficiency", ";p_{t};Eff", pt_bin, pt_min, pt_max);
  TEfficiency *H_PhoEff_pt = new TEfficiency("Pho reco efficiency", ";p_{t};Eff", pt_bin, pt_min, pt_max);
  TEfficiency *H_MLPhoEff_pt = new TEfficiency("MLPho reco efficiency", ";p_{t};Eff", pt_bin, pt_min, pt_max);
  TEfficiency *H_EleEff_eta = new TEfficiency("Ele reco efficiency", ";#eta;Eff", eta_bin, eta_min, eta_max);
  TEfficiency *H_PhoEff_eta = new TEfficiency("Pho reco efficiency", ";#eta;Eff", eta_bin, eta_min, eta_max);
  TEfficiency *H_MLPhoEff_eta = new TEfficiency("MLPho reco efficiency", ";#eta;Eff", eta_bin, eta_min, eta_max);

  TEfficiency *H_MLPho2Eff_pt = new TEfficiency("MLPho reco 2 pho efficiency", ";p{t};Eff", pt_bin, pt_min, pt_max);
  TEfficiency *H_MLPho2Eff_eta = new TEfficiency("MLPho reco 2 pho efficiency", ";#eta;Eff", eta_bin, eta_min, eta_max);
  
  /*
  for (UInt_t i = 0; i < n_dRs; i++)
  {
    H_EleNeigh[i] = new TH1F(Form("Ele Neighbors %.3f", dRs[i]), Form("Ele Neighbors %.3f", dRs[i]), n_EleNeigh + 1, -0.5, 0.5 + n_EleNeigh);
    H_PhoNeigh[i] = new TH1F(Form("Pho Neighbors %.3f", dRs[i]), Form("Pho Neighbors %.3f", dRs[i]), n_PhoNeigh + 1, -0.5, 0.5 + n_PhoNeigh);
  }
  */



  for (Long64_t jentry = 0; jentry < nentries; jentry++)
  {
    if (jentry % 1000 == 0)
      std::cout << jentry << std::endl;

    // if(jentry>1) break; //Uncomment to test on 1 event
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0)
      break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;

    // Arrays for electron match quantities
    Int_t *e_matchId = new Int_t[nElectron];
    Double_t *e_matchdR2 = new Double_t[nElectron];
    Double_t *e_ptGenSum = new Double_t[nElectron];
    Bool_t *e_match = new Bool_t[nElectron];
    Int_t *EleMatches = new Int_t[nElectron];
    //Double_t *EleNeigh = new Double_t[n_dRs];

    std::fill(e_matchId, e_matchId + nElectron, -10);
    std::fill(e_matchdR2, e_matchdR2 + nElectron, 100.);
    std::fill(e_ptGenSum, e_ptGenSum + nElectron, 0);
    std::fill(e_match, e_match + nElectron, 0);
    std::fill(EleMatches, EleMatches + nElectron, 0);
    //td::fill(EleNeigh, EleNeigh + n_dRs, -1);

    // Arrays for photon match quantities
    Int_t *ph_matchId = new Int_t[nPhoton];
    Double_t *ph_matchdR2 = new Double_t[nPhoton];
    Double_t *ph_ptGenSum = new Double_t[nPhoton];
    Bool_t *ph_match = new Bool_t[nPhoton];
    Int_t *PhoMatches = new Int_t[nPhoton];

    std::fill(ph_matchId, ph_matchId + nPhoton, -10);
    std::fill(ph_matchdR2, ph_matchdR2 + nPhoton, 100.);
    std::fill(ph_ptGenSum, ph_ptGenSum + nPhoton, 0);
    std::fill(ph_match, ph_match + nPhoton, 0);
    std::fill(PhoMatches, PhoMatches + nPhoton, 0);

    // Arrays for MLphoton match quantities
    Int_t *MLph_matchId = new Int_t[nMLPhoton];
    Double_t *MLph_matchdR2 = new Double_t[nMLPhoton];
    Double_t *MLph_ptGenSum = new Double_t[nMLPhoton];

    Bool_t *MLph_match = new Bool_t[nMLPhoton];
    std::fill(MLph_matchId, MLph_matchId + nMLPhoton, -10);
    std::fill(MLph_matchdR2, MLph_matchdR2 + nMLPhoton, 100.);
    std::fill(MLph_ptGenSum, MLph_ptGenSum + nMLPhoton, 0);
    std::fill(MLph_match, MLph_match + nMLPhoton, 0);

    //Lorentz vector with saved quantities for Eles
    LorentzVector p_GenEle1, p_GenEle2, p_GenEle3, p_GenEle4;
    LorentzVector p_RecEle1, p_RecEle2, p_RecPho1, p_RecPho2;
    LorentzVectorE p_RecEle1_E, p_RecEle2_E;
    //Did a match with an a (5_1) or an anti-a (5_1) occur?
    Bool_t p_GenEle5_1 = 0;
    Bool_t p_GenEle5__1 = 0;

    Bool_t p_RecEle5_1 = 0;
    Bool_t p_RecEle5__1 = 0;

    Bool_t p_RecPho5_1 = 0;
    Bool_t p_RecPho5__1 = 0;
    
    //Masses of the reco eles daughters of anti/a
    Double_t M1 = -1;
    Double_t M2 = -1;
    //Does the Reco mass pass the mass cut?
    Bool_t M1_grt_M = 1;
    Bool_t M2_grt_M = 1;

    for (UInt_t i = 0; i < nGenPart; i++) 
    {
      if (fabs(GenPart_eta[i]) > ETA_MAX) //In barrel
        continue;
      // loop to identify aa couples and find aa displacement
      if (GenPart_pdgId[i] == 5000001)
        for (UInt_t j = 0; j < nGenPart; j++)
          if (GenPart_pdgId[j] == -5000001 && fabs(GenPart_eta[j]) <= ETA_MAX)
            H_aa_dphideta->Fill(TVector2::Phi_mpi_pi(GenPart_phi[i] - GenPart_phi[j]), GenPart_eta[i] - GenPart_eta[j]);
      //No need to look at particles coming directly from collisions other than +-5000001
      if (GenPart_genPartIdxMother[i] == -1)
        continue;
      //We're only interested in electrons and photons (There are only electrons...)
      if (fabs(GenPart_pdgId[i]) != 11 && GenPart_pdgId[i] != 22)
        continue;
      //we only look at daughters of anti/a
      if (abs(GenPart_pdgId[GenPart_genPartIdxMother[i]]) == 5000001)
      {
        // loop on gen parts to reco true inv mass
        if (!p_GenEle5_1 || !p_GenEle5__1) //if at least one gen ee hasn't been found, keep looking
          for (UInt_t j = i + 1; j < nGenPart; j++) //no need to re check checked particles
          {
            if ((fabs(GenPart_eta[j]) > ETA_MAX) || (GenPart_pdgId[GenPart_genPartIdxMother[j]] != GenPart_pdgId[GenPart_genPartIdxMother[i]]) || fabs(GenPart_pdgId[j]) != 11)
              continue; //Same conditions as b4 + particles must have the same mother (anti/a)
            if (GenPart_pdgId[GenPart_genPartIdxMother[j]] == 5000001 && (!p_GenEle5_1))
            {//If this is a match, save 4 momenta etc
              p_GenEle1.SetCoordinates(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], e_mass);
              p_GenEle2.SetCoordinates(GenPart_pt[j], GenPart_eta[j], GenPart_phi[j], e_mass);

              if ((p_GenEle1 + p_GenEle2).M() * 1000 < M_MIN) //If this event doesn't pass the cut, don't save stuff
              {
                M1_grt_M = 0;
                break;
              }
              p_GenEle5_1 = 1; //We found a Gen couple of ee!
              //dR of gen diEle
              H_GenEledR->Fill(sqrt((GenPart_eta[i] - GenPart_eta[j]) * (GenPart_eta[i] - GenPart_eta[j]) + (GenPart_phi[i] - GenPart_phi[j]) * (GenPart_phi[i] - GenPart_phi[j])));
  
              break;//no point in looking for another daughter e
            }
            else if (GenPart_pdgId[GenPart_genPartIdxMother[j]] == -5000001 && (!p_GenEle5__1))
            {
              p_GenEle3.SetCoordinates(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], e_mass);
              p_GenEle4.SetCoordinates(GenPart_pt[j], GenPart_eta[j], GenPart_phi[j], e_mass);

              if ((p_GenEle3 + p_GenEle4).M() * 1000 < M_MIN)
              {
                M2_grt_M = 0;
                break;
              }
              p_GenEle5__1 = 1;

              H_GenEledR->Fill(sqrt((GenPart_eta[i] - GenPart_eta[j]) * (GenPart_eta[i] - GenPart_eta[j]) + (GenPart_phi[i] - GenPart_phi[j]) * (GenPart_phi[i] - GenPart_phi[j])));
              break;
            }
          }
        //If this gen part doesn't pass the cut, don't even proceed with the analysis
        if (GenPart_pdgId[GenPart_genPartIdxMother[i]] == 5000001 && !M1_grt_M)
          continue;
        if (GenPart_pdgId[GenPart_genPartIdxMother[i]] == -5000001 && !M2_grt_M)
          continue;

        Double_t eta = GenPart_eta[i];
        Double_t phi = GenPart_phi[i];

        //std::fill(EleNeigh, EleNeigh + n_dRs, 0); 

        Bool_t GenPart_EleMatch = 0; // Does the GenPart match with a Reco Ele?

        for (UInt_t j = 0; j < nElectron; j++)//loop for reco Ele
        {
          Double_t e_eta = Electron_eta[j];
          Double_t e_phi = Electron_phi[j];
          Double_t e_dR2 = (e_eta - eta) * (e_eta - eta) + (e_phi - phi) * (e_phi - phi);

          if (fabs(e_eta) < ETA_MAX && e_dR2 < dR2) //if matched
          {
            e_match[j] = 1;//this reco Ele has a gen match
            EleMatches[j]++;
            GenPart_EleMatch = 1;//this gen part has a reco match

            H_EleMass_matched->Fill(Electron_mass[j]); 

            if (e_dR2 < e_matchdR2[j])//find best match
            {
              e_matchdR2[j] = e_dR2;
              e_matchId[j] = i;

              if (GenPart_pdgId[GenPart_genPartIdxMother[i]] == 5000001)
              {//coords of reco eles
                p_RecEle1.SetCoordinates(Electron_pt[j], Electron_eta[j], Electron_phi[j], fabs(Electron_mass[j]));
                p_RecEle1_E.SetCoordinates(Electron_pt[j], Electron_eta[j], Electron_phi[j], 1 / (Electron_eInvMinusPInv[j] + 1 / (Electron_pt[j] * cosh(Electron_eta[j]))));
                p_RecEle5_1 = 1;
                M1 = Electron_mass[j];
              }
              if (GenPart_pdgId[GenPart_genPartIdxMother[i]] == -5000001)
              {
                p_RecEle2.SetCoordinates(Electron_pt[j], Electron_eta[j], Electron_phi[j], fabs(Electron_mass[j]));
                p_RecEle2_E.SetCoordinates(Electron_pt[j], Electron_eta[j], Electron_phi[j], 1 / (Electron_eInvMinusPInv[j] + 1 / (Electron_pt[j] * cosh(Electron_eta[j]))));
                p_RecEle5__1 = 1;
                M2 = Electron_mass[j];
              }
            }

            // Fill Reco if match
            H_phi_Ele->Fill(Electron_phi[j]);
            H_eta_Ele->Fill(Electron_eta[j]);
            H_pt_Ele->Fill(Electron_pt[j]);
            // Fill Gen if match
            H_phi_GenEle->Fill(GenPart_phi[i]);
            H_eta_GenEle->Fill(GenPart_eta[i]);
            H_pt_GenEle->Fill(GenPart_pt[i]);

            // Fill scatter if match
            H_phiEle_GenReco->Fill(GenPart_phi[i], Electron_phi[j]);
            H_etaEle_GenReco->Fill(GenPart_eta[i], Electron_eta[j]);

            e_ptGenSum[j] += GenPart_pt[i];
          }
          /*
          for (UInt_t k = 0; k < n_dRs; k++)
            if (e_dR2 < dRs[k])
              EleNeigh[k]++; // # ele in diff dR
              */
        }
        /*for (UInt_t j = 0; j < n_dRs; j++)
          H_EleNeigh[j]->Fill(EleNeigh[j]);*/

        H_EleEff_pt->Fill(GenPart_EleMatch, GenPart_pt[i]);
        H_EleEff_eta->Fill(GenPart_EleMatch, GenPart_eta[i]);

        Int_t PhoMatchRecoId = -10;
        Bool_t GenPart_PhoMatch = 0;

        for (UInt_t j = 0; j < nPhoton; j++)//Same stuff as b4
        {
          Double_t ph_eta = Photon_eta[j];
          Double_t ph_phi = Photon_phi[j];
          Double_t ph_dR2 = (ph_eta - eta) * (ph_eta - eta) + (ph_phi - phi) * (ph_phi - phi);
          if (fabs(ph_eta) < ETA_MAX && ph_dR2 < dR2)
          {
            ph_match[j] = 1;
            PhoMatches[j]++;
            GenPart_PhoMatch = 1;

            if (ph_dR2 < ph_matchdR2[j])
            {
              ph_matchdR2[j] = ph_dR2;
              ph_matchId[j] = i;
              PhoMatchRecoId = j;

              if (GenPart_pdgId[GenPart_genPartIdxMother[i]] == 5000001)
              {
                p_RecPho1.SetCoordinates(Photon_pt[j], Photon_eta[j], Photon_phi[j], fabs(Photon_mass[j]));
                p_RecPho5_1 = 1;
              }
              if (GenPart_pdgId[GenPart_genPartIdxMother[i]] == -5000001)
              {
                p_RecPho2.SetCoordinates(Photon_pt[j], Photon_eta[j], Photon_phi[j], fabs(Photon_mass[j]));
                p_RecPho5__1 = 1;
              }
            }
            // FIll Reco if match
            H_phi_Pho->Fill(Photon_phi[j]);
            H_eta_Pho->Fill(Photon_eta[j]);
            H_pt_Pho->Fill(Photon_pt[j]);
            // Fill Gen if match
            H_phi_GenPho->Fill(GenPart_phi[i]);
            H_eta_GenPho->Fill(GenPart_eta[i]);
            H_pt_GenPho->Fill(GenPart_pt[i]);

            // Fill scatter if match
            H_phiPho_GenReco->Fill(GenPart_phi[i], Photon_phi[j]);
            H_etaPho_GenReco->Fill(GenPart_eta[i], Photon_eta[j]);

            ph_ptGenSum[j] += GenPart_pt[i];
          }
        }
        H_PhoEff_pt->Fill(GenPart_PhoMatch, GenPart_pt[i]);
        H_PhoEff_eta->Fill(GenPart_PhoMatch, GenPart_eta[i]);

        // if(PhoMatchRecoId == -10) continue;//MLPhoton analysis only if match

        Int_t MLPhoMatchRecoId = -10;
        Bool_t GenPart_MLPhoMatch = 0;

        for (UInt_t j = 0; j < nMLPhoton; j++)
        {

          Double_t MLph_eta = MLPhoton_eta[j];
          Double_t MLph_phi = MLPhoton_phi[j];
          Double_t MLph_dR2 = (MLph_eta - eta) * (MLph_eta - eta) + (MLph_phi - phi) * (MLph_phi - phi);
          if (fabs(MLph_eta) < ETA_MAX && MLph_dR2 < dR2)
          {
            MLph_match[j] = 1;
            GenPart_MLPhoMatch = 1;
            if (MLph_dR2 < MLph_matchdR2[j])
            {
              MLph_matchdR2[j] = MLph_dR2;
              MLph_matchId[j] = i;
              MLPhoMatchRecoId = j;
            }
            // FIll Reco if match
            H_phi_MLPho->Fill(MLPhoton_phi[j]);
            H_eta_MLPho->Fill(MLPhoton_eta[j]);
            H_pt_MLPho->Fill(MLPhoton_pt[j]);
            // Fill Gen if match
            H_phi_GenMLPho->Fill(GenPart_phi[i]);
            H_eta_GenMLPho->Fill(GenPart_eta[i]);
            H_pt_GenMLPho->Fill(GenPart_pt[i]);

            // Fill scatter if match
            H_phiMLPho_GenReco->Fill(GenPart_phi[i], MLPhoton_phi[j]);
            H_etaMLPho_GenReco->Fill(GenPart_eta[i], MLPhoton_eta[j]);

            MLph_ptGenSum[j] += GenPart_pt[i];

            if (PhoMatchRecoId != -10)
            {
              H_etaPho_RecoVsML->Fill(Photon_eta[PhoMatchRecoId], MLPhoton_eta[j]);
              H_MLPho_dphi->Fill((Photon_phi[PhoMatchRecoId] - MLPhoton_phi[j]) / fabs(Photon_phi[PhoMatchRecoId]));
              H_MLPho_deta->Fill((Photon_eta[PhoMatchRecoId] - MLPhoton_eta[j]) / fabs(Photon_eta[PhoMatchRecoId]));
              H_MLPho_dpt->Fill((Photon_pt[PhoMatchRecoId] - MLPhoton_pt[j]) / fabs(Photon_pt[PhoMatchRecoId]));

              if (fabs(Photon_eta[PhoMatchRecoId] - MLPhoton_eta[j]) > 1.4)
              {
                Int_t j1 = ph_matchId[PhoMatchRecoId];
                Int_t j2 = MLph_matchId[MLPhoMatchRecoId];
                std::cout << jentry << std::endl;
                printf("%d\t%d\t%d\t%f\t%f\t%f\n", i, j1, j2, GenPart_eta[j1], Photon_eta[PhoMatchRecoId], MLPhoton_eta[j]);
              }
            }
          }
        }
        H_MLPhoEff_pt->Fill(GenPart_MLPhoMatch, GenPart_pt[i]);
        H_MLPhoEff_eta->Fill(GenPart_MLPhoMatch, GenPart_eta[i]);
        
      } // In barrel + isDaughter + stable if
    } // GenPart for


    if (p_GenEle5_1 && p_GenEle5__1) //If both aa are GenFound + pass cuts
    {
      H_s_4e->Fill((p_GenEle1 + p_GenEle2 + p_GenEle3 + p_GenEle4).M());
      H_s_ee->Fill((p_GenEle1 + p_GenEle2).M() * 1000);
      H_s_ee->Fill((p_GenEle3 + p_GenEle4).M() * 1000);
      H->Fill((p_GenEle1 + p_GenEle2).M() * 1000, (p_GenEle3 + p_GenEle4).M() * 1000);

      if ((p_GenEle1 + p_GenEle2).M() < (p_GenEle3 + p_GenEle4).M())//Print lowest/highest mass kinematic variables
      {
        H_pt5_1->Fill(p_GenEle1.pt());
        H_eta5_1->Fill(p_GenEle1.eta());
        H_phi5_1->Fill(p_GenEle1.phi());

        H_pt5_2->Fill(p_GenEle2.pt());
        H_eta5_2->Fill(p_GenEle2.eta());
        H_phi5_2->Fill(p_GenEle2.phi());

        H_pt5_3->Fill(p_GenEle3.pt());
        H_eta5_3->Fill(p_GenEle3.eta());
        H_phi5_3->Fill(p_GenEle3.phi());

        H_pt5_4->Fill(p_GenEle4.pt());
        H_eta5_4->Fill(p_GenEle4.eta());
        H_phi5_4->Fill(p_GenEle4.phi());

        H_eem_dphideta->Fill(TVector2::Phi_mpi_pi(p_GenEle1.phi() - p_GenEle2.phi()), p_GenEle1.eta() - p_GenEle2.eta());
        H_eeM_dphideta->Fill(TVector2::Phi_mpi_pi(p_GenEle3.phi() - p_GenEle4.phi()), p_GenEle3.eta() - p_GenEle4.eta());
      }
      else
      {
        H_pt5_1->Fill(p_GenEle3.pt());
        H_eta5_1->Fill(p_GenEle3.eta());
        H_phi5_1->Fill(p_GenEle3.phi());

        H_pt5_2->Fill(p_GenEle4.pt());
        H_eta5_2->Fill(p_GenEle4.eta());
        H_phi5_2->Fill(p_GenEle4.phi());

        H_pt5_3->Fill(p_GenEle1.pt());
        H_eta5_3->Fill(p_GenEle1.eta());
        H_phi5_3->Fill(p_GenEle1.phi());

        H_pt5_4->Fill(p_GenEle2.pt());
        H_eta5_4->Fill(p_GenEle2.eta());
        H_phi5_4->Fill(p_GenEle2.phi());

        H_eeM_dphideta->Fill(TVector2::Phi_mpi_pi(p_GenEle1.phi() - p_GenEle2.phi()), (p_GenEle1.eta() - p_GenEle2.eta()));
        H_eem_dphideta->Fill(TVector2::Phi_mpi_pi(p_GenEle3.phi() - p_GenEle4.phi()), (p_GenEle3.eta() - p_GenEle4.eta()));
      }
    }

    if (p_RecEle5_1 && p_RecEle5__1) //If both aa are Reco matched
    {
      H_s_eeReco->Fill(p_RecEle1.M() * 1000);
      H_s_eeReco->Fill(p_RecEle2.M() * 1000);
      H_s_4eReco->Fill((p_RecEle1 + p_RecEle2).M());
      H_s_eeReco_M->Fill(M1 * 1000);
      H_s_eeReco_M->Fill(M2 * 1000);
      H_s_eeReco_E->Fill(p_RecEle1_E.M() * 1000);
      H_s_eeReco_E->Fill(p_RecEle2_E.M() * 1000);
      H_s_4eReco_E->Fill((p_RecEle1_E + p_RecEle2_E).M());
      // std::cout<<p_RecEle1_E.M()<<std::endl;
    }

    if (p_RecPho5_1 && p_RecPho5__1)
    {
      H_s_ggReco->Fill(p_RecPho1.M() * 1000);
      H_s_ggReco->Fill(p_RecPho2.M() * 1000);
      H_s_4gReco->Fill((p_RecPho1 + p_RecPho2).M());
    }

    if (p_GenEle5_1 && p_GenEle5__1 && p_RecEle5_1 && p_RecEle5__1)
    {
      H_s_eeGenReco->Fill((p_GenEle1 + p_GenEle2).M() * 1000, p_RecEle1.M() * 1000);
      H_s_eeGenReco->Fill((p_GenEle3 + p_GenEle4).M() * 1000, p_RecEle2.M() * 1000);
    }
    for (UInt_t i = 0; i < nElectron; i++)
    {
      if (!e_match[i])
        continue;
      H_ptEle_GenReco->Fill(e_ptGenSum[i], Electron_pt[i]);
      H_pt_EleSum->Fill((Float_t)e_ptGenSum[i], 2);
      H_EleMatches->Fill(EleMatches[i]);
      H_dRMatchEle->Fill(sqrt(e_matchdR2[i]));
    }

    for (UInt_t i = 0; i < nPhoton; i++)
    {
      if (!ph_match[i])
        continue;
      H_ptPho_GenReco->Fill(ph_ptGenSum[i], Photon_pt[i]);
      H_pt_PhoSum->Fill((Float_t)ph_ptGenSum[i], 2);
      H_PhoMatches->Fill(PhoMatches[i]);
      H_dRMatchPho->Fill(sqrt(ph_matchdR2[i]));
    }

    for (UInt_t i = 0; i < nMLPhoton; i++)
    {
      if (!MLph_match[i])
        continue;
      H_ptMLPho_GenReco->Fill(MLph_ptGenSum[i], MLPhoton_pt[i]);
      H_pt_MLPhoSum->Fill((Float_t)MLph_ptGenSum[i], 2);
    }

    delete[] e_matchId;
    delete[] e_matchdR2;
    delete[] e_ptGenSum;
    delete[] e_match;
    delete[] EleMatches;
    //delete[] EleNeigh;

    delete[] ph_matchId;
    delete[] ph_matchdR2;
    delete[] ph_ptGenSum;
    delete[] ph_match;
    delete[] PhoMatches;

    delete[] MLph_matchId;
    delete[] MLph_matchdR2;
    delete[] MLph_ptGenSum;
  } // Ext For
  //Save stuff
  // Plot Gen quantities for Ele/Pho/aa
  HistoSave(H_phi_GenEle, "_.pdf", "#phi");
  HistoSave(H_phi_GenPho, "_.pdf", "#phi");

  HistoSave(H_eta_GenEle, "_.pdf", "#eta");
  HistoSave(H_eta_GenPho, "_.pdf", "#eta");

  HistoSave(H_pt_GenEle, "_.pdf", "pt [GeV]");
  HistoSave(H_pt_GenPho, "_.pdf", "pt [GeV]");

  Histo2Save(H_aa_dphideta, "Images/aa/dPhidEta.pdf", "#Delta#phi", "#Delta#eta");

  // Plot Reco quantities for Ele/Pho
  HistoSave(H_phi_Ele, "_.pdf", "", "#phi");
  HistoSave(H_phi_Pho, "_.pdf", "", "#phi");

  HistoSave(H_eta_Ele, "_.pdf", "", "#eta");
  HistoSave(H_eta_Pho, "_.pdf", "", "#eta");

  HistoSave(H_pt_Ele, "_.pdf", "", "pt");
  HistoSave(H_pt_Pho, "_.pdf", "", "pt");

  HistoSave(H_pt_EleSum, "_.pdf", "", "pt");
  HistoSave(H_pt_PhoSum, "_.pdf", "", "pt");

  // Plot ML quantities
  HistoSave(H_phi_MLPho, "_.pdf", "", "#phi");
  HistoSave(H_eta_MLPho, "_.pdf", "", "#eta");
  HistoSave(H_pt_MLPho, "_.pdf", "", "pt");
  HistoSave(H_pt_MLPhoSum, "_.pdf", "", "pt");

  // Plot comparison Gen/Reco
  CompareHisto(H_eta_Ele, H_eta_GenEle, 0.1, 0.75, 0.2, 0.85, "Reco", "Gen", "#eta", "", "Images/Ele/EtaEleComparison.pdf");
  CompareHisto(H_phi_Ele, H_phi_GenEle, 0.1, 0.75, 0.2, 0.85, "Reco", "Gen", "#phi", "", "Images/Ele/PhiEleComparison.pdf");
  CompareHisto(H_eta_Pho, H_eta_GenPho, 0.1, 0.75, 0.2, 0.85, "Reco", "Gen", "#eta", "", "Images/Pho/EtaPhoComparison.pdf");
  CompareHisto(H_phi_Pho, H_phi_GenPho, 0.1, 0.75, 0.2, 0.85, "Reco", "Gen", "#phi", "", "Images/Pho/PhiPhoComparison.pdf");
  CompareHisto(H_pt_Pho, H_pt_GenPho, H_pt_PhoSum, 0.1, 0.75, 0.3, 0.85, "Reco", "Gen", "Gen pt1+pt2", "pt [GeV]", "", "Images/Pho/PtPhoComparison.pdf");
  CompareHisto(H_pt_Ele, H_pt_GenEle, H_pt_EleSum, 0.1, 0.75, 0.3, 0.85, "Reco", "Gen", "Gen pt1+pt2", "pt [GeV]", "", "Images/Ele/PtEleComparison.pdf");

  // Plot scatters
  Histo2Save(H_phiEle_GenReco, "Images/Ele/phiGenReco.pdf", "Gen #phi", "Reco #phi");
  Histo2Save(H_phiPho_GenReco, "Images/Pho/phiGenReco.pdf", "Gen #phi", "Reco #phi");
  Histo2Save(H_phiMLPho_GenReco, "Images/MLPho/phiGenReco.pdf", "Gen #phi", "Reco #phi");

  Histo2Save(H_etaEle_GenReco, "Images/Ele/etaGenReco.pdf", "Gen #eta", "Reco #eta");
  Histo2Save(H_etaPho_GenReco, "Images/Pho/etaGenReco.pdf", "Gen #eta", "Reco #eta");
  Histo2Save(H_etaMLPho_GenReco, "Images/MLPho/etaGenReco.pdf", "Gen #eta", "Reco #eta");

  Histo2Save(H_ptEle_GenReco, "Images/Ele/ptGenReco.pdf", "Gen pt 1+2 [GeV]", "Reco pt [GeV]");
  Histo2Save(H_ptPho_GenReco, "Images/Pho/ptGenReco.pdf", "Gen pt 1+2 [GeV]", "Reco pt [GeV]");
  Histo2Save(H_ptMLPho_GenReco, "Images/MLPho/ptGenReco.pdf", "Gen pt 1+2 [GeV]", "Reco pt [GeV]");

  HistoSave(H_dRMatchEle, "Images/Ele/dRmatched.pdf", "dR", "Arbitrary units");
  HistoSave(H_dRMatchPho, "Images/Pho/dRmatched.pdf", "dR", "Arbitrary units");

  // Plot ML comparison Gen/Reco
  CompareHisto(H_eta_MLPho, H_eta_GenMLPho, 0.1, 0.75, 0.2, 0.85, "ML Reco", "Gen", "#eta", "", "Images/MLPho/EtaPhoComparison.pdf");
  CompareHisto(H_phi_MLPho, H_phi_GenMLPho, 0.1, 0.75, 0.2, 0.85, "ML Reco", "Gen", "#phi", "", "Images/MLPho/PhiPhoComparison.pdf");
  CompareHisto(H_pt_MLPho, H_pt_GenMLPho, H_pt_MLPhoSum, 0.1, 0.75, 0.3, 0.85, "ML Reco", "Gen", "Gen pt1+pt2", "pt [GeV]", "", "Images/MLPho/PtPhoComparison.pdf");

  // Plot # of matches per reco part
  HistoSave(H_EleMatches, "Images/Ele/EleMatches.pdf", "# Matches", "Arbitrary units");
  HistoSave(H_PhoMatches, "Images/Pho/PhoMatches.pdf", "# Matches", "Arbitrary units");

  // Plot invariant mass
  HistoSave(H_s_ee, "_.pdf", "#sqrt(s) [MeV]", "Arbitrary units");
  HistoSave(H_s_4e, "_.pdf", "#sqrt(s) [GeV]", "Arbitrary units");

  CompareHisto(H_s_eeReco, H_s_ee, H_s_eeReco_M, 0.1, 0.75, 0.2, 0.85, "Reco", "Gen", "Reco with M", "#sqrt(s) [MeV]", "ee #sqrt(s)", "Images/InvMass/InvMass_ee_M.pdf");
  CompareHisto(H_s_eeReco, H_s_ee, H_s_ggReco, 0.1, 0.75, 0.2, 0.85, "Reco ee", "Gen", "Reco #gamma#gamma", "#sqrt(s) [MeV]", "#sqrt(s)", "Images/InvMass/InvMass_eegg.pdf");

  CompareHisto(H_s_4eReco, H_s_4e, H_s_4gReco, 0.1, 0.75, 0.2, 0.85, "Reco", "Gen", "Reco #gamma#gamma", "#sqrt(s) [GeV]", "4e #sqrt(s)", "Images/InvMass/InvMass_4eg.pdf");

  CompareHisto(H_s_eeReco, H_s_ee, H_s_eeReco_E, 0.1, 0.75, 0.2, 0.85, "Reco", "Gen", "Reco with E", "#sqrt(s) [MeV]", "ee #sqrt(s)", "Images/InvMass/InvMass_ee_E.pdf");
  std::cout << H_s_eeReco_E->Integral();

  // Plot dR between Gen ee
  HistoSave(H_GenEledR, "Images/GendR_ee.pdf", "dR ee", "Arbitrary units");

  // Plot eta for ML Vs non Ml
  Histo2Save(H_etaPho_RecoVsML, "Images/MLPho/etaMLVsnonML.pdf", "Reco #eta", "ML Reco #eta");
  HistoSave(H_MLPho_dphi, "Images/MLpho/phiML_nonML.pdf", "(#phi_{nonML}-#phi_{ML}) / #phi_{nonML}", "Arbitrary units");
  HistoSave(H_MLPho_deta, "Images/MLpho/etaML_nonML.pdf", "(#eta_{nonML}-#eta_{ML}) / #eta_{nonML}", "Arbitrary units");
  HistoSave(H_MLPho_dpt, "Images/MLpho/ptML_nonML.pdf", "(#pt_{nonML}-#pt_{ML}) / pt_{non_ML}", "Arbitrary units");

  Histo2Save(H_s_eeGenReco, "Images/InvMass/InvMass_ee_GenReco.pdf", "Gen  sqrt{s} [MeV]", "Reco  sqrt{s} [MeV]");

  HistoSave(H_EleMass_matched, "Images/InvMass/EleMatchedMass.pdf", "m [GeV]", "Arbitrary units");

  // Plot aa quantities
  Histo2Save(H, "Images/aa/ma1_Vs_ma2.pdf", "#sqrt(s) a1 [MeV]", "#sqrt(s) a2 [MeV]");
  // CompareHisto(H_pt5_1, H_pt5_2, 0.1, 0.75, 0.2, 0.85, "1", "2", "pt of low m Eles", "pt [GeV]", "Images/aa/pt.pdf");
  // CompareHisto(H_eta5_1, H_eta5_2, 0.1, 0.75, 0.2, 0.85, "1", "2", "eta of low m Eles", "eta", "Images/aa/eta.pdf");
  // CompareHisto(H_phi5_1, H_phi5_2, 0.1, 0.75, 0.2, 0.85, "1", "2", "phi of low m Eles", "phi", "Images/aa/phi.pdf");
  CompareHisto(H_phi5_1, H_phi5_2, H_phi5_3, H_phi5_4, 0.1, 0.75, 0.25, 0.85, "Low Mass", "Low Mass", "High Mass", "High Mass", "#phi of 4e", "phi", "Images/aa/phi.pdf");
  CompareHisto(H_eta5_1, H_eta5_2, H_eta5_3, H_eta5_4, 0.1, 0.75, 0.2, 0.85, "Low Mass", "Low Mass", "High Mass", "High Mass", "#eta of 4e", "eta", "Images/aa/eta.pdf");
  CompareHisto(H_pt5_1, H_pt5_2, H_pt5_3, H_pt5_4, 0.4, 0.75, 0.55, 0.85, "Low Mass", "Low Mass", "High Mass", "High Mass", "p_{t} of 4e", "p_{t}", "Images/aa/pt.pdf");

  Histo2Save(H_eem_dphideta, "Images/aa/dPhidEta_lowM.pdf", "#Delta#phi", "#Delta#ele");
  Histo2Save(H_eeM_dphideta, "Images/aa/dPhidEta_hiM.pdf", "#Delta#phi", "#Delta#ele");

  TH2F *H_eeM_dphideta_cu = new TH2F("", "#Delta#phi Vs #Delta#eta ee high m crystal units", 5, -2.5, 2.5, 5, -2.5, 2.5);
  TH2F *H_eem_dphideta_cu = new TH2F("", "#Delta#phi Vs #Delta#eta ee low m crystal units", 5, -2.5, 2.5, 5, -2.5, 2.5);

  for (int i = 1; i <= H_eem_dphideta->GetNbinsX(); ++i)
    for (int j = 1; j <= H_eem_dphideta->GetNbinsY(); ++j)
    {
      double content = H_eem_dphideta->GetBinContent(i, j);
      double x = H_eem_dphideta->GetXaxis()->GetBinCenter(i);
      double y = H_eem_dphideta->GetYaxis()->GetBinCenter(j);

      // Fill the new rebinned histogram, adjusting bin width
      H_eem_dphideta_cu->Fill(x / dPhi, y / dEta, content);
    }

  for (int i = 1; i <= H_eeM_dphideta->GetNbinsX(); ++i)
    for (int j = 1; j <= H_eeM_dphideta->GetNbinsY(); ++j)
    {
      double content = H_eeM_dphideta->GetBinContent(i, j);
      double x = H_eeM_dphideta->GetXaxis()->GetBinCenter(i);
      double y = H_eeM_dphideta->GetYaxis()->GetBinCenter(j);

      // Fill the new rebinned histogram, adjusting bin width
      H_eeM_dphideta_cu->Fill(x / dPhi, y / dEta, content);
    }

  Histo2Save(H_eem_dphideta_cu, "Images/aa/dPhidEta_lowM_cu.pdf", "#Delta#phi [cu]", "#Delta#ele [cu]");
  Histo2Save(H_eeM_dphideta_cu, "Images/aa/dPhidEta_hiM_cu.pdf", "#Delta#phi [cu]", "#Delta#ele [cu]");

  EffSave(H_EleEff_pt, "Images/Ele/RecoEff_pt.pdf");
  EffSave(H_EleEff_eta, "Images/Ele/RecoEff_eta.pdf");
  EffSave(H_PhoEff_pt, "Images/Pho/RecoEff_pt.pdf");
  EffSave(H_PhoEff_eta, "Images/Pho/RecoEff_eta.pdf");
  EffSave(H_MLPhoEff_pt, "Images/MLPho/RecoEff_pt.pdf");
  EffSave(H_MLPhoEff_eta, "Images/MLPho/RecoEff_eta.pdf");

  /*
  //Plot #of neighbors for different dR
  TGraph **G_EleNeigh = new TGraph*[n_dRs];
  TGraph **G_PhoNeigh = new TGraph*[n_dRs];
  TCanvas *c = new TCanvas("c", "c");
  TLegend *l = new TLegend(0.7, 0.7, 0.9, 0.9);

  for(UInt_t i = 0; i < n_dRs; i++)
    {
      G_EleNeigh[i] = new TGraph();
      for(UInt_t j = 0; j <= n_EleNeigh; j++)
  {
    G_EleNeigh[i]->SetPoint(j, j, H_EleNeigh[i]->GetBinContent(j + 1));
  }
      G_EleNeigh[i]->SetMarkerStyle(20 + i);  // Different marker styles
      G_EleNeigh[i]->SetMarkerColor(TColor::GetColorPalette(i)); // Palette color
      G_EleNeigh[i]->SetLineColor(TColor::GetColorPalette(i));
      G_EleNeigh[i]->SetTitle(Form("dR = %.3f", dRs[i]));
      l->AddEntry(G_EleNeigh[i], Form("dR = %.3f", dRs[i]), "l");
    }

  G_EleNeigh[0]->Draw("APL");
  for(UInt_t i = 1; i < n_dRs; i++) G_EleNeigh[i]->Draw("PL SAME");
  l->Draw();
  c->Update();
  c->SaveAs("Images/Ele/Neighbors.pdf");

  HistoSave(H_EleNeigh[0], "Prova.pdf", "", "");

  for(UInt_t i = 0; i < n_dRs; i++)
    {
      G_PhoNeigh[i] = new TGraph();
      for(UInt_t j = 0; j < n_PhoNeigh; j++) G_PhoNeigh[i]->SetPoint(j, j, H_PhoNeigh[i]->GetBinContent(j));
      G_PhoNeigh[i]->SetTitle(Form("dR = %.3f", dRs[i]));
    }
  TCanvas *c1 = new TCanvas("", "");
  THStack *HS_EleNeigh = new THStack("hs", "Ele Neigh");
  for(UInt_t i = 0; i < n_dRs; i++) HS_EleNeigh->Add(H_EleNeigh[i]);
  HS_EleNeigh->Draw("HIST");
  c1->SaveAs("ProvaStacked.pdf");
  */

  std::cout << "Finished!\n";
}

void HistoSave(TH1F *histo, TString filename, TString XLabel, TString YLabel)
{
  TCanvas *c = new TCanvas("canvas", "c");
  histo->SetXTitle(XLabel);
  histo->SetYTitle(YLabel);
  histo->Scale(1. / histo->Integral(), "width");

  histo->SetMinimum(0);
  gStyle->SetOptStat("nuo");

  histo->Draw("HIST");
  c->SaveAs(filename);
  filename.Resize(filename.Length() - 3);
  filename += "png";
  c->SaveAs(filename);
  delete c;
}

void Histo2Save(TH2F *histo, TString filename, TString XLabel, TString YLabel)
{
  TCanvas *c = new TCanvas("canvas", "c");
  histo->SetXTitle(XLabel);
  histo->SetYTitle(YLabel);
  histo->SetMinimum(0);
  histo->Draw("COLZ");
  c->SaveAs(filename);
  filename.Resize(filename.Length() - 3);
  filename += "png";
  c->SaveAs(filename);
  delete c;
}

void CompareHisto(TH1F *histo1, TH1F *histo2, Float_t x1, Float_t y1, Float_t x2, Float_t y2, TString lab1, TString lab2, TString Xlab, TString title, TString path)
{
  TCanvas *c = new TCanvas("canvas", "c");
  histo1->SetTitle(title);

  histo1->SetXTitle(Xlab);

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2);

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2);

  histo1->Scale(1. / histo1->Integral(), "width");
  histo2->Scale(1. / histo2->Integral(), "width");

  Double_t max1 = histo1->GetMaximum();
  Double_t max2 = histo2->GetMaximum();

  if (max1 < max2)
  {
    max1 = max2;
  }

  histo1->SetMaximum(max1 * 1.1);
  histo1->SetMinimum(0);
  histo1->Draw("HIST");
  histo2->Draw("HIST SAME");

  histo1->SetYTitle("Arbitrary units");

  TLegend *leg = new TLegend(x1, y1, x2, y2);
  leg->AddEntry(histo1, lab1, "l");
  leg->AddEntry(histo2, lab2, "l");
  leg->Draw();

  histo1->SetStats();
  histo2->SetStats();

  TPaveStats *stats1 = (TPaveStats *)histo1->GetListOfFunctions()->FindObject("stats");
  TPaveStats *stats2 = (TPaveStats *)histo2->GetListOfFunctions()->FindObject("stats");

  if (stats1 && stats2)
  {
    // Move the second statistics box to avoid overlap
    stats1->SetX1NDC(0.7);  // Set new X1 position
    stats1->SetY1NDC(0.6);  // Set new Y1 position
    stats1->SetX2NDC(0.9);  // Set new X2 position
    stats1->SetY2NDC(0.75); // Set new Y2 position

    stats2->SetX1NDC(0.7);  // Set new X1 position
    stats2->SetY1NDC(0.75); // Set new Y1 position
    stats2->SetX2NDC(0.9);  // Set new X2 position
    stats2->SetY2NDC(0.9);  // Set new Y2 position

    // stats1->Draw();
    // stats2->Draw();
    histo1->SetStats(0);
    histo2->SetStats(0);
  }

  // Update the canvas to apply changes
  c->Update();

  c->SaveAs(path);
  path.Resize(path.Length() - 3);
  path += "png";
  c->SaveAs(path);

  delete c;
}

// Overload to compare 3
void CompareHisto(TH1F *histo1, TH1F *histo2, TH1F *histo3, Float_t x1, Float_t y1, Float_t x2, Float_t y2, TString lab1, TString lab2, TString lab3, TString Xlab, TString title, TString path)
{
  TCanvas *c = new TCanvas("canvas", "c");
  histo1->SetTitle(title);

  gStyle->SetOptStat(11);

  histo1->Scale(1. / histo1->Integral(), "width");
  histo2->Scale(1. / histo2->Integral(), "width");
  histo3->Scale(1. / histo3->Integral(), "width");

  Double_t max1 = histo1->GetMaximum();
  Double_t max2 = histo2->GetMaximum();
  Double_t max3 = histo3->GetMaximum();

  if (max1 < max2)
    max1 = max2;
  if (max1 < max3)
    max1 = max3;
  histo1->SetMaximum(max1 * 1.1);

  histo1->SetXTitle(Xlab);
  histo1->SetYTitle("Arbitrary units");

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2);
  histo1->SetLineStyle(0);
  histo1->SetMinimum(0);

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2);
  histo2->SetLineStyle(0);
  histo2->SetMinimum(0);

  histo3->SetLineColor(kGreen);
  histo3->SetLineWidth(2);
  histo3->SetLineStyle(0);
  histo3->SetMinimum(0);

  histo1->SetStats(kTRUE);
  histo2->SetStats(kTRUE);
  histo3->SetStats(kTRUE);

  histo1->Draw("HIST");
  histo2->Draw("HIST SAME");
  histo3->Draw("HIST SAME");

  histo1->SetYTitle("Entries");

  TLegend *leg = new TLegend(x1, y1, x2, y2);
  leg->AddEntry(histo1, lab1, "l");
  leg->AddEntry(histo2, lab2, "l");
  leg->AddEntry(histo3, lab3, "l");

  leg->SetTextSize(0.03);

  leg->Draw();

  gPad->Update();
  TPaveStats *stats1 = (TPaveStats *)histo1->GetListOfFunctions()->FindObject("stats");
  TPaveStats *stats2 = (TPaveStats *)histo2->GetListOfFunctions()->FindObject("stats");
  TPaveStats *stats3 = (TPaveStats *)histo3->GetListOfFunctions()->FindObject("stats");

  if (stats1 && stats2 && stats3)
  {
    // Move the second statistics box to avoid overlap
    stats1->SetX1NDC(0.7); // Set new X1 position
    stats1->SetY1NDC(0.8); // Set new Y1 position
    stats1->SetX2NDC(0.9); // Set new X2 position
    stats1->SetY2NDC(0.9); // Set new Y2 position
    stats1->SetOptStat(11);
    stats1->Draw();
    gPad->Update();

    stats2->SetX1NDC(0.7); // Set new X1 position
    stats2->SetY1NDC(0.7); // Set new Y1 position
    stats2->SetX2NDC(0.9); // Set new X2 position
    stats2->SetY2NDC(0.8); // Set new Y2 position
    stats2->SetOptStat(11);
    stats2->Draw();
    gPad->Update();

    stats3->SetX1NDC(0.7); // Set new X1 position
    stats3->SetY1NDC(0.6); // Set new Y1 position
    stats3->SetX2NDC(0.9); // Set new X2 position
    stats3->SetY2NDC(0.7); // Set new Y2 position
    stats3->SetOptStat(11);
    stats3->Draw();
    gPad->Update();

    // histo1->SetStats(0);
    // histo2->SetStats(0);
    // histo3->SetStats(0);
  }
  else
    std::cout << "no stats!" << stats1 << " " << stats2 << " " << stats3 << std::endl;

  // Update the canvas to apply changes
  c->Update();

  c->SaveAs(path);
  path.Resize(path.Length() - 3);
  path += "png";
  c->SaveAs(path);

  delete c;
}

// Overload to compare 4
void CompareHisto(TH1F *histo1, TH1F *histo2, TH1F *histo3, TH1F *histo4, Float_t x1, Float_t y1, Float_t x2, Float_t y2, TString lab1, TString lab2, TString lab3, TString lab4, TString Xlab, TString title, TString path)
{
  TCanvas *c = new TCanvas("canvas", "c");
  histo1->SetTitle(title);

  gStyle->SetOptStat(11);

  histo1->Scale(1. / histo1->Integral(), "width");
  histo2->Scale(1. / histo2->Integral(), "width");
  histo3->Scale(1. / histo3->Integral(), "width");
  histo4->Scale(1. / histo4->Integral(), "width");

  Double_t max1 = histo1->GetMaximum();
  Double_t max2 = histo2->GetMaximum();
  Double_t max3 = histo3->GetMaximum();
  Double_t max4 = histo4->GetMaximum();

  if (max1 < max2)
    max1 = max2;
  if (max1 < max3)
    max1 = max3;
  if (max1 < max4)
    max1 = max4;
  histo1->SetMaximum(max1 * 1.1);

  histo1->SetXTitle(Xlab);
  histo1->SetYTitle("Arbitrary units");

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2);
  histo1->SetLineStyle(0);
  histo1->SetMinimum(0);

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2);
  histo2->SetLineStyle(0);
  histo2->SetMinimum(0);

  histo3->SetLineColor(kGreen);
  histo3->SetLineWidth(2);
  histo3->SetLineStyle(0);
  histo3->SetMinimum(0);

  histo4->SetLineColor(kMagenta);
  histo4->SetLineWidth(2);
  histo4->SetLineStyle(0);
  histo4->SetMinimum(0);

  histo1->SetStats(kTRUE);
  histo2->SetStats(kTRUE);
  histo3->SetStats(kTRUE);
  histo4->SetStats(kTRUE);

  histo1->Draw("HIST");
  histo2->Draw("HIST SAME");
  histo3->Draw("HIST SAME");
  histo4->Draw("HIST SAME");

  histo1->SetYTitle("Entries");

  TLegend *leg = new TLegend(x1, y1, x2, y2);
  leg->AddEntry(histo1, lab1, "l");
  leg->AddEntry(histo2, lab2, "l");
  leg->AddEntry(histo3, lab3, "l");
  leg->AddEntry(histo4, lab4, "l");

  leg->SetTextSize(0.03);

  leg->Draw();

  gPad->Update();
  TPaveStats *stats1 = (TPaveStats *)histo1->GetListOfFunctions()->FindObject("stats");
  TPaveStats *stats2 = (TPaveStats *)histo2->GetListOfFunctions()->FindObject("stats");
  TPaveStats *stats3 = (TPaveStats *)histo3->GetListOfFunctions()->FindObject("stats");
  TPaveStats *stats4 = (TPaveStats *)histo4->GetListOfFunctions()->FindObject("stats");

  if (stats1 && stats2 && stats3)
  {
    // Move the second statistics box to avoid overlap
    stats1->SetX1NDC(0.7); // Set new X1 position
    stats1->SetY1NDC(0.8); // Set new Y1 position
    stats1->SetX2NDC(0.9); // Set new X2 position
    stats1->SetY2NDC(0.9); // Set new Y2 position
    stats1->SetOptStat(11);
    stats1->Draw();
    gPad->Update();

    stats2->SetX1NDC(0.7); // Set new X1 position
    stats2->SetY1NDC(0.7); // Set new Y1 position
    stats2->SetX2NDC(0.9); // Set new X2 position
    stats2->SetY2NDC(0.8); // Set new Y2 position
    stats2->SetOptStat(11);
    stats2->Draw();
    gPad->Update();

    stats3->SetX1NDC(0.7); // Set new X1 position
    stats3->SetY1NDC(0.6); // Set new Y1 position
    stats3->SetX2NDC(0.9); // Set new X2 position
    stats3->SetY2NDC(0.7); // Set new Y2 position
    stats3->SetOptStat(11);
    stats3->Draw();
    gPad->Update();

    stats4->SetX1NDC(0.7); // Set new X1 position
    stats4->SetY1NDC(0.6); // Set new Y1 position
    stats4->SetX2NDC(0.9); // Set new X2 position
    stats4->SetY2NDC(0.7); // Set new Y2 position
    stats4->SetOptStat(11);
    stats4->Draw();
    gPad->Update();

    // histo1->SetStats(0);
    // histo2->SetStats(0);
    // histo3->SetStats(0);
  }
  else
    std::cout << "no stats!" << stats1 << " " << stats2 << " " << stats3 << std::endl;

  // Update the canvas to apply changes
  c->Update();

  c->SaveAs(path);
  path.Resize(path.Length() - 3);
  path += "png";
  c->SaveAs(path);

  delete c;
}

void SaveHistoDiff(TH1F *hist1, TH1F *hist2, TString title, TString stats, TString path)
{
  TCanvas *c = new TCanvas("canvas", "c");

  TH1F histo_d = *hist2;
  histo_d.Add(hist1, -1);
  int nBin = histo_d.GetNbinsX();
  int Xmax = histo_d.GetXaxis()->GetXmax();
  int Xmin = histo_d.GetXaxis()->GetXmin();
  TH1F *H_d = new TH1F(title, stats, Xmax - Xmin + 5, Xmin - 2, Xmax + 2);

  for (int i = 1; i <= nBin; i++)
  {
    H_d->Fill(histo_d.GetBinContent(i));
  }
  H_d->Draw();
  c->SaveAs(path);
  delete c;

  return;
}

int isDaughter(Int_t *partPdgIds, Int_t *partMothIds, Int_t id)
{
  if (abs(partPdgIds[partMothIds[id]]) == 5000001)
  {
    return 1;
  }
  if (partMothIds[id] == -1)
  {
    return 0;
  }

  return isDaughter(partPdgIds, partMothIds, partMothIds[id]);
}

void CompareHistoZoom(TH1F *histo1, TH1F *histo2, Float_t x1, Float_t y1, Float_t x2, Float_t y2, TString lab1, TString lab2, TString Xlab, TString title, TString path, Float_t xz_1, Float_t yz_1, Float_t xz_2, Float_t yz_2, Float_t xmin, Float_t xmax)
{
  TCanvas *c = new TCanvas("canvas", "c");
  histo1->SetTitle(title);

  histo1->SetXTitle(Xlab);

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2);

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2);

  histo1->Scale(1. / histo1->Integral(), "width");
  histo2->Scale(1. / histo2->Integral(), "width");

  Double_t max1 = histo1->GetMaximum();
  Double_t max2 = histo2->GetMaximum();

  if (max1 < max2)
  {
    max1 = max2;
  }

  histo1->SetMaximum(max1 * 1.1);
  histo1->SetMinimum(0);
  histo1->Draw("HIST");
  histo2->Draw("HIST SAME");

  histo1->SetYTitle("Arbitrary units");

  TLegend *leg = new TLegend(x1, y1, x2, y2);
  leg->AddEntry(histo1, lab1, "l");
  leg->AddEntry(histo2, lab2, "l");
  leg->Draw();

  histo1->SetStats();
  histo2->SetStats();

  TPaveStats *stats1 = (TPaveStats *)histo1->GetListOfFunctions()->FindObject("stats");
  TPaveStats *stats2 = (TPaveStats *)histo2->GetListOfFunctions()->FindObject("stats");

  if (stats1 && stats2)
  {
    // Move the second statistics box to avoid overlap
    stats1->SetX1NDC(0.7);  // Set new X1 position
    stats1->SetY1NDC(0.6);  // Set new Y1 position
    stats1->SetX2NDC(0.9);  // Set new X2 position
    stats1->SetY2NDC(0.75); // Set new Y2 position

    stats2->SetX1NDC(0.7);  // Set new X1 position
    stats2->SetY1NDC(0.75); // Set new Y1 position
    stats2->SetX2NDC(0.9);  // Set new X2 position
    stats2->SetY2NDC(0.9);  // Set new Y2 position

    // stats1->Draw();
    // stats2->Draw();
    histo1->SetStats(0);
    histo2->SetStats(0);
  }

  // Update the canvas to apply changes
  c->Update();

  TPad *p = new TPad("p", "p", xz_1, yz_1, xz_2, yz_2);
  p->Draw();
  p->cd();
  TH1F *hc1 = (TH1F *)histo1->DrawCopy("HIST");
  TH1F *hc2 = (TH1F *)histo2->DrawCopy("HIST SAME");
  hc1->GetXaxis()->SetRangeUser(xmin, xmax);
  hc2->GetXaxis()->SetRangeUser(xmin, xmax);

  c->SaveAs(path);
  path.Resize(path.Length() - 3);
  path += "png";
  c->SaveAs(path);

  delete c;
}

void EffSave(TEfficiency *histo, TString filename)
{
  TCanvas *c = new TCanvas("canvas", "c");
  histo->Draw("AP");
  c->SaveAs(filename);
  filename.Resize(filename.Length() - 3);
  filename += "png";
  c->SaveAs(filename);
  delete c;
}
