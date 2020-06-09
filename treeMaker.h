//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed May  6 15:34:15 2020 by ROOT version 6.18/04
// from TTree treeMaker/tree
// found on file: NCUGlobalTuples_1.root
//////////////////////////////////////////////////////////

#ifndef treeMaker_h
#define treeMaker_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TClonesArray.h"
#include "vector"

class treeMaker {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         pu_nTrueInt;
   Int_t           pu_nPUVert;
   Bool_t          isData;
   ULong64_t       eventId;
   ULong64_t       runId;
   ULong64_t       lumiSection;
   ULong64_t       bunchXing;
   Int_t           nVtx;
   TClonesArray    *vertexP3;
   Float_t         pfMetCorrPt;
   Float_t         pfMetCorrPhi;
   Float_t         pfMetCorrSumEt;
   Float_t         pfMetCorrSig;
   vector<float>   *pfMetCorrUnc;
   Float_t         pfMetRawPt;
   Float_t         pfMetRawPhi;
   Float_t         pfMetRawSumEt;
   Float_t         pfMetRawCov00;
   Float_t         pfMetRawCov01;
   Float_t         pfMetRawCov10;
   Float_t         pfMetRawCov11;
   Float_t         pfmvaMetPt;
   Float_t         pfmvaMetPhi;
   Float_t         pfmvaMetSumEt;
   Float_t         pfmvaMetSig;
   Float_t         pfpuppiMETPt;
   Float_t         pfpuppiMETPhi;
   Float_t         pfpuppiMETSumEt;
   Float_t         pfpuppiMETSig;
   vector<float>   *pfpuppiMETUnc;
   Int_t           hlt_nTrigs;
   vector<bool>    *hlt_trigResult;
   vector<string>  *hlt_trigName;
   vector<int>     *hlt_trigPrescale;
   Int_t           hlt_nfilters;
   vector<bool>    *hlt_filterResult;
   Bool_t          hlt_hbhet;
   vector<string>  *hlt_filterName;
   Bool_t          hlt_filterbadChCandidate;
   Bool_t          hlt_filterbadPFMuon;
   Bool_t          hlt_filterbadGlobalMuon;
   Bool_t          hlt_filtercloneGlobalMuon;
   Float_t         genMET_true;
   Float_t         genMET_calo;
   Float_t         genMET_caloNonPrompt;
   Float_t         ptHat;
   Float_t         mcWeight;
   Float_t         HT;
   vector<float>   *pdf;
   Float_t         originalLHEweight;
   vector<float>   *pdfscaleSysWeights;
   Int_t           nGenPar;
   TClonesArray    *genParP4;
   vector<int>     *genParQ;
   vector<int>     *genParId;
   vector<int>     *genParSt;
   vector<int>     *genMomParId;
   vector<int>     *genParIndex;
   vector<int>     *genNMo;
   vector<int>     *genNDa;
   vector<int>     *genMo1;
   vector<int>     *genMo2;
   vector<int>     *genDa1;
   vector<int>     *genDa2;
   vector<int>     *genStFlag;
   Int_t           ak4nGenJet;
   TClonesArray    *ak4GenJetP4;
   Int_t           ak8nGenJet;
   TClonesArray    *ak8GenJetP4;
   Float_t         eleRho;
   Int_t           nEle;
   TClonesArray    *eleP4;
   vector<bool>    *eleInBarrel;
   vector<bool>    *eleInEndcap;
   vector<int>     *eleCharge;
   vector<int>     *eleChargeConsistent;
   vector<float>   *elecaloEnergy;
   vector<float>   *eleScEt;
   vector<float>   *eleScEn;
   vector<float>   *eleScPreEn;
   vector<float>   *eleScEta;
   vector<float>   *eleScPhi;
   vector<float>   *eleScRawEn;
   vector<float>   *eleScEtaWidth;
   vector<float>   *eleScPhiWidth;
   vector<float>   *eleR9;
   vector<float>   *eleHoverE;
   vector<float>   *eleD0;
   vector<float>   *eleDz;
   vector<float>   *eleEoverP;
   vector<float>   *eleBrem;
   vector<float>   *eledEtaAtVtx;
   vector<float>   *eledPhiAtVtx;
   vector<float>   *eleSigmaIEtaIEta;
   vector<float>   *eleSigmaIEtaIPhi;
   vector<float>   *eleSigmaIPhiIPhi;
   vector<bool>    *eleConvVeto;
   vector<int>     *eleMissHits;
   vector<float>   *eleEoverPInv;
   vector<float>   *eleEtaseedAtVtx;
   vector<float>   *eleE1x5;
   vector<float>   *eleE2x5;
   vector<float>   *eleE5x5;
   vector<float>   *eleSigmaIEtaIEtaFull5x5;
   vector<float>   *eleE1x5Full5x5;
   vector<float>   *eleE2x5Full5x5;
   vector<float>   *eleE5x5Full5x5;
   vector<float>   *eleR9Full5x5;
   vector<float>   *eleChHadIso;
   vector<float>   *eleNeHadIso;
   vector<float>   *eleGamIso;
   vector<float>   *elePUPt;
   vector<float>   *eleEcalPFClusterIso;
   vector<float>   *eleHcalPFClusterIso;
   vector<float>   *eleMiniIso_ch;
   vector<float>   *eleMiniIso_nh;
   vector<float>   *eleMiniIso_ph;
   vector<float>   *eleMiniIso_pu;
   vector<float>   *eleMiniIso_r;
   vector<float>   *eleMiniIsoBeta;
   vector<float>   *eleMiniIsoEA;
   vector<bool>    *eleEcalDrivenSeed;
   vector<bool>    *eleEcalDriven;
   vector<float>   *eleDr03EcalRecHitSumEt;
   vector<float>   *eleDr03HcalDepth1TowerSumEt;
   vector<float>   *eleDr03HcalDepth2TowerSumEt;
   vector<float>   *eleDr03HcalTowerSumEt;
   vector<float>   *eleDr03TkSumPt;
   vector<bool>    *eleIsPassVeto;
   vector<bool>    *eleIsPassLoose;
   vector<bool>    *eleIsPassMedium;
   vector<bool>    *eleIsPassTight;
   vector<bool>    *eleIsPassHEEP;
   vector<bool>    *eleIsPassVetoNoIso;
   vector<bool>    *eleIsPassLooseNoIso;
   vector<bool>    *eleIsPassMediumNoIso;
   vector<bool>    *eleIsPassTightNoIso;
   vector<bool>    *eleIsPassHEEPNoIso;
   vector<bool>    *eleIsPassMVAMedium;
   vector<bool>    *eleIsPassMVATight;
   vector<float>   *eleMVAValue;
   vector<int>     *eleMVACategory;
   Int_t           nMu;
   TClonesArray    *muP4;
   vector<int>     *muType;
   vector<int>     *muCharge;
   vector<bool>    *isGlobalMuon;
   vector<bool>    *isTrackerMuon;
   vector<bool>    *isPFMuon;
   vector<bool>    *isTightMuon;
   vector<bool>    *isLooseMuon;
   vector<bool>    *isMediumMuon;
   vector<bool>    *isSoftMuon;
   vector<bool>    *isHighPtMuon;
   vector<bool>    *isCustomTrackerMuon;
   vector<int>     *muITrkID;
   vector<int>     *muSegID;
   vector<int>     *muNSegs;
   vector<bool>    *muGood;
   vector<bool>    *muIsGood;
   vector<float>   *muTrkPt;
   vector<float>   *muTrkPtErr;
   vector<float>   *mudxy;
   vector<float>   *mudz;
   vector<float>   *musegmentCompatibility;
   vector<float>   *muchi2LocalPosition;
   vector<float>   *mutrkKink;
   vector<float>   *muInnerdxy;
   vector<float>   *muInnerdz;
   vector<int>     *muTrkLayers;
   vector<int>     *muPixelLayers;
   vector<int>     *muPixelHits;
   vector<float>   *muTrkQuality;
   vector<int>     *muHits;
   vector<float>   *muChi2NDF;
   vector<float>   *muInnervalidFraction;
   vector<int>     *muMatches;
   vector<float>   *muTrkIso;
   vector<float>   *muHcalIso;
   vector<float>   *muEcalIso;
   vector<float>   *muChHadIso;
   vector<float>   *muNeHadIso;
   vector<float>   *muGamIso;
   vector<float>   *muPUPt;
   vector<float>   *muInnerTrkPt;
   vector<float>   *muMiniIso_ch;
   vector<float>   *muMiniIso_nh;
   vector<float>   *muMiniIso_ph;
   vector<float>   *muMiniIso_pu;
   vector<float>   *muMiniIso_r;
   vector<float>   *muMiniIsoBeta;
   vector<float>   *muMiniIsoEA;
   Int_t           HPSTau_n;
   vector<float>   *taupt;
   TClonesArray    *HPSTau_4Momentum;
   TClonesArray    *HPSTau_Vposition;
   vector<bool>    *HPSTau_leadPFChargedHadrCand;
   vector<bool>    *HPSTau_leadPFChargedHadrCand_trackRef;
   vector<bool>    *disc_againstElectronLoose;
   vector<bool>    *disc_againstElectronMedium;
   vector<bool>    *disc_againstElectronTight;
   vector<bool>    *disc_againstElectronLooseMVA5;
   vector<bool>    *disc_againstElectronMediumMVA5;
   vector<bool>    *disc_againstElectronTightMVA5;
   vector<bool>    *disc_againstElectronVLooseMVA5;
   vector<bool>    *disc_againstElectronVTightMVA5;
   vector<bool>    *disc_againstMuonLoose;
   vector<bool>    *disc_againstMuonMedium;
   vector<bool>    *disc_againstMuonTight;
   vector<bool>    *disc_againstMuonLoose2;
   vector<bool>    *disc_againstMuonMedium2;
   vector<bool>    *disc_againstMuonTight2;
   vector<bool>    *disc_againstMuonLooseMVA;
   vector<bool>    *disc_againstMuonMediumMVA;
   vector<bool>    *disc_againstMuonTightMVA;
   vector<bool>    *disc_againstMuonLoose3;
   vector<bool>    *disc_againstMuonTight3;
   vector<bool>    *disc_byVLooseCombinedIsolationDeltaBetaCorr;
   vector<bool>    *disc_byLooseCombinedIsolationDeltaBetaCorr;
   vector<bool>    *disc_byMediumCombinedIsolationDeltaBetaCorr;
   vector<bool>    *disc_byTightCombinedIsolationDeltaBetaCorr;
   vector<bool>    *disc_byLooseIsolation;
   vector<bool>    *disc_byVLooseIsolationMVA3newDMwLT;
   vector<bool>    *disc_byLooseIsolationMVA3newDMwLT;
   vector<bool>    *disc_byMediumIsolationMVA3newDMwLT;
   vector<bool>    *disc_byTightIsolationMVA3newDMwLT;
   vector<bool>    *disc_byVTightIsolationMVA3newDMwLT;
   vector<bool>    *disc_byVVTightIsolationMVA3newDMwLT;
   vector<bool>    *disc_byVLooseIsolationMVA3newDMwoLT;
   vector<bool>    *disc_byLooseIsolationMVA3newDMwoLT;
   vector<bool>    *disc_byMediumIsolationMVA3newDMwoLT;
   vector<bool>    *disc_byTightIsolationMVA3newDMwoLT;
   vector<bool>    *disc_byVTightIsolationMVA3newDMwoLT;
   vector<bool>    *disc_byVVTightIsolationMVA3newDMwoLT;
   vector<bool>    *disc_byVLooseIsolationMVA3oldDMwLT;
   vector<bool>    *disc_byLooseIsolationMVA3oldDMwLT;
   vector<bool>    *disc_byMediumIsolationMVA3oldDMwLT;
   vector<bool>    *disc_byTightIsolationMVA3oldDMwLT;
   vector<bool>    *disc_byVTightIsolationMVA3oldDMwLT;
   vector<bool>    *disc_byVVTightIsolationMVA3oldDMwLT;
   vector<bool>    *disc_byVLooseIsolationMVA3oldDMwoLT;
   vector<bool>    *disc_byLooseIsolationMVA3oldDMwoLT;
   vector<bool>    *disc_byMediumIsolationMVA3oldDMwoLT;
   vector<bool>    *disc_byTightIsolationMVA3oldDMwoLT;
   vector<bool>    *disc_byVTightIsolationMVA3oldDMwoLT;
   vector<bool>    *disc_byVVTightIsolationMVA3oldDMwoLT;
   vector<bool>    *disc_byLooseCombinedIsolationDeltaBetaCorr3Hits;
   vector<bool>    *disc_byMediumCombinedIsolationDeltaBetaCorr3Hits;
   vector<bool>    *disc_byTightCombinedIsolationDeltaBetaCorr3Hits;
   vector<bool>    *disc_decayModeFinding;
   vector<bool>    *disc_decayModeFindingNewDMs;
   vector<float>   *disc_chargedIsoPtSum;
   vector<float>   *disc_neutralIsoPtSum;
   vector<float>   *disc_puCorrPtSum;
   vector<float>   *HPSTau_NewVz;
   vector<int>     *HPSTau_charge;
   Int_t           nPho;
   TClonesArray    *phoP4;
   vector<bool>    *phoIsPassTight;
   vector<bool>    *phoIsPassLoose;
   vector<bool>    *phoIsPassMedium;
   vector<float>   *phoIDMVA;
   vector<float>   *phoSCE;
   vector<float>   *phoSCRawE;
   vector<float>   *phoSCEta;
   vector<float>   *phoSCPhi;
   vector<float>   *phoSCEtaWidth;
   vector<float>   *phoSCPhiWidth;
   vector<float>   *phoSCBrem;
   vector<int>     *phohasPixelSeed;
   vector<int>     *phoEleVeto;
   vector<float>   *phoR9;
   vector<float>   *phoHoverE;
   vector<float>   *phoSigmaIEtaIEta;
   vector<float>   *phoSigmaIEtaIPhi;
   vector<float>   *phoSigmaIPhiIPhi;
   vector<float>   *phoSigmaIEtaIEtaFull5x5;
   vector<float>   *phoR9Full5x5;
   vector<float>   *phoPFChIso;
   vector<float>   *phoPFPhoIso;
   vector<float>   *phoPFNeuIso;
   Int_t           THINnJet;
   TClonesArray    *THINjetP4;
   Float_t         THINjetRho;
   Int_t           THINjetNPV;
   TClonesArray    *THINgenjetP4;
   vector<float>   *THINgenjetEM;
   vector<float>   *THINgenjetHAD;
   vector<float>   *THINgenjetINV;
   vector<float>   *THINgenjetAUX;
   vector<float>   *THINmatchedDR;
   vector<float>   *THINjetRawFactor;
   TClonesArray    *THINunCorrJetP4;
   vector<float>   *THINjetArea;
   vector<float>   *THINjetCorrUncUp;
   vector<float>   *THINjetCorrUncDown;
   vector<int>     *THINjetCharge;
   vector<int>     *THINjetPartonFlavor;
   vector<int>     *THINjetHadronFlavor;
   vector<bool>    *THINjetPassIDLoose;
   vector<bool>    *THINjetPassIDTight;
   vector<float>   *THINjetCEmEF;
   vector<float>   *THINjetCHadEF;
   vector<float>   *THINjetPhoEF;
   vector<float>   *THINjetNEmEF;
   vector<float>   *THINjetNHadEF;
   vector<float>   *THINjetEleEF;
   vector<float>   *THINjetMuoEF;
   vector<int>     *THINjetCMulti;
   vector<int>     *THINjetEleMulti;
   vector<int>     *THINjetMuoMulti;
   vector<float>   *THINjetSSV;
   vector<float>   *THINjetCSV;
   vector<float>   *THINjetSSVHE;
   vector<float>   *THINjetCISVV2;
   vector<float>   *THINjetTCHP;
   vector<float>   *THINjetTCHE;
   vector<float>   *THINjetJP;
   vector<float>   *THINjetJBP;
   vector<float>   *THINPUJetID;
   vector<bool>    *THINisPUJetIDLoose;
   vector<bool>    *THINisPUJetIDMedium;
   vector<bool>    *THINisPUJetIDTight;
   Int_t           FATnJet;
   TClonesArray    *FATjetP4;
   Float_t         FATjetRho;
   Int_t           FATjetNPV;
   TClonesArray    *FATgenjetP4;
   vector<float>   *FATgenjetEM;
   vector<float>   *FATgenjetHAD;
   vector<float>   *FATgenjetINV;
   vector<float>   *FATgenjetAUX;
   vector<float>   *FATmatchedDR;
   vector<float>   *FATjetRawFactor;
   TClonesArray    *FATunCorrJetP4;
   vector<float>   *FATjetArea;
   vector<float>   *FATjetCorrUncUp;
   vector<float>   *FATjetCorrUncDown;
   vector<int>     *FATjetCharge;
   vector<int>     *FATjetPartonFlavor;
   vector<int>     *FATjetHadronFlavor;
   vector<bool>    *FATjetPassIDLoose;
   vector<bool>    *FATjetPassIDTight;
   vector<float>   *FATjetCEmEF;
   vector<float>   *FATjetCHadEF;
   vector<float>   *FATjetPhoEF;
   vector<float>   *FATjetNEmEF;
   vector<float>   *FATjetNHadEF;
   vector<float>   *FATjetEleEF;
   vector<float>   *FATjetMuoEF;
   vector<int>     *FATjetCMulti;
   vector<int>     *FATjetEleMulti;
   vector<int>     *FATjetMuoMulti;
   vector<float>   *FATjetSSV;
   vector<float>   *FATjetCSV;
   vector<float>   *FATjetSSVHE;
   vector<float>   *FATjetCISVV2;
   vector<float>   *FATjetTCHP;
   vector<float>   *FATjetTCHE;
   vector<float>   *FATjetJP;
   vector<float>   *FATjetJBP;
   vector<float>   *FATjetTau1;
   vector<float>   *FATjetTau2;
   vector<float>   *FATjetTau3;
   vector<float>   *FATjetTau21;
   vector<float>   *FATjetSDmass;
   vector<float>   *FATjetGenSDmass;
   vector<int>     *FATnSubSDJet;
   vector<vector<int> > *FATsubjetSDFatJetIndex;
   vector<vector<float> > *FATsubjetSDPx;
   vector<vector<float> > *FATsubjetSDPy;
   vector<vector<float> > *FATsubjetSDPz;
   vector<vector<float> > *FATsubjetSDE;
   vector<vector<float> > *FATsubjetSDRawFactor;
   vector<vector<int> > *FATsubjetSDPartonFlavor;
   vector<vector<int> > *FATsubjetSDHadronFlavor;
   vector<vector<float> > *FATsubjetSDCSV;
   vector<float>   *FATjetPRmass;
   vector<float>   *FATjetPRmassL2L3Corr;
   vector<float>   *FATjetPuppiTau1;
   vector<float>   *FATjetPuppiTau2;
   vector<float>   *FATjetPuppiTau3;
   vector<float>   *FATjetPuppiSDmass;
   TClonesArray    *FATjetPuppiP4;
   TClonesArray    *FATjetPuppiSDRawP4;
   vector<int>     *FATnSubSDPuppiJet;
   vector<vector<int> > *FATsubjetSDPuppiFatJetIndex;
   vector<vector<float> > *FATsubjetSDPuppiPx;
   vector<vector<float> > *FATsubjetSDPuppiPy;
   vector<vector<float> > *FATsubjetSDPuppiPz;
   vector<vector<float> > *FATsubjetSDPuppiE;
   vector<vector<float> > *FATsubjetSDPuppiCSV;
   vector<float>   *FATjet_DoubleSV;
   vector<int>     *FATjet_nSV;
   vector<vector<float> > *FATjet_SVMass;
   Int_t           ADDnJet;
   TClonesArray    *ADDjetP4;
   vector<float>   *ADDjet_DoubleSV;
   vector<int>     *ADDjet_nSV;
   vector<vector<float> > *ADDjet_SVMass;
   Int_t           AK8PuppinJet;
   TClonesArray    *AK8PuppijetP4;
   Float_t         AK8PuppijetRho;
   Int_t           AK8PuppijetNPV;
   TClonesArray    *AK8PuppigenjetP4;
   vector<float>   *AK8PuppigenjetEM;
   vector<float>   *AK8PuppigenjetHAD;
   vector<float>   *AK8PuppigenjetINV;
   vector<float>   *AK8PuppigenjetAUX;
   vector<float>   *AK8PuppimatchedDR;
   vector<float>   *AK8PuppijetRawFactor;
   TClonesArray    *AK8PuppiunCorrJetP4;
   vector<float>   *AK8PuppijetArea;
   vector<float>   *AK8PuppijetCorrUncUp;
   vector<float>   *AK8PuppijetCorrUncDown;
   vector<int>     *AK8PuppijetCharge;
   vector<int>     *AK8PuppijetPartonFlavor;
   vector<int>     *AK8PuppijetHadronFlavor;
   vector<bool>    *AK8PuppijetPassIDLoose;
   vector<bool>    *AK8PuppijetPassIDTight;
   vector<float>   *AK8PuppijetCEmEF;
   vector<float>   *AK8PuppijetCHadEF;
   vector<float>   *AK8PuppijetPhoEF;
   vector<float>   *AK8PuppijetNEmEF;
   vector<float>   *AK8PuppijetNHadEF;
   vector<float>   *AK8PuppijetEleEF;
   vector<float>   *AK8PuppijetMuoEF;
   vector<int>     *AK8PuppijetCMulti;
   vector<int>     *AK8PuppijetEleMulti;
   vector<int>     *AK8PuppijetMuoMulti;
   vector<float>   *AK8PuppijetSSV;
   vector<float>   *AK8PuppijetCSV;
   vector<float>   *AK8PuppijetSSVHE;
   vector<float>   *AK8PuppijetCISVV2;
   vector<float>   *AK8PuppijetTCHP;
   vector<float>   *AK8PuppijetTCHE;
   vector<float>   *AK8PuppijetJP;
   vector<float>   *AK8PuppijetJBP;
   vector<float>   *AK8PuppijetTau1;
   vector<float>   *AK8PuppijetTau2;
   vector<float>   *AK8PuppijetTau3;
   vector<float>   *AK8PuppijetTau21;
   vector<float>   *AK8PuppijetSDmass;
   vector<float>   *AK8PuppijetGenSDmass;
   vector<int>     *AK8PuppinSubSDJet;
   vector<vector<int> > *AK8PuppisubjetSDFatJetIndex;
   vector<vector<float> > *AK8PuppisubjetSDPx;
   vector<vector<float> > *AK8PuppisubjetSDPy;
   vector<vector<float> > *AK8PuppisubjetSDPz;
   vector<vector<float> > *AK8PuppisubjetSDE;
   vector<vector<float> > *AK8PuppisubjetSDRawFactor;
   vector<vector<int> > *AK8PuppisubjetSDPartonFlavor;
   vector<vector<int> > *AK8PuppisubjetSDHadronFlavor;
   vector<vector<float> > *AK8PuppisubjetSDCSV;
   vector<float>   *AK8Puppijet_DoubleSV;
   vector<int>     *AK8Puppijet_nSV;
   vector<vector<float> > *AK8Puppijet_SVMass;
   Int_t           CA15PuppinJet;
   TClonesArray    *CA15PuppijetP4;
   Float_t         CA15PuppijetRho;
   Int_t           CA15PuppijetNPV;
   TClonesArray    *CA15PuppigenjetP4;
   vector<float>   *CA15PuppigenjetEM;
   vector<float>   *CA15PuppigenjetHAD;
   vector<float>   *CA15PuppigenjetINV;
   vector<float>   *CA15PuppigenjetAUX;
   vector<float>   *CA15PuppimatchedDR;
   vector<float>   *CA15PuppijetRawFactor;
   TClonesArray    *CA15PuppiunCorrJetP4;
   vector<float>   *CA15PuppijetArea;
   vector<float>   *CA15PuppijetCorrUncUp;
   vector<float>   *CA15PuppijetCorrUncDown;
   vector<int>     *CA15PuppijetCharge;
   vector<int>     *CA15PuppijetPartonFlavor;
   vector<int>     *CA15PuppijetHadronFlavor;
   vector<bool>    *CA15PuppijetPassIDLoose;
   vector<bool>    *CA15PuppijetPassIDTight;
   vector<float>   *CA15PuppijetCEmEF;
   vector<float>   *CA15PuppijetCHadEF;
   vector<float>   *CA15PuppijetPhoEF;
   vector<float>   *CA15PuppijetNEmEF;
   vector<float>   *CA15PuppijetNHadEF;
   vector<float>   *CA15PuppijetEleEF;
   vector<float>   *CA15PuppijetMuoEF;
   vector<int>     *CA15PuppijetCMulti;
   vector<int>     *CA15PuppijetEleMulti;
   vector<int>     *CA15PuppijetMuoMulti;
   vector<float>   *CA15PuppijetSSV;
   vector<float>   *CA15PuppijetCSV;
   vector<float>   *CA15PuppijetSSVHE;
   vector<float>   *CA15PuppijetCISVV2;
   vector<float>   *CA15PuppijetTCHP;
   vector<float>   *CA15PuppijetTCHE;
   vector<float>   *CA15PuppijetJP;
   vector<float>   *CA15PuppijetJBP;
   vector<float>   *CA15PuppijetTau1;
   vector<float>   *CA15PuppijetTau2;
   vector<float>   *CA15PuppijetTau3;
   vector<float>   *CA15PuppijetTau21;
   vector<float>   *CA15PuppijetSDmass;
   vector<float>   *CA15PuppijetGenSDmass;
   vector<int>     *CA15PuppinSubSDJet;
   vector<vector<int> > *CA15PuppisubjetSDFatJetIndex;
   vector<vector<float> > *CA15PuppisubjetSDPx;
   vector<vector<float> > *CA15PuppisubjetSDPy;
   vector<vector<float> > *CA15PuppisubjetSDPz;
   vector<vector<float> > *CA15PuppisubjetSDE;
   vector<vector<float> > *CA15PuppisubjetSDRawFactor;
   vector<vector<int> > *CA15PuppisubjetSDPartonFlavor;
   vector<vector<int> > *CA15PuppisubjetSDHadronFlavor;
   vector<vector<float> > *CA15PuppisubjetSDCSV;
   vector<float>   *CA15Puppi_doublebtag;
   vector<float>   *CA15PuppiECF_2_3_10;
   vector<float>   *CA15PuppiECF_1_2_10;
   vector<float>   *CA15Puppijet_DoubleSV;
   vector<int>     *CA15Puppijet_nSV;
   vector<vector<float> > *CA15Puppijet_SVMass;

   // List of branches
   TBranch        *b_pu_nTrueInt;   //!
   TBranch        *b_pu_nPUVert;   //!
   TBranch        *b_isData;   //!
   TBranch        *b_eventId;   //!
   TBranch        *b_runId;   //!
   TBranch        *b_lumiSection;   //!
   TBranch        *b_bunchXing;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_vertexP3;   //!
   TBranch        *b_pfMetCorrPt;   //!
   TBranch        *b_pfMetCorrPhi;   //!
   TBranch        *b_pfMetCorrSumEt;   //!
   TBranch        *b_pfMetCorrSig;   //!
   TBranch        *b_pfMetCorrUnc;   //!
   TBranch        *b_pfMetRawPt;   //!
   TBranch        *b_pfMetRawPhi;   //!
   TBranch        *b_pfMetRawSumEt;   //!
   TBranch        *b_pfMetRawCov00;   //!
   TBranch        *b_pfMetRawCov01;   //!
   TBranch        *b_pfMetRawCov10;   //!
   TBranch        *b_pfMetRawCov11;   //!
   TBranch        *b_pfmvaMetPt;   //!
   TBranch        *b_pfmvaMetPhi;   //!
   TBranch        *b_pfmvaMetSumEt;   //!
   TBranch        *b_pfmvaMetSig;   //!
   TBranch        *b_pfpuppiMETPt;   //!
   TBranch        *b_pfpuppiMETPhi;   //!
   TBranch        *b_pfpuppiMETSumEt;   //!
   TBranch        *b_pfpuppiMETSig;   //!
   TBranch        *b_pfpuppiMETUnc;   //!
   TBranch        *b_hlt_nTrigs;   //!
   TBranch        *b_hlt_trigResult;   //!
   TBranch        *b_hlt_trigName;   //!
   TBranch        *b_hlt_trigPrescale;   //!
   TBranch        *b_hlt_nfilters;   //!
   TBranch        *b_hlt_filterResult;   //!
   TBranch        *b_hlt_hbhet;   //!
   TBranch        *b_hlt_filterName;   //!
   TBranch        *b_hlt_filterbadChCandidate;   //!
   TBranch        *b_hlt_filterbadPFMuon;   //!
   TBranch        *b_hlt_filterbadGlobalMuon;   //!
   TBranch        *b_hlt_filtercloneGlobalMuon;   //!
   TBranch        *b_genMET_true;   //!
   TBranch        *b_genMET_calo;   //!
   TBranch        *b_genMET_caloNonPrompt;   //!
   TBranch        *b_ptHat;   //!
   TBranch        *b_mcWeight;   //!
   TBranch        *b_HT;   //!
   TBranch        *b_pdf;   //!
   TBranch        *b_originalLHEweight;   //!
   TBranch        *b_pdfscaleSysWeights;   //!
   TBranch        *b_nGenPar;   //!
   TBranch        *b_genParP4;   //!
   TBranch        *b_genParQ;   //!
   TBranch        *b_genParId;   //!
   TBranch        *b_genParSt;   //!
   TBranch        *b_genMomParId;   //!
   TBranch        *b_genParIndex;   //!
   TBranch        *b_genNMo;   //!
   TBranch        *b_genNDa;   //!
   TBranch        *b_genMo1;   //!
   TBranch        *b_genMo2;   //!
   TBranch        *b_genDa1;   //!
   TBranch        *b_genDa2;   //!
   TBranch        *b_genStFlag;   //!
   TBranch        *b_ak4nGenJet;   //!
   TBranch        *b_ak4GenJetP4;   //!
   TBranch        *b_ak8nGenJet;   //!
   TBranch        *b_ak8GenJetP4;   //!
   TBranch        *b_eleRho;   //!
   TBranch        *b_nEle;   //!
   TBranch        *b_eleP4;   //!
   TBranch        *b_eleInBarrel;   //!
   TBranch        *b_eleInEndcap;   //!
   TBranch        *b_eleCharge;   //!
   TBranch        *b_eleChargeConsistent;   //!
   TBranch        *b_elecaloEnergy;   //!
   TBranch        *b_eleScEt;   //!
   TBranch        *b_eleScEn;   //!
   TBranch        *b_eleScPreEn;   //!
   TBranch        *b_eleScEta;   //!
   TBranch        *b_eleScPhi;   //!
   TBranch        *b_eleScRawEn;   //!
   TBranch        *b_eleScEtaWidth;   //!
   TBranch        *b_eleScPhiWidth;   //!
   TBranch        *b_eleR9;   //!
   TBranch        *b_eleHoverE;   //!
   TBranch        *b_eleD0;   //!
   TBranch        *b_eleDz;   //!
   TBranch        *b_eleEoverP;   //!
   TBranch        *b_eleBrem;   //!
   TBranch        *b_eledEtaAtVtx;   //!
   TBranch        *b_eledPhiAtVtx;   //!
   TBranch        *b_eleSigmaIEtaIEta;   //!
   TBranch        *b_eleSigmaIEtaIPhi;   //!
   TBranch        *b_eleSigmaIPhiIPhi;   //!
   TBranch        *b_eleConvVeto;   //!
   TBranch        *b_eleMissHits;   //!
   TBranch        *b_eleEoverPInv;   //!
   TBranch        *b_eleEtaseedAtVtx;   //!
   TBranch        *b_eleE1x5;   //!
   TBranch        *b_eleE2x5;   //!
   TBranch        *b_eleE5x5;   //!
   TBranch        *b_eleSigmaIEtaIEtaFull5x5;   //!
   TBranch        *b_eleE1x5Full5x5;   //!
   TBranch        *b_eleE2x5Full5x5;   //!
   TBranch        *b_eleE5x5Full5x5;   //!
   TBranch        *b_eleR9Full5x5;   //!
   TBranch        *b_eleChHadIso;   //!
   TBranch        *b_eleNeHadIso;   //!
   TBranch        *b_eleGamIso;   //!
   TBranch        *b_elePUPt;   //!
   TBranch        *b_eleEcalPFClusterIso;   //!
   TBranch        *b_eleHcalPFClusterIso;   //!
   TBranch        *b_eleMiniIso_ch;   //!
   TBranch        *b_eleMiniIso_nh;   //!
   TBranch        *b_eleMiniIso_ph;   //!
   TBranch        *b_eleMiniIso_pu;   //!
   TBranch        *b_eleMiniIso_r;   //!
   TBranch        *b_eleMiniIsoBeta;   //!
   TBranch        *b_eleMiniIsoEA;   //!
   TBranch        *b_eleEcalDrivenSeed;   //!
   TBranch        *b_eleEcalDriven;   //!
   TBranch        *b_eleDr03EcalRecHitSumEt;   //!
   TBranch        *b_eleDr03HcalDepth1TowerSumEt;   //!
   TBranch        *b_eleDr03HcalDepth2TowerSumEt;   //!
   TBranch        *b_eleDr03HcalTowerSumEt;   //!
   TBranch        *b_eleDr03TkSumPt;   //!
   TBranch        *b_eleIsPassVeto;   //!
   TBranch        *b_eleIsPassLoose;   //!
   TBranch        *b_eleIsPassMedium;   //!
   TBranch        *b_eleIsPassTight;   //!
   TBranch        *b_eleIsPassHEEP;   //!
   TBranch        *b_eleIsPassVetoNoIso;   //!
   TBranch        *b_eleIsPassLooseNoIso;   //!
   TBranch        *b_eleIsPassMediumNoIso;   //!
   TBranch        *b_eleIsPassTightNoIso;   //!
   TBranch        *b_eleIsPassHEEPNoIso;   //!
   TBranch        *b_eleIsPassMVAMedium;   //!
   TBranch        *b_eleIsPassMVATight;   //!
   TBranch        *b_eleMVAValue;   //!
   TBranch        *b_eleMVACategory;   //!
   TBranch        *b_nMu;   //!
   TBranch        *b_muP4;   //!
   TBranch        *b_muType;   //!
   TBranch        *b_muCharge;   //!
   TBranch        *b_isGlobalMuon;   //!
   TBranch        *b_isTrackerMuon;   //!
   TBranch        *b_isPFMuon;   //!
   TBranch        *b_isTightMuon;   //!
   TBranch        *b_isLooseMuon;   //!
   TBranch        *b_isMediumMuon;   //!
   TBranch        *b_isSoftMuon;   //!
   TBranch        *b_isHighPtMuon;   //!
   TBranch        *b_isCustomTrackerMuon;   //!
   TBranch        *b_muITrkID;   //!
   TBranch        *b_muSegID;   //!
   TBranch        *b_muNSegs;   //!
   TBranch        *b_muGood;   //!
   TBranch        *b_muIsGood;   //!
   TBranch        *b_muTrkPt;   //!
   TBranch        *b_muTrkPtErr;   //!
   TBranch        *b_mudxy;   //!
   TBranch        *b_mudz;   //!
   TBranch        *b_musegmentCompatibility;   //!
   TBranch        *b_muchi2LocalPosition;   //!
   TBranch        *b_mutrkKink;   //!
   TBranch        *b_muInnerdxy;   //!
   TBranch        *b_muInnerdz;   //!
   TBranch        *b_muTrkLayers;   //!
   TBranch        *b_muPixelLayers;   //!
   TBranch        *b_muPixelHits;   //!
   TBranch        *b_muTrkQuality;   //!
   TBranch        *b_muHits;   //!
   TBranch        *b_muChi2NDF;   //!
   TBranch        *b_muInnervalidFraction;   //!
   TBranch        *b_muMatches;   //!
   TBranch        *b_muTrkIso;   //!
   TBranch        *b_muHcalIso;   //!
   TBranch        *b_muEcalIso;   //!
   TBranch        *b_muChHadIso;   //!
   TBranch        *b_muNeHadIso;   //!
   TBranch        *b_muGamIso;   //!
   TBranch        *b_muPUPt;   //!
   TBranch        *b_muInnerTrkPt;   //!
   TBranch        *b_muMiniIso_ch;   //!
   TBranch        *b_muMiniIso_nh;   //!
   TBranch        *b_muMiniIso_ph;   //!
   TBranch        *b_muMiniIso_pu;   //!
   TBranch        *b_muMiniIso_r;   //!
   TBranch        *b_muMiniIsoBeta;   //!
   TBranch        *b_muMiniIsoEA;   //!
   TBranch        *b_HPSTau_n;   //!
   TBranch        *b_taupt;   //!
   TBranch        *b_HPSTau_4Momentum;   //!
   TBranch        *b_HPSTau_Vposition;   //!
   TBranch        *b_HPSTau_leadPFChargedHadrCand;   //!
   TBranch        *b_HPSTau_leadPFChargedHadrCand_trackRef;   //!
   TBranch        *b_disc_againstElectronLoose;   //!
   TBranch        *b_disc_againstElectronMedium;   //!
   TBranch        *b_disc_againstElectronTight;   //!
   TBranch        *b_disc_againstElectronLooseMVA5;   //!
   TBranch        *b_disc_againstElectronMediumMVA5;   //!
   TBranch        *b_disc_againstElectronTightMVA5;   //!
   TBranch        *b_disc_againstElectronVLooseMVA5;   //!
   TBranch        *b_disc_againstElectronVTightMVA5;   //!
   TBranch        *b_disc_againstMuonLoose;   //!
   TBranch        *b_disc_againstMuonMedium;   //!
   TBranch        *b_disc_againstMuonTight;   //!
   TBranch        *b_disc_againstMuonLoose2;   //!
   TBranch        *b_disc_againstMuonMedium2;   //!
   TBranch        *b_disc_againstMuonTight2;   //!
   TBranch        *b_disc_againstMuonLooseMVA;   //!
   TBranch        *b_disc_againstMuonMediumMVA;   //!
   TBranch        *b_disc_againstMuonTightMVA;   //!
   TBranch        *b_disc_againstMuonLoose3;   //!
   TBranch        *b_disc_againstMuonTight3;   //!
   TBranch        *b_disc_byVLooseCombinedIsolationDeltaBetaCorr;   //!
   TBranch        *b_disc_byLooseCombinedIsolationDeltaBetaCorr;   //!
   TBranch        *b_disc_byMediumCombinedIsolationDeltaBetaCorr;   //!
   TBranch        *b_disc_byTightCombinedIsolationDeltaBetaCorr;   //!
   TBranch        *b_disc_byLooseIsolation;   //!
   TBranch        *b_disc_byVLooseIsolationMVA3newDMwLT;   //!
   TBranch        *b_disc_byLooseIsolationMVA3newDMwLT;   //!
   TBranch        *b_disc_byMediumIsolationMVA3newDMwLT;   //!
   TBranch        *b_disc_byTightIsolationMVA3newDMwLT;   //!
   TBranch        *b_disc_byVTightIsolationMVA3newDMwLT;   //!
   TBranch        *b_disc_byVVTightIsolationMVA3newDMwLT;   //!
   TBranch        *b_disc_byVLooseIsolationMVA3newDMwoLT;   //!
   TBranch        *b_disc_byLooseIsolationMVA3newDMwoLT;   //!
   TBranch        *b_disc_byMediumIsolationMVA3newDMwoLT;   //!
   TBranch        *b_disc_byTightIsolationMVA3newDMwoLT;   //!
   TBranch        *b_disc_byVTightIsolationMVA3newDMwoLT;   //!
   TBranch        *b_disc_byVVTightIsolationMVA3newDMwoLT;   //!
   TBranch        *b_disc_byVLooseIsolationMVA3oldDMwLT;   //!
   TBranch        *b_disc_byLooseIsolationMVA3oldDMwLT;   //!
   TBranch        *b_disc_byMediumIsolationMVA3oldDMwLT;   //!
   TBranch        *b_disc_byTightIsolationMVA3oldDMwLT;   //!
   TBranch        *b_disc_byVTightIsolationMVA3oldDMwLT;   //!
   TBranch        *b_disc_byVVTightIsolationMVA3oldDMwLT;   //!
   TBranch        *b_disc_byVLooseIsolationMVA3oldDMwoLT;   //!
   TBranch        *b_disc_byLooseIsolationMVA3oldDMwoLT;   //!
   TBranch        *b_disc_byMediumIsolationMVA3oldDMwoLT;   //!
   TBranch        *b_disc_byTightIsolationMVA3oldDMwoLT;   //!
   TBranch        *b_disc_byVTightIsolationMVA3oldDMwoLT;   //!
   TBranch        *b_disc_byVVTightIsolationMVA3oldDMwoLT;   //!
   TBranch        *b_disc_byLooseCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_disc_byMediumCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_disc_byTightCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_disc_decayModeFinding;   //!
   TBranch        *b_disc_decayModeFindingNewDMs;   //!
   TBranch        *b_disc_chargedIsoPtSum;   //!
   TBranch        *b_disc_neutralIsoPtSum;   //!
   TBranch        *b_disc_puCorrPtSum;   //!
   TBranch        *b_HPSTau_NewVz;   //!
   TBranch        *b_HPSTau_charge;   //!
   TBranch        *b_nPho;   //!
   TBranch        *b_phoP4;   //!
   TBranch        *b_phoIsPassTight;   //!
   TBranch        *b_phoIsPassLoose;   //!
   TBranch        *b_phoIsPassMedium;   //!
   TBranch        *b_phoIDMVA;   //!
   TBranch        *b_phoSCE;   //!
   TBranch        *b_phoSCRawE;   //!
   TBranch        *b_phoSCEta;   //!
   TBranch        *b_phoSCPhi;   //!
   TBranch        *b_phoSCEtaWidth;   //!
   TBranch        *b_phoSCPhiWidth;   //!
   TBranch        *b_phoSCBrem;   //!
   TBranch        *b_phohasPixelSeed;   //!
   TBranch        *b_phoEleVeto;   //!
   TBranch        *b_phoR9;   //!
   TBranch        *b_phoHoverE;   //!
   TBranch        *b_phoSigmaIEtaIEta;   //!
   TBranch        *b_phoSigmaIEtaIPhi;   //!
   TBranch        *b_phoSigmaIPhiIPhi;   //!
   TBranch        *b_phoSigmaIEtaIEtaFull5x5;   //!
   TBranch        *b_phoR9Full5x5;   //!
   TBranch        *b_phoPFChIso;   //!
   TBranch        *b_phoPFPhoIso;   //!
   TBranch        *b_phoPFNeuIso;   //!
   TBranch        *b_THINnJet;   //!
   TBranch        *b_THINjetP4;   //!
   TBranch        *b_THINjetRho;   //!
   TBranch        *b_THINjetNPV;   //!
   TBranch        *b_THINgenjetP4;   //!
   TBranch        *b_THINgenjetEM;   //!
   TBranch        *b_THINgenjetHAD;   //!
   TBranch        *b_THINgenjetINV;   //!
   TBranch        *b_THINgenjetAUX;   //!
   TBranch        *b_THINmatchedDR;   //!
   TBranch        *b_THINjetRawFactor;   //!
   TBranch        *b_THINunCorrJetP4;   //!
   TBranch        *b_THINjetArea;   //!
   TBranch        *b_THINjetCorrUncUp;   //!
   TBranch        *b_THINjetCorrUncDown;   //!
   TBranch        *b_THINjetCharge;   //!
   TBranch        *b_THINjetPartonFlavor;   //!
   TBranch        *b_THINjetHadronFlavor;   //!
   TBranch        *b_THINjetPassIDLoose;   //!
   TBranch        *b_THINjetPassIDTight;   //!
   TBranch        *b_THINjetCEmEF;   //!
   TBranch        *b_THINjetCHadEF;   //!
   TBranch        *b_THINjetPhoEF;   //!
   TBranch        *b_THINjetNEmEF;   //!
   TBranch        *b_THINjetNHadEF;   //!
   TBranch        *b_THINjetEleEF;   //!
   TBranch        *b_THINjetMuoEF;   //!
   TBranch        *b_THINjetCMulti;   //!
   TBranch        *b_THINjetEleMulti;   //!
   TBranch        *b_THINjetMuoMulti;   //!
   TBranch        *b_THINjetSSV;   //!
   TBranch        *b_THINjetCSV;   //!
   TBranch        *b_THINjetSSVHE;   //!
   TBranch        *b_THINjetCISVV2;   //!
   TBranch        *b_THINjetTCHP;   //!
   TBranch        *b_THINjetTCHE;   //!
   TBranch        *b_THINjetJP;   //!
   TBranch        *b_THINjetJBP;   //!
   TBranch        *b_THINPUJetID;   //!
   TBranch        *b_THINisPUJetIDLoose;   //!
   TBranch        *b_THINisPUJetIDMedium;   //!
   TBranch        *b_THINisPUJetIDTight;   //!
   TBranch        *b_FATnJet;   //!
   TBranch        *b_FATjetP4;   //!
   TBranch        *b_FATjetRho;   //!
   TBranch        *b_FATjetNPV;   //!
   TBranch        *b_FATgenjetP4;   //!
   TBranch        *b_FATgenjetEM;   //!
   TBranch        *b_FATgenjetHAD;   //!
   TBranch        *b_FATgenjetINV;   //!
   TBranch        *b_FATgenjetAUX;   //!
   TBranch        *b_FATmatchedDR;   //!
   TBranch        *b_FATjetRawFactor;   //!
   TBranch        *b_FATunCorrJetP4;   //!
   TBranch        *b_FATjetArea;   //!
   TBranch        *b_FATjetCorrUncUp;   //!
   TBranch        *b_FATjetCorrUncDown;   //!
   TBranch        *b_FATjetCharge;   //!
   TBranch        *b_FATjetPartonFlavor;   //!
   TBranch        *b_FATjetHadronFlavor;   //!
   TBranch        *b_FATjetPassIDLoose;   //!
   TBranch        *b_FATjetPassIDTight;   //!
   TBranch        *b_FATjetCEmEF;   //!
   TBranch        *b_FATjetCHadEF;   //!
   TBranch        *b_FATjetPhoEF;   //!
   TBranch        *b_FATjetNEmEF;   //!
   TBranch        *b_FATjetNHadEF;   //!
   TBranch        *b_FATjetEleEF;   //!
   TBranch        *b_FATjetMuoEF;   //!
   TBranch        *b_FATjetCMulti;   //!
   TBranch        *b_FATjetEleMulti;   //!
   TBranch        *b_FATjetMuoMulti;   //!
   TBranch        *b_FATjetSSV;   //!
   TBranch        *b_FATjetCSV;   //!
   TBranch        *b_FATjetSSVHE;   //!
   TBranch        *b_FATjetCISVV2;   //!
   TBranch        *b_FATjetTCHP;   //!
   TBranch        *b_FATjetTCHE;   //!
   TBranch        *b_FATjetJP;   //!
   TBranch        *b_FATjetJBP;   //!
   TBranch        *b_FATjetTau1;   //!
   TBranch        *b_FATjetTau2;   //!
   TBranch        *b_FATjetTau3;   //!
   TBranch        *b_FATjetTau21;   //!
   TBranch        *b_FATjetSDmass;   //!
   TBranch        *b_FATjetGenSDmass;   //!
   TBranch        *b_FATnSubSDJet;   //!
   TBranch        *b_FATsubjetSDFatJetIndex;   //!
   TBranch        *b_FATsubjetSDPx;   //!
   TBranch        *b_FATsubjetSDPy;   //!
   TBranch        *b_FATsubjetSDPz;   //!
   TBranch        *b_FATsubjetSDE;   //!
   TBranch        *b_FATsubjetSDRawFactor;   //!
   TBranch        *b_FATsubjetSDPartonFlavor;   //!
   TBranch        *b_FATsubjetSDHadronFlavor;   //!
   TBranch        *b_FATsubjetSDCSV;   //!
   TBranch        *b_FATjetPRmass;   //!
   TBranch        *b_FATjetPRmassL2L3Corr;   //!
   TBranch        *b_FATjetPuppiTau1;   //!
   TBranch        *b_FATjetPuppiTau2;   //!
   TBranch        *b_FATjetPuppiTau3;   //!
   TBranch        *b_FATjetPuppiSDmass;   //!
   TBranch        *b_FATjetPuppiP4;   //!
   TBranch        *b_FATjetPuppiSDRawP4;   //!
   TBranch        *b_FATnSubSDPuppiJet;   //!
   TBranch        *b_FATsubjetSDPuppiFatJetIndex;   //!
   TBranch        *b_FATsubjetSDPuppiPx;   //!
   TBranch        *b_FATsubjetSDPuppiPy;   //!
   TBranch        *b_FATsubjetSDPuppiPz;   //!
   TBranch        *b_FATsubjetSDPuppiE;   //!
   TBranch        *b_FATsubjetSDPuppiCSV;   //!
   TBranch        *b_FATjet_DoubleSV;   //!
   TBranch        *b_FATjet_nSV;   //!
   TBranch        *b_FATjet_SVMass;   //!
   TBranch        *b_ADDnJet;   //!
   TBranch        *b_ADDjetP4;   //!
   TBranch        *b_ADDjet_DoubleSV;   //!
   TBranch        *b_ADDjet_nSV;   //!
   TBranch        *b_ADDjet_SVMass;   //!
   TBranch        *b_AK8PuppinJet;   //!
   TBranch        *b_AK8PuppijetP4;   //!
   TBranch        *b_AK8PuppijetRho;   //!
   TBranch        *b_AK8PuppijetNPV;   //!
   TBranch        *b_AK8PuppigenjetP4;   //!
   TBranch        *b_AK8PuppigenjetEM;   //!
   TBranch        *b_AK8PuppigenjetHAD;   //!
   TBranch        *b_AK8PuppigenjetINV;   //!
   TBranch        *b_AK8PuppigenjetAUX;   //!
   TBranch        *b_AK8PuppimatchedDR;   //!
   TBranch        *b_AK8PuppijetRawFactor;   //!
   TBranch        *b_AK8PuppiunCorrJetP4;   //!
   TBranch        *b_AK8PuppijetArea;   //!
   TBranch        *b_AK8PuppijetCorrUncUp;   //!
   TBranch        *b_AK8PuppijetCorrUncDown;   //!
   TBranch        *b_AK8PuppijetCharge;   //!
   TBranch        *b_AK8PuppijetPartonFlavor;   //!
   TBranch        *b_AK8PuppijetHadronFlavor;   //!
   TBranch        *b_AK8PuppijetPassIDLoose;   //!
   TBranch        *b_AK8PuppijetPassIDTight;   //!
   TBranch        *b_AK8PuppijetCEmEF;   //!
   TBranch        *b_AK8PuppijetCHadEF;   //!
   TBranch        *b_AK8PuppijetPhoEF;   //!
   TBranch        *b_AK8PuppijetNEmEF;   //!
   TBranch        *b_AK8PuppijetNHadEF;   //!
   TBranch        *b_AK8PuppijetEleEF;   //!
   TBranch        *b_AK8PuppijetMuoEF;   //!
   TBranch        *b_AK8PuppijetCMulti;   //!
   TBranch        *b_AK8PuppijetEleMulti;   //!
   TBranch        *b_AK8PuppijetMuoMulti;   //!
   TBranch        *b_AK8PuppijetSSV;   //!
   TBranch        *b_AK8PuppijetCSV;   //!
   TBranch        *b_AK8PuppijetSSVHE;   //!
   TBranch        *b_AK8PuppijetCISVV2;   //!
   TBranch        *b_AK8PuppijetTCHP;   //!
   TBranch        *b_AK8PuppijetTCHE;   //!
   TBranch        *b_AK8PuppijetJP;   //!
   TBranch        *b_AK8PuppijetJBP;   //!
   TBranch        *b_AK8PuppijetTau1;   //!
   TBranch        *b_AK8PuppijetTau2;   //!
   TBranch        *b_AK8PuppijetTau3;   //!
   TBranch        *b_AK8PuppijetTau21;   //!
   TBranch        *b_AK8PuppijetSDmass;   //!
   TBranch        *b_AK8PuppijetGenSDmass;   //!
   TBranch        *b_AK8PuppinSubSDJet;   //!
   TBranch        *b_AK8PuppisubjetSDFatJetIndex;   //!
   TBranch        *b_AK8PuppisubjetSDPx;   //!
   TBranch        *b_AK8PuppisubjetSDPy;   //!
   TBranch        *b_AK8PuppisubjetSDPz;   //!
   TBranch        *b_AK8PuppisubjetSDE;   //!
   TBranch        *b_AK8PuppisubjetSDRawFactor;   //!
   TBranch        *b_AK8PuppisubjetSDPartonFlavor;   //!
   TBranch        *b_AK8PuppisubjetSDHadronFlavor;   //!
   TBranch        *b_AK8PuppisubjetSDCSV;   //!
   TBranch        *b_AK8Puppijet_DoubleSV;   //!
   TBranch        *b_AK8Puppijet_nSV;   //!
   TBranch        *b_AK8Puppijet_SVMass;   //!
   TBranch        *b_CA15PuppinJet;   //!
   TBranch        *b_CA15PuppijetP4;   //!
   TBranch        *b_CA15PuppijetRho;   //!
   TBranch        *b_CA15PuppijetNPV;   //!
   TBranch        *b_CA15PuppigenjetP4;   //!
   TBranch        *b_CA15PuppigenjetEM;   //!
   TBranch        *b_CA15PuppigenjetHAD;   //!
   TBranch        *b_CA15PuppigenjetINV;   //!
   TBranch        *b_CA15PuppigenjetAUX;   //!
   TBranch        *b_CA15PuppimatchedDR;   //!
   TBranch        *b_CA15PuppijetRawFactor;   //!
   TBranch        *b_CA15PuppiunCorrJetP4;   //!
   TBranch        *b_CA15PuppijetArea;   //!
   TBranch        *b_CA15PuppijetCorrUncUp;   //!
   TBranch        *b_CA15PuppijetCorrUncDown;   //!
   TBranch        *b_CA15PuppijetCharge;   //!
   TBranch        *b_CA15PuppijetPartonFlavor;   //!
   TBranch        *b_CA15PuppijetHadronFlavor;   //!
   TBranch        *b_CA15PuppijetPassIDLoose;   //!
   TBranch        *b_CA15PuppijetPassIDTight;   //!
   TBranch        *b_CA15PuppijetCEmEF;   //!
   TBranch        *b_CA15PuppijetCHadEF;   //!
   TBranch        *b_CA15PuppijetPhoEF;   //!
   TBranch        *b_CA15PuppijetNEmEF;   //!
   TBranch        *b_CA15PuppijetNHadEF;   //!
   TBranch        *b_CA15PuppijetEleEF;   //!
   TBranch        *b_CA15PuppijetMuoEF;   //!
   TBranch        *b_CA15PuppijetCMulti;   //!
   TBranch        *b_CA15PuppijetEleMulti;   //!
   TBranch        *b_CA15PuppijetMuoMulti;   //!
   TBranch        *b_CA15PuppijetSSV;   //!
   TBranch        *b_CA15PuppijetCSV;   //!
   TBranch        *b_CA15PuppijetSSVHE;   //!
   TBranch        *b_CA15PuppijetCISVV2;   //!
   TBranch        *b_CA15PuppijetTCHP;   //!
   TBranch        *b_CA15PuppijetTCHE;   //!
   TBranch        *b_CA15PuppijetJP;   //!
   TBranch        *b_CA15PuppijetJBP;   //!
   TBranch        *b_CA15PuppijetTau1;   //!
   TBranch        *b_CA15PuppijetTau2;   //!
   TBranch        *b_CA15PuppijetTau3;   //!
   TBranch        *b_CA15PuppijetTau21;   //!
   TBranch        *b_CA15PuppijetSDmass;   //!
   TBranch        *b_CA15PuppijetGenSDmass;   //!
   TBranch        *b_CA15PuppinSubSDJet;   //!
   TBranch        *b_CA15PuppisubjetSDFatJetIndex;   //!
   TBranch        *b_CA15PuppisubjetSDPx;   //!
   TBranch        *b_CA15PuppisubjetSDPy;   //!
   TBranch        *b_CA15PuppisubjetSDPz;   //!
   TBranch        *b_CA15PuppisubjetSDE;   //!
   TBranch        *b_CA15PuppisubjetSDRawFactor;   //!
   TBranch        *b_CA15PuppisubjetSDPartonFlavor;   //!
   TBranch        *b_CA15PuppisubjetSDHadronFlavor;   //!
   TBranch        *b_CA15PuppisubjetSDCSV;   //!
   TBranch        *b_CA15Puppi_doublebtag;   //!
   TBranch        *b_CA15PuppiECF_2_3_10;   //!
   TBranch        *b_CA15PuppiECF_1_2_10;   //!
   TBranch        *b_CA15Puppijet_DoubleSV;   //!
   TBranch        *b_CA15Puppijet_nSV;   //!
   TBranch        *b_CA15Puppijet_SVMass;   //!

   treeMaker();

   treeMaker(string inputFile,string outputFile, TTree *tree=0);
   virtual ~treeMaker();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(string outputFile);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

treeMaker::treeMaker() : fChain(0)
{

}

#ifdef treeMaker_cxx
treeMaker::treeMaker(string inputFile,string outputFile, TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(inputFile.data());
      if (!f || !f->IsOpen()) {
         f = new TFile(inputFile.data());
      }
      TDirectory * dir = (TDirectory*)f->Get(Form("%s:/tree",inputFile.data()));
      dir->GetObject("treeMaker",tree);

   }
   Init(tree);
   cout <<"load File"<<endl;
//   Loop(outputFile);
}

treeMaker::~treeMaker()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t treeMaker::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t treeMaker::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void treeMaker::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   vertexP3 = 0;
   pfMetCorrUnc = 0;
   pfpuppiMETUnc = 0;
   hlt_trigResult = 0;
   hlt_trigName = 0;
   hlt_trigPrescale = 0;
   hlt_filterResult = 0;
   hlt_filterName = 0;
   pdf = 0;
   pdfscaleSysWeights = 0;
   genParP4 = 0;
   genParQ = 0;
   genParId = 0;
   genParSt = 0;
   genMomParId = 0;
   genParIndex = 0;
   genNMo = 0;
   genNDa = 0;
   genMo1 = 0;
   genMo2 = 0;
   genDa1 = 0;
   genDa2 = 0;
   genStFlag = 0;
   ak4GenJetP4 = 0;
   ak8GenJetP4 = 0;
   eleP4 = 0;
   eleInBarrel = 0;
   eleInEndcap = 0;
   eleCharge = 0;
   eleChargeConsistent = 0;
   elecaloEnergy = 0;
   eleScEt = 0;
   eleScEn = 0;
   eleScPreEn = 0;
   eleScEta = 0;
   eleScPhi = 0;
   eleScRawEn = 0;
   eleScEtaWidth = 0;
   eleScPhiWidth = 0;
   eleR9 = 0;
   eleHoverE = 0;
   eleD0 = 0;
   eleDz = 0;
   eleEoverP = 0;
   eleBrem = 0;
   eledEtaAtVtx = 0;
   eledPhiAtVtx = 0;
   eleSigmaIEtaIEta = 0;
   eleSigmaIEtaIPhi = 0;
   eleSigmaIPhiIPhi = 0;
   eleConvVeto = 0;
   eleMissHits = 0;
   eleEoverPInv = 0;
   eleEtaseedAtVtx = 0;
   eleE1x5 = 0;
   eleE2x5 = 0;
   eleE5x5 = 0;
   eleSigmaIEtaIEtaFull5x5 = 0;
   eleE1x5Full5x5 = 0;
   eleE2x5Full5x5 = 0;
   eleE5x5Full5x5 = 0;
   eleR9Full5x5 = 0;
   eleChHadIso = 0;
   eleNeHadIso = 0;
   eleGamIso = 0;
   elePUPt = 0;
   eleEcalPFClusterIso = 0;
   eleHcalPFClusterIso = 0;
   eleMiniIso_ch = 0;
   eleMiniIso_nh = 0;
   eleMiniIso_ph = 0;
   eleMiniIso_pu = 0;
   eleMiniIso_r = 0;
   eleMiniIsoBeta = 0;
   eleMiniIsoEA = 0;
   eleEcalDrivenSeed = 0;
   eleEcalDriven = 0;
   eleDr03EcalRecHitSumEt = 0;
   eleDr03HcalDepth1TowerSumEt = 0;
   eleDr03HcalDepth2TowerSumEt = 0;
   eleDr03HcalTowerSumEt = 0;
   eleDr03TkSumPt = 0;
   eleIsPassVeto = 0;
   eleIsPassLoose = 0;
   eleIsPassMedium = 0;
   eleIsPassTight = 0;
   eleIsPassHEEP = 0;
   eleIsPassVetoNoIso = 0;
   eleIsPassLooseNoIso = 0;
   eleIsPassMediumNoIso = 0;
   eleIsPassTightNoIso = 0;
   eleIsPassHEEPNoIso = 0;
   eleIsPassMVAMedium = 0;
   eleIsPassMVATight = 0;
   eleMVAValue = 0;
   eleMVACategory = 0;
   muP4 = 0;
   muType = 0;
   muCharge = 0;
   isGlobalMuon = 0;
   isTrackerMuon = 0;
   isPFMuon = 0;
   isTightMuon = 0;
   isLooseMuon = 0;
   isMediumMuon = 0;
   isSoftMuon = 0;
   isHighPtMuon = 0;
   isCustomTrackerMuon = 0;
   muITrkID = 0;
   muSegID = 0;
   muNSegs = 0;
   muGood = 0;
   muIsGood = 0;
   muTrkPt = 0;
   muTrkPtErr = 0;
   mudxy = 0;
   mudz = 0;
   musegmentCompatibility = 0;
   muchi2LocalPosition = 0;
   mutrkKink = 0;
   muInnerdxy = 0;
   muInnerdz = 0;
   muTrkLayers = 0;
   muPixelLayers = 0;
   muPixelHits = 0;
   muTrkQuality = 0;
   muHits = 0;
   muChi2NDF = 0;
   muInnervalidFraction = 0;
   muMatches = 0;
   muTrkIso = 0;
   muHcalIso = 0;
   muEcalIso = 0;
   muChHadIso = 0;
   muNeHadIso = 0;
   muGamIso = 0;
   muPUPt = 0;
   muInnerTrkPt = 0;
   muMiniIso_ch = 0;
   muMiniIso_nh = 0;
   muMiniIso_ph = 0;
   muMiniIso_pu = 0;
   muMiniIso_r = 0;
   muMiniIsoBeta = 0;
   muMiniIsoEA = 0;
   taupt = 0;
   HPSTau_4Momentum = 0;
   HPSTau_Vposition = 0;
   HPSTau_leadPFChargedHadrCand = 0;
   HPSTau_leadPFChargedHadrCand_trackRef = 0;
   disc_againstElectronLoose = 0;
   disc_againstElectronMedium = 0;
   disc_againstElectronTight = 0;
   disc_againstElectronLooseMVA5 = 0;
   disc_againstElectronMediumMVA5 = 0;
   disc_againstElectronTightMVA5 = 0;
   disc_againstElectronVLooseMVA5 = 0;
   disc_againstElectronVTightMVA5 = 0;
   disc_againstMuonLoose = 0;
   disc_againstMuonMedium = 0;
   disc_againstMuonTight = 0;
   disc_againstMuonLoose2 = 0;
   disc_againstMuonMedium2 = 0;
   disc_againstMuonTight2 = 0;
   disc_againstMuonLooseMVA = 0;
   disc_againstMuonMediumMVA = 0;
   disc_againstMuonTightMVA = 0;
   disc_againstMuonLoose3 = 0;
   disc_againstMuonTight3 = 0;
   disc_byVLooseCombinedIsolationDeltaBetaCorr = 0;
   disc_byLooseCombinedIsolationDeltaBetaCorr = 0;
   disc_byMediumCombinedIsolationDeltaBetaCorr = 0;
   disc_byTightCombinedIsolationDeltaBetaCorr = 0;
   disc_byLooseIsolation = 0;
   disc_byVLooseIsolationMVA3newDMwLT = 0;
   disc_byLooseIsolationMVA3newDMwLT = 0;
   disc_byMediumIsolationMVA3newDMwLT = 0;
   disc_byTightIsolationMVA3newDMwLT = 0;
   disc_byVTightIsolationMVA3newDMwLT = 0;
   disc_byVVTightIsolationMVA3newDMwLT = 0;
   disc_byVLooseIsolationMVA3newDMwoLT = 0;
   disc_byLooseIsolationMVA3newDMwoLT = 0;
   disc_byMediumIsolationMVA3newDMwoLT = 0;
   disc_byTightIsolationMVA3newDMwoLT = 0;
   disc_byVTightIsolationMVA3newDMwoLT = 0;
   disc_byVVTightIsolationMVA3newDMwoLT = 0;
   disc_byVLooseIsolationMVA3oldDMwLT = 0;
   disc_byLooseIsolationMVA3oldDMwLT = 0;
   disc_byMediumIsolationMVA3oldDMwLT = 0;
   disc_byTightIsolationMVA3oldDMwLT = 0;
   disc_byVTightIsolationMVA3oldDMwLT = 0;
   disc_byVVTightIsolationMVA3oldDMwLT = 0;
   disc_byVLooseIsolationMVA3oldDMwoLT = 0;
   disc_byLooseIsolationMVA3oldDMwoLT = 0;
   disc_byMediumIsolationMVA3oldDMwoLT = 0;
   disc_byTightIsolationMVA3oldDMwoLT = 0;
   disc_byVTightIsolationMVA3oldDMwoLT = 0;
   disc_byVVTightIsolationMVA3oldDMwoLT = 0;
   disc_byLooseCombinedIsolationDeltaBetaCorr3Hits = 0;
   disc_byMediumCombinedIsolationDeltaBetaCorr3Hits = 0;
   disc_byTightCombinedIsolationDeltaBetaCorr3Hits = 0;
   disc_decayModeFinding = 0;
   disc_decayModeFindingNewDMs = 0;
   disc_chargedIsoPtSum = 0;
   disc_neutralIsoPtSum = 0;
   disc_puCorrPtSum = 0;
   HPSTau_NewVz = 0;
   HPSTau_charge = 0;
   phoP4 = 0;
   phoIsPassTight = 0;
   phoIsPassLoose = 0;
   phoIsPassMedium = 0;
   phoIDMVA = 0;
   phoSCE = 0;
   phoSCRawE = 0;
   phoSCEta = 0;
   phoSCPhi = 0;
   phoSCEtaWidth = 0;
   phoSCPhiWidth = 0;
   phoSCBrem = 0;
   phohasPixelSeed = 0;
   phoEleVeto = 0;
   phoR9 = 0;
   phoHoverE = 0;
   phoSigmaIEtaIEta = 0;
   phoSigmaIEtaIPhi = 0;
   phoSigmaIPhiIPhi = 0;
   phoSigmaIEtaIEtaFull5x5 = 0;
   phoR9Full5x5 = 0;
   phoPFChIso = 0;
   phoPFPhoIso = 0;
   phoPFNeuIso = 0;
   THINjetP4 = 0;
   THINgenjetP4 = 0;
   THINgenjetEM = 0;
   THINgenjetHAD = 0;
   THINgenjetINV = 0;
   THINgenjetAUX = 0;
   THINmatchedDR = 0;
   THINjetRawFactor = 0;
   THINunCorrJetP4 = 0;
   THINjetArea = 0;
   THINjetCorrUncUp = 0;
   THINjetCorrUncDown = 0;
   THINjetCharge = 0;
   THINjetPartonFlavor = 0;
   THINjetHadronFlavor = 0;
   THINjetPassIDLoose = 0;
   THINjetPassIDTight = 0;
   THINjetCEmEF = 0;
   THINjetCHadEF = 0;
   THINjetPhoEF = 0;
   THINjetNEmEF = 0;
   THINjetNHadEF = 0;
   THINjetEleEF = 0;
   THINjetMuoEF = 0;
   THINjetCMulti = 0;
   THINjetEleMulti = 0;
   THINjetMuoMulti = 0;
   THINjetSSV = 0;
   THINjetCSV = 0;
   THINjetSSVHE = 0;
   THINjetCISVV2 = 0;
   THINjetTCHP = 0;
   THINjetTCHE = 0;
   THINjetJP = 0;
   THINjetJBP = 0;
   THINPUJetID = 0;
   THINisPUJetIDLoose = 0;
   THINisPUJetIDMedium = 0;
   THINisPUJetIDTight = 0;
   FATjetP4 = 0;
   FATgenjetP4 = 0;
   FATgenjetEM = 0;
   FATgenjetHAD = 0;
   FATgenjetINV = 0;
   FATgenjetAUX = 0;
   FATmatchedDR = 0;
   FATjetRawFactor = 0;
   FATunCorrJetP4 = 0;
   FATjetArea = 0;
   FATjetCorrUncUp = 0;
   FATjetCorrUncDown = 0;
   FATjetCharge = 0;
   FATjetPartonFlavor = 0;
   FATjetHadronFlavor = 0;
   FATjetPassIDLoose = 0;
   FATjetPassIDTight = 0;
   FATjetCEmEF = 0;
   FATjetCHadEF = 0;
   FATjetPhoEF = 0;
   FATjetNEmEF = 0;
   FATjetNHadEF = 0;
   FATjetEleEF = 0;
   FATjetMuoEF = 0;
   FATjetCMulti = 0;
   FATjetEleMulti = 0;
   FATjetMuoMulti = 0;
   FATjetSSV = 0;
   FATjetCSV = 0;
   FATjetSSVHE = 0;
   FATjetCISVV2 = 0;
   FATjetTCHP = 0;
   FATjetTCHE = 0;
   FATjetJP = 0;
   FATjetJBP = 0;
   FATjetTau1 = 0;
   FATjetTau2 = 0;
   FATjetTau3 = 0;
   FATjetTau21 = 0;
   FATjetSDmass = 0;
   FATjetGenSDmass = 0;
   FATnSubSDJet = 0;
   FATsubjetSDFatJetIndex = 0;
   FATsubjetSDPx = 0;
   FATsubjetSDPy = 0;
   FATsubjetSDPz = 0;
   FATsubjetSDE = 0;
   FATsubjetSDRawFactor = 0;
   FATsubjetSDPartonFlavor = 0;
   FATsubjetSDHadronFlavor = 0;
   FATsubjetSDCSV = 0;
   FATjetPRmass = 0;
   FATjetPRmassL2L3Corr = 0;
   FATjetPuppiTau1 = 0;
   FATjetPuppiTau2 = 0;
   FATjetPuppiTau3 = 0;
   FATjetPuppiSDmass = 0;
   FATjetPuppiP4 = 0;
   FATjetPuppiSDRawP4 = 0;
   FATnSubSDPuppiJet = 0;
   FATsubjetSDPuppiFatJetIndex = 0;
   FATsubjetSDPuppiPx = 0;
   FATsubjetSDPuppiPy = 0;
   FATsubjetSDPuppiPz = 0;
   FATsubjetSDPuppiE = 0;
   FATsubjetSDPuppiCSV = 0;
   FATjet_DoubleSV = 0;
   FATjet_nSV = 0;
   FATjet_SVMass = 0;
   ADDjetP4 = 0;
   ADDjet_DoubleSV = 0;
   ADDjet_nSV = 0;
   ADDjet_SVMass = 0;
   AK8PuppijetP4 = 0;
   AK8PuppigenjetP4 = 0;
   AK8PuppigenjetEM = 0;
   AK8PuppigenjetHAD = 0;
   AK8PuppigenjetINV = 0;
   AK8PuppigenjetAUX = 0;
   AK8PuppimatchedDR = 0;
   AK8PuppijetRawFactor = 0;
   AK8PuppiunCorrJetP4 = 0;
   AK8PuppijetArea = 0;
   AK8PuppijetCorrUncUp = 0;
   AK8PuppijetCorrUncDown = 0;
   AK8PuppijetCharge = 0;
   AK8PuppijetPartonFlavor = 0;
   AK8PuppijetHadronFlavor = 0;
   AK8PuppijetPassIDLoose = 0;
   AK8PuppijetPassIDTight = 0;
   AK8PuppijetCEmEF = 0;
   AK8PuppijetCHadEF = 0;
   AK8PuppijetPhoEF = 0;
   AK8PuppijetNEmEF = 0;
   AK8PuppijetNHadEF = 0;
   AK8PuppijetEleEF = 0;
   AK8PuppijetMuoEF = 0;
   AK8PuppijetCMulti = 0;
   AK8PuppijetEleMulti = 0;
   AK8PuppijetMuoMulti = 0;
   AK8PuppijetSSV = 0;
   AK8PuppijetCSV = 0;
   AK8PuppijetSSVHE = 0;
   AK8PuppijetCISVV2 = 0;
   AK8PuppijetTCHP = 0;
   AK8PuppijetTCHE = 0;
   AK8PuppijetJP = 0;
   AK8PuppijetJBP = 0;
   AK8PuppijetTau1 = 0;
   AK8PuppijetTau2 = 0;
   AK8PuppijetTau3 = 0;
   AK8PuppijetTau21 = 0;
   AK8PuppijetSDmass = 0;
   AK8PuppijetGenSDmass = 0;
   AK8PuppinSubSDJet = 0;
   AK8PuppisubjetSDFatJetIndex = 0;
   AK8PuppisubjetSDPx = 0;
   AK8PuppisubjetSDPy = 0;
   AK8PuppisubjetSDPz = 0;
   AK8PuppisubjetSDE = 0;
   AK8PuppisubjetSDRawFactor = 0;
   AK8PuppisubjetSDPartonFlavor = 0;
   AK8PuppisubjetSDHadronFlavor = 0;
   AK8PuppisubjetSDCSV = 0;
   AK8Puppijet_DoubleSV = 0;
   AK8Puppijet_nSV = 0;
   AK8Puppijet_SVMass = 0;
   CA15PuppijetP4 = 0;
   CA15PuppigenjetP4 = 0;
   CA15PuppigenjetEM = 0;
   CA15PuppigenjetHAD = 0;
   CA15PuppigenjetINV = 0;
   CA15PuppigenjetAUX = 0;
   CA15PuppimatchedDR = 0;
   CA15PuppijetRawFactor = 0;
   CA15PuppiunCorrJetP4 = 0;
   CA15PuppijetArea = 0;
   CA15PuppijetCorrUncUp = 0;
   CA15PuppijetCorrUncDown = 0;
   CA15PuppijetCharge = 0;
   CA15PuppijetPartonFlavor = 0;
   CA15PuppijetHadronFlavor = 0;
   CA15PuppijetPassIDLoose = 0;
   CA15PuppijetPassIDTight = 0;
   CA15PuppijetCEmEF = 0;
   CA15PuppijetCHadEF = 0;
   CA15PuppijetPhoEF = 0;
   CA15PuppijetNEmEF = 0;
   CA15PuppijetNHadEF = 0;
   CA15PuppijetEleEF = 0;
   CA15PuppijetMuoEF = 0;
   CA15PuppijetCMulti = 0;
   CA15PuppijetEleMulti = 0;
   CA15PuppijetMuoMulti = 0;
   CA15PuppijetSSV = 0;
   CA15PuppijetCSV = 0;
   CA15PuppijetSSVHE = 0;
   CA15PuppijetCISVV2 = 0;
   CA15PuppijetTCHP = 0;
   CA15PuppijetTCHE = 0;
   CA15PuppijetJP = 0;
   CA15PuppijetJBP = 0;
   CA15PuppijetTau1 = 0;
   CA15PuppijetTau2 = 0;
   CA15PuppijetTau3 = 0;
   CA15PuppijetTau21 = 0;
   CA15PuppijetSDmass = 0;
   CA15PuppijetGenSDmass = 0;
   CA15PuppinSubSDJet = 0;
   CA15PuppisubjetSDFatJetIndex = 0;
   CA15PuppisubjetSDPx = 0;
   CA15PuppisubjetSDPy = 0;
   CA15PuppisubjetSDPz = 0;
   CA15PuppisubjetSDE = 0;
   CA15PuppisubjetSDRawFactor = 0;
   CA15PuppisubjetSDPartonFlavor = 0;
   CA15PuppisubjetSDHadronFlavor = 0;
   CA15PuppisubjetSDCSV = 0;
   CA15Puppi_doublebtag = 0;
   CA15PuppiECF_2_3_10 = 0;
   CA15PuppiECF_1_2_10 = 0;
   CA15Puppijet_DoubleSV = 0;
   CA15Puppijet_nSV = 0;
   CA15Puppijet_SVMass = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("pu_nTrueInt", &pu_nTrueInt, &b_pu_nTrueInt);
   fChain->SetBranchAddress("pu_nPUVert", &pu_nPUVert, &b_pu_nPUVert);
   fChain->SetBranchAddress("isData", &isData, &b_isData);
   fChain->SetBranchAddress("eventId", &eventId, &b_eventId);
   fChain->SetBranchAddress("runId", &runId, &b_runId);
   fChain->SetBranchAddress("lumiSection", &lumiSection, &b_lumiSection);
   fChain->SetBranchAddress("bunchXing", &bunchXing, &b_bunchXing);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("vertexP3", &vertexP3, &b_vertexP3);
   fChain->SetBranchAddress("pfMetCorrPt", &pfMetCorrPt, &b_pfMetCorrPt);
   fChain->SetBranchAddress("pfMetCorrPhi", &pfMetCorrPhi, &b_pfMetCorrPhi);
   fChain->SetBranchAddress("pfMetCorrSumEt", &pfMetCorrSumEt, &b_pfMetCorrSumEt);
   fChain->SetBranchAddress("pfMetCorrSig", &pfMetCorrSig, &b_pfMetCorrSig);
   fChain->SetBranchAddress("pfMetCorrUnc", &pfMetCorrUnc, &b_pfMetCorrUnc);
   fChain->SetBranchAddress("pfMetRawPt", &pfMetRawPt, &b_pfMetRawPt);
   fChain->SetBranchAddress("pfMetRawPhi", &pfMetRawPhi, &b_pfMetRawPhi);
   fChain->SetBranchAddress("pfMetRawSumEt", &pfMetRawSumEt, &b_pfMetRawSumEt);
   fChain->SetBranchAddress("pfMetRawCov00", &pfMetRawCov00, &b_pfMetRawCov00);
   fChain->SetBranchAddress("pfMetRawCov01", &pfMetRawCov01, &b_pfMetRawCov01);
   fChain->SetBranchAddress("pfMetRawCov10", &pfMetRawCov10, &b_pfMetRawCov10);
   fChain->SetBranchAddress("pfMetRawCov11", &pfMetRawCov11, &b_pfMetRawCov11);
   fChain->SetBranchAddress("pfmvaMetPt", &pfmvaMetPt, &b_pfmvaMetPt);
   fChain->SetBranchAddress("pfmvaMetPhi", &pfmvaMetPhi, &b_pfmvaMetPhi);
   fChain->SetBranchAddress("pfmvaMetSumEt", &pfmvaMetSumEt, &b_pfmvaMetSumEt);
   fChain->SetBranchAddress("pfmvaMetSig", &pfmvaMetSig, &b_pfmvaMetSig);
   fChain->SetBranchAddress("pfpuppiMETPt", &pfpuppiMETPt, &b_pfpuppiMETPt);
   fChain->SetBranchAddress("pfpuppiMETPhi", &pfpuppiMETPhi, &b_pfpuppiMETPhi);
   fChain->SetBranchAddress("pfpuppiMETSumEt", &pfpuppiMETSumEt, &b_pfpuppiMETSumEt);
   fChain->SetBranchAddress("pfpuppiMETSig", &pfpuppiMETSig, &b_pfpuppiMETSig);
   fChain->SetBranchAddress("pfpuppiMETUnc", &pfpuppiMETUnc, &b_pfpuppiMETUnc);
   fChain->SetBranchAddress("hlt_nTrigs", &hlt_nTrigs, &b_hlt_nTrigs);
   fChain->SetBranchAddress("hlt_trigResult", &hlt_trigResult, &b_hlt_trigResult);
   fChain->SetBranchAddress("hlt_trigName", &hlt_trigName, &b_hlt_trigName);
   fChain->SetBranchAddress("hlt_trigPrescale", &hlt_trigPrescale, &b_hlt_trigPrescale);
   fChain->SetBranchAddress("hlt_nfilters", &hlt_nfilters, &b_hlt_nfilters);
   fChain->SetBranchAddress("hlt_filterResult", &hlt_filterResult, &b_hlt_filterResult);
   fChain->SetBranchAddress("hlt_hbhet", &hlt_hbhet, &b_hlt_hbhet);
   fChain->SetBranchAddress("hlt_filterName", &hlt_filterName, &b_hlt_filterName);
   fChain->SetBranchAddress("hlt_filterbadChCandidate", &hlt_filterbadChCandidate, &b_hlt_filterbadChCandidate);
   fChain->SetBranchAddress("hlt_filterbadPFMuon", &hlt_filterbadPFMuon, &b_hlt_filterbadPFMuon);
   fChain->SetBranchAddress("hlt_filterbadGlobalMuon", &hlt_filterbadGlobalMuon, &b_hlt_filterbadGlobalMuon);
   fChain->SetBranchAddress("hlt_filtercloneGlobalMuon", &hlt_filtercloneGlobalMuon, &b_hlt_filtercloneGlobalMuon);
   fChain->SetBranchAddress("genMET_true", &genMET_true, &b_genMET_true);
   fChain->SetBranchAddress("genMET_calo", &genMET_calo, &b_genMET_calo);
   fChain->SetBranchAddress("genMET_caloNonPrompt", &genMET_caloNonPrompt, &b_genMET_caloNonPrompt);
   fChain->SetBranchAddress("ptHat", &ptHat, &b_ptHat);
   fChain->SetBranchAddress("mcWeight", &mcWeight, &b_mcWeight);
   fChain->SetBranchAddress("HT", &HT, &b_HT);
   fChain->SetBranchAddress("pdf", &pdf, &b_pdf);
   fChain->SetBranchAddress("originalLHEweight", &originalLHEweight, &b_originalLHEweight);
   fChain->SetBranchAddress("pdfscaleSysWeights", &pdfscaleSysWeights, &b_pdfscaleSysWeights);
   fChain->SetBranchAddress("nGenPar", &nGenPar, &b_nGenPar);
   fChain->SetBranchAddress("genParP4", &genParP4, &b_genParP4);
   fChain->SetBranchAddress("genParQ", &genParQ, &b_genParQ);
   fChain->SetBranchAddress("genParId", &genParId, &b_genParId);
   fChain->SetBranchAddress("genParSt", &genParSt, &b_genParSt);
   fChain->SetBranchAddress("genMomParId", &genMomParId, &b_genMomParId);
   fChain->SetBranchAddress("genParIndex", &genParIndex, &b_genParIndex);
   fChain->SetBranchAddress("genNMo", &genNMo, &b_genNMo);
   fChain->SetBranchAddress("genNDa", &genNDa, &b_genNDa);
   fChain->SetBranchAddress("genMo1", &genMo1, &b_genMo1);
   fChain->SetBranchAddress("genMo2", &genMo2, &b_genMo2);
   fChain->SetBranchAddress("genDa1", &genDa1, &b_genDa1);
   fChain->SetBranchAddress("genDa2", &genDa2, &b_genDa2);
   fChain->SetBranchAddress("genStFlag", &genStFlag, &b_genStFlag);
   fChain->SetBranchAddress("ak4nGenJet", &ak4nGenJet, &b_ak4nGenJet);
   fChain->SetBranchAddress("ak4GenJetP4", &ak4GenJetP4, &b_ak4GenJetP4);
   fChain->SetBranchAddress("ak8nGenJet", &ak8nGenJet, &b_ak8nGenJet);
   fChain->SetBranchAddress("ak8GenJetP4", &ak8GenJetP4, &b_ak8GenJetP4);
   fChain->SetBranchAddress("eleRho", &eleRho, &b_eleRho);
   fChain->SetBranchAddress("nEle", &nEle, &b_nEle);
   fChain->SetBranchAddress("eleP4", &eleP4, &b_eleP4);
   fChain->SetBranchAddress("eleInBarrel", &eleInBarrel, &b_eleInBarrel);
   fChain->SetBranchAddress("eleInEndcap", &eleInEndcap, &b_eleInEndcap);
   fChain->SetBranchAddress("eleCharge", &eleCharge, &b_eleCharge);
   fChain->SetBranchAddress("eleChargeConsistent", &eleChargeConsistent, &b_eleChargeConsistent);
   fChain->SetBranchAddress("elecaloEnergy", &elecaloEnergy, &b_elecaloEnergy);
   fChain->SetBranchAddress("eleScEt", &eleScEt, &b_eleScEt);
   fChain->SetBranchAddress("eleScEn", &eleScEn, &b_eleScEn);
   fChain->SetBranchAddress("eleScPreEn", &eleScPreEn, &b_eleScPreEn);
   fChain->SetBranchAddress("eleScEta", &eleScEta, &b_eleScEta);
   fChain->SetBranchAddress("eleScPhi", &eleScPhi, &b_eleScPhi);
   fChain->SetBranchAddress("eleScRawEn", &eleScRawEn, &b_eleScRawEn);
   fChain->SetBranchAddress("eleScEtaWidth", &eleScEtaWidth, &b_eleScEtaWidth);
   fChain->SetBranchAddress("eleScPhiWidth", &eleScPhiWidth, &b_eleScPhiWidth);
   fChain->SetBranchAddress("eleR9", &eleR9, &b_eleR9);
   fChain->SetBranchAddress("eleHoverE", &eleHoverE, &b_eleHoverE);
   fChain->SetBranchAddress("eleD0", &eleD0, &b_eleD0);
   fChain->SetBranchAddress("eleDz", &eleDz, &b_eleDz);
   fChain->SetBranchAddress("eleEoverP", &eleEoverP, &b_eleEoverP);
   fChain->SetBranchAddress("eleBrem", &eleBrem, &b_eleBrem);
   fChain->SetBranchAddress("eledEtaAtVtx", &eledEtaAtVtx, &b_eledEtaAtVtx);
   fChain->SetBranchAddress("eledPhiAtVtx", &eledPhiAtVtx, &b_eledPhiAtVtx);
   fChain->SetBranchAddress("eleSigmaIEtaIEta", &eleSigmaIEtaIEta, &b_eleSigmaIEtaIEta);
   fChain->SetBranchAddress("eleSigmaIEtaIPhi", &eleSigmaIEtaIPhi, &b_eleSigmaIEtaIPhi);
   fChain->SetBranchAddress("eleSigmaIPhiIPhi", &eleSigmaIPhiIPhi, &b_eleSigmaIPhiIPhi);
   fChain->SetBranchAddress("eleConvVeto", &eleConvVeto, &b_eleConvVeto);
   fChain->SetBranchAddress("eleMissHits", &eleMissHits, &b_eleMissHits);
   fChain->SetBranchAddress("eleEoverPInv", &eleEoverPInv, &b_eleEoverPInv);
   fChain->SetBranchAddress("eleEtaseedAtVtx", &eleEtaseedAtVtx, &b_eleEtaseedAtVtx);
   fChain->SetBranchAddress("eleE1x5", &eleE1x5, &b_eleE1x5);
   fChain->SetBranchAddress("eleE2x5", &eleE2x5, &b_eleE2x5);
   fChain->SetBranchAddress("eleE5x5", &eleE5x5, &b_eleE5x5);
   fChain->SetBranchAddress("eleSigmaIEtaIEtaFull5x5", &eleSigmaIEtaIEtaFull5x5, &b_eleSigmaIEtaIEtaFull5x5);
   fChain->SetBranchAddress("eleE1x5Full5x5", &eleE1x5Full5x5, &b_eleE1x5Full5x5);
   fChain->SetBranchAddress("eleE2x5Full5x5", &eleE2x5Full5x5, &b_eleE2x5Full5x5);
   fChain->SetBranchAddress("eleE5x5Full5x5", &eleE5x5Full5x5, &b_eleE5x5Full5x5);
   fChain->SetBranchAddress("eleR9Full5x5", &eleR9Full5x5, &b_eleR9Full5x5);
   fChain->SetBranchAddress("eleChHadIso", &eleChHadIso, &b_eleChHadIso);
   fChain->SetBranchAddress("eleNeHadIso", &eleNeHadIso, &b_eleNeHadIso);
   fChain->SetBranchAddress("eleGamIso", &eleGamIso, &b_eleGamIso);
   fChain->SetBranchAddress("elePUPt", &elePUPt, &b_elePUPt);
   fChain->SetBranchAddress("eleEcalPFClusterIso", &eleEcalPFClusterIso, &b_eleEcalPFClusterIso);
   fChain->SetBranchAddress("eleHcalPFClusterIso", &eleHcalPFClusterIso, &b_eleHcalPFClusterIso);
   fChain->SetBranchAddress("eleMiniIso_ch", &eleMiniIso_ch, &b_eleMiniIso_ch);
   fChain->SetBranchAddress("eleMiniIso_nh", &eleMiniIso_nh, &b_eleMiniIso_nh);
   fChain->SetBranchAddress("eleMiniIso_ph", &eleMiniIso_ph, &b_eleMiniIso_ph);
   fChain->SetBranchAddress("eleMiniIso_pu", &eleMiniIso_pu, &b_eleMiniIso_pu);
   fChain->SetBranchAddress("eleMiniIso_r", &eleMiniIso_r, &b_eleMiniIso_r);
   fChain->SetBranchAddress("eleMiniIsoBeta", &eleMiniIsoBeta, &b_eleMiniIsoBeta);
   fChain->SetBranchAddress("eleMiniIsoEA", &eleMiniIsoEA, &b_eleMiniIsoEA);
   fChain->SetBranchAddress("eleEcalDrivenSeed", &eleEcalDrivenSeed, &b_eleEcalDrivenSeed);
   fChain->SetBranchAddress("eleEcalDriven", &eleEcalDriven, &b_eleEcalDriven);
   fChain->SetBranchAddress("eleDr03EcalRecHitSumEt", &eleDr03EcalRecHitSumEt, &b_eleDr03EcalRecHitSumEt);
   fChain->SetBranchAddress("eleDr03HcalDepth1TowerSumEt", &eleDr03HcalDepth1TowerSumEt, &b_eleDr03HcalDepth1TowerSumEt);
   fChain->SetBranchAddress("eleDr03HcalDepth2TowerSumEt", &eleDr03HcalDepth2TowerSumEt, &b_eleDr03HcalDepth2TowerSumEt);
   fChain->SetBranchAddress("eleDr03HcalTowerSumEt", &eleDr03HcalTowerSumEt, &b_eleDr03HcalTowerSumEt);
   fChain->SetBranchAddress("eleDr03TkSumPt", &eleDr03TkSumPt, &b_eleDr03TkSumPt);
   fChain->SetBranchAddress("eleIsPassVeto", &eleIsPassVeto, &b_eleIsPassVeto);
   fChain->SetBranchAddress("eleIsPassLoose", &eleIsPassLoose, &b_eleIsPassLoose);
   fChain->SetBranchAddress("eleIsPassMedium", &eleIsPassMedium, &b_eleIsPassMedium);
   fChain->SetBranchAddress("eleIsPassTight", &eleIsPassTight, &b_eleIsPassTight);
   fChain->SetBranchAddress("eleIsPassHEEP", &eleIsPassHEEP, &b_eleIsPassHEEP);
   fChain->SetBranchAddress("eleIsPassVetoNoIso", &eleIsPassVetoNoIso, &b_eleIsPassVetoNoIso);
   fChain->SetBranchAddress("eleIsPassLooseNoIso", &eleIsPassLooseNoIso, &b_eleIsPassLooseNoIso);
   fChain->SetBranchAddress("eleIsPassMediumNoIso", &eleIsPassMediumNoIso, &b_eleIsPassMediumNoIso);
   fChain->SetBranchAddress("eleIsPassTightNoIso", &eleIsPassTightNoIso, &b_eleIsPassTightNoIso);
   fChain->SetBranchAddress("eleIsPassHEEPNoIso", &eleIsPassHEEPNoIso, &b_eleIsPassHEEPNoIso);
   fChain->SetBranchAddress("eleIsPassMVAMedium", &eleIsPassMVAMedium, &b_eleIsPassMVAMedium);
   fChain->SetBranchAddress("eleIsPassMVATight", &eleIsPassMVATight, &b_eleIsPassMVATight);
   fChain->SetBranchAddress("eleMVAValue", &eleMVAValue, &b_eleMVAValue);
   fChain->SetBranchAddress("eleMVACategory", &eleMVACategory, &b_eleMVACategory);
   fChain->SetBranchAddress("nMu", &nMu, &b_nMu);
   fChain->SetBranchAddress("muP4", &muP4, &b_muP4);
   fChain->SetBranchAddress("muType", &muType, &b_muType);
   fChain->SetBranchAddress("muCharge", &muCharge, &b_muCharge);
   fChain->SetBranchAddress("isGlobalMuon", &isGlobalMuon, &b_isGlobalMuon);
   fChain->SetBranchAddress("isTrackerMuon", &isTrackerMuon, &b_isTrackerMuon);
   fChain->SetBranchAddress("isPFMuon", &isPFMuon, &b_isPFMuon);
   fChain->SetBranchAddress("isTightMuon", &isTightMuon, &b_isTightMuon);
   fChain->SetBranchAddress("isLooseMuon", &isLooseMuon, &b_isLooseMuon);
   fChain->SetBranchAddress("isMediumMuon", &isMediumMuon, &b_isMediumMuon);
   fChain->SetBranchAddress("isSoftMuon", &isSoftMuon, &b_isSoftMuon);
   fChain->SetBranchAddress("isHighPtMuon", &isHighPtMuon, &b_isHighPtMuon);
   fChain->SetBranchAddress("isCustomTrackerMuon", &isCustomTrackerMuon, &b_isCustomTrackerMuon);
   fChain->SetBranchAddress("muITrkID", &muITrkID, &b_muITrkID);
   fChain->SetBranchAddress("muSegID", &muSegID, &b_muSegID);
   fChain->SetBranchAddress("muNSegs", &muNSegs, &b_muNSegs);
   fChain->SetBranchAddress("muGood", &muGood, &b_muGood);
   fChain->SetBranchAddress("muIsGood", &muIsGood, &b_muIsGood);
   fChain->SetBranchAddress("muTrkPt", &muTrkPt, &b_muTrkPt);
   fChain->SetBranchAddress("muTrkPtErr", &muTrkPtErr, &b_muTrkPtErr);
   fChain->SetBranchAddress("mudxy", &mudxy, &b_mudxy);
   fChain->SetBranchAddress("mudz", &mudz, &b_mudz);
   fChain->SetBranchAddress("musegmentCompatibility", &musegmentCompatibility, &b_musegmentCompatibility);
   fChain->SetBranchAddress("muchi2LocalPosition", &muchi2LocalPosition, &b_muchi2LocalPosition);
   fChain->SetBranchAddress("mutrkKink", &mutrkKink, &b_mutrkKink);
   fChain->SetBranchAddress("muInnerdxy", &muInnerdxy, &b_muInnerdxy);
   fChain->SetBranchAddress("muInnerdz", &muInnerdz, &b_muInnerdz);
   fChain->SetBranchAddress("muTrkLayers", &muTrkLayers, &b_muTrkLayers);
   fChain->SetBranchAddress("muPixelLayers", &muPixelLayers, &b_muPixelLayers);
   fChain->SetBranchAddress("muPixelHits", &muPixelHits, &b_muPixelHits);
   fChain->SetBranchAddress("muTrkQuality", &muTrkQuality, &b_muTrkQuality);
   fChain->SetBranchAddress("muHits", &muHits, &b_muHits);
   fChain->SetBranchAddress("muChi2NDF", &muChi2NDF, &b_muChi2NDF);
   fChain->SetBranchAddress("muInnervalidFraction", &muInnervalidFraction, &b_muInnervalidFraction);
   fChain->SetBranchAddress("muMatches", &muMatches, &b_muMatches);
   fChain->SetBranchAddress("muTrkIso", &muTrkIso, &b_muTrkIso);
   fChain->SetBranchAddress("muHcalIso", &muHcalIso, &b_muHcalIso);
   fChain->SetBranchAddress("muEcalIso", &muEcalIso, &b_muEcalIso);
   fChain->SetBranchAddress("muChHadIso", &muChHadIso, &b_muChHadIso);
   fChain->SetBranchAddress("muNeHadIso", &muNeHadIso, &b_muNeHadIso);
   fChain->SetBranchAddress("muGamIso", &muGamIso, &b_muGamIso);
   fChain->SetBranchAddress("muPUPt", &muPUPt, &b_muPUPt);
   fChain->SetBranchAddress("muInnerTrkPt", &muInnerTrkPt, &b_muInnerTrkPt);
   fChain->SetBranchAddress("muMiniIso_ch", &muMiniIso_ch, &b_muMiniIso_ch);
   fChain->SetBranchAddress("muMiniIso_nh", &muMiniIso_nh, &b_muMiniIso_nh);
   fChain->SetBranchAddress("muMiniIso_ph", &muMiniIso_ph, &b_muMiniIso_ph);
   fChain->SetBranchAddress("muMiniIso_pu", &muMiniIso_pu, &b_muMiniIso_pu);
   fChain->SetBranchAddress("muMiniIso_r", &muMiniIso_r, &b_muMiniIso_r);
   fChain->SetBranchAddress("muMiniIsoBeta", &muMiniIsoBeta, &b_muMiniIsoBeta);
   fChain->SetBranchAddress("muMiniIsoEA", &muMiniIsoEA, &b_muMiniIsoEA);
   fChain->SetBranchAddress("HPSTau_n", &HPSTau_n, &b_HPSTau_n);
   fChain->SetBranchAddress("taupt", &taupt, &b_taupt);
   fChain->SetBranchAddress("HPSTau_4Momentum", &HPSTau_4Momentum, &b_HPSTau_4Momentum);
   fChain->SetBranchAddress("HPSTau_Vposition", &HPSTau_Vposition, &b_HPSTau_Vposition);
   fChain->SetBranchAddress("HPSTau_leadPFChargedHadrCand", &HPSTau_leadPFChargedHadrCand, &b_HPSTau_leadPFChargedHadrCand);
   fChain->SetBranchAddress("HPSTau_leadPFChargedHadrCand_trackRef", &HPSTau_leadPFChargedHadrCand_trackRef, &b_HPSTau_leadPFChargedHadrCand_trackRef);
   fChain->SetBranchAddress("disc_againstElectronLoose", &disc_againstElectronLoose, &b_disc_againstElectronLoose);
   fChain->SetBranchAddress("disc_againstElectronMedium", &disc_againstElectronMedium, &b_disc_againstElectronMedium);
   fChain->SetBranchAddress("disc_againstElectronTight", &disc_againstElectronTight, &b_disc_againstElectronTight);
   fChain->SetBranchAddress("disc_againstElectronLooseMVA5", &disc_againstElectronLooseMVA5, &b_disc_againstElectronLooseMVA5);
   fChain->SetBranchAddress("disc_againstElectronMediumMVA5", &disc_againstElectronMediumMVA5, &b_disc_againstElectronMediumMVA5);
   fChain->SetBranchAddress("disc_againstElectronTightMVA5", &disc_againstElectronTightMVA5, &b_disc_againstElectronTightMVA5);
   fChain->SetBranchAddress("disc_againstElectronVLooseMVA5", &disc_againstElectronVLooseMVA5, &b_disc_againstElectronVLooseMVA5);
   fChain->SetBranchAddress("disc_againstElectronVTightMVA5", &disc_againstElectronVTightMVA5, &b_disc_againstElectronVTightMVA5);
   fChain->SetBranchAddress("disc_againstMuonLoose", &disc_againstMuonLoose, &b_disc_againstMuonLoose);
   fChain->SetBranchAddress("disc_againstMuonMedium", &disc_againstMuonMedium, &b_disc_againstMuonMedium);
   fChain->SetBranchAddress("disc_againstMuonTight", &disc_againstMuonTight, &b_disc_againstMuonTight);
   fChain->SetBranchAddress("disc_againstMuonLoose2", &disc_againstMuonLoose2, &b_disc_againstMuonLoose2);
   fChain->SetBranchAddress("disc_againstMuonMedium2", &disc_againstMuonMedium2, &b_disc_againstMuonMedium2);
   fChain->SetBranchAddress("disc_againstMuonTight2", &disc_againstMuonTight2, &b_disc_againstMuonTight2);
   fChain->SetBranchAddress("disc_againstMuonLooseMVA", &disc_againstMuonLooseMVA, &b_disc_againstMuonLooseMVA);
   fChain->SetBranchAddress("disc_againstMuonMediumMVA", &disc_againstMuonMediumMVA, &b_disc_againstMuonMediumMVA);
   fChain->SetBranchAddress("disc_againstMuonTightMVA", &disc_againstMuonTightMVA, &b_disc_againstMuonTightMVA);
   fChain->SetBranchAddress("disc_againstMuonLoose3", &disc_againstMuonLoose3, &b_disc_againstMuonLoose3);
   fChain->SetBranchAddress("disc_againstMuonTight3", &disc_againstMuonTight3, &b_disc_againstMuonTight3);
   fChain->SetBranchAddress("disc_byVLooseCombinedIsolationDeltaBetaCorr", &disc_byVLooseCombinedIsolationDeltaBetaCorr, &b_disc_byVLooseCombinedIsolationDeltaBetaCorr);
   fChain->SetBranchAddress("disc_byLooseCombinedIsolationDeltaBetaCorr", &disc_byLooseCombinedIsolationDeltaBetaCorr, &b_disc_byLooseCombinedIsolationDeltaBetaCorr);
   fChain->SetBranchAddress("disc_byMediumCombinedIsolationDeltaBetaCorr", &disc_byMediumCombinedIsolationDeltaBetaCorr, &b_disc_byMediumCombinedIsolationDeltaBetaCorr);
   fChain->SetBranchAddress("disc_byTightCombinedIsolationDeltaBetaCorr", &disc_byTightCombinedIsolationDeltaBetaCorr, &b_disc_byTightCombinedIsolationDeltaBetaCorr);
   fChain->SetBranchAddress("disc_byLooseIsolation", &disc_byLooseIsolation, &b_disc_byLooseIsolation);
   fChain->SetBranchAddress("disc_byVLooseIsolationMVA3newDMwLT", &disc_byVLooseIsolationMVA3newDMwLT, &b_disc_byVLooseIsolationMVA3newDMwLT);
   fChain->SetBranchAddress("disc_byLooseIsolationMVA3newDMwLT", &disc_byLooseIsolationMVA3newDMwLT, &b_disc_byLooseIsolationMVA3newDMwLT);
   fChain->SetBranchAddress("disc_byMediumIsolationMVA3newDMwLT", &disc_byMediumIsolationMVA3newDMwLT, &b_disc_byMediumIsolationMVA3newDMwLT);
   fChain->SetBranchAddress("disc_byTightIsolationMVA3newDMwLT", &disc_byTightIsolationMVA3newDMwLT, &b_disc_byTightIsolationMVA3newDMwLT);
   fChain->SetBranchAddress("disc_byVTightIsolationMVA3newDMwLT", &disc_byVTightIsolationMVA3newDMwLT, &b_disc_byVTightIsolationMVA3newDMwLT);
   fChain->SetBranchAddress("disc_byVVTightIsolationMVA3newDMwLT", &disc_byVVTightIsolationMVA3newDMwLT, &b_disc_byVVTightIsolationMVA3newDMwLT);
   fChain->SetBranchAddress("disc_byVLooseIsolationMVA3newDMwoLT", &disc_byVLooseIsolationMVA3newDMwoLT, &b_disc_byVLooseIsolationMVA3newDMwoLT);
   fChain->SetBranchAddress("disc_byLooseIsolationMVA3newDMwoLT", &disc_byLooseIsolationMVA3newDMwoLT, &b_disc_byLooseIsolationMVA3newDMwoLT);
   fChain->SetBranchAddress("disc_byMediumIsolationMVA3newDMwoLT", &disc_byMediumIsolationMVA3newDMwoLT, &b_disc_byMediumIsolationMVA3newDMwoLT);
   fChain->SetBranchAddress("disc_byTightIsolationMVA3newDMwoLT", &disc_byTightIsolationMVA3newDMwoLT, &b_disc_byTightIsolationMVA3newDMwoLT);
   fChain->SetBranchAddress("disc_byVTightIsolationMVA3newDMwoLT", &disc_byVTightIsolationMVA3newDMwoLT, &b_disc_byVTightIsolationMVA3newDMwoLT);
   fChain->SetBranchAddress("disc_byVVTightIsolationMVA3newDMwoLT", &disc_byVVTightIsolationMVA3newDMwoLT, &b_disc_byVVTightIsolationMVA3newDMwoLT);
   fChain->SetBranchAddress("disc_byVLooseIsolationMVA3oldDMwLT", &disc_byVLooseIsolationMVA3oldDMwLT, &b_disc_byVLooseIsolationMVA3oldDMwLT);
   fChain->SetBranchAddress("disc_byLooseIsolationMVA3oldDMwLT", &disc_byLooseIsolationMVA3oldDMwLT, &b_disc_byLooseIsolationMVA3oldDMwLT);
   fChain->SetBranchAddress("disc_byMediumIsolationMVA3oldDMwLT", &disc_byMediumIsolationMVA3oldDMwLT, &b_disc_byMediumIsolationMVA3oldDMwLT);
   fChain->SetBranchAddress("disc_byTightIsolationMVA3oldDMwLT", &disc_byTightIsolationMVA3oldDMwLT, &b_disc_byTightIsolationMVA3oldDMwLT);
   fChain->SetBranchAddress("disc_byVTightIsolationMVA3oldDMwLT", &disc_byVTightIsolationMVA3oldDMwLT, &b_disc_byVTightIsolationMVA3oldDMwLT);
   fChain->SetBranchAddress("disc_byVVTightIsolationMVA3oldDMwLT", &disc_byVVTightIsolationMVA3oldDMwLT, &b_disc_byVVTightIsolationMVA3oldDMwLT);
   fChain->SetBranchAddress("disc_byVLooseIsolationMVA3oldDMwoLT", &disc_byVLooseIsolationMVA3oldDMwoLT, &b_disc_byVLooseIsolationMVA3oldDMwoLT);
   fChain->SetBranchAddress("disc_byLooseIsolationMVA3oldDMwoLT", &disc_byLooseIsolationMVA3oldDMwoLT, &b_disc_byLooseIsolationMVA3oldDMwoLT);
   fChain->SetBranchAddress("disc_byMediumIsolationMVA3oldDMwoLT", &disc_byMediumIsolationMVA3oldDMwoLT, &b_disc_byMediumIsolationMVA3oldDMwoLT);
   fChain->SetBranchAddress("disc_byTightIsolationMVA3oldDMwoLT", &disc_byTightIsolationMVA3oldDMwoLT, &b_disc_byTightIsolationMVA3oldDMwoLT);
   fChain->SetBranchAddress("disc_byVTightIsolationMVA3oldDMwoLT", &disc_byVTightIsolationMVA3oldDMwoLT, &b_disc_byVTightIsolationMVA3oldDMwoLT);
   fChain->SetBranchAddress("disc_byVVTightIsolationMVA3oldDMwoLT", &disc_byVVTightIsolationMVA3oldDMwoLT, &b_disc_byVVTightIsolationMVA3oldDMwoLT);
   fChain->SetBranchAddress("disc_byLooseCombinedIsolationDeltaBetaCorr3Hits", &disc_byLooseCombinedIsolationDeltaBetaCorr3Hits, &b_disc_byLooseCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("disc_byMediumCombinedIsolationDeltaBetaCorr3Hits", &disc_byMediumCombinedIsolationDeltaBetaCorr3Hits, &b_disc_byMediumCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("disc_byTightCombinedIsolationDeltaBetaCorr3Hits", &disc_byTightCombinedIsolationDeltaBetaCorr3Hits, &b_disc_byTightCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("disc_decayModeFinding", &disc_decayModeFinding, &b_disc_decayModeFinding);
   fChain->SetBranchAddress("disc_decayModeFindingNewDMs", &disc_decayModeFindingNewDMs, &b_disc_decayModeFindingNewDMs);
   fChain->SetBranchAddress("disc_chargedIsoPtSum", &disc_chargedIsoPtSum, &b_disc_chargedIsoPtSum);
   fChain->SetBranchAddress("disc_neutralIsoPtSum", &disc_neutralIsoPtSum, &b_disc_neutralIsoPtSum);
   fChain->SetBranchAddress("disc_puCorrPtSum", &disc_puCorrPtSum, &b_disc_puCorrPtSum);
   fChain->SetBranchAddress("HPSTau_NewVz", &HPSTau_NewVz, &b_HPSTau_NewVz);
   fChain->SetBranchAddress("HPSTau_charge", &HPSTau_charge, &b_HPSTau_charge);
   fChain->SetBranchAddress("nPho", &nPho, &b_nPho);
   fChain->SetBranchAddress("phoP4", &phoP4, &b_phoP4);
   fChain->SetBranchAddress("phoIsPassTight", &phoIsPassTight, &b_phoIsPassTight);
   fChain->SetBranchAddress("phoIsPassLoose", &phoIsPassLoose, &b_phoIsPassLoose);
   fChain->SetBranchAddress("phoIsPassMedium", &phoIsPassMedium, &b_phoIsPassMedium);
   fChain->SetBranchAddress("phoIDMVA", &phoIDMVA, &b_phoIDMVA);
   fChain->SetBranchAddress("phoSCE", &phoSCE, &b_phoSCE);
   fChain->SetBranchAddress("phoSCRawE", &phoSCRawE, &b_phoSCRawE);
   fChain->SetBranchAddress("phoSCEta", &phoSCEta, &b_phoSCEta);
   fChain->SetBranchAddress("phoSCPhi", &phoSCPhi, &b_phoSCPhi);
   fChain->SetBranchAddress("phoSCEtaWidth", &phoSCEtaWidth, &b_phoSCEtaWidth);
   fChain->SetBranchAddress("phoSCPhiWidth", &phoSCPhiWidth, &b_phoSCPhiWidth);
   fChain->SetBranchAddress("phoSCBrem", &phoSCBrem, &b_phoSCBrem);
   fChain->SetBranchAddress("phohasPixelSeed", &phohasPixelSeed, &b_phohasPixelSeed);
   fChain->SetBranchAddress("phoEleVeto", &phoEleVeto, &b_phoEleVeto);
   fChain->SetBranchAddress("phoR9", &phoR9, &b_phoR9);
   fChain->SetBranchAddress("phoHoverE", &phoHoverE, &b_phoHoverE);
   fChain->SetBranchAddress("phoSigmaIEtaIEta", &phoSigmaIEtaIEta, &b_phoSigmaIEtaIEta);
   fChain->SetBranchAddress("phoSigmaIEtaIPhi", &phoSigmaIEtaIPhi, &b_phoSigmaIEtaIPhi);
   fChain->SetBranchAddress("phoSigmaIPhiIPhi", &phoSigmaIPhiIPhi, &b_phoSigmaIPhiIPhi);
   fChain->SetBranchAddress("phoSigmaIEtaIEtaFull5x5", &phoSigmaIEtaIEtaFull5x5, &b_phoSigmaIEtaIEtaFull5x5);
   fChain->SetBranchAddress("phoR9Full5x5", &phoR9Full5x5, &b_phoR9Full5x5);
   fChain->SetBranchAddress("phoPFChIso", &phoPFChIso, &b_phoPFChIso);
   fChain->SetBranchAddress("phoPFPhoIso", &phoPFPhoIso, &b_phoPFPhoIso);
   fChain->SetBranchAddress("phoPFNeuIso", &phoPFNeuIso, &b_phoPFNeuIso);
   fChain->SetBranchAddress("THINnJet", &THINnJet, &b_THINnJet);
   fChain->SetBranchAddress("THINjetP4", &THINjetP4, &b_THINjetP4);
   fChain->SetBranchAddress("THINjetRho", &THINjetRho, &b_THINjetRho);
   fChain->SetBranchAddress("THINjetNPV", &THINjetNPV, &b_THINjetNPV);
   fChain->SetBranchAddress("THINgenjetP4", &THINgenjetP4, &b_THINgenjetP4);
   fChain->SetBranchAddress("THINgenjetEM", &THINgenjetEM, &b_THINgenjetEM);
   fChain->SetBranchAddress("THINgenjetHAD", &THINgenjetHAD, &b_THINgenjetHAD);
   fChain->SetBranchAddress("THINgenjetINV", &THINgenjetINV, &b_THINgenjetINV);
   fChain->SetBranchAddress("THINgenjetAUX", &THINgenjetAUX, &b_THINgenjetAUX);
   fChain->SetBranchAddress("THINmatchedDR", &THINmatchedDR, &b_THINmatchedDR);
   fChain->SetBranchAddress("THINjetRawFactor", &THINjetRawFactor, &b_THINjetRawFactor);
   fChain->SetBranchAddress("THINunCorrJetP4", &THINunCorrJetP4, &b_THINunCorrJetP4);
   fChain->SetBranchAddress("THINjetArea", &THINjetArea, &b_THINjetArea);
   fChain->SetBranchAddress("THINjetCorrUncUp", &THINjetCorrUncUp, &b_THINjetCorrUncUp);
   fChain->SetBranchAddress("THINjetCorrUncDown", &THINjetCorrUncDown, &b_THINjetCorrUncDown);
   fChain->SetBranchAddress("THINjetCharge", &THINjetCharge, &b_THINjetCharge);
   fChain->SetBranchAddress("THINjetPartonFlavor", &THINjetPartonFlavor, &b_THINjetPartonFlavor);
   fChain->SetBranchAddress("THINjetHadronFlavor", &THINjetHadronFlavor, &b_THINjetHadronFlavor);
   fChain->SetBranchAddress("THINjetPassIDLoose", &THINjetPassIDLoose, &b_THINjetPassIDLoose);
   fChain->SetBranchAddress("THINjetPassIDTight", &THINjetPassIDTight, &b_THINjetPassIDTight);
   fChain->SetBranchAddress("THINjetCEmEF", &THINjetCEmEF, &b_THINjetCEmEF);
   fChain->SetBranchAddress("THINjetCHadEF", &THINjetCHadEF, &b_THINjetCHadEF);
   fChain->SetBranchAddress("THINjetPhoEF", &THINjetPhoEF, &b_THINjetPhoEF);
   fChain->SetBranchAddress("THINjetNEmEF", &THINjetNEmEF, &b_THINjetNEmEF);
   fChain->SetBranchAddress("THINjetNHadEF", &THINjetNHadEF, &b_THINjetNHadEF);
   fChain->SetBranchAddress("THINjetEleEF", &THINjetEleEF, &b_THINjetEleEF);
   fChain->SetBranchAddress("THINjetMuoEF", &THINjetMuoEF, &b_THINjetMuoEF);
   fChain->SetBranchAddress("THINjetCMulti", &THINjetCMulti, &b_THINjetCMulti);
   fChain->SetBranchAddress("THINjetEleMulti", &THINjetEleMulti, &b_THINjetEleMulti);
   fChain->SetBranchAddress("THINjetMuoMulti", &THINjetMuoMulti, &b_THINjetMuoMulti);
   fChain->SetBranchAddress("THINjetSSV", &THINjetSSV, &b_THINjetSSV);
   fChain->SetBranchAddress("THINjetCSV", &THINjetCSV, &b_THINjetCSV);
   fChain->SetBranchAddress("THINjetSSVHE", &THINjetSSVHE, &b_THINjetSSVHE);
   fChain->SetBranchAddress("THINjetCISVV2", &THINjetCISVV2, &b_THINjetCISVV2);
   fChain->SetBranchAddress("THINjetTCHP", &THINjetTCHP, &b_THINjetTCHP);
   fChain->SetBranchAddress("THINjetTCHE", &THINjetTCHE, &b_THINjetTCHE);
   fChain->SetBranchAddress("THINjetJP", &THINjetJP, &b_THINjetJP);
   fChain->SetBranchAddress("THINjetJBP", &THINjetJBP, &b_THINjetJBP);
   fChain->SetBranchAddress("THINPUJetID", &THINPUJetID, &b_THINPUJetID);
   fChain->SetBranchAddress("THINisPUJetIDLoose", &THINisPUJetIDLoose, &b_THINisPUJetIDLoose);
   fChain->SetBranchAddress("THINisPUJetIDMedium", &THINisPUJetIDMedium, &b_THINisPUJetIDMedium);
   fChain->SetBranchAddress("THINisPUJetIDTight", &THINisPUJetIDTight, &b_THINisPUJetIDTight);
   fChain->SetBranchAddress("FATnJet", &FATnJet, &b_FATnJet);
   fChain->SetBranchAddress("FATjetP4", &FATjetP4, &b_FATjetP4);
   fChain->SetBranchAddress("FATjetRho", &FATjetRho, &b_FATjetRho);
   fChain->SetBranchAddress("FATjetNPV", &FATjetNPV, &b_FATjetNPV);
   fChain->SetBranchAddress("FATgenjetP4", &FATgenjetP4, &b_FATgenjetP4);
   fChain->SetBranchAddress("FATgenjetEM", &FATgenjetEM, &b_FATgenjetEM);
   fChain->SetBranchAddress("FATgenjetHAD", &FATgenjetHAD, &b_FATgenjetHAD);
   fChain->SetBranchAddress("FATgenjetINV", &FATgenjetINV, &b_FATgenjetINV);
   fChain->SetBranchAddress("FATgenjetAUX", &FATgenjetAUX, &b_FATgenjetAUX);
   fChain->SetBranchAddress("FATmatchedDR", &FATmatchedDR, &b_FATmatchedDR);
   fChain->SetBranchAddress("FATjetRawFactor", &FATjetRawFactor, &b_FATjetRawFactor);
   fChain->SetBranchAddress("FATunCorrJetP4", &FATunCorrJetP4, &b_FATunCorrJetP4);
   fChain->SetBranchAddress("FATjetArea", &FATjetArea, &b_FATjetArea);
   fChain->SetBranchAddress("FATjetCorrUncUp", &FATjetCorrUncUp, &b_FATjetCorrUncUp);
   fChain->SetBranchAddress("FATjetCorrUncDown", &FATjetCorrUncDown, &b_FATjetCorrUncDown);
   fChain->SetBranchAddress("FATjetCharge", &FATjetCharge, &b_FATjetCharge);
   fChain->SetBranchAddress("FATjetPartonFlavor", &FATjetPartonFlavor, &b_FATjetPartonFlavor);
   fChain->SetBranchAddress("FATjetHadronFlavor", &FATjetHadronFlavor, &b_FATjetHadronFlavor);
   fChain->SetBranchAddress("FATjetPassIDLoose", &FATjetPassIDLoose, &b_FATjetPassIDLoose);
   fChain->SetBranchAddress("FATjetPassIDTight", &FATjetPassIDTight, &b_FATjetPassIDTight);
   fChain->SetBranchAddress("FATjetCEmEF", &FATjetCEmEF, &b_FATjetCEmEF);
   fChain->SetBranchAddress("FATjetCHadEF", &FATjetCHadEF, &b_FATjetCHadEF);
   fChain->SetBranchAddress("FATjetPhoEF", &FATjetPhoEF, &b_FATjetPhoEF);
   fChain->SetBranchAddress("FATjetNEmEF", &FATjetNEmEF, &b_FATjetNEmEF);
   fChain->SetBranchAddress("FATjetNHadEF", &FATjetNHadEF, &b_FATjetNHadEF);
   fChain->SetBranchAddress("FATjetEleEF", &FATjetEleEF, &b_FATjetEleEF);
   fChain->SetBranchAddress("FATjetMuoEF", &FATjetMuoEF, &b_FATjetMuoEF);
   fChain->SetBranchAddress("FATjetCMulti", &FATjetCMulti, &b_FATjetCMulti);
   fChain->SetBranchAddress("FATjetEleMulti", &FATjetEleMulti, &b_FATjetEleMulti);
   fChain->SetBranchAddress("FATjetMuoMulti", &FATjetMuoMulti, &b_FATjetMuoMulti);
   fChain->SetBranchAddress("FATjetSSV", &FATjetSSV, &b_FATjetSSV);
   fChain->SetBranchAddress("FATjetCSV", &FATjetCSV, &b_FATjetCSV);
   fChain->SetBranchAddress("FATjetSSVHE", &FATjetSSVHE, &b_FATjetSSVHE);
   fChain->SetBranchAddress("FATjetCISVV2", &FATjetCISVV2, &b_FATjetCISVV2);
   fChain->SetBranchAddress("FATjetTCHP", &FATjetTCHP, &b_FATjetTCHP);
   fChain->SetBranchAddress("FATjetTCHE", &FATjetTCHE, &b_FATjetTCHE);
   fChain->SetBranchAddress("FATjetJP", &FATjetJP, &b_FATjetJP);
   fChain->SetBranchAddress("FATjetJBP", &FATjetJBP, &b_FATjetJBP);
   fChain->SetBranchAddress("FATjetTau1", &FATjetTau1, &b_FATjetTau1);
   fChain->SetBranchAddress("FATjetTau2", &FATjetTau2, &b_FATjetTau2);
   fChain->SetBranchAddress("FATjetTau3", &FATjetTau3, &b_FATjetTau3);
   fChain->SetBranchAddress("FATjetTau21", &FATjetTau21, &b_FATjetTau21);
   fChain->SetBranchAddress("FATjetSDmass", &FATjetSDmass, &b_FATjetSDmass);
   fChain->SetBranchAddress("FATjetGenSDmass", &FATjetGenSDmass, &b_FATjetGenSDmass);
   fChain->SetBranchAddress("FATnSubSDJet", &FATnSubSDJet, &b_FATnSubSDJet);
   fChain->SetBranchAddress("FATsubjetSDFatJetIndex", &FATsubjetSDFatJetIndex, &b_FATsubjetSDFatJetIndex);
   fChain->SetBranchAddress("FATsubjetSDPx", &FATsubjetSDPx, &b_FATsubjetSDPx);
   fChain->SetBranchAddress("FATsubjetSDPy", &FATsubjetSDPy, &b_FATsubjetSDPy);
   fChain->SetBranchAddress("FATsubjetSDPz", &FATsubjetSDPz, &b_FATsubjetSDPz);
   fChain->SetBranchAddress("FATsubjetSDE", &FATsubjetSDE, &b_FATsubjetSDE);
   fChain->SetBranchAddress("FATsubjetSDRawFactor", &FATsubjetSDRawFactor, &b_FATsubjetSDRawFactor);
   fChain->SetBranchAddress("FATsubjetSDPartonFlavor", &FATsubjetSDPartonFlavor, &b_FATsubjetSDPartonFlavor);
   fChain->SetBranchAddress("FATsubjetSDHadronFlavor", &FATsubjetSDHadronFlavor, &b_FATsubjetSDHadronFlavor);
   fChain->SetBranchAddress("FATsubjetSDCSV", &FATsubjetSDCSV, &b_FATsubjetSDCSV);
   fChain->SetBranchAddress("FATjetPRmass", &FATjetPRmass, &b_FATjetPRmass);
   fChain->SetBranchAddress("FATjetPRmassL2L3Corr", &FATjetPRmassL2L3Corr, &b_FATjetPRmassL2L3Corr);
   fChain->SetBranchAddress("FATjetPuppiTau1", &FATjetPuppiTau1, &b_FATjetPuppiTau1);
   fChain->SetBranchAddress("FATjetPuppiTau2", &FATjetPuppiTau2, &b_FATjetPuppiTau2);
   fChain->SetBranchAddress("FATjetPuppiTau3", &FATjetPuppiTau3, &b_FATjetPuppiTau3);
   fChain->SetBranchAddress("FATjetPuppiSDmass", &FATjetPuppiSDmass, &b_FATjetPuppiSDmass);
   fChain->SetBranchAddress("FATjetPuppiP4", &FATjetPuppiP4, &b_FATjetPuppiP4);
   fChain->SetBranchAddress("FATjetPuppiSDRawP4", &FATjetPuppiSDRawP4, &b_FATjetPuppiSDRawP4);
   fChain->SetBranchAddress("FATnSubSDPuppiJet", &FATnSubSDPuppiJet, &b_FATnSubSDPuppiJet);
   fChain->SetBranchAddress("FATsubjetSDPuppiFatJetIndex", &FATsubjetSDPuppiFatJetIndex, &b_FATsubjetSDPuppiFatJetIndex);
   fChain->SetBranchAddress("FATsubjetSDPuppiPx", &FATsubjetSDPuppiPx, &b_FATsubjetSDPuppiPx);
   fChain->SetBranchAddress("FATsubjetSDPuppiPy", &FATsubjetSDPuppiPy, &b_FATsubjetSDPuppiPy);
   fChain->SetBranchAddress("FATsubjetSDPuppiPz", &FATsubjetSDPuppiPz, &b_FATsubjetSDPuppiPz);
   fChain->SetBranchAddress("FATsubjetSDPuppiE", &FATsubjetSDPuppiE, &b_FATsubjetSDPuppiE);
   fChain->SetBranchAddress("FATsubjetSDPuppiCSV", &FATsubjetSDPuppiCSV, &b_FATsubjetSDPuppiCSV);
   fChain->SetBranchAddress("FATjet_DoubleSV", &FATjet_DoubleSV, &b_FATjet_DoubleSV);
   fChain->SetBranchAddress("FATjet_nSV", &FATjet_nSV, &b_FATjet_nSV);
   fChain->SetBranchAddress("FATjet_SVMass", &FATjet_SVMass, &b_FATjet_SVMass);
   fChain->SetBranchAddress("ADDnJet", &ADDnJet, &b_ADDnJet);
   fChain->SetBranchAddress("ADDjetP4", &ADDjetP4, &b_ADDjetP4);
   fChain->SetBranchAddress("ADDjet_DoubleSV", &ADDjet_DoubleSV, &b_ADDjet_DoubleSV);
   fChain->SetBranchAddress("ADDjet_nSV", &ADDjet_nSV, &b_ADDjet_nSV);
   fChain->SetBranchAddress("ADDjet_SVMass", &ADDjet_SVMass, &b_ADDjet_SVMass);
   fChain->SetBranchAddress("AK8PuppinJet", &AK8PuppinJet, &b_AK8PuppinJet);
   fChain->SetBranchAddress("AK8PuppijetP4", &AK8PuppijetP4, &b_AK8PuppijetP4);
   fChain->SetBranchAddress("AK8PuppijetRho", &AK8PuppijetRho, &b_AK8PuppijetRho);
   fChain->SetBranchAddress("AK8PuppijetNPV", &AK8PuppijetNPV, &b_AK8PuppijetNPV);
   fChain->SetBranchAddress("AK8PuppigenjetP4", &AK8PuppigenjetP4, &b_AK8PuppigenjetP4);
   fChain->SetBranchAddress("AK8PuppigenjetEM", &AK8PuppigenjetEM, &b_AK8PuppigenjetEM);
   fChain->SetBranchAddress("AK8PuppigenjetHAD", &AK8PuppigenjetHAD, &b_AK8PuppigenjetHAD);
   fChain->SetBranchAddress("AK8PuppigenjetINV", &AK8PuppigenjetINV, &b_AK8PuppigenjetINV);
   fChain->SetBranchAddress("AK8PuppigenjetAUX", &AK8PuppigenjetAUX, &b_AK8PuppigenjetAUX);
   fChain->SetBranchAddress("AK8PuppimatchedDR", &AK8PuppimatchedDR, &b_AK8PuppimatchedDR);
   fChain->SetBranchAddress("AK8PuppijetRawFactor", &AK8PuppijetRawFactor, &b_AK8PuppijetRawFactor);
   fChain->SetBranchAddress("AK8PuppiunCorrJetP4", &AK8PuppiunCorrJetP4, &b_AK8PuppiunCorrJetP4);
   fChain->SetBranchAddress("AK8PuppijetArea", &AK8PuppijetArea, &b_AK8PuppijetArea);
   fChain->SetBranchAddress("AK8PuppijetCorrUncUp", &AK8PuppijetCorrUncUp, &b_AK8PuppijetCorrUncUp);
   fChain->SetBranchAddress("AK8PuppijetCorrUncDown", &AK8PuppijetCorrUncDown, &b_AK8PuppijetCorrUncDown);
   fChain->SetBranchAddress("AK8PuppijetCharge", &AK8PuppijetCharge, &b_AK8PuppijetCharge);
   fChain->SetBranchAddress("AK8PuppijetPartonFlavor", &AK8PuppijetPartonFlavor, &b_AK8PuppijetPartonFlavor);
   fChain->SetBranchAddress("AK8PuppijetHadronFlavor", &AK8PuppijetHadronFlavor, &b_AK8PuppijetHadronFlavor);
   fChain->SetBranchAddress("AK8PuppijetPassIDLoose", &AK8PuppijetPassIDLoose, &b_AK8PuppijetPassIDLoose);
   fChain->SetBranchAddress("AK8PuppijetPassIDTight", &AK8PuppijetPassIDTight, &b_AK8PuppijetPassIDTight);
   fChain->SetBranchAddress("AK8PuppijetCEmEF", &AK8PuppijetCEmEF, &b_AK8PuppijetCEmEF);
   fChain->SetBranchAddress("AK8PuppijetCHadEF", &AK8PuppijetCHadEF, &b_AK8PuppijetCHadEF);
   fChain->SetBranchAddress("AK8PuppijetPhoEF", &AK8PuppijetPhoEF, &b_AK8PuppijetPhoEF);
   fChain->SetBranchAddress("AK8PuppijetNEmEF", &AK8PuppijetNEmEF, &b_AK8PuppijetNEmEF);
   fChain->SetBranchAddress("AK8PuppijetNHadEF", &AK8PuppijetNHadEF, &b_AK8PuppijetNHadEF);
   fChain->SetBranchAddress("AK8PuppijetEleEF", &AK8PuppijetEleEF, &b_AK8PuppijetEleEF);
   fChain->SetBranchAddress("AK8PuppijetMuoEF", &AK8PuppijetMuoEF, &b_AK8PuppijetMuoEF);
   fChain->SetBranchAddress("AK8PuppijetCMulti", &AK8PuppijetCMulti, &b_AK8PuppijetCMulti);
   fChain->SetBranchAddress("AK8PuppijetEleMulti", &AK8PuppijetEleMulti, &b_AK8PuppijetEleMulti);
   fChain->SetBranchAddress("AK8PuppijetMuoMulti", &AK8PuppijetMuoMulti, &b_AK8PuppijetMuoMulti);
   fChain->SetBranchAddress("AK8PuppijetSSV", &AK8PuppijetSSV, &b_AK8PuppijetSSV);
   fChain->SetBranchAddress("AK8PuppijetCSV", &AK8PuppijetCSV, &b_AK8PuppijetCSV);
   fChain->SetBranchAddress("AK8PuppijetSSVHE", &AK8PuppijetSSVHE, &b_AK8PuppijetSSVHE);
   fChain->SetBranchAddress("AK8PuppijetCISVV2", &AK8PuppijetCISVV2, &b_AK8PuppijetCISVV2);
   fChain->SetBranchAddress("AK8PuppijetTCHP", &AK8PuppijetTCHP, &b_AK8PuppijetTCHP);
   fChain->SetBranchAddress("AK8PuppijetTCHE", &AK8PuppijetTCHE, &b_AK8PuppijetTCHE);
   fChain->SetBranchAddress("AK8PuppijetJP", &AK8PuppijetJP, &b_AK8PuppijetJP);
   fChain->SetBranchAddress("AK8PuppijetJBP", &AK8PuppijetJBP, &b_AK8PuppijetJBP);
   fChain->SetBranchAddress("AK8PuppijetTau1", &AK8PuppijetTau1, &b_AK8PuppijetTau1);
   fChain->SetBranchAddress("AK8PuppijetTau2", &AK8PuppijetTau2, &b_AK8PuppijetTau2);
   fChain->SetBranchAddress("AK8PuppijetTau3", &AK8PuppijetTau3, &b_AK8PuppijetTau3);
   fChain->SetBranchAddress("AK8PuppijetTau21", &AK8PuppijetTau21, &b_AK8PuppijetTau21);
   fChain->SetBranchAddress("AK8PuppijetSDmass", &AK8PuppijetSDmass, &b_AK8PuppijetSDmass);
   fChain->SetBranchAddress("AK8PuppijetGenSDmass", &AK8PuppijetGenSDmass, &b_AK8PuppijetGenSDmass);
   fChain->SetBranchAddress("AK8PuppinSubSDJet", &AK8PuppinSubSDJet, &b_AK8PuppinSubSDJet);
   fChain->SetBranchAddress("AK8PuppisubjetSDFatJetIndex", &AK8PuppisubjetSDFatJetIndex, &b_AK8PuppisubjetSDFatJetIndex);
   fChain->SetBranchAddress("AK8PuppisubjetSDPx", &AK8PuppisubjetSDPx, &b_AK8PuppisubjetSDPx);
   fChain->SetBranchAddress("AK8PuppisubjetSDPy", &AK8PuppisubjetSDPy, &b_AK8PuppisubjetSDPy);
   fChain->SetBranchAddress("AK8PuppisubjetSDPz", &AK8PuppisubjetSDPz, &b_AK8PuppisubjetSDPz);
   fChain->SetBranchAddress("AK8PuppisubjetSDE", &AK8PuppisubjetSDE, &b_AK8PuppisubjetSDE);
   fChain->SetBranchAddress("AK8PuppisubjetSDRawFactor", &AK8PuppisubjetSDRawFactor, &b_AK8PuppisubjetSDRawFactor);
   fChain->SetBranchAddress("AK8PuppisubjetSDPartonFlavor", &AK8PuppisubjetSDPartonFlavor, &b_AK8PuppisubjetSDPartonFlavor);
   fChain->SetBranchAddress("AK8PuppisubjetSDHadronFlavor", &AK8PuppisubjetSDHadronFlavor, &b_AK8PuppisubjetSDHadronFlavor);
   fChain->SetBranchAddress("AK8PuppisubjetSDCSV", &AK8PuppisubjetSDCSV, &b_AK8PuppisubjetSDCSV);
   fChain->SetBranchAddress("AK8Puppijet_DoubleSV", &AK8Puppijet_DoubleSV, &b_AK8Puppijet_DoubleSV);
   fChain->SetBranchAddress("AK8Puppijet_nSV", &AK8Puppijet_nSV, &b_AK8Puppijet_nSV);
   fChain->SetBranchAddress("AK8Puppijet_SVMass", &AK8Puppijet_SVMass, &b_AK8Puppijet_SVMass);
   fChain->SetBranchAddress("CA15PuppinJet", &CA15PuppinJet, &b_CA15PuppinJet);
   fChain->SetBranchAddress("CA15PuppijetP4", &CA15PuppijetP4, &b_CA15PuppijetP4);
   fChain->SetBranchAddress("CA15PuppijetRho", &CA15PuppijetRho, &b_CA15PuppijetRho);
   fChain->SetBranchAddress("CA15PuppijetNPV", &CA15PuppijetNPV, &b_CA15PuppijetNPV);
   fChain->SetBranchAddress("CA15PuppigenjetP4", &CA15PuppigenjetP4, &b_CA15PuppigenjetP4);
   fChain->SetBranchAddress("CA15PuppigenjetEM", &CA15PuppigenjetEM, &b_CA15PuppigenjetEM);
   fChain->SetBranchAddress("CA15PuppigenjetHAD", &CA15PuppigenjetHAD, &b_CA15PuppigenjetHAD);
   fChain->SetBranchAddress("CA15PuppigenjetINV", &CA15PuppigenjetINV, &b_CA15PuppigenjetINV);
   fChain->SetBranchAddress("CA15PuppigenjetAUX", &CA15PuppigenjetAUX, &b_CA15PuppigenjetAUX);
   fChain->SetBranchAddress("CA15PuppimatchedDR", &CA15PuppimatchedDR, &b_CA15PuppimatchedDR);
   fChain->SetBranchAddress("CA15PuppijetRawFactor", &CA15PuppijetRawFactor, &b_CA15PuppijetRawFactor);
   fChain->SetBranchAddress("CA15PuppiunCorrJetP4", &CA15PuppiunCorrJetP4, &b_CA15PuppiunCorrJetP4);
   fChain->SetBranchAddress("CA15PuppijetArea", &CA15PuppijetArea, &b_CA15PuppijetArea);
   fChain->SetBranchAddress("CA15PuppijetCorrUncUp", &CA15PuppijetCorrUncUp, &b_CA15PuppijetCorrUncUp);
   fChain->SetBranchAddress("CA15PuppijetCorrUncDown", &CA15PuppijetCorrUncDown, &b_CA15PuppijetCorrUncDown);
   fChain->SetBranchAddress("CA15PuppijetCharge", &CA15PuppijetCharge, &b_CA15PuppijetCharge);
   fChain->SetBranchAddress("CA15PuppijetPartonFlavor", &CA15PuppijetPartonFlavor, &b_CA15PuppijetPartonFlavor);
   fChain->SetBranchAddress("CA15PuppijetHadronFlavor", &CA15PuppijetHadronFlavor, &b_CA15PuppijetHadronFlavor);
   fChain->SetBranchAddress("CA15PuppijetPassIDLoose", &CA15PuppijetPassIDLoose, &b_CA15PuppijetPassIDLoose);
   fChain->SetBranchAddress("CA15PuppijetPassIDTight", &CA15PuppijetPassIDTight, &b_CA15PuppijetPassIDTight);
   fChain->SetBranchAddress("CA15PuppijetCEmEF", &CA15PuppijetCEmEF, &b_CA15PuppijetCEmEF);
   fChain->SetBranchAddress("CA15PuppijetCHadEF", &CA15PuppijetCHadEF, &b_CA15PuppijetCHadEF);
   fChain->SetBranchAddress("CA15PuppijetPhoEF", &CA15PuppijetPhoEF, &b_CA15PuppijetPhoEF);
   fChain->SetBranchAddress("CA15PuppijetNEmEF", &CA15PuppijetNEmEF, &b_CA15PuppijetNEmEF);
   fChain->SetBranchAddress("CA15PuppijetNHadEF", &CA15PuppijetNHadEF, &b_CA15PuppijetNHadEF);
   fChain->SetBranchAddress("CA15PuppijetEleEF", &CA15PuppijetEleEF, &b_CA15PuppijetEleEF);
   fChain->SetBranchAddress("CA15PuppijetMuoEF", &CA15PuppijetMuoEF, &b_CA15PuppijetMuoEF);
   fChain->SetBranchAddress("CA15PuppijetCMulti", &CA15PuppijetCMulti, &b_CA15PuppijetCMulti);
   fChain->SetBranchAddress("CA15PuppijetEleMulti", &CA15PuppijetEleMulti, &b_CA15PuppijetEleMulti);
   fChain->SetBranchAddress("CA15PuppijetMuoMulti", &CA15PuppijetMuoMulti, &b_CA15PuppijetMuoMulti);
   fChain->SetBranchAddress("CA15PuppijetSSV", &CA15PuppijetSSV, &b_CA15PuppijetSSV);
   fChain->SetBranchAddress("CA15PuppijetCSV", &CA15PuppijetCSV, &b_CA15PuppijetCSV);
   fChain->SetBranchAddress("CA15PuppijetSSVHE", &CA15PuppijetSSVHE, &b_CA15PuppijetSSVHE);
   fChain->SetBranchAddress("CA15PuppijetCISVV2", &CA15PuppijetCISVV2, &b_CA15PuppijetCISVV2);
   fChain->SetBranchAddress("CA15PuppijetTCHP", &CA15PuppijetTCHP, &b_CA15PuppijetTCHP);
   fChain->SetBranchAddress("CA15PuppijetTCHE", &CA15PuppijetTCHE, &b_CA15PuppijetTCHE);
   fChain->SetBranchAddress("CA15PuppijetJP", &CA15PuppijetJP, &b_CA15PuppijetJP);
   fChain->SetBranchAddress("CA15PuppijetJBP", &CA15PuppijetJBP, &b_CA15PuppijetJBP);
   fChain->SetBranchAddress("CA15PuppijetTau1", &CA15PuppijetTau1, &b_CA15PuppijetTau1);
   fChain->SetBranchAddress("CA15PuppijetTau2", &CA15PuppijetTau2, &b_CA15PuppijetTau2);
   fChain->SetBranchAddress("CA15PuppijetTau3", &CA15PuppijetTau3, &b_CA15PuppijetTau3);
   fChain->SetBranchAddress("CA15PuppijetTau21", &CA15PuppijetTau21, &b_CA15PuppijetTau21);
   fChain->SetBranchAddress("CA15PuppijetSDmass", &CA15PuppijetSDmass, &b_CA15PuppijetSDmass);
   fChain->SetBranchAddress("CA15PuppijetGenSDmass", &CA15PuppijetGenSDmass, &b_CA15PuppijetGenSDmass);
   fChain->SetBranchAddress("CA15PuppinSubSDJet", &CA15PuppinSubSDJet, &b_CA15PuppinSubSDJet);
   fChain->SetBranchAddress("CA15PuppisubjetSDFatJetIndex", &CA15PuppisubjetSDFatJetIndex, &b_CA15PuppisubjetSDFatJetIndex);
   fChain->SetBranchAddress("CA15PuppisubjetSDPx", &CA15PuppisubjetSDPx, &b_CA15PuppisubjetSDPx);
   fChain->SetBranchAddress("CA15PuppisubjetSDPy", &CA15PuppisubjetSDPy, &b_CA15PuppisubjetSDPy);
   fChain->SetBranchAddress("CA15PuppisubjetSDPz", &CA15PuppisubjetSDPz, &b_CA15PuppisubjetSDPz);
   fChain->SetBranchAddress("CA15PuppisubjetSDE", &CA15PuppisubjetSDE, &b_CA15PuppisubjetSDE);
   fChain->SetBranchAddress("CA15PuppisubjetSDRawFactor", &CA15PuppisubjetSDRawFactor, &b_CA15PuppisubjetSDRawFactor);
   fChain->SetBranchAddress("CA15PuppisubjetSDPartonFlavor", &CA15PuppisubjetSDPartonFlavor, &b_CA15PuppisubjetSDPartonFlavor);
   fChain->SetBranchAddress("CA15PuppisubjetSDHadronFlavor", &CA15PuppisubjetSDHadronFlavor, &b_CA15PuppisubjetSDHadronFlavor);
   fChain->SetBranchAddress("CA15PuppisubjetSDCSV", &CA15PuppisubjetSDCSV, &b_CA15PuppisubjetSDCSV);
   fChain->SetBranchAddress("CA15Puppi_doublebtag", &CA15Puppi_doublebtag, &b_CA15Puppi_doublebtag);
   fChain->SetBranchAddress("CA15PuppiECF_2_3_10", &CA15PuppiECF_2_3_10, &b_CA15PuppiECF_2_3_10);
   fChain->SetBranchAddress("CA15PuppiECF_1_2_10", &CA15PuppiECF_1_2_10, &b_CA15PuppiECF_1_2_10);
   fChain->SetBranchAddress("CA15Puppijet_DoubleSV", &CA15Puppijet_DoubleSV, &b_CA15Puppijet_DoubleSV);
   fChain->SetBranchAddress("CA15Puppijet_nSV", &CA15Puppijet_nSV, &b_CA15Puppijet_nSV);
   fChain->SetBranchAddress("CA15Puppijet_SVMass", &CA15Puppijet_SVMass, &b_CA15Puppijet_SVMass);
   Notify();
}

Bool_t treeMaker::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void treeMaker::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t treeMaker::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef treeMaker_cxx
