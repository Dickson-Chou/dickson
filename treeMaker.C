#define treeMaker_cxx
#include "treeMaker.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <iostream>

using namespace std;

void treeMaker::Loop(string outputFile)
{
//   In a ROOT session, you can do:
//      root> .L treeMaker.C
//      root> treeMaker t
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
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;
//   treeMaker::Loop(string outputFile);
  
//   string outputFile = "_filtered1.root";

  

// get the TTree from the file 
//
   TH1F* hm = new TH1F("hm"," ",400,0,400 );
   TH1F* heve = new TH1F("heve", "" , 1 , 0.5 , 1.5);
   TH1F* FATjetPt = new TH1F("FATpt","",400,0,400 );
   TH1F* FATjetEta = new TH1F("FATEta","",8 , -4 , 4);

   // clone tree 
 /*  TFile* newfile_data = new TFile(outputFile.data(),"recreate");                                                                
   TTree* newtree = fChain->CloneTree(0);
   newtree->SetMaxTreeSize(5000000000); */
   //cout << "Saving "  << outputFile       << " tree" << endl;
   
   // now open new root file                                                    
   TFile* newfile_data = new TFile(outputFile.data(),"recreate");
   TTree* data = new TTree("data", "data");
   TClonesArray* ptrFATJetP4 = new TClonesArray("TLorentzVector");
   TClonesArray &FATJetP4 = *ptrFATJetP4;
   ptrFATJetP4->BypassStreamer();
   Float_t FATjetTau21;
   data->Branch("FATJetP4",&FATJetP4);
   data->Branch("FATjetTau21",&FATjetTau21);
   

//   Float_t ;
//   data->Branch("",)


   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nPassEvt=0;
   Long64_t nPassTrig=0;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   //for (Long64_t jentry=0; jentry<100;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
//        cout << nVtx << endl;
      if (Cut(ientry) < 0) continue;
      if (jentry%100==0) 	printf("%lld/%lld, %4.1f%% done.\r",jentry,nentries,(float)jentry/(float)nentries*100.);      
      heve->Fill(1);
      // this example filter keeps events with at least two electrons 
      // or two muons that form an invariant mass between 60 and 120 GeV and 
      // with pt > 60 GeV

      //trigger
      bool passTrigger=false;
      for(int it=0; it<hlt_nTrigs; it++){
         std::string trigName = hlt_trigName->at(it);
         bool result = hlt_trigResult->at(it);
         if ((trigName.find("HLT_PFHT800_v")!= std::string::npos && result==1 )||
            (trigName.find("HLT_PFHT900_v")!= std::string::npos && result==1) ||
            (trigName.find("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v")!= std::string::npos && result==1) ||
            (trigName.find("HLT_AK8PFJet360_TrimMass30_v")!= std::string::npos && result==1) ||
            (trigName.find("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_v")!= std::string::npos && result==1) ||
            (trigName.find("HLT_AK8PFHT650_TrimR0p1PT0p03Mass50_v")!= std::string::npos && result==1) ||
            (trigName.find("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v")!= std::string::npos && result==1) 
            
            )
         {
            passTrigger=true;
            break;
         }
    }
         if (!passTrigger)continue;
         nPassTrig++;

      int HIndex[2]={-1,-1};
      unsigned int nGoodJets=0;
 //     cout << FATnJet << "\t"  << endl;
 //     cout<< FATjetTau21->size() << endl;
      

      for (int ij=0 ; ij < FATnJet; ij++)  {
 //        if(FATjetTau21[ij].size()==0)continue;
		   TLorentzVector* FATJet_ij = (TLorentzVector*)FATjetP4->At(ij);
		   if(FATJet_ij->Pt()<300)continue;
		   if(fabs(FATJet_ij->Eta())>2.4)continue;
         if(FATjetTau21<0.55)continue;  
         //ptrFATJet = FATJet_ij;

         //cout <<  FATnJet << endl;
         for (int i=0 ; i < ij ; i++){
            TLorentzVector* FATJet_i = (TLorentzVector*)FATjetP4->At(i);
            if(FATJet_i->Pt()<300)continue;
            if(fabs(FATJet_i->Eta())>2.4)continue;
            if(FATjetTau21<0.55)continue;
            if(fabs ( ( FATJet_ij->Eta () )-( FATJet_i->Eta () ))>1.3)continue;
            FATjetPt->Fill(FATJet_ij->Pt());
            FATjetEta->Fill(FATJet_ij->Eta());
            new(FATJetP4[0]) TLorentzVector(*FATJet_ij);
            new(FATJetP4[1]) TLorentzVector(*FATJet_i);
            
            HIndex[0] = ij;
            HIndex[1] = i; 
            nGoodJets++;
            FATJetP4.Clear();
            break;
            
	     }
      }

      if(nGoodJets<2)continue;
      nPassEvt++;
//      newtree->Fill();
      nPassEvt++;
      data->Fill();
   } // end of loop over events
  // TCanvas *c1 = new TCanvas("c1","c1", 900, 600);
 // c1->cd();



   //newtree->Print();
   heve->Write();
   FATjetPt->Write();
   FATjetEta->Write();                                                          
//   newtree->AutoSave();
   newfile_data->Write();
  // c1.update();
  // delete newfile_data;
   


   
   cout << "nentries = " << nentries << endl;
   cout << "Number of passed events = " << nPassEvt << endl;
   cout << "Number of events passing triggers =" << nPassTrig << endl;
   cout << "Reduction rate = " << (double)nPassEvt/(double)nentries << endl;


}
