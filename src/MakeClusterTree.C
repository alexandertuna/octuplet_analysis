///
///  \file   MakeClusterTree.C
///
///  \author Christopher Rogan
///          (crogan@cern.ch)
///
///  \date   2016 Oct
///

#include <iostream>
#include <stdlib.h>

#include "include/PDOToCharge.hh"
#include "include/TDOToTime.hh"
//#include "include/MMDataAnalysis.hh"
#include "include/MMPacmanAlgo.hh"
#include "include/GeoOctuplet.hh"
#include "include/MMDataAnalysis.hh"
#include "include/SimpleTrackFitter.hh"
#include "include/ScintillatorClusterFilterer.hh"

using namespace std;

int main(int argc, char* argv[]){

  char inputFileName[400];
  char outputFileName[400];
  char PDOFileName[400];
  char TDOFileName[400];
  
  if ( argc < 5 ){
    cout << "Error at Input: please specify input/output .root files ";
    cout << " and (optional) PDO/TDO calibration files" << endl;
    cout << "Example:   ./MakeClusterTree.x -i input.root -o output.root" << endl;
    cout << "Example:   ./MakeClusterTree.x -i input.root -o output.root";
    cout << " -p PDOcalib.root -t TDOcalib.root" << endl;
    return 0;
  }

  bool b_input = false;
  bool b_out   = false;
  bool b_pdo   = false;
  bool b_tdo   = false;
  for (int i=1;i<argc-1;i++){
    if (strncmp(argv[i],"-i",2)==0){
      sscanf(argv[i+1],"%s", inputFileName);
      b_input = true;
    }
    if (strncmp(argv[i],"-o",2)==0){
      sscanf(argv[i+1],"%s", outputFileName);
      b_out = true;
    }
    if (strncmp(argv[i],"-p",2)==0){
      sscanf(argv[i+1],"%s", PDOFileName);
      b_pdo = true;
    }
    if (strncmp(argv[i],"-t",2)==0){
      sscanf(argv[i+1],"%s", TDOFileName);
      b_tdo = true;
    }
  }

  if(!b_input){
    cout << "Error at Input: please specify input file (-i flag)" << endl;
    return 0;
  }

  if(!b_out){
    cout << "Error at Input: please specify output file (-o flag)" << endl;
    return 0;
  }

  PDOToCharge* PDOCalibrator;
  if(b_pdo)
    PDOCalibrator = new PDOToCharge(PDOFileName);
  else
    PDOCalibrator = new PDOToCharge();

  TDOToTime* TDOCalibrator;
  if(b_tdo)
    TDOCalibrator = new TDOToTime(TDOFileName);
  else
    TDOCalibrator = new TDOToTime();

  MMPacmanAlgo* PACMAN = new MMPacmanAlgo(5,2.,0.5);

  ScintillatorClusterFilterer* FILTERER = new ScintillatorClusterFilterer();

  GeoOctuplet* GEOMETRY = new GeoOctuplet();

  MMDataAnalysis* DATA;
  TFile* f = new TFile(inputFileName, "READ");
  if(!f){
    cout << "Error: unable to open input file " << inputFileName << endl;
    return false;
  }
  TTree* T = (TTree*) f->Get("COMB_data");
  if(!T){
    cout << "Error: cannot find tree COMB_data in " << inputFileName << endl;
    return false;
  }

  DATA = (MMDataAnalysis*) new MMDataAnalysis(T);
  DATA->SetTP(0);

  int Nevent = DATA->GetNEntries();
  int debug  = 0;
  
  TFile* fout = new TFile(outputFileName, "RECREATE");
  fout->cd();
  int N_clus;
  vector<int> clus_MMFE8;
  vector<int> clus_Index;
  vector<double> clus_Charge;
  vector<double> clus_Time;
  vector<double> clus_Channel;

  TTree* cluster_tree = new TTree("ClusterTree",
				  "ClusterTree");
  cluster_tree->Branch("N_clus", &N_clus);
  cluster_tree->Branch("clus_MMFE8", &clus_MMFE8);
  cluster_tree->Branch("clus_Index", &clus_Index);
  cluster_tree->Branch("clus_Charge", &clus_Charge);
  cluster_tree->Branch("clus_Time", &clus_Time);
  cluster_tree->Branch("clus_Channel", &clus_Channel);

  for(int evt = 0; evt < Nevent; evt++){
       DATA->GetEntry(evt);
    if(evt%10000 == 0)
      cout << "Processing event # " << evt << " | " << Nevent << endl;

    if(GEOMETRY->RunNumber() < 0){
      GEOMETRY->SetRunNumber(DATA->RunNum);
    }

    if(!DATA->sc_EventHits.IsGoodEvent())
      continue;
    
    // Calibrate PDO -> Charge
    PDOCalibrator->Calibrate(DATA->mm_EventHits);
    // Calibrate TDO -> Time
    TDOCalibrator->Calibrate(DATA->mm_EventHits);
  
    // initialize PACMAN info for this event
    PACMAN->SetEventTrigBCID(DATA->mm_trig_BCID);
    PACMAN->SetEventPadTime(0); // add this

    // initialize cluster maker
    FILTERER->SetRunNumber(DATA->RunNum);
    
    int Nboards = DATA->mm_EventHits.GetNBoards();
    
    vector<MMClusterList> all_clusters;
    for(int i = 0; i < Nboards; i++){
      if(DATA->mm_EventHits[i].GetNHits() == 0)
	continue;
      
      MMClusterList clusters = PACMAN->Cluster(DATA->mm_EventHits[i]);
      if(clusters.GetNCluster() > 0)
	all_clusters.push_back(clusters);
    }

    // preselection
    if (all_clusters.size() < 4)
      continue;

    N_clus = 0;
    clus_MMFE8.clear();
    clus_Index.clear();
    clus_Charge.clear();
    clus_Time.clear();
    clus_Channel.clear();

    MMClusterList clusters_all;
    MMClusterList clusters_fit;

    // flatten the list of lists
    for (auto clus_list: all_clusters)
      for (auto clus: clus_list)
        clusters_all.AddCluster(*clus);

    // filter clusters with scintillator roads
    // no need for sum(res2) selection since roads are already restrictive
    for (auto botpair: DATA->sc_EventHits.GetBotPair()){
      clusters_fit = FILTERER->FilterClustersScint(clusters_all, *GEOMETRY, botpair.first->Channel(), DATA->mm_EventNum, debug);

      N_clus = clusters_fit.GetNCluster();
      for (auto clus: clusters_fit){
        clus_MMFE8  .push_back(clus->MMFE8());
        clus_Index  .push_back(clus->MMFE8Index());
        clus_Charge .push_back(clus->Charge());
        clus_Time   .push_back(clus->Time());
        clus_Channel.push_back(clus->Channel());
      }
    }

    if (clusters_fit.GetNCluster() < 8)
      continue;

    cluster_tree->Fill();

  } // end event loop

  fout->cd();
  cluster_tree->Write();
  fout->Close();
    
}
