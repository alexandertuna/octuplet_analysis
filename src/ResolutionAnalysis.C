///
///  \file   ResolutionAnalysis.C
///
///  \author Tunaface
///          (tuna@cern.ch)
///
///  \date   2017 Jan
///

#include "TH1D.h"
#include "TH2D.h"
#include "TMultiGraph.h"
#include "TLine.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TF1.h"
#include <iostream>
#include <algorithm>
#include <stdlib.h>
#include <numeric>

#include "include/MMPlot.hh"
#include "include/PDOToCharge.hh"
#include "include/TDOToTime.hh"
#include "include/MMDataAnalysis.hh"
#include "include/MMClusterAlgo.hh"
#include "include/MMPacmanAlgo.hh"
#include "include/GeoOctuplet.hh"
#include "include/GeoPlane.hh"
#include "include/SimpleTrackFitter.hh"
#include "include/ScintillatorClusterFilterer.hh"

using namespace std;

void progress(double time_diff, int nprocessed, int ntotal){
  double rate = (double)(nprocessed+1)/time_diff;
  std::cout.precision(1);
  std::cout << "\r > " << nprocessed << " / " << ntotal 
            << " | "   << std::fixed << 100*(double)(nprocessed)/(double)(ntotal) << "%"
            << " | "   << std::fixed << rate << "Hz"
            << " | "   << std::fixed << time_diff/60 << "m elapsed"
            << " | "   << std::fixed << (double)(ntotal-nprocessed)/(rate*60) << "m remaining    "
            << std::flush;
  std::cout.precision(6);
}

double xpos(TPHit* hit, GeoOctuplet* geo){
  return geo->Get(hit->MMFE8Index()).LocalXatYbegin(hit->Channel())
    + geo->Get(hit->MMFE8Index()).Origin().X();
}

double xpos_end(TPHit* hit, GeoOctuplet* geo){
  return geo->Get(hit->MMFE8Index()).LocalXatYend(hit->Channel())
    + geo->Get(hit->MMFE8Index()).Origin().X();
}

double xpos(MMCluster* hit, GeoOctuplet* geo){
  return geo->Get(hit->MMFE8Index()).LocalXatYbegin(hit->Channel())
    + geo->Get(hit->MMFE8Index()).Origin().X();
}

double xpos(MMCluster hit, GeoOctuplet* geo){
  return geo->Get(hit.MMFE8Index()).LocalXatYbegin(hit.Channel())
    + geo->Get(hit.MMFE8Index()).Origin().X();
}

double zpos(TPHit* hit, GeoOctuplet* geo){
  return geo->Get(geo->Index(hit->MMFE8())).Origin().Z();
}

double zpos(MMCluster* hit, GeoOctuplet* geo){
  return geo->Get(geo->Index(hit->MMFE8())).Origin().Z();
}

double channel_from_x(double xpos, int board, GeoOctuplet* geo, int begin){
  int SignChannel = geo->Get(board).SignChannel();
  double xorigin  = geo->Get(board).Origin().X();
  double alpha    = (begin) ? geo->Get(board).StripAlpha() : 0.0;
  double channel  = ((xpos - xorigin - tan(alpha)*200) / (SignChannel*0.4)) + 256.5;
  return channel;
}

double calculate_slope(std::vector<double> xs, std::vector<double> zs){

  double slope  = 0.0;
  double c      = 0.0;
  double sum_z2 = 0.0;
  double sum_z  = std::accumulate(zs.begin(), zs.end(), 0.0);
  double avg_z  = sum_z / (double)(zs.size());
  for (auto z: zs)
    sum_z2 += (z*z);
  double denom = sum_z2 - (double)(zs.size())*(avg_z*avg_z);

  for (unsigned int it = 0; it < xs.size(); it++){
    c = (zs[it] - avg_z) / denom;
    slope += (c*xs[it]);
  }

  return slope;
}

int satisfy_road(std::vector<int> boards, std::vector<int> vmms, std::vector<int> chs, 
                  int roadsize=64, int neighbor_up=1, int neighbor_dn=1){

  std::vector<int> marked  = { 0,  0,  0,  0,  0,  0,  0,  0};
  std::vector<int> offsets = {64, 64, 58, 71, 58, 71, 64, 64};
  std::vector<int> flipped = {0, 3, 5, 6};
  int nroads  = 512/roadsize * 2;
  int nch     = (int)(chs.size());
  int ch_tp   = 0;
  int bo      = 0;
  std::vector<int> chs_tp = {};

  // convert to TP strips
  for (int i = 0; i < nch; i++){
    ch_tp = 64*vmms[i] + chs[i]-1;
    if (std::find(flipped.begin(), flipped.end(), boards[i]) != flipped.end())
      ch_tp = 511 - ch_tp;
    ch_tp += offsets[bo];
    chs_tp.push_back(ch_tp);
  }

  // check roads
  for (int road = 0; road < nroads; road++){
    marked = {0, 0, 0, 0, 0, 0, 0, 0};
    for (int i = 0; i < nch; i++){
      bo = boards[i];
      if                     (chs_tp[i] / roadsize == road)   marked[bo]++;
      else if (neighbor_up && chs_tp[i] / roadsize == road+1) marked[bo]++;
      else if (neighbor_dn && chs_tp[i] / roadsize == road-1) marked[bo]++;
    }
    if (marked[0] + marked[1] >= 1 &&
        marked[6] + marked[7] >= 1 &&
        marked[2] + marked[3] + marked[4] + marked[5] >= 2)
      return 1;
  }
  return 0;
}

TH1D* create_track_vs_time(MMDataAnalysis* DATA){
  double start_time = 0.0;
  double end_time   = 0.0;

  DATA->GetEntry(0);
  start_time = DATA->mm_EventHits.Time();
  DATA->GetEntry(DATA->GetNEntries() - 1);
  end_time = DATA->mm_EventHits.Time();

  int minutes = (int)((end_time - start_time)/60.0);
  TH1D* hist = new TH1D("track_vs_time", ";Time since start of run [minutes];", minutes, 0, minutes);

  return hist;
}

int main(int argc, char* argv[]){

  char inputFileName[400];
  char outputFileName[400];
  char PDOFileName[400];
  char TDOFileName[400];
  char AlignFileName[400];
  
  if ( argc < 5 ){
    cout << "Error at Input: please specify input/output .root files ";
    cout << " and (optional) PDO/TDO calibration files" << endl;
    cout << "Example:   ./RunMMAnalysis.x -i input.root -o output.root" << endl;
    cout << "Example:   ./RunMMAnalysis.x -i input.root -o output.root";
    cout << " -p PDOcalib.root -t TDOcalib.root" << endl;
    cout << "Other options:" << endl;
    cout << "   -p PDOcalib.root : PDO calibration file" << endl;
    cout << "   -t TDOcalib.root : TDO calibration file" << endl;
    cout << "   -a alignment.root : alignment file" << endl;
    return 0;
  }

  bool b_input = false;
  bool b_out   = false;
  bool b_pdo   = false;
  bool b_tdo   = false;
  bool b_align = false;
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
    if (strncmp(argv[i],"-a",2)==0){
      sscanf(argv[i+1],"%s", AlignFileName);
      b_align = true;
    }
  }

  if(!b_input) std::cout << "Error at Input: please specify  input file (-i flag)" << std::endl;
  if(!b_out)   std::cout << "Error at Input: please specify output file (-o flag)" << std::endl;
  if(!b_input || !b_out) return 0;

  const int nboards = 8;
  int nboardshit = 0;
  int i = 0, ibo = 0, ich = 0, test = 0;
  int debug = 0;
  int hella_clusters = 0;
  int do_trigger = 0;
  int EventNum4Hist = 0;

  // class defs
  PDOToCharge* PDOCalibrator;
  TDOToTime*   TDOCalibrator;
  MMPacmanAlgo*                PACMAN       = new MMPacmanAlgo(5,2.,0.5);
  GeoOctuplet*                 GEOMETRY     = new GeoOctuplet();
  SimpleTrackFitter*           FITTER       = new SimpleTrackFitter();
  ScintillatorClusterFilterer* FILTERER     = new ScintillatorClusterFilterer();
  MMDataAnalysis*    DATA;
  if(b_pdo) PDOCalibrator = new PDOToCharge(PDOFileName);
  else      PDOCalibrator = new PDOToCharge();
  if(b_tdo) TDOCalibrator = new TDOToTime(TDOFileName);
  else      TDOCalibrator = new TDOToTime();

  if(b_align)
    GEOMETRY->SetAlignment(AlignFileName);

  TFile* f = new TFile(inputFileName, "READ");
  if(!f){
    std::cout << "Error: unable to open input file " << inputFileName << std::endl;
    return 0;
  }
  TTree* T = (TTree*) f->Get("COMB_data");
  if(!T){
    std::cout << "Error: cannot find tree COMB_data in " << inputFileName << std::endl;
    return 0;
  }

  DATA = (MMDataAnalysis*) new MMDataAnalysis(T);
  int Nevent = DATA->GetNEntries();

  TH2D* scint_bot_vs_top;

  TH2D* track_hits_v_board;
  TH2D* track_param_x;
  TH2D* track_param_y;
  TH2D* track_angle_6;
  TH2D* track_angle_7;
  TH2D* track_angle_8;
  TH2D* track_angle_N_denom;
  TH2D* track_angle_N_numer;
  TH2D* track_hits_vs_evt;
  TH2D* track_hits_vs_evt_fid;
  TH2D* track_clusmult_vs_theta;

  TH2D* track_mx_vs_cx_0123_4567;
  TH2D* track_my_vs_cy_0123_4567;

  TH2D* track_mx_vs_cx_0236_1457;
  TH2D* track_my_vs_cy_0236_1457;

  TH2D* track_mx_vs_cx_0256_1347;
  TH2D* track_my_vs_cy_0256_1347;

  TH2D* track_slope_x_scint_19;
  TH2D* track_slope_x_scint_20;
  TH2D* track_slope_x_scint_21;

  TH2D* track_angle_x_scint_19;
  TH2D* track_angle_x_scint_20;
  TH2D* track_angle_x_scint_21;

  TH2D* track_scint_vs_time;

  TH2D* track_x_vs_y_0123_4567;
  TH2D* track_x_vs_y_0236_1457;
  TH2D* track_x_vs_y_0256_1347;

  TH2D* track_N1_board_vs_residual;
  TH2D* track_N1_board_vs_residual_10deg;
  std::vector<TH2D*> track_N1_theta_x_vs_residual;
  std::vector<TH2D*> track_N1_theta_y_vs_residual;
  TH2D* track_N1_board_vs_residual_art;

  TH2D* strip_position_vs_board;
  TH2D* strip_timediff_vs_pdo;
  std::vector<TH2D*> strip_charge_vs_channel;
  std::vector<TH2D*> strip_pdo_vs_channel;
  std::vector<TH2D*> strip_tdo_vs_channel;
  std::vector<TH2D*> strip_tdoc_vs_channel;
  std::vector<TH2D*> strip_dbc_vs_channel;
  std::vector<TH2D*> strip_time_vs_channel;
  std::vector<TH2D*> strip_zpos_vs_channel;
  std::vector<TH2D*> strip_event_vs_channel;
  TH2D* clus_vs_board;
  TH2D* hits_vs_board;
  TH2D* dups_vs_board;
  TH2D* dups_vs_channel;
  TH2D* hits_per_clus_vs_board;
  TH2D* timediff_vs_board;
  TH1D* timediff_track;
  TH2D* timediff_vs_charge;

  TH1D* track_vs_time;
  TH2D* track_exp_hit;
  TH2D* track_obs_hit;

  TH2D* trig_per_event;
  TH1D* trig_ambig;
  TH2D* trig_dtheta_vsNX;
  TH2D* trig_dtheta_vsNX_ok;
  TH1D* trig_mm;
  TH1D* trig_art;
  TH1D* trig_theta;
  TH1D* trig_eventmm;
  TH2D* trig_dbc_vs_N;
  TH2D* trig_dbc_vs_theta;
  TH2D* trig_dbc_vs_evt;
  TH2D* trig_art_vs_evt;
  TH1D* trig_dbc_pairs;

  TH1D* trig_dtheta_all_NX;
  TH1D* trig_dtheta_all_2X;
  TH1D* trig_dtheta_all_3X;
  TH1D* trig_dtheta_all_4X;
  TH1D* trig_dtheta_nearby_NX;
  TH1D* trig_dtheta_nearby_2X;
  TH1D* trig_dtheta_nearby_3X;
  TH1D* trig_dtheta_nearby_4X;
  TH1D* trig_dtheta_nearby_15deg_NX;
  TH1D* trig_dtheta_nearby_15deg_2X;
  TH1D* trig_dtheta_nearby_15deg_3X;
  TH1D* trig_dtheta_nearby_15deg_4X;

  TH1D* trig_dtheta_near64_NX; TH1D* trig_dtheta_near64_3X; TH1D* trig_dtheta_near64_4X;
  TH1D* trig_dtheta_near48_NX; TH1D* trig_dtheta_near48_3X; TH1D* trig_dtheta_near48_4X;
  TH1D* trig_dtheta_near36_NX; TH1D* trig_dtheta_near36_3X; TH1D* trig_dtheta_near36_4X;
  TH1D* trig_dtheta_near24_NX; TH1D* trig_dtheta_near24_3X; TH1D* trig_dtheta_near24_4X;
  TH1D* trig_dtheta_near16_NX; TH1D* trig_dtheta_near16_3X; TH1D* trig_dtheta_near16_4X;
  TH1D* trig_dtheta_near12_NX; TH1D* trig_dtheta_near12_3X; TH1D* trig_dtheta_near12_4X;
  TH1D* trig_dtheta_near08_NX; TH1D* trig_dtheta_near08_3X; TH1D* trig_dtheta_near08_4X;
  TH1D* trig_dtheta_near04_NX; TH1D* trig_dtheta_near04_3X; TH1D* trig_dtheta_near04_4X;

  TH1D* trig_dtheta_nearby_15deg_3X_not0;
  TH1D* trig_dtheta_nearby_15deg_3X_not1;
  TH1D* trig_dtheta_nearby_15deg_3X_not6;
  TH1D* trig_dtheta_nearby_15deg_3X_not7;

  TH2D* trig_dtheta_vs_theta_all;
  TH2D* trig_dtheta_vs_theta_near24;
  TH2D* trig_dtheta_vs_theta_near12;

  TH2D* trig_dx_vs_theta_0;
  TH2D* trig_dx_vs_theta_1;
  TH2D* trig_dx_vs_theta_6;
  TH2D* trig_dx_vs_theta_7;

  TH1D* trig_dtheta_0;
  TH1D* trig_dtheta_1;
  TH1D* trig_dtheta_2;
  TH1D* trig_dtheta_3;
  TH1D* trig_dtheta_4;
  TH1D* trig_dtheta_5;

  TH1D* trig_dtheta_0_else;
  TH1D* trig_dtheta_1_else;
  TH1D* trig_dtheta_2_else;
  TH1D* trig_dtheta_3_else;
  TH1D* trig_dtheta_4_else;
  TH1D* trig_dtheta_5_else;

  TH1D* trig_dtheta_2X;
  TH1D* trig_dtheta_3X;
  TH1D* trig_dtheta_4X;

  TH2D* trig_artpos_hit_vs_clussize;
  TH2D* trig_artpos_hit_vs_clussize_board0;
  TH2D* trig_artpos_hit_vs_clussize_board1;
  TH2D* trig_artpos_hit_vs_clussize_board2;
  TH2D* trig_artpos_hit_vs_clussize_board3;
  TH2D* trig_artpos_hit_vs_clussize_board4;
  TH2D* trig_artpos_hit_vs_clussize_board5;
  TH2D* trig_artpos_hit_vs_clussize_board6;
  TH2D* trig_artpos_hit_vs_clussize_board7;
  TH2D* trig_artpos_bci_vs_clussize;
  TH2D* trig_artpos_tdo_vs_clussize;
  TH2D* trig_artpos_pdo_vs_clussize;

  TH1D* trig_artwin_art;
  TH1D* trig_artwin_mmfe_mimic;
  TH1D* trig_artwin_mmfe_first;

  TH1D* mmfe_x;
  TH1D* mmfe_y;
  TH1D* mmfe_y_1u1v;
  TH2D* mmfe_xy_1u1v;
  TH1D* trig_x;
  TH1D* trig_y;
  TH1D* trig_vs_mmfe_x;
  TH1D* trig_vs_mmfe_y;
  TH1D* trig_vs_mmfe_y_nominal;
  TH1D* trig_vs_mmfe_y_bestxxx;
  TH1D* trig_vs_mmfe_y_bestpair;
  TH2D* trig_vs_mmfe_x_vs_theta;
  TH2D* trig_vs_mmfe_y_vs_theta;

  TH2D* trig_dtheta_vs_refmaxresidual;
  TH2D* trig_dtheta_vs_refdtheta;
  TH2D* trig_dtheta_vs_maxclustersize;
  TH2D* trig_dtheta_vs_dx_in_clus;

  TH1D* trig_nx_missing;
  TH1D* trig_dx_in_clus;
  TH1D* trig_refsumx2_all;
  TH1D* trig_refsumx2_bad;
  TH1D* trig_refmaxresidual;
  TH1D* trig_refdtheta;
  TH1D* trig_dtheta_unmatched_0;
  TH1D* trig_dtheta_unmatched_1;
  TH1D* trig_dtheta_unmatched_2;
  
  TH1D* trig_art_bcid;
  TH1D* mmfe_hit_bcid;
  TH1D* mmfe_fir_bcid;

  TH1D* trig_dbc_scint;

//   const std::vector<double> zboard = {  1.260,
//                                        10.520,
//                                        33.710,
//                                        42.740,
//                                       114.240,
//                                       124.260,
//                                       147.380,
//                                       152.900};
  const std::vector<double> zboard = {  0.0,
                                       11.2,
                                       32.4,
                                       43.6,
                                      113.6,
                                      124.8,
                                      146.0,
                                      157.2};

  scint_bot_vs_top = new TH2D("scint_bot_vs_top", ";SC bottom channel;SC top channel;Tracks", 8, -1.5, 6.5, 8, 14.5, 22.5);

  track_angle_x_scint_19 = new TH2D("track_angle_x_scint_19", ";x angle;scint. bottom channel;Tracks, scint. top channel 19", 100, -30, 30, 8, -1.5, 6.5);
  track_angle_x_scint_20 = new TH2D("track_angle_x_scint_20", ";x angle;scint. bottom channel;Tracks, scint. top channel 20", 100, -30, 30, 8, -1.5, 6.5);
  track_angle_x_scint_21 = new TH2D("track_angle_x_scint_21", ";x angle;scint. bottom channel;Tracks, scint. top channel 21", 100, -30, 30, 8, -1.5, 6.5);

  track_slope_x_scint_19 = new TH2D("track_slope_x_scint_19", ";x slope;scint. bottom channel;Tracks, scint. top channel 19", 100, -0.6, 0.6, 8, -1.5, 6.5);
  track_slope_x_scint_20 = new TH2D("track_slope_x_scint_20", ";x slope;scint. bottom channel;Tracks, scint. top channel 20", 100, -0.6, 0.6, 8, -1.5, 6.5);
  track_slope_x_scint_21 = new TH2D("track_slope_x_scint_21", ";x slope;scint. bottom channel;Tracks, scint. top channel 21", 100, -0.6, 0.6, 8, -1.5, 6.5);

  track_scint_vs_time = new TH2D("track_scint_vs_time", ";MM Event Number/1000;track angle, bottom scint. 0;Tracks", 1000, 0, 1000, 100, -30, 30);

  track_hits_v_board    = new TH2D("track_hits_v_board",    ";Board;N(Clusters);Tracks", 8, -0.5, 7.5, 6, -0.5, 5.5);
  track_hits_vs_evt     = new TH2D("track_hits_vs_evt",     ";MM Event Number;N(clusters in track);", 1000, 0, 1000, 10, -0.5, 9.5);
  track_hits_vs_evt_fid = new TH2D("track_hits_vs_evt_fid", ";MM Event Number;N(clusters in track);", 1000, 0, 1000, 10, -0.5, 9.5);

  track_param_x = new TH2D("track_param_x", ";x slope; x constant;Tracks",   100, -0.6, 0.6, 100,  -80, 280);
  track_param_y = new TH2D("track_param_y", ";y slope; y constant;Tracks",   100,   -4,   4, 100, -300, 500);
  track_angle_6 = new TH2D("track_angle_6", ";#theta#lower[0.5]{xz} [degrees];#theta#lower[0.5]{yz} [degrees];Tracks", 100, -35, 35, 100, -220, 220);
  track_angle_7 = new TH2D("track_angle_7", ";#theta#lower[0.5]{xz} [degrees];#theta#lower[0.5]{yz} [degrees];Tracks", 100, -35, 35, 100, -220, 220);
  track_angle_8 = new TH2D("track_angle_8", ";#theta#lower[0.5]{xz} [degrees];#theta#lower[0.5]{yz} [degrees];Tracks", 100, -35, 35, 100, -220, 220);
  track_angle_N_denom = new TH2D("track_angle_N_denom", ";#theta#lower[0.5]{xz} [degrees];#theta#lower[0.5]{yz} [degrees];Tracks", 700, -35, 35, 20, -220, 220);
  track_angle_N_numer = new TH2D("track_angle_N_numer", ";#theta#lower[0.5]{xz} [degrees];#theta#lower[0.5]{yz} [degrees];Tracks", 700, -35, 35, 20, -220, 220);

  track_clusmult_vs_theta = new TH2D("track_clusmult_vs_theta", ";#theta#lower[0.5]{xz} [degrees];hits in cluster", 100, -35, 35, 13, -0.5, 12.5);

  track_mx_vs_cx_0123_4567 = new TH2D("track_mx_vs_cx_0123_vs_4567", ";#Delta(x slope);#Delta(x constant);Events with 8 boards", 100, -0.6, 0.6, 400,   -40.0,   40.0);
  track_my_vs_cy_0123_4567 = new TH2D("track_my_vs_cy_0123_vs_4567", ";#Delta(y slope);#Delta(y constant);Events with 8 boards", 100,  -40,  40, 400, -2000.0, 2000.0);

  track_mx_vs_cx_0256_1347 = new TH2D("track_mx_vs_cx_0256_vs_1347", ";#Delta(x slope);#Delta(x constant);Events with 8 boards", 100, -0.08, 0.08, 400,    -8.0,    8.0);
  track_my_vs_cy_0256_1347 = new TH2D("track_my_vs_cy_0256_vs_1347", ";#Delta(y slope);#Delta(y constant);Events with 8 boards", 100,   -10,   10, 400,  -800.0,  800.0);

  track_mx_vs_cx_0236_1457 = new TH2D("track_mx_vs_cx_0236_vs_1457", ";#Delta(x slope);#Delta(x constant);Events with 8 boards", 100, -0.08, 0.08, 400,    -8.0,    8.0);
  track_my_vs_cy_0236_1457 = new TH2D("track_my_vs_cy_0236_vs_1457", ";#Delta(y slope);#Delta(y constant);Events with 8 boards", 100,   -20,   20, 400, -2000.0, 2000.0);

  // track_x_vs_y_0123_4567 = new TH2D("track_x_vs_y_0123_vs_4567", ";#Delta(x @ middle);#Delta(y @ middle);Events with 8 boards", 100, -40, 40, 100, -2000, 2000);
  // track_x_vs_y_0256_1347 = new TH2D("track_x_vs_y_0256_vs_1347", ";#Delta(x @ middle);#Delta(y @ middle);Events with 8 boards", 100,  -5,  5, 100,  -200,  200);
  // track_x_vs_y_0236_1457 = new TH2D("track_x_vs_y_0236_vs_1457", ";#Delta(x @ middle);#Delta(y @ middle);Events with 8 boards", 100,  -5,  5, 100,  -800,  800);

  track_x_vs_y_0123_4567 = new TH2D("track_x_vs_y_0123_vs_4567", ";#Delta(x @ middle);#Delta(y @ middle);Events with 8 boards", 400, -40, 40, 400, -2000, 2000);
  track_x_vs_y_0256_1347 = new TH2D("track_x_vs_y_0256_vs_1347", ";#Delta(x @ middle);#Delta(y @ middle);Events with 8 boards", 400, -40, 40, 400, -2000, 2000);
  track_x_vs_y_0236_1457 = new TH2D("track_x_vs_y_0236_vs_1457", ";#Delta(x @ middle);#Delta(y @ middle);Events with 8 boards", 400, -40, 40, 400, -2000, 2000);

  track_N1_board_vs_residual       = new TH2D("track_N1_board_vs_residual",       ";board;x_{cluster} - x_{track, proj.};Tracks", 8, -0.5, 7.5, 200, -5.0, 5.0);
  track_N1_board_vs_residual_10deg = new TH2D("track_N1_board_vs_residual_10deg", ";board;x_{cluster} - x_{track, proj.};Tracks", 8, -0.5, 7.5, 200, -5.0, 5.0);
  track_N1_board_vs_residual_art   = new TH2D("track_N1_board_vs_residual_art",   ";board;x_{cluster} - x_{track, proj.};Tracks", 8, -0.5, 7.5, 200, -5.0, 5.0);
  for (ibo = 0; ibo < nboards; ibo++){
    track_N1_theta_x_vs_residual.push_back(new TH2D(Form("track_N1_theta_x_vs_residual_%i", ibo), ";x theta;x_{cluster} - x_{track, proj.};Tracks", 100, -35, 35, 200, -5.0, 5.0));
    track_N1_theta_y_vs_residual.push_back(new TH2D(Form("track_N1_theta_y_vs_residual_%i", ibo), ";y theta;x_{cluster} - x_{track, proj.};Tracks", 100, -35, 35, 200, -5.0, 5.0));
  }

  strip_position_vs_board = new TH2D("strip_position_vs_board", ";strip number;MMFE number;charge [fC]", 512, 0.5, 512.5, 8, -0.5, 7.5);

  strip_timediff_vs_pdo = new TH2D("strip_timediff_vs_pdo", ";#DeltaBC;PDO;Strips", 60, -0.5, 59.5, 512, 0, 1024);

  for (ibo = 0; ibo < nboards; ibo++){
    strip_charge_vs_channel.push_back(new TH2D(Form("strip_charge_vs_channel_%i", ibo), ";strip number;Charge [fC];strip",    512, 0.5, 512.5,  512,   0,  128));
    strip_pdo_vs_channel.push_back   (new TH2D(Form("strip_pdo_vs_channel_%i",    ibo), ";strip number;PDO [counts];strip",   512, 0.5, 512.5,  512,   0, 2048));
    strip_tdo_vs_channel.push_back   (new TH2D(Form("strip_tdo_vs_channel_%i",    ibo), ";strip number;TDO [counts];strip",   512, 0.5, 512.5,  256,   0,  256));
    strip_tdoc_vs_channel.push_back  (new TH2D(Form("strip_tdoc_vs_channel_%i",   ibo), ";strip number;TDO corr. [ns];strip", 512, 0.5, 512.5,  110, -10,  100));
    strip_dbc_vs_channel.push_back   (new TH2D(Form("strip_dbc_vs_channel_%i",    ibo), ";strip number;#Delta BC;strip",      512, 0.5, 512.5,   64,   0,   64));
    strip_time_vs_channel.push_back  (new TH2D(Form("strip_time_vs_channel_%i",   ibo), ";strip number;Time [ns];strip",      512, 0.5, 512.5,  100, 300, 1000));
    strip_zpos_vs_channel.push_back  (new TH2D(Form("strip_zpos_vs_channel_%i",   ibo), ";strip number;z_{drift} [mm];strip", 512, 0.5, 512.5,  100,  10,   50));
    strip_event_vs_channel.push_back (new TH2D(Form("strip_event_vs_channel_%i",  ibo), ";strip number;Event number;strip",   512, 0.5, 512.5, 1000,   0, 1000));
  }
  dups_vs_channel = new TH2D("dups_vs_channel", ";strip number;MMFE number;Duplicates", 512, 0.5, 512.5, 8, -0.5, 7.5);

  clus_vs_board          = new TH2D("clus_vs_board",          ";MMFE number;clusters;Events",            8, -0.5, 7.5, 32, -0.5, 31.5);
  hits_vs_board          = new TH2D("hits_vs_board",          ";MMFE number;strips;Events",              8, -0.5, 7.5, 32, -0.5, 31.5);
  dups_vs_board          = new TH2D("dups_vs_board",          ";MMFE number;duplicate strips;Events",    8, -0.5, 7.5, 32, -0.5, 31.5);
  hits_per_clus_vs_board = new TH2D("hits_per_clus_vs_board", ";MMFE number;hits in a cluster;Clusters", 8, -0.5, 7.5, 13, -0.5, 12.5);
  timediff_vs_board      = new TH2D("timediff_vs_board",      ";MMFE number;#DeltaBCID;Events",          8, -0.5, 7.5, 100, -0.5, 99.5);
  timediff_track         = new TH1D("timediff_track",         ";#DeltaBCID;Strips on-track",             100, -0.5, 99.5);
  timediff_vs_charge     = new TH2D("timediff_vs_charge",     ";charge [fC];#DeltaBCID;Events",          200, 0, 200, 100, -0.5, 99.5);

  track_exp_hit = new TH2D("track_exp_hit", ";VMM;Board;Trigger", 10, -1.5, 8.5, 17, -0.75, 7.75);
  track_obs_hit = new TH2D("track_obs_hit", ";VMM;Board;Trigger", 10, -1.5, 8.5, 17, -0.75, 7.75);
  track_vs_time = create_track_vs_time(DATA);

  trig_per_event      = new TH2D("trig_per_event",      ";N(triggers);N(boards in MM track)", 20, -0.5, 19.5, 7, 2.5, 9.5);
  trig_ambig          = new TH1D("trig_ambig",          ";N(spurious ART, first) - N(spurious ART, other)", 13, -6.5, 6.5);
  trig_dtheta_vsNX    = new TH2D("trig_dtheta_vsNX",    ";#theta(TP)-#theta(MM);N(ART, X);Events", 1000, -500, 500, 5, 0.5, 5.5);
  trig_dtheta_vsNX_ok = new TH2D("trig_dtheta_vsNX_ok", ";#theta(TP)-#theta(MM);N(ART, X);Events", 1000, -500, 500, 5, 0.5, 5.5);
  //trig_dtheta_vsNX    = new TH2D("trig_dtheta_vsNX",    ";#theta(TP)-#theta(MM);N(ART, X);Events", 600, -100, 100, 5, 0.5, 5.5);
  //trig_dtheta_vsNX_ok = new TH2D("trig_dtheta_vsNX_ok", ";#theta(TP)-#theta(MM);N(ART, X);Events", 600, -100, 100, 5, 0.5, 5.5);
  trig_mm             = new TH1D("trig_mm",             ";N(MM clusters);",                        10, -0.5, 9.5);
  trig_art            = new TH1D("trig_art",            ";N(ART);",                                10, -0.5, 9.5);
  trig_theta          = new TH1D("trig_theta",          ";#theta [deg.];",                         100, -35, 35);
  trig_eventmm        = new TH1D("trig_eventmm",        ";MM Event Number;",                       1000, 0, 1000);
  trig_dbc_vs_N       = new TH2D("trig_dbc_vs_N",       ";#DeltaBC;N(ART)",                        10, -0.5, 9.5, 7, 2.5, 9.5);
  trig_dbc_vs_theta   = new TH2D("trig_dbc_vs_theta",   ";#DeltaBC;#theta [deg]",                  10, -0.5, 9.5, 500, -30, 20);
  trig_dbc_vs_evt     = new TH2D("trig_dbc_vs_evt",     ";MM Event Number;#DeltaBC",               1000, -100, 1000, 10, -0.5, 9.5);
  trig_dbc_pairs      = new TH1D("trig_dbc_pairs",      ";#DeltaBCID, pairs;",                     23, -11.5, 11.5);
  trig_art_vs_evt     = new TH2D("trig_art_vs_evt",     ";MM Event Number;N(ART in trigger)",      1000, -100, 1000, 10, -0.5, 9.5);

  trig_dtheta_all_NX          = new TH1D("trig_dtheta_all_NX",          ";#theta(TP)-#theta(MM);", 1000, -100, 100);
  trig_dtheta_all_2X          = new TH1D("trig_dtheta_all_2X",          ";#theta(TP)-#theta(MM);", 1000, -100, 100);
  trig_dtheta_all_3X          = new TH1D("trig_dtheta_all_3X",          ";#theta(TP)-#theta(MM);", 1000, -100, 100);
  trig_dtheta_all_4X          = new TH1D("trig_dtheta_all_4X",          ";#theta(TP)-#theta(MM);", 1000, -100, 100);
  trig_dtheta_nearby_NX       = new TH1D("trig_dtheta_nearby_NX",       ";#theta(TP)-#theta(MM);", 1000, -100, 100);
  trig_dtheta_nearby_2X       = new TH1D("trig_dtheta_nearby_2X",       ";#theta(TP)-#theta(MM);", 1000, -100, 100);
  trig_dtheta_nearby_3X       = new TH1D("trig_dtheta_nearby_3X",       ";#theta(TP)-#theta(MM);", 1000, -100, 100);
  trig_dtheta_nearby_4X       = new TH1D("trig_dtheta_nearby_4X",       ";#theta(TP)-#theta(MM);", 1000, -100, 100);
  trig_dtheta_nearby_15deg_NX = new TH1D("trig_dtheta_nearby_15deg_NX", ";#theta(TP)-#theta(MM);", 1000, -100, 100);
  trig_dtheta_nearby_15deg_2X = new TH1D("trig_dtheta_nearby_15deg_2X", ";#theta(TP)-#theta(MM);", 1000, -100, 100);
  trig_dtheta_nearby_15deg_3X = new TH1D("trig_dtheta_nearby_15deg_3X", ";#theta(TP)-#theta(MM);", 1000, -100, 100);
  trig_dtheta_nearby_15deg_4X = new TH1D("trig_dtheta_nearby_15deg_4X", ";#theta(TP)-#theta(MM);", 1000, -100, 100);

  trig_dtheta_nearby_15deg_3X_not0 = new TH1D("trig_dtheta_nearby_15deg_3X_not0", ";#theta(TP)-#theta(MM);", 1000, -100, 100);
  trig_dtheta_nearby_15deg_3X_not1 = new TH1D("trig_dtheta_nearby_15deg_3X_not1", ";#theta(TP)-#theta(MM);", 1000, -100, 100);
  trig_dtheta_nearby_15deg_3X_not6 = new TH1D("trig_dtheta_nearby_15deg_3X_not6", ";#theta(TP)-#theta(MM);", 1000, -100, 100);
  trig_dtheta_nearby_15deg_3X_not7 = new TH1D("trig_dtheta_nearby_15deg_3X_not7", ";#theta(TP)-#theta(MM);", 1000, -100, 100);

  trig_dtheta_near64_NX = new TH1D("trig_dtheta_near64_NX", ";#theta(TP)-#theta(MM);", 1000, -100, 100);
  trig_dtheta_near48_NX = new TH1D("trig_dtheta_near48_NX", ";#theta(TP)-#theta(MM);", 1000, -100, 100);
  trig_dtheta_near36_NX = new TH1D("trig_dtheta_near36_NX", ";#theta(TP)-#theta(MM);", 1000, -100, 100);
  trig_dtheta_near24_NX = new TH1D("trig_dtheta_near24_NX", ";#theta(TP)-#theta(MM);", 1000, -100, 100);
  trig_dtheta_near16_NX = new TH1D("trig_dtheta_near16_NX", ";#theta(TP)-#theta(MM);", 1000, -100, 100);
  trig_dtheta_near12_NX = new TH1D("trig_dtheta_near12_NX", ";#theta(TP)-#theta(MM);", 1000, -100, 100);
  trig_dtheta_near08_NX = new TH1D("trig_dtheta_near08_NX", ";#theta(TP)-#theta(MM);", 1000, -100, 100);
  trig_dtheta_near04_NX = new TH1D("trig_dtheta_near04_NX", ";#theta(TP)-#theta(MM);", 1000, -100, 100);

  trig_dtheta_near64_3X = new TH1D("trig_dtheta_near64_3X", ";#theta(TP)-#theta(MM);", 1000, -100, 100);
  trig_dtheta_near48_3X = new TH1D("trig_dtheta_near48_3X", ";#theta(TP)-#theta(MM);", 1000, -100, 100);
  trig_dtheta_near36_3X = new TH1D("trig_dtheta_near36_3X", ";#theta(TP)-#theta(MM);", 1000, -100, 100);
  trig_dtheta_near24_3X = new TH1D("trig_dtheta_near24_3X", ";#theta(TP)-#theta(MM);", 1000, -100, 100);
  trig_dtheta_near16_3X = new TH1D("trig_dtheta_near16_3X", ";#theta(TP)-#theta(MM);", 1000, -100, 100);
  trig_dtheta_near12_3X = new TH1D("trig_dtheta_near12_3X", ";#theta(TP)-#theta(MM);", 1000, -100, 100);
  trig_dtheta_near08_3X = new TH1D("trig_dtheta_near08_3X", ";#theta(TP)-#theta(MM);", 1000, -100, 100);
  trig_dtheta_near04_3X = new TH1D("trig_dtheta_near04_3X", ";#theta(TP)-#theta(MM);", 1000, -100, 100);

  trig_dtheta_near64_4X = new TH1D("trig_dtheta_near64_4X", ";#theta(TP)-#theta(MM);", 1000, -100, 100);
  trig_dtheta_near48_4X = new TH1D("trig_dtheta_near48_4X", ";#theta(TP)-#theta(MM);", 1000, -100, 100);
  trig_dtheta_near36_4X = new TH1D("trig_dtheta_near36_4X", ";#theta(TP)-#theta(MM);", 1000, -100, 100);
  trig_dtheta_near24_4X = new TH1D("trig_dtheta_near24_4X", ";#theta(TP)-#theta(MM);", 1000, -100, 100);
  trig_dtheta_near16_4X = new TH1D("trig_dtheta_near16_4X", ";#theta(TP)-#theta(MM);", 1000, -100, 100);
  trig_dtheta_near12_4X = new TH1D("trig_dtheta_near12_4X", ";#theta(TP)-#theta(MM);", 1000, -100, 100);
  trig_dtheta_near08_4X = new TH1D("trig_dtheta_near08_4X", ";#theta(TP)-#theta(MM);", 1000, -100, 100);
  trig_dtheta_near04_4X = new TH1D("trig_dtheta_near04_4X", ";#theta(TP)-#theta(MM);", 1000, -100, 100);

  trig_refsumx2_all = new TH1D("trig_refsumx2_all", "sum of squared residuals of ref.", 1000, 0, 2);
  trig_refsumx2_bad = new TH1D("trig_refsumx2_bad", "sum of squared residuals of ref.", 1000, 0, 2);

  trig_dtheta_vs_refmaxresidual = new TH2D("trig_dtheta_vs_refmaxresidual", "", 500, -100, 100, 500, 0, 8);
  trig_dtheta_vs_refdtheta      = new TH2D("trig_dtheta_vs_refdtheta",      "", 500, -100, 100, 500, 0, 8);
  trig_dtheta_vs_maxclustersize = new TH2D("trig_dtheta_vs_maxclustersize", "", 500, -100, 100,  30, -0.5, 29.5);
  trig_dtheta_vs_dx_in_clus     = new TH2D("trig_dtheta_vs_dx_in_clus",     "", 500, -100, 100,  30, -0.5, 29.5);

  trig_refmaxresidual = new TH1D("trig_refmaxresidual", "max residuals of ref.", 1000, 0, 8);
  trig_refdtheta      = new TH1D("trig_refdtheta",      "theta full #minus theta X", 1000, -40, 40);

  trig_dtheta_0 = new TH1D("trig_dtheta_0", "thetaTP #minus thetaMM", 1000, -80, 80);
  trig_dtheta_1 = new TH1D("trig_dtheta_1", "thetaTP #minus thetaMM", 1000, -80, 80);
  trig_dtheta_2 = new TH1D("trig_dtheta_2", "thetaTP #minus thetaMM", 1000, -80, 80);
  trig_dtheta_3 = new TH1D("trig_dtheta_3", "thetaTP #minus thetaMM", 1000, -80, 80);
  trig_dtheta_4 = new TH1D("trig_dtheta_4", "thetaTP #minus thetaMM", 1000, -80, 80);
  trig_dtheta_5 = new TH1D("trig_dtheta_5", "thetaTP #minus thetaMM", 1000, -80, 80);

  trig_dtheta_2X = new TH1D("trig_dtheta_2X", "thetaTP #minus thetaMM", 1000, -80, 80);
  trig_dtheta_3X = new TH1D("trig_dtheta_3X", "thetaTP #minus thetaMM", 1000, -80, 80);
  trig_dtheta_4X = new TH1D("trig_dtheta_4X", "thetaTP #minus thetaMM", 1000, -80, 80);

  trig_dtheta_0_else = new TH1D("trig_dtheta_0_else", "thetaTP #minus thetaMM", 1000, -80, 80);
  trig_dtheta_1_else = new TH1D("trig_dtheta_1_else", "thetaTP #minus thetaMM", 1000, -80, 80);
  trig_dtheta_2_else = new TH1D("trig_dtheta_2_else", "thetaTP #minus thetaMM", 1000, -80, 80);
  trig_dtheta_3_else = new TH1D("trig_dtheta_3_else", "thetaTP #minus thetaMM", 1000, -80, 80);
  trig_dtheta_4_else = new TH1D("trig_dtheta_4_else", "thetaTP #minus thetaMM", 1000, -80, 80);
  trig_dtheta_5_else = new TH1D("trig_dtheta_5_else", "thetaTP #minus thetaMM", 1000, -80, 80);

  trig_nx_missing = new TH1D("trig_nx_missing", "",  6, -2.5,  3.5);
  trig_dx_in_clus = new TH1D("trig_dx_in_clus", "", 30, -0.5, 29.5);

  trig_dtheta_unmatched_0 = new TH1D("trig_dtheta_unmatched_0", "#Delta#theta, 0 unmatched", 1000, -100, 100);
  trig_dtheta_unmatched_1 = new TH1D("trig_dtheta_unmatched_1", "#Delta#theta, 1 unmatched", 1000, -100, 100);
  trig_dtheta_unmatched_2 = new TH1D("trig_dtheta_unmatched_2", "#Delta#theta, 2 unmatched", 1000, -100, 100);

  trig_dtheta_vs_theta_all    = new TH2D("trig_dtheta_vs_theta_all",    ";#theta(MM);#theta(TP)-#theta(MM)", 100, -35, 35, 1000, -100, 100);
  trig_dtheta_vs_theta_near24 = new TH2D("trig_dtheta_vs_theta_near24", ";#theta(MM);#theta(TP)-#theta(MM)", 100, -35, 35, 1000, -100, 100);
  trig_dtheta_vs_theta_near12 = new TH2D("trig_dtheta_vs_theta_near12", ";#theta(MM);#theta(TP)-#theta(MM)", 100, -35, 35, 1000, -100, 100);

  trig_dx_vs_theta_0 = new TH2D("trig_dx_vs_theta_0", ";#theta(MM);#Deltax(ART, MMFE track) [mm]", 100, -35, 35, 1000, -16, 16);
  trig_dx_vs_theta_1 = new TH2D("trig_dx_vs_theta_1", ";#theta(MM);#Deltax(ART, MMFE track) [mm]", 100, -35, 35, 1000, -16, 16);
  trig_dx_vs_theta_6 = new TH2D("trig_dx_vs_theta_6", ";#theta(MM);#Deltax(ART, MMFE track) [mm]", 100, -35, 35, 1000, -16, 16);
  trig_dx_vs_theta_7 = new TH2D("trig_dx_vs_theta_7", ";#theta(MM);#Deltax(ART, MMFE track) [mm]", 100, -35, 35, 1000, -16, 16);

  trig_artpos_bci_vs_clussize = new TH2D("trig_artpos_bci_vs_clussize", ";ART index in cluster, by BCID;Cluster size",     8, -0.5, 7.5, 8, 0.5, 8.5);
  trig_artpos_tdo_vs_clussize = new TH2D("trig_artpos_tdo_vs_clussize", ";ART index in cluster, by BCID-TDO;Cluster size", 8, -0.5, 7.5, 8, 0.5, 8.5);
  trig_artpos_pdo_vs_clussize = new TH2D("trig_artpos_pdo_vs_clussize", ";ART index in cluster, by PDO;Cluster size",      8, -0.5, 7.5, 8, 0.5, 8.5);
  trig_artwin_art             = new TH1D("trig_artwin_art",        ";BC Window;Events", 13, -0.5, 12.5);
  trig_artwin_mmfe_mimic      = new TH1D("trig_artwin_mmfe_mimic", ";BC Window;Events", 13, -0.5, 12.5);
  trig_artwin_mmfe_first      = new TH1D("trig_artwin_mmfe_first", ";BC Window;Events", 13, -0.5, 12.5);

  trig_artpos_hit_vs_clussize         = new TH2D("trig_artpos_hit_vs_clussize",        ";ART index in cluster, by channel;Cluster size",  8, -0.5, 7.5, 8, 0.5, 8.5);
  trig_artpos_hit_vs_clussize_board0  = new TH2D("trig_artpos_hit_vs_clussize_board0", ";ART index in cluster, by channel;Cluster size",  8, -0.5, 7.5, 8, 0.5, 8.5);
  trig_artpos_hit_vs_clussize_board1  = new TH2D("trig_artpos_hit_vs_clussize_board1", ";ART index in cluster, by channel;Cluster size",  8, -0.5, 7.5, 8, 0.5, 8.5);
  trig_artpos_hit_vs_clussize_board2  = new TH2D("trig_artpos_hit_vs_clussize_board2", ";ART index in cluster, by channel;Cluster size",  8, -0.5, 7.5, 8, 0.5, 8.5);
  trig_artpos_hit_vs_clussize_board3  = new TH2D("trig_artpos_hit_vs_clussize_board3", ";ART index in cluster, by channel;Cluster size",  8, -0.5, 7.5, 8, 0.5, 8.5);
  trig_artpos_hit_vs_clussize_board4  = new TH2D("trig_artpos_hit_vs_clussize_board4", ";ART index in cluster, by channel;Cluster size",  8, -0.5, 7.5, 8, 0.5, 8.5);
  trig_artpos_hit_vs_clussize_board5  = new TH2D("trig_artpos_hit_vs_clussize_board5", ";ART index in cluster, by channel;Cluster size",  8, -0.5, 7.5, 8, 0.5, 8.5);
  trig_artpos_hit_vs_clussize_board6  = new TH2D("trig_artpos_hit_vs_clussize_board6", ";ART index in cluster, by channel;Cluster size",  8, -0.5, 7.5, 8, 0.5, 8.5);
  trig_artpos_hit_vs_clussize_board7  = new TH2D("trig_artpos_hit_vs_clussize_board7", ";ART index in cluster, by channel;Cluster size",  8, -0.5, 7.5, 8, 0.5, 8.5);

  trig_art_bcid = new TH1D("trig_art_bcid", ";BCID of ARTs;",                      11, -0.5, 10.5);
  mmfe_hit_bcid = new TH1D("mmfe_hit_bcid", ";BCID of MMFE strips;",               11, -0.5, 10.5);
  mmfe_fir_bcid = new TH1D("mmfe_fir_bcid", ";BCID of earliest strip in cluster;", 11, -0.5, 10.5);

  trig_dbc_scint = new TH1D("trig_dbc_scint", ";#DeltaBC(trig, scint);", 8191, -4095.5, 4095.5);

  trig_x = new TH1D("trig_x", ";x [mm];Events", 1000, -100, 250);
  trig_y = new TH1D("trig_y", ";y [mm];Events", 1000, -150, 300);
  mmfe_x = new TH1D("mmfe_x", ";x [mm];Events", 1000, -100, 250);
  mmfe_y = new TH1D("mmfe_y", ";y [mm];Events", 1000, -150, 300);
  mmfe_y_1u1v  = new TH1D("mmfe_y_1u1v",  ";y [mm];Events",       1000, -150, 300);
  mmfe_xy_1u1v = new TH2D("mmfe_xy_1u1v", ";x [mm];y [mm];Events", 200, -150, 300, 200, -150, 300);

  trig_vs_mmfe_x          = new TH1D("trig_vs_mmfe_x",          ";x(MMTP) - x(MMFE) [mm];Events",        1000, -8, 8);
  trig_vs_mmfe_y          = new TH1D("trig_vs_mmfe_y",          ";y(MMTP) - y(MMFE) [mm];Events",        1000, -300, 300);
  trig_vs_mmfe_x_vs_theta = new TH2D("trig_vs_mmfe_x_vs_theta", ";#theta;x(MMTP) - x(MMFE) [mm];Events", 400, -30, 20, 400, -8, 8);
  trig_vs_mmfe_y_vs_theta = new TH2D("trig_vs_mmfe_y_vs_theta", ";#theta;y(MMTP) - y(MMFE) [mm];Events", 400, -30, 20, 400, -300, 300);

  trig_vs_mmfe_y_nominal  = new TH1D("trig_vs_mmfe_y_nominal",  ";y(MMTP) - y(MMFE) [mm];Events", 1000, -300, 300);
  trig_vs_mmfe_y_bestxxx  = new TH1D("trig_vs_mmfe_y_bestxxx",  ";y(MMTP) - y(MMFE) [mm];Events", 1000, -300, 300);
  trig_vs_mmfe_y_bestpair = new TH1D("trig_vs_mmfe_y_bestpair", ";y(MMTP) - y(MMFE) [mm];Events", 1000, -300, 300);

  // cataloging duplicates
  MMFE8Hits duplicates;

  // collecting clusters and the nominal fit
  std::vector<MMClusterList> clusters_perboard;
  MMClusterList clusters_all;
  MMClusterList clusters_road;
  MMClusterList clusters_tp;
  MMClusterList clusters_x;
  MMTrack track;
  MMTrack track_tp;

  // N-1 fits and resolutions
  MMClusterList clusters_N1;
  MMTrack track_N1;

  // comparing two 4-hit tracks: different quads
  MMTrack track_0123;
  MMTrack track_4567;
  MMClusterList clus_0123;
  MMClusterList clus_4567;

  // comparing two 4-hit tracks: diff-quad Xs, diff-quad UV
  MMTrack track_0256;
  MMTrack track_1347;
  MMClusterList clus_0256;
  MMClusterList clus_1347;

  // comparing two 4-hit tracks: diff-quad Xs, same-quad UV
  MMTrack track_0236;
  MMTrack track_1457;
  MMClusterList clus_0236;
  MMClusterList clus_1457;

  // miscellaneous helpers
  std::vector<double> residuals = {};
  double residual = 0.0;
  TCanvas* can;
  GeoPlane plane;
  double z_middle = 0.0;
  for (ibo = 0; ibo < nboards; ibo++)
    z_middle += GEOMETRY->Get(ibo).Origin().Z();
  z_middle /= (double)(nboards);
  double track_time_start = 0.0;

  // output
  TFile* fout = new TFile(outputFileName, "RECREATE");
  fout->mkdir("event_displays");
  MMPlot();

  // progress bar
  std::chrono::time_point<std::chrono::system_clock> time_start;
  std::chrono::duration<double> elapsed_seconds;
  time_start = std::chrono::system_clock::now();

  // loop
  for(int evt = 0; evt < Nevent; evt++){

    // get the event
    DATA->GetEntry(evt);
    if(evt % 1000 == 0){
      elapsed_seconds = (std::chrono::system_clock::now() - time_start);
      progress(elapsed_seconds.count(), evt, Nevent);
    }

    if (evt == 0)
      track_time_start = DATA->mm_EventHits.Time();
    if (GEOMETRY->RunNumber() < 0)
      GEOMETRY->SetRunNumber(DATA->RunNum);
    if (!DATA->sc_EventHits.IsGoodEvent())
      continue;

    do_trigger = (DATA->RunNum == 3522 || DATA->RunNum == 3539);
    EventNum4Hist = DATA->mm_EventNum;

    // for treating 3518-3520 as one big run
    if (DATA->RunNum == 3519 && Nevent > 400000)
      EventNum4Hist = EventNum4Hist + 200000;
    else if (DATA->RunNum == 3520 && Nevent > 400000)
      EventNum4Hist = EventNum4Hist + 500000;

    if (evt > 3980800)
      break;

    // reset clusters and tracks
    track.Reset();
    track_tp.Reset();
    track_N1.Reset();
    clusters_all.Reset();
    clusters_road.Reset();
    clusters_tp.Reset();
    clusters_N1.Reset();
    for (auto clus_list: clusters_perboard)
      clus_list.Reset();
    clusters_perboard.clear();
    
    // calibrate
    PDOCalibrator->Calibrate(DATA->mm_EventHits);
    TDOCalibrator->Calibrate(DATA->mm_EventHits);
    PACMAN->SetEventTrigBCID(DATA->mm_trig_BCID);
    PACMAN->SetEventPadTime(0);

    FILTERER->SetRunNumber(DATA->RunNum);    

    // run pacman
    nboardshit = DATA->mm_EventHits.GetNBoards();
    for(i = 0; i < nboardshit; i++){
      if(DATA->mm_EventHits[i].GetNHits() == 0)
        continue;
      MMClusterList board_clusters = PACMAN->Cluster(DATA->mm_EventHits[i]);
      if(board_clusters.GetNCluster() > 0)
        clusters_perboard.push_back(board_clusters);

      // strips quantities vs channel
      if (true){
        for(ich = 0; ich < DATA->mm_EventHits[i].GetNHits(); ich++){
          ibo = GEOMETRY->Index(DATA->mm_EventHits[i].MMFE8());
          timediff_vs_board->Fill(ibo, DATA->mm_trig_BCID - DATA->mm_EventHits[i][ich].BCID());
          strip_timediff_vs_pdo->Fill(DATA->mm_trig_BCID - DATA->mm_EventHits[i][ich].BCID(), DATA->mm_EventHits[i][ich].PDO());
          if( !PACMAN->IsGoodHit(DATA->mm_EventHits[i][ich]) )
            continue;
          timediff_vs_charge->Fill(DATA->mm_EventHits[i][ich].Charge(), DATA->mm_trig_BCID - DATA->mm_EventHits[i][ich].BCID());
          strip_charge_vs_channel[ibo]->Fill(DATA->mm_EventHits[i][ich].Channel(), DATA->mm_EventHits[i][ich].Charge());
          strip_pdo_vs_channel[ibo]   ->Fill(DATA->mm_EventHits[i][ich].Channel(), DATA->mm_EventHits[i][ich].PDO());
          strip_tdo_vs_channel[ibo]   ->Fill(DATA->mm_EventHits[i][ich].Channel(), DATA->mm_EventHits[i][ich].TDO());
          strip_dbc_vs_channel[ibo]   ->Fill(DATA->mm_EventHits[i][ich].Channel(), DATA->mm_trig_BCID - DATA->mm_EventHits[i][ich].BCID());
          strip_event_vs_channel[ibo] ->Fill(DATA->mm_EventHits[i][ich].Channel(), EventNum4Hist/1000.0);

          //strip_tdoc_vs_channel[ibo]  ->Fill(DATA->mm_EventHits[i][ich].Channel(), DATA->mm_EventHits[i][ich].Time());
          //strip_time_vs_channel[ibo]  ->Fill(DATA->mm_EventHits[i][ich].Channel(),  (DATA->mm_trig_BCID - DATA->mm_EventHits[i][ich].BCID())*25 + DATA->mm_EventHits[i][ich].Time());
          //strip_zpos_vs_channel[ibo]  ->Fill(DATA->mm_EventHits[i][ich].Channel(), ((DATA->mm_trig_BCID - DATA->mm_EventHits[i][ich].BCID())*25 + DATA->mm_EventHits[i][ich].Time()) * vdrift);
        }
      }

      // finding a pseudo-event display
      if (hella_clusters < 0 && board_clusters.size() > 2){

        // require a nice number of strips
        int has3strips = 0;
        for (auto clus: board_clusters)
          if (clus->size() >= 3 && clus->NHoles() < 2)
            has3strips++;

        // require some spatial separation
        double min_sep = 1000.0;
        for (auto clus1: board_clusters)
          for (auto clus2: board_clusters)
            if (clus1->Channel() > clus2->Channel())
              min_sep = std::min(min_sep, std::fabs(clus1->Channel() - clus2->Channel()));

        // announce
        if (has3strips > 2 and min_sep > 40){
          std::cout << std::endl;
          std::cout << "Event " << DATA->mm_EventNum              << std::endl;
          std::cout << "Board " << board_clusters[0].MMFE8Index() << std::endl;
          std::cout << "N(cl) " << board_clusters.size()          << std::endl;
          for (auto clus: board_clusters)
            std::cout << " N(hits) = " << clus->size() << " @ " << (int)(clus->Channel()) << std::endl;
          std::cout << std::endl;
          hella_clusters++;
        }
      }

      // writing pseudo-event display for clustering visualization
      if (DATA->RunNum == 3522 && DATA->mm_EventNum == 6897 && board_clusters.size() > 0)
        for(ich = 0; ich < DATA->mm_EventHits[i].GetNHits(); ich++)
          if(PACMAN->IsGoodHit(DATA->mm_EventHits[i][ich]))
            strip_position_vs_board->Fill(DATA->mm_EventHits[i][ich].Channel(), GEOMETRY->Index(board_clusters[0].MMFE8()), DATA->mm_EventHits[i][ich].Charge());

    }

    // hits, duplicates, clusters per board
    // ------------------------------------
    for (ibo = 0; ibo < nboards; ibo++){
      test = -1;
      for (i = 0; i < (int)(clusters_perboard.size()); i++){
        // hits on this board!
        if (GEOMETRY->Index(DATA->mm_EventHits[i].MMFE8()) == ibo){
          test = i;
          clus_vs_board->Fill(ibo, clusters_perboard[i].size());
          hits_vs_board->Fill(ibo, DATA->mm_EventHits[i].GetNHits());
          dups_vs_board->Fill(ibo, DATA->mm_EventHits[i].GetNDuplicates());
          duplicates = DATA->mm_EventHits[i].GetDuplicates();
          for(ich = 0; ich < (int)(duplicates.GetNHits()); ich++){
            // do you want to fill this with a weight of duplicates[ich].GetNHits()-1?
            dups_vs_channel->Fill(duplicates[ich].Channel(), ibo);
          }
        }
      }

      // no hits on this board!
      if (test == -1){
        clus_vs_board->Fill(ibo, 0);
        hits_vs_board->Fill(ibo, 0);
        dups_vs_board->Fill(ibo, 0);
      }
    }
    
    for (auto clus_list: clusters_perboard)
      for (auto clus: clus_list)
        hits_per_clus_vs_board->Fill(GEOMETRY->Index(clus->MMFE8()), clus->GetNHits());
    
    // preselection quality to run tracking
    // require at least 4 boards hit
    if (clusters_perboard.size() < 4)
      continue;

    // flatten these lists
    for (auto clus_list: clusters_perboard)
      for (auto clus: clus_list)
        clusters_all.AddCluster(*clus);
    
    // fit it!
    for (auto botpair: DATA->sc_EventHits.GetBotPair()){
      clusters_road.Reset();
      clusters_road = FILTERER->FilterClustersScint(clusters_all, *GEOMETRY, botpair.first->Channel(), DATA->mm_EventNum, 0);
      track = FITTER->Fit(clusters_road, *GEOMETRY, DATA->mm_EventNum);
    }

    // triggers!
    if (!track.IsFit())
      continue;
    if (!track.IsTrigCand())
      continue;

    // angles
    if (clusters_road.size() == 6) track_angle_6->Fill(atan(track.SlopeX())*180/3.14159, atan(track.SlopeY())*180/3.14159);
    if (clusters_road.size() == 7) track_angle_7->Fill(atan(track.SlopeX())*180/3.14159, atan(track.SlopeY())*180/3.14159);
    if (clusters_road.size() == 8) track_angle_8->Fill(atan(track.SlopeX())*180/3.14159, atan(track.SlopeY())*180/3.14159);

    // clusters per track
    track_hits_vs_evt->Fill(EventNum4Hist/1000.0, clusters_road.size());
    if (GEOMETRY->IsFiducial(track))
      track_hits_vs_evt_fid->Fill(EventNum4Hist/1000.0, clusters_road.size());

    // rate of good tracks
    if (clusters_road.size() >= 6 && GEOMETRY->IsFiducial(track))
      track_vs_time->Fill( (DATA->mm_EventHits.Time() - track_time_start)/60.0 );

    // timing of on-track hits
    if (clusters_road.size() >= 7)
      for (auto clus: clusters_road)
        for (i = 0; i < (int)(clus->GetNHits()); i++)
          timediff_track->Fill(DATA->mm_trig_BCID - clus->Get(i).BCID());

    // some sanity plots
    track_param_x->Fill(track.SlopeX(), track.ConstX());
    track_param_y->Fill(track.SlopeY(), track.ConstY());

    // strip multiplicity vs angle
    for (auto clus: clusters_road)
      track_clusmult_vs_theta->Fill(atan(track.SlopeX())*180/3.14159, clus->size());

    // interlude: some scintillator plots
    for (auto toppair: DATA->sc_EventHits.GetTopPair()){
      for (auto botpair: DATA->sc_EventHits.GetBotPair()){
        scint_bot_vs_top->Fill(botpair.first->Channel(), toppair.first->Channel());
        if (toppair.first->Channel() == 19){
          track_angle_x_scint_19->Fill(atan(track.SlopeX())*180/3.14159, botpair.first->Channel());
          track_slope_x_scint_19->Fill(     track.SlopeX(),              botpair.first->Channel());
        }
        else if (toppair.first->Channel() == 20){
          track_angle_x_scint_20->Fill(atan(track.SlopeX())*180/3.14159, botpair.first->Channel());
          track_slope_x_scint_20->Fill(     track.SlopeX(),              botpair.first->Channel());
        }
        else if (toppair.first->Channel() == 21){
          track_angle_x_scint_21->Fill(atan(track.SlopeX())*180/3.14159, botpair.first->Channel());
          track_slope_x_scint_21->Fill(     track.SlopeX(),              botpair.first->Channel());
        }

        if (botpair.first->Channel() == 0) 
          track_scint_vs_time->Fill(EventNum4Hist/1000.0, atan(track.SlopeX())*180/3.14159);
      }
    }

    if (debug){
      can = Plot_Track2D(Form("track2D_%05d_road", DATA->mm_EventNum), track, *GEOMETRY, &clusters_road);
      fout->cd("event_displays");
      // can->Write();
      delete can;
    }

    //
    // VMM-level efficiency
    //
    double xproj = 0.0;
    int vmm_min  = 0, vmm_max  = 0;
    int chan_min = 0, chan_max = 0;
    for (ibo = 0; ibo < 8; ibo++){

      // make track from other boards
      clusters_N1.Reset();
      track_N1.Reset();
      for (auto clus: clusters_road)
        if (clus->MMFE8Index() != ibo)
          clusters_N1.AddCluster(*clus);
      track_N1 = FITTER->Fit(clusters_N1, *GEOMETRY, DATA->mm_EventNum);
      if (!track_N1.IsFit())
        continue;
      if (track_N1.NX() + track_N1.NU() + track_N1.NV() < 6)
        continue;

      // find track position at board and corresponding channel
      // find it at both ends of the board, to be safe for UV
      xproj    = track_N1.SlopeX()*GEOMETRY->Get(ibo).Origin().Z() + track_N1.ConstX();
      chan_min = (int)(channel_from_x(xproj, ibo, GEOMETRY, true));
      chan_max = (int)(channel_from_x(xproj, ibo, GEOMETRY, false));

      if (chan_min <= 0 || chan_max <= 0)
        continue;

      vmm_min = (chan_min-1)/64;
      vmm_max = (chan_max-1)/64;
      if (vmm_min != vmm_max)
        continue;

      chan_min = (chan_min-1)%64 - 1;
      chan_max = (chan_max-1)%64 - 1;

      // avoid edge cases
      if (chan_min==1 || chan_min==2 || chan_min==63 || chan_min==64)
        continue;
      if (chan_max==1 || chan_max==2 || chan_max==63 || chan_max==64)
        continue;

      track_exp_hit->Fill(vmm_min, ibo);
      for (auto clus: clusters_road)
        if (clus->MMFE8Index() == ibo)
          track_obs_hit->Fill(vmm_min, ibo);
    }

    //
    // N-1 histograms (residuals)
    //
    if (clusters_road.size() == 8){
      for (auto clus: clusters_road){
      
        clusters_N1.Reset();
        track_N1.Reset();

        for (auto clus_other: clusters_road)
          if (clus != clus_other)
            clusters_N1.AddCluster(*clus_other);
      
        track_N1 = FITTER->Fit(clusters_N1, *GEOMETRY, DATA->mm_EventNum);
        residual = GEOMETRY->GetResidualX(*clus, track_N1);
        ibo      = GEOMETRY->Index(clus->MMFE8());

        track_N1_board_vs_residual       ->Fill(ibo,            residual);
        track_N1_theta_x_vs_residual[ibo]->Fill(atan(track_N1.SlopeX())*180/3.14159, residual);
        track_N1_theta_y_vs_residual[ibo]->Fill(atan(track_N1.SlopeY())*180/3.14159, residual);
        if (std::fabs(atan(track_N1.SlopeX())*180/3.14159) < 10)
          track_N1_board_vs_residual_10deg->Fill(ibo, residual);
      }
    }

    //
    // 4-hit track vs. 4-hit track histograms
    //
    if (clusters_road.size() == 8){
      track_0123.Reset(); clus_0123.Reset();
      track_4567.Reset(); clus_4567.Reset();
      track_0256.Reset(); clus_0256.Reset();
      track_1347.Reset(); clus_1347.Reset();
      track_0236.Reset(); clus_0236.Reset();
      track_1457.Reset(); clus_1457.Reset();
    
      for (auto clus: clusters_road){
        ibo = GEOMETRY->Index(clus->MMFE8());
        if (ibo == 0 || ibo == 1 || ibo == 2 || ibo == 3) clus_0123.AddCluster(*clus);
        if (ibo == 4 || ibo == 5 || ibo == 6 || ibo == 7) clus_4567.AddCluster(*clus);
        if (ibo == 0 || ibo == 2 || ibo == 5 || ibo == 6) clus_0256.AddCluster(*clus);
        if (ibo == 1 || ibo == 3 || ibo == 4 || ibo == 7) clus_1347.AddCluster(*clus);
        if (ibo == 0 || ibo == 2 || ibo == 3 || ibo == 6) clus_0236.AddCluster(*clus);
        if (ibo == 1 || ibo == 4 || ibo == 5 || ibo == 7) clus_1457.AddCluster(*clus);
      }

      track_0123 = FITTER->Fit(clus_0123, *GEOMETRY, DATA->mm_EventNum);
      track_4567 = FITTER->Fit(clus_4567, *GEOMETRY, DATA->mm_EventNum);
      track_0256 = FITTER->Fit(clus_0256, *GEOMETRY, DATA->mm_EventNum);
      track_1347 = FITTER->Fit(clus_1347, *GEOMETRY, DATA->mm_EventNum);
      track_0236 = FITTER->Fit(clus_0236, *GEOMETRY, DATA->mm_EventNum);
      track_1457 = FITTER->Fit(clus_1457, *GEOMETRY, DATA->mm_EventNum);
    
      track_mx_vs_cx_0123_4567->Fill(track_0123.SlopeX() - track_4567.SlopeX(), track_0123.ConstX() - track_4567.ConstX());
      track_my_vs_cy_0123_4567->Fill(track_0123.SlopeY() - track_4567.SlopeY(), track_0123.ConstY() - track_4567.ConstY());

      track_mx_vs_cx_0256_1347->Fill(track_0256.SlopeX() - track_1347.SlopeX(), track_0256.ConstX() - track_1347.ConstX());
      track_my_vs_cy_0256_1347->Fill(track_0256.SlopeY() - track_1347.SlopeY(), track_0256.ConstY() - track_1347.ConstY());

      track_mx_vs_cx_0236_1457->Fill(track_0236.SlopeX() - track_1457.SlopeX(), track_0236.ConstX() - track_1457.ConstX());
      track_my_vs_cy_0236_1457->Fill(track_0236.SlopeY() - track_1457.SlopeY(), track_0236.ConstY() - track_1457.ConstY());    

      track_x_vs_y_0123_4567->Fill(track_0123.PointZ(z_middle).X() - track_4567.PointZ(z_middle).X(),
                                   track_0123.PointZ(z_middle).Y() - track_4567.PointZ(z_middle).Y());
      track_x_vs_y_0256_1347->Fill(track_0256.PointZ(z_middle).X() - track_1347.PointZ(z_middle).X(),
                                   track_0256.PointZ(z_middle).Y() - track_1347.PointZ(z_middle).Y());
      track_x_vs_y_0236_1457->Fill(track_0236.PointZ(z_middle).X() - track_1457.PointZ(z_middle).X(),
                                   track_0236.PointZ(z_middle).Y() - track_1457.PointZ(z_middle).Y());
    }

    //
    // TRIGGERS
    //
    if (!do_trigger)
      continue;

    // scintillator
    int dB = 0;
    if (DATA->RunNum == 3539)
      dB = 138;
    else
      dB = 1944;
    DATA->tp_EventTracks.SetSciBCID(DATA->sc_EventHits.TPBCID(), DATA->sc_EventHits.TPph());
    DATA->tp_EventTracks.SetSciOffset(dB);

    TPTrack* neo = DATA->tp_EventTracks.Highlander(clusters_road, true, 16);

    if (neo)
      trig_dbc_scint->Fill(DATA->tp_EventTracks.deltaBCID(neo->BCID()));

    if (!track.IsTrigCand())
      continue;

    // denominators of efficiencies go here
    track_angle_N_denom->Fill(atan(track.SlopeX())*180/3.14159, atan(track.SlopeY())*180/3.14159);

    if (!neo)
      continue;
    if (!neo->IsTrigCand())
      continue;

    // numerators
    track_angle_N_numer->Fill(atan(track.SlopeX())*180/3.14159, atan(track.SlopeY())*180/3.14159);

    bool badvmm = false;
    for (auto art: *neo)
      if (art->MMFE8Index()==0 && art->VMM()==3)
        badvmm = true;
    if (badvmm)
      continue;
    
    // triggers
    double thetaMM = atan(track.SlopeX());
    double thetaTP = atan(neo->MxLocal());
    double dtheta  = thetaTP - thetaMM;
    trig_per_event->Fill(DATA->tp_EventTracks.GetNTrack(), clusters_road.size());
    trig_dtheta_vsNX->Fill(dtheta*1000, neo->NX());
    if (neo->NMatch() == neo->size())
      trig_dtheta_vsNX_ok->Fill(dtheta*1000, neo->NX());

    // dtheta
    int nearby = 0; 
    int near64 = 0, near48 = 0, near36 = 0, near24 = 0;
    int near16 = 0, near12 = 0, near08 = 0, near04 = 0;
    int near06 = 0;
    int nearUV = 0;
    double xtrack = 0;

    // match in X and nearby
    nearby = 1;
    near64 = 1; near48 = 1; near36 = 1; near24 = 1; 
    near16 = 1; near12 = 1; near08 = 1; near04 = 1; 
    near06 = 1;
    for (auto art: *neo){
      if (!art->isX())
        continue;
      xtrack = track.SlopeX()*zboard[art->MMFE8Index()] + track.ConstX();
      nearby = nearby && (fabs(xtrack - xpos(art, GEOMETRY)) < 24*0.4);

      near64 = near64 && (fabs(xtrack - xpos(art, GEOMETRY)) < 64*0.4);
      near48 = near48 && (fabs(xtrack - xpos(art, GEOMETRY)) < 48*0.4);
      near36 = near36 && (fabs(xtrack - xpos(art, GEOMETRY)) < 36*0.4);
      near24 = near24 && (fabs(xtrack - xpos(art, GEOMETRY)) < 24*0.4);
      near16 = near16 && (fabs(xtrack - xpos(art, GEOMETRY)) < 16*0.4);
      near12 = near12 && (fabs(xtrack - xpos(art, GEOMETRY)) < 12*0.4);
      near08 = near08 && (fabs(xtrack - xpos(art, GEOMETRY)) <  8*0.4);
      near06 = near06 && (fabs(xtrack - xpos(art, GEOMETRY)) <  6*0.4);
      near04 = near04 && (fabs(xtrack - xpos(art, GEOMETRY)) <  4*0.4);
    }

    // match in UV
    nearUV = 1;
    double xposmid = 0;
    for (auto art: *neo){
      if (art->isX())
        continue;
      xposmid = (xpos(art, GEOMETRY) + xpos_end(art, GEOMETRY)) / 2.0;
      xtrack  = track.SlopeX()*zboard[art->MMFE8Index()] + track.ConstX();
      nearUV  = nearUV && (fabs(xtrack - xposmid) < 40*0.4);
    }

    if (true)   trig_dtheta_vs_theta_all   ->Fill(thetaMM*180/3.14159, dtheta*1000);
    if (near24) trig_dtheta_vs_theta_near24->Fill(thetaMM*180/3.14159, dtheta*1000);
    if (near12) trig_dtheta_vs_theta_near12->Fill(thetaMM*180/3.14159, dtheta*1000);

    if (near64 && neo->NX() >= 2) trig_dtheta_near64_NX->Fill(dtheta*1000);
    if (near48 && neo->NX() >= 2) trig_dtheta_near48_NX->Fill(dtheta*1000);
    if (near36 && neo->NX() >= 2) trig_dtheta_near36_NX->Fill(dtheta*1000);
    if (near24 && neo->NX() >= 2) trig_dtheta_near24_NX->Fill(dtheta*1000);
    if (near16 && neo->NX() >= 2) trig_dtheta_near16_NX->Fill(dtheta*1000);
    if (near12 && neo->NX() >= 2) trig_dtheta_near12_NX->Fill(dtheta*1000);
    if (near08 && neo->NX() >= 2) trig_dtheta_near08_NX->Fill(dtheta*1000);
    if (near04 && neo->NX() >= 2) trig_dtheta_near04_NX->Fill(dtheta*1000);

    if (near64 && neo->NX() == 3) trig_dtheta_near64_3X->Fill(dtheta*1000);
    if (near48 && neo->NX() == 3) trig_dtheta_near48_3X->Fill(dtheta*1000);
    if (near36 && neo->NX() == 3) trig_dtheta_near36_3X->Fill(dtheta*1000);
    if (near24 && neo->NX() == 3) trig_dtheta_near24_3X->Fill(dtheta*1000);
    if (near16 && neo->NX() == 3) trig_dtheta_near16_3X->Fill(dtheta*1000);
    if (near12 && neo->NX() == 3) trig_dtheta_near12_3X->Fill(dtheta*1000);
    if (near08 && neo->NX() == 3) trig_dtheta_near08_3X->Fill(dtheta*1000);
    if (near04 && neo->NX() == 3) trig_dtheta_near04_3X->Fill(dtheta*1000);

    if (near64 && neo->NX() == 4) trig_dtheta_near64_4X->Fill(dtheta*1000);
    if (near48 && neo->NX() == 4) trig_dtheta_near48_4X->Fill(dtheta*1000);
    if (near36 && neo->NX() == 4) trig_dtheta_near36_4X->Fill(dtheta*1000);
    if (near24 && neo->NX() == 4) trig_dtheta_near24_4X->Fill(dtheta*1000);
    if (near16 && neo->NX() == 4) trig_dtheta_near16_4X->Fill(dtheta*1000);
    if (near12 && neo->NX() == 4) trig_dtheta_near12_4X->Fill(dtheta*1000);
    if (near08 && neo->NX() == 4) trig_dtheta_near08_4X->Fill(dtheta*1000);
    if (near04 && neo->NX() == 4) trig_dtheta_near04_4X->Fill(dtheta*1000);

    if (neo->NX() >= 2) trig_dtheta_all_NX->Fill(dtheta*1000);
    if (neo->NX() == 2) trig_dtheta_all_2X->Fill(dtheta*1000);
    if (neo->NX() == 3) trig_dtheta_all_3X->Fill(dtheta*1000);
    if (neo->NX() == 4) trig_dtheta_all_4X->Fill(dtheta*1000);

    if (nearby) {
      if (neo->NX() >= 2) trig_dtheta_nearby_NX->Fill(dtheta*1000);
      if (neo->NX() == 2) trig_dtheta_nearby_2X->Fill(dtheta*1000);
      if (neo->NX() == 3) trig_dtheta_nearby_3X->Fill(dtheta*1000);
      if (neo->NX() == 4) trig_dtheta_nearby_4X->Fill(dtheta*1000);
    }
    
    if (nearby && fabs(thetaTP)*180/3.14159 > 15) {
      if (neo->NX() >= 2) trig_dtheta_nearby_15deg_NX->Fill(dtheta*1000);
      if (neo->NX() == 2) trig_dtheta_nearby_15deg_2X->Fill(dtheta*1000);
      if (neo->NX() == 3) trig_dtheta_nearby_15deg_3X->Fill(dtheta*1000);
      if (neo->NX() == 4) trig_dtheta_nearby_15deg_4X->Fill(dtheta*1000);

      std::vector<int> bos = {};
      for (auto art: *neo)
        if (art->isX())
          bos.push_back(art->MMFE8Index());

      if (neo->NX() == 3){
        if (std::find(bos.begin(), bos.end(), 0) == bos.end()) trig_dtheta_nearby_15deg_3X_not0->Fill(dtheta*1000);
        if (std::find(bos.begin(), bos.end(), 1) == bos.end()) trig_dtheta_nearby_15deg_3X_not1->Fill(dtheta*1000);
        if (std::find(bos.begin(), bos.end(), 6) == bos.end()) trig_dtheta_nearby_15deg_3X_not6->Fill(dtheta*1000);
        if (std::find(bos.begin(), bos.end(), 7) == bos.end()) trig_dtheta_nearby_15deg_3X_not7->Fill(dtheta*1000);
      }
    }

    for (auto art: *neo){
      if (!art->isX())
        continue;
      if (!near06)
        continue;
      xtrack = track.SlopeX()*zboard[art->MMFE8Index()] + track.ConstX();
      if (art->MMFE8Index() == 0) trig_dx_vs_theta_0->Fill(thetaMM*180/3.14159, xpos(art, GEOMETRY) - xtrack);
      if (art->MMFE8Index() == 1) trig_dx_vs_theta_1->Fill(thetaMM*180/3.14159, xpos(art, GEOMETRY) - xtrack);
      if (art->MMFE8Index() == 6) trig_dx_vs_theta_6->Fill(thetaMM*180/3.14159, xpos(art, GEOMETRY) - xtrack);
      if (art->MMFE8Index() == 7) trig_dx_vs_theta_7->Fill(thetaMM*180/3.14159, xpos(art, GEOMETRY) - xtrack);
    }
    
    // BCID window / collection time!
    int dbcid, dbcpair;
    dbcid = neo->BCIDWindow();

    if (dbcid >= 0){

      trig_dbc_vs_N  ->Fill(dbcid, neo->size());
      trig_dbc_vs_evt->Fill(EventNum4Hist/1000.0, dbcid);
      trig_dbc_vs_theta->Fill(dbcid, thetaTP*180/3.14159);
      trig_art_vs_evt->Fill(EventNum4Hist/1000.0, neo->size());

      for (auto art1: *neo){
        for (auto art2: *neo){
          if (art1->MMFE8Index() > art2->MMFE8Index()){
            dbcpair = art1->BCID() - art2->BCID();
            if (dbcpair >  4000) dbcpair -= 4096;
            if (dbcpair < -4000) dbcpair += 4096;
            trig_dbc_pairs->Fill(dbcpair);
          }
        }
      }
    }

    //
    // N-1 histograms (residuals)
    //
    if (near06 && clusters_road.size() >= 7){
      for (auto art: *neo){
      
        clusters_N1.Reset();
        track_N1.Reset();

        for (auto clus: clusters_road)
          if (art->MMFE8() != clus->MMFE8())
            clusters_N1.AddCluster(*clus);
      
        track_N1 = FITTER->Fit(clusters_N1, *GEOMETRY, DATA->mm_EventNum);
        ibo      = GEOMETRY->Index(art->MMFE8());

        // using z-half
        // residual = GEOMETRY->GetResidualX(MMCluster(MMHit(art->MMFE8(), art->VMM(), art->VMMChannel(), DATA->RunNum)), track_N1);

        // using z-board
        xtrack = track_N1.SlopeX()*zboard[art->MMFE8Index()] + track_N1.ConstX();
        residual = xpos(art, GEOMETRY) - xtrack;

        track_N1_board_vs_residual_art->Fill(ibo, residual);
      }
    }

    trig_art->Fill(neo->size());
    trig_mm->Fill(clusters_road.size());
    trig_eventmm->Fill(EventNum4Hist/1000.0);
    trig_theta->Fill(thetaTP*180/3.14159);

    // GBT event displays!
    //if (near24 && thetaTP*180/3.14159 < -10 && false)
    //  std::cout << "GBT Event " << neo->EventNumGBT() << " Win " << neo->BCIDWindow() << std::endl;

    // is the ART the first strip in the cluster?
    // ---------------------------------------------------------------
    std::string thatsme;
    int matched = 0;
    int artHIT = 0, artPDO = 0, artBCI = 0;
    ptrdiff_t artposHIT = 0, artposPDO = 0, artposBCI = 0;
    std::vector<int> hits = {}, pdos = {}, bcis = {};
    int allARTfirst = 1;

    if (near24 && thetaTP*180/3.14159 < -10){

      for (auto art: *neo){

        matched = 0;
        pdos.clear();
        hits.clear();
        bcis.clear();
        artHIT = -1;
        artPDO = -1;
        artBCI = -1;
        if (debug)
          std::cout << std::endl;

        for (auto clus: clusters_road){

          if (art->MMFE8Index() != clus->MMFE8Index())
            continue;

          if (debug){
            std::cout << "Evt " << evt << std::endl;
            std::cout << Form("ART Board %i VMM %i CH %2i | Theta %.1f", 
                              art->MMFE8Index(), art->VMM(), 
                              (int)(art->VMMChannel()), thetaTP*180/3.14159)
                      << std::endl;
          }

          // collect attributes of the cluster
          for (auto hit: *clus){
            hits.push_back((int)(hit->Channel()));
            pdos.push_back(hit->PDO());
            bcis.push_back(hit->BCID());
          }

          // match to ART
          for (auto hit: *clus){
            if (fabs(art->Channel() - hit->Channel()) < 0.1){
              matched = 1;
              artHIT = (int)(hit->Channel());
              artPDO = hit->PDO();
              artBCI = hit->BCID();
            }
          }

          // debug
          if (debug){
            ich = 0;
            for (auto hit: *clus){
              thatsme = (art->MMFE8()==clus->MMFE8() && fabs(art->Channel() - hit->Channel()) < 0.1) ? "<----- match ART" : "";
              std::cout << Form("MMF Board %i VMM %i CH %2i | Hit %i BC %4i TDO %2i PDO %2i %s", 
                                hit->MMFE8Index(), hit->VMM(), (int)(hit->VMMChannel()),
                                ich, hit->BCID(), hit->TDO(), hit->PDO(),
                                thatsme.c_str())
                        << std::endl;
              ich++;
            }
          }

          if (!matched)
            continue;

          std::sort(   bcis.begin(), bcis.end());
          std::reverse(pdos.begin(), pdos.end());

          // sort by strip number -- see my notebook for why the boards are as follows
          if (art->MMFE8Index() == 2 || art->MMFE8Index() == 3 || art->MMFE8Index() == 4 || art->MMFE8Index() == 5)
            std::sort(hits.begin(), hits.end());
          else
            std::reverse(hits.begin(), hits.end());

          // find the ART position within the cluster
          artposHIT = std::find(hits.begin(), hits.end(), artHIT) - hits.begin();
          artposPDO = std::find(pdos.begin(), pdos.end(), artPDO) - pdos.begin();
          artposBCI = std::find(bcis.begin(), bcis.end(), artBCI) - bcis.begin();

          if (debug){
            std::cout << "ART HIT position: " << artposHIT << std::endl;
            std::cout << "ART PDO position: " << artposPDO << std::endl;
            std::cout << "ART BC  position: " << artposBCI << std::endl;
          }

        }

        if (!matched)
          continue;

        trig_artpos_bci_vs_clussize->Fill(artposBCI, pdos.size());
        trig_artpos_pdo_vs_clussize->Fill(artposPDO, pdos.size());
        trig_artpos_hit_vs_clussize->Fill(artposHIT, hits.size());
        if (art->MMFE8Index() == 0) trig_artpos_hit_vs_clussize_board0->Fill(artposHIT, hits.size());
        if (art->MMFE8Index() == 1) trig_artpos_hit_vs_clussize_board1->Fill(artposHIT, hits.size());
        if (art->MMFE8Index() == 2) trig_artpos_hit_vs_clussize_board2->Fill(artposHIT, hits.size());
        if (art->MMFE8Index() == 3) trig_artpos_hit_vs_clussize_board3->Fill(artposHIT, hits.size());
        if (art->MMFE8Index() == 4) trig_artpos_hit_vs_clussize_board4->Fill(artposHIT, hits.size());
        if (art->MMFE8Index() == 5) trig_artpos_hit_vs_clussize_board5->Fill(artposHIT, hits.size());
        if (art->MMFE8Index() == 6) trig_artpos_hit_vs_clussize_board6->Fill(artposHIT, hits.size());
        if (art->MMFE8Index() == 7) trig_artpos_hit_vs_clussize_board7->Fill(artposHIT, hits.size());

        if (artHIT > 0)
          allARTfirst = 0;

      }
    }

    // the BCID distributions
    if (near24 && thetaTP*180/3.14159 < -10){

      // the BCID distributions: ART
      for (auto art: *neo){
        dbcid = art->BCID() - neo->BCIDEarliest();
        dbcid = (dbcid < 0) ? dbcid+4096 : dbcid;
        trig_art_bcid->Fill(dbcid);
      }

      // the BCID distributions: MMFE clusters
      int min_bcid = 0;
      std::vector<int> bcids = {};
      for (auto clus: clusters_road)
        for (auto hit: *clus)
          bcids.push_back(hit->BCID());
      min_bcid = *std::min_element(bcids.begin(), bcids.end());
      for (auto bc: bcids)
        mmfe_hit_bcid->Fill(bc - min_bcid);

      // the BCID distributions: MMFE earliest strip within cluster
      bcids.clear();
      int min_bcid_cluster = 9999;
      for (auto clus: clusters_road){
        min_bcid_cluster = 9999;
        for (auto hit: *clus)
          if (hit->BCID() < min_bcid_cluster)
            min_bcid_cluster = hit->BCID();
        bcids.push_back(min_bcid_cluster);
      }
      min_bcid = *std::min_element(bcids.begin(), bcids.end());
      for (auto bc: bcids)
        mmfe_fir_bcid->Fill(bc - min_bcid);
    }
    
    // X,Y resolution

    // collect channels
    std::vector<double> mm_x_x = {}, mm_x_u = {}, mm_x_v = {}, mm_z_xx = {}, mm_z_uv = {};
    std::vector<double> tp_x_x = {}, tp_x_u = {}, tp_x_v = {}, tp_z_xx = {}, tp_z_uv = {};
    for (auto art: *neo){
      if (art->isX()) tp_x_x.push_back(xpos(art, GEOMETRY));
      if (art->isU()) tp_x_u.push_back(xpos(art, GEOMETRY));
      if (art->isV()) tp_x_v.push_back(xpos(art, GEOMETRY));

      if (art->isX())               tp_z_xx.push_back(zpos(art, GEOMETRY));
      if (art->isU() || art->isV()) tp_z_uv.push_back(zpos(art, GEOMETRY));
    }        
    for (auto clus: clusters_road){
      if (clus->isX()) mm_x_x.push_back(xpos(clus, GEOMETRY));
      if (clus->isU()) mm_x_u.push_back(xpos(clus, GEOMETRY));
      if (clus->isV()) mm_x_v.push_back(xpos(clus, GEOMETRY));

      if (clus->isX())                mm_z_xx.push_back(zpos(clus, GEOMETRY));
      if (clus->isU() || clus->isV()) mm_z_uv.push_back(zpos(clus, GEOMETRY));
    }        
    
    double avg_tp_x = 0, avg_mm_x = 0;
    double avg_tp_u = 0, avg_mm_u = 0;
    double avg_tp_v = 0, avg_mm_v = 0;
    double cot_th = cos(1.5*3.14159/180)/sin(1.5*3.14159/180);

    double avg_tp_z_uv = 0, avg_tp_z_xx = 0;
    double avg_mm_z_uv = 0, avg_mm_z_xx = 0;

    double duv = 0, dz = 0;
    double tp_x = 0, mm_x = 0;
    double tp_y = 0, mm_y = 0;

    //if (tp_x_x.size() >= 2 && tp_x_u.size() == 2 && tp_x_v.size() == 2 &&
    //    mm_x_x.size() >= 2 && mm_x_u.size() == 2 && mm_x_v.size() == 2 && near24){
    if (near24){

      avg_tp_x = std::accumulate(tp_x_x.begin(), tp_x_x.end(), 0.0) / (double)(tp_x_x.size());
      avg_tp_u = std::accumulate(tp_x_u.begin(), tp_x_u.end(), 0.0) / (double)(tp_x_u.size());
      avg_tp_v = std::accumulate(tp_x_v.begin(), tp_x_v.end(), 0.0) / (double)(tp_x_v.size());

      avg_mm_x = std::accumulate(mm_x_x.begin(), mm_x_x.end(), 0.0) / (double)(mm_x_x.size());
      avg_mm_u = std::accumulate(mm_x_u.begin(), mm_x_u.end(), 0.0) / (double)(mm_x_u.size());
      avg_mm_v = std::accumulate(mm_x_v.begin(), mm_x_v.end(), 0.0) / (double)(mm_x_v.size());

      avg_tp_z_uv = std::accumulate(tp_z_uv.begin(), tp_z_uv.end(), 0.0) / (double)(tp_z_uv.size());
      avg_mm_z_uv = std::accumulate(mm_z_uv.begin(), mm_z_uv.end(), 0.0) / (double)(mm_z_uv.size());

      avg_tp_z_xx = std::accumulate(tp_z_xx.begin(), tp_z_xx.end(), 0.0) / (double)(tp_z_xx.size());
      avg_mm_z_xx = std::accumulate(mm_z_xx.begin(), mm_z_xx.end(), 0.0) / (double)(mm_z_xx.size());

      tp_x = avg_tp_x;
      tp_y = cot_th/2 * (avg_tp_v - avg_tp_u);
      mm_x = track.SlopeX()*avg_tp_z_xx + track.ConstX();
      mm_y = track.SlopeY()*avg_tp_z_uv + track.ConstY();

      //std::cout << Form("%6i | TP = %7.3f | MM@MM = %7.3f | MM@TP = %7.3f", evt, avg_tp_x,
      //                  track.SlopeX()*avg_mm_z_xx + track.ConstX(),
      //                  track.SlopeX()*avg_tp_z_xx + track.ConstX()) << std::endl;

      trig_x->Fill(tp_x);
      trig_y->Fill(tp_y);
      mmfe_x->Fill(mm_x);
      mmfe_y->Fill(mm_y);

      trig_vs_mmfe_x->Fill(tp_x - mm_x);
      trig_vs_mmfe_y->Fill(tp_y - mm_y);
      trig_vs_mmfe_x_vs_theta->Fill(thetaTP*180/3.14159, tp_x - mm_x);
      trig_vs_mmfe_y_vs_theta->Fill(thetaTP*180/3.14159, tp_y - mm_y);
    }

    // Y resolution
    if (near24 && nearUV && neo->NU() > 0 && neo->NV() > 0){

      double artU_x = 0, artU_z = 0, artV_x = 0, artV_z = 0;
      double dx_bestxxx  = 9999;
      double dy_bestxxx  = 9999;
      double dy_bestpair = 9999;
      double dy_nominal  = 0;
      double av_uv = 0;

      // in total
      artU_x = DATA->tp_EventTracks.AverageU( *neo, *GEOMETRY);
      artU_z = DATA->tp_EventTracks.AverageZU(*neo, *GEOMETRY);
      artV_x = DATA->tp_EventTracks.AverageV( *neo, *GEOMETRY);
      artV_z = DATA->tp_EventTracks.AverageZV(*neo, *GEOMETRY);
      duv = artU_x - artV_x;
      dz  = artV_z - artU_z;
      tp_y = -1 * cot_th * (duv + dz*neo->MxLocal()) / 2;
      tp_y += 217.9;
      mm_y = track.SlopeY()*(artU_z+artV_z)/2 + track.ConstY();
      dy_nominal = mm_y-tp_y;

      // pair-wise
      for (auto artU: *neo){
        for (auto artV: *neo){
          if (!artU->isU()) continue;
          if (!artV->isV()) continue;
          artU_x = xpos_end(artU, GEOMETRY);
          artV_x = xpos_end(artV, GEOMETRY);
          artU_z = zpos(artU, GEOMETRY);
          artV_z = zpos(artV, GEOMETRY);
          duv = artU_x - artV_x;
          dz  = artV_z - artU_z;
          tp_y = -1 * cot_th * (duv + dz*neo->MxLocal()) / 2;
          tp_y += 217.9;
          mm_y = track.SlopeY()*(artU_z+artV_z)/2 + track.ConstY();

          tp_x = track.SlopeX()*(artU_z+artV_z)/2 + track.ConstX();
          //tp_x = tp_x at the zavg of the uv pair;

          av_uv = (artU_x + artV_x) / 2.0;
          //std::cout << Form("%4i | tp_x = %7.1f | av_uv = %7.1f | mmy-tpy = %7.1f | nom = %7.1f", 
          //                  evt, tp_x, av_uv, mm_y-tp_y, dy_nominal) << std::endl;

          if (std::fabs(tp_x-av_uv) < std::fabs(dx_bestxxx)){
            dx_bestxxx = tp_x-av_uv;
            dy_bestxxx = mm_y-tp_y;
          }

          if (std::fabs(mm_y-tp_y) < std::fabs(dy_bestpair))
            dy_bestpair = mm_y-tp_y;

        }
      }
      trig_vs_mmfe_y_nominal ->Fill(dy_nominal);
      trig_vs_mmfe_y_bestxxx ->Fill(dy_bestxxx);
      trig_vs_mmfe_y_bestpair->Fill(dy_bestpair);
      //std::cout << Form("%4i | nom = %7.1f | best = %7.1f | xxx = %7.1f", evt, dy_nominal, dy_bestpair, dy_bestxxx) << std::endl;
    }

    double quad = GEOMETRY->GetQuadraticSumOfResidualsX(clusters_road, track, true);
    double maxresid = 0.0;
    if (near24)                     trig_refsumx2_all->Fill(quad);
    if (near24 && dtheta*1000 > 20) trig_refsumx2_bad->Fill(quad);

    // reference residuals
    residuals.clear();
    for (auto clus: clusters_road)
      residuals.push_back(std::fabs(GEOMETRY->GetResidualX(*clus, track)));
    maxresid = *std::max_element(residuals.begin(), residuals.end());

    // compare full cluster fit with only fitting X planes
    double slope_x   = calculate_slope(mm_x_x, mm_z_xx);
    double dtheta_mm = thetaMM - atan(slope_x);

    // guess at TP inefficiency
    int mmx_minus_tpx = (int)(mm_x_x.size()) - (int)(tp_x_x.size());

    // compare the ART hit with MM
    int board_in_track = 0;
    int maxclustersize = 0;
    int n_unmatched = 0;
    std::vector<double> strips = {};
    for (auto art: *neo){
      if (! art->isX())
        continue;
      board_in_track = 0;
      for (auto clus: clusters_road){
        if (art->MMFE8() != clus->MMFE8())
          continue;
        board_in_track = 1;
        strips.clear();
        for (auto str: *clus)
          strips.push_back(str->Channel());
        maxclustersize = std::max(maxclustersize, (int)(clus->size()));
        if (art->Channel() < *std::min_element(strips.begin(), strips.end())-5.1 ||
            art->Channel() > *std::max_element(strips.begin(), strips.end())+5.1){
          n_unmatched += 1;
        }
      }
      if (! board_in_track)
        n_unmatched += 1;
    }

    // how far is the ART from the strip closest to the track?
    double min_dx = 0.0;
    double min_ch = 0.0;
    double art_dx = 0.0;
    double max_dx = 0.0;
    for (auto art: *neo){
      if (! art->isX())
        continue;
      for (auto clus: clusters_road){
        if (art->MMFE8() != clus->MMFE8())
          continue;
        plane = GEOMETRY->Get(art->MMFE8Index());

        min_dx = 999.0;
        for (auto str: *clus){
          residual = std::fabs(plane.GetResidualX(str->Channel(), track));
          if (residual < min_dx){
            min_dx = residual;
            min_ch = str->Channel();
          }
        }
        art_dx = std::fabs(art->Channel() - min_ch);
        if (art_dx > max_dx)
          max_dx = art_dx;
      }
    }

    if (near24){
      trig_nx_missing->Fill(mmx_minus_tpx);
      trig_dx_in_clus->Fill(max_dx);
      trig_refdtheta->Fill(dtheta_mm*1000);
      trig_refmaxresidual->Fill(maxresid);
      trig_dtheta_vs_refmaxresidual->Fill(dtheta*1000, maxresid);
      trig_dtheta_vs_refdtheta     ->Fill(dtheta*1000, dtheta_mm*1000);
      trig_dtheta_vs_dx_in_clus    ->Fill(dtheta*1000, max_dx);
      trig_dtheta_vs_maxclustersize->Fill(dtheta*1000, maxclustersize);
    }

    if (near24 && maxresid < 1.2 && dtheta_mm*1000 < 3 && max_dx < 6 && mmx_minus_tpx == 0){
      if (n_unmatched == 0) trig_dtheta_unmatched_0->Fill(dtheta*1000);
      if (n_unmatched == 1) trig_dtheta_unmatched_1->Fill(dtheta*1000);
      if (n_unmatched >= 2) trig_dtheta_unmatched_2->Fill(dtheta*1000);
    }

    // summarize
    if (near24){
      trig_dtheta_0->Fill(dtheta*1000);
      if (maxresid < 1.2){
        trig_dtheta_1->Fill(dtheta*1000);
        if (dtheta_mm*1000 < 3){
          trig_dtheta_2->Fill(dtheta*1000);
          if (n_unmatched == 0){
            trig_dtheta_3->Fill(dtheta*1000);
            if (max_dx < 6){
              trig_dtheta_4->Fill(dtheta*1000);
              if (mmx_minus_tpx == 0){
                trig_dtheta_5->Fill(dtheta*1000);
                if      (neo->NX() == 2) trig_dtheta_2X->Fill(dtheta*1000);
                else if (neo->NX() == 3) trig_dtheta_3X->Fill(dtheta*1000);
                else if (neo->NX() == 4) trig_dtheta_4X->Fill(dtheta*1000);
              }
              else
                trig_dtheta_5_else->Fill(dtheta*1000);
            }
            else
              trig_dtheta_4_else->Fill(dtheta*1000);
          }
          else
            trig_dtheta_3_else->Fill(dtheta*1000);
        }
        else
          trig_dtheta_2_else->Fill(dtheta*1000);
      }
      else
        trig_dtheta_1_else->Fill(dtheta*1000);
    }

    // sadness?
    if (false &&

        // near24 &&
        // std::fabs(dtheta*1000) > 15 &&
        // maxresid < 1.2 &&
        // dtheta_mm*1000 < 3 &&
        // n_unmatched == 0 &&
        // max_dx < 6 &&
        // mmx_minus_tpx == 0 &&

        // near64 && nearUV && near04 &&
        // clusters_road.size() >= 7   &&
        // std::fabs(thetaTP*180/3.14159) > 10 && 
        // clusters_all.size()  >= 10  &&

        // evt < 2000

        std::fabs(thetaMM*180/3.14159) > 10
        ){

      // std::cout << Form("\"%05d\",", DATA->mm_EventNum) << std::endl;

      track_tp.Reset();
      clusters_tp.Reset();

      for (auto art: *neo)
        clusters_tp.AddCluster(MMCluster(MMHit(art->MMFE8(), art->VMM(), art->VMMChannel(), DATA->RunNum)));
      track_tp = FITTER->Fit(clusters_tp, *GEOMETRY, DATA->mm_EventNum);

      fout->cd("event_displays");
      can = Plot_Track2D(Form("track2D_%05d_MMall", DATA->mm_EventNum), track, *GEOMETRY, &clusters_all);
      can->Write(); 
      delete can;
      can = Plot_Track2D(Form("track2D_%05d_MMfit", DATA->mm_EventNum), track, *GEOMETRY, &clusters_road);
      can->Write(); 
      delete can;
      can = Plot_Track2D(Form("track2D_%05d_TPfit", DATA->mm_EventNum), track_tp, *GEOMETRY, &clusters_tp);
      can->Write(); 
      delete can;
    }


  }

  std::cout << std::endl;

  // write to file
  fout->cd();
  fout->mkdir("histograms");
  fout->cd("histograms");

  trig_refsumx2_all->Write();
  trig_refsumx2_bad->Write();
  trig_refmaxresidual->Write();
  trig_refdtheta->Write();
  trig_nx_missing->Write();
  trig_dx_in_clus->Write();

  trig_dtheta_0->Write();
  trig_dtheta_1->Write();
  trig_dtheta_2->Write();
  trig_dtheta_3->Write();
  trig_dtheta_4->Write();
  trig_dtheta_5->Write();

  trig_dtheta_2X->Write();
  trig_dtheta_3X->Write();
  trig_dtheta_4X->Write();

  trig_dtheta_0_else->Write();
  trig_dtheta_1_else->Write();
  trig_dtheta_2_else->Write();
  trig_dtheta_3_else->Write();
  trig_dtheta_4_else->Write();
  trig_dtheta_5_else->Write();

  trig_dtheta_vs_refmaxresidual->Write();
  trig_dtheta_vs_refdtheta->Write();
  trig_dtheta_vs_dx_in_clus->Write();
  trig_dtheta_vs_maxclustersize->Write();

  trig_dtheta_unmatched_0->Write();
  trig_dtheta_unmatched_1->Write();
  trig_dtheta_unmatched_2->Write();

  scint_bot_vs_top->Write();

  track_param_x->Write();
  track_param_y->Write();
  track_angle_6->Write();
  track_angle_7->Write();
  track_angle_8->Write();
  track_angle_N_denom->Write();
  track_angle_N_numer->Write();
  track_hits_vs_evt->Write();
  track_hits_vs_evt_fid->Write();
  track_clusmult_vs_theta->Write();

  track_exp_hit->Write();
  track_obs_hit->Write();
  track_vs_time->Write();

  track_mx_vs_cx_0123_4567->Write();
  track_my_vs_cy_0123_4567->Write();

  track_mx_vs_cx_0256_1347->Write();
  track_my_vs_cy_0256_1347->Write();

  track_mx_vs_cx_0236_1457->Write();
  track_my_vs_cy_0236_1457->Write();

  track_x_vs_y_0123_4567->Write();
  track_x_vs_y_0256_1347->Write();
  track_x_vs_y_0236_1457->Write();

  track_scint_vs_time->Write();

  track_slope_x_scint_19->Write();
  track_slope_x_scint_20->Write();
  track_slope_x_scint_21->Write();

  track_angle_x_scint_19->Write();
  track_angle_x_scint_20->Write();
  track_angle_x_scint_21->Write();

  for (auto hist: track_N1_theta_x_vs_residual)
    hist->Write();
  for (auto hist: track_N1_theta_y_vs_residual)
    hist->Write();
  track_N1_board_vs_residual->Write();
  track_N1_board_vs_residual_10deg->Write();
  track_N1_board_vs_residual_art->Write();

  strip_timediff_vs_pdo->Write();
  strip_position_vs_board->Write();
  for (auto hist: strip_dbc_vs_channel)    hist->Write();
  for (auto hist: strip_charge_vs_channel) hist->Write();
  for (auto hist: strip_pdo_vs_channel)    hist->Write();
  for (auto hist: strip_tdo_vs_channel)    hist->Write();
  for (auto hist: strip_event_vs_channel)  hist->Write();
  //for (auto hist: strip_tdoc_vs_channel)   hist->Write();
  //for (auto hist: strip_time_vs_channel)   hist->Write();
  //for (auto hist: strip_zpos_vs_channel)   hist->Write();
    
  clus_vs_board->Write();
  hits_vs_board->Write();
  dups_vs_board->Write();
  dups_vs_channel->Write();
  hits_per_clus_vs_board->Write();
  timediff_vs_board->Write();
  timediff_track->Write();
  timediff_vs_charge->Write();

  if (do_trigger){
    trig_per_event->Write();
    trig_ambig->Write();
    trig_dtheta_vsNX->Write();
    trig_dtheta_vsNX_ok->Write();
    trig_mm->Write();
    trig_art->Write();
    trig_eventmm->Write();
    trig_dbc_vs_N->Write();
    trig_dbc_vs_theta->Write();
    trig_dbc_vs_evt->Write();
    trig_dbc_pairs->Write();
    trig_art_vs_evt->Write();
    trig_theta->Write();

    trig_art_bcid->Write();
    mmfe_hit_bcid->Write();
    mmfe_fir_bcid->Write();

    trig_dbc_scint->Write();

    mmfe_x->Write();
    mmfe_y->Write();
    mmfe_xy_1u1v->Write();
    mmfe_y_1u1v->Write();
    trig_x->Write();
    trig_y->Write();
    trig_vs_mmfe_x->Write();
    trig_vs_mmfe_y->Write();
    trig_vs_mmfe_x_vs_theta->Write();
    trig_vs_mmfe_y_vs_theta->Write();
    trig_vs_mmfe_y_nominal->Write();
    trig_vs_mmfe_y_bestxxx->Write();
    trig_vs_mmfe_y_bestpair->Write();

    trig_artwin_art->Write();
    trig_artwin_mmfe_mimic->Write();
    trig_artwin_mmfe_first->Write();

    trig_artpos_hit_vs_clussize->Write();
    trig_artpos_hit_vs_clussize_board0->Write();
    trig_artpos_hit_vs_clussize_board1->Write();
    trig_artpos_hit_vs_clussize_board2->Write();
    trig_artpos_hit_vs_clussize_board3->Write();
    trig_artpos_hit_vs_clussize_board4->Write();
    trig_artpos_hit_vs_clussize_board5->Write();
    trig_artpos_hit_vs_clussize_board6->Write();
    trig_artpos_hit_vs_clussize_board7->Write();
    trig_artpos_bci_vs_clussize->Write();
    trig_artpos_tdo_vs_clussize->Write();
    trig_artpos_pdo_vs_clussize->Write();

    trig_dtheta_all_NX->Write();
    trig_dtheta_all_2X->Write();
    trig_dtheta_all_3X->Write();
    trig_dtheta_all_4X->Write();
    trig_dx_vs_theta_0->Write();
    trig_dx_vs_theta_1->Write();
    trig_dx_vs_theta_6->Write();
    trig_dx_vs_theta_7->Write();
    trig_dtheta_vs_theta_all   ->Write();
    trig_dtheta_vs_theta_near24->Write();
    trig_dtheta_vs_theta_near12->Write();
    trig_dtheta_nearby_NX->Write();
    trig_dtheta_nearby_2X->Write();
    trig_dtheta_nearby_3X->Write();
    trig_dtheta_nearby_4X->Write();
    trig_dtheta_nearby_15deg_NX->Write();
    trig_dtheta_nearby_15deg_2X->Write();
    trig_dtheta_nearby_15deg_3X->Write();
    trig_dtheta_nearby_15deg_4X->Write();
    trig_dtheta_nearby_15deg_3X_not0->Write();
    trig_dtheta_nearby_15deg_3X_not1->Write();
    trig_dtheta_nearby_15deg_3X_not6->Write();
    trig_dtheta_nearby_15deg_3X_not7->Write();
    trig_dtheta_near64_NX->Write();
    trig_dtheta_near48_NX->Write();
    trig_dtheta_near36_NX->Write();
    trig_dtheta_near24_NX->Write();
    trig_dtheta_near16_NX->Write();
    trig_dtheta_near12_NX->Write();
    trig_dtheta_near08_NX->Write();
    trig_dtheta_near04_NX->Write();
    trig_dtheta_near64_3X->Write();
    trig_dtheta_near48_3X->Write();
    trig_dtheta_near36_3X->Write();
    trig_dtheta_near24_3X->Write();
    trig_dtheta_near16_3X->Write();
    trig_dtheta_near12_3X->Write();
    trig_dtheta_near08_3X->Write();
    trig_dtheta_near04_3X->Write();
    trig_dtheta_near64_4X->Write();
    trig_dtheta_near48_4X->Write();
    trig_dtheta_near36_4X->Write();
    trig_dtheta_near24_4X->Write();
    trig_dtheta_near16_4X->Write();
    trig_dtheta_near12_4X->Write();
    trig_dtheta_near08_4X->Write();
    trig_dtheta_near04_4X->Write();
  }

  fout->Close();    
}


