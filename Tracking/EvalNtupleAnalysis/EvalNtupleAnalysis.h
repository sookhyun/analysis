#ifndef __EVALNTUPLEANALYSIS_H_
#define __EVALNTUPLEANALYSIS_H_

#include <TFile.h>
#include <TNtuple.h>
#include <TH1D.h>
#include <TEfficiency.h>

class TH1D;

class EvalNtupleAnalysis
{
  public:
    EvalNtupleAnalysis();
    virtual ~EvalNtupleAnalysis() {}

    int init();
    int open(std::string fname);
    int process_ntp_vertex();
    int close();
    int plot_data();
    int write_output(std::string oname);

    void select_gembed(float gembed){gembed_elect = gembed;}

  private:

    int set_branches_vertex(TNtuple* ntuple);
    int set_branches_gpoint(TNtuple* ntuple);
    int set_branches_track(TNtuple* ntuple);
    int set_branches_gtrack(TNtuple* ntuple);

    TFile* infile;
    TFile* ofile;
    TNtuple* ntp_vertex;
    TNtuple* ntp_gpoint;
    TNtuple* ntp_track;
    TNtuple* ntp_gtrack;

    TH1D* h_nevt_pileup;
    TH1D* h_nevt_pileup_maps; 
    TH1D* h_den;
    TH1D* h_num;
    TH1D* h_zvtx_res;

    std::map<unsigned int, unsigned int> npileup_den;
    std::map<unsigned int, unsigned int> npileup_num;

    float event;
    float gvz;
    float vz;
    float gntracksmaps;
    float ntracks;
    float gembed;
    float gembed_elect;
    float gnembed;

    unsigned int curevent;
    unsigned int nevents;
    unsigned int nevents_vtx;
    bool vertex_found;
    

};

#endif
