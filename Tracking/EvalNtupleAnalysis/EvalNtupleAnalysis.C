#include "EvalNtupleAnalysis.h"
#include <cstring>
#include <iostream>
#include <cmath>

#include <TH1D.h>

using namespace std;

EvalNtupleAnalysis::EvalNtupleAnalysis() :
	infile(0),
	ofile(0),
	ntp_vertex(0),
	ntp_gpoint(0),
	ntp_track(0),
	ntp_gtrack(0),
	h_nevt_pileup(0), h_nevt_pileup_maps(0),h_den(0), h_num(0), h_zvtx_res(0),
	npileup_den(), npileup_num(),	
	event(-999.), gvz(-999), vz(-999), gntracksmaps(-999), ntracks(-999), gembed(999), gembed_elect(999),
	gnembed(15), curevent(9999), nevents(0), nevents_vtx(0), vertex_found(false)
{
	init();

}

// init: open a new file, read ntuples, set branches
   
int EvalNtupleAnalysis::open(std::string fname ="test.root"){

    infile = new TFile(fname.c_str());
    if (!infile) {
	cout<<"file "<< fname.c_str()<<" does not exist!"<<endl;
	exit(1);
    }
    ntp_vertex = (TNtuple*) infile->Get("ntp_vertex");
    if (!ntp_vertex) {
	cout<<"ntuple vertex cannot be read! "<<endl;
	return 0;
    }
 
    ntp_gpoint = (TNtuple*) infile->Get("ntp_gpoint");
    if (!ntp_gpoint) {
        cout<<"ntuple gpoint cannot be read! "<<endl;
        return 0;
    }

    ntp_track = (TNtuple*) infile->Get("ntp_track");
    if (!ntp_track) {
        cout<<"ntuple track cannot be read! "<<endl;
        return 0;
    }

    ntp_gtrack = (TNtuple*) infile->Get("ntp_gtrack");
    if (!ntp_gtrack) {
        cout<<"ntuple gtrack cannot be read! "<<endl;
        return 0;
    }

    

  return 1;
}

int EvalNtupleAnalysis::process_ntp_vertex(){

    unsigned int totentries = ntp_vertex->GetEntries();
    cout<<"total number of entries in this file : "<<totentries<<endl;
    set_branches_vertex(ntp_vertex);

    unsigned int npileupevents = 0;
    unsigned int npileupevents_vtx =0;

      for (unsigned int ientry = 0; ientry<totentries; ++ientry){
        ntp_vertex->GetEntry(ientry);
//        cout<<"reading entry: "<<endl;
//	cout<<"event "<<event <<" gvz "<<gvz<<" vz "<< vz<<" gntracksmaps "<<gntracksmaps <<endl;  

        if (gvz!=gvz || gntracksmaps!=gntracksmaps) continue;
	if (gembed_elect >=0 ){
        	if (gembed!=gembed_elect || gntracksmaps<2 || fabs(gvz)>13.0 ) continue;
	        if (event != curevent){
        	        if (curevent!=9999) {
                	// process current event
                        cout<<"processing event "<<curevent<<endl;
                	++nevents;
                	if (vertex_found) ++nevents_vtx;
                	curevent = event;
                	vertex_found = false;
                	}else{
                	curevent = event;
                	}
        	}
		if (fabs(gvz-vz)<0.05) {
			vertex_found=true;
			h_zvtx_res->Fill(gvz-vz);
		}
	}else {
		if (gembed>=0 || gntracksmaps<2 || fabs(gvz)>13.0) continue;
                if (event != curevent){
                        if (curevent!=9999) {
                        // process current event
                        cout<<"processing event "<<curevent<<endl;
			cout<<"npileupevents "<<npileupevents<<" npileupevents_vtx "<<npileupevents_vtx<<endl;	
			cout<<"gnembed  "<<gnembed<<endl;
			npileup_den[npileupevents] += npileupevents;
			npileup_num[npileupevents] += npileupevents_vtx;			
			h_nevt_pileup->Fill(gnembed-1); // subtract 1 for in-time event
			h_nevt_pileup_maps->Fill(npileupevents);
		
                        curevent = event;
                        }else{
                        curevent = event;
                        }
			npileupevents=0;
			npileupevents_vtx=0;
                }
		
		vertex_found=false;
		if (fabs(gvz-vz)<0.05) vertex_found=true;
		++nevents;
		++npileupevents;	
		if (vertex_found){ 
			++nevents_vtx;
			++npileupevents_vtx; 
			h_zvtx_res->Fill(gvz-vz);
		}

		
	}
        //cout<<"cut event "<<curevent<<" event "<<event<<endl;

      }

      cout<<"processed events: "<<nevents<<" , events with vertex found:"<<nevents_vtx<<endl;	

    return 0;
}

int EvalNtupleAnalysis::close(){

    infile->Close();
    return 0; 
}


int EvalNtupleAnalysis::plot_data(){

	unsigned int inpileup =0.;

	for (std::map<unsigned int, unsigned int>::iterator iter = npileup_den.begin();	
		iter != npileup_den.end();
		++iter )
	{
		inpileup =  iter->first;
		h_den->Fill(inpileup, npileup_den[iter->first]);
		h_num->Fill(inpileup, npileup_num[iter->first]);
		cout<<"inpileup "<<inpileup<<" npileup "<<npileup_num[iter->first]<<endl;
	}

	
	return 0;
}

// end : write out histograms
int EvalNtupleAnalysis::write_output(std::string oname){

    TEfficiency* teff = new TEfficiency(*h_num, *h_den);
    teff->SetTitle("Vertexing efficiency for Out-of-time events");
    teff->Draw("q");

    ofile = new TFile(oname.c_str(),"recreate");
    ofile->cd();
    h_nevt_pileup->Write();
    h_nevt_pileup_maps->Write();
    h_den->Write();
    h_num->Write();
    h_zvtx_res->Write();
    teff->Write();

    ofile->Close();
    return 0;
}

int EvalNtupleAnalysis::init(){

        h_nevt_pileup = new TH1D("h_nevents_pileup","Number of pileup events per in-time event",100,0.5,100.5);
        h_nevt_pileup_maps = new TH1D("h_nevents_pileup_maps","Number of pileup events within MAPS acceptance per in-time event",30,0.5,30.5);
        h_den = new TH1D("den","denominator",20,0.5,20.5);
        h_num = new TH1D("num","numerator",20,0.5,20.5);
	h_den->GetXaxis()->SetTitle("number of pile-up events");
	h_zvtx_res = new TH1D("h_zvtx_res","z-vertex resolution",200, -0.05,0.05);

    return 0;
}

int EvalNtupleAnalysis::set_branches_vertex(TNtuple* ntuple){

    ntuple->SetBranchAddress("event",&event);
    ntuple->SetBranchAddress("gvz",&gvz);
    ntuple->SetBranchAddress("vz",&vz);
    ntuple->SetBranchAddress("gntracksmaps",&gntracksmaps);
    ntuple->SetBranchAddress("ntracks",&ntracks);
    ntuple->SetBranchAddress("gembed",&gembed);
    ntuple->SetBranchAddress("gnembed",&gnembed);

    return 0;
}

 
int EvalNtupleAnalysis::set_branches_gpoint(TNtuple* ntuple){

    return 0;
}

int EvalNtupleAnalysis::set_branches_track(TNtuple* ntuple){

    return 0;
}

int EvalNtupleAnalysis::set_branches_gtrack(TNtuple* ntuple){

    return 0;
}



