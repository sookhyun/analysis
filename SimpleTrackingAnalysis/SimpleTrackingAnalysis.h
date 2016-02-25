#ifndef __SIMPLETRACKINGANALYSIS_H__
#define __SIMPLETRACKINGANALYSIS_H__


// --- need to check all these includes...
#include <fun4all/SubsysReco.h>
#include <vector>

class PHCompositeNode;
//class PHG4HoughTransform;
class TH1D;
class TH2D;
//class TProfile;
class TProfile2D;
class RawTower;
class RawCluster;
class RawTowerContainer;
class RawClusterContainer;

class SimpleTrackingAnalysis: public SubsysReco
{

 public:

  SimpleTrackingAnalysis(const std::string &name="SimpleTrackingAnalysis");

  int Init(PHCompositeNode*);
  int process_event(PHCompositeNode*);
  int End(PHCompositeNode*);

  void set_nlayers(unsigned int x) {nlayers = x;}

  void set_verbosity(int x) {verbosity = x;}
  //void set_docalocuts(bool x) {docalocuts = x;}

  std::vector<RawTower*> get_ordered_towers(RawCluster*, RawTowerContainer*);
  std::vector<RawTower*> get_ordered_towers(const RawTowerContainer*);
  std::vector<RawCluster*> get_ordered_clusters(const RawClusterContainer*);
  void inspect_ordered_towers(const std::vector<RawTower*>&);
  void inspect_ordered_towers(const std::vector<RawTower*>&, int);
  void inspect_ordered_towers(const std::vector<RawTower*>&, double, int);


 private:

  // all values taken from sPHENIX pCDR
  // note the somewhat mysterious overlap
  // between the outer edge of the EMC
  // and the inner edge of the inner HCal
  // that shouldn't be a problem here, though
  const double emc_radius_inner =  90.0;
  const double hci_radius_inner = 115.7;
  const double hco_radius_inner = 182.0;
  const double emc_radius_outer = 116.1;
  const double hci_radius_outer = 137.0;
  const double hco_radius_outer = 268.5;

  // PHG4HoughTransform _hough; ERROR

  double magneticfield;

  int verbosity;
  //bool docalocuts;

  // event counter
  unsigned long long nevents;

  // number of layers
  unsigned int nlayers;

  // output histograms ---------------------------------------------------------

  TH2D* _truept_dptoverpt;                    // momentum resolution

  TH2D* _truept_dca;                          // dca resolution

  TH1D* _truept_particles_leaving7Hits;       // pattern reco eff baseline

  TH1D* _truept_particles_recoWithExactHits;  // pattern reco eff by nhits
  TH1D* _truept_particles_recoWithin1Hit;
  TH1D* _truept_particles_recoWithin2Hits;

  TH1D* _truept_particles_recoWithin3Percent; // parttern reco eff by momentum match
  TH1D* _truept_particles_recoWithin4Percent;
  TH1D* _truept_particles_recoWithin5Percent;

  TH1D* _recopt_tracks_all;                   // purity baseline (non-embedded particles)

  TH1D* _recopt_tracks_recoWithExactHits;     // purity by nhit match
  TH1D* _recopt_tracks_recoWithin1Hit;
  TH1D* _recopt_tracks_recoWithin2Hits;

  TH1D* _recopt_tracks_recoWithin3Percent;    // purity by momentum match
  TH1D* _recopt_tracks_recoWithin4Percent;
  TH1D* _recopt_tracks_recoWithin5Percent;

  TH1D* _recopt_tracks_withcalocuts_all;                   // purity baseline (non-embedded particles)

  TH1D* _recopt_tracks_withcalocuts_recoWithExactHits;     // purity by nhit match
  TH1D* _recopt_tracks_withcalocuts_recoWithin1Hit;
  TH1D* _recopt_tracks_withcalocuts_recoWithin2Hits;

  TH1D* _recopt_tracks_withcalocuts_recoWithin3Percent;    // purity by momentum match
  TH1D* _recopt_tracks_withcalocuts_recoWithin4Percent;
  TH1D* _recopt_tracks_withcalocuts_recoWithin5Percent;

  TH2D* _recopt_quality;                      // quality distributions
  TH2D* _truept_quality_particles_recoWithin4Percent;
  TH2D* _recopt_quality_tracks_all;
  TH2D* _recopt_quality_tracks_recoWithin4Percent;

  TH1D* _dx_vertex;                           // vertex resolution
  TH1D* _dy_vertex;
  TH1D* _dz_vertex;



  TH2D* _energy_ratio_emc;
  TH2D* _energy_ratio_hci;
  TH2D* _energy_ratio_hco;
  TH2D* _energy_ratio_hct;
  TH2D* _energy_ratio_tot_dumb;
  TH2D* _energy_ratio_tot_smart;

  TH2D* _energy_ratio_elb_emc;
  TH2D* _energy_ratio_elb_hci;
  TH2D* _energy_ratio_elb_hco;
  TH2D* _energy_ratio_elb_hct;
  TH2D* _energy_ratio_elb_tot_dumb;
  TH2D* _energy_ratio_elb_tot_smart;

  TH2D* _energy_ratio_eub_emc;
  TH2D* _energy_ratio_eub_hci;
  TH2D* _energy_ratio_eub_hco;
  TH2D* _energy_ratio_eub_hct;
  TH2D* _energy_ratio_eub_tot_dumb;
  TH2D* _energy_ratio_eub_tot_smart;

  TH2D* _energy_dphi_emc;
  TH2D* _energy_dphi_hci;
  TH2D* _energy_dphi_hco;

  TH2D* _energy_deta_emc;
  TH2D* _energy_deta_hci;
  TH2D* _energy_deta_hco;

  TProfile2D* _towers_3x3_emc;
  TProfile2D* _towers_5x5_emc;
  TProfile2D* _towers_7x7_emc;
  TProfile2D* _towers_9x9_emc;

  TH2D* _tower_energy_emc[10];
  TH2D* _towersum_energy_emc[10];

  TProfile2D* _towers_3x3_hci;
  TProfile2D* _towers_5x5_hci;
  TProfile2D* _towers_7x7_hci;
  TProfile2D* _towers_9x9_hci;

  TH2D* _tower_energy_hci[10];
  TH2D* _towersum_energy_hci[10];

  TProfile2D* _towers_3x3_hco;
  TProfile2D* _towers_5x5_hco;
  TProfile2D* _towers_7x7_hco;
  TProfile2D* _towers_9x9_hco;

  TH2D* _tower_energy_hco[10];
  TH2D* _towersum_energy_hco[10];

};

#endif // __SIMPLETRACKINGANALYSIS_H__
