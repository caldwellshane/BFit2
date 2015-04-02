// Always enclose header file contents with these ifndef/endif directives.
#ifndef _bdn_histograms_h
#define _bdn_histograms_h
#include "TH1.h"
#include "TH2.h"

#ifndef _bdn_histograms_cxx
#define EXTERNAL extern
#else
#define EXTERNAL
#endif

void book_histograms();

// Binning constants
//////////////////////////////////////////////////////////////////////////////////////////
// TDC Singles
static const Double_t	tMin			= -1000.0; // ns
static const Double_t	tMax			= 25500.0; // ns
static const Int_t		tBins			= 26500;   // bin = 1 ns
// ADC Singles
static const Double_t	aMin			= 0.0;     // ADC channel
static const Double_t	aMax			= 4096.0;  // ADC channel
static const Int_t		aBins			= 4096;    // bin = 1 ADC channel
// Time-Of-Flight (tof)
static const Double_t	TOFMin			= tMin; // ns
static const Double_t 	TOFMax			= tMax; // ns
static const Int_t		TOFBins			= 2*tBins; // bin = 0.5 ns due to averaging of 2 integers (t_dEa & t_dEb)
// MCP Maps (raw data on each axis)
static const Double_t	mapMin			= -1.1;
static const Double_t	mapMax			= 1.1;
static const Int_t		mapBins			= 1100;
// Ion Speed (v)
static const Double_t	vMin			= -1000.0;	// mm/us
static const Double_t	vMax			=  5000.0;	// mm/us
static const Int_t		vBins			=  6000;	// bin = 1 mm/us
// Inverse speed (vInv)
static const Double_t	vInvMin			= -0.1;  	// us/mm
static const Double_t	vInvMax			=  0.5;  	// us/mm
static const Int_t		vInvBins		=  6000; 	// bin = .001 us/mm
// Neutron energy (En)
static const Double_t	EnMin			= -1000.0;	// keV
static const Double_t	EnMax			=  5000.0;	// keV
static const Int_t		EnBins			=  6000;	// bin = 1 keV
// High-gain gamma energy
static const Double_t	eMin			= -100.0; // keV
static const Double_t	eMax			= 5000.0; // keV
static const Int_t		eBins			= 5100;   // bin = 1 kev
// Low-gain gamma energy
static const Double_t	eMinLG			= -100.0; // keV
static const Double_t	eMaxLG			= 10000.0; // keV
static const Int_t		eBinsLG			= 10100;   // bin = 1 kev
// Events vs Cycle Time
static const Double_t 	tCycMin 		= -1000.0;  // ms
static const Double_t 	tCycMax 		= 301000.0; // ms
static const Int_t		tCycBins		= 302000;   // bin = 1 ms
// Events vs RF phase
static const Double_t	rfPhaseMin		= -1.1; // #
static const Double_t	rfPhaseMax		= 1.1; // #
static const Int_t		rfPhaseBins		= 2200;
static const Double_t	wholeRFCycle	= 1.0; // range of rf_phase (not necessarily 2*pi)

// Histograms
//////////////////////////////////////////////////////////////////////////////////////////
// TDC Singles
EXTERNAL TH1I *ht_B_dEa;
EXTERNAL TH1I *ht_B_dEb;
EXTERNAL TH1I *ht_B_E;
EXTERNAL TH1I *ht_L_dEa;
EXTERNAL TH1I *ht_L_dEb;
EXTERNAL TH1I *ht_L_E;
EXTERNAL TH1I *ht_R_mcp;
EXTERNAL TH1I *ht_R_ge;
EXTERNAL TH1I *ht_T_mcp;
EXTERNAL TH1I *ht_T_ge;
EXTERNAL TH1I *ht_rf;
// Plastic timings
EXTERNAL TH1I *ht_B_dE;
EXTERNAL TH1I *ht_L_dE;
EXTERNAL TH1I *ht_B_dEdiff;
EXTERNAL TH1I *ht_L_dEdiff;
EXTERNAL TH1I *ht_B_dE_E;
EXTERNAL TH1I *ht_L_dE_E;
EXTERNAL TH1I *ht_B_dEmin;
EXTERNAL TH1I *ht_L_dEmin;
// ADC Singles
	//TH1I *ha_bkgd_R_ge;
	//TH1I *ha_bkgd_T_ge;
EXTERNAL TH1I *ha_B_dEa;
EXTERNAL TH1I *ha_B_dEb;
EXTERNAL TH1I *ha_B_dEsum;
EXTERNAL TH1I *ha_B_E;
EXTERNAL TH1I *ha_L_dEa;
EXTERNAL TH1I *ha_L_dEb;
EXTERNAL TH1I *ha_L_dEsum;
EXTERNAL TH1I *ha_L_E;
EXTERNAL TH1I *ha_R_mcpA;
EXTERNAL TH1I *ha_R_mcpB;
EXTERNAL TH1I *ha_R_mcpC;
EXTERNAL TH1I *ha_R_mcpD;
EXTERNAL TH1I *ha_R_mcpE;
EXTERNAL TH1I *ha_R_mcpSum;
EXTERNAL TH1I *ha_R_ge;
EXTERNAL TH1I *ha_R_ge_highE;
EXTERNAL TH1I *ha_T_mcpA;
EXTERNAL TH1I *ha_T_mcpB;
EXTERNAL TH1I *ha_T_mcpC;
EXTERNAL TH1I *ha_T_mcpD;
EXTERNAL TH1I *ha_T_mcpE;
EXTERNAL TH1I *ha_T_mcpSum;
EXTERNAL TH1I *ha_T_ge;
EXTERNAL TH1I *ha_T_ge_highE;
// Corrected MCP post spectra
EXTERNAL TH1I *ha_R_mcpA_corr;
EXTERNAL TH1I *ha_R_mcpB_corr;
EXTERNAL TH1I *ha_R_mcpC_corr;
EXTERNAL TH1I *ha_R_mcpD_corr;
EXTERNAL TH1I *ha_R_mcpE_corr;
EXTERNAL TH1I *ha_R_mcpSum_corr;
EXTERNAL TH1I *ha_T_mcpA_corr;
EXTERNAL TH1I *ha_T_mcpB_corr;
EXTERNAL TH1I *ha_T_mcpC_corr;
EXTERNAL TH1I *ha_T_mcpD_corr;
EXTERNAL TH1I *ha_T_mcpE_corr;
EXTERNAL TH1I *ha_T_mcpSum_corr;
// Time-Of-Flight (tof)
EXTERNAL TH1I *h_tof;
EXTERNAL TH1I *h_tof_LT;
EXTERNAL TH1I *h_tof_LR;
EXTERNAL TH1I *h_tof_BT;
EXTERNAL TH1I *h_tof_BR;
EXTERNAL TH1I *h_bkgd_tof;
EXTERNAL TH1I *h_bkgd_tof_LT;
EXTERNAL TH1I *h_bkgd_tof_LR;
EXTERNAL TH1I *h_bkgd_tof_BT;
EXTERNAL TH1I *h_bkgd_tof_BR;
EXTERNAL TH1I *h_E_tof;
EXTERNAL TH1I *h_E_tof_LT;
EXTERNAL TH1I *h_E_tof_LR;
EXTERNAL TH1I *h_E_tof_BT;
EXTERNAL TH1I *h_E_tof_BR;
EXTERNAL TH1I *h_bkgd_E_tof;
EXTERNAL TH1I *h_bkgd_E_tof_LT;
EXTERNAL TH1I *h_bkgd_E_tof_LR;
EXTERNAL TH1I *h_bkgd_E_tof_BT;
EXTERNAL TH1I *h_bkgd_E_tof_BR;
EXTERNAL TH1I *h_ge_tof;
EXTERNAL TH1I *h_ge_tof_RT;
EXTERNAL TH1I *h_ge_tof_RR;
EXTERNAL TH1I *h_ge_tof_TT;
EXTERNAL TH1I *h_ge_tof_TR;
EXTERNAL TH1I *h_dE_ge_tof;
EXTERNAL TH1I *h_dE_ge_tof_LT;
EXTERNAL TH1I *h_dE_ge_tof_LR;
EXTERNAL TH1I *h_dE_ge_tof_BT;
EXTERNAL TH1I *h_dE_ge_tof_BR;
EXTERNAL TH1I *h_dE_E_tof;
EXTERNAL TH1I *h_dE_E_tof_BB;
EXTERNAL TH1I *h_dE_E_tof_BL;
EXTERNAL TH1I *h_dE_E_tof_LB;
EXTERNAL TH1I *h_dE_E_tof_LL;
// Ion speed (v)
EXTERNAL TH1I *h_v;
EXTERNAL TH1I *h_v_LT;
EXTERNAL TH1I *h_v_LR;
EXTERNAL TH1I *h_v_BT;
EXTERNAL TH1I *h_v_BR;
EXTERNAL TH1I *h_bkgd_v;
EXTERNAL TH1I *h_bkgd_v_LT;
EXTERNAL TH1I *h_bkgd_v_LR;
EXTERNAL TH1I *h_bkgd_v_BT;
EXTERNAL TH1I *h_bkgd_v_BR;
// Inverse speed (vInv)
EXTERNAL TH1I *h_vInv;
EXTERNAL TH1I *h_vInv_LT;
EXTERNAL TH1I *h_vInv_LR;
EXTERNAL TH1I *h_vInv_BT;
EXTERNAL TH1I *h_vInv_BR;
EXTERNAL TH1I *h_bkgd_vInv;
EXTERNAL TH1I *h_bkgd_vInv_LT;
EXTERNAL TH1I *h_bkgd_vInv_LR;
EXTERNAL TH1I *h_bkgd_vInv_BT;
EXTERNAL TH1I *h_bkgd_vInv_BR;
// Neutron energy (En)
EXTERNAL TH1I *h_En;
EXTERNAL TH1I *h_En_LT;
EXTERNAL TH1I *h_En_LR;
EXTERNAL TH1I *h_En_BT;
EXTERNAL TH1I *h_En_BR;
EXTERNAL TH1I *h_bkgd_En;
EXTERNAL TH1I *h_bkgd_En_LT;
EXTERNAL TH1I *h_bkgd_En_LR;
EXTERNAL TH1I *h_bkgd_En_BT;
EXTERNAL TH1I *h_bkgd_En_BR;
// MCP maps
EXTERNAL TH1I *h_R_mcpX;
EXTERNAL TH1I *h_R_mcpY;
EXTERNAL TH1I *h_T_mcpX;
EXTERNAL TH1I *h_T_mcpY;
EXTERNAL TH2I *h_T_mcpMap;
EXTERNAL TH2I *h_R_mcpMap;
EXTERNAL TH2I *h_T_mcpMapPhys;
EXTERNAL TH2I *h_R_mcpMapPhys;
EXTERNAL TH2I *h_T_mcpMapPhysFidArea;
EXTERNAL TH2I *h_R_mcpMapPhysFidArea;
EXTERNAL TH2I *h_T_mcpMapPhys_3post;
EXTERNAL TH2I *h_R_mcpMapPhys_3post;
EXTERNAL TH2I *h_T_mcpMapPhysFidArea_3post;
EXTERNAL TH2I *h_R_mcpMapPhysFidArea_3post;
EXTERNAL TH2I *h_T_mcpMap_all;
EXTERNAL TH2I *h_R_mcpMap_all;
EXTERNAL TH2I *h_T_mcpMap_post50;
EXTERNAL TH2I *h_R_mcpMap_post50;
EXTERNAL TH2I *h_T_mcpMap_post100;
EXTERNAL TH2I *h_R_mcpMap_post100;
EXTERNAL TH2I *h_T_mcpMap_post200;
EXTERNAL TH2I *h_R_mcpMap_post200;
EXTERNAL TH2I *h_T_mcpMap_post250;
EXTERNAL TH2I *h_R_mcpMap_post250;
EXTERNAL TH2I *h_bkgd_T_mcpMap;
EXTERNAL TH2I *h_bkgd_R_mcpMap;
EXTERNAL TH2I *h_bkgd_T_mcpMapPhys;
EXTERNAL TH2I *h_bkgd_R_mcpMapPhys;
EXTERNAL TH2I *h_bkgd_T_mcpMapPhys_3post;
EXTERNAL TH2I *h_bkgd_R_mcpMapPhys_3post;
EXTERNAL TH2I *h_bkgd_T_mcpMap_all;
EXTERNAL TH2I *h_bkgd_R_mcpMap_all;
EXTERNAL TH2I *h_bkgd_T_mcpMap_post50;
EXTERNAL TH2I *h_bkgd_R_mcpMap_post50;
EXTERNAL TH2I *h_bkgd_T_mcpMap_post100;
EXTERNAL TH2I *h_bkgd_R_mcpMap_post100;
EXTERNAL TH2I *h_bkgd_T_mcpMap_post200;
EXTERNAL TH2I *h_bkgd_R_mcpMap_post200;
EXTERNAL TH2I *h_bkgd_T_mcpMap_post250;
EXTERNAL TH2I *h_bkgd_R_mcpMap_post250;
EXTERNAL TH2I *h_T_zero_mcpMap;
EXTERNAL TH2I *h_R_zero_mcpMap;
EXTERNAL TH2I *h_T_lowTOF_mcpMap;
EXTERNAL TH2I *h_R_lowTOF_mcpMap;
EXTERNAL TH2I *h_T_fast_mcpMap;
EXTERNAL TH2I *h_R_fast_mcpMap;
EXTERNAL TH2I *h_T_slow_mcpMap;
EXTERNAL TH2I *h_R_slow_mcpMap;
EXTERNAL TH2I *h_T_oops_mcpMap;
EXTERNAL TH2I *h_R_oops_mcpMap;
EXTERNAL TH2I *h_LT_zero_mcpMap;
EXTERNAL TH2I *h_LR_zero_mcpMap;
EXTERNAL TH2I *h_BT_zero_mcpMap;
EXTERNAL TH2I *h_BR_zero_mcpMap;
EXTERNAL TH2I *h_LT_lowTOF_mcpMap;
EXTERNAL TH2I *h_LR_lowTOF_mcpMap;
EXTERNAL TH2I *h_BT_lowTOF_mcpMap;
EXTERNAL TH2I *h_BR_lowTOF_mcpMap;
EXTERNAL TH2I *h_LT_fast_mcpMap;
EXTERNAL TH2I *h_LR_fast_mcpMap;
EXTERNAL TH2I *h_BT_fast_mcpMap;
EXTERNAL TH2I *h_BR_fast_mcpMap;
EXTERNAL TH2I *h_LT_slow_mcpMap;
EXTERNAL TH2I *h_LR_slow_mcpMap;
EXTERNAL TH2I *h_BT_slow_mcpMap;
EXTERNAL TH2I *h_BR_slow_mcpMap;
EXTERNAL TH2I *h_LT_oops_mcpMap;
EXTERNAL TH2I *h_LR_oops_mcpMap;
EXTERNAL TH2I *h_BT_oops_mcpMap;
EXTERNAL TH2I *h_BR_oops_mcpMap;
// Events vs Cycle Time
// Declared as TH1D so that they can be multiplied by a TH1D containing per-bin deadtime corrections
EXTERNAL TH2I *h_state_vs_cycle_time;
EXTERNAL TH1D *h_cycles_vs_cycle_time;
EXTERNAL TH1D *h_betas_vs_cycle_time_observed;
EXTERNAL TH1D *h_L_betas_vs_cycle_time_observed;
EXTERNAL TH1D *h_B_betas_vs_cycle_time_observed;
EXTERNAL TH1D *h_all_vs_cycle_time_observed;
EXTERNAL TH1D *h_T_zero_vs_cycle_time_observed;
EXTERNAL TH1D *h_R_zero_vs_cycle_time_observed;
EXTERNAL TH1D *h_zero_vs_cycle_time_observed;
EXTERNAL TH1D *h_T_lowTOF_vs_cycle_time_observed;
EXTERNAL TH1D *h_R_lowTOF_vs_cycle_time_observed;
EXTERNAL TH1D *h_lowTOF_vs_cycle_time_observed;
EXTERNAL TH1D *h_T_fast_vs_cycle_time_observed;
EXTERNAL TH1D *h_R_fast_vs_cycle_time_observed;
EXTERNAL TH1D *h_fast_vs_cycle_time_observed;
EXTERNAL TH1D *h_T_slow_vs_cycle_time_observed;
EXTERNAL TH1D *h_R_slow_vs_cycle_time_observed;
EXTERNAL TH1D *h_slow_vs_cycle_time_observed;
EXTERNAL TH1D *h_T_oops_vs_cycle_time_observed;
EXTERNAL TH1D *h_R_oops_vs_cycle_time_observed;
EXTERNAL TH1D *h_oops_vs_cycle_time_observed;
EXTERNAL TH1D *h_T_CE_vs_cycle_time_observed;
EXTERNAL TH1D *h_R_CE_vs_cycle_time_observed;
EXTERNAL TH1D *h_CE_vs_cycle_time_observed;
EXTERNAL TH1D *h_B_dEE_vs_cycle_time_observed;
EXTERNAL TH1D *h_L_dEE_vs_cycle_time_observed;
EXTERNAL TH1D *h_dEE_vs_cycle_time_observed;
EXTERNAL TH1D *h_LT_bg_vs_cycle_time_observed;
EXTERNAL TH1D *h_LR_bg_vs_cycle_time_observed;
EXTERNAL TH1D *h_BT_bg_vs_cycle_time_observed;
EXTERNAL TH1D *h_BR_bg_vs_cycle_time_observed;
EXTERNAL TH1D *h_bg_vs_cycle_time_observed;
EXTERNAL TH1D *h_LT_bg_gt2MeV_vs_cycle_time_observed;
EXTERNAL TH1D *h_LR_bg_gt2MeV_vs_cycle_time_observed;
EXTERNAL TH1D *h_BT_bg_gt2MeV_vs_cycle_time_observed;
EXTERNAL TH1D *h_BR_bg_gt2MeV_vs_cycle_time_observed;
EXTERNAL TH1D *h_bg_gt2MeV_vs_cycle_time_observed;
// Beta-Gamma, ADC
EXTERNAL TH1I *ha_bg_LT;
EXTERNAL TH1I *ha_bg_LR;
EXTERNAL TH1I *ha_bg_BT;
EXTERNAL TH1I *ha_bg_BR;
EXTERNAL TH1I *ha_sgnl_bg_LT;
EXTERNAL TH1I *ha_sgnl_bg_LR;
EXTERNAL TH1I *ha_sgnl_bg_BT;
EXTERNAL TH1I *ha_sgnl_bg_BR;
EXTERNAL TH1I *ha_bkgd_bg_LT;
EXTERNAL TH1I *ha_bkgd_bg_LR;
EXTERNAL TH1I *ha_bkgd_bg_BT;
EXTERNAL TH1I *ha_bkgd_bg_BR;
// Beta-Gamma, keV (calibrated)
EXTERNAL TH1I *he_bg_LT;
EXTERNAL TH1I *he_bg_LR;
EXTERNAL TH1I *he_bg_BT;
EXTERNAL TH1I *he_bg_BR;
EXTERNAL TH1I *he_bg;
EXTERNAL TH1I *he_sgnl_bg_LT;
EXTERNAL TH1I *he_sgnl_bg_LR;
EXTERNAL TH1I *he_sgnl_bg_BT;
EXTERNAL TH1I *he_sgnl_bg_BR;
EXTERNAL TH1I *he_sgnl_bg;
EXTERNAL TH1I *he_bkgd_bg_LT;
EXTERNAL TH1I *he_bkgd_bg_LR;
EXTERNAL TH1I *he_bkgd_bg_BT;
EXTERNAL TH1I *he_bkgd_bg_BR;
EXTERNAL TH1I *he_bkgd_bg;
// Gamma singles, ADC
//EXTERNAL TH1I *ha_R_ge;
//EXTERNAL TH1I *ha_T_ge;
EXTERNAL TH1I *ha_sgnl_R_ge;
EXTERNAL TH1I *ha_sgnl_T_ge;
EXTERNAL TH1I *ha_bkgd_R_ge;
EXTERNAL TH1I *ha_bkgd_T_ge;
// Gamma singles, keV (calibrated)
EXTERNAL TH1I *he_R_ge;
EXTERNAL TH1I *he_T_ge;
EXTERNAL TH1I *he_ge;
EXTERNAL TH1I *he_sgnl_R_ge;
EXTERNAL TH1I *he_sgnl_T_ge;
EXTERNAL TH1I *he_sgnl_ge;
EXTERNAL TH1I *he_bkgd_R_ge;
EXTERNAL TH1I *he_bkgd_T_ge;
EXTERNAL TH1I *he_bkgd_ge;
//EXTERNAL TH1I *he_R_ge_highE;
//EXTERNAL TH1I *he_T_ge_highE;
//EXTERNAL TH1I *he_ge_highE;
// Events vs RF phase
//EXTERNAL TH1D *ht_rf_phase_observed;
EXTERNAL TH1D *h_all_vs_rf_phase_observed;
EXTERNAL TH1D *h_slow_vs_rf_phase_observed;
EXTERNAL TH1D *h_LT_slow_vs_rf_phase_observed;
EXTERNAL TH1D *h_LR_slow_vs_rf_phase_observed;
EXTERNAL TH1D *h_BT_slow_vs_rf_phase_observed;
EXTERNAL TH1D *h_BR_slow_vs_rf_phase_observed;
EXTERNAL TH1D *h_oops_vs_rf_phase_observed;
EXTERNAL TH1D *h_LT_oops_vs_rf_phase_observed;
EXTERNAL TH1D *h_LR_oops_vs_rf_phase_observed;
EXTERNAL TH1D *h_BT_oops_vs_rf_phase_observed;
EXTERNAL TH1D *h_BR_oops_vs_rf_phase_observed;
EXTERNAL TH1D *h_bkgd_slow_vs_rf_phase_observed;
EXTERNAL TH1D *h_bkgd_LT_slow_vs_rf_phase_observed;
EXTERNAL TH1D *h_bkgd_LR_slow_vs_rf_phase_observed;
EXTERNAL TH1D *h_bkgd_BT_slow_vs_rf_phase_observed;
EXTERNAL TH1D *h_bkgd_BR_slow_vs_rf_phase_observed;
// Special diagnostics for dE-MCP instantaneous coinc peak
EXTERNAL TH1I *ht_B_dE_zero_time_singles;
EXTERNAL TH1I *ht_B_dEa_zero_time_singles;
EXTERNAL TH1I *ht_B_dEb_zero_time_singles;
EXTERNAL TH1I *ht_T_mcp_zero_time_singles;
EXTERNAL TH1I *h_bkgd_tof_dEmin;

#endif
