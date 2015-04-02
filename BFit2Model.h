// Always enclose header file contents with these ifndef/endif directives.
#ifndef _BFitModel_h
#define _BFitModel_h

#include "Rtypes.h"
#include "TString.h"

namespace BFitNamespace {
	
	enum ParIndex { nCyc, dt, DC, r1, r2, r3, p, rho, epsT, epsU, epsV, epsW, epsX, epsY, epsZ, gammaT1, gammaT2, gammaT3, gammaU1, gammaU2, gammaU3 };
	// a[dt]	= bin width in ms
	// a[DC]	= DC detection rate in cycles/ms
	// a[r1]	= Species 1 injection rate in cycles/ms
	// a[r2]	= Species 2 injection rate in cycles/ms
	// a[r3]	= Species 3 injection rate in cycles/ms
	// a[p]		= Fraction of injected ions that is succesfully trapped
	// a[rho]	= Fraction of already-trapped ions that is retained in trap at each capture
	// a[epsT]	= Detection efficiency for trapped pop T
// Note: epsU is multiplied on top of other untrapped eps's
	// a[epsU]	= Detection efficiency for all untrapped pops U
	// a[epsV]	= Detection efficiency for untrapped pop V
	// a[epsW]	= Detection efficiency for untrapped pop W
	// a[epsX]	= Detection efficiency for untrapped pop X
	// a[epsY]	= Detection efficiency for untrapped pop Y
	// a[epsZ]	= Detection efficiency for untrapped pop Z
	// a[gammaT1]	= Non-radioactive decay rate of trapped species 1 population (T1) (1/s)
	// a[gammaT2]	= Non-radioactive decay rate of trapped species 2 population (T2) (1/s)
	// a[gammaT3]	= Non-radioactive decay rate of trapped species 3 population (T3) (1/s)
	// a[gammaU1]	= Non-radioactive decay rate of untrapped species 1 population (U1) (1/s)
	// a[gammaU2]	= Non-radioactive decay rate of untrapped species 2 population (U2) (1/s)
	// a[gammaU3]	= Non-radioactive decay rate of untrapped species 3 population (U3) (1/s)
	// NB: the words are array indices, ie. numbers! If you put p instead of a[p] you will get 
	
// Function used to detect when the fitter changes the parameters
	bool CompareParArrays (const Double_t*, const Double_t*, size_t n, Double_t eps);
	void ComputeParameterDependentVars (Double_t*);
//	void ComputeTimeDependentVars (Double_t*, Double_t*);
//	void ComputePopulations (Double_t*, Double_t*);
	
// Ion populations by type -- Each takes species number as first argument
	Double_t Ttot		(Int_t, Double_t*, Double_t);
	Double_t Utot		(Int_t, Double_t*, Double_t);
	Double_t Vtot		(Int_t, Double_t*, Double_t);
	Double_t Wtot		(Int_t, Double_t*, Double_t);
	Double_t Ztot		(Int_t, Double_t*, Double_t);
	Double_t Xtot		(Int_t, Double_t*, Double_t);
	Double_t Ytot		(Int_t, Double_t*, Double_t);
// For untrapped pops, just the component that grows during capture
	Double_t Ucap		(Int_t, Double_t*, Double_t);
	Double_t Vcap		(Int_t, Double_t*, Double_t);
	Double_t Wcap		(Int_t, Double_t*, Double_t);
	Double_t Zcap		(Int_t, Double_t*, Double_t);
	Double_t Xcap		(Int_t, Double_t*, Double_t);
	Double_t Ycap		(Int_t, Double_t*, Double_t);
	Double_t Ybkgd		(Int_t, Double_t);
// Helper functions to calculate pops
	Double_t sigmaI			(Double_t, Double_t, Int_t);
	Double_t sigmaII		(Double_t, Double_t, Double_t, Int_t);
	Double_t sigmaIII		(Double_t, Double_t, Double_t, Double_t, Int_t);
	Double_t sigmaIV		(Double_t, Double_t, Double_t, Double_t, Double_t, Int_t);
//	Double_t SigmaT			(Double_t, Double_t, Int_t);
//	Double_t SigmaW			(Double_t, Double_t, Double_t, Int_t);
//	Double_t SigmaZ			(Double_t, Double_t, Double_t, Int_t);
//	Double_t SigmaY2		(Double_t, Double_t, Double_t, Double_t, Int_t);
//	Double_t SigmaY2T1		(Double_t, Double_t, Double_t, Double_t, Int_t);
	Double_t Y2Background	(Double_t, Double_t, Double_t);
	Double_t Y3Background	(Double_t, Double_t, Double_t, Double_t);
//	Double_t Y2InitialValue	(Double_t, Double_t*, Double_t, Double_t);
//	Double_t Y3InitialValue	(Double_t, Double_t*, Double_t, Double_t);
	
//	Double_t H_V1_Y2		(Double_t*, Double_t, Double_t);
//	Double_t H_W1_Y2		(Double_t*, Double_t, Double_t);
//	Double_t H_Z1_Y2_Zpart	(Double_t*, Double_t, Double_t);
//	Double_t H_Z1_Y2_Tpart	(Double_t*, Double_t, Double_t);
//	Double_t H_V2_Y3		(Double_t*, Double_t, Double_t);
//	Double_t H_W2_Y3		(Double_t*, Double_t, Double_t);
//	Double_t H_Z2_Y3_Zpart	(Double_t*, Double_t, Double_t);
//	Double_t H_Z2_Y3_Tpart	(Double_t*, Double_t, Double_t);
//	Double_t H_X2_Y3_Xpart	(Double_t*, Double_t, Double_t);
//	Double_t H_X2_Y3_Tpart	(Double_t*, Double_t, Double_t);
//	Double_t H_Y2_Y3		(Double_t*, Double_t, Double_t);
//	Double_t H_U1_U3		(Double_t*, Double_t, Double_t);
	
// Trapped and untrapped populations 1, 2, 3
	Double_t T1 (Double_t*, Double_t*);
/*	Double_t T2 (Double_t*, Double_t*);
	Double_t T3 (Double_t*, Double_t*);
	Double_t U1 (Double_t*, Double_t*);
	Double_t U2 (Double_t*, Double_t*);
	Double_t U3 (Double_t*, Double_t*);
	Double_t V1 (Double_t*, Double_t*);
	Double_t V2 (Double_t*, Double_t*);
	Double_t V3 (Double_t*, Double_t*);
	Double_t W1 (Double_t*, Double_t*);
	Double_t W2 (Double_t*, Double_t*);
	Double_t W3 (Double_t*, Double_t*);
	Double_t X2 (Double_t*, Double_t*);
	Double_t X3 (Double_t*, Double_t*);
	Double_t Y2 (Double_t*, Double_t*);
	Double_t Y3 (Double_t*, Double_t*);
	Double_t Z1 (Double_t*, Double_t*);
	Double_t Z2 (Double_t*, Double_t*);
	Double_t Z3 (Double_t*, Double_t*);
	Double_t All (Double_t*, Double_t*);*/
	
// Functions to compare to data: (obs. decay rate)x(bin dt) = counts by bin
//	Double_t yT (Int_t, Double_t*, Double_t);
//	Double_t yU (Int_t, Double_t*, Double_t);
//	Double_t yV (Int_t, Double_t*, Double_t);
//	Double_t yW (Int_t, Double_t*, Double_t);
//	Double_t yZ (Int_t, Double_t*, Double_t);
//	Double_t yX (Int_t, Double_t*, Double_t);
//	Double_t yY (Int_t, Double_t*, Double_t);
// These are kept as a layer having the signature that works easily in ROOT
// They just call the yT, yV, etc...
	Double_t yDC (Double_t*, Double_t*);
	Double_t yT1 (Double_t*, Double_t*);
	Double_t yT2 (Double_t*, Double_t*);
	Double_t yT3 (Double_t*, Double_t*);
	Double_t yU1 (Double_t*, Double_t*);
	Double_t yU2 (Double_t*, Double_t*);
	Double_t yU3 (Double_t*, Double_t*);
	Double_t yV1 (Double_t*, Double_t*);
	Double_t yV2 (Double_t*, Double_t*);
	Double_t yV3 (Double_t*, Double_t*);
	Double_t yW1 (Double_t*, Double_t*);
	Double_t yW2 (Double_t*, Double_t*);
	Double_t yW3 (Double_t*, Double_t*);
	Double_t yX2 (Double_t*, Double_t*);
	Double_t yX3 (Double_t*, Double_t*);
	Double_t yY2 (Double_t*, Double_t*);
	Double_t yY3 (Double_t*, Double_t*);
	Double_t yZ1 (Double_t*, Double_t*);
	Double_t yZ2 (Double_t*, Double_t*);
	Double_t yZ3 (Double_t*, Double_t*);
	Double_t yAll (Double_t*, Double_t*);
	
// Detection rates (/ms) for calculating integrals --> # of betas
	Double_t rDC (Double_t*, Double_t*);
	Double_t rT1 (Double_t*, Double_t*);
	Double_t rT2 (Double_t*, Double_t*);
	Double_t rT3 (Double_t*, Double_t*);
	Double_t rU1 (Double_t*, Double_t*);
	Double_t rU2 (Double_t*, Double_t*);
	Double_t rU3 (Double_t*, Double_t*);
	Double_t rV1 (Double_t*, Double_t*);
	Double_t rV2 (Double_t*, Double_t*);
	Double_t rV3 (Double_t*, Double_t*);
	Double_t rW1 (Double_t*, Double_t*);
	Double_t rW2 (Double_t*, Double_t*);
	Double_t rW3 (Double_t*, Double_t*);
	Double_t rX2 (Double_t*, Double_t*);
	Double_t rX3 (Double_t*, Double_t*);
	Double_t rY2 (Double_t*, Double_t*);
	Double_t rY3 (Double_t*, Double_t*);
	Double_t rZ1 (Double_t*, Double_t*);
	Double_t rZ2 (Double_t*, Double_t*);
	Double_t rZ3 (Double_t*, Double_t*);
	Double_t rAll (Double_t*, Double_t*);
	
// Offset functions to improve visualization: oTi = yTi + yDC
	Double_t oT1 (Double_t*, Double_t*);
	Double_t oT2 (Double_t*, Double_t*);
	Double_t oT3 (Double_t*, Double_t*);
	Double_t oU1 (Double_t*, Double_t*);
	Double_t oU2 (Double_t*, Double_t*);
	Double_t oU3 (Double_t*, Double_t*);
	
}

// Variables to hold integrals of functions
static Double_t T1_integral = 0.0;
static Double_t T2_integral = 0.0;
static Double_t T3_integral = 0.0;
static Double_t U1_integral = 0.0;
static Double_t U2_integral = 0.0;
static Double_t U3_integral = 0.0;
static Double_t DC_integral = 0.0;
static Double_t All_integral = 0.0;

static Double_t U1_integral_trap_empty = 0.0;
static Double_t U2_integral_trap_empty = 0.0;
static Double_t U3_integral_trap_empty = 0.0;
static Double_t U1_integral_trap_full = 0.0;
static Double_t U2_integral_trap_full = 0.0;
static Double_t U3_integral_trap_full = 0.0;

static Double_t Integral_sum = 0.0;

static Double_t T1_integral_error = 0.0;
static Double_t T2_integral_error = 0.0;
static Double_t T3_integral_error = 0.0;
static Double_t U1_integral_error = 0.0;
static Double_t U2_integral_error = 0.0;
static Double_t U3_integral_error = 0.0;
static Double_t DC_integral_error = 0.0;
static Double_t All_integral_error = 0.0;

static Double_t Integral_sum_error = 0.0;

#endif
