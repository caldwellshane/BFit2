// Always enclose header file contents with these ifndef/endif directives.
#ifndef _BFitModel_h
#define _BFitModel_h

#include "Rtypes.h"
#include "TString.h"
//#include "CSVtoStruct.h"

namespace BFitNamespace {
	
	// NB: these words are array indices, ie. integers!
	// *** If you put p instead of a[p] you will get the integer 6, not the trapping efficiency!
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
	
/////////////////////////////////////////////////////////////////////////////
// Global variables
/////////////////////////////////////////////////////////////////////////////
	extern char parNames[30][5];// = {"nCyc", "dt", "DC", "r1", "r2", "r3", "p", "rho", "epsT", "epsU", "epsV", "epsW", "epsX", "epsY", "epsZ", "gT1", "gT2", "gT3", "gU1", "gU2", "gU3"};
	
	// Structs
	//BDNCase_t	stBDNCases[FILE_ROWS_BDN];
	//BFitCase_t	stBFitCases[FILE_ROWS_BFit];
	//extern Int_t		iBDNCaseIndex, iBFitCaseIndex; // global index to identify case
	//extern Int_t		iNumStructs_BDN, iNumStructs_BFit;
	
	extern Double_t	rateScale; // scale down rates to help Hessian matrix inversion
	// These never change during fitting and are set in BFit()
	extern Int_t	nPars; // number of model parameters
	extern Double_t	iota; // tiny number for avoiding divide-by-0
	extern Double_t	tCap, tBac, tCyc; // cycle times
	extern Double_t	t1, t2, t3; // radioactive half-lives
	//extern Double_t	tZeroArg[1], tCycArg[1]; // specifc time values 1-element array for passing by reference
	extern Int_t	nCycles; // number of cycles in dataset -- don't confuse with parameter indec nCyc
	extern Int_t	nCapMax; // number of injections per cycle
	extern bool		b134sbFlag; // flag for 134sb cases which get special treatment
	extern bool		bEpsXEqualsEpsYFlag; // flag for forcing epsX := epsY
	// These change only when the parameters change during fitting and are set in yAll()
	extern Int_t 		nParChanges; // counts # of times pars have changed
	extern Double_t	*lastPar; // holds most recent paramter values for comparison
	extern Double_t	tT1, tT2, tT3, tU1, tU2, tU3; // modified lifetimes
	extern Double_t	aT1, aT2, aT3, aU1, aU2, aU3; // Decay factors for one capture interval tCap
	extern Double_t	eU1tCyc, eU2tCyc, eU3tCyc; // background decay at end of cycle: exp(-tCyc/tU1), ...
	extern Double_t	tT1U2, tT1U3, tU1U2, tU1U3, tT2U3, tU2U3; // Coefficients used in Y2 and Y3
	extern Double_t	cT1, cU1, cU2, cZT2, cZU2, cZU3, cXT1, cXU2, cXU3, cYU1, cYU2, cYU3, ThetaU, ThetaY; // various products of lifetimes
	//extern Double_t	ST1_1cap, SW11_1cap, SZ11_1cap, ST2_1cap, SW22_1cap, SZ12_1cap, SZ22_1cap; // specific values of the SigmaT, SigmaW, and SigmaZ functions
	extern Double_t	*sigmaT1, *sigmaT2, *sigmaT3;
	extern Double_t	*sigmaV1, *sigmaV2, *sigmaV3;
	extern Double_t	*sigmaW1, *sigmaW2, *sigmaW3;
	extern Double_t	*sigmaZ1, *sigmaZ2, *sigmaZ3;
	extern Double_t			  *sigmaX2, *sigmaX3;
	extern Double_t			  *sigmaY2, *sigmaY3;
	extern Double_t	*sY2v1, *sY2w1, *sY2z1;
	extern Double_t	*sY3v2, *sY3w2, *sY3z2, *sY3x2, *sY3v1, *sY3w1, *sY3z1;
	extern Double_t	*timeOfCapt; 
	extern Double_t	Gamma_T1_U2, Gamma_T2_U3;
	extern Double_t	Gamma_U1_U2, Gamma_U2_U3;
	extern Double_t	*I_V1_Y2, *I_V2_Y3;
	extern Double_t	*I_W1_Y2, *I_W2_Y3;
	extern Double_t	*I_Z1_Y2, *I_Z2_Y3;
	extern Double_t	*I_X1_Y2, *I_X2_Y3;
	extern Double_t			  *I_Y2_Y3;
	extern Double_t	*Y2InitialValues, *Y3InitialValues;
	extern Double_t	ampT1, ampT2, ampT3; // amplitudes of the T pops
	extern Double_t	ampV1, ampV2, ampV3; // amplitudes of the V pops
	extern Double_t	ampW1, ampW2, ampW3; // amplitudes of the W pops
	extern Double_t	ampZ1, ampZ2, ampZ3; // amplitudes of the Z pops
	extern Double_t		   ampX2, ampX3;  // amplitudes of the X pops
	extern Double_t		ampY2ptA, ampY2ptB; // amplitudes for Y2
	extern Double_t	ampY3fromV2, ampY3fromW2, ampY3fromZ2, ampY3fromX2, ampY3fromY2_ST1, ampY3fromY2_SW11, ampY3fromY2_SZ11; // amplitudes for Y3
	extern Double_t	V10, V20, V30, W10, W20, W30, Z10, Z20, Z30, X20, X30, Y20, Y30, U10, U20, U30; // initial value (t=0) for each pop
	// These change for every function call and are set in yAll()
	extern Int_t		nCap; // current injection number
	extern Double_t	eU10, eU20, eU30; // decay from t=0 for background parts: exp(-t/tT1), ...
	extern Double_t	eT1nCap, eT2nCap, eT3nCap, eU1nCap, eU2nCap, eU3nCap; // decay from last injection for capture parts: exp(-(t-tB-(nCap-1)*tA)/tT1), ...
	extern Double_t	T1val, T2val, T3val, U1val, U2val, U3val, V1val, V2val, V3val, W1val, W2val, W3val, Z1val, Z2val, Z3val, X2val, X3val, Y2val, Y3val; // value of each pop
	/////////////////////////////////////////////////////////////////////////////
	extern Int_t sema;	
	
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

static Double_t U1_integral_trap_empty_error = 0.0;
static Double_t U2_integral_trap_empty_error = 0.0;
static Double_t U3_integral_trap_empty_error = 0.0;
static Double_t U1_integral_trap_full_error  = 0.0;
static Double_t U2_integral_trap_full_error  = 0.0;
static Double_t U3_integral_trap_full_error  = 0.0;

#endif
