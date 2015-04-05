//////////////////////////////////////////////////////////////////////////
//
// BFit2Model by Shane Caldwell
// 2015-01-22
//
// This is an overhaul of BFit designed to make it faster. The model functions,
// at least to begin with, are the same. The key upgrade is that the yAll function
// precomputes dozens of special values that are used throughout the evaluations
// of the component functions. These values live in the global namespace and are
// declared above main() in BFit2.cxx.
//
// All of the precomputed values are parameter-dependent and independent of t. Computing
// them for every t and in every component was a huge waste of time. The fitting
// was just too slow. (Also made for ugly code, but made it easy to build up the 
// model one function at a time.) So the yAll function detects a change in the parameters
// made by the TH1::Fit routine, and then updates these parameter-dependent quantities.
//
// I also tried precomputing a number of time-dependent quantites that appear repeatedly:
// namely the injection number n and the various Exp decay factors. This causes a 
// problem with the way I calculate the initial values, eg. Y20 from Ycap(t=tCyc).
// I ended up going back and just doing these in  each function. It's not the fastest 
// possible way, but it is fast enough and allows an easy construction of the initial
// value in the ComputeParameterDependentVars(), where it belongs.
//
// I have tried to lay out the code in a logical way, but the precomputation of data
// means that the definitions of values are not near their use. If you are reading this
// code to understand its logic, I recommend that you begin by learning
// how the T and V populations work. These are conceptually very simple. The T pops
// have no background component, while V is like a T but with a background part.
// The other populations, except for Y, follow the same logic.
//
// The Y populations are another beast. They differ in that there is no attempt to
// put them into closed form, as I have done for the others using the SigmaT, SigmaW,
// and SigmaZ functions. Instead I use the YxInitialValue() functions to compute a
// Yx that has the right boundary value at each injection, then I sum a bunch of
// these in a 'for' loop in Ycap(). The background parts of the Y pops contain
// feeding from other untrapped pops, as is evidenced by the computations of
// Y20 and Y30 in yAll as well as in the functions YxBackground().
//
// My thesis contains an appendix that describes
// the model in detail, including derivations of all functions.
//
// May the betas be ever in your favor!
//   --Shane
//////////////////////////////////////////////////////////////////////////
// 2015-03-22
//	- Many changes since the above was written.
//	- Y2 and Y3 are now in closed form and appear to be exactly correct.
//	- (Tried several different ways to do this. Ended up with the sigmaI, sigmaII,
//	sigmaIII, sigmaIV algebra to hold the partial-integral information.)
//	- Today I am cleaning out unnecessary code and simplifying the relationships
//	between functions, including eliminating some intermediate functions that
//	weren't needed.
//	- Moved the population functions to their own file: BFit2Populations.cxx.
//	

#include "BFit2Model.h"
#include "CSVtoStruct.h"
#include "TMath.h"
#include "time.h"
//#include "string.h"
#include <iostream>
using namespace std;

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Fitting function -- sum of all components
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Double_t BFitNamespace::yAll(Double_t *t, Double_t *a) {
	using namespace BFitNamespace;
	// Global variables that are constant throughout the fit
	//extern char parNames[30][5];
	extern BFitCase_t stBFitCases[FILE_ROWS_BFit];
	extern Int_t iBFitCaseIndex;
	//extern Int_t nPars;
	//extern Double_t iota, t1, t2, t3;
	// Global variables that depend only on the parameters
	//extern Int_t nParChanges; // counts # of times pars have changed
	//extern Double_t *lastPar; // holds most recent paramter values for comparison
	// Global variables that depend on t
	//extern Double_t T1val, T2val, T3val, U1val, U2val, U3val, V1val, V2val, V3val, W1val, W2val, W3val, Z1val, Z2val, Z3val, X2val, X3val, Y2val, Y3val;
	// Local variables
	static Double_t f;
	static Int_t index; // index of parameter array
	//extern Int_t sema;
	
	Int_t hey = 0;
	while (sema) {
		if (hey==0) {
			printf("hey! blocking\nhey! blocking\nhey! blocking\nhey! blocking\nhey! blocking\nhey! blocking\n");
		}
		hey++;
		usleep(1000);
	}
	sema = 1;
	// When parameters change:
	if (!CompareParArrays(a,lastPar,sizeof(Double_t),iota)) {
		ComputeParameterDependentVars(a);
		memcpy(lastPar,a,nPars*sizeof(Double_t));
		nParChanges++;
		// Print updated paramters
		printf("New pars (%4d): ",nParChanges);
		for (index = 0; index < nPars; index++) {
			if (stBFitCases[iBFitCaseIndex].pbToggle[index]) printf("%s=%.4e ",parNames[index],a[index]);
		}
		printf("\n");
	}
	f = yDC(t,a) + yT1(t,a) + yT2(t,a) + yT3(t,a) + yU1(t,a) + yU2(t,a) + yU3(t,a);
	sema = 0;
	return f;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// CompareParArrays -- used to compare two parameter arrays to detect a change
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bool BFitNamespace::CompareParArrays (const Double_t *par1, const Double_t *par2, size_t n, Double_t eps) {
	if (sizeof(par1) != sizeof(par2)) cout << "Tried to compare to non-equal-length arrays! Results not guaranteed." << endl << endl;
	for (size_t i=0; i<n; i++) if (fabs(par1[i]-par2[i]) > eps) return false;
	return true;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ComputeParameterDependentVars -- when pars change, update parameter-dependent globals
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void BFitNamespace::ComputeParameterDependentVars (Double_t *a) {
	using namespace BFitNamespace;
	using namespace TMath;
// Global variables that are constant throughout the fit
	//extern Int_t nCapMax;
	//extern Double_t iota, tCap, tBac, tCyc, t1, t2, t3;
	//extern Double_t eU1tCyc, eU2tCyc, eU3tCyc;
	//extern bool b134sbFlag;
// Global variables that depend only on the parameters
	//extern Double_t tT1, tT2, tT3, tU1, tU2, tU3;
	//extern Double_t aT1, aT2, aT3, aU1, aU2, aU3;
	//extern Double_t tT1U2, tT1U3, tU1U2, tU1U3, tT2U3, tU2U3;
	//extern Double_t	*timeOfCapt;
	//extern Double_t	*sigmaT1, *sigmaT2, *sigmaT3;
	//extern Double_t	*sigmaV1, *sigmaV2, *sigmaV3;
	//extern Double_t	*sigmaW1, *sigmaW2, *sigmaW3;
	//extern Double_t	*sigmaZ1, *sigmaZ2, *sigmaZ3;
	//extern Double_t           *sigmaX2, *sigmaX3;
	//extern Double_t           *sigmaY2, *sigmaY3;
	//extern Double_t	*sY2v1, *sY2w1, *sY2z1;
	//extern Double_t	*sY3v2, *sY3w2, *sY3z2, *sY3x2, *sY3v1, *sY3w1, *sY3z1;
	//extern Double_t ampT1, ampT2, ampT3;
	//extern Double_t ampV1, ampV2, ampV3;
	//extern Double_t ampW1, ampW2, ampW3;
	//extern Double_t ampZ1, ampZ2, ampZ3;
	//extern Double_t        ampX2, ampX3;
	//extern Double_t V10, V20, V30;
	//extern Double_t W10, W20, W30;
	//extern Double_t Z10, Z20, Z30;
	//extern Double_t      X20, X30;
	//extern Double_t      Y20, Y30;
	//extern Double_t U10, U20, U30;
	static Int_t k;
	
// Special cases:
	if (b134sbFlag) a[gammaT3] = a[gammaT2];
// Offset gammaUi to avoid gammaTi == gammaUi (== 0)
	a[gammaU1] += 2*iota;
	a[gammaU2] += 2*iota;
	a[gammaU3] += 2*iota;
// Modified lifetimes:
	tT1 = 1.0 / ( 1.0/t1 + a[gammaT1]/1000.0 ); // net variable lifetime (1/e) in ms
	tT2 = 1.0 / ( 1.0/t2 + a[gammaT2]/1000.0 ); // net variable lifetime (1/e) in ms
	tT3 = 1.0 / ( 1.0/t3 + a[gammaT3]/1000.0 ); // net variable lifetime (1/e) in ms
	tU1 = 1.0 / ( 1.0/t1 + a[gammaU1]/1000.0 ); // net variable lifetime (1/e) in ms
	tU2 = 1.0 / ( 1.0/t2 + a[gammaU2]/1000.0 ); // net variable lifetime (1/e) in ms
	tU3 = 1.0 / ( 1.0/t3 + a[gammaU3]/1000.0 ); // net variable lifetime (1/e) in ms
	if (tT1==tU1 || tT2==tU2 || tT3==tU3) { // this probably won't be needed, since I put the 1000*iota offset in tUi
	// In the future I hope you won't need the variable lifetimes!
		printf("\n************************************************************");
		printf("\n*** WARNING: Found gamma_Ti = gamma_Ui for some i.");
		printf("\n*** Z_i population will wrongly evaluate to 0, and");
		printf("\n*** Y_{i+1} populations will be too small.");
		printf("\n*** Suggestion: Slightly change either gamma_Ti or gamma_Ui.");
		printf("\n************************************************************\n\n");
	}
	//printf("tT1=%f, tT2=%f, tT3=%f\ntU1=%f, tU2=%f, tU3=%f\n", tT1, tT2, tT3, tU1, tU2, tU3);
// Decay factors
	aT1 = Exp(-tCap/tT1);
	aT2 = Exp(-tCap/tT2);
	aT3 = Exp(-tCap/tT3);
	aU1 = Exp(-tCap/tU1);
	aU2 = Exp(-tCap/tU2);
	aU3 = Exp(-tCap/tU3);
// Special factors
	tT1U2	= tT1*tU2/(tU2-tT1);
	tT1U3	= tT1*tU3/(tU3-tT1);
	tU1U2	= tU1*tU2/(tU2-tU1);
	tU1U3	= tU1*tU3/(tU3-tU1);
	tT2U3	= tT2*tU3/(tU3-tT2);
	tU2U3	= tU2*tU3/(tU3-tU2);
// Amplitudes -- Ti -- factor of 0.001 is to put a[ri] from ions/sec to ions/ms
	ampT1		= 0.001 * a[r1] * tCap * a[p];
	ampT2		= 0.001 * a[r2] * tCap * a[p];
	ampT3		= 0.001 * a[r3] * tCap * a[p];
// Amplitudes -- Vi
	ampV1		= 0.001 * a[r1] * tCap * (1-a[p]);
	ampV2		= 0.001 * a[r2] * tCap * (1-a[p]);
	ampV3		= 0.001 * a[r3] * tCap * (1-a[p]);
// Amplitudes -- Wi
	ampW1		= 0.001 * a[r1] * tCap * a[p] * (1-a[rho]);
	ampW2		= 0.001 * a[r2] * tCap * a[p] * (1-a[rho]);
	ampW3		= 0.001 * a[r3] * tCap * a[p] * (1-a[rho]);
// Amplitudes -- Zi
	ampZ1		= 0.001 * a[r1] * tCap * a[p] * a[gammaT1]/(a[gammaT1]-a[gammaU1]);//(a[gammaT1]+iota)/((a[gammaT1]-a[gammaU1])+iota);
	ampZ2		= 0.001 * a[r2] * tCap * a[p] * a[gammaT2]/(a[gammaT2]-a[gammaU2]);
	ampZ3		= 0.001 * a[r3] * tCap * a[p] * a[gammaT3]/(a[gammaT3]-a[gammaU3]);
	//printf("ampZ1=%f, ampZ2=%f, ampZ3=%f\n", ampZ1, ampZ2, ampZ3);
// Amplitudes -- Xi
	ampX2		= 0.001 * a[r1] * tCap * a[p] * (1/t1) * (tT1*tU2/(tU2-tT1));
	ampX3		= 0.001 * a[r2] * tCap * a[p] * (1/t2) * (tT2*tU3/(tU3-tT2));
// Sigmas:
// See sigmas and _cap functions
	for (k=1; k<=nCapMax; k++) {
	// zero index of these set to zero in BFit()
		sigmaT1[k] = sigmaI(a[rho],aT1,k);
		sigmaT2[k] = sigmaI(a[rho],aT2,k);
		sigmaT3[k] = sigmaI(a[rho],aT3,k);
		sigmaV1[k] = sigmaI(1.0000,aU1,k);
		sigmaV2[k] = sigmaI(1.0000,aU2,k);
		sigmaV3[k] = sigmaI(1.0000,aU3,k);
		sigmaW1[k] = aT1 * sigmaII(a[rho],aT1,aU1,k);
		sigmaW2[k] = aT2 * sigmaII(a[rho],aT2,aU2,k);
		sigmaW3[k] = aT3 * sigmaII(a[rho],aT3,aU3,k);
		sigmaZ1[k] = (aU1-aT1) * sigmaII(a[rho],aT1,aU1,k) + sigmaT1[k]; // corresponds to $\sigma'_{Z1}$ in thesis
		sigmaZ2[k] = (aU2-aT2) * sigmaII(a[rho],aT2,aU2,k) + sigmaT2[k]; // corresponds to $\sigma'_{Z2}$ in thesis
		sigmaZ3[k] = (aU3-aT3) * sigmaII(a[rho],aT3,aU3,k) + sigmaT3[k]; // corresponds to $\sigma'_{Z3}$ in thesis
		sigmaX2[k] = (aU2-aT1) * sigmaII(a[rho],aT1,aU2,k) + sigmaT1[k]; // corresponds to $\sigma'_{X2}$ in thesis
		sigmaX3[k] = (aU3-aT2) * sigmaII(a[rho],aT2,aU3,k) + sigmaT2[k]; // corresponds to $\sigma'_{X3}$ in thesis
		// Components of sigmaY2:
		sY2v1[k]  = tU1U2/t1 * (aU2-aU1) *             sigmaII (1.0000,aU1,aU2,k);
		sY2w1[k]  = tU1U2/t1 * (aU2-aU1) *  aT1      * sigmaIII(a[rho],aT1,aU1,aU2,k);
		sY2z1[k]  = tU1U2/t1 * (aU2-aU1) * (aU1-aT1) * sigmaIII(a[rho],aT1,aU1,aU2,k) + ( tU1U2/t1 * (aU2-aU1) - tT1U2/t1 * (aU2-aT1) ) * sigmaII(a[rho],aT1,aU2,k);
		// Components of sigmaY3:
		sY3v2[k]  = tU2U3/t2 * (aU3-aU2) *             sigmaII (1.0000,aU2,aU3,k);
		sY3w2[k]  = tU2U3/t2 * (aU3-aU2) *  aT2      * sigmaIII(a[rho],aT2,aU2,aU3,k);
		sY3z2[k]  = tU2U3/t2 * (aU3-aU2) * (aU2-aT2) * sigmaIII(a[rho],aT2,aU2,aU3,k) + ( tU2U3/t2 * (aU3-aU2) - tT2U3/t2 * (aU3-aT2) ) * sigmaII(a[rho],aT2,aU3,k);
		sY3x2[k]  = tU2U3/t2 * (aU3-aU2) * (aU2-aT1) * sigmaIII(a[rho],aT1,aU2,aU3,k) + ( tU2U3/t2 * (aU3-aU2) - tT1U3/t2 * (aU3-aT1) ) * sigmaII(a[rho],aT1,aU3,k);
		sY3v1[k]  = tU1U2/t1 *       ( tU2U3/t2 * (aU3-aU2) * (aU2-aU1) * sigmaIII(1.0000,aU1,aU2,aU3,    k) + ( tU2U3/t2 * (aU3-aU2) - tU1U3/t2 * (aU3-aU1) ) * sigmaII (1.0000,aU1,aU3,    k) );
		sY3w1[k]  = tU1U2/t1 * aT1 * ( tU2U3/t2 * (aU3-aU2) * (aU2-aU1) * sigmaIV (a[rho],aT1,aU1,aU2,aU3,k) + ( tU2U3/t2 * (aU3-aU2) - tU1U3/t2 * (aU3-aU1) ) * sigmaIII(a[rho],aT1,aU1,aU3,k) );
		sY3z1[k]  = tU2U3/t2 * (aU3-aU2) * tU1U2/t1 * (aU2-aU1) * (aU1-aT1) * sigmaIV (a[rho],aT1,aU1,aU2,aU3,k)
				  //
				  + tU2U3/t2 * (aU3-aU2) * tU1U2/t1 * (aU2-aU1)             * sigmaIII(a[rho],aT1,aU2,aU3,    k)
				  - tU2U3/t2 * (aU3-aU2) * tT1U2/t1 * (aU2-aT1)             * sigmaIII(a[rho],aT1,aU2,aU3,    k)
				  //
				  + tU2U3/t2 * (aU3-aU2) * tU1U2/t1             * (aU1-aT1) * sigmaIII(a[rho],aT1,aU1,aU3,    k)
				  - tU1U3/t2 * (aU3-aU1) * tU1U2/t1             * (aU1-aT1) * sigmaIII(a[rho],aT1,aU1,aU3,    k)
				  //
				  + tU2U3/t2 * (aU3-aU2) * tU1U2/t1                         * sigmaII (a[rho],aT1,aU3,        k)
				  - tU1U3/t2 * (aU3-aU1) * tU1U2/t1                         * sigmaII (a[rho],aT1,aU3,        k)
				  //
				  - tU2U3/t2 * (aU3-aU2) * tT1U2/t1                         * sigmaII (a[rho],aT1,aU3,        k)
				  + tT1U3/t2 * (aU3-aT1) * tT1U2/t1                         * sigmaII (a[rho],aT1,aU3,        k);
		// sigmaY2 & sigmaY3 -- these have amplitudes and therefore are the full initial values of Y2 and Y3 for each tooth, unlike other sigmas
		sigmaY2[k] = ampV1*sY2v1[k] + ampW1*sY2w1[k] + ampZ1*sY2z1[k];
		sigmaY3[k] = ampV2*sY3v2[k] + ampW2*sY3w2[k] + ampZ2*sY3z2[k] + ampX2*sY3x2[k] + ampV1*sY3v1[k] + ampW1*sY3w1[k] + ampZ1*sY3z1[k];
	}
	
// End of cycle exp factors
	eU1tCyc	= Exp(-tCyc/tU1);
	eU2tCyc	= Exp(-tCyc/tU2);
	eU3tCyc	= Exp(-tCyc/tU3);
// Initial values of populations
	V10 = Vcap(1,a,tCyc) / (1-eU1tCyc);
	V20 = Vcap(2,a,tCyc) / (1-eU2tCyc);
	V30 = Vcap(3,a,tCyc) / (1-eU3tCyc);
	//////////////////////////////////////////
	W10 = Wcap(1,a,tCyc) / (1-eU1tCyc);
	W20 = Wcap(2,a,tCyc) / (1-eU2tCyc);
	W30 = Wcap(3,a,tCyc) / (1-eU3tCyc);
	//////////////////////////////////////////
	Z10 = Zcap(1,a,tCyc) / (1-eU1tCyc);
	Z20 = Zcap(2,a,tCyc) / (1-eU2tCyc);
	Z30 = Zcap(3,a,tCyc) / (1-eU3tCyc);
	//////////////////////////////////////////
	X20 = Xcap(2,a,tCyc) / (1-eU2tCyc);
	X30 = Xcap(3,a,tCyc) / (1-eU3tCyc);
	//////////////////////////////////////////
	U10 = ( V10 + W10 + Z10 );
	Y20 = ( Ycap(2,a,tCyc) + U10 * tU1/t1 * tU2/(tU2-tU1) * ( eU2tCyc - eU1tCyc ) ) / ( 1 - eU2tCyc );
	//////////////////////////////////////////
	U20 = ( V20 + W20 + Z20 + X20 + Y20 );
	Y30 = ( Ycap(3,a,tCyc)
		+ U20 * tU2/t2 * tU3/(tU3-tU2) * ( eU3tCyc - eU2tCyc )
		+ U10 * tU1/t1 * tU2/t2 * tU3/(tU3-tU2)/(tU3-tU1)/(tU2-tU1) * ( tU1 * (tU3-tU2) * eU1tCyc - tU2 * (tU3-tU1) * eU2tCyc + tU3 * (tU2-tU1) * eU3tCyc ) ) / ( 1 - eU3tCyc );
//	U30 = ( V30 + W30 + Z30 + X30 + Y30 );
//	printf("Parameter-dependent values computed.\n");
//	printf("Y20=%f, Y30=%f, Ycap(3,a,tCyc)=%f\n", Y20, Y30, Ycap(3,a,tCyc));
}

//////////////////////////////////////////////////////////////////////////
// "r" functions
// Instantaneous detection rate
//////////////////////////////////////////////////////////////////////////
Double_t BFitNamespace::rDC (Double_t *t, Double_t *a) {
	return a[nCyc]*a[DC]*0.001; // 1/sec to 1/ms
}
Double_t BFitNamespace::rT1 (Double_t *t, Double_t *a) {
	using namespace BFitNamespace;
	//extern Double_t t1;
	return a[nCyc]*a[epsT]*Ttot(1,a,t[0])/t1;
}
Double_t BFitNamespace::rT2 (Double_t *t, Double_t *a) {
	using namespace BFitNamespace;
	//extern Double_t t2;
	return a[nCyc]*a[epsT]*Ttot(2,a,t[0])/t2;
}
Double_t BFitNamespace::rT3 (Double_t *t, Double_t *a) {
	using namespace BFitNamespace;
	//extern Double_t t3;
	return a[nCyc]*a[epsT]*Ttot(3,a,t[0])/t3;
}
Double_t BFitNamespace::rU1 (Double_t *t, Double_t *a) {
	using namespace BFitNamespace;
	//extern Double_t t1;
	return rV1(t,a) + rW1(t,a) + rZ1(t,a);
	//a[nCyc]*a[epsU]*(a[epsV]*Vtot(1,a,t[0]) + a[epsW]*Wtot(1,a,t[0]) + a[epsZ]*Ztot(1,a,t[0]))/t1;
}
Double_t BFitNamespace::rU2 (Double_t *t, Double_t *a) {
	using namespace BFitNamespace;
	//extern Double_t t2;
	return rV2(t,a) + rW2(t,a) + rZ2(t,a) + rX2(t,a) + rY2(t,a);
	//a[nCyc]*a[epsU]*(a[epsV]*Vtot(2,a,t[0]) + a[epsW]*Wtot(2,a,t[0]) + a[epsZ]*Ztot(2,a,t[0]) + a[epsX]*Xtot(2,a,t[0]) + a[epsY]*Ytot(2,a,t[0]))/t2;
}
Double_t BFitNamespace::rU3 (Double_t *t, Double_t *a) {
	using namespace BFitNamespace;
	//extern Double_t t3;
	return rV3(t,a) + rW3(t,a) + rZ3(t,a) + rX3(t,a) + rY3(t,a);
	//a[nCyc]*a[epsU]*(a[epsV]*Vtot(3,a,t[0]) + a[epsW]*Wtot(3,a,t[0]) + a[epsZ]*Ztot(3,a,t[0]) + a[epsX]*Xtot(3,a,t[0]) + a[epsY]*Ytot(3,a,t[0]))/t3;
}
Double_t BFitNamespace::rV1 (Double_t *t, Double_t *a) {
	using namespace BFitNamespace;
	//extern Double_t t1;
	return a[nCyc]*a[epsU]*a[epsV]*Vtot(1,a,t[0])/t1;
}
Double_t BFitNamespace::rV2 (Double_t *t, Double_t *a) {
	using namespace BFitNamespace;
	//extern Double_t t2;
	return a[nCyc]*a[epsU]*a[epsV]*Vtot(2,a,t[0])/t2;
}
Double_t BFitNamespace::rV3 (Double_t *t, Double_t *a) {
	using namespace BFitNamespace;
	//extern Double_t t3;
	return a[nCyc]*a[epsU]*a[epsV]*Vtot(3,a,t[0])/t3;
}
Double_t BFitNamespace::rW1 (Double_t *t, Double_t *a) {
	using namespace BFitNamespace;
	//extern Double_t t1;
	return a[nCyc]*a[epsU]*a[epsW]*Wtot(1,a,t[0])/t1;
}
Double_t BFitNamespace::rW2 (Double_t *t, Double_t *a) {
	using namespace BFitNamespace;
	//extern Double_t t2;
	return a[nCyc]*a[epsU]*a[epsW]*Wtot(2,a,t[0])/t2;
}
Double_t BFitNamespace::rW3 (Double_t *t, Double_t *a) {
	using namespace BFitNamespace;
	//extern Double_t t3;
	return a[nCyc]*a[epsU]*a[epsW]*Wtot(3,a,t[0])/t3;
}
Double_t BFitNamespace::rZ1 (Double_t *t, Double_t *a) {
	using namespace BFitNamespace;
	//extern Double_t t1;
	return a[nCyc]*a[epsU]*a[epsZ]*Ztot(1,a,t[0])/t1;
}
Double_t BFitNamespace::rZ2 (Double_t *t, Double_t *a) {
	using namespace BFitNamespace;
	//extern Double_t t2;
	return a[nCyc]*a[epsU]*a[epsZ]*Ztot(2,a,t[0])/t2;
}
Double_t BFitNamespace::rZ3 (Double_t *t, Double_t *a) {
	using namespace BFitNamespace;
	//extern Double_t t3;
	return a[nCyc]*a[epsU]*a[epsZ]*Ztot(3,a,t[0])/t3;
}
Double_t BFitNamespace::rX2 (Double_t *t, Double_t *a) {
	using namespace BFitNamespace;
	//extern Double_t t2;
	return a[nCyc]*a[epsU]*a[epsX]*Xtot(2,a,t[0])/t2;
}
Double_t BFitNamespace::rX3 (Double_t *t, Double_t *a) {
	using namespace BFitNamespace;
	//extern Double_t t3;
	return a[nCyc]*a[epsU]*a[epsX]*Xtot(3,a,t[0])/t3;
}
Double_t BFitNamespace::rY2 (Double_t *t, Double_t *a) {
	using namespace BFitNamespace;
	//extern Double_t t2;
	return a[nCyc]*a[epsU]*a[epsY]*Ytot(2,a,t[0])/t2;
}
Double_t BFitNamespace::rY3 (Double_t *t, Double_t *a) {
	using namespace BFitNamespace;
	//extern Double_t t3;
	return a[nCyc]*a[epsU]*a[epsY]*Ytot(3,a,t[0])/t3;
}
Double_t BFitNamespace::rAll (Double_t *t, Double_t *a) {
	using namespace BFitNamespace;
	return rDC(t,a) + rT1(t,a) + rT2(t,a) + rT3(t,a) + rU1(t,a) + rU2(t,a) + rU3(t,a);
}

//////////////////////////////////////////////////////////////////////////
// "y" functions
// Functions to plot: (obs. decay rate)x(bin dt) = counts by bin = y
// -- these have the right signature for use in a TF1
// -- must run ComputerParameterDependentVars() before calling these!
/////////////////////////////////////////////////////////////////////
Double_t BFitNamespace::yDC (Double_t *t, Double_t *a) { return a[dt]*BFitNamespace::rDC(t,a); }
Double_t BFitNamespace::yT1 (Double_t *t, Double_t *a) { return a[dt]*BFitNamespace::rT1(t,a); }
Double_t BFitNamespace::yT2 (Double_t *t, Double_t *a) { return a[dt]*BFitNamespace::rT2(t,a); }
Double_t BFitNamespace::yT3 (Double_t *t, Double_t *a) { return a[dt]*BFitNamespace::rT3(t,a); }
Double_t BFitNamespace::yU1 (Double_t *t, Double_t *a) { return a[dt]*BFitNamespace::rU1(t,a); }
Double_t BFitNamespace::yU2 (Double_t *t, Double_t *a) { return a[dt]*BFitNamespace::rU2(t,a); }
Double_t BFitNamespace::yU3 (Double_t *t, Double_t *a) { return a[dt]*BFitNamespace::rU3(t,a); }
Double_t BFitNamespace::yV1 (Double_t *t, Double_t *a) { return a[dt]*BFitNamespace::rV1(t,a); }
Double_t BFitNamespace::yV2 (Double_t *t, Double_t *a) { return a[dt]*BFitNamespace::rV2(t,a); }
Double_t BFitNamespace::yV3 (Double_t *t, Double_t *a) { return a[dt]*BFitNamespace::rV3(t,a); }
Double_t BFitNamespace::yW1 (Double_t *t, Double_t *a) { return a[dt]*BFitNamespace::rW1(t,a); }
Double_t BFitNamespace::yW2 (Double_t *t, Double_t *a) { return a[dt]*BFitNamespace::rW2(t,a); }
Double_t BFitNamespace::yW3 (Double_t *t, Double_t *a) { return a[dt]*BFitNamespace::rW3(t,a); }
Double_t BFitNamespace::yZ1 (Double_t *t, Double_t *a) { return a[dt]*BFitNamespace::rZ1(t,a); }
Double_t BFitNamespace::yZ2 (Double_t *t, Double_t *a) { return a[dt]*BFitNamespace::rZ2(t,a); }
Double_t BFitNamespace::yZ3 (Double_t *t, Double_t *a) { return a[dt]*BFitNamespace::rZ3(t,a); }
Double_t BFitNamespace::yX2 (Double_t *t, Double_t *a) { return a[dt]*BFitNamespace::rX2(t,a); }
Double_t BFitNamespace::yX3 (Double_t *t, Double_t *a) { return a[dt]*BFitNamespace::rX3(t,a); }
Double_t BFitNamespace::yY2 (Double_t *t, Double_t *a) { return a[dt]*BFitNamespace::rY2(t,a); }
Double_t BFitNamespace::yY3 (Double_t *t, Double_t *a) { return a[dt]*BFitNamespace::rY3(t,a); }

//////////////////////////////////////////////////////////////////////////
// "o" functions
// Offset functions to improve visualization: oT1 = yT1 + yDC
//////////////////////////////////////////////////////////////////////////
Double_t BFitNamespace::oT1 (Double_t *t, Double_t *a) { return BFitNamespace::yDC(t,a) + BFitNamespace::yT1(t,a); }
Double_t BFitNamespace::oT2 (Double_t *t, Double_t *a) { return BFitNamespace::yDC(t,a) + BFitNamespace::yT2(t,a); }
Double_t BFitNamespace::oT3 (Double_t *t, Double_t *a) { return BFitNamespace::yDC(t,a) + BFitNamespace::yT3(t,a); }
Double_t BFitNamespace::oU1 (Double_t *t, Double_t *a) { return BFitNamespace::yDC(t,a) + BFitNamespace::yU1(t,a); }
Double_t BFitNamespace::oU2 (Double_t *t, Double_t *a) { return BFitNamespace::yDC(t,a) + BFitNamespace::yU2(t,a); }
Double_t BFitNamespace::oU3 (Double_t *t, Double_t *a) { return BFitNamespace::yDC(t,a) + BFitNamespace::yU3(t,a); }

/*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// sigmas
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Double_t BFitNamespace::sigmaI (Double_t r, Double_t a, Int_t n) {
	using namespace TMath;
	return (1-Power(r*a,n))/(1-r*a);
}
Double_t BFitNamespace::sigmaII (Double_t r, Double_t a1, Double_t a2, Int_t n) {
	using namespace TMath;
	//extern Double_t iota, tCap;
	return 1/(1-r*a1) * ( (1-Power(a2,n-1))/(1-a2) - r*a1 * (Power(a2,n-1)-Power(r*a1,n-1)+iota)/(a2-r*a1+iota) );
}
Double_t BFitNamespace::sigmaIII (Double_t r, Double_t aT1, Double_t aU1, Double_t aU2, Int_t n) {
	using namespace TMath;
	//extern Double_t iota, tCap;
	return 1/(1-r*aT1) * ( 1/(1-aU1) * ( (1-Power(aU2,n-1))/(1-aU2) - (Power(aU2,n-1)-Power(aU1,n-1))/(aU2-aU1) ) - (r*aT1+iota)/(aU1-r*aT1+iota) * ( (Power(aU2,n-1)-Power(aU1,n-1))/(aU2-aU1) - (Power(aU2,n-1)-Power(r*aT1,n-1))/(aU2-r*aT1) ) );
}
Double_t BFitNamespace::sigmaIV (Double_t r, Double_t aT1, Double_t aU1, Double_t aU2, Double_t aU3, Int_t n) {
	using namespace TMath;
	//extern Double_t iota;
	return 1/(1-r*aT1) * ( 1/(1-aU1) * ( 1/(1-aU2) * ( (1-Power(aU3,n-1))/(1-aU3) - (Power(aU3,n-1)-Power(aU2,n-1))/(aU3-aU2) ) - 1/(aU2-aU1) * ( (Power(aU3,n-1)-Power(aU2,n-1))/(aU3-aU2) - (Power(aU3,n-1)-Power(aU1,n-1))/(aU3-aU1) ) ) - (r*aT1+iota)/(aU1-r*aT1+iota) * ( ( 1/(aU2-aU1) * ( (Power(aU3,n-1)-Power(aU2,n-1))/(aU3-aU2) - (Power(aU3,n-1)-Power(aU1,n-1))/(aU3-aU1) ) ) - 1/(aU2-r*aT1) * ( (Power(aU3,n-1)-Power(aU2,n-1))/(aU3-aU2) - (Power(aU3,n-1)-Power(r*aT1,n-1))/(aU3-r*aT1) ) ) );
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// T & U populations
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Double_t BFitNamespace::T1 (Double_t *t, Double_t *a) {
	using namespace BFitNamespace;
	using namespace TMath;
	//extern Double_t tCap, tBac, tCyc, tT1;
	//extern Double_t ampT1, *sigmaT1;
	static Double_t f;
	static Int_t n;
	f = 0.0;
	n = Ceil((t[0]-tBac)/tCap);
	if (tBac <= t[0] && t[0] <= tCyc) {
		f = ampT1 * sigmaT1[n] * Exp(-(t[0]-tBac-(n-1)*tCap)/tT1);
	}
	return f;
}
Double_t BFitNamespace::Ttot (Int_t i, Double_t *a, Double_t tvar) {
	using namespace BFitNamespace;
	using namespace TMath;
	//extern Double_t tCap, tBac, tCyc, tT1, tT2, tT3;
	//extern Double_t ampT1, ampT2, ampT3, *sigmaT1, *sigmaT2, *sigmaT3;
	static Double_t f;
	static Int_t n;
	f = 0.0;
	n = Ceil((tvar-tBac)/tCap);
	if (tBac <= tvar && tvar <= tCyc) {
		if (i==1) f = ampT1 * sigmaT1[n] * Exp(-(tvar-tBac-(n-1)*tCap)/tT1);
		if (i==2) f = ampT2 * sigmaT2[n] * Exp(-(tvar-tBac-(n-1)*tCap)/tT2);
		if (i==3) f = ampT3 * sigmaT3[n] * Exp(-(tvar-tBac-(n-1)*tCap)/tT3);
	}
	return f;
}
Double_t BFitNamespace::Utot (Int_t i, Double_t *a, Double_t tvar) {
	using namespace BFitNamespace;
	using namespace TMath;
	//extern Double_t tBac, tCyc;
	static Double_t f;
	f = 0.0;
	if (tBac <= tvar && tvar <= tCyc) {
		if (i==1) f = Vtot(1,a,tvar) + Wtot(1,a,tvar) + Ztot(1,a,tvar);
		if (i==2) f = Vtot(2,a,tvar) + Wtot(2,a,tvar) + Ztot(2,a,tvar) + Xtot(2,a,tvar) + Ytot(2,a,tvar);
		if (i==3) f = Vtot(3,a,tvar) + Wtot(3,a,tvar) + Ztot(3,a,tvar) + Xtot(3,a,tvar) + Ytot(3,a,tvar);
	}
	return f;
}
Double_t BFitNamespace::Ucap (Int_t i, Double_t *a, Double_t tvar) {
	using namespace BFitNamespace;
	using namespace TMath;
	//extern Double_t tBac, tCyc;
	static Double_t f;
	f = 0.0;
	if (i==1) f = Vcap(1,a,tvar) + Wcap(1,a,tvar) + Zcap(1,a,tvar);
	if (i==2) f = Vcap(2,a,tvar) + Wcap(2,a,tvar) + Zcap(2,a,tvar) + Xcap(2,a,tvar) + Ycap(2,a,tvar);
	if (i==3) f = Vcap(3,a,tvar) + Wcap(3,a,tvar) + Zcap(3,a,tvar) + Xcap(3,a,tvar) + Ycap(3,a,tvar);
	return f;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// V populations
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Double_t BFitNamespace::Vcap (Int_t i, Double_t *a, Double_t tvar) {
	using namespace BFitNamespace;
	using namespace TMath;
	//extern Double_t tCap, tBac, tU1, tU2, tU3;
	//extern Double_t ampV1, ampV2, ampV3, *sigmaV1, *sigmaV2, *sigmaV3;
	static Double_t f;
	static Int_t n;
	f = 0.0;
	n = Ceil((tvar-tBac)/tCap);
	if (i==1) f = ampV1 * sigmaV1[n] * Exp(-(tvar-tBac-(n-1)*tCap)/tU1);//SigmaT(1,tU1,n) * Exp(-(tvar-tBac-(n-1)*tCap)/tU1);
	if (i==2) f = ampV2 * sigmaV2[n] * Exp(-(tvar-tBac-(n-1)*tCap)/tU2);//SigmaT(1,tU2,n) * Exp(-(tvar-tBac-(n-1)*tCap)/tU2);
	if (i==3) f = ampV3 * sigmaV3[n] * Exp(-(tvar-tBac-(n-1)*tCap)/tU3);//SigmaT(1,tU3,n) * Exp(-(tvar-tBac-(n-1)*tCap)/tU3);
	return f;
}
Double_t BFitNamespace::Vtot (Int_t i, Double_t *a, Double_t tvar) {
	//extern Double_t tBac, tCyc, tU1, tU2, tU3, V10, V20, V30;
	static Double_t f;
	f = 0.0; //catch bad values of tvar
	if (0 <= tvar && tvar <= tCyc) {
		if (i==1) f += V10 * TMath::Exp(-tvar/tU1);
		if (i==2) f += V20 * TMath::Exp(-tvar/tU2);
		if (i==3) f += V30 * TMath::Exp(-tvar/tU3);
	}
	if (tBac <= tvar && tvar <= tCyc)
		f += BFitNamespace::Vcap(i,a,tvar);
	return f;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// W populations
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Double_t BFitNamespace::Wcap (Int_t i, Double_t *a, Double_t tvar) {
	using namespace BFitNamespace;
	using namespace TMath;
	//extern Double_t tCap, tBac, tT1, tT2, tT3, tU1, tU2, tU3;
	//extern Double_t ampW1, ampW2, ampW3, *sigmaW1, *sigmaW2, *sigmaW3;
	static Double_t f;
	static Int_t n;
	f = 0.0;
	n = Ceil((tvar-tBac)/tCap);
//	if (i==1) f += ampW1 * SigmaW(a[rho],tT1,tU1,n) * Exp(-(tvar-tBac)/tU1);//* Exp(-(tvar-tBac-(n-1)*tCap)/tU1);
//	if (i==2) f += ampW2 * SigmaW(a[rho],tT2,tU2,n) * Exp(-(tvar-tBac)/tU2);//* Exp(-(tvar-tBac-(n-1)*tCap)/tU2);
//	if (i==3) f += ampW3 * SigmaW(a[rho],tT3,tU3,n) * Exp(-(tvar-tBac)/tU3);//* Exp(-(tvar-tBac-(n-1)*tCap)/tU3);
	if (i==1) f += ampW1 * sigmaW1[n] * Exp(-(tvar-tBac-(n-1)*tCap)/tU1);//* Exp(-(tvar-tBac-(n-1)*tCap)/tU1);
	if (i==2) f += ampW2 * sigmaW2[n] * Exp(-(tvar-tBac-(n-1)*tCap)/tU2);//* Exp(-(tvar-tBac-(n-1)*tCap)/tU2);
	if (i==3) f += ampW3 * sigmaW3[n] * Exp(-(tvar-tBac-(n-1)*tCap)/tU3);//* Exp(-(tvar-tBac-(n-1)*tCap)/tU3);
//	printf("Wcap: i=%d, tvar=%f, par=%f, f=%f\n", i, tvar, SigmaW[a[rho],tT3,tU3,n], f);
	return f;
}
Double_t BFitNamespace::Wtot (Int_t i, Double_t *a, Double_t tvar) {
	//extern Double_t tBac, tCyc, tU1, tU2, tU3, W10, W20, W30;
	static Double_t f;
	f = 0.0; //catch bad values of tvar
	if (0 <= tvar && tvar <= tCyc) {
		if (i==1) f += W10 * TMath::Exp(-tvar/tU1);
		if (i==2) f += W20 * TMath::Exp(-tvar/tU2);
		if (i==3) f += W30 * TMath::Exp(-tvar/tU3);
	}
	if (tBac <= tvar && tvar <= tCyc)
		f += BFitNamespace::Wcap(i,a,tvar);
//	printf("Wtot: i=%d, par=%f, tvar=%f, f=%f\n",i,W30,tvar,f);
	return f;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Z populations
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Double_t BFitNamespace::Zcap (Int_t i, Double_t *a, Double_t tvar) {
	using namespace BFitNamespace;
	using namespace TMath;
	//extern Double_t tCap, tBac, tT1, tT2, tT3, tU1, tU2, tU3;
	//extern Double_t ampT2, ampZ1, ampZ2, ampZ3, *sigmaT1, *sigmaT2, *sigmaT3, *sigmaZ1, *sigmaZ2, *sigmaZ3;
	static Double_t g, f;
	static Int_t n;
	f = 0.0;
	n = Ceil((tvar-tBac)/tCap);
	if (i==1) f += ampZ1 * ( sigmaZ1[n]*Exp(-(tvar-tBac-(n-1)*tCap)/tU1) - sigmaT1[n]*Exp(-(tvar-tBac-(n-1)*tCap)/tT1) );
	if (i==2) f += ampZ2 * ( sigmaZ2[n]*Exp(-(tvar-tBac-(n-1)*tCap)/tU2) - sigmaT2[n]*Exp(-(tvar-tBac-(n-1)*tCap)/tT2) );
	if (i==3) f += ampZ3 * ( sigmaZ3[n]*Exp(-(tvar-tBac-(n-1)*tCap)/tU3) - sigmaT3[n]*Exp(-(tvar-tBac-(n-1)*tCap)/tT3) );
	return f;
}
Double_t BFitNamespace::Ztot (Int_t i, Double_t *a, Double_t tvar) {
	//extern Double_t tBac, tCyc, tU1, tU2, tU3, Z10, Z20, Z30;
	static Double_t f;
	f = 0.0; //catch bad values of tvar
	if (0 <= tvar && tvar <= tCyc) {
		if (i==1) f += Z10 * TMath::Exp(-tvar/tU1);
		if (i==2) f += Z20 * TMath::Exp(-tvar/tU2);
		if (i==3) f += Z30 * TMath::Exp(-tvar/tU3);
	}
	if (tBac <= tvar && tvar <= tCyc)
		f += BFitNamespace::Zcap(i,a,tvar);
	return f;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// X populations
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Double_t BFitNamespace::Xcap (Int_t i, Double_t *a, Double_t tvar) {
	using namespace BFitNamespace;
	using namespace TMath;
	//extern Double_t tCap, tBac, tT1, tT2, tT3, tU1, tU2, tU3;
	//extern Double_t ampX2, ampX3, *sigmaT1, *sigmaT2, *sigmaX2, *sigmaX3;
	static Double_t f;
	static Int_t n;
	f = 0.0;
	n = Ceil((tvar-tBac)/tCap);
	if (i==2) f += ampX2 * ( sigmaX2[n] * Exp(-(tvar-tBac-(n-1)*tCap)/tU2) - sigmaT1[n] * Exp(-(tvar-tBac-(n-1)*tCap)/tT1) );
	if (i==3) f += ampX3 * ( sigmaX3[n] * Exp(-(tvar-tBac-(n-1)*tCap)/tU3) - sigmaT2[n] * Exp(-(tvar-tBac-(n-1)*tCap)/tT2) );
	return f;
}
Double_t BFitNamespace::Xtot (Int_t i, Double_t *a, Double_t tvar) {
	//extern Double_t tBac, tCyc, tU2, tU3, X20, X30;
	static Double_t f;
	f = 0.0; //catch bad values of tvar
	if (0 <= tvar && tvar <= tCyc) {
		if (i==2) f += X20 * TMath::Exp(-tvar/tU2);
		if (i==3) f += X30 * TMath::Exp(-tvar/tU3);
	}
	if (tBac <= tvar && tvar <= tCyc)
		f += BFitNamespace::Xcap(i,a,tvar);
	return f;
}



// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Y populations
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Double_t BFitNamespace::Ytot (Int_t i, Double_t *a, Double_t tvar) {
	//extern Int_t nCap;
	//extern Double_t tBac, tCyc, tU1, tU2, tU3, U10, U20, Y20, Y30;
	static Double_t f;
	f = 0.0;
	if (0 <= tvar && tvar <= tCyc) {
		f += BFitNamespace::Ybkgd(i,tvar);
	//	if (i==2) f += Ybkgd(2,tvar);//Y20 * Exp(-tvar/tU2) + U10 * tU1/t1 * tU2/(tU2-tU1) * ( Exp(-tvar/tU2) - Exp(-tvar/tU1) );
	//	if (i==3) f += Ybkgd(3,tvar);//Y3Background(U10,U20,Y30,tvar);
	}
	if (tBac < tvar && tvar <= tCyc)
		f += BFitNamespace::Ycap(i,a,tvar);
	return f;
}

Double_t BFitNamespace::Ybkgd (Int_t i, Double_t tvar) {
	using namespace TMath;
	//extern Double_t tBac, tCyc, t1, t2, tU1, tU2, tU3, U10, U20, Y20, Y30;
	Double_t f = 0.0;
	if (i==2) f += Y20 * Exp(-tvar/tU2)
				 + U10 * tU1/t1 * tU2/(tU2-tU1) * ( Exp(-tvar/tU2) - Exp(-tvar/tU1) );
	if (i==3) f += Y30 * Exp(-tvar/tU3)
				 + U20 * tU2/t2 * tU3/(tU3-tU2) * ( Exp(-tvar/tU3) - Exp(-tvar/tU2) )
				 + U10 * tU1/t1 * tU2/t2 * tU3/(tU3-tU2)/(tU3-tU1)/(tU2-tU1) * ( tU1*(tU3-tU2)*Exp(-tvar/tU1) - tU2*(tU3-tU1)*Exp(-tvar/tU2) + tU3*(tU2-tU1)*Exp(-tvar/tU3) );
	return f;
}

Double_t BFitNamespace::Ycap (Int_t i, Double_t *a, Double_t tvar) {
	using namespace BFitNamespace;
	using namespace TMath;
	//extern Double_t tCap, tBac, t1, t2, tT1, tT2, tU1, tU2, tU3, aT1, aU1, aU2, aU3, tU1U2, tT1U2, tU2U3, tT2U3, tU1U3, tT1U3;
	//extern Double_t ampV1, ampW1, ampZ1, ampV2, ampW2, ampZ2, ampX2, ampY2;
	//extern Double_t *sigmaT1, *sigmaV1, *sigmaW1, *sigmaZ1, *sigmaT2, *sigmaV2, *sigmaW2, *sigmaZ2, *sigmaX2, *sigmaY2, *sigmaY3, *sY2v1, *sY2w1, *sY2z1, *sY3w1;
	static Double_t tn, f;
	static Int_t n;
	f = 0.0;
	n = Ceil((tvar-tBac)/tCap);
	//if (tvar==tBac) n=1;
	tn = tvar-tBac-(n-1)*tCap;
	if (i==2) {
		f = sigmaY2[n] * Exp(-tn/tU2)
			+ tU1U2/t1 * ( Exp(-tn/tU2) - Exp(-tn/tU1) ) * ( ampV1*sigmaV1[n] + ampW1*sigmaW1[n] + ampZ1*sigmaZ1[n] )
			- tT1U2/t1 * ( Exp(-tn/tU2) - Exp(-tn/tT1) ) *   ampZ1*sigmaT1[n];
		if (IsNaN(f)) printf("t=%f, Y2 =%f\n", tvar, f);
	}
	if (i==3) {
		f = sigmaY3[n] * Exp(-tn/tU3)
		// Components from V2, W2, Z2, X2
			+ tU2U3/t2 * ( Exp(-tn/tU3) - Exp(-tn/tU2) ) * ( ampV2*sigmaV2[n] + ampW2*sigmaW2[n] + ampZ2*sigmaZ2[n] + ampX2*sigmaX2[n] )
			- tT2U3/t2 * ( Exp(-tn/tU3) - Exp(-tn/tT2) ) *   ampZ2*sigmaT2[n]
			- tT1U3/t2 * ( Exp(-tn/tU3) - Exp(-tn/tT1) ) *   ampX2*sigmaT1[n]
		// Components from Y2 <== V1, W1, Z1
			+ ampV1 * tU1U2/t1 * ( tU2U3/t2 * ( Exp(-tn/tU3) - Exp(-tn/tU2) ) * ( sigmaV1[n] + sY2v1[n]*t1/tU1U2 ) - tU1U3/t2 * ( Exp(-tn/tU3) - Exp(-tn/tU1) ) * sigmaV1[n] )
			//
			+ ampW1 * tU1U2/t1 * ( tU2U3/t2 * ( Exp(-tn/tU3) - Exp(-tn/tU2) ) * ( sigmaW1[n] + sY2w1[n]*t1/tU1U2 ) - tU1U3/t2 * ( Exp(-tn/tU3) - Exp(-tn/tU1) ) * sigmaW1[n] )
			//
			+ ampZ1 * tU2U3/t2 * ( Exp(-tn/tU3) - Exp(-tn/tU2) ) * sY2z1[n]
			//
			+ ampZ1 * tU2U3/t2 * ( Exp(-tn/tU3) - Exp(-tn/tU2) ) * tU1U2/t1 * sigmaZ1[n]
			- ampZ1 * tU1U3/t2 * ( Exp(-tn/tU3) - Exp(-tn/tU1) ) * tU1U2/t1 * sigmaZ1[n]
			//
			- ampZ1 * tU2U3/t2 * ( Exp(-tn/tU3) - Exp(-tn/tU2) ) * tT1U2/t1 * sigmaT1[n]
			+ ampZ1 * tT1U3/t2 * ( Exp(-tn/tU3) - Exp(-tn/tT1) ) * tT1U2/t1 * sigmaT1[n];
		//printf("t=%f, Y3=%f\n", tvar, f);
		if (IsNaN(f)) printf("t=%f, Y3 =%f\n", tvar, f);
	}
	
	return f;
}
*/
