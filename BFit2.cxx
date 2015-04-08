/*
2013-12-07 Shane Caldwell
	Made from fit_beta_singles.cxx.
	Will modify as necessary to validate against B_monte_carlo output.
2014-02-26
	Making a standalone C++ version of B_fit.cpp, which is a ROOT program.
	Big difference: this program does not output plots to screen. Instead it saves them to a
	subdirectory of the ROOT file that is given to it.
2014-06-29
	Just now noting this: Program now (as of like easrly March) uses the CSVtoStruct.h/cxx
*/

#include <unistd.h>
#include "stdio.h"
#include <iostream>
#include <iomanip>
#include <cstring>
#include "time.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFitResult.h"
#include "TMatrixDSym.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "Fit/FitConfig.h"
//#include "include/bdn.h"
//#include "include/sb135.h"
//#include "bdn_cases.h"
#include "bdnHistograms.h"
#include "CSVtoStruct.h"
#include "BFit2Model.h"
using namespace std;

//#include "TVirtualFitter.h"


/////////////////////////////////////////////////////////////////////////////
// Global variables
/////////////////////////////////////////////////////////////////////////////

// Structs
BDNCase_t	stBDNCases[FILE_ROWS_BDN];
BFitCase_t	stBFitCases[FILE_ROWS_BFit];
Int_t		iBDNCaseIndex, iBFitCaseIndex; // global index to identify case
Int_t		iNumStructs_BDN, iNumStructs_BFit;

namespace BFitNamespace {
char parNames[30][5] = {"nCyc", "dt", "DC", "r1", "r2", "r3", "p", "rho", "epsT", "epsU", "epsV", "epsW", "epsX", "epsY", "epsZ", "gT1", "gT2", "gT3", "gU1", "gU2", "gU3"};

Double_t	rateScale = 0.001; // scale down rates to help Hessian matrix inversion
// These never change during fitting and are set in BFit()
Int_t		nPars; // number of model parameters
Double_t	iota; // tiny number for avoiding divide-by-0
Double_t	tCap, tBac, tCyc; // cycle times
Double_t	t1, t2, t3; // radioactive half-lives
//Double_t	tZeroArg[1], tCycArg[1]; // specifc time values 1-element array for passing by reference
Int_t		nCycles; // number of cycles in dataset -- don't confuse with parameter indec nCyc
Int_t		nCapMax; // number of injections per cycle
bool		b134sbFlag = 0; // flag for 134sb cases which get special treatment
// These change only when the parameters change during fitting and are set in yAll()
Int_t 		nParChanges; // counts # of times pars have changed
Double_t	*lastPar; // holds most recent paramter values for comparison
Double_t	tT1, tT2, tT3, tU1, tU2, tU3; // modified lifetimes
Double_t	aT1, aT2, aT3, aU1, aU2, aU3; // Decay factors for one capture interval tCap
Double_t	eU1tCyc, eU2tCyc, eU3tCyc; // background decay at end of cycle: exp(-tCyc/tU1), ...
Double_t	tT1U2, tT1U3, tU1U2, tU1U3, tT2U3, tU2U3; // Coefficients used in Y2 and Y3
Double_t	cT1, cU1, cU2, cZT2, cZU2, cZU3, cXT1, cXU2, cXU3, cYU1, cYU2, cYU3, ThetaU, ThetaY; // various products of lifetimes
//Double_t	ST1_1cap, SW11_1cap, SZ11_1cap, ST2_1cap, SW22_1cap, SZ12_1cap, SZ22_1cap; // specific values of the SigmaT, SigmaW, and SigmaZ functions
Double_t	*sigmaT1, *sigmaT2, *sigmaT3;
Double_t	*sigmaV1, *sigmaV2, *sigmaV3;
Double_t	*sigmaW1, *sigmaW2, *sigmaW3;
Double_t	*sigmaZ1, *sigmaZ2, *sigmaZ3;
Double_t			  *sigmaX2, *sigmaX3;
Double_t			  *sigmaY2, *sigmaY3;
Double_t	*sY2v1, *sY2w1, *sY2z1;
Double_t	*sY3v2, *sY3w2, *sY3z2, *sY3x2, *sY3v1, *sY3w1, *sY3z1;
Double_t	*timeOfCapt; 
Double_t	Gamma_T1_U2, Gamma_T2_U3;
Double_t	Gamma_U1_U2, Gamma_U2_U3;
Double_t	*I_V1_Y2, *I_V2_Y3;
Double_t	*I_W1_Y2, *I_W2_Y3;
Double_t	*I_Z1_Y2, *I_Z2_Y3;
Double_t	*I_X1_Y2, *I_X2_Y3;
Double_t			  *I_Y2_Y3;
Double_t	*Y2InitialValues, *Y3InitialValues;
Double_t	ampT1, ampT2, ampT3; // amplitudes of the T pops
Double_t	ampV1, ampV2, ampV3; // amplitudes of the V pops
Double_t	ampW1, ampW2, ampW3; // amplitudes of the W pops
Double_t	ampZ1, ampZ2, ampZ3; // amplitudes of the Z pops
Double_t		   ampX2, ampX3;  // amplitudes of the X pops
Double_t		ampY2ptA, ampY2ptB; // amplitudes for Y2
Double_t	ampY3fromV2, ampY3fromW2, ampY3fromZ2, ampY3fromX2, ampY3fromY2_ST1, ampY3fromY2_SW11, ampY3fromY2_SZ11; // amplitudes for Y3
Double_t	V10, V20, V30, W10, W20, W30, Z10, Z20, Z30, X20, X30, Y20, Y30, U10, U20, U30; // initial value (t=0) for each pop
// These change for every function call and are set in yAll()
Int_t		nCap; // current injection number
Double_t	eU10, eU20, eU30; // decay from t=0 for background parts: exp(-t/tT1), ...
Double_t	eT1nCap, eT2nCap, eT3nCap, eU1nCap, eU2nCap, eU3nCap; // decay from last injection for capture parts: exp(-(t-tB-(nCap-1)*tA)/tT1), ...
Double_t	T1val, T2val, T3val, U1val, U2val, U3val, V1val, V2val, V3val, W1val, W2val, W3val, Z1val, Z2val, Z3val, X2val, X3val, Y2val, Y3val; // value of each pop
/////////////////////////////////////////////////////////////////////////////
Int_t sema = 0;
}

// Functions
Double_t intErr (TF1*, Double_t*, Double_t, Double_t);
void HistPrep (TH1*, Int_t, Int_t, char*);
void FuncPrep (TF1*, Double_t*, Int_t, Int_t, Int_t);
int BFit ();

// MAIN FUNCTION
int main (int argc, char *argv[]) {
	time_t tm = time(NULL);
	cout << "Timestamp = " << asctime(localtime(&tm));
	char *csvBDNCases, *csvBFitCases;
	csvBDNCases = "BDNCases.csv_transposed";
	csvBFitCases = "BFitCases.csv_transposed";
	cout << endl << "Importing metadata from CSV files..." << endl;
	iNumStructs_BDN  = CSVtoStruct_BDN  (csvBDNCases, stBDNCases);
	cout << "Imported " << iNumStructs_BDN << " BDN cases" << endl;
	iNumStructs_BFit = CSVtoStruct_BFit (csvBFitCases, stBFitCases);
	cout << "Imported " << iNumStructs_BFit << " BFit cases" << endl;
	iBDNCaseIndex  = FindStructIndex ( stBDNCases,  sizeof(BDNCase_t),  iNumStructs_BDN,  argv[1] );
	iBFitCaseIndex = FindStructIndex ( stBFitCases, sizeof(BFitCase_t), iNumStructs_BFit, argv[2] );
	if ( iBDNCaseIndex == -1 || iBFitCaseIndex == -1 )
	{ // One of the read-ins failed and already printed a message about it
		cout << "How to run this program:" << endl;
		cout << "'./BFit <BDN case code> <B_fit case code>'" << endl;
		cout << "where valid case codes are listed in the CSV files." << endl << endl;
		return -1; // error return
	}
	cout << "Performing BFit with BDN case " << stBDNCases[iBDNCaseIndex].pcsCaseCode << " and BFit case " << stBFitCases[iBFitCaseIndex].pcsCaseCode << endl << endl;
	return BFit(); // return status of BFit
}

int BFit () {
	
// Timer for debugging code
	clock_t timer, timerStart, timerStop;
	timerStart = clock();
	cout << "BFit started. Timer = "<< (Float_t)timerStart/CLOCKS_PER_SEC << " sec." << endl;
	int iReturn = SUCCESS;
	using namespace BFitNamespace;
	using namespace TMath;
	TFile *outfile = new TFile("BFit.root","recreate");
	
// Local copies of the relevant metadata structs
	BDNCase_t  stBDNCase = stBDNCases[iBDNCaseIndex];
	BFitCase_t stBFitCase = stBFitCases[iBFitCaseIndex];
	//printf("Bdn Case index = %d, BFit case index = %d\n", iBDNCaseIndex, iBFitCaseIndex);
	
// Assign global variables
	nPars	= stBFitCase.iNPars;
	tCap	= 1000.0 * stBDNCase.dCaptureTime;	// Time between BPT captures (ms)
	tBac	= 1000.0 * stBDNCase.dBackgroundTime; // Time spent in background measurment, per cycle (ms)
	tCyc	= 1000.0 * stBDNCase.dCycleTime;	// Time between BPT ejections (ms)
	t1		= 1000.0 * stBDNCase.dLifetime1[0]; // radioactive lifetime (1/e) in ms
	t2		= 1000.0 * stBDNCase.dLifetime2[0]; // radioactive lifetime (1/e) in ms
	t3		= 1000.0 * stBDNCase.dLifetime3[0]; // radioactive lifetime (1/e) in ms
	nCycles	= stBDNCase.nCycles; // Number of cycles in dataset
	iota	= 0.0000000001;
//	printf("          iota = %.4e\n",iota);
//	printf("          iota = %32.16f\n",iota);
//	printf("1000000 + iota = %32.16f\n",1000000+iota);
//	return 0;
	nCapMax	= Ceil((tCyc-tBac)/tCap);
	//printf("nCapMax=%d\n",nCapMax);
// Array to hold injection times
	timeOfCapt	= new Double_t [nCapMax+1];
	Int_t k;
	for (k=0; k<nCapMax+1; k++)
		timeOfCapt[k] = tBac + k*tCap;
// Arrays to hold the sigma values in each [index] tooth:
// Zero index not used (+1 element), for notational consistency
	sigmaT1	= new Double_t [nCapMax+1];
	sigmaT2	= new Double_t [nCapMax+1];
	sigmaT3	= new Double_t [nCapMax+1];
	sigmaV1	= new Double_t [nCapMax+1];
	sigmaV2	= new Double_t [nCapMax+1];
	sigmaV3	= new Double_t [nCapMax+1];
	sigmaW1	= new Double_t [nCapMax+1];
	sigmaW2	= new Double_t [nCapMax+1];
	sigmaW3	= new Double_t [nCapMax+1];
	sigmaZ1	= new Double_t [nCapMax+1];
	sigmaZ2	= new Double_t [nCapMax+1];
	sigmaZ3	= new Double_t [nCapMax+1];
	sigmaX2	= new Double_t [nCapMax+1];
	sigmaX3	= new Double_t [nCapMax+1];
	sigmaY2	= new Double_t [nCapMax+1];
	sigmaY3	= new Double_t [nCapMax+1];
	sY2v1	= new Double_t [nCapMax+1];
	sY2w1	= new Double_t [nCapMax+1];
	sY2z1	= new Double_t [nCapMax+1];
	sY3v2	= new Double_t [nCapMax+1];
	sY3w2	= new Double_t [nCapMax+1];
	sY3z2	= new Double_t [nCapMax+1];
	sY3x2	= new Double_t [nCapMax+1];
	sY3v1	= new Double_t [nCapMax+1];
	sY3w1	= new Double_t [nCapMax+1];
	sY3z1	= new Double_t [nCapMax+1];
	sigmaT1[0] = sigmaT2[0] = sigmaT3[0] = 0.0; // zero index not used -- set to zero for definiteness
	sigmaV1[0] = sigmaV2[0] = sigmaV3[0] = 0.0; // zero index not used -- set to zero for definiteness
	sigmaW1[0] = sigmaW2[0] = sigmaW3[0] = 0.0; // zero index not used -- set to zero for definiteness
	sigmaZ1[0] = sigmaZ2[0] = sigmaZ3[0] = 0.0; // zero index not used -- set to zero for definiteness
				sigmaX2[0] = sigmaX3[0] = 0.0; // zero index not used -- set to zero for definiteness
				sigmaY2[0] = sigmaY3[0] = 0.0; // zero index not used -- set to zero for definiteness
	sY2v1[0]  = sY2w1[0]  = sY2z1[0]  = 0.0; // zero index not used -- set to zero for definiteness
	sY3v2[0]  = sY3w2[0]  = sY3z2[0]  = sY3x2[0] = 0.0; // zero index not used -- set to zero for definiteness
	sY3v1[0]  = sY3w1[0]  = sY3z1[0]  = 0.0; // zero index not used -- set to zero for definiteness
	
// Print message
	TString separator = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
	cout << endl << separator << endl;
	cout << "BETA SINGLES MODEL" << endl << separator << endl;
	cout << "Case: " << stBDNCase.pcsCaseCode << endl;
	cout << "File: " << stBDNCase.pcsFilePath << endl;
	cout << "Histogram: " << stBFitCase.pcsHistName << endl;
	printf("Total cycle time\t= %10.3f s\n",	stBDNCase.dCycleTime);
	printf("Background time\t\t= %10.3f s\n",	stBDNCase.dBackgroundTime);
	printf("Capture cycle time\t= %10.3f s\n",	stBDNCase.dCaptureTime);
	printf("Last capture time\t= %10.3f s\n",	stBDNCase.dLastCaptureTime);
	cout << stBDNCase.pcsSpecies1Name; printf(" halflife\t\t= %10.3f s  (lifetime %10.3f s)\n", stBDNCase.dHalfLife1[0], stBDNCase.dLifetime1[0]);
	cout << stBDNCase.pcsSpecies2Name; printf(" halflife\t\t= %10.3f s  (lifetime %10.3f s)\n", stBDNCase.dHalfLife2[0], stBDNCase.dLifetime2[0]);
	cout << stBDNCase.pcsSpecies3Name; printf(" halflife\t\t= %10.3f s  (lifetime %10.3f s)\n", stBDNCase.dHalfLife3[0], stBDNCase.dLifetime3[0]);
	cout << separator << endl;
	
// Get histogram from ROOT file
	TFile *f = new TFile(stBDNCase.pcsFilePath);
	TH1D *h	= (TH1D*)f->Get(stBFitCase.pcsHistName);
	Double_t dBinWidth		= stBFitCase.pdSeed[dt];
	Double_t dNBins			= tCyc/dBinWidth;// # of bins covered by funtion  //h->GetNbinsX();
	Double_t pointsPerBin	= dBinWidth;//20;
	Double_t nPoints		= pointsPerBin * dNBins;
	Double_t dRebinFactor	= dBinWidth/(h->GetBinWidth(1));
//	printf("dBinWidth=%f, GetBinWidth=%f, dRebinFactor=%f\n", dBinWidth, h->GetBinWidth(1), dRebinFactor);
//	Double_t dRebinFactor = stBFitCase.iBinWidth/(h->GetBinWidth(1));
	Double_t binZero  = (1.0/dRebinFactor)*(0.0  - tCycMin) + 1.0;
	Double_t binCap  = (1.0/dRebinFactor)*(tBac - tCycMin) + 1.0;
	Double_t binCycle = (1.0/dRebinFactor)*(tCyc - tCycMin) + 1.0;
	TH1D *h1	= (TH1D*)h->Rebin(dRebinFactor,stBFitCase.pcsHistName);
	TH1D *h2	= (TH1D*)h->Rebin(dRebinFactor,stBFitCase.pcsHistName);
	cout << "HISTOGRAM REBINNED" << endl;
	cout << "   Number of hist bins = " << dNBins << endl;
	cout << "   Number of fn points = " << nPoints << endl;
	cout << "   Number of pts/bin   = " << pointsPerBin << endl;
	cout << separator << endl;
	if (stBFitCase.bDoFit) cout << "Doing fit with option string: " << stBFitCase.pcsOptions << endl;
	else cout << "Not fitting! Drawing functions using parameter seed values." << endl;
	cout << separator << endl;
	
// Define functions
//	Double_t tMax	= 1000.0 * stBDNCase.dCycleTime;
// Species populations
	TF1 *fyDC	= new TF1("fyDC", yDC, 0.0, tCyc, nPars);
	TF1 *fyT1	= new TF1("fyT1", yT1, 0.0, tCyc, nPars);
	TF1 *fyT2	= new TF1("fyT2", yT2, 0.0, tCyc, nPars);
	TF1 *fyT3	= new TF1("fyT3", yT3, 0.0, tCyc, nPars);
	TF1 *fyU1	= new TF1("fyU1", yU1, 0.0, tCyc, nPars);
	TF1 *fyU2	= new TF1("fyU2", yU2, 0.0, tCyc, nPars);
	TF1 *fyU3	= new TF1("fyU3", yU3, 0.0, tCyc, nPars);
//**************************************************************************
// This one fits the data:
	TF1 *fyAll	= new TF1("fyAll",yAll, 0.0, tCyc, nPars);
//**************************************************************************
// Beta rates to be used by TF1::Integral() and TF1::IntegralError()
	TF1 *frDC	= new TF1("frDC", rDC, 0.0, tCyc, nPars);
	TF1 *frT1	= new TF1("frT1", rT1, 0.0, tCyc, nPars);
	TF1 *frT2	= new TF1("frT2", rT2, 0.0, tCyc, nPars);
	TF1 *frT3	= new TF1("frT3", rT3, 0.0, tCyc, nPars);
	TF1 *frU1	= new TF1("frU1", rU1, 0.0, tCyc, nPars);
	TF1 *frU2	= new TF1("frU2", rU2, 0.0, tCyc, nPars);
	TF1 *frU3	= new TF1("frU3", rU3, 0.0, tCyc, nPars);
	TF1 *frAll	= new TF1("frAll",rAll, 0.0, tCyc, nPars);
// Offset functions for plotting
	TF1 *foT1	= new TF1("foT1", oT1, 0.0, tCyc, nPars);
	TF1 *foT2	= new TF1("foT2", oT2, 0.0, tCyc, nPars);
	TF1 *foT3	= new TF1("foT3", oT3, 0.0, tCyc, nPars);
	TF1 *foU1	= new TF1("foU1", oU1, 0.0, tCyc, nPars);
	TF1 *foU2	= new TF1("foU2", oU2, 0.0, tCyc, nPars);
	TF1 *foU3	= new TF1("foU3", oU3, 0.0, tCyc, nPars);
// Initiailize function to fit the data
	char pcsLifetime1ParName[100]; sprintf(pcsLifetime1ParName,"%s radioactive lifetime (1/e)", stBDNCase.pcsSpecies1Name);
	char pcsLifetime2ParName[100]; sprintf(pcsLifetime2ParName,"%s radioactive lifetime (1/e)", stBDNCase.pcsSpecies2Name);
	char pcsLifetime3ParName[100]; sprintf(pcsLifetime3ParName,"%s radioactive lifetime (1/e)", stBDNCase.pcsSpecies3Name);
	fyAll->SetParName(nCyc,"# Cyc's");
	fyAll->SetParName(DC,"DC rate");
	fyAll->SetParName(r1,"Rate 1");
	fyAll->SetParName(r2,"Rate 2");
	fyAll->SetParName(r3,"Rate 3");
	fyAll->SetParName(p,"Capt eff");
	fyAll->SetParName(rho,"Capt ret");
	fyAll->SetParName(epsT,"T eff");
	fyAll->SetParName(epsU,"U eff");
	fyAll->SetParName(epsV,"V eff");
	fyAll->SetParName(epsW,"W eff");
	fyAll->SetParName(epsX,"X eff");
	fyAll->SetParName(epsY,"Y eff");
	fyAll->SetParName(epsZ,"Z eff");
//	fyAll->SetParName(tau1,pcsLifetime1ParName);//"Lifetime 1");//pcsSpecies1ParName);
//	fyAll->SetParName(tau2,pcsLifetime2ParName);//"Lifetime 2");//pcsSpecies2ParName);
//	fyAll->SetParName(tau3,pcsLifetime3ParName);//"Lifetime 3");//pcsSpecies3ParName);
	fyAll->SetParName(gammaT1,"T1 loss");
	fyAll->SetParName(gammaT2,"T2 loss");
	fyAll->SetParName(gammaT3,"T3 loss");
	fyAll->SetParName(gammaU1,"U1 loss");
	fyAll->SetParName(gammaU2,"U2 loss");
	fyAll->SetParName(gammaU3,"U3 loss");
	fyAll->SetParName(dt,"Bin width");
	
// Initial parameter values and initial step sizes
// err contains initial step sizes now, will contain error estimates later.
	Int_t		*tog;
	Double_t	*par;
	Double_t	*err;
	tog	= stBFitCase.pbToggle;
	par = stBFitCase.pdSeed;
	err = stBFitCase.pdStep;
	Int_t index;
	par[nCyc] = nCycles;
// Special cases -- modifications to parameters -- catch right after param import
	if (!strcmp(stBDNCases[iBDNCaseIndex].pcsCaseCode,"134sb01") ||
		!strcmp(stBDNCases[iBDNCaseIndex].pcsCaseCode,"134sb02") ||
		!strcmp(stBDNCases[iBDNCaseIndex].pcsCaseCode,"134sb03") ||
		!strcmp(stBDNCases[iBDNCaseIndex].pcsCaseCode,"134sb0103"))
	{
		b134sbFlag = 1;
		par[gammaT3] = par[gammaT2];
	}
// Initialize global lastPar to values that guarantee computation of parameter-dependent variables in first eval of yAll
	lastPar	= new Double_t [nPars];
	memcpy(lastPar,par,nPars*sizeof(Double_t));
//	printf("\n\n[lastPar = %f, %f, %f, ...]\n\n",lastPar[0],lastPar[1],lastPar[2]);
//	for (index = 0; index < nPars; index++) lastPar[index] = -1.0;
// Initialize parameter-dependent variables
	BFitNamespace::ComputeParameterDependentVars(par);
// Place initial par values into fitting function
// Reduce rate params to agree in size with others, to help Hessian matrix inversion	
	for (index = 0; index < nPars; index++) {
		if (index == DC || index == r1 || index == r2 || index == r3) {
			fyAll->SetParameter(index,rateScale*par[index]);
			fyAll->SetParError(index,rateScale*err[index]);
		}
		else {
			fyAll->SetParameter(index,par[index]);
			fyAll->SetParError(index,err[index]);
		}
	}
	//fyAll->SetParameters(par);
	//fyAll->SetParErrors(err);
	//if (stBFitCase.pbToggle[rho]) fyAll->SetParLimits(rho,0,1);
	//fyAll->SetParLimits(gammaT1,-0.2,0.2);
	//fyAll->SetParLimits(r2,-100,1e4);
	//fyAll->SetParLimits(epsV,-1,1);
	//fyAll->SetParLimits(gammaT2,0,0.1);
// Print seed values that are assigned to the fit function
	cout << "PARAMETER SEED VALUES" << endl << separator << endl;
	cout << setw(14) << "Par name" << setw(10) << "Varying?" << "\t" << "Par init val and step" << endl << separator << endl;
	for (index = 0; index < nPars; index++) {
		cout << setw(14) << fyAll->GetParName(index) << setw(10) << tog[index] << "\t" << fyAll->GetParameter(index) << " +/- " << fyAll->GetParError(index) << endl;
	}
	cout << separator << endl;
//	if (!strcmp(stBDNCases[iBDNCaseIndex].pcsCaseCode,"134sb01") ||
//		!strcmp(stBDNCases[iBDNCaseIndex].pcsCaseCode,"134sb02") ||
//		!strcmp(stBDNCases[iBDNCaseIndex].pcsCaseCode,"134sb03") ||
//		!strcmp(stBDNCases[iBDNCaseIndex].pcsCaseCode,"134sb0103"))
//	{
//		b134sbFlag = 1;
//		fyAll->SetParameter(gammaT3, stBFitCase.pdSeed[gammaT2]);
//		par
//		cout << "134-Sb data detected. Forcing gammaT3 = gammaT2. YOU SHOULD GUARANTEE THAT X3 = Y3 = 0. You can set epsX = epsY = 0." << endl << endl;
//	}
	if (b134sbFlag) cout << "134-Sb data detected. Forcing gammaT3 = gammaT2. YOU SHOULD GUARANTEE THAT X3 = Y3 = 0. You can set epsX = epsY = 0." << endl << endl;
	else cout << endl;
// Initialize all functions to parameter seed values
	fyDC	-> SetParameters(par);
	fyT1	-> SetParameters(par);
	fyT2	-> SetParameters(par);
	fyT3	-> SetParameters(par);
	fyU1	-> SetParameters(par);
	fyU2	-> SetParameters(par);
	fyU3	-> SetParameters(par);
	frDC	-> SetParameters(par);
	frT1	-> SetParameters(par);
	frT2	-> SetParameters(par);
	frT3	-> SetParameters(par);
	frU1	-> SetParameters(par);
	frU2	-> SetParameters(par);
	frU3	-> SetParameters(par);
	frAll	-> SetParameters(par);
	foT1	-> SetParameters(par);
	foT2	-> SetParameters(par);
	foT3	-> SetParameters(par);
	foU1	-> SetParameters(par);
	foU2	-> SetParameters(par);
	foU3	-> SetParameters(par);
	
	if (stBFitCase.bDoFit) {
	// Fix the parameters that are supposed to be fixed
		for (index = 0; index < nPars; index++) {
			if (tog[index] == 0) fyAll->FixParameter(index, stBFitCase.pdSeed[index]);
		}
		if (b134sbFlag && tog[gammaT2]==0) fyAll->FixParameter(gammaT3, stBFitCase.pdSeed[gammaT2]);
	// Print initial parameters: (see similar code in yAll() in BFit2Model.cxx)
		if (nParChanges == 0) {
			printf("Ini pars (   0): ",nParChanges);
			for (index = 0; index < nPars; index++) {
				if (stBFitCases[iBFitCaseIndex].pbToggle[index]) printf("%s=%.4e ",parNames[index],par[index]);
			}
			printf("\n");
		}
	//	fyAll->SetParLimits(epsT,0,1);
	//	fyAll->SetParLimits(epsU,0,1);
	//	fyAll->SetParLimits(epsV,0,1);
	//	fyAll->SetParLimits(epsW,0,1);
	//	fyAll->SetParLimits(epsX,0,1);
	//	fyAll->SetParLimits(epsY,0,1);
	//	fyAll->SetParLimits(epsZ,0,1);
	// Do fit and get results
		timer = clock();
	//	TVirtualFitter::SetDefaultFitter("Minuit2");
	//	TVirtualFitter *fitter;
	//	
		TFitResultPtr fit = h1->Fit(fyAll,stBFitCase.pcsOptions);
		timer = clock() - timer;
		printf("\nFitting done in %d clicks (%f seconds).\n", timer, (Float_t)timer/CLOCKS_PER_SEC);
		printf("Fit status = %d\n",fit->Status());
		printf("Cov status = %d\n",fit->CovMatrixStatus());
		TMatrixDSym cov = fit->GetCovarianceMatrix();
	//	Double_t *cov2 = fit->GetCovarianceMatrix();
	//	Int_t ii=0, jj=0;
		Double_t *covArray = new Double_t [nPars*nPars];
		covArray = cov.GetMatrixArray();
		
	//	for (Int_t ii=0; ii<8; ii++) {
	//	//	printf("%.3e  ",cov.GetMatrixArray()[ii]);
	//		for (Int_t jj=0; jj<8; jj++) {
	//		//	printf("%.6e  ", cov[ii][jj]);
	//			printf("%.6e  ", covArray[ii+jj*nPars] );
	//		}
	//		cout << endl;
	//	}
		fit->Print("V");
	// Get parameter values from fit
		for (index = 0; index < nPars; index++) {
			par[index] = fyAll->GetParameter(index);
			err[index] = fyAll->GetParError(index);
		//	printf("%12s err = %f\n", fyAll->GetParName(index), err[index]);
		}
	// par is now up to date...
	// Set other functions to parameter values from fit
		BFitNamespace::ComputeParameterDependentVars(par);
		fyDC	-> SetParameters(par);
		fyT1	-> SetParameters(par);
		fyT2	-> SetParameters(par);
		fyT3	-> SetParameters(par);
		fyU1	-> SetParameters(par);
		fyU2	-> SetParameters(par);
		fyU3	-> SetParameters(par);
		frDC	-> SetParameters(par);
		frT1	-> SetParameters(par);
		frT2	-> SetParameters(par);
		frT3	-> SetParameters(par);
		frU1	-> SetParameters(par);
		frU2	-> SetParameters(par);
		frU3	-> SetParameters(par);
		frAll	-> SetParameters(par);
		foT1	-> SetParameters(par);
		foT2	-> SetParameters(par);
		foT3	-> SetParameters(par);
		foU1	-> SetParameters(par);
		foU2	-> SetParameters(par);
		foU3	-> SetParameters(par);
		
		if (0) {
			for (Int_t iPar = 0; iPar < nPars; iPar++) {
			//	grads[iPar] = frT2->GradientPar(iPar,par,0.0001);
				printf("%12s val = %f, err = %f\n", fyAll->GetParName(iPar), par[iPar], err[iPar]);
			}
			TF1 *fInt	= new TF1("frT2", rT2, 0.0, tCyc, nPars);
			sleep(1.5);
			Double_t grads[nPars];
			printf("\n");
			Double_t timeVar[1] = {0.0};
		//	fInt->GradientPar(timeVar,grads,0.001);
		//	fInt->GradientPar(par,grads,0.001);
			for (Int_t iPar = 0; iPar < nPars; iPar++) {
			//	grads[iPar] = frT2->GradientPar(iPar,par,0.0001);
				printf("%s grad = %.24f\n", parNames[iPar], grads[iPar]);
			}
			printf("\n");
			
			Double_t par1[nPars], par2[nPars];
			memcpy(par1,par,nPars*sizeof(Double_t));
			memcpy(par2,par,nPars*sizeof(Double_t));
			par2[r2] = par2[r2]*1.1;
			fInt->SetParameters(par1);
			Double_t fr_int_1 = fInt->Integral(0.0, tCyc);
			fInt->SetParameters(par2);
			Double_t fr_int_2 = fInt->Integral(0.0, tCyc);
			sleep(0.25);
			printf("int1 = %f, int2 = %f\n", fr_int_1, fr_int_2);
		}
		
	// Estimate # of betas detected and error
		T1_integral = frT1->Integral( 0.0, tCyc);
		T2_integral = frT2->Integral( 0.0, tCyc);
		T3_integral = frT3->Integral( 0.0, tCyc);
		U1_integral = frU1->Integral( 0.0, tCyc);
		U2_integral = frU2->Integral( 0.0, tCyc);
		U3_integral = frU3->Integral( 0.0, tCyc);
		DC_integral = frDC->Integral( 0.0, tCyc);
		All_integral = frAll->Integral( 0.0, tCyc);
		Integral_sum = DC_integral + T1_integral + T2_integral + T3_integral + U1_integral + U2_integral + U3_integral;
		
		timer = clock() - timer;
		printf("\nIntegrals computed in %d clicks (%f seconds).\n", timer, (Float_t)timer/CLOCKS_PER_SEC);
		
		if (stBFitCase.bComputeOtherIntegrals) {
			U1_integral_trap_empty = frU1->Integral( 0.0, tBac);
			U2_integral_trap_empty = frU2->Integral( 0.0, tBac);
			U3_integral_trap_empty = frU3->Integral( 0.0, tBac);
			U1_integral_trap_full  = frU1->Integral( tBac, tCyc);
			U2_integral_trap_full  = frU2->Integral( tBac, tCyc);
			U3_integral_trap_full  = frU3->Integral( tBac, tCyc);
			
			timer = clock() - timer;
			printf("Other integrals computed in %d clicks (%f seconds).\n", timer, (Float_t)timer/CLOCKS_PER_SEC);
		}
		
//		T1_integral_error = frT1->IntegralError( 0.0, tCyc, par, cov.GetMatrixArray() );
	//	printf("Getting rT2 integral error in 5 seconds...\n"); sleep(1);
	//	printf("Getting rT2 integral error in 4 seconds...\n"); sleep(1);
	//	printf("Getting rT2 integral error in 3 seconds...\n"); sleep(1);
	//	printf("Getting rT2 integral error in 2 seconds...\n"); sleep(1);
	//	printf("Getting rT2 integral error in 1 seconds...\n"); sleep(1);
		
	//	T2_integral_error = frT2->IntegralError( 0.0, tCyc, par, cov.GetMatrixArray() );
	//	T1_integral_error = intErr( frT1, covArray, 0.0, tCyc );
		
		if (stBFitCase.bUseMyErrorFunction) {
			T1_integral_error = intErr( frT1, covArray, 0.0, tCyc );
			T2_integral_error = intErr( frT2, covArray, 0.0, tCyc );
			T3_integral_error = intErr( frT3, covArray, 0.0, tCyc );
			U1_integral_error = intErr( frU1, covArray, 0.0, tCyc );
			U2_integral_error = intErr( frU2, covArray, 0.0, tCyc );
			U3_integral_error = intErr( frU3, covArray, 0.0, tCyc );
			DC_integral_error = intErr( frDC, covArray, 0.0, tCyc );
			All_integral_error = intErr( frAll, covArray, 0.0, tCyc );
		}
		else {
			T1_integral_error = frT1->IntegralError( 0.0, tCyc, par, cov.GetMatrixArray() );
			T2_integral_error = frT2->IntegralError( 0.0, tCyc, par, cov.GetMatrixArray() );
			T3_integral_error = frT3->IntegralError( 0.0, tCyc, par, cov.GetMatrixArray() );
			U1_integral_error = frU1->IntegralError( 0.0, tCyc, par, cov.GetMatrixArray() );
			U2_integral_error = frU2->IntegralError( 0.0, tCyc, par, cov.GetMatrixArray() );
			U3_integral_error = frU3->IntegralError( 0.0, tCyc, par, cov.GetMatrixArray() );
			All_integral_error = frAll->IntegralError( 0.0, tCyc, par, cov.GetMatrixArray() );
		}
//		T3_integral_error = frT3->IntegralError( 0.0, tCyc, par, cov.GetMatrixArray() );
//		U1_integral_error = frU1->IntegralError( 0.0, tCyc, par, cov.GetMatrixArray() );
//		U2_integral_error = frU2->IntegralError( 0.0, tCyc, par, cov.GetMatrixArray() );
//		U3_integral_error = frU3->IntegralError( 0.0, tCyc, par, cov.GetMatrixArray() );
//		DC_integral_error = frDC->IntegralError( 0.0, tCyc, par, cov.GetMatrixArray() );
//		All_integral_error = frAll->IntegralError( 0.0, tCyc, par, cov.GetMatrixArray() );
		
	//	printf("tCyc=%f\n",tCyc);
	//	T1_integral_error = frT1->IntegralError( 0.0, 1000*tCyc);
	//	T2_integral_error = frT2->IntegralError( 0.0, tCyc);
	//	T3_integral_error = frT3->IntegralError( 0.0, tCyc);
	//	U1_integral_error = frU1->IntegralError( 0.0, tCyc);
	//	U2_integral_error = frU2->IntegralError( 0.0, tCyc);
	//	U3_integral_error = frU3->IntegralError( 0.0, tCyc);
	//	DC_integral_error = frDC->IntegralError( 0.0, tCyc);
	//	All_integral_error = frAll->IntegralError( 0.0, tCyc);
		
		Integral_sum_error = Sqrt( Power(DC_integral_error,2.0) + Power(T1_integral_error,2.0) + Power(T2_integral_error,2.0) + Power(T3_integral_error,2.0) + Power(U1_integral_error,2.0) + Power(U2_integral_error,2.0) + Power(U3_integral_error,2.0) );
		
		timer = clock() - timer;
		printf("Errors computed in %d clicks (%f seconds).\n", timer, (Float_t)timer/CLOCKS_PER_SEC);
		
		cout << endl << separator << endl;
		printf("NUMBER OF BETAS DETECTED, by population:\n");
		cout << separator << endl;
		printf("T1 integral = %.1f +/- %.1f\n", T1_integral, T1_integral_error);
		printf("U1 integral = %.1f +/- %.1f\n", U1_integral, U1_integral_error);
		printf("T2 integral = %.1f +/- %.1f\n", T2_integral, T2_integral_error);
		printf("U2 integral = %.1f +/- %.1f\n", U2_integral, U2_integral_error);
		printf("T3 integral = %.1f +/- %.1f\n", T3_integral, T3_integral_error);
		printf("U3 integral = %.1f +/- %.1f\n", U3_integral, U3_integral_error);
		printf("DC integral = %.1f +/- %.1f\n", DC_integral, DC_integral_error);
		cout << separator << endl;
		printf("Sum of above = %.1f +/- %.1f <-- no cov in unc\n", Integral_sum, Integral_sum_error);
		printf("All integral = %.1f +/- %.1f\n", All_integral, All_integral_error);
		if (stBFitCase.bComputeOtherIntegrals) {
			cout << separator << endl;
			printf("U1 with trap emtpy = %.1f; trap full = %.1f\n", U1_integral_trap_empty, U1_integral_trap_full);
			printf("U2 with trap emtpy = %.1f; trap full = %.1f\n", U2_integral_trap_empty, U2_integral_trap_full);
			printf("U3 with trap emtpy = %.1f; trap full = %.1f\n", U3_integral_trap_empty, U3_integral_trap_full);
//			printf("All untrapped with trap empty = %f (bin %f to bin %f)\n", h1->Integral(binZero, binCap-1) - 0.001*par[DC]*tBac, binZero, binCap);
//			printf("All data area in histogram = %f\n", h1->Integral(binZero, binCycle-1));
//			printf("All data area in histogram = %f\n", h1->Integral());
		}
		cout << separator << endl << endl;
		
	} // end if (stBFitCase.bDoFit)
	
// Draw
	Double_t yMin, yMax, yRange;
	
	fyAll->SetLineColor(kBlack);
	fyDC->SetLineColor(kBlack);
	foT1->SetLineColor(kGreen);
	foT2->SetLineColor(kBlue);
	foT3->SetLineColor(kRed);
	foU1->SetLineColor(kGreen);
	foU2->SetLineColor(kBlue);
	foU3->SetLineColor(kRed);
	h1->SetLineColor(16);
//	h2->SetLineColor(16);
	
	fyDC->SetLineStyle(2);
	foU1->SetLineStyle(2);
	foU2->SetLineStyle(2);
	foU3->SetLineStyle(2);
	
//	Double_t nPoints = 2 * (1000.0 * stBDNCase.dCycleTime / dRebinFactor); // leading interger is number of points per bin
// Now defined above
	fyAll->SetNpx(nPoints);
	foT1->SetNpx(nPoints);
	foT2->SetNpx(nPoints);
	foT3->SetNpx(nPoints);
	foU1->SetNpx(nPoints);
	foU2->SetNpx(nPoints);
	foU3->SetNpx(nPoints);
	
	FuncPrep(fyAll,par,nPoints,kBlack,1);
	FuncPrep(foT1,par,nPoints,kGreen,1);
	FuncPrep(foT2,par,nPoints,kBlue,1);
	FuncPrep(foT3,par,nPoints,kRed,1);
	FuncPrep(foU1,par,nPoints,kGreen,7);
	FuncPrep(foU2,par,nPoints,kBlue,7);
	FuncPrep(foU3,par,nPoints,kRed,7);
//	Float_t wdth = 2;
//	fyAll->SetLineWidth(wdth);
//	foT1->SetLineWidth(wdth);
//	foT2->SetLineWidth(1.5*wdth);
//	foT3->SetLineWidth(wdth);
//	foU1->SetLineWidth(wdth);
//	foU2->SetLineWidth(wdth);
//	foU3->SetLineWidth(wdth);
	
	Int_t bsmPlot = 0;
	Float_t scale = 1.0;
	//if (bsmPlot) scale = 2.0;
	TCanvas *c_BFit	= new TCanvas("c_BFit","Beta singles fit",scale*945,scale*600);
	gStyle->SetOptStat("e");
	gStyle->SetOptFit(1111);
	h1 ->Draw("HIST");
	h1->SetXTitle("Cycle time (ms)");
	char yTitle[100];
	sprintf(yTitle, "Detections / %d ms", Nint(dRebinFactor));
	h1->SetYTitle(yTitle);
	h1->GetYaxis()->SetTitleOffset(1.4);
	fyAll ->Draw("SAME");
	foU1->Draw("SAME");
	foT1->Draw("SAME");
	foU3->Draw("SAME");
	foT3->Draw("SAME");
	foU2->Draw("SAME");
	foT2->Draw("SAME");
	fyDC->Draw("SAME");
	//c_B_fit->Update();
	
	//Double_t xMin = h1->GetXaxis()->GetXmin();
	//Double_t xMax = h1->GetXaxis()->GetXmax();
	yMin = par[nCyc]*(0.001*par[DC])*par[dt];
	yMax = h1->GetMaximum();
//	printf("\nY Range is %f to %f\n\n",yMin,yMax);
	yRange = yMax - yMin;
	yMin = yMin - 0.05*yRange;
	yMax = yMax + 0.30*yRange;
//	yMin = 0.0;
//	yMax = 14000.0;
	//printf("Range = (%f,%f)\n", yMin, yMax);
	h1->GetYaxis()->SetRangeUser(yMin,yMax);
//	h2->GetYaxis()->SetRangeUser(yMin,yMax);
	h1->GetXaxis()->SetRangeUser(-1000,tCyc+1000);
//	h2->GetXaxis()->SetRangeUser(-1000,tCyc+1000);
	
	if (bsmPlot) { // special case for a plot in Beta Singles Model appendix in my thesis
		h1->SetXTitle("Cycle time (ms)");
		h1->SetYTitle("Detections / 500 ms");
		h1->SetTitle("Simulation vs Model, All populations");
		h1->GetXaxis()->SetRangeUser(0,tCyc);
		h1->GetYaxis()->SetTitleOffset(1.0);
		TLegend *leg_1 = new TLegend(0.15, 0.62, 0.40, 0.87);
		leg_1->AddEntry(h1 , "Data");
		leg_1->AddEntry(fyDC, "DC");
		char pcsT1Name[STRING_SIZE], pcsT2Name[STRING_SIZE], pcsT3Name[STRING_SIZE];
		char pcsU1Name[STRING_SIZE], pcsU2Name[STRING_SIZE], pcsU3Name[STRING_SIZE];
		sprintf(pcsT1Name,"Species 1, trapped");
		sprintf(pcsT2Name,"Species 2, trapped");
		sprintf(pcsT3Name,"Species 3, trapped");
		sprintf(pcsU1Name,"Species 1, untrapped");
		sprintf(pcsU2Name,"Species 2, untrapped");
		sprintf(pcsU3Name,"Species 3, untrapped");
		leg_1->AddEntry(foT1, pcsT1Name);
		leg_1->AddEntry(foU1, pcsU1Name);
		leg_1->AddEntry(foT2, pcsT2Name);
		leg_1->AddEntry(foU2, pcsU2Name);
		leg_1->AddEntry(foT3, pcsT3Name);
		leg_1->AddEntry(foU3, pcsU3Name);
		leg_1->AddEntry(fyAll , "Fit function (all species)");
		leg_1->SetFillColor(0);
		leg_1->Draw();
	}
	else { // normal case
		TLegend *leg_1 = new TLegend(0.13, 0.69, 0.35, 0.94);
		leg_1->AddEntry(h1 , "Data");
		leg_1->AddEntry(fyDC, "DC");
		char pcsT1Name[STRING_SIZE], pcsT2Name[STRING_SIZE], pcsT3Name[STRING_SIZE];
		char pcsU1Name[STRING_SIZE], pcsU2Name[STRING_SIZE], pcsU3Name[STRING_SIZE];
		strcpy(pcsT1Name, stBDNCase.pcsSpecies1Name); strcat(pcsT1Name, " trapped");
		strcpy(pcsT2Name, stBDNCase.pcsSpecies2Name); strcat(pcsT2Name, " trapped");
		strcpy(pcsT3Name, stBDNCase.pcsSpecies3Name); strcat(pcsT3Name, " trapped");
		strcpy(pcsU1Name, stBDNCase.pcsSpecies1Name); strcat(pcsU1Name, " untrapped");
		strcpy(pcsU2Name, stBDNCase.pcsSpecies2Name); strcat(pcsU2Name, " untrapped");
		strcpy(pcsU3Name, stBDNCase.pcsSpecies3Name); strcat(pcsU3Name, " untrapped");
		leg_1->AddEntry(foT1, pcsT1Name);
		leg_1->AddEntry(foU1, pcsU1Name);
		leg_1->AddEntry(foT2, pcsT2Name);
		leg_1->AddEntry(foU2, pcsU2Name);
		leg_1->AddEntry(foT3, pcsT3Name);
		leg_1->AddEntry(foU3, pcsU3Name);
		leg_1->AddEntry(fyAll , "Fit function (all species)");
		leg_1->SetFillColor(0);
		leg_1->Draw();
	}
	
	TPaveText *pt = new TPaveText(.70, .85, .91, .94, "NDC");
	pt->AddText(stBFitCase.pcsCaseCode);
	pt->SetTextFont(82);
	pt->Draw();
	
	if (stBFitCase.bDoFit) {
		char int_DC[100];
		char int_T1[100], int_T2[100], int_T3[100];
		char int_U1[100], int_U2[100], int_U3[100];
		char int_All[100];
//		sprintf(int1,"T1 = %8.0f +/- %6.0f; U1 = %8.0f +/- %6.0f", T1_integral, T1_integral_error, U1_integral, U1_integral_error);
//		sprintf(int2,"T2 = %8.0f +/- %6.0f; U2 = %8.0f +/- %6.0f", T2_integral, T2_integral_error, U2_integral, U2_integral_error);
//		sprintf(int3,"T3 = %8.0f +/- %6.0f; U3 = %8.0f +/- %6.0f", T3_integral, T3_integral_error, U3_integral, U3_integral_error);
		sprintf(int_DC,  "A=%9.0f #pm%9.0f", DC_integral, DC_integral_error);
		sprintf(int_T1,  "A=%9.0f #pm%9.0f", T1_integral, T1_integral_error);
		sprintf(int_U1,  "A=%9.0f #pm%9.0f", U1_integral, U1_integral_error);
		sprintf(int_T2,  "A=%9.0f #pm%9.0f", T2_integral, T2_integral_error);
		sprintf(int_U2,  "A=%9.0f #pm%9.0f", U2_integral, U2_integral_error);
		sprintf(int_T3,  "A=%9.0f #pm%9.0f", T3_integral, T3_integral_error);
		sprintf(int_U3,  "A=%9.0f #pm%9.0f", U3_integral, U3_integral_error);
		sprintf(int_All, "A=%9.0f #pm%9.0f", All_integral, All_integral_error);
		TLegend *calc = new TLegend(.13, 0.50, 0.35, 0.69);
		calc->AddEntry(fyDC, int_DC);
		calc->AddEntry(foT1, int_T1);
		calc->AddEntry(foU1, int_U1);
		calc->AddEntry(foT2, int_T2);
		calc->AddEntry(foU2, int_U2);
		calc->AddEntry(foT3, int_T3);
		calc->AddEntry(foU3, int_U3);
		calc->AddEntry(fyAll, int_All);
		calc->SetTextFont(82);
		calc->SetTextSize(0.02);
		calc->SetFillColor(0);
		calc->Draw();
	} // end if (stBFitCase.bDoFit)
	
	gPad->Update();
	TPaveStats *stats_1 = (TPaveStats*)h1->FindObject("stats");
	stats_1->SetTextFont(82);
	//stats_1->SetTextSize(0.02);
	stats_1->SetX1NDC(.35);
	stats_1->SetX2NDC(.60);
	stats_1->SetY1NDC(.69);
	stats_1->SetY2NDC(.94);
	gPad->Update();
	c_BFit->Modified();
	outfile->WriteTObject(c_BFit);
	outfile->WriteTObject(h1);
	
	if (stBFitCase.bMonteCarlo) {
		
		TH1D *h_DDC	= (TH1D*)f->Get("h_DDC_cyctime");
		TH1D *h_DT1	= (TH1D*)f->Get("h_DT1_cyctime");
		TH1D *h_DT2	= (TH1D*)f->Get("h_DT2_cyctime");
		TH1D *h_DT3	= (TH1D*)f->Get("h_DT3_cyctime");
		TH1D *h_DU1	= (TH1D*)f->Get("h_DU1_cyctime");
		TH1D *h_DU2	= (TH1D*)f->Get("h_DU2_cyctime");
		TH1D *h_DU3	= (TH1D*)f->Get("h_DU3_cyctime");
		
		h2->SetLineColor(kBlack);
		if (stBFitCase.bHasDDC) h_DDC->SetLineColor(kBlack);
//		h_DT1->SetLineColor(kGreen);
//		h_DU1->SetLineColor(kGreen);
//		h_DT2->SetLineColor(kBlue);
//		h_DU2->SetLineColor(kBlue);
//		h_DT3->SetLineColor(kRed);
//		h_DU3->SetLineColor(kRed);
		
		if (stBFitCase.bHasDDC) h_DDC->SetLineStyle(7);
//		h_DU1->SetLineStyle(7);
//		h_DU2->SetLineStyle(7);
//		h_DU3->SetLineStyle(7);
		
		if (stBFitCase.bHasDDC) h_DDC->Rebin(dRebinFactor);
	//	h_DT1->Rebin(dRebinFactor);
	//	h_DT2->Rebin(dRebinFactor);
	//	h_DT3->Rebin(dRebinFactor);
	//	h_DU1->Rebin(dRebinFactor);
	//	h_DU2->Rebin(dRebinFactor);
	//	h_DU3->Rebin(dRebinFactor);
	//	
	//	fyT1->SetLineColor(kGreen);
	//	fyT2->SetLineColor(kBlue);
	//	fyT3->SetLineColor(kRed);
		fyU1->SetLineColor(kGreen);
		fyU2->SetLineColor(kBlue);
		fyU3->SetLineColor(kRed);
		
//		fyU1->SetLineStyle(7);
//		fyU2->SetLineStyle(7);
//		fyU3->SetLineStyle(7);
		
	//	fyT1->SetNpx(nPoints);
	//	fyT2->SetNpx(nPoints);
	//	fyT3->SetNpx(nPoints);
		fyU1->SetNpx(nPoints);
		fyU2->SetNpx(nPoints);
		fyU3->SetNpx(nPoints);
		
		h2->GetYaxis()->SetRangeUser(0,yMax);
		h2->GetXaxis()->SetRangeUser(-1000,tCyc+1000);
		
//		TCanvas *c_decays_cyctime = new TCanvas("c_decays_cyctime","Decays versus cycle time",945,600);
//		h2->Draw();
//		fyAll->Draw("SAME");
//		if (stBFitCase.bHasDDC) fyDC->Draw("SAME");
//		fyT1->Draw("SAME");
//		fyT2->Draw("SAME");
//		fyT3->Draw("SAME");
//		fyU1->Draw("SAME");
//		fyU2->Draw("SAME");
//		fyU3->Draw("SAME");
//		if (stBFitCase.bHasDDC) h_DDC->Draw("SAME");
//		h_DT1->Draw("SAME");
//		h_DT2->Draw("SAME");
//		h_DT3->Draw("SAME");
//		h_DU1->Draw("SAME");
//		h_DU2->Draw("SAME");
//		h_DU3->Draw("SAME");
//		gPad->Update();
//		TPaveStats *stats_2 = (TPaveStats*)h2->FindObject("stats");
//		stats_2->SetX1NDC(.13);
//		stats_2->SetX2NDC(.32);
//		stats_2->SetY1NDC(.80);
//		stats_2->SetY2NDC(.88);
//		gPad->Update();
//		c_decays_cyctime->Modified();
//		outfile->WriteTObject(c_decays_cyctime);
		
		printf("T1 entries:    %10d\n",(Int_t)h_DT1->GetEntries());
		printf("U1 entries:    %10d\n",(Int_t)h_DU1->GetEntries());
		printf("T2 entries:    %10d\n",(Int_t)h_DT2->GetEntries());
		printf("U2 entries:    %10d\n",(Int_t)h_DU2->GetEntries());
		printf("T3 entries:    %10d\n",(Int_t)h_DT3->GetEntries());
		printf("U3 entries:    %10d\n",(Int_t)h_DU3->GetEntries());
		if (stBFitCase.bHasDDC) printf("DC entries:    %10d\n",(Int_t)h_DDC->GetEntries());
		printf("Total entries: %10d\n",(Int_t)h1->GetEntries());
		
		if (stBFitCase.bHasVWXY) { // used for Monte Carlo data where V, W, X, Y pops are known
			
			TH1D *h_DV1	= (TH1D*)f->Get("h_DV1_cyctime");
			TH1D *h_DV2	= (TH1D*)f->Get("h_DV2_cyctime");
			TH1D *h_DV3	= (TH1D*)f->Get("h_DV3_cyctime");
			TH1D *h_DW1	= (TH1D*)f->Get("h_DW1_cyctime");
			TH1D *h_DW2	= (TH1D*)f->Get("h_DW2_cyctime");
			TH1D *h_DW3	= (TH1D*)f->Get("h_DW3_cyctime");
			TH1D *h_DZ1	= (TH1D*)f->Get("h_DZ1_cyctime");
			TH1D *h_DZ2	= (TH1D*)f->Get("h_DZ2_cyctime");
			TH1D *h_DZ3	= (TH1D*)f->Get("h_DZ3_cyctime");
			TH1D *h_DX2	= (TH1D*)f->Get("h_DX2_cyctime");
			TH1D *h_DX3	= (TH1D*)f->Get("h_DX3_cyctime");
			TH1D *h_DY2	= (TH1D*)f->Get("h_DY2_cyctime");
			TH1D *h_DY3	= (TH1D*)f->Get("h_DY3_cyctime");
			
			HistPrep(h_DT1,dRebinFactor,dBinWidth,"T");
			HistPrep(h_DT2,dRebinFactor,dBinWidth,"T");
			HistPrep(h_DT3,dRebinFactor,dBinWidth,"T");
			HistPrep(h_DV1,dRebinFactor,dBinWidth,"V");
			HistPrep(h_DV2,dRebinFactor,dBinWidth,"V");
			HistPrep(h_DV3,dRebinFactor,dBinWidth,"V");
			HistPrep(h_DW1,dRebinFactor,dBinWidth,"W");
			HistPrep(h_DW2,dRebinFactor,dBinWidth,"W");
			HistPrep(h_DW3,dRebinFactor,dBinWidth,"W");
			HistPrep(h_DZ1,dRebinFactor,dBinWidth,"Z");
			HistPrep(h_DZ2,dRebinFactor,dBinWidth,"Z");
			HistPrep(h_DZ3,dRebinFactor,dBinWidth,"Z");
			HistPrep(h_DX2,dRebinFactor,dBinWidth,"X");
			HistPrep(h_DX3,dRebinFactor,dBinWidth,"X");
			HistPrep(h_DY2,dRebinFactor,dBinWidth,"Y");
			HistPrep(h_DY3,dRebinFactor,dBinWidth,"Y");
			
/*			h_DV1 ->Rebin(dRebinFactor);
			h_DV2 ->Rebin(dRebinFactor);
			h_DV3 ->Rebin(dRebinFactor);
			h_DW1 ->Rebin(dRebinFactor);
			h_DW2 ->Rebin(dRebinFactor);
			h_DW3 ->Rebin(dRebinFactor);
			h_DZ1 ->Rebin(dRebinFactor);
			h_DZ2 ->Rebin(dRebinFactor);
			h_DZ3 ->Rebin(dRebinFactor);
			h_DX2 ->Rebin(dRebinFactor);
			h_DX3 ->Rebin(dRebinFactor);
			h_DY2 ->Rebin(dRebinFactor);
			h_DY3 ->Rebin(dRebinFactor);
			
			h_DV1->SetLineColor(kBlack);
			h_DV2->SetLineColor(kBlack);
			h_DV3->SetLineColor(kBlack);
			h_DX2->SetLineColor(kBlack);
			h_DY2->SetLineColor(kBlack);
			h_DX3->SetLineColor(kBlack);
			h_DY3->SetLineColor(kBlack);
*/			
//			h_DX2->SetLineStyle(8);
//			h_DY2->SetLineStyle(3);
//			h_DX3->SetLineStyle(8);
//			h_DY3->SetLineStyle(3);
			
			TH1D *h_DU2_1	= (TH1D*)h_DU2->Rebin(1,"h_DU2_1");
			TH1D *h_DU3_1	= (TH1D*)h_DU3->Rebin(1,"h_DU3_1");
			
			h_DU2_1->SetLineStyle(1);
			h_DU3_1->SetLineStyle(1);
			h_DU2_1->SetLineColor(kBlack);
			h_DU3_1->SetLineColor(kBlack);
			
			TF1 *fyV1	= new TF1("fyV1", yV1, 0.0, tCyc, nPars);
			TF1 *fyV2	= new TF1("fyV2", yV2, 0.0, tCyc, nPars);
			TF1 *fyV3	= new TF1("fyV3", yV3, 0.0, tCyc, nPars);
			TF1 *fyW1	= new TF1("fyW1", yW1, 0.0, tCyc, nPars);
			TF1 *fyW2	= new TF1("fyW2", yW2, 0.0, tCyc, nPars);
			TF1 *fyW3	= new TF1("fyW3", yW3, 0.0, tCyc, nPars);
			TF1 *fyZ1	= new TF1("fyZ1", yZ1, 0.0, tCyc, nPars);
			TF1 *fyZ2	= new TF1("fyZ2", yZ2, 0.0, tCyc, nPars);
			TF1 *fyZ3	= new TF1("fyZ3", yZ3, 0.0, tCyc, nPars);
			TF1 *fyX2	= new TF1("fyX2", yX2, 0.0, tCyc, nPars);
			TF1 *fyX3	= new TF1("fyX3", yX3, 0.0, tCyc, nPars);
			TF1 *fyY2	= new TF1("fyY2", yY2, 0.0, tCyc, nPars);
			TF1 *fyY3	= new TF1("fyY3", yY3, 0.0, tCyc, nPars);
			
			FuncPrep(fyT1,par,nPoints,kGreen,2);//9);
			FuncPrep(fyT2,par,nPoints,kBlue,2);//9);
			FuncPrep(fyT3,par,nPoints,kRed,2);//9);
			FuncPrep(fyV1,par,nPoints,kGreen,2);//9);
			FuncPrep(fyV2,par,nPoints,kBlue,2);//9);
			FuncPrep(fyV3,par,nPoints,kRed,2);//9);
			FuncPrep(fyW1,par,nPoints,kGreen,2);//2);
			FuncPrep(fyW2,par,nPoints,kBlue,2);//2);
			FuncPrep(fyW3,par,nPoints,kRed,2);//2);
			FuncPrep(fyZ1,par,nPoints,kGreen,2);//8);
			FuncPrep(fyZ2,par,nPoints,kBlue,2);//8);
			FuncPrep(fyZ3,par,nPoints,kRed,2);//8);
			FuncPrep(fyX2,par,nPoints,kBlue,2);//5);
			FuncPrep(fyX3,par,nPoints,kRed,2);//5);
			FuncPrep(fyY2,par,nPoints,kBlue,2);
			FuncPrep(fyY3,par,nPoints,kRed,2);
/*			fyV1->SetParameters(par);
			fyV2->SetParameters(par);
			fyV3->SetParameters(par);
			fyW1->SetParameters(par);
			fyW2->SetParameters(par);
			fyW3->SetParameters(par);
			fyZ1->SetParameters(par);
			fyZ2->SetParameters(par);
			fyZ3->SetParameters(par);
			fyX2->SetParameters(par);
			fyX3->SetParameters(par);
			fyY2->SetParameters(par);
			fyY3->SetParameters(par);
			
			fyV1->SetNpx(nPoints);
			fyV2->SetNpx(nPoints);
			fyV3->SetNpx(nPoints);
			fyW1->SetNpx(nPoints);
			fyW2->SetNpx(nPoints);
			fyW3->SetNpx(nPoints);
			fyZ1->SetNpx(nPoints);
			fyZ2->SetNpx(nPoints);
			fyZ3->SetNpx(nPoints);
			fyX2->SetNpx(nPoints);
			fyX3->SetNpx(nPoints);
			fyY2->SetNpx(nPoints);
			fyY3->SetNpx(nPoints);
			
			fyV1->SetLineColor(kGreen);
			fyV2->SetLineColor(kBlue);
			fyV3->SetLineColor(kRed);
			fyW1->SetLineColor(kGreen);
			fyW2->SetLineColor(kBlue);
			fyW3->SetLineColor(kRed);
			fyZ1->SetLineColor(kGreen);
			fyZ2->SetLineColor(kBlue);
			fyZ3->SetLineColor(kRed);
			fyX2->SetLineColor(kBlue);
			fyX3->SetLineColor(kRed);
			fyY2->SetLineColor(kBlue);
			fyY3->SetLineColor(kRed);
*/			
			
			TCanvas *c_feeding = new TCanvas("c_feeding","Decays from pops X and Y versus cycle time",945,600);
			h1->Draw();
//			h_DU3_1->Draw("SAME");
//			fyU2->Draw("SAME");
//			fyU3->Draw("SAME");
			h_DV3->Draw("SAME");
			fyV3->Draw("SAME");
			h_DV1->Draw("SAME");
			h_DV2->Draw("SAME");
			fyV1->Draw("SAME");
			fyV2->Draw("SAME");
//			h_DX2->Draw("SAME");
//			h_DX3->Draw("SAME");
//			fyX2->Draw("SAME");
//			fyX3->Draw("SAME");
//			h_DY2->Draw("SAME");
//			h_DY3->Draw("SAME");
//			fyY2->Draw("SAME");
//			fyY3->Draw("SAME");
//			h_DU2_1->Draw("SAME");
			outfile->WriteTObject(c_feeding);
			
			TCanvas *c_Ti = new TCanvas("c_Ti","Decays from T pops versus cycle time",945,600);
			h_DT2->Draw();
			h_DT1->Draw("SAME");
			h_DT3->Draw("SAME");
			fyT2->Draw("SAME");
			fyT3->Draw("SAME");
			fyT1->Draw("SAME");
			yMax = Max(Max(h_DT1->GetMaximum(),h_DT2->GetMaximum()),h_DT3->GetMaximum());
			yMin = 0.0;//-0.05 * yMax;
			yMax =  1.05 * yMax;
			h_DT2->GetYaxis()->SetRangeUser(yMin,yMax);
			outfile->WriteTObject(c_Ti);
			
			TCanvas *c_Vi = new TCanvas("c_Vi","Decays from V pops versus cycle time",945,600);
			h_DV2->Draw();
			h_DV1->Draw("SAME");
			h_DV3->Draw("SAME");
			fyV2->Draw("SAME");
			fyV3->Draw("SAME");
			fyV1->Draw("SAME");
			yMax = Max(Max(h_DV1->GetMaximum(),h_DV2->GetMaximum()),h_DV3->GetMaximum());
			yMin = -20.0;//-0.05 * yMax;
			yMax =  1.05 * yMax;
			h_DV2->GetYaxis()->SetRangeUser(yMin,yMax);
			outfile->WriteTObject(c_Vi);
			
			TCanvas *c_Wi = new TCanvas("c_Wi","Decays from W pops versus cycle time",945,600);
			h_DW2->Draw();
			h_DW1->Draw("SAME");
			h_DW3->Draw("SAME");
			fyW2->Draw("SAME");
			fyW3->Draw("SAME");
			fyW1->Draw("SAME");
			yMax = Max(Max(h_DW1->GetMaximum(),h_DW2->GetMaximum()),h_DW3->GetMaximum());
			yMin = -10.0;//-0.05 * yMax;
			yMax =  1.05 * yMax;
			h_DW2->GetYaxis()->SetRangeUser(yMin,yMax);
			outfile->WriteTObject(c_Wi);
			
			TCanvas *c_Zi = new TCanvas("c_Zi","Decays from Z pops versus cycle time",945,600);
			h_DZ2->Draw();
			h_DZ1->Draw("SAME");
			h_DZ3->Draw("SAME");
			fyZ2->Draw("SAME");
			fyZ3->Draw("SAME");
			fyZ1->Draw("SAME");
			yMax = Max(Max(h_DZ1->GetMaximum(),h_DZ2->GetMaximum()),h_DZ3->GetMaximum());
			yMin = -10.0;//-0.05 * yMax;
			yMax =  1.05 * yMax;
			h_DZ2->GetYaxis()->SetRangeUser(yMin,yMax);
			outfile->WriteTObject(c_Zi);
			
			TCanvas *c_Xi = new TCanvas("c_Xi","Decays from X pops versus cycle time",945,600);
			h_DX2->Draw();
			h_DX3->Draw("SAME");
			fyX3->Draw("SAME");
			fyX2->Draw("SAME");
			yMax = Max(h_DX2->GetMaximum(),h_DX3->GetMaximum());
			yMin = 0.0;//-0.05 * yMax;
			yMax =  1.05 * yMax;
			h_DX2->GetYaxis()->SetRangeUser(yMin,yMax);
			outfile->WriteTObject(c_Xi);
			
			TCanvas *c_Yi = new TCanvas("c_Yi","Decays from Y pops versus cycle time",945,600);
			h_DY2->Draw();
			h_DY3->Draw("SAME");
			fyY3->Draw("SAME");
			fyY2->Draw("SAME");
			yMax = Max(h_DY2->GetMaximum(),h_DY3->GetMaximum());
			yMax = Max(yMax,fyY2->GetMaximum(0.0, tCyc, 0.1, 20));
			yMax = Max(yMax,fyY3->GetMaximum(0.0, tCyc, 0.1, 20));
			yMin = 0.0;//-0.05 * yMax;
			yMax =  1.05 * yMax;
			h_DY2->GetYaxis()->SetRangeUser(yMin,yMax);
			outfile->WriteTObject(c_Yi);
			
			TCanvas *c_Ui = new TCanvas("c_Ui","Decays from U pops versus cycle time",945,600);
			h_DU2->Draw();
			h_DU1->Draw("SAME");
			h_DU3->Draw("SAME");
			fyU2->Draw("SAME");
			fyU3->Draw("SAME");
			fyU1->Draw("SAME");
			yMax = Max(Max(h_DU1->GetMaximum(),h_DU2->GetMaximum()),h_DU3->GetMaximum());
			yMin = -0.05 * yMax;
			yMax =  1.05 * yMax;
			h_DU2->GetYaxis()->SetRangeUser(yMin,yMax);
			outfile->WriteTObject(c_Ui);
			
//			TF1 *fyT1copy;//	= new TF1("fyT1F", yT1, 0.0, tCyc, nPars);
//			fyT1->Copy(fyTF1copy)
//			FuncPrep(fyT1copy,par,nPoints,kBlack);
			
			
			TCanvas *c_T1V1W1Z1X2 = new TCanvas("c_T1V1W1Z1X2","Feeding from T1",945,600);
			h_DT1->Draw();
			h_DV1->Draw("SAME");
			h_DW1->Draw("SAME");
			h_DZ1->Draw("SAME");
			h_DX2->Draw("SAME");
//			fyT1->SetLineAttributes(kBlack,2,1);
			fyT1->Draw("SAME");
			yMax = 0.0;
			yMax = Max(yMax,h_DT1->GetMaximum());
			yMax = Max(yMax,h_DW1->GetMaximum());
			yMin = -0.05 * yMax;
			yMax =  1.05 * yMax;
			h_DT1->GetYaxis()->SetRangeUser(yMin,yMax);
			outfile->WriteTObject(c_T1V1W1Z1X2);
			
			TCanvas *c_T2V2W2Z2X3Y3 = new TCanvas("c_T2V2W2Z2X3Y3","Feeding from T2",945,600);
			h_DT2->Draw();
			h_DV2->Draw("SAME");
			h_DW2->Draw("SAME");
			h_DZ2->Draw("SAME");
			h_DX3->Draw("SAME");
			h_DY3->Draw("SAME");
			fyT2->Draw("SAME");
			fyV2->Draw("SAME");
			fyW2->Draw("SAME");
			fyZ2->Draw("SAME");
			fyX3->Draw("SAME");
			fyY3->Draw("SAME");
			yMax = 0.0;
			yMax = Max(yMax,h_DT2->GetMaximum());
			yMax = Max(yMax,h_DV2->GetMaximum());
			yMax = Max(yMax,h_DW2->GetMaximum());
			yMax = Max(yMax,h_DZ2->GetMaximum());
			yMax = Max(yMax,h_DX3->GetMaximum());
			yMax = Max(yMax,h_DY3->GetMaximum());
			yMax = Max(yMax,fyT2->GetMaximum(tBac,tCyc,10.0,20));
			yMax = Max(yMax,fyV2->GetMaximum(tBac,tCyc,10.0,20));
			yMax = Max(yMax,fyW2->GetMaximum(tBac,tCyc,10.0,20));
			yMax = Max(yMax,fyZ2->GetMaximum(tBac,tCyc,10.0,20));
			yMax = Max(yMax,fyX3->GetMaximum(tBac,tCyc,10.0,20));
			yMax = Max(yMax,fyY3->GetMaximum(tBac,tCyc,10.0,20));
			yMin = -0.05 * yMax;
			yMax =  1.05 * yMax;
			h_DT2->GetYaxis()->SetRangeUser(yMin,yMax);
			outfile->WriteTObject(c_T2V2W2Z2X3Y3);
			
			TCanvas *c_T3V3W3Z3 = new TCanvas("c_T3V3W3Z3","Feeding from T3",945,600);
			h_DT3->Draw();
			h_DV3->Draw("SAME");
			h_DW3->Draw("SAME");
			h_DZ3->Draw("SAME");
			fyT3->Draw("SAME");
			fyV3->Draw("SAME");
			fyW3->Draw("SAME");
			fyZ3->Draw("SAME");
			yMax = 0.0;
			yMax = Max(yMax,h_DT3->GetMaximum());
			yMax = Max(yMax,h_DV3->GetMaximum());
			yMax = Max(yMax,h_DW3->GetMaximum());
			yMax = Max(yMax,h_DZ3->GetMaximum());
			yMin = -0.05 * yMax;
			yMax =  1.05 * yMax;
			h_DT3->GetYaxis()->SetRangeUser(yMin,yMax);
			outfile->WriteTObject(c_T3V3W3Z3);
			
			Double_t tval1 = tCyc;
			Double_t tvaln[1] = {tval1};
			printf("f(%f) = %f\n", tval1, Ycap(3,par,tval1));
			//printf("f(%f) = %f\n", tval1, yAll(tvaln,par));
			
		} // end if (has_VWXY == 1)
		
	} // end if (draw_monte_carlo == 1)
	
//	TFile *outfile = new TFile("BFit.root","recreate");
//	outfile->WriteTObject(c_BFit);
//	if (stBFitCase.bMonteCarlo)	outfile->WriteTObject(c_decays_cyctime);
//	if (stBFitCase.bHasVXWY)	outfile->WriteTObject(c_feeding);
	outfile->Close();
	
	timerStop = clock();
	cout << endl;
	cout << "BFit2 done. Timer = " << (Float_t)timerStop/CLOCKS_PER_SEC << " sec." << endl;
	cout << "Elapsed time = " << (Float_t)(timerStop-timerStart)/CLOCKS_PER_SEC << " sec." << endl << endl;
	
	delete [] lastPar;
	delete [] timeOfCapt;
	delete [] sigmaT1;
	delete [] sigmaT2;
	delete [] sigmaT3;
	delete [] sigmaV1;
	delete [] sigmaV2;
	delete [] sigmaV3;
	delete [] sigmaW1;
	delete [] sigmaW2;
	delete [] sigmaW3;
	delete [] sigmaZ1;
	delete [] sigmaZ2;
	delete [] sigmaZ3;
	delete [] sigmaX2;
	delete [] sigmaX3;
	delete [] sigmaY2;
	delete [] sigmaY3;
	delete [] sY2v1;
	delete [] sY2w1;
	delete [] sY2z1;
	delete [] sY3v2;
	delete [] sY3w2;
	delete [] sY3z2;
	delete [] sY3x2;
	delete [] sY3v1;
	delete [] sY3w1;
	delete [] sY3z1;
	
	return iReturn;
}

void HistPrep (TH1 *h, Int_t rebin, Int_t binWidth, char* pop) {
	h->Rebin(rebin);
	h->SetLineColor(kBlack);
	
	//printf("rebin=%d\n",rebin);
	char xTitle[100], yTitle[100], Title[100];
	sprintf(xTitle,"Cycle time (ms)");
	sprintf(yTitle,"Detections / %d ms", binWidth);
	sprintf(Title,"Simulation vs Model, %s populations", pop);
	h->SetXTitle(xTitle);
	h->SetYTitle(yTitle);
	h->SetTitle(Title);
	h->GetXaxis()->SetRangeUser(0,BFitNamespace::tCyc);
	h->GetYaxis()->SetTitleOffset(1.4);
}

void FuncPrep (TF1 *f, Double_t *pars, Int_t nPoints, Int_t color, Int_t style) {
	f->SetParameters(pars);
	f->SetNpx(nPoints);
	f->SetLineColor(color);
	f->SetLineStyle(style);
}

Double_t intErr (TF1 *fn, Double_t *cov, Double_t t1, Double_t t2) {
	using namespace TMath;
//	extern BFitCase_t stBFitCase;
	Double_t dp, I0 = 0, I1 = 0, I2 = 0, variance = 0;
	Double_t *p0, *p1, *p2;
	Int_t i, j; // param indices
	Int_t		nPar	= fn->GetNpar();
				p0		= fn->GetParameters(); // This won't be changed; should be best-fit values
				p1		= new Double_t [nPar];
				p2		= new Double_t [nPar];
	Double_t	step = 1; // Fraction of 1-sigma in parameter to use as parameter step size
	Double_t	dIdp[nPar]; // to hold gradient of integral wrt each param
	for (i=0; i<nPar; i++) dIdp[i] = 0; // initialize so we can skip some
	// Three lines: Update PDV's; Ensure fn params are right; Get integral at best-fit values
//	BFitNamespace::ComputeParameterDependentVars(p0);
	I0=fn->Integral(t1,t2,p0);
	
	// Loop over params and get derivative of integral wrt each param
	for (i=0; i<nPar; i++) {
		if (cov[i+i*nPar]>0) {// skip fixed params, which have cov(i,i) = 0
			dp = step * Sqrt(cov[i+i*nPar]); // Parameter step size = step * (sigma of i^th param)
			memcpy(p1,p0,nPar*sizeof(Double_t)); // Initialize p1
			memcpy(p2,p0,nPar*sizeof(Double_t)); // Initialize p2
			p1[i] = p0[i] - 0.5*dp; // p1: change i^th param by -dp/2
			p2[i] = p0[i] + 0.5*dp; // p2: change i^th param by +dp/2
			// Get integrals for + and - changes:
			BFitNamespace::ComputeParameterDependentVars(p1);
			I1 = fn->Integral(t1,t2,p1);
			BFitNamespace::ComputeParameterDependentVars(p2);
			I2 = fn->Integral(t1,t2,p2);
			// Get derivative from finite difference about central value
			dIdp[i] = (I2-I1) / dp;
		//	printf("par %2d:  dp=%.4e, p0=%+.4e, p1=%+.4e, p2=%+.4e, I0=%.6e, I1=%.6e, I2=%.6e, dIdp=%+.6e\n", i, dp, p0[i], p1[i], p2[i], I0, I1, I2, dIdp[i]);
		//	for (j=0; j<nPar; j++)
		//		printf("p0[%d]=%+.4e, p1[%d]=%+.4e, p2[%d]=%+.4e\n", j, p0[j], j, p1[j], j, p2[j]);
		}
	}
	// Now get variances:
	for (i=0; i<nPar; i++) {
		if (cov[i+i*nPar]>0) { // skip fixed params, which have cov(i,i) = 0
			for (j=i; j<nPar; j++) {
				if (i==j)	variance +=     dIdp[i] * dIdp[j] * cov[i+j*nPar]; // Update variance --     diagonal elements
				else		variance += 2 * dIdp[i] * dIdp[j] * cov[i+j*nPar]; // Update variance -- off-diagonal elements
			}
		}
	}
	return Sqrt(variance);
}

//int find_struct_index (void *p, char* pcsSearchString, int numStructs, int struct_size ) {
//	// We use a dirty trick here!
//	// We tell the co
//	int struct_index = 0;
//	for (struct_index = 0; struct_index < numStructs; struct_index++ )
//	{
//		printf("%s\n",(char*)p);
//		
//		if (strcmp((char*)p, pcsSearchString) == 0) break; // treat struct as char to get casecode of any struct
//		// pcsCaseCode must be the first member of the struct!!
//		//printf("%d",struct_index);
//		p += struct_size;
//	}
//	if (struct_index == numStructs)
//	{
//		cout << "No match found for case code entered. Execute the program like this:" << endl;
//		cout << "'$ ./PROGRAM <BDN case code>" << endl << endl;
//		return 0;
//	}
//	return struct_index;
//}

//////////////////////////////////////////////////////////////////////////////////////////////
//	B_stBFitCase_t *m;
//	stBDNCase_index = 0;
//	for (m = stBDNCases; stBDNCase_index < iNumStructs_bdn; ++m, ++stBDNCase_index )
//	{
//		if (strcmp(m->pcsCaseCode, argv[2]) == 0) break;
//		printf("%d",stBDNCase_index);		
//	}
//	if (stBDNCase_index == iNumStructs_bdn)
//	{
//		cout << "No match found for casecode entered. Execute the program like this:" << endl;
//		cout << "'$ ./PROGRAM 137i07' for 137-I run 7" << endl << endl;
//		return 0;
//	}
///////////////////////////////////////////////////////////////////////////////////////////////
//	stBDNCase_t *m;
//	for (m = stBDNCases, gdx = 0; m->code != 0 && strcmp(m->code, argv[1]) != 0; ++m, ++gdx);
//	if (m->code == 0) {
//		cout << "No match found for casecode entered. Execute the program like this:" << endl;
//		cout << "'$ ./PROGRAM 137i07' for 137-I run 7" << endl << endl;
//		return 0;
//	}
//	B_fit(m->code, gdx);
//	return 0;
///////////////////////////////////////////////////////////////////////////////////////////////
