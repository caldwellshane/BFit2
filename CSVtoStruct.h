//Function Name:  ReadCSVtoStruct
//
//      Directory Name: 
//
//      Author:         C. Peters
//      Date:           28-Feb-2014
//
//      Installation:   Argonne National Laboratory
//
//      To Compile:     g++ -g -O0 ReadCSVtoStruct.cxx -o ReadCSVtoStruct
//
//      Revision(s):
//
//      Note(s):        This program dynamically allocates a number of values based on NPARAMS
//                      which is impored from the CSV file.  These memory locations are NOT
//                      'freed' inside the program, since a different function might still need
//                      them.  It is ~assumed~ the kernel will free the memory on program exit.
//      
//  			Also, there is VERY little error checking!
//
// 2014-02-28 Shane Caldwell
//	Adding bdn_case_t to these definitions.
//	Notation: every struct member starts with a code indicating the variable type"
//		- p  = pointer
//		- cs = C string
//		- d  = double
//		- i  = int
//		- b  = bool, but these are actually C strings. Usage: if (bStringVar)
//	EVERY STRUCT MUST HAVE THE CASE CODE AS A char AS ITS FIRST MEMBER

// Configuration defines
#define SUCCESS 0
#define ERROR   1
//#define DEBUG
#define STRING_SIZE     4096
#define FILE_ROWS_BDN     68 // number of structs in stBDNCases
#define FILE_ROWS_BFit    48 // number of structs in stBFitCases

struct BDNCase_t
{
///////////////////////////////////////////////////////////////////////////////////////////////
// IMPORTED VALUES
///////////////////////////////////////////////////////////////////////////////////////////////
	char pcsCaseCode[STRING_SIZE]; // Case code
	char pcsFilePath[STRING_SIZE]; // Full path of ROOT file containing summed data
	char pcsFileName[STRING_SIZE]; // ROOT file containing summed data
	char pcsIsotopeName[STRING_SIZE]; // BDN precursor of interest
	char pcsExperimentDates[STRING_SIZE]; // Dates of experiment
	char pcsLogBookPages[STRING_SIZE]; // Log book pages for the experiment
	char pcsRunFiles[STRING_SIZE]; // Run files for the experiment
// Species names
	char pcsPrecursorName[STRING_SIZE]; // Species 1: Parent of nucleus of interest
	char pcsEmitterName[STRING_SIZE]; // Species 2: Nucleus of interest
	char pcsDaughterName[STRING_SIZE]; // Species 3: Daugter of nucleus of interest
// Nearest isobars
	char pcsSpecies1Name[STRING_SIZE]; // Species 1: Parent of nucleus of interest
	char pcsSpecies2Name[STRING_SIZE]; // Species 2: Nucleus of interest
	char pcsSpecies3Name[STRING_SIZE]; // Species 3: Daugter of nucleus of interest
// Masses; array holding {value, sigma}
	double dPrecursorMassAMU; // Atomic mass of BDN precursor (AMU)
    double dEmitterMassAMU; // Atomic mass of BDN emitter (AMU)
    double dDaughterMassAMU; // Atomic mass of BDN daughter (AMU)
// Half-lives; array holding {value, sigma}
	double dHalfLife1[2]; // Two-element array holding half-life [0] and uncertainty [1] for species 1
    double dHalfLife2[2]; // Two-element array holding half-life [0] and uncertainty [1] for species 2
    double dHalfLife3[2]; // Two-element array holding half-life [0] and uncertainty [1] for species 3
// MCP stuff
	double dRightMCPWidth; // Width of MCP (mm)
	double dRightGridDistance; // Trap center to MCP grid (mm)
	double dRightMCPDistance; // Trap center to MCP surface (mm)
	double dRightMCPBiasKV; // MCP bias in kV
	double dTopMCPWidth; // Width of MCP (mm)
	double dTopGridDistance; // Trap center to MCP grid (mm)
	double dTopMCPDistance; // Trap center to MCP surface (mm)
	double dTopMCPBiasKV; // MCP bias in kV
// Neutron energy threshold
	double dNeutronEnergyThresholdKeV; // Neutron energy threshold
// Timing constants describing the experiment, in ms
	double dCycleTime; // BPT cycle time (background + trapping)
	double dBackgroundTime; // Duration of background measurement
	double dCaptureTime; // Time between captures
	double dLastCaptureTime; // Time between last capture and ejection
// Timing constants affecting deadtime {value, sigma}
	double dCaptVetoDurn; // Overall duration of capture pulse veto (us)
	double dCaptVetoOver; // Overlap of capture pulse veto with the ensuing capture cycle (us)
	double dEvtDeadtime[2]; // Deadtime per event and uncertainty (us)
// Timing constants affecting deadtime {value, sigma}
	double dRFFrequencyHz; // Trap RF frequency (kHz)
	double dRFMeterAmplitude; // Trap RF amplitude from precision meter ("Vac")
	double dRFFnGenAmplitude; // Trap RF amplitude on function generator (Vpp)
	double dRFAmpPower; // Trap RF power from amplifier (Watts)
// Case metadata
	double dBeamRateCPTSpecies1;
	double dBeamRateCPTSpecies2;
	double dBeamRateCPTSpecies3;
	double dRuntimeSec;
	double nRunFiles;
	double nEvents;
	double nCycles;
// Deadtime corrections
	double dAvgTriggerRate;
	double dAvgDeadtimeCorrectionFactor;
	double dAvgDeadtimeCorrectionFactorBackground;
	double dAvgDeadtimeCorrectionFactorTrapping;
///////////////////////////////////////////////////////////////////////////////////////////////
// CALCULATED VALUES
///////////////////////////////////////////////////////////////////////////////////////////////
// 1/e lifetimes; array holding {value, sigma}
	double dLifetime1[2]; // Two-element array holding (1/e) lifetime [0] and uncertainty [1] for species 1
	double dLifetime2[2]; // Two-element array holding (1/e) lifetime [0] and uncertainty [1] for species 2
	double dLifetime3[2]; // Two-element array holding (1/e) lifetime [0] and uncertainty [1] for species 3
// Masses in keV
	double dPrecursorMassKeV;
	double dEmitterMassKeV;
	double dDaughterMassKeV;
// Kinematics
	double dFastIonMassKeV;
	double dNeutronEnergyMassFactorKeV;
	double dQBetaKeV;
	double dNeutronSeparationEnergyKeV;
	double dQBetaNeutronKeV;
	double dMaxNeutronEnergyKeV;
	double dMinFastIonEnergyKeV;
	double dMaxFastIonEnergyKeV;
	double dMinFastIonSpeed;
	double dMinFastIonSpeedS;
	double dMinFastIonSpeedZ;
	double dMaxFastIonSpeed;
	double dRightGridAcceleration;
	double dTopGridAcceleration;
	double dRightMCPMinFastIonTOF;
	double dRightMCPMaxFastIonTOF;
	double dTopMCPMinFastIonTOF;
	double dTopMCPMaxFastIonTOF;
};

struct BFitCase_t
{
	char pcsCaseCode[STRING_SIZE];
	char pcsHistName[STRING_SIZE];
	int iBinWidth;
	int iNPars;
	int *pbToggle;
	double *pdSeed;
	double *pdStep;
	char pcsOptions[STRING_SIZE];
	char bDoFit;
	char bUseMyErrorFunction;
	char bMonteCarlo;
	char bHasDDC;
	char bHasVWXY;
	char bComputeOtherIntegrals;
};

//Public function prototypes
int FindStructIndex  ( void *p, int iStructSize, int iNumStructs, char* pcsSearchString );
int CSVtoStruct_BDN  ( char *pcsFileName, BDNCase_t  *pstStruct );
int CSVtoStruct_BFit ( char *pcsFileName, BFitCase_t *pstStruct );
