//	Function Name:	CSVtoStruct
//
//	Directory Name:	
//
//	Author:		C. Peters
//	Date:		28-Feb-2014
//
//	Installation:	Argonne National Laboratory
//
//	To Compile:	g++ -g -O0 CSVtoStruct.cxx -o CSVtoStruct
//
//	Revision(s):
//
//	Note(s):	This program dynamically allocates a number of values based on NPARAMS
//			which is impored from the CSV file.  These memory locations are NOT
//			'freed' inside the program, since a different function might still need
//			them.  It is ~assumed~ the kernel will free the memory on program exit.
//			
//			Also, there is VERY little error checking!
//	
//			
// 2014-03-02 Shane Caldwell
//	Modified Chris Peters's original code.
//	- Changed some names of things
//	- Removed the main function from this file
//	- Made a new version of CSVtoStruct and ParseToStruct for each type of stuct we need to use
//	- Added FindStructIndex, which Chris also wrote for me 
// Include Files
#include "CSVtoStruct.h"
#include <string>
#include <iostream>
#include <fstream>
#include <cstring>
#include <stdlib.h>
#include "TMath.h"
using namespace std;

// Private Functions
int FindStructIndex ( void*, int, int, char* );
int ParseToStruct_BFit	( char*, BFitCase_t* );
int ParseToStruct_BDN	( char*, BDNCase_t* );

// Begin code....

//int main(int argc, char *argv[])
//{
//	// THIS IS WHAT SHANES PROGRAM HAS TO DO
//
//	int iReturn = SUCCESS;
//	B_fit_case_t stB_fit_case_t[FILE_ROWS];
//
//	int iNumStructs = ReadCSVtoStruct(argv[1], stB_fit_case_t);
//	
//	printf("Sucessfully imported %i rows.\n", iNumStructs);
//	return iReturn;
//}

int FindStructIndex ( void *p, int iStructSize, int iNumStructs, char* pcsSearchString ) {
	///////////////////////////////////////////////////////////////////////////////////////
	// This function exists to turn a case code into the corresponding array index.
	//
	// Inputs:
	//	void *p - Pointer of undefined (void) type
	//			- User must guarantee that this points to an array of structs (eg. BDNCase_t's)
	//			- User must guarantee that the first member of each struct is a case code
	//	int iStructSize - The sizeof() each struct in the array that p points to
	//	int iNumStructs - The number of elements in the array that p points to
	//	char* pcsSearchString - The case code for the data to be found
	//
	// Output:
	//	int iStructIndex - The index of the data for the case code specified
	///////////////////////////////////////////////////////////////////////////////////////
//	cout << "FindStructIndex called." << endl;
	int iStructIndex = 0;
	// Loop over the array of structs
	for (iStructIndex = 0; iStructIndex < iNumStructs; iStructIndex++ )
	{// Pick off first struct member **as a C string** and compare it to our case code
//		printf("%s\n",(char*)p); // diagnostic: see what that C string actually holds
		if (strcmp((char*)p, pcsSearchString) == 0) break; // exit loop with iStructIndex at current value
		p += iStructSize; // advance pointer to next struct
	}
	if (iStructIndex == iNumStructs)
	{// We have reached the end of the struct array without matching the case code
		cout << "No match found for case code entered." << endl;
		return -1; // value for error catching by calling program
	}
	return iStructIndex;
}

//////////////////////////////////////////////////////////////////////////
// Read BDN struct
//////////////////////////////////////////////////////////////////////////
int CSVtoStruct_BDN (char *pcsFileName, BDNCase_t *pstStruct) // change per struct type
{
	int iReturn = SUCCESS;
	char   	*pcsLine = NULL;
	int	iStructIndex = 0;
	FILE *file = fopen(pcsFileName, "r");
	size_t LineLength = 0;
	
	// Allocate the dynamic line buffer
	pcsLine = (char*) malloc(STRING_SIZE);

	// Check if the file opened	
	if (file)
	{
		// Throw away the first line
		getline(&pcsLine, &LineLength, file);

		// Read the rest of the file
		while ( getline(&pcsLine, &LineLength, file) && iStructIndex < FILE_ROWS_BDN) // change per struct type
		{
#ifdef DEBUG		
			printf("\n\nLine[%i]: <%s>\n\n", iStructIndex, pcsLine);
#endif		
			// Call the parse function
			iReturn = ParseToStruct_BDN (pcsLine, &pstStruct[iStructIndex]); // change per struct type

			// Increment the index
			iStructIndex++;
		}
		
		free (pcsLine);
		fclose(file);
	}
	else
	{
		printf("ERROR!  Could not find filename %s in local directory!\n", pcsFileName);
		iReturn = 0;
	}

	return iStructIndex;
}

int ParseToStruct_BDN (char *pcsLine, BDNCase_t *pstStruct) // change per struct type
{
	int iReturn = SUCCESS;
	
	char *pcsResult;
	int iNumParams = 0;
	int iParamIndex;
	
	const double ln2		= 0.69314718056;
	
// Get case code
	pcsResult = strtok(pcsLine,",");
	strcpy(pstStruct->pcsCaseCode, pcsResult);
#ifdef DEBUG
	printf("CaseCode = %s\n", pcsResult);
#endif

// Get file path
	pcsResult = strtok(NULL,",");
	strcpy(pstStruct->pcsFilePath, pcsResult);
#ifdef DEBUG
	printf("FilePath = %s\n", pcsResult);
#endif

// Get file name
	pcsResult = strtok(NULL,",");
	strcpy(pstStruct->pcsFileName, pcsResult);
#ifdef DEBUG
	printf("FileName = %s\n", pcsResult);
#endif

// Get isotope name
	pcsResult = strtok(NULL,",");
	strcpy(pstStruct->pcsIsotopeName, pcsResult);
#ifdef DEBUG
	printf("IsotopeName = %s\n", pcsResult);
#endif

// Get experiment dates
	pcsResult = strtok(NULL,",");
	strcpy(pstStruct->pcsExperimentDates, pcsResult);
#ifdef DEBUG
	printf("ExperimentDates = %s\n", pcsResult);
#endif

// Get log book pages
	pcsResult = strtok(NULL,",");
	strcpy(pstStruct->pcsLogBookPages, pcsResult);
#ifdef DEBUG
	printf("LogBookPages = %s\n", pcsResult);
#endif

// Get run files
	pcsResult = strtok(NULL,",");
	strcpy(pstStruct->pcsRunFiles, pcsResult);
#ifdef DEBUG
	printf("RunFiles = %s\n", pcsResult);
#endif

// Get Precursor name
	pcsResult = strtok(NULL,",");
	strcpy(pstStruct->pcsPrecursorName, pcsResult);
#ifdef DEBUG
	printf("PrecursorName = %s\n", pcsResult);
#endif

// Get Emitter name
	pcsResult = strtok(NULL,",");
	strcpy(pstStruct->pcsEmitterName, pcsResult);
#ifdef DEBUG
	printf("EmitterName = %s\n", pcsResult);
#endif

// Get Daughter name
	pcsResult = strtok(NULL,",");
	strcpy(pstStruct->pcsDaughterName, pcsResult);
#ifdef DEBUG
	printf("DaughterName = %s\n", pcsResult);
#endif

// Get species 1 name
	pcsResult = strtok(NULL,",");
	strcpy(pstStruct->pcsSpecies1Name, pcsResult);
#ifdef DEBUG
	printf("Species1Name = %s\n", pstStruct->pcsSpecies1Name);
#endif

// Get species 2 name
	pcsResult = strtok(NULL,",");
	strcpy(pstStruct->pcsSpecies2Name, pcsResult);
//	strcpy(pstStruct->pcsSpecies2Name, pstStruct->pcsPrecursorName);
#ifdef DEBUG
	printf("Species2Name = %s\n", pstStruct->pcsSpecies2Name);
#endif

// Get species 3 name
	pcsResult = strtok(NULL,",");
	strcpy(pstStruct->pcsSpecies3Name, pcsResult);
//	strcpy(pstStruct->pcsSpecies3Name, pstStruct->pcsEmitterName);
#ifdef DEBUG
	printf("Species3Name = %s\n", pstStruct->pcsSpecies3Name);
#endif

// Get Precursor Mass
	pcsResult = strtok(NULL,",");
	pstStruct->dPrecursorMassAMU = 0.000001*atof(pcsResult);
#ifdef DEBUG
		printf("Precursor Mass = %.15g\n", pstStruct->dPrecursorMassAMU);
#endif

// Get Emitter Mass
	pcsResult = strtok(NULL,",");
	pstStruct->dEmitterMassAMU = 0.000001*atof(pcsResult);
#ifdef DEBUG
		printf("Emitter Mass = %.15g\n", pstStruct->dEmitterMassAMU);
#endif

// Get Daughter Mass
	pcsResult = strtok(NULL,",");
	pstStruct->dDaughterMassAMU = 0.000001*atof(pcsResult);
#ifdef DEBUG
		printf("Daughter Mass = %.15g\n", pstStruct->dDaughterMassAMU);
#endif

// Get species 1 half-life and uncertainty (s)
	pcsResult = strtok(NULL,",");
	pstStruct->dHalfLife1[0] = atof(pcsResult);
	pcsResult = strtok(NULL,",");
	pstStruct->dHalfLife1[1] = atof(pcsResult);
#ifdef DEBUG
		printf("HalfLife1 = %.15g +/- %.15g\n", pstStruct->dHalfLife1[0], pstStruct->dHalfLife1[1]);
#endif

// Get species 2 half-life and uncertainty (s)
	pcsResult = strtok(NULL,",");
	pstStruct->dHalfLife2[0] = atof(pcsResult);
	pcsResult = strtok(NULL,",");
	pstStruct->dHalfLife2[1] = atof(pcsResult);
#ifdef DEBUG
		printf("HalfLife2 = %.15g +/- %.15g\n", pstStruct->dHalfLife2[0], pstStruct->dHalfLife2[1]);
#endif

// Get species 3 half-life and uncertainty (s)
	pcsResult = strtok(NULL,",");
	pstStruct->dHalfLife3[0] = atof(pcsResult);
	pcsResult = strtok(NULL,",");
	pstStruct->dHalfLife3[1] = atof(pcsResult);
#ifdef DEBUG
		printf("HalfLife3 = %.15g +/- %.15g\n", pstStruct->dHalfLife3[0], pstStruct->dHalfLife3[1]);
#endif
	
// Get species 1 lifetime and uncertainty (s)
//	pcsResult = strtok(NULL,",");
//	pstStruct->dLifetime1[0] = atof(pcsResult);
//	pcsResult = strtok(NULL,",");
//	pstStruct->dLifetime1[1] = atof(pcsResult);
	pstStruct->dLifetime1[0] = (pstStruct->dHalfLife1[0])/ln2;
	pstStruct->dLifetime1[1] = (pstStruct->dHalfLife1[1])/ln2;
#ifdef DEBUG
		printf("Lifetime1 = %.15g +/- %.15g\n", pstStruct->dLifetime1[0], pstStruct->dLifetime1[1]);
#endif

// Get species 2 lifetime and uncertainty (s)
//	pcsResult = strtok(NULL,",");
//	pstStruct->dLifetime2[0] = atof(pcsResult);
//	pcsResult = strtok(NULL,",");
//	pstStruct->dLifetime2[1] = atof(pcsResult);
	pstStruct->dLifetime2[0] = (pstStruct->dHalfLife2[0])/ln2;
	pstStruct->dLifetime2[1] = (pstStruct->dHalfLife2[1])/ln2;
#ifdef DEBUG
		printf("Lifetime2 = %.15g +/- %.15g\n", pstStruct->dLifetime2[0], pstStruct->dLifetime2[1]);
#endif

// Get species 3 lifetime and uncertainty (s)
//	pcsResult = strtok(NULL,",");
//	pstStruct->dLifetime3[0] = atof(pcsResult);
//	pcsResult = strtok(NULL,",");
//	pstStruct->dLifetime3[1] = atof(pcsResult);
	pstStruct->dLifetime3[0] = (pstStruct->dHalfLife3[0])/ln2;
	pstStruct->dLifetime3[1] = (pstStruct->dHalfLife3[1])/ln2;
#ifdef DEBUG
		printf("Lifetime3 = %.15g +/- %.15g\n", pstStruct->dLifetime3[0], pstStruct->dLifetime3[1]);
#endif

// Get Right MCP width
	pcsResult = strtok(NULL,",");
	pstStruct->dRightMCPWidth = atof(pcsResult);
#ifdef DEBUG
		printf("Right MCP Width = %.15g\n", pstStruct->dRightMCPWidth);
#endif

// Get Distance to Right MCP Grid
	pcsResult = strtok(NULL,",");
	pstStruct->dRightGridDistance = atof(pcsResult);
#ifdef DEBUG
		printf("Distance to Right MCP Grid = %.15g\n", pstStruct->dRightGridDistance);
#endif

// Get Distance to Right MCP Surface
	pcsResult = strtok(NULL,",");
	pstStruct->dRightMCPDistance = atof(pcsResult);
#ifdef DEBUG
		printf("Distance to Right MCP Surface = %.15g\n", pstStruct->dRightMCPDistance);
#endif

// Get Right MCP bias
	pcsResult = strtok(NULL,",");
	pstStruct->dRightMCPBiasKV = atof(pcsResult);
#ifdef DEBUG
		printf("RightMCP Bias = %.15g\n", pstStruct->dRightMCPBiasKV);
#endif

// Get Top MCP width
	pcsResult = strtok(NULL,",");
	pstStruct->dTopMCPWidth = atof(pcsResult);
#ifdef DEBUG
		printf("Top MCP Width = %.15g\n", pstStruct->dTopMCPWidth);
#endif

// Get Distance to Top MCP Grid
	pcsResult = strtok(NULL,",");
	pstStruct->dTopGridDistance = atof(pcsResult);
#ifdef DEBUG
		printf("Distance to Top MCP Grid = %.15g\n", pstStruct->dTopGridDistance);
#endif

// Get Distance to Top MCP Surface
	pcsResult = strtok(NULL,",");
	pstStruct->dTopMCPDistance = atof(pcsResult);
#ifdef DEBUG
		printf("Distance to Top MCP Surface = %.15g\n", pstStruct->dTopMCPDistance);
#endif

// Get Top MCP bias
	pcsResult = strtok(NULL,",");
	pstStruct->dTopMCPBiasKV = atof(pcsResult);
#ifdef DEBUG
		printf("Top MCP Bias = %.15g\n", pstStruct->dTopMCPBiasKV);
#endif

// Get neutron energy threshold
	pcsResult = strtok(NULL,",");
	pstStruct->dNeutronEnergyThresholdKeV = atof(pcsResult);
#ifdef DEBUG
		printf("Neutron Energy threshold (keV) = %.15g\n", pstStruct->dNeutronEnergyThresholdKeV);
#endif

// Get cycle time
	pcsResult = strtok(NULL,",");
	pstStruct->dCycleTime = atof(pcsResult);
#ifdef DEBUG
		printf("Cycle time = %.15g\n", pstStruct->dCycleTime);
#endif

// Get background time
	pcsResult = strtok(NULL,",");
	pstStruct->dBackgroundTime = atof(pcsResult);
#ifdef DEBUG
		printf("Background time = %.15g\n", pstStruct->dBackgroundTime);
#endif

// Get capture time
	pcsResult = strtok(NULL,",");
	pstStruct->dCaptureTime = atof(pcsResult);
#ifdef DEBUG
		printf("Capture time = %.15g\n", pstStruct->dCaptureTime);
#endif

// Get last capture time
	pcsResult = strtok(NULL,",");
	pstStruct->dLastCaptureTime = atof(pcsResult);
#ifdef DEBUG
		printf("Last capture time = %.15g\n", pstStruct->dLastCaptureTime);
#endif

// Get capture veto duration
	pcsResult = strtok(NULL,",");
	pstStruct->dCaptVetoDurn = atof(pcsResult);
#ifdef DEBUG
		printf("Capture veto duration = %.15g\n", pstStruct->dCaptVetoDurn);
#endif

// Get capture veto overlap
	pcsResult = strtok(NULL,",");
	pstStruct->dCaptVetoOver = atof(pcsResult);
#ifdef DEBUG
		printf("Capture veto overlap with ensuing capture cycle = %.15g\n", pstStruct->dCaptVetoOver);
#endif

// Get per-event deadtime and uncertainty
	pcsResult = strtok(NULL,",");
	pstStruct->dEvtDeadtime[0] = atof(pcsResult);
	pcsResult = strtok(NULL,",");
	pstStruct->dEvtDeadtime[1] = atof(pcsResult);
#ifdef DEBUG
		printf("Per-event deadtime = %.15g +/- %.15g\n", pstStruct->dEvtDeadtime[0], pstStruct->dEvtDeadtime[1]);
#endif

// Get RF frequency
	pcsResult = strtok(NULL,",");
	pstStruct->dRFFrequencyHz = atof(pcsResult);
#ifdef DEBUG
		printf("RF Frequency (Hz) = %.15g\n", pstStruct->dRFFrequencyHz);
#endif

// Get RF meter amplitude
	pcsResult = strtok(NULL,",");
	pstStruct->dRFMeterAmplitude = atof(pcsResult);
#ifdef DEBUG
		printf("RF Meter Amplitude (\"Vac\") = %.15g\n", pstStruct->dRFMeterAmplitude);
#endif

// Get RF function generator amplitude
	pcsResult = strtok(NULL,",");
	pstStruct->dRFFnGenAmplitude = atof(pcsResult);
#ifdef DEBUG
		printf("RF Function Generator Amplitude (Vpp) = %.15g\n", pstStruct->dRFFnGenAmplitude);
#endif

// Get RF amplifier power
	pcsResult = strtok(NULL,",");
	pstStruct->dRFAmpPower = atof(pcsResult);
#ifdef DEBUG
		printf("RF Amplifier Power (Watts) = %.15g\n", pstStruct->dRFAmpPower);
#endif

// Get Species 1 beam rate to CPT
	pcsResult = strtok(NULL,",");
	pstStruct->dBeamRateCPTSpecies1 = atof(pcsResult);
#ifdef DEBUG
		printf("Beam rate of species 1 to CPT = %.15g\n", pstStruct->dBeamRateCPTSpecies1);
#endif

// Get Species 2 beam rate to CPT
	pcsResult = strtok(NULL,",");
	pstStruct->dBeamRateCPTSpecies2 = atof(pcsResult);
#ifdef DEBUG
		printf("Beam rate of species 2 to CPT = %.15g\n", pstStruct->dBeamRateCPTSpecies2);
#endif

// Get Species 3 beam rate to CPT
	pcsResult = strtok(NULL,",");
	pstStruct->dBeamRateCPTSpecies3 = atof(pcsResult);
#ifdef DEBUG
		printf("Beam rate of species 3 to CPT = %.15g\n", pstStruct->dBeamRateCPTSpecies3);
#endif

// Get runtime in seconds
	pcsResult = strtok(NULL,",");
	pstStruct->dRuntimeSec = 3600.0 * atof(pcsResult);
#ifdef DEBUG
		printf("Run time in seconds = %.15g\n", pstStruct->dRuntimeSec);
#endif

// Get number of run files
	pcsResult = strtok(NULL,",");
	pstStruct->nRunFiles = atof(pcsResult);
#ifdef DEBUG
		printf("Number of run files = %.15g\n", pstStruct->nRunFiles);
#endif

// Get number of events
	pcsResult = strtok(NULL,",");
	pstStruct->nEvents = atof(pcsResult);
#ifdef DEBUG
		printf("Number of events = %.15g\n", pstStruct->nEvents);
#endif

// Get number of cycles
	pcsResult = strtok(NULL,",");
	pstStruct->nCycles = atof(pcsResult);
#ifdef DEBUG
		printf("Number of cycles = %.15g\n", pstStruct->nCycles);
#endif

// Get trigger rate
	pcsResult = strtok(NULL,",");
	pstStruct->dAvgTriggerRate = atof(pcsResult);
#ifdef DEBUG
		printf("Average trigger rate -- whole cycle = %.15g\n", pstStruct->dAvgTriggerRate);
#endif

// Get deadtime correction
	pcsResult = strtok(NULL,",");
	pstStruct->dAvgDeadtimeCorrectionFactor = atof(pcsResult);
#ifdef DEBUG
		printf("Average deadtime correction -- whole cycle = %.15g\n", pstStruct->dAvgDeadtimeCorrectionFactor);
#endif

// Get deadtime correction
	pcsResult = strtok(NULL,",");
	pstStruct->dAvgDeadtimeCorrectionFactorBackground = atof(pcsResult);
#ifdef DEBUG
		printf("Average deadtime correction -- background = %.15g\n", pstStruct->dAvgDeadtimeCorrectionFactorBackground);
#endif

// Get deadtime correction
	pcsResult = strtok(NULL,",");
	pstStruct->dAvgDeadtimeCorrectionFactorTrapping = atof(pcsResult);
#ifdef DEBUG
		printf("Average deadtime correction -- trapping = %.15g\n", pstStruct->dAvgDeadtimeCorrectionFactorTrapping);
#endif
	
	using namespace TMath;
	const double c						= 299792.458;	// mm/us
	const double keVperAMU				= 931494.0023;
	const double dElectronMassAMU		= 0.0005485794411963045; // amu
	const double dNeutronMassAMU		= 1.00866491585; // neutron mass in amu
	const double dNeutronMassKeV		= keVperAMU*dNeutronMassAMU;
	const double dElectronMassKeV		= keVperAMU*dElectronMassAMU;
	const double charge					= 2;
//	const double m_ion					= m_atom - *m_electron; // 134Te2+ mass in amu
//	const double mc2_ion				= m_ion*keVperAMU; //keV/c^2
//	const double mc2_factor				= 0.5*m_ion*(m_ion/m_neutron)*keVperAMU; // keV/c^2
	
	pstStruct->dPrecursorMassKeV			= keVperAMU*pstStruct->dPrecursorMassAMU;
	pstStruct->dEmitterMassKeV				= keVperAMU*pstStruct->dEmitterMassAMU;
	pstStruct->dDaughterMassKeV				= keVperAMU*pstStruct->dDaughterMassAMU;
	pstStruct->dFastIonMassKeV				= pstStruct->dDaughterMassKeV - charge*dElectronMassKeV;
	pstStruct->dNeutronEnergyMassFactorKeV	= (pstStruct->dFastIonMassKeV/c/c)*(pstStruct->dFastIonMassKeV/dNeutronMassKeV);
	pstStruct->dQBetaKeV					= pstStruct->dPrecursorMassKeV - pstStruct->dEmitterMassKeV;
	pstStruct->dNeutronSeparationEnergyKeV	= pstStruct->dDaughterMassKeV + dNeutronMassKeV - pstStruct->dEmitterMassKeV;
	pstStruct->dQBetaNeutronKeV				= pstStruct->dQBetaKeV - pstStruct->dNeutronSeparationEnergyKeV;
	pstStruct->dMaxNeutronEnergyKeV			= pstStruct->dQBetaNeutronKeV/(1+(dNeutronMassKeV/pstStruct->dFastIonMassKeV));
	pstStruct->dMaxFastIonEnergyKeV			= pstStruct->dQBetaNeutronKeV/(1+(pstStruct->dFastIonMassKeV/dNeutronMassKeV));
	pstStruct->dMinFastIonEnergyKeV			= pstStruct->dNeutronEnergyThresholdKeV/(1+(pstStruct->dFastIonMassKeV/dNeutronMassKeV));
	pstStruct->dMinFastIonSpeed				= c*Sqrt(2*pstStruct->dMinFastIonEnergyKeV/pstStruct->dFastIonMassKeV);
	pstStruct->dMaxFastIonSpeed				= c*Sqrt(2*pstStruct->dMaxFastIonEnergyKeV/pstStruct->dFastIonMassKeV);
	
// Ballpark estimates -- the outermost location of the slowest Fast Ion: what is the grid position of this ion?
// Used to find the max TOF
// Does not need to be done separately for each MCP, should be close enough
// I get a 1.4mm shift toward the center for the slowest Fast Ion
	const double dMaxMCPTransverse			= Sqrt(2.0 * Power(pstStruct->dRightMCPWidth/2.0, 2.0));
	const double dMaxGridTransverse			= dMaxMCPTransverse -1.4; //-1.4mm to correct for deflection by MCP field
	const double theta						= ATan(dMaxGridTransverse/pstStruct->dRightGridDistance);
	
	pstStruct->dMinFastIonSpeedS			= pstStruct->dMinFastIonSpeed * Sin(theta);
	pstStruct->dMinFastIonSpeedZ			= pstStruct->dMinFastIonSpeed * Cos(theta);
	
	double vMin		= pstStruct->dMinFastIonSpeedZ;
	double vMax		= pstStruct->dMaxFastIonSpeed;
	
// Right MCP calculations
	const double dRightGridGap				= pstStruct->dRightMCPDistance - pstStruct->dRightGridDistance;
	pstStruct->dRightGridAcceleration		= charge*pstStruct->dRightMCPBiasKV/dRightGridGap/(pstStruct->dFastIonMassKeV/c/c);
	const double dMaxRightGridDistance		= Sqrt(Power(pstStruct->dRightGridDistance,2.0) + Power(dMaxGridTransverse,2.0));
	const double dMaxRightMCPDistance		= Sqrt(Power(pstStruct->dRightMCPDistance, 2.0) + Power(dMaxMCPTransverse,2.0));
	double dGridRight						= pstStruct->dRightGridDistance;
	double accelRight						= pstStruct->dRightGridAcceleration;
	//printf("%s: var=%f\n",pstStruct->pcsCaseCode, pstStruct->dNeutronSeparationEnergyKeV);//pstStruct->dMaxFastIonEnergyKeV);
	pstStruct->dRightMCPMinFastIonTOF		= (dGridRight/vMax) + (vMax/accelRight)*(Sqrt(2*accelRight*dRightGridGap/vMax/vMax + 1) - 1);
	pstStruct->dRightMCPMaxFastIonTOF		= (dGridRight/vMin) + (vMin/accelRight)*(Sqrt(2*accelRight*dRightGridGap/vMin/vMin + 1) - 1);
	
// Top MCP calculations
	const double dTopGridGap				= pstStruct->dTopMCPDistance - pstStruct->dTopGridDistance;
	pstStruct->dTopGridAcceleration			= charge*pstStruct->dTopMCPBiasKV/dTopGridGap/(pstStruct->dFastIonMassKeV/c/c);
	const double dMaxTopGridDistance		= Sqrt(Power(pstStruct->dTopGridDistance,2.0) + Power(dMaxGridTransverse,2.0));
	const double dMaxTopMCPDistance			= Sqrt(Power(pstStruct->dTopMCPDistance, 2.0) + Power(dMaxMCPTransverse,2.0));
	double dGridTop							= pstStruct->dTopGridDistance;
	double accelTop							= pstStruct->dTopGridAcceleration;
	pstStruct->dTopMCPMinFastIonTOF			= (dGridTop/vMax) + (vMax/accelTop)*(Sqrt(2*accelTop*dTopGridGap/vMax/vMax + 1) - 1);
	pstStruct->dTopMCPMaxFastIonTOF			= (dGridTop/vMin) + (vMin/accelTop)*(Sqrt(2*accelTop*dTopGridGap/vMin/vMin + 1) - 1);
	
#ifdef DEBUG
	printf("Precursor mass (keV)            = %f\n", pstStruct->dPrecursorMassKeV);
	printf("Emitter mass (keV)              = %f\n", pstStruct->dEmitterMassKeV);
	printf("Daughter mass (keV)             = %f\n", pstStruct->dDaughterMassKeV);
	printf("Fast Ion mass (keV)             = %f\n", pstStruct->dFastIonMassKeV);
	printf("Neutron energy mass factor (keV)= %f\n", pstStruct->dNeutronEnergyMassFactorKeV);
	printf("Electron mass (keV)             = %f\n", dElectronMassKeV);
	printf("Neutron mass (keV)              = %f\n", dNeutronMassKeV);
	printf("Q-Beta (keV)                    = %f\n", pstStruct->dQBetaKeV);
	printf("Neutron separation energy (keV) = %f\n", pstStruct->dNeutronSeparationEnergyKeV);
	printf("Q-Beta-Neutron (keV)            = %f\n", pstStruct->dQBetaNeutronKeV);
	printf("Min Neutron Energy (keV)        = %f\n", pstStruct->dNeutronEnergyThresholdKeV);
	printf("Max Neutron Energy (keV)        = %f\n", pstStruct->dMaxNeutronEnergyKeV);
	printf("Min Fast Ion Energy (keV)       = %f\n", pstStruct->dMinFastIonEnergyKeV);
	printf("Max Fast Ion Energy (keV)       = %f\n", pstStruct->dMaxFastIonEnergyKeV);
	printf("Min Fast Ion Speed              = %f\n", pstStruct->dMinFastIonSpeed);
	printf("Min Fast Ion Speed - S comp.    = %f\n", pstStruct->dMinFastIonSpeedS);
	printf("Min Fast Ion Speed - Z comp.    = %f\n", pstStruct->dMinFastIonSpeedZ);
	printf("Max Fast Ion Speed              = %f\n", pstStruct->dMaxFastIonSpeed);
// Right MCP
	printf("(Right MCP) - Max transverse grid pos. (mm) = %f\n", dMaxGridTransverse);
	printf("(Right MCP) - Cos(Max angle)                = %f\n", Cos(dMaxGridTransverse/pstStruct->dRightGridDistance));
	printf("(Right MCP) - Theta for slowest fast ion    = %f\n", theta);
	printf("(Right MCP) - Transverse shift with min energy & max angle = -%fmm\n", (vMin/accelRight)*(Sqrt(2*accelRight*dRightGridGap/vMin/vMin + 1) - 1)*c*Sqrt(2*pstStruct->dMinFastIonEnergyKeV/pstStruct->dFastIonMassKeV)*Sin(theta));
	printf("Right MCP - Grid acceleration (mm/us/us)    = %f\n", pstStruct->dRightGridAcceleration);
	printf("Right MCP - Min Fast Ion TOF                = %f\n", pstStruct->dRightMCPMinFastIonTOF);
	printf("Right MCP - Max Fast Ion TOF                = %f\n", pstStruct->dRightMCPMaxFastIonTOF);
// Top MCP
	printf("Top MCP - Grid acceleration (mm/us/us)      = %f\n", pstStruct->dTopGridAcceleration);
	printf("Top MCP - Min Fast Ion TOF                  = %f\n", pstStruct->dTopMCPMinFastIonTOF);
	printf("Top MCP - Max Fast Ion TOF                  = %f\n", pstStruct->dTopMCPMaxFastIonTOF);
#endif
	
	return iReturn;
}

//////////////////////////////////////////////////////////////////////////
// Read B_Fit struct
//////////////////////////////////////////////////////////////////////////

int CSVtoStruct_BFit (char *pcsFileName, BFitCase_t *pstStruct) // change per struct type
{
	int iReturn = SUCCESS;
	char   	*pcsLine = NULL;
	int	iStructIndex = 0;
	FILE *file = fopen(pcsFileName, "r");
	size_t LineLength = 0;

	// Allocate the dynamic line buffer
	pcsLine = (char*) malloc(STRING_SIZE);

	// Check if the file opened	
	if (file)
	{
		// Throw away the first line
		getline(&pcsLine, &LineLength, file);

		// Read the rest of the file
		while ( getline(&pcsLine, &LineLength, file) && iStructIndex < FILE_ROWS_BFit) // change per struct type
		{
#ifdef DEBUG		
			printf("\n\nLine[%i]: <%s>\n\n", iStructIndex, pcsLine);
#endif		
			// Call the parse function
			iReturn = ParseToStruct_BFit (pcsLine, &pstStruct[iStructIndex]); // change per struct type

			// Increment the index
			iStructIndex++;
		}
		
		free (pcsLine);
		fclose(file);
	}
	else
	{
		printf("ERROR!  Could not find filename %s in local directory!\n", pcsFileName);
		iReturn = 0;
	}

	return iStructIndex;
}


int ParseToStruct_BFit (char *pcsLine, BFitCase_t *pstStruct) // change per struct type
{
	int iReturn = SUCCESS;

	char *pcsResult;
	int iNumParams = 0;
	int iParamIndex;

	// Get the Casecode
	pcsResult = strtok(pcsLine,",");
	strcpy(pstStruct->pcsCaseCode, pcsResult);
#ifdef DEBUG
	printf("CaseCode = %s\n", pcsResult);
#endif

	// Get the Histname
	pcsResult = strtok(NULL,",");
	strcpy(pstStruct->pcsHistName, pcsResult);
#ifdef DEBUG
	printf("Histname = %s\n", pstStruct->pcsHistName);
#endif

// Get the bin width -- No longer used!
	pcsResult = strtok(NULL,",");
    pstStruct->iBinWidth = atoi(pcsResult);
#ifdef DEBUG
	printf("Bin width is: %d\n", pstStruct->iBinWidth);
#endif

	// Get the number of parameters to follow
	pcsResult = strtok(NULL,",");
    iNumParams = atoi(pcsResult); // used below for looping over parameters
    pstStruct->iNPars = iNumParams;
#ifdef DEBUG
	printf("Total number of parameters is: %d\n", iNumParams);
#endif

	// Allocate that many triplet parameters in the struct
	// NOTE: THESE ARE NOT DE-ALLOCATED ANYWHERE, AND WOULD BE A MEMORY LEAK
	// TO A PROGRAM THAT RUNS FOR A LONG TIME!
	pstStruct->pbToggle = (int*)malloc(iNumParams * sizeof(int));
	pstStruct->pdSeed  =  (double*)malloc(iNumParams * sizeof(double));
	pstStruct->pdStep  =  (double*)malloc(iNumParams * sizeof(double));

	// Now, loop through the rest of the triplet of tokens that many times
	for (iParamIndex = 0; iParamIndex < iNumParams; iParamIndex++)
	{
		// Do the toggle bool
		pcsResult = strtok(NULL,",");
		pstStruct->pbToggle[iParamIndex] = atoi(pcsResult);
#ifdef DEBUG
		printf("bool Toggle[%i] = %i\n", iParamIndex, pstStruct->pbToggle[iParamIndex]);
#endif

		// Do the Seed double
		pcsResult = strtok(NULL,",");
		pstStruct->pdSeed[iParamIndex] = atof(pcsResult);
		//sscanf(pcsResult, "%d", &pstStruct->pdSeed[iParamIndex]);
#ifdef DEBUG
		printf("double Seed[%i]%s = %.15g\n", iParamIndex, pcsResult, pstStruct->pdSeed[iParamIndex]);
#endif

		// Do the Step double
		pcsResult = strtok(NULL,",");
		pstStruct->pdStep[iParamIndex] = atof(pcsResult);
#ifdef DEBUG
		printf("double Step[%i]%s = %.15g\n", iParamIndex, pcsResult, pstStruct->pdStep[iParamIndex]);
#endif

	}

	// Now we get the last few items	
	
	// Options
	pcsResult = strtok(NULL,",");
	strncpy(pstStruct->pcsOptions, pcsResult, STRING_SIZE);
	
	// DoFit
	pcsResult = strtok(NULL,",");
	pstStruct->bDoFit = atoi(pcsResult);
	
	// DoFit
	pcsResult = strtok(NULL,",");
	pstStruct->bUseMyErrorFunction = atoi(pcsResult);
	
	// MonteCarlo
	pcsResult = strtok(NULL,",");
	pstStruct->bMonteCarlo = atoi(pcsResult);
	
	//HasDDC
	pcsResult = strtok(NULL,",");
	pstStruct->bHasDDC = atoi(pcsResult);
	
	// HasCrazyAcronym
	pcsResult = strtok(NULL,",");
	pstStruct->bHasVWXY = atoi(pcsResult);
	
	// computeOtherIntegrals
	pcsResult = strtok(NULL,",");
	pstStruct->bComputeOtherIntegrals = atoi(pcsResult);

#ifdef DEBUG
	printf("Options = %s\nDoFit = %i\nMonteCarlo = %i\nHasDDC = %i\nHasVWXY = %i\nComputeOtherIntegrals = %i\n", 
		pstStruct->pcsOptions, pstStruct->bDoFit, pstStruct->bMonteCarlo, 
		pstStruct->bHasDDC, pstStruct->bHasVWXY, pstStruct->bComputeOtherIntegrals);
#endif

	return iReturn;
}
