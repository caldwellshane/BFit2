#include "BFit2Model.h"
#include "TMath.h"
#include <iostream>
using namespace std;

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// sigmas
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Double_t BFitNamespace::sigmaI (Double_t r, Double_t a, Int_t n) {
	using namespace TMath;
	return (1-Power(r*a,n))/(1-r*a);
}
Double_t BFitNamespace::sigmaII (Double_t r, Double_t a1, Double_t a2, Int_t n) {
	using namespace TMath;
	extern Double_t iota, tCap;
	return 1/(1-r*a1) * ( (1-Power(a2,n-1))/(1-a2) - r*a1 * (Power(a2,n-1)-Power(r*a1,n-1)+iota)/(a2-r*a1+iota) );
}
Double_t BFitNamespace::sigmaIII (Double_t r, Double_t aT1, Double_t aU1, Double_t aU2, Int_t n) {
	using namespace TMath;
	extern Double_t iota, tCap;
	return 1/(1-r*aT1) * ( 1/(1-aU1) * ( (1-Power(aU2,n-1))/(1-aU2) - (Power(aU2,n-1)-Power(aU1,n-1))/(aU2-aU1) ) - (r*aT1+iota)/(aU1-r*aT1+iota) * ( (Power(aU2,n-1)-Power(aU1,n-1))/(aU2-aU1) - (Power(aU2,n-1)-Power(r*aT1,n-1))/(aU2-r*aT1) ) );
}
Double_t BFitNamespace::sigmaIV (Double_t r, Double_t aT1, Double_t aU1, Double_t aU2, Double_t aU3, Int_t n) {
	using namespace TMath;
	extern Double_t iota;
	return 1/(1-r*aT1) * ( 1/(1-aU1) * ( 1/(1-aU2) * ( (1-Power(aU3,n-1))/(1-aU3) - (Power(aU3,n-1)-Power(aU2,n-1))/(aU3-aU2) ) - 1/(aU2-aU1) * ( (Power(aU3,n-1)-Power(aU2,n-1))/(aU3-aU2) - (Power(aU3,n-1)-Power(aU1,n-1))/(aU3-aU1) ) ) - (r*aT1+iota)/(aU1-r*aT1+iota) * ( ( 1/(aU2-aU1) * ( (Power(aU3,n-1)-Power(aU2,n-1))/(aU3-aU2) - (Power(aU3,n-1)-Power(aU1,n-1))/(aU3-aU1) ) ) - 1/(aU2-r*aT1) * ( (Power(aU3,n-1)-Power(aU2,n-1))/(aU3-aU2) - (Power(aU3,n-1)-Power(r*aT1,n-1))/(aU3-r*aT1) ) ) );
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// T populations
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Double_t BFitNamespace::Ttot (Int_t i, Double_t *a, Double_t tvar) {
	using namespace BFitNamespace;
	using namespace TMath;
	extern Double_t tCap, tBac, tCyc, tT1, tT2, tT3;
	extern Double_t ampT1, ampT2, ampT3, *sigmaT1, *sigmaT2, *sigmaT3;
	static Double_t f;
	static Int_t n;
	f = 0.0;
	n = Ceil((tvar-tBac)/tCap);
	if (tvar==tBac) n=1;
	if (tBac <= tvar && tvar <= tCyc) {
		if (i==1) f = ampT1 * sigmaT1[n] * Exp(-(tvar-tBac-(n-1)*tCap)/tT1);
		if (i==2) f = ampT2 * sigmaT2[n] * Exp(-(tvar-tBac-(n-1)*tCap)/tT2);
		if (i==3) f = ampT3 * sigmaT3[n] * Exp(-(tvar-tBac-(n-1)*tCap)/tT3);
	}
	return f;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// U populations -- made of V, W, Z, X, Y
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Double_t BFitNamespace::Utot (Int_t i, Double_t *a, Double_t tvar) {
	using namespace BFitNamespace;
	using namespace TMath;
	extern Double_t tBac, tCyc;
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
	extern Double_t tBac, tCyc;
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
Double_t BFitNamespace::Vtot (Int_t i, Double_t *a, Double_t tvar) {
	extern Double_t tBac, tCyc, tU1, tU2, tU3, V10, V20, V30;
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
Double_t BFitNamespace::Vcap (Int_t i, Double_t *a, Double_t tvar) {
	using namespace BFitNamespace;
	using namespace TMath;
	extern Double_t tCap, tBac, tU1, tU2, tU3;
	extern Double_t ampV1, ampV2, ampV3, *sigmaV1, *sigmaV2, *sigmaV3;
	static Double_t f;
	static Int_t n;
	f = 0.0;
	n = Ceil((tvar-tBac)/tCap);
	if (i==1) f = ampV1 * sigmaV1[n] * Exp(-(tvar-tBac-(n-1)*tCap)/tU1);
	if (i==2) f = ampV2 * sigmaV2[n] * Exp(-(tvar-tBac-(n-1)*tCap)/tU2);
	if (i==3) f = ampV3 * sigmaV3[n] * Exp(-(tvar-tBac-(n-1)*tCap)/tU3);
	return f;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// W populations
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Double_t BFitNamespace::Wtot (Int_t i, Double_t *a, Double_t tvar) {
	extern Double_t tBac, tCyc, tU1, tU2, tU3, W10, W20, W30;
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
Double_t BFitNamespace::Wcap (Int_t i, Double_t *a, Double_t tvar) {
	using namespace BFitNamespace;
	using namespace TMath;
	extern Double_t tCap, tBac, tT1, tT2, tT3, tU1, tU2, tU3;
	extern Double_t ampW1, ampW2, ampW3, *sigmaW1, *sigmaW2, *sigmaW3;
	static Double_t f;
	static Int_t n;
	f = 0.0;
	n = Ceil((tvar-tBac)/tCap);
	if (i==1) f += ampW1 * sigmaW1[n] * Exp(-(tvar-tBac-(n-1)*tCap)/tU1);
	if (i==2) f += ampW2 * sigmaW2[n] * Exp(-(tvar-tBac-(n-1)*tCap)/tU2);
	if (i==3) f += ampW3 * sigmaW3[n] * Exp(-(tvar-tBac-(n-1)*tCap)/tU3);
//	printf("Wcap: i=%d, tvar=%f, par=%f, f=%f\n", i, tvar, SigmaW[a[rho],tT3,tU3,n], f);
	return f;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Z populations
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Double_t BFitNamespace::Ztot (Int_t i, Double_t *a, Double_t tvar) {
	extern Double_t tBac, tCyc, tU1, tU2, tU3, Z10, Z20, Z30;
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
Double_t BFitNamespace::Zcap (Int_t i, Double_t *a, Double_t tvar) {
	using namespace BFitNamespace;
	using namespace TMath;
	extern Double_t tCap, tBac, tT1, tT2, tT3, tU1, tU2, tU3;
	extern Double_t ampZ1, ampZ2, ampZ3, *sigmaT1, *sigmaT2, *sigmaT3, *sigmaZ1, *sigmaZ2, *sigmaZ3;
	static Double_t f;
	static Int_t n;
	f = 0.0;
	n = Ceil((tvar-tBac)/tCap);
	if (i==1) f += ampZ1 * ( sigmaZ1[n]*Exp(-(tvar-tBac-(n-1)*tCap)/tU1) - sigmaT1[n]*Exp(-(tvar-tBac-(n-1)*tCap)/tT1) );
	if (i==2) f += ampZ2 * ( sigmaZ2[n]*Exp(-(tvar-tBac-(n-1)*tCap)/tU2) - sigmaT2[n]*Exp(-(tvar-tBac-(n-1)*tCap)/tT2) );
	if (i==3) f += ampZ3 * ( sigmaZ3[n]*Exp(-(tvar-tBac-(n-1)*tCap)/tU3) - sigmaT3[n]*Exp(-(tvar-tBac-(n-1)*tCap)/tT3) );
	return f;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// X populations
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Double_t BFitNamespace::Xtot (Int_t i, Double_t *a, Double_t tvar) {
	extern Double_t tBac, tCyc, tU2, tU3, X20, X30;
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
Double_t BFitNamespace::Xcap (Int_t i, Double_t *a, Double_t tvar) {
	using namespace BFitNamespace;
	using namespace TMath;
	extern Double_t tCap, tBac, tT1, tT2, tU2, tU3;
	extern Double_t ampX2, ampX3, *sigmaT1, *sigmaT2, *sigmaX2, *sigmaX3;
	static Double_t f;
	static Int_t n;
	f = 0.0;
	n = Ceil((tvar-tBac)/tCap);
	if (i==2) f += ampX2 * ( sigmaX2[n] * Exp(-(tvar-tBac-(n-1)*tCap)/tU2) - sigmaT1[n] * Exp(-(tvar-tBac-(n-1)*tCap)/tT1) );
	if (i==3) f += ampX3 * ( sigmaX3[n] * Exp(-(tvar-tBac-(n-1)*tCap)/tU3) - sigmaT2[n] * Exp(-(tvar-tBac-(n-1)*tCap)/tT2) );
	return f;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Y populations
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Double_t BFitNamespace::Ytot (Int_t i, Double_t *a, Double_t tvar) {
	extern Int_t nCap;
	extern Double_t tBac, tCyc, tU1, tU2, tU3, U10, U20, Y20, Y30;
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
	extern Double_t tBac, tCyc, t1, t2, tU1, tU2, tU3, U10, U20, Y20, Y30;
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
	extern Double_t tCap, tBac, t1, t2, tT1, tT2, tU1, tU2, tU3, aT1, aU1, aU2, aU3, tU1U2, tT1U2, tU2U3, tT2U3, tU1U3, tT1U3;
	extern Double_t ampV1, ampW1, ampZ1, ampV2, ampW2, ampZ2, ampX2, ampY2;
	extern Double_t *sigmaT1, *sigmaV1, *sigmaW1, *sigmaZ1, *sigmaT2, *sigmaV2, *sigmaW2, *sigmaZ2, *sigmaX2, *sigmaY2, *sigmaY3, *sY2v1, *sY2w1, *sY2z1, *sY3w1;
	static Double_t tn, expT1, expT2, expU1, expU2, expU3, f;
	static Int_t n;
	f = 0.0;
	n = Ceil((tvar-tBac)/tCap);
	//if (tvar==tBac) n=1;
	tn = tvar-tBac-(n-1)*tCap;
	if (i==2) {
		expT1 = Exp(-tn/tT1);
		expU1 = Exp(-tn/tU1);
		expU2 = Exp(-tn/tU2);
		f = sigmaY2[n] * expU2
			+ ampV1 * sigmaV1[n] * tU1U2/t1 * ( expU2 - expU1 )
			+ ampW1 * sigmaW1[n] * tU1U2/t1 * ( expU2 - expU1 )
			+ ampZ1 * sigmaZ1[n] * tU1U2/t1 * ( expU2 - expU1 )
			- ampZ1 * sigmaT1[n] * tT1U2/t1 * ( expU2 - expT1 );
		if (IsNaN(f) && tvar==tBac) printf("t=%f, Y2 =%f\n", tvar, f);
	}
	if (i==3) {
		expT1 = Exp(-tn/tT1);
		expT2 = Exp(-tn/tT2);
		expU1 = Exp(-tn/tU1);
		expU2 = Exp(-tn/tU2);
		expU3 = Exp(-tn/tU3);
		f = sigmaY3[n] * expU3
				// Components from V2, W2, Z2, X2
			+ ampV2 * sigmaV2[n] * tU2U3/t2 * ( expU3 - expU2 )
			+ ampW2 * sigmaW2[n] * tU2U3/t2 * ( expU3 - expU2 )
			+ ampZ2 * sigmaZ2[n] * tU2U3/t2 * ( expU3 - expU2 )
			- ampZ2 * sigmaT2[n] * tT2U3/t2 * ( expU3 - expT2 )
			+ ampX2 * sigmaX2[n] * tU2U3/t2 * ( expU3 - expU2 )
			- ampX2 * sigmaT1[n] * tT1U3/t2 * ( expU3 - expT1 )
				// Components from Y2 <== V1, W1, Z1
			+ ampV1 *   sY2v1[n] *            tU2U3/t2 * ( expU3 - expU2 )
			+ ampV1 * sigmaV1[n] * tU1U2/t1 * tU2U3/t2 * ( expU3 - expU2 )
			- ampV1 * sigmaV1[n] * tU1U2/t1 * tU1U3/t2 * ( expU3 - expU1 )
			+ ampW1 *   sY2w1[n] *            tU2U3/t2 * ( expU3 - expU2 )
			+ ampW1 * sigmaW1[n] * tU1U2/t1 * tU2U3/t2 * ( expU3 - expU2 )
			- ampW1 * sigmaW1[n] * tU1U2/t1 * tU1U3/t2 * ( expU3 - expU1 )
			+ ampZ1 *   sY2z1[n] *            tU2U3/t2 * ( expU3 - expU2 )
			+ ampZ1 * sigmaZ1[n] * tU1U2/t1 * tU2U3/t2 * ( expU3 - expU2 )
			- ampZ1 * sigmaZ1[n] * tU1U2/t1 * tU1U3/t2 * ( expU3 - expU1 )
			- ampZ1 * sigmaT1[n] * tT1U2/t1 * tU2U3/t2 * ( expU3 - expU2 )
			+ ampZ1 * sigmaT1[n] * tT1U2/t1 * tT1U3/t2 * ( expU3 - expT1 );
		//printf("t=%f, Y3=%f\n", tvar, f);
		if (IsNaN(f) && tvar==tBac) printf("t=%f, Y3 =%f\n", tvar, f);
	}
	return f;
}

