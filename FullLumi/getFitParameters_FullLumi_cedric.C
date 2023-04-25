#include <cstdio>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include "TFile.h"
#include <TCanvas.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TMinuit.h>
#include "TH1.h"
#include "TF1.h"
#include "Riostream.h"
#include "TGraphErrors.h"
#include "TROOT.h"
#include "TStyle.h"
#include <TH2.h>
#include <TH2F.h>
#include <TH3.h>
#include <TH3F.h>
#include "TF2.h"
#include "TGraph2D.h"
#include "TMath.h"
#include "TFitter.h"
#include <string>
#include <iomanip>


void getResults(int NTOY, TString cuts, 
                Double_t *Kbest_fit, Double_t *DKbest_fit,
                Double_t *Mbest_fit, Double_t *DMbest_fit, 
                Double_t *corrKM_fit, Double_t *DcorrKM_fit);

//Reference values for the input parameters
const Double_t Kref = 0.13726; const Double_t dKu = 0.00081;
const Double_t Mref = 0.0525;  const Double_t dMu = 0.003;

//range and step for the parameters grid (they depend on the specific job settings)
const Int_t sigmaLim  = 5; const Int_t sigmaStep = 4;
const Int_t ngrid = sigmaLim*sigmaStep*2 + 1;

Double_t oneSigma_MS    =  0.02;//1.0;//************************** CHANGE THIS VALUE ACCORDING TO THE +/-1sigma SHIFT USED FOR THE MULTIPLE SCATTERING MODELIZATION
Double_t offset_MS      =  1.01;//0.0;//************************** CHANGE THIS VALUE ACCORDING TO THE NOMINAL MODELIZATION OF MULTIPLE SCATTERING EFFECTS

TGraph *g_nLL; TCanvas *c_nLL;
TH2D *h_NormNuisance;
TGraphErrors *NormNuisanceM[ngrid]; TCanvas *cM_NormNuisance;
TGraphErrors *NormNuisanceK[ngrid]; TCanvas *cK_NormNuisance;
TH2D *h_MultipleScatteringNuisance; TCanvas *c_MultipleScatteringNuisance;
TGraphErrors *MSNuisanceM[ngrid]; TCanvas *cM_MSNuisance;
TGraphErrors *MSNuisanceK[ngrid]; TCanvas *cK_MSNuisance;
TH2D *h_SingleHitResNuisance; TCanvas *c_SingleHitResNuisance;

//interpolateNLL graphs
TH2D *h_nLL2D; TCanvas *c_nLL2D;
TGraphErrors *resultsM[ngrid]; TCanvas *cM; TF1 *fitfunction[ngrid];
TGraphErrors *corr_KM; TCanvas *cKM;
TGraphErrors *resultsK; TCanvas *cK;


//ŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋ    FUNCTIONS TO CALCULATE aµHLO    ŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋ

const double alpha = 1./137.0359991;
const double mmu   = 0.1056583715;//GeV
const double me    = 0.51099906e-3;
const double Emu   = 150;//GeV
//  t(x)
Double_t t_vs_x(Double_t x) {
        Double_t mu2=mmu*mmu;
        return mu2*(x)*(x)/(x-1);
}

//  x(t)
Double_t x_vs_t(Double_t x) {
        double mu2=mmu*mmu;
         Double_t t = x;
         return t/(2*mmu*mmu)*(1 - sqrt(1 - 4*mmu*mmu/t));
}

//  Ee(theta_e)
double Eevsth(double thetae) {
	double r = sqrt(Emu*Emu-mmu*mmu)/(Emu+me);
        return  me*(1.+r*r*cos(thetae)*cos(thetae))/(1-r*r*cos(thetae)*cos(thetae));
}

//  t(theta_e) 
double t_vs_the(double thetae) {
	double Ee = Eevsth(thetae);
        return 2*me*me - 2*me*Ee;
}

TF1 *vphad_K;
Double_t fvphad_K(Double_t *x, Double_t *par) {
	Double_t K = *x;
        Double_t t = par[0];
        Double_t M = par[1];
        Double_t squ = TMath::Sqrt(1. - 4.*M/t);
        return K*M*( - 5./9. - 4./3.*M/t + 2*(4./3.*M*M/t/t + M/(3.*t) - 1./6.)/squ*TMath::Log(TMath::Abs((1. - squ)/(1. + squ))) );
}

TF1* vphad_M;
Double_t fvphad_M(Double_t *x, Double_t *par) {
	Double_t M = *x;
        Double_t t = par[0];
        Double_t K = par[1];
        Double_t squ = TMath::Sqrt(1. - 4.*M/t);
        return K*M*( - 5./9. - 4./3.*M/t + 2*(4./3.*M*M/t/t + M/(3.*t) - 1./6.)/squ*TMath::Log(TMath::Abs((1. - squ)/(1. + squ))) );
}

Double_t amuhlo(Double_t *x, Double_t *par) {
        Double_t t = t_vs_x(*x);
        Double_t K = par[0];
        Double_t M = par[1];
        Double_t squ = TMath::Sqrt(1. - 4.*M/t);
        Double_t vphad = K*M*( - 5./9. - 4./3.*M/t + 2*(4./3.*M*M/t/t + M/(3.*t) - 1./6.)/squ*TMath::Log(TMath::Abs((1. - squ)/(1. + squ))) );
        return alpha/TMath::Pi()*(1 - *x)*vphad;
}

Double_t error_amuhlo(Double_t *x, Double_t *par) {
        Double_t t = t_vs_x(*x);
        Double_t K      = par[0];
        Double_t M      = par[1];
        Double_t DK     = par[2];
        Double_t DM     = par[3];
        Double_t corrKM = par[4];

        Double_t squ = TMath::Sqrt(1. - 4.*M/t);

        Double_t params[2] = {t, M};
        Double_t Dvphad_DK = vphad_K->Derivative(K, params);
        params[1] = K;
        Double_t Dvphad_DM = vphad_M->Derivative(M, params);

        Double_t Dvphad = TMath::Sqrt( TMath::Power(Dvphad_DK*DK, 2) + TMath::Power(Dvphad_DM*DM, 2) + 2*corrKM*DK*DM*Dvphad_DK*Dvphad_DM );
        return alpha/TMath::Pi()*(1 - *x)*Dvphad;
}

//ŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋŋ




void getFitParameters_FullLumi(int NTOY= 1, TString cuts = "thmu0.2_the32") {

	Double_t KBest  = 0, DKBest  = 0;
	Double_t MBest  = 0, DMBest  = 0;
	Double_t corrKM = 0, DcorrKM = 0;
	getResults(NTOY, cuts, &KBest, &DKBest, &MBest, &DMBest, &corrKM, &DcorrKM);

	//find KBest_index
	Double_t K0 = Kref - sigmaLim*dKu;//Kreal for iK = 0
	double KBest_index = -1;
	double fractional_part = modf((KBest - K0)/(dKu/sigmaStep), &KBest_index);
	if(fractional_part > 0.5) KBest_index++;
	if(KBest_index < 0) KBest_index = 0;
	if(KBest_index >= ngrid) KBest_index = ngrid-1;

	//find KBest_1sigma_index
	double KBest_1sigma_index = -1;
	fractional_part = modf(((KBest + DKBest) - K0)/(dKu/sigmaStep), &KBest_1sigma_index);
	if(fractional_part > 0.5) KBest_1sigma_index++;
	if(KBest_1sigma_index < 0) KBest_1sigma_index = 0;
	if(KBest_1sigma_index >= ngrid) KBest_1sigma_index = ngrid-1;

	//find MBest_index
	Double_t M0 = Mref - sigmaLim*dMu;//Mreal for iM = 0
	double MBest_index = -1;
	fractional_part = modf((MBest - M0)/(dMu/sigmaStep), &MBest_index);
	if(fractional_part > 0.5) MBest_index++;
	if(MBest_index < 0) MBest_index = 0;
	if(MBest_index >= ngrid) MBest_index = ngrid-1;

	//find MBest_1sigma_index
	double MBest_1sigma_index = -1;
	fractional_part = modf(((MBest + DMBest) - M0)/(dMu/sigmaStep), &MBest_1sigma_index);
	if(fractional_part > 0.5) MBest_1sigma_index++;
	if(MBest_1sigma_index < 0) MBest_1sigma_index = 0;
	if(MBest_1sigma_index >= ngrid) MBest_1sigma_index = ngrid-1;

	cout<<"KBest_index = "<<KBest_index<<", MBest_index = "<<MBest_index<<endl;
	cout<<"KBest_1sigma_index = "<<KBest_1sigma_index<<", MBest_1sigma_index = "<<MBest_1sigma_index<<endl;

	cout<<"KBest  = "<<KBest<<" +/- "<<DKBest
	  <<"\nMBest  = "<<MBest<<" +/- "<<DMBest
	  <<"\ncorrKM = "<<corrKM<<" +/- "<<DcorrKM<<endl;

	//calculate amuHLO expected
	Double_t t_max = t_vs_the(0);//GeV2
        Double_t x_max = x_vs_t(t_max);
        cout<<"tmax = "<<t_max<<" GeV2; xmax = "<<x_max<<endl;

        TF1 *func_amuhlo = new TF1("func_amuhlo", amuhlo, 0, 1, 2);
        func_amuhlo->SetParameters(Kref, Mref);
        func_amuhlo->SetNpx(1e5);
        Double_t amuhlo_100_ref = func_amuhlo->Integral(0, 1);
        Double_t amuhlo_87_ref  = func_amuhlo->Integral(0, x_max);
        Double_t amuhlo_13_ref  = func_amuhlo->Integral(x_max, 1);
        Double_t amuhlo_err = 3e-10;

        cout<<"REFERENCE VALUES:"<<endl;
        cout<<"amuHLO = "<<amuhlo_100_ref<<endl;
        cout<<"amuHLO(87%) = "<<amuhlo_87_ref<<endl;
        cout<<"amuHLO(13%) = "<<amuhlo_13_ref<<endl;

	//calculate amuHLO measured
	func_amuhlo->SetParameters(KBest, MBest);
        Double_t amuhlo_100 = func_amuhlo->Integral(0, 1);
        Double_t amuhlo_87  = func_amuhlo->Integral(0, x_max);
        Double_t amuhlo_13  = func_amuhlo->Integral(x_max, 1);

        //calculate error on amuHLO measured
        TF1 *func_error_amuhlo = new TF1("func_error_amuhlo", error_amuhlo, 0, 1, 5);
        func_error_amuhlo->SetParameters(KBest, MBest, DKBest, DMBest, corrKM);
        func_error_amuhlo->SetNpx(1e5);

        vphad_K = new TF1("vphad_K", fvphad_K, Kref - 5*dKu, Kref + 5*dKu, 2);
        vphad_K->SetParameters(t_max, MBest);
        vphad_K->SetNpx(1e5);

        vphad_M = new TF1("vphad_M", fvphad_M, Mref - 5*dMu, Mref + 5*dMu, 2);
        vphad_M->SetParameters(t_max, KBest);
        vphad_M->SetNpx(1e5);

        Double_t Damuhlo_100 = func_error_amuhlo->Integral(0, 1);
        Double_t Damuhlo_87  = func_error_amuhlo->Integral(0, x_max);
        Double_t Damuhlo_13  = func_error_amuhlo->Integral(x_max, 1);

        cout<<"\nMEASURED VALUES:"<<endl;
        cout<<"amuHLO = "<<amuhlo_100<<" +/- "<<Damuhlo_100<<endl;
        cout<<"amuHLO(87%) = "<<amuhlo_87<<" +/- "<<Damuhlo_87<<endl;
        cout<<"amuHLO(13%) = "<<amuhlo_13<<" +/- "<<Damuhlo_13<<endl;

	//EVALUATE ERROR ON AMUHLO
        TF2 *pdf = new TF2("pdf", "bigaus");
        pdf->SetRange(KBest - 10*DKBest, MBest - 10*DMBest, KBest + 10*DKBest, MBest + 10*DMBest);
        pdf->SetParameters(10, KBest, DKBest, MBest, DMBest, corrKM);
        pdf->SetNpx(1e3);
        pdf->SetNpy(1e3);
        TH1D *h_amuhlo100 = new TH1D("h_amuhlo100", "a_{#mu}^{HLO} distribution; a_{#mu}^{HLO}; Entries", 100, amuhlo_100 - 10e-10, amuhlo_100 + 10e-10);
        TH1D *h_amuhlo87  = new TH1D("h_amuhlo87",  "a_{#mu}^{HLO}(87%) distribution; a_{#mu}^{HLO}(87%); Entries", 100, amuhlo_87 - 10e-10, amuhlo_87 + 10e-10);
        TH1D *h_amuhlo13  = new TH1D("h_amuhlo13",  "a_{#mu}^{HLO}(13%) distribution; a_{#mu}^{HLO}(13%); Entries", 100, amuhlo_13 - 10e-11, amuhlo_13 + 10e-11);

        for(int i = 0; i < 1e5; i++) {
                double K = 0, M = 0;
                pdf->GetRandom2(K, M);
                func_amuhlo->SetParameters(K, M);
                h_amuhlo100->Fill(func_amuhlo->Integral(0, 1));
                h_amuhlo87->Fill(func_amuhlo->Integral(0, x_max));
                h_amuhlo13->Fill(func_amuhlo->Integral(x_max, 1));
        }

        TFitResultPtr fitresult = h_amuhlo100->Fit("gaus", "QS");
        amuhlo_100 = fitresult->Parameter(1);
        Damuhlo_100 = fitresult->Parameter(2);

        fitresult = h_amuhlo87->Fit("gaus", "QS");
        amuhlo_87 = fitresult->Parameter(1);
        Damuhlo_87 = fitresult->Parameter(2);

        fitresult = h_amuhlo13->Fit("gaus", "QS");
        amuhlo_13 = fitresult->Parameter(1);
        Damuhlo_13 = fitresult->Parameter(2);

        cout<<"\nMEASURED VALUES:"<<endl;
        cout<<"amuHLO = "<<amuhlo_100<<" +/- "<<Damuhlo_100<<endl;
        cout<<"amuHLO(87%) = "<<amuhlo_87<<" +/- "<<Damuhlo_87<<endl;
        cout<<"amuHLO(13%) = "<<amuhlo_13<<" +/- "<<Damuhlo_13<<endl;



	TF1 *f_NormNuisance = new TF1("f_Norm", "pol1", -2*sigmaLim, 2*sigmaLim);
	f_NormNuisance->SetNpx(1e3);
	f_NormNuisance->SetParameters(1, 1);

	NormNuisanceK[static_cast<int>(MBest_index)]->Fit(f_NormNuisance, "0");
	Double_t NormNuisance = f_NormNuisance->Eval( (KBest - Kref)/dKu );
	Double_t NormNuisance_1sigma_K = f_NormNuisance->Eval( (KBest + DKBest - Kref)/dKu );

	TF1 *f_NormNuisance_M = new TF1("f_Norm_M", "pol2", -2*sigmaLim, 2*sigmaLim);
	f_NormNuisance_M->SetNpx(1e3);
	f_NormNuisance_M->SetParameters(1, -1, 1);
	NormNuisanceM[static_cast<int>(KBest_index)]->Fit(f_NormNuisance_M, "0");
	Double_t NormNuisance_1sigma_M = f_NormNuisance_M->Eval( (MBest + DMBest - Mref)/dMu );


	TF1 *f_MultipleScatteringNuisance = new TF1("f_MultipleScattering", "pol1", -2*sigmaLim, 2*sigmaLim);
	f_MultipleScatteringNuisance->SetNpx(1e3);
	f_MultipleScatteringNuisance->SetParameters(1, 1);
	MSNuisanceK[static_cast<int>(MBest_index)]->Fit(f_MultipleScatteringNuisance, "0");
	Double_t MultipleScatteringNuisance = f_MultipleScatteringNuisance->Eval( (KBest - Kref)/dKu );
	Double_t MultipleScatteringNuisance_1sigma_K = f_MultipleScatteringNuisance->Eval( (KBest + DKBest - Kref)/dKu );

	TF1 *f_MultipleScatteringNuisance_M = new TF1("f_MultipleScattering_M", "pol2", -2*sigmaLim, 2*sigmaLim);
	f_MultipleScatteringNuisance_M->SetNpx(1e3);
	f_MultipleScatteringNuisance_M->SetParameters(1, -1, 1);
	MSNuisanceM[static_cast<int>(KBest_index)]->Fit(f_MultipleScatteringNuisance_M, "0");
	Double_t MultipleScatteringNuisance_1sigma_M = f_MultipleScatteringNuisance_M->Eval( (MBest + DMBest - Mref)/dMu );



/*	Double_t NormNuisance = h_NormNuisance->GetBinContent(KBest_index+1, MBest_index+1);
	Double_t MultipleScatteringNuisance = h_MultipleScatteringNuisance->GetBinContent(KBest_index+1, MBest_index+1);
	Double_t NormNuisance_1sigma_K = h_NormNuisance->GetBinContent(KBest_1sigma_index+1, MBest_index+1);
	Double_t MultipleScatteringNuisance_1sigma_K = h_MultipleScatteringNuisance->GetBinContent(KBest_1sigma_index+1, MBest_index+1);
	Double_t NormNuisance_1sigma_M = h_NormNuisance->GetBinContent(KBest_index+1, MBest_1sigma_index+1);
	Double_t MultipleScatteringNuisance_1sigma_M = h_MultipleScatteringNuisance->GetBinContent(KBest_index+1, MBest_1sigma_index+1);
*/



        //build the full covariance matrix V:
        //
        //            [ A   B ]
        // H = V^-1 = [ BT  C ]
        //
        //
        //            [ Vnuisance       Vcorrelation ]
        // V        = [ VcorrelationT       Vkm      ]
	//
	//
	// SIGMA    = [ sigmaK    0   ]
	//            [   0    sigmaM ]
	//
	//
        // A^-1 = NuisanceCovMatrix
        //
	// B = - A*P*SIGMA
	//
	// C = Vkm^-1 + BT*A^-1*B
	//

        //get the NuisanceCovMatrix for K = KBest
        TFile *NuisanceCovMatrix_InputFile = new TFile("./HesseMatrix/robustHesse"+cuts+Form("_iK%1.0f_iM%1.0f_TOY%i.root", KBest_index, MBest_index, NTOY));
        if(!NuisanceCovMatrix_InputFile->IsOpen()) { cout<<"Cannot find NuisanceCovMatrix input file"<<endl; return ; }
        TH2D *hCovMatrix = (TH2D*) NuisanceCovMatrix_InputFile->Get("h_covariance");
        //convert TH2D to a TMatrixD
        TMatrixD InverseA(hCovMatrix->GetNbinsX(), hCovMatrix->GetNbinsY());
        for(int i = 0; i < hCovMatrix->GetNbinsX(); i++) {
                for(int j = 0; j < hCovMatrix->GetNbinsY(); j++) {
                        InverseA[i][j] = hCovMatrix->GetBinContent(i+1, j+1);
                }
        }
        cout<<"\n##### NuisanceCovMatrix(KBest) #####"<<endl;
        InverseA.Print();

	TMatrixD A = InverseA; A.Invert();

	TMatrixD Vkm(2, 2);
	Vkm[0][0] = DKBest*DKBest;        Vkm[0][1] = DKBest*DMBest*corrKM;
	Vkm[1][0] = DKBest*DMBest*corrKM; Vkm[1][1] = DMBest*DMBest;

	TMatrixD InverseVkm = Vkm; InverseVkm.Invert();

	TMatrixD SIGMA(2, 2);
	SIGMA[0][0] = DKBest; SIGMA[0][1] = 0;
	SIGMA[1][0] = 0;      SIGMA[1][1] = DMBest;

	TMatrixD P(hCovMatrix->GetNbinsX(), 1);
	//sigmaK
	P[0][0] = - (NormNuisance_1sigma_K - NormNuisance);
	//P[1][0] = - (MultipleScatteringNuisance_1sigma_K - MultipleScatteringNuisance);
//	P[2][0] = - (SingleHitResNuisance_1sigma_K - SingleHitResNuisance);
	//sigmaM
	P[0][1] = - (NormNuisance_1sigma_M - NormNuisance);
	//P[1][1] = - (MultipleScatteringNuisance_1sigma_M - MultipleScatteringNuisance);
//	P[2][1] = - (SingleHitResNuisance_1sigma_M - SingleHitResNuisance);

        //Step 1: find B
	TMatrixD B = A*P*SIGMA;

        //Transpose B
        TMatrixD BT(2, hCovMatrix->GetNbinsX());
	BT.Transpose(B);

        //Step 2: find C
        TMatrixD C = InverseVkm + BT*InverseA*B;

        //Step 3: print the full covariance matrix (...starting from the hessian H, which is its inverse...)
        TMatrixD Hessian(hCovMatrix->GetNbinsX()+2, hCovMatrix->GetNbinsX()+2);
        for(int i = 0; i < hCovMatrix->GetNbinsX(); i++) {
                for(int j = 0; j < hCovMatrix->GetNbinsX(); j++) {
                        Hessian[i][j] = A[i][j];
                }
                Hessian[i][hCovMatrix->GetNbinsX()]   = B[i][0];
                Hessian[i][hCovMatrix->GetNbinsX()+1] = B[i][1];

                Hessian[hCovMatrix->GetNbinsX()][i]   = BT[0][i];
                Hessian[hCovMatrix->GetNbinsX()+1][i] = BT[1][i];
        }
	Hessian[hCovMatrix->GetNbinsX()][hCovMatrix->GetNbinsX()]     = C[0][0];
	Hessian[hCovMatrix->GetNbinsX()][hCovMatrix->GetNbinsX()+1]   = C[0][1];
	Hessian[hCovMatrix->GetNbinsX()+1][hCovMatrix->GetNbinsX()]   = C[1][0];
	Hessian[hCovMatrix->GetNbinsX()+1][hCovMatrix->GetNbinsX()+1] = C[1][1];

	TMatrixD FullCovMatrix = Hessian; FullCovMatrix.Invert();
        cout<<"\n#### Final Full Covariance Matrix: #####"<<endl;
        FullCovMatrix.Print();

	//get all the useful quantities
        DKBest   = TMath::Sqrt(FullCovMatrix[hCovMatrix->GetNbinsX()][hCovMatrix->GetNbinsX()]);
        DMBest   = TMath::Sqrt(FullCovMatrix[hCovMatrix->GetNbinsX()+1][hCovMatrix->GetNbinsX()+1]);
	corrKM   = FullCovMatrix[hCovMatrix->GetNbinsX()+1][hCovMatrix->GetNbinsX()]/(DKBest*DMBest);

	Double_t DNormNuisance = TMath::Sqrt(FullCovMatrix[0][0]);
        Double_t corr_K_NormNuisance = FullCovMatrix[0][hCovMatrix->GetNbinsX()]/(DKBest*DNormNuisance);
        Double_t corr_M_NormNuisance = FullCovMatrix[0][hCovMatrix->GetNbinsX()+1]/(DMBest*DNormNuisance);

	Double_t DMultipleScatteringNuisance = TMath::Sqrt(FullCovMatrix[1][1]);
        Double_t corr_K_MultipleScatteringNuisance = FullCovMatrix[1][hCovMatrix->GetNbinsX()]/(DKBest*DMultipleScatteringNuisance);
        Double_t corr_M_MultipleScatteringNuisance = FullCovMatrix[1][hCovMatrix->GetNbinsX()+1]/(DMBest*DMultipleScatteringNuisance);

//	Double_t DSingleHitResNuisance = TMath::Sqrt(FullCovMatrix[2][2]);
//	Double_t corr_K_SingleHitResNuisance = FullCovMatrix[2][hCovMatrix->GetNbinsX()]/(DKBest*DSingleHitResNuisance);
//	Double_t corr_M_SingleHitResNuisance = FullCovMatrix[2][hCovMatrix->GetNbinsX()+1]/(DMBest*DSingleHitResNuisance);

	Double_t corr_NormNuisance_MultipleScatteringNuisance = FullCovMatrix[0][1]/(DNormNuisance*DMultipleScatteringNuisance);
//	Double_t corr_NormNuisance_SingleHitResNuisance       = FullCovMatrix[0][2]/(DNormNuisance*DSingleHitResNuisance);
//	Double_t corr_MultipleScatteringNuisance_SingleHitResNuisance = FullCovMatrix[1][2]/(DMultipleScatteringNuisance*DSingleHitResNuisance);

	cout<<"\n****** FIT RESULTS ******"<<endl;
	cout<<"K = "<<KBest<<" +/- "<<DKBest<<endl;
	cout<<"M = "<<MBest<<" +/- "<<DMBest<<endl;
	cout<<"corrKM = "<<corrKM<<endl;
	cout<<"NuNorm = "<<NormNuisance<<" +/- "<<DNormNuisance<<endl;
	cout<<"MuMS   = "<<MultipleScatteringNuisance<<" +/- "<<DMultipleScatteringNuisance<<endl;
//	cout<<"MuIntr = "<<SingleHitResNuisance<<" +/- "<<DSingleHitResNuisance<<endl;
	cout<<"\ncorr(K, NuNorm) = "<<corr_K_NormNuisance<<endl;
	cout<<  "corr(K, MuMS)   = "<<corr_K_MultipleScatteringNuisance<<endl;
//	cout<<  "corr(K, MuIntr) = "<<corr_K_SingleHitResNuisance<<endl;
	cout<<  "corr(M, NuNorm) = "<<corr_M_NormNuisance<<endl;
	cout<<  "corr(M, MuMS)   = "<<corr_M_MultipleScatteringNuisance<<endl;
//	cout<<  "corr(M, MuIntr) = "<<corr_M_SingleHitResNuisance<<endl;
	cout<<  "corr(NuNorm, MuMS)   = "<<corr_NormNuisance_MultipleScatteringNuisance<<endl;
//	cout<<  "corr(NuNorm, MuIntr) = "<<corr_NormNuisance_SingleHitResNuisance<<endl;
//	cout<<  "corr(MuMS, MuIntr)   = "<<corr_MultipleScatteringNuisance_SingleHitResNuisance<<endl;



	TString superFileName = "resultsFit2D_" + cuts + ".root";
	TFile *superFile = new TFile("fitParameters_" + cuts + ".root", "UPDATE");

	FullCovMatrix.Write("FullCovMatrix", TObject::kOverwrite);	

	h_NormNuisance->Write("", TObject::kOverwrite);
	for(int iM = 0; iM < ngrid; iM++) NormNuisanceM[iM]->Write("", TObject::kOverwrite);
	cM_NormNuisance->Write("", TObject::kOverwrite);
	for(int iK = 0; iK < ngrid; iK++) NormNuisanceK[iK]->Write("", TObject::kOverwrite);
	cK_NormNuisance->Write("", TObject::kOverwrite);

	h_MultipleScatteringNuisance->Write("", TObject::kOverwrite);
	for(int iM = 0; iM < ngrid; iM++) MSNuisanceM[iM]->Write("", TObject::kOverwrite);
	cM_MSNuisance->Write("", TObject::kOverwrite);
	for(int iK = 0; iK < ngrid; iK++) MSNuisanceK[iK]->Write("", TObject::kOverwrite);
	cK_MSNuisance->Write("", TObject::kOverwrite);

//		g_SingleHitResNuisance->Write("", TObject::kOverwrite);
//		c_SingleHitResNuisance->Write("", TObject::kOverwrite);

	for(int iM = 0; iM < ngrid; iM++) resultsM[iM]->Write("", TObject::kOverwrite);
	cM->Write("", TObject::kOverwrite);
	corr_KM->Write("", TObject::kOverwrite);
	cKM->Write("", TObject::kOverwrite);
	resultsK->Write("", TObject::kOverwrite);
	cK->Write("", TObject::kOverwrite);
	h_nLL2D->Write("", TObject::kOverwrite);
	c_nLL2D->Write("", TObject::kOverwrite);


	TTree *results = new TTree("results", "tree with full fit results");
	results->Branch("K", &KBest, "K/D");
	results->Branch("DK", &DKBest, "DK/D");
	results->Branch("M", &MBest, "M/D");
	results->Branch("DM", &DMBest, "DM/D");
	results->Branch("corrKM", &corrKM, "corrKM/D");
	results->Branch("NormNuisance", &NormNuisance, "NormNuisance/D");
	results->Branch("DNormNuisance", &DNormNuisance, "DNormNuisance/D");
	results->Branch("MSNuisance", &MultipleScatteringNuisance, "MSNuisance/D");
	results->Branch("DMSNuisance", &DMultipleScatteringNuisance, "DMSNuisance/D");
//	results->Branch("IntrResNuisance", &SingleHitResNuisance, "IntrResNuisance/D");
//	results->Branch("DIntrResNuisance", &DSingleHitResNuisance, "DIntrResNuisance/D");
	results->Branch("corrKNormNuisance", &corr_K_NormNuisance, "corrKNormNuisance/D");
	results->Branch("corrKMSNuisance", &corr_K_MultipleScatteringNuisance, "corrKMSNuisance/D");
//	results->Branch("corrKIntrResNuisance", &corr_K_SingleHitResNuisance, "corrKIntrResNuisance/D");
	results->Branch("corrMNormNuisance", &corr_M_NormNuisance, "corrMNormNuisance/D");
	results->Branch("corrMMSNuisance", &corr_M_MultipleScatteringNuisance, "corrMMSNuisance/D");
//	results->Branch("corrMIntrResNuisance", &corr_M_SingleHitResNuisance, "corrMIntrResNuisance/D");
	results->Branch("corrNormMS", &corr_NormNuisance_MultipleScatteringNuisance, "corrNormMS/D");
//	results->Branch("corrNormIntrRes", &corr_NormNuisance_SingleHitResNuisance, "corrNormIntrRes/D");
//	results->Branch("corrMSIntrRes", &corr_MultipleScatteringNuisance_SingleHitResNuisance, "corrMSIntrRes/D");
	results->Branch("amuhlo100", &amuhlo_100, "amuhlo100/D");
	results->Branch("Damuhlo100", &Damuhlo_100, "Damuhlo100/D");
	results->Branch("amuhlo87", &amuhlo_87, "amuhlo87/D");
	results->Branch("Damuhlo87", &Damuhlo_87, "Damuhlo87/D");
	results->Branch("amuhlo13", &amuhlo_13, "amuhlo13/D");
	results->Branch("Damuhlo13", &Damuhlo_13, "Damuhlo13/D");
	results->Fill();
	results->Write("", TObject::kOverwrite);


return ;
}





void getResults(int NTOY, TString cuts, Double_t *Kbest_fit, Double_t *DKbest_fit, Double_t *Mbest_fit, Double_t *DMbest_fit, Double_t *CorrKM, Double_t *DCorrKM) {


        //TH2 likelihood grid
        h_nLL2D = new TH2D("h_nLL2D", "likelihood distribution for the (K, M) template grid; #frac{K_{template} - K_{ref}}{#sigma_{K}}; #frac{M_{template} - M_{ref}}{#sigma_{M}}",
                            ngrid, -sigmaLim, sigmaLim, ngrid, -sigmaLim, sigmaLim);

	cM = new TCanvas("cM", "Likelihood as a function of M, K fixed", 1080*1.3, 720*1.3);
	cM->Divide(8, 5);


	//TH2 containing the best fit value of the normalization nuisance
	h_NormNuisance = new TH2D("h_NormNuisance", "#nu normalization nuisance parameter grid, for K and M templates; #frac{K_{template} - K_{ref}}{#sigma_{K}}; #frac{M_{template} - M_{ref}}{#sigma_{M}}", 
				  ngrid, -sigmaLim, sigmaLim, ngrid, -sigmaLim, sigmaLim);
        cM_NormNuisance = new TCanvas("cM_NormNuisance", "#nu nuisance as a function of M, K fixed", 1080*1.3, 720*1.3);
        cM_NormNuisance->Divide(8, 5);
        cK_NormNuisance = new TCanvas("cK_NormNuisance", "#nu nuisance as a function of K, M fixed", 1080*1.3, 720*1.3);
        cK_NormNuisance->Divide(8, 5);

	//TH2 containing the best fit value of the shape MS nuisance
	h_MultipleScatteringNuisance = new TH2D("h_MultipleScatteringNuisance", "Multiple scattering nuisance parameter grid, for K and M templates; #frac{K_{template} - K_{ref}}{#sigma_{K}}; #frac{M_{template} - M_{ref}}{#sigma_{M}}",
						ngrid, -sigmaLim, sigmaLim, ngrid, -sigmaLim, sigmaLim);
        cM_MSNuisance = new TCanvas("cM_MSNuisance", "#theta_{MS} nuisance as a function of M, K fixed", 1080*1.3, 720*1.3);
        cM_MSNuisance->Divide(8, 5);
        cK_MSNuisance = new TCanvas("cK_MSNuisance", "#theta_{MS} nuisance as a function of K, M fixed", 1080*1.3, 720*1.3);
        cK_MSNuisance->Divide(8, 5);

	std::map<int, std::pair<double, double>> Mbest_giveniK;//key: iK. value: first = Mfit giveniK, second = errorMfit_giveniK
	std::map<double, int> iK_givenNLL;//key: nLL. value: iK. in questo modo ordino sui nLL e seleziono solo i punti attorno al minimo
        TFile *results_file;
	double nLL_min0 = 0;
	const int MINPONTS_M = 5;
	for(int iK = 0; iK < ngrid; iK++) {

		int nuisance_counter = 0;
		NormNuisanceM[iK] = new TGraphErrors();
		NormNuisanceM[iK]->SetTitle(Form("iK = %i; M [#sigma_{M} units]; #nu", iK));
		NormNuisanceM[iK]->SetMarkerStyle(8);
		MSNuisanceM[iK] = new TGraphErrors();
		MSNuisanceM[iK]->SetTitle(Form("iK = %i; M [#sigma_{M} units]; #theta_{MS}", iK));
		MSNuisanceM[iK]->SetMarkerStyle(8);

		std::map<double, double> points_MnLL_giveniK;//key: nLL, value: M_template
		for(int iM = 0; iM < ngrid; iM++) {
			TString file_dir = "./OutputCombine/higgsCombine2D_" + cuts + Form("_iK%i_iM%i_TOY%i.root", iK, iM, NTOY);
			results_file = new TFile(file_dir);
			if(results_file->IsOpen()) {
				TTree *results_tree = (TTree*) results_file->Get("limit");
				if(results_tree != nullptr) {
					Double_t nLL;
					Double_t nLL0;
					results_tree->SetBranchAddress("nll", &nLL);
					results_tree->SetBranchAddress("nll0", &nLL0);

					Float_t LumiNuisance;
					Float_t LumiNuisance_error;
					results_tree->SetBranchAddress("trackedParam_provalnN_error", &LumiNuisance);
					results_tree->SetBranchAddress("trackedError_provalnN_error", &LumiNuisance_error);

					Float_t MS;
					Float_t MS_error;
					results_tree->SetBranchAddress("trackedParam_MultipleScattering", &MS);
					results_tree->SetBranchAddress("trackedError_MultipleScattering", &MS_error);

					results_tree->GetEntry(0);

					Double_t sumNLL = 2*static_cast<double>(nLL + nLL0);
					h_nLL2D->SetBinContent(iK+1, iM+1, sumNLL);
					points_MnLL_giveniK[sumNLL] = - sigmaLim + static_cast<double>(iM)/sigmaStep;//points in dMu units within +/- sigmaLim

					h_NormNuisance->SetBinContent(iK+1, iM+1, LumiNuisance);//lnNuisance);
					NormNuisanceM[iK]->SetPoint(nuisance_counter, - sigmaLim + static_cast<double>(iM)/sigmaStep, LumiNuisance);
					NormNuisanceM[iK]->SetPointError(nuisance_counter, 0, LumiNuisance_error);

					h_MultipleScatteringNuisance->SetBinContent(iK+1, iM+1, MS);
					MSNuisanceM[iK]->SetPoint(nuisance_counter, - sigmaLim + static_cast<double>(iM)/sigmaStep, MS);
					MSNuisanceM[iK]->SetPointError(nuisance_counter, 0,  MS_error);

					nuisance_counter++;
					results_file->Close();
				}
			}
                }

		cM_NormNuisance->cd(iK+1);
		NormNuisanceM[iK]->Draw("APE");
		cM_MSNuisance->cd(iK+1);
		MSNuisanceM[iK]->Draw("APE");

		int counter = 0;
		resultsM[iK] = new TGraphErrors();
		resultsM[iK]->SetName(Form("resultsM_K%i", iK));
                resultsM[iK]->SetTitle(Form("iK = %i; M [#sigma_{M} units]", iK));
		resultsM[iK]->SetMarkerStyle(8);
		//fill resultsM with points (M, nLL) around the minimum nLL and within nLL = nLLmin+1
		double nLL_relativeMinimum = points_MnLL_giveniK.begin()->first;
//		if(iK == 0) nLL_min0 = nLL_relativeMinimum;
                nLL_min0 += nLL_relativeMinimum;
		for(auto& element : points_MnLL_giveniK) {
			if(element.first <= nLL_relativeMinimum+1) {
				resultsM[iK]->SetPoint(counter, element.second, element.first);
				//resultsM[iK]->SetPointError(counter, 0, element.first*1e-6);
				counter++;
			}
		}
		if(resultsM[iK]->GetN() < MINPONTS_M) {
			counter = 0;
			auto it = points_MnLL_giveniK.begin();
			while(counter <= MINPONTS_M-1) {
				resultsM[iK]->SetPoint(counter, it->second, it->first);
				counter++;
				++it;
			}
		}
	}


	for(int iM = 0; iM < ngrid; iM++) {
                int nuisance_counter = 0;
                NormNuisanceK[iM] = new TGraphErrors();
                NormNuisanceK[iM]->SetTitle(Form("iM = %i; K [#sigma_{K} units]; #nu", iM));
                NormNuisanceK[iM]->SetMarkerStyle(8);
                MSNuisanceK[iM] = new TGraphErrors();
                MSNuisanceK[iM]->SetTitle(Form("iM = %i; K [#sigma_{K} units]; #theta_{MS}", iM));
                MSNuisanceK[iM]->SetMarkerStyle(8);
                for(int iK = 0; iK < ngrid; iK++) {
			TString file_dir = "./OutputCombine/higgsCombine2D_" + cuts + Form("_iK%i_iM%i_TOY%i.root", iK, iM, NTOY);
			results_file = new TFile(file_dir);
                        if(results_file->IsOpen()) {
                        	TTree *results_tree = (TTree*) results_file->Get("limit");
				if(results_tree != nullptr) {
					Float_t LumiNuisance;
					Float_t LumiNuisance_error;
					results_tree->SetBranchAddress("trackedParam_provalnN_error", &LumiNuisance);
					results_tree->SetBranchAddress("trackedError_provalnN_error", &LumiNuisance_error);

					Float_t MS;
					Float_t MS_error;
					results_tree->SetBranchAddress("trackedParam_MultipleScattering", &MS);
					results_tree->SetBranchAddress("trackedError_MultipleScattering", &MS_error);

					results_tree->GetEntry(0);

					NormNuisanceK[iM]->SetPoint(nuisance_counter, - sigmaLim + static_cast<double>(iK)/sigmaStep, LumiNuisance);
					NormNuisanceK[iM]->SetPointError(nuisance_counter, 0, LumiNuisance_error);
					MSNuisanceK[iM]->SetPoint(nuisance_counter, - sigmaLim + static_cast<double>(iK)/sigmaStep, MS);
					MSNuisanceK[iM]->SetPointError(nuisance_counter, 0, MS_error);

					nuisance_counter++;
					results_file->Close();
				}
			}
		}
		cK_NormNuisance->cd(iM+1);
		NormNuisanceK[iM]->Draw("APE");
		cK_MSNuisance->cd(iM+1);
		MSNuisanceK[iM]->Draw("APE");
	}



	nLL_min0 = nLL_min0/ngrid;

	for(int iK = 0; iK < ngrid; iK++) {
                for(int iM = 0; iM < resultsM[iK]->GetN(); iM++) {
                        resultsM[iK]->SetPoint(iM, resultsM[iK]->GetX()[iM], resultsM[iK]->GetY()[iM] - nLL_min0);
                }

                for(int iM = 0; iM < ngrid; iM++) h_nLL2D->SetBinContent(iK+1, iM+1, h_nLL2D->GetBinContent(iK+1, iM+1) - nLL_min0);

		fitfunction[iK] = new TF1(Form("fitfunction_K%i", iK), "pol2", -5, 5);
		fitfunction[iK]->SetNpx(1e3);
		cM->cd(iK+1);
		resultsM[iK]->Fit(fitfunction[iK], "QS");
		resultsM[iK]->Draw("APE");
		fitfunction[iK]->Draw("SAME");
		Double_t nLLMin = fitfunction[iK]->GetMinimum();
		Double_t Mfit   = fitfunction[iK]->GetMinimumX();//in unit of dMu
		Double_t ErrorMfit = fitfunction[iK]->GetX(nLLMin+1, Mfit, Mfit+3, 1e-12, 1e6);//in unit of dMu
		Mfit = Mref + Mfit*dMu;//in GeV2
		ErrorMfit = Mref + ErrorMfit*dMu - Mfit;
		std::pair<double, double> newpair = {Mfit, ErrorMfit};
		Mbest_giveniK[iK] = newpair;
		iK_givenNLL[nLLMin] = iK;
        }

	resultsK = new TGraphErrors();
	resultsK->SetName("resultsK");
	resultsK->SetTitle("K template fit; K_{template}; -2lnL");
	resultsK->SetMarkerStyle(8);

	int counter = 0;
	double nLLAbsMinimum = iK_givenNLL.begin()->first;
	for(auto& element : iK_givenNLL) {
		if(element.first <= nLLAbsMinimum+1) {
			double K = - sigmaLim + static_cast<double>(element.second)/sigmaStep;//in units of dKu, within +/- sigmaLim
			resultsK->SetPoint(counter, K, element.first);
			counter++;
		}
	}
	if(resultsK->GetN() < MINPONTS_M) {
		counter = 0;
		auto it = iK_givenNLL.begin();
		while(counter <= MINPONTS_M-1) {
			double K = - sigmaLim + static_cast<double>(it->second)/sigmaStep;
			resultsK->SetPoint(counter, K, it->first);
			counter++;
			++it;
		}
	}

        c_nLL2D = new TCanvas("c_nLL2D", "", 1080*1.3, 720*1.3);
        h_nLL2D->Draw("ZCOL");

	//find KBest, DKBest
	cK = new TCanvas("cK", "Likelihood as a function of K", 1080, 720);
	TF1 *funcK = new TF1("funcK", "pol2", -5, 5);
	resultsK->Fit("funcK");
	resultsK->Draw("APE");
	funcK->Draw("SAME");
	nLLAbsMinimum = funcK->GetMinimum();
	*Kbest_fit  = funcK->GetMinimumX();//in units of dKu
	*DKbest_fit = funcK->GetX(nLLAbsMinimum+1, *Kbest_fit, *Kbest_fit+3, 1e-12, 1e6);//in units of dKu
	*Kbest_fit  = Kref + (*Kbest_fit)*dKu;
	*DKbest_fit = Kref + (*DKbest_fit)*dKu - (*Kbest_fit);

	//find MBest, DMBest, corrKM
	corr_KM = new TGraphErrors();
	corr_KM->SetName("corr_KM");
	corr_KM->SetTitle("correlation between (K_{template}, M_{fit}(K_{template})); K_{template}; M_{fit}");
	corr_KM->SetMarkerStyle(8);
	for(int k = 0; k < resultsK->GetN(); k++) {
		Double_t k_template = resultsK->GetX()[k];
		k_template = Kref + k_template*dKu;
		Double_t nLL_ktemplate = resultsK->GetY()[k];
		corr_KM->SetPoint(k, k_template, Mbest_giveniK[iK_givenNLL[nLL_ktemplate]].first);
	}
	Double_t sum = 0;
	for(auto& p : Mbest_giveniK) sum += p.second.second;
	*DMbest_fit = sum/Mbest_giveniK.size();

	cKM = new TCanvas("cKM", "", 1080, 720);
	TF1 *func_corrKM = new TF1("func_corrKM", "pol1", Kref - 5*dKu, Kref + 5*dKu);
	corr_KM->Fit("func_corrKM");
	corr_KM->Draw("APE");
	func_corrKM->Draw("SAME");

	*Mbest_fit = func_corrKM->Eval(*Kbest_fit);

	Double_t b  = func_corrKM->GetParameter(1);
	Double_t Db = func_corrKM->GetParError(1);
	*CorrKM  = b*(*DKbest_fit)/(*DMbest_fit)/TMath::Sqrt(1 + TMath::Power(b*(*DKbest_fit)/(*DMbest_fit), 2));
	*DCorrKM = Db/b*TMath::Power((*DMbest_fit)/(*DKbest_fit)/b, 2)/TMath::Power(TMath::Power((*DMbest_fit)/(*DKbest_fit)/b, 2) + 1, 3./2.);

	delete results_file;
	delete funcK;
	delete func_corrKM;

return ;
}




