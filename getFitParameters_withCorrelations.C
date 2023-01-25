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


#define SAVE 1


//this script makes the parabolic interpolation to get the best fit value of K and interpolates the nuisance parameters values from the fits performed by combine
//to calculate the nuisance parameters values corresponding to the best fit value of K.
//all the uncertainties are calculated considering the covariance matrix computed by combine.
//this script also calculates the final covariance matrix between K and the nuisance parameters.


void getResults(int NTOY, TString cuts, Double_t *Kbest_fit, Double_t *DKbest_fit);

TGraph *g_nLL; TCanvas *c_nLL;
TGraphErrors *g_NormNuisance; TCanvas *c_NormNuisance;
TGraphErrors *g_MultipleScatteringNuisance; TCanvas *c_MultipleScatteringNuisance;
//TGraphErrors *g_SingleHitResNuisance; TCanvas *c_SingleHitResNuisance;
//TGraphErrors *g_EbeamNuisance; TCanvas *c_EbeamNuisance;

//Reference values for the input parameters
//Double_t Kref   = 0.13726; Double_t dKu    = 0.045;
Double_t Kref   = 0.13724; Double_t dKu    = 0.045;
Double_t Mref   = 0.0525;  Double_t dMu    = 0.;
//range and step for the parameters grid (they depend on the specific job settings)
Int_t sigmaLim  = 5, sigmaStep = 4;
Int_t ngrid = sigmaLim*sigmaStep*2 + 1;


//Double_t starting_Ebeam = 149.950, ending_Ebeam = 150.050, step_Ebeam = 0.005;
//Int_t nSteps = (ending_Ebeam - starting_Ebeam)/step_Ebeam + 1;
//Double_t offset = 150;

Double_t oneSigma_MS    =  1.0;//************************** CHANGE THIS VALUE ACCORDING TO THE +/-1sigma SHIFT USED FOR THE MULTIPLE SCATTERING MODELIZATION
Double_t offset_MS      =  0.0;//0.6;//************************** CHANGE THIS VALUE ACCORDING TO THE NOMINAL MODELIZATION OF MULTIPLE SCATTERING EFFECTS
//Double_t oneSigma_Intr  =  1.0;
//Double_t oneSigma_Ebeam =  20.;
//Double_t offset_Intr    =  5.0;
//Double_t offset_Ebeam   =  0.0;




void getFitParameters_withCorrelations(int NTOY = 1, TString cuts = "thmu0.2_the32") {

	Double_t KBest = 0, DKBest = 0;
	//perform the parabolic interpolation to calculate the best fit value of K and interpolate the nuisance parameters to get their value for K = K_best_fit
	getResults(NTOY, cuts, &KBest, &DKBest);
	
	//find the template which is closest to the best fit value of K
	double KBest_index = -1;
	double fractional_part = modf( (KBest - g_nLL->GetX()[0])/(dKu/sigmaStep), &KBest_index);
	if(fractional_part > 0.5) KBest_index++;
	if(KBest_index < 0) KBest_index = 0;
	if(KBest_index >= ngrid) KBest_index = ngrid-1;
	

        TF1 *f_NormNuisance = new TF1("f_Norm", "pol1", Kref - 10*dKu, Kref + 10*dKu);
        f_NormNuisance->SetNpx(1e3);
        f_NormNuisance->SetParameters(1, 1);
	g_NormNuisance->Fit(f_NormNuisance, "Q0");
	
	TF1 *f_MultipleScatteringNuisance = new TF1("f_MS", "pol1", Kref - 10*dKu, Kref + 10*dKu);
	f_MultipleScatteringNuisance->SetNpx(1e3);
	f_MultipleScatteringNuisance->SetParameters(1, 1);
	g_MultipleScatteringNuisance->Fit(f_MultipleScatteringNuisance, "Q0");

/*        TF1 *f_SingleHitResNuisance = new TF1("f_1HitRes", "pol1", Kref - 10*dKu, Kref + 10*dKu);
        f_SingleHitResNuisance->SetNpx(1e3);
        f_SingleHitResNuisance->SetParameters(1, 1);
	g_SingleHitResNuisance->Fit(f_SingleHitResNuisance, "Q0");

        TF1 *f_EbeamNuisance = new TF1("f_Ebeam", "pol1", Kref - 10*dKu, Kref + 10*dKu);
        f_EbeamNuisance->SetNpx(1e3);
        f_EbeamNuisance->SetParameters(1, 1);
	g_EbeamNuisance->Fit(f_EbeamNuisance, "Q0");*/


        //get the nuisance covariance matrix for K = K_best_fit
        TFile *NuisanceCovMatrix_InputFile = new TFile("./HesseMatrix/robustHesse"+cuts+Form("_iK%1.0f_TOY%i.root", KBest_index, NTOY));
        if(!NuisanceCovMatrix_InputFile->IsOpen()) { cout<<"Cannot find NuisanceCovMatrix input file"<<endl; return ; }
        TH2D *hCovMatrix = (TH2D*) NuisanceCovMatrix_InputFile->Get("h_covariance");
        //convert TH2D into a TMatrixD
        TMatrixD NuisanceCovMatrix(hCovMatrix->GetNbinsX(), hCovMatrix->GetNbinsY());
        for(int i = 0; i < hCovMatrix->GetNbinsX(); i++) {
                for(int j = 0; j < hCovMatrix->GetNbinsY(); j++) {
                        NuisanceCovMatrix[i][j] = hCovMatrix->GetBinContent(i+1, j+1);
                }
        }
        cout<<"\n##### NuisanceCovMatrix(K_best_fit) #####"<<endl;
        NuisanceCovMatrix.Print();

        //build the full covariance matrix V:
        //
        //        [ A   b ]
        // V^-1 = [ bT  c ]
        //
        //
        //        [ Vnuisance       Vcorrelation ]
        // V    = [ VcorrelationT       Vkk      ]
        //

        //Step 1: find b
        //NuisanceCovMatrix*b = p/sigmaK
	TVectorD p(hCovMatrix->GetNbinsX());
	p[0] = - (f_NormNuisance->Eval(KBest + DKBest) - f_NormNuisance->Eval(KBest))/DKBest;
	p[1] = - (f_MultipleScatteringNuisance->Eval(KBest + DKBest) - f_MultipleScatteringNuisance->Eval(KBest))/DKBest;
//	p[2] = - (f_SingleHitResNuisance->Eval(KBest + DKBest) - f_SingleHitResNuisance->Eval(KBest))/DKBest;
//	p[3] = - (f_EbeamNuisance->Eval(KBest + DKBest) - f_EbeamNuisance->Eval(KBest))/DKBest;

	TVectorD b = NormalEqn(NuisanceCovMatrix, p);

        //Step 2: find c
        //Transpose b (still don't know how to do that with root...)
        TMatrixD bT(1, hCovMatrix->GetNbinsX());
        for(int i = 0; i < hCovMatrix->GetNbinsX(); i++) bT[0][i] = b[i];
        Double_t c = 1./DKBest/DKBest + (bT*NuisanceCovMatrix*b)[0];

        //Step 3: print the full covariance matrix (...starting from its inverse...)
        TMatrixD InverseNuisanceCovMatrix = NuisanceCovMatrix;
        InverseNuisanceCovMatrix.Invert();
        TMatrixD FullCovMatrix(hCovMatrix->GetNbinsX()+1, hCovMatrix->GetNbinsX()+1);
        for(int i = 0; i < hCovMatrix->GetNbinsX(); i++) {
                for(int j = 0; j < hCovMatrix->GetNbinsX(); j++) {
                        FullCovMatrix[i][j] = InverseNuisanceCovMatrix[i][j];
                }
                FullCovMatrix[i][hCovMatrix->GetNbinsX()] = b[i];
                FullCovMatrix[hCovMatrix->GetNbinsX()][i] = b[i];
        }
        FullCovMatrix[hCovMatrix->GetNbinsX()][hCovMatrix->GetNbinsX()] = c;
        FullCovMatrix.Invert();
        cout<<"\n#### Final Full Covariance Matrix: #####"<<endl;
        FullCovMatrix.Print();

        DKBest   = TMath::Sqrt(FullCovMatrix[hCovMatrix->GetNbinsX()][hCovMatrix->GetNbinsX()]);

	Double_t  NormNuisance = f_NormNuisance->Eval(KBest);
	Double_t DNormNuisance = TMath::Sqrt(FullCovMatrix[0][0]);
	Double_t  MultipleScatteringNuisance = f_MultipleScatteringNuisance->Eval(KBest);
	Double_t DMultipleScatteringNuisance = TMath::Sqrt(FullCovMatrix[1][1]);
/*	Double_t  SingleHitResNuisance = f_SingleHitResNuisance->Eval(KBest);
	Double_t DSingleHitResNuisance = TMath::Sqrt(FullCovMatrix[2][2]);
	Double_t  EbeamNuisance = f_EbeamNuisance->Eval(KBest);
	Double_t DEbeamNuisance = TMath::Sqrt(FullCovMatrix[3][3]);*/
        
        Double_t corr_K_NU = FullCovMatrix[0][hCovMatrix->GetNbinsX()]/(DKBest*DNormNuisance);
	Double_t corr_K_MultipleScatteringNuisance = FullCovMatrix[1][hCovMatrix->GetNbinsX()]/(DKBest*DMultipleScatteringNuisance);
/*        Double_t corr_K_SingleHitResNuisance = FullCovMatrix[1][hCovMatrix->GetNbinsX()]/(DKBest*DSingleHitResNuisance);
	Double_t corr_K_EbeamNuisance = FullCovMatrix[3][hCovMatrix->GetNbinsX()]/(DKBest*DMultipleScatteringNuisance);*/

	Double_t corr_NU_MultipleScatteringNuisance = FullCovMatrix[0][1]/(DNormNuisance*DMultipleScatteringNuisance);
//	Double_t corr_NU_SingleHitResNuisance       = FullCovMatrix[0][1]/(DNormNuisance*DSingleHitResNuisance);
//	Double_t corr_NU_EbeamNuisance = FullCovMatrix[0][3]/(DNormNuisance*DEbeamNuisance);

//	Double_t corr_MultipleScatteringNuisance_SingleHitResNuisance = FullCovMatrix[1][2]/(DMultipleScatteringNuisance*DSingleHitResNuisance);
//	Double_t corr_MultipleScatteringNuisance_EbeamNuisance = FullCovMatrix[1][3]/(DMultipleScatteringNuisance*DEbeamNuisance);

//	Double_t corr_SingleHitResNuisance_EbeamNuisance = FullCovMatrix[2][3]/(DSingleHitResNuisance*DEbeamNuisance);


	 MultipleScatteringNuisance =  MultipleScatteringNuisance*oneSigma_MS + offset_MS;
	DMultipleScatteringNuisance = DMultipleScatteringNuisance*oneSigma_MS;
/*	 SingleHitResNuisance =  SingleHitResNuisance*oneSigma_Intr + offset_Intr;
	DSingleHitResNuisance = DSingleHitResNuisance*oneSigma_Intr;
	 EbeamNuisance =  EbeamNuisance*oneSigma_Ebeam + offset_Ebeam;
	DEbeamNuisance = DEbeamNuisance*oneSigma_Ebeam;*/


	cout<<"****** FIT RESULTS ******"<<endl;
	cout<<"K = "<<KBest<<" +/- "<<DKBest<<endl;
	cout<<"NuNorm = "<<NormNuisance<<" +/- "<<DNormNuisance<<endl;
	cout<<"MuMS   = "<<MultipleScatteringNuisance<<" +/- "<<DMultipleScatteringNuisance<<endl;
//	cout<<"MuIntr = "<<SingleHitResNuisance<<" +/- "<<DSingleHitResNuisance<<endl;
//	cout<<"MuEbeam = "<<EbeamNuisance<<" +/- "<<DEbeamNuisance<<endl;
	cout<<"\ncorr(K, NuNorm) = "<<corr_K_NU<<endl;
	cout<<  "corr(K, MuMS)   = "<<corr_K_MultipleScatteringNuisance<<endl;
//	cout<<  "corr(K, MuIntr) = "<<corr_K_SingleHitResNuisance<<endl;
//	cout<<  "corr(K, MuEbeam) = "<<corr_K_EbeamNuisance<<endl;
	cout<<  "corr(NuNorm, MuMS)   = "<<corr_NU_MultipleScatteringNuisance<<endl;
//	cout<<  "corr(NuNorm, MuIntr) = "<<corr_NU_SingleHitResNuisance<<endl;
//	cout<<  "corr(NuNorm, MuEbeam) = "<<corr_NU_EbeamNuisance<<endl;
//	cout<<  "corr(MuMS, MuIntr)   = "<<corr_MultipleScatteringNuisance_SingleHitResNuisance<<endl;
//	cout<<  "corr(MuMS, MuEbeam)   = "<<corr_MultipleScatteringNuisance_EbeamNuisance<<endl;
//	cout<<  "corr(MuIntr, MuEbeam)   = "<<corr_SingleHitResNuisance_EbeamNuisance<<endl;


	
	
	
	
	TString superFileName = "resultsFit2D_" + cuts + ".root";
	TFile *superFile = new TFile(superFileName, "UPDATE");
	TH1D *h_chi2_nll            = (TH1D*) superFile->Get("h_chi2_nll");
	TH1D *h_K_2D                = (TH1D*) superFile->Get("h_K_2D");
	TH1D *h_K_errors            = (TH1D*) superFile->Get("h_K_errors");
	TH1D *h_LumiNuisance        = (TH1D*) superFile->Get("h_LumiNuisance");
	TH1D *h_LumiNuisance_errors = (TH1D*) superFile->Get("h_LumiNuisance_errors");
	TH1D *h_MS                  = (TH1D*) superFile->Get("h_MultipleScattering");
	TH1D *h_MS_errors           = (TH1D*) superFile->Get("h_MultipleScattering_errors");
/*	TH1D *h_Intr                = (TH1D*) superFile->Get("h_SingleHitRes");
	TH1D *h_Intr_errors         = (TH1D*) superFile->Get("h_SingleHitRes_errors");
	TH1D *h_Ebeam               = (TH1D*) superFile->Get("h_Ebeam");
	TH1D *h_Ebeam_errors        = (TH1D*) superFile->Get("h_Ebeam_errors");*/
	
	TH1D *h_corrK_NU            = (TH1D*) superFile->Get("h_corrK_NU");
	TH1D *h_corrK_MUMS              = (TH1D*) superFile->Get("h_corrK_MUMS");
//	TH1D *h_corrK_MUIntr        = (TH1D*) superFile->Get("h_corrK_MUIntr");
//	TH1D *h_corrK_MUEbeam           = (TH1D*) superFile->Get("h_corrK_MUEbeam");
	
	TH1D *h_corrNU_MUMS             = (TH1D*) superFile->Get("h_corrNU_MUMS");
//	TH1D *h_corrNU_MUIntr           = (TH1D*) superFile->Get("h_corrNU_MUIntr");
//	TH1D *h_corrNU_MUEbeam          = (TH1D*) superFile->Get("h_corrNU_MUEbeam");
	
//	TH1D *h_corrMUMS_MUIntr         = (TH1D*) superFile->Get("h_corrMUMS_MUIntr");
//	TH1D *h_corrMUMS_MUEbeam        = (TH1D*) superFile->Get("h_corrMUMS_MUEbeam");
	
//	TH1D *h_corrMUIntr_MUEbeam      = (TH1D*) superFile->Get("h_corrMUIntr_MUEbeam");



	if(h_chi2_nll == nullptr) h_chi2_nll = new TH1D("h_chi2_nll", "#Chi^{2} of the nLL parabolic fits; #Chi^2/ndof; Entries", 200, 0, 10);

	if(h_K_2D == nullptr) h_K_2D = new TH1D("h_K_2D", "Best fit results K from toys. fit 2D; K; Entries", 200, Kref - 10*dKu, Kref + 10*dKu);
        if(h_K_errors == nullptr) h_K_errors = new TH1D("h_K_errors", "Distribution of #DeltaK from toys. fit 2D; #DeltaK; Entries", 300, 0.001, 2*dKu);

	if(h_LumiNuisance == nullptr) h_LumiNuisance = new TH1D("h_LumiNuisance", "Distribution of #nu normalization nuisance from toys. fit 2D; #nu; Entries", 200, -2, 2);
	if(h_LumiNuisance_errors == nullptr) h_LumiNuisance_errors = new TH1D("h_LumiNuisance_errors", "Distribution of #Delta#nu from toys. fit 2D; #Delta#nu; Entries", 200, 0, 1);

	if(h_MS == nullptr) h_MS = new TH1D("h_MultipleScattering", "Distribution of #mu_{MS} shape nuisance from toys. fit 2D; #mu_{MS}; Entries", 200, -2, 2);
	if(h_MS_errors == nullptr) h_MS_errors = new TH1D("h_MultipleScattering_errors", "Distribution of #Delta#mu_{MS} from toys. fit 2D; #Delta#mu_{MS}; Entries", 200, 0, 1);

/*	if(h_Intr == nullptr) h_Intr = new TH1D("h_SingleHitRes", "Distribution of #mu_{Intr} shape nuisance from toys. fit 2D; #mu_{Intr}; Entries", 200, 3, 6);
	if(h_Intr_errors == nullptr) h_Intr_errors = new TH1D("h_SingleHitRes_errors", "Distribution of #Delta#mu_{Intr} from toys. fit 2D; #Delta#mu_{Intr}; Entries", 200, 0, 1);

	if(h_Ebeam == nullptr) h_Ebeam = new TH1D("h_Ebeam", "Distribution of #mu_{Ebeam} shape nuisance from toys. fit 2D; #mu_{Ebeam}; Entries", 200, 4, 8);
	if(h_Ebeam_errors == nullptr) h_Ebeam_errors = new TH1D("h_Ebeam_errors", "Distribution of #Delta#mu_{Ebeam} from toys. fit 2D; #Delta#mu_{Ebeam}; Entries", 200, 0, 1);
*/	
	
	if(h_corrK_NU == nullptr) h_corrK_NU = new TH1D("h_corrK_NU", "Distribution of #rho(K, #nu) from toys. fit 2D; #rho(K, #nu); Entries", 200, -1, 1);
	if(h_corrK_MUMS == nullptr) h_corrK_MUMS = new TH1D("h_corrK_MUMS", "Distribution of #rho(K, #mu_{MS}) from toys. fit 2D; #rho(K, #mu_{MS}); Entries", 200, -1, 1);
//	if(h_corrK_MUIntr == nullptr) h_corrK_MUIntr = new TH1D("h_corrK_MUIntr", "Distribution of #rho(K, #mu_{Intr}) from toys. fit 2D; #rho(K, #mu_{Intr}); Entries", 200, -1, 1);
//	if(h_corrK_MUEbeam == nullptr) h_corrK_MUEbeam = new TH1D("h_corrK_MUEbeam", "Distribution of #rho(K, #mu_{E_{beam}}) from toys. fit 2D; #rho(K, #mu_{E_{beam}}); Entries", 200, -1, 1);

	if(h_corrNU_MUMS == nullptr) h_corrNU_MUMS = new TH1D("h_corrNU_MUMS", "Distribution of #rho(#nu, #mu_{MS}) from toys. fit 2D; #rho(#nu, #mu_{MS}); Entries", 200, -1, 1);
//	if(h_corrNU_MUIntr == nullptr) h_corrNU_MUIntr = new TH1D("h_corrNU_MUIntr", "Distribution of #rho(#nu, #mu_{Intr}) from toys. fit 2D; #rho(#nu, #mu_{Intr}); Entries", 200, -1, 1);
//	if(h_corrNU_MUEbeam == nullptr) h_corrNU_MUEbeam = new TH1D("h_corrNU_MUEbeam", "Distribution of #rho(#nu, #mu_{E_{beam}}) from toys. fit 2D; #rho(#nu, #mu_{E_{beam}}); Entries", 200, -1, 1);

//	if(h_corrMUMS_MUIntr == nullptr) h_corrMUMS_MUIntr = new TH1D("h_corrMUMS_MUIntr", "Distribution of #rho(#mu_{MS}, #mu_{Intr}) from toys. fit 2D; #rho(#mu_{MS}, #mu_{Intr}); Entries", 200, -1, 1);
//	if(h_corrMUMS_MUEbeam == nullptr) h_corrMUMS_MUEbeam = new TH1D("h_corrMUMS_MUEbeam", "Distribution of #rho(#mu_{MS}, #mu_{E_{beam}}) from toys. fit 2D; #rho(#mu_{MS}, #mu_{E_{beam}}); Entries", 200, -1, 1);

//	if(h_corrMUIntr_MUEbeam == nullptr) h_corrMUIntr_MUEbeam = new TH1D("h_corrMUIntr_MUEbeam", "Distribution of #rho(#mu_{Intr}. #mu_{E_{beam}}) from toys. fit 2D; #rho(#mu_{Intr}. #mu_{E_{beam}}); Entries", 200, -1, 1);


	h_chi2_nll->Fill(g_nLL->GetFunction("f_nLL")->GetChisquare()/g_nLL->GetFunction("f_nLL")->GetNDF());

	h_K_2D->Fill(KBest);
	h_K_errors->Fill(DKBest);
	h_LumiNuisance->Fill(NormNuisance);
	h_LumiNuisance_errors->Fill(DNormNuisance);
	h_MS->Fill(MultipleScatteringNuisance);
	h_MS_errors->Fill(DMultipleScatteringNuisance);
/*	h_Intr->Fill(SingleHitResNuisance);
	h_Intr_errors->Fill(DSingleHitResNuisance);
	h_Ebeam->Fill(EbeamNuisance);
	h_Ebeam_errors->Fill(DEbeamNuisance);*/

	h_corrK_NU->Fill(corr_K_NU);
	h_corrK_MUMS->Fill(corr_K_MultipleScatteringNuisance);
//	h_corrK_MUIntr->Fill(corr_K_SingleHitResNuisance);
//	h_corrK_MUEbeam->Fill(corr_K_EbeamNuisance);

	h_corrNU_MUMS->Fill(corr_NU_MultipleScatteringNuisance);
//	h_corrNU_MUIntr->Fill(corr_NU_SingleHitResNuisance);
//	h_corrNU_MUEbeam->Fill(corr_NU_EbeamNuisance);

//	h_corrMUMS_MUIntr->Fill(corr_MultipleScatteringNuisance_SingleHitResNuisance);
//	h_corrMUMS_MUEbeam->Fill(corr_MultipleScatteringNuisance_EbeamNuisance);

//	h_corrMUIntr_MUEbeam->Fill(corr_SingleHitResNuisance_EbeamNuisance);


	h_chi2_nll->Write("", TObject::kOverwrite);

	h_K_2D->Write("", TObject::kOverwrite);
	h_K_errors->Write("", TObject::kOverwrite);
	h_LumiNuisance->Write("", TObject::kOverwrite);
	h_LumiNuisance_errors->Write("", TObject::kOverwrite);
	h_MS->Write("", TObject::kOverwrite);
	h_MS_errors->Write("", TObject::kOverwrite);
/*	h_Intr->Write("", TObject::kOverwrite);
	h_Intr_errors->Write("", TObject::kOverwrite);
	h_Ebeam->Write("", TObject::kOverwrite);
	h_Ebeam_errors->Write("", TObject::kOverwrite);*/

	h_corrK_NU->Write("", TObject::kOverwrite);
	h_corrK_MUMS->Write("", TObject::kOverwrite);
//	h_corrK_MUIntr->Write("", TObject::kOverwrite);
//	h_corrK_MUEbeam->Write("", TObject::kOverwrite);

	h_corrNU_MUMS->Write("", TObject::kOverwrite);
//	h_corrNU_MUIntr->Write("", TObject::kOverwrite);
//	h_corrNU_MUEbeam->Write("", TObject::kOverwrite);

//	h_corrMUMS_MUIntr->Write("", TObject::kOverwrite);
//	h_corrMUMS_MUEbeam->Write("", TObject::kOverwrite);

//	h_corrMUIntr_MUEbeam->Write("", TObject::kOverwrite);

	if(NTOY < 10) {
                if(!superFile->GetDirectory("plots_lnL_toys")) superFile->mkdir("plots_lnL_toys");
                superFile->cd("plots_lnL_toys");
		g_nLL->Write(Form("lnL_toy%i", NTOY), TObject::kOverwrite);

                if(!superFile->GetDirectory("plots_nuisance_toys")) superFile->mkdir("plots_nuisance_toys");
                superFile->cd("plots_nuisance_toys");
		g_NormNuisance->Write(Form("normNuisance_toy%i", NTOY), TObject::kOverwrite);
		g_MultipleScatteringNuisance->Write(Form("MultipleScattering_toy%i", NTOY), TObject::kOverwrite);
//		g_SingleHitResNuisance->Write(Form("SingleHitRes_toy%i", NTOY), TObject::kOverwrite);
//		g_EbeamNuisance->Write(Form("Ebeam_toy%i", NTOY), TObject::kOverwrite);
	}
	superFile->Close();


return ;
}






void getResults(int NTOY, TString cuts, Double_t *Kbest_fit, Double_t *DKbest_fit) {

        //TGraph containing the best fit value of the normalization nuisance
	g_NormNuisance = new TGraphErrors();
        g_NormNuisance->SetName("g_NormNuisance");
        g_NormNuisance->SetTitle("#nu normalization nuisance parameter vs K_{template}; K_{template}; #nu normalization nuisance");
        g_NormNuisance->SetMarkerStyle(8);
        g_NormNuisance->SetMarkerColor(kBlue);

        //TGraph containing the best fit value of the shape MS nuisance
	g_MultipleScatteringNuisance = new TGraphErrors();
        g_MultipleScatteringNuisance->SetName("g_MultipleScatteringNuisance");
        g_MultipleScatteringNuisance->SetTitle("Multiple scattering nuisance parameter vs K_{template}; K_{template}; Multiple Scattering nuisance");
        g_MultipleScatteringNuisance->SetMarkerStyle(8);
        g_MultipleScatteringNuisance->SetMarkerColor(kBlue);

        //TGraph containing the best fit value of the shape SingleHitRes nuisance
/*	g_SingleHitResNuisance = new TGraphErrors();
        g_SingleHitResNuisance->SetName("g_SingleHitResNuisance");
        g_SingleHitResNuisance->SetTitle("Single hit resolution nuisance parameter vs K_{template}; K_{template}; Single hit resolution nuisance");
        g_SingleHitResNuisance->SetMarkerStyle(8);
        g_SingleHitResNuisance->SetMarkerColor(kBlue);*/

        //TGraph containing the best fit value of the shape Ebeam nuisance
/*	g_EbeamNuisance = new TGraphErrors();
        g_EbeamNuisance->SetName("g_EbeamNuisance");
        g_EbeamNuisance->SetTitle("Beam energy nuisance parameter vs K_{template}; K_{template}; Ebeam nuisance");
        g_EbeamNuisance->SetMarkerStyle(8);
        g_EbeamNuisance->SetMarkerColor(kBlue);*/

        //TGraph containing the negative LogLikelihood at minimum
        g_nLL = new TGraph();
        g_nLL->SetName("g_nLL");
        g_nLL->SetTitle("2(nLL+nLL0) vs K_{template}; K; 2(nLL+nLL0)");
        g_nLL->SetMarkerStyle(8);
        g_nLL->SetMarkerColor(kBlue);

	Double_t sumNLL = 0;

        TFile *results_file;

	TString resultsFileName = "./OutputCombine/higgsCombine2D_" + cuts;
        int j = 0;
        for(int iK = 0; iK < ngrid; iK++) {
        	TString file_dir = resultsFileName + Form("_iK%i_TOY%i.root", iK, NTOY);
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

/*				Float_t SigmaIntr;
				Float_t SigmaIntr_error;
				results_tree->SetBranchAddress("trackedParam_SingleHitRes", &SigmaIntr);
				results_tree->SetBranchAddress("trackedError_SingleHitRes", &SigmaIntr_error);

				Float_t Ebeam;
				Float_t Ebeam_error;
				results_tree->SetBranchAddress("trackedParam_Ebeam", &Ebeam);
				results_tree->SetBranchAddress("trackedError_Ebeam", &Ebeam_error);*/

				results_tree->GetEntry(0);

				Double_t idK = iK - sigmaLim*sigmaStep;
				Double_t Ktemplate = Kref + idK*dKu/sigmaStep;

				sumNLL = static_cast<double>(nLL + nLL0);
				sumNLL = 2*sumNLL;

				g_nLL->SetPoint(j, Ktemplate, sumNLL);
				g_NormNuisance->SetPoint(j, Ktemplate, LumiNuisance);
				g_NormNuisance->SetPointError(j, 0, LumiNuisance_error);
				g_MultipleScatteringNuisance->SetPoint(j, Ktemplate, MS);
				g_MultipleScatteringNuisance->SetPointError(j, 0, MS_error);
/*				g_SingleHitResNuisance->SetPoint(j, Ktemplate, SigmaIntr);
				g_SingleHitResNuisance->SetPointError(j, 0, SigmaIntr_error);
				g_EbeamNuisance->SetPoint(j, Ktemplate, Ebeam);
				g_EbeamNuisance->SetPointError(j, 0, Ebeam_error);*/

				if(NTOY < 10 || NTOY%100 == 0) {
					cout<<iK<<") K = "<<Ktemplate<<" GeV. 2(nLL+nLL0) = "<<sumNLL
					    	<<"\tLumiNuisance = "<<LumiNuisance<<" +/- "<<LumiNuisance_error
					    	<<"\tMultipleScattering = "<<MS*oneSigma_MS + offset_MS<<" +/- "<<MS_error*oneSigma_MS
//					    	<<"\tSingleHitRes = "<<SigmaIntr*oneSigma_Intr + offset_Intr<<" +/- "<<SigmaIntr_error*oneSigma_Intr
//					    	<<"\tEbeam = "<<Ebeam*oneSigma_Ebeam + offset_Ebeam<<" +/- "<<Ebeam_error*oneSigma_Ebeam
					    	<<endl;
				}

				j++;
			}
		}
		results_file->Close();
	}

	Double_t min_lnL = TMath::MinElement(g_nLL->GetN(), g_nLL->GetY());
	Int_t min_index = -1;
	for(int i = 0; i < g_nLL->GetN(); i++) {
		g_nLL->SetPoint(i, g_nLL->GetX()[i], g_nLL->GetY()[i] - min_lnL);
		if(g_nLL->GetY()[i] == 0) min_index = i;
	}

	Double_t fit_minlim = 0;
	Double_t deltaNLL_fitThreshold = 60;
	for(int i = min_index; i >= 0; i--) if(g_nLL->GetY()[i] <= deltaNLL_fitThreshold) fit_minlim = g_nLL->GetX()[i]*(1-0.01);
	Double_t fit_maxlim = 0;
	for(int i = min_index; i < g_nLL->GetN(); i++) if(g_nLL->GetY()[i] <= deltaNLL_fitThreshold) fit_maxlim = g_nLL->GetX()[i]*(1+0.01);

	TF1 *f_nLL = new TF1("f_nLL", "pol2", Kref - 10*dKu, Kref + 10*dKu);
	f_nLL->SetNpx(1e3);

	c_nLL = new TCanvas("c_nLL", "", 1080, 720);
	c_nLL->SetGrid();
	c_nLL->Draw();
	g_nLL->Fit("f_nLL", "", "", fit_minlim, fit_maxlim);
	g_nLL->Draw("AP");
	f_nLL->Draw("SAME");
	g_nLL->SetMinimum(-1);

	Double_t K_best = f_nLL->GetMinimumX();
	Double_t min_nLL = f_nLL->GetMinimum();
	Double_t K_error = 1./TMath::Sqrt(f_nLL->GetParameter(2));

	TLine *l_1sigma = new TLine(g_nLL->GetXaxis()->GetXmin(), min_nLL+1, g_nLL->GetXaxis()->GetXmax(), min_nLL+1);
	l_1sigma->SetLineStyle(5);
	l_1sigma->SetLineColor(kGreen+1);
	l_1sigma->SetLineWidth(4);
	l_1sigma->Draw("SAME");

	if(NTOY < 10 || NTOY%100 == 0) cout<<NTOY<<") K_best = "<<K_best<<" +/- "<<K_error<<endl;

	*Kbest_fit = K_best;
	*DKbest_fit = K_error;

        c_NormNuisance = new TCanvas("c_NormNuisance", "", 1080, 720);
        c_NormNuisance->SetGrid();
        c_NormNuisance->Draw();
        g_NormNuisance->Draw("APE");

	c_MultipleScatteringNuisance = new TCanvas("c_MultipleScatteringNuisance", "", 1080, 720);
        c_MultipleScatteringNuisance->SetGrid();
        c_MultipleScatteringNuisance->Draw();
        g_MultipleScatteringNuisance->Draw("APE");

/*        c_SingleHitResNuisance = new TCanvas("c_SingleHitResNuisance", "", 1080, 720);
        c_SingleHitResNuisance->SetGrid();
        c_SingleHitResNuisance->Draw();
        g_SingleHitResNuisance->Draw("APE");

        c_EbeamNuisance = new TCanvas("c_EbeamNuisance", "", 1080, 720);
        c_EbeamNuisance->SetGrid();
        c_EbeamNuisance->Draw();
        g_EbeamNuisance->Draw("APE");
*/
        
        
        
        delete results_file;
	delete f_nLL;
return ;
}
