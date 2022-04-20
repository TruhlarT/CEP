 // c++ headers
#include <iostream>
#include <utility>
#include <sstream> 
#include <algorithm> 
#include <stdio.h> 
#include <stdlib.h> 
#include <vector> 
#include <fstream> 
#include <cmath> 
#include <cstdlib>
#include <unistd.h>

// ROOT headers
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include <TH2.h> 
#include <TF1.h> 
#include <TF2.h> 
#include <THStack.h> 
#include <TStyle.h> 
#include <TGraph.h> 
#include <TGraph2D.h> 
#include <TGraphErrors.h> 
#include <TCanvas.h> 
#include <TLegend.h> 
#include <TGaxis.h> 
#include <TString.h> 
#include <TColor.h> 
#include <TLine.h> 
#include <TExec.h> 
#include <TFitResultPtr.h> 
#include <TFitResult.h> 
#include <TLatex.h> 
#include <TMath.h>
#include <TLorentzVector.h>
#include "TCanvas.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TDirectory.h"


using namespace std;

const double speedOfLight = 299792458; // m/s
const double pionMass = 0.13957; // GeV /c^2
const double protonMass = 0.93827; // GeV /c^2
const double convertToDegree = 57.2957795;

enum SIDE {E = 0, East = 0, W = 1, West = 1, nSides};
enum PARTICLES {Pion = 0, Kaon = 1, Proton = 2, nParticles};
enum COMBINATIONS {ElInel = 0, El = 1, Inel = 2, nCombination};
enum RP_ID {E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D, nRomanPots};

TString particleLables[nParticles] = { TString("Pion"), TString("Kaon"), TString("Proton")};
TString combinationLabel[nCombination] = { TString("El+Inel"), TString("El"), TString("Inel")};
TString sideLabel[nSides] = { TString("East"), TString("West")};

TFile* fout;
TFile* TPCeff[6];
TFile* data;
TFile* graniitti;

TTree* tree;
TTree* treeBack;

TH1D* hInvMassGran[nCombination][nParticles];
TH1F* hInvMass[4];
TH1D* hDeltaPhi[3];
TH2D* hPxPy[3];

 // 0 = pi- 1 = K- 2 = pbar
TH3F* hTPCeff[6]; // 3 = pi+ 4 = K+ 5 = p
TH1D* hInvMassCorr[nParticles][2][3]; // 0 - signal, 1 - Background
                                    //  - El + Inel, 1 - El, 2 - Inel 
TH1D* hInvMassUncorr[nParticles][2][3]; // 0 - signal, 1 - Background

Int_t nTracks, totalCharge, nTofTrks; 
UInt_t runNumber;
Double_t VPDSumEast, VPDSumWest, VPDTimeDiff;


Double_t chiPair[nParticles]; 
Double_t invMass[nParticles];
Double_t missingPt, deltaTOF, mSquared;
Double_t deltaDeltaTOF[nParticles];
Double_t deltaTOFExpected[nParticles];

/////////////////////////////////

Bool_t elastic, fourPiState;
Bool_t trigger[17];

UInt_t BBCSmall[nSides]; // BBC truncated sum, small tiles
UInt_t BBCLarge[nSides]; // BBC truncated sum, large tiles 
Double_t xCorrelationsRp[nSides];
Double_t yCorrelationsRp[nSides];
Double_t thetaRp[nSides];
Double_t phiRp[nSides];
Double_t timeRp[nSides];
Double_t pRp[nSides];
Double_t ptRp[nSides];
Double_t etaRp[nSides];
Double_t rpX[nSides];
Double_t rpZ[nSides];
Double_t rpY[nSides];
Double_t t[nSides];
Double_t xi[nSides];

Double_t ADC[nRomanPots][2]; // RP_ID   0,    1,    2,   3,   4,   5,   6, 7
Double_t TAC[nRomanPots][2]; // RP_name E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D

Double_t dEdx[4];
Double_t momentum[4];
Double_t transMomentum[4];
Double_t TOFtime[4];
Double_t TOFlength[4];
Double_t charge[4];
Double_t nSigmaTPC[nParticles][4];
Double_t vertexId[4];
Double_t vertexesZ[4];
Double_t DcaXY[4];
Double_t DcaZ[4];
Double_t NhitsFit[4];
Double_t NhitsDEdx[4];
Double_t Eta[4];
Double_t Phi[4];
Double_t Chi2[4];

const int xMax = 20;
const int yMax = 20;
const int zMax = 24;

const double textSize = 0.04;
const double labelSize = 0.05;
const int fontStyle = 42;

void Init();


int protonsInside, protonsTotal;

void preliminaryPlots()
{
    
    TString output = "/home/truhlar/Desktop/STAR/CEP/Analysis/MissingPt/output.root";
    TString input = "/home/truhlar/Desktop/STAR/CEP/Analysis/MissingPt/ana0305.root";

    data = TFile::Open(input, "read");
    if (!data)
    {
        cout<<"Error: cannot open "<<input<<endl;
        return;
    }


    fout = new TFile(output,"RECREATE");
    Init(); // Preparing histograms 

    gStyle->SetLineWidth(2);      //axis line
    gStyle->SetFrameLineWidth(2); //frame line


    TH1I* hMissingSig = (TH1I*)data -> Get("hMissPt_All_0");
    TH1I* hMissingBcg = (TH1I*)data -> Get("hMissPt_All_1");

    TCanvas* newCanvas = new TCanvas("newCanvas","newCanvas",900,500);
    //TCanvas* newCanvas = new TCanvas("newCanvas","newCanvas",800,700);
    gPad->SetMargin(0.09,0.02,0.11,0.05); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
    gPad->SetTickx();
    gPad->SetTicky(); 
    gPad->SetLogy(0);

    TGaxis::SetMaxDigits(3);
    hMissingSig->SetTitle("; p_{T}^{miss} [GeV] ; Events / 10 MeV");
    hMissingSig->SetStats(false);
    hMissingSig->GetXaxis()->SetTitleFont(fontStyle);
    hMissingSig->GetXaxis()->SetTitleFont(fontStyle);
    hMissingSig->GetXaxis()->SetLabelFont(fontStyle);
    hMissingSig->GetYaxis()->SetLabelFont(fontStyle);
    hMissingSig->GetXaxis()->SetLabelSize(labelSize);
    hMissingSig->GetYaxis()->SetLabelSize(labelSize);
    hMissingSig->GetXaxis()->SetTitleSize(labelSize);
    hMissingSig->GetYaxis()->SetTitleSize(labelSize);
    hMissingSig->GetXaxis()->SetTitleOffset(1.0);
    hMissingSig->GetYaxis()->SetTitleOffset(0.9);
    //hMissingSig->GetYaxis()->SetRangeUser(0, 0.082);   
    hMissingSig->SetMarkerColor(1);
    hMissingSig->SetMarkerSize(1);
    hMissingSig->SetMarkerStyle(20);
    hMissingSig->SetLineColor(1);
    hMissingSig->SetLineStyle(1);
    hMissingSig->SetLineWidth(1);
    hMissingSig->Draw("E");

    hMissingBcg->SetMarkerColor(2);
    hMissingBcg->SetMarkerSize(1);
    hMissingBcg->SetMarkerStyle(22);
    hMissingBcg->SetLineColor(2);
    hMissingBcg->SetLineStyle(1);
    hMissingBcg->SetLineWidth(1);
    hMissingBcg->Draw("same E");

    TPaveText *textPub = new TPaveText(0.2,0.82,0.7,0.92,"brNDC");
    textPub -> SetTextSize(textSize+0.01);
    textPub -> SetTextAlign(22);
    textPub -> SetFillColor(0);
    textPub -> SetTextFont(42);
    textPub -> SetTextAlign(12);
    //textPub -> AddText("p + p #rightarrow p + h^{+}h^{-} + p     #sqrt{s} = 510 GeV");
    textPub -> AddText("p + p #rightarrow p + h^{+}h^{-} + p                  #sqrt{s} = 510 GeV");
    textPub -> Draw("same");

    TPaveText *text;
    text = new TPaveText(0.27,0.69,0.52,0.84,"brNDC");
    text -> SetTextSize(textSize);
    text -> SetFillColor(0);
    text -> SetTextFont(42);
    text -> SetTextAlign(32);
    text -> AddText("h^{+}, h^{-} kinematics:");   
    text -> AddText("p_{T} > 0.2 GeV");
    text -> AddText("|#eta| < 0.7");
    text -> Draw("same");

    text = new TPaveText(0.54,0.63,0.8,0.84,"brNDC");
    text -> SetTextSize(textSize);
    text -> SetFillColor(0);
    text -> SetTextFont(42);
    text -> SetTextAlign(12);
    text -> AddText("Forward proton kinematics:");
    text -> AddText("(p_{x} + 0.6)^{2} + p_{y}^{2} < 1.25 GeV^{2}");
    text -> AddText("0.4 GeV < |p_{y}| < 0.8 GeV");
    text -> AddText("p_{x} > -0.27 GeV");
    text -> Draw("same");

    TLine *line = new TLine(0.1,0,0.1,hMissingSig->GetMaximum());
    line->SetLineStyle(10);
    line->SetLineColor(1);
    line->SetLineWidth(4);
    line->Draw("same");


    TLegend *leg1 = new TLegend(0.57, 0.49, 0.9, 0.59);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(textSize);
    leg1->SetTextFont(42);
    leg1->SetMargin(0.1);
    leg1->AddEntry(hMissingSig, "Data (unlike-sign pairs)","ple");
    leg1->AddEntry(hMissingBcg, "Data (like-sign pairs)","ple");
    leg1->Draw("same");


    fout->cd();
    newCanvas->Update();
    newCanvas->Write("missingPt");
    newCanvas->Close();

    fout->Write();
    fout->Close();
    
}



void Init(){

    for (int i = 0; i < 3; ++i)
    {                          
        hInvMassCorr[0][0][i]  = new TH1D("corrInvMass" + particleLables[0] + combinationLabel[i] + "Sig", "Corrected inv. mass " + particleLables[0] + combinationLabel[i], 64, 0.3, 3.5);
        hInvMassCorr[1][0][i]  = new TH1D("corrInvMass" + particleLables[1] + combinationLabel[i] + "Sig", "Corrected inv. mass " + particleLables[1] + combinationLabel[i], 44, 0.8, 3);
        hInvMassCorr[2][0][i]  = new TH1D("corrInvMass" + particleLables[2] + combinationLabel[i] + "Sig", "Corrected inv. mass " + particleLables[2] + combinationLabel[i], 24, 1.6, 4);
        hInvMassCorr[0][1][i]  = new TH1D("corrInvMass" + particleLables[0] + combinationLabel[i] + "Bcg", "Corrected inv. mass " + particleLables[0] + combinationLabel[i], 64, 0.3, 3.5);
        hInvMassCorr[1][1][i]  = new TH1D("corrInvMass" + particleLables[1] + combinationLabel[i] + "Bcg", "Corrected inv. mass " + particleLables[1] + combinationLabel[i], 44, 0.8, 3);
        hInvMassCorr[2][1][i]  = new TH1D("corrInvMass" + particleLables[2] + combinationLabel[i] + "Bcg", "Corrected inv. mass " + particleLables[2] + combinationLabel[i], 24, 1.6, 4);

        hInvMassUncorr[0][0][i]  = new TH1D("uncorrInvMass" + particleLables[0] + combinationLabel[i] + "Sig", "Uncorrected inv. mass " + particleLables[0] + combinationLabel[i], 64, 0.3, 3.5);
        hInvMassUncorr[1][0][i]  = new TH1D("uncorrInvMass" + particleLables[1] + combinationLabel[i] + "Sig", "Uncorrected inv. mass " + particleLables[1] + combinationLabel[i], 44, 0.8, 3);
        hInvMassUncorr[2][0][i]  = new TH1D("uncorrInvMass" + particleLables[2] + combinationLabel[i] + "Sig", "Uncorrected inv. mass " + particleLables[2] + combinationLabel[i], 24, 1.6, 4);
        hInvMassUncorr[0][1][i]  = new TH1D("uncorrInvMass" + particleLables[0] + combinationLabel[i] + "Bcg", "Uncorrected inv. mass " + particleLables[0] + combinationLabel[i], 64, 0.3, 3.5);
        hInvMassUncorr[1][1][i]  = new TH1D("uncorrInvMass" + particleLables[1] + combinationLabel[i] + "Bcg", "Uncorrected inv. mass " + particleLables[1] + combinationLabel[i], 44, 0.8, 3);
        hInvMassUncorr[2][1][i]  = new TH1D("uncorrInvMass" + particleLables[2] + combinationLabel[i] + "Bcg", "Uncorrected inv. mass " + particleLables[2] + combinationLabel[i], 24, 1.6, 4);
    }

    for (int i = 0; i < nSides; ++i)
        hPxPy[i] = new TH2D("pxpy" + sideLabel[i], "px py on " + sideLabel[i], 220, -1.1, 1.1, 220, -1.1, 1.1);
    hPxPy[2] = new TH2D("pxpyEast+West", "px py on East+West", 220, -1.1, 1.1, 220, -1.1, 1.1);
       
    hDeltaPhi[0] = new TH1D("deltaPhi", "deltaPhi", 180, -0.5, 179.5);
    hDeltaPhi[1] = new TH1D("deltaPhiEl", "deltaPhiEl", 180, -0.5, 179.5);
    hDeltaPhi[2] = new TH1D("deltaPhiInel", "deltaPhiinel", 180, -0.5, 179.5);   
}

