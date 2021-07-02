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
void ConnectInput(TTree* tree);
void PlotPionsPlot();
bool ProtonFiducial();
void RunGraniitti();

int protonsInside, protonsTotal;

void rafal()
{

    TString output = "/home/truhlar/Desktop/STAR/CEP/Analysis/Outputs/rafalNew.root";
    TString input = "/home/truhlar/Desktop/STAR/CEP/Analysis/Data/rafal.root";
    //TString graniitti_input = "/home/truhlar/Desktop/STAR/CEP/Analysis/Graniitti/2pi_100k_1.root";
    TString graniitti_input = "/home/truhlar/Desktop/STAR/Graniitti_new/GRANIITTI/output/RootFiles/200/200new.root";

    data = TFile::Open(input, "read");
    if (!data)
    {
        cout<<"Error: cannot open "<<input<<endl;
        return;
    }

    graniitti = TFile::Open(graniitti_input, "read");
    if (!graniitti)
    {
        cout<<"Error: cannot open "<<graniitti_input<<endl;
        return;
    }

    fout = new TFile(output,"RECREATE");
    Init(); // Preparing histograms 

    gStyle->SetLineWidth(2);      //axis line
    gStyle->SetFrameLineWidth(2); //frame line

    protonsInside = 0;
    protonsTotal = 0;

    RunGraniitti();
    //PlotRPPlot();
    fout->cd();
    PlotPionsPlot();

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

void ConnectInput(TTree* tree)
{

// PID and some quality event info
    tree->SetBranchAddress("mSquared", &mSquared);
    tree->SetBranchAddress("elastic", &elastic);
    tree->SetBranchAddress("nSigTrk1Pion", &nSigmaTPC[Pion][0]);
    tree->SetBranchAddress("nSigTrk2Pion", &nSigmaTPC[Pion][1]);
    tree->SetBranchAddress("nSigTrk3Pion", &nSigmaTPC[Pion][2]);
    tree->SetBranchAddress("nSigTrk4Pion", &nSigmaTPC[Pion][3]);
    tree->SetBranchAddress("mSquared", &mSquared);
    for (int iPart = 0; iPart < nParticles; ++iPart)
    {
        tree->SetBranchAddress("invMass" + particleLables[iPart], &invMass[iPart]);
        tree->SetBranchAddress("chiPair" + particleLables[iPart], &chiPair[iPart]);
    }

// RP track info  
    for (int i = 0; i < nSides; ++i)
    {
        tree->SetBranchAddress("t" + sideLabel[i], &t[i]);
        tree->SetBranchAddress("phiRp" + sideLabel[i], &phiRp[i]);
        tree->SetBranchAddress("xCorrelationsRp" + sideLabel[i], &xCorrelationsRp[i]);
        tree->SetBranchAddress("yCorrelationsRp" + sideLabel[i], &yCorrelationsRp[i]);
    }
// Vertex info
    tree->SetBranchAddress("vertexZ", &vertexesZ[0]);

// Central track info
    for (int i = 0; i < 4; ++i)
    {
        tree->SetBranchAddress(Form("momentum%i",i), &momentum[i]);
        tree->SetBranchAddress(Form("transMomentum%i",i), &transMomentum[i]);
        tree->SetBranchAddress(Form("charge%i",i), &charge[i]);
        tree->SetBranchAddress(Form("DcaXY%i",i), &DcaXY[i]);
        tree->SetBranchAddress(Form("DcaZ%i",i), &DcaZ[i]);
        tree->SetBranchAddress(Form("NhitsFit%i",i), &NhitsFit[i]);
        tree->SetBranchAddress(Form("NhitsDEdx%i",i), &NhitsDEdx[i]);
        tree->SetBranchAddress(Form("Eta%i",i), &Eta[i]);
        tree->SetBranchAddress(Form("Phi%i",i), &Phi[i]);

    }

// event info
    tree->SetBranchAddress("fourPiState", &fourPiState);

}

void PlotPionsPlot()
{
    fout->cd();
    TCanvas* newCanvas = new TCanvas("newCanvas","newCanvas",800,700);
    gPad->SetMargin(0.13,0.03,0.105,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
    gPad->SetTickx();
    gPad->SetTicky(); 
    gPad->SetLogy(0);

    TFile* rafalData = TFile::Open("/home/truhlar/Desktop/STAR/CEP/Analysis/Data/rafal.root", "read");
    fout->cd();
    TH1D *hPipi = (TH1D*)rafalData->Get("Figure 8/Hist1D_y1");
    TH1D *hPiInel = (TH1D*)rafalData->Get("Figure 12 (top, left)/Hist1D_y1");
    TH1D *hPiEl = (TH1D*)rafalData->Get("Figure 12 (top, right)/Hist1D_y1");


    TString stateLabel = "#pi^{+}#pi^{-}";

    TH1D *hist, *histCompare, *histGraniitti;
    TString yLabel = "Probability per event / 50 MeV";


    hist = (TH1D*)hPipi->Clone("hist"); 
    histGraniitti = (TH1D*)hInvMassGran[ElInel][Pion]->Clone("histGraniitti");
    cout<<"Here histGraniitti with "<< hInvMassGran[ElInel][Pion]->GetEntries()<<endl;
    Double_t scaleFactor; 
    scaleFactor =   1 /hist->Integral();
    hist->Scale(scaleFactor*2);
    //cout<<"Normalizing to "<<hist->Integral(12,64)<<endl;
    scaleFactor =   1 /histGraniitti->Integral();
    histGraniitti->Scale(scaleFactor);

    hist->SetTitle(" ; m(" + stateLabel + ") [GeV]; " + yLabel);
    hist->SetStats(false);
    hist->GetXaxis()->SetTitleFont(fontStyle);
    hist->GetXaxis()->SetTitleFont(fontStyle);
    hist->GetXaxis()->SetLabelFont(fontStyle);
    hist->GetYaxis()->SetLabelFont(fontStyle);
    hist->GetXaxis()->SetLabelSize(labelSize);
    hist->GetYaxis()->SetLabelSize(labelSize);
    hist->GetXaxis()->SetTitleSize(labelSize);
    hist->GetYaxis()->SetTitleSize(labelSize);
    hist->GetXaxis()->SetTitleOffset(0.9);
    hist->GetYaxis()->SetTitleOffset(1.3);
    hist->GetYaxis()->SetRangeUser(0, 0.082);   
    hist->SetMarkerColor(1);
    hist->SetMarkerSize(1);
    hist->SetMarkerStyle(20);
    hist->SetLineColor(1);
    hist->SetLineStyle(1);
    hist->SetLineWidth(1);
    hist->Draw("hist");

    gStyle->SetOptStat("");
    gStyle->SetPalette(1);


    histGraniitti->SetMarkerColor(4);
    histGraniitti->SetMarkerSize(1);
    histGraniitti->SetMarkerStyle(22);
    histGraniitti->SetLineColor(4);
    histGraniitti->SetLineStyle(1);
    histGraniitti->SetLineWidth(1);
    histGraniitti->Draw("same hist");


    TPaveText *textSTAR;
    textSTAR = new TPaveText(0.17,0.89,0.33,0.95,"brNDC");
    textSTAR -> SetTextSize(textSize+0.02);
    textSTAR -> SetFillColor(0);
    textSTAR -> SetTextFont(72);
    textSTAR -> AddText("STAR");
    textSTAR -> Draw("same");
    textSTAR = new TPaveText(0.17,0.84,0.33,0.89,"brNDC");
    textSTAR -> SetTextSize(textSize);
    textSTAR -> SetTextFont(52);
    textSTAR -> SetFillColor(0);
    textSTAR -> AddText("Internal");
    textSTAR -> Draw("same");

    TPaveText *textPub = new TPaveText(0.35,0.88,0.88,0.95,"brNDC");
    textPub -> SetTextSize(textSize);
    textPub -> SetTextAlign(22);
    textPub -> SetFillColor(0);
    textPub -> SetTextFont(42);
    textPub -> SetTextAlign(12);
    textPub -> AddText("p + p #rightarrow p + " + stateLabel + " + p       #sqrt{s} = 200 GeV");
    //textPub -> AddText("#sqrt{s} = 200 GeV");
    textPub -> Draw("same");


    TPaveText *text;
    text = new TPaveText(0.27,0.69,0.52,0.84,"brNDC");
    text -> SetTextSize(textSize);
    text -> SetFillColor(0);
    text -> SetTextFont(42);
    text -> SetTextAlign(32);
    text -> AddText("#pi^{+}, #pi^{-} kinematics:");   
    text -> AddText("p_{T} > 0.2 GeV");
    text -> AddText("|#eta| < 0.7");
    text -> Draw("same");

    text = new TPaveText(0.54,0.63,0.88,0.84,"brNDC");
    text -> SetTextSize(textSize);
    text -> SetFillColor(0);
    text -> SetTextFont(42);
    text -> SetTextAlign(12);
    text -> AddText("Forward proton kinematics:");
    text -> AddText("(p_{x} + 0.6)^{2} + p_{y}^{2} < 1.25 GeV^{2}");
    text -> AddText("0.4 GeV < |p_{y}| < 0.8 GeV");
    text -> AddText("p_{x} > -0.27 GeV");
    text -> Draw("same");


    TLegend *leg1 = new TLegend(0.54, 0.44, 0.87, 0.59);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(textSize);
    leg1->SetTextFont(42);
    leg1->SetMargin(0.1);
    leg1->AddEntry(hist, "Data (unlike-sign pairs)","ple");
    leg1->AddEntry(histGraniitti, "Graniitti","ple");
    leg1->Draw("same");

    text = new TPaveText(0.51,0.22,0.93,0.39,"brNDC");
    text -> SetTextSize(textSize);
    text -> SetFillColor(0);
    text -> SetTextFont(42);
    text -> SetTextAlign(22);
    text -> AddText("Statistical errors only");
    text -> AddText("Acceptance corrected");
    text -> AddText("Not background subtracted");
    text -> Draw("same");

    newCanvas->Update();
    newCanvas->Write("pions");
   // newCanvas->Close(); 
/////////////   El + Inel /////////////////////////////////////
    TH1D* histGraniittiEl;
    hist = (TH1D*)hPiInel->Clone("hist"); 
    histCompare = (TH1D*)hPiEl->Clone("histCompare");
    histGraniitti = (TH1D*)hInvMassGran[Inel][Pion]->Clone("histGraniittiInel");
    cout<<"Here histGraniittiInel with "<< hInvMassGran[Inel][Pion]->GetEntries()<<endl;
    histGraniittiEl = (TH1D*)hInvMassGran[El][Pion]->Clone("histGraniittiEl");
    cout<<"Here histGraniittiEl with "<< hInvMassGran[El][Pion]->GetEntries()<<endl;
    cout<<"_________________"<<endl;
    cout<<"Ratio-data (inel/el): "<<hist->Integral()<<" / "<<histCompare->Integral()<<" = "<<hist->Integral() / histCompare->Integral() <<endl;
    cout<<"Ratio-graniitti (inel/el): "<<histGraniitti->Integral()<<" / "<<histGraniittiEl->Integral()<<" = "<<histGraniitti->Integral() / histGraniittiEl->Integral() <<endl;
    cout<<"_________________"<<endl;

    scaleFactor =   1 /hist->Integral();
    hist->Scale(scaleFactor*2);
    scaleFactor =   1 /histCompare->Integral();
    histCompare->Scale(scaleFactor*2);
    scaleFactor =   1 /histGraniitti->Integral();
    histGraniitti->Scale(scaleFactor);
    scaleFactor =   1 /histGraniittiEl->Integral();
    histGraniittiEl->Scale(scaleFactor);

    hist->SetTitle(" ; m(" + stateLabel + ") [GeV]; " + yLabel);
    hist->SetStats(false);
    hist->GetXaxis()->SetTitleFont(fontStyle);
    hist->GetXaxis()->SetTitleFont(fontStyle);
    hist->GetXaxis()->SetLabelFont(fontStyle);
    hist->GetYaxis()->SetLabelFont(fontStyle);
    hist->GetXaxis()->SetLabelSize(labelSize);
    hist->GetYaxis()->SetLabelSize(labelSize);
    hist->GetXaxis()->SetTitleSize(labelSize);
    hist->GetYaxis()->SetTitleSize(labelSize);
    hist->GetXaxis()->SetTitleOffset(0.9);
    hist->GetYaxis()->SetTitleOffset(1.3);
    hist->GetYaxis()->SetRangeUser(0, 0.12);   
    hist->SetMarkerColor(4);
    hist->SetMarkerSize(1);
    hist->SetMarkerStyle(20);
    hist->SetLineColor(4);
    hist->SetLineStyle(1);
    hist->SetLineWidth(1);
    hist->Draw("hist");

    gStyle->SetOptStat("");
    gStyle->SetPalette(1);

    histCompare->SetMarkerColor(1);
    histCompare->SetMarkerSize(1);
    histCompare->SetMarkerStyle(21);
    histCompare->SetLineColor(1);
    histCompare->SetLineStyle(1);
    histCompare->SetLineWidth(1);
    histCompare->Draw("SAME hist");

    histGraniitti->SetMarkerColor(2);
    histGraniitti->SetMarkerSize(1);
    histGraniitti->SetMarkerStyle(21);
    histGraniitti->SetLineColor(2);
    histGraniitti->SetLineStyle(1);
    histGraniitti->SetLineWidth(1);
    histGraniitti->Draw("SAME hist");

    histGraniittiEl->SetMarkerColor(6);
    histGraniittiEl->SetMarkerSize(1);
    histGraniittiEl->SetMarkerStyle(21);
    histGraniittiEl->SetLineColor(6);
    histGraniittiEl->SetLineStyle(1);
    histGraniittiEl->SetLineWidth(1);
    histGraniittiEl->Draw("SAME hist");

    textSTAR = new TPaveText(0.17,0.89,0.33,0.95,"brNDC");
    textSTAR -> SetTextSize(textSize+0.02);
    textSTAR -> SetFillColor(0);
    textSTAR -> SetTextFont(72);
    textSTAR -> AddText("STAR");
    textSTAR -> Draw("same");
    textSTAR = new TPaveText(0.17,0.84,0.33,0.89,"brNDC");
    textSTAR -> SetTextSize(textSize);
    textSTAR -> SetTextFont(52);
    textSTAR -> SetFillColor(0);
    textSTAR -> AddText("Internal");
    textSTAR -> Draw("same");

    textPub = new TPaveText(0.35,0.88,0.88,0.95,"brNDC");
    textPub -> SetTextSize(textSize);
    textPub -> SetTextAlign(22);
    textPub -> SetFillColor(0);
    textPub -> SetTextFont(42);
    textPub -> SetTextAlign(12);
    textPub -> AddText("p + p #rightarrow p + " + stateLabel + " + p       #sqrt{s} = 200 GeV");
    //textPub -> AddText("#sqrt{s} = 200 GeV");
    textPub -> Draw("same");

    text = new TPaveText(0.27,0.68,0.52,0.83,"brNDC");
    text -> SetTextSize(textSize);
    text -> SetFillColor(0);
    text -> SetTextFont(42);
    text -> SetTextAlign(32);
    text -> AddText("#pi^{+}, #pi^{-} kinematics:");   
    text -> AddText("p_{T} > 0.2 GeV");
    text -> AddText("|#eta| < 0.7");
    text -> Draw("same");

    text = new TPaveText(0.54,0.62,0.88,0.83,"brNDC");
    text -> SetTextSize(textSize);
    text -> SetFillColor(0);
    text -> SetTextFont(42);
    text -> SetTextAlign(12);
    text -> AddText("Forward proton kinematics:");
    text -> AddText("(p_{x} + 0.6)^{2} + p_{y}^{2} < 1.25 GeV^{2}");
    text -> AddText("0.4 GeV < |p_{y}| < 0.8 GeV");
    text -> AddText("p_{x} > -0.27 GeV");
    text -> Draw("same");


    leg1 = new TLegend(0.45, 0.38, 0.78, 0.58);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(textSize);
    leg1->SetTextFont(42);
    leg1->SetMargin(0.1);
    leg1->AddEntry(hist, "Data, #Delta#varphi < 90^{#circ} (unlike-sign pairs)","ple");
    leg1->AddEntry(histCompare, "Data, #Delta#varphi > 90^{#circ} (unlike-sign pairs)","ple");
    leg1->AddEntry(histGraniitti, "Graniitti, #Delta#varphi < 90^{#circ}","ple");
    leg1->AddEntry(histGraniittiEl, "Graniitti, #Delta#varphi > 90^{#circ}","ple");
    leg1->Draw("same");

    text = new TPaveText(0.71,0.27,0.88,0.44,"brNDC");
    text -> SetTextSize(textSize);
    text -> SetFillColor(0);
    text -> SetTextFont(42);
    text -> SetTextAlign(22);
    text -> AddText("Statistical errors only");
    text -> AddText("Acceptance corrected");
    text -> AddText("Not background subtracted");
    text -> Draw("same");

    newCanvas->Update();
    newCanvas->Write("pionsEl+Inel");
///////////////////////////// Inelastic ///////////////////////////////////////////
    hist->SetTitle(" ; m(" + stateLabel + ") [GeV]; " + yLabel);
    hist->SetStats(false);
    hist->GetXaxis()->SetTitleFont(fontStyle);
    hist->GetXaxis()->SetTitleFont(fontStyle);
    hist->GetXaxis()->SetLabelFont(fontStyle);
    hist->GetYaxis()->SetLabelFont(fontStyle);
    hist->GetXaxis()->SetLabelSize(labelSize);
    hist->GetYaxis()->SetLabelSize(labelSize);
    hist->GetXaxis()->SetTitleSize(labelSize);
    hist->GetYaxis()->SetTitleSize(labelSize);
    hist->GetXaxis()->SetTitleOffset(0.9);
    hist->GetYaxis()->SetTitleOffset(1.3);
    hist->GetYaxis()->SetRangeUser(0, 0.12);   
    hist->SetMarkerColor(4);
    hist->SetMarkerSize(1);
    hist->SetMarkerStyle(20);
    hist->SetLineColor(4);
    hist->SetLineStyle(1);
    hist->SetLineWidth(1);
    hist->Draw("hist");

    gStyle->SetOptStat("");
    gStyle->SetPalette(1);


    histGraniitti->SetMarkerColor(2);
    histGraniitti->SetMarkerSize(1);
    histGraniitti->SetMarkerStyle(21);
    histGraniitti->SetLineColor(2);
    histGraniitti->SetLineStyle(1);
    histGraniitti->SetLineWidth(1);
    histGraniitti->Draw("SAME hist");


    textSTAR = new TPaveText(0.17,0.89,0.33,0.95,"brNDC");
    textSTAR -> SetTextSize(textSize+0.02);
    textSTAR -> SetFillColor(0);
    textSTAR -> SetTextFont(72);
    textSTAR -> AddText("STAR");
    textSTAR -> Draw("same");
    textSTAR = new TPaveText(0.17,0.84,0.33,0.89,"brNDC");
    textSTAR -> SetTextSize(textSize);
    textSTAR -> SetTextFont(52);
    textSTAR -> SetFillColor(0);
    textSTAR -> AddText("Internal");
    textSTAR -> Draw("same");

    textPub = new TPaveText(0.35,0.88,0.88,0.95,"brNDC");
    textPub -> SetTextSize(textSize);
    textPub -> SetTextAlign(22);
    textPub -> SetFillColor(0);
    textPub -> SetTextFont(42);
    textPub -> SetTextAlign(12);
    textPub -> AddText("p + p #rightarrow p + " + stateLabel + " + p       #sqrt{s} = 200 GeV");
    //textPub -> AddText("#sqrt{s} = 200 GeV");
    textPub -> Draw("same");

    text = new TPaveText(0.27,0.68,0.52,0.83,"brNDC");
    text -> SetTextSize(textSize);
    text -> SetFillColor(0);
    text -> SetTextFont(42);
    text -> SetTextAlign(32);
    text -> AddText("#pi^{+}, #pi^{-} kinematics:");   
    text -> AddText("p_{T} > 0.2 GeV");
    text -> AddText("|#eta| < 0.7");
    text -> Draw("same");

    text = new TPaveText(0.54,0.62,0.88,0.83,"brNDC");
    text -> SetTextSize(textSize);
    text -> SetFillColor(0);
    text -> SetTextFont(42);
    text -> SetTextAlign(12);
    text -> AddText("Forward proton kinematics:");
    text -> AddText("(p_{x} + 0.6)^{2} + p_{y}^{2} < 1.25 GeV^{2}");
    text -> AddText("0.4 GeV < |p_{y}| < 0.8 GeV");
    text -> AddText("p_{x} > -0.27 GeV");
    text -> Draw("same");


    leg1 = new TLegend(0.45, 0.38, 0.78, 0.58);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(textSize);
    leg1->SetTextFont(42);
    leg1->SetMargin(0.1);
    leg1->AddEntry(hist, "Data, #Delta#varphi < 90^{#circ} (unlike-sign pairs)","ple");
    leg1->AddEntry(histGraniitti, "Graniitti, #Delta#varphi < 90^{#circ}","ple");
    leg1->Draw("same");

    text = new TPaveText(0.71,0.27,0.88,0.44,"brNDC");
    text -> SetTextSize(textSize);
    text -> SetFillColor(0);
    text -> SetTextFont(42);
    text -> SetTextAlign(22);
    text -> AddText("Statistical errors only");
    text -> AddText("Acceptance corrected");
    text -> AddText("Not background subtracted");
    text -> Draw("same");

    newCanvas->Update();
    newCanvas->Write("pionsInel");
//////////////////////////////////////////// Elastic ///////////////////////////////

    histCompare->SetTitle(" ; m(" + stateLabel + ") [GeV]; " + yLabel);
    histCompare->SetStats(false);
    histCompare->GetXaxis()->SetTitleFont(fontStyle);
    histCompare->GetXaxis()->SetTitleFont(fontStyle);
    histCompare->GetXaxis()->SetLabelFont(fontStyle);
    histCompare->GetYaxis()->SetLabelFont(fontStyle);
    histCompare->GetXaxis()->SetLabelSize(labelSize);
    histCompare->GetYaxis()->SetLabelSize(labelSize);
    histCompare->GetXaxis()->SetTitleSize(labelSize);
    histCompare->GetYaxis()->SetTitleSize(labelSize);
    histCompare->GetXaxis()->SetTitleOffset(0.9);
    histCompare->GetYaxis()->SetTitleOffset(1.3);
    histCompare->GetYaxis()->SetRangeUser(0, 0.12);   
    histCompare->SetMarkerColor(1);
    histCompare->SetMarkerSize(1);
    histCompare->SetMarkerStyle(21);
    histCompare->SetLineColor(1);
    histCompare->SetLineStyle(1);
    histCompare->SetLineWidth(1);
    histCompare->Draw("hist");
    histCompare->Draw("SAME hist");
    gStyle->SetOptStat("");
    gStyle->SetPalette(1);



    histGraniittiEl->SetMarkerColor(6);
    histGraniittiEl->SetMarkerSize(1);
    histGraniittiEl->SetMarkerStyle(21);
    histGraniittiEl->SetLineColor(6);
    histGraniittiEl->SetLineStyle(1);
    histGraniittiEl->SetLineWidth(1);
    histGraniittiEl->Draw("SAME hist");

    textSTAR = new TPaveText(0.17,0.89,0.33,0.95,"brNDC");
    textSTAR -> SetTextSize(textSize+0.02);
    textSTAR -> SetFillColor(0);
    textSTAR -> SetTextFont(72);
    textSTAR -> AddText("STAR");
    textSTAR -> Draw("same");
    textSTAR = new TPaveText(0.17,0.84,0.33,0.89,"brNDC");
    textSTAR -> SetTextSize(textSize);
    textSTAR -> SetTextFont(52);
    textSTAR -> SetFillColor(0);
    textSTAR -> AddText("Internal");
    textSTAR -> Draw("same");

    textPub = new TPaveText(0.35,0.88,0.88,0.95,"brNDC");
    textPub -> SetTextSize(textSize);
    textPub -> SetTextAlign(22);
    textPub -> SetFillColor(0);
    textPub -> SetTextFont(42);
    textPub -> SetTextAlign(12);
    textPub -> AddText("p + p #rightarrow p + " + stateLabel + " + p       #sqrt{s} = 200 GeV");
    //textPub -> AddText("#sqrt{s} = 200 GeV");
    textPub -> Draw("same");

    text = new TPaveText(0.27,0.68,0.52,0.83,"brNDC");
    text -> SetTextSize(textSize);
    text -> SetFillColor(0);
    text -> SetTextFont(42);
    text -> SetTextAlign(32);
    text -> AddText("#pi^{+}, #pi^{-} kinematics:");   
    text -> AddText("p_{T} > 0.2 GeV");
    text -> AddText("|#eta| < 0.7");
    text -> Draw("same");

    text = new TPaveText(0.54,0.62,0.88,0.83,"brNDC");
    text -> SetTextSize(textSize);
    text -> SetFillColor(0);
    text -> SetTextFont(42);
    text -> SetTextAlign(12);
    text -> AddText("Forward proton kinematics:");
    text -> AddText("(p_{x} + 0.6)^{2} + p_{y}^{2} < 1.25 GeV^{2}");
    text -> AddText("0.4 GeV < |p_{y}| < 0.8 GeV");
    text -> AddText("p_{x} > -0.27 GeV");
    text -> Draw("same");


    leg1 = new TLegend(0.45, 0.38, 0.78, 0.58);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(textSize);
    leg1->SetTextFont(42);
    leg1->SetMargin(0.1);
    leg1->AddEntry(histCompare, "Data, #Delta#varphi > 90^{#circ} (unlike-sign pairs)","ple");
    leg1->AddEntry(histGraniittiEl, "Graniitti, #Delta#varphi > 90^{#circ}","ple");
    leg1->Draw("same");

    text = new TPaveText(0.71,0.27,0.88,0.44,"brNDC");
    text -> SetTextSize(textSize);
    text -> SetFillColor(0);
    text -> SetTextFont(42);
    text -> SetTextAlign(22);
    text -> AddText("Statistical errors only");
    text -> AddText("Acceptance corrected");
    text -> AddText("Not background subtracted");
    text -> Draw("same");

    newCanvas->Update();
    newCanvas->Write("pionsEl");
}


bool ProtonFiducial()
{
    double x,y;
    bool protonInRange[nSides] = {false, false};
    for (int i = 0; i < nSides; ++i)
    {
        x = xCorrelationsRp[i];
        y = yCorrelationsRp[i];
        if( abs(y) < 0.8 && abs(y) > 0.4 && x > -0.27 && (x + 0.6)*(x + 0.6) + y*y < 1.25 )
            protonInRange[i] = true;
    }

    if(protonInRange[0] && protonInRange[1])
        return true;
    else
        return false;
}


void RunGraniitti()
{

    TTree* tree[nParticles]; // 3 = 4PI state
    tree[Pion] = dynamic_cast<TTree*>( graniitti->Get("pionTree") );
//    tree[Kaon] = dynamic_cast<TTree*>( graniitti->Get("kaonTree") );
//    tree[Proton] = dynamic_cast<TTree*>( graniitti->Get("protonTree") );
//    tree[FourPions] = dynamic_cast<TTree*>( graniitti->Get("4pionTree") );

    if (!tree[Pion] ){//|| !tree[Kaon] || !tree[Proton]){
        cout<<"Error: cannot open one of the TTree in Graniitti"<<endl;
        return;
    }

    TString granForwardCuts =   TString("px_proton1 > - 0.2 && px_proton2 > - 0.2 &&")
                                + TString("abs(py_proton1) < 0.4 && abs(py_proton1) > 0.2 && abs(py_proton2) < 0.4 && abs(py_proton2) > 0. &&")
                                + TString("(px_proton1 + 0.3)*(px_proton1 + 0.3) +py_proton1*py_proton1 < 0.25 &&")
                                + TString("(px_proton2 + 0.3)*(px_proton2 + 0.3) +py_proton2*py_proton2 < 0.25");
    TString graniittiCuts = "eta_part1 > -0.7 && eta_part1 < 0.7 && eta_part2 > -0.7 && eta_part2 < 0.7 && t_proton1 < -0.02 && t_proton2 < -0.02 && t_proton1 > -0.24  && t_proton2 > -0.24 && " + granForwardCuts;
    TString partCuts[nParticles];
    partCuts[Pion] = "pT_part1 > 0.2 && pT_part2 > 0.2";
    partCuts[Kaon] = "pT_part1 > 0.3 && pT_part2 > 0.3 && TMath::Min(pT_part1,pT_part2) < 0.7";
    partCuts[Proton] = "pT_part1 > 0.4 && pT_part2 > 0.4 && TMath::Min(pT_part1,pT_part2) < 1.1";
    //partCuts[FourPions] = "pT_part1 > 0.2 && pT_part2 > 0.2 && pT_part3 > 0.2 && pT_part4 > 0.2";

    Double_t binning[4][3] =  {{64, 0.3, 3.5},{44, 0.8, 3}, {24, 1.6, 4}, {50,0.5,4.5}};

    TString variable, usedCuts; 
    TString combCuts[] = { TString(""), TString("&& TMath::Abs(deltaphi_pp) > 1.570796327"), TString("&& TMath::Abs(deltaphi_pp) < 1.570796327")}; // 57.2957795 = 1 rad => 1.6 = 90 cca
    TString stateLabel[] = {"#pi^{+}#pi^{-}","K^{+}K^{-}","p#bar{p}"};

    Int_t nBins;
    Float_t min, max;

    for (int iComb = 0; iComb < nCombination; ++iComb)
    {
        usedCuts= graniittiCuts + combCuts[iComb];
        for (int iPart = 0; iPart < 1; ++iPart) // nParticles; ++iPart)
        {
            usedCuts+= " && " + partCuts[iPart];    
            variable = "invMass_state";
            nBins = binning[iPart][0];
            min = binning[iPart][1];
            max = binning[iPart][2];

            tree[iPart]->Draw(variable +">>" + variable +"Sig1(" + nBins + "," + min + "," + max + ")",usedCuts);
            hInvMassGran[iComb][iPart] = (TH1D*)gPad->GetPrimitive(variable +"Sig1")->Clone(Form("granInvMass_%i_%i",iComb, iPart)); 
            cout<<"Combination: "<<iComb<<" Particle: "<<iPart<<" Entries: "<< hInvMassGran[iComb][iPart]->GetEntries()<< endl;  
        }
    }

}