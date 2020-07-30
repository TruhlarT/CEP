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

TTree* tree;
TTree* treeBack;

TH1F* hInvMass[4];
TH1D* hDeltaPhi[3];
TH2D* hPxPy[3];


TH1D* hZvertex[2];   
TH1D* hT[2]; 
TH1D* hNhits[2];   
TH1D* hNhitsFit[2];   
TH1D* hDCAxy[2]; 
TH1D* hDCAz[2]; 
TH1D* hEta[2]; 

 // 0 = pi- 1 = K- 2 = pbar
TH3F* hTPCeff[6]; // 3 = pi+ 4 = K+ 5 = p
TH1D* hInvMassCorr[2]; // 0 - signal, 1 - Background
                                    //  - El + Inel, 1 - El, 2 - Inel 
TH1D* hInvMassUncorr[2]; // 0 - signal, 1 - Background

TH1D* hRapidityCorr[nParticles][2][3];  
TH1D* hRapidityUncorr[nParticles][2][3]; 

TH1D* hDeltaPhiCorr[nParticles][2][3];  
TH1D* hDeltaPhiUncorr[nParticles][2][3]; 

Int_t nTracks, totalCharge, nTofTrks; 
UInt_t runNumber;
Double_t VPDSumEast, VPDSumWest, VPDTimeDiff;


Double_t chiPair[nParticles]; 
Double_t invMass[nParticles];
Double_t missingPt, deltaTOF, mSquared, pairRapidity;
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
const double labelSize = 0.045;

void Init();
void ConnectInput(TTree* tree);
void Make(int signal);


void PlotPlots();
void PlotCutsFlow();

void FourPiPlot()
{

    TString output = "/home/truhlar/Desktop/STAR/CEP/Analysis/Outputs/4piPlots.root";
    TString input = "/home/truhlar/Desktop/STAR/CEP/Analysis/Data/4pi.root";

    TString TPCeffInput[6];
    TPCeffInput[0] = "/home/truhlar/Desktop/STAR/CEP/Analysis/Data/NewEff/effPionsM.root";
    TPCeffInput[1] = "/home/truhlar/Desktop/STAR/CEP/Analysis/Data/NewEff/effKaonsM.root";
    TPCeffInput[2] = "/home/truhlar/Desktop/STAR/CEP/Analysis/Data/NewEff/effProtonsM.root";
    TPCeffInput[3] = "/home/truhlar/Desktop/STAR/CEP/Analysis/Data/NewEff/effPionsP.root";
    TPCeffInput[4] = "/home/truhlar/Desktop/STAR/CEP/Analysis/Data/NewEff/effKaonsP.root";
    TPCeffInput[5] = "/home/truhlar/Desktop/STAR/CEP/Analysis/Data/NewEff/effProtonsP.root";

    for (int i = 0; i < 6; ++i)
    {
        TPCeff[i] = TFile::Open(TPCeffInput[i], "read");
        if (!TPCeff[i])
        {
            cout<<"Error: cannot open "<<TPCeffInput[i]<<endl;
            return 4 +i;
        }
        hTPCeff[i] = (TH3F*)TPCeff[i] -> Get("effRafal"); 
    }


    data = TFile::Open(input, "read");
    if (!data)
    {
        cout<<"Error: cannot open "<<input<<endl;
        return;
    }


    fout = new TFile(output,"RECREATE");
    Init(); // Preparing histograms 

    tree = dynamic_cast<TTree*>( data->Get("recTree") );
    treeBack = dynamic_cast<TTree*>( data->Get("Background") );

    gStyle->SetLineWidth(2);      //axis line
    gStyle->SetFrameLineWidth(2); //frame line


    ConnectInput(tree);
    for(Long64_t iev=0; iev<tree->GetEntries(); ++iev)
    { //get the event
        tree->GetEntry(iev); 
        Make(0);
    } 

    ConnectInput(treeBack);
    for(Long64_t iev=0; iev<treeBack->GetEntries(); ++iev)
    { //get the event
        treeBack->GetEntry(iev); 
        Make(1);
    }


    PlotPlots();
    PlotCutsFlow();

    fout->Write();
    fout->Close();
    
}




void Make(int signal)
{
    double effTotal, effTPC, effTOF;
    unsigned int PID;
   // cout<< vertexesZ[0] <<" "<< NhitsFit[0]<<" "<<NhitsFit[1] <<" "<< NhitsDEdx[0]<<" "<<NhitsDEdx[1] <<" "<<DcaZ[0] <<" "<<DcaZ[1] <<" "<<DcaXY[0] <<" "<<DcaXY[1] <<" "<<Eta[0] <<" "<<Eta[1] <<" "<< !fourPiState<<endl;
            
    if(fourPiState)
    {
        hZvertex[signal]->Fill(vertexesZ[0]);  
        hT[signal]->Fill(t[0]); 
        hT[signal]->Fill(t[1]); 
        for (int iTrack = 0; iTrack < 4; ++iTrack)
        {
            hNhits[signal]->Fill(NhitsDEdx[iTrack]);   
            hNhitsFit[signal]->Fill(NhitsFit[iTrack]);   
            hDCAxy[signal]->Fill(DcaXY[iTrack]);  
            hDCAz[signal]->Fill(DcaZ[iTrack]); 
            hEta[signal]->Fill(Eta[iTrack]);

        }
    }


    if(vertexesZ[0] < 80 && vertexesZ[0] > -80 && NhitsFit[0] >=25 && NhitsFit[1] >= 25 && NhitsDEdx[0] >= 15 && 
    NhitsDEdx[1] >= 15 && DcaZ[0] < 1 && DcaZ[0] > -1 && DcaZ[1] < 1 && DcaZ[1] > -1 && DcaXY[0] < 1.5 && 
    DcaXY[1] < 1.5 && Eta[0] > -0.7 && Eta[0] < 0.7 && Eta[1] > -0.7 && Eta[1] < 0.7 &&  t[0] < -0.12 && t[1] < -0.12 && 
    t[0] > -1.0  && t[1] > -1.0 && fourPiState && NhitsFit[2] >=25 && NhitsFit[3] >= 25 && NhitsDEdx[2] >= 15 && 
    NhitsDEdx[3] >= 15 && DcaZ[2] < 1 && DcaZ[2] > -1 && DcaZ[3] < 1 && DcaZ[3] > -1 && DcaXY[2] < 1.5 && DcaXY[3] < 1.5 && 
    Eta[2] > -0.7 && Eta[2] < 0.7 && Eta[3] > -0.7 && Eta[3] < 0.7)
    {
 
        effTotal = 1;
        effTOF = 1;
        if(nSigmaTPC[Pion][0] < 3 && nSigmaTPC[Pion][1] < 3 && nSigmaTPC[Pion][2] < 3 && nSigmaTPC[Pion][3] < 3) 
        {
            for (int iTrack = 0; iTrack < 4; ++iTrack)
            {
                PID = 0;
                if(charge[iTrack] > 0)
                    PID = 3;
                
                Double_t phiToMC = Phi[iTrack];
                if( phiToMC < 0)
                    phiToMC = 2*3.14159265359 + phiToMC;
                effTPC = hTPCeff[0 + PID]->GetBinContent( hTPCeff[0 + PID]->GetXaxis()->FindBin(phiToMC), hTPCeff[0 + PID]->GetYaxis()->FindBin(transMomentum[iTrack]), hTPCeff[0 + PID]->GetZaxis()->FindBin(Eta[iTrack])); 
                                
                //cout<<hTPCeff[0 + PID]->GetXaxis()->FindBin(vertexesZ[iTrack])<<" | "<< hTPCeff[0 + PID]->GetYaxis()->FindBin(transMomentum[iTrack])<<" | "<< hTPCeff[0 + PID]->GetZaxis()->FindBin(Eta[iTrack])<<endl;

                effTotal = effTotal*effTPC*effTOF;
            }
            if(effTotal != 0 && transMomentum[0] > 0.2 && transMomentum[1] > 0.2 && transMomentum[2] > 0.2 && transMomentum[3] > 0.2)
            {
                hInvMassCorr[signal]->Fill(invMass[Pion], 1/effTotal);
                hInvMassCorr[signal]->Fill(invMass[Pion], 1/effTotal);


            }
            if(transMomentum[0] > 0.2 && transMomentum[1] > 0.2 && transMomentum[2] > 0.2 && transMomentum[3] > 0.2)
            {
                hInvMassUncorr[signal]->Fill(invMass[Pion]);
                hInvMassUncorr[signal]->Fill(invMass[Pion]);
            }

        }

    }
}


void Init()
{

   
    hInvMassCorr[0]  = new TH1D("corrInvMassSig", "Corrected inv. mass ", 50, 0.5, 4.5);
    hInvMassCorr[1]  = new TH1D("corrInvMassBcg", "Corrected inv. mass ", 50, 0.5, 4.5);

    hInvMassUncorr[0]  = new TH1D("uncorrInvMassSig", "Uncorrected inv. mass ", 50, 0.5, 4.5);
    hInvMassUncorr[1]  = new TH1D("uncorrInvMassBcq", "Uncorrected inv. mass ", 50, 0.5, 4.5);

    hZvertex[0]  = new TH1D("ZvertexSig", "Corrected inv. mass ", 100, -200, 250);
    hZvertex[1]  = new TH1D("ZvertexBcg", "Corrected inv. mass ", 100, -200, 250);

    hNhits[0]  = new TH1D("nHitsSig", "Corrected inv. mass ", 61, 0.0, 60);
    hNhits[1]  = new TH1D("nHitsBcg", "Corrected inv. mass ", 61, 0.0, 60);

    hNhitsFit[0]  = new TH1D("nHitsFitSig", "Corrected inv. mass ", 61, 0.0, 60);
    hNhitsFit[1]  = new TH1D("nHitsFitBcg", "Corrected inv. mass ", 61, 0.0, 60);

    hDCAxy[0]  = new TH1D("DCAxySig", "Corrected inv. mass ", 100, 0.0, 3.5);
    hDCAxy[1]  = new TH1D("DCAxyBcg", "Corrected inv. mass ", 100, 0.0, 3.5);

    hDCAz[0]  = new TH1D("DCAzSig", "Corrected inv. mass ", 100, -1.5, 1.5);
    hDCAz[1]  = new TH1D("DCAzBcg", "Corrected inv. mass ", 100, -1.5, 1.5);

    hEta[0]  = new TH1D("etaSig", "Corrected inv. mass ", 100, -2.0, 3.5);
    hEta[1]  = new TH1D("etaBcg", "Corrected inv. mass ", 100, -2.0, 3.5);

    hT[0]  = new TH1D("tSig", "Corrected inv. mass ", 100, -2.0, 0.0);
    hT[1]  = new TH1D("tBcg", "Corrected inv. mass ", 100, -2.0, 0.0);

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
    tree->SetBranchAddress("pairRapidity", &pairRapidity);
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


void PlotPlots()
{

    TCanvas *cCanvas = new TCanvas("cCanvas","cCanvas",800,700);
    gPad->SetMargin(0.12,0.03,0.1,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
    gPad->SetTickx();
    gPad->SetTicky();  
    gStyle->SetOptStat("");
    gStyle->SetPalette(1);
    gStyle->SetLineWidth(2);      //axis line
    gStyle->SetFrameLineWidth(2); //frame line

    TH1D* histSignal;
    TH1D* histBackground;

    histSignal = (TH1D*)hZvertex[0]->Clone("histSignal");
    histBackground = (TH1D*)hZvertex[1]->Clone("histBackground");

    histSignal->SetTitle(" ; Z_{vrtx} [cm]; Number of events");
    histSignal->SetStats(false);
    histSignal->GetXaxis()->SetTitleFont(42);
    histSignal->GetYaxis()->SetTitleFont(42);
    histSignal->GetXaxis()->SetLabelFont(42);
    histSignal->GetYaxis()->SetLabelFont(42);
    histSignal->GetXaxis()->SetLabelSize(labelSize);
    histSignal->GetYaxis()->SetLabelSize(labelSize);
    histSignal->GetXaxis()->SetTitleSize(labelSize);
    histSignal->GetYaxis()->SetTitleSize(labelSize);
    histSignal->GetXaxis()->SetTitleOffset(1.0);
    histSignal->GetYaxis()->SetTitleOffset(1.3);  
    histSignal->GetYaxis()->SetRangeUser(0.0, TMath::Max(histSignal->GetMaximum(),histBackground->GetMaximum())*1.2);
    histSignal->SetMarkerColor(1);
    histSignal->SetMarkerSize(1);
    histSignal->SetMarkerStyle(20);
    histSignal->SetLineColor(1);
    histSignal->SetLineStyle(1);
    histSignal->SetLineWidth(1);
    histSignal->Draw("E");
    
    gStyle->SetOptStat("");
    gStyle->SetPalette(1);

    histBackground->SetMarkerColor(2);
    histBackground->SetMarkerSize(1);
    histBackground->SetMarkerStyle(22);
    histBackground->SetLineColor(2);
    histBackground->SetLineStyle(1);
    histBackground->SetLineWidth(1);
    histBackground->Draw("ESAME");


    TPaveText *textPub = new TPaveText(0.65,0.75,0.92,0.88,"brNDC");
    textPub -> SetTextSize(textSize);
    textPub -> SetTextAlign(22);
    textPub -> SetFillColor(0);
    textPub -> SetTextFont(42);
    textPub -> AddText("p + p #rightarrow p + #pi^{+}#pi^{+}#pi^{-}#pi^{-} + p");
    textPub -> AddText("#sqrt{s} = 510 GeV");
    textPub -> Draw("same");

    TPaveText *textSTAR;
    textSTAR = new TPaveText(0.65,0.89,0.9,0.95,"brNDC");
    textSTAR -> SetTextSize(textSize);
    textSTAR -> SetTextAlign(22);
    textSTAR -> SetFillColor(0);
    textSTAR -> SetTextFont(62);
    textSTAR->AddText("THIS THESIS");
    textSTAR -> Draw("same");

    TLegend* leg1 = new TLegend(0.55, 0.65, 0.78, 0.74);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(textSize);
    leg1->SetTextFont(42);
    leg1->AddEntry(histSignal,"Data (unlike-sign pairs)","pe");
    leg1->AddEntry(histBackground,"Data (like-sign pairs)","pe");
    leg1->Draw("same");


    TLine *left02 = new TLine(-80,0.0,-80,histSignal->GetMaximum()*0.6);
    left02->SetLineStyle(10);
    left02->SetLineColor(1);
    left02->SetLineWidth(4);
    left02->Draw("same");

    left02 = new TLine(80,0.0,80,histSignal->GetMaximum()*0.6);
    left02->SetLineStyle(10);
    left02->SetLineColor(1);
    left02->SetLineWidth(4);
    left02->Draw("same");

    cCanvas->Update();
    cCanvas->Write("zVertex");


    histSignal = (TH1D*)hNhits[0]->Clone("histSignal");
    histBackground = (TH1D*)hNhits[1]->Clone("histBackground");

    histSignal->SetTitle(" ; N^{dE/dx}_{hits} ; Number of tracks");
    histSignal->SetStats(false);
    histSignal->GetXaxis()->SetTitleFont(42);
    histSignal->GetYaxis()->SetTitleFont(42);
    histSignal->GetXaxis()->SetLabelFont(42);
    histSignal->GetYaxis()->SetLabelFont(42);
    histSignal->GetXaxis()->SetLabelSize(labelSize);
    histSignal->GetYaxis()->SetLabelSize(labelSize);
    histSignal->GetXaxis()->SetTitleSize(labelSize);
    histSignal->GetYaxis()->SetTitleSize(labelSize);
    histSignal->GetXaxis()->SetTitleOffset(1.0);
    histSignal->GetYaxis()->SetTitleOffset(1.3);  
    histSignal->GetYaxis()->SetRangeUser(0.0, TMath::Max(histSignal->GetMaximum(),histBackground->GetMaximum())*1.2);
    histSignal->SetMarkerColor(1);
    histSignal->SetMarkerSize(1);
    histSignal->SetMarkerStyle(20);
    histSignal->SetLineColor(1);
    histSignal->SetLineStyle(1);
    histSignal->SetLineWidth(1);
    histSignal->Draw("E");
    
    gStyle->SetOptStat("");
    gStyle->SetPalette(1);

    histBackground->SetMarkerColor(2);
    histBackground->SetMarkerSize(1);
    histBackground->SetMarkerStyle(22);
    histBackground->SetLineColor(2);
    histBackground->SetLineStyle(1);
    histBackground->SetLineWidth(1);
    histBackground->Draw("ESAME");


    textPub = new TPaveText(0.65,0.75,0.92,0.88,"brNDC");
    textPub -> SetTextSize(textSize);
    textPub -> SetTextAlign(22);
    textPub -> SetFillColor(0);
    textPub -> SetTextFont(42);
    textPub -> AddText("p + p #rightarrow p + #pi^{+}#pi^{+}#pi^{-}#pi^{-} + p");
    textPub -> AddText("#sqrt{s} = 510 GeV");
    textPub -> Draw("same");

    textSTAR = new TPaveText(0.65,0.89,0.9,0.95,"brNDC");
    textSTAR -> SetTextSize(textSize);
    textSTAR -> SetTextAlign(22);
    textSTAR -> SetFillColor(0);
    textSTAR -> SetTextFont(62);
    textSTAR->AddText("THIS THESIS");
    textSTAR -> Draw("same");

    leg1 = new TLegend(0.55, 0.65, 0.78, 0.74);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(textSize);
    leg1->SetTextFont(42);
    leg1->AddEntry(histSignal,"Data (unlike-sign pairs)","pe");
    leg1->AddEntry(histBackground,"Data (like-sign pairs)","pe");
    leg1->Draw("same");

    left02 = new TLine(15,0.0,15,histSignal->GetMaximum()*0.6);
    left02->SetLineStyle(10);
    left02->SetLineColor(1);
    left02->SetLineWidth(4);
    left02->Draw("same");

    cCanvas->Update();
    cCanvas->Write("nHits");


    histSignal = (TH1D*)hNhitsFit[0]->Clone("histSignal");
    histBackground = (TH1D*)hNhitsFit[1]->Clone("histBackground");

    histSignal->SetTitle(" ; N^{fit}_{hits} ; Number of tracks");
    histSignal->SetStats(false);
    histSignal->GetXaxis()->SetTitleFont(42);
    histSignal->GetYaxis()->SetTitleFont(42);
    histSignal->GetXaxis()->SetLabelFont(42);
    histSignal->GetYaxis()->SetLabelFont(42);
    histSignal->GetXaxis()->SetLabelSize(labelSize);
    histSignal->GetYaxis()->SetLabelSize(labelSize);
    histSignal->GetXaxis()->SetTitleSize(labelSize);
    histSignal->GetYaxis()->SetTitleSize(labelSize);
    histSignal->GetXaxis()->SetTitleOffset(1.0);
    histSignal->GetYaxis()->SetTitleOffset(1.3);  
    histSignal->GetYaxis()->SetRangeUser(0.0, TMath::Max(histSignal->GetMaximum(),histBackground->GetMaximum())*1.2);
    histSignal->SetMarkerColor(1);
    histSignal->SetMarkerSize(1);
    histSignal->SetMarkerStyle(20);
    histSignal->SetLineColor(1);
    histSignal->SetLineStyle(1);
    histSignal->SetLineWidth(1);
    histSignal->Draw("E");
    
    gStyle->SetOptStat("");
    gStyle->SetPalette(1);

    histBackground->SetMarkerColor(2);
    histBackground->SetMarkerSize(1);
    histBackground->SetMarkerStyle(22);
    histBackground->SetLineColor(2);
    histBackground->SetLineStyle(1);
    histBackground->SetLineWidth(1);
    histBackground->Draw("ESAME");


    textPub = new TPaveText(0.65,0.75,0.92,0.88,"brNDC");
    textPub -> SetTextSize(textSize);
    textPub -> SetTextAlign(22);
    textPub -> SetFillColor(0);
    textPub -> SetTextFont(42);
    textPub -> AddText("p + p #rightarrow p + #pi^{+}#pi^{+}#pi^{-}#pi^{-} + p");
    textPub -> AddText("#sqrt{s} = 510 GeV");
    textPub -> Draw("same");

    textSTAR = new TPaveText(0.65,0.89,0.9,0.95,"brNDC");
    textSTAR -> SetTextSize(textSize);
    textSTAR -> SetTextAlign(22);
    textSTAR -> SetFillColor(0);
    textSTAR -> SetTextFont(62);
    textSTAR->AddText("THIS THESIS");
    textSTAR -> Draw("same");

    leg1 = new TLegend(0.55, 0.65, 0.78, 0.74);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(textSize);
    leg1->SetTextFont(42);
    leg1->AddEntry(histSignal,"Data (unlike-sign pairs)","pe");
    leg1->AddEntry(histBackground,"Data (like-sign pairs)","pe");
    leg1->Draw("same");

    left02 = new TLine(25,0.0,25,histSignal->GetMaximum()*0.6);
    left02->SetLineStyle(10);
    left02->SetLineColor(1);
    left02->SetLineWidth(4);
    left02->Draw("same");

    cCanvas->Update();
    cCanvas->Write("nHitsFit");


    histSignal = (TH1D*)hDCAxy[0]->Clone("histSignal");
    histBackground = (TH1D*)hDCAxy[1]->Clone("histBackground");

    histSignal->SetTitle(" ; DCA_{xy} [cm]; Number of tracks");
    histSignal->SetStats(false);
    histSignal->GetXaxis()->SetTitleFont(42);
    histSignal->GetYaxis()->SetTitleFont(42);
    histSignal->GetXaxis()->SetLabelFont(42);
    histSignal->GetYaxis()->SetLabelFont(42);
    histSignal->GetXaxis()->SetLabelSize(labelSize);
    histSignal->GetYaxis()->SetLabelSize(labelSize);
    histSignal->GetXaxis()->SetTitleSize(labelSize);
    histSignal->GetYaxis()->SetTitleSize(labelSize);
    histSignal->GetXaxis()->SetTitleOffset(1.0);
    histSignal->GetYaxis()->SetTitleOffset(1.3);  
    histSignal->GetYaxis()->SetRangeUser(0.0, TMath::Max(histSignal->GetMaximum(),histBackground->GetMaximum())*1.2);
    histSignal->SetMarkerColor(1);
    histSignal->SetMarkerSize(1);
    histSignal->SetMarkerStyle(20);
    histSignal->SetLineColor(1);
    histSignal->SetLineStyle(1);
    histSignal->SetLineWidth(1);
    histSignal->Draw("E");
    
    gStyle->SetOptStat("");
    gStyle->SetPalette(1);

    histBackground->SetMarkerColor(2);
    histBackground->SetMarkerSize(1);
    histBackground->SetMarkerStyle(22);
    histBackground->SetLineColor(2);
    histBackground->SetLineStyle(1);
    histBackground->SetLineWidth(1);
    histBackground->Draw("ESAME");


    textPub = new TPaveText(0.65,0.75,0.92,0.88,"brNDC");
    textPub -> SetTextSize(textSize);
    textPub -> SetTextAlign(22);
    textPub -> SetFillColor(0);
    textPub -> SetTextFont(42);
    textPub -> AddText("p + p #rightarrow p + #pi^{+}#pi^{+}#pi^{-}#pi^{-} + p");
    textPub -> AddText("#sqrt{s} = 510 GeV");
    textPub -> Draw("same");

    textSTAR = new TPaveText(0.65,0.89,0.9,0.95,"brNDC");
    textSTAR -> SetTextSize(textSize);
    textSTAR -> SetTextAlign(22);
    textSTAR -> SetFillColor(0);
    textSTAR -> SetTextFont(62);
    textSTAR->AddText("THIS THESIS");
    textSTAR -> Draw("same");

    leg1 = new TLegend(0.55, 0.65, 0.78, 0.74);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(textSize);
    leg1->SetTextFont(42);
    leg1->AddEntry(histSignal,"Data (unlike-sign pairs)","pe");
    leg1->AddEntry(histBackground,"Data (like-sign pairs)","pe");
    leg1->Draw("same");

    left02 = new TLine(1.5,0.0,1.5,histSignal->GetMaximum()*0.6);
    left02->SetLineStyle(10);
    left02->SetLineColor(1);
    left02->SetLineWidth(4);
    left02->Draw("same");

    cCanvas->Update();
    cCanvas->Write("dcaXY");

    histSignal = (TH1D*)hDCAz[0]->Clone("histSignal");
    histBackground = (TH1D*)hDCAz[1]->Clone("histBackground");

    histSignal->SetTitle(" ; DCA_{z} [cm]; Number of tracks");
    histSignal->SetStats(false);
    histSignal->GetXaxis()->SetTitleFont(42);
    histSignal->GetYaxis()->SetTitleFont(42);
    histSignal->GetXaxis()->SetLabelFont(42);
    histSignal->GetYaxis()->SetLabelFont(42);
    histSignal->GetXaxis()->SetLabelSize(labelSize);
    histSignal->GetYaxis()->SetLabelSize(labelSize);
    histSignal->GetXaxis()->SetTitleSize(labelSize);
    histSignal->GetYaxis()->SetTitleSize(labelSize);
    histSignal->GetXaxis()->SetTitleOffset(1.0);
    histSignal->GetYaxis()->SetTitleOffset(1.3);  
    histSignal->GetYaxis()->SetRangeUser(0.0, TMath::Max(histSignal->GetMaximum(),histBackground->GetMaximum())*1.2);
    histSignal->SetMarkerColor(1);
    histSignal->SetMarkerSize(1);
    histSignal->SetMarkerStyle(20);
    histSignal->SetLineColor(1);
    histSignal->SetLineStyle(1);
    histSignal->SetLineWidth(1);
    histSignal->Draw("E");
    
    gStyle->SetOptStat("");
    gStyle->SetPalette(1);

    histBackground->SetMarkerColor(2);
    histBackground->SetMarkerSize(1);
    histBackground->SetMarkerStyle(22);
    histBackground->SetLineColor(2);
    histBackground->SetLineStyle(1);
    histBackground->SetLineWidth(1);
    histBackground->Draw("ESAME");


    textPub = new TPaveText(0.65,0.75,0.92,0.88,"brNDC");
    textPub -> SetTextSize(textSize);
    textPub -> SetTextAlign(22);
    textPub -> SetFillColor(0);
    textPub -> SetTextFont(42);
    textPub -> AddText("p + p #rightarrow p + #pi^{+}#pi^{+}#pi^{-}#pi^{-} + p");
    textPub -> AddText("#sqrt{s} = 510 GeV");
    textPub -> Draw("same");

    textSTAR = new TPaveText(0.65,0.89,0.9,0.95,"brNDC");
    textSTAR -> SetTextSize(textSize);
    textSTAR -> SetTextAlign(22);
    textSTAR -> SetFillColor(0);
    textSTAR -> SetTextFont(62);
    textSTAR->AddText("THIS THESIS");
    textSTAR -> Draw("same");

    leg1 = new TLegend(0.55, 0.65, 0.78, 0.74);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(textSize);
    leg1->SetTextFont(42);
    leg1->AddEntry(histSignal,"Data (unlike-sign pairs)","pe");
    leg1->AddEntry(histBackground,"Data (like-sign pairs)","pe");
    leg1->Draw("same");

    left02 = new TLine(-1.0,0.0,-1.0,histSignal->GetMaximum()*0.6);
    left02->SetLineStyle(10);
    left02->SetLineColor(1);
    left02->SetLineWidth(4);
    left02->Draw("same");

    left02 = new TLine(1.0,0.0,1.0,histSignal->GetMaximum()*0.6);
    left02->SetLineStyle(10);
    left02->SetLineColor(1);
    left02->SetLineWidth(4);
    left02->Draw("same");

    cCanvas->Update();
    cCanvas->Write("dcaZ");


    histSignal = (TH1D*)hEta[0]->Clone("histSignal");
    histBackground = (TH1D*)hEta[1]->Clone("histBackground");

    histSignal->SetTitle(" ; #eta ; Number of tracks");
    histSignal->SetStats(false);
    histSignal->GetXaxis()->SetTitleFont(42);
    histSignal->GetYaxis()->SetTitleFont(42);
    histSignal->GetXaxis()->SetLabelFont(42);
    histSignal->GetYaxis()->SetLabelFont(42);
    histSignal->GetXaxis()->SetLabelSize(labelSize);
    histSignal->GetYaxis()->SetLabelSize(labelSize);
    histSignal->GetXaxis()->SetTitleSize(labelSize);
    histSignal->GetYaxis()->SetTitleSize(labelSize);
    histSignal->GetXaxis()->SetTitleOffset(1.0);
    histSignal->GetYaxis()->SetTitleOffset(1.3);  
    histSignal->GetYaxis()->SetRangeUser(0.0, TMath::Max(histSignal->GetMaximum(),histBackground->GetMaximum())*1.2);
    histSignal->SetMarkerColor(1);
    histSignal->SetMarkerSize(1);
    histSignal->SetMarkerStyle(20);
    histSignal->SetLineColor(1);
    histSignal->SetLineStyle(1);
    histSignal->SetLineWidth(1);
    histSignal->Draw("E");
    
    gStyle->SetOptStat("");
    gStyle->SetPalette(1);

    histBackground->SetMarkerColor(2);
    histBackground->SetMarkerSize(1);
    histBackground->SetMarkerStyle(22);
    histBackground->SetLineColor(2);
    histBackground->SetLineStyle(1);
    histBackground->SetLineWidth(1);
    histBackground->Draw("ESAME");


    textPub = new TPaveText(0.65,0.75,0.92,0.88,"brNDC");
    textPub -> SetTextSize(textSize);
    textPub -> SetTextAlign(22);
    textPub -> SetFillColor(0);
    textPub -> SetTextFont(42);
    textPub -> AddText("p + p #rightarrow p + #pi^{+}#pi^{+}#pi^{-}#pi^{-} + p");
    textPub -> AddText("#sqrt{s} = 510 GeV");
    textPub -> Draw("same");

    textSTAR = new TPaveText(0.65,0.89,0.9,0.95,"brNDC");
    textSTAR -> SetTextSize(textSize);
    textSTAR -> SetTextAlign(22);
    textSTAR -> SetFillColor(0);
    textSTAR -> SetTextFont(62);
    textSTAR->AddText("THIS THESIS");
    textSTAR -> Draw("same");

    leg1 = new TLegend(0.55, 0.65, 0.78, 0.74);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(textSize);
    leg1->SetTextFont(42);
    leg1->AddEntry(histSignal,"Data (unlike-sign pairs)","pe");
    leg1->AddEntry(histBackground,"Data (like-sign pairs)","pe");
    leg1->Draw("same");

    left02 = new TLine(-0.7,0.0,-0.7,histSignal->GetMaximum()*0.6);
    left02->SetLineStyle(10);
    left02->SetLineColor(1);
    left02->SetLineWidth(4);
    left02->Draw("same");

    left02 = new TLine(0.7,0.0,0.7,histSignal->GetMaximum()*0.6);
    left02->SetLineStyle(10);
    left02->SetLineColor(1);
    left02->SetLineWidth(4);
    left02->Draw("same");

    cCanvas->Update();
    cCanvas->Write("eta");

    histSignal = (TH1D*)hT[0]->Clone("histSignal");
    histBackground = (TH1D*)hT[1]->Clone("histBackground");

    histSignal->SetTitle(" ; t [GeV^{2}]; Number of tracks");
    histSignal->SetStats(false);
    histSignal->GetXaxis()->SetTitleFont(42);
    histSignal->GetYaxis()->SetTitleFont(42);
    histSignal->GetXaxis()->SetLabelFont(42);
    histSignal->GetYaxis()->SetLabelFont(42);
    histSignal->GetXaxis()->SetLabelSize(labelSize);
    histSignal->GetYaxis()->SetLabelSize(labelSize);
    histSignal->GetXaxis()->SetTitleSize(labelSize);
    histSignal->GetYaxis()->SetTitleSize(labelSize);
    histSignal->GetXaxis()->SetTitleOffset(1.0);
    histSignal->GetYaxis()->SetTitleOffset(1.3);  
    histSignal->GetYaxis()->SetRangeUser(0.0, TMath::Max(histSignal->GetMaximum(),histBackground->GetMaximum())*1.2);
    histSignal->SetMarkerColor(1);
    histSignal->SetMarkerSize(1);
    histSignal->SetMarkerStyle(20);
    histSignal->SetLineColor(1);
    histSignal->SetLineStyle(1);
    histSignal->SetLineWidth(1);
    histSignal->Draw("E");
    
    gStyle->SetOptStat("");
    gStyle->SetPalette(1);

    histBackground->SetMarkerColor(2);
    histBackground->SetMarkerSize(1);
    histBackground->SetMarkerStyle(22);
    histBackground->SetLineColor(2);
    histBackground->SetLineStyle(1);
    histBackground->SetLineWidth(1);
    histBackground->Draw("ESAME");


    textPub = new TPaveText(0.65,0.75,0.92,0.88,"brNDC");
    textPub -> SetTextSize(textSize);
    textPub -> SetTextAlign(22);
    textPub -> SetFillColor(0);
    textPub -> SetTextFont(42);
    textPub -> AddText("p + p #rightarrow p + #pi^{+}#pi^{+}#pi^{-}#pi^{-} + p");
    textPub -> AddText("#sqrt{s} = 510 GeV");
    textPub -> Draw("same");

    textSTAR = new TPaveText(0.65,0.89,0.9,0.95,"brNDC");
    textSTAR -> SetTextSize(textSize);
    textSTAR -> SetTextAlign(22);
    textSTAR -> SetFillColor(0);
    textSTAR -> SetTextFont(62);
    textSTAR->AddText("THIS THESIS");
    textSTAR -> Draw("same");

    leg1 = new TLegend(0.55, 0.65, 0.78, 0.74);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(textSize);
    leg1->SetTextFont(42);
    leg1->AddEntry(histSignal,"Data (unlike-sign pairs)","pe");
    leg1->AddEntry(histBackground,"Data (like-sign pairs)","pe");
    leg1->Draw("same");

    left02 = new TLine(-1.0,0.0,-1.0,histSignal->GetMaximum()*0.6);
    left02->SetLineStyle(10);
    left02->SetLineColor(1);
    left02->SetLineWidth(4);
    left02->Draw("same");

    left02 = new TLine(-0.12,0.0,-0.12,histSignal->GetMaximum()*0.6);
    left02->SetLineStyle(10);
    left02->SetLineColor(1);
    left02->SetLineWidth(4);
    left02->Draw("same");

    cCanvas->Update();
    cCanvas->Write("t");

    cCanvas->Close();

}


void PlotCutsFlow()
{
    TFile* anaData = TFile::Open("/home/truhlar/Desktop/STAR/CEP/Analysis/Data/4pi.root", "read");
    TH1I* hCutsData = (TH1I*)anaData -> Get("AnalFlow2");
    TH1I* hCuts = new TH1I("AnalysisFlow", "CutsFlow", 16, 1, 17);
// //////////////////////////////////////////////////////////
// Plot Cuts Flow
    TString Labels[] = { TString("All"), TString("CEP trigger"), TString("2 RP tracks"), TString("Fiducial RP cut"), 
                      TString("4 TPC-TOF tracks"), TString("1 vertex"), TString("Tot. charge 0"), 
                      TString("|z_{vrtx}| < 80 cm"), TString("N_{hits}^{fit} #geq 25"), TString("N_{hits}^{dE/dx} #geq 15"),
                      TString("|DCA(z)| < 1 cm"), TString("DCA(xy) < 1.5 cm"), TString("|#eta| < 0.7"), 
                      TString("0.12< -t < 1.0 GeV^{2}"), TString("p_{T}^{miss} < 0.1 GeV"), TString("PID")};
                      //TString("-1.0 GeV^{2} < t < - 0.12 GeV^{2}")
    //gPad->SetMargin(0.9,0.02,0.1,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
    TCanvas* newCanvas = new TCanvas("newCanvas","newCanvas",950,600);
    gPad->SetMargin(0.08,0.03,0.1,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
    gPad->SetTickx();
    gPad->SetTicky(); 
    gPad->SetLogy();

    for(int iBin = 1; iBin < hCuts->GetNbinsX(); ++iBin)
    {
        hCuts->SetBinContent(iBin,hCutsData->GetBinContent(iBin));
        hCuts->GetXaxis()->SetBinLabel(iBin, Labels[iBin-1]);
    }

    double PIDbinVal = 0;
    PIDbinVal += hInvMassUncorr[0]->GetEntries();

    hCuts->SetBinContent(hCuts->GetNbinsX(),PIDbinVal);
    hCuts->GetXaxis()->SetBinLabel(hCuts->GetNbinsX(), Labels[hCuts->GetNbinsX()-1]);

    hCuts->SetTitle("; ; Number of events");
    hCuts->SetStats(false);
    hCuts->GetXaxis()->SetTitleFont(42);
    hCuts->GetYaxis()->SetTitleFont(42);
    hCuts->GetXaxis()->SetLabelFont(62);
    hCuts->GetYaxis()->SetLabelFont(42);
    hCuts->GetXaxis()->SetLabelSize(labelSize-0.01);
    hCuts->GetYaxis()->SetLabelSize(labelSize);
    hCuts->GetXaxis()->SetTitleSize(labelSize);
    hCuts->GetYaxis()->SetTitleSize(labelSize);
    hCuts->GetXaxis()->SetTitleOffset(0.9);
    hCuts->GetYaxis()->SetTitleOffset(0.85); 
    //hCuts->GetYaxis()->SetRangeUser(4000, 3000000000); 
    hCuts->SetMarkerColor(1);
    hCuts->SetMarkerSize(1.5);
    hCuts->SetMarkerStyle(20);
    hCuts->SetLineColor(1);
    hCuts->SetLineStyle(1);
    hCuts->SetLineWidth(2);
    hCuts->LabelsOption("d");
    hCuts->Draw("");
    gStyle->SetPaintTextFormat("1.2g");
    hCuts->Draw("TEXT15 same");

    TPaveText *textPub = new TPaveText(0.35,0.85,0.7,0.95,"brNDC");
    textPub -> SetTextSize(textSize+0.01);
    textPub -> SetTextAlign(22);
    textPub -> SetFillColor(0);
    textPub -> SetTextFont(42);
    textPub -> SetTextAlign(12);
    textPub -> AddText("p + p #rightarrow p + #pi^{+}#pi^{+}#pi^{-}#pi^{-} + p     #sqrt{s} = 510 GeV");
    textPub -> Draw("same");

    fout->cd();
    newCanvas->Update();
    newCanvas->Write("CutsFlow");
    newCanvas->Close();
}