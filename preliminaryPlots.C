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

const double textSize = 0.05;

void Init();
void ConnectInput(TTree* tree);
void Make(int signal);
void PlotRPPlot();
void PlotPionsPlot();
void PlotKaonsPlot();
void PlotProtonsPlot();
bool ProtonFiducial();

void preliminaryPlots()
{

    TString output = "/home/truhlar/Desktop/STAR/CEP/Analysis/Outputs/preliminaryPlot.root";
    TString input = "/home/truhlar/Desktop/STAR/CEP/Analysis/Data/ppRun17.root";

    TString TPCeffInput[6];
    TPCeffInput[0] = "/home/truhlar/Desktop/STAR/CEP/Analysis/Data/effPionsM.root";
    TPCeffInput[1] = "/home/truhlar/Desktop/STAR/CEP/Analysis/Data/effKaonsM.root";
    TPCeffInput[2] = "/home/truhlar/Desktop/STAR/CEP/Analysis/Data/effProtonsM.root";
    TPCeffInput[3] = "/home/truhlar/Desktop/STAR/CEP/Analysis/Data/effPionsP.root";
    TPCeffInput[4] = "/home/truhlar/Desktop/STAR/CEP/Analysis/Data/effKaonsP.root";
    TPCeffInput[5] = "/home/truhlar/Desktop/STAR/CEP/Analysis/Data/effProtonsP.root";

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

    gStyle->SetLineWidth(2);      //axis line
    gStyle->SetFrameLineWidth(2); //frame line

    PlotRPPlot();
/*    PlotPionsPlot();
    PlotKaonsPlot();
    PlotProtonsPlot();
*/
    gStyle->SetLineWidth(2);      //axis line
    gStyle->SetFrameLineWidth(2); //frame line

    fout->Write();
    fout->Close();
    
}




void Make(int signal)
{
    double effTotal, effTPC, effTOF;
    unsigned int PID;
   // cout<< vertexesZ[0] <<" "<< NhitsFit[0]<<" "<<NhitsFit[1] <<" "<< NhitsDEdx[0]<<" "<<NhitsDEdx[1] <<" "<<DcaZ[0] <<" "<<DcaZ[1] <<" "<<DcaXY[0] <<" "<<DcaXY[1] <<" "<<Eta[0] <<" "<<Eta[1] <<" "<< !fourPiState<<endl;
            

    if(vertexesZ[0] < 80 && vertexesZ[0] > -80 && NhitsFit[0] >=25 && NhitsFit[1] >= 25 && NhitsDEdx[0] >= 15 && NhitsDEdx[1] >= 15 && DcaZ[0] < 1 && DcaZ[0] > -1 && DcaZ[1] < 1 && DcaZ[1] > -1 && DcaXY[0] < 1.5 && DcaXY[1] < 1.5 && Eta[0] > -0.7 && Eta[0] < 0.7 && Eta[1] > -0.7 && Eta[1] < 0.7 &&  t[0] < -0.12 && t[1] < -0.12 && t[0] > -1.0  && t[1] > -1.0 && !fourPiState)
    {
        for (int i = 0; i < nSides; ++i)
        {
            hPxPy[i]->Fill(xCorrelationsRp[i], yCorrelationsRp[i]);
            hPxPy[2]->Fill(xCorrelationsRp[i], yCorrelationsRp[i]);
            //cout<<xCorrelationsRp[i] <<" : "<< yCorrelationsRp[i]<<endl;
        }   

        if(!ProtonFiducial())
            return;

        int combination = 2;
        if(elastic)
            combination = 1;
        effTotal = 1;
        effTOF = 1;
        if(chiPair[Pion] > 9 && chiPair[Kaon] > 9 && chiPair[Proton] < 9 && mSquared > 0.6) // it is... proton!
        {
            for (int iTrack = 0; iTrack < 2; ++iTrack)
            {
                PID = 0;
                if(charge[iTrack] > 0)
                    PID = 3;
                Double_t phiToMC = Phi[iTrack];
                if( phiToMC < 0)
                    phiToMC = 2*3.14159265359 + phiToMC;
                effTPC = hTPCeff[2 + PID]->GetBinContent( hTPCeff[2 + PID]->GetXaxis()->FindBin(phiToMC), hTPCeff[2 + PID]->GetYaxis()->FindBin(transMomentum[iTrack]), hTPCeff[2 + PID]->GetZaxis()->FindBin(Eta[iTrack])); 
                
                effTotal = effTotal*effTPC*effTOF;
            }
            if(effTotal != 0 && transMomentum[0] > 0.4 && transMomentum[1] > 0.4 && (transMomentum[0] < 1.1 || transMomentum[1] < 1.1) )
            {
                hInvMassCorr[Proton][signal][0]->Fill(invMass[Proton], 1/effTotal);
                hInvMassCorr[Proton][signal][combination]->Fill(invMass[Proton], 1/effTotal);
            }
            if(transMomentum[0] > 0.4 && transMomentum[1] > 0.4 && (transMomentum[0] < 1.1 || transMomentum[1] < 1.1) )
            {
                hInvMassUncorr[Proton][signal][0]->Fill(invMass[Proton]);
                hInvMassUncorr[Proton][signal][combination]->Fill(invMass[Proton]);
            }

        }
        else if(chiPair[Pion] > 9 && chiPair[Kaon] < 9 && chiPair[Proton] > 9 && mSquared > 0.15) // it is... kaon!
        {
            for (int iTrack = 0; iTrack < 2; ++iTrack)
            {
                PID = 0;
                if(charge[iTrack] > 0)
                    PID = 3;
                Double_t phiToMC = Phi[iTrack];
                if( phiToMC < 0)
                    phiToMC = 2*3.14159265359 + phiToMC;
                effTPC = hTPCeff[1 + PID]->GetBinContent( hTPCeff[1 + PID]->GetXaxis()->FindBin(phiToMC), hTPCeff[1 + PID]->GetYaxis()->FindBin(transMomentum[iTrack]), hTPCeff[1 + PID]->GetZaxis()->FindBin(Eta[iTrack])); 
                
                effTotal = effTotal*effTPC*effTOF;
            }
            if(effTotal != 0 && transMomentum[0] > 0.3 && transMomentum[1] > 0.3 && (transMomentum[0] < 0.7 || transMomentum[1] < 0.7) )
            {
                hInvMassCorr[Kaon][signal][0]->Fill(invMass[Kaon], 1/effTotal);
                hInvMassCorr[Kaon][signal][combination]->Fill(invMass[Kaon], 1/effTotal);
            }
            if(transMomentum[0] > 0.3 && transMomentum[1] > 0.3 && (transMomentum[0] < 0.7 || transMomentum[1] < 0.7) )
            {
                hInvMassUncorr[Kaon][signal][0]->Fill(invMass[Kaon]);
                hInvMassUncorr[Kaon][signal][combination]->Fill(invMass[Kaon]);
            }
        }
        else if( chiPair[Pion] < 12) // it is... pion!
        {
            for (int iTrack = 0; iTrack < 2; ++iTrack)
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
            if(effTotal != 0 && transMomentum[0] > 0.2 && transMomentum[1] > 0.2)
            {
                hInvMassCorr[Pion][signal][0]->Fill(invMass[Pion], 1/effTotal);
                hInvMassCorr[Pion][signal][combination]->Fill(invMass[Pion], 1/effTotal);
            }
            if(transMomentum[0] > 0.2 && transMomentum[1] > 0.2)
            {
                hInvMassUncorr[Pion][signal][0]->Fill(invMass[Pion]);
                hInvMassUncorr[Pion][signal][combination]->Fill(invMass[Pion]);
            }

        }

    }
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
        hPxPy[i] = new TH2D("pxpy" + sideLabel[i], "px py on " + sideLabel[i], 100, -1.4, 1.4, 100, -1.4, 1.4);
    hPxPy[2] = new TH2D("pxpyEast+West", "px py on East+West", 100, -1.4, 1.4, 100, -1.4, 1.4);
       
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

void PlotRPPlot()
{
    TCanvas *cCanvas2D = new TCanvas("cCanvas2D","cCanvas2D",800,700);
    gPad->SetMargin(0.09,0.13,0.1,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
    gStyle->SetPalette(1);
    gPad->SetTickx();
    gPad->SetTicky();
    gPad->SetLogz(); 


    hPxPy[2]->SetTitle(" ; p_{x} [GeV/c]; p_{y} [GeV/c]");;
    hPxPy[2]->Draw("colz");

    const Int_t n = 300;
    Double_t x[n], y[n];
    Double_t tmp;
    for(int i = 0; i < n; ++i){
        x[i] = 0.194 + (0.354*i)/n;
        //tmp = (x[i] -1.163)*(x[i]-1.163) - 0.464*0.464;
        tmp = (x[i] +0.6)*(x[i]+0.6) - 1.2*1.2;
        y[i] = -sqrt(abs(tmp));
    }
    TGraph* gr = new TGraph(n,x,y);
    gr->SetLineWidth(4);
    gr->Draw("same");


    TLine *left02 = new TLine(-0.3,-0.35,0.548,-0.35);
    left02->SetLineStyle(1);
    left02->SetLineColor(1);
    left02->SetLineWidth(4);
    left02->Draw("same");


    TLine *left01 = new TLine(-0.3,-0.9,-0.3,-0.35);
    left01->SetLineStyle(1);
    left01->SetLineColor(1);
    left01->SetLineWidth(4);
    left01->Draw("same");

    left01 = new TLine(-0.3,-0.9,0.194,-0.9);
    left01->SetLineStyle(1);
    left01->SetLineColor(1);
    left01->SetLineWidth(4);
    left01->Draw("same");          
// UP
    left02 = new TLine(-0.3,0.35,0.548,0.35);
    left02->SetLineStyle(1);
    left02->SetLineColor(1);
    left02->SetLineWidth(4);
    left02->Draw("same");


    left01 = new TLine(-0.3,0.35,-0.3,0.9);
    left01->SetLineStyle(1);
    left01->SetLineColor(1);
    left01->SetLineWidth(4);
    left01->Draw("same");

    left01 = new TLine(-0.3,0.9,0.194,0.9);
    left01->SetLineStyle(1);
    left01->SetLineColor(1);
    left01->SetLineWidth(4);
    left01->Draw("same");

    for(int i = 0; i < n; ++i){
        x[i] = 0.194 + (0.354*i)/n;
        //tmp = (x[i] -1.31)*(x[i]-1.31) - 0.725*0.725;
        tmp = (x[i] +0.6)*(x[i]+0.6) - 1.2*1.2;
        y[i] = sqrt(abs(tmp));
    }
    gr = new TGraph(n,x,y);
    gr->SetLineWidth(4);
    gr->Draw("same");

    cCanvas2D->Update();
    cCanvas2D->Write("hPxPyEast+West");


    cCanvas2D->Close();
}

void PlotPionsPlot()
{
    TCanvas* newCanvas = new TCanvas("newCanvas","newCanvas",800,700);
    gPad->SetMargin(0.11,0.02,0.105,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
    gPad->SetTickx();
    gPad->SetTicky(); 
    gPad->SetLogy(0);

    TString stateLabel = "#pi^{+}#pi^{-}";

    TH1D *hist, *histCompare;
    TString yLabel = "Probability per event / 50 MeV";


    hist = (TH1D*)hInvMassCorr[Pion][0][ElInel]->Clone("hist"); 
    histCompare = (TH1D*)hInvMassCorr[Pion][1][ElInel]->Clone("histCompare");
    
    Double_t scaleFactor; 
    scaleFactor =   1 /hist->Integral();
    hist->Scale(scaleFactor);
    histCompare->Scale(scaleFactor);

    hist->SetTitle(" ; m(" + stateLabel + ") [GeV]; " + yLabel);
    hist->SetStats(false);
    hist->GetXaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetTitleFont(42);
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelFont(42);
    hist->GetXaxis()->SetTitleSize(0.045);
    hist->GetYaxis()->SetTitleSize(0.045);
    hist->GetXaxis()->SetTitleOffset(0.9);
    hist->GetYaxis()->SetTitleOffset(1.1);
    hist->GetYaxis()->SetRangeUser(0, 0.075);   
    hist->SetMarkerColor(1);
    hist->SetMarkerSize(1);
    hist->SetMarkerStyle(20);
    hist->SetLineColor(1);
    hist->SetLineStyle(1);
    hist->SetLineWidth(1);
    hist->Draw("E");

    gStyle->SetOptStat("");
    gStyle->SetPalette(1);
    gStyle->SetLineWidth(2);      //axis line
    gStyle->SetFrameLineWidth(2); //frame line

    histCompare->SetMarkerColor(2);
    histCompare->SetMarkerSize(1);
    histCompare->SetMarkerStyle(22);
    histCompare->SetLineColor(2);
    histCompare->SetLineStyle(1);
    histCompare->SetLineWidth(1);
    histCompare->Draw("ESAME");

    TPaveText *textPub = new TPaveText(0.5,0.9,0.88,0.94,"brNDC");
    textPub -> SetTextSize(textSize);
    textPub -> SetTextAlign(22);
    textPub -> SetFillColor(0);
    textPub -> SetTextFont(42);
    textPub -> AddText("p + p #rightarrow p' + " + stateLabel + " + p'");
    textPub -> AddText("#sqrt{s} = 510 GeV");
    textPub -> Draw("same");

    TPaveText *textSTAR;
    textSTAR = new TPaveText(0.16,0.9,0.5,0.945,"brNDC");
    textSTAR -> SetTextSize(textSize);
    textSTAR -> SetFillColor(0);
    textSTAR -> SetTextFont(72);
    textSTAR -> AddText("STAR Preliminary");
    textSTAR -> Draw("same");

    TPaveText *text;
    text = new TPaveText(0.52,0.8,0.88,0.88,"brNDC");
    text -> SetTextSize(textSize);
    text -> SetFillColor(0);
    text -> SetTextFont(42);
    text -> SetTextAlign(12);
    text -> AddText("p':   0.12 < - t < 1.0 GeV^{2}");
    text -> Draw("same");

    text = new TPaveText(0.52,0.67,0.88,0.78,"brNDC");
    text -> SetTextSize(textSize);
    text -> SetFillColor(0);
    text -> SetTextFont(42);
    text -> SetTextAlign(12);
    text -> AddText("#pi^{+}, #pi^{-}:   p_{T} > 0.2 GeV");
    text -> AddText("            |#eta| < 0.7");
    text -> Draw("same");

    TLegend* leg1 = new TLegend(0.52, 0.54, 0.78, 0.66);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(textSize);
    leg1->SetTextFont(42);
    leg1->AddEntry(hist, "Data (unlike-sign pairs)","pe");
    leg1->AddEntry(histCompare, "Data (like-sign pairs)","pe");
    leg1->Draw("same");

    text = new TPaveText(0.52,0.44,0.88,0.52,"brNDC");
    text -> SetTextSize(textSize);
    text -> SetFillColor(0);
    text -> SetTextFont(42);
    text -> SetTextAlign(12);
    text -> AddText("non-exclusive background not-subtracted");
    text -> Draw("same");

    newCanvas->Update();
    newCanvas->Write("pions");
   // newCanvas->Close(); 

}

void PlotKaonsPlot()
{
    TCanvas* newCanvas = new TCanvas("newCanvas","newCanvas",800,700);
    gPad->SetMargin(0.11,0.02,0.105,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
    gPad->SetTickx();
    gPad->SetTicky(); 
    gPad->SetLogy(0);

    TString stateLabel = "K^{+}K^{-}";

    TH1D *hist, *histCompare;
    TString yLabel = "Probability per event / 50 MeV";


    hist = (TH1D*)hInvMassCorr[Kaon][0][ElInel]->Clone("hist"); 
    histCompare = (TH1D*)hInvMassCorr[Kaon][1][ElInel]->Clone("histCompare");
    
    Double_t scaleFactor; 
    scaleFactor =   1 /hist->Integral();
    hist->Scale(scaleFactor);
    histCompare->Scale(scaleFactor);

    hist->SetTitle(" ; m(" + stateLabel + ") [GeV]; " + yLabel);
    hist->SetStats(false);
    hist->GetXaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetTitleFont(42);
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelFont(42);
    hist->GetXaxis()->SetTitleSize(0.045);
    hist->GetYaxis()->SetTitleSize(0.045);
    hist->GetXaxis()->SetTitleOffset(0.9);
    hist->GetYaxis()->SetTitleOffset(1.1);
    hist->GetYaxis()->SetRangeUser(0, 0.205);   
    hist->SetMarkerColor(1);
    hist->SetMarkerSize(1);
    hist->SetMarkerStyle(20);
    hist->SetLineColor(1);
    hist->SetLineStyle(1);
    hist->SetLineWidth(1);
    hist->Draw("E");

    gStyle->SetOptStat("");
    gStyle->SetPalette(1);
    gStyle->SetLineWidth(2);      //axis line
    gStyle->SetFrameLineWidth(2); //frame line

    histCompare->SetMarkerColor(2);
    histCompare->SetMarkerSize(1);
    histCompare->SetMarkerStyle(22);
    histCompare->SetLineColor(2);
    histCompare->SetLineStyle(1);
    histCompare->SetLineWidth(1);
    histCompare->Draw("ESAME");

    TPaveText *textPub = new TPaveText(0.5,0.9,0.88,0.94,"brNDC");
    textPub -> SetTextSize(textSize);
    textPub -> SetTextAlign(22);
    textPub -> SetFillColor(0);
    textPub -> SetTextFont(42);
    textPub -> AddText("p + p #rightarrow p' + " + stateLabel + " + p'    #sqrt{s} = 510 GeV");
    textPub -> Draw("same");

    TPaveText *textSTAR;
    textSTAR = new TPaveText(0.15,0.9,0.4,0.94,"brNDC");
    textSTAR -> SetTextSize(textSize);
    textSTAR -> SetFillColor(0);
    textSTAR -> SetTextFont(72);
    textSTAR -> AddText("STAR Preliminary");
    textSTAR -> Draw("same");

    TPaveText *text;
    text = new TPaveText(0.52,0.62,0.88,0.8,"brNDC");
    text -> SetTextSize(textSize);
    text -> SetFillColor(0);
    text -> SetTextFont(42);
    text -> SetTextAlign(12);
    text -> AddText("K^{+}, K^{-}:   p_{T} > 0.3 GeV");
    text -> AddText("             min(p_{T}^{+}, p_{T}^{-}) < 0.7 GeV");
    text -> AddText("             |#eta| < 0.7");
    text -> Draw("same");


    text = new TPaveText(0.52,0.8,0.88,0.88,"brNDC");
    text -> SetTextSize(textSize);
    text -> SetFillColor(0);
    text -> SetTextFont(42);
    text -> SetTextAlign(12);
    text -> AddText("p':   0.12 < - t < 1.0 GeV^{2}");
    text -> Draw("same");

    TLegend* leg1 = new TLegend(0.52, 0.5, 0.76, 0.62);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(textSize);
    leg1->SetTextFont(42);
    leg1->AddEntry(hist, "Data (unlike-sign pairs)","pe");
    leg1->AddEntry(histCompare, "Data (like-sign pairs)","pe");
    leg1->Draw("same");

    newCanvas->Update();
    newCanvas->Write("kaons");
   // newCanvas->Close(); 

}

void PlotProtonsPlot()
{
    TCanvas* newCanvas = new TCanvas("newCanvas","newCanvas",800,700);
    gPad->SetMargin(0.11,0.02,0.105,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
    gPad->SetTickx();
    gPad->SetTicky(); 
    gPad->SetLogy(0);

    TString stateLabel = "p#bar{p}";

    TH1D *hist, *histCompare;
    TString yLabel = "Probability per event / 100 MeV";


    hist = (TH1D*)hInvMassCorr[Proton][0][ElInel]->Clone("hist"); 
    histCompare = (TH1D*)hInvMassCorr[Proton][1][ElInel]->Clone("histCompare");
    
    Double_t scaleFactor; 
    scaleFactor =   1 /hist->Integral();
    hist->Scale(scaleFactor);
    histCompare->Scale(scaleFactor);

    hist->SetTitle(" ; m(" + stateLabel + ") [GeV]; " + yLabel);
    hist->SetStats(false);
    hist->GetXaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetTitleFont(42);
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelFont(42);
    hist->GetXaxis()->SetTitleSize(0.045);
    hist->GetYaxis()->SetTitleSize(0.045);
    hist->GetXaxis()->SetTitleOffset(0.9);
    hist->GetYaxis()->SetTitleOffset(1.1);
    hist->GetYaxis()->SetRangeUser(0, 0.29);   
    hist->SetMarkerColor(1);
    hist->SetMarkerSize(1);
    hist->SetMarkerStyle(20);
    hist->SetLineColor(1);
    hist->SetLineStyle(1);
    hist->SetLineWidth(1);
    hist->Draw("E");

    gStyle->SetOptStat("");
    gStyle->SetPalette(1);
    gStyle->SetLineWidth(2);      //axis line
    gStyle->SetFrameLineWidth(2); //frame line

    histCompare->SetMarkerColor(2);
    histCompare->SetMarkerSize(1);
    histCompare->SetMarkerStyle(22);
    histCompare->SetLineColor(2);
    histCompare->SetLineStyle(1);
    histCompare->SetLineWidth(1);
    histCompare->Draw("ESAME");

    TPaveText *textPub = new TPaveText(0.5,0.9,0.88,0.94,"brNDC");
    textPub -> SetTextSize(textSize);
    textPub -> SetTextAlign(22);
    textPub -> SetFillColor(0);
    textPub -> SetTextFont(42);
    textPub -> AddText("p + p #rightarrow p' + " + stateLabel + " + p'    #sqrt{s} = 510 GeV");
    textPub -> Draw("same");

    TPaveText *textSTAR;
    textSTAR = new TPaveText(0.15,0.9,0.4,0.94,"brNDC");
    textSTAR -> SetTextSize(textSize);
    textSTAR -> SetFillColor(0);
    textSTAR -> SetTextFont(72);
    textSTAR -> AddText("STAR Preliminary");
    textSTAR -> Draw("same");

    TPaveText *text;
    text = new TPaveText(0.52,0.62,0.88,0.8,"brNDC");
    text -> SetTextSize(textSize);
    text -> SetFillColor(0);
    text -> SetTextFont(42);
    text -> SetTextAlign(12);
    text -> AddText("p, #bar{p}:^{ }   p_{T} > 0.4 GeV");
    text -> AddText("           min(p_{T}^{+}, p_{T}^{-}) < 1.1 GeV");
    text -> AddText("           |#eta| < 0.7");
    text -> Draw("same");


    text = new TPaveText(0.52,0.8,0.88,0.88,"brNDC");
    text -> SetTextSize(textSize);
    text -> SetFillColor(0);
    text -> SetTextFont(42);
    text -> SetTextAlign(12);
    text -> AddText("p':   0.12 < - t < 1.0 GeV^{2}");
    text -> Draw("same");

    TLegend* leg1 = new TLegend(0.52, 0.5, 0.76, 0.62);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(textSize);
    leg1->SetTextFont(42);
    leg1->AddEntry(hist, "Data (unlike-sign pairs)","pe");
    leg1->AddEntry(histCompare, "Data (like-sign pairs)","pe");
    leg1->Draw("same");

    newCanvas->Update();
    newCanvas->Write("protons");
    newCanvas->Close(); 

}



bool ProtonFiducial()
{
    return true;
    for (int i = 0; i < nSides; ++i)
    {
        if(abs(yCorrelationsRp[i]) < 0.9 && abs(yCorrelationsRp[i]) > 0.35  && xCorrelationsRp[i] > -0.3 )
            return true;
    }

}




/*
void PlotRPPlot()
{
    TCanvas *cCanvas2D = new TCanvas("cCanvas2D","cCanvas2D",800,700);
    gPad->SetMargin(0.09,0.13,0.1,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
    gStyle->SetPalette(1);
    gPad->SetTickx();
    gPad->SetTicky();
    gPad->SetLogz(); 

    for (int i = 0; i < nSides; ++i)
    {
        hPxPy[i]->SetTitle(" ; p_{x} [GeV/c]; p_{y} [GeV/c]");;
        hPxPy[i]->Draw("colz");

        if(i == East)
        {
            // Down....
            const Int_t n = 300;
            Double_t x[n], y[n];
            Double_t tmp;
            for(int i = 0; i < n; ++i){
                x[i] = 0.153 + (0.3*i)/n;
                tmp = (x[i] -1.424)*(x[i]-1.424) - 0.88*0.88;
                y[i] = -sqrt(abs(tmp));
            }
            TGraph* gr = new TGraph(n,x,y);
            gr->SetLineWidth(4);
            gr->Draw("same");


            TLine *left02 = new TLine(-0.28,-0.413,0.453,-0.413);
            left02->SetLineStyle(1);
            left02->SetLineColor(1);
            left02->SetLineWidth(4);
            left02->Draw("same");


            TLine *left01 = new TLine(-0.28,-0.9,-0.28,-0.413);
            left01->SetLineStyle(1);
            left01->SetLineColor(1);
            left01->SetLineWidth(4);
            left01->Draw("same");

            left01 = new TLine(-0.28,-0.917,0.153,-0.917);
            left01->SetLineStyle(1);
            left01->SetLineColor(1);
            left01->SetLineWidth(4);
            left01->Draw("same");

    // Up....

            left02 = new TLine(-0.3,0.36,0.5,0.36);
            left02->SetLineStyle(1);
            left02->SetLineColor(1);
            left02->SetLineWidth(4);
            left02->Draw("same");


            left01 = new TLine(-0.3,0.36,-0.3,0.84);
            left01->SetLineStyle(1);
            left01->SetLineColor(1);
            left01->SetLineWidth(4);
            left01->Draw("same");

            left01 = new TLine(-0.3,0.84,0.2,0.84);
            left01->SetLineStyle(1);
            left01->SetLineColor(1);
            left01->SetLineWidth(4);
            left01->Draw("same");

            for(int i = 0; i < n; ++i){
                x[i] = 0.2 + (0.3*i)/n;
                tmp = (x[i] -1.31)*(x[i]-1.31) - 0.725*0.725;
                y[i] = sqrt(abs(tmp));
            }
            gr = new TGraph(n,x,y);
            gr->SetLineWidth(4);
            gr->Draw("same");
        }else if(i == West)
        {
            // Down....
            const Int_t n = 300;
            Double_t x[n], y[n];
            Double_t tmp;
            for(int i = 0; i < n; ++i){
                x[i] = 0.15 + (0.4*i)/n;
                tmp = (x[i] -1.163)*(x[i]-1.163) - 0.464*0.464;
                y[i] = -sqrt(abs(tmp));
            }
            TGraph* gr = new TGraph(n,x,y);
            gr->SetLineWidth(4);
            gr->Draw("same");


            TLine *left02 = new TLine(-0.25,-0.4,0.55,-0.4);
            left02->SetLineStyle(1);
            left02->SetLineColor(1);
            left02->SetLineWidth(4);
            left02->Draw("same");


            TLine *left01 = new TLine(-0.25,-0.9,-0.25,-0.4);
            left01->SetLineStyle(1);
            left01->SetLineColor(1);
            left01->SetLineWidth(4);
            left01->Draw("same");

            left01 = new TLine(-0.25,-0.9,0.15,-0.9);
            left01->SetLineStyle(1);
            left01->SetLineColor(1);
            left01->SetLineWidth(4);
            left01->Draw("same");

    // Up....

            left02 = new TLine(-0.22,0.37,0.44,0.37);
            left02->SetLineStyle(1);
            left02->SetLineColor(1);
            left02->SetLineWidth(4);
            left02->Draw("same");


            left01 = new TLine(-0.22,0.37,-0.22,0.8);
            left01->SetLineStyle(1);
            left01->SetLineColor(1);
            left01->SetLineWidth(4);
            left01->Draw("same");

            left01 = new TLine(-0.22,0.8,0.2,0.8);
            left01->SetLineStyle(1);
            left01->SetLineColor(1);
            left01->SetLineWidth(4);
            left01->Draw("same");

            for(int i = 0; i < n; ++i){
                x[i] = 0.2 + (0.24*i)/n;
                tmp = (x[i] -1.368)*(x[i]-1.368) - 0.851*0.851;
                y[i] = sqrt(abs(tmp));
            }
            gr = new TGraph(n,x,y);
            gr->SetLineWidth(4);
            gr->Draw("same");            
        }

        cCanvas2D->Update();
        cCanvas2D->Write("hPxPy"+sideLabel[i]);
    }

    hPxPy[2]->SetTitle(" ; p_{x} [GeV/c]; p_{y} [GeV/c]");;
    hPxPy[2]->Draw("colz");

    const Int_t n = 300;
    Double_t x[n], y[n];
    Double_t tmp;
    for(int i = 0; i < n; ++i){
        x[i] = 0.15 + (0.4*i)/n;
        //tmp = (x[i] -1.163)*(x[i]-1.163) - 0.464*0.464;
        tmp = (x[i] +0.6)*(x[i]+0.6) - 1.2*1.2;
        y[i] = -sqrt(abs(tmp));
    }
    TGraph* gr = new TGraph(n,x,y);
    gr->SetLineWidth(4);
    gr->Draw("same");


    TLine *left02 = new TLine(-0.28,-0.4,0.55,-0.4);
    left02->SetLineStyle(1);
    left02->SetLineColor(1);
    left02->SetLineWidth(4);
    left02->Draw("same");


    TLine *left01 = new TLine(-0.28,-0.9,-0.28,-0.4);
    left01->SetLineStyle(1);
    left01->SetLineColor(1);
    left01->SetLineWidth(4);
    left01->Draw("same");

    left01 = new TLine(-0.28,-0.9,0.15,-0.9);
    left01->SetLineStyle(1);
    left01->SetLineColor(1);
    left01->SetLineWidth(4);
    left01->Draw("same");          
// UP
    left02 = new TLine(-0.3,0.36,0.5,0.36);
    left02->SetLineStyle(1);
    left02->SetLineColor(1);
    left02->SetLineWidth(4);
    left02->Draw("same");


    left01 = new TLine(-0.3,0.36,-0.3,0.84);
    left01->SetLineStyle(1);
    left01->SetLineColor(1);
    left01->SetLineWidth(4);
    left01->Draw("same");

    left01 = new TLine(-0.3,0.84,0.2,0.84);
    left01->SetLineStyle(1);
    left01->SetLineColor(1);
    left01->SetLineWidth(4);
    left01->Draw("same");

    for(int i = 0; i < n; ++i){
        x[i] = 0.2 + (0.3*i)/n;
        //tmp = (x[i] -1.31)*(x[i]-1.31) - 0.725*0.725;
        tmp = (x[i] +0.6)*(x[i]+0.6) - 1.2*1.2;
        y[i] = sqrt(abs(tmp));
    }
    gr = new TGraph(n,x,y);
    gr->SetLineWidth(4);
    gr->Draw("same");

    cCanvas2D->Update();
    cCanvas2D->Write("hPxPyEast+West");


    cCanvas2D->Close();
} */