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

const bool OldProduction = false;
TString label;

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
void Make(int signal);
void PlotCutsFlow();
void PlotRPPlot();
void PlotPionsPlot();
void PlotKaonsPlot();
void PlotProtonsPlot();
bool ProtonFiducial();
void RunGraniitti();

void SetGraphStyle(TH1* hist);
void SetTextStyle(TPaveText* text);
void SetLegendStyle(TLegend* leg1);
void SetMarkerStyle(TH1* hist);
void SetLineStyle(TLine* line);


int protonsInside, protonsTotal;

void Graniitti()
{
    TString output;
    TString input;
    if(OldProduction){
        label = "Preliminary";
        output = "/home/truhlar/Desktop/STAR/CEP/Analysis/Outputs/ToBeShown.root"; 
        input = "/home/truhlar/Desktop/STAR/CEP/Analysis/Data/ppRun17.root";
    }else{
        label = "Internal";
        output = "/home/truhlar/Desktop/STAR/CEP/Analysis/Outputs/graniitti.root";
        input = "/home/truhlar/Desktop/STAR/CEP/Analysis/Data/P20ic.root";
    }

    TString graniitti_input = "/home/truhlar/Desktop/STAR/Graniitti_new/GRANIITTI/output/RootFiles/510/510.root";

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

    graniitti = TFile::Open(graniitti_input, "read");
    if (!graniitti)
    {
        cout<<"Error: cannot open "<<graniitti_input<<endl;
        return;
    }

    fout = new TFile(output,"RECREATE");
    Init(); // Preparing histograms 

    tree = dynamic_cast<TTree*>( data->Get("recTree") );
    treeBack = dynamic_cast<TTree*>( data->Get("Background") );

    gStyle->SetLineWidth(2);      //axis line
    gStyle->SetFrameLineWidth(2); //frame line

    protonsInside = 0;
    protonsTotal = 0;
    ConnectInput(tree);
    for(Long64_t iev=0; iev<tree->GetEntries(); ++iev)
    { //get the event
        tree->GetEntry(iev); 
        Make(0);
    } 
    cout<<"Protons: "<<protonsInside<<" / "<<protonsTotal<<endl;
    ConnectInput(treeBack);
    for(Long64_t iev=0; iev<treeBack->GetEntries(); ++iev)
    { //get the event
        treeBack->GetEntry(iev); 
        Make(1);
    }

    RunGraniitti();
    PlotRPPlot();

    PlotPionsPlot();
    PlotKaonsPlot();
    PlotProtonsPlot();

    PlotCutsFlow();

    int a;
    cin>>a;

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

        double deltaPhi = TMath::Abs(phiRp[East] - phiRp[West])*convertToDegree;
        //cout<< deltaPhi<<" = "<< phiRp[East] << " - " << phiRp[West]<<endl;
        if(deltaPhi > 180)
            deltaPhi = 360 - deltaPhi;

        hDeltaPhi[0]->Fill(deltaPhi);
        if(elastic)
            hDeltaPhi[1]->Fill(deltaPhi);
        else
            hDeltaPhi[2]->Fill(deltaPhi);

        protonsTotal++;
        if(!ProtonFiducial())
            return;
        protonsInside++;
        int combination = 2;
        //if(elastic)
        if(deltaPhi > 90)
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

void PlotRPPlot()
{
    TCanvas *cCanvas2D = new TCanvas("cCanvas2D","cCanvas2D",800,700);
    gPad->SetMargin(0.09,0.16,0.11,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
    gStyle->SetPalette(1);
    gPad->SetTickx();
    gPad->SetTicky();
    gPad->SetLogz(); 

    gStyle->SetOptStat("");
    gStyle->SetPalette(1);

    hPxPy[2]->SetTitle(" ; p_{x} [GeV]; p_{y} [GeV]; Events / (10 MeV #times 10 MeV)");
    hPxPy[2]->GetXaxis()->SetTitleFont(fontStyle);
    hPxPy[2]->GetXaxis()->SetTitleFont(fontStyle);
    hPxPy[2]->GetZaxis()->SetTitleFont(42);
    hPxPy[2]->GetXaxis()->SetLabelFont(fontStyle);
    hPxPy[2]->GetYaxis()->SetLabelFont(fontStyle);
    hPxPy[2]->GetZaxis()->SetLabelFont(42);
    hPxPy[2]->GetXaxis()->SetLabelSize(labelSize);
    hPxPy[2]->GetYaxis()->SetLabelSize(labelSize);
    hPxPy[2]->GetZaxis()->SetLabelSize(labelSize);
    hPxPy[2]->GetXaxis()->SetTitleSize(labelSize);
    hPxPy[2]->GetYaxis()->SetTitleSize(labelSize);
    hPxPy[2]->GetZaxis()->SetTitleSize(labelSize);
    hPxPy[2]->GetXaxis()->SetTitleOffset(1.0);
    hPxPy[2]->GetYaxis()->SetTitleOffset(0.7);
    hPxPy[2]->GetZaxis()->SetTitleOffset(1.1);
    //hPxPy[2]->GetYaxis()->SetRangeUser(-1.1, 1.4); 
    //hPxPy[2]->GetXaxis()->SetRangeUser(-0.8, 0.8); 
    hPxPy[2]->SetStats(false);
    hPxPy[2]->SetLineColor(1);
    hPxPy[2]->SetLineWidth(4);
    hPxPy[2]->Draw("colz");


    TPaveText *textPub = new TPaveText(0.12,0.9,0.32,0.95,"brNDC");
    SetTextStyle(textPub);
    textPub -> SetTextAlign(12);
    textPub -> AddText("p+p #rightarrow p+h^{+}h^{-}+p");
    textPub -> Draw("same");
    textPub = new TPaveText(0.62,0.9,0.79,0.95,"brNDC");
    SetTextStyle(textPub);
    textPub -> SetTextAlign(12);
    textPub -> AddText("#sqrt{s}=510 GeV");
    textPub -> Draw("same");


    TLegend* leg1 = new TLegend(0.15, 0.52, 0.7, 0.57);
    SetLegendStyle(leg1);
    leg1->SetMargin(0.05);
    leg1->AddEntry(hPxPy[2], "Forward proton fiducial region","l");
    leg1->Draw("same");

    const Int_t n = 100;
    Double_t x[n], y[n];
    Double_t tmp;
    for(int i = 0; i < n; ++i){
        x[i] = 0.175 + (0.27*i)/n;
        //tmp = (x[i] -1.163)*(x[i]-1.163) - 0.464*0.464;
        tmp = (x[i] +0.6)*(x[i]+0.6) - 1.25;
        y[i] = -sqrt(abs(tmp));
    }
    TGraph* gr = new TGraph(n,x,y);
    gr->SetLineWidth(4);
    gr->Draw("same");


    TLine *left02 = new TLine(-0.27,-0.4,0.445,-0.4);
    SetLineStyle(left02);
    left02->Draw("same");


    TLine *left01 = new TLine(-0.27,-0.8,-0.27,-0.4);
    SetLineStyle(left01);
    left01->Draw("same");

    left01 = new TLine(-0.27,-0.8,0.185,-0.8);
    SetLineStyle(left01);
    left01->Draw("same");          
// UP
    left02 = new TLine(-0.27,0.4,0.445,0.4);
    SetLineStyle(left02);
    left02->Draw("same");


    left01 = new TLine(-0.27,0.4,-0.27,0.8);
    SetLineStyle(left01);
    left01->Draw("same");

    left01 = new TLine(-0.27,0.8,0.185,0.8);
    SetLineStyle(left01);
    left01->Draw("same");

    for(int i = 0; i < n; ++i){
        x[i] = 0.175 + (0.27*i)/n;
        //tmp = (x[i] -1.31)*(x[i]-1.31) - 0.725*0.725;
        tmp = (x[i] +0.6)*(x[i]+0.6) - 1.25;
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
    TCanvas* newCanvas = new TCanvas("pions","pions",800,700);
    gPad->SetMargin(0.13,0.03,0.105,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
    gPad->SetTickx();
    gPad->SetTicky(); 
    gPad->SetLogy(0);

    TString stateLabel = "#pi^{+}#pi^{-}";

    TH1D *hist, *histCompare, *histGraniitti;
    TString yLabel = "Probability per event / 50 MeV";


    hist = (TH1D*)hInvMassCorr[Pion][0][ElInel]->Clone("hist"); 
    histCompare = (TH1D*)hInvMassCorr[Pion][1][ElInel]->Clone("histCompare");


    TH1D* histGraniittiEl, *histGraniittiInel;
    histGraniittiInel = (TH1D*)hInvMassGran[Inel][Pion]->Clone("histGraniittiInel");
    histGraniittiEl = (TH1D*)hInvMassGran[El][Pion]->Clone("histGraniittiEl");

    Double_t inel = hInvMassCorr[Pion][0][Inel]->Integral(); 
    Double_t el = hInvMassCorr[Pion][0][El]->Integral();
    //cout<<"Elastic / inelastic = "<< el<< " / "<< inel << "\n";
    histGraniittiInel->Scale(inel/((el+inel)*histGraniittiInel->Integral()));
    histGraniittiEl->Scale(el/((el+inel)*histGraniittiEl->Integral()));
    //cout<<"Scaled Elastic / inelastic = "<< histGraniittiEl->Integral() << " / "<< histGraniittiInel->Integral() << "\n";
    
    histGraniitti = (TH1D*)histGraniittiInel->Clone("histGraniitti");
    histGraniitti->Add(histGraniittiEl, 1.0);
    //cout<<"Sum: "<< histGraniitti->Integral() << "\n";

    Double_t scaleFactor; 
    scaleFactor =   1 /hist->Integral();
    hist->Scale(scaleFactor);
    histCompare->Scale(scaleFactor);
    cout<<"Graniitti integral: "<<histGraniitti->Integral()<<endl;

    hist->SetTitle(" ; m(" + stateLabel + ") [GeV]; " + yLabel);
    SetGraphStyle(hist);          
    hist->GetYaxis()->SetRangeUser(0, 0.092);   
    hist->Draw("hist E");

    gStyle->SetOptStat("");
    gStyle->SetPalette(1);

    SetMarkerStyle(histCompare);
    histCompare->SetMarkerColor(2);
    histCompare->SetMarkerStyle(22);
    histCompare->SetLineColor(2);
    histCompare->Draw("same hist E");

    SetMarkerStyle(histGraniitti);
    histGraniitti->Draw("same hist E");


    TPaveText *textSTAR;
    textSTAR = new TPaveText(0.17,0.89,0.33,0.95,"brNDC");
    SetTextStyle(textSTAR);
    textSTAR -> SetTextSize(textSize+0.02);
    textSTAR -> SetTextFont(72);
    textSTAR -> AddText("STAR");
    textSTAR -> Draw("same");
    textSTAR = new TPaveText(0.17,0.84,0.33,0.89,"brNDC");
    SetTextStyle(textSTAR);
    textSTAR -> SetTextFont(52);
    textSTAR -> AddText(label);
    textSTAR -> Draw("same");

    TPaveText *textPub = new TPaveText(0.35,0.88,0.88,0.95,"brNDC");
    SetTextStyle(textPub);
    textPub -> SetTextAlign(12);
    textPub -> AddText("p + p #rightarrow p + " + stateLabel + " + p       #sqrt{s} = 510 GeV");
    //textPub -> AddText("#sqrt{s} = 510 GeV");
    textPub -> Draw("same");


    TPaveText *text;
    text = new TPaveText(0.27,0.69,0.52,0.84,"brNDC");
    SetTextStyle(text);
    text -> AddText("#pi^{+}, #pi^{-} kinematics:");   
    text -> AddText("p_{T} > 0.2 GeV");
    text -> AddText("|#eta| < 0.7");
    text -> Draw("same");

    text = new TPaveText(0.54,0.63,0.88,0.84,"brNDC");
    SetTextStyle(text);
    text -> SetTextAlign(12);
    text -> AddText("Forward proton kinematics:");
    text -> AddText("(p_{x} + 0.6)^{2} + p_{y}^{2} < 1.25 GeV^{2}");
    text -> AddText("0.4 GeV < |p_{y}| < 0.8 GeV");
    text -> AddText("p_{x} > -0.27 GeV");
    text -> Draw("same");


    TLegend *leg1 = new TLegend(0.54, 0.44, 0.87, 0.59);
    SetLegendStyle(leg1);
    leg1->AddEntry(hist, "Data (unlike-sign pairs)","ple");
    leg1->AddEntry(histCompare, "Data (like-sign pairs)","ple");
    leg1->AddEntry(histGraniitti, "Graniitti","ple");
    leg1->Draw("same");

    text = new TPaveText(0.51,0.17,0.93,0.39,"brNDC");
    SetTextStyle(text);
    text -> AddText("Statistical errors only");
    text -> AddText("Acceptance corrected");
    text -> AddText("Not efficiency corrected");
    text -> AddText("Not background subtracted");
    text -> Draw("same");

    newCanvas->Update();
    newCanvas->Write("pions");
    //newCanvas->Close(); 
/////////////   El + Inel /////////////////////////////////////
    newCanvas = new TCanvas("pionsInel","pionsInel",800,700);
    gPad->SetMargin(0.13,0.03,0.105,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
    gPad->SetTickx();
    gPad->SetTicky(); 
    gPad->SetLogy(0);
    hist = (TH1D*)hInvMassCorr[Pion][0][Inel]->Clone("hist"); 
    histCompare = (TH1D*)hInvMassCorr[Pion][0][El]->Clone("histCompare");
    histGraniitti = (TH1D*)hInvMassGran[Inel][Pion]->Clone("histGraniittiInel");
    cout<<"Here histGraniittiInel with "<< hInvMassGran[Inel][Pion]->GetEntries()<<endl;
    histGraniittiEl = (TH1D*)hInvMassGran[El][Pion]->Clone("histGraniittiEl");
    cout<<"Here histGraniittiEl with "<< hInvMassGran[El][Pion]->GetEntries()<<endl;
    cout<<"_________________"<<endl;
    cout<<"Ratio-data (inel/el): "<<hist->Integral()<<" / "<<histCompare->Integral()<<" = "<<hist->Integral() / histCompare->Integral() <<endl;
    cout<<"Ratio-graniitti (inel/el): "<<histGraniitti->Integral()<<" / "<<histGraniittiEl->Integral()<<" = "<<histGraniitti->Integral() / histGraniittiEl->Integral() <<endl;
    cout<<"_________________"<<endl;

    scaleFactor =   1 /hist->Integral();
    hist->Scale(scaleFactor);
    scaleFactor =   1 /histCompare->Integral();
    histCompare->Scale(scaleFactor);
    scaleFactor =   1 /histGraniitti->Integral();
    histGraniitti->Scale(scaleFactor);
    scaleFactor =   1 /histGraniittiEl->Integral();
    histGraniittiEl->Scale(scaleFactor);

///////////////////////////// Inelastic ///////////////////////////////////////////
    hist->SetTitle(" ; m(" + stateLabel + ") [GeV]; " + yLabel);
    SetGraphStyle(hist);
    hist->GetYaxis()->SetRangeUser(0, 0.17);   
    hist->Draw("hist E");

    gStyle->SetOptStat("");
    gStyle->SetPalette(1);


    SetMarkerStyle(histGraniitti);
    histGraniitti->Draw("SAME hist E");


    textSTAR = new TPaveText(0.17,0.89,0.33,0.95,"brNDC");
    SetTextStyle(textSTAR);
    textSTAR -> SetTextSize(textSize+0.02);
    textSTAR -> SetTextFont(72);
    textSTAR -> AddText("STAR");
    textSTAR -> Draw("same");
    textSTAR = new TPaveText(0.17,0.84,0.33,0.89,"brNDC");
    SetTextStyle(textSTAR);
    textSTAR -> SetTextFont(52);
    textSTAR -> AddText(label);
    textSTAR -> Draw("same");

    textPub = new TPaveText(0.35,0.88,0.88,0.95,"brNDC");
    SetTextStyle(textPub);
    textPub -> SetTextAlign(12);
    textPub -> AddText("p + p #rightarrow p + " + stateLabel + " + p       #sqrt{s} = 510 GeV");
    //textPub -> AddText("#sqrt{s} = 510 GeV");
    textPub -> Draw("same");

    text = new TPaveText(0.27,0.68,0.52,0.83,"brNDC");
    SetTextStyle(text);
    text -> SetTextAlign(32);
    text -> AddText("#pi^{+}, #pi^{-} kinematics:");   
    text -> AddText("p_{T} > 0.2 GeV");
    text -> AddText("|#eta| < 0.7");
    text -> Draw("same");

    text = new TPaveText(0.54,0.62,0.88,0.83,"brNDC");
    SetTextStyle(text);
    text -> SetTextAlign(12);
    text -> AddText("Forward proton kinematics:");
    text -> AddText("(p_{x} + 0.6)^{2} + p_{y}^{2} < 1.25 GeV^{2}");
    text -> AddText("0.4 GeV < |p_{y}| < 0.8 GeV");
    text -> AddText("p_{x} > -0.27 GeV");
    text -> Draw("same");


    leg1 = new TLegend(0.45, 0.48, 0.78, 0.58);
    SetLegendStyle(leg1);
    leg1->AddEntry(hist, "Data, #Delta#varphi < 90^{#circ} (unlike-sign pairs)","ple");
    leg1->AddEntry(histGraniitti, "Graniitti, #Delta#varphi < 90^{#circ}","ple");
    leg1->Draw("same");

    text = new TPaveText(0.6,0.22,0.88,0.44,"brNDC");
    SetTextStyle(text);
    text -> AddText("Statistical errors only");
    text -> AddText("Acceptance corrected");
    text -> AddText("Not efficiency corrected");
    text -> AddText("Not background subtracted");
    text -> Draw("same");

    newCanvas->Update();
    newCanvas->Write("pionsInel");
    //newCanvas->Close();
//////////////////////////////////////////// Elastic ///////////////////////////////
    newCanvas = new TCanvas("pionsEl","pionsEl",800,700);
    histCompare->SetTitle(" ; m(" + stateLabel + ") [GeV]; " + yLabel);
    SetGraphStyle(histCompare);
    histCompare->GetYaxis()->SetRangeUser(0, 0.12);   
    histCompare->Draw("hist E");
    gStyle->SetOptStat("");
    gStyle->SetPalette(1);

    SetMarkerStyle(histGraniittiEl);
    histGraniittiEl->Draw("SAME hist E");

    textSTAR = new TPaveText(0.17,0.89,0.33,0.95,"brNDC");
    SetTextStyle(textSTAR);
    textSTAR -> SetTextSize(textSize+0.02);
    textSTAR -> SetTextFont(72);
    textSTAR -> AddText("STAR");
    textSTAR -> Draw("same");
    textSTAR = new TPaveText(0.17,0.84,0.33,0.89,"brNDC");
    SetTextStyle(textSTAR);
    textSTAR -> SetTextFont(52);
    textSTAR -> AddText(label);
    textSTAR -> Draw("same");

    textPub = new TPaveText(0.35,0.88,0.88,0.95,"brNDC");
    SetTextStyle(textPub);
    textPub -> SetTextAlign(12);
    textPub -> AddText("p + p #rightarrow p + " + stateLabel + " + p       #sqrt{s} = 510 GeV");
    //textPub -> AddText("#sqrt{s} = 510 GeV");
    textPub -> Draw("same");

    text = new TPaveText(0.27,0.68,0.52,0.83,"brNDC");
    SetTextStyle(text);
    text -> SetTextAlign(32);
    text -> AddText("#pi^{+}, #pi^{-} kinematics:");   
    text -> AddText("p_{T} > 0.2 GeV");
    text -> AddText("|#eta| < 0.7");
    text -> Draw("same");

    text = new TPaveText(0.54,0.62,0.88,0.83,"brNDC");
    SetTextStyle(text);
    text -> SetTextAlign(12);
    text -> AddText("Forward proton kinematics:");
    text -> AddText("(p_{x} + 0.6)^{2} + p_{y}^{2} < 1.25 GeV^{2}");
    text -> AddText("0.4 GeV < |p_{y}| < 0.8 GeV");
    text -> AddText("p_{x} > -0.27 GeV");
    text -> Draw("same");


    leg1 = new TLegend(0.45, 0.48, 0.78, 0.58);
    SetLegendStyle(leg1);
    leg1->AddEntry(histCompare, "Data, #Delta#varphi > 90^{#circ} (unlike-sign pairs)","ple");
    leg1->AddEntry(histGraniittiEl, "Graniitti, #Delta#varphi > 90^{#circ}","ple");
    leg1->Draw("same");

    text = new TPaveText(0.6,0.22,0.88,0.44,"brNDC");
    SetTextStyle(text);
    text -> AddText("Statistical errors only");
    text -> AddText("Acceptance corrected");
    text -> AddText("Not efficiency corrected");
    text -> AddText("Not background subtracted");
    text -> Draw("same");

    newCanvas->Update();
    newCanvas->Write("pionsEl");
    //newCanvas->Close();
}

void PlotKaonsPlot()
{
    TCanvas* newCanvas = new TCanvas("kaons","kaons",800,700);
    gPad->SetMargin(0.13,0.02,0.105,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
    gPad->SetTickx();
    gPad->SetTicky(); 
    gPad->SetLogy(0);

    TString stateLabel = "K^{+}K^{-}";

    TH1D *hist, *histCompare;
    TString yLabel = "Probability per event / 50 MeV";

    hist = (TH1D*)hInvMassCorr[Kaon][0][ElInel]->Clone("hist"); 
    histCompare = (TH1D*)hInvMassCorr[Kaon][1][ElInel]->Clone("histCompare");

    TH1D* histGraniittiEl, *histGraniittiInel, *histGraniitti;
    histGraniittiInel = (TH1D*)hInvMassGran[Inel][Kaon]->Clone("histGraniittiInel");
    histGraniittiEl = (TH1D*)hInvMassGran[El][Kaon]->Clone("histGraniittiEl");

    Double_t inel = hInvMassCorr[Kaon][0][Inel]->Integral(); 
    Double_t el = hInvMassCorr[Kaon][0][El]->Integral();

    histGraniittiInel->Scale(inel/((el+inel)*histGraniittiInel->Integral()));
    histGraniittiEl->Scale(el/((el+inel)*histGraniittiEl->Integral()));

    histGraniitti = (TH1D*)histGraniittiInel->Clone("histGraniitti");
    histGraniitti->Add(histGraniittiEl, 1.0);

 
    Double_t scaleFactor; 
    scaleFactor =   1 /hist->Integral();
    hist->Scale(scaleFactor);
    histCompare->Scale(scaleFactor);

    hist->SetTitle(" ; m(" + stateLabel + ") [GeV]; " + yLabel);
    SetGraphStyle(hist);
    hist->GetYaxis()->SetRangeUser(0, 0.355);   
    hist->Draw("hist E");

    gStyle->SetOptStat("");
    gStyle->SetPalette(1);

    SetMarkerStyle(histCompare);
    histCompare->SetMarkerColor(2);
    histCompare->SetMarkerStyle(22);
    histCompare->SetLineColor(2);
    histCompare->Draw("SAME hist E");

    SetMarkerStyle(histGraniitti);
    histGraniitti->Draw("same hist E");

    TPaveText *textSTAR;
    textSTAR = new TPaveText(0.17,0.89,0.33,0.95,"brNDC");
    SetTextStyle(textSTAR);
    textSTAR -> SetTextSize(textSize+0.02);
    textSTAR -> SetTextFont(72);
    textSTAR -> AddText("STAR");
    textSTAR -> Draw("same");
    textSTAR = new TPaveText(0.17,0.84,0.33,0.89,"brNDC");
    SetTextStyle(textSTAR);
    textSTAR -> SetTextFont(52);
    textSTAR -> AddText(label);
    textSTAR -> Draw("same");

    TPaveText *textPub = new TPaveText(0.35,0.88,0.88,0.95,"brNDC");
    SetTextStyle(textPub);
    textPub -> SetTextAlign(12);
    textPub -> AddText("p + p #rightarrow p + " + stateLabel + " + p       #sqrt{s} = 510 GeV");
    //textPub -> AddText("#sqrt{s} = 510 GeV");
    textPub -> Draw("same");


    TPaveText *text;
    text = new TPaveText(0.26,0.62,0.55,0.81,"brNDC");
    SetTextStyle(text);
    text -> SetTextAlign(31);
    text -> AddText("K^{+}, K^{-} kinematics:");
    text -> AddText("min(p_{T}^{+},p_{T}^{-}) < 0.7 GeV");   
    text -> AddText("p_{T} > 0.3 GeV");
    text -> AddText("|#eta| < 0.7");
    text -> Draw("same");

    text = new TPaveText(0.55,0.62,0.89,0.81,"brNDC");
    SetTextStyle(text);
    text -> SetTextAlign(11);
    text -> AddText("Forward proton kinematics:");
    text -> AddText("(p_{x} + 0.6)^{2} + p_{y}^{2} < 1.25 GeV^{2}");
    text -> AddText("0.4 GeV < |p_{y}| < 0.8 GeV");
    text -> AddText("p_{x} > -0.27 GeV");
    text -> Draw("same");


    TLegend* leg1 = new TLegend(0.54, 0.44, 0.87, 0.59);
    SetLegendStyle(leg1);
    leg1->AddEntry(hist, "Data (unlike-sign pairs)","ple");
    leg1->AddEntry(histCompare, "Data (like-sign pairs)","ple");
    leg1->AddEntry(histGraniitti, "Graniitti","ple");
    leg1->Draw("same");

    text = new TPaveText(0.51,0.18,0.93,0.40,"brNDC");
    SetTextStyle(text);
    text -> AddText("Statistical errors only");
    text -> AddText("Acceptance corrected");
    text -> AddText("Not efficiency corrected");
    text -> AddText("Not background subtracted");
    text -> Draw("same");

    newCanvas->Update();
    newCanvas->Write("kaons");
    //newCanvas->Close(); 
///////////////////////////// Inelastic ///////////////////////////////////////////
    newCanvas = new TCanvas("kaonsInel","kaonsInel",800,700);
    hist = (TH1D*)hInvMassCorr[Kaon][0][Inel]->Clone("hist"); 
    histCompare = (TH1D*)hInvMassCorr[Kaon][0][El]->Clone("histCompare");
    histGraniitti = (TH1D*)hInvMassGran[Inel][Kaon]->Clone("histGraniittiInel");
    histGraniittiEl = (TH1D*)hInvMassGran[El][Kaon]->Clone("histGraniittiEl");

    scaleFactor =   1 /hist->Integral();
    hist->Scale(scaleFactor);
    scaleFactor =   1 /histCompare->Integral();
    histCompare->Scale(scaleFactor);
    scaleFactor =   1 /histGraniitti->Integral();
    histGraniitti->Scale(scaleFactor);
    scaleFactor =   1 /histGraniittiEl->Integral();
    histGraniittiEl->Scale(scaleFactor);

    hist->SetTitle(" ; m(" + stateLabel + ") [GeV]; " + yLabel);
    SetGraphStyle(hist);
    hist->GetYaxis()->SetRangeUser(0, 0.615);   
    hist->Draw("hist E");

    gStyle->SetOptStat("");
    gStyle->SetPalette(1);


    SetMarkerStyle(histGraniitti);
    histGraniitti->Draw("SAME hist E");


    textSTAR = new TPaveText(0.17,0.89,0.33,0.95,"brNDC");
    SetTextStyle(textSTAR);
    textSTAR -> SetTextSize(textSize+0.02);
    textSTAR -> SetTextFont(72);
    textSTAR -> AddText("STAR");
    textSTAR -> Draw("same");
    textSTAR = new TPaveText(0.17,0.84,0.33,0.89,"brNDC");
    SetTextStyle(textSTAR);
    textSTAR -> SetTextFont(52);
    textSTAR -> AddText(label);
    textSTAR -> Draw("same");

    textPub = new TPaveText(0.35,0.88,0.88,0.95,"brNDC");
    SetTextStyle(textPub);
    textPub -> SetTextAlign(12);
    textPub -> AddText("p + p #rightarrow p + " + stateLabel + " + p       #sqrt{s} = 510 GeV");
    //textPub -> AddText("#sqrt{s} = 510 GeV");
    textPub -> Draw("same");

    text = new TPaveText(0.26,0.62,0.55,0.81,"brNDC");
    SetTextStyle(text);
    text -> SetTextAlign(31);
    text -> AddText("K^{+}, K^{-} kinematics:");
    text -> AddText("min(p_{T}^{+},p_{T}^{-}) < 0.7 GeV");   
    text -> AddText("p_{T} > 0.3 GeV");
    text -> AddText("|#eta| < 0.7");
    text -> Draw("same");

    text = new TPaveText(0.54,0.62,0.88,0.83,"brNDC");
    SetTextStyle(text);
    text -> SetTextAlign(12);
    text -> AddText("Forward proton kinematics:");
    text -> AddText("(p_{x} + 0.6)^{2} + p_{y}^{2} < 1.25 GeV^{2}");
    text -> AddText("0.4 GeV < |p_{y}| < 0.8 GeV");
    text -> AddText("p_{x} > -0.27 GeV");
    text -> Draw("same");


    leg1 = new TLegend(0.45, 0.48, 0.78, 0.58);
    SetLegendStyle(leg1);
    leg1->AddEntry(hist, "Data, #Delta#varphi < 90^{#circ} (unlike-sign pairs)","ple");
    leg1->AddEntry(histGraniitti, "Graniitti, #Delta#varphi < 90^{#circ}","ple");
    leg1->Draw("same");

    text = new TPaveText(0.51,0.18,0.93,0.40,"brNDC");
    SetTextStyle(text);
    text -> AddText("Statistical errors only");
    text -> AddText("Acceptance corrected");
    text -> AddText("Not efficiency corrected");
    text -> AddText("Not background subtracted");
    text -> Draw("same");

    newCanvas->Update();
    newCanvas->Write("kaonsInel");
    //newCanvas->Close();
//////////////////////////////////////////// Elastic ///////////////////////////////
    newCanvas = new TCanvas("kaonsEl","kaonsEl",800,700);
    histCompare->SetTitle(" ; m(" + stateLabel + ") [GeV]; " + yLabel);
    SetGraphStyle(histCompare);
    histCompare->GetYaxis()->SetRangeUser(0, 0.315);   
    histCompare->Draw("hist E");

    gStyle->SetOptStat("");
    gStyle->SetPalette(1);


    SetMarkerStyle(histGraniittiEl);
    histGraniittiEl->Draw("SAME hist E");

    textSTAR = new TPaveText(0.17,0.89,0.33,0.95,"brNDC");
    SetTextStyle(textSTAR);
    textSTAR -> SetTextSize(textSize+0.02);
    textSTAR -> SetTextFont(72);
    textSTAR -> AddText("STAR");
    textSTAR -> Draw("same");
    textSTAR = new TPaveText(0.17,0.84,0.33,0.89,"brNDC");
    SetTextStyle(textSTAR);
    textSTAR -> SetTextSize(textSize);
    textSTAR -> SetTextFont(52);
    textSTAR -> AddText(label);
    textSTAR -> Draw("same");

    textPub = new TPaveText(0.35,0.88,0.88,0.95,"brNDC");
    SetTextStyle(textPub);
    textPub -> SetTextAlign(12);
    textPub -> AddText("p + p #rightarrow p + " + stateLabel + " + p       #sqrt{s} = 510 GeV");
    textPub -> Draw("same");

    text = new TPaveText(0.16,0.64,0.45,0.83,"brNDC");
    SetTextStyle(text);
    text -> SetTextAlign(12);
    text -> AddText("K^{+}, K^{-} kinematics:");
    text -> AddText("min(p_{T}^{+},p_{T}^{-}) < 0.7 GeV");   
    text -> AddText("p_{T} > 0.3 GeV");
    text -> AddText("|#eta| < 0.7");
    text -> Draw("same");

    text = new TPaveText(0.54,0.62,0.88,0.83,"brNDC");
    SetTextStyle(text);
    text -> SetTextAlign(12);
    text -> AddText("Forward proton kinematics:");
    text -> AddText("(p_{x} + 0.6)^{2} + p_{y}^{2} < 1.25 GeV^{2}");
    text -> AddText("0.4 GeV < |p_{y}| < 0.8 GeV");
    text -> AddText("p_{x} > -0.27 GeV");
    text -> Draw("same");


    leg1 = new TLegend(0.45, 0.5, 0.78, 0.6);
    SetLegendStyle(leg1);
    leg1->AddEntry(histCompare, "Data, #Delta#varphi > 90^{#circ} (unlike-sign pairs)","ple");
    leg1->AddEntry(histGraniittiEl, "Graniitti, #Delta#varphi > 90^{#circ}","ple");
    leg1->Draw("same");

    text = new TPaveText(0.61,0.2,0.88,0.44,"brNDC");
    SetTextStyle(text);
    text -> AddText("Statistical errors only");
    text -> AddText("Acceptance corrected");
    text -> AddText("Not efficiency corrected");
    text -> AddText("Not background subtracted");
    text -> Draw("same");

    newCanvas->Update();
    newCanvas->Write("kaonsEl");
    //newCanvas->Close();
}

void PlotProtonsPlot()
{
    TCanvas* newCanvas = new TCanvas("protons","protons",800,700);
    gPad->SetMargin(0.13,0.02,0.105,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
    gPad->SetTickx();
    gPad->SetTicky(); 
    gPad->SetLogy(0);

    TString stateLabel = "p#bar{p}";

    TH1D *hist, *histCompare;
    TString yLabel = "Probability per event / 100 MeV";


    hist = (TH1D*)hInvMassCorr[Proton][0][ElInel]->Clone("hist"); 
    histCompare = (TH1D*)hInvMassCorr[Proton][1][ElInel]->Clone("histCompare");
    
    TH1D* histGraniittiEl, *histGraniittiInel, *histGraniitti;
    histGraniittiInel = (TH1D*)hInvMassGran[Inel][Proton]->Clone("histGraniittiInel");
    histGraniittiEl = (TH1D*)hInvMassGran[El][Proton]->Clone("histGraniittiEl");

    Double_t inel = hInvMassCorr[Proton][0][Inel]->Integral(); 
    Double_t el = hInvMassCorr[Proton][0][El]->Integral();
    //cout<<"Elastic / inelastic = "<< el<< " / "<< inel << "\n";
    histGraniittiInel->Scale(inel/((el+inel)*histGraniittiInel->Integral()));
    histGraniittiEl->Scale(el/((el+inel)*histGraniittiEl->Integral()));
    //cout<<"Scaled Elastic / inelastic = "<< histGraniittiEl->Integral() << " / "<< histGraniittiInel->Integral() << "\n";
    
    histGraniitti = (TH1D*)histGraniittiInel->Clone("histGraniitti");
    histGraniitti->Add(histGraniittiEl, 1.0);
    //cout<<"Sum: "<< histGraniitti->Integral() << "\n";

    Double_t scaleFactor; 
    scaleFactor =   1 /hist->Integral();
    hist->Scale(scaleFactor);
    histCompare->Scale(scaleFactor);

    hist->SetTitle(" ; m(" + stateLabel + ") [GeV]; " + yLabel);
    SetGraphStyle(hist);
    hist->GetXaxis()->SetRangeUser(1.8, 4.2); 
    hist->GetYaxis()->SetRangeUser(0, 0.4);   
    hist->Draw("hist E");

    gStyle->SetOptStat("");
    gStyle->SetPalette(1);

    SetMarkerStyle(histCompare);
    histCompare->SetMarkerColor(2);
    histCompare->SetMarkerStyle(22);
    histCompare->SetLineColor(2);
    histCompare->Draw("same hist E");

    SetMarkerStyle(histGraniitti);
    histGraniitti->Draw("same hist E");


    TPaveText *textSTAR;
    textSTAR = new TPaveText(0.17,0.89,0.33,0.95,"brNDC");
    SetTextStyle(textSTAR);
    textSTAR -> SetTextSize(textSize+0.02);
    textSTAR -> SetTextFont(72);
    textSTAR -> AddText("STAR");
    textSTAR -> Draw("same");
    textSTAR = new TPaveText(0.17,0.84,0.33,0.89,"brNDC");
    SetTextStyle(textSTAR);
    textSTAR -> SetTextSize(textSize);
    textSTAR -> SetTextFont(52);
    textSTAR -> AddText(label);
    textSTAR -> Draw("same");

    TPaveText *textPub = new TPaveText(0.35,0.88,0.88,0.95,"brNDC");
    SetTextStyle(textPub);
    textPub -> SetTextAlign(12);
    textPub -> AddText("p + p #rightarrow p + " + stateLabel + " + p       #sqrt{s} = 510 GeV");
    //textPub -> AddText("#sqrt{s} = 510 GeV");
    textPub -> Draw("same");


    TPaveText *text;
    text = new TPaveText(0.26,0.64,0.55,0.85,"brNDC");
    SetTextStyle(text);
    text -> SetTextAlign(32);
    text -> AddText("p, #bar{p} kinematics:");   
    text -> AddText("min(p_{T}^{+},p_{T}^{-}) < 1.1 GeV");   
    text -> AddText("p_{T} > 0.4 GeV");
    text -> AddText("|#eta| < 0.7");
    text -> Draw("same");

    text = new TPaveText(0.54,0.64,0.88,0.85,"brNDC");
    SetTextStyle(text);
    text -> SetTextAlign(12);
    text -> AddText("Forward proton kinematics:");
    text -> AddText("(p_{x} + 0.6)^{2} + p_{y}^{2} < 1.25 GeV^{2}");
    text -> AddText("0.4 GeV < |p_{y}| < 0.8 GeV");
    text -> AddText("p_{x} > -0.27 GeV");
    text -> Draw("same");

    text = new TPaveText(0.51,0.22,0.93,0.44,"brNDC");
    SetTextStyle(text);
    text -> AddText("Statistical errors only");
    text -> AddText("Acceptance corrected");
    text -> AddText("Not efficiency corrected");
    text -> AddText("Not background subtracted");
    text -> Draw("same");

    TLegend* leg1 = new TLegend(0.54, 0.48, 0.87, 0.63);
    SetLegendStyle(leg1);
    leg1->AddEntry(hist, "Data (unlike-sign pairs)","ple");
    leg1->AddEntry(histCompare, "Data (like-sign pairs)","ple");
    leg1->AddEntry(histGraniitti, "Graniitti","ple");
    leg1->Draw("same");

    newCanvas->Update();
    newCanvas->Write("protons");
    //newCanvas->Close();

///////////////////////////// Inelastic ///////////////////////////////////////////
    newCanvas = new TCanvas("protonsInel","protonsInel",800,700);   
 
    hist = (TH1D*)hInvMassCorr[Proton][0][Inel]->Clone("hist"); 
    histCompare = (TH1D*)hInvMassCorr[Proton][0][El]->Clone("histCompare");
    histGraniitti = (TH1D*)hInvMassGran[Inel][Proton]->Clone("histGraniittiInel");
    histGraniittiEl = (TH1D*)hInvMassGran[El][Proton]->Clone("histGraniittiEl");

    scaleFactor =   1 /hist->Integral();
    hist->Scale(scaleFactor);
    scaleFactor =   1 /histCompare->Integral();
    histCompare->Scale(scaleFactor);
    scaleFactor =   1 /histGraniitti->Integral();
    histGraniitti->Scale(scaleFactor);
    scaleFactor =   1 /histGraniittiEl->Integral();
    histGraniittiEl->Scale(scaleFactor);

    hist->SetTitle(" ; m(" + stateLabel + ") [GeV]; " + yLabel);
    SetGraphStyle(hist);
    hist->GetXaxis()->SetRangeUser(1.8, 4.2); 
    hist->GetYaxis()->SetRangeUser(0, 0.37);   
    hist->Draw("hist E");

    gStyle->SetOptStat("");
    gStyle->SetPalette(1);


    SetMarkerStyle(histGraniitti);
    histGraniitti->Draw("SAME hist E");

    textSTAR = new TPaveText(0.17,0.89,0.33,0.95,"brNDC");
    SetTextStyle(textSTAR);
    textSTAR -> SetTextSize(textSize+0.02);
    textSTAR -> SetTextFont(72);
    textSTAR -> AddText("STAR");
    textSTAR -> Draw("same");
    textSTAR = new TPaveText(0.17,0.84,0.33,0.89,"brNDC");
    SetTextStyle(textSTAR);
    textSTAR -> SetTextFont(52);
    textSTAR -> SetFillColor(0);
    textSTAR -> AddText(label);
    textSTAR -> Draw("same");

    textPub = new TPaveText(0.35,0.88,0.88,0.95,"brNDC");
    SetTextStyle(textPub);
    textPub -> SetTextAlign(12);
    textPub -> AddText("p + p #rightarrow p + " + stateLabel + " + p       #sqrt{s} = 510 GeV");
    //textPub -> AddText("#sqrt{s} = 510 GeV");
    textPub -> Draw("same");

    text = new TPaveText(0.26,0.64,0.55,0.85,"brNDC");
    SetTextStyle(text);
    text -> SetTextAlign(32);
    text -> AddText("p, #bar{p} kinematics:");   
    text -> AddText("min(p_{T}^{+},p_{T}^{-}) < 1.1 GeV");   
    text -> AddText("p_{T} > 0.4 GeV");
    text -> AddText("|#eta| < 0.7");
    text -> Draw("same");

    text = new TPaveText(0.54,0.64,0.88,0.85,"brNDC");
    SetTextStyle(text);
    text -> SetTextAlign(12);
    text -> AddText("Forward proton kinematics:");
    text -> AddText("(p_{x} + 0.6)^{2} + p_{y}^{2} < 1.25 GeV^{2}");
    text -> AddText("0.4 GeV < |p_{y}| < 0.8 GeV");
    text -> AddText("p_{x} > -0.27 GeV");
    text -> Draw("same");


    leg1 = new TLegend(0.45, 0.51, 0.78, 0.61);
    SetLegendStyle(leg1);
    leg1->AddEntry(hist, "Data, #Delta#varphi < 90^{#circ} (unlike-sign pairs)","ple");
    leg1->AddEntry(histGraniitti, "Graniitti, #Delta#varphi < 90^{#circ}","ple");
    leg1->Draw("same");

    text = new TPaveText(0.51,0.21,0.93,0.43,"brNDC");
    SetTextStyle(text);
    text -> AddText("Statistical errors only");
    text -> AddText("Acceptance corrected");
    text -> AddText("Not efficiency corrected");
    text -> AddText("Not background subtracted");
    text -> Draw("same");

    newCanvas->Update();
    newCanvas->Write("protonsInel");
    //newCanvas->Close();
//////////////////////////////////////////// Elastic ///////////////////////////////
    newCanvas = new TCanvas("protonsEl","protonsEl",800,700);
    histCompare->SetTitle(" ; m(" + stateLabel + ") [GeV]; " + yLabel);
    SetGraphStyle(histCompare);
    histCompare->GetXaxis()->SetRangeUser(1.8, 4.2); 
    histCompare->GetYaxis()->SetRangeUser(0, 0.44);   
    histCompare->Draw("hist E");
    gStyle->SetOptStat("");
    gStyle->SetPalette(1);

    SetMarkerStyle(histGraniittiEl);
    histGraniittiEl->Draw("SAME hist E");

    textSTAR = new TPaveText(0.17,0.89,0.33,0.95,"brNDC");
    SetTextStyle(textSTAR);
    textSTAR -> SetTextSize(textSize+0.02);
    textSTAR -> SetTextFont(72);
    textSTAR -> AddText("STAR");
    textSTAR -> Draw("same");

    textSTAR = new TPaveText(0.17,0.84,0.33,0.89,"brNDC");
    SetTextStyle(textSTAR);
    textSTAR -> SetTextFont(52);
    textSTAR -> AddText(label);
    textSTAR -> Draw("same");

    textPub = new TPaveText(0.35,0.88,0.88,0.95,"brNDC");
    SetTextStyle(textPub);
    textPub -> SetTextAlign(12);
    textPub -> AddText("p + p #rightarrow p + " + stateLabel + " + p       #sqrt{s} = 510 GeV");
    //textPub -> AddText("#sqrt{s} = 510 GeV");
    textPub -> Draw("same");

    text = new TPaveText(0.26,0.64,0.55,0.85,"brNDC");
    SetTextStyle(text);
    text -> SetTextAlign(32);
    text -> AddText("p, #bar{p} kinematics:");   
    text -> AddText("min(p_{T}^{+},p_{T}^{-}) < 1.1 GeV");   
    text -> AddText("p_{T} > 0.4 GeV");
    text -> AddText("|#eta| < 0.7");
    text -> Draw("same");

    text = new TPaveText(0.54,0.64,0.88,0.85,"brNDC");
    SetTextStyle(text);
    text -> SetTextAlign(12);
    text -> AddText("Forward proton kinematics:");
    text -> AddText("(p_{x} + 0.6)^{2} + p_{y}^{2} < 1.25 GeV^{2}");
    text -> AddText("0.4 GeV < |p_{y}| < 0.8 GeV");
    text -> AddText("p_{x} > -0.27 GeV");
    text -> Draw("same");;


    leg1 = new TLegend(0.45, 0.48, 0.78, 0.58);
    SetLegendStyle(leg1);
    leg1->AddEntry(hist, "Data, #Delta#varphi < 90^{#circ} (unlike-sign pairs)","ple");
    leg1->AddEntry(histGraniitti, "Graniitti, #Delta#varphi < 90^{#circ}","ple");
    leg1->Draw("same");

    text = new TPaveText(0.51,0.18,0.93,0.40,"brNDC");
    SetTextStyle(text);
    text -> AddText("Statistical errors only");
    text -> AddText("Acceptance corrected");
    text -> AddText("Not efficiency corrected");
    text -> AddText("Not background subtracted");
    text -> Draw("same");

    newCanvas->Update();
    newCanvas->Write("protonsEl");
    //newCanvas->Close();
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

void PlotCutsFlow()
{
    
    TFile* anaData = TFile::Open("/home/truhlar/Desktop/STAR/CEP/Analysis/Data/anaFlow2.root", "read");
    TH1I* hCutsData = (TH1I*)anaData -> Get("AnalFlow");
    TH1I* hCuts = new TH1I("AnalysisFlow", "CutsFlow", 16, 1, 17);
// //////////////////////////////////////////////////////////
// Plot Cuts Flow
    TString Labels[] = { TString("All"), TString("CEP trigger"), TString("2 RP tracks"), TString("Fiducial RP cut"), 
                      TString("2 TPC-TOF tracks"), TString("1 vertex"), TString("Tot. charge 0"), 
                      TString("|z_{vrtx}| < 80 cm"), TString("N_{hits}^{fit} #geq 25"), TString("N_{hits}^{dE/dx} #geq 15"),
                      TString("|DCA(z)| < 1 cm"), TString("DCA(xy) < 1.5 cm"), TString("|#eta| < 0.7"), 
                      TString("0.12< -t < 1.0 GeV^{2}"), TString("p_{T}^{miss} < 0.1 GeV"), TString("PID")};
                      //TString("-1.0 GeV^{2} < t < - 0.12 GeV^{2}")
    //gPad->SetMargin(0.9,0.02,0.1,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
    TCanvas* newCanvas = new TCanvas("newCanvas","newCanvas",950,600);

    for(int iBin = 1; iBin < hCuts->GetNbinsX(); ++iBin)
    {
        hCuts->SetBinContent(iBin,hCutsData->GetBinContent(iBin));
        hCuts->GetXaxis()->SetBinLabel(iBin, Labels[iBin-1]);
    }

    double PIDbinVal = 0;
    for (int i = 0; i < nParticles; ++i)
        PIDbinVal += hInvMassUncorr[i][0][0]->GetEntries();

    hCuts->SetBinContent(hCuts->GetNbinsX(),PIDbinVal);
    hCuts->GetXaxis()->SetBinLabel(hCuts->GetNbinsX(), Labels[hCuts->GetNbinsX()-1]);

    hCuts->SetTitle("; ; Number of events");
    SetGraphStyle(hCuts);
    hCuts->GetXaxis()->SetLabelFont(62);
    hCuts->GetXaxis()->SetLabelSize(labelSize-0.01);
    hCuts->SetMarkerSize(1.5);
    hCuts->SetLineWidth(2);    
    gPad->SetMargin(0.08,0.03,0.1,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
    gPad->SetTickx();
    gPad->SetTicky(); 
    gPad->SetLogy();
    hCuts->GetYaxis()->SetTitleOffset(0.85); 
    hCuts->GetYaxis()->SetRangeUser(4000, 3000000000); 
    hCuts->LabelsOption("d");
    hCuts->Draw("");
    gStyle->SetPaintTextFormat("1.2g");
    hCuts->Draw("TEXT15 same");

    TPaveText *textPub = new TPaveText(0.4,0.85,0.7,0.95,"brNDC");
    SetTextStyle(textPub);
    textPub -> AddText("p + p #rightarrow p + h^{+}h^{-} + p       #sqrt{s} = 510 GeV");
    textPub -> Draw("same");

    fout->cd();
    newCanvas->Update();
    newCanvas->Write("CutsFlow");
    //newCanvas->Close();
}

void RunGraniitti()
{


    TTree* tree[nParticles]; // 3 = 4PI state
    tree[Pion] = dynamic_cast<TTree*>( graniitti->Get("pionTree") );
    tree[Kaon] = dynamic_cast<TTree*>( graniitti->Get("kaonTree") );
    tree[Proton] = dynamic_cast<TTree*>( graniitti->Get("protonTree") );
//    tree[FourPions] = dynamic_cast<TTree*>( graniitti->Get("4pionTree") );

    if (!tree[Pion] || !tree[Kaon] || !tree[Proton]){
        cout<<"Error: cannot open one of the TTree in Graniitti"<<endl;
        return;
    }

    TString granForwardCuts =   TString("px_proton1 > - 0.3 && px_proton1 < 0.5 && abs(py_proton1) < 0.9 && abs(py_proton1) > 0.3 &&")
                                + TString("abs(py_proton1) < -px_proton1+1.05 && px_proton2 > - 0.3 && px_proton2 < 0.5 &&")
                                + TString(" abs(py_proton2) < 0.9 && abs(py_proton2) > 0.3 && abs(py_proton2) < -px_proton2+1.05");
    TString graniittiCuts = "eta_part1 > -0.7 && eta_part1 < 0.7 && eta_part2 > -0.7 && eta_part2 < 0.7 && t_proton1 < -0.12 && t_proton2 < -0.12 && t_proton1 > -1.0  && t_proton2 > -1.0 && " + granForwardCuts;
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
        for (int iPart = 0; iPart < 3; ++iPart) // nParticles; ++iPart)
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


void SetGraphStyle(TH1* hist)
{
    gPad->SetMargin(0.13,0.03,0.105,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
    gPad->SetTickx();
    gPad->SetTicky(); 
    gPad->SetLogy(0);
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
    hist->SetMarkerColor(1);
    hist->SetMarkerSize(1.5);
    hist->SetMarkerStyle(29);
    hist->SetLineColor(1);
    hist->SetLineStyle(1);
    hist->SetLineWidth(1);

}//SetGraphStyle

void SetTextStyle(TPaveText* text)
{
    text -> SetTextSize(textSize);
    text -> SetTextAlign(22);
    text -> SetFillColor(0);
    text -> SetFillStyle(0);
    text -> SetTextFont(fontStyle);
    text -> SetBorderSize(0);
}//SetTextStyle

void SetLegendStyle(TLegend* leg1)
{
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(textSize);
    leg1->SetTextFont(fontStyle);
    leg1->SetMargin(0.1);
}//SetLegendStyle

void SetMarkerStyle(TH1* hist)
{
    hist->SetMarkerColor(4);
    hist->SetMarkerSize(1);
    hist->SetMarkerStyle(20);
    hist->SetLineColor(4);
    hist->SetLineStyle(1);
    hist->SetLineWidth(1);

}//SetMarkerStyle

 
void SetLineStyle(TLine* line)
{
    line->SetLineStyle(1);
    line->SetLineColor(1);
    line->SetLineWidth(4);
}//SetLineStyle
