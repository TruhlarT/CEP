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

TH1F* hInvMass[4];
TH1D* hDeltaPhi[3];
TH2D* hPxPy[3];

 // 0 = pi- 1 = K- 2 = pbar
TH3F* hTPCeff[6]; // 3 = pi+ 4 = K+ 5 = p
TH1D* hInvMassCorr[nParticles][2][3]; // 0 - signal, 1 - Background
                                    //  - El + Inel, 1 - El, 2 - Inel 
TH1D* hInvMassUncorr[nParticles][2][3]; // 0 - signal, 1 - Background

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


void PlotDeltaPhi();
void PlotRap();
void PlotKaonsPlot();
void PlotProtonsPlot();


void thesisPlots()
{

    TString output = "/home/truhlar/Desktop/STAR/CEP/Analysis/Outputs/thesisPlot.root";
    TString input = "/home/truhlar/Desktop/STAR/CEP/Analysis/Data/anaFlow2.root";


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


    PlotDeltaPhi();
    PlotRap();
    PlotKaonsPlot();
    PlotProtonsPlot();

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
 
        double deltaPhi = TMath::Abs(phiRp[East] - phiRp[West])*convertToDegree;
        if(deltaPhi > 180)
            deltaPhi = 360 - deltaPhi;

        hDeltaPhi[0]->Fill(deltaPhi);
        if(elastic)
            hDeltaPhi[1]->Fill(deltaPhi);
        else
            hDeltaPhi[2]->Fill(deltaPhi);

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

                hRapidityCorr[Proton][signal][0]->Fill(pairRapidity, 1/effTotal);
                hRapidityCorr[Proton][signal][combination]->Fill(pairRapidity, 1/effTotal);

                hDeltaPhiCorr[Proton][signal][0]->Fill(deltaPhi, 1/effTotal);
                hDeltaPhiCorr[Proton][signal][combination]->Fill(deltaPhi, 1/effTotal);
            }
            if(transMomentum[0] > 0.4 && transMomentum[1] > 0.4 && (transMomentum[0] < 1.1 || transMomentum[1] < 1.1) )
            {
                hInvMassUncorr[Proton][signal][0]->Fill(invMass[Proton]);
                hInvMassUncorr[Proton][signal][combination]->Fill(invMass[Proton]);

                hRapidityUncorr[Proton][signal][0]->Fill(pairRapidity);
                hRapidityUncorr[Proton][signal][combination]->Fill(pairRapidity);

                hDeltaPhiUncorr[Proton][signal][0]->Fill(deltaPhi);
                hDeltaPhiUncorr[Proton][signal][combination]->Fill(deltaPhi);
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

                hRapidityCorr[Kaon][signal][0]->Fill(pairRapidity, 1/effTotal);
                hRapidityCorr[Kaon][signal][combination]->Fill(pairRapidity, 1/effTotal);

                hDeltaPhiCorr[Kaon][signal][0]->Fill(deltaPhi, 1/effTotal);
                hDeltaPhiCorr[Kaon][signal][combination]->Fill(deltaPhi, 1/effTotal);
            }
            if(transMomentum[0] > 0.3 && transMomentum[1] > 0.3 && (transMomentum[0] < 0.7 || transMomentum[1] < 0.7) )
            {
                hInvMassUncorr[Kaon][signal][0]->Fill(invMass[Kaon]);
                hInvMassUncorr[Kaon][signal][combination]->Fill(invMass[Kaon]);

                hRapidityUncorr[Kaon][signal][0]->Fill(pairRapidity);
                hRapidityUncorr[Pion][signal][combination]->Fill(pairRapidity);

                hDeltaPhiUncorr[Kaon][signal][0]->Fill(deltaPhi);
                hDeltaPhiUncorr[Kaon][signal][combination]->Fill(deltaPhi);
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

                hRapidityCorr[Pion][signal][0]->Fill(pairRapidity, 1/effTotal);
                hRapidityCorr[Pion][signal][combination]->Fill(pairRapidity, 1/effTotal);

                hDeltaPhiCorr[Pion][signal][0]->Fill(deltaPhi, 1/effTotal);
                hDeltaPhiCorr[Pion][signal][combination]->Fill(deltaPhi, 1/effTotal);
            }
            if(transMomentum[0] > 0.2 && transMomentum[1] > 0.2)
            {
                hInvMassUncorr[Pion][signal][0]->Fill(invMass[Pion]);
                hInvMassUncorr[Pion][signal][combination]->Fill(invMass[Pion]);

                hRapidityUncorr[Pion][signal][0]->Fill(pairRapidity);
                hRapidityUncorr[Pion][signal][combination]->Fill(pairRapidity);

                hDeltaPhiUncorr[Pion][signal][0]->Fill(deltaPhi);
                hDeltaPhiUncorr[Pion][signal][combination]->Fill(deltaPhi);
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
    
        hRapidityCorr[0][0][i]  = new TH1D("corrRapidity" + particleLables[0] + combinationLabel[i] + "Sig", "Corrected rapidity " + particleLables[0] + combinationLabel[i], 16, -0.8, 0.8);
        hRapidityCorr[1][0][i]  = new TH1D("corrRapidity" + particleLables[1] + combinationLabel[i] + "Sig", "Corrected rapidity " + particleLables[1] + combinationLabel[i], 10, -0.8, 0.8);
        hRapidityCorr[2][0][i]  = new TH1D("corrRapidity" + particleLables[2] + combinationLabel[i] + "Sig", "Corrected rapidity " + particleLables[2] + combinationLabel[i], 8, -0.8, 0.8);
        hRapidityCorr[0][1][i]  = new TH1D("corrRapidity" + particleLables[0] + combinationLabel[i] + "Bcg", "Corrected rapidity " + particleLables[0] + combinationLabel[i], 16, -0.8, 0.8);
        hRapidityCorr[1][1][i]  = new TH1D("corrRapidity" + particleLables[1] + combinationLabel[i] + "Bcg", "Corrected rapidity " + particleLables[1] + combinationLabel[i], 10, -0.8, 0.8);
        hRapidityCorr[2][1][i]  = new TH1D("corrRapidity" + particleLables[2] + combinationLabel[i] + "Bcg", "Corrected rapidity " + particleLables[2] + combinationLabel[i], 8, -0.8, 0.8);

        hRapidityUncorr[0][0][i]  = new TH1D("uncorrRapidity" + particleLables[0] + combinationLabel[i] + "Sig", "Uncorrected rapidity " + particleLables[0] + combinationLabel[i], 16, -0.8, 0.8);
        hRapidityUncorr[1][0][i]  = new TH1D("uncorrRapidity" + particleLables[1] + combinationLabel[i] + "Sig", "Uncorrected rapidity " + particleLables[1] + combinationLabel[i], 10, -0.8, 0.8);
        hRapidityUncorr[2][0][i]  = new TH1D("uncorrRapidity" + particleLables[2] + combinationLabel[i] + "Sig", "Uncorrected rapidity " + particleLables[2] + combinationLabel[i], 8, -0.8, 0.8);
        hRapidityUncorr[0][1][i]  = new TH1D("uncorrRapidity" + particleLables[0] + combinationLabel[i] + "Bcg", "Uncorrected rapidity " + particleLables[0] + combinationLabel[i], 16, -0.8, 0.8);
        hRapidityUncorr[1][1][i]  = new TH1D("uncorrRapidity" + particleLables[1] + combinationLabel[i] + "Bcg", "Uncorrected rapidity " + particleLables[1] + combinationLabel[i], 10, -0.8, 0.8);
        hRapidityUncorr[2][1][i]  = new TH1D("uncorrRapidity" + particleLables[2] + combinationLabel[i] + "Bcg", "Uncorrected rapidity " + particleLables[2] + combinationLabel[i], 8, -0.8, 0.8);

        hDeltaPhiCorr[0][0][i]  = new TH1D("corrDeltaPhi" + particleLables[0] + combinationLabel[i] + "Sig", "Corrected delta phi " + particleLables[0] + combinationLabel[i], 18, 0.0, 180);
        hDeltaPhiCorr[1][0][i]  = new TH1D("corrDeltaPhi" + particleLables[1] + combinationLabel[i] + "Sig", "Corrected delta phi " + particleLables[1] + combinationLabel[i], 18, 0.0, 180);
        hDeltaPhiCorr[2][0][i]  = new TH1D("corrDeltaPhi" + particleLables[2] + combinationLabel[i] + "Sig", "Corrected delta phi " + particleLables[2] + combinationLabel[i], 9, 0.0, 180);
        hDeltaPhiCorr[0][1][i]  = new TH1D("corrDeltaPhi" + particleLables[0] + combinationLabel[i] + "Bcg", "Corrected delta phi " + particleLables[0] + combinationLabel[i], 18, 0.0, 180);
        hDeltaPhiCorr[1][1][i]  = new TH1D("corrDeltaPhi" + particleLables[1] + combinationLabel[i] + "Bcg", "Corrected delta phi " + particleLables[1] + combinationLabel[i], 18, 0.0, 180);
        hDeltaPhiCorr[2][1][i]  = new TH1D("corrDeltaPhi" + particleLables[2] + combinationLabel[i] + "Bcg", "Corrected delta phi " + particleLables[2] + combinationLabel[i], 9, 0.0, 180);

        hDeltaPhiUncorr[0][0][i]  = new TH1D("uncorrDeltaPhi" + particleLables[0] + combinationLabel[i] + "Sig", "Uncorrected delta phi " + particleLables[0] + combinationLabel[i], 18, 0.0, 180);
        hDeltaPhiUncorr[1][0][i]  = new TH1D("uncorrDeltaPhi" + particleLables[1] + combinationLabel[i] + "Sig", "Uncorrected delta phi " + particleLables[1] + combinationLabel[i], 18, 0.0, 180);
        hDeltaPhiUncorr[2][0][i]  = new TH1D("uncorrDeltaPhi" + particleLables[2] + combinationLabel[i] + "Sig", "Uncorrected delta phi " + particleLables[2] + combinationLabel[i], 9, 0.0, 180);
        hDeltaPhiUncorr[0][1][i]  = new TH1D("uncorrDeltaPhi" + particleLables[0] + combinationLabel[i] + "Bcg", "Uncorrected delta phi " + particleLables[0] + combinationLabel[i], 18, 0.0, 180);
        hDeltaPhiUncorr[1][1][i]  = new TH1D("uncorrDeltaPhi" + particleLables[1] + combinationLabel[i] + "Bcg", "Uncorrected delta phi " + particleLables[1] + combinationLabel[i], 18, 0.0, 180);
        hDeltaPhiUncorr[2][1][i]  = new TH1D("uncorrDeltaPhi" + particleLables[2] + combinationLabel[i] + "Bcg", "Uncorrected delta phi " + particleLables[2] + combinationLabel[i], 9, 0.0, 180);
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



void PlotDeltaPhi()
{
    TCanvas* newCanvas = new TCanvas("newCanvas","newCanvas",800,700);
    gPad->SetMargin(0.13,0.04,0.105,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
    gPad->SetTickx();
    gPad->SetTicky(); 
    gPad->SetLogy(0);
    gStyle->SetOptStat("");
    gStyle->SetPalette(1);
    gStyle->SetLineWidth(2);      //axis line
    gStyle->SetFrameLineWidth(2); //frame line


    TString stateLabel[] = {"#pi^{+}#pi^{-}","K^{+}K^{-}","p#bar{p}"};

    TH1D *hist, *histCompare;
    TString yLabel = "Probability per event";
    for (int comb = 0; comb < nCombination; ++comb)
    {
        for (int iPart = 0; iPart < nParticles; ++iPart)
        {
            

            hist = (TH1D*)hDeltaPhiCorr[iPart][0][comb]->Clone("hist"); 
            histCompare = (TH1D*)hDeltaPhiCorr[iPart][1][comb]->Clone("histCompare");

            Double_t scaleFactor; 
            scaleFactor =   1 /hist->Integral();
            hist->Scale(scaleFactor);
            histCompare->Scale(scaleFactor);

            hist->SetTitle(" ; #Delta#varphi [deg]; " + yLabel);
            hist->SetStats(false);
            hist->GetXaxis()->SetTitleFont(42);
            hist->GetYaxis()->SetTitleFont(42);
            hist->GetXaxis()->SetLabelFont(42);
            hist->GetYaxis()->SetLabelFont(42);
            hist->GetXaxis()->SetLabelSize(labelSize);
            hist->GetYaxis()->SetLabelSize(labelSize);
            hist->GetXaxis()->SetTitleSize(labelSize);
            hist->GetYaxis()->SetTitleSize(labelSize);
            hist->GetXaxis()->SetTitleOffset(1.0);
            hist->GetYaxis()->SetTitleOffset(1.4);
            //hist->GetYaxis()->SetRangeUser(0, 0.235);   
            hist->SetMarkerColor(1);
            hist->SetMarkerSize(1);
            hist->SetMarkerStyle(20);
            hist->SetLineColor(1);
            hist->SetLineStyle(1);
            hist->SetLineWidth(1);
            hist->Draw("E");

            gStyle->SetOptStat("");
            gStyle->SetPalette(1);

            histCompare->SetMarkerColor(2);
            histCompare->SetMarkerSize(1);
            histCompare->SetMarkerStyle(22);
            histCompare->SetLineColor(2);
            histCompare->SetLineStyle(1);
            histCompare->SetLineWidth(1);
            histCompare->Draw("ESAME");

            TPaveText *textSTAR;
            textSTAR = new TPaveText(0.3,0.9,0.7,0.95,"brNDC");
            textSTAR -> SetTextSize(textSize);
            textSTAR -> SetFillColor(0);
            textSTAR -> SetTextAlign(22);
            textSTAR -> SetTextFont(62);
            textSTAR->AddText("THIS THESIS");
            textSTAR -> Draw("same");

            TPaveText *textPub = new TPaveText(0.3,0.79,0.7,0.89,"brNDC");
            textPub -> SetTextSize(textSize);
            textPub -> SetTextAlign(22);
            textPub -> SetFillColor(0);
            textPub -> SetTextFont(42);
            textPub -> SetTextAlign(22);
            textPub -> AddText("p + p #rightarrow p + " + stateLabel[iPart] + " + p");
            textPub -> AddText("#sqrt{s} = 510 GeV");
            textPub -> Draw("same");

            if(comb != 0)
            {
                TString phiLabel = "#Delta#varphi > 90^{#circ}";
                if(comb == Inel)
                    phiLabel = "#Delta#varphi < 90^{#circ}";
                textPub = new TPaveText(0.3,0.6,0.7,0.66,"brNDC");
                textPub -> SetTextSize(textSize);
                textPub -> SetTextAlign(22);
                textPub -> SetFillColor(0);
                textPub -> SetTextFont(42);
                textPub -> AddText(phiLabel);
                textPub -> Draw("same");
            }

            TLegend* leg1 = new TLegend(0.3, 0.68, 0.7, 0.78);
            leg1->SetFillStyle(0);
            leg1->SetBorderSize(0);
            leg1->SetTextSize(textSize);
            leg1->SetTextFont(42);
            leg1 -> SetTextAlign(22);
            leg1->SetMargin(0.1);
            leg1->AddEntry(hist, "Data (unlike-sign pairs)","pe");
            leg1->AddEntry(histCompare, "Data (like-sign pairs)","pe");
            leg1->Draw("same");

            newCanvas->Update();
            newCanvas->Write("deltaPhi" + particleLables[iPart] + combinationLabel[comb]);
        }
    }
    newCanvas->Close();

}

void PlotRap()
{
    TCanvas* newCanvas = new TCanvas("newCanvas","newCanvas",800,700);
    gPad->SetMargin(0.13,0.04,0.105,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
    gPad->SetTickx();
    gPad->SetTicky(); 
    gPad->SetLogy(0);
    gStyle->SetOptStat("");
    gStyle->SetPalette(1);
    gStyle->SetLineWidth(2);      //axis line
    gStyle->SetFrameLineWidth(2); //frame line


    TString stateLabel[] = {"#pi^{+}#pi^{-}","K^{+}K^{-}","p#bar{p}"};

    TH1D *hist, *histCompare;
    TString yLabel = "Probability per event";
    for (int comb = 0; comb < nCombination; ++comb)
    {
        for (int iPart = 0; iPart < nParticles; ++iPart)
        {
            

            hist = (TH1D*)hRapidityCorr[iPart][0][comb]->Clone("hist"); 
            histCompare = (TH1D*)hRapidityCorr[iPart][1][comb]->Clone("histCompare");

            Double_t scaleFactor; 
            scaleFactor =   1 /hist->Integral();
            hist->Scale(scaleFactor);
            histCompare->Scale(scaleFactor);

            hist->SetTitle(" ; y(" + stateLabel[iPart] + "); " + yLabel);
            hist->SetStats(false);
            hist->GetXaxis()->SetTitleFont(42);
            hist->GetYaxis()->SetTitleFont(42);
            hist->GetXaxis()->SetLabelFont(42);
            hist->GetYaxis()->SetLabelFont(42);
            hist->GetXaxis()->SetLabelSize(labelSize);
            hist->GetYaxis()->SetLabelSize(labelSize);
            hist->GetXaxis()->SetTitleSize(labelSize);
            hist->GetYaxis()->SetTitleSize(labelSize);
            hist->GetXaxis()->SetTitleOffset(1.0);
            hist->GetYaxis()->SetTitleOffset(1.4);
            hist->GetYaxis()->SetRangeUser(0, hist->GetMaximum()*1.5);
            //if(iPart == Proton)
              // hist->GetYaxis()->SetRangeUser(0, hist->GetMaximum()*1.8);
            hist->SetMarkerColor(1);
            hist->SetMarkerSize(1);
            hist->SetMarkerStyle(20);
            hist->SetLineColor(1);
            hist->SetLineStyle(1);
            hist->SetLineWidth(1);
            hist->Draw("E");

            gStyle->SetOptStat("");
            gStyle->SetPalette(1);

            histCompare->SetMarkerColor(2);
            histCompare->SetMarkerSize(1);
            histCompare->SetMarkerStyle(22);
            histCompare->SetLineColor(2);
            histCompare->SetLineStyle(1);
            histCompare->SetLineWidth(1);
            histCompare->Draw("ESAME");

            TPaveText *textSTAR;
            textSTAR = new TPaveText(0.15,0.9,0.3,0.95,"brNDC");
            textSTAR -> SetTextSize(textSize);
            textSTAR -> SetFillColor(0);
            textSTAR -> SetTextAlign(12);
            textSTAR -> SetTextFont(62);
            textSTAR->AddText("THIS THESIS");
            textSTAR -> Draw("same");

            TPaveText *textPub = new TPaveText(0.4,0.9,0.95,0.95,"brNDC");
            textPub -> SetTextSize(textSize);
            textPub -> SetTextAlign(32);
            textPub -> SetFillColor(0);
            textPub -> SetTextFont(42);
            textPub -> AddText("p + p #rightarrow p + " + stateLabel[iPart] + " + p      #sqrt{s} = 510 GeV");
            //textPub -> AddText("#sqrt{s} = 510 GeV");
            textPub -> Draw("same");

            if(comb != 0)
            {
                TString phiLabel = "#Delta#varphi > 90^{#circ}";
                if(comb == Inel)
                    phiLabel = "#Delta#varphi < 90^{#circ}";
                textPub = new TPaveText(0.15,0.82,0.95,0.88,"brNDC");
                textPub -> SetTextSize(textSize);
                textPub -> SetTextAlign(32);
                textPub -> SetFillColor(0);
                textPub -> SetTextFont(42);
                textPub -> AddText(phiLabel);
                textPub -> Draw("same");
            }

            TLegend* leg1 = new TLegend(0.15, 0.78, 0.7, 0.88);
            leg1->SetFillStyle(0);
            leg1->SetBorderSize(0);
            leg1->SetTextSize(textSize);
            leg1->SetTextFont(42);
            leg1 -> SetTextAlign(12);
            leg1->SetMargin(0.1);
            leg1->AddEntry(hist, "Data (unlike-sign pairs)","pe");
            leg1->AddEntry(histCompare, "Data (like-sign pairs)","pe");
            leg1->Draw("same");

            newCanvas->Update();
            newCanvas->Write("pairRap" + particleLables[iPart] + combinationLabel[comb]);
        }
    }
    newCanvas->Close();
}

void PlotKaonsPlot()
{

    TCanvas* newCanvas = new TCanvas("newCanvas","newCanvas",800,700);
    gPad->SetMargin(0.13,0.025,0.105,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
    gPad->SetTickx();
    gPad->SetTicky(); 
    gPad->SetLogy(0);
    gStyle->SetOptStat("");
    gStyle->SetPalette(1);
    gStyle->SetLineWidth(2);      //axis line
    gStyle->SetFrameLineWidth(2); //frame line

    TString stateLabel = "K^{+}K^{-}";

    TH1D *hist, *histCompare;
    TString yLabel = "Probability per event / 50 MeV";


    hist = (TH1D*)hInvMassCorr[Kaon][0][Inel]->Clone("hist"); 
    histCompare = (TH1D*)hInvMassCorr[Kaon][0][El]->Clone("histCompare");
    

    Double_t scaleFactor =   1 /hist->Integral();
    hist->Scale(scaleFactor);
    histCompare->Scale(scaleFactor);

    hist->SetTitle(" ; m(" + stateLabel + ") [GeV]; " + yLabel);
    hist->SetStats(false);
    hist->GetXaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetTitleFont(42);
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelFont(42);
    hist->GetXaxis()->SetLabelSize(labelSize);
    hist->GetYaxis()->SetLabelSize(labelSize);
    hist->GetXaxis()->SetTitleSize(labelSize);
    hist->GetYaxis()->SetTitleSize(labelSize);
    hist->GetXaxis()->SetTitleOffset(0.9);
    hist->GetYaxis()->SetTitleOffset(1.3);
    hist->GetYaxis()->SetRangeUser(0, 0.37);   
    hist->SetMarkerColor(4);
    hist->SetMarkerSize(1);
    hist->SetMarkerStyle(20);
    hist->SetLineColor(4);
    hist->SetLineStyle(1);
    hist->SetLineWidth(1);
    hist->Draw("E");

    gStyle->SetOptStat("");
    gStyle->SetPalette(1);

    histCompare->SetMarkerColor(1);
    histCompare->SetMarkerSize(1);
    histCompare->SetMarkerStyle(21);
    histCompare->SetLineColor(1);
    histCompare->SetLineStyle(1);
    histCompare->SetLineWidth(1);
    histCompare->Draw("ESAME");


    TPaveText *textSTAR;
    textSTAR = new TPaveText(0.15,0.9,0.3,0.95,"brNDC");
    textSTAR -> SetTextSize(textSize);
    textSTAR -> SetFillColor(0);
    textSTAR -> SetTextAlign(12);
    textSTAR -> SetTextFont(62);
    textSTAR->AddText("THIS THESIS");
    textSTAR -> Draw("same");

    TPaveText *textPub = new TPaveText(0.4,0.9,0.95,0.95,"brNDC");
    textPub -> SetTextSize(textSize);
    textPub -> SetTextAlign(32);
    textPub -> SetFillColor(0);
    textPub -> SetTextFont(42);
    textPub -> AddText("p + p #rightarrow p + " + stateLabel + " + p      #sqrt{s} = 510 GeV");
    //textPub -> AddText("#sqrt{s} = 510 GeV");
    textPub -> Draw("same");


    TLegend* leg1 = new TLegend(0.3, 0.78, 0.95, 0.88);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(textSize);
    leg1->SetTextFont(42);
    leg1 -> SetTextAlign(32);
    leg1->SetMargin(0.44);
    leg1->AddEntry(hist, "Data, #Delta#varphi < 90^{#circ} (unlike-sign pairs)","pe");
    leg1->AddEntry(histCompare, "Data, #Delta#varphi > 90^{#circ} (unlike-sign pairs)","pe");
    leg1->Draw("same");

    newCanvas->Update();
    newCanvas->Write("kaonsEl+Inel");
    newCanvas->Close();
}

void PlotProtonsPlot()
{

    TCanvas* newCanvas = new TCanvas("newCanvas","newCanvas",800,700);
    gPad->SetMargin(0.13,0.025,0.105,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
    gPad->SetTickx();
    gPad->SetTicky(); 
    gPad->SetLogy(0);
    gStyle->SetOptStat("");
    gStyle->SetPalette(1);
    gStyle->SetLineWidth(2);      //axis line
    gStyle->SetFrameLineWidth(2); //frame line

    TString stateLabel = "p#bar{p}";

    TH1D *hist, *histCompare;
    TString yLabel = "Probability per event / 100 MeV";


    hist = (TH1D*)hInvMassCorr[Proton][0][Inel]->Clone("hist"); 
    histCompare = (TH1D*)hInvMassCorr[Proton][0][El]->Clone("histCompare");
    

    Double_t scaleFactor =   1 /hist->Integral();
    hist->Scale(scaleFactor);
    histCompare->Scale(scaleFactor);

    hist->SetTitle(" ; m(" + stateLabel + ") [GeV]; " + yLabel);
    hist->SetStats(false);
    hist->GetXaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetTitleFont(42);
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelFont(42);
    hist->GetXaxis()->SetLabelSize(labelSize);
    hist->GetYaxis()->SetLabelSize(labelSize);
    hist->GetXaxis()->SetTitleSize(labelSize);
    hist->GetYaxis()->SetTitleSize(labelSize);
    hist->GetXaxis()->SetTitleOffset(0.9);
    hist->GetYaxis()->SetTitleOffset(1.3);
    hist->GetYaxis()->SetRangeUser(0, histCompare->GetMaximum()*1.45);   
    hist->SetMarkerColor(4);
    hist->SetMarkerSize(1);
    hist->SetMarkerStyle(20);
    hist->SetLineColor(4);
    hist->SetLineStyle(1);
    hist->SetLineWidth(1);
    hist->Draw("E");

    gStyle->SetOptStat("");
    gStyle->SetPalette(1);

    histCompare->SetMarkerColor(1);
    histCompare->SetMarkerSize(1);
    histCompare->SetMarkerStyle(21);
    histCompare->SetLineColor(1);
    histCompare->SetLineStyle(1);
    histCompare->SetLineWidth(1);
    histCompare->Draw("ESAME");


    TPaveText *textSTAR;
    textSTAR = new TPaveText(0.15,0.9,0.3,0.95,"brNDC");
    textSTAR -> SetTextSize(textSize);
    textSTAR -> SetFillColor(0);
    textSTAR -> SetTextAlign(12);
    textSTAR -> SetTextFont(62);
    textSTAR->AddText("THIS THESIS");
    textSTAR -> Draw("same");

    TPaveText *textPub = new TPaveText(0.4,0.9,0.95,0.95,"brNDC");
    textPub -> SetTextSize(textSize);
    textPub -> SetTextAlign(32);
    textPub -> SetFillColor(0);
    textPub -> SetTextFont(42);
    textPub -> AddText("p + p #rightarrow p + " + stateLabel + " + p      #sqrt{s} = 510 GeV");
    //textPub -> AddText("#sqrt{s} = 510 GeV");
    textPub -> Draw("same");


    TLegend* leg1 = new TLegend(0.3, 0.78, 0.94, 0.88);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(textSize);
    leg1->SetTextFont(42);
    leg1 -> SetTextAlign(32);
    leg1->SetMargin(0.44);
    leg1->AddEntry(hist, "Data, #Delta#varphi < 90^{#circ} (unlike-sign pairs)","pe");
    leg1->AddEntry(histCompare, "Data, #Delta#varphi > 90^{#circ} (unlike-sign pairs)","pe");
    leg1->Draw("same");

    newCanvas->Update();
    newCanvas->Write("protonsEl+Inel");
    newCanvas->Close();
}
