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

TString combinationLabel[nCombination] = { TString("El+Inel"), TString("El"), TString("Inel")};
TString sideLabel[nSides] = { TString("East"), TString("West")};

TFile* fout;
TFile* TPCeff[6];
TFile* data;
TFile* graniitti;

TTree* tree;
TTree* treeBack;
TTree* granitiTree[4];

TH1D* hInvMass[4];
TH1D* hDeltaPhi[3];
TH2D* hPxPy[3];

TH1D* gDeltaPhi[4];
TH1D* gRap[4];
TH1D* gInvMass[4][3];

 // 0 = pi- 1 = K- 2 = pbar
TH3F* hTPCeff[6]; // 3 = pi+ 4 = K+ 5 = p
TH1D* hInvMassCorr[nParticles+1][2][3]; // 0 - signal, 1 - Background
                                    //  - El + Inel, 1 - El, 2 - Inel 
 // 0 - signal, 1 - Background

TH1D* hRapidityCorr[nParticles+1][2][3];  

TH1D* hDeltaPhiCorr[nParticles+1][2][3];  

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


TString particleLables[nParticles+1] = { TString("Pion"), TString("Kaon"), TString("Proton"), TString("4Pion")};
TString stateLabel[] = {"#pi^{+}#pi^{-}","K^{+}K^{-}","p#bar{p}","#pi^{+}#pi^{+}#pi^{-}#pi^{-}"};

const int xMax = 20;
const int yMax = 20;
const int zMax = 24;

const double textSize = 0.04;
const double labelSize = 0.045;

void Init();
void ConnectInput(TTree* tree);
void Make(int signal);

void RunGraniti();

void PlotDeltaPhi();
void PlotRap();
void PlotMassPlots();


void granitiPlots()
{

    TString output = "/home/truhlar/Desktop/STAR/CEP/Analysis/Outputs/granitiPlot.root";
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

    RunGraniti();

    PlotDeltaPhi();
    PlotRap();
    PlotMassPlots();

    fout->Write();
    fout->Close();
    
}


void RunGraniti()
{

    TString granInput = "/home/truhlar/Downloads/graniitti/output/Graniitti.root";
    graniitti = TFile::Open(granInput, "read");
    if (!graniitti){
        cout<<"Error: cannot open: "<<granInput<<endl;
        return 6;
    }

    granitiTree[Pion] = dynamic_cast<TTree*>( graniitti->Get("pionTree") );
    granitiTree[Kaon] = dynamic_cast<TTree*>( graniitti->Get("kaonTree") );
    granitiTree[Proton] = dynamic_cast<TTree*>( graniitti->Get("protonTree") );
    granitiTree[3] = dynamic_cast<TTree*>( graniitti->Get("4pionTree") );

    if (!granitiTree[Pion] || !granitiTree[Kaon] || !granitiTree[Proton] || !granitiTree[3]){
        cout<<"Error: cannot open one of the TTree in Graniitti"<<endl;
        return;
    }
    fout->cd();
    TString granForwardCuts =   TString("px_proton1 > - 0.27 && abs(py_proton1) < 0.8 && abs(py_proton1) > 0.4 && ")
                            + TString("px_proton2 > - 0.27 && abs(py_proton2) < 0.8 && abs(py_proton2) > 0.4 && ")
                            + TString("(px_proton1 + 0.6)*(px_proton1 + 0.6) + py_proton1*py_proton1 < 1.25 && ")
                            + TString("(px_proton2 + 0.6)*(px_proton2 + 0.6) + py_proton2*py_proton2 < 1.25");

    TString graniittiCuts = "eta_part1 > -0.7 && eta_part1 < 0.7 && eta_part2 > -0.7 && eta_part2 < 0.7 && t_proton1 < -0.12 && t_proton2 < -0.12 && t_proton1 > -1.0  && t_proton2 > -1.0 && " + granForwardCuts;
    TString partCuts[] = { TString("pT_part1 > 0.2 && pT_part2 > 0.2 "), 
                        TString("pT_part1 > 0.3 && pT_part2 > 0.3 && TMath::Min(pT_part1,pT_part2) < 0.7 "),
                        TString("pT_part1 > 0.4 && pT_part2 > 0.4 && TMath::Min(pT_part1,pT_part2) < 1.1 "),
                        TString("pT_part1 > 0.2 && pT_part2 > 0.2 && pT_part3 > 0.2 && pT_part4 > 0.2 && eta_part3 > -0.7 && eta_part3 < 0.7 && eta_part4 > -0.7 && eta_part4 < 0.7 ")}; 

    Double_t binning[4][3] =  {{64, 0.3, 3.5},{44, 0.8, 3}, {24, 1.6, 4}, {50,0.5,4.5}};
    Int_t rapBins[4] = {16, 10, 8, 8};
    TString combCuts[] = { TString(""), TString("&& TMath::Abs(deltaphi_pp) > 1.570796327"), TString("&& TMath::Abs(deltaphi_pp) < 1.570796327")}; // 57.2957795 = 1 rad => 1.6 = 90 cca

    TString usedCuts = graniittiCuts;
    TString variable;
    Int_t nBins;
    Double_t min, max;
    TH1D* hist;
    for (int iPart = 0; iPart < nParticles + 1; ++iPart)
    {
        // Plot t 
        usedCuts+= " && " + partCuts[iPart];
        //cout<<usedCuts<<endl; 

        //////////////////////////////////////////////////
        // Plot phi 
        variable = "phi";
        nBins = 18;
        if(iPart > 1)
            nBins = 9;
        min = 0.0;
        max = 180.0;

        granitiTree[iPart]->Draw("TMath::Abs(deltaphi_pp*57.2957795)>>" + variable +"Sig1(" + nBins + "," + min + "," + max + ")",usedCuts);
        hist = (TH1D*)gPad->GetPrimitive(variable +"Sig1");   
        gDeltaPhi[iPart] = (TH1D*)hist->Clone("gDeltaPhi" + particleLables[iPart]);


        //////////////////////////////////////////////
        // Plot inv mass
        variable = "invMass_state";
        nBins = binning[iPart][0];
        min = binning[iPart][1];
        max = binning[iPart][2];

        granitiTree[iPart]->Draw(variable +">>" + variable +"Sig1(" + nBins + "," + min + "," + max + ")",usedCuts);
        hist = (TH1D*)gPad->GetPrimitive(variable +"Sig1");   
        gInvMass[iPart][0] = (TH1D*)hist->Clone("gInvMass" + particleLables[iPart]);
        //cout<<"Part "<<iPart<<" "<<gInvMass[iPart][0]<<endl;
        if(iPart < 3)
        {
            for (int comb = 0; comb < 2; ++comb)
            {
                TString cutsPhi = usedCuts + combCuts[comb+1];
                granitiTree[iPart]->Draw(variable +">>" + variable +"Sig1(" + nBins + "," + min + "," + max + ")",cutsPhi);
                hist = (TH1D*)gPad->GetPrimitive(variable +"Sig1"); 
                gInvMass[iPart][comb+1] = (TH1D*)hist->Clone("gDeltaPhi" + particleLables[iPart] + combinationLabel[comb+1]);
            }
        }    

        //////////////////////////////////////////////
        // Plot pair rapidity
        variable = "rap_state";
        nBins = rapBins[iPart];
        min = -0.8;
        max = 0.8;

        granitiTree[iPart]->Draw(variable +">>" + variable +"Sig1(" + nBins + "," + min + "," + max + ")",usedCuts);
        hist = (TH1D*)gPad->GetPrimitive(variable +"Sig1");   
        gRap[iPart] = (TH1D*)hist->Clone("gRap" + particleLables[iPart]);

    }

}

void Make(int signal)
{
    double effTotal, effTPC, effTOF;
    unsigned int PID;
   // cout<< vertexesZ[0] <<" "<< NhitsFit[0]<<" "<<NhitsFit[1] <<" "<< NhitsDEdx[0]<<" "<<NhitsDEdx[1] <<" "<<DcaZ[0] <<" "<<DcaZ[1] <<" "<<DcaXY[0] <<" "<<DcaXY[1] <<" "<<Eta[0] <<" "<<Eta[1] <<" "<< !fourPiState<<endl;
            
    double deltaPhi = TMath::Abs(phiRp[East] - phiRp[West])*convertToDegree;
    if(deltaPhi > 180)
        deltaPhi = 360 - deltaPhi;
    hDeltaPhi[0]->Fill(deltaPhi);
    if(elastic)
        hDeltaPhi[1]->Fill(deltaPhi);
    else
        hDeltaPhi[2]->Fill(deltaPhi);

    if(vertexesZ[0] < 80 && vertexesZ[0] > -80 && NhitsFit[0] >=25 && NhitsFit[1] >= 25 && NhitsDEdx[0] >= 15 && NhitsDEdx[1] >= 15 && DcaZ[0] < 1 && DcaZ[0] > -1 && DcaZ[1] < 1 && DcaZ[1] > -1 && DcaXY[0] < 1.5 && DcaXY[1] < 1.5 && Eta[0] > -0.7 && Eta[0] < 0.7 && Eta[1] > -0.7 && Eta[1] < 0.7 &&  t[0] < -0.12 && t[1] < -0.12 && t[0] > -1.0  && t[1] > -1.0 && !fourPiState)
    {
 


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
                hInvMassCorr[3][signal][0]->Fill(invMass[Pion], 1/effTotal);

                hRapidityCorr[3][signal][0]->Fill(pairRapidity, 1/effTotal);

                hDeltaPhiCorr[3][signal][0]->Fill(deltaPhi, 1/effTotal);
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
        hInvMassCorr[3][0][i]  = new TH1D("corrInvMass" + particleLables[3] + combinationLabel[i] + "Sig", "Corrected inv. mass " + particleLables[3] + combinationLabel[i], 50, 0.5, 4.5);
        hInvMassCorr[0][1][i]  = new TH1D("corrInvMass" + particleLables[0] + combinationLabel[i] + "Bcg", "Corrected inv. mass " + particleLables[0] + combinationLabel[i], 64, 0.3, 3.5);
        hInvMassCorr[1][1][i]  = new TH1D("corrInvMass" + particleLables[1] + combinationLabel[i] + "Bcg", "Corrected inv. mass " + particleLables[1] + combinationLabel[i], 44, 0.8, 3);
        hInvMassCorr[2][1][i]  = new TH1D("corrInvMass" + particleLables[2] + combinationLabel[i] + "Bcg", "Corrected inv. mass " + particleLables[2] + combinationLabel[i], 24, 1.6, 4);
        hInvMassCorr[3][1][i]  = new TH1D("corrInvMass" + particleLables[3] + combinationLabel[i] + "Bcg", "Corrected inv. mass " + particleLables[3] + combinationLabel[i], 50, 0.5, 4.5);

        hRapidityCorr[0][0][i]  = new TH1D("corrRapidity" + particleLables[0] + combinationLabel[i] + "Sig", "Corrected rapidity " + particleLables[0] + combinationLabel[i], 16, -0.8, 0.8);
        hRapidityCorr[1][0][i]  = new TH1D("corrRapidity" + particleLables[1] + combinationLabel[i] + "Sig", "Corrected rapidity " + particleLables[1] + combinationLabel[i], 10, -0.8, 0.8);
        hRapidityCorr[2][0][i]  = new TH1D("corrRapidity" + particleLables[2] + combinationLabel[i] + "Sig", "Corrected rapidity " + particleLables[2] + combinationLabel[i], 8, -0.8, 0.8);
        hRapidityCorr[3][0][i]  = new TH1D("corrRapidity" + particleLables[3] + combinationLabel[i] + "Sig", "Corrected rapidity " + particleLables[3] + combinationLabel[i], 8, -0.8, 0.8);
        hRapidityCorr[0][1][i]  = new TH1D("corrRapidity" + particleLables[0] + combinationLabel[i] + "Bcg", "Corrected rapidity " + particleLables[0] + combinationLabel[i], 16, -0.8, 0.8);
        hRapidityCorr[1][1][i]  = new TH1D("corrRapidity" + particleLables[1] + combinationLabel[i] + "Bcg", "Corrected rapidity " + particleLables[1] + combinationLabel[i], 10, -0.8, 0.8);
        hRapidityCorr[2][1][i]  = new TH1D("corrRapidity" + particleLables[2] + combinationLabel[i] + "Bcg", "Corrected rapidity " + particleLables[2] + combinationLabel[i], 8, -0.8, 0.8);
        hRapidityCorr[3][1][i]  = new TH1D("corrRapidity" + particleLables[3] + combinationLabel[i] + "Bcg", "Corrected rapidity " + particleLables[3] + combinationLabel[i], 8, -0.8, 0.8);

        hDeltaPhiCorr[0][0][i]  = new TH1D("corrDeltaPhi" + particleLables[0] + combinationLabel[i] + "Sig", "Corrected delta phi " + particleLables[0] + combinationLabel[i], 18, 0.0, 180);
        hDeltaPhiCorr[1][0][i]  = new TH1D("corrDeltaPhi" + particleLables[1] + combinationLabel[i] + "Sig", "Corrected delta phi " + particleLables[1] + combinationLabel[i], 18, 0.0, 180);
        hDeltaPhiCorr[2][0][i]  = new TH1D("corrDeltaPhi" + particleLables[2] + combinationLabel[i] + "Sig", "Corrected delta phi " + particleLables[2] + combinationLabel[i], 9, 0.0, 180);
        hDeltaPhiCorr[3][0][i]  = new TH1D("corrDeltaPhi" + particleLables[3] + combinationLabel[i] + "Sig", "Corrected delta phi " + particleLables[3] + combinationLabel[i], 9, 0.0, 180);
        hDeltaPhiCorr[0][1][i]  = new TH1D("corrDeltaPhi" + particleLables[0] + combinationLabel[i] + "Bcg", "Corrected delta phi " + particleLables[0] + combinationLabel[i], 18, 0.0, 180);
        hDeltaPhiCorr[1][1][i]  = new TH1D("corrDeltaPhi" + particleLables[1] + combinationLabel[i] + "Bcg", "Corrected delta phi " + particleLables[1] + combinationLabel[i], 18, 0.0, 180);
        hDeltaPhiCorr[2][1][i]  = new TH1D("corrDeltaPhi" + particleLables[2] + combinationLabel[i] + "Bcg", "Corrected delta phi " + particleLables[2] + combinationLabel[i], 9, 0.0, 180);
        hDeltaPhiCorr[3][1][i]  = new TH1D("corrDeltaPhi" + particleLables[3] + combinationLabel[i] + "Bcg", "Corrected delta phi " + particleLables[3] + combinationLabel[i], 9, 0.0, 180);

     }

/*
    for (int i = 0; i < 4; ++i)
    {
        gDeltaPhi[i] = new TH1D("granDeltaPhi" + particleLables[i], "granDeltaPhi");
        gRap[i] = new TH1D("granRap" + particleLables[i], "granRap");
        for (int j = 0; j < 2; ++j)
            gInvMassPhi[i][j] = new TH1D("granInvMass" + particleLables[i] + combinationLabel[j+1], "granInvMass " + particleLables[i] + combinationLabel[j+1]);

        gInvMass[i] = new TH1D("granInvMass" + particleLables[i], "granInvMass " + particleLables[i]);
    }
*/
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
    gPad->SetMargin(0.15,0.04,0.105,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
    gPad->SetTickx();
    gPad->SetTicky(); 
    gPad->SetLogy(0);
    gStyle->SetOptStat("");
    gStyle->SetPalette(1);
    gStyle->SetLineWidth(2);      //axis line
    gStyle->SetTickLength(0.04,"Y");
    gStyle->SetTickLength(0.04,"X");
    gStyle->SetFrameLineWidth(2); //frame line

    gPad->SetMargin(0.15,0.04,0.3,0.02);
    gStyle->SetPadTickY(1);
    // Upper plot will be in pad1
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.295, 1, 1.0);
    pad1->SetTickx();
    pad1->SetTicky();
    pad1->SetTopMargin(0.04);
    pad1->SetRightMargin(0.04);
    pad1->SetLeftMargin(0.12);
    pad1->SetBottomMargin(0); // Upper and lower plot are joined
    gStyle->SetFrameLineWidth(2); //frame line
    gStyle->SetTickLength(0.04,"Y");
    gStyle->SetTickLength(0.04,"X");
    pad1->Draw();             // Draw the upper pad: pad1

    newCanvas->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.00, 1, 0.2985);
    pad2->SetTickx();
    pad2->SetTicky();
    pad2->SetTopMargin(0);
    pad2->SetRightMargin(0.04);
    pad2->SetLeftMargin(0.12);
    pad2->SetBottomMargin(0.3);
    pad2->Draw();

    TH1D *hist, *histCompare;
    TString name = "phiPlot";
    TString yLabelDef = "Probability per event";
    TString yLabel;
    Double_t nBins[4] = {10, 10, 20, 20};
    Double_t binning[4][2] =  {{0.0, 180.0},{0.0, 180.0}, {0.0, 180.0}, {0.0, 180.0}};
    Double_t limits[4][2] = { {0.01, 2.4}, {0.01, 2.4}, {0.01, 2.4}, {0.01, 2.4}}; 
    for (int i = 0; i < nParticles +1; ++i)
    { 
        int comb = 0;
        yLabel = yLabelDef + " / " + nBins[i] + " deq";

        pad1->cd(); 
        hist = (TH1D*)hDeltaPhiCorr[i][0][comb]->Clone("hist");
        hist -> Add(hDeltaPhiCorr[i][1][comb], -1); 
         //cout<<"Part FURTHER "<<i<<" "<<gInvMass[i][0]<<endl;
        histCompare = (TH1D*)gDeltaPhi[i]->Clone("histCompare");
        // Getting histogram created from TTree making warning when use Divide

        hist->SetMarkerColor(4);
        hist->SetMarkerSize(1);
        hist->SetMarkerStyle(20);
        hist->SetLineColor(4);
        hist->SetLineStyle(1);
        hist->SetLineWidth(1);
        hist->SetStats(false);

        histCompare->SetMarkerColor(1);
        histCompare->SetMarkerSize(1);
        histCompare->SetMarkerStyle(21);
        histCompare->SetLineColor(1);
        histCompare->SetLineStyle(1);
        histCompare->SetLineWidth(1);
        histCompare->SetStats(false);

        hist->GetXaxis()->SetTitleFont(43);
        hist->GetYaxis()->SetTitleFont(43);
        hist->GetXaxis()->SetTitleSize(30);
        hist->GetYaxis()->SetTitleSize(30);
        hist->GetXaxis()->SetTitleOffset(1.0);
        hist->GetYaxis()->SetTitleOffset(1.4);
        hist->GetYaxis()->SetLabelFont(43);
        hist->GetYaxis()->SetLabelSize(30);

        //hist->GetYaxis()->SetTickLength(0.04);
        hist->GetXaxis()->SetTickLength(0.04);
        //histCompare->GetYaxis()->SetTickLength(0.04);
        histCompare->GetXaxis()->SetTickLength(0.04);

        Double_t scaleFactor =   1 /hist->Integral();
        hist->Scale(scaleFactor);
        scaleFactor =   1 /histCompare->Integral();
        histCompare->Scale(scaleFactor);

        hist->GetYaxis()->SetRangeUser(0.0001, max(hist->GetMaximum(),histCompare->GetMaximum())*1.2);

        hist->SetTitle(" ; #Delta#varphi [deg] ; " + yLabel);
        hist->Draw();              
        histCompare->Draw("same");         

        TPaveText *textSTAR = new TPaveText(0.7,0.85,0.85,0.91,"brNDC"); 
        textSTAR -> SetTextSize(0.055);
        textSTAR -> SetFillColor(0);
        textSTAR -> SetTextFont(62);
        textSTAR -> AddText("THIS THESIS");
        textSTAR -> Draw("same");

        textSTAR = new TPaveText(0.7,0.65,0.85,0.8,"brNDC"); 
        textSTAR -> SetTextSize(30);
        textSTAR -> SetFillColor(0);
        textSTAR -> SetTextFont(43);
        textSTAR -> AddText("p + p #rightarrow p + " + stateLabel[i] +" + p");
        textSTAR -> AddText("#sqrt{s} = 510 GeV");
        textSTAR -> Draw("same");

        TLegend *legendPID = new TLegend(0.68,0.45,0.92,0.62,"","brNDC");
        legendPID->SetFillStyle(0);
        legendPID->SetBorderSize(0);
        legendPID -> SetTextSize(30);
        legendPID -> SetTextFont(43);
        legendPID -> AddEntry(hist, "Data", "p");
        legendPID -> AddEntry(histCompare, "Graniitti", "p");
        legendPID -> Draw("same");

        textSTAR = new TPaveText(0.68,0.32,0.92,0.44,"brNDC"); 
        textSTAR -> SetTextSize(30);
        textSTAR -> SetFillColor(0);
        textSTAR -> SetTextFont(43);
        textSTAR -> AddText("Like-sign background");
        textSTAR -> AddText("subtracted");
        textSTAR -> Draw("same");

        if(comb != 0)
        {
            TString phiLabel = "#Delta#varphi > 90^{#circ}";
            if(comb == Inel)
                phiLabel = "#Delta#varphi < 90^{#circ}";
            textSTAR = new TPaveText(0.68,0.35,0.92,0.45,"brNDC");
            textSTAR -> SetTextSize(30);
            textSTAR -> SetTextAlign(22);
            textSTAR -> SetFillColor(0);
            textSTAR -> SetTextFont(43);
            textSTAR -> AddText(phiLabel);
            textSTAR -> Draw("same");
        }

        pad2->cd();       

        // Define the ratio plot
        TH1D *h3 = (TH1D*)hist->Clone("h3");
        h3->GetYaxis()->SetLabelSize(0.);
        h3->SetMinimum(limits[i][0]);  // Define Y ..
        h3->SetMaximum(limits[i][1]); // .. range
        h3->Divide(histCompare);
        h3->SetMarkerColor(1);
        h3->SetMarkerSize(1);
        h3->SetMarkerStyle(20);
        h3->SetLineColor(1);
        h3->SetLineStyle(1);
        h3->SetLineWidth(1);

        // Ratio plot (h3) settings
        h3->SetTitle(""); // Remove the ratio title

        // Y axis ratio plot settings
        h3->GetYaxis()->SetTitle("Data/ Graniitti");
        h3->GetYaxis()->SetNdivisions(504);
        h3->GetYaxis()->SetTitleSize(30);
        h3->GetYaxis()->SetTitleFont(43);
        h3->GetYaxis()->SetTitleOffset(1.4);
        h3->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
        h3->GetYaxis()->SetLabelSize(0.0);

        // X axis ratio plot settings
        h3->GetXaxis()->SetTitleSize(30);
        h3->GetXaxis()->SetTitleFont(43);
        h3->GetXaxis()->SetTitleOffset(3.2);
        h3->GetXaxis()->SetLabelOffset(0.025);
        h3->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
        h3->GetXaxis()->SetLabelSize(30);
        h3->GetYaxis()->SetTickLength(0.03);
        h3->GetXaxis()->SetTickLength(0.08);
        h3->Draw("ep");       // Draw the ratio plot

        TGaxis *axis = new TGaxis(binning[i][0], limits[i][0], binning[i][0], limits[i][1], limits[i][0], limits[i][1], 504, "");
        axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
        axis->SetLabelSize(30);
        axis->SetLabelOffset(0.01);
        axis->SetTickLength(0.08);
        gStyle->SetLineWidth(1);      //axis line
        gStyle->SetFrameLineWidth(2); //frame line
        axis->Draw("same");

        TLine *unity = new TLine(binning[i][0], 1., binning[i][1] , 1.);
        unity->SetLineColor(kBlack);
        unity->SetLineWidth(2);
        unity->Draw();

        newCanvas->Update();
        newCanvas->Write(name + particleLables[i] + combinationLabel[comb]);
        
    } 

}

void PlotRap()
{
    TCanvas* newCanvas = new TCanvas("newCanvas","newCanvas",800,700);
    gPad->SetMargin(0.15,0.04,0.105,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
    gPad->SetTickx();
    gPad->SetTicky(); 
    gPad->SetLogy(0);
    gStyle->SetOptStat("");
    gStyle->SetPalette(1);
    gStyle->SetLineWidth(2);      //axis line
    gStyle->SetTickLength(0.04,"Y");
    gStyle->SetTickLength(0.04,"X");
    gStyle->SetFrameLineWidth(2); //frame line

    gPad->SetMargin(0.15,0.04,0.3,0.02);
    gStyle->SetPadTickY(1);
    // Upper plot will be in pad1
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.295, 1, 1.0);
    pad1->SetTickx();
    pad1->SetTicky();
    pad1->SetTopMargin(0.04);
    pad1->SetRightMargin(0.04);
    pad1->SetLeftMargin(0.12);
    pad1->SetBottomMargin(0); // Upper and lower plot are joined
    gStyle->SetFrameLineWidth(2); //frame line
    gStyle->SetTickLength(0.04,"Y");
    gStyle->SetTickLength(0.04,"X");
    pad1->Draw();             // Draw the upper pad: pad1

    newCanvas->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.00, 1, 0.2985);
    pad2->SetTickx();
    pad2->SetTicky();
    pad2->SetTopMargin(0);
    pad2->SetRightMargin(0.04);
    pad2->SetLeftMargin(0.12);
    pad2->SetBottomMargin(0.3);
    pad2->Draw();




    TH1D *hist, *histCompare;
    TString name = "rapPlot";
    TString yLabelDef = "Probability per event";
    TString yLabel;
    Double_t nBins[4] = {0.1, 0.16, 0.2, 0.2};
    Double_t binning[4][2] =  {{-0.8, 0.8}, {-0.8, 0.8}, {-0.8, 0.8}, {-0.8, 0.8}};
    Double_t limits[4][2] = { {0.01, 2.4}, {0.01, 2.4}, {0.01, 2.4}, {0.01, 2.4}}; 
    for (int i = 0; i < nParticles +1; ++i)
    { 
        yLabel = yLabelDef + " / " + nBins[i];
        int comb = 0;
        if(i == nParticles && comb > 0)
            continue;

        pad1->cd();
        hist = (TH1D*)hRapidityCorr[i][0][comb]->Clone("hist");
        hist -> Add(hRapidityCorr[i][1][comb], -1); 
         //cout<<"Part FURTHER "<<i<<" "<<gInvMass[i][0]<<endl;
        histCompare = (TH1D*)gRap[i]->Clone("histCompare");
        // Getting histogram created from TTree making warning when use Divide

        hist->SetMarkerColor(4);
        hist->SetMarkerSize(1);
        hist->SetMarkerStyle(20);
        hist->SetLineColor(4);
        hist->SetLineStyle(1);
        hist->SetLineWidth(1);
        hist->SetStats(false);

        histCompare->SetMarkerColor(1);
        histCompare->SetMarkerSize(1);
        histCompare->SetMarkerStyle(21);
        histCompare->SetLineColor(1);
        histCompare->SetLineStyle(1);
        histCompare->SetLineWidth(1);
        histCompare->SetStats(false);

        hist->GetXaxis()->SetTitleFont(43);
        hist->GetYaxis()->SetTitleFont(43);
        hist->GetXaxis()->SetTitleSize(30);
        hist->GetYaxis()->SetTitleSize(30);
        hist->GetXaxis()->SetTitleOffset(1.0);
        hist->GetYaxis()->SetTitleOffset(1.4);
        hist->GetYaxis()->SetLabelFont(43);
        hist->GetYaxis()->SetLabelSize(30);

        //hist->GetYaxis()->SetTickLength(0.04);
        hist->GetXaxis()->SetTickLength(0.04);
        //histCompare->GetYaxis()->SetTickLength(0.04);
        histCompare->GetXaxis()->SetTickLength(0.04);

        Double_t scaleFactor =   1 /hist->Integral();
        hist->Scale(scaleFactor);
        scaleFactor =   1 /histCompare->Integral();
        histCompare->Scale(scaleFactor);

        hist->GetYaxis()->SetRangeUser(0.0001, max(hist->GetMaximum(),histCompare->GetMaximum())*1.2);

        hist->SetTitle(" ; y(" + stateLabel[i] + "); " + yLabel);
        hist->Draw();              
        histCompare->Draw("same");         

        TPaveText *textSTAR = new TPaveText(0.7,0.85,0.85,0.91,"brNDC"); 
        textSTAR -> SetTextSize(0.055);
        textSTAR -> SetFillColor(0);
        textSTAR -> SetTextFont(62);
        textSTAR -> AddText("THIS THESIS");
        textSTAR -> Draw("same");

        textSTAR = new TPaveText(0.7,0.65,0.85,0.8,"brNDC"); 
        textSTAR -> SetTextSize(30);
        textSTAR -> SetFillColor(0);
        textSTAR -> SetTextFont(43);
        textSTAR -> AddText("p + p #rightarrow p + " + stateLabel[i] +" + p");
        textSTAR -> AddText("#sqrt{s} = 510 GeV");
        textSTAR -> Draw("same");

        TLegend *legendPID = new TLegend(0.68,0.45,0.92,0.62,"","brNDC");
        legendPID->SetFillStyle(0);
        legendPID->SetBorderSize(0);
        legendPID -> SetTextSize(30);
        legendPID -> SetTextFont(43);
        legendPID -> AddEntry(hist, "Data", "p");
        legendPID -> AddEntry(histCompare, "Graniitti", "p");
        legendPID -> Draw("same");

        textSTAR = new TPaveText(0.68,0.32,0.92,0.44,"brNDC"); 
        textSTAR -> SetTextSize(30);
        textSTAR -> SetFillColor(0);
        textSTAR -> SetTextFont(43);
        textSTAR -> AddText("Like-sign background");
        textSTAR -> AddText("subtracted");
        textSTAR -> Draw("same");

        if(comb != 0)
        {
            TString phiLabel = "#Delta#varphi > 90^{#circ}";
            if(comb == Inel)
                phiLabel = "#Delta#varphi < 90^{#circ}";
            textSTAR = new TPaveText(0.68,0.35,0.92,0.45,"brNDC");
            textSTAR -> SetTextSize(30);
            textSTAR -> SetTextAlign(22);
            textSTAR -> SetFillColor(0);
            textSTAR -> SetTextFont(43);
            textSTAR -> AddText(phiLabel);
            textSTAR -> Draw("same");
        }

        pad2->cd();       

        // Define the ratio plot
        TH1D *h3 = (TH1D*)hist->Clone("h3");
        h3->GetYaxis()->SetLabelSize(0.);
        h3->SetMinimum(limits[i][0]);  // Define Y ..
        h3->SetMaximum(limits[i][1]); // .. range
        h3->Divide(histCompare);
        h3->SetMarkerColor(1);
        h3->SetMarkerSize(1);
        h3->SetMarkerStyle(20);
        h3->SetLineColor(1);
        h3->SetLineStyle(1);
        h3->SetLineWidth(1);

        // Ratio plot (h3) settings
        h3->SetTitle(""); // Remove the ratio title

        // Y axis ratio plot settings
        h3->GetYaxis()->SetTitle("Data/ Graniitti");
        h3->GetYaxis()->SetNdivisions(504);
        h3->GetYaxis()->SetTitleSize(30);
        h3->GetYaxis()->SetTitleFont(43);
        h3->GetYaxis()->SetTitleOffset(1.4);
        h3->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
        h3->GetYaxis()->SetLabelSize(0.0);

        // X axis ratio plot settings
        h3->GetXaxis()->SetTitleSize(30);
        h3->GetXaxis()->SetTitleFont(43);
        h3->GetXaxis()->SetTitleOffset(3.2);
        h3->GetXaxis()->SetLabelOffset(0.025);
        h3->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
        h3->GetXaxis()->SetLabelSize(30);
        h3->GetYaxis()->SetTickLength(0.03);
        h3->GetXaxis()->SetTickLength(0.08);
        h3->Draw("ep");       // Draw the ratio plot

        TGaxis *axis = new TGaxis(binning[i][0], limits[i][0], binning[i][0], limits[i][1], limits[i][0], limits[i][1], 504, "");
        axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
        axis->SetLabelSize(30);
        axis->SetLabelOffset(0.01);
        axis->SetTickLength(0.08);
        gStyle->SetLineWidth(1);      //axis line
        gStyle->SetFrameLineWidth(2); //frame line
        axis->Draw("same");

        TLine *unity = new TLine(binning[i][0], 1., binning[i][1] , 1.);
        unity->SetLineColor(kBlack);
        unity->SetLineWidth(2);
        unity->Draw();

        newCanvas->Update();
        newCanvas->Write(name + particleLables[i] + combinationLabel[comb]);
        
    } 
}


void PlotMassPlots()
{

    TCanvas* newCanvas = new TCanvas("newCanvas","newCanvas",800,700);
    gPad->SetMargin(0.15,0.04,0.105,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
    gPad->SetTickx();
    gPad->SetTicky(); 
    gPad->SetLogy(0);
    gStyle->SetOptStat("");
    gStyle->SetPalette(1);
    gStyle->SetLineWidth(2);      //axis line
    gStyle->SetTickLength(0.04,"Y");
    gStyle->SetTickLength(0.04,"X");
    gStyle->SetFrameLineWidth(2); //frame line

    gPad->SetMargin(0.15,0.04,0.3,0.02);
    gStyle->SetPadTickY(1);
    // Upper plot will be in pad1
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.295, 1, 1.0);
    pad1->SetTickx();
    pad1->SetTicky();
    pad1->SetTopMargin(0.04);
    pad1->SetRightMargin(0.04);
    pad1->SetLeftMargin(0.12);
    pad1->SetBottomMargin(0); // Upper and lower plot are joined
    gStyle->SetFrameLineWidth(2); //frame line
    gStyle->SetTickLength(0.04,"Y");
    gStyle->SetTickLength(0.04,"X");
    pad1->Draw();             // Draw the upper pad: pad1

    newCanvas->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.00, 1, 0.2985);
    pad2->SetTickx();
    pad2->SetTicky();
    pad2->SetTopMargin(0);
    pad2->SetRightMargin(0.04);
    pad2->SetLeftMargin(0.12);
    pad2->SetBottomMargin(0.3);
    pad2->Draw();

    TH1D *hist, *histCompare;
    TString name = "massPlot";
    TString yLabelDef = "Probability per event";
    TString yLabel;
    Double_t nBins[4] = {50, 50, 100, 80};
    Double_t binning[4][2] =  {{0.3, 3.5},{0.8, 3}, {1.6, 4}, {0.5,4.5}};
    Double_t limits[4][2] = { {0.01, 2.4}, {0.01, 2.4}, {0.01, 2.4}, {0.01, 2.4}}; 
    for (int i = 0; i < nParticles +1; ++i)
    { 
        yLabel = yLabelDef + " / " + nBins[i] + " MeV";
        for (int comb = 0; comb < nCombination; ++comb)
        { 
            if(i == nParticles && comb > 0)
                continue;

            pad1->cd();
            hist = (TH1D*)hInvMassCorr[i][0][comb]->Clone("hist");
            hist -> Add(hInvMassCorr[i][1][comb], -1); 
             //cout<<"Part FURTHER "<<i<<" "<<gInvMass[i][0]<<endl;
            histCompare = (TH1D*)gInvMass[i][comb]->Clone("histCompare");
            // Getting histogram created from TTree making warning when use Divide

            hist->SetMarkerColor(4);
            hist->SetMarkerSize(1);
            hist->SetMarkerStyle(20);
            hist->SetLineColor(4);
            hist->SetLineStyle(1);
            hist->SetLineWidth(1);
            hist->SetStats(false);

            histCompare->SetMarkerColor(1);
            histCompare->SetMarkerSize(1);
            histCompare->SetMarkerStyle(21);
            histCompare->SetLineColor(1);
            histCompare->SetLineStyle(1);
            histCompare->SetLineWidth(1);
            histCompare->SetStats(false);

            hist->GetXaxis()->SetTitleFont(43);
            hist->GetYaxis()->SetTitleFont(43);
            hist->GetXaxis()->SetTitleSize(30);
            hist->GetYaxis()->SetTitleSize(30);
            hist->GetXaxis()->SetTitleOffset(1.0);
            hist->GetYaxis()->SetTitleOffset(1.4);
            hist->GetYaxis()->SetLabelFont(43);
            hist->GetYaxis()->SetLabelSize(30);

            //hist->GetYaxis()->SetTickLength(0.04);
            hist->GetXaxis()->SetTickLength(0.04);
            //histCompare->GetYaxis()->SetTickLength(0.04);
            histCompare->GetXaxis()->SetTickLength(0.04);

            Double_t scaleFactor =   1 /hist->Integral();
            hist->Scale(scaleFactor);
            scaleFactor =   1 /histCompare->Integral();
            histCompare->Scale(scaleFactor);

            hist->GetYaxis()->SetRangeUser(0.0001, max(hist->GetMaximum(),histCompare->GetMaximum())*1.2);

            hist->SetTitle(" ; m(" + stateLabel[i] + ") [GeV] ; " + yLabel);
            hist->Draw();              
            histCompare->Draw("same");         

            TPaveText *textSTAR = new TPaveText(0.7,0.85,0.85,0.91,"brNDC"); 
            textSTAR -> SetTextSize(0.055);
            textSTAR -> SetFillColor(0);
            textSTAR -> SetTextFont(62);
            textSTAR -> AddText("THIS THESIS");
            textSTAR -> Draw("same");

            textSTAR = new TPaveText(0.7,0.65,0.85,0.8,"brNDC"); 
            textSTAR -> SetTextSize(30);
            textSTAR -> SetFillColor(0);
            textSTAR -> SetTextFont(43);
            textSTAR -> AddText("p + p #rightarrow p + " + stateLabel[i] +" + p");
            textSTAR -> AddText("#sqrt{s} = 510 GeV");
            textSTAR -> Draw("same");

            textSTAR = new TPaveText(0.68,0.32,0.92,0.44,"brNDC"); 
            textSTAR -> SetTextSize(30);
            textSTAR -> SetFillColor(0);
            textSTAR -> SetTextFont(43);
            textSTAR -> AddText("Like-sign background");
            textSTAR -> AddText("subtracted");
            textSTAR -> Draw("same");

            TLegend *legendPID = new TLegend(0.68,0.45,0.92,0.62,"","brNDC");
            legendPID->SetFillStyle(0);
            legendPID->SetBorderSize(0);
            legendPID -> SetTextSize(30);
            legendPID -> SetTextFont(43);
            legendPID -> AddEntry(hist, "Data", "p");
            legendPID -> AddEntry(histCompare, "Graniitti", "p");
            legendPID -> Draw("same");

            if(comb != 0)
            {
                TString phiLabel = "#Delta#varphi > 90^{#circ}";
                if(comb == Inel)
                    phiLabel = "#Delta#varphi < 90^{#circ}";
                textSTAR = new TPaveText(0.68,0.35,0.92,0.45,"brNDC");
                textSTAR -> SetTextSize(30);
                textSTAR -> SetTextAlign(22);
                textSTAR -> SetFillColor(0);
                textSTAR -> SetTextFont(43);
                textSTAR -> AddText(phiLabel);
                textSTAR -> Draw("same");
            }

            pad2->cd();       

            // Define the ratio plot
            TH1D *h3 = (TH1D*)hist->Clone("h3");
            h3->GetYaxis()->SetLabelSize(0.);
            h3->SetMinimum(limits[i][0]);  // Define Y ..
            h3->SetMaximum(limits[i][1]); // .. range
            h3->Divide(histCompare);
            h3->SetMarkerColor(1);
            h3->SetMarkerSize(1);
            h3->SetMarkerStyle(20);
            h3->SetLineColor(1);
            h3->SetLineStyle(1);
            h3->SetLineWidth(1);

            // Ratio plot (h3) settings
            h3->SetTitle(""); // Remove the ratio title

            // Y axis ratio plot settings
            h3->GetYaxis()->SetTitle("Data/ Graniitti");
            h3->GetYaxis()->SetNdivisions(504);
            h3->GetYaxis()->SetTitleSize(30);
            h3->GetYaxis()->SetTitleFont(43);
            h3->GetYaxis()->SetTitleOffset(1.4);
            h3->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
            h3->GetYaxis()->SetLabelSize(0.0);

            // X axis ratio plot settings
            h3->GetXaxis()->SetTitleSize(30);
            h3->GetXaxis()->SetTitleFont(43);
            h3->GetXaxis()->SetTitleOffset(3.2);
            h3->GetXaxis()->SetLabelOffset(0.025);
            h3->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
            h3->GetXaxis()->SetLabelSize(30);
            h3->GetYaxis()->SetTickLength(0.03);
            h3->GetXaxis()->SetTickLength(0.08);
            h3->Draw("ep");       // Draw the ratio plot

            TGaxis *axis = new TGaxis(binning[i][0], limits[i][0], binning[i][0], limits[i][1], limits[i][0], limits[i][1], 504, "");
            axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
            axis->SetLabelSize(30);
            axis->SetLabelOffset(0.01);
            axis->SetTickLength(0.08);
            gStyle->SetLineWidth(1);      //axis line
            gStyle->SetFrameLineWidth(2); //frame line
            axis->Draw("same");

            TLine *unity = new TLine(binning[i][0], 1., binning[i][1] , 1.);
            unity->SetLineColor(kBlack);
            unity->SetLineWidth(2);
            unity->Draw();

            newCanvas->Update();
            newCanvas->Write(name + particleLables[i] + combinationLabel[comb]);
        }
    } 


}
