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
enum RP_ID {E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D, nRomanPots};

TString particleLables[nParticles] = { TString("Pion"), TString("Kaon"), TString("Proton")};
TString sideLabel[nSides] = { TString("East"), TString("West")};
TString rpNames[nRomanPots] = { TString("E1U"), TString("E1D"), TString("E2U"), TString("E2D"), TString("W1U"), TString("W1D"), TString("W2U"), TString("W2D")};

TFile* data;
TFile* fout;
TFile* TPCeff;
TFile* TOFeff;

TTree* recTree;

TH3D* hTOFeff[6]; // 0 = pi- 1 = K- 2 = pbar
TH3D* hTPCeff[6]; // 3 = pi+ 4 = K+ 5 = p

TH1D* hInvMass[nParticles];
TH1D* hDeltaTime;
TH1D* hvertexZ;
TH1D* hvertexRP;
TH1D* hvertexRPGolden;
TH1D* hDeltaVertex;

TH1D* hEtaDif;
TH1D* hLFactor;

TH1D* pTGoldenPions;

TH1D* hEtaDist;
TH1D* hEta[20];
TH1D* hEtaDistCalib;

TH1D* hTime[2];
TH1D* hTofLength[2];
TH1D* hRunNumber[2];

TH2D* hSlewing[16];
TH2D* hPosition[nSides];
TH2D* hdEdx;
TH2D* hMomentum[nParticles];
TH2D* hTransMomentum[nParticles];

TH1D* hchiPairPion[nParticles];
TH1D* hchiPairKaon[nParticles];
TH1D* hchiPairProton[nParticles];
TH1D* hMSquered[nParticles];
TH1D* hPtCharged[2];
TH1D* hRatio;

TH2D* hDEdxRafal[nParticles], *hDEdxTomas[nParticles], *hDEdxTomas2[nParticles], *hDEdxDaniel[nParticles];


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

///////////////////
UInt_t	nTOFBadTracks, etaTotal;
UInt_t	nTOFGoodTracks;
UInt_t	nWierdTOFInfo, triggerOn, triggerOff;

void Init();
void ConnectInput();
void Make();

void analysis()
{
	TString input = "/home/truhlar/Desktop/STAR/CEP/Analysis/Data/deltaPhi.root";
    TString TPCeffInput = "/home/truhlar/Desktop/STAR/CEP/Analysis/Data/etaPhiEfficiency_16_01_19_delta015_twoRuns.root";
    TString TOFeffInput = "/home/truhlar/Desktop/STAR/CEP/Analysis/Data/effWithBinningForSystematics.root";
	TString output = "/home/truhlar/Desktop/STAR/CEP/Analysis/Outputs/anaOutput.root";

    TPCeff = TFile::Open(TPCeffInput, "read");
    if (!TPCeff)
    {
        cout<<"Error: cannot open "<<TPCeffInput<<endl;
        return;
    }

    TOFeff = TFile::Open(TOFeffInput, "read");
    if (!TOFeff)
    {
        cout<<"Error: cannot open "<<TOFeffInput<<endl;
        return;
    }

	data = TFile::Open(input, "read");
	if (!data)
	{
		cout<<"Error: cannot open "<<input<<endl;
		return;
	}


    for (int i = 0; i < 6; ++i)
    {    
        hTPCeff[i] = (TH3D*)TPCeff -> Get(Form("hTPCEffiCD%i121",i));
        hTOFeff[i] = (TH3D*)TOFeff -> Get(Form("hTOFEffiCD%i12",i)); 
    }


	fout = new TFile(output,"RECREATE");
	Init(); // Preparing histograms 
	ConnectInput(); // Connecting input

    TCanvas *cCanvas = new TCanvas("cCanvas","cCanvas",800,700);
    gPad->SetMargin(0.12,0.02,0.1,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top) 
    gStyle->SetOptStat("");

    TH1F* histSignal;
    TH1F* histBackground;

    TString variable = "deltaPhi";
    TString cuts = "!fourPiState && vertexZ<80 && vertexZ > -80 && NhitsFit1 >=25 && NhitsFit0 >= 25 && NhitsDEdx1 >= 15 && NhitsDEdx0 >= 15 && DcaZ1 < 1 && DcaZ1 > -1 && DcaZ0 < 1 && DcaZ0 > -1 && DcaXY1 < 1.5 && DcaXY0 < 1.5 && Eta1 > -0.7 && Eta1 < 0.7 && Eta0 > -0.7 && Eta0 < 0.7 &&  tEast < -0.12 && tWest < -0.12 && tEast > -1.0  && tWest > -1.0";
    TString usedCuts = cuts + "&& !elastic";
    int nBins = 100;
    double min = 0.0;
    double max = 200.0;
    recTree->Draw("TMath::Abs( " + variable +"*57.2957795)>>" + variable +"Bcg1(" + nBins + "," + min + "," + max + ")",usedCuts);
    histBackground = (TH1F*)gPad->GetPrimitive(variable +"Bcg1");
    histBackground->SetMarkerColor(2);
    histBackground->SetMarkerSize(1);
    histBackground->SetMarkerStyle(20);


    usedCuts = cuts + "&& elastic";
    recTree->Draw("TMath::Abs( " + variable +"*57.2957795)>>" + variable +"Sig1(" + nBins + "," + min + "," + max + ")",usedCuts);
    histSignal = (TH1F*)gPad->GetPrimitive(variable +"Sig1");   
    histSignal->SetTitle(" ; |#Delta #varphi| [deg]; Number of events");
    histSignal->GetXaxis()->SetTitleFont(42);
    histSignal->GetYaxis()->SetTitleFont(42);
    histSignal->GetXaxis()->SetLabelFont(42);
    histSignal->GetYaxis()->SetLabelFont(42);
    histSignal->GetXaxis()->SetTitleSize(0.045);
    histSignal->GetYaxis()->SetTitleSize(0.045);
    histSignal->GetXaxis()->SetTitleOffset(0.9);
    histSignal->GetYaxis()->SetTitleOffset(1.3);
    histSignal->GetYaxis()->SetRangeUser(5.0, TMath::Max(histSignal->GetMaximum(),histBackground->GetMaximum())*1.1);
    histSignal->SetMarkerColor(4);
    histSignal->SetMarkerSize(1);
    histSignal->SetMarkerStyle(20);
    histSignal->SetLineColor(4);
    histSignal->SetLineStyle(1);
    histSignal->SetLineWidth(1);
    histSignal->Draw("E");
    histBackground->Draw("ESAME");

    TPaveText *textPub = new TPaveText(0.7,0.75,0.92,0.88,"brNDC");
    textPub -> SetTextSize(0.04);
    textPub -> SetTextAlign(22);
    textPub -> SetFillColor(0);
    textPub -> SetTextFont(42);
    textPub -> AddText("p + p #rightarrow p + X + p");
    textPub -> AddText("#sqrt{s} = 510 GeV");
    textPub -> Draw("same");

    TPaveText *textSTAR;
    textSTAR = new TPaveText(0.75,0.89,0.9,0.95,"brNDC");
    textSTAR -> SetTextSize(0.04);
    textSTAR -> SetFillColor(0);
    textSTAR -> SetTextFont(62);
    textSTAR->AddText("STAR Internal");
    textSTAR -> Draw("same");

    TLegend* leg1 = new TLegend(0.6, 0.65, 0.78, 0.74);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.04);
    leg1->SetTextFont(42);
    leg1->AddEntry(histSignal,"Elastic","p");
    leg1->AddEntry(histBackground,"Inelastic","p");
    leg1->Draw("same");


    cCanvas->Update();
    cCanvas->Write(variable); 


    TH1D* hEta0;
/*    for(int i = 0; i < 21; ++i)
    {
        for(int j = 0; j < 21; ++j)
        {
            hEta0 = hTOFeff[0]->ProjectionZ("",i,i,j,j);
            hEta0->SetTitle( Form(" TOF eff %d, %d, %d, %d ; #eta ; ",i,i,j,j));
            hEta0->Draw("");  
            cCanvas->Update();
            cCanvas->SaveAs( Form("etaEffTOF/%d%d%d%d.png",i,i,j,j));
        }
    }
    hEta0 = hTPCeff[0]->ProjectionZ("",11,11,10,10);
    //TH1D* hEta = hTOFeff[3]->ProjectionZ();
    hEta0->SetTitle(" TPC eff ; #eta ; ");
    hEta0->Draw("HIST");  
    return;

	////////////// Making histograms
	nTOFBadTracks = 0;
	nTOFGoodTracks = 0;
	nWierdTOFInfo = 0;
	Long64_t nev = recTree->GetEntries();
	cout<<"Proccesing "<<nev<<" events"<<endl;
	//event loop
	//nev = 1000;
    etaTotal = 0;
	for(Long64_t iev=0; iev<nev; ++iev) 
	{ //get the event
		recTree->GetEntry(iev); 
		Make();
	} 

    for(int i = 0; i < 20; ++i)
    {
        hLFactor->SetBinContent(i+1, hEta[i]->GetEntries()/etaTotal);
        hLFactor->GetXaxis()->SetBinLabel(i+1, Form("%.2f",0.6 + i*0.02));
    }
    gStyle->SetOptStat("");
    hLFactor->SetTitle("; |#eta| cut; f_{Lost}");

	Int_t totalTracks = nTOFGoodTracks+nTOFBadTracks+nWierdTOFInfo;
	cout<<"Total number of tracks: "<< totalTracks <<endl;
	cout<<"Tracks with good TOF inof: "<<nTOFGoodTracks<<" "<< (nTOFGoodTracks*100.0)/totalTracks<<endl;
	cout<<"Tracks with bad TOF inof: "<<nTOFBadTracks<<" "<< (nTOFBadTracks*100.0)/totalTracks<<endl;
	cout<<"Tracks with just one bad TOF inof: "<<nWierdTOFInfo<<" "<< (nWierdTOFInfo*100.0)/totalTracks<<endl;
    cout<<"Trigger on: "<<triggerOn<<" and off: "<<triggerOff<<endl;
    fout->cd();


    hDEdxTomas[0]->SetTitle(" ; #frac{q}{e} #times p_{T} [GeV/c] ;dE/dx [keV/cm]");
    hDEdxTomas[0]->SetStats(0);
    hDEdxTomas[0]->SetMarkerColor(2);
    hDEdxTomas[0]->SetMarkerSize(1);
    hDEdxTomas[0]->SetMarkerStyle(29);
    hDEdxTomas[0]->Draw("SCAT");
    hDEdxTomas[1]->SetMarkerColor(3);
    hDEdxTomas[1]->SetMarkerSize(1);
    hDEdxTomas[1]->SetMarkerStyle(29);
    hDEdxTomas[1]->Draw("SAME0");
    hDEdxTomas[2]->SetMarkerColor(4);
    hDEdxTomas[2]->SetMarkerSize(1);
    hDEdxTomas[2]->SetMarkerStyle(29);
    hDEdxTomas[2]->Draw("SAME0");
    TLegend *legendPID = new TLegend(0.72,0.55,0.85,0.80,"","brNDC");
    legendPID -> AddEntry(hDEdxTomas[0], "#pi^{+} + #pi^{-}", "p");
    legendPID -> AddEntry(hDEdxTomas[1], "K^{+} + K^{-}", "p");
    legendPID -> AddEntry(hDEdxTomas[2], "p + #bar{p}", "p");
    legendPID -> Draw("same");

    cCanvas->Update();
    cCanvas->Write("dEdxByTomas");

    hDEdxRafal[0]->SetStats(0);
    hDEdxRafal[0]->SetTitle(" ; #frac{q}{e} #times p_{T} [GeV/c] ;dE/dx [keV/cm]");
    hDEdxRafal[0]->SetMarkerColor(2);
    hDEdxRafal[0]->SetMarkerSize(1);
    hDEdxRafal[0]->SetMarkerStyle(29);
    hDEdxRafal[0]->Draw("SCAT");
    hDEdxRafal[1]->SetMarkerColor(3);
    hDEdxRafal[1]->SetMarkerSize(1);
    hDEdxRafal[1]->SetMarkerStyle(29);
    hDEdxRafal[1]->Draw("SAME0");
    hDEdxRafal[2]->SetMarkerColor(4);
    hDEdxRafal[2]->SetMarkerSize(1);
    hDEdxRafal[2]->SetMarkerStyle(29);
    hDEdxRafal[2]->Draw("SAME0");
    legendPID = new TLegend(0.72,0.55,0.85,0.80,"","brNDC");
    legendPID -> AddEntry(hDEdxRafal[0], "#pi^{+} + #pi^{-}", "p");
    legendPID -> AddEntry(hDEdxRafal[1], "K^{+} + K^{-}", "p");
    legendPID -> AddEntry(hDEdxRafal[2], "p + #bar{p}", "p");
    legendPID -> Draw("same");

    cCanvas->Update();
    cCanvas->Write("dEdxByRafal");

    hDEdxDaniel[0]->SetStats(0);
    hDEdxDaniel[0]->SetTitle(" ; #frac{q}{e} #times p_{T} [GeV/c] ;dE/dx [keV/cm]");
    hDEdxDaniel[0]->SetMarkerColor(2);
    hDEdxDaniel[0]->SetMarkerSize(1);
    hDEdxDaniel[0]->SetMarkerStyle(29);
    hDEdxDaniel[0]->Draw("SCAT");
    hDEdxDaniel[1]->SetMarkerColor(3);
    hDEdxDaniel[1]->SetMarkerSize(1);
    hDEdxDaniel[1]->SetMarkerStyle(29);
    hDEdxDaniel[1]->Draw("SAME0");
    hDEdxDaniel[2]->SetMarkerColor(4);
    hDEdxDaniel[2]->SetMarkerSize(1);
    hDEdxDaniel[2]->SetMarkerStyle(29);
    hDEdxDaniel[2]->Draw("SAME0");
    legendPID = new TLegend(0.72,0.55,0.85,0.80,"","brNDC");
    legendPID -> AddEntry(hDEdxDaniel[0], "#pi^{+} + #pi^{-}", "p");
    legendPID -> AddEntry(hDEdxDaniel[1], "K^{+} + K^{-}", "p");
    legendPID -> AddEntry(hDEdxDaniel[2], "p + #bar{p}", "p");
    legendPID -> Draw("same");

    cCanvas->Update();
    cCanvas->Write("dEdxByDaniel");

    hDEdxTomas2[0]->SetStats(0);
    hDEdxTomas2[0]->SetTitle(" ; #frac{q}{e} #times p_{T} [GeV/c] ;dE/dx [keV/cm]");
    hDEdxTomas2[0]->SetMarkerColor(2);
    hDEdxTomas2[0]->SetMarkerSize(1);
    hDEdxTomas2[0]->SetMarkerStyle(29);
    hDEdxTomas2[0]->Draw("SCAT");
    hDEdxTomas2[1]->SetMarkerColor(3);
    hDEdxTomas2[1]->SetMarkerSize(1);
    hDEdxTomas2[1]->SetMarkerStyle(29);
    hDEdxTomas2[1]->Draw("SAME0");
    hDEdxTomas2[2]->SetMarkerColor(4);
    hDEdxTomas2[2]->SetMarkerSize(1);
    hDEdxTomas2[2]->SetMarkerStyle(29);
    hDEdxTomas2[2]->Draw("SAME0");
    legendPID = new TLegend(0.72,0.55,0.85,0.80,"","brNDC");
    legendPID -> AddEntry(hDEdxTomas2[0], "#pi^{+} + #pi^{-}", "p");
    legendPID -> AddEntry(hDEdxTomas2[1], "K^{+} + K^{-}", "p");
    legendPID -> AddEntry(hDEdxTomas2[2], "p + #bar{p}", "p");
    legendPID -> Draw("same");

    cCanvas->Update();
    cCanvas->Write("dEdxByTomas2");

    for (int i = 1; i < 100; ++i)
    {
        if(hPtCharged[0]->GetBinContent(i) != 0)    
            hRatio->SetBinContent(i,hPtCharged[1]->GetBinContent(i)/hPtCharged[0]->GetBinContent(i));
        else
            hRatio->SetBinContent(i,0);
    }

    hEtaDistCalib->SetTitle(" ; #eta ; Efficiency corrected number of tracks");

    hEtaDistCalib->GetXaxis()->SetTitleFont(42);
    hEtaDistCalib->GetYaxis()->SetTitleFont(42);
    hEtaDistCalib->GetXaxis()->SetLabelFont(42);
    hEtaDistCalib->GetYaxis()->SetLabelFont(42);
    hEtaDistCalib->GetXaxis()->SetTitleSize(0.045);
    hEtaDistCalib->GetYaxis()->SetTitleSize(0.045);
    hEtaDistCalib->GetXaxis()->SetTitleOffset(0.9);
    hEtaDistCalib->GetYaxis()->SetTitleOffset(1.3);
    hEtaDistCalib->GetXaxis()->SetRangeUser(-1.0,3.5);
    hEtaDistCalib->SetMarkerColor(4);
    hEtaDistCalib->SetMarkerSize(1);
    hEtaDistCalib->SetMarkerStyle(20);
    hEtaDistCalib->SetLineColor(4);
    hEtaDistCalib->SetLineStyle(1);
    hEtaDistCalib->SetLineWidth(1);
    hEtaDistCalib->Draw("E");

    TPaveText *textPub = new TPaveText(0.7,0.75,0.92,0.88,"brNDC");
    textPub -> SetTextSize(0.04);
    textPub -> SetTextAlign(22);
    textPub -> SetFillColor(0);
    textPub -> SetTextFont(42);
    textPub -> AddText("p + p #rightarrow p + X + p");
    textPub -> AddText("#sqrt{s} = 510 GeV");
    int NentriesEl = hEtaDistCalib->GetEntries();
    TString tileIdStrEl; tileIdStrEl.Form("%i h_{cand}^{Exc}",NentriesEl);
    textPub -> AddText(tileIdStrEl);
    textPub -> Draw("same");

    TPaveText *textSTAR;
    textSTAR = new TPaveText(0.75,0.89,0.9,0.95,"brNDC");
    textSTAR -> SetTextSize(0.04);
    textSTAR -> SetFillColor(0);
    textSTAR -> SetTextFont(62);
    textSTAR->AddText("STAR Internal");
    textSTAR -> Draw("same");

    TLegend* leg1 = new TLegend(0.6, 0.65, 0.8, 0.74);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.04);
    leg1->SetTextFont(42);
    leg1->AddEntry(hEtaDistCalib,"In+El (unlike-sign pairs)","p");
    leg1->Draw("same");

    cCanvas->Update();
    cCanvas->SaveAs("etaDistCalib.png");   
*/
    cCanvas->Close();
	//fout->Write();
	fout->Close();
    
}




void Make()
{
    double effTotal, effTPC, effTOF;
    unsigned int PID;
    for (int i = 0; i < 17; ++i)
    {
        if(trigger[i])
            triggerOn++;
        else
            triggerOff++;
    }
	
    hdEdx->Fill(log10(momentum[0]),log10(dEdx[0]));
    hdEdx->Fill(log10(momentum[1]),log10(dEdx[1]));

/*
    if(Eta[0] < -0.7 || Eta[0] > 0.7 || Eta[1] < -0.7 || Eta[1] > 0.7 )
        hEtaDif->Fill(Eta[0] - Eta[1]); 

    hPosition[East]->Fill(rpX[East], rpY[East]);
    hPosition[West]->Fill(rpX[West], rpY[West]);

	Double_t sqrtWest = sqrt(1 + (protonMass*protonMass)/(pRp[West]*pRp[West]));
    Double_t sqrtEast = sqrt(1 + (protonMass*protonMass)/(pRp[East]*pRp[East]));
    Double_t deltaTRP = (timeRp[East] - timeRp[West]);
   // cout<<"RP time: "<<timeRp[East]*speedOfLight<<"  "<< timeRp[West]*speedOfLight<< " "<<sqrtWest << " "<<rpZ[East]<< " "<<sqrtEast<<endl;
    Double_t z0 = ((deltaTRP*speedOfLight + rpZ[West]*sqrtWest + rpZ[East]*sqrtEast)/(sqrtWest + sqrtEast))*100; // covert from m to cm 
    Double_t z01 = ((deltaTRP*speedOfLight)/(2))*100; // covert from m to cm 

    hvertexZ->Fill(vertexesZ[0]);
    hvertexRP->Fill(z0 - z01);
*/
  //  cout<<vertexesZ[0] <<" "<< NhitsFit[0]<<" "<< NhitsFit[1]<<" "<< NhitsDEdx[0]<<" "<< NhitsDEdx[1]<<" "<< DcaZ[0]<<" "<<DcaZ[1] <<" "<<DcaXY[0] <<" "<<DcaXY[1] <<" "<<Eta[0] <<" "<<Eta[1]<<endl;
	//cout<<fourPiState<<endl;


    if(vertexesZ[0]<80 && vertexesZ[0] > -80 && NhitsFit[0] >=25 && NhitsFit[1] >= 25 && NhitsDEdx[0] >= 15 && NhitsDEdx[1] >= 15 && DcaZ[0] < 1 && DcaZ[0] > -1 && DcaZ[1] < 1 && DcaZ[1] > -1 && DcaXY[0] < 1.5 && DcaXY[1] < 1.5  && !fourPiState){ // && Eta[0] > -0.7 && Eta[0] < 0.7 && Eta[1] > -0.7 && Eta[1] < 0.7

        ++etaTotal;
        ++etaTotal;
        for(int i = 0; i < 20; ++i)
        {
            if(Eta[0] > (-0.6 - i*0.02) && Eta[0] < (0.6 + i*0.02) && Eta[1] > (-0.6 - i*0.02) && Eta[1] < (0.6 + i*0.02))
            {
                hEta[i]->Fill(Eta[0]);
                hEta[i]->Fill(Eta[1]);
            }
        }

        if(!(Eta[0] > -0.7 && Eta[0] < 0.7 && Eta[1] > -0.7 && Eta[1] < 0.7))
            return;


        for (int iTrack = 0; iTrack < 2; ++iTrack)
        {
            hEtaDist->Fill(Eta[iTrack]);
            if( Eta[iTrack] < - 1.0 || Eta[iTrack] > 1.0)
                continue;
            effTPC = hTPCeff[0]->GetBinContent( hTPCeff[0]->GetXaxis()->FindBin(vertexesZ[iTrack]), hTPCeff[0]->GetYaxis()->FindBin(transMomentum[iTrack]), hTPCeff[0]->GetZaxis()->FindBin(Eta[iTrack]));
            effTOF = hTOFeff[0]->GetBinContent( hTOFeff[0]->GetXaxis()->FindBin(vertexesZ[iTrack]), hTOFeff[0]->GetYaxis()->FindBin(transMomentum[iTrack]), hTOFeff[0]->GetZaxis()->FindBin(Eta[iTrack]));  
            effTotal = effTPC*effTOF;
            if(effTotal)
            {
                hEtaDistCalib->Fill(Eta[iTrack], 1/effTotal);
                if(elastic && nSigmaTPC[Pion][0] < 3 && nSigmaTPC[Pion][0] > -3 && nSigmaTPC[Pion][1] > -3 && nSigmaTPC[Pion][1] < 3 && (chiPair[0] < 3 || chiPair[1] > 3 || chiPair[2] < 3 || mSquared < 0.2 || mSquared > 0.32) && (chiPair[0] < 3 || chiPair[1] < 3 || chiPair[2] > 3 || mSquared < 0.7 || mSquared > 1.1))
                    pTGoldenPions->Fill(transMomentum[iTrack], 1/effTotal);
            }
        }

/*
		hDeltaVertex->Fill(vertexesZ[0] - z0);
		hvertexRPGolden->Fill(z0);

        if(charge[0] == 1)
            hPtCharged[1]->Fill(momentum[0]);
        else
            hPtCharged[0]->Fill(momentum[0]);
        if(charge[1] == 1)
            hPtCharged[1]->Fill(momentum[1]);
        else
            hPtCharged[0]->Fill(momentum[1]);
        
//        for (int i = 0; i < 16; ++i)
//        {
//            if(ADC[i/2][i%2] != 0 && TAC[i/2][i%2] != 0)
//                hSlewing[i]->Fill(ADC[i/2][i%2], TAC[i/2][i%2]);
//        }

		if(TOFtime[1] < 0 || TOFlength[0] < 0 || TOFlength[1] < 0 || TOFtime[0] < 0){
			hRunNumber[1]->Fill(runNumber);
			if(TOFtime[1] == TOFlength[1] && TOFtime[1] == -999)
			{
				nTOFBadTracks++;
				hTime[1]->Fill(TOFtime[1]);
				hTofLength[1]->Fill(TOFlength[1]);
			}else
			{
				nTOFGoodTracks++;
				hDeltaTime->Fill(TOFtime[0] - TOFtime[1]);
				hTime[0]->Fill(TOFtime[1]);
				hTofLength[0]->Fill(TOFlength[1]);
			}
			if(TOFtime[0] == TOFlength[0] && TOFtime[0] == -999)
			{
				nTOFBadTracks++;
				hTime[1]->Fill(TOFtime[0]);
				hTofLength[1]->Fill(TOFlength[0]);
			}else
			{
				nTOFGoodTracks++;
				hDeltaTime->Fill(TOFtime[0] - TOFtime[1]);
				hTime[0]->Fill(TOFtime[0]);
				hTofLength[0]->Fill(TOFlength[0]);
			}

			if(TOFtime[0]*TOFlength[0] < 0 || TOFtime[1]*TOFlength[1] < 0)
				nWierdTOFInfo++;
		}else
		{
			nTOFGoodTracks++;
			nTOFGoodTracks++;
			hDeltaTime->Fill(TOFtime[0] - TOFtime[1]);
			hTime[0]->Fill(TOFtime[0]);
			hTime[0]->Fill(TOFtime[1]);
			hTofLength[0]->Fill(TOFlength[0]);
			hTofLength[0]->Fill(TOFlength[1]);
			hRunNumber[0]->Fill(runNumber);

            if(deltaDeltaTOF[Proton] > -0.5 && deltaDeltaTOF[Proton] < 0.5 )
            {
                hchiPairPion[Proton]->Fill(chiPair[Pion]);
                hchiPairKaon[Proton]->Fill(chiPair[Kaon]);
                hchiPairProton[Proton]->Fill(chiPair[Proton]);
                hMSquered[Proton]->Fill(mSquared); 
            }else if(deltaDeltaTOF[Kaon] > -0.5 && deltaDeltaTOF[Kaon] < 0.5 )
            {
                hchiPairPion[Kaon]->Fill(chiPair[Pion]);
                hchiPairKaon[Kaon]->Fill(chiPair[Kaon]);
                hchiPairProton[Kaon]->Fill(chiPair[Proton]);
                hMSquered[Kaon]->Fill(mSquared);   
            }else if(deltaDeltaTOF[Pion] > -0.5 && deltaDeltaTOF[Pion] < 0.5 )
            {
                hchiPairPion[Pion]->Fill(chiPair[Pion]);
                hchiPairKaon[Pion]->Fill(chiPair[Kaon]);
                hchiPairProton[Pion]->Fill(chiPair[Proton]);
                hMSquered[Pion]->Fill(mSquared);
            }

            effTotal = 1;
            PID = 0;
            if(chiPair[Pion] > 3 && chiPair[Kaon] > 3 && chiPair[Proton] < 3 && mSquared > 0.7 && mSquared < 1.1)
            {
                for (int iTrack = 0; iTrack < 2; ++iTrack)
                {
                    hDEdxRafal[Proton]->Fill(momentum[iTrack]*charge[iTrack],dEdx[iTrack]); 
                    if(charge[iTrack] > 0)
                        PID = 3;
                    effTPC = hTPCeff[2 + PID]->GetBinContent( hTPCeff[2 + PID]->GetXaxis()->FindBin(vertexesZ[iTrack]), hTPCeff[2 + PID]->GetYaxis()->FindBin(transMomentum[iTrack]), hTPCeff[2 + PID]->GetZaxis()->FindBin(Eta[iTrack]));
                    effTOF = hTOFeff[2 + PID]->GetBinContent( hTOFeff[2 + PID]->GetXaxis()->FindBin(vertexesZ[iTrack]), hTOFeff[2 + PID]->GetYaxis()->FindBin(transMomentum[iTrack]), hTOFeff[2 + PID]->GetZaxis()->FindBin(Eta[iTrack]));  
                    effTotal = effTotal*effTPC*effTOF;
                }
                if(effTotal)
                    hInvMass[Proton]->Fill(invMass[Proton], 1/effTotal);
            }
            else if(chiPair[Pion] > 3 && chiPair[Kaon] < 3 && chiPair[Proton] > 3 && mSquared > 0.2 && mSquared < 0.32)
            {
                for (int iTrack = 0; iTrack < 2; ++iTrack)
                {
                    hDEdxRafal[Kaon]->Fill(momentum[iTrack]*charge[iTrack],dEdx[iTrack]);
                    if(charge[iTrack] > 0)
                        PID = 3;
                    effTPC = hTPCeff[1 + PID]->GetBinContent( hTPCeff[1 + PID]->GetXaxis()->FindBin(vertexesZ[iTrack]), hTPCeff[1 + PID]->GetYaxis()->FindBin(transMomentum[iTrack]), hTPCeff[1 + PID]->GetZaxis()->FindBin(Eta[iTrack]));
                    effTOF = hTOFeff[1 + PID]->GetBinContent( hTOFeff[1 + PID]->GetXaxis()->FindBin(vertexesZ[iTrack]), hTOFeff[1 + PID]->GetYaxis()->FindBin(transMomentum[iTrack]), hTOFeff[1 + PID]->GetZaxis()->FindBin(Eta[iTrack]));  
                    effTotal = effTotal*effTPC*effTOF;
                }
                if(effTotal)
                    hInvMass[Kaon]->Fill(invMass[Kaon], 1/effTotal);
            }
            else if( nSigmaTPC[Pion][0] < 3 && nSigmaTPC[Pion][0] > -3 && nSigmaTPC[Pion][1] > -3 && nSigmaTPC[Pion][1] < 3)
            {
                for (int iTrack = 0; iTrack < 2; ++iTrack)
                {
                    hDEdxRafal[Pion]->Fill(momentum[iTrack]*charge[iTrack],dEdx[iTrack]);
                    if(charge[iTrack] > 0)
                        PID = 3;
                    effTPC = hTPCeff[0 + PID]->GetBinContent( hTPCeff[0 + PID]->GetXaxis()->FindBin(vertexesZ[iTrack]), hTPCeff[0 + PID]->GetYaxis()->FindBin(transMomentum[iTrack]), hTPCeff[0 + PID]->GetZaxis()->FindBin(Eta[iTrack]));
                    effTOF = hTOFeff[0 + PID]->GetBinContent( hTOFeff[0 + PID]->GetXaxis()->FindBin(vertexesZ[iTrack]), hTOFeff[0 + PID]->GetYaxis()->FindBin(transMomentum[iTrack]), hTOFeff[0 + PID]->GetZaxis()->FindBin(Eta[iTrack])); 
                    effTotal = effTotal*effTPC*effTOF;
                }
                if(effTotal)
                    hInvMass[Pion]->Fill(invMass[Pion], 1/effTotal);
            }

            if(chiPair[Proton] < 3 && mSquared > 0.7 && mSquared < 1.1)
                for (int iTrack = 0; iTrack < 2; ++iTrack)
                    hDEdxTomas[Proton]->Fill(momentum[iTrack]*charge[iTrack],dEdx[iTrack]);
            else if(chiPair[Kaon] < 3 && mSquared > 0.2 && mSquared < 0.32)
                for (int iTrack = 0; iTrack < 2; ++iTrack)
                    hDEdxTomas[Kaon]->Fill(momentum[iTrack]*charge[iTrack],dEdx[iTrack]);
            else if( nSigmaTPC[Pion][0] < 3 && nSigmaTPC[Pion][0] > -3 && nSigmaTPC[Pion][1] > -3 && nSigmaTPC[Pion][1] < 3)
                for (int iTrack = 0; iTrack < 2; ++iTrack)
                    hDEdxTomas[Pion]->Fill(momentum[iTrack]*charge[iTrack],dEdx[iTrack]);

            if(chiPair[Proton] < 3 && ((mSquared > 0.7 && mSquared < 1.1) || (deltaDeltaTOF[Proton] > -0.5 && deltaDeltaTOF[Proton] < 0.5)))
                for (int iTrack = 0; iTrack < 2; ++iTrack)
                    hDEdxTomas2[Proton]->Fill(momentum[iTrack]*charge[iTrack],dEdx[iTrack]);
            else if(chiPair[Kaon] < 3 && ((mSquared > 0.2 && mSquared < 0.32) || (deltaDeltaTOF[Kaon] > -0.5 && deltaDeltaTOF[Kaon] < 0.5)))
                for (int iTrack = 0; iTrack < 2; ++iTrack)
                    hDEdxTomas2[Kaon]->Fill(momentum[iTrack]*charge[iTrack],dEdx[iTrack]);
            else if( nSigmaTPC[Pion][0] < 3 && nSigmaTPC[Pion][0] > -3 && nSigmaTPC[Pion][1] > -3 && nSigmaTPC[Pion][1] < 3)
                for (int iTrack = 0; iTrack < 2; ++iTrack)
                    hDEdxTomas2[Pion]->Fill(momentum[iTrack]*charge[iTrack],dEdx[iTrack]);

            if(chiPair[Proton] < 3 && deltaDeltaTOF[Proton] > -0.5 && deltaDeltaTOF[Proton] < 0.5)
                for (int iTrack = 0; iTrack < 2; ++iTrack)
                    hDEdxDaniel[Proton]->Fill(momentum[iTrack]*charge[iTrack],dEdx[iTrack]);
            else if(chiPair[Kaon] < 3 && deltaDeltaTOF[Kaon] > -0.5 && deltaDeltaTOF[Kaon] < 0.5)
                for (int iTrack = 0; iTrack < 2; ++iTrack)
                    hDEdxDaniel[Kaon]->Fill(momentum[iTrack]*charge[iTrack],dEdx[iTrack]);
            else if( nSigmaTPC[Pion][0] < 3 && nSigmaTPC[Pion][0] > -3 && nSigmaTPC[Pion][1] > -3 && nSigmaTPC[Pion][1] < 3)
                for (int iTrack = 0; iTrack < 2; ++iTrack)
                    hDEdxDaniel[Pion]->Fill(momentum[iTrack]*charge[iTrack],dEdx[iTrack]);
        }
        */
	}


}


void Init(){

    hdEdx = new TH2D("dEdx", "dEdx", 100,-1.5,2,100,0,2.5);
    hEtaDist = new TH1D("etaDist", "etaDist", 100, -2,2);
    hEtaDist->SetTitle(" ; #eta ; ");
    hEtaDistCalib = new TH1D("etaDistCalib", "etaDistCalib", 100, -2,2);
    hEtaDistCalib->SetTitle(" ; #eta ; ");

    pTGoldenPions = new TH1D("pTGoldenPi", "pTGoldenPi", 200, 0,4);

    hLFactor =  new TH1D("LFactor", "Loosing factor",20,-0.5,19.5);
    for(int i = 0; i < 20; ++i)
        hEta[i] = new TH1D(Form("etaDist_%d",i), "etaDist", 100, -2,2);

    hEtaDif = new TH1D("etaDif", "etaDif", 100, -2,2);
	hvertexZ = new TH1D("vertexZ", "vertexZ", 200, -200, 200);
	hvertexRP = new TH1D("vertexRP","vertexRP", 200, -200, 200);
	hvertexRPGolden = new TH1D("vertexRPGolden", "vertexRP for golden events", 200, -200, 200);
	hDeltaVertex = new TH1D("DeltaVertex", "vertexZ - vertexRP for golden events", 200, -200, 200);

    hPosition[East] = new TH2D("PositionE", "Position East", 100,-0.1,0.1,100,-0.1,0.1);
    hPosition[West] = new TH2D("PositionW", "Position West", 100,-0.1,0.1,100,-0.1,0.1);

    hPtCharged[0] = new TH1D("ptPositive", "Pt of negative primary tracks", 100, 0 ,10);
    hPtCharged[1] = new TH1D("ptNegative", "Pt of positive primary tracks", 100, 0 ,10);

    hRatio = new TH1D("ratio", "Ratio of positive and negative primary tracks", 100, 0 ,10);
    for (int i = 0; i < 16; ++i)
        hSlewing[i] = new TH2D("Slewing_" + rpNames[i/2] + Form("_%i",i%2), "Slewing_" + rpNames[i/2] + Form("_%i",i%2), 100,0,800,200,0,1600);


    for (int i = 0; i < nParticles; ++i)
    {
        hDEdxRafal[i] = new TH2D("dEdxRafal" + particleLables[i], "dEdx" + particleLables[i], 120,-3,3,80,0,10);
        hDEdxTomas[i] = new TH2D("dEdxTomas" + particleLables[i], "dEdx" + particleLables[i], 120,-3,3,80,0,10);
        hDEdxDaniel[i] = new TH2D("dEdxDaniel" + particleLables[i], "dEdx" + particleLables[i], 120,-3,3,80,0,10);
        hDEdxTomas2[i] = new TH2D("dEdxTomas2" + particleLables[i], "dEdx" + particleLables[i], 120,-3,3,80,0,10);
        hchiPairPion[i]  = new TH1D("chiPairPionFor"  + particleLables[i], "chiPairPion For" + particleLables[i], 100, 0, 35);
        hchiPairKaon[i]  = new TH1D("chiPairKaonFor" + particleLables[i], "chiPairKaon For" + particleLables[i], 100, 0, 35);
        hchiPairProton[i]  = new TH1D("chiPairProtonFor" + particleLables[i], "chiPairProton For" + particleLables[i], 100, 0, 35);
        hMSquered[i]  = new TH1D("mSquaredFor" + particleLables[i], "mSquared for " + particleLables[i] , 200, -0.5, 2.5);
        hInvMass[i]  = new TH1D("invMass" + particleLables[i], "inv. mass " + particleLables[i] , 100, 0, 3.5);
        hMomentum[i] = new TH2D("momentumCorr" + particleLables[i], "momentumCorr" + particleLables[i], 100,-3,3,100,-3,3);
        hTransMomentum[i] = new TH2D("transMomentumCorr" + particleLables[i], "transMomentumCorr" + particleLables[i], 100,-3,3,100,-3,3);
    }

	fout->mkdir("GoodTOFinfo")->cd();
	hDeltaTime = new TH1D("DeltaTime", "DeltaTime", 200, -50, 50);

	for (int i = 0; i < 2; ++i)
	{
		if(i==1)
			fout->mkdir("BadTOFinfo")->cd();
		hTime[i] = new TH1D("Time", "Time", 2000, -1000, 50000);
		hTofLength[i] = new TH1D("TOF length", "TOF length", 200, -1010.0, 500.0);
		hRunNumber[i] = new TH1D("runNumber", "run number", 71025, 18078045, 18149069);
	}


}

void ConnectInput(){
	recTree = dynamic_cast<TTree*>( data->Get("recTree") );
	if (!recTree)
	{
		cout<<"Error: cannot get recTree"<<endl;
		return;
	}

// PID and some quality event info
    recTree->SetBranchAddress("missingPt", &missingPt);
    recTree->SetBranchAddress("deltaTOF", &deltaTOF);
    recTree->SetBranchAddress("mSquared", &mSquared); 
    recTree->SetBranchAddress("nSigTrk1Pion", &nSigmaTPC[Pion][0]);
    recTree->SetBranchAddress("nSigTrk2Pion", &nSigmaTPC[Pion][1]);
    for (int iPart = 0; iPart < nParticles; ++iPart)
    {
        recTree->SetBranchAddress("invMass" + particleLables[iPart], &invMass[iPart]);
        recTree->SetBranchAddress("chiPair" + particleLables[iPart], &chiPair[iPart]);
        recTree->SetBranchAddress("deltaTOFExpected" + particleLables[iPart], &deltaTOFExpected[iPart]);
        recTree->SetBranchAddress("deltaDeltaTOF" + particleLables[iPart], &deltaDeltaTOF[iPart]);  
    }


// Vertex info
    recTree->SetBranchAddress("vertexZ", &vertexesZ[0]);

// Central track info
    for (int i = 0; i < 4; ++i)
    {
        recTree->SetBranchAddress(Form("dEdx%i",i), &dEdx[i]);
        recTree->SetBranchAddress(Form("momentum%i",i), &momentum[i]);
        recTree->SetBranchAddress(Form("transMomentum%i",i), &transMomentum[i]);
        recTree->SetBranchAddress(Form("charge%i",i), &charge[i]);
        recTree->SetBranchAddress(Form("TOFtime%i",i), &TOFtime[i]);
        recTree->SetBranchAddress(Form("TOFlength%i",i), &TOFlength[i]);
        recTree->SetBranchAddress(Form("DcaXY%i",i), &DcaXY[i]);
        recTree->SetBranchAddress(Form("DcaZ%i",i), &DcaZ[i]);
        recTree->SetBranchAddress(Form("NhitsFit%i",i), &NhitsFit[i]);
        recTree->SetBranchAddress(Form("NhitsDEdx%i",i), &NhitsDEdx[i]);
        recTree->SetBranchAddress(Form("Eta%i",i), &Eta[i]);
        recTree->SetBranchAddress(Form("Phi%i",i), &Phi[i]);
        recTree->SetBranchAddress(Form("Chi2%i",i), &Chi2[i]);
    }
// RP track info  
    for (int i = 0; i < nSides; ++i)
    {
        recTree->SetBranchAddress("rpX" + sideLabel[i], &rpX[i]);
        recTree->SetBranchAddress("rpY" + sideLabel[i], &rpY[i]);
        recTree->SetBranchAddress("rpZ" + sideLabel[i], &rpZ[i]);
        recTree->SetBranchAddress("thetaRp" + sideLabel[i], &thetaRp[i]);
        recTree->SetBranchAddress("phiRp" + sideLabel[i], &phiRp[i]);
        recTree->SetBranchAddress("timeRp" + sideLabel[i], &timeRp[i]);
        recTree->SetBranchAddress("pRp" + sideLabel[i], &pRp[i]);
        recTree->SetBranchAddress("ptRp" + sideLabel[i], &ptRp[i]);
        recTree->SetBranchAddress("etaRp" + sideLabel[i], &etaRp[i]);
        recTree->SetBranchAddress("xCorrelationsRp" + sideLabel[i], &xCorrelationsRp[i]);
        recTree->SetBranchAddress("yCorrelationsRp" + sideLabel[i], &yCorrelationsRp[i]);
        recTree->SetBranchAddress("t" + sideLabel[i], &t[i]);
        recTree->SetBranchAddress("xi" + sideLabel[i], &xi[i]);
    }
// RP event info
    for (int i = 0; i < nRomanPots; ++i)
    {
        recTree->SetBranchAddress("ADC_" + rpNames[i] + "V", &ADC[i][0]);
        recTree->SetBranchAddress("ADC_" + rpNames[i] + "H", &ADC[i][1]);
        recTree->SetBranchAddress("TAC_" + rpNames[i] + "V", &TAC[i][0]);
        recTree->SetBranchAddress("TAC_" + rpNames[i] + "H", &TAC[i][1]);
    }

// event info
    recTree->SetBranchAddress("elastic", &elastic);
    recTree->SetBranchAddress("fourPiState", &fourPiState);
    recTree->SetBranchAddress("runNumber", &runNumber);
    recTree->SetBranchAddress("VPDTimeDiff", &VPDTimeDiff);
    recTree->SetBranchAddress("VPDSumWest", &VPDSumWest);
    recTree->SetBranchAddress("VPDSumEast", &VPDSumEast);
    for (int i = 0; i < nSides; ++i)
    {
        recTree->SetBranchAddress("BBCSmall" + sideLabel[i], &BBCSmall[i]);
        recTree->SetBranchAddress("BBCLarge" + sideLabel[i], &BBCLarge[i]); 
    }
    recTree->SetBranchAddress("RP_CPT2_570701", &trigger[3]);
    recTree->SetBranchAddress("RP_CPT2noBBCL_570705", &trigger[7]);
    recTree->SetBranchAddress("RP_CPT2_570711", &trigger[9]);
    recTree->SetBranchAddress("RP_CPT2_590701", &trigger[12]);
    recTree->SetBranchAddress("RP_CPT2noBBCL_590705", &trigger[14]);
    recTree->SetBranchAddress("RP_CPTnoBBCL_590708", &trigger[15]);
}


