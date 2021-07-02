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
TTree* bcgTree;

TH1D* hvertexZ;

/////////////////////////////////

Bool_t elastic, fourPiState;
Bool_t trigger[17];

Double_t chiPair[nParticles]; 
Double_t invMass[nParticles];
Double_t missingPt, deltaTOF, mSquared;

Double_t t[nSides];
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


void Init();
void ConnectInput();
void Make();

void zvertex()
{
	TString input = "/home/truhlar/Desktop/STAR/CEP/Analysis/Data/anaFlow2.root";
	data = TFile::Open(input, "read");
	if (!data)
	{
		cout<<"Error: cannot open "<<input<<endl;
		return;
	}

    recTree = dynamic_cast<TTree*>( data->Get("recTree") );
    bcgTree = dynamic_cast<TTree*>( data->Get("Background") );
    if (!recTree || !bcgTree)
    {
        cout<<"Error: cannot get recTree or bcgTree"<<endl;
        return;
    }

	Init(); // Preparing histograms 
	ConnectInput(); // Connecting input

    TCanvas *cCanvas = new TCanvas("cCanvas","cCanvas",800,700);
    gPad->SetMargin(0.12,0.02,0.1,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top) 
    //gStyle->SetOptStat("");

    TH1F* histSignal;
    TH1F* histBackground;

    TString variable = "vertexZ";
    TString cuts = "!fourPiState && NhitsFit1 >=25 && NhitsFit0 >= 25 && NhitsDEdx1 >= 15 && NhitsDEdx0 >= 15 && DcaZ1 < 1 && DcaZ1 > -1 && DcaZ0 < 1 && DcaZ0 > -1 && DcaXY1 < 1.5 && DcaXY0 < 1.5 && Eta1 > -0.7 && Eta1 < 0.7 && Eta0 > -0.7 && Eta0 < 0.7 &&  tEast < -0.12 && tWest < -0.12 && tEast > -1.0  && tWest > -1.0";
    TString usedCuts = cuts;
    int nBins = 104;
    double min = -220.0;
    double max = 240.0;
    bcgTree->Draw(variable+">>" + variable +"Bcg1(" + nBins + "," + min + "," + max + ")",usedCuts);
    histBackground = (TH1F*)gPad->GetPrimitive(variable +"Bcg1");
    histBackground->SetMarkerColor(2);
    histBackground->SetMarkerSize(1);
    histBackground->SetMarkerStyle(20);


    recTree->Draw(variable + ">>" + variable +"Sig1(" + nBins + "," + min + "," + max + ")",usedCuts);
    histSignal = (TH1F*)gPad->GetPrimitive(variable +"Sig1");   
    histSignal->SetTitle(" ; Z_{vrtx} [cm]; Number of events");
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
    textSTAR->AddText("THIS THESIS");
    textSTAR -> Draw("same");

    TLegend* leg1 = new TLegend(0.6, 0.65, 0.78, 0.74);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.04);
    leg1->SetTextFont(42);
    leg1->AddEntry(histSignal,"El+Inel (unlike-sign pairs)","p");
    leg1->AddEntry(histBackground,"El+Inel (like-sign pairs)","p");
    leg1->Draw("same");


    cCanvas->Update();
    //cCanvas->Write(variable); 
    //cCanvas->SaveAs(path + name[i] + ".pdf");
    TF1 * f1 = new TF1("f1","gaus");
    gStyle->SetOptStat("1111");
    gStyle->SetOptFit(1);    
    histSignal -> Add(histBackground, -1);
    histSignal -> Fit("f1");
    cout<<f1->GetParameter(2)<<endl;
    histSignal->Draw("E");
    f1->Draw("same");
    textPub = new TPaveText(0.7,0.72,0.92,0.88,"brNDC");
    textPub -> SetTextSize(0.04);
    textPub -> SetTextAlign(22);
    textPub -> SetFillColor(0);
    textPub -> SetTextFont(42);
    textPub -> AddText("p + p #rightarrow p + X + p");
    textPub -> AddText("#sqrt{s} = 510 GeV");
    textPub -> AddText("Like-sign pairs subtracted");
    textPub -> Draw("same");

    textSTAR -> Draw("same");

    gPad->SetTickx();
    gPad->SetTicky(); 
    gStyle->SetLineWidth(2);      //axis line
    gStyle->SetFrameLineWidth(2); //frame line


    leg1 = new TLegend(0.7, 0.65, 0.78, 0.74);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.04);
    leg1->SetTextFont(42);
    leg1->AddEntry(histSignal,"Data","p");
    leg1->AddEntry(f1,"Gaussian fit","l");
    leg1->Draw("same");
    cCanvas->Update();








    //cCanvas->Close();
    
}




void Make()
{

	
}


void Init(){

	hvertexZ = new TH1D("vertexZ", "vertexZ", 200, -200, 200);

}

void ConnectInput(){

// PID and some quality event info
    recTree->SetBranchAddress("missingPt", &missingPt);
    recTree->SetBranchAddress("deltaTOF", &deltaTOF);
    recTree->SetBranchAddress("mSquared", &mSquared); 
    recTree->SetBranchAddress("nSigTrk1Pion", &nSigmaTPC[Pion][0]);
    recTree->SetBranchAddress("nSigTrk2Pion", &nSigmaTPC[Pion][1]);
    for (int iPart = 0; iPart < nParticles; ++iPart)
    {
        recTree->SetBranchAddress("invMass" + particleLables[iPart], &invMass[iPart]);
        recTree->SetBranchAddress("chiPair" + particleLables[iPart], &chiPair[iPart]);}


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


// event info
    recTree->SetBranchAddress("elastic", &elastic);
    recTree->SetBranchAddress("fourPiState", &fourPiState);

    recTree->SetBranchAddress("RP_CPT2_570701", &trigger[3]);
    recTree->SetBranchAddress("RP_CPT2noBBCL_570705", &trigger[7]);
    recTree->SetBranchAddress("RP_CPT2_570711", &trigger[9]);
    recTree->SetBranchAddress("RP_CPT2_590701", &trigger[12]);
    recTree->SetBranchAddress("RP_CPT2noBBCL_590705", &trigger[14]);
    recTree->SetBranchAddress("RP_CPTnoBBCL_590708", &trigger[15]);
}


