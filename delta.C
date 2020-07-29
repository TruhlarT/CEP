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


void delta()
{
	TString input = "/home/truhlar/Desktop/STAR/CEP/Analysis/Data/NewEff/effPionsP.root";
	data = TFile::Open(input, "read");
	if (!data)
	{
		cout<<"Error: cannot open "<<input<<endl;
		return;
	}


    TCanvas *cCanvas = new TCanvas("cCanvas","cCanvas",800,700);
    gPad->SetMargin(0.1,0.04,0.1,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top) 
    gStyle->SetOptStat("");
    gStyle->SetPalette(1);
    gStyle->SetLineWidth(2);      //axis line
    gStyle->SetFrameLineWidth(2); //frame line


    TH1F* hDelta = (TH1F*)data -> Get("delta");
    gPad->SetLogy();
    hDelta->SetTitle(" ; #delta^{2} ;Number of tracks");
    hDelta->GetXaxis()->SetTitleFont(42);
    hDelta->GetYaxis()->SetTitleFont(42);
    hDelta->GetXaxis()->SetLabelFont(42);
    hDelta->GetYaxis()->SetLabelFont(42);
    hDelta->GetXaxis()->SetTitleSize(0.045);
    hDelta->GetYaxis()->SetTitleSize(0.045);
    hDelta->GetXaxis()->SetTitleOffset(0.9);
    hDelta->GetYaxis()->SetTitleOffset(1.0);
    hDelta->GetXaxis()->SetRangeUser(0.0,0.05);
    hDelta->SetMarkerColor(4);
    hDelta->SetMarkerSize(1);
    hDelta->SetMarkerStyle(20);
    hDelta->SetLineColor(4);
    hDelta->SetLineStyle(1);
    hDelta->SetLineWidth(1);
    hDelta->Draw("E");

    TLine *left02 = new TLine(0.001,0.0,0.001,hDelta->GetMaximum()/2);
    left02->SetLineStyle(10);
    left02->SetLineColor(1);
    left02->SetLineWidth(4);
    left02->Draw("same");


    TLegend* leg1 = new TLegend(0.6, 0.65, 0.78, 0.74);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.04);
    leg1->SetTextFont(42);
    leg1->AddEntry(hDelta,"STARsim","p");
    leg1->Draw("same");

    cCanvas->Update();



    //cCanvas->Close();
    
}



