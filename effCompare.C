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

TFile* fout;
TFile* TPCeff;
TFile* MyEff;

TH1F* hMyEff;
TH1F* hRafalEff;
TH1F* hComparison;

TH3F* hTPCeff;
TH3F* hMyEff3D;

void Init();
void ConnectInput();
void Make();

void effCompare()
{
    TString rafalEff = "/home/truhlar/Desktop/STAR/CEP/Analysis/Data/etaPhiEfficiency_16_01_19_delta015_twoRuns.root";
    TString myEff = "/home/truhlar/Desktop/STAR/CEP/Analysis/Data/myEffNew.root";
	TString output = "/home/truhlar/Desktop/STAR/CEP/Analysis/Outputs/effCompare.root";

    TPCeff = TFile::Open(rafalEff, "read");
    if (!TPCeff)
    {
        cout<<"Error: cannot open "<<rafalEff<<endl;
        return;
    }

    MyEff = TFile::Open(myEff, "read");
    if (!MyEff)
    {
        cout<<"Error: cannot open "<<myEff<<endl;
        return;
    }

    hTPCeff = (TH3F*)TPCeff -> Get(Form("hTPCEffiCD%i120",3));
/*
    cout<<"Scanning X-axis"<<endl;
    TH1* hX = hTPCeff->Project3D("x");
    for (int i = 1; i <= hX->GetNbinsX(); ++i)
        cout<<"Bin "<<i<<" start at "<< hX->GetBinLowEdge(i)<<" w Width: "<< hX->GetBinWidth(i)<<endl;
    
    cout<<"Scanning Y-axis"<<endl;
    TH1* hY = hTPCeff->Project3D("y");
    for (int i = 1; i <= hY->GetNbinsX(); ++i)
        cout<<"Bin "<<i<<" start at "<< hY->GetBinLowEdge(i)<<" w Width: "<< hY->GetBinWidth(i)<<endl;
    
    cout<<"Scanning Z-axis"<<endl;
    TH1* hZ = hTPCeff->Project3D("z");
    for (int i = 1; i <= hZ->GetNbinsX(); ++i)
        cout<<"Bin "<<i<<" start at "<< hZ->GetBinLowEdge(i)<<" w Width: "<< hZ->GetBinWidth(i)<<endl;
    return;
*/

    hMyEff3D = (TH3F*)MyEff -> Get("effRafal");


	fout = new TFile(output,"RECREATE");
	Init(); // Preparing histograms 
    Make();

	fout->Write();
	fout->Close();
    
}




void Make()
{
    Double_t ratio, rafal, mine;
    Int_t x, y, z;
    Int_t xMax, yMax, zMax;
    xMax = 20;
    yMax = 20;
    zMax = 24;
    x = 1;
    y = 1;
    z = 1;
    for (int i = 1; i <= xMax*yMax*zMax; ++i)
    {
        rafal = hTPCeff->GetBinContent( x, y, z);
        mine = hMyEff3D->GetBinContent( x, y, z);

        hMyEff->SetBinContent(i, mine);
        hRafalEff->SetBinContent(i, rafal);
        
        if(mine)
        {
            ratio = rafal / mine;
            hComparison->SetBinContent(i, ratio);
        }

        x++;
        if(x == xMax+1)
        {
            x = 1;
            y++;
            if(y == yMax+1)
            {
                y = 1;
                z++;
            }
        }    
    }

}


void Init(){

    hMyEff = new TH1F("myEff", "my efficiencies", 9600, -0.5, 9599.5);
    hRafalEff = new TH1F("rafalEff", "Rafal efficiencies", 9600, -0.5, 9599.5);
    hComparison = new TH1F("comparison", "Rafal / mine", 9600, -0.5, 9599.5);

}

