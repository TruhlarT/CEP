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

enum SIDE {E = 0, East = 0, W = 1, West = 1, nSides};
enum RP_ID {E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D, nRomanPots};
enum BRANCH_ID { EU, ED, WU, WD, nBranches };
enum ARM_ID { EU_WU, ED_WD, nArms };
enum STATION_ID { E1, E2, W1, W2, nStations };
enum RP_CONFIGURATION {EUD, EDU, IUU, IDD, nConfiguration}; 

const int nTriggers = 17;
const int triggerID[] = { 570209, 570219, 570229, 570701, 570702, 570703, 570704, 570705, 
                  570709, 570711, 570712, 570719, 590701, 590703, 590705, 590708, 590709};
// 570702 RP_UPC // 570712 RP_UPC // 570703 RP_SDT // 570709 RP_ET // 570719 RP_ET // 570701 RP_CPT2 // 570711 RP_CPT2 // 570705 RP_CPT2noBBCL // 570704 RP_Zerobias // 590703 RP_SDT // 590709 RP_ET // 590701 RP_CPT2 // 590705 RP_CPT2noBBCL // 590708 RP_CPTnoBBCL // 570209 JPsi*HTTP // 570219 JPsi*HTTP // 570229 JPsi*HTTP

TString rpNames[nRomanPots] = { TString("E1U"), TString("E1D"), TString("E2U"), TString("E2D"), TString("W1U"), TString("W1D"), TString("W2U"), TString("W2D")};

TFile* data;
TFile* fout;

TH2D* hPosition;
TH2D* hSlewing[12][2]; // 12 + 2xPMTs
// East - Up - 3x x-region, East - Down - 3x x-region
// West - Up - 3x x-region, West - Down - 3x x-region 

TTree *recTree[nRomanPots];


Double_t ADC[2]; 
Double_t TAC[2];
Double_t rpX, rpY, rpZ; 

Double_t y[nSides][4] = 
{
    { -0.066, -0.024, 0.024, 0.062},
    { -0.057, -0.024, 0.026, 0.072}
};
Double_t x[nSides][4] = 
{
    { -0.02, 0.0, 0.02, 0.045},
    { -0.016, 0.0, 0.02, 0.042}
};

Int_t rpID;

void Init();
void ConnectInput();
void Make();

void slewing()
{
	TString input = "/home/truhlar/Desktop/STAR/CEP/Analysis/Data/SlewRHICfNew.root";
	TString output = "/home/truhlar/Desktop/STAR/CEP/Analysis/Outputs/Slew.root";

	data = TFile::Open(input, "read");
	if (!data)
	{
		cout<<"Error: cannot open "<<input<<endl;
		return;
	}


	fout = new TFile(output,"RECREATE");
	Init(); // Preparing histograms 
	ConnectInput(); // Connecting input


	////////////// Making histograms
	for (int i = 0; i < nRomanPots; ++i)
	{ //get the event   
        rpID = i; 
    	Long64_t nev = recTree[i]->GetEntries();
    	cout<<"Proccesing RP"<<i<<" with "<<nev<<" events"<<endl;
        for(Long64_t iev=0; iev<nev; ++iev) 
        {
            recTree[i]->GetEntry(iev); 
            Make();
        }	   
    } 
    fout->cd();
	fout->Write();
	fout->Close();
    
}




void Make()
{
    if(rpX != -999 && rpY != -999)
        hPosition->Fill(rpX,rpY);

    if(rpY < y[rpID/4][2 - (rpID%2)*2] || rpY > y[rpID/4][3 - (rpID%2)*2])
    {
        return;
    }


    for (int i = 0; i < 3; ++i)
    {
        if(rpX > x[rpID/4][i] && rpX < x[rpID/4][i+1])
        {
            hSlewing[(rpID/4)*6 + (rpID%2)*3 + i][0]->Fill(ADC[0],TAC[0]);
            hSlewing[(rpID/4)*6 + (rpID%2)*3 + i][1]->Fill(ADC[1],TAC[1]);
        }
    }

}


void Init()
{

    hPosition = new TH2D("Position", "Position", 100,-0.1,0.1,100,-0.1,0.1);

    for (int i = 0; i < 12; ++i)
    {
        hSlewing[i][0] = new TH2D(Form("Slewing%iV",i), Form("Slewing%iV",i), 100,0,800,100,0,1400);
    }

    for (int i = 0; i < 12; ++i)
    {
        hSlewing[i][1] = new TH2D(Form("Slewing%iH",i), Form("Slewing%iH",i), 100,0,800,100,0,1400);
        
    }
}

void ConnectInput(){

    //standard reconstructed tree
    for (int i = 0; i < nRomanPots; ++i)
    {
        recTree[i] = new TTree(rpNames[i], rpNames[i]);
        recTree[i] = dynamic_cast<TTree*>( data->Get(rpNames[i]) );
        if (!recTree[i])
        {
            cout<<"Error: cannot get recTree"<<endl;
            return;
        }
        recTree[i]->SetBranchAddress("ADC_V", &ADC[0]);
        recTree[i]->SetBranchAddress("ADC_H", &ADC[1]);
        recTree[i]->SetBranchAddress("TAC_V", &TAC[0]);
        recTree[i]->SetBranchAddress("TAC_H", &TAC[1]);
        recTree[i]->SetBranchAddress("rpX", &rpX);
        recTree[i]->SetBranchAddress("rpY", &rpY);
        recTree[i]->SetBranchAddress("rpZ", &rpZ);
    }
}


