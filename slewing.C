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
TString side[nSides] = { TString("East"), TString("West")};
TString position[2] = { TString("Up"), TString("Down")};
TString pmtID[2] = { TString("V"), TString("H")};
TString branchNames[nBranches] = { TString("EU"), TString("ED"), TString("WU"), TString("WD")};

TFile* data;
TFile* fout;

TH2D* hPosition;
TH2D* hSlewing[24][2]; // 12 + 2xPMTs
// East - Up - 3x x-region, East - Down - 6x x-region
// West - Up - 3x x-region, West - Down - 6x x-region 

const int nIterations = 4;
TH2D* hDeltaTac[16][6][nIterations];
Float_t TACaverage[100][nIterations];


TTree *recTree[nRomanPots];


Double_t ADC[2][2]; 
Double_t TAC[2][2];
Double_t rpX[2], rpY[2], rpZ[2]; 

Double_t y[nSides][4] = 
{
    { -0.066, -0.024, 0.024, 0.062},
    { -0.057, -0.024, 0.026, 0.072}
};
Double_t x[nSides][7] = 
{
    { -0.02, -0.01, 0.0, 0.01, 0.02, 0.03, 0.045},
    { -0.016, -0.01, 0.0, 0.01, 0.02, 0.03, 0.042}
};

Int_t branchID;

void Init();
void ConnectInput();
void Make();

void MakeIteration(int iter);

void slewing()
{
	TString input = "/home/truhlar/Desktop/STAR/CEP/bigTestSlewing.root";
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

    for (int i = 0; i < 100; ++i)
    {
        for (int j = 0; j < nIterations; ++j)
        {
            
            TACaverage[i][j] = 0;
        }
    }

	////////////// Making histograms
	for (int i = 0; i < nBranches; ++i)
	{ //get the event   
        branchID = i; 
    	Long64_t nev = recTree[i]->GetEntries();
    	cout<<"Proccesing branch " + branchNames[i] + " with "<<nev<<" events"<<endl;
        for(Long64_t iev=0; iev<nev; ++iev) 
        {
            recTree[i]->GetEntry(iev); 
            Make();
            MakeIteration(0);
        }	   
    } 

    double average, norm;
    for (int iter = 1; iter < nIterations; ++iter)
    {
        for (int i = 0; i < 16; ++i)
        {
            for (int j = 0; j < 6; ++j)
            { 
                for(int adcBin = 0; adcBin < 100; ++adcBin)
                {
                    average = 0;
                    norm = 0;
                    for(int dTACbin = 0; dTACbin < 100; ++dTACbin)
                    {
                        average += (20*dTACbin - 995)*hDeltaTac[i][j][iter - 1]->GetBinContent(adcBin,dTACbin);
                        norm += hDeltaTac[i][j][iter - 1]->GetBinContent(adcBin,dTACbin);
                    }
                    if(!norm)
                        continue;
                    for (int ii = iter; ii < nIterations; ++ii)
                    {
                        TACaverage[adcBin][ii] += pow(-1,iter)*(average / norm);
                    }
                }
            }
        }
        for (int i = 0; i < nBranches; ++i)
        { //get the event   
            branchID = i; 
            Long64_t nev = recTree[i]->GetEntries();
            cout<<"Proccesing branch " + branchNames[i] + " with "<<iter<<" iteration"<<endl;
            for(Long64_t iev=0; iev<nev; ++iev) 
            {
                recTree[i]->GetEntry(iev); 
                MakeIteration(iter);
            }      
        }
    }
    fout->cd();
	fout->Write();
	fout->Close();
    
}


void MakeIteration(int iter)
{

    for (int i = 0; i < 4; ++i)
    {
        for (int k = 0; k < 6; ++k)
        {
            if(rpX[0] > x[branchID/2][k] && rpX[0] < x[branchID/2][k+1] && rpX[1] > x[branchID/2][k] && rpX[1] < x[branchID/2][k+1])
            {
                hDeltaTac[4*branchID + i][k][iter]->Fill(ADC[i%2][i/2], TAC[1][i/2] - TAC[0][i/2] + TACaverage[int(ADC[i%2][i/2])/8][iter]);
            }
        }
    }

}

void Make()
{
    for (int i = 0; i < 2; ++i)
    {
        hPosition->Fill(rpX[i],rpY[i]);
        
        if(rpY[i] > y[branchID/2][2 - (branchID%2)*2] && rpY[i] < y[branchID/2][3 - (branchID%2)*2])
        {
            for (int k = 0; k < 6; ++k)
            {
                if(rpX[i] > x[branchID/2][k] && rpX[i] < x[branchID/2][k+1])
                {
                    hSlewing[(branchID/2)*12 + (branchID%2)*6 + k][0]->Fill(ADC[i][0],TAC[i][0]);
                    hSlewing[(branchID/2)*12 + (branchID%2)*6 + k][1]->Fill(ADC[i][1],TAC[i][1]);
                }
            }
        }
    }

}


void Init()
{

    hPosition = new TH2D("Position", "Position", 100,-0.1,0.1,100,-0.1,0.1);

    for (int i = 0; i < 24; ++i)
    {
        hSlewing[i][0] = new TH2D("Slewing_" + side[i/12] + position[(i/6)%2] + Form("(%.2lf-%.2lf)",x[i/12][(i%6)],x[i/12][(i%6) + 1]) + "_V", "Slewing " + side[i/12] + position[(i/6)%2] + Form("(%.2lf-%.2lf)",x[i/12][(i%6)],x[(i/6)%2][(i%6) + 1]) + "_V", 100,0,800,100,0,1400);
    }

    for (int i = 0; i < 24; ++i)
    {
        hSlewing[i][1] = new TH2D("Slewing_" + side[i/12] + position[(i/6)%2] + Form("(%.2lf-%.2lf)",x[i/12][(i%6)],x[i/12][(i%6) + 1]) + "_H", "Slewing " + side[i/12] + position[(i/6)%2] + Form("(%.2lf-%.2lf)",x[i/12][(i%6)],x[(i/6)%2][(i%6) + 1]) + "_H", 100,0,800,100,0,1400);
        
    }

    for (int iter = 0; iter < nIterations; ++iter)
    {

        fout->mkdir(Form("Iter%i", iter))->cd();

        for (int i = 0; i < 16; ++i)
        {
            for (int j = 0; j < 6; ++j)
            {
                hDeltaTac[i][j][iter] = new TH2D("DeltaTac_" + rpNames[(i/8)*4 + (i/4)%2] + "-" + rpNames[(i/8)*4 + (i/4)%2 + 2] +  Form("(%.2lf-%.2lf)",x[i/8][j],x[i/8][j+1]) + pmtID[(i%4)/2] + Form("%i",(i%4)%2 +1), "DeltaTac_" + rpNames[(i/8)*4 + (i/4)%2] + "-" + rpNames[(i/8)*4 + (i/4)%2 + 2] +  Form("(%.2lf-%.2lf)",x[i/8][j],x[i/8][j+1]) + pmtID[(i%4)/2] + Form("%i",(i%4)%2 +1), 100,0,800,100,-1000,1000);
            }                                            
        }
    }
}



void ConnectInput(){

    //standard reconstructed tree
    for (int i = 0; i < nBranches; ++i)
    {
        recTree[i] = new TTree(branchNames[i], branchNames[i]);
        recTree[i] = dynamic_cast<TTree*>( data->Get(branchNames[i]) );
        if (!recTree[i])
        {
            cout<<"Error: cannot get recTree"<<endl;
            return;
        }
        for (int j = 0; j < 2; ++j)
        {
            recTree[i]->SetBranchAddress(Form("ADC_V_%i",j+1), &ADC[j][0]);
            recTree[i]->SetBranchAddress(Form("ADC_H_%i",j+1), &ADC[j][1]);
            recTree[i]->SetBranchAddress(Form("TAC_V_%i",j+1), &TAC[j][0]);
            recTree[i]->SetBranchAddress(Form("TAC_H_%i",j+1), &TAC[j][1]);
            recTree[i]->SetBranchAddress(Form("rpX_%i",j+1), &rpX[j]);
            recTree[i]->SetBranchAddress(Form("rpY_%i",j+1), &rpY[j]);
            recTree[i]->SetBranchAddress(Form("rpZ_%i",j+1), &rpZ[j]);
        }
    }
}


