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

TString rpNames[nRomanPots] = { TString("E1U"), TString("E1D"), TString("W1U"), TString("W1D"),  TString("E2U"), TString("E2D"), TString("W2U"), TString("W2D")};
TString side[nSides] = { TString("East"), TString("West")};
TString position[2] = { TString("Up"), TString("Down")};
TString pmtID[2] = { TString("V"), TString("H")};
TString branchNames[nBranches] = { TString("EU"), TString("ED"), TString("WU"), TString("WD")};

TFile* data;
TFile* fout;

TH1D* hCorrection[16];


const int nIterations = 4;
TH2D* hDeltaTac[8][nIterations];
double tacAvrg[8][100][nIterations];
/// meaning of 8:    0 - channel1 EU  2 - channel1 ED ....
///                  1 - channel2 EU  3 - channel2 ED ....

TTree *recTree[nRomanPots];

int counter[nIterations];

Double_t ADC[2][2]; 
Double_t TAC[2][2];

Int_t branchID;

void Init();
void ConnectInput();
void Make();

void MakeIteration(int iter);
double avrg(int bin, int iteration, int ihist);
double tac(int iter,int station,int channel);

void slewing()
{
	TString input = "/home/truhlar/Desktop/STAR/CEP/Analysis/Data/slewingElastic.root";
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

    for (int iCh = 0; iCh < 8; ++iCh)
        for (int i = 0; i < 100; ++i)
            for (int j = 0; j < nIterations; ++j)    
                tacAvrg[iCh][i][j] = 0;

    for (int i = 0; i < nIterations; ++i)
    {
        counter[i] = 0;
    }
	////////////// Making histograms

    double average, norm;
    for (int iter = 0; iter < nIterations; ++iter)
    {
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
            for (int iCh = 0; iCh < 2; ++iCh)
            {
                //cout<<"____________Making fit__________"<<endl;
                TF1 *gaus = new TF1("gaus","gaus(0)",hDeltaTac[branchID*2 + iCh][iter]->GetMean(2) - 200, hDeltaTac[branchID*2 + iCh][iter]->GetMean(2) + 200);
                hDeltaTac[branchID*2 + iCh][iter]->FitSlicesY(gaus,0,-1,0,"QR");
                TH1D *h2_1 = (TH1D*)gDirectory->Get("DeltaTac_" + rpNames[i + (iter%2)*4] + Form("Channel_%i",iCh) + "_1");
                if(!h2_1)
                    cout<<"shit"<<endl;
                for (int iBin = 0; iBin < 100; ++iBin)
                {
                    tacAvrg[branchID*2 + iCh][iBin][iter] = h2_1->GetBinContent(iBin+1);
                    //cout<<"Mean by fit: "<< h2_1->GetBinContent(iBin+1)<<" for bin: "<< iBin<<endl;
                }
                    //tacAvrg[branchID*2 + iCh][iBin][iter] = avrg(iBin+1,iter,branchID*2 + iCh );
            }
            
        }
    }
    for (int iBranch = 0; iBranch < 4; ++iBranch)
    {
        for (int iRp = 0; iRp < 2; ++iRp)
        {
            for (int iChnl = 0; iChnl < 2; ++iChnl)
            {
                for (int iBin = 1; iBin < 101; ++iBin)
                {
                    double tmp = 0;
                    for (int i = iRp; i < nIterations; i+=2)
                    {
                        tmp =  tmp - tacAvrg[2*iBranch + iChnl][iBin-1][i];
                    }
                    hCorrection[iBranch*2 + iRp*8 + iChnl]->SetBinContent(iBin, tmp);
                    
                }
                
            }
        }
    }


    fout->cd();
	fout->Write();
	fout->Close();
    
}

double avrg(int bin, int iteration, int ihist)
{
    //cout <<  << endl;
    TF1 *f = new TF1("f","gaus(0)");
    hDeltaTac[ihist][iteration]->FitSlicesY();
    cout<< f->GetParameter(0) <<endl;
    return f->GetParameter(0);
}
 
void MakeIteration(int iter)
{
    int activRp = iter%2;
    int fixRp = (iter+1)%2;
    double deltaT[2]; // TAC [station][chanell]
    for (int iChnl = 0; iChnl < 2; ++iChnl)
    {
        deltaT[iChnl] = tac(iter, activRp, iChnl) - ((tac(iter, fixRp, 0) + tac(iter, fixRp, 1))/2);
        //cout<<deltaT[iChnl]<<" = "<<tac(iter, activRp, iChnl)<<" - "<<((tac(iter, fixRp, 0) + tac(iter, fixRp, 1))/2) <<endl;
        hDeltaTac[2*branchID + iChnl][iter] -> Fill(ADC[activRp][iChnl], deltaT[iChnl]); 
        if(deltaT[iChnl] < 1)
            counter[iter]++;              
    }

}

double tac(int iter, int station, int channel)
{
    double tmp = TAC[station][channel];
    for (int i = station; i < iter; i+=2)
    {
        tmp =  tmp - tacAvrg[2*branchID + channel][hDeltaTac[0][0]->GetXaxis()->FindBin(ADC[station][channel] -1)][i];
    }
    return tmp; 
}

void Make()
{


}


void Init()
{
    for (int i = 0; i < 16; ++i)
    {
        hCorrection[i] = new TH1D("corr" + rpNames[i/2] + Form("Channel_%i",i%2), "corr" + rpNames[i/2] + Form("Channel_%i",i%2),100,0,100);
    }

    for (int iter = 0; iter < nIterations; ++iter)
    {

        fout->mkdir(Form("Iter%i", iter))->cd();

        for (int i = 0; i < 8; ++i)
        {                      
                hDeltaTac[i][iter] = new TH2D("DeltaTac_" + rpNames[i/2 + (iter%2)*4] + Form("Channel_%i",i%2), "DeltaTac_" + rpNames[i/2 + (iter%2)*4] + Form("Channel_%i",i%2), 100,0,800,100,-1000,1000);                                           
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
        }
    }
}


