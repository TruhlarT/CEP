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
#include <TH3.h> 
#include <THStack.h> 
#include <TStyle.h> 
#include <TGraph.h> 
#include <TGraph2D.h> 
#include <TGraphErrors.h> 
#include <TCanvas.h> 
#include <TLegend.h> 
#include <TPad.h>
#include <TMinuit.h>
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

#include "BasicPlots.h"
#include "PID.h"
#include "trackQuality.h"
#include "Plot.h"
#include "FourPi.h"
#include "RPplots.h"

using namespace std;

enum PARTICLES {Pion = 0, Kaon = 1, Proton = 2, nParticles};
enum COMBINATIONS {ElInEl = 0, El = 1, InEl = 2, nCombination};

TString particleLables[nParticles] = { TString("Pion"), TString("Kaon"), TString("Proton")};
TString combinationLabel[nCombination] = { TString("El+InEl"), TString("El"), TString("InEl")};

TFile* TPCeff;
TFile* TOFeff;
TFile* fout;

TTree* tree;
TTree* treeBack;

TH3D* hTOFeff[6]; // 0 = pi- 1 = K- 2 = pbar
TH3D* hTPCeff[6]; // 3 = pi+ 4 = K+ 5 = p
TH1D* hInvMassCorr[nParticles][2][3]; // 0 - signal, 1 - Background
                                    //  - El + InEl, 1 - El, 2 - InEl 
TH1D* hInvMassUncorr[nParticles][2][3]; // 0 - signal, 1 - Background

TH2D* hMomentum[nParticles];
TH2D* hTransMomentum[nParticles];

Double_t chiPair[nParticles]; 
Double_t invMass[nParticles];
Double_t mSquared;
Double_t momentum[4];
Double_t transMomentum[4];
Double_t charge[4];
Double_t nSigmaTPC[nParticles][4];
Double_t vertexesZ[4];
Double_t DcaXY[4];
Double_t DcaZ[4];
Double_t NhitsFit[4];
Double_t NhitsDEdx[4];
Double_t Eta[4];
Double_t Phi[4];

Bool_t fourPiState;
Bool_t elastic;

void Init();
void ConnectInput(TTree* tree);
void Make(int signal);
void PlotMoneyPlot();

//_____________________________________________________________________________
int main(int argc, char** argv) {
	//open output file

	bool showText = false;
	bool showCutsLine = true;

	if(argc != 2){
		cout<<"Error: wrong input"<<endl;
		cout<<"You should do: ./Analysis DataSet"<<endl;
		return 1;
	}

	//const string dataName = "StRP_production_0000.root";
	const string& DataSet = argv[1];
	TString dataName = "/home/truhlar/Desktop/STAR/CEP/Analysis/Data/" + DataSet + ".root";
	TString output = "/home/truhlar/Desktop/STAR/CEP/Analysis/Outputs/" + DataSet + "/";
	

	TString cutsOption[] = { TString("vertexZ<80 && vertexZ > -80"), TString("NhitsFit1 >=25 && NhitsFit0 >= 25"), TString("NhitsDEdx1 >= 15 && NhitsDEdx0 >= 15"), TString("DcaZ1 < 1 && DcaZ1 > -1 && DcaZ0 < 1 && DcaZ0 > -1"), TString("DcaXY1 < 1.5 && DcaXY0 < 1.5"), TString("Eta1 > -0.7 && Eta1 < 0.7 && Eta0 > -0.7 && Eta0 < 0.7"), TString("!fourPiState")};
	TString cutsLabels[] = { TString("|z_{vtx}| < 80 cm"), TString("N_{hits}^{fit} #geq 25"), TString("N_{hits}^{dE/dx} #geq 15"), TString("|DCA(z)| < 1 cm"), TString("DCA(XY) < 1.5 cm"), TString("|#eta| < 0.7"), TString("!fourPiState")  };

	TFile* data = TFile::Open(dataName, "read");
	if (!data){
		cout<<"Error: cannot open: "<<dataName<<endl;
		return 2;
    }

	tree = dynamic_cast<TTree*>( data->Get("recTree") );
	treeBack = dynamic_cast<TTree*>( data->Get("Background") );
	
	if (!tree || !treeBack){
		cout<<"Error: cannot open one of the TTree"<<endl;
		return 3;
    }

    TString TPCeffInput = "/home/truhlar/Desktop/STAR/CEP/Analysis/Data/etaPhiEfficiency_16_01_19_delta015_twoRuns.root";
    TString TOFeffInput = "/home/truhlar/Desktop/STAR/CEP/Analysis/Data/effWithBinningForSystematics.root";

    TPCeff = TFile::Open(TPCeffInput, "read");
    if (!TPCeff)
    {
        cout<<"Error: cannot open "<<TPCeffInput<<endl;
        return 4;
    }

    TOFeff = TFile::Open(TOFeffInput, "read");
    if (!TOFeff)
    {
        cout<<"Error: cannot open "<<TOFeffInput<<endl;
        return 5;
    }

    for (int i = 0; i < 6; ++i)
    {    
        hTPCeff[i] = (TH3D*)TPCeff -> Get(Form("hTPCEffiCD%i121",i)); // hTPCEffiCD%i120" for dead sector 19
        hTOFeff[i] = (TH3D*)TOFeff -> Get(Form("hTOFEffiCD%i12",i)); 
    }

	fout = new TFile(output +"StRP.root","RECREATE");

	Plot tool;

//////////////////////////////////////////////////////////////////////
//					No cuts applied 
//////////////////////////////////////////////////////////////////////
	
    fout->mkdir("BasicPlots")->cd();
	BasicPlots Plots(data, fout, output, showText);
	Plots.PlotHistogram();

    TDirectory* curDir = fout->mkdir("El+InEl");
    curDir->cd();
    TString usedCuts = "";
    for (int i = 0; i < nCombination; ++i)
    {
        if(i == El)
        {
            curDir = fout->mkdir("El");
            curDir->cd();
            usedCuts = "elastic";
        }else if( i == InEl)
        {
            curDir = fout->mkdir("InEl");
            curDir->cd();
            usedCuts = "!elastic";
        }

        curDir->mkdir("PID")->cd();
        PID PIDPlotsWithCuts(data, fout, output, showCutsLine, usedCuts);
        PIDPlotsWithCuts.PlotHistogram();

        curDir->mkdir("trackQuality")->cd();
        trackQuality TrackPlotsWithCuts(data, fout, output, showCutsLine, usedCuts);
        TrackPlotsWithCuts.PlotHistogram();

        curDir->mkdir("RPplots")->cd();
        RPplots RPplotsWithCuts(data, fout, output, showCutsLine, usedCuts);
        RPplotsWithCuts.PlotHistogram();
    }

//////////////////////////////////////////////////////////////////////
//				All cuts applied
//////////////////////////////////////////////////////////////////////

	TDirectory* cutsDir = fout->mkdir("Cuts");
	cutsDir->cd();

	TString cuts = "";
	int size = (sizeof(cutsOption)/sizeof(*cutsOption));

	TCanvas* newCanvas = new TCanvas("newCanvas","newCanvas",800,700);
	gPad->SetMargin(0.11,0.02,0.105,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
	gPad->SetTickx();
	gPad->SetTicky(); 
	gPad->SetLogy();
	gStyle->SetOptStat("");
	gStyle->SetPalette(1);
	gStyle->SetLineWidth(2);      //axis line
	gStyle->SetFrameLineWidth(2); //frame line

	TH1F* hCutsFlow = new TH1F("cuts", "cuts", size, 0, size-1);
	
	tree->Draw("invMassPion>>invMassPionSignal");
	TH1F *tmpHist = (TH1F*)gPad->GetPrimitive("invMassPionSignal");
	hCutsFlow->SetBinContent(1,tmpHist->GetEntries());
	hCutsFlow->GetXaxis()->SetBinLabel(1, "excl. events");
	
	for (int i = 0; i < size; ++i)
	{
		cuts += cutsOption[i]; 
        cout<<"Applying cuts: "<<cuts<<endl;
		tree->Draw("invMassPion>>invMassPionSignal",cuts);
		TH1F *tmpHist = (TH1F*)gPad->GetPrimitive("invMassPionSignal");
		if(i != size -1)
        {
            hCutsFlow->GetXaxis()->SetBinLabel(i+2, cutsLabels[i]);
            hCutsFlow->SetBinContent(i+2,tmpHist->GetEntries());
            cout<<"Number of entries: "<<tmpHist->GetEntries()<<endl;
        }

        cutsDir->mkdir(Form("trckQ_%i",i))->cd();
        trackQuality pokus(data, fout, output, showCutsLine, cuts);
        pokus.PlotHistogram();

        if(i != size -1)
            cuts += " && ";
	}
    cutsDir->cd();
	cout<<cuts<<endl;

	hCutsFlow->SetTitle("; ; Number of events");
	tool.SetGraphStyle(hCutsFlow);
	tool.SetMarkerStyle(hCutsFlow);
	//hCutsFlow->GetYaxis()->SetRangeUser(800,300000000);
	hCutsFlow->GetXaxis()->SetLabelSize(0.045);
	hCutsFlow->Draw();
	tool.DrawTextStar(hCutsFlow,2);

	TPaveText *textPub = new TPaveText(0.75,0.74,0.9,0.9,"brNDC");
	tool.SetTextStyle(textPub);
	textPub -> AddText("p + p #rightarrow p + X + p");
	textPub -> AddText("#sqrt{s} = 510 GeV");
	textPub -> AddText("Cuts flow");
	textPub -> Draw("same");

	newCanvas->Update();
	//newCanvas->SaveAs(output + "BasicPlots/Cuts.png");
	newCanvas->Write("CutsFlow");

    TDirectory* goldenDir = fout->mkdir("GoldenEvents");
    goldenDir->cd();
    TDirectory* currentDir = goldenDir->mkdir("El+InEl");
    currentDir->cd();

    usedCuts = cuts;
    for (int i = 0; i < nCombination; ++i)
    {
        if(i == El)
        {
            currentDir = goldenDir->mkdir("El");
            currentDir->cd();
            usedCuts = cuts + "&& elastic";
        }else if( i == InEl)
        {
            currentDir = goldenDir->mkdir("InEl");
            currentDir->cd();
            usedCuts = cuts + "&& !elastic";
        }

    	currentDir->mkdir("PID")->cd();
    	PID PIDPlotsWithCuts(data, fout, output, showCutsLine, usedCuts);
    	PIDPlotsWithCuts.PlotHistogram();

        currentDir->mkdir("trackQuality")->cd();
        trackQuality TrackPlotsWithCuts(data, fout, output, showCutsLine, usedCuts);
        TrackPlotsWithCuts.PlotHistogram();

        currentDir->mkdir("RPplots")->cd();
        RPplots RPplotsWithCuts(data, fout, output, showCutsLine, usedCuts);
        RPplotsWithCuts.PlotHistogram();
    }


///////////////////////////////////////////////////////////////
//    TDirectory* TomasDir = cutsDir->mkdir("PID_Tomas");
//    TomasDir->cd();
//    TomasDir->mkdir("Pions")->cd();
//    PID PIDPlotsWithCuts1(data, fout, output, showCutsLine, cuts + "&& nSigTrk1Pion < 3 && nSigTrk1Pion > -3 && nSigTrk2Pion > -3 && nSigTrk2Pion < 3 && (chiPairKaon > 3 || mSquared < 0.2 || mSquared > 0.32) && (chiPairProton > 3 || mSquared < 0.7 || mSquared > 1.1)");
//    PIDPlotsWithCuts1.PlotHistogram();
//  
//    TomasDir->mkdir("Kaons")->cd();
//    PID PIDPlotsWithCuts2(data, fout, output, showCutsLine, cuts + "&& chiPairKaon < 3 && mSquared > 0.2 && mSquared < 0.32 && (chiPairProton > 3 || mSquared < 0.7 || mSquared > 1.1)");
//    PIDPlotsWithCuts2.PlotHistogram();
//  
//    TomasDir->mkdir("Protons")->cd();
//    PID PIDPlotsWithCuts3(data, fout, output, showCutsLine, cuts + "&& chiPairProton < 3 && mSquared > 0.7 && mSquared < 1.1");
//    PIDPlotsWithCuts3.PlotHistogram();
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//    TDirectory* TomasDir2 = cutsDir->mkdir("PID_Tomas2");
//    TomasDir2->cd();
//    TomasDir2->mkdir("Pions")->cd();
//    PID PIDPlotsWithCuts4(data, fout, output, showCutsLine, cuts + "&& chiPairPion < 3.46 && (chiPairPion < 3 || chiPairKaon > 3 || mSquared < 0.15) && (chiPairPion < 3 || chiPairKaon < 3 || chiPairProton > 3 || mSquared < 0.6)");
//    PIDPlotsWithCuts4.PlotHistogram();
//
//    TomasDir2->mkdir("Kaons")->cd();
//    PID PIDPlotsWithCuts5(data, fout, output, showCutsLine, cuts + "&& chiPairPion > 3 && chiPairKaon < 3 && mSquared > 0.15 && (chiPairPion < 3 || chiPairKaon < 3 || chiPairProton > 3 || mSquared < 0.6)");
//    PIDPlotsWithCuts5.PlotHistogram();
//
//    TomasDir2->mkdir("Protons")->cd();
//    PID PIDPlotsWithCuts6(data, fout, output, showCutsLine, cuts + "&& chiPairPion > 3 && chiPairKaon > 3 && chiPairProton < 3 && mSquared > 0.6");
//    PIDPlotsWithCuts6.PlotHistogram();
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//    TDirectory* RafalDir = cutsDir->mkdir("PID_Rafal");
//    RafalDir->cd();
//    RafalDir->mkdir("Pions")->cd();
//    PID PIDPlotsWithCuts7(data, fout, output, showCutsLine, cuts + "&& nSigTrk1Pion < 3 && nSigTrk1Pion > -3 && nSigTrk2Pion > -3 && nSigTrk2Pion < 3 && (chiPairPion < 3 || chiPairKaon > 3 || chiPairProton < 3 || mSquared < 0.2 || mSquared > 0.32) && (chiPairPion < 3 || chiPairKaon < 3 || chiPairProton > 3 || mSquared < 0.7 || mSquared > 1.1)");
//    PIDPlotsWithCuts7.PlotHistogram();
//
//    RafalDir->mkdir("Kaons")->cd();
//    PID PIDPlotsWithCuts8(data, fout, output, showCutsLine, cuts + "&& chiPairPion > 3 && chiPairKaon < 3 && chiPairProton > 3 && mSquared > 0.2 && mSquared < 0.32 && (chiPairPion < 3 || chiPairKaon < 3 || chiPairProton > 3 || mSquared < 0.7 || mSquared > 1.1)");
//    PIDPlotsWithCuts8.PlotHistogram();
//
//    RafalDir->mkdir("Protons")->cd();
//    PID PIDPlotsWithCuts9(data, fout, output, showCutsLine, cuts + "&& chiPairPion > 3 && chiPairKaon > 3 && chiPairProton < 3 && mSquared > 0.6 && mSquared < 1.1");
//    PIDPlotsWithCuts9.PlotHistogram();
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//    TDirectory* DanielDir = cutsDir->mkdir("PID_Daniel");
//    DanielDir->cd();
//    DanielDir->mkdir("Pions")->cd();
//    PID PIDPlotsWithCuts10(data, fout, output, showCutsLine, cuts + "&& nSigTrk1Pion < 3 && nSigTrk1Pion > -3 && nSigTrk2Pion > -3 && nSigTrk2Pion < 3 && (chiPairKaon > 3 || deltaDeltaTOFKaon < -0.5 || deltaDeltaTOFKaon > 0.5) && (chiPairProton > 3 || deltaDeltaTOFProton < -0.5 || deltaDeltaTOFProton > 0.5)");
//    PIDPlotsWithCuts10.PlotHistogram();
//
//    DanielDir->mkdir("Kaons")->cd();
//    PID PIDPlotsWithCuts11(data, fout, output, showCutsLine, cuts + "&& chiPairKaon < 3 && deltaDeltaTOFKaon > -0.5 && deltaDeltaTOFKaon < 0.5  && (chiPairProton > 3 || deltaDeltaTOFProton < -0.5 || deltaDeltaTOFProton > 0.5)");
//    PIDPlotsWithCuts11.PlotHistogram();
//
//    DanielDir->mkdir("Protons")->cd();
//    PID PIDPlotsWithCuts12(data, fout, output, showCutsLine, cuts + "&& chiPairProton < 3 && deltaDeltaTOFProton > -0.5 && deltaDeltaTOFProton < 0.5");
//    PIDPlotsWithCuts12.PlotHistogram();
///////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////
//				PID cuts applied
//////////////////////////////////////////////////////////////////////

    Init();
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

    PlotMoneyPlot();



    TCanvas *cCanvas2D = new TCanvas("cCanvas2D","cCanvas2D",800,700);
    gPad->SetMargin(0.09,0.13,0.1,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
    gStyle->SetPalette(1);
    for (int i = 0; i < nParticles; ++i)
    {
        hMomentum[i]->Draw("colz");
        cCanvas2D->Write("momentum" + particleLables[i]);  
        hTransMomentum[i]->Draw("colz");
        cCanvas2D->Write("transMomentum" + particleLables[i]);
    }


//////////////////////////////////////////////////////////////////////
//              4 pions state
//////////////////////////////////////////////////////////////////////
 
    TDirectory* fourPiDir = fout->mkdir("4pions");
    fourPiDir->cd();

    FourPi fourPiPlots(data, fout, output, showText);
    fourPiPlots.PlotHistogram();
    
//////////////////////////////////////////////////////////////////////

	fout->Close();
	data->Close();

	cout<<"Ending Analysis... GOOD BYE!"<<endl;
	return 0;
}//main


void Init()
{
    for (int i = 0; i < nParticles; ++i)
    {
        hMomentum[i] = new TH2D("momentumCorr" + particleLables[i], "momentumCorr" + particleLables[i], 100,0,3,100,0,3);
        hTransMomentum[i] = new TH2D("transMomentumCorr" + particleLables[i], "transMomentumCorr" + particleLables[i], 100,0,3,100,0,3);
    }
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


}


void ConnectInput(TTree* tree)
{

// PID and some quality event info
    tree->SetBranchAddress("mSquared", &mSquared);
    tree->SetBranchAddress("elastic", &elastic);
    tree->SetBranchAddress("nSigTrk1Pion", &nSigmaTPC[Pion][0]);
    tree->SetBranchAddress("nSigTrk2Pion", &nSigmaTPC[Pion][1]);
    tree->SetBranchAddress("mSquared", &mSquared);
    for (int iPart = 0; iPart < nParticles; ++iPart)
    {
        tree->SetBranchAddress("invMass" + particleLables[iPart], &invMass[iPart]);
        tree->SetBranchAddress("chiPair" + particleLables[iPart], &chiPair[iPart]);
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

void Make(int signal)
{
    double effTotal, effTPC, effTOF;
    unsigned int PID;
   // cout<< vertexesZ[0] <<" "<< NhitsFit[0]<<" "<<NhitsFit[1] <<" "<< NhitsDEdx[0]<<" "<<NhitsDEdx[1] <<" "<<DcaZ[0] <<" "<<DcaZ[1] <<" "<<DcaXY[0] <<" "<<DcaXY[1] <<" "<<Eta[0] <<" "<<Eta[1] <<" "<< !fourPiState<<endl; 
    if(vertexesZ[0] < 80 && vertexesZ[0] > -80 && NhitsFit[0] >=25 && NhitsFit[1] >= 25 && NhitsDEdx[0] >= 15 && NhitsDEdx[1] >= 15 && DcaZ[0] < 1 && DcaZ[0] > -1 && DcaZ[1] < 1 && DcaZ[1] > -1 && DcaXY[0] < 1.5 && DcaXY[1] < 1.5 && Eta[0] > -0.7 && Eta[0] < 0.7 && Eta[1] > -0.7 && Eta[1] < 0.7  && !fourPiState)
    {
        int combination = 2;
        if(elastic)
            combination = 1;
        effTotal = 1;
        if(chiPair[Pion] > 9 && chiPair[Kaon] > 9 && chiPair[Proton] < 9 && mSquared > 0.6) // it is... proton!
        {
            for (int iTrack = 0; iTrack < 2; ++iTrack)
            {
                PID = 0;
                if(charge[iTrack] > 0)
                    PID = 3;
                effTPC = hTPCeff[2 + PID]->GetBinContent( hTPCeff[2 + PID]->GetXaxis()->FindBin(vertexesZ[iTrack]), hTPCeff[2 + PID]->GetYaxis()->FindBin(transMomentum[iTrack]), hTPCeff[2 + PID]->GetZaxis()->FindBin(Eta[iTrack]));
                effTOF = hTOFeff[2 + PID]->GetBinContent( hTOFeff[2 + PID]->GetXaxis()->FindBin(vertexesZ[iTrack]), hTOFeff[2 + PID]->GetYaxis()->FindBin(transMomentum[iTrack]), hTOFeff[2 + PID]->GetZaxis()->FindBin(Eta[iTrack]));  
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
                effTPC = hTPCeff[1 + PID]->GetBinContent( hTPCeff[1 + PID]->GetXaxis()->FindBin(vertexesZ[iTrack]), hTPCeff[1 + PID]->GetYaxis()->FindBin(transMomentum[iTrack]), hTPCeff[1 + PID]->GetZaxis()->FindBin(Eta[iTrack]));
                effTOF = hTOFeff[1 + PID]->GetBinContent( hTOFeff[1 + PID]->GetXaxis()->FindBin(vertexesZ[iTrack]), hTOFeff[1 + PID]->GetYaxis()->FindBin(transMomentum[iTrack]), hTOFeff[1 + PID]->GetZaxis()->FindBin(Eta[iTrack]));  
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
                effTPC = hTPCeff[0 + PID]->GetBinContent( hTPCeff[0 + PID]->GetXaxis()->FindBin(vertexesZ[iTrack]), hTPCeff[0 + PID]->GetYaxis()->FindBin(transMomentum[iTrack]), hTPCeff[0 + PID]->GetZaxis()->FindBin(Eta[iTrack]));
                effTOF = hTOFeff[0 + PID]->GetBinContent( hTOFeff[0 + PID]->GetXaxis()->FindBin(vertexesZ[iTrack]), hTOFeff[0 + PID]->GetYaxis()->FindBin(transMomentum[iTrack]), hTOFeff[0 + PID]->GetZaxis()->FindBin(Eta[iTrack])); 
                
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


void PlotMoneyPlot()
{
    TCanvas* newCanvas = new TCanvas("newCanvas","newCanvas",800,700);
    gPad->SetMargin(0.11,0.02,0.105,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
    gPad->SetTickx();
    gPad->SetTicky(); 
    gPad->SetLogy(0);
    gStyle->SetOptStat("");
    gStyle->SetPalette(1);
    gStyle->SetLineWidth(2);      //axis line
    gStyle->SetFrameLineWidth(2); //frame line

    TString stateLabel[] = {"#pi^{+}#pi^{-}","K^{+}K^{-}","p#bar{p}"}; 
    Plot tool;

    TH1D *hist, *histCompare;
    TString yLabel, name;
    yLabel = "Number of events";
    name = "uncorrInvMass";


    TDirectory* moneyDir = fout->mkdir("MoneyPlots");
    moneyDir->cd();
    moneyDir->mkdir("Uncorrected")->cd();

    for (int corrIndex = 0; corrIndex < 2; ++corrIndex) // 0 - uncorrected, 1 - corrected
    {
        if(corrIndex == 1)
        {
            yLabel = "Corrected counts";
            name = "corrInvMass";
            moneyDir->mkdir("Corrected")->cd();
        }

        for (int i = 0; i < nParticles; ++i)
        { 
            for (int comb = 0; comb < nCombination; ++comb)
            { 
                hist = hInvMassUncorr[i][0][comb];
                histCompare = hInvMassUncorr[i][1][comb];
                hist->SetTitle(" ; m(" + stateLabel[i] + ") [GeV/c^{2}]; " + yLabel);
                tool.SetGraphStyle(hist,4,20,1,4,1,1,0.9,1.3);
                tool.SetMarkerStyle(hist);
                hist->Draw("E");
                tool.DrawText(hist, i+1, true);
                tool.DrawTextStar(hist, 1);
                tool.SetMarkerStyle(histCompare,2,20,1,2,1,1);
                histCompare->Draw("ESAME");

                TLegend* leg1 = new TLegend(0.58, 0.64, 0.78, 0.74);
                tool.SetLegendStyle(leg1);
                leg1->AddEntry(hist, combinationLabel[comb] + " (unlike-sign pairs)","p");
                leg1->AddEntry(histCompare, combinationLabel[comb] + " (like-sign pairs)","p");
                leg1->Draw("same");

                newCanvas->Update();
                newCanvas->Write(name + particleLables[i] + combinationLabel[comb]);
            }
        }
    }

    moneyDir->mkdir("RatioPlots")->cd();

    gPad->SetMargin(0.9,0.02,0.3,0.02);
    gStyle->SetPadTickY(1);
    gStyle->SetTickLength(0.02,"Y");
    gStyle->SetTickLength(0.03,"X");
    // Upper plot will be in pad1
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.295, 1, 1.0);
    pad1->SetTopMargin(0.04);
    pad1->SetRightMargin(0.02);
    pad1->SetBottomMargin(0); // Upper and lower plot are joined
    pad1->Draw();             // Draw the upper pad: pad1
             
    newCanvas->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.00, 1, 0.2985);
    pad2->SetTopMargin(0);
    pad2->SetRightMargin(0.02);
    pad2->SetBottomMargin(0.3);
    pad2->Draw();

    name = "ratioPlot";
    yLabel = "Normalized counts";
    for (int i = 0; i < nParticles; ++i)
    { 
        for (int comb = 0; comb < nCombination; ++comb)
        { 
            pad1->cd();
            hist = hInvMassCorr[i][0][comb];
            histCompare = hInvMassUncorr[i][0][comb];
            
            tool.SetMarkerStyle(hist);
            tool.SetMarkerStyle(histCompare,2,20,1,2,1,1);

            hist->GetXaxis()->SetTitleFont(43);
            hist->GetYaxis()->SetTitleFont(43);
            hist->GetXaxis()->SetTitleSize(25);
            hist->GetYaxis()->SetTitleSize(25);
            hist->GetXaxis()->SetTitleOffset(0.9);
            hist->GetYaxis()->SetTitleOffset(1.5);
            hist->GetYaxis()->SetLabelFont(43);
            hist->GetYaxis()->SetLabelSize(20);

            Double_t scaleFactor =   1 /hist->Integral();
            hist->Scale(scaleFactor);
            scaleFactor =   1 /histCompare->Integral();
            histCompare->Scale(scaleFactor);

            hist->GetYaxis()->SetRangeUser(0.001, hist->GetMaximum()*1.2);

            hist->SetTitle(" ; m(" + stateLabel[i] + ") [GeV/c^{2}] ; " + yLabel);
            hist->Draw();              
            histCompare->Draw("same");         

            TPaveText *textSTAR = new TPaveText(0.75,0.89,0.85,0.95,"brNDC"); 
            textSTAR -> SetTextSize(0.05);
            textSTAR -> SetFillColor(0);
            textSTAR -> SetTextFont(62);
            textSTAR -> AddText("STAR Internal");
            textSTAR -> Draw("same");

            textSTAR = new TPaveText(0.75,0.74,0.85,0.88,"brNDC"); 
            textSTAR -> SetTextSize(25);
            textSTAR -> SetFillColor(0);
            textSTAR -> SetTextFont(43);
            textSTAR -> AddText("p + p #rightarrow p + " + stateLabel[i] +" + p");
            textSTAR -> AddText("#sqrt{s} = 510 GeV");
            textSTAR -> Draw("same");

            TLegend *legendPID = new TLegend(0.68,0.58,0.92,0.73,"","brNDC");
            tool.SetLegendStyle(legendPID);
            legendPID -> SetTextSize(25);
            legendPID -> SetTextFont(43);
            legendPID -> AddEntry(hist, "Corrected", "p");
            legendPID -> AddEntry(histCompare, "Uncorrected", "p");
            legendPID -> Draw("same");

            pad2->cd();       

            // Define the ratio plot
            TH1D *h3 = (TH1D*)hist->Clone("h3");
            h3->GetYaxis()->SetLabelSize(0.);
            h3->SetMarkerColor(kBlack);
            h3->SetMinimum(0.2);  // Define Y ..
            h3->SetMaximum(2.2); // .. range
            h3->Divide(histCompare);
            h3->SetMarkerStyle(20);
            h3->Draw("ep");       // Draw the ratio plot

            // Ratio plot (h3) settings
            h3->SetTitle(""); // Remove the ratio title

            // Y axis ratio plot settings
            h3->GetYaxis()->SetTitle("Ratio");
            h3->GetYaxis()->SetNdivisions(505);
            h3->GetYaxis()->SetTitleSize(25);
            h3->GetYaxis()->SetTitleFont(43);
            h3->GetYaxis()->SetTitleOffset(1.5);
            h3->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
            h3->GetYaxis()->SetLabelSize(0.0);

            // X axis ratio plot settings
            h3->GetXaxis()->SetTitleSize(25);
            h3->GetXaxis()->SetTitleFont(43);
            h3->GetXaxis()->SetTitleOffset(3.2);
            h3->GetXaxis()->SetLabelOffset(0.025);
            h3->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
            h3->GetXaxis()->SetLabelSize(20);

            TGaxis *axis = new TGaxis(h3->GetXaxis()->GetXmin(), 0.2, h3->GetXaxis()->GetXmin(), 2.2, 0.2, 2.2, 5, "");
            axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
            axis->SetLabelSize(20);
            axis->Draw("same");

            TLine *unity = new TLine(h3->GetXaxis()->GetXmin(), 1., h3->GetXaxis()->GetXmax() , 1.);
            unity->SetLineColor(kBlack);
            unity->SetLineWidth(2);
            unity->Draw();

            newCanvas->Update();
            newCanvas->Write(name + particleLables[i] + combinationLabel[comb]);
        }
    }

}