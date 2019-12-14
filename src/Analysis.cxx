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

using namespace std;

enum PARTICLES {Pion = 0, Kaon = 1, Proton = 2, nParticles};
TString particleLables[nParticles] = { TString("Pion"), TString("Kaon"), TString("Proton")};

TFile* TPCeff;
TFile* TOFeff;

TTree* tree;
TTree* treeBack;

TH3D* hTOFeff[6]; // 0 = pi- 1 = K- 2 = pbar
TH3D* hTPCeff[6]; // 3 = pi+ 4 = K+ 5 = p
TH1D* hInvMassCorr[nParticles][2]; // 0 - signal, 1 - Background
TH1D* hInvMassUncorr[nParticles][2]; // 0 - signal, 1 - Background

Double_t nSigPair[nParticles]; 
Double_t invMass[nParticles];
Double_t mSquared;
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

void Init();
void ConnectInput(TTree* tree);
void Make(int signal);

//_____________________________________________________________________________
int main(int argc, char** argv) {
	//open output file

	bool showText = false;
	bool showCutsLine = false;

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
        hTPCeff[i] = (TH3D*)TPCeff -> Get(Form("hTPCEffiCD%i120",i));
        hTOFeff[i] = (TH3D*)TOFeff -> Get(Form("hTOFEffiCD%i12",i)); 
    }


	TFile* fout = new TFile(output +"StRP.root","RECREATE");

	Plot tool;

//////////////////////////////////////////////////////////////////////
//					No cuts applied 
//////////////////////////////////////////////////////////////////////
	fout->mkdir("BasicPlots")->cd();
	BasicPlots Plots(data, fout, output, showText);
	Plots.PlotHistogram();

	fout->mkdir("PID")->cd();
	PID PIDPlots(data, fout, output, showCutsLine);
	PIDPlots.PlotHistogram();

	fout->mkdir("trackQuality")->cd();
	trackQuality TrackPlots(data, fout, output, showCutsLine);
	TrackPlots.PlotHistogram();

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

	TH1F* hCutsFlow = new TH1F("cuts", "cuts", size+1, 0, size);
	
	tree->Draw("invMassPion>>invMassPionSignal");
	TH1F *tmpHist = (TH1F*)gPad->GetPrimitive("invMassPionSignal");
	hCutsFlow->SetBinContent(1,tmpHist->GetEntries());
	hCutsFlow->GetXaxis()->SetBinLabel(1, "excl. events");
	
	for (int i = 0; i < size; ++i)
	{
		cuts += cutsOption[i]; 

		tree->Draw("invMassPion>>invMassPionSignal",cuts);
		TH1F *tmpHist = (TH1F*)gPad->GetPrimitive("invMassPionSignal");
		hCutsFlow->GetXaxis()->SetBinLabel(i+2, cutsLabels[i]);
		hCutsFlow->SetBinContent(i+2,tmpHist->GetEntries());

		if(i != size -1)
			cuts += " && ";
	}

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
	textPub -> AddText("p + p #rightarrow p' + X + p'");
	textPub -> AddText("#sqrt{s} = 510 GeV");
	textPub -> AddText("Cuts flow");
	textPub -> Draw("same");

	newCanvas->Update();
	//newCanvas->SaveAs(output + "BasicPlots/Cuts.png");
	newCanvas->Write("CutsFlow");

	cutsDir->mkdir("PID")->cd();
	PID PIDPlotsWithCuts(data, fout, output, showCutsLine, cuts);
	PIDPlotsWithCuts.PlotHistogram();

    cutsDir->mkdir("trackQuality")->cd();
    trackQuality TrackPlotsWithCuts(data, fout, output, showCutsLine, cuts);
    TrackPlotsWithCuts.PlotHistogram();

///////////////////////////////////////////////////////////////
    TDirectory* TomasDir = cutsDir->mkdir("PID_Tomas");
    TomasDir->cd();
    TomasDir->mkdir("Pions")->cd();
    PID PIDPlotsWithCuts1(data, fout, output, showCutsLine, cuts + "&& nSigTrk1Pion < 3 && nSigTrk1Pion > -3 && nSigTrk2Pion > -3 && nSigTrk2Pion < 3 && (nSigPairKaon > 3 || mSquared < 0.2 || mSquared > 0.32) && (nSigPairProton > 3 || mSquared < 0.7 || mSquared > 1.1)");
    PIDPlotsWithCuts1.PlotHistogram();

    TomasDir->mkdir("Kaons")->cd();
    PID PIDPlotsWithCuts2(data, fout, output, showCutsLine, cuts + "&& nSigPairKaon < 3 && mSquared > 0.2 && mSquared < 0.32 && (nSigPairProton > 3 || mSquared < 0.7 || mSquared > 1.1)");
    PIDPlotsWithCuts2.PlotHistogram();

    TomasDir->mkdir("Protons")->cd();
    PID PIDPlotsWithCuts3(data, fout, output, showCutsLine, cuts + "&& nSigPairProton < 3 && mSquared > 0.7 && mSquared < 1.1");
    PIDPlotsWithCuts3.PlotHistogram();
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
    TDirectory* TomasDir2 = cutsDir->mkdir("PID_Tomas2");
    TomasDir2->cd();
    TomasDir2->mkdir("Pions")->cd();
    PID PIDPlotsWithCuts4(data, fout, output, showCutsLine, cuts + "&& nSigPairPion < 3.46 && (nSigPairPion < 3 || nSigPairKaon > 3 || mSquared < 0.15) && (nSigPairPion < 3 || nSigPairKaon < 3 || nSigPairProton > 3 || mSquared < 0.6)");
    PIDPlotsWithCuts4.PlotHistogram();

    TomasDir2->mkdir("Kaons")->cd();
    PID PIDPlotsWithCuts5(data, fout, output, showCutsLine, cuts + "&& nSigPairPion > 3 && nSigPairKaon < 3 && mSquared > 0.15 && (nSigPairPion < 3 || nSigPairKaon < 3 || nSigPairProton > 3 || mSquared < 0.6)");
    PIDPlotsWithCuts5.PlotHistogram();

    TomasDir2->mkdir("Protons")->cd();
    PID PIDPlotsWithCuts6(data, fout, output, showCutsLine, cuts + "&& nSigPairPion > 3 && nSigPairKaon > 3 && nSigPairProton < 3 && mSquared > 0.6");
    PIDPlotsWithCuts6.PlotHistogram();
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
    TDirectory* RafalDir = cutsDir->mkdir("PID_Rafal");
    RafalDir->cd();
    RafalDir->mkdir("Pions")->cd();
    PID PIDPlotsWithCuts7(data, fout, output, showCutsLine, cuts + "&& nSigTrk1Pion < 3 && nSigTrk1Pion > -3 && nSigTrk2Pion > -3 && nSigTrk2Pion < 3 && (nSigPairPion < 3 || nSigPairKaon > 3 || nSigPairProton < 3 || mSquared < 0.2 || mSquared > 0.32) && (nSigPairPion < 3 || nSigPairKaon < 3 || nSigPairProton > 3 || mSquared < 0.7 || mSquared > 1.1)");
    PIDPlotsWithCuts7.PlotHistogram();

    RafalDir->mkdir("Kaons")->cd();
    PID PIDPlotsWithCuts8(data, fout, output, showCutsLine, cuts + "&& nSigPairPion > 3 && nSigPairKaon < 3 && nSigPairProton > 3 && mSquared > 0.2 && mSquared < 0.32 && (nSigPairPion < 3 || nSigPairKaon < 3 || nSigPairProton > 3 || mSquared < 0.7 || mSquared > 1.1)");
    PIDPlotsWithCuts8.PlotHistogram();

    RafalDir->mkdir("Protons")->cd();
    PID PIDPlotsWithCuts9(data, fout, output, showCutsLine, cuts + "&& nSigPairPion > 3 && nSigPairKaon > 3 && nSigPairProton < 3 && mSquared > 0.6 && mSquared < 1.1");
    PIDPlotsWithCuts9.PlotHistogram();
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
    TDirectory* DanielDir = cutsDir->mkdir("PID_Daniel");
    DanielDir->cd();
    DanielDir->mkdir("Pions")->cd();
    PID PIDPlotsWithCuts10(data, fout, output, showCutsLine, cuts + "&& nSigTrk1Pion < 3 && nSigTrk1Pion > -3 && nSigTrk2Pion > -3 && nSigTrk2Pion < 3 && (nSigPairKaon > 3 || deltaDeltaTOFKaon < -0.5 || deltaDeltaTOFKaon > 0.5) && (nSigPairProton > 3 || deltaDeltaTOFProton < -0.5 || deltaDeltaTOFProton > 0.5)");
    PIDPlotsWithCuts10.PlotHistogram();

    DanielDir->mkdir("Kaons")->cd();
    PID PIDPlotsWithCuts11(data, fout, output, showCutsLine, cuts + "&& nSigPairKaon < 3 && deltaDeltaTOFKaon > -0.5 && deltaDeltaTOFKaon < 0.5  && (nSigPairProton > 3 || deltaDeltaTOFProton < -0.5 || deltaDeltaTOFProton > 0.5)");
    PIDPlotsWithCuts11.PlotHistogram();

    DanielDir->mkdir("Protons")->cd();
    PID PIDPlotsWithCuts12(data, fout, output, showCutsLine, cuts + "&& nSigPairProton < 3 && deltaDeltaTOFProton > -0.5 && deltaDeltaTOFProton < 0.5");
    PIDPlotsWithCuts12.PlotHistogram();
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

    gPad->SetLogy(0);
    TDirectory* uncorrDir = fout->mkdir("FinalPlots-uncorrected");
    uncorrDir->cd();

	hInvMassUncorr[Pion][0]->SetTitle(" ; m(#pi^{+}#pi^{-}) [GeV/c^{2}]; Number of events");
	tool.SetGraphStyle(hInvMassUncorr[Pion][0],4,20,1,4,1,1,0.9,1.3);
	tool.SetMarkerStyle(hInvMassUncorr[Pion][0]);
	hInvMassUncorr[Pion][0]->Draw("E");
	tool.DrawText(hInvMassUncorr[Pion][0], 1, true);
	tool.DrawTextStar(hInvMassUncorr[Pion][0]);
    tool.SetMarkerStyle(hInvMassUncorr[Pion][1],2,20,1,2,1,1);
    hInvMassUncorr[Pion][1]->Draw("ESAME");

    TLegend* leg1 = new TLegend(0.58, 0.7, 0.78, 0.8);
    tool.SetLegendStyle(leg1);
    leg1->AddEntry(hInvMassUncorr[Pion][0],"In+El (unlike-sign pairs)","p");
    leg1->AddEntry(hInvMassUncorr[Pion][1],"In+El (like-sign pairs)","p");
    leg1->Draw("same");

	newCanvas->Update();
	newCanvas->Write("uncorrInvMassPion");

	hInvMassUncorr[Kaon][0]->SetTitle(" ; m(K^{+}K^{-}) [GeV/c^{2}]; Number of events");
	tool.SetGraphStyle(hInvMassUncorr[Kaon][0],4,20,1,4,1,1,0.9,1.3);
	tool.SetMarkerStyle(hInvMassUncorr[Kaon][0]);
	hInvMassUncorr[Kaon][0]->Draw("E");
	tool.DrawText(hInvMassUncorr[Kaon][0], 2, true);
	tool.DrawTextStar(hInvMassUncorr[Kaon][0]);
    tool.SetMarkerStyle(hInvMassUncorr[Kaon][1],2,20,1,2,1,1);
    hInvMassUncorr[Kaon][1]->Draw("ESAME");

    leg1 = new TLegend(0.58, 0.7, 0.78, 0.8);
    tool.SetLegendStyle(leg1);
    leg1->AddEntry(hInvMassUncorr[Kaon][0],"In+El (unlike-sign pairs)","p");
    leg1->AddEntry(hInvMassUncorr[Kaon][1],"In+El (like-sign pairs)","p");
    leg1->Draw("same");

	newCanvas->Update();
	newCanvas->Write("uncorrInvMassKaon");

	hInvMassUncorr[Proton][0]->SetTitle(" ; m(p#bar{p}) [GeV/c^{2}]; Number of events");
	tool.SetGraphStyle(hInvMassUncorr[Proton][0],4,20,1,4,1,1,0.9,1.3);
	tool.SetMarkerStyle(hInvMassUncorr[Proton][0]);
	hInvMassUncorr[Proton][0]->Draw("E");
	tool.DrawText(hInvMassUncorr[Proton][0],3, true);
	tool.DrawTextStar(hInvMassUncorr[Proton][0]);
    tool.SetMarkerStyle(hInvMassUncorr[Proton][1],2,20,1,2,1,1);
    hInvMassUncorr[Proton][1]->Draw("ESAME");   

    leg1 = new TLegend(0.58, 0.7, 0.78, 0.8);
    tool.SetLegendStyle(leg1);
    leg1->AddEntry(hInvMassUncorr[Proton][0],"In+El (unlike-sign pairs)","p");
    leg1->AddEntry(hInvMassUncorr[Proton][1],"In+El (like-sign pairs)","p");
    leg1->Draw("same");

	newCanvas->Update();
	newCanvas->Write("uncorrInvMassProton");

    gPad->SetLogy(0);
    TDirectory* corrDir = fout->mkdir("FinalPlots-corrected");
    corrDir->cd();
 
    hInvMassCorr[Pion][0]->SetTitle(" ; m(#pi^{+}#pi^{-}) [GeV/c^{2}]; Corrected counts");
    tool.SetGraphStyle(hInvMassCorr[Pion][0],4,20,1,4,1,1,0.9,1.3);
    tool.SetMarkerStyle(hInvMassCorr[Pion][0]);
    hInvMassCorr[Pion][0]->Draw("E");
    tool.DrawText(hInvMassCorr[Pion][0], 1, true);
    tool.DrawTextStar(hInvMassCorr[Pion][0]);
    tool.SetMarkerStyle(hInvMassCorr[Pion][1],2,20,1,2,1,1);
    hInvMassCorr[Pion][1]->Draw("ESAME");

    leg1 = new TLegend(0.58, 0.7, 0.78, 0.8);
    tool.SetLegendStyle(leg1);
    leg1->AddEntry(hInvMassCorr[Pion][0],"In+El (unlike-sign pairs)","p");
    leg1->AddEntry(hInvMassCorr[Pion][1],"In+El (like-sign pairs)","p");
    leg1->Draw("same");

    newCanvas->Update();
    newCanvas->Write("corrInvMassPion");

    hInvMassCorr[Kaon][0]->SetTitle(" ; m(K^{+}K^{-}) [GeV/c^{2}]; Corrected counts");
    tool.SetGraphStyle(hInvMassCorr[Kaon][0],4,20,1,4,1,1,0.9,1.3);
    tool.SetMarkerStyle(hInvMassCorr[Kaon][0]);
    hInvMassCorr[Kaon][0]->Draw("E");
    tool.DrawText(hInvMassCorr[Kaon][0], 2, true);
    tool.DrawTextStar(hInvMassCorr[Kaon][0]);
    tool.SetMarkerStyle(hInvMassCorr[Kaon][1],2,20,1,2,1,1);
    hInvMassCorr[Kaon][1]->Draw("ESAME");

    leg1 = new TLegend(0.58, 0.7, 0.78, 0.8);
    tool.SetLegendStyle(leg1);
    leg1->AddEntry(hInvMassCorr[Kaon][0],"In+El (unlike-sign pairs)","p");
    leg1->AddEntry(hInvMassCorr[Kaon][1],"In+El (like-sign pairs)","p");
    leg1->Draw("same");

    newCanvas->Update();
    newCanvas->Write("corrInvMassKaon");

    hInvMassCorr[Proton][0]->SetTitle(" ; m(p#bar{p}) [GeV/c^{2}]; Corrected counts");
    tool.SetGraphStyle(hInvMassCorr[Proton][0],4,20,1,4,1,1,0.9,1.3);
    tool.SetMarkerStyle(hInvMassCorr[Proton][0]);
    hInvMassCorr[Proton][0]->Draw("E");
    tool.DrawText(hInvMassCorr[Proton][0],3, true);
    tool.DrawTextStar(hInvMassCorr[Proton][0]);
    tool.SetMarkerStyle(hInvMassCorr[Proton][1],2,20,1,2,1,1);
    hInvMassCorr[Proton][1]->Draw("ESAME");   

    leg1 = new TLegend(0.58, 0.7, 0.78, 0.8);
    tool.SetLegendStyle(leg1);
    leg1->AddEntry(hInvMassCorr[Proton][0],"In+El (unlike-sign pairs)","p");
    leg1->AddEntry(hInvMassCorr[Proton][1],"In+El (like-sign pairs)","p");
    leg1->Draw("same");

    newCanvas->Update();
    newCanvas->Write("corrInvMassProton");
	//hCutsSum->Draw();
//	textCut->Draw("same");
//	newCanvas->Update();
//	newCanvas->Write("CutsSummary");

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
    
    hInvMassCorr[0][0]  = new TH1D("corrInvMass" + particleLables[0] + "Sig", "Corrected inv. mass " + particleLables[0] , 64, 0.3, 3.5);
    hInvMassCorr[1][0]  = new TH1D("corrInvMass" + particleLables[1] + "Sig", "Corrected inv. mass " + particleLables[1] , 44, 0.8, 3);
    hInvMassCorr[2][0]  = new TH1D("corrInvMass" + particleLables[2] + "Sig", "Corrected inv. mass " + particleLables[2] , 24, 1.6, 4);
    hInvMassCorr[0][1]  = new TH1D("corrInvMass" + particleLables[0] + "Bcg", "Corrected inv. mass " + particleLables[0] , 64, 0.3, 3.5);
    hInvMassCorr[1][1]  = new TH1D("corrInvMass" + particleLables[1] + "Bcg", "Corrected inv. mass " + particleLables[1] , 44, 0.8, 3);
    hInvMassCorr[2][1]  = new TH1D("corrInvMass" + particleLables[2] + "Bcg", "Corrected inv. mass " + particleLables[2] , 24, 1.6, 4);

    hInvMassUncorr[0][0]  = new TH1D("uncorrInvMass" + particleLables[0] + "Sig", "Uncorrected inv. mass " + particleLables[0] , 64, 0.3, 3.5);
    hInvMassUncorr[1][0]  = new TH1D("uncorrInvMass" + particleLables[1] + "Sig", "Uncorrected inv. mass " + particleLables[1] , 44, 0.8, 3);
    hInvMassUncorr[2][0]  = new TH1D("uncorrInvMass" + particleLables[2] + "Sig", "Uncorrected inv. mass " + particleLables[2] , 24, 1.6, 4);
    hInvMassUncorr[0][1]  = new TH1D("uncorrInvMass" + particleLables[0] + "Bcg", "Uncorrected inv. mass " + particleLables[0] , 64, 0.3, 3.5);
    hInvMassUncorr[1][1]  = new TH1D("uncorrInvMass" + particleLables[1] + "Bcg", "Uncorrected inv. mass " + particleLables[1] , 44, 0.8, 3);
    hInvMassUncorr[2][1]  = new TH1D("uncorrInvMass" + particleLables[2] + "Bcg", "Uncorrected inv. mass " + particleLables[2] , 24, 1.6, 4);
}


void ConnectInput(TTree* tree)
{

// PID and some quality event info
    tree->SetBranchAddress("mSquared", &mSquared); 
    tree->SetBranchAddress("nSigTrk1Pion", &nSigmaTPC[Pion][0]);
    tree->SetBranchAddress("nSigTrk2Pion", &nSigmaTPC[Pion][1]);
    for (int iPart = 0; iPart < nParticles; ++iPart)
    {
        tree->SetBranchAddress("invMass" + particleLables[iPart], &invMass[iPart]);
        tree->SetBranchAddress("nSigPair" + particleLables[iPart], &nSigPair[iPart]);
    }


// Vertex info
    tree->SetBranchAddress("vertexZ", &vertexesZ[0]);

// Central track info
    for (int i = 0; i < 4; ++i)
    {
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

        effTotal = 1;
        if(nSigPair[Pion] > 3 && nSigPair[Kaon] > 3 && nSigPair[Proton] < 3 && mSquared > 0.6) // it is... proton!
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
                hInvMassCorr[Proton][signal]->Fill(invMass[Proton], 1/effTotal);
            if(transMomentum[0] > 0.4 && transMomentum[1] > 0.4 && (transMomentum[0] < 1.1 || transMomentum[1] < 1.1) )
                hInvMassUncorr[Proton][signal]->Fill(invMass[Proton]);
        }
        else if(nSigPair[Pion] > 3 && nSigPair[Kaon] < 3 && nSigPair[Proton] > 3 && mSquared > 0.15) // it is... kaon!
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
                hInvMassCorr[Kaon][signal]->Fill(invMass[Kaon], 1/effTotal);
            if(transMomentum[0] > 0.3 && transMomentum[1] > 0.3 && (transMomentum[0] < 0.7 || transMomentum[1] < 0.7) )
                hInvMassUncorr[Kaon][signal]->Fill(invMass[Kaon]);
        }
        else if( nSigPair[Pion] < sqrt(12)) // it is... pion!
        {
            for (int iTrack = 0; iTrack < 2; ++iTrack)
            {
                PID = 0;
                if(charge[iTrack] > 0)
                    PID = 3;
                effTPC = hTPCeff[0 + PID]->GetBinContent( hTPCeff[0 + PID]->GetXaxis()->FindBin(vertexesZ[iTrack]), hTPCeff[0 + PID]->GetYaxis()->FindBin(transMomentum[iTrack]), hTPCeff[0 + PID]->GetZaxis()->FindBin(Eta[iTrack]));
                effTOF = hTOFeff[0 + PID]->GetBinContent( hTOFeff[0 + PID]->GetXaxis()->FindBin(vertexesZ[iTrack]), hTOFeff[0 + PID]->GetYaxis()->FindBin(transMomentum[iTrack]), hTOFeff[0 + PID]->GetZaxis()->FindBin(Eta[iTrack])); 
                effTotal = effTotal*effTPC*effTOF;
            }
            if(effTotal != 0 && transMomentum[0] > 0.2 && transMomentum[1] > 0.2)
                hInvMassCorr[Pion][signal]->Fill(invMass[Pion], 1/effTotal);
            if(transMomentum[0] > 0.2 && transMomentum[1] > 0.2)
                hInvMassUncorr[Pion][signal]->Fill(invMass[Pion]);
        }

    }

}