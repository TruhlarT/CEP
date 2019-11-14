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

#include "BasicPlots.h"
#include "PID.h"
#include "trackQuality.h"
#include "Plot.h"
#include "FourPi.h"

using namespace std;

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
	

	TString cutsOption[] = { TString("vertexZ<80 && vertexZ > -80"), TString("NhitsFit1 >=25 && NhitsFit0 >= 25"), TString("NhitsDEdx1 >= 15 && NhitsDEdx0 >= 15"), TString("DcaZ1 < 1 && DcaZ1 > -1 && DcaZ0 < 1 && DcaZ0 > -1"), TString("DcaXY1 < 1.5 && DcaXY0 < 1.5"), TString("Eta1 > -0.8 && Eta1 < 0.8 && Eta0 > -0.8 && Eta0 < 0.8"), TString("!fourPiState")};
	TString cutsLabels[] = { TString("|z_{vtx}| < 80 cm"), TString("N_{hits}^{fit} #geq 25"), TString("N_{hits}^{dE/dx} #geq 15"), TString("|DCA(z)| < 1 cm"), TString("DCA(XY) < 1.5 cm"), TString("|#eta| < 0.8"), TString("!fourPiState")  };

	TFile* data = TFile::Open(dataName, "read");
	if (!data){
		cout<<"Error: cannot open: "<<dataName<<endl;
		return 2;
    }

	TTree* tree = dynamic_cast<TTree*>( data->Get("recTree") );
	TTree* treeBack = dynamic_cast<TTree*>( data->Get("Background") );
	
	if (!tree || !treeBack){
		cout<<"Error: cannot open one of the TTree"<<endl;
		return 3;
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
    PID PIDPlotsWithCuts4(data, fout, output, showCutsLine, cuts + "&& nSigTrk1Pion < 3 && nSigTrk1Pion > -3 && nSigTrk2Pion > -3 && nSigTrk2Pion < 3 && (nSigPairKaon > 3 || ((mSquared < 0.2 || mSquared > 0.32) && (deltaDeltaTOFKaon < -0.5 || deltaDeltaTOFKaon > 0.5))) && (nSigPairProton > 3 || ((mSquared < 0.7 || mSquared > 1.1) && (deltaDeltaTOFProton < -0.5 || deltaDeltaTOFProton > 0.5)))");
    PIDPlotsWithCuts4.PlotHistogram();

    TomasDir2->mkdir("Kaons")->cd();
    PID PIDPlotsWithCuts5(data, fout, output, showCutsLine, cuts + "&& nSigPairKaon < 3 && ((mSquared > 0.2 && mSquared < 0.32) || (deltaDeltaTOFKaon > -0.5 && deltaDeltaTOFKaon < 0.5)) && (nSigPairProton > 3 || ((mSquared < 0.7 || mSquared > 1.1) && (deltaDeltaTOFProton < -0.5 || deltaDeltaTOFProton > 0.5)))");
    PIDPlotsWithCuts5.PlotHistogram();

    TomasDir2->mkdir("Protons")->cd();
    PID PIDPlotsWithCuts6(data, fout, output, showCutsLine, cuts + "&& nSigPairProton < 3 && ((mSquared > 0.7 && mSquared < 1.1) || (deltaDeltaTOFProton > -0.5 && deltaDeltaTOFProton < 0.5))");
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
    PID PIDPlotsWithCuts9(data, fout, output, showCutsLine, cuts + "&& nSigPairPion > 3 && nSigPairKaon > 3 && nSigPairProton < 3 && mSquared > 0.7 && mSquared < 1.1");
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
	gPad->SetLogy(0);
	TDirectory* PIDDir = fout->mkdir("AppliedPID");
	PIDDir->cd();

	TH1F* hCutsSum = new TH1F("cutsSummary", "cuts", size+1, 0, size);
	TPaveText *textCut = new TPaveText(0.1,0.1,0.95,0.9,"brNDC");
	tool.SetTextStyle(textCut);

	TString appliedCuts = "nSigTrk1Pion < 3 && nSigTrk1Pion > -3 && nSigTrk2Pion > -3 && nSigTrk2Pion < 3 && (nSigPairKaon > 3 || mSquared < 0.2 || mSquared > 0.32) && (nSigPairProton > 3 || mSquared < 0.7 || mSquared > 1.1) && " + cuts;
//	textCut -> AddText("All track quality cuts applied");
//	textCut -> AddText("#pi #pi: (nSigPairProton >= 3 || nSigPairKaon <= 3 || nSigPairPion <= 3)");
//	textCut -> AddText("&& (nSigPairKaon >= 3 || nSigPairProton >= 3 || nSigPairPion <= 3)");

	treeBack->Draw("invMassPion>>invMassPionBackground(50,0,2.5)",appliedCuts);
	TH1F* tmpHist2 = (TH1F*)gPad->GetPrimitive("invMassPionBackground");
	tool.SetMarkerStyle(tmpHist2,2,20,1,2,1,1);

	tree->Draw("invMassPion>>invMassPionSignal(50,0,2.5)",appliedCuts);
	tmpHist = (TH1F*)gPad->GetPrimitive("invMassPionSignal");
	tmpHist->SetTitle(" ; m(#pi^{+}#pi^{-}) [GeV/c^{2}]; Number of events");
	//tmpHist->GetXaxis()->SetRangeUser(0,2.5);
	tool.SetGraphStyle(tmpHist,4,20,1,4,1,1,0.9,1.3);
	tool.SetMarkerStyle(tmpHist);
	tmpHist->Draw("E");
	tool.DrawText(tmpHist, 1, true);
	tool.DrawTextStar(tmpHist);
	tmpHist2->Draw("ESAME");

	TLegend* leg1 = new TLegend(0.58, 0.7, 0.78, 0.8);
	tool.SetLegendStyle(leg1);
	leg1->AddEntry(tmpHist,"In+El (unlike-sign pairs)","p");
	leg1->AddEntry(tmpHist2,"In+El (like-sign pairs)","p");
	leg1->Draw("same");

	newCanvas->Update();
	newCanvas->Write("invMassPion");

	appliedCuts = "nSigPairKaon < 3 && mSquared > 0.2 && mSquared < 0.32 && (nSigPairProton > 3 || mSquared < 0.7 || mSquared > 1.1) && " + cuts;
	textCut -> AddText("KK: nSigPairKaon < 3 && nSigPairProton > 3 && nSigPairPion > 3");

	treeBack->Draw("invMassKaon>>invMassKaonBackground(50,0.5,3.5)",appliedCuts);
	tmpHist2 = (TH1F*)gPad->GetPrimitive("invMassKaonBackground");
	tool.SetMarkerStyle(tmpHist2,2,20,1,2,1,1);

	tree->Draw("invMassKaon>>invMassKaonSignal(50,0.5,3.5)",appliedCuts);
	tmpHist = (TH1F*)gPad->GetPrimitive("invMassKaonSignal");
	tmpHist->SetTitle(" ; m(K^{+}K^{-}) [GeV/c^{2}]; Number of events");
	//tmpHist->GetXaxis()->SetRangeUser(0,2.5);
	tool.SetGraphStyle(tmpHist,4,20,1,4,1,1,0.9,1.3);
	tool.SetMarkerStyle(tmpHist);
	tmpHist->Draw("E");
	tool.DrawText(tmpHist, 2, true);
	tool.DrawTextStar(tmpHist);
	tmpHist2->Draw("ESAME");

	leg1 = new TLegend(0.58, 0.7, 0.78, 0.8);
	tool.SetLegendStyle(leg1);
	leg1->AddEntry(tmpHist,"In+El (unlike-sign pairs)","p");
	leg1->AddEntry(tmpHist2,"In+El (like-sign pairs)","p");
	leg1->Draw("same");

	newCanvas->Update();
	newCanvas->Write("invMassKaon");

	appliedCuts = "nSigPairProton < 3 && mSquared > 0.7 && mSquared < 1.13 && " + cuts;
	textCut -> AddText("pp: nSigPairProton < 3 && nSigPairKaon > 3 && nSigPairPion > 3");

	treeBack->Draw("invMassProton>>invMassProtonBackground(50,1.0,4.5)",appliedCuts);
	tmpHist2 = (TH1F*)gPad->GetPrimitive("invMassProtonBackground");
	tool.SetMarkerStyle(tmpHist2,2,20,1,2,1,1);

	tree->Draw("invMassProton>>invMassProtonSignal(50,1.0,4.5)",appliedCuts);
	tmpHist = (TH1F*)gPad->GetPrimitive("invMassProtonSignal");
	tmpHist->SetTitle(" ; m(p#bar{p}) [GeV/c^{2}]; Number of events");
	//tmpHist->GetXaxis()->SetRangeUser(0,2.5);
	tool.SetGraphStyle(tmpHist,4,20,1,4,1,1,0.9,1.3);
	tool.SetMarkerStyle(tmpHist);
	tmpHist->Draw("E");
	tool.DrawText(tmpHist,3, true);
	tool.DrawTextStar(tmpHist);
	tmpHist2->Draw("ESAME");

	leg1 = new TLegend(0.58, 0.7, 0.78, 0.8);
	tool.SetLegendStyle(leg1);
	leg1->AddEntry(tmpHist,"In+El (unlike-sign pairs)","p");
	leg1->AddEntry(tmpHist2,"In+El (like-sign pairs)","p");
	leg1->Draw("same");

	newCanvas->Update();
	newCanvas->Write("invMassProton");

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

