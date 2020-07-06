#include "BasicPlots.h"
#include "Plot.h"


//_____________________________________________________________________________
BasicPlots::BasicPlots(TFile* dataInput, TFile* fileOut, TString outnam, bool text, TString inputCuts)
{
	//constructor
	TEXT = text;
	bcground = true;

	output = outnam;
	cuts = inputCuts;
	if(inputCuts != "")
		cutsWithPrefix = " && " + inputCuts;
	else
		cutsWithPrefix="";

	data = dataInput;
	fout = fileOut;

	siz = 0.045;
	font = 42;

	cout << "BasicPlots::BasicPlots() called" << endl;

}//BasicPlots

//_____________________________________________________________________________
BasicPlots::~BasicPlots()
{
  //destructor

  cout << "BasicPlots::~BasicPlots() destructor called" << endl;

}//~BasicPlots


//_____________________________________________________________________________
void BasicPlots::PlotHistogram() {

	if (!data){
		cout<<"Error: cannot open input file"<<endl;
		return;
    }

	TTree* tree = dynamic_cast<TTree*>( data->Get("recTree") );
	TTree* treeBack = dynamic_cast<TTree*>( data->Get("Background") );
	
	if (!tree || !treeBack){
		cout<<"Error: cannot open one of the TTree"<<endl;
		return;
    }

	if (!fout){
		cout<<"Error: cannot open output file"<<endl;
		return;
    }

	TCanvas *cCanvas = new TCanvas("cCanvas","cCanvas",800,700);
	gPad->SetMargin(0.11,0.02,0.105,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
	gPad->SetTickx();
	gPad->SetTicky();  
	gStyle->SetOptStat("");
	gStyle->SetPalette(1);
	gStyle->SetLineWidth(2);      //axis line
	gStyle->SetFrameLineWidth(2); //frame line
	//gStyle->SetErrorX(0); // suppres X errors

	Plot tool;
//////////////////////////////////////////

	hCuts = (TH1F*)data -> Get("AnalysisFlow");

// //////////////////////////////////////////////////////////
// Plot Cuts Flow
	TString Labels[] = { TString("All"), TString("2 TPC-TOF tracks"), 
	                  TString("Same vertex"), TString("TotCharge 0"), TString(""),};
	//gPad->SetMargin(0.9,0.02,0.1,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
	gPad->SetLogy();

	TH1F* hCutsFlow = new TH1F("cuts", "cuts", 8, 1, 9);
	int i = 1;
	for(int tb=1; tb<5; ++tb){ 
		hCutsFlow->GetXaxis()->SetBinLabel(tb, Labels[tb-1]);
		hCutsFlow->SetBinContent(tb,hCuts->GetBinContent(i));
		i++;
	}
//	hCutsFlow->GetXaxis()->SetBinLabel(8, "4 #pi");
//	tree->Draw("nSigPPion>>nSigPPion","fourPiState");
//	hCutsFlow->SetBinContent(8,((TH1F*)gPad->GetPrimitive("nSigPPion"))->GetEntries());
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

	
	cCanvas->Update();
	//cCanvas->SaveAs(output + "BasicPlots/Cuts.png");
	cCanvas->Write("CutsFlow");
  
}//BasicPlots::PlotHistogram

