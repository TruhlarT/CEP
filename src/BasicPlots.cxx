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
	hMissingPtTPC = (TH1D*)data -> Get("All/MissingPt_TPC2t_Combi2part");
	hMissingPtTOF = (TH1D*)data -> Get("All/MissingPt_TOF2trk_Combi2part");
	hMissingPtQ0 = (TH1D*)data -> Get("All/MissingPt_Q0_Combi2part");
	hMissingPtExc = (TH1D*)data -> Get("All/MissingPt_Excl_Combi2part");

// //////////////////////////////////////////////////////////
// Plot Cuts Flow
	TString Labels[] = { TString("All"), TString("CPT2noBBCL"), TString("El+Inel"), TString("2 TPC-TOF tracks"), 
	                  TString("Same vertex"), TString("TotCharge 0"), TString("p_{T}^{miss} < 0.1 GeV/c"), TString(""),};
	//gPad->SetMargin(0.9,0.02,0.1,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
	gPad->SetLogy();

	TH1F* hCutsFlow = new TH1F("cuts", "cuts", 8, 1, 9);
	int i = 1;
	for(int tb=1; tb<8; ++tb){ 
		if(tb==3){
			hCutsFlow->GetXaxis()->SetBinLabel(tb, Labels[tb-1]);
			hCutsFlow->SetBinContent(tb,hCuts->GetBinContent(i)+hCuts->GetBinContent((++i)++));
			continue;
		}
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
//////////////////////////////////////////
//////////////////////////////////////////
////////////// 
	gPad->SetLogy();
	hMissingPtTOF->SetStats(0);
	hMissingPtTOF->SetTitle(" ; p_{T}^{miss} [GeV];Number of events");
	tool.SetGraphStyle(hMissingPtTOF);
	//tool.SetMarkerStyle(hMissingPtTPC);
	hMissingPtTPC->SetMinimum(10);

	//hMissingPtTPC->Draw();	
    hMissingPtTOF->SetFillColorAlpha(7, 0.2);
    hMissingPtTOF->SetLineColor(28);
    hMissingPtTOF->SetFillStyle(1001);
    hMissingPtTOF->Draw(); 
	hMissingPtQ0->SetFillColor(4);
	hMissingPtQ0->SetLineColor(4);
	hMissingPtQ0->SetFillStyle(1001);
	hMissingPtQ0->Draw("same");
	hMissingPtExc->SetFillColorAlpha(2, 0.3);
	hMissingPtExc->SetLineColor(2);
	hMissingPtExc->SetFillStyle(1001);
	hMissingPtExc->Draw("same");
	cCanvas->cd();


	leg1 = new TLegend(0.66,0.73,0.95,0.9);
	tool.SetLegendStyle(leg1);
	leg1->SetMargin(0.1);
	//leg1 -> AddEntry(hMissingPtTPC, "2 TPC tracks", "l");
	leg1 -> AddEntry(hMissingPtTOF, "2 TPC-TOF tracks", "fl");
	leg1 -> AddEntry(hMissingPtQ0, "Total charge 0", "fl");
	leg1 -> AddEntry(hMissingPtExc, "Exclusive", "fl");
	leg1->Draw("same");

	TPaveText *textPub1 = new TPaveText(0.25,0.75,0.35,0.75,"brNDC");
	tool.SetTextStyle(textPub1);
	if(TEXT)
		textPub1 -> AddText("Exclusive peak");
	textPub1 -> Draw("same");

    textPub = new TPaveText(0.34,0.9,0.95,0.95,"brNDC");
    textPub -> SetTextSize(0.04);
    textPub -> SetTextAlign(12);
    textPub -> SetFillColor(0);
    textPub -> SetTextFont(42);
    textPub -> AddText("p + p #rightarrow p + #pi^{+}#pi^{+}#pi^{-}#pi^{-} + p      #sqrt{s} = 510 GeV");
    //textPub -> AddText("#sqrt{s} = 510 GeV");
    textPub -> Draw("same");

    TPaveText *textSTAR = new TPaveText(0.13,0.9,0.25,0.95,"brNDC");
    textSTAR -> SetTextSize(0.04);
    textSTAR -> SetTextAlign(12);
    textSTAR -> SetFillColor(0);
    textSTAR -> SetTextFont(62);
    textSTAR -> AddText("THIS THESIS");
    textSTAR -> Draw("same");

	TLine *left = new TLine(0.1,0,0.1,hMissingPtTOF->GetMaximum()*1.1);
	tool.SetLineStyle(left,10,1,4);
    left->Draw("same");

	cCanvas->Update();
	cCanvas-> Write("hMissingPt");
	//cCanvas->SaveAs(output + "BasicPlots/cMissingPt.png");
	cCanvas->Close();

  
}//BasicPlots::PlotHistogram



/*
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

*/