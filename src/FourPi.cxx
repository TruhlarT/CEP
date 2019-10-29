#include "FourPi.h"
#include "Plot.h"


//_____________________________________________________________________________
FourPi::FourPi(TFile* dataInput, TFile* fileOut, TString outnam, bool text, TString inputCuts)
{
	//constructor
	output = outnam;

	if(inputCuts != "")
	{
		cuts = inputCuts + " && fourPiState";
	}
	else
	{
		cuts = "fourPiState";
	}



	data = dataInput;
	fout = fileOut;

	TEXT = text;

	cout << "trackQuality::trackQuality() called" << endl;
}//FourPi

//_____________________________________________________________________________
FourPi::~FourPi()
{
  //destructor

  cout << "FourPi::~FourPi() destructor called" << endl;

}//~FourPi


//_____________________________________________________________________________
void FourPi::PlotHistogram() {

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
	gPad->SetMargin(0.12,0.02,0.1,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
	gPad->SetTickx();
	gPad->SetTicky();  
	gStyle->SetOptStat("");
	gStyle->SetPalette(1);
	gStyle->SetLineWidth(2);      //axis line
	gStyle->SetFrameLineWidth(2); //frame line

	tree->UseCurrentStyle();
	treeBack->UseCurrentStyle();

	Plot tool;
	TString variable, cutsOption;
	Int_t nBins, nInputs;
	Float_t min, max;
//////////////////////////////////////////
	// Plot Cuts Flow for 4 pion
	tmpHist = (TH1F*)data -> Get("AnalysisFlow2");
	TString Labels2[] = { TString("All"), TString("CPT trigger"), TString("El+InEl"), TString("4 TPC-TOF tracks"), 
	                  TString("Same vertex"), TString("TotCharge 0"), TString("p_{T}^{miss} < 0.1 GeV/c"), TString(""),};
	gPad->SetLogy();

	TH1F* hCutsFlow2 = new TH1F("cuts2", "cuts", 8, 1, 9);
	int i = 1;
	for(int tb=1; tb<8; ++tb){ 
		if(tb==3){
			hCutsFlow2->GetXaxis()->SetBinLabel(tb, Labels2[tb-1]);
			hCutsFlow2->SetBinContent(tb,tmpHist->GetBinContent(i)+tmpHist->GetBinContent((++i)++));
			continue;
		}
		hCutsFlow2->GetXaxis()->SetBinLabel(tb, Labels2[tb-1]);
		hCutsFlow2->SetBinContent(tb,tmpHist->GetBinContent(i));
		i++;
	}
	hCutsFlow2->SetTitle("; ; Number of events");
	tool.SetGraphStyle(hCutsFlow2);
	tool.SetMarkerStyle(hCutsFlow2);
	//hCutsFlow->GetYaxis()->SetRangeUser(800,300000000);
	hCutsFlow2->GetXaxis()->SetLabelSize(0.045);
	hCutsFlow2->Draw();
	tool.DrawTextStar(hCutsFlow2,2);

	TPaveText* textPub = new TPaveText(0.75,0.74,0.9,0.9,"brNDC");
	tool.SetTextStyle(textPub);
	textPub -> AddText("p + p #rightarrow p' + #pi^{+}#pi^{+}#pi^{-}#pi^{-} + p'");
	textPub -> AddText("#sqrt{s} = 510 GeV");
	textPub -> AddText("Cuts flow");
	textPub -> Draw("same");

	
	cCanvas->Update();
	cCanvas->Write("CutsFlow");
	gPad->SetLogy(0);
//////////////////////////////////////////////////////////////
	TH1F* hMissingPtTPC = (TH1F*)data -> Get("All/MissingPt_TPC2t_Combi4part");
	TH1F* hMissingPtQ0 = (TH1F*)data -> Get("All/MissingPt_Q0_Combi4part");
	TH1F* hMissingPtExc = (TH1F*)data -> Get("All/MissingPt_Excl_Combi4part");
	gPad->SetLogy();
	hMissingPtTPC->SetStats(0);
	hMissingPtTPC->SetTitle(" ; p_{T}^{miss} [GeV/c];Number of events");
	tool.SetGraphStyle(hMissingPtTPC);
	hMissingPtTPC->SetMinimum(10);

	hMissingPtTPC->Draw();	
	hMissingPtQ0->SetFillColor(4);
	hMissingPtQ0->SetLineColor(4);
	hMissingPtQ0->SetFillStyle(1001);
	hMissingPtQ0->Draw("same");
	hMissingPtExc->SetFillColorAlpha(2, 0.5);
	hMissingPtExc->SetLineColor(2);
	hMissingPtExc->SetFillStyle(1001);
	hMissingPtExc->Draw("same");
	cCanvas->cd();


	leg1 = new TLegend(0.3,0.8,0.5,0.95);
	tool.SetLegendStyle(leg1);
	leg1 -> AddEntry(hMissingPtTPC, "4 TPC-TOF tracks", "l");
	leg1 -> AddEntry(hMissingPtQ0, "Total charge 0", "fl");
	leg1 -> AddEntry(hMissingPtExc, "Exclusive", "fl");
	leg1->Draw("same");

	TPaveText *textPub1 = new TPaveText(0.25,0.75,0.35,0.75,"brNDC");
	tool.SetTextStyle(textPub1);
	if(TEXT)
		textPub1 -> AddText("Exclusive peak");
	textPub1 -> Draw("same");

	TPaveText *textPub2 = new TPaveText(0.75,0.79,0.9,0.9,"brNDC");
	tool.SetTextStyle(textPub2);
	textPub2 -> AddText("p + p #rightarrow p' + #pi^{+}#pi^{+}#pi^{-}#pi^{-} + p'");
	textPub2 -> AddText("#sqrt{s} = 510 GeV");
	textPub2 -> Draw("same");
	tool.DrawTextStar(tmpHist,2);

	TLine *left = new TLine(0.1,0,0.1,600000);
	tool.SetLineStyle(left,10,1,4);
    left->Draw("same");

	cCanvas->Update();
	cCanvas-> Write("hMissingPt");
	gPad->SetLogy(0);
///////////////////////////////////////////////////////////////
///////////// Plot Inv Mass Inelastic + Elastic
	variable = "invMassPion";
	cutsOption = cuts;
	nBins = 50;
	min = 0.5;
	max = 4.5;

	treeBack->Draw(variable + ">>" + variable + "Bcg(" + nBins + "," + min + "," + max + ")",cutsOption);
	histBackground = (TH1F*)gPad->GetPrimitive(variable + "Bcg");
	tool.SetMarkerStyle(histBackground,2,20,1,2,1,1);

	tree->Draw(variable + ">>" + variable + "Sig(" + nBins + "," + min + "," + max + ")",cutsOption);
	histSignal = (TH1F*)gPad->GetPrimitive(variable + "Sig");
	histSignal->SetTitle(" ; m(#pi^{+}#pi^{+}#pi^{-}#pi^{-}) [GeV/c^{2}]; Number of events");
	tool.SetGraphStyle(histSignal,4,20,1,4,1,1,0.9,1.3);
	tool.SetMarkerStyle(histSignal);
	histSignal->Draw("E");
	tool.DrawText(histSignal,4,true);
	tool.DrawTextStar(histSignal);
	histBackground->Draw("ESAME");

	leg1 = new TLegend(0.58, 0.7, 0.78, 0.8);
	tool.SetLegendStyle(leg1);
	leg1->AddEntry(histSignal,"In+El (unlike-sign pairs)","p");
	leg1->AddEntry(histBackground,"In+El (like-sign pairs)","p");
	leg1->Draw("same");
	
	if(TEXT){
		TPaveText *textPub3 = new TPaveText(0.6,0.6,0.7,0.7,"brNDC");
		tool.SetTextStyle(textPub3);
		textPub3 -> AddText("Drop at about 1 GeV");
		textPub3 -> Draw("same");

		TPaveText *textPub4 = new TPaveText(0.75,0.5,0.8,0.6,"brNDC");
		tool.SetTextStyle(textPub4);
		textPub4 -> AddText("Peak at about 1270 MeV");
		textPub4 -> Draw("same");	
	}


	cCanvas->Update();
	cCanvas->Write(variable);
//////////////////////////////////////////
// Plot Inv Mass Inelastic

	variable = "invMassPion";
	cutsOption = cuts + "&& !elastic";
	nBins = 50;
	min = 0.5;
	max = 4.5;

	treeBack->Draw(variable + ">>" + variable + "Bcg(" + nBins + "," + min + "," + max + ")",cutsOption);
	histBackground = (TH1F*)gPad->GetPrimitive(variable + "Bcg");
	tool.SetMarkerStyle(histBackground,2,20,1,2,1,1);

	tree->Draw(variable + ">>" + variable + "Sig(" + nBins + "," + min + "," + max + ")",cutsOption);
	histSignal = (TH1F*)gPad->GetPrimitive(variable + "Sig");
	histSignal->SetTitle(" ; m(#pi^{+}#pi^{+}#pi^{-}#pi^{-} [GeV/c^{2}]; Number of events");
	tool.SetGraphStyle(histSignal,4,20,1,4,1,1,0.9,1.3);
	tool.SetMarkerStyle(histSignal);
	histSignal->Draw("E");
	tool.DrawText(histSignal,4,true);
	tool.DrawTextStar(histSignal);
	histBackground->Draw("ESAME");

	leg1 = new TLegend(0.58, 0.7, 0.78, 0.8);
	tool.SetLegendStyle(leg1);
	leg1->AddEntry(histSignal,"In (unlike-sign pairs)","p");
	leg1->AddEntry(histBackground,"In (like-sign pairs)","p");
	leg1->Draw("same");

	cCanvas->Update();
	cCanvas->Write(variable + "In");
//////////////////////////////////////////
// Plot Inv Mass Elastic
	variable = "invMassPion";
	cutsOption = cuts + "&& elastic";
	nBins = 50;
	min = 0.5;
	max = 4.5;

	treeBack->Draw(variable + ">>" + variable + "Bcg(" + nBins + "," + min + "," + max + ")",cutsOption);
	histBackground = (TH1F*)gPad->GetPrimitive(variable + "Bcg");
	tool.SetMarkerStyle(histBackground,2,20,1,2,1,1);

	tree->Draw(variable + ">>" + variable + "Sig(" + nBins + "," + min + "," + max + ")",cutsOption);
	histSignal = (TH1F*)gPad->GetPrimitive(variable + "Sig");
	histSignal->SetTitle(" ; m(#pi^{+}#pi^{+}#pi^{-}#pi^{-}) [GeV/c^{2}]; Number of events");
	tool.SetGraphStyle(histSignal,4,20,1,4,1,1,0.9,1.3);
	tool.SetMarkerStyle(histSignal);
	histSignal->Draw("E");
	tool.DrawText(histSignal,4,true);
	tool.DrawTextStar(histSignal);
	histBackground->Draw("ESAME");

	leg1 = new TLegend(0.58, 0.7, 0.78, 0.8);
	tool.SetLegendStyle(leg1);
	leg1->AddEntry(histSignal,"El (unlike-sign pairs)","p");
	leg1->AddEntry(histBackground,"El (like-sign pairs)","p");
	leg1->Draw("same");
	
	cCanvas->Update();
	cCanvas->Write(variable + "El");
//////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////	
// Plot z vertex
	variable = "vertexZ";
	cutsOption = cuts;
	nBins = 100;
	min = -200;
	max = 200;

	treeBack->Draw(variable + ">>" + variable + "Bcg(" + nBins + "," + min + "," + max + ")",cutsOption);
	histBackground = (TH1F*)gPad->GetPrimitive(variable + "Bcg");
	tool.SetMarkerStyle(histBackground,2,20,1,2,1,1);

	tree->Draw(variable + ">>" + variable + "Sig(" + nBins + "," + min + "," + max + ")",cutsOption);
	histSignal = (TH1F*)gPad->GetPrimitive(variable + "Sig");
	histSignal->SetTitle(" ; Z_{vrtx} [cm]; Number of events");
	tool.SetGraphStyle(histSignal,4,20,1,4,1,1,0.9,1.3);
	tool.SetMarkerStyle(histSignal);
	histSignal->Draw("E");
	tool.DrawText(histSignal,4,true,0.68, 0.75, 0.9, 0.88);
	tool.DrawTextStar(histSignal,2);
	histBackground->Draw("ESAME");

	leg1 = new TLegend(0.15,0.83,0.28,0.93); 
	tool.SetLegendStyle(leg1);
	leg1->AddEntry(histSignal,"In+El (unlike-sign pairs)","p");
	leg1->AddEntry(histBackground,"In+El (like-sign pairs)","p");
	leg1->Draw("same");
	if(TEXT){
		TLine *left4 = new TLine(-80,0,-80,140);
		tool.SetLineStyle(left4,10,1,4);
	    left4->Draw("same");

		TLine *left5 = new TLine(80,0,80,140);
		tool.SetLineStyle(left5,10,1,4);
	    left5->Draw("same");		
	}
	
	cCanvas->Update();
	cCanvas->Write(variable);
//////////////////////////////////////////
// Plot DCAXY vertex
	variable = "DcaXY";
	cutsOption = cuts;
	nInputs = 4;
	nBins = 100;
	min = 0;
	max = 3.5;

	treeBack->Draw(variable + "1>>" + variable + "Bcg1(" + nBins + "," + min + "," + max + ")",cutsOption);
	histBackground = (TH1F*)gPad->GetPrimitive(variable +"Bcg1");
	for(int i = 2; i <= nInputs; ++i)
	{	
		treeBack->Draw(variable + i + ">>" + variable + "Bcg" + i + "(" + nBins + "," + min + "," + max + ")",cutsOption);
		tmpHist = (TH1F*)gPad->GetPrimitive(variable +"Bcg" + i);
		histBackground->Add(tmpHist);
	}
	tool.SetMarkerStyle(histBackground,2,20,1,2,1,1);

	tree->Draw(variable + "1>>" + variable + "Sig1(" + nBins + "," + min + "," + max + ")",cutsOption);
	histSignal = (TH1F*)gPad->GetPrimitive(variable +"Sig1");
	for(int i = 2; i <= nInputs; ++i)
	{	
		tree->Draw(variable + i + ">>" + variable + "Sig" + i + "(" + nBins + "," + min + "," + max + ")",cutsOption);
		tmpHist = (TH1F*)gPad->GetPrimitive(variable +"Sig" + i);
		histSignal->Add(tmpHist);
	}
	histSignal->SetTitle(" ; DCA_{xy} [cm]; Number of tracks");
	tool.SetGraphStyle(histSignal,4,20,1,4,1,1,0.9,1.4);
	tool.SetMarkerStyle(histSignal);
	histSignal->Draw("E");
	tool.DrawText(histSignal,4,false,0.68, 0.75, 0.9, 0.88);
	tool.DrawTextStar(histSignal,2);
	histBackground->Draw("ESAME");

	leg1 = new TLegend(0.6, 0.65, 0.78, 0.74);
	tool.SetLegendStyle(leg1);
	leg1->AddEntry(histSignal,"In+El (unlike-sign pairs)","p");
	leg1->AddEntry(histBackground,"In+El (like-sign pairs)","p");
	leg1->Draw("same");
	if(TEXT){
		TLine *left6 = new TLine(1.5,0,1.5,700);
		tool.SetLineStyle(left6,10,1,4);
	    left6->Draw("same");
	}
	cCanvas->Update();
	cCanvas->Write(variable);
//////////////////////////////////////////
// Plot DCAZ vertex
	variable = "DcaZ";
	cutsOption = cuts;
	nInputs = 4;
	nBins = 100;
	min = -1.5;
	max = 1.5;

	treeBack->Draw(variable + "1>>" + variable + "Bcg1(" + nBins + "," + min + "," + max + ")",cutsOption);
	histBackground = (TH1F*)gPad->GetPrimitive(variable +"Bcg1");
	for(int i = 2; i <= nInputs; ++i)
	{	
		treeBack->Draw(variable + i + ">>" + variable + "Bcg" + i + "(" + nBins + "," + min + "," + max + ")",cutsOption);
		tmpHist = (TH1F*)gPad->GetPrimitive(variable +"Bcg" + i);
		histBackground->Add(tmpHist);
	}
	tool.SetMarkerStyle(histBackground,2,20,1,2,1,1);

	tree->Draw(variable + "1>>" + variable + "Sig1(" + nBins + "," + min + "," + max + ")",cutsOption);
	histSignal = (TH1F*)gPad->GetPrimitive(variable +"Sig1");
	for(int i = 2; i <= nInputs; ++i)
	{	
		tree->Draw(variable + i + ">>" + variable + "Sig" + i + "(" + nBins + "," + min + "," + max + ")",cutsOption);
		tmpHist = (TH1F*)gPad->GetPrimitive(variable +"Sig" + i);
		histSignal->Add(tmpHist);
	}
	histSignal->SetTitle(" ; DCA_{z} [cm]; Number of tracks");
	tool.SetGraphStyle(histSignal,4,20,1,4,1,1,0.9,1.4);
	tool.SetMarkerStyle(histSignal);
	histSignal->Draw("E");
	tool.DrawText(histSignal,4,false,0.68, 0.75, 0.9, 0.88);
	tool.DrawTextStar(histSignal,2);
	histBackground->Draw("ESAME");

	leg1 = new TLegend(0.6, 0.65, 0.78, 0.74);
	tool.SetLegendStyle(leg1);
	leg1->AddEntry(histSignal,"In+El (unlike-sign pairs)","p");
	leg1->AddEntry(histBackground,"In+El (like-sign pairs)","p");
	leg1->Draw("same");
	if(TEXT){
		TLine *left7 = new TLine(1,0,1,5000);
		tool.SetLineStyle(left7,10,1,4);
	    left7->Draw("same");

		TLine *left8 = new TLine(-1,0,-1,5000);
		tool.SetLineStyle(left8,10,1,4);
	    left8->Draw("same");
	}
	cCanvas->Update();
	cCanvas->Write(variable);
//////////////////////////////////////////
// Plot NhitsDEdx vertex
	variable = "NhitsDEdx";
	cutsOption = cuts;
	nInputs = 4;
	nBins = 61;
	min = 0;
	max = 60;

	treeBack->Draw(variable + "1>>" + variable + "Bcg1(" + nBins + "," + min + "," + max + ")",cutsOption);
	histBackground = (TH1F*)gPad->GetPrimitive(variable +"Bcg1");
	for(int i = 2; i <= nInputs; ++i)
	{	
		treeBack->Draw(variable + i + ">>" + variable + "Bcg" + i + "(" + nBins + "," + min + "," + max + ")",cutsOption);
		tmpHist = (TH1F*)gPad->GetPrimitive(variable +"Bcg" + i);
		histBackground->Add(tmpHist);
	}
	tool.SetMarkerStyle(histBackground,2,20,1,2,1,1);

	tree->Draw(variable + "1>>" + variable + "Sig1(" + nBins + "," + min + "," + max + ")",cutsOption);
	histSignal = (TH1F*)gPad->GetPrimitive(variable +"Sig1");
	for(int i = 2; i <= nInputs; ++i)
	{	
		tree->Draw(variable + i + ">>" + variable + "Sig" + i + "(" + nBins + "," + min + "," + max + ")",cutsOption);
		tmpHist = (TH1F*)gPad->GetPrimitive(variable +"Sig" + i);
		histSignal->Add(tmpHist);
	}
	histSignal->SetTitle(" ; N^{dE/dx}_{hits} [-]; Number of tracks");
	tool.SetGraphStyle(histSignal,4,20,1,4,1,1,0.9,1.4);
	tool.SetMarkerStyle(histSignal);
	histSignal->Draw("E");
	tool.DrawText(histSignal,4,false,0.68, 0.75, 0.9, 0.88);
	tool.DrawTextStar(histSignal,2);
	histBackground->Draw("ESAME");

	leg1 = new TLegend(0.6, 0.65, 0.78, 0.74);
	tool.SetLegendStyle(leg1);
	leg1->AddEntry(histSignal,"In+El (unlike-sign pairs)","p");
	leg1->AddEntry(histBackground,"In+El (like-sign pairs)","p");
	leg1->Draw("same");

	if(TEXT){
		TLine *left9 = new TLine(15,0,15,1200);
		tool.SetLineStyle(left9,10,1,4);
	    left9->Draw("same");
	}

	cCanvas->Update();
	cCanvas->Write(variable);
//////////////////////////////////////////
// Plot Eta vertex
	variable = "Eta";
	cutsOption = cuts;
	nInputs = 4;
	nBins = 100;
	min = -2;
	max = 3.5;

	treeBack->Draw(variable + "1>>" + variable + "Bcg1(" + nBins + "," + min + "," + max + ")",cutsOption);
	histBackground = (TH1F*)gPad->GetPrimitive(variable +"Bcg1");
	for(int i = 2; i <= nInputs; ++i)
	{	
		treeBack->Draw(variable + i + ">>" + variable + "Bcg" + i + "(" + nBins + "," + min + "," + max + ")",cutsOption);
		tmpHist = (TH1F*)gPad->GetPrimitive(variable +"Bcg" + i);
		histBackground->Add(tmpHist);
	}
	tool.SetMarkerStyle(histBackground,2,20,1,2,1,1);

	tree->Draw(variable + "1>>" + variable + "Sig1(" + nBins + "," + min + "," + max + ")",cutsOption);
	histSignal = (TH1F*)gPad->GetPrimitive(variable +"Sig1");
	for(int i = 2; i <= nInputs; ++i)
	{	
		tree->Draw(variable + i + ">>" + variable + "Sig" + i + "(" + nBins + "," + min + "," + max + ")",cutsOption);
		tmpHist = (TH1F*)gPad->GetPrimitive(variable +"Sig" + i);
		histSignal->Add(tmpHist);
	}
	histSignal->SetTitle(" ; #eta [-]; Number of tracks");
	tool.SetGraphStyle(histSignal,4,20,1,4,1,1,0.9,1.4);
	tool.SetMarkerStyle(histSignal);
	histSignal->Draw("E");
	tool.DrawText(histSignal,4,false,0.68, 0.75, 0.9, 0.88);
	tool.DrawTextStar(histSignal,2);
	histBackground->Draw("ESAME");

	leg1 = new TLegend(0.6, 0.65, 0.78, 0.74);
	tool.SetLegendStyle(leg1);
	leg1->AddEntry(histSignal,"In+El (unlike-sign pairs)","p");
	leg1->AddEntry(histBackground,"In+El (like-sign pairs)","p");
	leg1->Draw("same");

	if(TEXT){
		TLine *left2 = new TLine(-0.8,0,-0.8,1000);
		tool.SetLineStyle(left2,10,1,4);
	    left2->Draw("same");

		TLine *left3 = new TLine(0.8,0,0.8,800);
		tool.SetLineStyle(left3,10,1,4);
	    left3->Draw("same");
	}
	
	cCanvas->Update();
	cCanvas->Write(variable);
//////////////////////////////////////////
	cCanvas->Close();
}//FourPi::PlotHistogram


