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

	tree->UseCurrentStyle();
	treeBack->UseCurrentStyle();

	Plot tool;
// Plot Inv Mass Inelastic + Elastic
	treeBack->Draw("invMass>>invMassBackground(50,0,2.5)","nSigPPion<3"+cutsWithPrefix);
	TH1F *tmpHist2 = (TH1F*)gPad->GetPrimitive("invMassBackground");
	tool.SetMarkerStyle(tmpHist2,2,20,1,2,1,1);

	tree->Draw("invMass>>invMassSignal(50,0,2.5)","nSigPPion<3"+cutsWithPrefix);
	TH1F *tmpHist = (TH1F*)gPad->GetPrimitive("invMassSignal");
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
	////cCanvas->SaveAs(output + "BasicPlots/invMass.png");
	cCanvas->Write("invMass");
//////////////////////////////////////////
// Plot Inv Mass Inelastic
	treeBack->Draw("invMass>>invMassBackground(50,0,2.5)","nSigPPion<3 && !elastic"+cutsWithPrefix);
	tmpHist2 = (TH1F*)gPad->GetPrimitive("invMassBackground");
	tool.SetMarkerStyle(tmpHist2,2,20,1,2,1,1);

	tree->Draw("invMass>>invMassSignal(50,0,2.5)","nSigPPion<3 && !elastic"+cutsWithPrefix);
	tmpHist = (TH1F*)gPad->GetPrimitive("invMassSignal");
	tmpHist->SetTitle(" ; m(#pi^{+}#pi^{-}) [GeV/c^{2}]; Number of events");
	//tmpHist->GetXaxis()->SetRangeUser(0,2.5);
	tool.SetGraphStyle(tmpHist,4,20,1,4,1,1,0.9,1.3);
	tool.SetMarkerStyle(tmpHist);
	tmpHist->Draw("E");
	tool.DrawText(tmpHist, 1, true);
	tool.DrawTextStar(tmpHist);
	tmpHist2->Draw("ESAME");

	leg1 = new TLegend(0.58, 0.7, 0.78, 0.8);
	tool.SetLegendStyle(leg1);
	leg1->AddEntry(tmpHist,"In (unlike-sign pairs)","p");
	leg1->AddEntry(tmpHist2,"In (like-sign pairs)","p");
	leg1->Draw("same");

	cCanvas->Update();
	////cCanvas->SaveAs(output + "BasicPlots/invMassIn.png");
	cCanvas->Write("invMassIn");
//////////////////////////////////////////
// Plot Inv Mass Elastic
	treeBack->Draw("invMass>>invMassBackground(50,0,2.5)","nSigPPion<3 && elastic"+cutsWithPrefix);
	tmpHist2 = (TH1F*)gPad->GetPrimitive("invMassBackground");
	tool.SetMarkerStyle(tmpHist2,2,20,1,2,1,1);

	tree->Draw("invMass>>invMassSignal(50,0,2.5)","nSigPPion<3 && elastic");
	tmpHist = (TH1F*)gPad->GetPrimitive("invMassSignal");
	tmpHist->SetTitle(" ; m(#pi^{+}#pi^{-}) [GeV/c^{2}]; Number of events");
	//tmpHist->GetXaxis()->SetRangeUser(0,2.5);
	tool.SetGraphStyle(tmpHist,4,20,1,4,1,1,0.9,1.3);
	tool.SetMarkerStyle(tmpHist);
	tmpHist->Draw("E");
	tool.DrawText(tmpHist, 1, true);
	tool.DrawTextStar(tmpHist);
	tmpHist2->Draw("ESAME");

	leg1 = new TLegend(0.58, 0.7, 0.78, 0.8);
	tool.SetLegendStyle(leg1);
	leg1->AddEntry(tmpHist,"El (unlike-sign pairs)","p");
	leg1->AddEntry(tmpHist2,"El (like-sign pairs)","p");
	leg1->Draw("same");

	cCanvas->Update();
	////cCanvas->SaveAs(output + "BasicPlots/invMassEl.png");
	cCanvas->Write("invMassEl");
//////////////////////////////////////////

	hCuts = (TH1F*)data -> Get("AnalysisFlow");
	hMissingPtTPC = (TH1D*)data -> Get("All/MissingPt_TPC2t_Combi");
	hMissingPtTOF = (TH1D*)data -> Get("All/MissingPt_TOF2trk_Combi");
	hMissingPtQ0 = (TH1D*)data -> Get("All/MissingPt_Q0_Combi");
	hMissingPtExc = (TH1D*)data -> Get("All/MissingPt_Excl_Combi");

// //////////////////////////////////////////////////////////
// Plot Cuts Flow
	TString Labels[] = { TString("All"), TString("CPT trigger"), TString("El+InEl"), TString("2 TPC-TOF tracks"), 
	                  TString("Same vertex"), TString("TotCharge 0"), TString("p_{T}^{miss} < 0.1 GeV/c"), TString(""),};
	//gPad->SetMargin(0.9,0.02,0.1,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
	gPad->SetLogy();

	TH1F* hCutsFlow = new TH1F("cuts", "cuts", 8, 1, 9);
	int i = 1;
	for(int tb=1; tb<8; ++tb){ 
		if(tb==3){
			hCutsFlow->GetXaxis()->SetBinLabel(tb, Labels[tb-1]);
			hCutsFlow->SetBinContent(tb,hCuts->GetBinContent(i)+hCuts->GetBinContent((++i)++));
			cout<<hCutsFlow->GetBinContent(tb)<<endl;
			continue;
		}
		hCutsFlow->GetXaxis()->SetBinLabel(tb, Labels[tb-1]);
		hCutsFlow->SetBinContent(tb,hCuts->GetBinContent(i));
		cout<<hCutsFlow->GetBinContent(tb)<<endl;
		i++;
	}
	//hCutsFlow->GetXaxis()->SetBinLabel(8, "PID");
	//tree->Draw("nSigPPion>>nSigPPion","nSigPPion<3");
	//hCutsFlow->SetBinContent(8,((TH1F*)gPad->GetPrimitive("nSigPPion"))->GetEntries());
	hCutsFlow->SetTitle("; ; Number of events");
	tool.SetGraphStyle(hCutsFlow);
	tool.SetMarkerStyle(hCutsFlow);
	//hCutsFlow->GetYaxis()->SetRangeUser(800,300000000);
	hCutsFlow->GetXaxis()->SetLabelSize(0.045);
	hCutsFlow->Draw();
	tool.DrawTextStar(tmpHist,2);

	TPaveText *textPub = new TPaveText(0.75,0.74,0.9,0.9,"brNDC");
	tool.SetTextStyle(textPub);
	textPub -> AddText("p + p #rightarrow p' + X + p'");
	textPub -> AddText("#sqrt{s} = 510 GeV");
	textPub -> AddText("Cuts flow");
	textPub -> Draw("same");

	
	cCanvas->Update();
	//cCanvas->SaveAs(output + "BasicPlots/Cuts.png");
	cCanvas->Write("CutsFlow");
//////////////////////////////////////////
////////////// 
	gPad->SetLogy();
	hMissingPtTPC->SetStats(0);
	hMissingPtTPC->SetTitle(" ; p_{T}^{miss} [GeV/c];Number of events");
	tool.SetGraphStyle(hMissingPtTPC);
	//tool.SetMarkerStyle(hMissingPtTPC);
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
	leg1 -> AddEntry(hMissingPtTPC, "2 TPC-TOF tracks", "l");
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
	textPub2 -> AddText("p + p #rightarrow p' + X + p'");
	textPub2 -> AddText("#sqrt{s} = 510 GeV");
	textPub2 -> Draw("same");
	tool.DrawTextStar(tmpHist,2);

	TLine *left = new TLine(0.1,0,0.1,600000);
	tool.SetLineStyle(left,10,1,4);
    left->Draw("same");

	cCanvas->Update();
	cCanvas-> Write("hMissingPt");
	//cCanvas->SaveAs(output + "BasicPlots/cMissingPt.png");
	cCanvas->Close();

  
}//BasicPlots::PlotHistogram

