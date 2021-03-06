#include "trackQuality.h"
#include "Plot.h"

//_____________________________________________________________________________
trackQuality::trackQuality(TFile* dataInput, TFile* fileOut, TString outnam, bool text, TString inputCuts)
{
	//constructor
	output = outnam;

	cuts = inputCuts;
	if(inputCuts != "")
		cutsWithPrefix = " && " + inputCuts;
	else
		cutsWithPrefix="";

	data = dataInput;
	fout = fileOut;

	TEXT = text;

	dataLabel = "Data";
	/*
	if(inputCuts.Contains("!elastic"))
	{
        dataLabel = "Inel";
	}
    else if(inputCuts.Contains("elastic"))
    {
        dataLabel = "El";
    }
    else
    {
        dataLabel = "El+Inel";
    }
	*/
	cout << "trackQuality::trackQuality() called" << endl;

}//trackQuality

//_____________________________________________________________________________
trackQuality::~trackQuality()
{
  //destructor

  cout << "trackQuality::~trackQuality() destructor called" << endl;

}//~trackQuality


//_____________________________________________________________________________
void trackQuality::PlotHistogram(){

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

	Int_t nBins, nBins2, nInputs;
	Float_t min, max, min2, max2;

	Plot tool;
// Plot z vertex
	treeBack->Draw("vertexZ>>vertexZBcg",cuts);
	tmpHist2 = (TH1F*)gPad->GetPrimitive("vertexZBcg");
	tool.SetMarkerStyle(tmpHist2,2,20,1,2,1,1);
	
	tree->Draw("vertexZ>>vertexZSig",cuts);
	tmpHist = (TH1F*)gPad->GetPrimitive("vertexZSig");
	tmpHist->SetTitle(" ; Z_{vrtx} [cm]; Number of events");
	//tmpHist->GetXaxis()->SetRangeUser(0,2.5);
	tool.SetGraphStyle(tmpHist,4,20,1,4,1,1,0.9,1.3);
	tool.SetMarkerStyle(tmpHist);
	tmpHist->Draw("E");
	tool.DrawText(tmpHist,0,true,0.68, 0.75, 0.9, 0.88);
	tool.DrawTextStar(tmpHist,2);
	tmpHist2->Draw("ESAME");
	TLegend* leg1 = new TLegend(0.15,0.83,0.28,0.93); 
	tool.SetLegendStyle(leg1);
	leg1->AddEntry(tmpHist,dataLabel + " (unlike-sign pairs)","p");
	leg1->AddEntry(tmpHist2,dataLabel + " (like-sign pairs)","p");
	leg1->Draw("same");
	if(TEXT){
		TLine *left4 = new TLine(-80,0,-80,tmpHist->GetMaximum()/2);
		tool.SetLineStyle(left4,10,1,4);
	   left4->Draw("same");

		TLine *left5 = new TLine(80,0,80,tmpHist->GetMaximum()/2);
		tool.SetLineStyle(left5,10,1,4);
	   left5->Draw("same");		
	}
	cCanvas->Update();
	//cCanvas->SaveAs( output + "trackQuality/zVertex.png");
	cCanvas->Write("zVertex");
//////////////////////////////////////////
// Plot DCAXY vertex
	treeBack->Draw("DcaXY0>>DcaXY1Bcg(100,0,3.5)",cuts);
	tmpHist2 = (TH1F*)gPad->GetPrimitive("DcaXY1Bcg");
	treeBack->Draw("DcaXY1>>DcaXY2Bcg(100,0,3.5)",cuts);
	tmpHist = (TH1F*)gPad->GetPrimitive("DcaXY2Bcg");
	tmpHist2->Add(tmpHist);
	tool.SetMarkerStyle(tmpHist2,2,20,1,2,1,1);

	tree->Draw("DcaXY0>>DcaXY1Sig(100,0,3.5)",cuts);
	tmpHist = (TH1F*)gPad->GetPrimitive("DcaXY1Sig");
	tree->Draw("DcaXY1>>DcaXY2Sig(100,0,3.5)",cuts);
	tmpHist3 = (TH1F*)gPad->GetPrimitive("DcaXY2Sig");
	tmpHist->Add(tmpHist3);
	tmpHist->SetTitle(" ; DCA_{xy} [cm]; Number of tracks");
	tool.SetGraphStyle(tmpHist,4,20,1,4,1,1,0.9,1.4);
	tool.SetMarkerStyle(tmpHist);
	tmpHist->Draw("E");
	tool.DrawText(tmpHist,0,false,0.68, 0.75, 0.9, 0.88);
	tool.DrawTextStar(tmpHist,2);
	tmpHist2->Draw("ESAME");

	leg1 = new TLegend(0.6, 0.65, 0.78, 0.74);
	tool.SetLegendStyle(leg1);
	leg1->AddEntry(tmpHist,dataLabel + " (unlike-sign pairs)","p");
	leg1->AddEntry(tmpHist2,dataLabel + " (like-sign pairs)","p");
	leg1->Draw("same");
	if(TEXT){
		TLine *left6 = new TLine(1.5,0,1.5,tmpHist->GetMaximum()/2);
		tool.SetLineStyle(left6,10,1,4);
	   left6->Draw("same");
	}
	cCanvas->Update();
	//cCanvas->SaveAs( output + "trackQuality/DcaXY.png");
	cCanvas->Write("DcaXY");
//////////////////////////////////////////
// Plot DCAZ vertex
	treeBack->Draw("DcaZ0>>DcaZ1Bcg(100,-1.5,1.5)",cuts);
	tmpHist2 = (TH1F*)gPad->GetPrimitive("DcaZ1Bcg");
	treeBack->Draw("DcaZ1>>DcaZ2Bcg(100,-1.5,1.5)",cuts);
	tmpHist = (TH1F*)gPad->GetPrimitive("DcaZ2Bcg");
	tmpHist2->Add(tmpHist);
	tool.SetMarkerStyle(tmpHist2,2,20,1,2,1,1);

	tree->Draw("DcaZ1>>DcaZ1Sig(100,-1.5,1.5)",cuts);
	tmpHist = (TH1F*)gPad->GetPrimitive("DcaZ1Sig");
	tree->Draw("DcaZ0>>DcaZ2Sig(100,-1.5,1.5)",cuts);
	tmpHist3 = (TH1F*)gPad->GetPrimitive("DcaZ2Sig");
	tmpHist->Add(tmpHist3);
	tmpHist->SetTitle(" ; DCA_{z} [cm]; Number of tracks");
	tool.SetGraphStyle(tmpHist,4,20,1,4,1,1,0.9,1.4);
	tool.SetMarkerStyle(tmpHist);
	tmpHist->Draw("E");
	tool.DrawText(tmpHist,0,false,0.68, 0.75, 0.9, 0.88);
	tool.DrawTextStar(tmpHist,2);
	tmpHist2->Draw("ESAME");
	leg1 = new TLegend(0.6, 0.65, 0.78, 0.74);
	tool.SetLegendStyle(leg1);
	leg1->AddEntry(tmpHist,dataLabel + " (unlike-sign pairs)","p");
	leg1->AddEntry(tmpHist2,dataLabel + " (like-sign pairs)","p");
	leg1->Draw("same");
	if(TEXT){
		TLine *left7 = new TLine(1,0,1,tmpHist->GetMaximum()/2);
		tool.SetLineStyle(left7,10,1,4);
	   left7->Draw("same");

		TLine *left8 = new TLine(-1,0,-1,tmpHist->GetMaximum()/2);
		tool.SetLineStyle(left8,10,1,4);
	   left8->Draw("same");
	}
	cCanvas->Update();
	//cCanvas->SaveAs( output + "trackQuality/DcaZ.png");
	cCanvas->Write("DcaZ");
//////////////////////////////////////////
	// Plot NhitsFit 
	TString variable = "NhitsFit";
	treeBack->Draw(variable +"1>>" + variable +"Bcg1(61,0,60)",cuts);
	tmpHist = (TH1F*)gPad->GetPrimitive(variable +"Bcg1");
	treeBack->Draw(variable +"0>>" + variable +"Bcg2(61,0,60)",cuts);
	tmpHist2 = (TH1F*)gPad->GetPrimitive(variable +"Bcg2");
	tmpHist->Add(tmpHist2);
	tool.SetMarkerStyle(tmpHist2,2,20,1,2,1,1);

	tree->Draw(variable +"1>>" + variable +"Sig1(61,0,60)",cuts);
	tmpHist = (TH1F*)gPad->GetPrimitive(variable +"Sig1");
	tree->Draw(variable +"0>>" + variable +"Sig2(61,0,60)",cuts);
	tmpHist3 = (TH1F*)gPad->GetPrimitive(variable +"Sig2");
	tmpHist->Add(tmpHist3);
	tmpHist->SetTitle(" ; N^{fit}_{hits} ; Number of tracks");
	tool.SetGraphStyle(tmpHist,4,20,1,4,1,1,0.9,1.4);
	tool.SetMarkerStyle(tmpHist);
	tmpHist->GetXaxis()->SetRangeUser(10,60);
	tmpHist->Draw("E");
	tool.DrawText(tmpHist,0,false,0.68, 0.75, 0.9, 0.88);
	tool.DrawTextStar(tmpHist,2);
	tmpHist2->Draw("ESAME");


	leg1 = new TLegend(0.6, 0.65, 0.78, 0.74);
	tool.SetLegendStyle(leg1);
	leg1->AddEntry(tmpHist,dataLabel + " (unlike-sign pairs)","p");
	leg1->AddEntry(tmpHist2,dataLabel + " (like-sign pairs)","p");
	leg1->Draw("same");

	if(TEXT){
		TLine *left9 = new TLine(25,0,25,tmpHist->GetMaximum()/2);
		tool.SetLineStyle(left9,10,1,4);
	   left9->Draw("same");
	}

	cCanvas->Update();
	//cCanvas->SaveAs( output + "trackQuality/" + variable + ".png");
	cCanvas->Write(variable);


/////////////////////////////////////////////
	// Plot NhitsDEdx 
	variable = "NhitsDEdx";

	treeBack->Draw(variable +"1>>" + variable +"Bcg1(61,0,60)",cuts);
	tmpHist = (TH1F*)gPad->GetPrimitive(variable +"Bcg1");
	treeBack->Draw(variable +"0>>" + variable +"Bcg2(61,0,60)",cuts);
	tmpHist2 = (TH1F*)gPad->GetPrimitive(variable +"Bcg2");
	tmpHist->Add(tmpHist2);
	tool.SetMarkerStyle(tmpHist2,2,20,1,2,1,1);

	tree->Draw(variable +"1>>" + variable +"Sig1(61,0,60)",cuts);
	tmpHist = (TH1F*)gPad->GetPrimitive(variable +"Sig1");
	tree->Draw(variable +"0>>" + variable +"Sig2(61,0,60)",cuts);
	tmpHist3 = (TH1F*)gPad->GetPrimitive(variable +"Sig2");
	tmpHist->Add(tmpHist3);
	tmpHist->SetTitle(" ; N^{dE/dx}_{hits} ; Number of tracks");
	tool.SetGraphStyle(tmpHist,4,20,1,4,1,1,0.9,1.4);
	tool.SetMarkerStyle(tmpHist);
	tmpHist->GetXaxis()->SetRangeUser(5,60);
	tmpHist->Draw("E");
	tool.DrawText(tmpHist,0,false,0.68, 0.75, 0.9, 0.88);
	tool.DrawTextStar(tmpHist,2);
	tmpHist2->Draw("ESAME");


	leg1 = new TLegend(0.6, 0.65, 0.78, 0.74);
	tool.SetLegendStyle(leg1);
	leg1->AddEntry(tmpHist,dataLabel + " (unlike-sign pairs)","p");
	leg1->AddEntry(tmpHist2,dataLabel + " (like-sign pairs)","p");
	leg1->Draw("same");

	if(TEXT){
		TLine *left9 = new TLine(15,0,15,tmpHist->GetMaximum()/2);
		tool.SetLineStyle(left9,10,1,4);
	   left9->Draw("same");
	}

	cCanvas->Update();
	//cCanvas->SaveAs( output + "trackQuality/" + variable + ".png");
	cCanvas->Write(variable);	
//////////////////////////////////////////
// Plot Eta 
	treeBack->Draw("Eta1>>Eta1Bcg(100,-2,3.5)",cuts);
	tmpHist2 = (TH1F*)gPad->GetPrimitive("Eta1Bcg");
	treeBack->Draw("Eta0>>Eta2Bcg(100,-2,3.5)",cuts);
	tmpHist = (TH1F*)gPad->GetPrimitive("Eta2Bcg");
	tmpHist2->Add(tmpHist);
	tool.SetMarkerStyle(tmpHist2,2,20,1,2,1,1);

	tree->Draw("Eta1>>Eta1Sig(100,-2,3.5)",cuts);
	tmpHist = (TH1F*)gPad->GetPrimitive("Eta1Sig");
	tree->Draw("Eta0>>Eta2Sig(100,-2,3.5)",cuts);
	tmpHist3 = (TH1F*)gPad->GetPrimitive("Eta2Sig");
	tmpHist->Add(tmpHist3);
	tmpHist->SetTitle(" ; #eta ; Number of tracks");
	tool.SetGraphStyle(tmpHist,4,20,1,4,1,1,0.9,1.3);
	tool.SetMarkerStyle(tmpHist);
	tmpHist->Draw("E");
	tool.DrawText(tmpHist,0,false,0.68, 0.75, 0.9, 0.88);
	tool.DrawTextStar(tmpHist,2);
	tmpHist2->Draw("ESAME");

	leg1 = new TLegend(0.6, 0.65, 0.8, 0.74);
	tool.SetLegendStyle(leg1);
	leg1->AddEntry(tmpHist,dataLabel + " (unlike-sign pairs)","p");
	leg1->AddEntry(tmpHist2,dataLabel + " (like-sign pairs)","p");
	leg1->Draw("same");
	if(TEXT){
		TLine *left2 = new TLine(-0.7,0,-0.7,tmpHist->GetMaximum()/2);
		tool.SetLineStyle(left2,10,1,4);
	   left2->Draw("same");

		TLine *left3 = new TLine(0.7,0,0.7,tmpHist->GetMaximum()/2);
		tool.SetLineStyle(left3,10,1,4);
	   left3->Draw("same");
	}
	cCanvas->Update();
	//cCanvas->SaveAs( output + "trackQuality/Eta.png");
	cCanvas->Write("Eta");
//////////////////////////////////////////
// Plot t 
    variable = "t";
    nBins = 100;
    min = -2.0;
    max = 0.0;
    treeBack->Draw(variable +"East>>" + variable +"Bcg1(" + nBins + "," + min + "," + max + ")",cuts);
    histBackground = (TH1F*)gPad->GetPrimitive(variable +"Bcg1");
    treeBack->Draw(variable +"West>>" + variable +"Bcg2(" + nBins + "," + min + "," + max + ")",cuts);
    tmpHist = (TH1F*)gPad->GetPrimitive(variable +"Bcg2");
    histBackground->Add(tmpHist);
    tool.SetMarkerStyle(histBackground,2,20,1,2,1,1);


    tree->Draw(variable +"East>>" + variable +"Sig1(" + nBins + "," + min + "," + max + ")",cuts);
    histSignal = (TH1F*)gPad->GetPrimitive(variable +"Sig1");
    tree->Draw(variable +"West>>" + variable +"Sig2(" + nBins + "," + min + "," + max + ")",cuts);
    tmpHist = (TH1F*)gPad->GetPrimitive(variable +"Sig2");
    histSignal->Add(tmpHist);   
    histSignal->SetTitle(" ; t; Number of tracks");
    tool.SetGraphStyle(histSignal,4,20,1,4,1,1,0.9,1.4);
    tool.SetMarkerStyle(histSignal);
    histSignal->Draw("E");
    tool.DrawText(histSignal,0,false,0.15,0.74,0.28,0.9,12);
    tool.DrawTextStar(histSignal,3);
    histBackground->Draw("ESAME");

    if(TEXT){
        TLine *left4 = new TLine(-1.0,0,-1.0,tmpHist->GetMaximum()/2);
        tool.SetLineStyle(left4,10,1,4);
        left4->Draw("same");

        TLine *left5 = new TLine(-0.12,0,-0.12  ,tmpHist->GetMaximum()/2);
        tool.SetLineStyle(left5,10,1,4);
        left5->Draw("same");     
    }
                                
    leg1 = new TLegend(0.15, 0.65, 0.28, 0.74);
    tool.SetLegendStyle(leg1);
    leg1->AddEntry(histSignal,dataLabel + " (unlike-sign pairs)","p");
    leg1->AddEntry(histBackground,dataLabel + " (like-sign pairs)","p");
    leg1->Draw("same");


    cCanvas->Update();
    cCanvas->Write(variable);
///////////////////////////////////////////////////////////////////////

	TCanvas *cCanvas2D = new TCanvas("cCanvas2D","cCanvas2D",800,700);
	gPad->SetMargin(0.09,0.13,0.1,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
	gStyle->SetPalette(1);
	gPad->SetTickx();
	gPad->SetTicky(); 

//////////////////////////////////////////
// Plot ZvrtxVsEta vertex
	tree->Draw("Eta0:vertexZ>>ZvrtxVsEta2Sig(100,-200,200,100,-1,1.5)",cuts,"colz");
	tmp2DHist2 = (TH2F*)gPad->GetPrimitive("ZvrtxVsEta2Sig");
	tree->Draw("Eta1:vertexZ>>ZvrtxVsEtaSig(100,-200,200,100,-1,1.5)",cuts,"colz");
	tmp2DHist = (TH2F*)gPad->GetPrimitive("ZvrtxVsEtaSig");
	tmp2DHist->Add(tmp2DHist2);
	tmp2DHist->SetTitle(" ; Z_{vrtx} [cm]; #eta");
	tool.SetGraphStyle(tmp2DHist,4,20,1,4,1,1,0.9,1.0);
	tmp2DHist->Draw("colz");
	tool.DrawText(tmp2DHist,0,true,0.61,0.75,0.76,0.9);
	tool.DrawTextStar(tmp2DHist,1);
	if(TEXT){
		TLine *left = new TLine(-80,-0.7,-80,0.7);
		tool.SetLineStyle(left);
	   left->Draw("same");

		TLine *left1 = new TLine(80,-0.7,80,0.7);
		tool.SetLineStyle(left1);
	   left1->Draw("same");

		TLine *line = new TLine(-80,0.7,80,0.7);
		tool.SetLineStyle(line);
	   line->Draw("same");

		TLine *line2 = new TLine(-80,-0.7,80,-0.7);
		tool.SetLineStyle(line2);
	   line2->Draw("same");	
	}
	cCanvas2D->Update();
	//cCanvas2D->SaveAs( output + "trackQuality/ZvrtxVsEta.png");
	cCanvas2D->Write("ZvrtxVsEta");
//////////////////////////////////////////
// Plot PhiVsDCAXY vertex
	tree->Draw("Phi0:DcaXY0>>PhiVsDCAXY2Sig(50,0,4,50,-3.5,3.5)",cuts,"colz");
	tmp2DHist2 = (TH2F*)gPad->GetPrimitive("PhiVsDCAXY2Sig");
	tree->Draw("Phi1:DcaXY1>>PhiVsDCAXYSig(50,0,4,50,-3.5,3.5)",cuts,"colz");
	tmp2DHist = (TH2F*)gPad->GetPrimitive("PhiVsDCAXYSig");
	tmp2DHist->Add(tmp2DHist2);
	tmp2DHist->SetTitle(" ; DCA_{xy} [cm] ; #phi");
	tool.SetGraphStyle(tmp2DHist,4,20,1,4,1,1,0.9,1.0);
	tmp2DHist->Draw("colz");
	tool.DrawText(tmp2DHist,0,true,0.61,0.75,0.76,0.9);
	tool.DrawTextStar(tmp2DHist,1);

	cCanvas2D->Update();
	//cCanvas2D->SaveAs( output + "trackQuality/PhiVsDCAXY.png");
	cCanvas2D->Write("PhiVsDCAXY");
//////////////////////////////////////////
	cCanvas->Close();
	cCanvas2D->Close();
}

