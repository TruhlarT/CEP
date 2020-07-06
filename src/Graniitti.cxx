#include "Graniitti.h"
#include "Plot.h"

//_____________________________________________________________________________
Graniitti::Graniitti(TFile* dataInput, TFile* graniittiInput, TFile* fileOut, TString outnam, bool text, TString inputCuts)
{
	//constructor
	output = outnam;

	data = dataInput;
	graniitti = graniittiInput;
	fout = fileOut;

	cuts = inputCuts;
	if(inputCuts != "")
		cutsWithPrefix = " && " + inputCuts;
	else
		cutsWithPrefix="";

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

	TEXT = text;
	cout << "Graniitti::Graniitti() called" << endl;

}//Graniitti

//_____________________________________________________________________________
Graniitti::~Graniitti()
{
  //destructor

  cout << "Graniitti::~Graniitti() destructor called" << endl;

}//~Graniitti


//_____________________________________________________________________________
void Graniitti::PlotHistogram(){

	if (!data){
		cout<<"Error: cannot open input file"<<endl;
		return;
    }

	TTree* treeSig = dynamic_cast<TTree*>( data->Get("recTree") );
	TTree* treeBack = dynamic_cast<TTree*>( data->Get("Background") );
	
	if (!treeSig || !treeBack){
		cout<<"Error: cannot open one of the TTree"<<endl;
		return;
    }

	if (!fout){
		cout<<"Error: cannot open output file"<<endl;
		return;
    }

    TTree* tree[nParticles+1]; // 3 = 4PI state
	tree[Pion] = dynamic_cast<TTree*>( graniitti->Get("pionTree") );
	tree[Kaon] = dynamic_cast<TTree*>( graniitti->Get("kaonTree") );
	tree[Proton] = dynamic_cast<TTree*>( graniitti->Get("protonTree") );
	tree[3] = dynamic_cast<TTree*>( graniitti->Get("4pionTree") );

	if (!tree[Pion] || !tree[Kaon] || !tree[Proton] || !tree[3]){
		cout<<"Error: cannot open one of the TTree in Graniitti"<<endl;
		return;
    }

	TString granForwardCuts = 	TString("px_proton1 > - 0.3 && px_proton1 < 0.5 && abs(py_proton1) < 0.9 && abs(py_proton1) > 0.3 &&")
								+ TString("abs(py_proton1) < -px_proton1+1.05 && px_proton2 > - 0.3 && px_proton2 < 0.5 &&")
								+ TString(" abs(py_proton2) < 0.9 && abs(py_proton2) > 0.3 && abs(py_proton2) < -px_proton2+1.05");
//	TString granForwardCuts = 	TString("((theta_proton1 > 0.005/((phi_proton1+0.7)^5) + 0.0011 + 0.0001*phi_proton1*phi_proton1 &&")
//								+ TString("theta_proton1 < 0.005/((phi_proton1+0.6)^4) + 0.0032 + 0.0001*phi_proton1*phi_proton1 &&")
//								+ TString("theta_proton1 < 0.02/((phi_proton1-0.2)^3) -0.0008) ||")
//								+ TString("(theta_proton1 > 0.005/((-phi_proton1+0.7)^5) + 0.0011 + 0.0001*phi_proton1*phi_proton1 &&")
//								+ TString("theta_proton1 < 0.005/((-phi_proton1+0.55)^4) + 0.0034 + 0.0003*phi_proton1*phi_proton1 &&")
//								+ TString("theta_proton1 < 0.02/((-phi_proton1-0.2)^3) -0.0008)) && ") 
//								+ TString("((theta_proton2 > 0.005/((phi_proton2+0.7)^5) + 0.0011 + 0.0001*phi_proton2*phi_proton2 &&")
//								+ TString("theta_proton2 < 0.005/((phi_proton2+0.6)^4) + 0.0032 + 0.0001*phi_proton2*phi_proton2 &&")
//								+ TString("theta_proton2 < 0.02/((phi_proton2-0.2)^3) -0.0008) ||")
//								+ TString("(theta_proton2 > 0.005/((-phi_proton2+0.7)^5) + 0.0011 + 0.0001*phi_proton2*phi_proton2 &&")
//								+ TString("theta_proton2 < 0.005/((-phi_proton2+0.55)^4) + 0.0034 + 0.0003*phi_proton2*phi_proton2 &&")
//								+ TString("theta_proton2 < 0.02/((-phi_proton2-0.2)^3) -0.0008))");
	TString graniittiCuts = "eta_part1 > -0.7 && eta_part1 < 0.7 && eta_part2 > -0.7 && eta_part2 < 0.7 && t_proton1 < -0.12 && t_proton2 < -0.12 && t_proton1 > -1.0  && t_proton2 > -1.0 && " + granForwardCuts;
	TString partCuts[nParticles+1];
	partCuts[Pion] = "pT_part1 > 0.2 && pT_part2 > 0.2";
	partCuts[Kaon] = "pT_part1 > 0.3 && pT_part2 > 0.3 && TMath::Min(pT_part1,pT_part2) < 0.7";
	partCuts[Proton] = "pT_part1 > 0.4 && pT_part2 > 0.4 && TMath::Min(pT_part1,pT_part2) < 1.1";
	partCuts[3] = "pT_part1 > 0.2 && pT_part2 > 0.2 && pT_part3 > 0.2 && pT_part4 > 0.2";


	Double_t binning[4][3] =  {{64, 0.3, 3.5},{44, 0.8, 3}, {24, 1.6, 4}, {50,0.5,4.5}};

	TCanvas *cCanvas = new TCanvas("cCanvas","cCanvas",800,700);
	gPad->SetMargin(0.12,0.02,0.1,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
	gPad->SetTickx();
	gPad->SetTicky();  
	gStyle->SetOptStat("");
	gStyle->SetPalette(1);
	gStyle->SetLineWidth(2);      //axis line
	gStyle->SetFrameLineWidth(2); //frame line

	Plot tool;

	TString variable, usedCuts;	
	TString combLabels[] = { TString("El+Inel"), TString("El"), TString("Inel")};
	TString combCuts[] = { TString(""), TString("&& TMath::Abs(deltaphi_pp) > 1.570796327"), TString("&& TMath::Abs(deltaphi_pp) < 1.570796327")}; // 57.2957795 = 1 rad => 1.6 = 90 cca
	TString particleLables[nParticles+1] = { TString("Pion"), TString("Kaon"), TString("Proton"), TString("4Pion")};
	TString stateLabel[] = {"#pi^{+}#pi^{-}","K^{+}K^{-}","p#bar{p}","#pi^{+}#pi^{+}#pi^{-}#pi^{-}"};

	Int_t nBins;
	Float_t min, max;

	TDirectory* currDir = fout->GetDirectory("Graniitti");
	if(!currDir){
		cout<<"Graniitti directory not found!!"<<endl;
		cout<<"--------Error-----"<<endl;
		cout<<"--------Error-----"<<endl;
		cout<<"--------Error-----"<<endl;
		cout<<"--------Error-----"<<endl;
		cout<<"------------------"<<endl;
		return;
	}

	for (int iComb = 0; iComb < 3; ++iComb)
	{
		dataLabel = combLabels[iComb];
		currDir->mkdir(combLabels[iComb])->cd();
		usedCuts= graniittiCuts + combCuts[iComb];
		for (int iPart = 0; iPart < nParticles + 1; ++iPart)
		{
			// Plot t 
			usedCuts+= " && " + partCuts[iPart];
			//cout<<usedCuts<<endl;	
			variable = "t";
			nBins = 100;
			min = -2.0;
			max = 0.0;
			tree[iPart]->Draw(variable +"_proton1>>" + variable +"Sig1(" + nBins + "," + min + "," + max + ")",usedCuts);
			tmpHist = (TH1F*)gPad->GetPrimitive(variable +"Sig1");
			tree[iPart]->Draw(variable +"_proton2>>" + variable +"Sig2(" + nBins + "," + min + "," + max + ")",usedCuts);
			tmpHist2 = (TH1F*)gPad->GetPrimitive(variable +"Sig2");
			tmpHist->Add(tmpHist2);   
			tmpHist->SetTitle(" ; t; Number of protons");
			tool.SetGraphStyle(tmpHist,4,20,1,4,1,1,0.9,1.4);
			tool.SetMarkerStyle(tmpHist);
			tmpHist->Draw("E");
			tool.DrawText(tmpHist,0,false,0.15,0.74,0.28,0.9,12, true);
			tool.DrawTextStar(tmpHist,3);
			//cout<<tmpHist->GetEntries()<<endl;
			                        
			TLegend* leg1 = new TLegend(0.15, 0.65, 0.28, 0.74);
			tool.SetLegendStyle(leg1);
			leg1->AddEntry(tmpHist,dataLabel + " Graniitti","p");
			leg1->Draw("same");


			cCanvas->Update();
			cCanvas->Write(variable + particleLables[iPart]);
			//////////////////////////////////////
			// Plot tSum
			variable = "t";
			nBins = 100;
			min = 0.0;
			max = 4.0;
			tree[iPart]->Draw("TMath::Abs( " + variable +"_proton1 + " + variable +"_proton2)>>" + variable +"Sig1(" + nBins + "," + min + "," + max + ")",usedCuts);
			tmpHist = (TH1F*)gPad->GetPrimitive(variable +"Sig1");   
			tmpHist->SetTitle(" ; |t_{1} + t_{2}| [GeV^{2}]; Number of events");
			tool.SetGraphStyle(tmpHist,4,20,1,4,1,1,0.9,1.4);
			tool.SetMarkerStyle(tmpHist);
			tmpHist->Draw("E");
			tool.DrawText(tmpHist,0,false,0.68, 0.75, 0.9, 0.88, 22, true);
			tool.DrawTextStar(tmpHist,2);


			leg1 = new TLegend(0.6, 0.65, 0.78, 0.74);
			tool.SetLegendStyle(leg1);
			leg1->AddEntry(tmpHist,dataLabel + " Graniitti","p");
			leg1->Draw("same");

			cCanvas->Update();
			cCanvas->Write(variable + "Sum" + particleLables[iPart]);    
			//////////////////////////////////////////////////
			// Plot phi 
			variable = "phi";
			nBins = 100;
			min = 0.0;
			max = 350.0;

			tree[iPart]->Draw("TMath::Abs(deltaphi_pp*57.2957795)>>" + variable +"Sig1(" + nBins + "," + min + "," + max + ")",usedCuts);
			tmpHist = (TH1F*)gPad->GetPrimitive(variable +"Sig1");   
			tmpHist->SetTitle(" ; |#Delta#varphi| [deg]; Number of events");
			tool.SetGraphStyle(tmpHist,4,20,1,4,1,1,0.9,1.4);
			tool.SetMarkerStyle(tmpHist);
			tmpHist->Draw("E");
			tool.DrawText(tmpHist,0,false,0.68, 0.75, 0.9, 0.88, 22, true);
			tool.DrawTextStar(tmpHist,2);


			leg1 = new TLegend(0.6, 0.65, 0.78, 0.74);
			tool.SetLegendStyle(leg1);
			leg1->AddEntry(tmpHist,dataLabel + " Graniitti","p");
			leg1->Draw("same");


			cCanvas->Update();
			cCanvas->Write(variable + particleLables[iPart]); 
			//////////////////////////////////////////////
			// Plot inv mass
			variable = "invMass_state";
			nBins = binning[iPart][0];
			min = binning[iPart][1];
			max = binning[iPart][2];

			tree[iPart]->Draw(variable +">>" + variable +"Sig1(" + nBins + "," + min + "," + max + ")",usedCuts);
			tmpHist = (TH1F*)gPad->GetPrimitive(variable +"Sig1");   
			tmpHist->SetTitle(" ; m(" + stateLabel[iPart] + ") [GeV/c^{2}]; Number of events");
			tool.SetGraphStyle(tmpHist,4,20,1,4,1,1,0.9,1.3);
			tool.SetMarkerStyle(tmpHist);
			tmpHist->Draw("E");
			tool.DrawText(tmpHist, iPart+1, true,  0.68, 0.75, 0.9, 0.88, 22, true);
			tool.DrawTextStar(tmpHist, 1);
			tmpHist->Write("invMass" + particleLables[iPart]);

			leg1 = new TLegend(0.58, 0.64, 0.78, 0.74);
			tool.SetLegendStyle(leg1);
			leg1->AddEntry(tmpHist, combLabels[iComb] + " Graniitti","p");
			leg1->Draw("same");

			cCanvas->Update();
			cCanvas->Write(variable + particleLables[iPart]);
			//////////////////////////////////////////////
			// Plot pair rapidity
			variable = "rap_state";
			nBins = 35;
			min = -0.7;
			max = 1.1;

			tree[iPart]->Draw(variable +">>" + variable +"Sig1(" + nBins + "," + min + "," + max + ")",usedCuts);
			tmpHist = (TH1F*)gPad->GetPrimitive(variable +"Sig1");   
			tmpHist->SetTitle(" ; y(" + stateLabel[iPart] + "); Number of events");
			tool.SetGraphStyle(tmpHist,4,20,1,4,1,1,0.9,1.3);
			tool.SetMarkerStyle(tmpHist);
			tmpHist->Draw("E");
			tool.DrawText(tmpHist, iPart+1, true,  0.68, 0.75, 0.9, 0.88, 22, true);
			tool.DrawTextStar(tmpHist, 1);

			leg1 = new TLegend(0.58, 0.64, 0.78, 0.74);
			tool.SetLegendStyle(leg1);
			leg1->AddEntry(tmpHist, combLabels[iComb] + " Graniitti","p");
			leg1->Draw("same");

			cCanvas->Update();
			cCanvas->Write(variable + particleLables[iPart]);
		}
	}

/////////////////////////////////////////////////////
	cCanvas->Close();

}


//_____________________________________________________________________________
void Graniitti::PlotComparison(){


}