#include <TStyle.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TCanvas.h>
#include <TNtupleD.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <RooRealVar.h>
#include <RooConstVar.h>
#include <RooDataSet.h>
#include <RooGaussian.h>
#include <RooChebychev.h>
#include <RooBernstein.h>
#include <RooExponential.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include "myloop.h"
#include "plotDressing.h"
#include <iostream>

using namespace RooFit;

// General fitting options
#define NUMBER_OF_CPU       1
#define DO_MINOS            kTRUE
// 0 - w/o DISPLAY
// 1 - w/  DISPLAY
#define DISPLAY             1

#define SOURCE              "myloop_reference.root"

//-----------------------------------------------------------------
// Definition of channel #
// channel = 1: B+ -> J/psi K+
// channel = 2: B0 -> J/psi K*
// channel = 3: B0 -> J/psi Ks
// channel = 4: Bs -> J/psi phi
// channel = 5: Jpsi + pipi
// channel = 6: Lambda_b -> Jpsi + Lambda

void correctedYield(int channel = 1, int beamSpotErrEstimate = 1)
{
    double ptBins[] = {10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 100.};
    int ptBinsSize = (sizeof(ptBins)/sizeof(*ptBins))-1;
    cout<<"pT bin size = "<<ptBinsSize<<'\t'<<ptBins[0]<<endl;
    
    TH1D *bPtHist = new TH1D("bPtHist", "B-meson p_{T} distribution", ptBinsSize, &ptBins[0]);
    
    double mass_min, mass_max, mass_peak;
    int nbins;
    TString ntuple_name = "", xaxis_title = "";
    
    switch (channel) {
        case 1:
            mass_min = 5.0; mass_max = 6.0;
            mass_peak = BP_MASS;
            nbins = 50;
            ntuple_name = "ntkp";
            xaxis_title = "M_{J/#psi K^{#pm}} [GeV]";
            break;
        case 2:
            mass_min = 5.0; mass_max = 6.0;
            mass_peak = B0_MASS;
            nbins = 50;
            ntuple_name = "ntkstar";
            xaxis_title = "M_{J/#psi K^{#pm}#pi^{#mp}} [GeV]";
            break;
        case 3:
            mass_min = 5.0; mass_max = 6.0;
            mass_peak = B0_MASS;
            nbins = 50;
            ntuple_name = "ntks";
            xaxis_title = "M_{J/#psi K^{0}_{S}} [GeV]";
            break;
        case 4:
            mass_min = 5.0; mass_max = 6.0;
            mass_peak = BS_MASS;
            nbins = 50;
            ntuple_name = "ntphi";
            xaxis_title = "M_{J/#psi K^{#pm}K^{#mp}} [GeV]";
            break;
        case 5:
            mass_min = 3.6; mass_max = 4.0;
            mass_peak = PSI2S_MASS;
            nbins = 80;
            ntuple_name = "ntmix";
            xaxis_title = "M_{J/#psi #pi^{#pm}#pi^{#mp}} [GeV]";
            break;
        case 6:
            mass_min = 5.3; mass_max = 6.3;
            mass_peak = LAMBDAB_MASS;
            nbins = 50;
            ntuple_name = "ntlambda";
            xaxis_title = "M_{J/#psi #Lambda} [GeV]";
            break;
    }
    
    TFile *fin = new TFile(SOURCE);
    TTree *tin = (TTree*)fin->Get(ntuple_name);
    
    ReducedBranches br;
    br.setbranchadd(tin);
    
    int n_br_queued = 0;
    ReducedBranches br_queue[32];
    
    for (int evt=0;evt<tin->GetEntries();evt++) {
        tin->GetEntry(evt);
        
        if (channel==1) { // cuts for B+ -> J/psi K+
            switch (beamSpotErrEstimate) {
                case 1:
                    if (br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v2]!=1) continue;
                    break;
                case 2:
                    if (br.hltbook[HLT_Dimuon16_Jpsi_v1]!=1) continue;
                    break;
                case 3:
                    if (br.hltbook[HLT_Dimuon0_Jpsi_Muon_v1]!=1) continue;
                    break;
                case 4:
                    if (br.hltbook[HLT_Dimuon10_Jpsi_Barrel_v1]!=1) continue;
                    break;
            }
            if (br.vtxprob<=0.1) continue;
            if (br.tk1pt<=1.0) continue;
            if (br.lxy/br.errxy<=3.0) continue;
            if (br.cosalpha2d<=0.99) continue;
            
            if (br.mass >= 5.16 && br.mass <= 5.365)
                bPtHist->Fill(br.pt);
            
        }else
        if (channel==2) { // cuts for B0 -> J/psi K*
            if (br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v2]!=1) continue;
            if (br.vtxprob<=0.1) continue;
            if (br.lxy/br.errxy<=3.0) continue;
            if (br.cosalpha2d<=0.99) continue;
            if (fabs(br.tktkmass-KSTAR_MASS)>=0.05) continue;
            
            TLorentzVector v4_tk1, v4_tk2;
            v4_tk1.SetPtEtaPhiM(br.tk1pt,br.tk1eta,br.tk1phi,KAON_MASS);
            v4_tk2.SetPtEtaPhiM(br.tk2pt,br.tk2eta,br.tk2phi,KAON_MASS);
            if (fabs((v4_tk1+v4_tk2).Mag()-PHI_MASS)<=0.01) continue;
            
            if (n_br_queued==0) {
                memcpy(&br_queue[n_br_queued],&br,sizeof(ReducedBranches));
                n_br_queued++;
            }else
            if (br.run == br_queue[n_br_queued-1].run && br.event == br_queue[n_br_queued-1].event) { // same event
                memcpy(&br_queue[n_br_queued],&br,sizeof(ReducedBranches));
                n_br_queued++;
                if (n_br_queued>=32) printf("Warning: maximum queued branches reached.\n");
            }
            
            if (br.run != br_queue[n_br_queued-1].run || br.event != br_queue[n_br_queued-1].event || evt==tin->GetEntries()-1) {
                for (int i=0; i<n_br_queued; i++) {
                    
                    bool isBestKstarMass = true;
                    for (int j=0; j<n_br_queued; j++) {
                        if (j==i) continue;
                        if (br_queue[i].mu1idx==br_queue[j].mu1idx &&
                            br_queue[i].mu2idx==br_queue[j].mu2idx &&
                            br_queue[i].tk1idx==br_queue[j].tk1idx &&
                            br_queue[i].tk2idx==br_queue[j].tk2idx) {
                        
                            if (fabs(br_queue[j].tktkmass-KSTAR_MASS)<fabs(br_queue[i].tktkmass-KSTAR_MASS)) {
                                isBestKstarMass = false;
                                continue;
                            }
                        }
                    }
                                 
                  //  if (isBestKstarMass) _nt->Fill(&br_queue[i].mass);
                }
                
                n_br_queued = 0;
                memcpy(&br_queue[n_br_queued],&br,sizeof(ReducedBranches));
                n_br_queued++;
            }
        }else
        if (channel==3) { // cuts for B0 -> J/psi Ks
            if (br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v2]!=1) continue;
            if (br.vtxprob<=0.1) continue;
            if (br.lxy/br.errxy<=3.0) continue;
            if (br.tktkblxy/br.tktkberrxy<=3.0) continue;
            if (br.cosalpha2d<=0.99) continue;
            if (fabs(br.tktkmass-KSHORT_MASS)>=0.015) continue;
                
           // _nt->Fill(&br.mass);
        }else
        if (channel==4) { // cuts for Bs -> J/psi phi
            if (br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v2]!=1) continue;
            if (br.vtxprob<=0.1) continue;
            if (br.lxy/br.errxy<=3.0) continue;
            if (br.cosalpha2d<=0.99) continue;
            if (fabs(br.tktkmass-PHI_MASS)>=0.010) continue;
            
            TLorentzVector v4_tk1, v4_tk2;
            v4_tk1.SetPtEtaPhiM(br.tk1pt,br.tk1eta,br.tk1phi,KAON_MASS);
            v4_tk2.SetPtEtaPhiM(br.tk2pt,br.tk2eta,br.tk2phi,PION_MASS);
            if (fabs((v4_tk1+v4_tk2).Mag()-KSTAR_MASS)<=0.05) continue;
            v4_tk1.SetPtEtaPhiM(br.tk1pt,br.tk1eta,br.tk1phi,PION_MASS);
            v4_tk2.SetPtEtaPhiM(br.tk2pt,br.tk2eta,br.tk2phi,KAON_MASS);
            if (fabs((v4_tk1+v4_tk2).Mag()-KSTAR_MASS)<=0.05) continue;
                
          //  _nt->Fill(&br.mass);
        }else
        if (channel==5) { // cuts for psi(2S)/X(3872) -> J/psi pipi
            if (br.vtxprob<=0.2) continue;
            if (fabs(br.tk1eta)>=1.6) continue;
            if (fabs(br.tk2eta)>=1.6) continue;
            
           // _nt->Fill(&br.mass);
        }else
        if (channel==6) {
            if (br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v2]!=1) continue;
            if (br.vtxprob<=0.1) continue;
            if (br.lxy/br.errxy<=3.0) continue;
            if (br.tktkblxy/br.tktkberrxy<=3.0) continue;
            if (br.cosalpha2d<=0.99) continue;
            if (fabs(br.tktkmass-LAMBDA_MASS)>=0.015) continue;
            
            TLorentzVector v4_tk1, v4_tk2;
            v4_tk1.SetPtEtaPhiM(br.tk1pt,br.tk1eta,br.tk1phi,PION_MASS);
            v4_tk2.SetPtEtaPhiM(br.tk2pt,br.tk2eta,br.tk2phi,PION_MASS);
            if (fabs((v4_tk1+v4_tk2).Mag()-KSHORT_MASS)<=0.015) continue;
            
            
          //  _nt->Fill(&br.mass);
        }
    }
    fin->Close();
    
    //------- Reading acceptance * efficiency numbers -----------------
    
    TFile *fileEffAcc = new TFile("./EffAcc.root","r");
    TH1D *h_ptEffAcc = (TH1D*)fileEffAcc->Get("h_bp_pt_AccEff");
    cout<<"Number of bins in acceptance*efficiency = "<<h_ptEffAcc->GetNbinsX()<<endl;
    
    //------- Published 7 TeV numbers----------------------------------
    
    double cross_section7TeV[] = {4.07, 1.47, 0.412, 0.181, 0.042};
   // double cross_section7TeVE[] = {0.47, 0.13, 0.041, 0.015, 0.007}; //Statistical errors
    double cross_section7TeVE[] = {0.31, 0.09, 0.026, 0.012, 0.004}; //Systematic errors
    double ptBins_7TeV[] = {5., 10., 13., 17., 24., 30.};
    int ptBinsSize_7TeV = (sizeof(ptBins_7TeV)/sizeof(*ptBins_7TeV))-1;
    cout<<"pT bin size (7 TeV) = "<<ptBinsSize_7TeV<<'\t'<<ptBins_7TeV[0]<<endl;
    
    TH1D *bPtHist_7TeV = new TH1D("bPtHist_7TeV", "7 TeV B-meson p_{T} distribution", ptBinsSize_7TeV, &ptBins_7TeV[0]);
    
    for(int j = 1; j <= bPtHist_7TeV->GetNbinsX(); j++)
    {
        cout<<j<<'\t'<<cross_section7TeV[j-1]<<endl;
        bPtHist_7TeV->SetBinContent(j,cross_section7TeV[j-1]);
        bPtHist_7TeV->SetBinError(j,cross_section7TeVE[j-1]);
    }
    
    //-----------------------------------------------------------------
    
    double etaMax = 2.4;
    double etaMin = -2.4;
    double etaWid = etaMax - etaMin;
    double nEvents = 739847.;
    double BR = (1.027e-03) * (5.961e-02);
    double Lint = 47; //47 (pb)-1
    
    bPtHist->Sumw2();
    bPtHist->Scale(1/nEvents);
    
    for(int i = 1; i <= bPtHist->GetNbinsX(); i++)
    {
        double binWid = bPtHist->GetBinWidth(i);
        double accEff = h_ptEffAcc->GetBinContent(i);
        
        bPtHist->SetBinContent(i, bPtHist->GetBinContent(i)/2.0/binWid/etaWid/BR/accEff/Lint);
        bPtHist->SetBinError(i, bPtHist->GetBinError(i)/2.0/binWid/etaWid/BR/accEff/Lint);
    }
    
#if DISPLAY
    TCanvas *c1 = canvasDressing("c1");
    c1->SetLogy();
    
    bPtHist->SetTitle("");
    bPtHist->GetXaxis()->SetTitle("p_{T}[GeV]");
    bPtHist->GetXaxis()->SetLabelFont(42);
    bPtHist->GetXaxis()->SetLabelOffset(0.01);
    bPtHist->GetXaxis()->SetTitleSize(0.06);
    bPtHist->GetXaxis()->SetTitleOffset(1.09);
    bPtHist->GetXaxis()->SetLabelFont(42);
    bPtHist->GetXaxis()->SetLabelSize(0.055);
    bPtHist->GetXaxis()->SetTitleFont(42);
    bPtHist->GetYaxis()->SetTitle("d#sigma/dp_{T} (pp #rightarrow B^{+}X; |y|<2.4) [#mub/GeV]");
    bPtHist->GetYaxis()->SetLabelFont(42);
    bPtHist->GetYaxis()->SetLabelOffset(0.01);
    bPtHist->GetYaxis()->SetTitleOffset(1.14);
    bPtHist->GetYaxis()->SetTitleSize(0.06);
    bPtHist->GetYaxis()->SetTitleFont(42);
    bPtHist->GetYaxis()->SetLabelFont(42);
    bPtHist->GetYaxis()->SetLabelSize(0.055);
    bPtHist->SetMarkerStyle(20);
    bPtHist->SetMarkerColor(kBlack);
    bPtHist->SetLineColor(kRed-7);
    bPtHist->SetLineWidth(4);
   // bPtHist->SetMaximum(2.0);
    bPtHist->Draw();
    
    bPtHist_7TeV->SetMarkerStyle(21);
    bPtHist_7TeV->SetMarkerColor(kBlack);
    bPtHist_7TeV->SetLineWidth(4);
    bPtHist_7TeV->SetLineColor(kYellow);
  //  bPtHist_7TeV->Draw("same");
  //  LegendpTSpectrum();
    
#endif
    
}
