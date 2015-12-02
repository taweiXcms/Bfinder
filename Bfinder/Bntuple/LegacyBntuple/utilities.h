#include <iostream>
#include <utility>
#include <TRandom.h>
#include <TRandom3.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TProfile.h>
#include <TCut.h>
#include <TH2F.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TBox.h>
#include <TF1.h>
#include <TDirectory.h>
#include <TDirectoryFile.h>
#include <TGraph.h>
#include <TMath.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TNtuple.h>
#include <TGraphAsymmErrors.h>

static const int nColor = 5;
static const int colorCode[nColor] = {
	1, 2, kGreen+1, 4, 6
};

// Algos
static const int nAlgos = 8;
//static const char *algoName[nAlgos] = { "", "icPu5", "akPu2PF", "akPu3PF", "akPu4PF", "akPu2Calo", "akPu3Calo", "akPu4Calo" };
//static const char *algoNamePP[nAlgos] = { "", "icPu5", "ak2PF", "ak3PF", "ak4PF", "ak2Calo", "ak3Calo", "ak4Calo" };
//static const char *algoNameGen[nAlgos] = { "", "icPu5", "akPu2PF", "akPu3PF", "akPu4PF", "akPu2PF", "akPu3PF", "akPu4PF" };

// Centrality binning
const int nbins_cent=	     6;
Double_t boundaries_cent[nbins_cent+1] = {   0,2,4,12,20,28,36 };
Double_t ncoll[nbins_cent] = { 1660, 1310, 745, 251, 62.8, 10.8 };

// Track pT binning
const int nbins_trackPt=	     22;
Double_t boundaries_trackPt[nbins_trackPt+1] = {   0.5,0.6,0.7,0.8,0.9, 1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9, 2.0, 4.0, 8.0, 20.0, 30.0, 50.0, 100.0 , 200.0 };



// Missing pT binning
const int nbins_MPT=	     6;
Double_t boundaries_MPT[nbins_MPT+1] = {   0,2,4,12,20,28,36   };

// Remove error 
void removeError(TH1F *h)
{
	for (int i=1;i<=h->GetNbinsX();i++)
	{
		h->SetBinError(i,0);
	}   
	
}


class JetData
{
   public:
   int hiBin;
   float hiHFplus;
   float leadingJetPt;  
   float leadingJetEta;  
   float leadingJetPhi;  
   float subleadingJetPt;
   float subleadingJetEta;
   float subleadingJetPhi;
   float genleadingJetPt;  
   float genleadingJetEta;  
   float genleadingJetPhi;  
   float gensubleadingJetPt;
   float gensubleadingJetEta;
   float gensubleadingJetPhi;
   float mpt;
   float cormpt;
   float cormpt2;
   float genPMpt;
   float genMpt;
   float v2;
   float phi0;
   float n;
   float corV2;
   float corPhi0;
   float corN;
   int    leadingJetIt;  
   int    subleadingJetIt;
   
   
   JetData(TTree *t, bool setBranch=1) {
      if (setBranch){
      t->Branch("hiBin",&hiBin,"hiBin/I");
      t->Branch("hiHFplus",&hiHFplus,"hiHFplus/F");
      t->Branch("v2",&v2,"v2/F");
      t->Branch("phi0",&phi0,"phi0/F");
      t->Branch("n",&n,"n/F");
      t->Branch("corV2",&corV2,"corV2/F");
      t->Branch("corPhi0",&corPhi0,"corPhi0/F");
      t->Branch("corN",&corN,"corN/F");
      t->Branch("leadingJetPt",&leadingJetPt,"leadingJetPt/F");
      t->Branch("leadingJetPhi",&leadingJetPhi,"leadingJetPhi/F");
      t->Branch("leadingJetEta",&leadingJetEta,"leadingJetEta/F");
      t->Branch("subleadingJetPt",&subleadingJetPt,"subleadingJetPt/F");
      t->Branch("subleadingJetPhi",&subleadingJetPhi,"subleadingJetPhi/F");
      t->Branch("subleadingJetEta",&subleadingJetEta,"subleadingJetEta/F");
      t->Branch("genleadingJetPt",&genleadingJetPt,"genleadingJetPt/F");
      t->Branch("genleadingJetPhi",&genleadingJetPhi,"genleadingJetPhi/F");
      t->Branch("genleadingJetEta",&genleadingJetEta,"genleadingJetEta/F");
      t->Branch("gensubleadingJetPt",&gensubleadingJetPt,"gensubleadingJetPt/F");
      t->Branch("gensubleadingJetPhi",&gensubleadingJetPhi,"gensubleadingJetPhi/F");
      t->Branch("gensubleadingJetEta",&gensubleadingJetEta,"gensubleadingJetEta/F");
      t->Branch("leadingJetIt",&leadingJetIt,"leadingJetIt/I");
      t->Branch("subleadingJetIt",&subleadingJetIt,"subleadingJetIt/I");
      t->Branch("mpt",&mpt,"mpt/F");
      t->Branch("cormpt",&cormpt,"cormpt/F");
      t->Branch("cormpt2",&cormpt2,"cormpt2/F");
      t->Branch("genMpt",&genMpt,"genMpt/F");
      t->Branch("genPMpt",&genPMpt,"genPMpt/F");
      t->SetAlias("Aj","(leadingJetPt-subleadingJetPt)/(leadingJetPt+subleadingJetPt)");
      t->SetAlias("genAj","(genleadingJetPt-gensubleadingJetPt)/(genleadingJetPt+gensubleadingJetPt)");
      }
   };
};

class HistoData
{
   public:
   TH1F *hRecoPt;
   TH1F *hCorrectedPt;
   TH1F *hGenPt;
   TH1F *hClosurePt;   

   TH1F *hRecoEta;
   TH1F *hCorrectedEta;
   TH1F *hGenEta;
   TH1F *hClosureEta;   

   TH1F *hRecoPhi;
   TH1F *hCorrectedPhi;
   TH1F *hGenPhi;
   TH1F *hClosurePhi;   

   TH1F *hRecoDR;
   TH1F *hCorrectedDR;
   TH1F *hGenDR;
   TH1F *hClosureDR;   


   char *histoName;
   HistoData(char *title) {
      hRecoPt       = new TH1F (Form("hRecoPt_%s",title),";Reconstructed Tracks p_{T} (GeV/c); Entries",nbins_trackPt,boundaries_trackPt);
      hGenPt        = new TH1F (Form("hGenPt_%s",title),";SimTracks p_{T} (GeV/c); Entries",nbins_trackPt,boundaries_trackPt);
      hCorrectedPt  = new TH1F (Form("hCorrectedPt_%s",title),";Eff Corrected Tracks p_{T} (GeV/c); Entries",nbins_trackPt,boundaries_trackPt);
      hRecoPt->Sumw2();
      hGenPt->Sumw2();
      hCorrectedPt->Sumw2();

      hRecoEta       = new TH1F (Form("hRecoEta_%s",title),";Reconstructed Tracks #eta; Entries",48,-2.4,2.4);
      hGenEta        = new TH1F (Form("hGenEta_%s",title),";SimTracks #eta; Entries",48,-2.4,2.4);
      hCorrectedEta  = new TH1F (Form("hCorrectedEta_%s",title),";Eff Corrected Tracks #eta; Entries",48,-2.4,2.4);
      hRecoEta->Sumw2();
      hGenEta->Sumw2();
      hCorrectedEta->Sumw2();

      hRecoPhi       = new TH1F (Form("hRecoPhi_%s",title),";Reconstructed Tracks #Phi; Entries",20,-3.1416,3.1416);
      hGenPhi        = new TH1F (Form("hGenPhi_%s",title),";SimTracks #Phi; Entries",20,-3.1416,3.1416);
      hCorrectedPhi  = new TH1F (Form("hCorrectedPhi_%s",title),";Eff Corrected Tracks #Phi; Entries",20,-3.1416,3.1416);
      hRecoPhi->Sumw2();
      hGenPhi->Sumw2();
      hCorrectedPhi->Sumw2();

      hRecoDR       = new TH1F (Form("hRecoDR_%s",title),";Reconstructed Tracks #DR; Entries",50,0,5);
      hGenDR        = new TH1F (Form("hGenDR_%s",title),";SimTracks #DR; Entries",50,0,5);
      hCorrectedDR  = new TH1F (Form("hCorrectedDR_%s",title),";Eff Corrected Tracks #DR; Entries",50,0,5);
      hRecoDR->Sumw2();
      hGenDR->Sumw2();
      hCorrectedDR->Sumw2();


      histoName = title;
   } ;  
   
   void calcEff() {
      hClosurePt = (TH1F*) hCorrectedPt->Clone(Form("hClosurePt_%s",histoName));
      removeError(hGenPt);
      hClosurePt->Divide(hGenPt);
      hClosurePt->SetYTitle("Ratio");
      hClosureEta = (TH1F*) hCorrectedEta->Clone(Form("hClosureEta_%s",histoName));
      removeError(hGenEta);
      hClosureEta->Divide(hGenEta);
      hClosureEta->SetYTitle("Ratio");
      hClosurePhi = (TH1F*) hCorrectedPhi->Clone(Form("hClosurePhi_%s",histoName));
      removeError(hGenPhi);
      hClosurePhi->Divide(hGenPhi);
      hClosurePhi->SetYTitle("Ratio");
      hClosureDR = (TH1F*) hCorrectedDR->Clone(Form("hClosureDR_%s",histoName));
      removeError(hGenDR);
      hClosureDR->Divide(hGenDR);
      hClosureDR->SetYTitle("Ratio");
   }
};



// make a histogram from TF1 function
TH1F *functionHist(TF1 *f, TH1F* h,char *fHistname)
{
	TH1F *hF = (TH1F*)h->Clone(fHistname);
	for (int i=1;i<=h->GetNbinsX();i++)
	{
		double var = f->Integral(h->GetBinLowEdge(i),h->GetBinLowEdge(i+1))/h->GetBinWidth(i);
		hF->SetBinContent(i,var);
		hF->SetBinError(i,0);
	}
	return hF;
}

TLegend *myLegend(double x1,double y1,double x2, double y2)
{
	TLegend *leg = new TLegend(x1,y1,x2,y2);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	return leg; 
	
}


// draw envelope using a systematic uncertainty histogram
TH1F* drawEnvelope(TH1F *h,char *opt,int color = kGray,int fillStyle = 0, int fillColor = 0,double shift = 0)
{
	TH1F *hClone = (TH1F*) h->Clone(Form("%s_mirror",h->GetTitle()));
	TH1F *hMirror = (TH1F*) h->Clone(Form("%s_mirror2",h->GetTitle()));
	for (int i=1;i<=h->GetNbinsX();i++)
	{
		double val = h->GetBinContent(i);
		hMirror->SetBinContent(i,1-fabs(val-1)+shift);
		hClone->SetBinContent(i,val+shift);
	}
	
	//   hMirror->SetLineStyle(2);
	//   h->SetLineStyle(2);
	hMirror->SetLineColor(color);
	hMirror->SetFillColor(fillColor);
	hMirror->SetFillStyle(fillStyle);
	hClone->SetLineColor(color);
	hClone->SetFillColor(fillColor);
	hClone->SetFillStyle(fillStyle);
	hClone->Draw(opt);
	hMirror->Draw(opt);
	return hMirror;
}

void makeHistTitle(TH1 *h,char *title, char *xTitle, char *yTitle, int color = -1, bool centerTitle = 1)
{
	h->SetTitle(title);
	h->SetXTitle(xTitle);
	h->SetYTitle(yTitle);
	
	if (centerTitle) {
		h->GetXaxis()->CenterTitle();
		h->GetYaxis()->CenterTitle();
		
	}
	
	if (color!=-1) {
		h->SetLineColor(color);
		h->SetMarkerColor(color);
	}
	
	h->GetYaxis()->SetNdivisions(610); 
	
	h->GetYaxis()->SetLabelFont(43);
	h->GetYaxis()->SetTitleFont(43);
	h->GetYaxis()->SetLabelSize(20);
	h->GetYaxis()->SetTitleSize(22);
	h->GetYaxis()->SetTitleOffset(1.6);
	
	h->GetXaxis()->SetLabelFont(43);
	h->GetXaxis()->SetTitleFont(43);
	h->GetXaxis()->SetLabelSize(20);
	h->GetXaxis()->SetTitleSize(22);
	h->GetXaxis()->SetTitleOffset(3.1);
	
	h->GetXaxis()->SetNoExponent();
	h->GetXaxis()->SetMoreLogLabels();
	
	h->GetXaxis()->SetTitleOffset(1.4);
	
}


// Remove bins with error > central value
void cleanup(TH1F *h)
{
	for (int i=1;i<=h->GetNbinsX();i++)
	{
		double val1 = h->GetBinContent(i);
		double valErr1 = h->GetBinError(i);
		if (valErr1>=val1) {
			h->SetBinContent(i,0);
			h->SetBinError(i,0);
		}
	}   
	
}

/*
// rebin the spectra
TH1F *rebin(TH1F *h, char *histName)
{
	TH1F *hRebin = new TH1F(Form("%s_rebin",h->GetName()),Form("rebin %s",h->GetTitle()),nbins_recrebin,boundaries_recrebin);
	for (int i=1;i<=h->GetNbinsX();i++)
	{
		double val=h->GetBinContent(i);
		double valErr=h->GetBinError(i);
		int binNum = hRebin->FindBin(h->GetBinCenter(i));
		double val1 = hRebin->GetBinContent(binNum);
		double valErr1 = hRebin->GetBinError(binNum);
		hRebin->SetBinContent(binNum,val+val1);
		hRebin->SetBinError(binNum,sqrt(valErr1*valErr1+valErr*valErr));
	}
	cleanup(hRebin);
	hRebin->SetName(histName);
	return hRebin;
}*/

// divide by bin width
void divideBinWidth(TH1 *h)
{
	h->Sumw2();
	for (int i=0;i<=h->GetNbinsX();i++)
	{
		Float_t val = h->GetBinContent(i);
		Float_t valErr = h->GetBinError(i);
		val/=h->GetBinWidth(i);
		valErr/=h->GetBinWidth(i);
		h->SetBinContent(i,val);
		h->SetBinError(i,valErr);
	}
	h->GetXaxis()->CenterTitle();
	h->GetYaxis()->CenterTitle();
}


// make systematic histogram
void checkMaximumSys(TH1F *hSys, TH1F *h, int opt=0)
{
	if (h->GetNbinsX()!=hSys->GetNbinsX()) {
		cout <<"ERROR! Different NBins in subroutine checkMaximumSys!"<<endl;
	} else {
		double val = 1;
		for (int i=1;i<=h->GetNbinsX();i++) {
			//cout <<i<<" "<<val<<" "<<hSys->GetBinContent(i)<<" "<<h->GetBinContent(i)<<endl;
			if (h->GetBinContent(i)==0) continue;
			if (fabs(hSys->GetBinContent(i))>val) val = fabs(hSys->GetBinContent(i));
			if (fabs(h->GetBinContent(i)-1)+1>val) val=fabs(h->GetBinContent(i)-1)+1;
			if (opt) hSys->SetBinContent(i,val);
			else     hSys->SetBinContent(i,fabs(h->GetBinContent(i)-1)+1);
		}
	}
}



void makeMultiPanelCanvasWithGap(TCanvas*& canv,
								 const Int_t columns,
								 const Int_t rows,
								 const Float_t leftOffset,
								 const Float_t bottomOffset,
								 const Float_t leftMargin,
								 const Float_t bottomMargin,
								 const Float_t edge, const Float_t asyoffset) {
	if (canv==0) {
		Error("makeMultiPanelCanvas","Got null canvas.");
		return;
	}
	canv->Clear();
	
	TPad* pad[columns][rows];
	
	Float_t Xlow[columns];
	Float_t Xup[columns];
	Float_t Ylow[rows];
	Float_t Yup[rows];
	
	Float_t PadWidth =
	(1.0-leftOffset)/((1.0/(1.0-leftMargin)) +
					  (1.0/(1.0-edge))+(Float_t)columns-2.0);
	Float_t PadHeight =
	(1.0-bottomOffset)/((1.0/(1.0-bottomMargin)) +
						(1.0/(1.0-edge))+(Float_t)rows-2.0);
	
	//PadHeight = 0.5*PadHeight;
	
	Xlow[0] = leftOffset;
	Xup[0] = leftOffset + PadWidth/(1.0-leftMargin);
	Xup[columns-1] = 1;
	Xlow[columns-1] = 1.0-PadWidth/(1.0-edge);
	
	Yup[0] = 1;
	Ylow[0] = 1.0-PadHeight/(1.0-edge);
	Ylow[rows-1] = bottomOffset;
	Yup[rows-1] = bottomOffset + PadHeight/(1.0-bottomMargin);
	
	for(Int_t i=1;i<columns-1;i++) {
		Xlow[i] = Xup[0] + (i-1)*PadWidth;
		Xup[i] = Xup[0] + (i)*PadWidth;
	}
	Int_t ct = 0;
	for(Int_t i=rows-2;i>0;i--) {
		if(i==rows-2){
			Ylow[i] = Yup[rows-1] + ct*PadHeight;
			Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
		}else{
			Ylow[i] = Yup[rows-1] + ct*PadHeight;
			Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
			//Yup[i] = 0.2*Yup[i];
		}
		ct++;
	}
	
	TString padName;
	for(Int_t i=0;i<columns;i++) {
		for(Int_t j=0;j<rows;j++) {
			canv->cd();
			padName = Form("p_%d_%d",i,j);
			//pad[i][j] = new TPad(padName.Data(),padName.Data(),
			//Xlow[i],Ylow[j],Xup[i],Yup[j]);
			// this is hacked version to create aysmmetric pads around low 
			if(j==0){
				pad[i][j] = new TPad(padName.Data(),padName.Data(),
									 Xlow[i],Ylow[j]-asyoffset,Xup[i],Yup[j]);
			}else{
				pad[i][j] = new TPad(padName.Data(),padName.Data(),
									 Xlow[i],Ylow[j],Xup[i],Yup[j]-asyoffset);
			}
			
			
			if(i==0) pad[i][j]->SetLeftMargin(leftMargin);
			else pad[i][j]->SetLeftMargin(0);
			
			if(i==(columns-1)) pad[i][j]->SetRightMargin(edge);
			else pad[i][j]->SetRightMargin(0);
			
			if(j==0) pad[i][j]->SetTopMargin(edge);
			//else pad[i][j]->SetTopMargin(0.01);
			else pad[i][j]->SetTopMargin(0.02);
			
			//if(j==0) pad[i][j]->SetTopMargin(edge*3.5);
			//else pad[i][j]->SetTopMargin(0.0);
			
			if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin);
			else pad[i][j]->SetBottomMargin(0.15);
			
			pad[i][j]->Draw();
			pad[i][j]->cd();
			pad[i][j]->SetNumber(columns*j+i+1);
		}
	}
}

void putCMSPrel(double x, double y, double size){
	TLatex *tex=0;
	tex = new TLatex(x,y,"CMS Preliminary");
	tex->SetTextSize(size);
	tex->SetLineWidth(2);
	tex->SetNDC();
	tex->Draw();
}
void drawText(const char *text, float xp, float yp, int size){
	TLatex *tex = new TLatex(xp,yp,text);
	tex->SetTextFont(63);
	tex->SetTextSize(size);
	tex->SetTextColor(kBlack);
	tex->SetLineWidth(1);
	//tex->SetTextFont(42);
	tex->SetNDC();
	tex->Draw();
}


TGraphAsymmErrors *getEfficiency(TH1 *hPass,TH1 *hTotal)
{
   TGraphAsymmErrors *g = new TGraphAsymmErrors;
   g->BayesDivide(hPass,hTotal);
   return g;
}

TGraphAsymmErrors *getEfficiency(TTree *t,char *variable,int nBin, double binL, double binH, TCut selection, TCut preselection, TCut weight = "1")
{
   TH1F *hPass = new TH1F("hPass","",nBin,binL,binH);
   TH1F *hTotal = new TH1F("hTotal","",nBin,binL,binH);
   hPass->Sumw2();
   hTotal->Sumw2();
   cout <<Form("%s*(%s)",weight.GetTitle(),(selection && preselection).GetTitle())<<endl;
   t->Project("hPass",variable, Form("%s*(%s)",weight.GetTitle(),(selection && preselection).GetTitle()));
   t->Project("hTotal",variable,Form("%s",preselection.GetTitle()) );
   
   TGraphAsymmErrors *g = getEfficiency(hPass,hTotal);
   
//   hPass->Divide(hPass,hTotal,1,1,"B");
   
//   TGraphAsymmErrors *g = new TGraphAsymmErrors(hPass);
   delete hPass;
   delete hTotal;
   return g;
}
