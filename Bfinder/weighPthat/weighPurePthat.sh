#!/bin/bash

###

###pp
#Dzero prompt
#sed '1iconst int nBins=9; float pthatBin[nBins]={0,5,10,15,30,50,80,120,170}; float crosssec[nBins+1]={6.885e+09,1.516e+09,2.113e+08,5.524e+07,4.411e+06,5.986e+05,8.142e+04,1.329e+04,2.404e+03,0.}; int genSignal[2]={1,2};' weighPurePthat.C > weighPurePthat_tmp.C
#INPUTFILE="ntD_EvtBase_20160513_DfinderMC_pp_20160502_dPt0tkPt0p5_D0Dstar_prompt_Dpt2Dy1p1tkPt0p7tkEta2Decay2p9Dalpha0p14Skim.root"
#OUTPUTFILE="ntD_EvtBase_20160513_DfinderMC_pp_20160502_dPt0tkPt0p5_D0Dstar_prompt_Dpt2Dy1p1tkPt0p7tkEta2Decay2p9Dalpha0p14Skim_pthatweight.root"

#Dzero nonprompt
#sed '1iconst int nBins=9; float pthatBin[nBins]={0,5,10,15,30,50,80,120,170}; float crosssec[nBins+1]={1.895e+08,1.393e+08,3.798e+07,1.206e+07,1.231e+06,1.853e+05,2.725e+04,4.566e+03,8.569e+02,0.}; int genSignal[2]={1,2};' weighPurePthat.C > weighPurePthat_tmp.C
#INPUTFILE="/data/wangj/MC2015/Dntuple/pp/revised/ntD_pp_Dzero_kpi_nonprompt/ntD_EvtBase_20160303_Dfinder_20160302_pp_Pythia8_nonprompt_D0_dPt0tkPt0p5_noweight.root"
#OUTPUTFILE="/data/wangj/MC2015/Dntuple/pp/revised/ntD_pp_Dzero_kpi_nonprompt/ntD_EvtBase_20160303_Dfinder_20160302_pp_Pythia8_nonprompt_D0_dPt0tkPt0p5_pthatweight.root"

###PbPb
#Dzero prompt
sed '1iconst int nBins=9; float pthatBin[nBins]={0,5,10,15,30,50,80,120,170}; float crosssec[nBins+1]={6.885e+09,1.516e+09,2.113e+08,5.524e+07,4.411e+06,5.986e+05,8.142e+04,1.329e+04,2.404e+03,0.}; int genSignal[2]={1,2};' weighPurePthat.C > weighPurePthat_tmp.C
#INPUTFILE="/data/wangj/MC2015/Dntuple/PbPb/revised/ntD_PbPb_Dzero_kpi_prompt/ntD_EvtBase_20160330_Dfinder_20160329_PbPb_Pythia8_prompt_D0_dPt1tkPt0p5_noweight.root"
INPUTFILE="ntD_EvtBase_20160513_DfinderMC_PbPb_20160502_dPt1tkPt0p5_D0_prompt_Dpt2Dy1p1tkPt0p7tkEta2Decay2p9Dalpha0p14Skim.root"
#OUTPUTFILE="/data/wangj/MC2015/Dntuple/PbPb/revised/ntD_PbPb_Dzero_kpi_prompt/ntD_EvtBase_20160330_Dfinder_20160329_PbPb_Pythia8_prompt_D0_dPt1tkPt0p5_pthatweight.root"
OUTPUTFILE="ntD_EvtBase_20160513_DfinderMC_PbPb_20160502_dPt1tkPt0p5_D0_prompt_Dpt2Dy1p1tkPt0p7tkEta2Decay2p9Dalpha0p14Skim_pthatweight.root"

#Dzero nonprompt
#sed '1iconst int nBins=9; float pthatBin[nBins]={0,5,10,15,30,50,80,120,170}; float crosssec[nBins+1]={1.895e+08,1.393e+08,3.798e+07,1.206e+07,1.231e+06,1.853e+05,2.725e+04,4.566e+03,8.569e+02,0.}; int genSignal[2]={1,2};' weighPurePthat.C > weighPurePthat_tmp.C
#INPUTFILE="/data/wangj/MC2015/Dntuple/PbPb/revised/ntD_PbPb_Dzero_kpi_nonprompt/ntD_EvtBase_20160330_Dfinder_20160329_PbPb_Pythia8_nonprompt_D0_dPt1tkPt0p5_noweight.root"
#OUTPUTFILE="/data/wangj/MC2015/Dntuple/PbPb/revised/ntD_PbPb_Dzero_kpi_nonprompt/ntD_EvtBase_20160330_Dfinder_20160329_PbPb_Pythia8_nonprompt_D0_dPt1tkPt0p5_pthatweight.root"

###

cp "$INPUTFILE" "$OUTPUTFILE"
g++ weighPurePthat_tmp.C $(root-config --cflags --libs) -g -o weighPurePthat_tmp.exe 
./weighPurePthat_tmp.exe "$INPUTFILE" "$OUTPUTFILE"
rm weighPurePthat_tmp.exe

rm weighPurePthat_tmp.C
