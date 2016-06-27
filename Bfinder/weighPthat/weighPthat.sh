#!/bin/bash

###

###pp
#Bplus
sed '1iconst int nBins=4; float pthatBin[nBins]={5,15,30,50}; float crosssec[nBins+1]={3.845e+07,1.317e+06,1.406e+05,1.446e+04,0.}; int genSignal[2]={1,1};' weighPthat.C > weighPthat_tmp.C
#INPUTFILE="Bntuple20160624_Pythia8_BuToJpsiK.root"
#OUTPUTFILE="Bntuple20160624_Pythia8_BuToJpsiK_pthatweight.root"
INPUTFILE="Bntuple20160624_pp_Pythia8_BuToJpsiK.root"
OUTPUTFILE="Bntuple20160624_pp_Pythia8_BuToJpsiK_pthatweight.root"

###

cp "$INPUTFILE" "$OUTPUTFILE"
g++ weighPthat_tmp.C $(root-config --cflags --libs) -g -o weighPthat_tmp.exe 
./weighPthat_tmp.exe "$INPUTFILE" "$OUTPUTFILE"
rm weighPthat_tmp.exe

rm weighPthat_tmp.C
