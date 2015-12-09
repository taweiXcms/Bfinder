#!/bin/bash
#path="/eos/cms/tier0/store/hidata/HIRun2015/HIOniaL1DoubleMu0/AOD/PromptReco-v1/000/262"
#path="/eos/cms/tier0/store/hidata/HIRun2015/HIHardProbes/AOD/PromptReco-v1/000/262/"
#path="/eos/cms/tier0/store/hidata/HIRun2015/HIMinimumBias1/AOD/PromptReco-v1/000/262/"
eos="/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select"
dir=`$eos ls $path`
#echo $dir 
for d in $dir
do
	subdir=`$eos ls $path/$d`
	#echo $d
	for sd in $subdir
	do
		#echo $sd
		file=`$eos ls $path/$d/$sd`
		for f in $file
		do
			echo "root://cms-xrd-global.cern.ch/$path/$d/$sd/$f"
		done 
	done 
done

#num="548 563 566 570 694 695 703 726 735"
