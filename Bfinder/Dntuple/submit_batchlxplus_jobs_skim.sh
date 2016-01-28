####REMEMBER TO MOVE TO PP OR PbPb MODE
DIRECTORYOUTPUT="/afs/cern.ch/work/g/ginnocen/HeavyFlavour/ProductionsBfinder/CMSSW_7_5_8/src/Bfinder/Bfinder/Dntuple/output"
NAME="loop.C"
g++ $NAME $(root-config --cflags --libs) -Wall -O2 -o "${NAME/%.C/}.exe"
NAMESKIM="skim.C"
g++ $NAMESKIM $(root-config --cflags --libs) -Wall -O2 -o "${NAMESKIM/%.C/}.exe"

rm -rf $DIRECTORYOUTPUT
mkdir $DIRECTORYOUTPUT
rm  mylistfinal.txt

eos ls /store/group/phys_heavyions/HeavyFlavourRun2/DfinderRun2/HIHardProbes/crab_DfinderData_PbPb_20160126_dPt0tkPt2p5_D0Dstar3p5p_FINALJSON_v6/160126_215426/0000/  | awk '{print "root://eoscms//eos/cms/store/group/phys_heavyions/HeavyFlavourRun2/DfinderRun2/HIHardProbes/crab_DfinderData_PbPb_20160126_dPt0tkPt2p5_D0Dstar3p5p_FINALJSON_v6/160126_215426/0000/" $0}' >> mylistfinal.txt
eos ls /store/group/phys_heavyions/HeavyFlavourRun2/DfinderRun2/HIHardProbes/crab_DfinderData_PbPb_20160126_dPt0tkPt2p5_D0Dstar3p5p_FINALJSON_v6/160126_215426/0001/  | awk '{print "root://eoscms//eos/cms/store/group/phys_heavyions/HeavyFlavourRun2/DfinderRun2/HIHardProbes/crab_DfinderData_PbPb_20160126_dPt0tkPt2p5_D0Dstar3p5p_FINALJSON_v6/160126_215426/0001/" $0}' >> mylistfinal.txt
eos ls /store/group/phys_heavyions/HeavyFlavourRun2/DfinderRun2/HIHardProbes/crab_DfinderData_PbPb_20160126_dPt0tkPt2p5_D0Dstar3p5p_FINALJSON_v6/160126_215426/0002/  | awk '{print "root://eoscms//eos/cms/store/group/phys_heavyions/HeavyFlavourRun2/DfinderRun2/HIHardProbes/crab_DfinderData_PbPb_20160126_dPt0tkPt2p5_D0Dstar3p5p_FINALJSON_v6/160126_215426/0002/" $0}' >> mylistfinal.txt

count=0 ; for i in `cat mylistfinal.txt` ; do echo cd $PWD/ > run_$count.sh ; echo 'export X509_USER_PROXY=~/x509_user_proxy'>> run_$count.sh ; 
echo eval \`scram runtime -sh\` >> run_$count.sh ; echo cd - >> run_$count.sh ; echo $PWD/loop.exe $i $DIRECTORYOUTPUT/ntuple_$count.root  >> run_$count.sh ;  echo $PWD/skim.exe $DIRECTORYOUTPUT/ntuple_$count.root $DIRECTORYOUTPUT/ntuple_skimmed_$count.root  >> run_$count.sh ; chmod +x run_$count.sh; count=$((count+1)) ; done
for i in `ls run_*` ; do bsub -q cmscaf1nd $i ; done

#cmscaf1nd
