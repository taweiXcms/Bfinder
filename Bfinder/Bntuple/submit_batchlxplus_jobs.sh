####REMEMBER TO MOVE TO PP OR PbPb MODE
DIRECTORYOUTPUT="/afs/cern.ch/work/g/ginnocen/HeavyFlavourRun2/BfinderData_pp_20151218_bPt0jpsiPt0tkPt0p5_v3"
NAME="loop.C"
g++ $NAME $(root-config --cflags --libs) -Wall -O2 -o "${NAME/%.C/}.exe"

rm  mylistfinal.txt

eos ls /store/group/phys_heavyions/HeavyFlavourRun2/BfinderRun2/DoubleMu/crab_BfinderData_pp_20151218_bPt0jpsiPt0tkPt0p5_v3/151218_121602/0000/  | awk '{print "root://eoscms//eos/cms/store/group/phys_heavyions/HeavyFlavourRun2/BfinderRun2/DoubleMu/crab_BfinderData_pp_20151218_bPt0jpsiPt0tkPt0p5_v3/151218_121602/0000/" $0}' >> mylistfinal.txt
eos ls /store/group/phys_heavyions/HeavyFlavourRun2/BfinderRun2/DoubleMu/crab_BfinderData_pp_20151218_bPt0jpsiPt0tkPt0p5_v3/151218_121602/0001/  | awk '{print "root://eoscms//eos/cms/store/group/phys_heavyions/HeavyFlavourRun2/BfinderRun2/DoubleMu/crab_BfinderData_pp_20151218_bPt0jpsiPt0tkPt0p5_v3/151218_121602/0001/" $0}' >> mylistfinal.txt
eos ls /store/group/phys_heavyions/HeavyFlavourRun2/BfinderRun2/DoubleMu/crab_BfinderData_pp_20151218_bPt0jpsiPt0tkPt0p5_v3/151218_121602/0002/  | awk '{print "root://eoscms//eos/cms/store/group/phys_heavyions/HeavyFlavourRun2/BfinderRun2/DoubleMu/crab_BfinderData_pp_20151218_bPt0jpsiPt0tkPt0p5_v3/151218_121602/0002/" $0}' >> mylistfinal.txt

count=0 ; for i in `cat mylistfinal.txt` ; do echo cd $PWD/ > run_$count.sh ; echo 'export X509_USER_PROXY=~/x509_user_proxy'>> run_$count.sh ; 
echo eval \`scram runtime -sh\` >> run_$count.sh ; echo cd - >> run_$count.sh ; echo $PWD/loop.exe $i $DIRECTORYOUTPUT/ntuple_$count.root  >> run_$count.sh ; chmod +x run_$count.sh; count=$((count+1)) ; done
for i in `ls run_*` ; do bsub -q cmscaf1nd $i ; done

#cmscaf1nd
