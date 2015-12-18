####REMEMBER TO MOVE TO PP OR PbPb MODE
DIRECTORYOUTPUT="/afs/cern.ch/work/g/ginnocen/HeavyFlavourRun2/DfinderData_pp_20151218_dPt0tkPt1_D0Dstar3p5p"
NAME="loop.C"
g++ $NAME $(root-config --cflags --libs) -Wall -O2 -o "${NAME/%.C/}.exe"

rm  mylistfinal.txt

eos ls /store/group/phys_heavyions/HeavyFlavourRun2/DfinderRun2/HeavyFlavor/crab_DfinderData_pp_20151218_dPt0tkPt1_D0Dstar3p5p/151218_092138/0000/  | awk '{print "root://eoscms//eos/cms/store/group/phys_heavyions/HeavyFlavourRun2/DfinderRun2/HeavyFlavor/crab_DfinderData_pp_20151218_dPt0tkPt1_D0Dstar3p5p/151218_092138/0000/" $0}' >> mylistfinal.txt

count=0 ; for i in `cat mylistfinal.txt` ; do echo cd $PWD/ > run_$count.sh ; echo 'export X509_USER_PROXY=~/x509_user_proxy'>> run_$count.sh ; 
echo eval \`scram runtime -sh\` >> run_$count.sh ; echo cd - >> run_$count.sh ; echo $PWD/loop.exe $i $DIRECTORYOUTPUT/ntuple_$count.root  >> run_$count.sh ; chmod +x run_$count.sh; count=$((count+1)) ; done
for i in `ls run_*` ; do bsub -q cmscaf1nd $i ; done

#cmscaf1nd
