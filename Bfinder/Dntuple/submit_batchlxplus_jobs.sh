DIRECTORYOUTPUT="root://eoscms//eos/cms//store/group/phys_heavyions/wangj/DntupleData/ntD_EvtBase_20160425_HighPtJet80_DfinderData_pp_20160329_dPt0tkPt0p5_D0Dstar"
NAME="loop.C"

rm mylistfinal.txt
rm run_*.sh
rm -rf LSFJOB_*

mkdir -p $DIRECTORYOUTPUT

g++ $NAME $(root-config --cflags --libs) -Wall -O2 -o "${NAME/%.C/}.exe"

eos ls /store/group/phys_heavyions/HeavyFlavourRun2/DfinderRun2/HighPtJet80/DfinderData_pp_20160329_dPt0tkPt0p5_D0Dstar/  | awk '{print "root://eoscms//eos/cms//store/group/phys_heavyions/HeavyFlavourRun2/DfinderRun2/HighPtJet80/DfinderData_pp_20160329_dPt0tkPt0p5_D0Dstar/" $0}' >> mylistfinal.txt

count=0 ; for i in `cat mylistfinal.txt` ; do echo cd $PWD/ > run_$count.sh ; echo 'export X509_USER_PROXY=~/x509_user_proxy'>> run_$count.sh ; 
echo eval \`scram runtime -sh\` >> run_$count.sh ; echo cd - >> run_$count.sh ; echo $PWD/loop.exe $i $DIRECTORYOUTPUT/ntuple_$count.root  >> run_$count.sh ; chmod +x run_$count.sh; count=$((count+1)) ; done
#for i in `ls run_*` ; do bsub -q cmscaf1nd $i ; done
for i in `ls run_*` ; do bsub -q 1nh $i ; done
