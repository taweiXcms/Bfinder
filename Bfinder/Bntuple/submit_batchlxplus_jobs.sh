DIRECTORYOUTPUT="root://eoscms//eos/cms//store/group/phys_heavyions/wangj/BntupleData/ntB_EvtBase_20160407_BfinderData_PbPb_20160406_bPt5jpsiPt0tkPt0p8_BpB0BsX_sub1"
NAME="loop.C"

rm mylistfinal.txt
rm run_*.sh
rm -rf LSFJOB_*

g++ $NAME $(root-config --cflags --libs) -Wall -O2 -o "${NAME/%.C/}.exe"

eos ls /store/user/twang/HI_Bfinder/Data/HIOniaL1DoubleMu0/BfinderData_PbPb_20160406_bPt5jpsiPt0tkPt0p8_BpB0BsX_sub1/  | awk '{print "root://eoscms//eos/cms//store/user/twang/HI_Bfinder/Data/HIOniaL1DoubleMu0/BfinderData_PbPb_20160406_bPt5jpsiPt0tkPt0p8_BpB0BsX_sub1/" $0}' >> mylistfinal.txt

count=0 ; for i in `cat mylistfinal.txt` ; do echo cd $PWD/ > run_$count.sh ; echo 'export X509_USER_PROXY=~/x509_user_proxy'>> run_$count.sh ; 
echo eval \`scram runtime -sh\` >> run_$count.sh ; echo cd - >> run_$count.sh ; echo $PWD/loop.exe $i $DIRECTORYOUTPUT/ntuple_$count.root  >> run_$count.sh ; chmod +x run_$count.sh; count=$((count+1)) ; done
#for i in `ls run_*` ; do bsub -q cmscaf1nd $i ; done
for i in `ls run_*` ; do bsub -q 1nh $i ; done