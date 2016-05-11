DIRECTORYOUTPUT="/store/group/phys_heavyions/HeavyFlavourRun2/DntupleData/ntD_EvtBase_20160511_HIHardProbes_DfinderData_PbPb_20160402_dPt0tkPt2p5_D0Dstar3p5p_skimhltbranches/"
INPUTDIR="/store/group/phys_heavyions/wangj/DntupleData/ntD_EvtBase_20160405_HIHardProbes_DfinderData_PbPb_20160402_dPt0tkPt2p5_D0Dstar3p5p_FINALJSON/"

LOGFOLDER="/afs/cern.ch/user/t/twang/scratch0/stdout_merge/"
NAME="skimTree_Dzero.C"
EOS="root://eoscms//eos/cms/"
isPbPb=1

rm mylistfinal.txt
rm run_*.sh
rm -rf LSFJOB_*
rm -rf $LOGFOLDER
mkdir $LOGFOLDER

g++ $NAME $(root-config --cflags --libs) -Wall -O2 -o "${NAME/%.C/}.exe"

eos ls $INPUTDIR  | awk '{print "" $0}' >> mylistfinal.txt

maxf=5000 ;
count=0 ; 
for i in `cat mylistfinal.txt` ; do 
if [ $count -lt $maxf ]; then
ifexist=`eos ls $DIRECTORYOUTPUT/$i`
if [ -z $ifexist ]; then
infn=`echo $i | awk -F "." '{print $1}'` ;
echo $infn ;
echo cd $PWD/ > run_$infn.sh ; 
echo 'export X509_USER_PROXY=~/x509_user_proxy'>> run_$infn.sh ; 
echo eval \`scram runtime -sh\` >> run_$infn.sh ; 
#echo cmsenv >> run_$infn.sh
echo cd - >> run_$infn.sh ; 
echo $PWD/skimTree_Dzero.exe $EOS$INPUTDIR$i $EOS$DIRECTORYOUTPUT/$infn.root $isPbPb  >> run_$infn.sh ; 
chmod +x run_$infn.sh; 
fi
count=$((count+1)) ; 
fi
done

for i in `ls run_*` ; do 
#echo test
#for i in `ls run_*` ; do bsub -q cmscaf1nd $i ; done
#for i in `ls run_*` ; do bsub -q 1nh $i ; done
bsub -o $LOGFOLDER -q 1nh $i; 
done

