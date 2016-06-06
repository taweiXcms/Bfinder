DIRECTORYOUTPUT="/store/group/phys_heavyions/HeavyFlavourRun2/BntupleData/DoubleMu/Bntuple_BfinderData_pp_20160419_bPt0jpsiPt0tkPt0p5/"
INPUTDIR="/store/group/phys_heavyions/HeavyFlavourRun2/BfinderRun2/DoubleMu/BfinderData_pp_20160419_bPt0jpsiPt0tkPt0p5/"

LOGFOLDER="/afs/cern.ch/user/t/twang/scratch0/stdout_merge/"
NAME="loop.C"
EOS="root://eoscms//eos/cms/"

rm mylistfinal.txt
rm run_*.sh
rm -rf LSFJOB_*
rm -rf $LOGFOLDER
mkdir $LOGFOLDER

g++ $NAME $(root-config --cflags --libs) -Wall -O2 -o "${NAME/%.C/}.exe"

eos ls $INPUTDIR | grep root | awk '{print "" $0}' >> mylistfinal.txt

maxf=1000 ;
count=0 ; 
for i in `cat mylistfinal.txt` ; do 
if [ $count -lt $maxf ]; then
ifexist=`eos ls $DIRECTORYOUTPUT/ntuple_$i`
if [ -z $ifexist ]; then
infn=`echo $i | awk -F "." '{print $1}'` ;
echo $infn ;
echo cd $PWD/ > run_$infn.sh ; 
echo 'export X509_USER_PROXY=~/x509_user_proxy'>> run_$infn.sh ; 
echo eval \`scram runtime -sh\` >> run_$infn.sh ; 
echo cd - >> run_$infn.sh ; 
echo $PWD/loop.exe $EOS$INPUTDIR$i $EOS$DIRECTORYOUTPUT/ntuple_$infn.root  >> run_$infn.sh ; 
chmod +x run_$infn.sh; 
fi
count=$((count+1)) ; 
fi
done

for i in `ls run_*` ; do 
#for i in `ls run_*` ; do bsub -q cmscaf1nd $i ; done
#for i in `ls run_*` ; do bsub -q 1nh $i ; done
bsub -o $LOGFOLDER -q 1nh $i; 
done

