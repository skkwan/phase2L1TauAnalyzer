jobName="2019_Aug29-ggHtautau_noPU-try2"
#
j=0
for i in {0..11}
do
    cat test-Analyzer.py > SUB-Analyzer-${i}.py
    #echo "process.source.skipEvents = cms.untracked.uint32(${j})" >> SUB-Analyzer-${i}.py
    j=$(( $j + 100))
    cat submit.py >> SUB-Analyzer-${i}.py
    cat submit-$i.py >> SUB-Analyzer-${i}.py


    mkdir -p /nfs_scratch/skkwan/${jobName}/SUB-ggHtautau-SUBPhase-$i/dags/daginputs

    #removing  --vsize-limit=7000
    farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=inputFileListAODSIM.txt --output-dir=/hdfs/store/user/skkwan/${jobName} --submit-dir=/nfs_scratch/skkwan/${jobName}/SUB-ggHtautau-SUBPhase-$i/submit --output-dag-file=/nfs_scratch/skkwan/${jobName}/SUB-ggHtautau-SUBPhase-$i/dags/dag  ${jobName}-ggHtautau-noPU  $CMSSW_BASE  $CMSSW_BASE/src/L1Trigger/phase2L1TauAnalyzer/test/aodsim-two-file-batch/SUB-Analyzer-$i.py   &
    
    if [ "$i" -eq "15" ]; then
	wait;
    fi

    if [ "$i" -eq "30" ]; then
	wait;
    fi

    if [ "$i" -eq "45" ]; then
	wait;
    fi

    if [ "$i" -eq "60" ]; then
	wait;
    fi


done
