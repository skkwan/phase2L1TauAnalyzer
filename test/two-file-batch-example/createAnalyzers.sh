jobName="2019_Sep16-DYToLL_NoPU_withL1Tracks_PFCand-MINI"
#
j=0
for i in {0..1}
do
    cat test-Analyzer.py > SUB-Analyzer-${i}.py
    #echo "process.source.skipEvents = cms.untracked.uint32(${j})" >> SUB-Analyzer-${i}.py
    j=$(( $j + 100))
    cat submit.py >> SUB-Analyzer-${i}.py
    cat submit-$i.py >> SUB-Analyzer-${i}.py

    mkdir -p /nfs_scratch/skkwan/${jobName}/SUB-DYToLL-M-50-NoPU-14TeV-SUBPhase-$i/dags/daginputs

    farmoutAnalysisJobs --vsize-limit=7000 --assume-input-files-exist  --input-file-list=inputFileList-dummy.txt --output-dir=/hdfs/store/user/skkwan/${jobName} --submit-dir=/nfs_scratch/skkwan/${jobName}/SUB-DYToLL-M50-NoPU-SUBPhase-$i/submit --output-dag-file=/nfs_scratch/skkwan/${jobName}/SUB-DYToLL-M50-NoPU-SUBPhase-$i/dags/dag  ${jobName}-DYToLL-M50-NoPU  $CMSSW_BASE  $CMSSW_BASE/src/L1Trigger/phase2L1TauAnalyzer/test/two-file-batch-example/SUB-Analyzer-$i.py     &
    
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
