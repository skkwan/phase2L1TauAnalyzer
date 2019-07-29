To run the two file solution using farmOutAnalysisJobs:
1. First place your analyzer in this folder and call it test-Analyzer.py
2. Next put the Secondary File list in inputFileList.txt and the mini-AOD list in inputFileList-MINI.txt
3. Create the submit footers, submit-x.py using the command ```python dasqueryMiniFEVT.py inputFileList.txt inputFileList-MINI.txt```

**For the next step PLEASE run a test job to see that it runs and outputs a file on /hdfs/ before submitting all jobs.** You can do this by changing the bash script to run on only two files like: ```for i in {0..1}```

4. Create and submit the analyzers ```bash createAnalyzers.sh```