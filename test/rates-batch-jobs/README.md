# To run the one-file solution (FEVT) using farmOutAnalysisJobs:
1. First place your analyzer in this folder and call it test-Analyzer-rates.py
2. Next put the FEVT file lists in the `samples` folder, with .txt file names specified in ```userInputSubmitPhaseIITau.sh```.
3. The command is ```bash userInputSubmitPhaseIITau.sh [jobName]```. For example,
   ```
   nohup bash userInputSubmitPhaseIITau.sh 2019_Sep25 > logSubmit_Sep25 &
   tail -f logSubmit_Sept25
   ```
