#import das_client
from sys import argv
import os

inputFile = argv[1] #contains a list of files that you're running over, FEVT version
inputFileMini = argv[2] #contains the list of files in mini format
j = 0
with open(inputFile,'r') as f:
    for x in f:
        file = open('submit-'+str(j)+'.py','w')
        j = j+1
        x = x.rstrip()
        output='\nprocess.source.secondaryFileNames = cms.untracked.vstring( \n "'+x+'") \nprocess.source.lumisToProcess = cms.untracked.VLuminosityBlockRange("1:'
        x = x.split('/store',1)[-1] 
        x = '/store'+x #gets rid of any xrootd info in the filename
        das_query_string = 'dasgoclient --query="lumi file=%s" --limit=0' % x
        das_output = os.popen(das_query_string).read() #query DAS via bash
        das_output = das_output.replace("root","root',")#commas between each secondary file
        das_output = das_output.replace("/store","'root://cmsxrootd.fnal.gov//store") #need xrootd info in secondary filename
        das_output = das_output.replace(",","\",\"1:") #don't want newline in string
        das_output = das_output.replace("\n","") #don't want newline in string
        das_output = das_output.replace("[","") #don't want in string
        das_output = das_output.replace("]","") #don't want in string

        output2='\nprocess.source.fileNames = cms.untracked.vstring( \n '
        das_query_string = 'dasgoclient --query="child file=%s" --limit=0' % x
        das_output_2 = os.popen(das_query_string).read() #query DAS via bash
        if not das_output_2:
            print("das_output_2 is empty, inputting the full miniaod list instead")
            with open(inputFileMini,"r") as f2:
                das_output_2 = f2.read()
        das_output_2 = das_output_2.replace("/store","'root://cmsxrootd.fnal.gov//store") #need xrootd info in secondary filename
        das_output_2 = das_output_2.replace("\n","") #don't want newline in string
        das_output_2 = das_output_2.replace("[","") #don't want in string
        das_output_2 = das_output_2.replace("]","") #don't want in string
        das_output_2 = das_output_2.replace("root'root","root','root")#commas between each secondary file

        #print das_output
        output = output+das_output+'")\n\n'+output2+das_output_2+"')\n\n"
        file.write(output)
        file.close()
        print output
