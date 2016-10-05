import os
import shutil




class SlurmScript():
    def __init__(self,
                 totalJobTime = None,
                 dlpolyJobs   = [],
                 jobName      = 'genericJob',
                 account      = 'Seed2016m',
                 partition    = 'compute-12',
                 ntasks       = 12,
                 scriptName   = None
                ):


    self.dlpolyJobs = dlpolyJobs

    if not totalJobTime:
        self.totalJobTime = sum([x.jobTime for x in dlpolyJobs])
    else:
        self.totalJobTime = totalJobTime

    self.jobName      = jobName
    self.account      = account
    self.partition    = partition
    self.ntasks       = ntasks
    if not scriptName:
        self.scriptName   = jobName + ".sh"
    else:
        self.scriptName   = scriptName

    #@classmethod later
    def makeFiles(self):
        for j in self.dlpolyJobs:
            os.makedirs(j.workingDirectory)
#            for x in contentsToCopy:
#                shutil.copy(x, newFolder)
        return

    def slurmJobString(self):
        import datetime
        outString  = "#!/bin/bash -l\n"
        outString += "#SBATCH --job-name=" + self.jobName + "\n"
        outString += "#SBATCH --time=" + str(datetime.timedelta(seconds=self.time)) + "\n"
        outString += "#SBATCH --account=" + self.account + "\n"
        outString += "#SBATCH --partition=" + self.partition + "\n"
        outString += "#SBATCH --ntasks=" + str(self.ntasks) + "\n"

        return outString

    def createJobScript(self):
        with open(self.scriptName, 'w') as outf:
            outf.write(self.slurmJobString())

    def submitJobScript(self):
        return



currentWorkingDirectory = os.getcwd()
materialName            = 'LTP'

def setupFolder(newFolder, *contentsToCopy):
    os.makedirs(newFolder)
    for x in contentsToCopy:
        shutil.copy(x, newFolder)
    return True

def writeSlurmScript():
    return True
