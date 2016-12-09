import subprocess, threading
from subprocess import PIPE, STDOUT
import tempfile

class RunCommand():
    ''' Run with threading as python2 doesn't have timeout
        code from http://stackoverflow.com/questions/1191374/using-module-subprocess-with-timeout '''
    def __init__(self, cmd=None, inString=None):
        self.cmd = cmd
        self.inString = inString
        self.process = None

    def run(self, timeout = 1e10):
        def target():
#            tempfile.tempdir = None
#            tempFile = tempfile.TemporaryFile(bufsize = 0)
#            tempFile = tempfile.NamedTemporaryFile(bufsize = 0)
            self.process = subprocess.Popen(self.cmd, shell=True, stdout=PIPE, stdin=PIPE, stderr=STDOUT)#, bufsize = 4096)
#            self.process = subprocess.Popen(self.cmd, shell=True, stdout=PIPE, stdin=inStr, stderr=STDOUT)
#            self.process = subprocess.Popen(self.cmd, shell=True, stdout=tempFile, stdin=PIPE, stderr=tempFile)
#            self.process = subprocess.Popen(self.cmd, shell=True, stdout=PIPE, stdin=inStr, stderr=STDOUT)
            self.output = self.process.communicate(input=self.inString)[0]
#            print 'tempFile', tempFile.read()#;exit()
#            print self.output

#        thread = threading.Thread(target=target(self.inString))
        thread = threading.Thread(target=target)
        thread.start()

        thread.join(timeout)
        if thread.is_alive():
            self.process.terminate()
            thread.join()
        return self.process.returncode

class RunCommandNew():
    ''' Run with threading as python2 doesn't have timeout
        code from http://stackoverflow.com/questions/1191374/using-module-subprocess-with-timeout '''
    def __init__(self, cmd=None, inString=None):
        self.cmd = cmd
        self.inString = inString
        self.process = None

    def run(self, timeout = 1e10):
        def target():
            tempFile = tempfile.NamedTemporaryFile()
            self.process = subprocess.Popen(self.cmd, shell=True, stdout=tempFile, stdin=PIPE, stderr=tempFile)
#            self.process = subprocess.Popen(self.cmd, shell=True, stdout=PIPE, stdin=inStr, stderr=STDOUT)
            self.process.communicate(input=self.inString)#;exit()
#            self.output = self.process.communicate(input=self.inString)[0]
            with open(tempFile.name, 'r') as outf:
                self.output = outf.read()

#        thread = threading.Thread(target=target(self.inString))
        thread = threading.Thread(target=target)
        thread.start()

        thread.join(timeout)
        if thread.is_alive():
            self.process.terminate()
            thread.join()
        return self.process.returncode

#try:
#    from hardcode import GULP_EXE
#    runner = RunCommand(GULP_EXE, open('LTP.gin').read())
#    runner.run(timeout = 3)
#    print runner.output
#except:
#    print "somsomt"


class RunCommandSafe():
    ''' Crap but shouldn't hit problem with size
        sort out files and popen another time '''

    def __init__(self, cmd=None, inString=None, baseFileName='TEMP'):
        self.cmd = cmd
        self.inString = inString
        self.process = None
        self.baseFileName = baseFileName

    def run(self):
        import os 
        print "Warning - DO NOT RUN IN PARALLEL"

        with open(self.baseFileName + 'INPUT', 'w') as outf:
            outf.write(self.inString)
#        os.system(self.cmd + ' < ' + 'TEMPINPUT > TEMPOUTPUT') 
        os.system(self.cmd + ' < ' + self.baseFileName + 'INPUT > ' + self.baseFileName + 'OUTPUT') 
        with open(self.baseFileName + 'OUTPUT', 'r') as outf:
            self.output = outf.read()
        os.remove(self.baseFileName + 'INPUT')
        os.remove(self.baseFileName + 'OUTPUT')
