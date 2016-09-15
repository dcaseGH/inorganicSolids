import subprocess, threading
from subprocess import PIPE, STDOUT

class RunCommand():
    ''' Run with threading as python2 doesn't have timeout
        code from http://stackoverflow.com/questions/1191374/using-module-subprocess-with-timeout '''
    def __init__(self, cmd=None, inString=None):
        self.cmd = cmd
        self.inString = inString
        self.process = None

    def run(self, timeout = 1e10):
        def target():
            self.process = subprocess.Popen(self.cmd, shell=True, stdout=PIPE, stdin=PIPE, stderr=STDOUT)
            self.output = self.process.communicate(input=self.inString)[0]

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
