import os
import subprocess
os.chdir('C:\Documents and Settings\Administrator\My Documents\dcProgs\dcprogs')
process = subprocess.Popen('hjcfit.exe', shell=True, stdout=subprocess.PIPE)
process.wait()
print process.returncode