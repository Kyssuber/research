#!/usr/bin/env python

'''
GOAL:
* This will create a bash script to run multiple serial jobs in parallel using slurms array option

USAGE:
* when logging into grawp
  module load Python3

* Run this in the output directory that has a subfolder for each galaxy.
  /mnt/astrophysics/rfinn/muchogalfit-output

* Try running this once and check the output of the JOB_{}.sh script.
* If it looks good, run again using the --submit flag to submit the script to slurm.

'''

import os
import subprocess
import sys
import glob

HOME = os.getenv("HOME")


def write_output(script_id, input_file, narray=1000, data_dir=None, wavelength=None, submit=False):
    ''' copying python from Matt Bellis, commands from Ryan Decker '''
    output = ""
    output += "#!/bin/bash\n"
    output += "\n"
    output += "# Set Job Name\n"
    output += "#SBATCH -J job\n"
    output += "\n"
    output += "# Set file to capture standard out and standard error and append the jobID (%j)\n"
    output += "#SBATCH -o job.out.%j\n"
    output += "\n"
    output += "#SBATCH --partition=normal\n"
    output += "\n"
    # running in array mode, rather than spawning narray independent processes
    output += "# for testing\n"
    output += f"#SBATCH --array=1-{narray}\n"
    output += "\n"
    output += "#Set the number of nodes\n"
    output += "#SBATCH -N 1\n"
    output += "#SBATCH --ntasks=1\n"
    output += "\n"
    output += "#Set the time limit for the job - 10.5 hour is specified here\n"
    output += "#SBATCH --time=10:30:00\n"
    output += "\n"
    output += "#SBATCH --cpus-per-task=1\n"
    output += "\n"
    output += "# Load any environmental modules needed\n"
    output += "module load Python3\n"
    # print version of typing-extensions
    #output += "pip3 list |grep typing"
    #output += "echo"
    #output += "printenv"
    #output += "yes Y | pip3 uninstall typing-extensions\n"
    output += "pip3 install typing-extensions==4.0.1\n"
    output += "pip3 list |grep typing \n"    
    output += "module load gnu9\n"
    output += "\n"
    output += "# Move to the directory needed - defaults to the submission directory\n"
    output += "\n"
    output += "# perform calculation\n"
    output += "#\n"
    if data_dir is not None:
        s = f'LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p {data_dir}/{input_file})\n'
    else:
        print('please provide a valid data directory')
        return
    if wavelength is None:
        print('please provide a valid wavelength as input to slurm constructor function ')
        return
        
    output += s
    output += "#\n"
    output += f"python3 {HOME}/github/mucho-galfit/parallel/run1galfit.py $LINE {wavelength}\n"

    outfname = f"JOB_{script_id}.sh"

    outfile = open(outfname, 'w')
    outfile.write(output)
    outfile.close()

    cmds = ['sbatch', outfname]
    if submit:
        print(f"Submitting job to process {input_file}")
        process = subprocess.Popen(cmds, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        print(stdout.decode())
        print(stderr.decode())
        pass


###########################
###########################
import argparse

parser = argparse.ArgumentParser(
    description='Program to create bash script to run galfit in parallel.  This uses the array option in slurm so that the processes are associated with the same process id, rather than submitting a bunch of individual jobs.  The script can be submitted by setting the --submit flag.')
parser.add_argument('--wavelength',
                    dest='wavelength',
                    default='W3',
                    #options = ['W3'],
                    help='Wavelength of images to analyze')

parser.add_argument('--testsample',
                    dest='testsample',
                    default=False,
                    action='store_true',
                    help='Set this option to submit a script with a test sample of 10 galaxies to slurm. Default is false. You can check output JOB_{scriptid}.sh.  If everything looks good, then add submit flag.')
parser.add_argument('--submit',
                    dest='submit',
                    default=False,
                    action='store_true',
                    help='Set this option to submit the script to slurm. Default is false. You can check output JOB_{scriptid}.sh.  If everything looks good, then add submit flag.')
args = parser.parse_args()

###########################################################
cwd = os.getcwd()

# what is this directory???
data_dir = f"{HOME}/research/wisesize/"
# this is the directory that has a subdirectory for each galaxy
# the following is the directory that grawp sees (remove rfinn if running on virgo)
data_dir = "/mnt/astrophysics/rfinn/muchogalfit-output/"


# this is the name of the shell script that will be created
script_id = f"VFIDall-{args.wavelength}"

print('data_dir = ', data_dir)
print()
print('script_id = ', script_id) 
print()
os.chdir(data_dir)

# this assumes that the data directory has a file called Dirs.txt
#
# Dirs.txt contains one line for each galaxy that will be analyzed.
#

outfile = "Dirs.txt"
os.system(f"ls -d VFID???? > {outfile}")
# count lines - need to give number of galaxies to slurm b/c this is equal to number of processes
infile = open(outfile, 'r')
nfiles = (len(infile.readlines()))
infile.close()

print(f"\nthe number of lines in Dirs.txt = {nfiles}\n")
# write out files and submit jobs
# for d in dirlist:

input_file = "Dirs.txt"

# set sample size to 10 galaxies for testing

if args.testsample:
    nfiles=20
    print(f"\nSetting number of jobs to {nfiles} for testing.\n")
    
write_output(script_id, input_file, narray=nfiles, data_dir=data_dir, \
                 wavelength=args.wavelength,submit=args.submit)

os.chdir(cwd)
