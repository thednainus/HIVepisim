import os, glob, subprocess

## This program simulates DNA sequence alignments using the program Seq-Gen

os.chdir('/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/R_projects/HIVepisim/Analyses/Results/VirusTreeSimulator_trees/FINAL')
my_files = glob.glob('*.tre')

output_dir = '/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/R_projects/HIVepisim/Analyses/Results/VirusTreeSimulator_trees/FINAL/Alignment/'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)


for filename in my_files:
    #seq-gen -mHKY -t0.5 -fe -l10000 -n1 < filename > filename+'.txt'
    # values for generating sequence alignment for HIV pol as described in van der Kyl and Berkhout, 2012
    # values of -s0.0028 based on paper by Patino-Galindo and Gonzalez-Candelas 2017
    # values of -t8.75 based on Chen et al 2004
    # Updated these values based on sequences for San Diego + background sequences
    cmd = ['/Applications/Seq-Gen-1.3.4/source/seq-gen', '-mHKY', '-l3012', '-t8.75', '-f0.389,0.165,0.228,0.218', '-s0.0028', '-n1', '-of', filename]
    run = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    stdout,stderr = run.communicate()
    #print stdout

    print stderr
    print filename



    if run.returncode != 0:
        raise RuntimeError("Seq-gen error: %s" % stderr)

    results = stdout


    OutFileName = output_dir+filename+'.fasta'


    OutFile = open(OutFileName,'a')

    OutFile.write(results)
    OutFile.close()