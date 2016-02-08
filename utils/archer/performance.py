import sys
import os
import argparse
import pprint
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

def main():
    #Check input args
    inputfile, outputfile = get_args()

    #Check if input file exists
    try:
        logfile = open(inputfile, 'r')
    except IOError:
        print ("Input file couldn't be opened. Exiting...")
        return

    #Elaborate data for output plot 
    job_id = get_job_id(inputfile)
    performance = get_performance(logfile)

    #Show output
    make_plot(job_id, performance)
    pprint.pprint(performance)    

def get_args():
    #Reads input arguments for input and output files
    parser = argparse.ArgumentParser(
           description='Extract performance data from given archer log file.')
    parser.add_argument('-i', dest="ifile", required=True, metavar="<inputfile>",
           help="Path to input archer log file. Tipically with extension '.o<JOBID>'.")
    parser.add_argument('-o', dest="ofile", required=True, metavar="<outputfile>",
           help="Path to output file containing the performance data.")
    args = parser.parse_args()
    return args.ifile, args.ofile
    
def get_job_id(filename):
    #Get Job ID from file's basename
    basename = os.path.basename(filename)
    return basename.split('.')[-1][1:8]

def get_performance(logfile):
    #Loop through file and save { thread_num: wall time } for each logged execution
    performance = {}
    for line in logfile:
        if "thread" in line:
            #Get 4th element in sentence
            thread_num = int(line.split(" ")[3])
	if "Wallclock" in line:
            #Get second element in value pair and remove last parenthesis
            wallclock_time = float(line.replace(")","").split(', ')[1])
            performance[thread_num] = wallclock_time
    return performance

def make_plot(title, a):
    #Plot performance in graph
    fig = plt.figure()
    plt.plot(a.keys(), a.values())     
    fig.savefig(title + ".png")


if __name__ == '__main__':
    main()


