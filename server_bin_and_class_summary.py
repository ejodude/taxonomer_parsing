import os, sys, argparse
from Bio import SeqIO
import copy
import itertools
import multiprocessing as mp #allows for parallelization of the classification to speed up script.
import numpy as np
import timeit
import time
import re

#########
#globals
#########

parser=argparse.ArgumentParser(description="This script will generate binner and classifier counts for the default taxonomer databases")
parser.add_argument("class_file", type=str, help="classifier output file")
parser.add_argument("output", type=str, help="file path for count output files")
arg=parser.parse_args()


class_file = arg.class_file
file_name = arg.output


#troubleshooting files 
# class_file = "MS_2696.std_classifier.txt"
# file_name = "test_summary.txt"


#let's precompile the dictionary to store info:
count_dict = { 
"ambiguous": np.array( [0,0] ),
"bacterial": np.array( [0,0] ),
"fungal": np.array( [0,0] ),
"human": np.array( [0,0] ),
"phage": np.array( [0,0] ),
"phix": np.array( [0,0] ),
"unclassified": np.array( [0,0] ),
"viral": np.array( [0,0] ) }
#the array stores binner and classified counts respectively. 




#################
#functions
#################


		
def write_count_tables( count_dict, file_name ):
	"""
	function writes the contents of count data binner and classifier
	counts to an output file

	Needs:
		count_array- array of DBs x 2 size; with [:,0] as binner counts and [:,1] as classifier counts
		db_array - text info on database names 
		file_name - variable for writing the files should be root file name ie. "file" 
	"""
	
	#open our two output files. 
	out = open( file_name, 'w' )

	#write header lines. 
	header = ["database", "count"]
	out.write("%s\n" %( "\t".join( header ) ) )


	for key in count_dict:
		dat = count_dict[ key ] 
		out.write("%s\t%s\n" %( key + "_binner", str( dat[0] ) ) )
		out.write("%s\t%s\n" %( key + "_classifier", str( dat[1] ) ) )

	out.close()
		


#################
#MAIN
#################

#for turning into import and functional call:
if 1:
	#open file and iterate, while appending to dictioanry the count data. 
	for line in open( class_file, 'r'):
		line = line.split("\t")
		if line[1] == "C":
			count_dict[ line[0] ] = count_dict[ line[0] ] + 1  #append to bin and classify counters
		else:
			count_dict[ line[0] ] = count_dict[ line[0] ] + [1,0] #only append to binner counter obj[0]


	#write output to file. 
	write_count_tables( count_dict, file_name )












