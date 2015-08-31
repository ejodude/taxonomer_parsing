"""
script pareses taxonomer output for processing through DESeq
"""

import argparse #std python imports
import os
import sys
import re
import glob
import numpy as np

###############
#GLOBALs
###############



parser=argparse.ArgumentParser(description="script pareses taxonomer output into experiment count matrix")
parser.add_argument("dir", type=str, help="directory with count tables")
parser.add_argument("out", type=str, help="output file for final combined count table")
parser.add_argument("file_flag", type=str, help="file flag for globbing the correct count files")
arg=parser.parse_args()


"""
EXAMPLE
dir = "/data2/eosborne/strep_proj/strep_fastqs/taxonomer_analysis/strepDB_classifier_results"
out = "/data2/eosborne/strep_proj/strep_fastqs/taxonomer_analysis/strepDB_classifier_results_table-FULL.txt"
file_flag = "classifier-strepDB"
"""
dir = arg.dir 
out = arg.out 
file_flag = arg.file_flag 


###############
#Functions/
###############


class tax_result:
	def __init__(self, line):
		line = line.strip().split("\t")
		self.taxonomy = line[0]
		self.level = line[1]
		self.count = int( line[2] )


#this loop processes a file into a data structure. 
def process_count_file( file ):
	file_dict = {} #store results for single file. 

	f = open( file, 'r' )
	for line in f:
		if line.split("\t")[0] != "taxa": #skip header
			result = tax_result( line ) 
			file_dict[ result.taxonomy ] = [ result.count, result.level ]

	f.close()
	return( file_dict )
	

#this loop processes the file results and incorporates them int the data struct 
#that summarizes the overall results. 
def combine_results( count_dict, sample_list, file_results ):

	#iterate over taxIDs in count dictionary (one's we've seen) already
	for key in count_dict: 

		if key in file_results: #do we have results for gene ?
			cdat = file_results[ key ][0]
			count_dict[ key ] = count_dict[ key ] + [ cdat ]

		else: #we dont have results for this taxon in the current file. 
			cdat = [0] #we have zero counts for taxID
			count_dict[ key ] = count_dict[ key ] + cdat 

	#now go through file_results. can use loop to instatiate new taxIDs.
	for k in file_results:
		if k not in count_dict: #if we haven't seen tax id before..
			cdat = file_results[ k ][0]
			empty_dat = [0] * ( len( sample_list ) -1 )
			count_dict[ k ] = empty_dat + [ cdat ]
		else:
			pass

	return( count_dict )

#function writes the counts to a 2-D array

def write_results( count_dict, file_results, out ):
	o = open(out, 'w') #open output

	tkeys = count_dict.keys() #grab dictionary keys
	tkeys.sort() #sort them on tax ID; should 
	data = [] #empty data list for 2-D count array
	levels = [] #keep track of level data. 
	for tk in tkeys:
		data.append( count_dict[ tk ] )
		levels.append( file_results[ tk ][1] )
	
	#turn 2D list data into an array. 
	dat_array = np.array(data)

	#now write this all to file. 
	taxIDs = np.array( [ tkeys ] ).T #convert taxIDs for rownames
	levels = np.array( [levels]  ).T #convert levels to another row dataset.
	col_meta_data = np.hstack((taxIDs, levels)) #join column data

	head = ["taxa", "level"] + sample_list #header line with samples. 
	head = "\t".join( head )
	np.savetxt(o, np.hstack((col_meta_data, dat_array)), delimiter='\t', fmt='%s', header= head )


###############
#MAIN
###############


count_dict = {} #initialize empty dict structure. 
sample_list =[] #sample IDs

files = glob.glob( dir + "/*%s*.txt" %( file_flag ) ) 

#iterate over the files. 
for file in files:
	sample_id = file.split("/")[-1].split(".")[0] #sample id
	sample_list.append( sample_id ) #append sample id to list.
	#can use len(sample_list) to initialize new taxIDs to the count dict

	file_results = process_count_file( file )
	count_dict = combine_results( count_dict, sample_list, file_results ) #update count dict

#write the data out. 
write_results( count_dict, file_results, out )

















