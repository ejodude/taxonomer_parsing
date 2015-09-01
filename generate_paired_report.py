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



parser=argparse.ArgumentParser(description="script pareses taxonomer output into experiment count matrix\
 for a paired exprimental-control sample")
parser.add_argument("expt_file", type=str, help="sample count tables; should be full count result")
parser.add_argument("ctl_file", type=str, help="control count table; should be full count result")
parser.add_argument("out", type=str, help="output file for final combined count table")
parser.add_argument("--method", type=str, help="normalization method: one of [ 'level_proportion' ] ")
parser.add_argument("--rd_cutoff", type=int, help="number of reads for taxon to be reported in the output file")
parser.set_defaults(rd_cutoff=10) #set default count cutoff to 10 if arg not supplied.
parser.set_defaults(method="level_proportion") #set default normalization method 
arg=parser.parse_args()


"""
EXAMPLE
expt_file = "/data2/eosborne/taxonomer_projs/schlaberg_CSF/taxonomer_results/bacterial/CSF-1_cat.fa.taxonomer_v2_full_cnt.txt"
ctl_file= "/data2/eosborne/taxonomer_projs/schlaberg_CSF/taxonomer_results/bacterial/CSF-control_cat.fa.taxonomer_v2.out.screened.v2_full_cnt.txt"
out = "test_matrix.txt"
method = "level_proportion"
cutoff = 10

"""
expt_file = arg.expt_file
ctl_file= arg.ctl_file
out = arg.out 
method = arg.method 
cutoff = arg.rd_cutoff

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
	

def normalize_counts( array, levels=None, method="level_proportion" ):
	if method == "level_proportion":
		norm_array = level_prop_norm( array, levels )
	else:
		raise Exception( "not a supported normalization method. try running with -h to see supported opts" )
	return( norm_array )



def level_prop_norm( array, levels ):
	"""
	function converts counts to proportions based on the level of the taxonomic assignment
	and the total read counts assigned to the taxonomic level. 

	proportion value = counts[ level N ]  / total reads classified[ level N] 
	"""
	norm_array = np.zeros( array.shape, dtype='f' ) #generate empty array

	for level in np.unique( levels ):
		#print level
		bools = levels == level
		tot_rds = np.sum( array[ bools ], axis= 0, dtype='f' ) #need a float array to get proportions
		norm_array[bools] = array[ bools ] / tot_rds

	return( norm_array ) 


def apply_count_filter( expt_array, cutoff = 10 ):
	"""
	function looks at experimental array for any taxa that have
	counts above certain cutoff in at least one of the samples. 
	
	Can be modified for more complex filtering... 
	
	returns boolean array of samples passing filtering 
	"""
	#grab any genes where at least one expt are greater-equal to cutoff. 
	filter_bools = np.any( expt_array >= cutoff, axis=1)
	return( filter_bools )
	


#function takes the full count matrices, filters and processes
#them for normalization.
def process_results( expt_dict, ctl_dict, cutoff, method ):
	"""
	#function takes the full count matrices, filters and processes
	#them for normalization.

	#does more logic on the experimental samples for filtering, which 
	#is why expt and control and seperated at first. 

	"""
	RV = {} #stores return dataset

	tkeys = expt_dict.keys() #grab dictionary keys
	expt_data = [] #empty data list for 2-D count array for experimentals
	ctl_data = [] #empty data list for 2-D count array for control samples
	levels = [] #keep track of level data. 
	for tk in tkeys:
		expt_data.append( expt_dict[ tk ] )
		ctl_data.append( ctl_dict[ tk ] )
		levels.append( file_results[ tk ][1] )

	#turn 2D list data into an array. 
	expt_array = np.array(expt_data)
	ctl_array = np.array(ctl_data)

	#now write this all to file. 
	taxIDs = np.array( tkeys ).T #convert taxIDs for rownames
	levels = np.array( levels  ).T #convert levels to another row dataset.

	#perform normalizations based on method supplied..
	norm_expt = normalize_counts( expt_array, levels, method=method )
	norm_ctl = normalize_counts( ctl_array, levels, method=method )

	#now filter based on read counts...
	filter_bools = apply_count_filter( expt_array, cutoff = cutoff)

	#generate return data structure. 
	RV[ 'expt' ] = norm_expt[ filter_bools ]
	RV[ 'ctl' ] = norm_ctl[ filter_bools ]
	RV[ 'taxIDs' ] = taxIDs[ filter_bools ] 
	RV[ 'levels' ] = levels[ filter_bools ] 

	return( RV )


#function writes the counts to a 2-D array
def write_results( RV, expt_list, ctl_list, out ):
	"""
	function takes in the normalized counts for expt and control and 
	combines them into a mega matrix, with the corresponding sample IDS

	Sorts the output based on taxonomic level and taxonomy name for
	easy interpretation
	"""


	o = open(out, 'w') #open output
	taxIDs = RV[ 'taxIDs'] 
	levels = RV[ 'levels'] 
	expt = RV[ 'expt'] 
	ctl = RV[ 'ctl'] 

	#combine expt and ctl matrix by row. 
	count_mat = np.concatenate( (expt, ctl ), axis=1 ) 
	sample_list = expt_list + ctl_list #combine sample names. 

	if len( sample_list ) != count_mat.shape[1]:
		raise Exception("mismatch in dimensions between sample and count matrix")


	#combine taxids and levels for sorting. sorts by level and then tax id.
	#returns indexes of sort. 
	sort_order = np.lexsort( ( taxIDs, levels ) )



	#now sort all of the data and write to output file... 
	sort_taxIDs = np.array( taxIDs[ sort_order ]  ) #convert taxIDs for rownames
	sort_levels = np.array( levels[ sort_order ]  ) #convert levels to another row dataset.
	sort_count = count_mat[ sort_order, : ] 
	col_meta_data = np.vstack((sort_taxIDs, sort_levels)).T #join column data

	#use numpy to write everything to output file. 
	head = ["taxa", "level"] + sample_list #header line with samples. 
	head = "\t".join( head )
	np.savetxt(o, np.hstack((col_meta_data, sort_count)), delimiter='\t', fmt='%s', header= head )






###############
#MAIN
###############

if 1:
	#process the experimental data
	
	expt_dict = {} #initialize empty dict structure. 
	expt_list =[] #sample IDs


	#process the expt file(s)
	files = expt_file.split(",")


	#iterate over the files. 
	for file in files:
		sample_id = file.split("/")[-1].split(".")[0] #sample id
		expt_list.append( sample_id ) #append sample id to list.
		#can use len(sample_list) to initialize new taxIDs to the count dict

		file_results = process_count_file( file )
		expt_dict = combine_results( expt_dict, expt_list, file_results ) #update count dict


	#process the control samples seperately
	ctl_dict = {} #initialize empty dict structure. 
	ctl_list =[] #sample IDs

	#process the expt file(s)
	files = ctl_file.split(",")

	#iterate over the files. 
	for file in files:
		sample_id = file.split("/")[-1].split(".")[0] #sample id
		ctl_list.append( sample_id ) #append sample id to list.
		#can use len(sample_list) to initialize new taxIDs to the count dict

		file_results = process_count_file( file )
		ctl_dict = combine_results( ctl_dict, ctl_list, file_results ) #update count dict


	#normalize the datasets and apply count filters to reduce search space. 
	count_dataset = process_results( expt_dict, ctl_dict, cutoff, method )

	#finally combine and write the data to output file. 
	write_results( count_dataset, expt_list, ctl_list, out )










