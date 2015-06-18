import os, sys, argparse
from Bio import SeqIO
import copy
import itertools
import multiprocessing as mp #allows for parallelization of the classification to speed up script.
import numpy as np
import timeit
import re

#########
#globals
#########

parser=argparse.ArgumentParser(description="run taxonomy count parsing")
parser.add_argument("key_file", type=str, help="taxonomer key file")
parser.add_argument("class_file", type=str, help="classifier output file")
parser.add_argument("output", type=str, help="file path for count output files")
parser.add_argument("--db", type=str, help="database to parse for server output.")
parser.add_argument("--server_out", dest='server_out', action='store_true', help ="boolean for whether results are server database (if custom db, do not use arg)")
parser.set_defaults(server_out=False)
parser.add_argument("--proc", type=str, help="optional arg for number of proceses to run with classifier; higher thread number will increase speed of classifier; defaults to 1", default = 1)
arg=parser.parse_args()


# key_file = "/data2/eosborne/taxonomer_dbs/strep_MLST/strep_pneumo-gene/strep_pneumo-gene_MLST_key.txt"
# server_output = False
# class_file = "/data2/eosborne/taxonomer_dbs/strep_MLST/mark_subtraction_db/aroE_classifier/test_binner-aroE-binner-aroE_concat_classify-ties.txt"

# key_file = "/data2/eosborne/taxonomer_dbs/strep_db/strepDB_key.txt"
# server_output = False
# class_file = "/data2/eosborne/strep_proj/strep_fastqs/10847R_reads/10847X31_S30722.pe-combo.strepDB_classifier.txt"
# 
# 
# key_file = "/data2/eosborne/taxonomer_dbs/db_KAn/hg38.out_key.txt"
# server_output = False
# #class_file = "/data2/eosborne/taxonomer_projs/MS_RNA_seq/fastqs/11168X9_S043-HC_150129_D00294_0156_AC6E82ANXX_6-pe-combo.hg38.k31_classifier.txt"
# class_file = "./11168X9_S043-HC_150129_D00294_0156_AC6E82ANXX_6-pe-combo.hg38.k31_classifier-HEAD.txt"
# file_name = "./test"
# 
# 
# key_file = "/data2/eosborne/taxonomer_dbs/taxonomer_taxid_keys/bacterial_16S_key.txt"
# server_output = True
# #class_file = "/data2/eosborne/taxonomer_projs/MS_RNA_seq/fastqs/11168X9_S043-HC_150129_D00294_0156_AC6E82ANXX_6-pe_combo.std_classifier.txt"
# class_file = "11168X9_S043-HC_150129_D00294_0156_AC6E82ANXX_6-pe_combo.std_classifier-HEAD.txt"
# database = "bacterial"
# file_name = "./test_server"

key_file = arg.key_file
class_file = arg.class_file
file_name = arg.output
server_output = arg.server_out
database = arg.db






#################
#functions
#################


#load in the taxID number and data relationship
def load_keys( key_file ):
	"""
	key_file - file with key id number and taxonomy information
	returns dictionary with 
	d[taxid] = [ [ parent, level, taxonomy ] ,  counter ]
	"""
	kfile = open( key_file, 'r' ) #open file

	key_dict = {} #create dict structure for look up
	for line in kfile: #iterate over lines. 
		line = line.strip().split("\t")
		k = line[0] #taxID
		data = [ line[1:], 0 ] # [ [taxonomic info], counter ]
		key_dict[ k ] = data

	kfile.close()
	return( key_dict ) 
	
	
#class object for processing VCF files. 
class classify_results:
	"""
	makes an iterator to process results to custom classifier database
	chunksize argument controls the numbers of classified sequences to
	chunk at a given time for post-processing in other functions. 
	"""
	def __init__(self,file, chunksize=10000):
		self.f = open(file,'r') #open file handle upon init
		self.line = self.f.readline() #read first line of the file. 
		self.chunksize = chunksize
	def __iter__(self):
		return self
	
	def next(self):
		cnt = 0 #boolean for chunking.
		return_array = [] #initialize empty array to store line data. 
		#check here if we are currently on last line, and raise StopIteration to exit next()
		if len(self.line) == 0: 
			raise StopIteration
		while cnt < self.chunksize:
			line = self.line.strip()
			if len( line ) == 0:
				return( return_array ) 
				break #break out of loop because we are at last line in file. 
			elif line[0] == "U":
				pass
			else: #we have a classified read. 
				return_array.append( line ) 
				cnt += 1
			self.line = self.f.readline()
		return( return_array )


#class object for processing VCF files. 
class classify_results_server:
	"""
	makes an iterator to process results to custom classifier database
	chunksize argument controls the numbers of classified sequences to
	chunk at a given time for post-processing in other functions. 
	"""
	def __init__(self,file, database, chunksize=10000 ):
		self.f = open(file,'r') #open file handle upon init
		self.line = self.f.readline() #read first line of the file. 
		self.chunksize = chunksize
		self.db = database
		self.charlen = len(database)
	def __iter__(self):
		return self
	
	def next(self):
		cnt = 0 #boolean for chunking.
		return_array = [] #initialize empty array to store line data. 
		#check here if we are currently on last line, and raise StopIteration to exit next()
		if len(self.line) == 0: 
			raise StopIteration
		while cnt < self.chunksize:
			line = self.line.strip()
			if len( line ) == 0:
				return( return_array ) 
				break #break out of loop because we are at last line in file. 
			if database + "\tC" in line: #gather classified reads to the particular DB
				return_array.append( line ) 
				cnt += 1
			else: #we have a classified read. 
				pass
			self.line = self.f.readline()
		return( return_array )





class individual_read:
	"""
	class individual_read standardized object for a line of classifier output
	stores the relevant data for a classifier line 
	"""
	def __init__(self,data):
		self.taxid = data[2] #taxid- numerical id for seq
		self.name =  data[3] #seqid - text id for sequence. 
		self.tax_level = data[5] #taxonomic level of classification
		self.read_id = data[1] #sequence read ID. 


class individual_server_read:
	"""
	class individual_read standardized object for a line of classifier output
	stores the relevant data for a classifier line 
	"""
	def __init__(self,data):
		self.taxid = data[2] #taxid- numerical id for seq
		self.name =  data[3] #seqid - text id for sequence. 
		self.tax_level = data[5] #taxonomic level of classification
		self.read_id = data[1] #sequence read ID. 
		self.db = data[0]



def process_classify_file( data ): 
	"""
	#iterates over file and counts the number of times a read is 
	#classfied to each unique ID number. This loop parses output from single
	#classfier database run. 
	
	data - taxonomer object data structure with returned classified lines in list format
	id_dict - parsed key file results with counter for read count tracking
	"""
	#print data
	#pull results from itertools
	results = data[0]
	id_dict = data[1]
	
	#process the object and append counters
	RV = copy.deepcopy(id_dict) #dictionary of returned results
	for result in results: #iterate over lines. 
		result = result.split("\t")
		taxid = result[3]
		RV[ taxid ][1] += 1

	return( RV ) #return modified dictioanry with updated counts





def process_individual_classify_file( data ): 
	"""
	#iterates over file and counts the number of times a read is 
	#classfied to each unique ID number. This loop parses output from single
	#classfier database run. 
	
	data - taxonomer object data structure with returned classified lines in list format
	id_dict - parsed key file results with counter for read count tracking
	"""
	#print data
	#pull results from itertools
	results = data[0]
	id_dict = data[1]
	
	#process the object and append counters
	RV = copy.deepcopy(id_dict) #dictionary of returned results
	for result in results: #iterate over lines. 
		result = result.split("\t")
		taxid = result[2]
		RV[ taxid ][1] += 1

	return( RV ) #return modified dictioanry with updated counts




def proc_counts( count_data, r):
	#might want to remove this bit of logic?
	"""
	function takes in the count data and iterator results object r
	and appends the reads classified in r to the count data. 
	essentially creates a running sum for the count data as we process
	reads without locking things up with a global counter variable
	
	"""
	if count_data[0].shape[0] != len( r.keys() ):
		raise Exception('key lengths differ btw count dict and the processed reads')
		
	counts = np.array( [cnt[1] for cnt in r.values()] ) #totals for given chunk we processed. 

	#append the counts in r to the count dataset
	count_data[1] = count_data[1] + counts

	return( count_data ) #return the tallies. 



def parse_taxonomy_counts( count_data, id_dict):
	"""
	function takes in count data and id_dict
	and uses parent-child relationships to perpetuate
	counts through the taxonomy. 
	count_data[0] -tax ids
	count_data[1] - classifier counts

	returns count data with new entry:
	count_data[2] which is the perpetuated count matrix.
	"""

	taxa = count_data[0]
	counts = count_data[1]
	#make deep copy to store our perpetuated counts
	appened_counts = copy.deepcopy( count_data[1] )
 
	#iterate over values in count_data
	for i in xrange( taxa.shape[0]): 
		tax = taxa[i]
		count = counts[i]
		tax_info = id_dict[tax][0] #[parent, tax level, name] 
	
		if count == 0: #we have no reads to perpetuate through taxonomy. 
			continue #so skip to the next iterator
		
		#perpetuate through the taxonomy until you hit root.
		current_tax = tax 
		while current_tax != '1':
			#print  current_tax
			#print count
			parent = tax_info[0] 
			#when you find parent, append the root child count
			appened_counts[ taxa == parent ] += count
	
			#update the current tax and info for next iter. 
			current_tax = parent
			tax_info = id_dict[ current_tax ][0]

	#add perpetuated counts to count data
	count_data.append( appened_counts )
	return( count_data ) 
		
		
		
def write_count_tables( count_data, id_dict, file_name ):
	"""
	function writes the contents of count data to two tables
		1). raw classifier counts "raw_cnt.txt"
		2). perpetuated count table "full_cnt.txt"
	for raw and perpetuated counts respectively

	Needs:
		count_data - counts and associated info
		id_dict - contain taxonomy information 
		file_name - variable for writing the files should be root file name ie. "file" 
	"""
	taxa = count_data[0]
	counts = count_data[1]
	appended_counts = count_data[2]

	#open our two output files. 
	raw = open( file_name + "_raw_cnt.txt", 'w' )
	appended = open( file_name + "_full_cnt.txt", 'w' )

	#write header lines. 
	header = ["taxa", "taxon_level", "count"]
	raw.write("%s\n" %( "\t".join( header ) ) )
	appended.write("%s\n" %( "\t".join( header ) ) )

	#iterate over numerically sorted taxon list. 
	for tax in sorted(taxa, key=lambda item: int(item)):
		tax_dat = id_dict[tax][0]
		level =  tax_dat[1]  #extra shit to turn this into easy int for sorting/filteirng.
		name = tax_dat[2]
		print_list = [ name, level ] 
	
		#now send the results to outfile. 
		tax_bool = taxa == tax #bool array for where tax id resided in data fields. 
		cnt = str( counts[ tax_bool ][0] )
		raw.write("%s\n" %(  "\t".join( print_list + [cnt] ) ) )
		cnt = str( appended_counts[ tax_bool ][0] )
		appended.write("%s\n" %(  "\t".join( print_list + [cnt] ) ) )



#################
#process ties
#################

	
#class object for processing VCF files. 
class tie_result_chunk:
	"""
	class tie_result generates an iterator for looping through ties in a taxonomer output 
	file. returns a list of ties for a given sequence read. 
	"""
	def __init__(self,file, chunksize=10000):
		self.f = open(file,'r')
		self.line = self.f.readline() #read first line of the file. 
		self.chunk = chunksize
	def __iter__(self):
		return self
	
	def next(self):
		read_array = []
		
		#check here if we are currently on last line, and raise StopIteration to exit next()
		if len(self.line) == 0: #
			raise StopIteration
		counter = 0	
		while len(read_array) < self.chunk:
			#print counter ##r
			counter += 1 ##r
			read = "" #initiate read on each iteration
			read_boolean = True
			line = self.line
			if len( line ) == 0: #checking for last line data
				return( read_array ) 
				break #break out of loop because we are at last line in file. 

			return_array = [] #initialize empty array to store chunked line data. 

			while read_boolean:
				line = self.line
				line_data = line.strip().split("\t")
				#print line_data
			
				if len( line ) == 0: #checking for last line data
					if len(return_array) != 0: #if we have data hanging around we need to write before ending
						read_array.append( return_array )
					return( read_array ) 
					break #break out of loop because we are at last line in file. 
				elif line_data[0] == "U":
					current_read = line_data[1] 
					pass #just pass and go to next line for unclassified reads
				else:
					current_read = line_data[1] 
				
				if read == "": #itialize the read variable with real data
					read = current_read #assign sequnce name to read var
					if line_data[0] == "C":
						return_array.append( line_data ) 
				elif current_read == read: #we are still processing tie for same read. 
					return_array.append( line_data )
				else: 
					read_boolean = False
					if line_data[0] == "C":
						#print "bitchin"
						#print return_array
						if len(return_array) == 0: #this is special case for reads without tie.
							return_array.append( line_data ) 
						read_array.append( return_array ) 
					break #break out of loop because we are at last line in file. 

				self.line = self.f.readline()
			
		return( read_array )


def process_individual_classify_file_tie( data ): 
	"""
	#iterates over file and counts the number of times a read is 
	#classfied to each unique ID number. This loop parses output from single
	#classfier database run. 
	
	data - taxonomer object data structure with returned classified lines in list format
	id_dict - parsed key file results with counter for read count tracking
	"""
	#print data
	#pull results from itertools
	results = data[0]
	id_dict = data[1]
	
	#process the object and append counters
	RV = copy.deepcopy(id_dict) #dictionary of returned results
	for result in results: #iterate over lines. 
		result = result.split("\t")
		if len(result) == 1: #we are not working with tied output 
			#rd = individual_read( result[0] ) #load classified read object
			#RV[ rd.taxid ][1] += 1 
			taxid = result[0][2]
			RV[ taxid ][1] += 1
		else: #we are working with tied output and need to append counts to all ties accordingly 
			for rds in result: #iterate over ties...
				taxid = result[0][2]
				RV[ taxid ][1] += 1
				#rd = individual_read( rds ) #load classified read object
				#RV[ rd.taxid ][1] += 1 #append counter

	return( RV ) #return modified dictioanry with updated counts






#################
#MAIN
#################

#for turning into import and functional call:
#def main(key_file, classify_file, proc_num = None, server_output = False, db = None )
if 1:

	#load in the tax IDs
	id_dict = load_keys( key_file )
	#load in the classify data into iterator object
	if server_output:
		classify_data = classify_results_server(class_file, database, chunksize=10000 )
	else:
		classify_data = classify_results(class_file, chunksize=10000 )


	#parse the number of processes to enact
	if arg.proc == None:
	 	proc_num = 1
	else:
		proc_num = arg.proc
	#proc_num = 1 ##r
	###
	#setup multiprocessing for the classification of SVs
	###
	
	p = mp.Pool( processes = proc_num )
	sys.stderr.write("...running parent process with job id %d \n can use this ID to exit \n" %(os.getpid() ) )
	"""
	now we process the classifier file. need to choose which function based on whether server or single
	database output (important variables are in different places). this allows for us to process the files
	extremely quickly using multiprocessing 
	"""
	if server_output:
		class_results = itertools.imap(process_classify_file, itertools.izip( classify_data, itertools.repeat(id_dict) ) )
	else:
		class_results = itertools.imap(process_individual_classify_file, itertools.izip( classify_data, itertools.repeat(id_dict) ) )
		#class_results = p.apply_async(process_individual_classify_file, itertools.izip( classify_data, itertools.repeat(id_dict) ) )
	


	#create a simple count data for piecing everything back together. 
	#count_dict = dict( zip( id_dict.keys() , [cnt[1] for cnt in id_dict.values()] ) )
	#count_dict[ tax id ] = counter 

	count_data = [ np.array( id_dict.keys() ) , np.array( [cnt[1] for cnt in id_dict.values()] )  ]
	#count_data = [ array(dictionary keys), array( counts ) ] 
	#allows for rapid boolean lookup and counter appends.



	start_time = timeit.default_timer() #for looking at performance

	#iterate over the iterator results sum across all processed dictionaries. 
	data_results = []
	for r in class_results:
		#print r
		#iterate over the object and parse the counts into count data 
		count_data = proc_counts( count_data, r )
	
	elapsed = timeit.default_timer() - start_time #for looking at performance

	p.close()


	#now finish by creating a copy of the counts that have been perpetuated
	#through the taxonomic levels (counts classified deeper in taxonomy are
	#collated and added to their parent nodes)
	count_data = parse_taxonomy_counts( count_data, id_dict)

	#and write the output to a file:
	write_count_tables( count_data, id_dict, file_name )


	#final output to std err that the run has finished. 
	sys.stderr.write("...counts processed for %s \n" %( class_file )  )















