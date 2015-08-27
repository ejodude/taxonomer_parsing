import os, sys, argparse
from Bio import SeqIO
import copy
import itertools
import multiprocessing as mp #allows for parallelization of the classification to speed up script.
import numpy as np
import timeit
import re
import gzip



parser=argparse.ArgumentParser(description="script grabs the reads from fastq file that classify to a taxonomic node and its children")
parser.add_argument("key_file", type=str, help="taxonomer key file")
parser.add_argument("class_file", type=str, help="classifier output file")
parser.add_argument("fastq_file", type=str, help="classifier output file")
parser.add_argument("--taxid", type=str, help="read classified to taxonomer id for parsing")
parser.add_argument("--db", type=str, help="database to parse for server output.")
parser.add_argument("--server_out", dest='server_out', action='store_true', help ="boolean for whether results are server database (if custom db, do not use arg)")
parser.set_defaults(server_out=False)
parser.add_argument("--proc", type=str, help="optional arg for number of proceses to run with classifier; higher thread number will increase speed of classifier; defaults to 1", default = 1)
arg=parser.parse_args()



"""

key_file = "/home/eosborne/tools/taxonomer_parsing/gg_99_key-EDIT.txt"
fastq_file = "/data2/eosborne/taxonomer_projs/MS_RNA_seq/fastqs/11168X9_S043-HC_150129_D00294_0156_AC6E82ANXX_6-pe_combo.fa.gz"
class_file = "/data2/eosborne/taxonomer_projs/MS_RNA_seq/fastqs/11168X9_S043-HC_150129_D00294_0156_AC6E82ANXX_6-pe_combo.std_classifier.txt"
taxid = '21'
output = "test_seq.fq"
database = "bacterial"
server_output = True

"""

key_file = arg.key_file
class_file = arg.class_file
fastq_file = arg.fastq_file
server_output = arg.server_out
database = arg.db
taxid = arg.taxid




#-------------------------------------------------------------------------------------#
#FUNCTIONS
#-------------------------------------------------------------------------------------#



#function needs to recreate the taxonomy by stepping through the lineage and information in the STI
#can store the lineage as we iterate so know the full taxonomy of the particular sequence as well as 
#the taoxnomic rank value (length of taxonomy) 
def parse_children( tax_id, id_dict, children_dict ):
	"""
	taxid - sequence ID from TRI file. unique ID for sequnce
	id_dict - populated dataset from key file.
	children_dict - dictionary to store children relationships. 
	"""
	#test 
	#tax_id = '460830'


	tax_names = []
	tax_numbers = [ ]
	is_A_list  = [] #store is_A relationships

	while tax_id != '1': #forces traversal through the tri file until we get to the root of taxonomy
		#print parent
		if tax_id == '0': #need this to process the root in the tri file. 
			break

		#print tax_id
	
		is_A_list = [tax_id] + is_A_list
		tax_numbers = [tax_id] +  tax_numbers
		tax_name = id_dict[ tax_id ][0][2] 
		tax_names = [tax_name] + tax_names #append tax_id to front of list
	
		#finally append info for the root of the taxonomy. 
		if tax_id not in children_dict: #create child entry. 
			children_dict[ tax_id ] = set( tax_numbers )
		else:
			t_num = children_dict[ tax_id ]
			children_dict[ tax_id ] = t_num  | set( tax_numbers ) #join children with set operator
		
		#now process to next lineage
		tax_id  = id_dict[ tax_id ][0][0]

	#finally append info for the root of the taxonomy. 
	is_A_list = [tax_id] + is_A_list
	tax_numbers = [tax_id] +  tax_numbers
	tax_name = id_dict[ tax_id ][0][2] 
	tax_names = [tax_name] + tax_names #append tax_id to front of list

	if tax_id not in children_dict: #create child entry. 
		children_dict[ tax_id ] = set( tax_numbers )
	else:
		t_num = children_dict[ tax_id ]
		children_dict[ tax_id ] = t_num  | set( tax_numbers ) #join children with set operator

	return( children_dict ) 



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
	
#function needs to recreate the taxonomy by stepping through the lineage and information in the STI
#can store the lineage as we iterate so know the full taxonomy of the particular sequence as well as 
#the taoxnomic rank value (length of taxonomy) 
def parse_children_of_taxa( tax_id, id_dict, children_dict, taxid ):
	"""
	taxid - sequence ID from TRI file. unique ID for sequnce
	id_dict - populated dataset from key file.
	children_dict - dictionary to store children relationships. 
	taxid - focal tax id to parse
	"""
	#test 
	#tax_id = '460830'


	tax_names = []
	tax_numbers = [ ]
	is_A_list  = [] #store is_A relationships

	while tax_id != '1': #forces traversal through the tri file until we get to the root of taxonomy
		#print parent
		if tax_id == '0': #need this to process the root in the tri file. 
			break

		#print tax_id
	
		is_A_list = [tax_id] + is_A_list
		tax_numbers = [tax_id] +  tax_numbers
		tax_name = id_dict[ tax_id ][0][2] 
		tax_names = [tax_name] + tax_names #append tax_id to front of list
	
		#finally append info for the root of the taxonomy. 
		if tax_id == taxid: #create child entry. 
			t_num = children_dict[ tax_id ]
			children_dict[ tax_id ] = t_num  +  tax_numbers #join children with set operator
			break
		#now process to next lineage
		tax_id  = id_dict[ tax_id ][0][0]
	
	return( children_dict ) 



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
	children - taxids that we want to retrieve sequences for
	"""
	#print data
	#pull results from itertools
	results = data[0]
	children = data[1]
	
	#process the object and append counters
	read_list = []
	for result in results: #iterate over lines. 
		result = result.split("\t")
		taxid = result[3]
		if taxid in children:
			#print result[2]
			read_list.append( result[2] )
	return( read_list ) #return modified dictioanry with updated counts



def process_individual_classify_file( data ): 
	"""
	#iterates over file and counts the number of times a read is 
	#classfied to each unique ID number. This loop parses output from single
	#classfier database run. 
	
	data - taxonomer object data structure with returned classified lines in list format
	children - taxids that we want to retrieve sequences for
	"""
	#print data
	#pull results from itertools
	results = data[0]
	children = data[1]
	
	#process the object and append counters
	for result in results: #iterate over lines. 
		result = result.split("\t")
		taxid = result[2]
		if taxid in children:
			read_list.append( result[1] )

	return( read_list ) #return modified dictioanry with updated counts




#-------------------------------------------------------------------------------------#
#MAIN
#-------------------------------------------------------------------------------------#
	

########
# load and process the tax ids and children of tax id
########

#load in the tax IDs from key file
id_dict = load_keys( key_file )



#figure out the parent - child relationships from key file...
children_dict = {}
children_dict[ taxid ] = [] #dict to store children taxa of a given taxID
for tax_id in id_dict.keys():
	#parse and update taxonomy dictionary. 
	children_dict = parse_children_of_taxa( tax_id, id_dict, children_dict, taxid ) 

#uniquify the children list. 
children = set(children_dict[taxid])


########
# process iterator for taxonomer output and retrieve the fastq ids that
# classify to desired sequences. 
########


#now we need to process the classifier output. 
#load in the classify data into iterator object
if server_output:
	classify_data = classify_results_server(class_file, database, chunksize=10000 )
else:
	classify_data = classify_results(class_file, chunksize=10000 )


#parse the number of processes to enact
if arg.proc == None:
	proc_num = 1
else:
	proc_num = int( arg.proc )


p = mp.Pool( processes = proc_num )
sys.stderr.write("...running parent process with job id %d \n can use this ID to exit \n" %(os.getpid() ) )
"""
now we process the classifier file. need to choose which function based on whether server or single
database output (important variables are in different places). this allows for us to process the files
extremely quickly using multiprocessing 
"""
if server_output:
	class_results = itertools.imap(process_classify_file, itertools.izip( classify_data, itertools.repeat(children) ) )
else:
	class_results = itertools.imap(process_individual_classify_file, itertools.izip( classify_data, itertools.repeat(children) ) )


#iterate over the iterator results sum across all processed dictionaries. 
fastq_ids = []
for r in class_results:
	#print r
	#iterate over the object and parse the counts into count data 
	fastq_ids = fastq_ids + r 

#sort the fastq ids.. will be easier to do lookup this way. 
fastq_ids = sorted( fastq_ids )



########
# process the input file into a single line. 
########

with gzip.open(fastq_file, 'rb') as f:
    for line in f:
    	if line[1:-1] in fastq_ids:
    		print line.strip()
    		print f.readline().strip()
    		
    		










