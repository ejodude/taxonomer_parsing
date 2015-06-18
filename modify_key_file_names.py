
"""
script cleans up the taxonomer results for the server
databases. translates sti number into the actual taxonomic 
path
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



parser=argparse.ArgumentParser(description="script cleans up the taxonomer results for the server \
databases. translates sti number into the actual taxonomic path. prints results to std out")
parser.add_argument("key_file", type=str, help="file with taxonomer key")
parser.add_argument("tax_file", type=str, help="file with taxonomic information in human readable format")
arg=parser.parse_args()



#for troubleshooting:
# key_file = "/data2/eosborne/taxonomer_projs/tax_parsing/gg_99_key_v2.txt"
# tax_file = "/data/Pneumonia_Samples/greengenes/gg_13_5_otus/taxonomy/99_otu_taxonomy.txt"

key_file = arg.key_file
tax_file = arg.tax_file




#####################
#functions.
#####################




#need to load the tax file into dictionary object with tax id 
#number as the key

def load_tax( tax_file):
	"""
	load in the tax file that dumps taxonomy for STI entries. 
	"""
	tax_dict = {}
	for line in open( tax_file, 'r' ):
		line = line.strip().split("\t")
		tax_dict[ line[0] ] = line[1]
	return( tax_dict )




def parse_taxid( string, tax_dict ):
	if string == "0__root": #parse root
		seqid = "root"
	else: #need to process lookup from sti number to taxa
		seqid = string.split(";")[-1] #grab terminal taxa info
		seqid = seqid.split("__")[-1]
	return( seqid ) 
	

def populate_keys( key_file, tax_dict, final_tax_dict ):
	"""
	function looks for matches between key file and tax dictionary 
	and populates the taxonomy for a given entry up to the root. 
	This way, taxonomic levels not represented by a STI id can still
	be translated into human readable format with dictionary lookup.
	
	returns a dictionary that will translate the terminal taxa of the
	key phylogeny into the human readable taxonomy:
		final_tax_dict[ '7__327390'] =
		'root;k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__; g__; s__'
	"""

	def translate_taxonomy( taxonomy, ktaxonomy, final_tax_dict):
		"""
		iterate over the taxonomy and create dictionary entires for each
		unique tri level and the create the full text based phylogeny
		ie. 
		{'0__root': 'root',
		'1__2': 'root;k__Bacteria',
		'2__10': 'root;k__Bacteria; p__Proteobacteria'}
		"""
		for i in xrange( len(ktaxonomy) ):
			key = ktaxonomy[i] #this will be the dict key (last numerical entry in key taxonomy
			phylogeny = ";".join( taxonomy[:i+1] ) #this will be the summation of the whole
			#human readable taxonomy
			#print key 
			#print phylogeny
			if key not in final_tax_dict:
				final_tax_dict[ key ] = phylogeny
		return( final_tax_dict )


	#iterate over the key file
	for line in open( key_file, 'r' ):
		line = line.strip().split("\t")
		key_taxonomy = line[-1] #key-based taxonomy of tri numbers
		seq_id = parse_taxid( key_taxonomy, tax_dict )

		#check if the seq id is in the tax_dictionary. if so, we can process
		#the seq_id to retrieve the taxonomy info (we only have info for terminal
		#taxa in the STI file, but can derive the full taxonomy from this info. 
		if re.search("_ID", seq_id ): #we have entry derived from an STI ID for lookup
			#process "_ID" from seq_id
			seq_id = seq_id.split("_ID")[0]
			tdat = tax_dict[seq_id] 
			taxonomy = tdat.split(";") #text-based taxonomy
			taxonomy = ["root"] + taxonomy #append root to taxonomy since it's not included
			ktaxonomy = key_taxonomy.split(";") #grab the numerical tri taxonomy from key file
			#now process the two taxonomies and update the final_tax_dict
			final_tax_dict = translate_taxonomy( taxonomy, ktaxonomy, final_tax_dict)

	return( final_tax_dict )



def print_key_update( key_file, final_tax_dict ):
	"""
	function prints the translated results with the 
	same ordering as the original key file
	"""
	#iterate over the key file
	for line in open( key_file, 'r' ):
		line = line.strip().split("\t")
		taxonomy = line[-1].split(";")
		tax_data = final_tax_dict[ taxonomy[ -1 ] ]
		line[-1] = tax_data
		print "%s" %( "\t".join(line) )




#####################
#MAIN
#####################




if 1:
	tax_dict = load_tax( tax_file ) #load the human readable file from database dump 

	#dictionary to store translated keys.
	final_tax_dict = {} 
	final_tax_dict = populate_keys( key_file, tax_dict, final_tax_dict )

	#finally finish by updating the key file 
	print_key_update( key_file, final_tax_dict )












