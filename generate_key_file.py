import os, sys, argparse
from Bio import SeqIO


parser = argparse.ArgumentParser(description="Creates the taxid map and seq_id map to be used by taxonomer")
parser.add_argument("sti", type=str,  help="sti file" )
parser.add_argument("tri", type=str, help="tri file" )
args = parser.parse_args()

#set the file args from argparse
sti = args.sti
tri = args.tri


#dummy seqs for troubleshooting
#sti = "/data2/eosborne/taxonomer_dbs/strep_MLST/strep_pneumo-gene/strep_pneumo-gene_MLST.sti"
#tri = "/data2/eosborne/taxonomer_dbs/strep_MLST/strep_pneumo-gene/strep_pneumo-gene_MLST.tri"

######################################
#functions
######################################

def process_names( names ):
	"""
	makes the names in format similar to Keith's taxonomy method. 
	"""
	p_list = []
	for i in xrange( len( names ) ):
		#print i
		p_list.append( str(i) + "__" + names[i] )

	RV = ";".join(p_list)
	return( RV )


#function needs to recreate the taxonomy by stepping through the lineage and information in the STI
#can store the lineage as we iterate so know the full taxonomy of the particular sequence as well as 
#the taoxnomic rank value (length of taxonomy) 
def parse_taxonomy( seq_id, lineage, key_dictionary ):
	"""
	seq_id - sequence ID from TRI file. unique ID for sequnce
	lineage - the lineage value (next level up in taxonomy) for sequence
	key_dictionary - dictionary to store various seqIDs and their relevant taxonomic info. 
	"""
	if seq_id in sti_dict:
		tax_id = sti_dict[ seq_id ]
		tax_names = [ tax_id ] #list of taxon names
	else:
		tax_id = str( seq_id )
		tax_names = [ tax_id ] #list of taxon names
	tax_numbers = [ seq_id ]
	is_A_list  = [] #store is_A relationships

	while lineage != '1': #forces traversal through the tri file until we get to the root of taxonomy
		#print lineage
		if lineage == '0': #need this to process the root in the tri file. 
			break
		is_A_list = [lineage] + is_A_list
		tax_numbers = [lineage] +  tax_numbers
		if lineage in sti_dict: #we have the next taxonomic representative in the sti file
			tax_id = sti_dict[ lineage ]
			tax_names = [tax_id] + tax_names #append tax_id to front of list
		else: #the taxon does not have a sequence representative. 
			tax_id = str( lineage ) 
			tax_names = [tax_id] + tax_names
		#now process to next lineage
		lineage  = tri_dict[ lineage ] 


	tax_names = ['root'] + tax_names #append tax_id to front of list
	tax_numbers = [lineage] + tax_numbers
	is_A_list = ['0'] + [lineage] + is_A_list

	#now append all of these reuslts to the final dictionary, which will be keyed 
	#off of the tax_numbers list (unique IDs for each taxonomic level.

	for i in xrange( len( tax_numbers ) ):
		id = tax_numbers[i]
		if id in key_dictionary:
			pass
		else:
			parent = is_A_list[i]
			level = i #taxonomic level (how far down in levels are we?)
			names = process_names( tax_names[:i+1] )
			key_dictionary[ id  ] = [ parent, level, names ]

	return( key_dictionary ) 



######################################
#MAIN
######################################


#open and parse sti and tri files using Biopython into dictionary structures. 
#1. STI
sti_dict = {} #dict to store sti results
handle = open(sti, "rU")
for record in SeqIO.parse(handle, "fasta") :
    sti_dict[ str(record.seq) ] =  record.id + "_ID" #allows tracking to specific sti seq and not tri ID
    #especially for numerical STI files for huge databaes. 
handle.close()

#2. TRI
tri_dict = {} #dict to store sti results
handle = open(tri, "rU")
for record in SeqIO.parse(handle, "fasta") :
    tri_dict[ record.id ] =  str(record.seq)
handle.close()



#now we need to parse the tri file relationships into the nested heirarchy oft he taxonomy. 
key_dictionary = {} #dict to store results for key file
for seq_id, lineage in tri_dict.iteritems():
	#print seq_id
	#print lineage
	#print "\n"
	#parse and update taxonomy dictionary. 
	key_dictionary = parse_taxonomy( seq_id, lineage, key_dictionary ) 




keyz = sorted( [int(i) for i in key_dictionary.keys()], key = int ) #numerical sort tax ids
keyz = [ str(k) for k in keyz] #convert back to strings

for k in keyz: #iterate and write to std out
	if k != '0':
		data = key_dictionary[k]
		data = [ str(item) for item in data ] 
		print "%s\t%s" %( k, "\t".join(data) )


	


