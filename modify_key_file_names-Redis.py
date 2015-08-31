
"""
script cleans up the taxonomer results for the server
databases. translates sti number into the actual taxonomic 
path

#need to be logged in as redis...
sudo -u redis /bin/bash

Start the server in screen or in another window:
sudo -u redis -H ~redis/bin/redis-cli -s ~redis/var/run/redis.sock

Run the python application
"""



import argparse #std python imports
import os
import sys
import re
import glob
import numpy as np
import redis



###############
#GLOBALs
###############



parser=argparse.ArgumentParser(description="script cleans up the taxonomer results for the server \
databases. translates sti number into the actual taxonomic path. prints results to std out")
parser.add_argument("key_file", type=str, help="file with taxonomer key")
parser.add_argument("--db", type=str, help="Redis database to query")
arg=parser.parse_args()


"""
#for troubleshooting:
key_file = "/data2/eosborne/taxonomer_projs/taxonomer_keys/uniprot90_viral_30.key"
db = "viral"
"""


key_file = arg.key_file
db = arg.db




#####################
#functions.
#####################


def parse_taxonomy( taxonomy ):

	def grab_redis_info(db, seqid, tax_list):
		"""
		Some quick sample redis commands:
		r = redis.Redis(unix_socket_path='/home/redis/var/run/redis.sock')
		
		r.get("viral:name:%s" %() )
		"""
		#perform redis db search for the tax id. 
		name = r.get("%s:name:%s" %( db, seqid ) )
		level = r.get("%s:category:%s" %( db, seqid ) )
		if name == None: #if we have no results
			name = seqid #set id to seqid
			level = len(tax_list) #set level to be the length of tax list
		return( name, level )


	tax_list = []
	for seqid in taxonomy: #iterate over the numerical taxonomy
		seqid = seqid.split("__")[-1]
	
		if seqid == "root": #parse root
			name = "root"
			level = ""
			
		else: #need to process lookup from sti number to taxa
			if seqid not in tax_dict:
				name, level = grab_redis_info(db, seqid, tax_list)
				#print level
				#print name
				tax_dict[ seqid ] = [ name, level ]
			else:
				name = tax_dict[ seqid ][0]
				level = tax_dict[ seqid ][1]
	
		tax_list.append( name )

	return( ";".join(tax_list) ) 



def print_key_update( key_file ):
	"""
	function prints the translated results with the 
	same ordering as the original key file
	"""


	
	#iterate over the key file
	for line in open( key_file, 'r' ):
		line = line.strip().split("\t")
		taxonomy = line[-1].split(";")
		new_taxonomy = parse_taxonomy( taxonomy )
		line[-1] = new_taxonomy #reassign to new taxonomy
		print "%s" %( "\t".join( line ) ) #print to std out




#####################
#MAIN
#####################




if 1:
	#stores the taxonomer dictionary 
	tax_dict = {}
	#open Redis db connection
	r = redis.Redis(unix_socket_path='/home/redis/var/run/redis.sock')
	#finally finish by updating the key file 
	print_key_update( key_file )







