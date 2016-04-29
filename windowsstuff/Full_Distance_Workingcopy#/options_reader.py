# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 16:43:47 2015

@author: csharp
"""
from collections import defaultdict
from optionsclass import OPTIONS
from trail_holder import analyse
#get current datetime and make folder for results/log_file to be placed in
from time import gmtime, strftime
date=strftime("%Y_%m_%d_%H_%M", gmtime())
import os
if os.path.exists(date)==True:
    date=date+"_1"
os.system(str("mkdir "+date))
log_file=open(str(str(date)+"/log_file.txt"),"w")

string=str("""

===============================================================================
                Log_file:"""+date+"""
===============================================================================
""")
log_file.write(string)

#reads in options_file.txt and extracts options
file_handle=open("options_file.txt","r")
options_dict=defaultdict(list)
file_dict={}
accepted_profiles=[]
for line in file_handle:
    if line[0]=="#":
        log_file.write(line)
        if "#NAME:" in line:
            parts=line.replace("\n","").split("\t")
            name=parts[1]
        if "#LEFT:" in line:
            j=OPTIONS()
            options_dict[name]=j
            parts=line.replace("\n","").split("\t")
            j.left=parts[1]
            accepted_profiles.append(parts[1])
	if "#DATABASE" in line:
	    parts=line.replace("\n","").split("\t")
	    database=parts[1]	
        if "#RIGHT:" in line:
            parts=line.replace("\n","").split("\t")
            print parts[1:][:]
            for K in parts[1:][:]:    
                accepted_profiles.append(str(K))
                j.right.append(str(K))
        if "#MIN_DIST" in line:
            parts=line.replace("\n","").split("\t")
            j.min_dist=int(parts[1])
        if "#MAX_DIST" in line:
            parts=line.replace("\n","").split("\t")
            j.max_dist=int(parts[1])
        if "#STRAND" in line:
            parts=line.replace("\n","").split("\t")
            j.strand=parts[1]
        if "#INPUT_FILE:" in line:
            parts=line.replace("\n","").split("\t")
            file_dict[parts[1]]=""
            filename=parts[1]
        if "#OUTPUT_FILE:" in line:
            parts=line.replace("\n","").split("\t")
            file_dict[filename]=parts[1]
print options_dict

if database=="BIGS" or database=="SANGER":
	for i in file_dict.keys():
	    file_input=i
	    file_output=file_dict[i]
	    file_handle=open(file_input,"r") #Hmmer output file
	    file_handleagain=open(file_input,"r") #Hmmer output file
	    file_towrite=open(file_output,"w") # file to write
	    file_handle2=open("iterated_dataset.csv","r") # file for id to species
	    Isolate2species={}
	    for line in file_handle2:
		things=line.split(",")
		Isolate2species[things[0].replace("_A0_contigs.fa","")]=things[5].replace(" ","_")
	    
	    #file_handle2.close()
	    querydict=[]    
	    file_towrite.write(str("Target_name,accession,tlen,query_name,species,query,qlen,E-value,score,bias,#,of,cEvalue,iE-value,score,bias,from,to,from,to,from,to,file"))
	    file_towrite.write("\n")
	    counter=0
	    for lines in file_handleagain:
		if "# Query file:      " in lines:
		    
		    querydict.append(lines.strip("# Query file:      ").strip("\n"))
	    for line in file_handle:
		
		if "# Query file:      " in line:
		    counter+=1
		parts=line.split(",")    
		#if "DNase" in parts[0] or "IMM_good" in parts[0] or "Col_E3_" in parts[0] or "imm" in parts[0] or "Cytotoxic" in parts[0] or "Cloacin_immun" in parts[0] or "_D" in parts[0]: # names of hmm profile used
		split_by_space=parts[0].split(" ")
		if split_by_space[0] in accepted_profiles:
		    party=[]
		    party2=parts[0].split(" ")
		    contig=party2[18]
		    for i in range(0,len(party2)):
		        if party2[i] != "":
		            party.append(party2[i]) 
		    
		    part=party[:3]
		    contig_length=party[5]
		    contig=party[3]
		    evalue=party[6]
		    HNHstart=party[17]
		    HNHstop=party[18]

		    stuff=party[3].split("_")


		    if len(stuff)==6:
		        part.append(str(stuff[0]+"_"+stuff[1]))
		    if len(stuff)==7:
		        part.append(str(stuff[0]+"_"+stuff[1]+"_"+stuff[2]))
		    if len(stuff)==4:
			part.append(str(stuff[0]+"_"+stuff[1]+"_"+stuff[2]))
			frame=stuff[3]
			
		    try:
			try:
		        	part.append(Isolate2species[str(stuff[0]+"_"+stuff[1])])
			except:
				part.append('NO_SPECIESDATA')
		        part.append(contig)
		        part.append(contig_length)
		        #print parts[0]
		        #print part[3], Isolate2species[parts[3]]   
		        
		        
		       # part.append(str(stuff[0]+"_"+stuff[1]+"_"+stuff[2]+"_"+stuff[3]+"_"+stuff[4]+"_"+stuff[5]))  
		       # part.append(querydict[counter])
		        
		        #print part
		        for i in range(5,len(party)):
		            part.append(party[i])
		       
		        part[7]=evalue
		        part[18]=HNHstart
		        part[19]=HNHstop
			if database=="SANGER":
				part[5]=str(querydict[counter].replace("_6frame.fa","")+"_"+frame)
		        for i in range(len(part)):
		            if part[i] != "":
		                file_towrite.write(part[i].strip("\n"))
		                file_towrite.write(",")
		        file_towrite.write("\n")
			
		    except:
		        pass
	    file_handle.close() #Hmmer output file
	    file_handleagain.close() #Hmmer output file
	    file_towrite.close()# file to write
	    file_handle2.close() # file for id to species   
else:
	profile_string=""
	for i in options_dict.values():
            profile_string=profile_string+str(i.left+" "+i.right[0]+" ")
	for i in file_dict.keys():
	    parse_command=str("./dist/practice_parser "+ i +" "+  file_dict[i] +" "+profile_string)
	    print parse_command
	    os.system(parse_command)
if database=="SANGER":
	database="PATRIC"
object_array=[]
for i in file_dict.values():
    object_array=analyse(object_array,i, database, options_dict)

log_file.write(str("Number of hits to profiles in parsed hmmer file:\t"+str(len(object_array))))
    
from trail_holder import plot_distance
for i in options_dict.keys():
    print i
    [forward,reverse]=plot_distance(object_array, i, date)

from function_holder import addBIGSspecies
object_array=addBIGSspecies(object_array)

from trail_holder import analyse_dist
for i in options_dict:
    if i=="0":
        continue
    log_file.write(str("\nFor the profile:\t"+i))
    
    print "For the profile:",i
    starts=[]
    stops=[]
    response=0
    response=raw_input("How many regions do you want to accept?")
    log_file.write(str("Acecepted ranges:"+str(response)))
    for j in range(int(response)):
        starts.append(raw_input("Define start index:"))
        stops.append(raw_input("Define stop index:"))
        
    
    [Results_per_region,matched_array] = analyse_dist(object_array,i, starts, stops, options_dict[i].strand)
    for i in range(len(Results_per_region)):
        log_file.write(str("\n\nRegion:"+str(i+1)+"\n"+"Start_index:"+str(starts[i])+"\n"+"Stop_index:"+str(stops[i])+"\n"+"Results in region:"+str(Results_per_region[i])+"\n"))

#remove duplicate results
matched_dictionary=defaultdict(list)
for i in matched_array:
    string=str(i.identifier+"_"+str(i.HNHstart)+"_"+str(i.IMMstart)+i.HNHcontig)
    matched_dictionary[string].append(i)

clean_array=[]
for i in matched_dictionary:
    clean_array.append(matched_dictionary[i][0])
matched_array=clean_array        
log_file.write(str("Total number of hits found:"+str(len(matched_array))))
        
from trail_holder import sequence_finder5
sequence_finder5(matched_array)

log_file.write("\nUsed getorf (EMBOSS) to find protein sequences within open reading frames")
counter_in=0
counter_out=0
counter_nofile=0
for i in matched_array:
    if i.finder=="FORCED":
        counter_out+=1
    elif i.finder=="":
        counter_in+=1
    elif i.finder=="NO_FILE":
        counter_nofile+=1
log_file.write("\nNumber of hits found in ORFs:"+str(counter_in))
log_file.write("\nNumber of hits found outside of ORFs:"+str(counter_out))
log_file.write("\nNumber of hits with no file:"+str(counter_nofile))

filename=raw_input("Create fasta files of hits for left profile, in file:")
filename=str(str(date)+"/"+filename)
file_handle=open(filename,"w")
counter=0
for i in matched_array:
     if i.HNHsequence!="":
         string=str(">"+i.identifier+"_"+str(counter)+"_"+str(i.HNHstart)+"\n"+i.HNHsequence+"\n")
         file_handle.write(string)
         counter+=1
 
file_handle.close()
import os 
if os.path.exists(filename)==True:
    print """
===============================================================================

        """, filename ,"""created!
        
===============================================================================
"""
else:
    """

===============================================================================

            Error there was a problem making your file

==============================================================================="""

filename=raw_input("Create fasta files of hits for right profile, in file:")
filename=str(str(date)+"/"+filename)
file_handle=open(filename,"w")
counter=0
for i in matched_array:
     if i.HNHsequence!="":
         string=str(">"+i.identifier+"_"+str(counter)+"_"+str(i.HNHstart)+"\n"+i.IMMsequence+"\n")
         file_handle.write(string)
         counter+=1
 
file_handle.close()

if os.path.exists(filename)==True:
    print """
===============================================================================

        """, filename ,"""created!
        
===============================================================================
"""
else:
    """

===============================================================================

            Error there was a problem making your file

==============================================================================="""

#make a fasta file of sequences and unique headers
from trail_holder import make_addprofiles
make_addprofiles(str(date+"/HNHsequences.fa"), matched_array)
#run fasta file against PFAM-A using HMMscan

import cPickle






from trail_holder import search_profiles
search_profiles(str(date+"/HNHsequences.fa"), str(date+"/HNHsequences.txt"))
#add the profiles to matched_array2
from trail_holder import addprofiles
addprofiles(str(date+"/HNHsequences.txt"), matched_array)
#opens a gui to allow the user to select the allowed domains

for i in matched_array:
    i.profiles=[]
from trail_holder import addprofiles
addprofiles(str(date+"/HNHsequences.txt"), matched_array)
from trail_holder import choose_profiles
profile_list=choose_profiles(matched_array)

pickle_file_binary=str(date+"/results_binary.pickle")
pickle_file=str(date+"/results.pickle")
with open(pickle_file_binary, "wb") as output_file:
    cPickle.dump(matched_array, output_file)
with open(pickle_file, "w") as output_file:
    cPickle.dump(matched_array, output_file)

log_file.write("\nPFAM profiles detected in protein and passed as allowed: ")
for i in profile_list:
    log_file.write(str("\n"+i))
    
match_it=[]
for i in matched_array:
    i.MATCH="TRUE"
    for j in i.profiles:
        if j not in profile_list:
            i.MATCH="FALSE"
seqdict=defaultdict(list)
speciesdict=defaultdict(int)
isolatedict=defaultdict(int)
for i in matched_array:
    if i.MATCH=="TRUE" and len(i.HNHsequence)>350 and len(i.HNHsequence)<950 and i.finder!="FORCED" and i.database=="BIGS":
        seqdict[i.HNHsequence].append(i.profiles)
        speciesdict[i.Species]+=1
        isolatedict[i.isolate]+=1
        match_it.append(i)
log_file.write(str("\n\nNumber of HNH sequences:"+str(len(seqdict))))
log_file.write(str("\n\nNumber of Species:"+str(len(speciesdict))))
log_file.write(str("\n\nNumber of isolates:"+str(len(isolatedict))))
print len(seqdict), "Number of HNH sequences"
print len(speciesdict), "Number of Species"
print len(isolatedict), "Number of isolates"    
    

log_file.close()
