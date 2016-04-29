# -*- coding: utf-8 -*-
"""
Created on Thu May 01 13:30:26 2014

@author: kell3284
"""

from Bio.Seq import Seq
#from Bio.SeqUtils import GC
#from Bio.Alphabet import generic_dna
#from Bio.Alphabet import IUPAC
#from Bio.Data import CodonTable
#from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
#from Bio.Blast import NCBIWWW
#from Bio.Blast import NCBIXML
from collections import defaultdict
from Bio.SeqRecord import SeqRecord

import os

############################################################################
    
############################################################################ 

def musclealign(Spec):
    import os
    command=str("./../Programs/MUSCLE/muscle -in "+Spec+".fa -out " +Spec+".afa")
    os.system(command)
    command_tree=str("./../Programs/MUSCLE/muscle -maketree -in "+Spec+".afa -out "+ Spec+".phy -cluster neighborjoining")
    os.system(command_tree)
    

############################################################################
    
############################################################################ 
     
#[isolate_seq, isolate_start]=BIGS_reader("BIGSdb-1.fasta", "BIGSdbtable-1.csv")
#[isolate_seqimm, isolate_startimm]=BIGS_reader("BIGSdb2-1.fasta", "BIGSdbtable2-1.csv")
#comparestarts(isolate_seq, isolate_start,isolate_seqimm, isolate_startimm,"sequencefile2.fasta")


#example maplab style bar chart
def make_bar_chart(Totals, matches, labels, title, yaxis):
    import numpy as np
    import matplotlib.pyplot as plt
    N=len(Totals)
    
    ind=np.arange(N)
    width=0.45
    fig, ax = plt.subplots()
    rects1 = ax.bar(ind, Totals, width, color='r') 
    rects2 = ax.bar(ind+width, matches, width, color='y')
    ax.set_ylabel(yaxis)
    ax.set_title(title)
    ax.set_xticks(ind+width)
    ax.set_xticklabels( labels, fontsize='10' )
    
    for label in ax.xaxis._ticklabels():
        label.set_rotation(90)
    
    ax.legend( (rects1[0], rects2[0]), ('Total isolates', 'Isolates containing a colicin') )
    plt.show()

############################################################################
    
############################################################################ 

def convert_for_tree(tree_file, tree_output, tree_data_output, Total_species, Matching_species ):
    
    
    from Bio import Phylo
    tree=Phylo.read(tree_file,"newick")
    tree=tree.as_phyloxml()    
  
    for clade in tree.get_terminals():
        clade.name=clade.name.replace(" ", "_")
        parts=clade.name.split("_")
        print parts
        clade.name=str(parts[0]+"_"+parts[1])
        
    
    Phylo.write(tree, tree_output, 'phyloxml')
    file_handle=open(tree_data_output, "w")
    counter=0
    file_handle.write("LABELS,mylabel_1,mylabel_2")
    file_handle.write("\n")
    file_handle.write("COLORS,#ff0000,#00ff00")
    file_handle.write("\n")
    for i in Total_species.keys():
        if i == "1986":
            continue
        name=str(i)
        number=((Matching_species[str(name)]/Total_species[str(name)])*500)
        #print(name, Matching_species[name])
        i=i.replace(" ", "_")
        parts=i.split("_")
        print parts
        try:
            i=str(parts[1]+"_"+parts[3])
        except:
            pass
        line= str(i+","+str(Total_species[name])+","+str(number))
        file_handle.write(line)
        file_handle.write("\n")
    #    if counter > 30:
    #        break
        counter+=1
    file_handle.close()
    return tree
    
    
    
############################################################################
    
############################################################################    
       
def sequence_finder2(start, immposition, immcontig,  typ,database):
#will return translated sequences of proteins found in the start variable from the files stored on harddrive. If files are not on harddrive directory must be changed-found in function holder.
    proteindict =defaultdict(list) #returns protein dictionary
    import sys
    seqdict=defaultdict(list)
    framesss=defaultdict(list)
    if database=="PATRIC":
        directory="/run/media/csharp/Dell_USB_Portable_HDD/ftp.patricbrc.org/patric2/genomes/genomes"
        
    else:
        directory= "/run/media/csharp/Dell_USB_Portable_HDD/contigs_ISOLATES_RefSet4"
        
    transeq="./../Emboss/EMBOSS-6.6.0/emboss/transeq -sequence" # runs transeq on the command line
    transeqend=" -outseq $p'_6frame.fa' -clean -table 11 -frame 6"   
    
    getorf="./../Emboss/EMBOSS-6.6.0/emboss/getorf -sequence"
    getorfend=" -outseq 'tester.txt' -table 11"
    
    
    if typ=="HNH":
        numb=0
    else:
        numb=1
    counter=0
    for things in start.keys():
        
        isolate_search=things
        bounter=0
        while bounter<99:
            try:
                endpos=int(immposition[str(isolate_search+" "+ str(start[isolate_search][bounter][numb]))][0])
            
                #print HNHposition[str(isolate_search+" "+ str(start[isolate_search][bounter][0]))]
                startpos=int(start[isolate_search][bounter ][numb])
                
                contig=immcontig[str(isolate_search+" "+ str(start[isolate_search][bounter][numb]))][0]
                if typ=="HNH":
                    temp=endpos
                    endpos=startpos
                    startpos=temp
                #print startpos, endpos, contig
                
                if database==("PATRIC"):
                   
                   parts=contig.split("/")
                   filename=str(parts[7])
                   filename2=str(parts[7]+"/"+parts[7]+".fna")
                   
                else:
                    parts=contig.split("_")
            
                    filename2=str(parts[0]+"_"+parts[1]+"_"+parts[2]+"_contigs.fa")
                    filename=str(parts[0]+"_"+parts[1]+"_"+parts[2]+"_contigs.fa")
                
                collect=""
                frame=int(contig[-1:])
                
                bounter+=1
                
                print filename
                for file in os.listdir(str(directory)):
                    
                    if file== filename:
                        
                        file_handle=open(str(directory+"/" + filename2),"r")
                        
                        for item in file_handle:
                            collect=collect+item
                           
                        parts=collect.split(">") 
                        for i in parts:
                           
                                
                            if database != "PATRIC":
                                if i[:len(contig[:-2])]==contig[:-2]:
                                    #print i[-1:len(contig[:-2]):-1]
                                    string_safe=str(i[len(contig):])
                                    print "HHHHHHHH" 
                                   
                                    break
                            else:
                                
                                if isolate_search in i and database == "PATRIC":
                                    j=i.split("]\n")
                                    
                                    #print i[-1:len(contig[:-2]):-1]
                                    string_safe=str(j[1])
                        print string_safe[:10]
                        file_handle.close()
                        file_handle8=open("seqtest.fa","w")
                        file_handle8.write(str(">Seqtester\n"))
                        file_handle8.write(string_safe)
                        file_handle8.close()
                        
                        command=str(transeq+" "+directory+"/" + filename2+transeqend)
                        #command_getorf=str(getorf+" "+directory+ "/" + filename2+getorfend)
                        command_getorf=str(getorf+" "+"seqtest.fa"+getorfend)

                        os.system(command_getorf)
                        os.system(command)

                        file_handle=open(str("_6frame.fa"),"r")
                        
                        for item in file_handle:
                            collect=collect+item
                            
                        parts=collect.split(">") 
                        
                        for i in parts[1:]:
                            
                                
                            if database != "PATRIC":
                                if i[:len(contig)]==contig:
                                    #print i[-1:len(contig[:-2]):-1]
                                    string=str(i[len(contig):]).replace("\n","")
                                    string=string.replace(">","")
                                    string=string.replace(" ","")
                                    
                                    break
                            else:
                                
                                
                                if str(i[:len(isolate_search)+2])==str(isolate_search+"_"+str(frame)) and database == "PATRIC":
                                    j=i.split("]\n")
                                    string=str(j[1]).replace("\n","")
                                    string=string.replace("\t","")
                                    
                                    break
                        file_handle.close()            
                        
                        
                        if frame<100:
                            seqdict[things].append(contig)
                           
                            #print str(protein[startpos-1:endpos+1])
                            framesss[things].append(str(string[startpos-700:endpos]))
                            counter+=1
                            print counter
                            
                            collect2=""
                            file_handle5=open("tester.txt",'r')
                            
                            for items in file_handle5:
                                collect2=collect2+items
                            parts2=collect2.split(">")
                                
                            for k in parts2[1:]:
                                
                                l=k.replace("\n","")
                                if str(string[startpos:endpos]) in l:
                                    N=l.replace("(REVERSE SENSE)","")
                                    M=N.split("]")
                                    M=M[1].replace("\n","")
                                    
                                    proteindict[things].append(M)
                                    y2="rm rm -f 'tester.txt'"
                                    os.system(y2)
                                    file_handle5.close()
                                    break
                                    
                        else:
                            break
                        
                        y="rm rm -f '_6frame.fa'"
                       
                        os.system(y)
                        
            except:
                sys.exc_info()[0]
                break
            
    return proteindict, seqdict, framesss

############################################################################
    
############################################################################ 


def convert_for_tree2(tree_file, tree_output, tree_data_output, Total_species, Matching_species1, Matching_species2, Both):
    
    
    from Bio import Phylo
    tree=Phylo.read(tree_file,"newick")
    tree=tree.as_phyloxml()    
    
    
    


    
    
    
    for clade in tree.get_terminals():
        clade.name=clade.name.replace(" ", "_")
        parts=clade.name.split("_")
        print parts
        clade.name=str(parts[0]+"_"+parts[1])
        
    
    Phylo.write(tree, tree_output, 'phyloxml')
    file_handle=open(tree_data_output, "w")
    counter=0
    file_handle.write("LABELS,mylabel_1,mylabel_2, mylabel_3")
    file_handle.write("\n")
    file_handle.write("COLORS,#ff0000,#00ff00,#0000ff")
    file_handle.write("\n")
    for i in Total_species.keys():
        if i == "1986":
            continue
        name=str(i)
        number=((Matching_species1[str(name)]/Total_species[str(name)])*100)
        number2=((Matching_species2[str(name)]/Total_species[str(name)])*100)
        number3=((Both[str(name)]/Total_species[str(name)])*100)
        #print(name, Matching_species[name])
        i=i.replace(" ", "_")
        parts=i.split("_")
        #print parts
        try:
            i=str(parts[1]+"_"+parts[3])
        except:
            pass
        line= str(i+","+str(str(number)+","+str(number2)+","+str(number3)))
        file_handle.write(line)
        file_handle.write("\n")
    #    if counter > 30:
    #        break
        counter+=1
    file_handle.close()
    return tree
############################################################################
    
############################################################################     

def make_hist(data, bino):
    import seaborn as sns
    import matplotlib.pyplot as plt
    
    sns.set(style="white", palette="muted")
    f, ax1 = plt.subplots()
    #sns.despine(left=True)
    
    rs = data
    
    b, g, r, p = sns.color_palette("muted", 4)
    
    d = rs
    
    sns.distplot(d, kde=False, bins=bino, color=b, ax=ax1)
    
    
    plt.setp(ax1)
    plt.tight_layout()

############################################################################
    
############################################################################ 

def sequence_finder3(object_array):
#will return translated sequences of proteins found in the start variable from the files stored on harddrive. If files are not on
# harddrive directory must be changed, returns protein translations from STOP to STOP sites and the nucleotide translations are
#included as contig
    import sys
    import os
    
    for search in object_array:
        database=search.database                    #different syntax between PATRIC and Bigs
        isolate_search=search.identifier
        print isolate_search
        if database=="PATRIC":
            directory="/run/media/csharp/Dell_USB_Portable_HDD/ftp.patricbrc.org/patric2/genomes/genomes"
            
        else:
            directory= "/run/media/csharp/Dell_USB_Portable_HDD/contigs_ISOLATES_connor"
            
        transeq="./../Emboss/EMBOSS-6.6.0/emboss/transeq -sequence" # runs transeq on the command line
        transeqend=" -outseq $p'_6frame.fa' -clean -table 11 -frame 6"   
        
        getorf="./../Emboss/EMBOSS-6.6.0/emboss/getorf -sequence"
        getorfend=" -outseq 'tester.txt' -table 11 -find 0" #produces STOP-STOP nucleotide translations
        getorfend2=" -outseq 'tester2.txt' -table 11 -find 2 " # produces nucleotide translations between STOP-STOP
        
        
        endposHNH=search.HNHstop
    
        endposIMM=search.IMMstop
        startposHNH=search.HNHstart
        startposIMM=search.IMMstart
        
        
        HNHcontig=search.HNHcontig
        IMMcontig=search.IMMcontig
        if database==("PATRIC"):
           
           parts=HNHcontig.split("/")
           filename=str(parts[9])
           filename2=str(parts[9]+"/"+parts[9]+".fna")
           
        else:
            parts=HNHcontig.split("_")
    
            filename2=str(parts[0]+"_"+parts[1]+"_"+parts[2]+"_contigs.fa")
            filename=str(parts[0]+"_"+parts[1]+"_"+parts[2]+"_contigs.fa")
        
        collect=""
        frameHNH=int(HNHcontig[-1:])
        frameIMM=int(IMMcontig[-1:])
        
        
        print filename
        for file in os.listdir(str(directory)):
            
            if file== filename:
                
                file_handle=open(str(directory+"/" + filename2),"r")
                
                for item in file_handle:
                    collect=collect+item
                   
                parts=collect.split(">") 
                for i in parts:
                   
                        
                    if database != "PATRIC":
                        if i[:len(HNHcontig[:-2])]==HNHcontig[:-2]:
                           
                            string_safe=str(i[len(HNHcontig):])
                            
                           
                            break
                    else:
                        
                        if isolate_search in i[:50] and database == "PATRIC":
                            print i[:50], isolate_search
                            j=i.split("]\n")
                            
                            string_safe=str(j[1])
                file_handle.close()
                file_handle8=open("seqtest.fa","w")
                file_handle8.write(str(">Seqtester\n"))
                file_handle8.write(string_safe)
                file_handle8.close()
                
                command=str(transeq+" "+directory+"/" + filename2+transeqend)
                command_getorf=str(getorf+" "+"seqtest.fa"+getorfend)
                command_getorf2=str(getorf+" "+"seqtest.fa"+getorfend2)                
                
                os.system(command_getorf)
                os.system(command)
                os.system(command_getorf2)

                file_handle=open(str("_6frame.fa"),"r")
                
                for item in file_handle:
                    collect=collect+item
                    
                parts=collect.split(">") 
                
                for i in parts[1:]:
                    
                        
                    if database != "PATRIC":
                        if i[:len(HNHcontig)]==HNHcontig:
                            #print i[-1:len(contig[:-2]):-1]
                            string=str(i[len(HNHcontig):]).replace("\n","")
                            string=string.replace(">","")
                            string=string.replace(" ","")
                            
                        if i[:len(IMMcontig)]==IMMcontig:
                            #print i[-1:len(contig[:-2]):-1]
                            stringIMM=str(i[len(IMMcontig):]).replace("\n","")
                            stringIMM=stringIMM.replace(">","")
                            stringIMM=stringIMM.replace(" ","")
                            
                    else:
                        
                        test="fail"
                        test2="fail"
                        if str(i[:len(isolate_search)+2])==str(isolate_search+"_"+str(frameHNH)) and database == "PATRIC":
                            j=i.split("]\n")
                            
                            string=str(j[1]).replace("\n","")
                            string=string.replace("\t","")
                            test="PASS"
                            
                            
                        if str(i[:len(isolate_search)+2])==str(isolate_search+"_"+str(frameIMM)) and database == "PATRIC":
                            j=i.split("]\n")
                            stringIMM=str(j[1]).replace("\n","")
                            stringIMM=stringIMM.replace("\t","") 
                            test2="PASS"
                        if test=="PASS" and test2=="PASS":
                        
                            break
                        
                file_handle.close()            
                
                
                
                    
                collect2=""
                file_handle5=open("tester.txt",'r')
                
                for items in file_handle5:
                    collect2=collect2+items
                parts2=collect2.split(">")
                #parts2.sort(key=len, reverse=True)
                test="FAIL"
                test2="FAIL"
                for k in parts2[1:]:
                    if k !="":
                        l=k.replace("\n","")
                        if str(string[startposHNH+5:endposHNH-5]) in l:
                             
                            N=l.replace("(REVERSE SENSE)","")
                            M=N.split("]")
                            positions=M[0].split(" ")
                            print M
                            start=int(positions[1].replace("[",""))
                            end=int(positions[3])
                            search.contig_start=end
                            search.contig_end=start
                            file_tester2=open("tester2.txt","r")
                            tester2_string=str(M[0])
                            collect4=""
                            for line in file_tester2:
                               collect4+=line
                               
                            parts_prom=collect4.split(">")
                            for next_iter in parts_prom:
                                if tester2_string in next_iter:
                                    parts_prom_again=next_iter.split("]")
                                    string456=parts_prom_again[1].replace("\n","")
                                    search.contig=string456.replace("(REVERSE SENSE)","")
                                    break    
                            file_tester2.close()
                            M=M[1].replace("\n","")
                            test="PASS"
                            search.HNHsequence=M
                        if str(stringIMM[startposIMM+5:endposIMM-5]) in l:
                            N=l.replace("(REVERSE SENSE)","")
                            M=N.split("]")
                            #print M, str(stringIMM[startposIMM+5:endposIMM-5])
                            M=M[1].replace("\n","")
                            test2="PASS"
                            search.IMMsequence=M
                        if test=="PASS" and test2=="PASS":
                           # y2="rm rm -f 'tester.txt'"
                            #os.system(y2)
                            file_handle5.close()
                            break
                        
                if search.HNHsequence=="":
                    search.HNHsequence=str(string[startposHNH-500:endposHNH+4])
                    search.finder="FORCED"
                    if search.HNHsequence=="":
                         search.HNHsequence=str(string[startposHNH-200:endposHNH+4])
                         search.finder="FORCED"
                if search.IMMsequence=="":
                    search.IMMsequence=str(string[startposIMM-40:endposHNH+50])
                
                y="rm rm -f '_6frame.fa'"
               
                os.system(y)

############################################################################
    
############################################################################ 


def object_writer(object_array, filename):
    file_handle=open(filename,'w')
    for i in object_array:
        string=str(">"+ i.identifier+"_"+str(i.HNHstart)+"\n"+i.HNHsequence+"\n")
        file_handle.write(string)
    file_handle.close()

############################################################################
    
############################################################################ 

def addprofiles(filename, OBJ_array):
    matched_array=OBJ_array
    file_handle=open(filename,"r")
    count=0
    for line in file_handle:
        count+=1
        if count>3:
            parts2=[]
            parts=line.split(" ")
            for i in parts:
                if i !="":
                    parts2.append(i)
            
            if float(parts2[4])<float(0.00005):
                
                stuff=parts2[2].split("_")
                #print stuff[:-1]
                start=stuff[-1]
                
                if len(stuff)==3:
                    string=str(stuff[0]+"_"+stuff[1])
                else:
                    string=str(stuff[0])
                print start, string    
                for j in matched_array:
                    if j.HNHstart== int(start) and j.identifier== string:
                        j.profiles.append(parts2[0])
    file_handle.close()

############################################################################
    
############################################################################ 

def analyse(object_array,filename, database):
    from classex1 import ISOLATE
    
    print """\r

===============================================================================

            Analysing """,filename,""" for protein co-ordinates...
            
===============================================================================

"""
    data= open(filename, 'r')
    
    
    e_values=[]  # list of e_values 
    counter=0 # counter to miss first line
     #everything in table matched or un matched with e_value over 0.000005
    
    Species={}#dictionary of Species
    for item in data:
        counter+=1
        parts=item.split(",")
        if counter>2 and float(parts[7])<float(0.000005):
            y=ISOLATE()
            y.database=database

            e_values.append(float(parts[7]))
            if database=="PATRIC":
                parts[4]=(parts[4]).replace("/run/media/csharp/Dell_USB_Portable_HDD/ftp.patricbrc.org/patric2/genomes/genomes/","/home/csharp/PATRIC/ftp.patricbrc.org/patric2/genomes/")
                shortened_isolate=parts[4].split("/")
                #print shortened_isolate
                short_isolate=shortened_isolate[7]
                Species[parts[3]]=short_isolate
                y.identifier=parts[3]
                y.isolate=short_isolate
                #print y.isolate
                #y.speciesWriter()
            else:
                Species[parts[3]]=parts[4]
                y.isolate=parts[3]
                y.identifier=parts[3]
                y.Species=(parts[4])
            if "Colicin-DNase" == parts[0]:
                y.HNHstart=int(parts[18])
                y.HNHstop=int(parts[19])
                y.HNHcontig=(parts[5])
                y.coltype="DNase"
                y.contiglength=int(parts[6])
                y.MATCH="HNH"
            if "short-DNases" in parts[0]:
                y.HNHstart=int(parts[18])
                y.HNHstop=int(parts[19])
                y.HNHcontig=(parts[5])
                y.coltype="DNase"
                y.contiglength=int(parts[6])
                y.MATCH="HNH"
            if "DNase_bias_profiles" == parts[0]:
                y.HNHstart=int(parts[18])
                y.HNHstop=int(parts[19])
                y.HNHcontig=(parts[5])
                y.coltype="BIAS"
                y.contiglength=int(parts[6])
                y.MATCH="HNH"
            

            elif "Colicin_D" in parts[0]:
                y.HNHstart=int(parts[18])
                y.HNHstop=int(parts[19])
                y.HNHcontig=(parts[5])
                y.contiglength=int(parts[6])
                y.coltype="D"
                y.MATCH="HNH"
                y.contiglength=int(parts[6])
            elif "Cytotoxic" in parts[0]:
                y.HNHstart=int(parts[18])
                y.HNHstop=int(parts[19])
                y.HNHcontig=(parts[5])
                y.contiglength=int(parts[6])
                y.coltype="Cyto"
                y.MATCH="HNH"
                y.contiglength=int(parts[6])
            elif "Cloacin_immun" in parts[0]:
                y.IMMstart=int(parts[18])
                y.IMMstop=int(parts[19])
                y.IMMcontig=(parts[5])
                y.contiglength=int(parts[6])
                y.coltype="Cyto"
                y.MATCH="IMM"
            
            elif "Colicin_im" in parts[0]:
                y.IMMstart=int(parts[18])
                y.IMMstop=int(parts[19])
                y.IMMcontig=(parts[5])
                y.contiglength=int(parts[6])
                y.coltype="D" 
                y.MATCH="IMM"
            elif "Col_Dimm" in parts[0]:
                y.IMMstart=int(parts[18])
                y.IMMstop=int(parts[19])
                y.IMMcontig=(parts[5])
                y.contiglength=int(parts[6])
                y.coltype="Dnnnnnn" 
                y.MATCH="IMM"
            elif "IMM_good_profiles" in parts[0]:
                y.IMMstart=int(parts[18])
                y.IMMstop=int(parts[19])
                y.IMMcontig=(parts[5])
                y.contiglength=int(parts[6])
                y.coltype="DNase" 
                y.MATCH="IMM"
            elif "Colicin_Pyocin" in parts[0]:
                y.IMMstart=int(parts[18])
                y.IMMstop=int(parts[19])
                y.IMMcontig=(parts[5])
                y.contiglength=int(parts[6])
                y.coltype="DNase"
                y.MATCH="IMM"
            y.speciesWriter()
            object_array.append(y)
    
    data.close()
    print """

===============================================================================

            Finished analysing """,filename,"""
            
===============================================================================

"""
    return object_array

############################################################################
    #elif "poreimm" in parts[0] or " Microcin" in parts[0] or "Colicin_im " in parts[0] or "ImmE5" in parts[0]:
############################################################################ 

def analyse_dist2(object_array):
    from collections import defaultdict
    object_dict=defaultdict(list)
    for i in object_array:
        object_dict[i.identifier].append(i)
        
    counter=0
    for K in object_dict.keys():
        counter+=1
        print counter
        for i in object_dict[K]:
            for j in object_dict[K]:
                if i.coltype==j.coltype and i.database==j.database:
                    if i.coltype=="BIAS" or i.coltype=="GOOD" or i.coltype=="D":
                        if i.MATCH=="IMM" and j.MATCH=="HNH":
                            if i.IMMstart-j.HNHstop <30 and i.IMMstart-j.HNHstop>0:
                                i.HNHstart=j.HNHstart
                                i.HNHstop=j.HNHstop
                                i.HNHcontig=j.HNHcontig
                                i.MATCH="TRUE"
                    elif i.coltype=="Pore":
                        if i.MATCH=="IMM" and j.MATCH=="HNH":
                            if i.contiglength-i.IMMstop-j.HNHstop <50 and i.IMMstart-j.HNHstop>28:
                            
                                i.HNHstart=j.HNHstart
                                i.HNHstop=j.HNHstop
                                i.HNHcontig=j.HNHcontig
                                i.MATCH="TRUE"
    counter=0  
    matched_array2=[]              
    for i in object_dict.keys():
        for j in object_dict[i]:
         if j.MATCH=="TRUE":
             matched_array2.append(j)
    return matched_array2

################################################################################################

################################################################################################


def plot_distance(object_array, coltype):
    import matplotlib.pyplot as pl
    import seaborn as sns
    from collections import defaultdict
    object_dict=defaultdict(list)
    for i in object_array:
        object_dict[i.identifier].append(i)
    keys=object_dict.keys()
    
    #finds the isolates containing a HNH and an Imm gene next to each other but with 
    #imm after the hnh returns a dict of lists with all matches
    
    isoforgraphs=[]
    isographs=[]
    reverse1=[]
    reversegraph=[]
    for M in range(0,80,1): 
        counter=0
        counter2=0
        reverse=0
        reverse2=0
        for k in keys:
            for j in object_dict[k]:
                for i in object_dict[k]:
                        if i.coltype==j.coltype and i.database==j.database and i.coltype==coltype:
                            if i.MATCH=="IMM" and j.MATCH=="HNH":
                                if i.IMMstop-j.HNHstart==M:
                                    counter+=1
                                if i.IMMstart-j.HNHstop == M:
                                    counter2+=1
                                if i.contiglength- i.IMMstop- j.HNHstop < M and i.IMMstart-j.HNHstop>0:
                            
                                    reverse+=1
                                if i.contiglength- i.IMMstop -j.HNHstop == M: 
                                    reverse2+=1
                        
        isographs.append(counter)
        isoforgraphs.append(counter2)
        reverse1.append(reverse)
        reversegraph.append(reverse2)
        

    fig, ax1 = pl.subplots()
    font = {'family' : 'normal','weight' : 'bold','size'   : 22}
    pl.rc('font', **font)
    ax1.plot(range(0,80,1),isoforgraphs,'b',range(0,80,1),isographs,'r',range(0,80,1),reversegraph,'g')
#    ax1.plot(isoforgraphs,linewidth=3)
#    ax1.plot(isographs,linewidth=3)
#    ax1.plot(reversegraph, linewidth=3)
    pl.xlabel('Intergenic region (residues)')
    pl.ylabel('Hits')
    labels=["left->right", "right->left", "reverse strand"]
    ax1.legend(labels)
    pl.show()
    return isographs, isoforgraphs
#######################################################################################################################

#######################################################################################################################

def sequence_finder4(object_array):
#will return translated sequences of proteins found in the start variable from the files stored on harddrive. If files are not on harddrive directory must be changed-found in function holder.
    import sys
    import os
    for search in object_array:
        database=search.database
        isolate_search=search.identifier
        print isolate_search
        if database=="PATRIC":
            directory="/run/media/csharp/Dell_USB_Portable_HDD/ftp.patricbrc.org/patric2/genomes/genomes"
            
        else:
            directory= "/run/media/csharp/Dell_USB_Portable_HDD/contigs_ISOLATES_connor"
            
        transeq="./../Emboss/EMBOSS-6.6.0/emboss/transeq -sequence" # runs transeq on the command line
        transeqend=" -outseq $p'_6frame.fa' -clean -table 11 -frame 6"   
        
        getorf="./../Programs/Prodigal-2.6.1/prodigal -i"
        getorfend=" -a tester.faa -d tester2.fa -f gbk -o output_file_waste.txt -s output_prod_genes.txt -p meta "
        
        
        
        endposHNH=search.HNHstop
    
        endposIMM=search.IMMstop
        startposHNH=search.HNHstart
        startposIMM=search.IMMstart
        
        
        HNHcontig=search.HNHcontig
        IMMcontig=search.IMMcontig
        if database==("PATRIC"):
           
           parts=HNHcontig.split("/")
           filename=str(parts[9])
           filename2=str(parts[9]+"/"+parts[9]+".fna")
           
        else:
            parts=HNHcontig.split("_")
    
            filename2=str(parts[0]+"_"+parts[1]+"_"+parts[2]+"_contigs.fa")
            filename=str(parts[0]+"_"+parts[1]+"_"+parts[2]+"_contigs.fa")
        
        collect=""
        frameHNH=int(HNHcontig[-1:])
        frameIMM=int(IMMcontig[-1:])
        
        
        print filename
        for file in os.listdir(str(directory)):
            
            if file== filename:
                
                file_handle=open(str(directory+"/" + filename2),"r")
                
                for item in file_handle:
                    collect=collect+item
                   
                parts=collect.split(">") 
                for i in parts:
                   
                        
                    if database != "PATRIC":
                        if i[:len(HNHcontig[:-2])]==HNHcontig[:-2]:
                           
                            string_safe=str(i[len(HNHcontig):])
                            
                           
                            break
                    else:
                        
                        if isolate_search in i[:50] and database == "PATRIC":
                            print i[:50], isolate_search
                            j=i.split("]\n")
                            
                            string_safe=str(j[1])
                file_handle.close()
                file_handle8=open("seqtest.fa","w")
                file_handle8.write(str(">Seqtester\n"))
                file_handle8.write(string_safe)
                file_handle8.close()
                
                command=str(transeq+" "+directory+"/" + filename2+transeqend)
                command_getorf=str(getorf+" "+"seqtest.fa"+getorfend)
                os.system(command_getorf)
                os.system(command)

                file_handle=open(str("_6frame.fa"),"r")
                
                for item in file_handle:
                    collect=collect+item
                    
                parts=collect.split(">") 
                
                for i in parts[1:]:
                    
                        
                    if database != "PATRIC":
                        if i[:len(HNHcontig)]==HNHcontig:
                            #print i[-1:len(contig[:-2]):-1]
                            string=str(i[len(HNHcontig):]).replace("\n","")
                            string=string.replace(">","")
                            string=string.replace(" ","")
                            
                        if i[:len(IMMcontig)]==IMMcontig:
                            #print i[-1:len(contig[:-2]):-1]
                            stringIMM=str(i[len(IMMcontig):]).replace("\n","")
                            stringIMM=stringIMM.replace(">","")
                            stringIMM=stringIMM.replace(" ","")
                            
                    else:
                        
                        test="fail"
                        test2="fail"
                        if str(i[:len(isolate_search)+2])==str(isolate_search+"_"+str(frameHNH)) and database == "PATRIC":
                            j=i.split("]\n")
                            
                            string=str(j[1]).replace("\n","")
                            string=string.replace("\t","")
                            test="PASS"
                            
                            
                        if str(i[:len(isolate_search)+2])==str(isolate_search+"_"+str(frameIMM)) and database == "PATRIC":
                            j=i.split("]\n")
                            stringIMM=str(j[1]).replace("\n","")
                            stringIMM=stringIMM.replace("\t","") 
                            test2="PASS"
                        if test=="PASS" and test2=="PASS":
                        
                            break
                        
                file_handle.close()            
                
                
                
                    
                collect2=""
                file_handle5=open("tester.faa",'r')
                
                for items in file_handle5:
                    collect2=collect2+items
                parts2=collect2.split(">")
                #parts2.sort(key=len, reverse=True)
                test="FAIL"
                test2="FAIL"
                for k in parts2[1:]:
                    if k !="":
                        l=k.replace("\n","")
                        if str(string[startposHNH+5:endposHNH-5]) in l:
                             
                            N=l.replace("(REVERSE SENSE)","")
                            M=N.split("gc_cont=0.")
                            positions=M[0].split(" ")
                            start=int(positions[2])
                            end=int(positions[4])
                            search.contig_start=end
                            search.contig_end=start
                            tester2_string=positions[8]
                            file_tester2=open("tester2.fa","r")
                            
                            collect4=""
                            for line in file_tester2:
                               collect4+=line
                               
                            parts_prom=collect4.split(">")
                            for next_iter in parts_prom:
                                if tester2_string in next_iter  :
                                    print tester2_string
                                    parts_prom_again=next_iter.split(tester2_string)
                                    search.contig=parts_prom_again[1][13:].replace("\n","")
                                    break    
                            file_tester2.close()
                            M=M[1][3:].replace("\n","")
                            test="PASS"
                            search.HNHsequence=M
                        if str(stringIMM[startposIMM+5:endposIMM-5]) in l:
                            N=l.replace("(REVERSE SENSE)","")
                            M=N.split("gc_cont=0.")
                            #print M, str(stringIMM[startposIMM+5:endposIMM-5])
                            M=M[1].replace("\n","")
                            test2="PASS"
                            search.IMMsequence=M
                        if test=="PASS" and test2=="PASS":
                           # y2="rm rm -f 'tester.txt'"
                            #os.system(y2)
                            file_handle5.close()
                            break
                        
                if search.HNHsequence=="":
                    search.HNHsequence=str(string[startposHNH-500:endposHNH+4])
                    search.finder="FORCED"
                    if search.HNHsequence=="":
                         search.HNHsequence=str(string[startposHNH-200:endposHNH+4])
                         search.finder="FORCED"
                if search.IMMsequence=="":
                    search.IMMsequence=str(string[startposIMM-40:endposHNH+50])
                
                y="rm rm -f '_6frame.fa'"
               
                os.system(y)

############################################################################
    
############################################################################

def blast(sequence):
    import os
    blast="./../Programs/Blast/ncbi-blast-2.2.30+/bin/blastp -query "
    blastend=" -db ../Programs/Blast/ncbi-blast-2.2.30+/bin/db/colicin_db -task blastp-fast -outfmt \"7  stitle evalue \" -max_target_seqs 4 -out blast_results.txt"
    command=str(blast+sequence+blastend)
    os.system(command)

############################################################################

############################################################################
def sequence_finder5(object_array):
#will return translated sequences of proteins found in the start variable from the files stored on harddrive. If files are not on harddrive directory must be changed-found in function holder.
    import sys
    import os
    for search in object_array:

        database=search.database
        print database
        isolate_search=search.identifier
        print isolate_search
        if database=="PATRIC":
            directory="/run/media/csharp/Dell_USB_Portable_HDD/ftp.patricbrc.org/patric2/genomes/genomes"
            
        else:
            directory= "/run/media/csharp/Dell_USB_Portable_HDD/All_contigs"
            
        transeq="./../Emboss/EMBOSS-6.6.0/emboss/transeq -sequence" # runs transeq on the command line
        transeqend=" -outseq $p'_6frame.fa' -clean -table 11 -frame 6"   
        
        getorf="./../Emboss/EMBOSS-6.6.0/emboss/getorf -sequence"
        getorfend=" -outseq 'tester.txt' -table 11 -find 1 " #produces START-STOP nucleotide translations(could change to zero and produce protein translations??)
        getorfend2=" -outseq 'tester2.txt' -table 11 -find 4 -flanking 400" # finds regions flanking start codons
        
        
        endposHNH=search.HNHstop
    
        endposIMM=search.IMMstop
        startposHNH=search.HNHstart
        startposIMM=search.IMMstart
        
        
        HNHcontig=search.HNHcontig
        
        IMMcontig=search.IMMcontig
        if database==("PATRIC"):
           
           parts=HNHcontig.split("/")
           filename=str(parts[9])
           filename2=str(parts[9]+"/"+parts[9]+".fna")
           
        else:
            parts=HNHcontig.split("_")
    
            filename2=str(parts[0]+"_"+parts[1]+"_"+parts[2]+"_contigs.fa")
            filename=str(parts[0]+"_"+parts[1]+"_"+parts[2]+"_contigs.fa")
        
        collect=""
        frameHNH=int(HNHcontig[-1:])
        frameIMM=int(IMMcontig[-1:])
        
        
        print filename
        for file in os.listdir(str(directory)):
            
            if file== filename:
                
                file_handle=open(str(directory+"/" + filename2),"r")
               
                for item in file_handle:
                    collect=collect+item
                   
                parts=collect.split(">") 
                for i in parts:
                       
                    if database != "PATRIC":
                        
                        if i[:len(HNHcontig[:-2])]==HNHcontig[:-2]:
                           
                            string_safe=str(i[len(HNHcontig):])
                            
                           
                            break
                    else:
                        
                        if isolate_search in i[:50] and database == "PATRIC":
                            print i[:50], isolate_search
                            j=i.split("]\n")
                            
                            string_safe=str(j[1])
                file_handle.close()
                file_handle8=open("seqtest.fa","w")
                file_handle8.write(str(">Seqtester\n"))
                file_handle8.write(string_safe)
                file_handle8.close()
                
                command=str(transeq+" "+directory+"/" + filename2+transeqend)
                command_getorf=str(getorf+" "+"seqtest.fa"+getorfend)
                command_getorf2=str(getorf+" "+"seqtest.fa"+getorfend2)                
                
                os.system(command_getorf)
                os.system(command)
                os.system(command_getorf2)

                file_handle=open(str("_6frame.fa"),"r")
                
                for item in file_handle:
                    collect=collect+item
                    
                parts=collect.split(">") 
                
                for i in parts[1:]:
                    
                        
                    if database != "PATRIC":
                        if i[:len(HNHcontig)]==HNHcontig:
                            #print i[-1:len(contig[:-2]):-1]
                            string=str(i[len(HNHcontig):]).replace("\n","")
                            string=string.replace(">","")
                            string=string.replace(" ","")
                            
                        if i[:len(IMMcontig)]==IMMcontig:
                            #print i[-1:len(contig[:-2]):-1]
                            stringIMM=str(i[len(IMMcontig):]).replace("\n","")
                            stringIMM=stringIMM.replace(">","")
                            stringIMM=stringIMM.replace(" ","")
                            
                    else:
                        
                        test="fail"
                        test2="fail"
                        if str(i[:len(isolate_search)+2])==str(isolate_search+"_"+str(frameHNH)) and database == "PATRIC":
                            j=i.split("]\n")
                            
                            string=str(j[1]).replace("\n","")
                            string=string.replace("\t","")
                            test="PASS"
                            
                            
                        if str(i[:len(isolate_search)+2])==str(isolate_search+"_"+str(frameIMM)) and database == "PATRIC":
                            j=i.split("]\n")
                            stringIMM=str(j[1]).replace("\n","")
                            stringIMM=stringIMM.replace("\t","") 
                            test2="PASS"
                        if test=="PASS" and test2=="PASS":
                        
                            break
                        
                file_handle.close()            
                
                
                    
                collect2=""
                file_handle5=open("tester.txt",'r')
                
                for items in file_handle5:
                    collect2=collect2+items
                parts2=collect2.split(">")
                #parts2.sort(key=len, reverse=True)
                test="FAIL"
                test2="FAIL"
                for k in parts2[1:]:
                    if k !="":
                        l=k.replace("\n","")
                        if str(string[startposHNH+10:endposHNH-10]) in l:
                             
                            N=l.replace("(REVERSE SENSE)","")
                            M=N.split("]")
                            positions=M[0].split(" ")
                            start=int(positions[1].replace("[",""))
                            end=int(positions[3])
                            search.contig_start=end
                            search.contig_end=start
                            file_tester2=open("tester2.txt","r")
                            tester2_string=str("Around codon at "+str(search.contig_start)+".")
                            tester2_string2=str("Around codon at "+str(search.contig_end)+".")
                            collect4=""
                            for line in file_tester2:
                               collect4+=line
                               
                            parts_prom=collect4.split(">")
                            for next_iter in parts_prom:
                                if tester2_string in next_iter or tester2_string2 in next_iter:
                                    parts_prom_again=next_iter.split(".")
                                    search.contig=parts_prom_again[1].replace("\n","")
                                    break    
                            file_tester2.close()
                            M=M[1].replace("\n","")
                            test="PASS"
                            search.HNHsequence=M
                        if str(stringIMM[startposIMM+5:endposIMM-5]) in l:
                            N=l.replace("(REVERSE SENSE)","")
                            M=N.split("]")
                            #print M, str(stringIMM[startposIMM+5:endposIMM-5])
                            M=M[1].replace("\n","")
                            test2="PASS"
                            search.IMMsequence=M
                        if test=="PASS" and test2=="PASS":
                           # y2="rm rm -f 'tester.txt'"
                            #os.system(y2)
                            file_handle5.close()
                            break
                        
                if search.HNHsequence=="":
                    search.HNHsequence=str(string[startposHNH-500:endposHNH+4])
                    search.finder="FORCED"
                    if search.HNHsequence=="":
                         search.HNHsequence=str(string[startposHNH-200:endposHNH+4])
                         search.finder="FORCED"
                if search.IMMsequence=="":
                    search.IMMsequence=str(string[startposIMM-40:endposHNH+50])
                
                y="rm rm -f '_6frame.fa'"
               
                os.system(y)

#####################################################################################################

########################################################################################################

def THE_PICKLER(filename, object_name, save_load):
    import cPickle as pickle
    if save_load.lower=="save":
        with open(filename,"wb") as f:
            pickle.dump(object_name,f)
    if save_load.lower=="load":
        file_handle=open(filename,"r")
        object_name=pickle.load(file_handle)
        return object_name
        
#########################################################################################
        
#######################################################################################################

def profiles_finder(Profile, object_name):
    counter=0
    for i in object_name:
        if Profile in i.profiles:
            counter+=1
            print i.identifier, i.Species, i.HNHsequence, len(i.HNHsequence), i.coltype, i.profiles, i.comments, i.finder      
    print "I found this many:" , counter
    
    
    
    
    
    
    
    
allowed_profiles=["HNH", "SHOCT", "Pyocin_S", "Colicin", "Cytotoxic", "Colicin_D", "E2R135", "IncFII_repA", "Colicin_Ia", "DUF1851"]

################################################################################################

##################################################################################################

def blast(sequence, matched_array):
    import os
    blast="./../Programs/Blast/ncbi-blast-2.2.30+/bin/blastp -query "
    blastend=" -db ../Programs/Blast/ncbi-blast-2.2.30+/bin/db/colicin_db -task blastp-fast -outfmt \"7  stitle \" -max_target_seqs 4 -out blast_results.txt"
    command=str(blast+sequence+blastend)
    os.system(command)
    blast(sequence)
    collect="" 
    file_handle=open("blast_results.txt","r")
    for line in file_handle:
        collect=collect+line
    parts=collect.split("# Query:")
    for i in parts:
        smaller_parts=i.split("\n")
        if "# 0 hits found" in smaller_parts:
            continue
            #print smaller_parts[0]
        else:
            for j in smaller_parts[4:7]:
                if "hypothetical" not in j and "predicted" not in j:
                    j=j.replace("[","|")
                    tiny_parts=j.split("|")
                    seqdic[smaller_parts[0].replace(" ","")].append(tiny_parts[4])
                
    for i in matched_array:
        if seqdic[i.HNHsequence.replace(" ","")]!=[]:
            print seqdic[i.HNHsequence.replace(" ","")][0]
            i.comments=seqdic[i.HNHsequence.replace(" ","")][0] 
    
    return matched_array

#################################################################################################

#################################################################################################

def sequence_finder7(object_array):
#will return translated sequences of proteins found in the start variable from the files stored on harddrive. If files are not on harddrive directory must be changed-found in function holder.
    import sys
    import os
    import shutil
    for search in object_array:
        database=search.database
        isolate_search=search.identifier
        print isolate_search
        
        #databases are stored in different directories
        if database=="PATRIC":
            directory="/run/media/csharp/Dell_USB_Portable_HDD/ftp.patricbrc.org/patric2/genomes/genomes"
            
        else:
            directory= "/run/media/csharp/Dell_USB_Portable_HDD/All_contigs"
            
        transeq="./../Emboss/EMBOSS-6.6.0/emboss/transeq -sequence"             # runs transeq on the command line
        transeqend=" -outseq $p'_6frame.fa' -clean -table 11 -frame 6"          #clean for * as stop codons, bacterial table and all 6 reading frames
        
        #Applies Prokka orf finder and annotation software
        getorf="./../Programs/prokka-1.10/bin/prokka "
        getorfend=" --force --outdir ~/windowsstuff/prokk_results --prefix tester --cpus 4 --norrna --notrna --quiet"
        
        
        #finds the start and stop positions of HNH and IMM as given by HMMer
        endposHNH=search.HNHstop
    
        endposIMM=search.IMMstop
        startposHNH=search.HNHstart
        startposIMM=search.IMMstart
        
        #finds contig filename, different naming conventions between bigs and PATRIC
        HNHcontig=search.HNHcontig
        IMMcontig=search.IMMcontig
        if database==("PATRIC"):
           
           parts=HNHcontig.split("/")
           filename=str(parts[9])
           filename2=str(parts[9]+"/"+parts[9]+".fna")
           
        else:
            parts=HNHcontig.split("_")
    
            filename2=str(parts[0]+"_"+parts[1]+"_"+parts[2]+"_contigs.fa")
            filename=str(parts[0]+"_"+parts[1]+"_"+parts[2]+"_contigs.fa")
            
        # frame is last character of contig
        collect=""
        frameHNH=int(HNHcontig[-1:])
        frameIMM=int(IMMcontig[-1:])
        
        print filename , search.HNHcontig
        for file in os.listdir(str(directory)):
            
            if file== filename:
                
                file_handle=open(str(directory+"/" + filename2),"r")
                # read file in as one long string collect
                for item in file_handle:
                    collect=collect+item
                #split into contigs by ">"   
                parts=collect.split(">") 
                for i in parts:
                   
                    # finds contig that shares name with HNHcontig    
                    if database != "PATRIC":
                        if i[:len(HNHcontig[:-2])]==HNHcontig[:-2]:
                           
                            string_safe=str(i[len(HNHcontig):])
                            
                           
                            break
                    else:
                        
                        if isolate_search in i[:50] and database == "PATRIC":
                            print i[:50], isolate_search
                            
                            # take all the data in contig and write to file
                            j=i.split("]\n")
                            string_safe=str(j[1])
                file_handle.close()
                file_handle8=open("seqtest.fa","w")
                #write contig to file
                file_handle8.write(str(">Seqtester\n"))
                file_handle8.write(string_safe)
                file_handle8.close()
                
                # finds the 6 frame translation of the contig    
                command=str(transeq+" "+directory+"/" + filename2+transeqend)
                
                #performs the Prokka pipline to the contig
                command_getorf=str(getorf+" "+"~/windowsstuff/seqtest.fa"+getorfend)
                print command_getorf
                os.system(command)
                string=str("/run/media/csharp/Dell_USB_Portable_HDD/gff_files/"+str(search.identifier+"_"+str(search.HNHstart)+"_"+str(search.IMMstart))+".gbk")
                print os.path.exists(string)
                if os.path.exists(string)==True:
                    print "Already exists"
                    temp_var="True"
                else:
                    print "Does not exit creating GBK file"
                    os.system(command_getorf)
                    temp_var="False"

                file_handle=open(str("_6frame.fa"),"r")
                
                for item in file_handle:
                    collect=collect+item
                    
                parts=collect.split(">") 
                
                for i in parts[1:]:
                    #from the 6 frame translation contig in the right reading frame and puts them as string and styring IMM
                        
                    if database != "PATRIC":
                        if i[:len(HNHcontig)]==HNHcontig:
                            #print i[-1:len(contig[:-2]):-1]
                            string=str(i[len(HNHcontig):]).replace("\n","")
                            string=string.replace(">","")
                            string=string.replace(" ","")
                            
                        if i[:len(IMMcontig)]==IMMcontig:
                            #print i[-1:len(contig[:-2]):-1]
                            stringIMM=str(i[len(IMMcontig):]).replace("\n","")
                            stringIMM=stringIMM.replace(">","")
                            stringIMM=stringIMM.replace(" ","")
                            
                    else:
                        
                        test="fail"
                        test2="fail"
                        if str(i[:len(isolate_search)+2])==str(isolate_search+"_"+str(frameHNH)) and database == "PATRIC":
                            j=i.split("]\n")
                            
                            string=str(j[1]).replace("\n","")
                            string=string.replace("\t","")
                            test="PASS"
                            
                            
                        if str(i[:len(isolate_search)+2])==str(isolate_search+"_"+str(frameIMM)) and database == "PATRIC":
                            j=i.split("]\n")
                            stringIMM=str(j[1]).replace("\n","")
                            stringIMM=stringIMM.replace("\t","") 
                            test2="PASS"
                        if test=="PASS" and test2=="PASS":
                        
                            break
                        
                file_handle.close()            
                
                
                
                    
                if os.path.getsize('prokk_results/tester.faa')==0:
                    pass
                else:
                    dict_of_residues={}
                    dict_of_nuc={}
                    file_handle5=open("prokk_results/tester.faa",'r')
                    for items in file_handle5:
                        if ">" in items:
                            line=items
                            dict_of_residues[line]=""
                            
                        else:
                            dict_of_residues[line]=str(dict_of_residues[line]+items.replace("\n",""))
                    
                    file_handle6=open("prokk_results/tester.faa","r")
                    for thing in file_handle6:
                        if ">" in thing:
                            stu=thing
                            dict_of_nuc[stu]=""
                    else:
                        dict_of_nuc[stu]=str(dict_of_nuc[stu]+thing.replace("\n",""))   
                        
                    for k in dict_of_residues:
                        if k !="":
                            l=dict_of_residues[k]
                            if str(string[startposHNH+5:endposHNH-5]) in l:
                                
                                search.HNHsequence=l
                                search.HNHnucleotide=dict_of_nuc[k]
                                break
                    for next_iter in dict_of_residues:
                        if next_iter !="":
                            l=dict_of_residues[next_iter]
                            if str(string[startposIMM+5:endposIMM-5]) in l:
                                
                                search.IMMsequence=l
                                search.IMMnucleotide=dict_of_nuc[k]
                                break
                
                
                #parts2.sort(key=len, reverse=True)
                            
                            
                           
                               
                # if the previous mathod doesnt work it puts it as forced and just adds the hmmer co-ordinates plus a preset amount            
                        
                if search.HNHsequence=="":
                    search.HNHsequence=str(string[startposHNH-500:endposHNH+4])
                    search.finder="FORCED"
                    if search.HNHsequence=="":
                         search.HNHsequence=str(string[startposHNH-200:endposHNH+4])
                         search.finder="FORCED"
                if search.IMMsequence=="":
                    search.IMMsequence=str(string[startposIMM-40:endposHNH+50])
                
                #removes 6 frame translation and copies gff file with unique identifier to hard drive
                y="rm rm -f '_6frame.fa'"
                string=str("/run/media/csharp/Dell_USB_Portable_HDD/gff_files/"+str(search.identifier+"_"+str(search.HNHstart)+"_"+str(search.IMMstart))+".gff")
                string1=str("/run/media/csharp/Dell_USB_Portable_HDD/gff_files/"+str(search.identifier+"_"+str(search.HNHstart)+"_"+str(search.IMMstart))+".gbk")
                if temp_var=="False":
                    shutil.move('prokk_results/tester.gff', string)
                    shutil.move('prokk_results/tester.gbk', string1)
                
                os.system(y)


################################################################################################################################################################################

#################################################################################################################################################################################

def sequence_finder8(object_array):
#Being used to add nucleotide and protein sequences to existing array
    print "Begining sequence finder..."
    import sys
    import os
    import shutil
    for search in object_array:
        database=search.database
        isolate_search=search.identifier
        print isolate_search
        
        #databases are stored in different directories
        if database=="PATRIC":
            directory="/run/media/csharp/Dell_USB_Portable_HDD/ftp.patricbrc.org/patric2/genomes/genomes"
            
        else:
            directory= "/run/media/csharp/Dell_USB_Portable_HDD/All_contigs_iter"
            
        transeq="./../Emboss/EMBOSS-6.6.0/emboss/transeq -sequence"             # runs transeq on the command line
        transeqend=" -outseq $p'_6frame.fa' -clean -table 11 -frame 6"          #clean for * as stop codons, bacterial table and all 6 reading frames
        
        #Applies Prokka orf finder and annotation software
        getorf="./../Programs/prokka-1.10/bin/prokka "
        getorfend=" --force --outdir ~/windowsstuff/prokka_results --prefix tester --cpus 4 --norrna --notrna --quiet"
        
        
        #finds the start and stop positions of HNH and IMM as given by HMMer
        endposHNH=search.HNHstop
    
        endposIMM=search.IMMstop
        startposHNH=search.HNHstart
        startposIMM=search.IMMstart
        
        #finds contig filename, different naming conventions between bigs and PATRIC
        HNHcontig=search.HNHcontig
        IMMcontig=search.IMMcontig
        if database==("PATRIC"):
           
           parts=HNHcontig.split("/")
           filename=str(parts[9])
           filename2=str(parts[9]+"/"+parts[9]+".fna")
           
        else:
            parts=HNHcontig.split("_")
    
            filename2=str(parts[0]+"_"+parts[1]+"_"+parts[2]+"_contigs.fa")
            filename=str(parts[0]+"_"+parts[1]+"_"+parts[2]+"_contigs.fa")
            
        # frame is last character of contig
        collect=""
        frameHNH=int(HNHcontig[-1:])
        frameIMM=int(IMMcontig[-1:])
        
        print "Opening ",filename ,"Finding", search.HNHcontig
        #check if file exists    
        if os.path.exists(str(directory+"/" + filename2))==True:
            
            file_handle=open(str(directory+"/" + filename2),"r")
            # read file in as one long string collect
            for item in file_handle:
                collect=collect+item
            #split into contigs by ">"   
            parts=collect.split(">") 
            for i in parts:
               
                # finds contig that shares name with HNHcontig    
                if database != "PATRIC":
                    if i[:len(HNHcontig[:-2])]==HNHcontig[:-2]:
                       
                        string_safe=str(i[len(HNHcontig):])
                        
                       
                        break
                else:
                    
                    if isolate_search in i[:50] and database == "PATRIC":
                        
                        # take all the data in contig and write to file
                        j=i.split("]\n")
                        string_safe=str(j[1])
            file_handle.close()
            file_handle8=open("seqtest.fa","w")
            #write contig to file
            file_handle8.write(str(">Seqtester\n"))
            file_handle8.write(string_safe)         #string_safe is the name of the contig
            file_handle8.close()
            
            # finds the 6 frame translation of the contig    
            command=str(transeq+" "+directory+"/" + filename2+transeqend)
            
            #performs the Prokka pipline to the contig
            
            
            os.system(command)
            string=str("/run/media/csharp/Dell_USB_Portable_HDD/gff_files/"+str(search.identifier+"_"+str(search.HNHstart)+"_"+str(search.IMMstart))+".gbk")
            if os.path.exists(string)==True:
                command_getorf=str(getorf+" "+"~/windowsstuff/seqtest.fa"+getorfend+" --fast")
                temp_var="True"
            else:
                #command_getorf=str(getorf+" "+"~/windowsstuff/seqtest.fa"+getorfend+" --compliant")
                command_getorf=str(getorf+" "+"~/windowsstuff/seqtest.fa"+getorfend)
                temp_var="False"
                #temp_var="True"
            print temp_var, command_getorf
            os.system(command_getorf)
#                if os.path.exists(string)==True:
#                    print "Already exists"
#                    temp_var="True"
#                else:
#                    print "Does not exit creating GBK file"
#                    os.system(command_getorf)
#                    temp_var="False"

            file_handle=open(str("_6frame.fa"),"r")
            
            for item in file_handle:
                collect=collect+item
                
            parts=collect.split(">") 
            
            for i in parts[1:]:
                #from the 6 frame translation contig in the right reading frame and puts them as string and styring IMM
                    
                if database != "PATRIC":
                    if i[:len(HNHcontig)]==HNHcontig:
                        #print i[-1:len(contig[:-2]):-1]
                        string=str(i[len(HNHcontig):]).replace("\n","")
                        string=string.replace(">","")
                        string=string.replace(" ","")
                        
                    if i[:len(IMMcontig)]==IMMcontig:
                        #print i[-1:len(contig[:-2]):-1]
                        stringIMM=str(i[len(IMMcontig):]).replace("\n","")
                        stringIMM=stringIMM.replace(">","")
                        stringIMM=stringIMM.replace(" ","")
                        
                else:
                    
                    test="fail"
                    test2="fail"
                    if str(i[:len(isolate_search)+2])==str(isolate_search+"_"+str(frameHNH)) and database == "PATRIC":
                        j=i.split("]\n")
                        
                        string=str(j[1]).replace("\n","")
                        string=string.replace("\t","")
                        test="PASS"
                        
                        
                    if str(i[:len(isolate_search)+2])==str(isolate_search+"_"+str(frameIMM)) and database == "PATRIC":
                        j=i.split("]\n")
                        stringIMM=str(j[1]).replace("\n","")
                        stringIMM=stringIMM.replace("\t","") 
                        test2="PASS"
                    if test=="PASS" and test2=="PASS":
                    
                        break
                    
            file_handle.close()            
            
            
            
                
            if os.path.getsize('prokka_results/tester.faa')==0:
                pass
            else:
                dict_of_residues={}
                dict_of_nuc={}
                file_handle5=open("prokka_results/tester.faa",'r')
                for items in file_handle5:
                    if ">" in items:
                        line=items
                        dict_of_residues[line]=""
                        
                    else:
                        dict_of_residues[line]=str(dict_of_residues[line]+items.replace("\n",""))
                
                file_handle6=open("prokka_results/tester.ffn","r")
                for thing in file_handle6:
                    if ">" in thing:
                        stu=thing
                        dict_of_nuc[stu]=""
                    else:
                        dict_of_nuc[stu]=str(dict_of_nuc[stu]+thing.replace("\n",""))      
                for k in dict_of_residues:
                    if k !="":
                        l=dict_of_residues[k]
                        if str(string[startposHNH+5:endposHNH-5]) in l:
                            print l
                            search.HNHsequence=l
                            search.HNHnucleotide=dict_of_nuc[k]
                            
                            break
                for next_iter in dict_of_residues:
                    if next_iter !="":
                        l=dict_of_residues[next_iter]
                        if str(string[startposIMM+5:endposIMM-5]) in l:
                            
                            search.IMMsequence=l
                            search.IMMnucleotide=dict_of_nuc[next_iter]
                            break
                file_handle9=open("prokka_results/tester.fsa","r")
                search.contig_nucleotides=""
                counter_of_file=0
                for line in file_handle9:
                    if counter_of_file==0:
                        counter_of_file+=1
                    else:
                        search.contig_nucleotides=str(search.contig_nucleotides+line.replace("\n",""))
                try:  
                    HNHpromoter_pos=search.contig_nucleotides.lower().rfind(search.HNHnucleotide[3:100].replace("\n","").lower())
                    IMMpromoter=search.contig_nucleotides.lower().rfind(search.HNHnucleotide[-100:-50].replace("\n","").lower())
                    search.HNHpromoter = search.contig_nucleotides[HNHpromoter_pos-500:HNHpromoter_pos].upper()
                    search.IMMpromoter = search.contig_nucleotides[IMMpromoter:IMMpromoter+200].upper()
                    pass
                except:
                    print "Could not find nucleotide sequence"
                search.contig_nucleotides=""
            #parts2.sort(key=len, reverse=True)
                        
                        
                       
                           
            # if the previous mathod doesnt work it puts it as forced and just adds the hmmer co-ordinates plus a preset amount            
                    
            if search.HNHsequence=="":
                search.HNHsequence=str(string[startposHNH-500:endposHNH+4])
                search.finder="FORCED"
                if search.HNHsequence=="":
                     search.HNHsequence=str(string[startposHNH-200:endposHNH+4])
                     search.finder="FORCED"
            if search.IMMsequence=="":
                search.IMMsequence=str(string[startposIMM-40:endposHNH+50])
            
            #removes 6 frame translation and copies gff file with unique identifier to hard drive
            y="rm rm -f '_6frame.fa'"
            string=str("/run/media/csharp/Dell_USB_Portable_HDD/gff_files/"+str(search.identifier+"_"+str(search.HNHstart)+"_"+str(search.IMMstart))+".gff")
            string1=str("/run/media/csharp/Dell_USB_Portable_HDD/gff_files/"+str(search.identifier+"_"+str(search.HNHstart)+"_"+str(search.IMMstart))+".gbk")
            if temp_var=="False":
                shutil.move('prokka_results/tester.gbk', string1)
            
            os.system(y)

##############################################################################################

################################################################################################

def gbk_maker(object_array):
#will return translated sequences of proteins found in the start variable from the files stored on harddrive. If files are not on harddrive directory must be changed-found in function holder.
    import sys
    import os
    import shutil
    for search in object_array:
        database=search.database
        isolate_search=search.identifier
        print isolate_search
        
        #databases are stored in different directories
        if database=="PATRIC":
            directory="/run/media/csharp/Dell_USB_Portable_HDD/ftp.patricbrc.org/patric2/genomes/genomes"
            
        else:
            directory= "/run/media/csharp/Dell_USB_Portable_HDD/contigs_ISOLATES_connor"
            
        transeq="./../Emboss/EMBOSS-6.6.0/emboss/transeq -sequence"             # runs transeq on the command line
        transeqend=" -outseq $p'_6frame.fa' -clean -table 11 -frame 6"          #clean for * as stop codons, bacterial table and all 6 reading frames
        
        #Applies Prokka orf finder and annotation software
        getorf="./../Programs/prokka-1.10/bin/prokka "
        getorfend=" --force --outdir ~/windowsstuff/prokk_results --prefix tester --cpus 3 --norrna --notrna --quiet"
        
        
        #finds the start and stop positions of HNH and IMM as given by HMMer
        endposHNH=search.HNHstop
    
        endposIMM=search.IMMstop
        startposHNH=search.HNHstart
        startposIMM=search.IMMstart
        
        #finds contig filename, different naming conventions between bigs and PATRIC
        HNHcontig=search.HNHcontig
        IMMcontig=search.IMMcontig
        if database==("PATRIC"):
           
           parts=HNHcontig.split("/")
           filename=str(parts[9])
           filename2=str(parts[9]+"/"+parts[9]+".fna")
           
        else:
            parts=HNHcontig.split("_")
    
            filename2=str(parts[0]+"_"+parts[1]+"_"+parts[2]+"_contigs.fa")
            filename=str(parts[0]+"_"+parts[1]+"_"+parts[2]+"_contigs.fa")
            
        # frame is last character of contig
        collect=""
        frameHNH=int(HNHcontig[-1:])
        frameIMM=int(IMMcontig[-1:])
        
        print filename , search.HNHcontig
        for file in os.listdir(str(directory)):
            
            if file== filename:
                
                file_handle=open(str(directory+"/" + filename2),"r")
                # read file in as one long string collect
                for item in file_handle:
                    collect=collect+item
                #split into contigs by ">"   
                parts=collect.split(">") 
                for i in parts:
                   
                    # finds contig that shares name with HNHcontig    
                    if database != "PATRIC":
                        if i[:len(HNHcontig[:-2])]==HNHcontig[:-2]:
                           
                            string_safe=str(i[len(HNHcontig):])
                            
                           
                            break
                    else:
                        
                        if isolate_search in i[:50] and database == "PATRIC":
                            print i[:50], isolate_search
                            
                            # take all the data in contig and write to file
                            j=i.split("]\n")
                            string_safe=str(j[1])
                file_handle.close()
                file_handle8=open("seqtest.fa","w")
                #write contig to file
                file_handle8.write(str(">Seqtester\n"))
                file_handle8.write(string_safe)
                file_handle8.close()
                
                # finds the 6 frame translation of the contig    
                
                #performs the Prokka pipline to the contig
                command_getorf=str(getorf+" "+"~/windowsstuff/seqtest.fa"+getorfend)
                print command_getorf
                string=str("/run/media/csharp/Dell_USB_Portable_HDD/gff_files/"+str(search.identifier+"_"+str(search.HNHstart)+"_"+str(search.IMMstart))+".gbk")
                print os.path.exists(string)
                if os.path.exists(string)==True:
                    print "Already exists"
                else:
                    print "Does not exit creating GBK file"
                    os.system(command_getorf)
                    shutil.move('prokk_results/tester.gbk', string)
                    
###################################################################################################################

#####################################################################################################################

def get_GBKfiles(object_vec, number):
    # for each object in object vec, searches for the corresponding genbank file, creates i.record then finds the HNH sequence and identifier the genomic environment found either side by distance number
    for i in object_vec:
        i.record=""
        i.env=[]
        i.env_products=[]
    for i in object_vec:
        try:
            string=str("/run/media/csharp/Dell_USB_Portable_HDD/gff_files/"+str(i.identifier+"_"+str(i.HNHstart)+"_"+str(i.IMMstart))+".gbk")
            i.record = SeqIO.read(string, "genbank")
            print "found"
            counter=0
            feature_counted=[]
            wanted_list=[]
            features=[x for x in i.record.features if x.type == "CDS" or x.type=="gene"]
            test="FAIL"
            for feature in features:
                counter+=1
                    
                if i.HNHsequence[int(len(i.HNHsequence)/3):-int(len(i.HNHsequence)/3)] in feature.qualifiers["translation"][0]:
                    #print feature.qualifiers[feature.qualifiers.keys()[-1]],feature.qualifiers["translation"][0], i.identifier, counter
                    test="PASS"
                    feature_counted.append(counter)
            feature_counted=set(feature_counted)
            for j in feature_counted:
                for k in range(-int(number),int(number)):
                    try:
                        wanted_list.append(features[j+k])
                        i.env=wanted_list
                        i.env_products.append(features[j+k].qualifiers["product"][0])
                    except:
                        continue
            i.record=""
        except:
            print "Didnt find"
            pass
    return object_vec

##################################################################################################################

####################################################################################################################

def analyse_dist3(object_array):
    #tests the second distribution found in DNases
    from collections import defaultdict
    object_dict=defaultdict(list)
    for i in object_array:
        object_dict[i.identifier].append(i)
        
    counter=0
    for K in object_dict.keys():
        counter+=1
        for i in object_dict[K]:
            for j in object_dict[K]:
                if i.coltype==j.coltype and i.database==j.database:
                    if  i.coltype=="D" or i.coltype=="Cyto" or i.coltype=="DNase":
                        if i.MATCH=="IMM" and j.MATCH=="HNH":
                            if i.IMMstart-j.HNHstop <12 and i.IMMstart-j.HNHstop>0:
                                i.HNHstart=j.HNHstart
                                i.HNHstop=j.HNHstop
                                i.HNHcontig=j.HNHcontig
                                i.MATCH="TRUE"
                   
    counter=0  
    matched_array=[]              
    for i in object_dict.keys():
        for j in object_dict[i]:
         if j.MATCH=="TRUE":
             matched_array.append(j)
    return matched_array
###############################################################################    
def plasmid_finder(object_array):
#uses the same contig finder from previous seq_finders passes it to cBar for plasmid prediction
    import subprocess
    for search in object_array:
        database=search.database
        isolate_search=search.identifier
        print isolate_search
        
        #databases are stored in different directories
        if database=="PATRIC":
            directory="/run/media/csharp/Dell_USB_Portable_HDD/ftp.patricbrc.org/patric2/genomes/genomes"
            
        else:
            directory= "/run/media/csharp/Dell_USB_Portable_HDD/contigs_ISOLATES_connor"
            
        plasmid_find="./../Programs/cBar.1.2/cBar.pl seqtest.fa seqtest.txt"
        
        
        #finds the start and stop positions of HNH and IMM as given by HMMer

        
        #finds contig filename, different naming conventions between bigs and PATRIC
        HNHcontig=search.HNHcontig
        if database==("PATRIC"):
           
           parts=HNHcontig.split("/")
           filename=str(parts[9])
           filename2=str(parts[9]+"/"+parts[9]+".fna")
           
        else:
            parts=HNHcontig.split("_")
    
            filename2=str(parts[0]+"_"+parts[1]+"_"+parts[2]+"_contigs.fa")
            filename=str(parts[0]+"_"+parts[1]+"_"+parts[2]+"_contigs.fa")
            
        # frame is last character of contig
        collect=""
        
        print filename , search.HNHcontig
        for file in os.listdir(str(directory)):
            
            if file== filename:
                
                file_handle=open(str(directory+"/" + filename2),"r")
                # read file in as one long string collect
                for item in file_handle:
                    collect=collect+item
                #split into contigs by ">"   
                parts=collect.split(">") 
                for i in parts:
                   
                    # finds contig that shares name with HNHcontig    
                    if database != "PATRIC":
                        if i[:len(HNHcontig[:-2])]==HNHcontig[:-2]:
                           
                            string_safe=str(i[len(HNHcontig):])
                            
                           
                            break
                    else:
                        
                        if isolate_search in i[:50] and database == "PATRIC":
                            print i[:50], isolate_search
                            
                            # take all the data in contig and write to file
                            j=i.split("]\n")
                            string_safe=str(j[1])
                file_handle.close()
                file_handle8=open("seqtest.fa","w")
                #write contig to file
                file_handle8.write(str(">Seqtester\n"))
                file_handle8.write(string_safe)
                file_handle8.close()
        #os.system(plasmid_find)
        result = subprocess.check_output(plasmid_find, shell=True)
        if result[-10]=="1":
            search.plasmid="Plasmid"
        else:
            search.plasmid="Chromosome"
            
###############################################################################
            
###############################################################################
            
def make_addprofiles(filename, OBJ_array):
    file_handle=open(filename,"w")
    for m in OBJ_array:
     if m.HNHsequence!="":
         string=str(">"+m.identifier+"_"+str(m.HNHstart)+"\n"+m.HNHsequence+"\n")
         file_handle.write(string)
    file_handle.close()


##################################################################################

#################################################################################

def addBIGSspecies(object_array):
    # finds the species information for a BIGS hit from a csv file and adds it to i.species
    Isolate2species=defaultdict(list)
    Species2isolate=defaultdict(list)
    file_handle2=open("IsolatesSpecies2.csv","r")
    for line in file_handle2:
        things=line.split(",")
        #print things[1]
        Isolate2species[things[0]].append(things[1].replace(" \n",""))
        Species2isolate[things[0]].append(things[2].replace(" \n",""))
    file_handle2.close()
    for i in object_array:
        if i.database=="BIGS":
            number=i.identifier.replace("ISOLATE_","")
            if Isolate2species[number] != [] :
                i.Species=(Species2isolate[number][0]).replace(" ","_")
                string=(Isolate2species[number][0]).replace(" ","")
                i.isolate=str((Species2isolate[number][0]).replace(" ","_")+"_" +string.replace("\n",""))
    return object_array          
###############################################################################
                
###############################################################################

def clean_array(object_array):
    from collections import defaultdict
    matched_dictionary=defaultdict(list)
    for i in object_array:
        string=str(i.identifier+"_"+str(i.HNHstart)+"_"+str(i.IMMstart)+i.HNHcontig)
        matched_dictionary[string].append(i)
    
        
    clean_array=[]
    for i in matched_dictionary:
        clean_array.append(matched_dictionary[i][0])
    object_array=clean_array
    return object_array

###############################################################################

###############################################################################
def choose_profiles(matched_array):
    #given an arry of hits this function will set up a GUI to select the allowed profiles and return them as strings in a vector
    old_profiles=[]
    for i in matched_array:
        for j in i.profiles:
            if j not in old_profiles:
                old_profiles.append(j)
    global trofile_list
    trofile_list=[]
    #old_profiles=["DUF12","DUF3456","DUF555"]
    
    label_dict={}
    import Tkinter as tk
    master=tk.Tk()
    
    def ammend_list():
        for i in label_dict.keys():
            if label_dict[i].get()==1:
               trofile_list.append(i) 
        
                
    tk.Label(master, text="Profile selection!").grid(row=0, sticky=tk.W)
    for i in range(len(old_profiles)):
        label_dict[old_profiles[i]]=tk.IntVar()
        tk.Checkbutton(master, text=old_profiles[i], variable=label_dict[old_profiles[i]]).grid(row=i+1, sticky=tk.W)
    #def ammend_list():
    #    for i in old_profiles:
    tk.Button(master, text="Select", command=ammend_list).grid(row=i+2, stick=tk.W)
    tk.Button(master, text='Quit', command=master.quit).grid(row=i+3, sticky=tk.W, pady=4)
    
    master.mainloop()
    profile_list=trofile_list
    return profile_list

###############################################################################

###############################################################################


        
        
        
    
    












