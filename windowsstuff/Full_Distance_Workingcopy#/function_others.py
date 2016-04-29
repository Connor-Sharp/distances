# -*- coding: utf-8 -*-
"""
Created on Mon Jan 19 12:44:36 2015

@author: csharp
"""

# -*- coding: utf-8 -*-
"""
Created on Thu May 01 13:30:26 2014

@author: kell3284
"""


from collections import defaultdict


import os


def analyse(object_array,filename, database, first, second):
    from classex1 import ISOLATE
    
    print """\r

===============================================================================

            Analysing """,filename,""" for protein co-ordinates...
            
===============================================================================

\r"""
    data= open(filename, 'r')
    
    
    e_values=[]  # list of e_values 
    counter=0 # counter to miss first line
     #everything in table matched or un matched with e_value over 0.000005
    
    Species={}#dictionary of Species
    for item in data:
        counter+=1
        parts=item.split(",")
        

        if counter>2 and float(parts[7])<float(0.0000005):
            y=ISOLATE()
            y.database=database

            e_values.append(float(parts[7]))
            if database=="PATRIC":
                shortened_isolate=parts[4].split("/")

                y.identifier=parts[3]
                y.isolate=parts[4]
		y.Species=shortened_isolate[:2]

            else:
                Species[parts[3]]=parts[4]
                y.isolate=parts[3]
                y.identifier=parts[3]
                y.Species=(parts[4])
            if first == parts[0]:
                y.HNHstart=int(parts[18])
                y.HNHstop=int(parts[19])
                y.HNHcontig=(parts[5])
                y.coltype="DNase"
                y.contiglength=int(parts[6])
                y.MATCH="HNH"
            
            
            elif second == parts[0]:
                y.IMMstart=int(parts[18])
                y.IMMstop=int(parts[19])
                y.IMMcontig=(parts[5])
                y.contiglength=int(parts[6])
                y.coltype="DNase"
                y.MATCH="IMM"
            
            y.speciesWriter()
            object_array.append(y)
    
    data.close()
    print """\r

===============================================================================

                Finished analysing """,filename,"""
            
===============================================================================

"""
    return object_array


def analyse_dist2(object_array, starts, stops, strand):
    from collections import defaultdict
    import sys
    object_dict=defaultdict(list)
    for i in object_array:
        object_dict[i.identifier].append(i)
    left="""
===============================================================================

                        Analysing distances...
            """
    right="""
    
===============================================================================\r"""
    counter=0
    print left
    maxi=len(object_dict)*len(starts)
    for region in range(len(starts)):
        start=int(starts[region])
        stop=int(stops[region])
        for K in object_dict.keys():
            counter+=1
            progress="====="*int((float(counter)/maxi)*10)
            prog2="     "*(10-int((float(counter)/maxi)*10))
            string="\tProgress ["+progress+prog2+"]"
            sys.stdout.write(string+"\r")
            sys.stdout.flush()
            for i in object_dict[K]:
                for j in object_dict[K]:
                    if i.coltype==j.coltype and i.database==j.database:
                        if strand == "POS":
                            if i.MATCH=="IMM" and j.MATCH=="HNH":
                                if i.IMMstart-j.HNHstop <stop and i.IMMstart-j.HNHstop>start:
                                    i.HNHstart=j.HNHstart
                                    i.HNHstop=j.HNHstop
                                    i.HNHcontig=j.HNHcontig
                                    i.MATCH="TRUE"
                        elif strand == "NEG":
                            if i.MATCH=="IMM" and j.MATCH=="HNH":
                                if i.contiglength-i.IMMstop-j.HNHstop <stop and i.IMMstart-j.HNHstop>start:
                                
                                    i.HNHstart=j.HNHstart
                                    i.HNHstop=j.HNHstop
                                    i.HNHcontig=j.HNHcontig
                                    i.MATCH="TRUE"
    counter=0  
    matched_array2=[]              
    for i in object_dict.keys():
        for j in object_dict[i]:
         if j.MATCH=="TRUE":
             j.region_no=region
             matched_array2.append(j)
    print """


                        Distances analysed
                    
==============================================================================="""
    return matched_array2

################################################################################################

################################################################################################


def plot_distance(object_array, coltype):
    import matplotlib.pyplot as pl
    #import seaborn as sns
    from collections import defaultdict
    object_dict=defaultdict(list)
    for i in object_array:
        object_dict[i.identifier].append(i)
    keys=object_dict.keys()
    
    #finds the isolates containing a HNH and an Imm gene next to each other but with 
    #imm after the hnh returns a dict of lists with all matches
    width=int(raw_input("Maximum distance to search:"))
    isoforgraphs=[]
    isographs=[]
    for M in range(0,width,1): 
        counter=0
        counter2=0
        for k in keys:
            for j in object_dict[k]:
                for i in object_dict[k]:
                    if coltype!="Pore":
                        if i.coltype==j.coltype and i.database==j.database and i.coltype==coltype:
                            if i.MATCH=="IMM" and j.MATCH=="HNH":
                                if i.IMMstart-j.HNHstop <M and i.IMMstart-j.HNHstop>0:
                                    counter+=1
                                if i.IMMstart-j.HNHstop == M:
                                    counter2+=1
                    elif coltype=="Pore":
                        if i.MATCH=="IMM" and j.MATCH=="HNH":
                            if i.contiglength- i.IMMstop- j.HNHstop < M and i.IMMstart-j.HNHstop>0:
                            
                                counter+=1
                            if i.contiglength- i.IMMstop -j.HNHstop == M: 
                                counter2+=1
                        
        isographs.append(counter)
        isoforgraphs.append(counter2)
        

    fig, ax1 = pl.subplots()
    font = {'family' : 'normal','weight' : 'bold','size'   : 22}
    pl.rc('font', **font)
    ax1.plot(range(0,width,1),isoforgraphs,'b',range(0,width,1),isographs,'r')
    ax1.plot(isoforgraphs,linewidth=3)
    ax1.plot(isographs,linewidth=3)
    pl.xlabel('Intergenic region (residues)')
    pl.ylabel('Hits')
    labels=["Distribution", "Cumulative"]
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
            
        transeq="transeq -sequence" # runs transeq on the command line
        transeqend=" -outseq $p'_6frame.fa' -clean -table 11 -frame 6"   
        
        getorf="prodigal -i"
        getorfend=" -a tester.faa -d tester2.fa -f gbk -o output_file_waste.txt -s output_prod_genes.txt -p meta"
        
        
        
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


############################################################################

############################################################################
def sequence_finder5(object_array):
#will return translated sequences of proteins found in the start variable from the files stored on harddrive. If files are not on harddrive directory must be changed-found in function holder.
    import sys
    import os

    left="""\r
    
===============================================================================

        Finding sequences from hits, this may take a while...
        
    
    Progress ["""
    
    right ="""
    
===============================================================================
"""



    

    
    counter=0
    for search in object_array:
        counter+=1
        progress="=" * int((float(counter)/float(len(object_array)))*float(50))
        prog2=" "*(50-int((float(counter)/float(len(object_array)))*float(50)))
        string="\tProgress ["+progress+prog2+"]"
        #sys.stdout.write(string+"\r")
        
        database=search.database
        isolate_search=search.identifier
        sys.stdout.write( left+progress+prog2+"]"+"\n"+isolate_search)

        transeq="transeq -sequence" # runs transeq on the command line
        transeqend=" -outseq $p'_6frame.fa' -clean -table 11 -frame 6"   
        
        getorf="getorf -sequence"
        getorfend=" -outseq 'tester.txt' -table 11 -find 1" #produces START-STOP nucleotide translations(could change to zero and produce protein translations??)
        getorfend2=" -outseq 'tester2.txt' -table 11 -find 4 -flanking 400" # finds regions flanking start codons
        getorfend3=" -outseq 'tester3.txt' -table 11 -find 2"
        
        
        endposHNH=search.HNHstop
    
        endposIMM=search.IMMstop
        startposHNH=search.HNHstart
        startposIMM=search.IMMstart
        
        
        HNHcontig=search.HNHcontig
        IMMcontig=search.IMMcontig

        collect=""
        frameHNH=int(HNHcontig[-1:])
        frameIMM=int(IMMcontig[-1:])
        
        if database !="BIGS":
            genome_file=str(HNHcontig[:-2]+".fna")
            file_handle=open(str(HNHcontig[:-2]+".fna"),"r")
	    print 'Opening', HNHcontig
        else:
            parts=HNHcontig.split("_")
            genome_file=str("/run/media/csharp/Dell_USB_Portable_HDD/All_contigs_iter/"+parts[0]+"_"+parts[1]+"_"+parts[2]+"_contigs.fa")
            if os.path.exists(genome_file)==False:
                print "File doesnt exit!!"
                search.finder="NO_FILE"
                continue
            file_handle=open(genome_file,"r")

        if database == "NCBI":
            string_safe=""
            counter=0
            for line in file_handle:
                if line[0]==">":
			print line
		else:
                    string_safe=string_safe+line
                    
                
            
                    
			
        for item in file_handle:
            collect=collect+item
           
        parts=collect.split(">") 
        for i in parts:
  
            if database == "BIGS":
                if i[:len(HNHcontig[:-2])]==HNHcontig[:-2]:
                   
                    string_safe=str(i[len(HNHcontig):])
                    
                   
                    break
            elif database == "PATRIC":
                
                if isolate_search in i[:50]:
                    print i[:50], isolate_search
		    if "]\n" in i:
                    	j=i.split("]\n")
		    else:
			j=i.split(str(isolate_search+"\n"))
                    
                    string_safe=str(j[1])

			

        file_handle.close()
        file_handle8=open("seqtest.fa","w")
        file_handle8.write(str(">Seqtester\n"))
        file_handle8.write(string_safe)
        file_handle8.close()
        
        command=str(transeq+" "+genome_file+transeqend)
        command_getorf=str(getorf+" "+"seqtest.fa"+getorfend)
        command_getorf2=str(getorf+" "+"seqtest.fa"+getorfend2)                
	command_getorf2=str(getorf+" "+"seqtest.fa"+getorfend3)         

        os.system(command_getorf)
        os.system(command)
        os.system(command_getorf2)

        file_handle=open(str("_6frame.fa"),"r")
        
        for item in file_handle:
            collect=collect+item
            
        parts=collect.split(">") 
        
        for i in parts[1:]:
            
                
            if database != "PATRIC" and database != "NCBI":
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
                    
            elif database == "PATRIC":
                
                test="fail"
                test2="fail"
                if str(i[:len(isolate_search)+2])==str(isolate_search+"_"+str(frameHNH)) and database != "BIGS":
                    j=i.split(str(str("]")+"\n"))
                    
                    string=str(j[1]).replace("\n","")
                    string=string.replace("\t","")
                    test="PASS"
                    
                    
                if str(i[:len(isolate_search)+2])==str(isolate_search+"_"+str(frameIMM)) and database != "BIGS":
                    j=i.split(str(str("]")+"\n"))
                    stringIMM=str(j[1]).replace("\n","")
                    stringIMM=stringIMM.replace("\t","") 
                    test2="PASS"
                if test=="PASS" and test2=="PASS":
                
                    break
            elif database=="NCBI":
                if i[:len(search.identifier)+2]==str(search.identifier+"_"+str(frameHNH)):
                    #print i[-1:len(contig[:-2]):-1]
                        string=str(i[len(search.identifier):]).replace("\n","")
                        string=string.replace(">","")
                        string=string.replace(" ","")
                    
                if i[:len(search.identifier)+2]==str(search.identifier+"_"+str(frameIMM)):
                    #print i[-1:len(contig[:-2]):-1]
                    stringIMM=str(i[len(search.identifier):]).replace("\n","")
                    stringIMM=stringIMM.replace(">","")
                    stringIMM=stringIMM.replace(" ","")
		
                
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
                    print M
                    positions=M[0].split(" ")
                    start=int(positions[1].replace("[",""))
                    end=int(positions[3])
                    search.contig_start=start
                    search.contig_end=end
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

                    file_tester3=open("tester3.txt","r")
                    nucleotide_end=str("- "+str(search.contig_end)+"]")
                    nuc_sequences=""
                    for line in file_tester3:
                       nuc_sequences+=line
                       if nucleotide_end in line:
                           print line
                       
                    nucletide_seqs=nuc_sequences.split(">")
                    for nuc_fasta in nucletide_seqs:
                        if nucleotide_end in nuc_fasta:
                            nuc_fasta=nuc_fasta.replace("(REVERSE SENSE)","")
                            nuc_fasta_parts=nuc_fasta.split("]")
                            search.HNHnucleotide=nuc_fasta_parts[1].replace("\n","")
                            break    
                    file_tester3.close()

                    M=M[1].replace("\n","")
                    test="PASS"
                    search.HNHsequence=M
                if str(stringIMM[startposIMM+5:endposIMM-5]) in l:
                    N=l.replace("(REVERSE SENSE)","")
                    M=N.split("]")
                    print M
                    positions=M[0].split(" ")
                    start=int(positions[1].replace("[",""))
                    end=int(positions[3])
                    search.contig_start=start
                    search.contig_end=end
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

                    file_tester3=open("tester3.txt","r")
                    nucleotide_end=str("- "+str(search.contig_end)+"]")
                    print nucleotide_end
                    nuc_sequences=""
                    for line in file_tester3:
                       nuc_sequences+=line
                       if nucleotide_end in line:
                           print line
                       
                    nucletide_seqs=nuc_sequences.split(">")
                    for nuc_fasta in nucletide_seqs:
                        if nucleotide_end in nuc_fasta:
                            nuc_fasta=nuc_fasta.replace("(REVERSE SENSE)","")
                            nuc_fasta_parts=nuc_fasta.split("]")
                            search.IMMnucleotide=nuc_fasta_parts[1].replace("\n","")
                            break    
                    file_tester3.close()
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
        print right
        os.system(y)

          
