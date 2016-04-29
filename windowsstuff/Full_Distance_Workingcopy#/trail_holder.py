# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 18:00:54 2015

@author: csharp
"""

def analyse(object_array,filename, database, options):
    from classex1 import ISOLATE
    options_dict=options
    print """\r

===============================================================================

            Analysing """,filename,""" for protein co-ordinates... E_value 5E-5
            
===============================================================================

"""
    data= open(filename, 'r')
    
    E_counter=0
    good_counter=0
    total_counter=0
    e_values=[]  # list of e_values 
    counter=0 # counter to miss first line
     #everything in table matched or un matched with e_value over 0.000005
    
    Species={}#dictionary of Species
    for item in data:
        counter+=1
        parts=item.split(",")
        total_counter+=1
        if counter>2 and float(parts[7])<float(0.01):
            good_counter+=1
            y=ISOLATE()
            y.found=[]
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
#
#            e_values.append(float(parts[7]))
#            if database=="PATRIC":
#                parts[4]=(parts[4]).replace("/run/media/csharp/Dell_USB_Portable_HDD/ftp.patricbrc.org/patric2/genomes/genomes/","/home/csharp/PATRIC/ftp.patricbrc.org/patric2/genomes/")
#                shortened_isolate=parts[4].split("/")
#                short_isolate=shortened_isolate[7]
#                Species[parts[3]]=short_isolate
#                y.identifier=parts[3]
#                y.isolate=short_isolate
            for names in options_dict.keys():
		
                if str(parts[0].replace("\n","")) == options_dict[names].left:
                    y.HNHstart=int(parts[18])
                    y.HNHstop=int(parts[19])
                    y.HNHcontig=(parts[5])
                    y.coltype=names
                    y.found.append(parts[0].replace("\n",""))
                    
                    y.contiglength=int(parts[6])
                    y.MATCH="HNH"
                elif str(parts[0].replace("\n","")) in  options_dict[names].right[0]:
                    y.IMMstart=int(parts[18])
                    y.IMMstop=int(parts[19])
                    y.IMMcontig=(parts[5])
                    y.coltype=names
                    y.contiglength=int(parts[6])
                    y.MATCH="IMM"  
                    y.found.append(parts[0].replace("\n",""))
            y.speciesWriter()
            object_array.append(y)
        else:
            E_counter+=1
    print "E-counter=", E_counter
    print "total_counter=", total_counter
    print "good_counter=", good_counter
    print "difference=", int(total_counter-(E_counter+good_counter))
    data.close()
    print """\r

===============================================================================

            Finished analysing """,filename,"""
            
===============================================================================

"""
    return object_array  



def parser(file_input, file_output):

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
        if "MotA_ExbB" in parts[0] or "TonB_C" in parts[0] or "TolB_N" in parts[0] or "ExbD" in parts[0] or "TolA" in parts[0]:
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
            try:
                part.append(Isolate2species[str(stuff[0]+"_"+stuff[1])])
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
                for i in range(len(part)):
                    if part[i] != "":
                        file_towrite.write(part[i].strip("\n"))
                        file_towrite.write(",")
                file_towrite.write("\n")
            except:
                pass
            print part
    file_handle.close() #Hmmer output file
    file_handleagain.close() #Hmmer output file
    file_towrite.close()# file to write
    file_handle2.close() # file for id to species
###############################################################################
    
###############################################################################

def plot_distance(object_array, coltype, date):
    import matplotlib.pyplot as pl
    import seaborn as sns
    from collections import defaultdict
    object_dict=defaultdict(list)
    for i in object_array:
        object_dict[i.identifier].append(i)
    keys=object_dict.keys()
    
    #finds the isolates containing a HNH and an Imm gene next to each other but with 
    #imm after the hnh returns a dict of lists with all matches
    sns.set_style("white")
    sns.set_context("notebook", font_scale=1.5, rc={"lines.linewidth": 2.5})
    width=int(raw_input("Maximum distance to search:"))
    isoforgraphs=[]
    isographs=[]
    reverse1=[]
    reversegraph=[]
    for M in range(0,width,1): 
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
    font = {'family' : 'normal','weight' : 'bold','size'   : 28}
    pl.rc('font', **font)
    ax1.plot(range(0,3*width,3),isoforgraphs,'black',range(0,3*width,3),isographs,'--',range(0,3*width,3),reversegraph,'-')
#    ax1.plot(isoforgraphs,linewidth=3)
#    ax1.plot(isographs,linewidth=3)
#    ax1.plot(reversegraph, linewidth=3)
    pl.xlabel('Intergenic region (bp)')
    pl.ylabel('Profile pair hits')
    labels=["col-imm", "imm-col", "reverse strand Imm"]
    #ax1.legend(labels)
    pl.show()
    figure_file=str(date+"/distance_figure"+coltype+".png")
    fig.savefig(figure_file, dpi=fig.dpi)
    return isographs, isoforgraphs 
    
    
#==============================================================================

#==============================================================================

def analyse_dist(object_array,name, starts, stops, strand):
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
    old_length=0
    new_length=0
    Results_in_region=[]
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
                    if i.coltype==j.coltype and i.database==j.database and i.coltype==name:
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
        new_length=len(matched_array2)
        Results_in_region.append(new_length-old_length)
        old_length=new_length
    print """


                        Distances analysed
                    
==============================================================================="""
    return Results_in_region, matched_array2
    
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
        getorfend=" -outseq 'tester.txt' -table 11 -find 1 " #produces START-STOP nucleotide translations(could change to zero and produce protein translations??)
        getorfend2=" -outseq 'tester2.txt' -table 11 -find 4 -flanking 400" # finds regions flanking start codons
        
        
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
                    j=i.split("]\n")
                    
                    string_safe=str(j[1])

			

        file_handle.close()
        file_handle8=open("seqtest.fa","w")
        file_handle8.write(str(">Seqtester\n"))
        file_handle8.write(string_safe)
        file_handle8.close()
        
        command=str(transeq+" "+genome_file+transeqend)
        command_getorf=str(getorf+" "+"seqtest.fa"+getorfend)
        command_getorf2=str(getorf+" "+"seqtest.fa"+getorfend2)                
        
        os.system(command_getorf)
        os.system(command)
        os.system(command_getorf2)

        file_handle=open(str("_6frame.fa"),"r")
        for item in file_handle:
            counter=0
            string=""
            stringIMM=""
            if search.identifier in item:
                counter+=1
            if counter>0:
                string=string+item.replace("\n","")
                stringIMM=stringIMM+item.replace("\n","")
                
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
                    j=i.split("]\n")
                    
                    string=str(j[1]).replace("\n","")
                    string=string.replace("\t","")
                    test="PASS"
                    
                    
                if str(i[:len(isolate_search)+2])==str(isolate_search+"_"+str(frameIMM)) and database != "BIGS":
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
        print right
        os.system(y)


###############################################################################
            
###############################################################################
            
            
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

def search_profiles(file_in, file_out):
    import os
    command=str("hmmscan --tblout "+file_out+" /home/csharp/hmmerfiles/Pfam-A.hmm "+file_in)
    os.system(command)
###############################################################################
###############################################################################

    
def addprofiles(filename, OBJ_array):
    matched_array=OBJ_array
    file_handle=open(filename,"r")
    count=0
    for line in file_handle:
        count+=1
        if count>3 and line[0]!="#":
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
