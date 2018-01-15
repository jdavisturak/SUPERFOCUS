#SUPER-FOCUS: A tool for agile functional analysis of shotgun metagenomic data || version 0.22.X.1
#-------------------------------------------------------------------------------------------------------------
#(c)            Silva, G. G. Z., Green K., B. E. Dutilh, and R. A. Edwards: 
#			   SUPER-FOCUS: A tool for agile functional analysis of shotgun metagenomic data. (under review)
#website: 	   https://edwards.sdsu.edu/SUPERFOCUS

import os,sys,numpy as np
import subprocess

options= "SUPER-FOCUS: A tool for agile functional analysis of shotgun metagenomic data\n"\
      "-------------------------------------------------------------------------------------------------------------\n"\
      "Options:\n"\
      "         -h     ------: print help\n"\
      "         -q     string: query file (FASTA or FASTQ format)\n"\
      "         -dir   string: output directory\n"\
      "         -o     string: project name (default 'my_project')\n"\
      "         -mi    float:  minimum identity (default 60 %)\n"\
      "         -ml    int:    minimum alignment (default 15 amino acids)\n"\
      "         -focus int:    runs FOCUS; 1 does run; 0 does not run: default 0 \n"\
      "         -t     int:    number of threads (default 8)\n"\
      "         -e     float:  e-value (default 0.00001)\n"\
      "         -fast  int:    0 (False)/ 1 (True) (default 1)\n"\
      "         -db    string: database (DB_90, DB_95, DB_98, or DB_100; default DB_98)\n"\
      "         -p     int:    amino acid input; 0 nucleotides; 1 amino acids (default 0)\n"\
      "         -a     string: aligner choice (rapsearch, blast, diamond; default rapsearch)\n"\
      "         -n     int:    normalize results based on counts; 0 does not normalize; 1 does normalize (default 1)\n"\
      "         -r     string: reduction table used in the SUPER-FOCUS reduction: rast / ncbi(default)\n"\
      "-------------------------------------------------------------------------------------------------------------\n"\
      "example> python superfocus.py -q query.fasta -dir myOutputdirectory"

#hash with all the default parameters
myproject={'-o':"my_project",'-dir':"","-q":"","-mi":'60',
           "-ml":"15","-n":"1","-focus":"0","-a":"rapsearch",
           "-t":"8","-e":"0.00001","-db":"DB_98","-p":"0","-fast":"1","-r":"ncbi"}

#gets the user parameters and add in the hash
def setParameters():
    if "/" in sys.argv[0]:
        myproject["dir"]="/".join(sys.argv[0].split("/")[:-1])+"/"
    else:
        myproject["dir"]=""
        
    userParameters=sys.argv[1:]

    if "-h" in sys.argv:
        pass
    else:
        for i in range(0,len(userParameters),2):
            try:
                myproject[userParameters[i]]=userParameters[i+1]
            except:
                if userParameters[i] in myproject:
                    print "Please inform a value for "+userParameters[i]
                else:
                    print userParameters[i]+" is not a valid parameter"
setParameters()

############################################# SUPER-FOCUS parameters #############################################
project_name=myproject["-o"]
project_output=myproject["-dir"]+"/"

#project_output=""
if project_output != "":
    #creates the folder to output the files
    os.system("mkdir -p "+project_output)
    dbLocation =  project_output + "/db"
    os.system("mkdir -p " + dbLocation)
else:
    os.system("mkdir -p "+project_name)
    project_output=project_name# a folder with the project name is created in the SF directory
    dbLocation = myproject["dir"] + "/db"
systemDbLocation = myproject["dir"] + "/db"

print "User's output folder is {}. Database location is {}".format(project_output, dbLocation)
query=myproject["-q"]
mi=float(myproject["-mi"])
ml=int(myproject["-ml"])
normalize=int(myproject["-n"])
f=str(myproject["-focus"]) #1 run focus | 0 not run focus

#Aligner parameters
aligner=myproject["-a"].lower()
T=myproject["-t"]#number of threads
evalue=myproject["-e"]
try:mydb=myproject["-db"].split("_")[1]
except: pass
my_alignments=project_name+"__alignments"
fast_mode=myproject["-fast"]#fast[T] or sensitive[F]
proteins=int(myproject["-p"]) #if this parameter is 1, the FOCUS reduction cannot be applied

###################################################################################################################

#returns the path for a given program name
def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

#This function runs FOCUS and returns the list of organisms present in a metagenome
def runNParseFocus():
    #run
    os.system("python "+myproject["dir"]+"focus/focus.py -m "+str(1./(10**10))+" -q "+query)

    #parse focus output
    organisms=[]
    f=open(query+"__FOCUS_output.txt")
    check=0
    for line in f:
        if "Genus Level" in line: check = 1
        elif check == 1: check+=1
        elif check==2:
            line=line.split()
            if len(line)==3:
                organisms.append(line[1].split("_")[0])
            else:
                break
    f.close()

    os.system("mv "+query+"__FOCUS_output.txt "+myproject["-dir"])
    return organisms

#This function gets the list of organisms present in a metagenome and formatted the DB according to the selected aligner
def formatDb(organisms):
    def sumColumn(matrix):
        return np.sum(matrix, axis=0)

    if myproject["-r"]=="rast":
        f=open(systemDbLocation+"/organisms2subsystem__RAST.txt")
    else:
        f=open(systemDbLocation+"/organisms2subsystem.txt")
    f.readline()
    matrix=[]

    for line in f:
        line=line.split()
        #if the organism is in the prediction, we keep the subsystems present on it
        if line[0] in organisms:
            matrix.append([int(x) for x in line[1:]])
    f.close()

    matrix=sumColumn(matrix)

    #return all the subsystems IDs present int the predicted organisms
    subsystems=[str(i+1) for i in range(0,1290) if matrix[i]!=0]

    path=" "+systemDbLocation+"/clusters/"+mydb+"_clusters/"
    print "\n\nFormatting DB: "+str(len(subsystems))+" are going to be used in the DB"
    dbname="DB__"+query.split("/")[-1]

    os.system("cat "+path+path.join([value+"_cluster.faa" for value in subsystems])+" > "+dbLocation+"/focus_reduction/"+dbname)
    
    
    #Here the database is formatted based on your option of aligner
    if aligner=="rapsearch":
        os.system("prerapsearch -d "+dbLocation+"/focus_reduction/"+dbname+" -n "+dbLocation+"/focus_reduction/rapsearch2/"+dbname+".db")
    elif aligner=="diamond":
        os.system("diamond makedb --in  "+dbLocation+"/focus_reduction/"+dbname+" --db "+dbLocation+"/focus_reduction/diamond/"+dbname+".db")
    else:
        os.system("makeblastdb -in "+dbLocation+"/focus_reduction/"+dbname+" -dbtype prot -out "+dbLocation+"/focus_reduction/blast/"+dbname+".db -title "+dbname)

    #delete fasta with sequences
    os.system("rm "+dbLocation+"/focus_reduction/"+dbname)
    
    return dbname+".db"


#This function calls the user aligner choice and aligns it to one of the 4 SUPER-FOCUS database
#INPUT: it the user query
#OUPUT: alignments in tabular format
def align(mydb):
    aligner=myproject["-a"]
    databaseMode="static"
    currentDbLocation = systemDbLocation+"/static"
    if "DB__" in mydb:#it run FOCUS
        databaseMode="focus_reduction"
        currentDbLocation=dbLocation+"/focus_reduction/"
        os.system("mkdir -p "+currentDbLocation+aligner+('2' if aligner == 'rapsearch' else ''))
    else:
        mydb=mydb+".db"
        
    if aligner=="rapsearch":
        os.system("rapsearch -a "+fast_mode+" -q "+query+" -d "+currentDbLocation+"/rapsearch2/"+mydb+" -o "+project_output+"/"+project_name+"__alignments -v 250 -z "+T+" -e "+evalue+" -b 0 -s f")

    elif aligner=="diamond":
        blast="blastx"
        if proteins==1:#we have proteins as input
            blast="blastp"
        dia = str("diamond "+blast+" -d "+currentDbLocation+"/diamond/"+mydb+".dmnd -q "+query+ " -a "+project_output+"/"+project_name+"__alignments.daa"+" -k 250 -p "+T+" -e "+evalue)
        if fast_mode==0:#fast mode
            dia += " --sensitive"
        print(dia)
        os.system(dia)
        print ("running diamond view")
        #dia2 = str("diamond view -a "+project_output+"/"+project_name+"__alignments.daa -o "+project_output+"/"+project_name+"__alignments.m8")
        dia2 = ["diamond",  "view", "-a", project_output+"/"+project_name+"__alignments.daa", "-o", project_output+"/"+project_name+"__alignments.m8"]
        a=subprocess.Popen(dia2)
        a.wait()
        #os.system(dia2)
        print(dia2)
        # else:#sensitive mode
        #     os.system("diamond "+blast+" -d "+dbLocation+"/"+databaseMode+"/diamond/"+mydb+" -q "+query+" -o "+project_output+"/"+project_name+"__alignments.m8 -k 250  -p "+T+" -e "+evalue+" --sensitive")
            
    elif aligner=="blast":
        if proteins==0:#we have nucleotides as input
            os.system("blastx -db "+currentDbLocation+"/blast/"+mydb+" -query "+query+" -out "+project_output+"/"+project_name+"__alignments.m8 -outfmt 6 -evalue "+evalue+" -max_target_seqs 250 -num_threads "+T)
        else:#if the user wants to input proteins, but the FOCUS prediction would be done
            os.system("blastp -db "+currentDbLocation+"/blast/"+mydb+" -query "+query+" -out "+project_output+"/"+project_name+"__alignments.m8 -outfmt 6 -evalue "+evalue+" -max_target_seqs 250 -num_threads "+T)

    #if databaseMode=="focus_reduction":
    #    [os.system("rm "+dbLocation+"/focus_reduction/"+aligner+"/"+mydb+'* 2> /dev/null') for aligner in ["blast","diamond","rapsearch2"]]
        

#This function parses the alignments of BLAST, RAPSEARCH2 or DIAMOND
#OUTPUT: a hash with the counts of hits into subsystems
def parse_alignments():
    subsystems={};functions={}
    #reads the alignments and normalizes it if needed
    def alignments_count(alignments):
        L=1.0;L2=1.0
        for sequenceID in alignments:
            #It gets how many subsystems were associated to a sequence to normalize the abundance later
            if normalize==1:
                L=len(alignments[sequenceID][1])*1.
                L2=len(alignments[sequenceID][2])*1.
            
            for subsystem in alignments[sequenceID][1]:
                try:#normalizes by the number of subsystem assigned to a sequence (if normalize = 1)
                    subsystems[subsystem]+=1/L
                except:#first time adding
                    subsystems[subsystem]=1/L

            for function in alignments[sequenceID][2]:
                try:#normalizes by the number of functions assigned to a sequence (if normalize = 1)
                    functions[function]+=1/L
                except:#first time adding
                    functions[function]=1/L

    alignments={} 
    f=open(project_output+"/"+my_alignments+".m8")

    if aligner=="rapsearch":
        [f.readline() for line in xrange(5)] #reads the file useless lines and avoids checking for "#"

    ID=""
    #parses the tabular output from BLAST, DIAMOND or RAPSearch2    
    for line in f:
        line=line.split("\t")

        #Only selects hits with min identity 60% and alignments of min 15 bp
        if (ID!=line[0] or ID=="") and float(line[2])>=mi and int(line[3])>=ml:            
            E=float(line[10])#E-value
            SS=line[1].split("__")[1]#Subsystem PK on SF database
            myfunction=line[1].split("__")[2]#Function

            #Saves other hits that had the same e-value than the best hit
            try:
                if alignments[line[0]][0]==E:
                    alignments[line[0]][1]+=[SS]
                    alignments[line[0]][2]+=[myfunction]
                else:#check that an evalue was diff from the best hit
                     #it means all the other hits for a sequenceX are not best hits
                     #it might look odd, but it reduces the number of times we get into the loop
                     #and it means less parsing time for large files
                     ID=line[0]
            except:
                #when we see a new sequence we clean the hash with the subsystems to save memory for big files
                #but before we see the subsystems and normalize if needed by the number of hits to a specific sequence
                alignments_count(alignments)
                alignments={line[0]:[E,[SS],[myfunction]]}
    f.close()
    alignments_count(alignments)#save the last sequence alignments results
    
    return [subsystems,functions]

#This function writes the SUPER-FOCUS results 
#INPUT: Hash with the subsystems assignments
#OUTPUT: Results in the levels 1, 2, and 3
def write_results(subsystems_assignments):
    DBpk={}

    #loads the subsystems PK (primary keys) 
    f=open(systemDbLocation+"/database_PKs.txt");f.readline()
    for line in f:
        line=line.replace("\n","").replace("\r","").split("\t")
        DBpk[line[0]]=line[1:]        
    f.close()

    #results
    for level in [0,1,2]:
        results={}
        for subsystem in subsystems_assignments:
            assigned_level=DBpk[subsystem][level]
            counts=subsystems_assignments[subsystem]

            try:results[assigned_level]+=counts
            except:results[assigned_level]=counts
            
        #writes the results for each subsystem
        K=results.keys();K.sort();S=sum(results.values())*1.
        o=open(project_output+"/"+project_name+"__results_subsystem_level_"+str(level+1)+".xls","w+")
        o.write("Query: "+query+"\nDB: DB_"+str(mydb)+"\n\nSubsystem Level "+str(level+1)+"\tNumber of assignments\tRelative Abundance (%)\n")
        for subsystem in K:
            o.write(subsystem+"\t"+str(results[subsystem])+"\t"+str((results[subsystem]/S)*100)+"\n")
        o.close()


#This function writes the SUPER-FOCUS functions assignments 
#INPUT: Hash with the functions assignments
#OUTPUT: Results in the function level
def write_functions(function_assignments):
        S=sum(function_assignments.values())*1.
        o=open(project_output+"/"+project_name+"__results_function_level.xls","w+")
        o.write("Query: "+query+"\nDB: DB_"+str(mydb)+"\n\nSEED Function\tNumber of assignments\tRelative Abundance (%)\n")
        for assignment in function_assignments:
            o.write(assignment+"\t"+str(function_assignments[assignment])+"\t"+str((function_assignments[assignment]/S)*100)+"\n")
        o.close()

def main2():
    p=0;a=1
    print "SUPER-FOCUS: A tool for agile functional analysis of shotgun metagenomic data"
    if myproject["-q"]=="":#check if the query is valid
        print "     -Invalid query file [-q]. Please select a FASTA/FASTQ file";p+=1
    if myproject["-db"].upper() not in ["DB_90","DB_95","DB_98","DB_100",]:#check if the DB is valid
        print "     -Invalid database [-db]. Please select DB_90","DB_95","DB_98 or","DB_100";p+=1
    if aligner not in ["blast","diamond","rapsearch"]:
         print "     -Invalid aligner [-a]. Please select blast","diamond"," or rapsearch";p+=1;a=0
    if a==1:#check if the aligner is installed
        a=myproject["-a"].lower()
        if "blast" in a:a=a+"x"
        if which(a)==None:#check if the aligner is valid
             print "     -"+aligner+" is not installed";p+=1
    if myproject["-n"] not in ["0","1"]:
        print "     -Invalid normalization [-n]. Please select 0 or 1";p+=1
    if myproject["-focus"] not in ["0","1"]:
        print "     -Invalid FOCUS selection [-focus]. Please select 0 or 1";p+=1
    if myproject["-p"] not in ["0","1"]:
        print "     -Invalid [-p]. Please select 0 or 1";p+=1
    if myproject["-dir"]=="":
        print "     -Invalid [-dir] output directory.";p+=1

    if p==0:# if it is 0, we can run SUPER-FOCUS
        if f=='1' and proteins==0:#Run FOCUS; otherwise, uses the default db
            organisms=runNParseFocus()
            mydb=formatDb(organisms)
        else:
            mydb=myproject["-db"].split("_")[-1]

        print "1) Aligning sequences in "+query+" against DB_"+mydb+" using "+aligner
        test=align(mydb)

        #test=1
        if test!=0:
            print "2) Alignment complete with success"
            #read alignments
            print "3) Reading alignments"
            parse_myproject=parse_alignments()
            subsystems_assignments=parse_myproject[0]
            functions_assignments=parse_myproject[1]

            print "4) Read alignments"
            #write results
            print "5) Writing results"
            write_results(subsystems_assignments)
            write_functions(functions_assignments)
            print "Done :) Please check the '"+project_output+project_name+"' folder"
            #os.system("rm "+project_output+"/"+my_alignments+".m8")
        else:
            print "There was a problem in the alignment process"

try:
    version=float(".".join((sys.version).split()[0].split(".")[:2]))#version of python installed
    if "-h" in sys.argv:
        print options
    elif (2.6<=version<3.0)==False:
        print "SUPER-FOCUS requires Python >= 2.6.X or < Python 3.Z"
    else:
        main2()
except:
    import traceback
    print traceback.format_exc()
    print options
