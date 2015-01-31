#################################
##           CCH  v1           ##
##        October  2014        ##
## A. P. Levine and D. P. Gale ##
#################################
#  Acknowledgments J. M. Levine

#Based on cCHv1.45 and Dv1.47

#Release Version 1.1
#06/10/14

#Import packages
import sys
#import csv, numpy
import time
import itertools
import math
from operator import itemgetter
import argparse

#Print the program name, version and current time
print
print "     #################################"
print "     ##            CCH v1           ##"
print "     ##         October  2014       ##"
print "     ## A. P. Levine and D. P. Gale ##"
print "     #################################"
print
print "Combinatorial Conflicting Homozygosity"
print

execution_time = time.localtime()
print time.strftime("%A %d/%m/%Y %H:%M:%S", execution_time) +"\n"

##see http://docs.python.org/dev/library/argparse.html and http://docs.python.org/dev/howto/argparse.html#combining-positional-and-optional-arguments
parser = argparse.ArgumentParser(prog="CCH release v1",description="")
parser.add_argument("-snp",help="SNP file [required]")
parser.add_argument("-ind",help="individuals file")
parser.add_argument("-name",help="analysis name [default OUTPUT]")
parser.add_argument("-number",help="number of individuals to analyse [default ALL]")
parser.add_argument("-length_snp","-ls",help="SNP no CH length threshold [default 100]")
parser.add_argument("-length_cm","-lc",help="cM no CH length threshold [default 2]")
parser.add_argument("-chr",help="chromosome(s) (e.g. 1,2,3) or all [default ALL]")
parser.add_argument("-progress",help="progress report off [default ON]",action="store_true")
parser.add_argument("-duplicates",help="duplicate removal off [default ON]",action="store_true")
parser.add_argument("-plot",help="plot results off [default ON]",action="store_true")
args = parser.parse_args()
#required, default
#SNP file
if args.snp==None:
	print("**ERROR** missing SNP file!\nStopping...\nSee usage instructions below:\n")
	parser.print_help()
	exit()
if args.snp!=None:
	print("SNP file:\t\t\t\t"+args.snp)
#Individuals file
if args.ind==None:
	print("Individuals file:\t\t\t"+"NOT PROVIDED")
if args.ind!=None:
	print("Individuals file:\t\t\t"+args.ind)
#Output file name
if args.name==None:
	print("Output name:\t\t\t\t"+"OUTPUT [by default]")
	log="output"
if args.name!=None:
	print("Output name:\t\t\t\t"+args.name)
	log=args.name

#individuals_position is a dictionary that stores the column number of each individual in the data file
individuals_positions=dict([])

#Load the data file
data_file=open(args.snp,"r")
#all_individuals stores the names of all individuals in the data file
all_individuals=dict([])
count=0
for i in data_file.next().rstrip("\n").split("\t")[6:]: #6 because of the preceding information columns
	all_individuals[i]=count
	count=count+1
print "Number of individuals in SNP file:\t"+str(len(all_individuals))

if args.ind==None:
	print("Number of individuals provided:\t\t"+str(len(all_individuals))+" [all by default]")

if args.ind!=None:
	#Load the list of individuals
	individuals_file=open(args.ind,"r")
	#analysis_individuals stores the names of the individuals to analyse
	analysis_individuals=[]
	#Add each individual to analysis_individuals
	for line in individuals_file:
		split=line.rstrip("\n").split("\t")
		analysis_individuals.extend(split)
	#Number of individuals in input file
	print("Number of individuals provided:\t\t"+str(len(analysis_individuals)))

if args.ind==None:
	analysis_individuals=all_individuals

#Number of individuals to choose
if args.number==None:
	print("Number of individuals to choose:\t"+str(len(analysis_individuals))+" [all by default]")
	choose=len(analysis_individuals)
if args.number!=None:
	print("Number of individuals to choose:\t"+args.number)
	choose=int(args.number)
#Minimum length of CH region (SNPs)
if args.length_snp==None:
	print("Minimum length of no CH region (SNPs):\t100 [by default]")
	threshold_snp=100
if args.length_snp!=None:
	print("Minimum length of no CH region (SNPs):\t"+args.length_snp)
	threshold_snp=int(args.length_snp)
#Minimum length of CH region (cM)
if args.length_cm==None:
	print("Minimum length of no CH region (cM):\t4 [by default]")
	threshold_cm=4
if args.length_cm!=None:
	print("Minimum length of no CH region (cM):\t"+args.length_cm)
	threshold_cm=float(args.length_cm)
#Chromosome(s)
if args.chr==None:
	print("Chromosome(s):\t\t\t\tALL [by default]")
	chr=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22']
if args.chr!=None:
	print("Chromosome(s):\t\t\t\t"+args.chr)
	if args.chr=="all" or args.chr=="ALL":
		chr=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22']
	else:
		chr=[]
		for j in args.chr.split(","):
			chr.append(j)

#Progress report
if args.progress==False:
	print("Print progress report:\t\t\tON [by default]")
	print_progress=1
if args.progress==True:
	print("Print progress report:\t\t\tOFF")
	print_progress=0
#Duplicate removal
if args.duplicates==False:
	print("Duplicate removal:\t\t\tON [by default]")
	duplicate_removal=1
if args.duplicates==True:
	print("Duplicate removal:\t\t\tOFF")
	duplicate_removal=0
#Plotting
if args.plot==False:
	print("Plot results:\t\t\t\tON [by default]")
	to_plot=1
	import matplotlib
	import matplotlib.pyplot as plt
if args.plot==True:
	print("Plot results:\t\t\t\tOFF")
	to_plot=0

print
if print_progress==1:
	print "Checking SNP file for individuals..."
for i in analysis_individuals:
	temp=all_individuals.get(i,"missing")
	if temp=="missing":
		print "**ERROR** there are no data for individual "+i+" !\nStopping..." #individual from analysis_individuals not found in all_individuals
		exit()
	if temp!="missing":
		individuals_positions[i]=temp+6
if print_progress==1:
	print "...Done!\nAll individuals found\n"

if print_progress==1:
	print "Loading SNP data to memory..."
	if len(analysis_individuals)>30:
		print "There are a lot of individuals in the SNP data file: this may take a while or the program might crash!"
#array will store the contents of the data file
array=[]
this_line=[] #this is a temporary field used to build the array line by line
#Loop through each line of the data_file
count=0
stored_count=0 #length of stored array
for line in data_file:
	#Split (tab delimited) and remove end of line marker
	split=line.rstrip("\n").split("\t")
	this_chr=split[0]
	if this_chr in chr:
		#Keep the first four columns
		this_line.append(split[0])
		this_line.append(split[1])
		this_line.append(split[2])
		this_line.append(split[3])
		#Count the progress
		stored_count=stored_count+1
		#Keep the columns that the individuals to be analysed are in (as above)
		for i in individuals_positions.values():
			this_line.append(split[i])
		#Store this information in array
		array.append(this_line)
		this_line=[]
	count=count+1
	if count%10000==0 and print_progress==1:
		sys.stdout.write(".")
		sys.stdout.flush()
if print_progress==1:
	print "..Done!"
	print "Total number of SNPs:\t"+str(stored_count)+"\n"

#new_positions is a dictionary which stores the position of each individual in the array (which is different as the array only includes the individuals to be kept for the analysis)
new_positions=dict([])
position = 4 #4 because of the preceeding information columns
#Loop through each individual and define the column they are in
for i in individuals_positions.keys():
	new_positions[i]=position
	position = position + 1

total=len(individuals_positions) #total is the total number of individuals in the individuals file
#choose is the number of individuals to choose, e.g. 17C10 total=17, choose=10
combinations = itertools.combinations(individuals_positions, choose) #combinations is a list of all possible combinations of choose individuals
number_of_analyses = math.factorial(total)/(math.factorial(choose)*(math.factorial(total-choose))) #calculate the number of analyses to be run

#Print the number of analyses to be run
if print_progress==1:
	print "Total number of combinations: " + str(number_of_analyses)
	print

#Open a log file to which the results will be written
logfile = open(log+".txt", "w")
#Prepare and save the header
ID=[]
for i in range(1,choose+1):
	ID.append(str("ID"+str(i))) #ID1 through to ID[choose]
logfile.write("Len_SNP\tStart_chr\tStart_bp\tStart_cM\tStart_SNP\tEnd_chr\tEnd_bp\tEnd_cM\tEnd_SNP\tLen_bp\tLen_cM\t"+"\t".join(ID)+"\n")

previous_line=['1'] #Set the previous line to '1' for the first loop
count=0 #count stores the number of preceding SNPs with no CH
region=[] #region temporarily stores the details of a region of no CH
analysis_count=0 #this counts the number of analyses that have been done
for each_combination in combinations: #Loop through each combination of choose individuals
	analysis_count=analysis_count+1 #Increase the analysis count
	if print_progress==1:
		print "Running combination", analysis_count, "of", number_of_analyses #Print the current analysis number
	columns=[]
	for i in each_combination:
		columns.append(new_positions[i]) #Identify the columns of the individuals to be analysed for the current selection
	index=0 #index is the row number in the array
	for i in array:
		index=index+1 #index=1 is the first line
		previous_count = count #Define the previous_count as the count from the last row of the loop through the array
		line_info = i[:4] #line_info stores the SNP details (name, position)
		to_analyse = []
		for j in columns: #Loop through each of the columns of relevance (as defined earlier)
			to_analyse.append(i[j]) #Append the data from the coumns of relevance to to-analyse
		occurrences = (to_analyse.count("AA"), to_analyse.count("BB")) #store the number of occurrences of AA and BB
		if max(occurrences[0], occurrences[1])==float(occurrences[0]+occurrences[1]): #If the maximum of the number of occurrences of AA and BB equals the total number of occurences of AA and BB then there is no CH
			score=1 #There is no CH so score=1
		else: #If the maximum number of occurrences of AA and BB does not equal the total number of occurrences of AA and BB there must be CH
			score=0 #There is CH so score=0
		if score==1: #If there is no CH keep counting
			count = count + 1
		if score<1:
			count = 0 #If there is CH reset the count
		if (previous_line[0] != line_info[0]): #If there is a chromosome transition (the current chromosome does not equal the previous row's chromosome) reset the count
			count = 0
		if index == stored_count: #Last SNP in map
			count = 0
		if index==3: #This is the first row of real data
			start = i #start stores the beginning of a region of no CH (at the beginning start it as the first row)
		if count == 0 and index>2: #If there is CH (and we are further than the second row, i.e. real data) store information from the start of the region of no CH and the previous line (the end of the region of no CH) in the temp variable
			temp=previous_count, start[0], start[1], start[2], start[3], previous_line[0], previous_line[1], previous_line[2], previous_line[3],round((float(previous_line[1])-float(start[1])),2),round((float(previous_line[2])-float(start[2])),2)
			for k in temp: #Append each element of the temp variable (beginning and end of the last region of no CH) to the region variable
				region.append(k)
			for k in each_combination: #Append the names of the individuals currently being analysed to the region variable
				region.append(k)
			start = i #The current line is the beginning of a new region (the run of no CH ended in the previous row)
		if region!=[]: #If we are at the end of a region of no CH (as there are data in the region variable) (this could have been included in the previous if test)...
			if region[0]>threshold_snp and region[10]>threshold_cm: #If the length of no CH in SNPs exceed the user-defined cut off paramater...
				logfile.write("\t".join(map(str,region))) #Write the region of no CH to the output
				logfile.write("\n") #Move to a new line in the output file
		previous_line = i #Set previous_line to the current line (i.e. the current line will be the previous_line in the next row)
		region=[] #Reset the region variable (if it had data in it then it would have just been written)
		if index%10000==0 and print_progress==1:
			sys.stdout.write(".")
			sys.stdout.flush()
	if print_progress==1:
		print "..Done!"

logfile.close()
if print_progress==1:
	print

if duplicate_removal==1:
	data_file=open(log+".txt", "r") #Load file for duplicate removal
	logfile = open(log+"_dr.txt", "w") #Output file (duplicates removed)
	logfile.write("Len_SNP\tStart_chr\tStart_bp\tStart_cM\tStart_SNP\tEnd_chr\tEnd_bp\tEnd_cM\tEnd_SNP\tLen_bp\tLen_cM\n") #Header
	if print_progress==1:
		print "Removing duplicates..."
	for each in chr: #Loop through each chromosome
		data_file.seek(0) #Go to the start of the data file
		this_chr=str(each) #this_chr is the first chromosome
		chr=[]
		for i in data_file: #Extract data for this_chr
			s=i.rstrip("\n").split()
			if s[1]==this_chr:
				chr.append(s) 
		chr_1=sorted(chr,key = lambda x:(int(x[2]),-int(x[6])),reverse=False) #Sort this_chr by start position
		previous_start=0 #Remove stage 1 overlaps
		chr_2=[]
		for i in chr_1: 
			if int(i[2])>previous_start:
				chr_2.append(i[0:11])
			previous_start=int(i[2])
		chr_3=sorted(chr_2,key = lambda x:(int(x[6]),int(x[2])),reverse=False) #Sort this_chr by end position
		previous_end=0 #Remove stage 2 overlaps
		chr_4=[]
		for i in chr_3:
			if int(i[6])>previous_end:
				chr_4.append(i[0:11])
			previous_end=int(i[2])
		for i in chr_4: #Output results for this_chr
			logfile.write("\t".join(map(str,i))+"\n")
		if print_progress==1:
			sys.stdout.write(".")
			sys.stdout.flush()
	logfile.close()
	if print_progress==1:
		print "..Done!"
		print

chr=[280,266,227,210,212,194,192,170,164,183,157,176,133,129,138,135,140,124,115,106,70,76] #Chromosome lengths (cM)
chr_sum=[0,280,546,773,983,1195,1389,1581,1751,1915,2098,2255,2431,2564,2693,2831,2966,3106,3230,3345,3451,3521] #Summated chromosome lengths (cM)
chr_mid=[140,413,659.5,878,1089,1292,1485,1666,1833,2006.5,2176.5,2343,2497.5,2628.5,2762,2898.5,3036,3168,3287.5,3398,3486,3559] #Summated chromosome mid-points (cM)
total=3597 #Total length (cM)

if to_plot==1:
	if print_progress==1:
		print "Plotting graph..."
	for i in chr_sum:
		plt.axvline(x=i, color="gray",linestyle = "--")
	plt.axhline(y=threshold_cm,color="green")
	data=open(log+"_dr.txt", "r")
	first=0
	y_max=[]
	count=0
	for i in data:
		s=i.rstrip("\n").split()
		if first>0:
			chr=int(s[1]) #Chromosome
			start=float(s[3]) #Start position (cM)
			end=float(s[7]) #End position (cM)
			length=float(s[10]) #Length (cM)
			y_max.append(length) #List of lengths (for y axis limit)
			if chr%2 == 0: #even
				col="cyan"
			if chr%2 == 1: #odd
				col="blue"
			plt.hlines(y=length,xmin=chr_sum[chr-1]+start,xmax=chr_sum[chr-1]+end,color=col,lw=1) #Horizontal line
			plt.vlines(x=chr_sum[chr-1]+start,ymin=0,ymax=length,color=col,lw=1) #Vertical line at start
			plt.vlines(x=chr_sum[chr-1]+end,ymin=0,ymax=length,color=col,lw=1) #Vertical line at end
		first=1
		count=count+1
		if count%100==0 and print_progress==1:
			sys.stdout.write(".")
			sys.stdout.flush()
	plt.xlim([0,total]) #X-axis limits (whole genome)
	if len(y_max)>0:
		plt.ylim([0,max(y_max)+1]) #Y-axis limits
	plt.xticks(chr_mid, range(1,23),fontsize=8,rotation=15) #X-axis ticks
	plt.ylabel("Length of no conflicting homozygosity (cM)") #Y-axis label
	plt.xlabel("Chromosome") #X-axis label
	plt.title(log+"\n"+time.strftime("%d/%m/%Y", execution_time))
	plt.figtext(0.9,0.04,"CCH v1",color="grey")
	#plt.show()
	plt.savefig(log+".pdf", papertype="a4", format="pdf",orientation="landscape")
	if print_progress==1:
		print "..Done!\n"

print "Finished!"

execution_time = time.localtime()
print time.strftime("%A %d/%m/%Y %H:%M:%S", execution_time)
