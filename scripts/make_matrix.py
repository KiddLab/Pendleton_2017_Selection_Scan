# make-matrix.py
# A. Pendleton
# runs commands to format files for and generating the matrix files for genotypes
# This script requires:
# 1. Input file of sample names (example sample file is provided below in script)
# 2. Path for input/output files. This script requires there to be a file with the same 
#	   filehandle root name to have an eigenstrat genotype file AND a plink .bim file 
#	   (e.g. FILE.bim and FILE.eigenstratgeno will become FILE_forMatrix.txt)
# 3. The FILE_forMatrix.txt file generated here will be the input file for matrix2png	
import glob
import numpy as np
import os
import sys
from optparse import  OptionParser

###############################################################################
USAGE = """
python make-matrix.py 		--samples < file with sample data > 
							--directory < Directory to write output files >

samples == sample name file that you want to appear on the left of your matrix2png file. 
directory == Directory with eigenstrat genotype (.eigenstratgeno) and .bim (PLINK) files, AND where output files will be written
"""

parser = OptionParser(USAGE)
parser.add_option('--samples',dest='samples', help = 'file with sample name data')
parser.add_option('--directory',dest='directory', help = 'Directory for input and output files')
(options, args) = parser.parse_args()

if options.directory is None:
    parser.error('output directory not given')
if options.samples is None:
    parser.error('sample data file not given')

###############################################################################
## EXAMPLE OF SAMPLE FILE
#The sample file input is simply the name that you want to appear on your matrix2png 
#	output plot. If you want the samples in a particular order, we have found it is best
#	to provide a number in front of your sample name / identifier, that can be used to 
#	sort the matrix after this step is complete. File structure is sample name + '\n' 
"""
04_Saluki_1233
09_AfghanHound_1735
38_Sub-Saharan_1756
33_LabradorRetriever_2972
14_Xolo_4669
18_GreatDane_6610
29_ScottishTerrier_8542
28_ToyPoodle_10442
15_Chihuahua_13131
08_SiberianHusky_14529
"""
###############################################################################


print 'Reading in results from the following directory', options.directory 

#Generating the row column IDs
sampleFile = open(options.samples, 'r')
sample_array = [] #sample array
#first line in your output file MUST be snp id. 
#So we add it in the line below as a "sample" for simplification in matrix construction
sample_array.append('SNP_ID') 

for line in sampleFile: #read through the sample file
	line = line.split()
	if '#' in line[0]:
		continue
	sampleID = line[0]
	sample_array.append(sampleID) #add sample ID to sample array
		
windowCount = 0  #keeps track of window count 


fHD_array = [] #filehandle array
#This script assumes that the input files are all in the directory provided in options.directory
#Within options.directory, you must have files that are in eigenstrat genotype AND .bim (PLINK) file format for your region of interest. 
#The below function (glob) will search for any file that ends with 'eigenstratgeno' 
for file in glob.glob(options.directory+"*.eigenstratgeno"):
	print file 

	name=file.replace(".eigenstratgeno","")
	FH=open(file, 'r')
	windowCount+=1
    
	outfile = name + '_forMatrix.txt'
	outFile = open(outfile, 'w')
	print 'Writing new matrix to', outfile
	
	bimfile = name + '.bim'
	bimFile = open(bimfile, 'r')
	snpList = []
	for line in bimFile:
		line = line.split()
		snpList.append(line[1])
		lenSNP = len(snpList)
		snp = '\n'.join(map(str,snpList))
	outFile.write("#SampleID\t"+"\t".join(map(str,snpList))+"\n")
	
	lineNum = 0
	array=[]
	for line in FH:
		lineNum += 1
		line=list(line.rstrip())
		for i in range(0,len(line)):
			if '9' in line[i]: #9 = missing data identifier for matrix2png. If you dont want missing data, then use plink to remove missing sites
				line[i] = '-'
		array.append(map(str,line))
	
	print 'Generated transposed matrix for %i SNPs' % (lineNum)
	matrix=np.array(array)
	print 'matrix dimensions', matrix.shape
	print matrix
	new_matrix=matrix.T
	print 'new matrix dimensions', new_matrix.shape
	for i in range(len(new_matrix)):
		outFile.write(sample_array[i+1]+"\t"+"\t".join(map(str,new_matrix[i]))+"\n")