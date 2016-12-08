# RNASeq Analysis Pipeline v 0.0.1  #
# 03/08/2016 Serdar Turkarslan      #
# Institute for Systems Biology     #

import glob, sys, os, string

# Create STAR genome index
def genomeIndex(genomeDir,starDir,genomeFasta):
    indexCmd = '%s --runMode genomeGenerate --runThreadN 4 --genomeDir %s --sjdbGTFfile %s --sjdbOverhang 74 --genomeFastaFiles %s --genomeSAindexNbases 9' %(starDir, genomeDir, genomeGFF, genomeFasta)
    print indexCmd
    print "\033[34m Indexing genome... \033[0m"
    os.system(indexCmd)

# -o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o- #

# Function to run and control pipeline
def runPipeline():
    fastqFolders = glob.glob('%s/*/' %(dataDir)) # Modify for location of your read folders
    print
    print 'FASTQ Folders: %s' %fastqFolders
    print
    # Start processing each folder
    for folder in fastqFolders:
        sampleFolder = folder.split("/")[1] # Folder containing reads
        fastqFilesFirst = glob.glob('%s/*R1*.fastq.gz' %(folder)) # Get list of 1st pair files
        for file in fastqFilesFirst:
            print
            print "\033[33m Processing %s \033[0m" %(file)
            fileName = file.split("_R1")[0] # Read file name
            
            print "FileName: %s" %fileName
            sampleResultsDir = resultsDir+ '/'+fileName.split("/")[1] # folder name in the results folder
            
            print "sampleResultsDir: %s" %sampleResultsDir
            sampleTitle = fileName.split("/")[2]
            
            libraryName = sampleTitle.split("_")[0] 
            
            lane = fileName.split("/")[2].split("_")[2]
            
            sampleName = fileName.split("/")[2].split("_")[1]
            
            sampleId = fileName.split("/")[2].split("_")[0] # uncomment for test run
            print "FileName: %s Lane: %s sample: %s ID: %s Library name: %s" %(fileName, lane, sampleName, sampleId, libraryName)
            print
            
            # Create read names
            firstPair = fileName + "_R1_001.fastq.gz" # or "_R1_001.fastq.gz"
            secondPair = fileName + "_R2_001.fastq.gz"
            
            print "First Pair: %s, Second Pair: %s" %(firstPair, secondPair)
            print

            ############ Run Functions ##############
            # Run trimmomatic
            firstPairTrimmedPaired, secondPairTrimmedPaired = runTrim(firstPair, secondPair)
        # Run STAR
        runStar(firstPairTrimmedPaired, secondPairTrimmedPaired, libraryName, sampleResultsDir, folder)
        
        #Run HTSeq count
        runHTseq(libraryName, sampleResultsDir)

        # Run cufflinks
        runCufflinks(libraryName, sampleResultsDir)
        

    return firstPair, secondPair, libraryName, sampleResultsDir, folder
  
# -o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o- #      


# Trim reads with Trimmomatic
def runTrim(firstPair, secondPair):
    print
    print "\033[34m Running Read Trimming... \033[0m"
    
    # Program Parameters
    illuminaClip = "ILLUMINACLIP:/users/sturkars/Trimmomatic-0.35/adapters/TruSeq3-PE.fa:2:30:10" #Remove adapters 
    leading = "LEADING:3" #Remove leading low quality or N bases
    trailing = "TRAILING:3" #Remove trailing low quality or N bases
    slidingWindow = "SLIDINGWINDOW:4:20" #Scan the read with a 4-base wide sliding window, clipping average quality per base < 20
    minLen = "MINLEN:36"
    paired = "PE" # or SE for single end
    threads = "-phred33 -threads 8"
      
    # define result files
    filesFolder = firstPair.split('/')[0] + "/" + firstPair.split('/')[1]
    firstPairTrimmedPaired = filesFolder+"/trimmed/"+firstPair.split('.fastq.gz')[0].split("/")[2] + "_paired_trimmed.fastq.gz"
    secondPairTrimmedPaired = filesFolder+"/trimmed/"+secondPair.split('.fastq.gz')[0].split("/")[2] + "_paired_trimmed.fastq.gz"
    firstPairTrimmedUnpaired = filesFolder+"/trimmed/"+firstPair.split('.fastq.gz')[0].split("/")[2] + "_unpaired_trimmed.fastq.gz"
    secondPairTrimmedUnpaired = filesFolder+"/trimmed/"+secondPair.split('.fastq.gz')[0].split("/")[2] + "_unpaired_trimmed.fastq.gz"
    trimDir = filesFolder+"/trimmed/"
    
    # create trim folder
    if not os.path.exists('%s' %(trimDir)):
        print '\033[31m %s directory does NOT exists.Creating one! \033[0m' %(trimDir)  
        os.makedirs('%s' %(trimDir))
    else:
        print '\033[31m %s directory exists. Not creating. \033[0m' %(trimDir)
    
    # define command
    cmd = '/users/sturkars/java/bin/java -Xmx128m -jar %s %s %s %s %s %s %s %s %s %s %s %s %s %s' %(trimmomaticPath, paired, threads, firstPair, secondPair, firstPairTrimmedPaired, firstPairTrimmedUnpaired, secondPairTrimmedPaired, secondPairTrimmedUnpaired, illuminaClip, leading, trailing, slidingWindow, minLen)
    print "Trimmomatic Command: ", cmd
    
    # Run trimmomatic
    os.system(cmd)
    return firstPairTrimmedPaired, secondPairTrimmedPaired
    print
    print
# -o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o- #


# Run STAR alignment
def runStar(firstPairTrimmedPaired, secondPairTrimmedPaired, libraryName, sampleResultsDir, folder):
    print
    print "\033[34m Running STAR alignment... \033[0m"
    # collect list of recalibrated bam files
    print folder
    # We need to comma seperate all firt read pairs and then second read pairs for ReadFilesIn
    readFilesInFirst = glob.glob('%strimmed/*_R1_*_paired_trimmed.fastq.gz' %(folder))
    readFilesInSecond = []
    # Get the first read pair filename and replace R1 with R2 to create second read file name
    for file in readFilesInFirst:
        secondPair = file.replace("_R1_", "_R2_")
        readFilesInSecond.append(secondPair)
    readFilesIn1stJoined = ",".join(readFilesInFirst)
    readFilesIn2ndJoined = ",".join(readFilesInSecond)
    
    # create results folder
    if not os.path.exists('%s' %(sampleResultsDir)):
        print '\033[31m %s directory does NOT exists.Creating one! \033[0m' %(sampleResultsDir)  
        os.makedirs('%s' %(sampleResultsDir))
    else:
        print '\033[31m %s directory exists. Not creating. \033[0m' %(sampleResultsDir)
    
    starOptions ="--runThreadN 8 --outSAMattributes All --genomeLoad NoSharedMemory --outFilterType Normal --alignSJoverhangMin 5 --alignSJDBoverhangMin 3 --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.4 --outFilterScoreMinOverLread 0.4 --quantMode TranscriptomeSAM GeneCounts"
    outFileNamePrefix = sampleResultsDir+"/"+libraryName+"_star_"
    
    # STAR run command
    cmd = '%s --genomeDir %s %s --readFilesIn %s %s --outFileNamePrefix %s --sjdbOverhang 74' %(starPath, genomeDir, starOptions, readFilesIn1stJoined, readFilesIn2ndJoined, outFileNamePrefix)
    # resort bam files for cufflinks
    cmd2 = '%s sort %s/%s_star_Aligned.sortedByCoord.out.bam %s/%s_Aligned.Sorted' %(samtoolsPath, sampleResultsDir, libraryName, sampleResultsDir, libraryName)
    # index sorted bam file
    cmd3 = '%s index %s/%s_Aligned.Sorted.bam' %(samtoolsPath, sampleResultsDir, libraryName)
    print "Star Command: ", cmd
    print cmd2
    print cmd3
    
    # Run star
    os.system(cmd)
    os.system(cmd2)
    os.system(cmd3)
 
# -o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o- #

# Run htseq counts
def runHTseq(libraryName, sampleResultsDir):
    print "\033[34m Running HTSeq counts... \033[0m"
    
    # create results folder
    htseqFolder = '%s/htseq-counts-exon' %(resultsDir)
    if not os.path.exists(htseqFolder):
        print '\033[31m %s directory does NOT exists.Creating one! \033[0m' %(htseqFolder)  
        os.makedirs(htseqFolder)
    else:
        print '\033[31m %s directory exists. Not creating. \033[0m' %(htseqFolder)
    
    htseqInputFile = sampleResultsDir+"/"+libraryName+"_star_Aligned.sortedByCoord.out.bam"
    cmd = '/users/sturkars/anaconda2/bin/python -m HTSeq.scripts.count -s "reverse" -t "exon" -r pos -f bam %s %s > %s/%s_htseqcounts-exon.txt' %(htseqInputFile,genomeGFF,htseqFolder,libraryName) 
    print cmd
    
    # htseq run command
    os.system(cmd)

# -o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o- #

# Run cufflinks
def runCufflinks(libraryName, sampleResultsDir):
    print "\033[34m Running Cufflinks... \033[0m"
    
    # create results folder
    cufflinksFolder = '%s/cufflinks' %(resultsDir)
    if not os.path.exists(cufflinksFolder):
        print '\033[31m %s directory does NOT exists.Creating one! \033[0m' %(cufflinksFolder)  
        os.makedirs(cufflinksFolder)
    else:
        print '\033[31m %s directory exists. Not creating. \033[0m' %(cufflinksFolder)
    
    cufflinksInputFile = sampleResultsDir+"/"+libraryName+"_Aligned.Sorted.bam"
    cufflinksOutputDir = sampleResultsDir+"/cufflinks"

    # Cufflinks run command
    cmd = '%s --library-type fr-firststrand -p 8 -o %s -G %s %s' %(cufflinksPath, cufflinksOutputDir, genomeGFF, cufflinksInputFile)
    
    # Cuffquant run command
    cmd2 = '%s -p 8 -o %s %s %s' %(cuffquantPath, cufflinksOutputDir, genomeGFF, cufflinksInputFile)

    print cmd
    print cmd2
    
    # Run cuffdiff and cuffquant
    os.system(cmd)
    os.system(cmd2)     

# -o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o- #

# Run cuffmerge
def runCuffmerge():
    print "\033[34m Running Cuffmerge... \033[0m"
    cufflinkFiles = glob.glob('%s/*/cufflinks/transcripts.gtf' %(resultsDir))
    outputFile = '%s/assembly_GTF_list.txt' %(resultsDir)
    with open(outputFile, 'w') as g:
        for file in cufflinkFiles:
            g.write('%s\n' %file)
    
    cmd = '%s -p 8 -g %s -o %s/cuffmerge_stats %s' %(cuffmergePath, genomeGFF, resultsDir, outputFile)
    print cmd
    
    # Cuffmerge command
    os.system(cmd)  

# -o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o- #

# Run cuffdiff
def runCuffdiff(sampleResultsDir):
    print "\033[34m Running Cuffdiff... \033[0m"
    cxbFiles = glob.glob('%s/*/cufflinks/*.cxb' %(resultsDir))
    cxbManifestFile = '%s/all_cxb_manifest.txt' %(resultsDir)
    cxbManifestLabelFile = '%s/all_cxb_manifest_labels.txt' %(resultsDir)
    # with open(cxbManifestFile, 'w') as h:
    with open(cxbManifestFile, 'r') as s:
        for line in s:
            cxbLine = line.split('\n')[0]
    # create labels for each sample
    with open(cxbManifestLabelFile, 'r') as t:
        for line in t:
            cxbLabel = line.split('\n')[0]        
            
    transcriptGTF = '%s/cuffmerge_stats/merged.gtf' %resultsDir
    outputDir = '%s/cuffdiff_output' %(resultsDir)
    cmd = '%s -p 8 -library-type fr-firststrand -L %s -o %s %s %s' %(cuffdiffPath, cxbLabel, outputDir, genomeGFF, cxbLine)
    print cmd
    
    # run cuffdiff
    os.system(cmd)            
    
# -o-o-o-o-o-o-o-o-o-o-o-o-  Data and Results folders  -o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o- # 

# Input files    
dataDir = "dvh-rnaseq-27882859" 
genomeDir = "reference-03" 
resultsDir = "results-03"
fastqcDir = '%s/fastqc' %(resultsDir)
cuffmergeDir = '%s/cuffmerge' %(resultsDir)
genomeFasta = glob.glob('%s/*GCA_000195755.1.30.dna.genome.fasta' %(genomeDir))
genomeFasta = ' '.join(genomeFasta)
genomeGFF = glob.glob('%s/*GCA_000195755.1.30.gtf' %(genomeDir))
genomeGFF = ' '.join(genomeGFF)
genomeMask = glob.glob('%s/*GCA_000195755.1.30.gtf' %(genomeDir))
genomeMask = ' '.join(genomeMask)  
tophatGenome = '%s/desulfovibrio' %(genomeDir)
rsemGenome = 'reference-rsem-dvh/Desulfovibrio-rsem'

#Create Results directory
if not os.path.exists('%s' %(resultsDir)):
    print '\033[34m %s directory does NOT exists.Creating one! \033[0m' %(resultsDir)  
    os.makedirs('%s' %(resultsDir))
else:
    print '\033[34m %s directory exists. Not creating. \033[0m' %(resultsDir)


# -o-o-o-o-o-o-o-o-o-o-o-o-  Program paths  -o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o- # 

starPath = "/users/sturkars/STAR-STAR_2.4.2a/bin/Linux_x86_64/STAR" 
trimmomaticPath = "/users/sturkars/Trimmomatic-0.35/trimmomatic-0.35.jar"
fastqc = "/users/sturkars/FastQC/fastqc" 
cufflinksPath = "/users/sturkars/cufflinks-2.2.1/cufflinks"
cuffmergePath = "/users/sturkars/cufflinks-2.2.1/cuffmerge"
cuffquantPath = "/users/sturkars/cufflinks-2.2.1/cuffquant"
cuffdiffPath = "/users/sturkars/cufflinks-2.2.1/cuffdiff"
samtoolsPath = "/users/sturkars/samtools-1.2/samtools"

# -o-o-o-o-o-o-o-o-o-o-o-o-  Run functions  -o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o- # 
# Genome index
genomeIndex(genomeDir,starPath,genomeFasta)

# Pipeline
firstPair, secondPair, libraryName, sampleResultsDir, folder = runPipeline()

# Cuffmerge
runCuffmerge()

# cuffdiff
runCuffdiff(sampleResultsDir)
