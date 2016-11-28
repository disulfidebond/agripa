#!/usr/bin/python

import sys
import subprocess
import argparse
import re
import time # used for debugging

def debugPrint(l):
    for v in l:
        kv = v
        for k,v in kv.items():
            print k,v
    print ""
    print "######"
    print ""

def preliminaryCheck(lfile, infile, outfile, cmd):
    if len(outfile) != 1:
	print "Error, only 1 output filename is allowed!  Please check your input command."
	return 0
    if cmd == 'format' and len(infile) > 1:
        print "Error, only 1 input file is allowed if reformatting an annotation file!  Please check your input command."
        return 0
    if cmd == 'merge':
	failcheck = 0
	if len(infile) != 4:
            print "Error!  It looks like you do not have the correct number of files for input"
            print "\nThe merge command requires, in this order, each with the \"-i\" option:\n\nThe fasta file with extension .fasta\nThe annotation file with extension .gtf"
            print "The fasta file with extension .fasta for the genome to be merged\nThe annotation file with extension .gtf for the genome to be merged.\n\n"
            return 0
        else: # sourceSeq == inputfile1, mergeSeq == inputfile2
            checkSourceSeqs_fasta = infile[0].split('.')
            try:
                checkSourceSeqs_fasta = checkSourceSeqs_fasta[1][0:2]
            except IndexError:
                    print "Please check file and file type of"
                    print infile[0]
                    print "before proceeding"
                    return 0
            checkSourceSeqs_gtf = infile[1].split('.')
            try:
                checkSourceSeqs_gtf = checkSourceSeqs_gtf[1][0]
            except IndexError:
                print "Please check file and file type of"
                print infile[1]
                print "before proceeding"
                return 0
            checkMergeSeqs_fasta = infile[2].split('.')
            try:
                checkMergeSeqs_fasta = checkMergeSeqs_fasta[1][0:2]
            except IndexError:
                print "Please check file and file type of"
                print infile[2]
                print "before proceeding"
                return 0
            checkMergeSeqs_gtf = infile[3].split('.')
            try:
		checkMergeSeqs_gtf = checkMergeSeqs_gtf[1][0]
            except IndexError:
                print "Please check file and file type of"
                print infile[3]
                print "before proceeding"
                return 0
            if checkSourceSeqs_fasta != "fa":
                print "Please check file "
                print infile[0]
                print "before proceeding"
                failcheck = 1
            if checkSourceSeqs_gtf != "g":
                print "Please check file "
                print infile[1]
                print "before proceeding"
                failcheck = 1
            if checkMergeSeqs_fasta != "fa":
                print "Please check file "
                print infile[2]
                print "before proceeding"
                failcheck = 1
            if checkMergeSeqs_gtf != "g":
                print "please check file "
                print infile[3]
                print "before proceeding"
                failcheck = 1
            if failcheck:
                return 0
    return 1
	

def parseCol8(c):
    rtpl = []
    tsplit = c.split(';')
    for i in tsplit:
        tmp = i.split("=")
        rtpl.append((tmp[0], tmp[1]))
    return rtpl
    
def parseCDS(c):
    return (c[3], c[4])  
    
# parse FASTA file
def fastaOpenAndReturnList(f):
    s = []
    p = re.compile(">")
    with open(f) as dfile:
        for l in dfile:
            m = re.match(p, l)
            if m:
                continue
            else:
                lAdd = l.rstrip('\n')
                s.append(lAdd)
        return s

def fastaListToString(f):
    s = ""
    for i in f:
        s += i
    return s
def gtfOpenAndParse(l): #        # note: only compatible with GTF
    r = []
    with open(l) as inputdatafile:
        for dline in inputdatafile:
            dline = dline.rstrip('\n')
            if len(dline) != 0:
                t = dline[0]
                if t != "#":
                    itmSplit = dline.split('\t')
                    cdsTuple = parseCDS(itmSplit)
                    parsedColumn = parseCol8(itmSplit[8])
                    d = dict()
                    try:
                        parsedColumn = parseCol8(itmSplit[8])
                        d = dict((x,y) for x,y in parsedColumn)
                    except TypeError:
                        d_tmp = itmSplit.pop()
                        d_tmp_kv = d_tmp.split("=")
                        d[d_tmp_kv[0]] = d_tmp_kv[1]
                    d["startCDS"] = cdsTuple[0]
                    d["stopCDS"] = cdsTuple[1]
                    r.append(d)
        return r

def pyBLAST(cdsStart, cdsStop, fasta, dbID):
    bsplit = []
    bsplit_res = []
    outfmtString = "7 qseqid sseqid evalue sstart send"
    # default options for remote BLAST
    # blastn -query testquery.fasta -db human -gapopen 5 -gapextend 2 -reward 2 -penalty -3 -word_size 11
    DBBLAST = dbID
    blastString = fasta[cdsStart:cdsStop]
    fastaBLAST = open('queryBLAST.fasta', 'w')
    queryString = ">query\n"
    queryString += blastString
    fastaBLAST.write(queryString)
    fastaBLAST.close()
    # v = subprocess.check_call(["echo", queryString], stdout=fastaBLAST)
    # if v != 0:
    #    print "Error, did not write query"
    blastout = subprocess.check_output(['blastn', '-query', 'queryBLAST.fasta', '-db', DBBLAST, '-gapopen', '5', '-gapextend', '2', '-reward', '2', '-penalty', '-3', '-word_size', '11', '-outfmt', outfmtString])
    # alternatively can return list and parse later
    # list valuse are 0:3 constant, 4: no. hits, 5: result queryID subj ID evalue start end 6: no. queries processed
    # return blastout
    blast_res = []
    try:
        bsplit = blastout.split("#")
        try:
            bsplit_res = bsplit[5].split("\n")
        except IndexError:
    		# no BLAST results
            return (-1, -1)
    except IndexError:
    	# error with BLAST
        return [(-1, -1)]
    # get no. of hits
    blastResultNumber = 0
    try:
        blastCountCheck = bsplit_res[0]
        blastCount = blastCountCheck.split(' ')
        blastResultNumber = int(blastCount[1])
    except ValueError:
        return [(-1, -1)]
    
    if not blastResultNumber:
        return [(-1, -1)] # error with BLAST result, or no results
    else:
        if blastResultNumber > 5:
            blastResultNumber = 5 # take only top 5 results if more than 5 exist
        else:
            blastResultNumber = blastResultNumber+1
        for i in xrange(1, blastResultNumber): # top 5 results, or less 
            bsplit_tmp = bsplit_res[i]
            bsplit_top_res = bsplit_tmp.split("\t")
            bsplit_tmpChk = float(bsplit_top_res[2])
            if (bsplit_tmpChk < 0.001):
                blast_res.append((bsplit_top_res[3], bsplit_top_res[4]))
        if len(blast_res) == 0:
            blast_res.append((-1, -1))
    return blast_res


def main():
    parser = argparse.ArgumentParser(description='AGRIPA: Annotation and Genomic Reference fIle PArsing tool')
    cmd = ""
    
    # other possibility is to use parser.add_argument('something', choices=['a', 'b', 'c'])
    # parser.add_args(['a']) tests
    
    # parser_merge = subparser.addparser('m', help='merge help')
    # parser_merge
    
    
    # inputFiles
    parser.add_argument('-i', '--input', dest='input', action='append', help='Input file, if apporopriate, provide multiple files with -i FILE -i FILE')
    
    # output files
    parser.add_argument('-o', '--output', dest='output', action='append', help='Output file, default filename is outputfile.txt, it is a very good idea to change this to something else.')
    
    # action taken
    parser.add_argument('-c', '--cmd', dest='cmd', choices=['merge', 'homology', 'format'], help='command can be merge, homology, or format')
    
    # blast ID
    parser.add_argument('-b', '--blastdb', dest='blast', action='append', help='blast database prefix id, one per fasta file is required')
    # arguments: merge, format, homology
    
    args = parser.parse_args()
    
    infile = args.input
    outfile = args.output
    cmd = args.cmd
    blastDB_ID = args.blast
    
    # prerun check
    lenList = len(infile)
    
    if cmd == 'merge': 
        if len(blastDB_ID) != 2:
            print "Error, you must provide two BLAST database prefix identifiers with the merge command"
            print "Use the syntax \"-b identifier -b identifier\" "
            print "Be ABSOLUTELY certain the order of the database identifiers matches the order of other fasta and annotation files"

    if not (preliminaryCheck(lenList, infile, outfile, cmd)):
        sys.exit()
    
    if cmd == 'format':
        print 'Sorry, this has not been implemented yet'
        sys.exit()
    elif cmd == 'homology':
        print 'Please select another option'
    	sys.exit()
    elif cmd == 'merge': # sourceSeq == inputfile 1, mergeSeq == inputfile2
        sourceSeqGtf = gtfOpenAndParse(infile[1])
        sourceSeqFasta = fastaOpenAndReturnList(infile[0])
        sourceSeqFastaString = fastaListToString(sourceSeqFasta)
        mergeSeqGtf = gtfOpenAndParse(infile[3])
        mergeSeqFasta = fastaOpenAndReturnList(infile[2])
        mergeSeqFastaString = fastaListToString(mergeSeqFasta)
        parseLen = len(sourceSeqGtf)
    	mergedDictList = []
    	for x in xrange(0, parseLen):
            # print "Iteration " + str(x)
            # debugPrint(mergedDictList)
            sourceDict = sourceSeqGtf[x]
            addMerged = 0
            sourceDB = blastDB_ID[0]
            if "gene" in sourceDict:
            # first option, attempt to match gene
                for y in mergeSeqGtf:
                    comparisonDict = y
                    if 'gene' in comparisonDict:
                        if comparisonDict['gene'] == sourceDict['gene']:
                            d = sourceDict
    	                    # mergedDict["gene"] = entryDict["gene"]
                            mergedEntry = sourceDict['gene'] + "_" + comparisonDict['gene']
                            d["mergedEntry"] = mergedEntry
                            mergedDictList.append(d)
                            addMerged = 1
                            break
                if addMerged == 1:
                    continue
            else:
                sourceDict_StartCDS = sourceDict["startCDS"]
                sourceDict_StopCDS = sourceDict["stopCDS"]
                sourceDict_startCDS_int = int(sourceDict_StartCDS)
                sourceDict_stopCDS_int = int(sourceDict_StopCDS)
                for v in mergeSeqGtf:
                    cmpDB = blastDB_ID[1]
                    comparisonDict = v
                    cmpDict_StartCDS = comparisonDict["startCDS"]
                    cmpDict_StopCDS = comparisonDict["stopCDS"]
                    cmp_startCDS_int = int(cmpDict_StartCDS)
                    cmp_stopCDS_int = int(cmpDict_StopCDS)                	
                    # compare startCDS and stopCDS
                    chkStart = abs(sourceDict_startCDS_int - cmp_startCDS_int)
                    chkStop = abs(sourceDict_stopCDS_int - cmp_stopCDS_int)
                    # chkLength = abs(source_stopCDS_int - source_startCDS_int)  
                    if chkStart < 250 and chkStop < 250 :
                        mergedEntry = sourceDict["ID"] + "CDS_Merge" + "_" + comparisonDict["ID"]
                        addMerged = 1
                        d = sourceDict
                        d["mergedEntry"] = mergedEntry
                        mergedDictList.append(d)
                        break
                if addMerged == 1:
                    continue
                sourceDictBLASTResult = pyBLAST(sourceDict_startCDS_int, sourceDict_stopCDS_int, mergeSeqFastaString, cmpDB)
                mEntry = ""
                if len(sourceDictBLASTResult) > 1:
                    mEntry = "_#"
                if sourceDictBLASTResult[0][0] == -1:
                    d = sourceDict
                    mergedEntry = d["ID"]
                    mergedEntry += "_" + "NOBLASTRESULT_NOMERGE" + mEntry
                    d["mergedEntry"] = mergedEntry
                    mergedDictList.append(d)
                    continue # no BLAST result
                test_cmp_startCDS_int = int(sourceDictBLASTResult[0][0])
                test_cmp_stopCDS_int = int(sourceDictBLASTResult[0][1])
                compareDictBLASTResult = pyBLAST(test_cmp_startCDS_int, test_cmp_stopCDS_int, sourceSeqFastaString, sourceDB)
                if compareDictBLASTResult[0][0] == -1:
                    d = sourceDict
                    mergedEntry = d["ID"]
                    mergedEntry += "_" + "BLASTRESULT_NOMERGE" + "_" + cmpDB + mEntry
                    d["mergedEntry"] = mergedEntry
                    mergedDictList.append(d)
                    continue # no second BLAST result
                else:
                    d = sourceDict
                    mergedEntry = d["ID"]
                    mergedEntry += "_" + "BLASTMERGE" + "_" + "CDS" + "_" + compareDictBLASTResult[0][0] + "_" + compareDictBLASTResult[0][1] + "_" + "CDS" + "_" + compareDictBLASTResult[0][0] + "_" + compareDictBLASTResult[0][1] + mEntry
                    d["mergedEntry"] = mergedEntry
                    mergedDictList.append(d)
                    continue
        outputfileItem = outfile[0] 
        for entry in mergedDictList:
            entryDict = entry
            outputfile = open(outputfileItem, 'a')
            outList = []
            stringOut = ""
            for k,v in entryDict.items():
                if k == "ID":
                    s = k + "=" + v + "\t"
                    outList.append((1, s))
                elif k == "gene":
                    s = k + "=" + v + "\t"
                    outList.append((3, s))
                elif k == "locus_tag":
                    s = k + "=" + v + "\t"
                    outList.append((2, s))
                elif k == "product":
                    s = k + "=" + v + "\t"
                    outList.append((4, s))
                else:
                    s = k + "=" + v + "\t"
                    outList.append((4, s))
            outList = sorted(outList)
            for i in outList:
                stringOut += i[1] + "\t"
            # print s
            outputfile.write(stringOut)
            outputfile.write('\n')
        outputfile.close()
    else:
        print 'Error, command not recognized!'
        sys.exit()
    
    
if __name__ == "__main__":
	main()
