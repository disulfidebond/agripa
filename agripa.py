#!/usr/bin/python

import sys
import subprocess
import argparse
import re
import time # used for debugging
import md5

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
    if cdsStart > cdsStop:
        tmp = cdsStop
        cdsStop = cdsStart
        cdsStart = tmp
    blastString = fasta[cdsStart:cdsStop]
    fastaBLAST = open('queryBLAST.fasta', 'w')
    queryString = ">query\n"
    queryString += blastString
    if len(queryString) == 0:
        print "ERROR"
        print queryString
        sys.exit()
    fastaBLAST.write(queryString)
    fastaBLAST.close()
    blastout = subprocess.check_output(['blastn', '-query', 'queryBLAST.fasta', '-db', DBBLAST, '-gapopen', '5', '-gapextend', '2', '-reward', '2', '-penalty', '-3', '-word_size', '11', '-outfmt', outfmtString])
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

def parseFastaHeader(h):
    pList = []
    multiP = 0
    try:
        parsedHeader = h.split(' ')
        pList = parsedHeader
    except AttributeError:
        print "Error with"
        print h
        sys.exit()
    try:
        checkHeader = h.split(')')
        if len(checkHeader) > 2:
            multiP = 1
    except AttributeError:
        print "Error parsing header line:"
        print h
        sys.exit()
    p = re.compile("\(")
    returnRes = ""
    if not multiP:
        for i in pList:
            m = re.match(p, i)
            if m:
                resHeader = i.split(')')
                resHeader = resHeader[0]
                returnRes = resHeader[1:]
                returnRes = returnRes.upper()
                return returnRes
    else:
        resHeader = h.split(')')
        n = len(resHeader)
        nPop = n-2 # n-(x+1) where x = item index, +1 for zero indexed
        resHeaderItem = resHeader.pop(nPop)
        resHeaderParsed = resHeaderItem.split('(')
        returnRes = resHeaderParsed.pop()
        returnRes = returnRes.upper()
        return returnRes
    return returnRes

def createFastaHashTable(f):
    # return values are
    # fastaDict -> {md5Key: fastaEntry, md5Key: fastaEntry}
    # tupleDictList -> [(iCt, {md5Key: fastaEntry}), (iCt, {md5Key: fastaEntry})]
    tupleDictList = []
    p = re.compile(">")
    iCt = 0
    with open(f) as fasta:
        tupleDict = dict()
        fastaDict = dict()
        for i in fasta:
            m = re.match(p, i)
            i = i.rstrip('\n')
#            iParsed = parseFastaHeader(i)
            if m:
                iParsed = parseFastaHeader(i) # patch fix: moved to prevent unnecessary iterations
                dictListKey = md5.md5(i).hexdigest()
                tupleDict[dictListKey] = i
                fastaDict[dictListKey] = i
                tupleDictList.append((iCt, tupleDict))
                tupleDict = dict()
                fastaDict[dictListKey] = iParsed
                iCt+=1
            else:
                continue
        return (fastaDict, tupleDictList)
def retrieveStringFromFastaWithIndex(c):
    p = re.compile(">")
    iCt = 0
    s = ""
    matchedHeader = 0
    f = 'GCF_000002035.5_GRCz10_rna.fna'
    with open(f) as fasta:
        for i in fasta:
            m = re.match(p, i)
            i = i.rstrip('\n')
            if m:
                if iCt == c:
                    matchedHeader = 1
                    iCt += 1
                    continue
                else:
                    iCt += 1
                    matchedHeader = 0
                    continue
            else:
                if matchedHeader:
                    s += i
                else:
                    continue
        return s

def listOfGeneNames(l):
    filedata = []
    fileDict = dict()
    print l
    with open(l) as inputdatafile:
        for dline in inputdatafile:
            filedata.append(dline.rstrip())
    for i in filedata:
        isplit = i.split(',')
        k = isplit[0]
        v = isplit[1]
        v = v.upper()
        fileDict[k] = v
    return fileDict


def parseBLASTresults(b):
    # blastParsedList = []
    blastOutList = b.split('\n')
    resultBLAST = blastOutList[0]
    resultBLASTList = resultBLAST.split(',')
    try:
	bitscoreOfBLAST = resultBLASTList.pop()
        evalueOfBLAST = resultBLASTList.pop()
        try:
            bitscoreOfBLAST = float(bitscoreOfBLAST)
            if bitscoreOfBLAST > 40.0:
                # print "# evalue" + evalueOfBLAST
                matchSubjectString = resultBLASTList[0]
                matchSubjectID = resultBLASTList[1:]
                matchSubjectID = ','.join(matchSubjectID)
                # print matchSubjectID + "," + matchSubjectString
                return (matchSubjectID, matchSubjectString)
        except ValueError:
            print "# Warning: Bitscore for " + resultBLAST + "unreadable!"
            try:
                evalueOfBLAST = float(evalueOfBLAST)
                if evalueOfBLAST < 0.001:
                    matchSubjectString = resultBLASTList[0]
                    matchSubjectID = resultBLASTList[1:]
                    matchSubjectID = ','.join(matchSubjectID)               
                    return (matchSubjectID, matchSubjectString)
            except ValueError:
                return (-1, -1)
    except IndexError:
        return (-1, -2)
    return (-1, 0)

def geneDictLookup(geneDict, q):
    s = ""
    for k,v in geneDict.items():
        if v == q:
            s = k
            return s
    return s
def geneListMatch(g, fastaTuple):
    for ct in fastaTuple:
        pval = ct[1]
        if g in pval:
            return (ct[0], pval[g])
    return (-1,-1)
      
def pyBLAST_recursive(geneList, queryGeneDict, queryTupleDictList, subjectDB, queryDB):
    resList = []
    for k,v in geneList.items():
        print "Checking Entry"
        print k + "\t" + v
        geneMD5Key_query = geneDictLookup(queryGeneDict, v)
        if len(geneMD5Key_query) == 0:
            print "Error, not found in lookup"
            geneMD5Key_query = 'NULL'
            print "EntryID Gene BLASTGene BLAST_Homologue R_BLAST_Homologue"
            print "# " + k + "\t" + v + "\t" + "ERROR" + "\t" + "LOOKUPERROR" + "\t" + "NONE"
            continue # no match
        matchResult = geneListMatch(geneMD5Key_query, queryTupleDictList)
        if matchResult[-1] == -1:
            # no match in fasta tuple due to error
            print "Error in hash"
            print geneMD5Key_query
            sys.exit()
        fastaString = retrieveStringFromFastaWithIndex(matchResult[0])
        with open('queryBLAST.fasta', 'w') as ffile:
            ffile.write(fastaString)
        outfmtString = "10 sseq stitle sseqid evalue bitscore"
        DBBLAST = subjectDB
        blastoutput = subprocess.check_output(['blastn', '-query', 'queryBLAST.fasta', '-db', DBBLAST, '-gapopen', '5', '-gapextend', '2', '-reward', '2', '-penalty', '-3', '-word_size', '11', '-outfmt', outfmtString])
        bsubjectParsed = parseBLASTresults(blastoutput)
        if bsubjectParsed[0] == -1:
            print "error with parsing"
            print "EntryID Gene BLASTGene BLAST_Homologue R_BLAST_Homologue"
            print "# " + k + "\t" + v + "\t" + "ERROR" + "\t" + "PARSERROR" + "\t" + "NONE"
            continue
        pResult = parseFastaHeader(bsubjectParsed[0])
        print "Subject BLAST\n###\n"
        print bsubjectParsed
        if pResult:
            DBBLAST = queryDB
            with open('queryBLAST.fasta', 'w') as ffile:
                ffile.write(bsubjectParsed[1])
            blastoutput = subprocess.check_output(['blastn', '-query', 'queryBLAST.fasta', '-db', DBBLAST, '-gapopen', '5', '-gapextend', '2', '-reward', '2', '-penalty', '-3', '-word_size', '11', '-outfmt', outfmtString])
            bqueryParsed = parseBLASTresults(blastoutput)
            # print pResult
            # print bsubjectParsed[0]
            if bqueryParsed[0] == -1:
                print "Error or No match with query BLAST"
                print "EntryID Gene BLASTGene BLAST_Homologue R_BLAST_Homologue"
                print "# " + k + "\t" + v + "\t" + pResult + "\t" + bsubjectParsed[0] + "\t" + "BLASTERROR_R_BLAST_NONE"
                d = dict()
                d["ID"] = k
                d["gene"] = pResult
                d["catfish_gene"] = bsubjectParsed[0]
                d["species_homologous"] = "NONE"
                resList.append(d)
                continue
            bQueryResult = parseFastaHeader(bqueryParsed[0])
            bQueryResult = bQueryResult.upper()
            print "SPECIES HOMOLOGUE " + bQueryResult
            print bqueryParsed[0]
            print bqueryParsed[1]
            if pResult == bQueryResult:
                if bQueryResult == v:
                    print "EntryID Gene BLASTGene BLAST_Homologue R_BLAST_Homologue"
                    print "# " + k + "\t" + v + "\t" + pResult + "\t" + bsubjectParsed[0] + "\t" + bqueryParsed[0] 
                    d = dict()
                    d["ID"] = k
                    d["gene"] = pResult
                    d["catfish_gene"] = bsubjectParsed[0]
                    d["species_homologous"] = bqueryParsed[0]
                    resList.append(d)
                    continue
                else:
                    print "EntryID Gene BLASTGene BLAST_Homologue R_BLAST_Homologue"
                    print "# " + k + "\t" + v + "\t" + pResult + "\t" + bsubjectParsed[0] + "\t" + "NOT_HOMOLOGOUS"
                    continue
            else:
                # no match, update results
                print "EntryID Gene BLASTGene BLAST_Homologue R_BLAST_Homologue"  
                print "# " + k + "\t" + v + "\t" + pResult + "\t" + bsubjectParsed[0] + "\t" + "BLASTSPECIES_NONE"
                d = dict()
                d["ID"] = k
                d["gene"] = pResult
                d["catfish_gene"] = bsubjectParsed[0]
                d["species_homologous"] = "NONE"
                resList.append(d)
                continue
        else:
            print "Error parsing!"
            print bsubjectParsed
            print blastoutput
            print matchResult
            sys.exit()
    return resList


def main():
    parser = argparse.ArgumentParser(description='AGRIPA: Annotation and Genomic Reference fIle PArsing tool')
    cmd = ""
    # inputFiles
    parser.add_argument('-i', '--input', dest='input', action='append', help='Input file, if apporopriate, provide multiple files with -i FILE -i FILE.  See Help docs for more information.')
    
    # fasta file for homology
    parser.add_argument('-f', '--fasta', dest='ffile', action='append', help='Fasta file with homology command, otherwise ignored.')
    
    # output files
    parser.add_argument('-o', '--output', dest='output', action='append', help='Output file, default filename is outputfile.txt, it is a very good idea to change this to something else.')
    
    # action taken
    parser.add_argument('-c', '--cmd', dest='cmd', choices=['merge', 'homology', 'format'], help='command can be merge, homology, or format')
    
    # blast ID
    parser.add_argument('-b', '--blastdb', dest='blast', action='append', help='blast database prefix id, one per fasta file is required')
    # arguments: merge, format, homology
    
    # blast homology names
    parser.add_argument('-sid', '--subjectID', dest='sid', action='append', help='if using homology command, this option is the name of the BLAST database that is queried against.  The BLAST db must be in the current working directory.  If homology is selected this option is required, otherwise this option is ignored.')
    parser.add_argument('-qid', '--queryID', dest='qid', action='append', help='if using homology command, this option is the name of the BLAST database that will be BLASTED against the subject.  The BLAST db must be in the current working directory.  If homology is selected this option is required, otherwise this option is ignored.')
    
    # formatted list of gene names, for homology parsing only
    parser.add_argument('-gList', '--genelist', dest='glist', action='append', help='if using homology command, a comma-separated list of gene names that will be queried is required, otherwise this option is ignored')
    
    args = parser.parse_args()
    
    infile = args.input
    outfile = args.output
    cmd = args.cmd
    blastDB_ID = args.blast
    homology_subjectID = args.sid
    homology_queryID = args.qid
    hfastaFile = args.ffile
    gListIn = args.glist
    
    homologyFastaList = []
    
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
        if not homology_subjectID:
            print 'You must provide a BLAST homology subject name that matches a BLAST database in the current working directory'
            sys.exit()
        if not homology_queryID:
            print 'You must provide a BLAST homology query name that matches a BLAST database in the current working directory'
            sys.exit()
        homologyFastaList = createFastaHashTable(hfastaFile)
        gList = listOfGeneNames(gListIn)
        subjectDB = homology_subjectID
        queryDB = homology_queryID
        l_dict_query = homologyFastaList[0]
        l_tupleDict_query = homologyFastaList[1]
        l_dict_subject = ""
        l_tuple_subject = ""
        resList = pyBLAST_recursive(gList, l_dict_query, l_tupleDict_query, subjectDB, queryDB)
        for i in resList:
            print i
        
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
                            mergedEntry = sourceDict['gene'] + "_" + comparisonDict['gene']
                            d["mergedEntry"] = "mergedGene" + "_" + mergedEntry
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
            # check for existence or absence of keys
            # the keys ID startCDS stopCDS mergedEntry will always be present
            if not 'locus_tag' in entryDict:
                entryDict['locus_tag'] = "None"
            if not 'gene' in entryDict:
                entryDict['gene'] = "None"
            if not 'product' in entryDict:
                entryDict['product'] = "None"
            for k,v in entryDict.items():
                if k == "ID":
                    s = k + "=" + v + "\t"
                    outList.append((1, s))
                if k == "gene":
                    s = k + "=" + v + "\t"
                    outList.append((3, s))
                if k == "locus_tag":
                    s = k + "=" + v + "\t"
                    outList.append((2, s))
                if k == "product":
                    s = k + "=" + v + "\t"
                    outList.append((6, s))
                if k == "startCDS":
                    s = k + "=" + v + "\t"
                    outList.append((4, s))
                if k == "stopCDS":
                    s = k + "=" + v + "\t"
                    outList.append((5, s))
                if k == "mergedEntry":
                    s = k + "=" + v + "\t"
                    outList.append((7, s))
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
