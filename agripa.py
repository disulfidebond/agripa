#!/usr/bin/python

import re
import subprocess
import sys
import md5
import argparse

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
    # NOTE: temporary patch is assign f to subjectDB fastafile
    # since this method is only used with to retrieve the fasta string from the 
    # subject fasta index, however this should be properly updated
    # reminder: 
    # subject -> human, zebrafish, etc. that is BLASTed first
    # query -> catfish, the query for the subject
    #     the result from the subject BLAST 
    #     is BLASTed as a reciprocal search against catfish
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
                    # print matchSubjectID + "," + matchSubjectString
                    return (matchSubjectID, matchSubjectString)
            except ValueError:
                # print "NULL" # error with BLAST result output
                return (-1, -1)
    except IndexError:
        # print "NULL" # error with BLAST result
        return (-1, -2)
    # return blastParsedList
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
    # query == what is being BLASTED
    # subject == what is being BLASTed against with the query
    resList = []
    for k,v in geneList.items():
        # Step 1: match geneList name to parsed Catfish (query) name
        #
        # return md5hash from dict of k,v
        # where k == md5hash of fasta header, v == parsed fasta header
        print "Checking Entry"
        print k + "\t" + v
        geneMD5Key_query = geneDictLookup(queryGeneDict, v)
        if len(geneMD5Key_query) == 0:
            print "Error, not found in lookup"
            geneMD5Key_query = 'NULL'
            print "EntryID Gene BLASTGene BLAST_Homologue R_BLAST_Homologue"
            print "# " + k + "\t" + v + "\t" + "ERROR" + "\t" + "LOOKUPERROR" + "\t" + "NONE"
            continue # no match
#        matchResult = geneListMatch(geneMD5Key_query, subjectTupleDictList)
        # Step 2: retrieve gene name and fasta string, then print to tmp file for BLAST
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

# todo: test run
# check if not above

def main():
    # subject
    # zebrafish GCF_000002035.5_GRCz10_rna.fna
    # human human_rna.fna
# NOTE: Query fasta must be parsed, but subject fasta is already parsed by BLAST
# no need to parse it again
#    l = createFastaHashTable('human_rna.fna')
#    l_dict_subject = l[0]
#    l_tupleDict_subject = l[1]
    l_dict_subject = ""
    l_tuple_subject = ""
    # query
    l = createFastaHashTable('GCF_001660625.1_IpCoco_1.2_rna.fna')
    l_dict_query = l[0]
    l_tupleDict_query = l[1]
    parser = argparse.ArgumentParser(description='AGRIPA: Annotation and Genomic Reference fIle PArsing tool')
    parser.add_argument('-i', '--input', dest='input', help='Input gene list')
    args = parser.parse_args()
    gFileIn = args.input
#    gList = listOfGeneNames('geneList_Catfish.csv')
    gList = listOfGeneNames(gFileIn)
    subjectDB = 'Danio_rerio'
#    subjectDB = 'Human'
    queryDB = 'CatfishGene'
    resList = pyBLAST_recursive(gList, l_dict_query, l_tupleDict_query, subjectDB, queryDB)
    xLen = len(resList)
    print "Length of resList is"
    print xLen
#    s = ""
#    print "OUTPUT\n\n###\n\n"

#    with open('homologyoutputFile.txt', 'a') as hFileOut:
#        for x in xrange(0, xLen):
#            xVal = resList[x]
#            s += xVal["ID"] + "\t" + xVal["gene"] + "\t" + xVal["catfish_gene"] + "\t" +xVal["human_homologous"]
#            # print s
#            hFileOut.write(s)

if __name__ == "__main__":
	main()

