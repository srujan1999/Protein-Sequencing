"""
Protein Sequencing Project
Name:
Roll Number:
"""

import hw6_protein_tests as test

project = "Protein" # don't edit this

### WEEK 1 ###

'''
readFile(filename)
#1 [Check6-1]
Parameters: str
Returns: str
'''
def readFile(filename):
    read_file=open(filename,"r").read().splitlines()
    string="".join(read_file)
    return string


'''
dnaToRna(dna, startIndex)
#2 [Check6-1]
Parameters: str ; int
Returns: list of strs
'''
def dnaToRna(dna, startIndex):
    replace_dna=dna[startIndex:].replace('T', 'U')
    Rna_list=[]
    for i in range(0,len(replace_dna),3):
        index=replace_dna[i:i+3]
        Rna_list.append(index)
    for i in Rna_list:
        if i=="UAA" or i=="UAG" or i=="UGA":
            word_index=Rna_list.index(i)
            return Rna_list[:word_index+1]
    return Rna_list
'''
makeCodonDictionary(filename)
#3 [Check6-1]
Parameters: str
Returns: dict mapping strs to strs
'''
def makeCodonDictionary(filename):
    import json
    open_file=open(filename,"r")
    data=json.load(open_file)
    codon_dict={}
    for i,j in data.items():
        for k in j:
            replace=k.replace('T','U')
            codon_dict[replace]=i
    return codon_dict


'''
generateProtein(codons, codonD)
#4 [Check6-1]
Parameters: list of strs ; dict mapping strs to strs
Returns: list of strs
'''
def generateProtein(codons, codonD):
    amino_list=[]
    for i in codons:
        if len(amino_list)==0 and i=="AUG":
            amino_list.append("Start")
        else:
            amino_list.append(codonD[i])
    return amino_list


'''
synthesizeProteins(dnaFilename, codonFilename)
#5 [Check6-1]
Parameters: str ; str
Returns: 2D list of strs
'''
def synthesizeProteins(dnaFilename, codonFilename):
    overall_protein=[]
    unused_base=0
    read_file=readFile(dnaFilename)
    codon_dict=makeCodonDictionary(codonFilename)
    startIndex=0
    while len(read_file)>startIndex:
        find_ATG=read_file[startIndex:startIndex+3]
        if find_ATG!="ATG":
            startIndex+=1
            unused_base+=1
        else:
            Rna_seq=dnaToRna(read_file, startIndex)
            protein_data=generateProtein(Rna_seq, codon_dict)
            overall_protein.append(protein_data)
            startIndex+=3*len(overall_protein[-1])
    print("total number of bases:",len(read_file),"unused-base count:",unused_base,"total number of proteins synthesized:",len(overall_protein))
    return overall_protein


def runWeek1():
    print("Human DNA")
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    print("Elephant DNA")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")


### WEEK 2 ###

'''
commonProteins(proteinList1, proteinList2)
#1 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs
Returns: 2D list of strs
'''
def commonProteins(proteinList1, proteinList2):
    uniq_protein=[ i for i in proteinList1 if i in proteinList2]
    return uniq_protein


'''
combineProteins(proteinList)
#2 [Check6-2]
Parameters: 2D list of strs
Returns: list of strs
'''
def combineProteins(proteinList):
    protein_list=[j for i in proteinList for j in i ]
    return protein_list


'''
aminoAcidDictionary(aaList)
#3 [Check6-2]
Parameters: list of strs
Returns: dict mapping strs to ints
'''
def aminoAcidDictionary(aaList):
    amino_dicts={}
    for i in aaList:
        amino_dicts[i]=aaList.count(i)
    return amino_dicts


'''
findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
#4 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs ; float
Returns: 2D list of values
'''
def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
    protein_list1 = combineProteins(proteinList1)
    protein_list2 = combineProteins(proteinList2)
    amino_dict1 = aminoAcidDictionary(protein_list1)
    amino_dict2 = aminoAcidDictionary(protein_list2)
    amino_acids= []
    for i in amino_dict1:
        if i not in amino_acids and i != "Start" and i != "Stop":
            amino_acids.append(i)
    for i in amino_dict2:
        if i not in amino_acids and i != "Start" and i != "Stop":
            amino_acids.append(i)    
    amino_freqs= []     
    for acid in amino_acids:  
        freq1 = 0
        freq2 = 0
        if acid in protein_list1:
            freq1 = (amino_dict1[acid]/ len(protein_list1))
        if acid in protein_list2:
            freq2 = (amino_dict2[acid]/len(protein_list2))
        diff = freq2 -freq1
        if (diff > cutoff) or (diff < -cutoff):
            amino_freqs.append([acid, freq1, freq2])
    return amino_freqs


'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):
    return


def runWeek2():
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")

    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)


### WEEK 3 ###

'''
makeAminoAcidLabels(proteinList1, proteinList2)
#2 [Hw6]
Parameters: 2D list of strs ; 2D list of strs
Returns: list of strs
'''
def makeAminoAcidLabels(proteinList1, proteinList2):
    return


'''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):
    return


'''
createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
#4 [Hw6] & #5 [Hw6]
Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
Returns: None
'''
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):
    import matplotlib.pyplot as plt
    return


'''
makeEdgeList(labels, biggestDiffs)
#5 [Hw6]
Parameters: list of strs ; 2D list of values
Returns: list of strs
'''
def makeEdgeList(labels, biggestDiffs):
    return


'''
runFullProgram()
#6 [Hw6]
Parameters: no parameters
Returns: None
'''
def runFullProgram():
    return


### RUN CODE ###

# This code runs the test cases to check your work
if __name__ == "__main__":
    # print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
    # test.week1Tests()
    # print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
    #runWeek1()
    #test.testReadFile()
    #test.testDnaToRna()
    #test.testMakeCodonDictionary()
    #test.testGenerateProtein()
    #test.testSynthesizeProteins()
    ## Uncomment these for Week 2 ##
    #test.testCommonProteins()
    #test.testCombineProteins()
    #test.testAminoAcidDictionary()
    test.testFindAminoAcidDifferences()
    """
    print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    test.week2Tests()
    print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    runWeek2()
    """

    ## Uncomment these for Week 3 ##
    """
    print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    test.week3Tests()
    print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()
    """
