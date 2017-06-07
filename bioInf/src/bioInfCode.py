'''
Created on Jun 6, 2017

@author: weesh
'''

import re

def revtrans(str):
# Reverse transcribe sequence
    str2 = str[::-1]
    str2 = str2.replace('A', 't').replace('T', 'a').replace('G', 'c').replace('C', 'g')
    str2 = str2.upper()
    return str2

def translate(str):
# Translate sequence
    str2 = str.replace('T', 'U')
    table = {'UU' : 'F', 'CU' : 'L', 'AUG':'M', 'AU':'I', 'GU':'V', 'UC':'S', 'CC' : 'P', 'AC' : 'T', 'GC' : 'A', 'UA' : 'Y', 'UAA' : '*', 'UAG' : '*', 'CA' : 'H', 'CAA' : 'Q', 'CAG' : 'Q', 'AA' : 'N', 'AAA' : 'K', 'AAG' : 'K', 'GA' : 'D', 'GAA' : 'E', 'GAG' : 'E', 'UG' : 'C', 'UGA' : '*', 'UGG' : 'W', 'CG' : 'R', 'AG' : 'R', 'GG' : 'G'}
    ans = ""
    for i in range(0, len(str), 3):
        slice = str2[i:i+3]
        if slice[:3] in table:
            ans+=table[slice[:3]]
        else:
            ans+=table[slice[:2]]
    return ans

def fastaRead():
# Fasta Reader

    file = open('fasta.in', 'r')
    lib = {}
    str = "" 
    title = ""
    for line in file:
        line = line.rstrip('\n')
        if line[0]=='>':
            if len(str)>0:
                lib[title] = str
            title = line[1:]
            str=""
        else:
            str+=line
    lib[title] = str
    return lib
def probeRNAP(str):   
# RNAse P Finder 
    
#     str = (str[::-1][:500][::-1]) # remove me

    searchObj1 = re.findall('AG.{4}AT.{10,60}ACA.AA.{4}G.{2}TA.{0,600}GAGGAA.{2}TC.{4}C', str,)
    searchObj2 = re.findall('GAGGAA.{2}TC.{4}C.{0,600}AG.{4}AT.{10,60}ACA.AA.{4}G.{2}TA', str,)
#     searchObj3 = re.search('.*GAGGAA.{2}TC.{4}C.{0,600}AG.{4}AT.{10,60}AC.*', str)
#     searchObj4 = re.search('.*AG.{4}AT.{10,60}AC.*', str,)
#     print(len(str),searchObj1, searchObj2)
#     print(len(str),len(searchObj1), len(searchObj2))
    if(str[0]=='N'):
        for str in searchObj2:
            print(str)
    elif(str[0]=='R'):
         for str in searchObj2:
            print(revtrans(str))
#     print(len(searchObj2))

    return (len(searchObj1)+len(searchObj2))
def processFasta():
# processFasta and probes for RNAP
    lib = fastaRead()
#     print(lib)
    for key in lib:
#         print("CODING ",len(lib[key]))
        if probeRNAP(lib[key])>0:
            print(key.rstrip(), "CONTAINS RNAse P", probeRNAP('N'+lib[key]), "time(s)")
        elif probeRNAP(revtrans(lib[key]))>0:
            print(key.rstrip(), "CONTAINS RNAse P", probeRNAP('R'+revtrans(lib[key])), "time(s)")
        else: 
            print(key.rstrip(), "DOES NOT CONTAIN RNAse P")

def merge(arr, str):
# merges sequences, array of possible additions, string

#     print(arr)
    if(len(arr)==0):
        return str
    sign = 0
    maxSim = -1
    next = ""
    for txt in arr:
        sim1 = sim(str, txt)
        if(abs(sim1)>maxSim):
            next=txt
            maxSim=abs(sim1)
            if sim1>0:
                sign=-1
            else:
                sign=1
    
    maxSim *= sign
    
#     if(len(str)<100):
#         print(maxSim, next, str)
    arr.remove(next)
    if not maxSim==len(next):
        if maxSim>0:
            return merge(arr, str+next[maxSim:])
        else:
#             print(next[:len(next)-(-1*maxSim)])
            return merge(arr, next[:len(next)-(-1*maxSim)]+str)
    else:
        return merge(arr, str)
def sim(src, txt):
# returns the number of shared characters

#     print (src, txt)
    if txt in src:
        return len(txt)
    else:
        for i in range(len(txt)-1, -1, -1):
            if(src[:i]==txt[len(txt)-i:]):
                return i
            elif(src[len(src)-i:]==txt[:i]):
                return -1*i
def jigsaw():
# merges a jigsaw of sequences
    arr = []
    file = open('jigsaw.in', 'r')
    for line in file:
        for word in line.split():
            arr.append(word[0:len(word)-1])
    return merge(arr, '')


''' main method below '''

# str = input("Enter a DNA sequence: ")
# print("The reverse transcribed sequence is:",revtrans(str))
# print("The translated sequence is:",translate(str))

# processFasta()

# print(jigsaw())



