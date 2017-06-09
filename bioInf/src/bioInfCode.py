'''
Created on Jun 6, 2017

@author: weesh
'''

import re
from builtins import str

spacing = 100
"BLAH BLAH BLAH"

def matches(char1, char2):
    char1 = char1.replace('A', 't').replace('C', 'g').replace('G', 'c').replace('T', 'a').upper()
    return (char1 == char2) or ((char1=='C' and char2 == 'T') or (char1 == 'A' and char2 == 'G'))

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
def probeRNAP(string):   
# RNAse P Finder 
    
#     str = (str[::-1][:500][::-1]) # remove me

#     searchObj1 = re.findall('.{100}AG.{4}AT.{10,60}ACA.AA.{4}G.{2}TA.{0,600}GAGGAA.{2}TC.{4}C.{100}', str,)
    searchObj2 = re.findall('.{'+str(spacing)+'}GAGGAA.{2}TC.{4}C.{0,600}AG.{4}AT.{10,60}ACA.AA.{4}G.{2}TA.{'+str(spacing)+'}', string,)
#     searchObj3 = re.search('.*GAGGAA.{2}TC.{4}C.{0,600}AG.{4}AT.{10,60}AC.*', str)
#     searchObj4 = re.search('.*AG.{4}AT.{10,60}AC.*', str,)
#     print(len(str),searchObj1, searchObj2)
#     print(len(str),len(searchObj1), len(searchObj2))
#     if(str[0]=='N'):
#         for str in searchObj2:
#             print(str)
#     elif(str[0]=='R'):
#          for str in searchObj2:
#             print((str))
    
    return searchObj2
#     print(len(searchObj2))
#     return (len(searchObj1)+len(searchObj2))

def processFasta():
# processFasta and probes for RNAP
    lib = fastaRead()
#     print(lib)
    for key in lib:
#         print("CODING ",len(lib[key]))
#         if len(probeRNAP(lib[key]))>0:
#             print(key.rstrip(), "CONTAINS RNAse P", probeRNAP('N'+lib[key]), "time(s)")
        arr = probeRNAP(revtrans(lib[key]))
        
        if len(arr)>0:
            print(key.rstrip(), "CONTAINS RNAse P", len(arr), "time(s)")
            for string in arr:
                print(string)
                for num in range(2, 10):
                    ans = formatPin(string[0:spacing+1], num)
                    print(num,ans.count('X'), ans)
#                 print(formatPin(string[0:spacing+1], 10))                 
                
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

def hairpin(s, loopSize):
# finds hairpin containing end of string
#     print(s)
#     print(len(s))
#     print(s.find('TTCGG'))
#     print(s[76]+" "+s[83])
    length = len(s)
    arr = [[0 for x in range(length)] for y in range(length)]
    
    for dif in range(1, length):
        for start in range(0, length-dif):
#             print(dif)
            end = dif+start
            b = matches(s[start], s[end])
            if dif<=loopSize-1:
#                 if b:
#                     arr[start][end] = 0
#                 else:
#                      arr[start][end] = 1
                arr[start][end] = 0
            else:
                if b:
                    arr[start][end] = arr[start+1][end-1]
                else:
                    arr[start][end] = min(min(arr[start][end-1], arr[start+1][end])+1, arr[start+1][end-1]+2)
    
#     print(arr[77][82],arr[76][83], matches(s[76], s[83]), arr[72][87], matches(s[72], s[87]))
    start = 0
    end = length-1
    minP = 1.000000000000
    for count in range(0, end-2*loopSize-10):
        tempP = (arr[count][end] * 1.00000)/((end-count)*1.00000000)
        if tempP < minP:
            start = count
            minP = tempP
#     start = 59
#     print("BOOGIE ", start, minP)
#     print("DARN", 59, (arr[59][end]*1.00000)/((end-59)*1.000000))
#     print(s[start:int(start+(end-start)/2)])
    
    sArr = [x+60 for x in range(1, 40)]
#     print('BB', sArr)
#     for arrray in range(56, 80):
#         print(arrray, arr[arrray][61:])
    list = []
#     print s[start:end+1]
    bob = start
    
    while not abs(start-end)<=loopSize-1:
        
        if not matches(s[start], s[end]):
            if arr[start][end-1] > arr[start+1][end]:
                list.append(start)
                start+=1
            else:
                list.append(end)
                end-=1
        else:
            start+=1
            end-=1
    
    list.append(start)
    list.append(end)
    list.append(bob)
#     print(list)
    return list

def formatPin(string, loopSize):
    list = hairpin(string, loopSize)
#                 print(string[0:spacing+1])
    start = list.pop()

    intE = list.pop()
    intS = list.pop()
#     print(string[intS:intE+1])
#     print(list)
    
#     print("*",string[start:len(string)])
    ans = ""
    for ind in range(start, len(string)):
        if(not ind in list and not (ind<=intE and ind>=intS)):
            ans+=(string[ind])
#     print("|",ans,"W/O ERRORS")

    topL = ans[:int(len(ans)/2)]
    botL = ans[int(len(ans)/2):] 
    botL = botL[::-1]
    
    
    tupList = []
    for ind in range(len(topL)):
        tupList.append((topL[ind], botL[ind]))
        
    topStr = ""
    botStr = ""
     
    for ind in tupList:
        topStr+=ind[0]
        botStr+=ind[1]
     
#     print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
#     print(topStr)
#     print(botStr)
#     print(topL)
#     print(botL)
#     print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
    
    extra = 0
    count=0
    flip = False
    region = string[start:len(string)]
#     print("XXX", region)
#     print("YYY", ans)
#     print(string[start:intS],string[intE+1:])

    for ind in range(start, len(string)):
        if(ind in range(intS, intE+1)):
            continue
#         if(flip and int(len(ans)-count) < len(ans)):
#             print(ans[int(len(ans)-count)],count, ind, string[ind])
#         else:
#             print(ans[count], count, ind, string[ind],'.')
#         if(flip):
#             print(ind, count)
#             print(ind, string[ind], int(len(ans)-count), ans[int(len(ans)-count)])
#         else:
#             print(ind, string[ind], count, ans[count])
#         if(flip):
#             print (len(ans)-count, len(ans), ind, len(string), count, len(ans))
#         print(not (count==0 and flip))
        if((not (count==0 and flip)) and ((string[ind] == ans[count] and not flip) or (flip and string[ind] == ans[int(len(ans)-count)]))):
            if(count < len(ans)/2 and not flip):
                count+=1
            else:
                count-=1
                flip = True
        elif not flip:
#             print(string[ind], count, "BOOK")
            tupList.insert(count+extra, (string[ind], 'X'))
            extra+=1
            topStr = ""
            botStr = ""
            for ind in tupList:
                topStr+=ind[0]
                botStr+=ind[1]
             
#             print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
#             print(topStr)
#             print(botStr)
#         #     print(topL)
#         #     print(botL)
#             print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
        else:
#             print(string[ind], count, "BOOK")
            tupList.insert((count+1)+extra, ('X', string[ind]))    
            extra+=1
            topStr = ""
            botStr = ""
            for ind in tupList:
                topStr+=ind[0]
                botStr+=ind[1]
             
#             print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
#             print(topStr)
#             print(botStr)
#         #     print(topL)
#         #     print(botL)
#             print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
    
    topStr = ""
    botStr = ""
    
    for ind in tupList:
        topStr+=ind[0]
        botStr+=ind[1]
    
    return (topStr+' ['+string[intS:intE+1]+ '] '+botStr[::-1])
    
    
#     
#     count = 0
#     for ind in range(len(string)-1, start-1, -1):
#         topStr = ""
#         botStr = ""
#         
#         for ined in tupList:
#             topStr+=ined[0]
#             botStr+=ined[1]
#         print("XXXXXXXXXXXXXXXXXXXDDDDDDDDDDDDDDDDDDDDDDDXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
#         print(topStr)
#         print(botStr)
#         print(string[len(string)-1:start-1:-1])
#         print("XXXXXXXXXXXXXXXXXXXDDDDDDDDDDDDDDDDDDDDDDDXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
# #         print("CIJC", botL)
#         print(string[ind], ind-start, botStr[len(botStr)-1-count], (len(botStr)-1-count), "REKTREKTREKTREKTREKT")
#         if(count<len(botStr)):
#             if(string[ind]==botStr[len(botStr)-1-count]):
#                 count+=1
#             else:
#                 tupList.insert(count, (' ', string[ind]))
#                 count+=1
#         elif(count<2*len(botStr)):
#             if(string[ind]==topStr[count-len(botStr)]):
#                 count+=1
#             else:
#                 tupList.insert(count-len(botStr), (string[ind], ' '))
#                 count+=1
#         else:
#             tupList.append((string[ind], ' '))
#       
#     topStr = ""
#     botStr = ""
#     
#     for ind in tupList:
#         topStr+=ind[0]
#         botStr+=ind[1]
#     
#     print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
#     print(topStr)
#     print(botStr)
#     print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
# #     
# #     for ch in ans[0:int(len(ans)/2)+1]:
# #         if(count in list):
# #             topL.append(string[count])
# #             count+=1
# #         topL.append(ch)
# #         count+=1
# #     topL = topL[::-1]
# #     
# #     for ch in ans[int(len(ans)):int(len(ans)/2)]:
# #         if(count in list):
# #             botL.append(string[count])
# #             count+=1
# #         while not count>=2*len(topL) and ch!=topL[count-len(topL)]:
# #             botL.append(' ')
# #             count+=1
# #         botL.append(ch)
# #         count+=1
#     print("TOP", topL)
#     print("BOT", botL)
    
''' main method below '''

# str = input("Enter a DNA sequence: ")
# print("The reverse transcribed sequence is:",revtrans(str))
# print("The translated sequence is:",translate(str))

processFasta()


# print(revtrans('GAGGAAAGTCCGGGCTC'))
# print(jigsaw())



