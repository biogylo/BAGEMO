#
#
#               BAGEMO : Basic Genetics Module
#
#   This is  a basic genetics module that includes some functions
#
#   CountBases(string):   Returns as a string the amount of A,C,G,T in a
#                           string of DNA as 4 integers separated by a
#                           space (' ') character.
#
#                    example:
#                           CountBases("AAAACCCGGT") -> "4 3 2 1"
#
#
import math
yes = ('S','s','Y','y','1','')
BaseComplements = {'A':'T','T':'A','G':'C','C':'G'} #The base pairs

def ncr(n, k):
    if 0 <= k <= n:
        ntok = 1
        ktok = 1
        for t in range(1, min(k, n - k) + 1):
            ntok *= n
            ktok *= t
            n -= 1
        return ntok // ktok
    else:
        return 0
def LoadCodons(): #Loads codons from file
    file = open("codons.txt",'r')
    data = file.read()
    file.close()
    data = list(data)
    dict  = {}
    for i in range(0,len(data),6):
        dict.update({(data[i]+data[i+1]+data[i+2]):data[i+4]})
    return dict
CodonDictionary = LoadCodons() #stores codons in a constant
def CountBases(string):
    # Gets a string of DNA and tells you, how much of each base is there
    a, c, g, t = 0, 0, 0, 0
    for char in string:
        if char == 'A':
            a+=1
        elif char == 'C':
            c+=1
        elif char == 'G':
            g+=1
        elif char == 'T' or char == 'U':
            t+=1
    return str(a)+" "+str(c)+" "+str(g)+" "+str(t)
def DNAtoRNA(string):
    #translates a string of DNA into a string of RNA (Tymine to Uracyl)
    string = list(string)
    for i in range(len(string)):
        if string[i] == 'T':
            string[i] = 'U'
    return "".join(string)
def ReverseComplement(string):
    string = list(string[::-1])
    for i in range(len(string)):
        string[i] = BaseComplements[string[i]]
    return "".join(string)
def RNAtoProtein(string):
    protein = []
    for i in range(0,len(string),3):
            protein.append(CodonDictionary[string[i:i+3:]])

    return "".join(protein)
def FindMotif(string,motif):
    subsize = len(motif)
    places = []
    print(range(len(string)-subsize+1))
    for i in range(len(string)-subsize+1):
        result = string.find(motif,i,i+subsize) + 1
        if result != 0:
            places.append(str(result)+' ')
    return "".join(places)
def MendelFirst(k,m,n):
    #for k homozygote dom , m heterozygote ,  n recessiv
    #return probability of having a dominant allele in offspring of any 2
    couples = ncr(k+m+n,2)
    doms = ncr(k,2) + k*m + k*n + m*n*0.5 + ncr(m,2)*0.75
    return doms/couples
def FibRab(n,k):
    if n == 1:
        return 1
    elif n == 2:
        return 1
    else:
        return FibRab(n-2,k)*k + FibRab(n-1,k)
def MortFibRab(n,m):
    if n > 1:
        return MortFibRab(n-2,m) + MortFibRab(n-1,m) - MortFibRab(n-m,m)
    elif n in [1,0]:
        return 1
    else:
        return 0
