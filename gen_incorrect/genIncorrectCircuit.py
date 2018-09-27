import random
import sys
import math


basicGates = ['00 1', '11 1', '10 1', '01 1', '01 0 \n'+'10 0', '10 0', '01 0']
types = {0: 'm', 1: 'z', 2: 'z'}
global ctrl

def basicPattern(mode, n, k):
    if not mode:
        return 'a'+ str(n) if k else 'b'+ str(n)
    else:
        return 'a_'+ str(n) + '_' if k else 'b_'+ str(n) + '_'

def outputPattern(mode, n):
    if not mode:
        return 'm0'+ str(n)
    else:
        return 'z_' + str(n) + '_'

def readOriginal(filename, mode):
    f = open(filename)
    start = 0
    names = []
    gates = []
    outputs = {}
    maxNode = 0
    for line in f.readlines():
        if not start and line[1:5] == "name":
            start = 1
        if start:
            if line[1:5] == "name":
                names.append(line.strip().split(" "))
                if names[-1][-1][0] == types[mode]:
                    outputs[names[-1][-1]] = len(names) - 1
                if [int(n[1:]) for n in names[-1][1:] if n[0]=='n'] != []:
                    maxNode = max(max([int(n[1:]) for n in names[-1][1:] if n[0]=='n']), maxNode)
            else:
                gates.append(line.strip())
    f.close()
    return names, gates, outputs, maxNode


def genBackModel(size, initialN, binaryList, whichOut, mode):
    backModel = ""
    basicModel = ""
    j = 0
    initial2 = initialN
    initial1 = initialN
    for i in range(size):
        if i<size/2:
            backModel += '.names ' + basicPattern(mode,2*i, 1) + ' '+ basicPattern(mode,2*i+1, 1) + ' n'+str(initialN) + "\n"
            if ctrl == 'easy': backModel += binaryList[j:j+2] + " 1"+ "\n"
            else: backModel += basicGates[random.randint(0, 6)] + '\n'
        else:
            backModel += '.names ' + basicPattern(mode,2*(i-4), 0) + ' ' + basicPattern(mode,2*(i-4)+1, 0) + ' n'+str(initialN) + "\n"
            if ctrl == 'easy': backModel += binaryList[j:j+2] + " 1" + "\n"
            else: backModel += basicGates[random.randint(0, 6)] + '\n'
        j += 2
        if i==size-1:
            basicModel += '.names ' + 'n'+str(initial1-1) + ' '+ 'n' +str(initialN+7) + ' '+whichOut+ "\n"
            basicModel += "01 1" + "\n" + "10 1" +"\n"
        else:
            basicModel += '.names ' + 'n'+str(initial2) + ' n'+str(initial2+1) + ' n'+str(initialN+8) + "\n"
            basicModel += "11 1" + "\n"
        initialN += 1
        initial2 += 2
    return backModel + basicModel, initialN+7

def genExtra(size, whichOuts, outputs, names, gates, maxNode, mode):
    newline = ""
    for whichOut in whichOuts:
        newline += ".names "
        preGate = names[outputs[whichOut]]
        for pre in preGate[1:-1]:
            newline += pre + " "
        newline += "n"+str(maxNode+1) + "\n"
        newline += gates[outputs[whichOut]] + "\n"
        k = random.randrange(0, pow(2, 2*size-1))
        k = str(bin(k))[2:]
        k = '0'*(2*size-len(k)) + k
        maxNode += 2
        lines, maxNode = genBackModel(size, maxNode, k, whichOut, mode)
        maxNode -= 1
        newline += lines
    return newline + ".end"

def writeNew(filename, newline, whichOuts):
    fr = open(filename)
    name = filename[0:filename.index('.')]
    fw = open(name+ 'in' + ".blif", 'w')
    stop = 1
    lines = fr.readlines()
    i = 0
    m = 0
    while i < len(lines):
        line = lines[i]
        if line[0:6]!= '.names' and not m:
            fw.write(line)
            i += 1
        else:
            m = 1
        if m:
            b = 1
            if line[0:3]==".en":
                break
            if len(line) > 7 and line[0:2]=='.n':
                for whichOut in whichOuts:
                    if whichOut in line.split()[-1]:
                        i += 2
                        b = 0
            if b :
                fw.write(line)
                i += 1
    fw.write(newline)
    fw.close()
    fr.close()


if __name__ == "__main__":
    # process the input
    filename = sys.argv[1]
    ctrl = sys.argv[2]
    s = 0
    if not filename or not ctrl:
        print("Please use arguments as follows:")
        print(".py(execution name)  .blif(multiplier.filename)  easy or hard (which mode)")
        s = 1
    if ctrl not in ['easy', 'hard']:
        print("Please use mode from the two kinds: easy or hard")
        s = 1
    type = ''
    size = ''
    for i in filename:
        if i.isalpha():type += i
        if i.isdigit():size += i
        if i == '.': break
    size = int(size)
    if type not in ['multi', 'mas', 'mont']:
        print("Please input correct type of multiplier!")
        s = 1
    if 2 ** math.log(size, 2) != size:
        print("Please input multiplier to the base of 2!")
        s = 1

    if not s:
        if type=='multi': mode = 0
        if type=='mas': mode = 1
        if type=='mont': mode = 2
        seed = random.randint(1,size)
        whichOut = []
        for _ in range(seed):
            n = random.randint(0,size-1)
            while outputPattern(mode,n) in whichOut:
                n = random.randint(0,size-1)
            whichOut.append(outputPattern(mode,n))

        names, gates, outputs, maxNode = readOriginal(filename, mode)
        newline = genExtra(size, whichOut, outputs, names, gates, maxNode, mode)
        writeNew(filename, newline, whichOut)
    else:
        pass
