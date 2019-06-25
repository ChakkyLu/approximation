from collections import defaultdict
import time
import sys
sys.path.append("..")
from _base.Circuit import Circuit
from _base.Implement import Implement
class T:
    def __init__(self):
        pass
    def cons(self):
        return defaultdict(list)

    def cons2(self):
        return defaultdict(int)

def method1():
    fr = open("kaka.csv")
    lines = fr.readlines()
    fr2 = open("haha.csv")
    lines2 = fr2.readlines()
    dict1, dict2 = defaultdict(str), defaultdict(str)

    for l1,l2 in zip(lines, lines2):
        l11 = l1.split(",")
        l22 = l2.split(",")

        dict1[l11[0]] += l1
        dict2[l22[0]] += l2

    fw1 = open("k1.csv", "w")
    fw2 = open("k2.csv", "w")

    for key, value in dict1.items():
        fw1.write(value)
        fw2.write(dict2[key])

    fr.close()
    fr2.close()
    fw1.close()
    fw2.close()

def method2(filename):
    ti = T()
    fr = open("haha.csv")
    lines = fr.readlines()
    dict2 = defaultdict(list)

    for l1 in lines:
        l11 = l1.split(",")
        dict2[l11[0]].append(l1.strip())


    spec_cir = Circuit()
    spec_cir.readBlif(filename)
    spec_cir.GetEachNodePriOut()

    fw = open("heihei.csv", "w")

    for key, value in dict2.items():
        pos = list(spec_cir.nodesDict[key])
        # for key2, value2 in value.items():
        # pos = list( spec_cir.nodesDict[key2] | pos1 )
        pos.sort(key=lambda x: int(x[1:]) if x[1] !='0' else int(x[2]))
        value2 = value
        newstr = "".join([v[:-1]+","+",".join(pos)+"\n" for v in value2])
        fw.write(newstr)

    fw.close()

def method3(filename):
    ti = T()
    fr = open("haha.csv")
    lines = fr.readlines()

    f2 = open("m6_2G.csv")
    lines2 = f2.readlines()

    nodeValues = {}
    twoNodeVal = defaultdict(ti.cons)
    twoNodeKey = defaultdict(ti.cons)
    prev_out = {}
    two_out = {}

    spec_cir = Circuit()
    spec_cir.readBlif(filename)
    spec_cir.GetEachNodePriOut()

    fw = open("kakak2.csv", "w")

    for l1 in lines:
        l11 = l1.split(",")
        nodeValues[l11[0]] = float(l11[3])
        # prev_out[l11[0]] = spec_cir.nodesDict[l11[0]]

    for l2 in lines2:
        l22 = l2.split(",")
        twoNodeVal[l22[0]][l22[2]].append(l2.strip())

    for key, value in twoNodeVal.items():
        pos1 = spec_cir.nodesDict[key]
        for key2, value2 in value.items():
            pos = list( spec_cir.nodesDict[key2] | pos1 )
            pos.sort(key=lambda x: int(x[1:]) if x[1] !='0' else int(x[2]))
            newstr = "".join([v+","+str(int(float(v.split(",")[-2])> nodeValues[v.split(",")[0]]))+ ","+",".join(pos)+"\n" for v in value2])
            fw.write(newstr)

def method4(filename):
    ti = T()
    fr = open("heihei.csv")
    lines = fr.readlines()

    nodeVal = {}

    fw = open("kakak2.csv", "w")

    spec_cir = Circuit()
    spec_cir.readBlif(filename)
    spec_cir.GetEachNodePriOut()

    for l1 in lines:
        l11 = l1.split(",")
        if l11[0] not in nodeVal:
            nodeVal[l11[0]] = [l11[2], str(len(spec_cir.nodesDict[l11[0]])), l11[4]]

    for k,v in nodeVal.items():
        fw.write(k+","+",".join(v)+"\n")

    fw.close()

def method5(filename1, filename2):
    f = open(filename1)
    cur_time = time.time()
    nodeValues = {}

    BESTCHOICE = {"AND": 0, "OR": 1, "XOR": 2}

    for l1 in f.readlines():
        l11 = l1.split(",")
        if l11[0] not in nodeValues:
            nodeValues[l11[0]] = {"err": float(l11[1]), "pk": int(l11[3]), "level": int(l11[2])}


    spec_cir = Circuit()
    spec_cir.readBlif(filename2)
    # nodeValues = spec_cir.GetBasicPattern()

    impl_cir = spec_cir.copy()

    nodes = spec_cir.nodes

    target_err = 0.001

    cur_err = 0

    cur_nodes = []

    i = 0

    f = 0
    while cur_err < target_err:
        last_err = cur_err
        # print(cur_nodes, cur_err)
        m = 1
        cur_target = ""
        if cur_nodes == []:
            m = 1
            cur_target = ""
            for node in nodes:
                k = nodeValues[node.name]["err"]
                if k < m:
                    m = k
                    cur_target = node
            cur_nodes.append(cur_target.name)
            cur_err += m
            nodes.remove(cur_target)
        else:
            choices = []
            m = 1
            cur_target = ""
            for node in nodes:
                k = max([nodeValues[node.name]["err"]+nodeValues[cur_node]["err"] for cur_node in cur_nodes])
                if k <= m:
                    if cur_target:
                        m = k
                        cur_target = node
                    else:
                        if not cur_target:
                            m = k
                            cur_target = node
            cur_err += m
            nodes.remove(cur_target)
            cur_nodes.append(cur_target.name)

        if cur_target.name == "n140":
            impl_cir.ModifyCir(cur_target.name, 2)
        else:
            impl_cir.ModifyCir(cur_target.name, BESTCHOICE[cur_target.type])

        sys.stdout.write("\rProgress: %d / %s" % (i, cur_err))
        i += 1


        # impl_cir.simplify()
    cur_time = time.time() - cur_time

    print("Exsim..." + "...cos..." + str(cur_time) + "s")
    # Im = Implement()
    # Im.DoExSimCir(spec_cir, impl_cir)
    impl_cir.writeBlif("le4.blif")
    print(spec_cir.getSize())
    print(impl_cir.getSize())

def method6(filename1):
    spec_cir = Circuit()
    spec_cir.readBlif(filename1)
    ap_err = spec_cir.GetBasicPattern()

    fw = open(spec_cir.model+".csv", "w")

    for key,value in ap_err.items():
        fw.write(key+","+str(value["err"])+","+str(value["level"])+","+str(value["pk"])+"\n")

    fw.close()

def method7(filename1,filename2):
    spec_cir = Circuit()
    spec_cir.readBlif(filename1)
    impl_cir = Circuit()
    impl_cir.readBlif(filename2)
    # ap_err = spec_cir.GetBasicPattern()
    print(spec_cir.getSize())
    print(impl_cir.getSize())

def method8(filename1):
    fw = open("data.csv", "w")

    fr1 = open(filename1)

    dict = defaultdict(list)
    dict2 = defaultdict(str)

    for line in fr1.readlines():
        line1 = line.strip().split(",")
        dict[line1[0]].append(line1[4])
        dict2[line1[0]]+= line

    for key,value in dict.items():
        if max(value)!=min(value):
            continue
        fw.write(dict2[key])

    fw.close()

def method9(filename1):
    fw = open("kk.txt", "w")
    t = 3
    for i in range(2,40):
        fw.write('p'+'0'*(2-len(str(i))) + str(i) + " ")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(">>>>Please use as follows:")
        print("simplify.py filename")
    else:
        # method6(sys.argv[1])
        method9(sys.argv[1])
