from .Node import Node
from collections import defaultdict
import sys

def sameInput(l1, l2):
    if len(l1) != len(l2): return False
    ll1 = [l1[0].name, l1[1].name]
    ll2 = [l2[0].name, l2[1].name]
    if ll1 == ll2 or ll1 == ll2[::-1]:
        return True
    else: return False

class Circuit:
    def __init__(self):
        self.nodes = []
        self.model = ""
        self.outputNodes = []
        self.nodeMaps = {}
        self.maps = {}
        self.size = 0

    def readBlif(self, filename):
        fr = open(filename)
        lines = fr.readlines()
        i = 0
        while ".model" not in lines[i]:
            i += 1
        self.model = lines[i].strip().split()[1]
        inputs = ""
        outputs = ""
        while ".inputs" not in lines[i]:
            i += 1
        j = i
        while ".outputs" not in lines[j]:
            j += 1
        k = j

        while ".names" not in lines[k]:
            k += 1
        inputs = " ".join([lines[m].strip().replace("\\", "") for m in range(i,j)])
        self.inputs = inputs.split()[1:]
        # print(self.inputs)
        outputs = " ".join([lines[m].strip().replace("\\", "") for m in range(j,k)])
        self.outputs = outputs.split()[1:]
        while k < len(lines) - 1:
            line = lines[k]
            if '.end' in line:
                break
            if '.names' in line:
                if '.names' in lines[k+2] or '.end' in lines[k+2]:
                    nodes = line.split()[1:]
                    gates = lines[k+1].strip()
                    node = Node(nodes[-1], nodes[:-1], gates)
                    self.nodes.append(node)
                    k += 2
                else:
                    if '.names' in lines[k+3] or '.end' in lines[k+3]:
                        nodes = line.split()[1:]
                        gates = lines[k+1].strip() + '\n' + lines[k+2].strip()
                        node = Node(nodes[-1],nodes[:-1], gates)
                        self.nodes.append(node)
                        k += 3
                    else:
                        print("BLIF NODE NOT SUPPORT")
                        return 0
        # for node in self.nodes:
        #     print(node.name)
        self.size = int(len(self.inputs) / 2)
        self.genNodes()
        self.genLevel()
        self.genOutput()

    def getSize(self):
        m = 1
        for node in self.nodes:
            if node.inputs and len(node.inputs) == 2:
                m += 1
        return m

    def writeBlif(self, filename):
        fw = open(filename, 'w')
        fw.write('.model ' + self.model + '\n')

        fw.write('.inputs ' + ' '.join(self.inputs) + '\n')

        fw.write('.outputs ' + ' '.join(self.outputs)+ '\n')

        for node in self.nodes:

            if len(node.inputs) == 2:
                fw.write('.names ' + node.inputs[0].name + ' ' + node.inputs[1].name + ' ' + node.name)
            if len(node.inputs) == 1:
                fw.write('.names ' + node.inputs[0].name + ' ' + node.name)
            if len(node.inputs) == 0:
                fw.write('.names ' + node.name)

            fw.write('\n')
            fw.write(node.returnType())
            fw.write('\n')

        fw.write('.end')
        fw.close()

    def getStrOfCir(self):
        content = []

        for node in self.nodes:
            newStr = ""
            if len(node.inputs) == 2:
                i1 = node.inputs[0].name.replace("m","k").replace("n","j")
                i2 = node.inputs[1].name.replace("m","k").replace("n","j")
                o1 = node.name.replace("m","k").replace("n","j")
                newStr += '.names ' + i1 + ' ' + i2 + ' ' + o1
            if len(node.inputs) == 1:
                i1 = node.inputs[0].name.replace("m","k").replace("n","j")
                o1 = node.name.replace("m","k").replace("n","j")
                newStr += '.names ' + i1 + ' ' + o1
            if len(node.inputs) == 0:
                o1 = node.name.replace("m","k").replace("n","j")
                newStr += '.names ' + o1

            newStr += "\n"
            newStr += node.returnType()
            newStr += "\n"
            # tmp = [newStr]
            content.append(newStr)

        return content

    def InNodes(self, name):
        for node in self.nodes:
            if node.name == name:
                return node
        newNode = Node(name, None, None)
        return newNode

    def genNodes(self):
        n = len(self.nodes)
        for i in range(n):
            node = self.nodes[i]
            if node.name in self.inputs: self.nodes[i].is_input = True
            if node.name in self.outputs: self.nodes[i].is_output = True
            if node.inputs:
                for input in node.inputs:
                    tmp = self.InNodes(input)
                    self.nodes[i].updateInputs(tmp)
            self.nodeMaps[node.name] = self.nodes[i]

    def genLevel(self):
        self.nodes.sort(key=lambda x:int(x.name[1:].replace("_","")) if not x.is_output else 1000000+int(x.name[1:].replace("_","")))
        n = len(self.nodes)
        for i in range(n):
            self.nodeMaps[self.nodes[i].name] = self.nodes[i]
            if self.nodes[i].is_input:
                self.nodes[i].level = 0
                for j in range(len(self.nodes[i].outputs)):
                    self.nodes[i].outputs[j].level = 1
            else:
                if self.nodes[i].inputs:
                    self.nodes[i].level = max([node.level for node in self.nodes[i].inputs]) + 1
                else:
                    self.nodes[i].level = 0

    def sort(self):
        self.nodes.sort(key=lambda x:x.level)

    def genOutput(self):
        n = len(self.nodes)
        for i in range(n):
            node = self.nodes[i]
            if node.inputs:
                for input in node.inputs:
                    if input.name in self.inputs:
                        input.is_input = True
                    if not input.is_input:
                        self.nodes[self.find_node(input.name)].outputs.add(node)
            if node.name in self.outputs:
                self.outputNodes.append(node)

    def copy(self):
        new_cir = Circuit()
        new_cir.nodes = []
        for node in self.nodes:
            new_cir.nodes.append(node.copy())
        # new_cir.nodes = self.nodes.copy()

        new_cir.inputs = self.inputs
        new_cir.model = self.model
        new_cir.outputs = self.outputs
        return new_cir

    def changeType(self, name, node_type):
        nodes = self.nodes
        for i in range(len(nodes)):
            if nodes[i].name == name:
                # print(nodes[i].inputs[0].name, nodes[i].inputs[1].name)
                if node_type in ["BUF_0", "BUF_1", "NOT_0", "NOT_1"]:
                    if node_type[0:3] == "BUF":
                        replace_node = circuit.nodeExit(nodes[i].inputs[int(node_type[-1])].name)
                        for j in range(i+1,len(nodes)):
                            if nodes[j].inputs[0].name == name:
                                nodes[j].inputs[0] = replace_node
                            if nodes[j].inputs[1].name == name:
                                nodes[j].inputs[1] = replace_node
                    else:
                        replace_node = circuit.nodeExit(nodes[i].inputs[int(node_type[-1])].name)
                        for j in range(i+1,len(nodes)):
                            if nodes[j].inputs[0].name == name:
                                nodes[j].inputs[0] = replace_node
                                tmp = nodes[j].returnType()
                                tmp = bool(int(not int(tmp[0]))) + tmp[1:]
                                if nodes[j].type == "XOR" or nodes[j].type == "XNOR":
                                    nodes[j].type = nodes[j].inverse()
                                else:
                                    nodes[j].typeofGate(tmp)
                            if nodes[j].inputs[1].name == name:
                                nodes[j].inputs[1] = replace_node
                                tmp = nodes[j].returnType()
                                tmp = tmp[0] + bool(int(not int(tmp[1]))) + tmp[2:]
                                if nodes[j].type == "XOR" or nodes[j].type == "XNOR":
                                    nodes[j].type = nodes[j].inverse()
                                else:
                                    nodes[j].typeofGate(tmp)

                    nodes.pop(i)
                else:
                    node[i].type = node_type
                break
        self.nodes = nodes

    def nodeExit(self, name):
        for node in self.nodes:
            if node.name == name:
                return node
        return None

    def find_node(self, name):
        for i in range(len(self.nodes)):
            if self.nodes[i].name == name:
                return i
        print(name)
        print("NOT FOUND")
        return -1

    def all_nodes_map(self, index):
        all_nodes = []
        for node in self.nodes[index+1:]:
            for input in node.inputs:
                i = input.name
                if i not in self.inputs and i not in all_nodes:
                    all_nodes.append(i)
        return all_nodes

    def simplify(self):
        self.nodeMaps = {}

        i = len(self.nodes) - 1

        while i > -1 :
            if self.nodes[i].is_output:
                self.nodeMaps[self.nodes[i].name] = self.nodes[i]
                if self.nodes[i].inputs:
                    for node in self.nodes[i].inputs:
                        self.nodeMaps[node.name] = node
            else:
                if self.nodes[i].name in self.nodeMaps:
                    if self.nodes[i].inputs:
                        for node in self.nodes[i].inputs:
                            self.nodeMaps[node.name] = node
                else:
                    self.nodes.pop(i)
                    i += 1
            i -= 1
        # self.genOutput()
        self.sort()

    def findNodeByInput(self, inputs0, inputs1):
        inputs = [inputs0, inputs1]
        for i in range(len(self.nodes)):
            if len(self.nodes[i].inputs) == 2 and inputs == [self.nodes[i].inputs[0].name, self.nodes[i].inputs[1].name] or [self.nodes[i].inputs[0].name, self.nodes[i].inputs[1].name] == inputs[::-1]:
                return (i, self.nodes[i].name)
        return (-1,-1)

    def AIGtoXOR(self):
        print("======Turn To XOR=====")
        i = 0
        while i < len(self.nodes)-1:
            j = i + 1
            cur_node = self.nodes[i]
            while j < len(self.nodes):
                next_node = self.nodes[j]
                if sameInput(next_node.inputs, cur_node.inputs):
                    k = j + 1
                    for m in range(k,len(self.nodes)):
                        if [self.nodes[m].inputs[0].name, self.nodes[m].inputs[1].name] == [cur_node.name, next_node.name] \
                        or [self.nodes[m].inputs[1].name, self.nodes[m].inputs[0].name] == [cur_node.name, next_node.name]:

                            if (self.nodes[m].type == "NOR" and cur_node.type in ["AND2_NP", "AND2_PN"]) \
                            or (self.nodes[m].type == "OR" and cur_node.type in ["AND", "NOR"]):
                                self.nodes[m].inputs = self.nodes[i].inputs
                                self.nodes[m].type = "XNOR"
                            else:
                                if (self.nodes[m].type == "NOR" and cur_node.type in ["AND", "NOR"]) \
                                or (self.nodes[m].type == "OR" and cur_node.type in ["AND2_NP", "AND2_PN"]):
                                    self.nodes[m].inputs = self.nodes[i].inputs
                                    self.nodes[m].type = "XOR"
                                else:
                                    print(self.nodes[m].type, cur_node.type, cur_node.name)
                                    continue
                            break

                    break
                j += 1

            i += 1
        self.simplify()

    # def simulate(self, inputs):
    #
    #     value_dict = defaultdict(int)
    #
    #     circuit_size = len(self.inputs)
    #     single_input_size = int(circuit_size/2)
    #
    #     inputs_a = inputs[:single_input_size][::-1]
    #     inputs_b = inputs[single_input_size:][::-1]
    #     for i in range(int(len(self.inputs)/2)):
    #         value_dict[self.inputs[i]] = int(inputs_a[i])
    #         value_dict[self.inputs[i+single_input_size]] = int(inputs_b[i])
    #     i = 0
    #
    #     for node in self.nodes:
    #         if node.type == 'XOR':
    #             value_dict[node.name] = int(value_dict[node.inputs[0].name]!=value_dict[node.inputs[1].name])
    #             continue
    #         if node.type == 'XNOR':
    #             value_dict[node.name] = int(value_dict[node.inputs[0].name]==value_dict[node.inputs[1].name])
    #             continue
    #         if node.type == "ZERO":
    #             value_dict[node.name] = 0
    #             continue
    #         if node.type == "ONE":
    #             value_dict[node.name] = 1
    #             continue
    #         if node.type == "BUF":
    #             value_dict[node.name] = int(value_dict[node.inputs[0].name])
    #             continue
    #         if node.type == "NOT":
    #             value_dict[node.name] = int(not value_dict[node.inputs[0].name])
    #             continue
    #         type = node.returnType()
    #         type = type.replace(' ', '')
    #
    #
    #         if (value_dict[node.inputs[0].name] == int(type[0])) and (value_dict[node.inputs[1].name] == int(type[1])):
    #             value_dict[node.name] = int(type[-1])
    #         else:
    #             value_dict[node.name] = int(not int(type[-1]))
    #
    #     return ''.join([str(value_dict[output]) for output in self.outputs][::-1])

    def simulate(self, inputs):

        value_dict = defaultdict(int)

        circuit_size = len(self.inputs)
        single_input_size = int(circuit_size/2)

        inputs_a = inputs[:single_input_size][::-1]
        inputs_b = inputs[single_input_size:][::-1]
        for i in range(int(len(self.inputs)/2)):
            value_dict[self.inputs[i]] = int(inputs_a[i])
            value_dict[self.inputs[i+single_input_size]] = int(inputs_b[i])
        i = 0

        for node in self.nodes:
            if node.type == 'XOR':
                value_dict[node.name] = int(value_dict[node.inputs[0].name]!=value_dict[node.inputs[1].name])
                continue
            if node.type == 'XNOR':
                value_dict[node.name] = int(value_dict[node.inputs[0].name]==value_dict[node.inputs[1].name])
                continue
            if node.type == "ZERO":
                value_dict[node.name] = 0
                continue
            if node.type == "ONE":
                value_dict[node.name] = 1
                continue
            if node.type == "BUF":
                value_dict[node.name] = int(value_dict[node.inputs[0].name])
                continue
            if node.type == "NOT":
                value_dict[node.name] = int(not value_dict[node.inputs[0].name])
                continue
            type = node.returnType()
            type = type.replace(' ', '')


            if (value_dict[node.inputs[0].name] == int(type[0])) and (value_dict[node.inputs[1].name] == int(type[1])):
                value_dict[node.name] = int(type[-1])
            else:
                value_dict[node.name] = int(not int(type[-1]))

        return ''.join([str(value_dict[output]) for output in self.outputs][::-1]), value_dict

    def construct(self):
        return [0]*2

    def simulateEachNode(self, node_name="", node_vector=[], node_pattern=[], fixed_nodes=[], fixed_pattern=""):
        if node_name:
            target_node = self.nodeExit(node_name)
            if not target_node:
                print("Given Node does not exist")
                return 0
            else:
                new_node_vector = []
                self.simplify(target_node, new_node_vector)
                self.nodes = new_node_vector
        else:
            tota_dict = {}

            for node in self.nodes:
                if node.name not in self.inputs:
                    tota_dict[node.name] = [0]*2

            p_pro = 0
            n_pro = 0
            j = 0

            circuit_size = len(self.inputs)

            cater_patterns = [[] for _ in range(len(node_pattern))]


            for inputs in range(2**circuit_size):

                inputs = str(bin(inputs))[2:]
                inputs = '0'* (circuit_size - len(inputs)) + inputs

                outputs, dict = self.simulate(inputs)

                for key in tota_dict.keys():
                    tota_dict[key][dict[key]] += 1

                if node_vector and node_pattern:
                    for i in range(len(node_pattern)):
                        node_p = node_pattern[i]
                        tmp = "".join([str(dict[p]) for p in node_vector])
                        if tmp == node_p:
                            if fixed_nodes:
                                if "".join([str(dict[p]) for p in fixed_nodes]) == fixed_pattern:
                                    cater_patterns[i].append(inputs)
                            else:
                                cater_patterns[i].append(inputs)

                sys.stdout.write("\rProgress: %d / %d" % (j+1, 2**circuit_size))
                j += 1

            print("\n")

            if node_vector and node_pattern:
                fw = open("answer.csv", "w")
                fw.write("Required Node:    " + ",".join(node_vector) + "\n")
                for i in range(len(node_pattern)):
                    node_p = node_pattern[i]
                    fw.write("Required Pattern:    " + node_p + "\n")
                    fw.write("Total Account:    " + str(len(cater_patterns[i])) + "\n")

    def genNodeVectorPattern(self):
        self.nodeVP = []
        for node in self.nodes:
            if node.inputs[0].name not in self.inputs and node.inputs[1].name not in self.inputs:
                n1, n2 = node.inputs[0], node.inputs[1]
                if n1.inputs[0].name != n2.inputs[0].name and n1.inputs[0].name != n2.inputs[1].name:
                    tmp = [n1.inputs[0].name, n1.inputs[1].name, n2.inputs[0].name, n2.inputs[1].name]
                    out = node.name
                    aP = AnaPat(tmp, out)
                    self.nodeVP.append(aP)
                    nn1 = n1.inputs[0]

    def analyzeCircuit(self, node_name, node_vector=[]):
        if not node_name:
            target_node = self.nodeExit(node_name)
            if not target_node:
                print("Given Node does not exist")
                return 0
            else:
                new_node_vector = []
                self.simplify(target_node, new_node_vector)
                self.nodes = new_node_vector
        else:
            tota_dict = {}

            for node in self.nodes:
                if node.name not in self.inputs:
                    tota_dict[node.name] = [0]*2


            p_pro = 0
            n_pro = 0
            j = 0

            circuit_size = len(self.inputs)

            cater_patterns = defaultdict(list)


            for inputs in range(2**circuit_size):

                inputs = str(bin(inputs))[2:]
                inputs = '0'* (circuit_size - len(inputs)) + inputs

                outputs, dict = self.simulate(inputs)

                for key in tota_dict.keys():
                    tota_dict[key][dict[key]] += 1
                tmp = "".join([str(dict[p]) for p in node_vector])

                if tmp not in cater_patterns[dict[node_name]]:
                    cater_patterns[dict[node_name]].append(tmp)

                sys.stdout.write("\rProgress: %d / %d" % (j+1, 2**circuit_size))
                j += 1

            print("\n")

            fw = open("concrete.csv", "w")
            fw.write("Target Node: %s" % node_name + "\n")
            fw.write("Local Input Node: " +  " ".join(node_vector) + "\n")

            for key, value in cater_patterns.items():
                fw.write("Target Value: %s \n" % key)
                fw.write("\n".join(value) + "\n")

            fw.close()

    def simplify2(self, target_node, new_node_vector):
        if target_node:
            new_node_vector.insert(0,target_node)
            for node in target_node.inputs:
                self.simplify(self.nodeExit(node.name), new_node_vector)

    def GetValueOfSpecificNode(self, target_nodes, target_patterns, target_pin):

        circuit_size = len(self.inputs)
        k = [0]*2

        for inputs in range(2**circuit_size):
            inputs = str(bin(inputs))[2:]
            inputs = '0'* (circuit_size - len(inputs)) + inputs

            outputs, dict = self.simulate(inputs)


            tmp = "".join([str(dict[p]) for p in target_nodes])
            if tmp == target_patterns:
                k[int(dict[target_pin])] += 1

        print(k)

    def simplifyXNOR(self):
        n = len(self.nodes)
        for i in range(n):
            if self.nodes[i].type == "AND2_NP":
                self.nodes[i].type = "AND"
                self.nodes[self.find_node(self.nodes[i].inputs[0].name)].type = "OR"
                for node in self.nodes[i].inputs[0].outputs:
                    if node.type == "XNOR":
                        self.nodes[self.find_node(node.name)].type = "XOR"
                        break

    def genLocalSuspect(self, target_node, target, total, level):
        if target_node.name in self.inputs:
            if target_node.name[0] == "a":
                total[level][int(target_node.name[1:])] = target
            else:
                total[level][int(target_node.name[1:]) + int( len(self.inputs)/2)] = target
        else:
            pVal = self.getLocalInput(target_node, int(target))
            i = 0
            for p in pVal:
                m = target_node.inputs[0]
                n = target_node.inputs[1]
                self.genLocalSuspect(m, int(p[0]), total, level+1+i)
                self.genLocalSuspect(n, int(p[1]), total, level+1+i)
                i += 1

    def kkkk(self):
        return [0]*len(self.inputs)

    def callLocal(self):
        total = defaultdict(self.kkkk)
        target_node = self.nodeExit("n36")
        pVal = self.GetSpecDict("AND", "ZERO")
        for p in pVal:
            m = target_node.inputs[0]
            n = target_node.inputs[1]
            self.genLocalSuspect(m, int(p[0]), total, 0)
            self.genLocalSuspect(n, int(p[1]), total, 0)
        print(total)

    def kk(self):
        return defaultdict(list)

    def construct2(self):
        return [0]*15

    def construct3(self):
        return set([])

    def GetEachNodePriOut(self):
        n = len(self.nodes)
        self.outputsInfo = defaultdict(int)
        self.nodesDict = defaultdict(self.construct3)
        for node in self.outputNodes:
            m = node.name
            if node.inputs:
                for n2 in node.inputs:
                    self.nodesDict[n2.name].add(m)

        for i in range(n-1,-1,-1):
            if self.nodes[i].name not in self.outputs:
                node = self.nodes[i]
                if node.inputs:
                    for n2 in node.inputs:
                        self.nodesDict[n2.name] =  self.nodesDict[n2.name] | self.nodesDict[node.name]

        self.poInfo = {}

        for node in self.nodes:
            key = node.name
            self.outputsInfo[key] = len(self.nodesDict[key])
            # print(key)
            if not node.is_output:
                self.poInfo[key] = min([int(v[1:]) for v in self.nodesDict[key]])
            else:
                self.poInfo[key] = int(key[1:])


    def GetSpecDict(self,prev,cur):
        if prev == "AND":
            if cur=="ZERO":
                return ["11 1"]
            if cur=="ONE":
                return ["01 0", "10 0", "00 0"]
            if cur=="BUF0":
                return ["10 0"]
            if cur=="BUF1":
                return ["01 0"]
            if cur=="NOT0":
                return ["00 0", "11 1", "01 0"]
            if cur=="NOT1":
                return ["00 0", "11 1", "10 0"]
        if prev == "OR":
            #01 1 10 1 00 0 11 1
            if cur=="ZERO":
                return ["11 1", "10 1", "01 1"]
            if cur=="ONE":
                return ["00 0"]
            if cur=="BUF0":
                return ["01 0",]
            if cur=="BUF1":
                return ["10 0"]
            if cur=="NOT0":
                return ["10 1", "00 0", "11 1"]
            if cur=="NOT1":
                return ["00 0", "11 1", "01 0"]

        if prev == "OR":
            #01 1 10 1 00 0 11 1
            if cur=="ZERO":
                return ["11 1", "10 1", "01 1"]
            if cur=="ONE":
                return ["00 0"]
            if cur=="BUF0":
                return ["01 0",]
            if cur=="BUF1":
                return ["10 0"]
            if cur=="NOT0":
                return ["10 1", "00 0", "11 1"]
            if cur=="NOT1":
                return ["00 0", "11 1", "01 0"]

        if prev == "XOR":
            #10 1 01 1 11 0 00 0
            if cur=="ZERO":
                return ["10 1", "01 1"]
            if cur=="ONE":
                return ["00 0", "11 0"]
            if cur=="BUF0":
                return ["01 1", "11 0"]
            if cur=="BUF1":
                return ["10 1", "11 0"]
            if cur=="NOT0":
                return ["10 1", "00 0"]
            if cur=="NOT1":
                return ["01 1", "00 0"]

    def getLocalInput(self, target_node, target_value):
        if target_value:
            if target_node.type == "AND":
                return ["11"]
            if target_node.type == "OR":
                return ["01", "10", "11"]
            if target_node.type == "XOR":
                return ["10", "01"]
            if target_node.type == "XNOR":
                return ["00", "11"]
        else:
            if target_node.type == "AND":
                return ["01", "10", "00"]
            if target_node.type == "OR":
                return ["00"]
            if target_node.type == "XOR":
                return ["11", "00"]
            if target_node.type == "XNOR":
                return ["01", "10"]

    def ModifyCir(self, node_name, type):
        index = self.find_node(node_name)
        node = self.nodes[index]
        # print(node.type)
        if type == 0:
            while len(node.inputs):
                node.inputs.pop(0)
            node.type = "ZERO"
        if type == 1:
            while len(node.inputs):
                node.inputs.pop(0)
            node.type = "ONE"
        if type == 2:
            node.inputs.pop(1)
            node.type = "BUF"
        if type == 3:
            node.inputs.pop(0)
            node.type = "BUF"
        if type == 4:
            node.inputs.pop(1)
            node.type = "NOT"
        if type == 5:
            node.inputs.pop(0)
            node.type = "NOT"
        self.nodes[index] = node

    def GetBasicPattern(self):
        self.GetEachNodePriOut()
        approximate_error = {}
        for node in self.nodes:
            if node.is_output:
                approximate_error[node.name] = {"err": 1/ 2**(len(self.outputs)-int(node.name[1:])), "pk": self.outputsInfo[node.name], "level": node.level}
            else:
                approximate_error[node.name] = {"err": 1/ 2**self.outputsInfo[node.name], "pk": self.outputsInfo[node.name], "level": node.level}
        return approximate_error



if __name__ == "__main__":
    filename = "m6s.blif"
    circuit = Circuit()
    circuit.readBlif(filename)
    circuit.callLocal()
