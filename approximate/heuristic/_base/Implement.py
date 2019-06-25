from .Circuit import Circuit
from collections import defaultdict
import time
import sys
sys.path.append("..")
from sat.genSAT import p5
import os
from multiprocessing import Pool
import math

class Implement:
    def __init__(self):
        self.maps = {}
        pass

    def DoExSim(self, filename1, filename2):
        spec_cir = Circuit()
        spec_cir.readBlif(filename1)
        impl_cir = Circuit()
        impl_cir.readBlif(filename2)

        circuit_size = len(spec_cir.inputs)

        j = 0
        WM = 0
        WCAE = float("-inf")
        MAE = 0
        pattern = []
        o1, o2 = [], []

        for inputs in range(2**circuit_size):

            inputs = str(bin(inputs))[2:]
            inputs = '0'* (circuit_size - len(inputs)) + inputs

            outputs = spec_cir.simulate(inputs)
            outputs2 = impl_cir.simulate(inputs)
            if outputs != outputs2:
                WM += 1

            cur_WCAE = float(abs(int(outputs,2)-int(outputs2,2))/(2**len(spec_cir.outputs)))


            if cur_WCAE >= 0.000108:
                print(inputs)
                print(int(outputs,2))
                print(int(outputs2,2))
                pattern.append(inputs)
                o1.append(outputs)
                o2.append(outputs2)

            if len(pattern) > 5: break

            WCAE = max(WCAE, cur_WCAE)
            MAE += cur_WCAE
            sys.stdout.write("\rProgress: %d / %d" % (j+1, 2**circuit_size))
            j += 1

        MAE = float( MAE / (2**circuit_size) )

        print("DIFFERENT MINTERMS: %s " % str(WM))
        print("worst case absolute error rate : %s " % str(WCAE))
        print("mean absolute error rate  %s " % str(MAE))
        print("================")
        print(pattern, o1, o2)

    def SingleCore(self, spec_cir, input_patterns_s, input_patterns_e):
        circuit_size = len(spec_cir.inputs)

        for inputs in range(input_patterns_s, input_patterns_e):

            inputs = str(bin(inputs))[2:]
            inputs = '0'* (circuit_size - len(inputs)) + inputs

            outputs = spec_cir.simulate(inputs)

            self.maps[inputs] = outputs
        return self.maps

    def GetMulDict(self, filename):
        spec_cir = Circuit()
        spec_cir.readBlif(filename)
        circuit_size = len(spec_cir.inputs)

        kdict = {}

        for i in range(2**32):
            kdict[i] = "0"*1

    def HeuForBooth(self, filename):
        print("Please input status:")
        print(">>> 123 : generate nodemap for booth multipliers")
        print(">>> 1 : generate approximate booth multipliers")

        status = int(input())

        if status == 123:

            datapath = os.path.abspath(os.path.join(os.getcwd(), "..")) + "/database/"
            satpath = os.path.abspath(os.path.join(os.getcwd(), "..")) + "/sat/"

            spec_cir = Circuit()
            spec_cir.readBlif(filename)
            circuit_size = len(spec_cir.inputs)
            size = spec_cir.size

            nodesMap_name = []
            nodesMap_type = []

            type = 1

            if type == 1:
                model = "model"
                newfile = str(size)
            else:
                model = "newmodel"
                newfile = "wt"+str(size)

            model16 = open(datapath+model+str(size)+".blif").readlines()
            sub16 = open(datapath+"sub"+str(size)+".blif").readlines()

            sub16 = "".join(sub16)
            model16 = "".join(model16)

            content = spec_cir.getStrOfCir()
            spec_cir.GetEachNodePriOut()

            fwrite = open(newfile+"base.csv", 'w')

            fixed_file = "for"+str(size)+"c.blif"

            for i in range(len(spec_cir.nodes)):
                #const zero
                node = spec_cir.nodes[i].copy()
                spec_cir.nodes[i].badstatus = 0
                if node.name in spec_cir.outputs:
                    fwrite.write(node.name+","+str(float(2**int(node.name[1:])/(2**len(spec_cir.outputs))))+"\n")
                else:
                    if node.inputs[0].name in spec_cir.inputs and node.inputs[1].name in spec_cir.inputs and  node.inputs[0].name[0] == node.inputs[1].name[0]:
                        spec_cir.nodes[i].badstatus = 1
                        print(node.name)
                        continue
                    if node.inputs[1].badstatus == 1 and node.inputs[0].name not in spec_cir.inputs:
                        print(node.name)
                        spec_cir.nodes[i].badstatus = 1
                        continue
                    err = 2 ** (spec_cir.poInfo[node.name]-len(spec_cir.outputs))
                    f = 1
                    fw = open(satpath+fixed_file, 'w')
                    node.inputs = []
                    node.type = "ZERO"
                    fw.write(model16)
                    newcontent = content.copy()
                    newcontent[i] = node.singleMSG()
                    # print(newcontent[i])
                    fw.write("".join(newcontent))
                    fw.write(sub16)
                    fw.close()
                    os.system('cd ../sat; abc -c \" read '+fixed_file+'; sat; \" > log1')
                    flog = open(satpath+"log1").readlines()[-1]
                    # print(flog)
                    if "UNSATISFIABLE" in flog:
                        # nodesMap_name.append(node.name)
                        # nodesMap_type.append('ZERO')
                        # f = 0
                        # nline = node.name+","+"0"+","+str(err)+"\n"
                        print(node.name+","+"0")
                        fwrite.write(node.name+","+"0"+","+str(err)+"\n")
                        continue
                    #const ONE
                    node = spec_cir.nodes[i].copy()
                    # print(spec_cir.poInfo[node.name])
                    fw = open(satpath+fixed_file, 'w')
                    node.inputs = []
                    node.type = "ONE"
                    fw.write(model16)
                    newcontent = content.copy()
                    newcontent[i] = node.singleMSG()
                    # print(newcontent[i])
                    fw.write("".join(newcontent))
                    fw.write(sub16)
                    fw.close()
                    os.system('cd ../sat; abc -c \" read '+fixed_file+'; sat; \" > log2')
                    flog = open(satpath+"log2").readlines()[-1]
                    # print(flog)
                    if "UNSATISFIABLE" in flog:
                        # if f==1:
                        # nodesMap_name.append(node.name)
                        # nodesMap_type.append('ONE')
                        print(node.name+","+"1")
                        fwrite.write(node.name+","+"1"+","+str(err)+"\n")

                    # sys.stdout.write("\rProgress: %d / %d" % (i+1, len(spec_cir.nodes)))
                    # sys.stdout.flush()

            fwrite.close()

        if status == 1:
            spec_cir = Circuit()
            spec_cir.readBlif(filename)
            size = spec_cir.size

            basefile = str(size)+"base.csv"

            datapath = os.path.abspath(os.path.join(os.getcwd(), "..")) + "/database/"
            newADD = ""

            if "wt" in spec_cir.model:
                basefile = "wt" + basefile
                newADD += "wt"

            ap_path = os.path.abspath(os.path.join(os.getcwd(), "..")) + "/approximate_circuit/"

            fbase = open(basefile).readlines()


            content = spec_cir.getStrOfCir()
            spec_cir.GetEachNodePriOut()

            # print(spec_cir.outputsInfo["n1055"])

            nodeMap = [(fcontent.split(",")[0], fcontent.split(",")[1], float(fcontent.split(",")[2])) for fcontent in fbase]
            nodeMap = sorted(nodeMap, key=lambda x:x[2])


            print("Please input target error: ")

            actual_target = float(input())
            target_err = actual_target
            cur_err = 0
            selected = []
            nnwithvalue = {}


            model16 = open(datapath+newADD+"model"+str(size)+".blif").readlines()
            sub16 = open(datapath+newADD+"sub"+str(size)+".blif").readlines()
            newsub16 = open(datapath+"newsub"+str(size)+".blif").readlines()

            sub16 = "".join(sub16)
            model16 = "".join(model16)
            newsub16 = "".join(newsub16)

            where = math.ceil(math.log(actual_target*(2**len(spec_cir.outputs)),2))
            sat_clause = p5(where,len(spec_cir.outputs), "")

            cur_time = time.time()
            content = spec_cir.getStrOfCir()
            content = "".join(content)

            satpath = os.path.abspath(os.path.join(os.getcwd(), "..")) + "/sat/"

            for nmap in nodeMap:
                nnwithvalue[nmap[0]] = nmap[2]
                node = spec_cir.nodeMaps[nmap[0]].copy()
                if (node.inputs[0].name in spec_cir.inputs) or (node.inputs[0].name in selected and node.inputs[1].name in selected):
                     # or \
                    # (node.inputs[0].name in spec_cir.inputs and node.inputs[1].name in selected):
                    if (cur_err+nmap[2]) <= target_err:
                        # print(nmap[0])
                        if spec_cir.nodeMaps[nmap[0]].inputs[0].name in selected:
                            cur_err -= nnwithvalue[spec_cir.nodeMaps[nmap[0]].inputs[0].name]

                        else:
                            if spec_cir.nodeMaps[nmap[0]].inputs[1].name in selected:
                                cur_err -= nnwithvalue[spec_cir.nodeMaps[nmap[0]].inputs[1].name]

                        cur_err += nmap[2]
                        spec_cir.ModifyCir(nmap[0], int(nmap[1]))

                        selected.append(nmap[0])


                        # if nmap[0] == "n1216":
                        #     break
                        content = spec_cir.getStrOfCir()
                        content = "".join(content)

                        fw = open(satpath+"forsat.blif", 'w')

                        fw.write(model16+content+sub16)

                        fw.close()

                        os.system('cd ../sat; abc -c \" read forsat.blif; resyn2; sat; \" > log')
                        flog = open(satpath+"log").readlines()[-1]

                        # print(flog)
                        if "UNSATISFIABLE" in flog:
                            selected.append(nmap[0])
                            pass
                        else:
                            spec_cir.nodes[spec_cir.find_node(node.name)] = node.copy()
                            print(nmap[0], int(nmap[1]), cur_err)
                            cur_err -= nmap[2]




            # print(cur_err)
            print(">>>It takes %s \'s " % str(time.time()-cur_time))
            # print(len(selected))


            # spec_cir.writeBlif(ap_path+spec_cir.model.lower()+"ap.blif");
            content = spec_cir.getStrOfCir()
            content = "".join(content)

            satpath = os.path.abspath(os.path.join(os.getcwd(), "..")) + "/sat/"
            fw = open(satpath+"forsat.blif", 'w')

            fw.write(model16+content+newsub16+sat_clause+".end")

            fw.close()

            os.system('cd ../sat; abc -c \" read forsat.blif; resyn2; sat; \" > log')
            flog = open(satpath+"log").readlines()[-1]

            print(flog)
            if "UNSATISFIABLE" in flog:
                print("見事にpassしました、bravo!")
            else:
                print("残念ながら、失敗しました！まだやり直しようよ！")

            spec_cir.writeBlif(ap_path+spec_cir.model.lower()+"ap.blif")

    def DoExSimCir(self, spec_cir, impl_cir):
        circuit_size = len(spec_cir.inputs)

        WM = 0
        WCAE = float("-inf")
        MAE = 0
        j = 0

        for inputs in range(2**circuit_size):

            inputs = str(bin(inputs))[2:]
            inputs = '0'* (circuit_size - len(inputs)) + inputs

            outputs = spec_cir.simulate(inputs)
            outputs2 = impl_cir.simulate(inputs)
            if outputs != outputs2:
                WM += 1

            WCAE = max(WCAE, float(abs(int(outputs,2)-int(outputs2,2))/(2**len(spec_cir.outputs))))
            MAE += float(abs(int(outputs,2)-int(outputs2,2))/ (2**len(spec_cir.outputs)))
            sys.stdout.write("\rProgress: %d / %d" % (j+1, 2**circuit_size))
            j += 1
        MAE = float( MAE / (2**circuit_size) )

        print("DIFFERENT MINTERMS: %s " % str(WM))
        print("worst case absolute error rate : %s " % str(WCAE))
        print("mean absolute error rate  %s " % str(MAE))

    def SingleChange(self, filename1):
        fw = open("result.csv", "w")
        spec_cir = Circuit()
        spec_cir.readBlif(filename1)
        impl_cir = spec_cir.copy()

        circuit_size = len(spec_cir.inputs)

        n = len(spec_cir.nodes)

        def basic(type, node1):
            node = node1.copy()
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
            return node

        # impl_cirs = [spec_cir.copy()] * 6
        for i in range(n):
            node = spec_cir.nodes[i]
            impl_cirs = [spec_cir.copy()] * 6
            if node.inputs:
                for j in range(6):
                    impl_cir = spec_cir.copy()
                    impl_cir.nodes[i] = basic(j, impl_cir.nodes[i])
                    print(impl_cir.nodes[i].type, spec_cir.nodes[i].type)

                    WM = 0
                    WCAE = float("-inf")
                    MAE = 0

                    for inputs in range(2**circuit_size):

                        inputs = str(bin(inputs))[2:]
                        inputs = '0'* (circuit_size - len(inputs)) + inputs

                        outputs = spec_cir.simulate(inputs)
                        outputs2 = impl_cir.simulate(inputs)
                        if outputs != outputs2:
                            WM += 1
                        WCAE = max(WCAE, float(abs(int(outputs,2)-int(outputs2,2))/(2**len(spec_cir.outputs))))
                        MAE += float(abs(int(outputs,2)-int(outputs2,2))/ (2**len(spec_cir.outputs)))
                    MAE = float( MAE / (2**circuit_size) )

                    print("DIFFERENT MINTERMS: %s " % str(WM))
                    print("worst case absolute error rate : %s " % str(WCAE))
                    print("mean absolute error rate  %s " % str(MAE))

                    fw.write(node.name +"," + str(j) + "," + str(WM) + "," + str(WCAE) + "," + str(MAE) )

    def CheckPattern(self, filename1, filename2, pp):
        spec_cir = Circuit()
        spec_cir.readBlif(filename1)
        impl_cir = Circuit()
        impl_cir.readBlif(filename2)

        circuit_size = len(spec_cir.inputs)

        WM = 0
        WCAE = float("-inf")
        MAE = 0

        for inputs in range(2**circuit_size):

            inputs = str(bin(inputs))[2:]
            inputs = '0'* (circuit_size - len(inputs)) + inputs

            outputs = spec_cir.simulate(inputs)
            outputs2 = impl_cir.simulate(inputs)

            if str(dict[pp[0]])+str(dict[pp[1]]) == "00" and str(dict2[pp[0]])+str(dict2[pp[1]]) == "11":
                print("YES")

    def simplfy(self, filename1):
        spec_cir = Circuit()
        spec_cir.readBlif(filename1)
        spec_cir.AIGtoXOR()
        # spec_cir.simplifyXNOR()
        save_path = os.path.abspath(os.path.join(os.getcwd(), "..")) + "/database/" + spec_cir.model.lower()+"c.blif"
        spec_cir.writeBlif(save_path)

    def Heuristic(self, filename1, filename2):
        start_time = int(time.time())
        spec_cir = Circuit()
        spec_cir.readBlif(filename1)
        impl_cir = Circuit()
        impl_cir.readBlif(filename2)

        circuit_size = len(spec_cir.inputs)
        size = int(circuit_size/2)

        j = 0
        WM = 0
        WCAE = float("-inf")
        MAE = 0
        pattern = []
        o1, o2 = [], []
        val = 0
        #
        print("Simulation.....")
        status = 1

        if status == 0:

            inputs_s = [0]*size + [0]*size

            inputs_s[int(0)] = 1
            inputs_s[size+int(7)] = 1
            # inputs_s[size+int(1)] = 1
            inputs_s = inputs_s[:size][::-1] + inputs_s[size:][::-1]

            inputs_s = "".join( list(map(str,inputs_s)) )

            print(inputs_s)

            for inputs in ["1"*(size*2)]:
                # inputs += "01"+""
                # inputs = str(bin(inputs))[2:]
                # inputs = "0"*(size-1-len(inputs)) + inputs
                # # inputs = "0"*size + inputs
                # inputs = inputs[:size-2] + "1" + inputs[:size-1] + inputs + "1"


                outputs = spec_cir.simulate(inputs)

                outputs2 = impl_cir.simulate(inputs)



                # outputs_c = str(bin(int("0"*1+"1"*15,2) * int("1"*16,2)))[2:]
                # outputs_c = int(inputs[:8],2) + int(inputs[8:],2)

                # print(outputs,outputs2, str(bin(outputs_c)))
                # print(int(outputs,2),int(outputs2,2))

                # print(outputs,outputs2)


                cur_WCAE = float(abs(int(outputs,2)-int(outputs2,2))/(2**len(spec_cir.outputs)))
                WCAE = max(cur_WCAE, WCAE)


                # print(inputs[0]=="1")
                if cur_WCAE > 0.9:
                    break

            print("worst case asbsolute error rate : %s " % str(WCAE))

        if status == 1:

            tmp = 0
            for inputs in range(2**circuit_size):

                inputs = str(bin(inputs))[2:]

                inputs = '0'* (circuit_size - len(inputs)) + inputs

                # if inputs[:size] == "0"*size or inputs[size:] == "0"*size:
                #     continue

                outputs, dict1 = spec_cir.simulate(inputs)
                # outputs = str(bin(abs(int(inputs[:size],2) * int(inputs[size:],2))))[2:]
                # outputs = '0'*(circuit_size-len(outputs)) + outputs
                outputs2, dict2 = impl_cir.simulate(inputs)
                # print(int(inputs[:int(len(inputs)/2)],2), int(inputs[int(len(inputs)/2):],2), int(outputs, 2), int(outputs2, 2))

                if outputs != outputs2:
                    WM += 1
                    if  int(outputs,2)-int(outputs2,2) == 25:
                        print(inputs, int(outputs,2), int(outputs2,2), int(outputs,2)-int(outputs2,2))
                        print(outputs, outputs2)
                    # break

                cur_WCAE = float(abs(int(outputs[:size],2)-int(outputs2[:size],2))/(2**size))

                # if cur_WCAE > 0.9:
                #     print(inputs[:spec_cir.size],inputs[spec_cir.size:])
                # print(inputs, int(outputs[:size],2), abs(int(outputs2[:size],2)))
                # if int(outputs[:size],2) < int(outputs2[:size],2):
                #     val += 1
                #
                #     # if (inputs[5]== "1" and inputs[-3]=="1"):
                #         # print(inputs[5], inputs[-3])
                #         # print(inputs, int(outputs[:size],2), abs(int(outputs2[:size],2)))
                #         # tmp += 1
                #     # if outputs[:4] == outputs2[:4]:
                #         # print(inputs,outputs[4:],outputs2[4:])
                #     pattern.append(inputs)
                #     o1.append(outputs)
                #     o2.append(outputs2)


                WCAE = max(WCAE, cur_WCAE)
                MAE += cur_WCAE
                # sys.stdout.write("\rProgress: %d / %d" % (j+1, 2**circuit_size))
                j += 1

                MAE = float( MAE / (2**circuit_size) )

            print("DIFFERENT MINTERMS: %s " % str(WM))
            # print(val,tmp)
            print("worst case absolute error rate : %s " % str(WCAE))
            print("worst case absolute error %s " % str(WCAE*(2**circuit_size)))
            print("================")
            # print(len(pattern))
            # print(pattern,o1,o2)
            print("Program Cost: %s 's " % str( float(time.time() - start_time)) )

    def getApproCir(self, filename):
        cur_time = time.time()
        spec_cir = Circuit()
        spec_cir.readBlif(filename)

        ap_err = spec_cir.GetBasicPattern()


        BESTCHOICE = {"AND": 0, "OR": 1, "XOR": 2, "XNOR":4, "NOR": 0, "AND2_PN":0, "AND2_NP":0, "NAND2_PN":1, "NAND2_NP":1}
        impl_cir = spec_cir.copy()
        nodes = spec_cir.nodes

        target_err = 0.01
        cur_err = 0
        cur_nodes = []
        i = 0

        while i < 0.5*len(spec_cir.nodes) and cur_err<target_err:
            m = 1
            cur_target = ""
            if cur_nodes == []:
                m = 1
                cur_target = ""
                for node in nodes:
                    k = ap_err[node.name]["err"]
                    if k < m:
                        m = k
                        cur_target = node
                cur_nodes.append(cur_target.name)
                cur_err += m
                cur_max = cur_err
                nodes.remove(cur_target)
            else:
                choices = []
                m = 1
                cur_target = ""
                for node in nodes:
                    k = cur_max +  ap_err[node.name]["err"]
                    # k = max([ap_err[node.name]["err"]+ap_err[cur_node]["err"] for cur_node in cur_nodes])
                    if k <= m:
                        m = k
                        cur_target = node

                cur_err += m
                nodes.remove(cur_target)
                cur_nodes.append(cur_target.name)
                if ap_err[cur_target.name]["err"] > cur_max: cur_max = ap_err[cur_target.name]["err"]

            print(cur_target.name)
            impl_cir.ModifyCir(cur_target.name, BESTCHOICE[cur_target.type])

            # sys.stdout.write("\rProgress: %d / %s" % (i, cur_err))
            i += 1


        cur_time = time.time() - cur_time

        save_path = os.path.abspath(os.path.join(os.getcwd(), "..")) + "/approximate_circuit/" + spec_cir.model.lower()+"_ap.blif"
        print("\n"+"FINISHED! COST TIME: " + str(cur_time) + "s")
        print("WCAE estimation : " + str(cur_err/10))
        impl_cir.writeBlif(save_path)
        self.DoExSimCir(spec_cir, impl_cir)

    def getUnit(self, filename):
        spec_cir = Circuit()
        spec_cir.readBlif(filename)

        ap_err = spec_cir.GetBasicPattern()

        unit_path = os.path.abspath(os.path.join(os.getcwd(), "..")) + "/database/" + spec_cir.model.lower()+".csv"

        fw = open(unit_path, "w")

        for key,value in ap_err.items():
            fw.write(key+","+str(value["level"])+","+str(value["pk"])+","+str(value["err"])+"\n")

        fw.close()

    def newMethod(self, filename):
        cur_time = time.time()
        spec_cir = Circuit()
        spec_cir.readBlif(filename)
        circuit_size = int(len(spec_cir.inputs)/2)
        output_size = len(spec_cir.outputs)
        cur_err = 0

        target_err = 1e-8
        #1e-5

        answers = spec_cir.copy()
        impl_cir = spec_cir.copy()

        print("Start......")



        for node in spec_cir.nodes:
            inputs = [1]*circuit_size + [0]*circuit_size
            if node.inputs[0].name in spec_cir.inputs and node.inputs[1].name in spec_cir.inputs:
                if node.type == "AND":
                    inputs[int(node.inputs[0].name[1:])] = 1
                    inputs[circuit_size+int(node.inputs[1].name[1:])] = 1
                    inputs = inputs[:circuit_size][::-1] +inputs[circuit_size:][::-1]
                    real_input = "".join(list(map(str,inputs)))
                    # print(real_input)
                    impl_cir.ModifyCir(node.name, 0)
                    outputs = impl_cir.simulate(real_input)
                    outputs2 = int(real_input[:circuit_size],2) + int(real_input[circuit_size:],2)
                    this_err = (abs(int(outputs,2)-outputs2)/2**output_size)
                    # print(node.name, this_err, real_input)
                    if (this_err + cur_err) < target_err:
                        # print(outputs,outputs2)
                        # print(node.name, this_err, real_input)
                        answers.ModifyCir(node.name, 0)
                        cur_err += this_err
                        impl_cir.nodes[impl_cir.find_node(node.name)] = node
                if node.type == "XOR":
                    inputs[int(node.inputs[0].name[1:])] = 1
                    inputs[circuit_size+int(node.inputs[1].name[1:])] = 0
                    inputs = inputs[:circuit_size][::-1] +inputs[circuit_size:][::-1]
                    real_input = "".join(list(map(str,inputs)))
                    # print(real_input)
                    impl_cir.ModifyCir(node.name, 0)
                    outputs = impl_cir.simulate(real_input)
                    outputs2 = int(real_input[:circuit_size],2) + int(real_input[circuit_size:],2)
                    this_err = (abs(int(outputs,2)-outputs2)/2**output_size)
                    # print(node.name, this_err, real_input)
                    if (this_err + cur_err) < target_err:
                        # print(outputs,outputs2)
                        # print(node.name, this_err, real_input)
                        answers.ModifyCir(node.name, 0)
                        cur_err += this_err
                        impl_cir.nodes[impl_cir.find_node(node.name)] = node
                if node.type == "XNOR":
                    inputs[int(node.inputs[0].name[1:])] = 1
                    inputs[circuit_size+int(node.inputs[1].name[1:])] = 1
                    inputs = inputs[:circuit_size][::-1] +inputs[circuit_size:][::-1]
                    real_input = "".join(list(map(str,inputs)))
                    # print(real_input)
                    impl_cir.ModifyCir(node.name, 0)
                    outputs = impl_cir.simulate(real_input)
                    outputs2 = int(real_input[:circuit_size],2) + int(real_input[circuit_size:],2)
                    this_err = (abs(int(outputs,2)-outputs2)/2**output_size)
                    # print(node.name, this_err, real_input)
                    if (this_err + cur_err) < target_err:
                        # print(outputs,outputs2)
                        # print(node.name, this_err, real_input)
                        answers.ModifyCir(node.name, 0)
                        cur_err += this_err
                        impl_cir.nodes[impl_cir.find_node(node.name)] = node


        save_path = os.path.abspath(os.path.join(os.getcwd(), "..")) + "/approximate_circuit/" + spec_cir.model.lower()+"_ap.blif"
        answers.writeBlif(save_path)
