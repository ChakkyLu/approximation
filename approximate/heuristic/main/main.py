import sys
sys.path.append("..")
from _base.Circuit import Circuit
from _base.Implement import Implement


class Solution:
    def __init__(self):
        self.Im = Implement()
        pass

    def ex_sim(self, filename1, filename2):
        try:
            open(filename1)
            open(filename2)
        except Exception as e:
            print(">>>>>>Can not open files!")
            return 0

        self.Im.DoExSim(filename1, filename2)

    def eq_single(self, filename1):
        try:
            open(filename1)
        except Exception as e:
            print(">>>>>>Can not open files!")
            return 0

        self.Im.SingleChange(filename1)

    def an_eachnode(self, filename1, filename2):
        try:
            open(filename1)
        except Exception as e:
            print(">>>>>>Can not open files!")
            return 0

        self.Im.CheckPattern(filename1, filename2, ["m05", "m06"])

    def change(self, filename1):
        try:
            open(filename1)
        except Exception as e:
            print(">>>>>>Can not open files!")
            return 0

        fr = open(filename1)
        fw = open("nn.blif","w")

        import re
        for line in fr.readlines():
            line = line.replace("P[","m")
            line = line.replace("IN1","a").replace("[","").replace("]","")
            line = line.replace("IN2","b")
            fw.write(line)

        fw.close()





if __name__ == "__main__":
    C = Solution()
    Im = Implement()

    if len(sys.argv) <= 2:
        print(">>>>Please use as follows:")
        print("simplify.py filename")
    else:
        command = sys.argv[1]

        if command == "exsim":
            if len(sys.argv) != 3:
                print(">>>>Please use as follows")
                print("simeq spec.blif impl.blif")
            else:
                print("=========")
                Im.GetMulDict(sys.argv[2])

        if command == "simeq":
            if len(sys.argv) != 4:
                print(">>>>Please use as follows")
                print("simeq spec.blif impl.blif")
            else:
                print("=========")
                C.ex_sim(sys.argv[2], sys.argv[3])
        if command == "ecsingle":
            if len(sys.argv) != 3:
                print(">>>>Please use as follows")
                print("ecsingle spec.blif")
            else:
                print("=========")
                C.eq_single(sys.argv[2])
        if command == "heu":
            if len(sys.argv) != 4:
                print(">>>>Please use as follows")
                print("ecsingle spec.blif")
            else:
                print("=========")
                Im.Heuristic(sys.argv[2], sys.argv[3])
        if command == "simplify":
            if len(sys.argv) != 3:
                print(">>>>Please use as follows")
                print("ecsingle spec.blif")
            else:
                print("=========")
                Im.simplfy(sys.argv[2])
        if command == "xor":
            if len(sys.argv) != 3:
                print(">>>>Please use as follows")
                print("ecsingle spec.blif")
            else:
                print(">>>Turning original ciruit to XOR format........")
                Im.turnXOR(sys.argv[2])
        if command == "approximate":
            if len(sys.argv) != 3:
                print(">>>>Please use as follows")
                print("ecsingle spec.blif")
            else:
                print(">>>Generating approximate circuit............")
                Im.getApproCir(sys.argv[2])
        if command == "newm":
            if len(sys.argv) != 3:
                print(">>>>Please use as follows")
                print("ecsingle spec.blif")
            else:
                print(">>>>Generating approximate circuit........")
                Im.newMethod(sys.argv[2])
        if command == "changeType":
            if len(sys.argv) != 3:
                print(">>>>Please use as follows")
                print("ecsingle spec.blif")
            else:
                print(">>>>Generating the preliminary analysis........")
        if command == "booth":
            if len(sys.argv) != 3:
                print(">>>>Please use as follows")
                print("ecsingle spec.blif")
            else:
                print(">>>>Generating the preliminary analysis........")
                Im.HeuForBooth(sys.argv[2])
