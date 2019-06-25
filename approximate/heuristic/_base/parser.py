class parser:
    def __init__(self):
        pass

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


    def genLocalSuspect(self, target_node, target, total, level):
        if target_node.name in self.inputs:
            total[level][int(target_node.name[1:])] = target
        pval = self.GetSpecDict("BUF0")
        for p in pVal:
            m = target_node.inputs[0]
            n = target_node.inputs[1]
            self.genLocalSuspect(m, p[0], total, level+1)
            self.genLocalSuspect(n, p[1], total, level+1)
