class Node:

    def __init__(self, name, inputs, type="BLIF"):
        self.name = name
        self.inputs = inputs if inputs else None
        self.outputs = set([])
        self.level = 0
        self.badstatus = 0
        if type:
            if type != "BLIF":
                self.typeofGate(type)
            else:
                self.type = type
        else: self.type = "Primary Input"
        self.is_input = False
        self.is_output = False

    def copy(self):
        node = Node(self.name, None, "")
        node.inputs = self.inputs.copy()
        node.outputs = self.outputs
        node.type = self.type
        return node

    def updateInputs(self, node):
        if self.inputs:
            if isinstance(self.inputs[0], str):
                self.inputs = []
                self.inputs.append(node)
            else:
                if isinstance(self.inputs[0], Node):
                    self.inputs.append(node)
        else:
            self.inputs = []
            self.inputs.append(node)

    def typeofGate(self, type):
        self.type = ""
        if type == "0":
            self.type = "ZERO"
        if type == "1":
            self.type = "ONE"
        if type == "1 1":
            self.type = "BUF"
        if type == "1 0":
            self.type = "NOT"
        if type == "11 1":
            self.type = "AND"
        if type == "11 0":
            self.type = "NAND"
        if type == "00 0":
            self.type = "OR"
        if type == "00 1":
            self.type = "NOR"
        if type == "10 1\n01 1" or type == "01 1\n10 1" or type=="11 0\n00 0" or type=="00 0\n11 0":
            self.type = "XOR"
        if type == "00 1\n11 1" or type == "11 1\n00 1" or type == "10 0\n01 0" or type == "01 0\n10 0":
            self.type = "XNOR"
        if type == "01 1":
            self.type = "AND2_NP"
        if type == "01 0":
            self.type = "NAND2_NP"
        if type == "10 1":
            self.type = "AND2_PN"
        if type == "10 0":
            self.type = "NAND2_PN"
        if not self.type:
            print(type)
            print("Gate Type Not Supported!")

    def returnType(self):
        if self.type == "AND":
            return "11 1"
        if self.type == "NAND":
            return "11 0"
        if self.type == "OR":
            return "00 0"
        if self.type == "NOR":
            return "00 1"
        if self.type == "XOR":
            return "10 1\n01 1"
        if self.type == "XNOR":
            return "00 1\n11 1"
        if self.type == "AND2_NP":
            return "01 1"
        if self.type == "NAND2_NP":
            return "01 0"
        if self.type == "AND2_PN":
            return "10 1"
        if self.type == "NAND2_PN":
            return "10 0"
        if self.type == "BUF":
            return "1 1"
        if self.type == "NOT":
            return "0 1"
        if self.type == "ONE":
            return "1"
        if self.type == "ZERO":
            return "0"
        return "Gate Type Not Exist"

    def inverse(self):
        if self.type == "ZERO":
            return "ONE"
        if self.type == "ONE":
            return "ZERO"
        if self.type == "BUF":
            return "NOT"
        if self.type == "NOT":
            return "BUF"
        if self.type == "AND":
            return "NAND"
        if self.type == "NAND":
            return "AND"
        if self.type == "OR":
            return "NOR"
        if self.type == "NOR":
            return "OR"
        if self.type == "XOR":
            return "XNOR"
        if self.type == "XNOR":
            return "XOR"
        if self.type == "AND2_NP":
            return "NAND2_NP"
        if self.type == "NAND2_NP":
            return "AND2_NP"
        if self.type == "AND2_PN":
            return "NAND2_PN"
        if self.type == "NAND2_PN":
            return "AND2_PN"
        return "CANNOT FIND THE INVERTER"

    def singleMSG(self):
        node = self
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
        return newStr
