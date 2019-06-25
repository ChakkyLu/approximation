def kd(i):
    size = 16
    if size <100:
        return str(i) if i>9 else '0'+str(i)
    else:
        if i<10:
            return '00'+str(i)
        if i<100:
            return '0'+str(i)
        return str(i)

def main():
    fw = open("buma.txt", 'w')
    size = 16
    newstr = ""



    def p1(newstr):
        start_c = "k"
        dict = {'in':'z10', 'carry':"c10", 'real': "l"}

        for i in range(size):
            newstr += ".names " + start_c + kd(i) + " " + dict["in"]+kd(i)+"\n" + "1 0\n"

        newstr += ".names " + dict["in"]+kd(0) +" " + dict["carry"]+kd(0)+"\n1 1\n"
        newstr += ".names " + dict["in"]+kd(0) +" " + dict["real"]+kd(0)+"\n1 0\n"

        for i in range(1,size):
            newstr += ".names " + dict["in"]+kd(i)+" " +  dict["carry"]+kd(i-1) + " " + dict["real"]+kd(i)+ "\n10 1\n01 1\n"
            newstr += ".names " + dict["in"]+kd(i)+" " +  dict["carry"]+kd(i-1) + " " + dict["carry"]+kd(i)+ "\n11 1\n"
        return newstr

    def p2(newstr):
        newstr += ".names " + "l" + " l".join([kd(i) for i in range(size)]) + " " + "l"+kd(size)+"\n" + "0"*size + " 0\n"
        newstr += ".names " + "m" + kd(size) + "\n0\n"
        fr = open('a'+str(size+1)+".blif")
        newstr += "".join(fr.readlines())
        return newstr

    def p3(newstr):
        start_c = "s"
        dict = {'in':'z20', 'carry':"c20", 'real': "r"}

        for i in range(size):
            newstr += ".names " + start_c + kd(i) + " " + dict["in"]+kd(i)+"\n" + "1 0\n"

        newstr += ".names " + dict["in"]+kd(0) +" " + dict["carry"]+kd(0)+"\n1 1\n"
        newstr += ".names " + dict["in"]+kd(0) +" " + dict["real"]+kd(0)+"\n1 0\n"

        for i in range(1,size):
            newstr += ".names " + dict["in"]+kd(i)+" " +  dict["carry"]+kd(i-1) + " " + dict["real"]+kd(i)+ "\n10 1\n01 1\n"
            newstr += ".names " + dict["in"]+kd(i)+" " +  dict["carry"]+kd(i-1) + " " + dict["carry"]+kd(i)+ "\n11 1\n"
        return newstr

    def p4(pos, newstr):
        size = 25
        a = size - pos
        newstr += ".names " + "s" + " s".join([kd(i) for i in range(pos,size)]) +" " + "h00"+"\n" + "0"*a+" 0\n"
        newstr += ".names " + "r" + " r".join([kd(i) for i in range(pos,size)]) +" " + "h01"+"\n" + "0"*a+" 0\n"
        newstr += ".names h00 h01 "+"s"+kd(size)+" out\n" + "1-0 1\n-11 1"
        return newstr


def p5(pos, size, newstr):
    # size = 16
    a = size - pos
    newstr += ".names " + "s" + " s".join([kd(i) for i in range(pos,size)]) +" " + "h00"+"\n" + "0"*a+" 0\n"
    newstr += ".names " + "r" + " r".join([kd(i) for i in range(pos,size)]) +" " + "h01"+"\n" + "0"*a+" 0\n"
    newstr += ".names h00 h01 "+"s"+kd(size)+" out\n" + "1-0 1\n-11 1\n"
    return newstr


if __name__ == "__main__":
    newstr = p5(20, 32, "")
    print(newstr)

# newstr = p1(newstr)
# newstr = p2(newstr)
# newstr = p3(newstr)
# newstr = p4(14,newstr)
# fw.write(newstr)
