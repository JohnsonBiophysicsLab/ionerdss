from angles import angles
import math

def mag(x):
    return math.sqrt(sum(i ** 2 for i in x))

def string2coordlist(input_string):
    # replace the continue spaces with space
    newstring = ""
    store_space = False
    for oneChar in input_string:
        if oneChar==" ":
            if store_space == False:
                pass
            else:
                newstring += oneChar
                store_space = False
        else:
            newstring += oneChar
            store_space = True
    input_string = newstring
    coord1 = float(input_string[:input_string.find(" ")])
    temp_string = input_string[input_string.find(" ")+1:]
    coord2 = float(temp_string[:temp_string.find(" ")])
    temp_string = temp_string[temp_string.find(" ")+1:]
    coord3 = float(temp_string)
    return [coord1, coord2, coord3]


# read in coords from file
inputpath = ("./coordR25d10.txt")
my_input = open(inputpath, 'r', encoding = 'utf-8')
string_list_input = my_input.readlines()
my_input.close()

COMA = string2coordlist(string_list_input[1])
NormA= string2coordlist(string_list_input[2])
COMB = string2coordlist(string_list_input[10])
NormB= string2coordlist(string_list_input[11])
COMC = string2coordlist(string_list_input[19])
NormC= string2coordlist(string_list_input[20])
COMD = string2coordlist(string_list_input[28])
NormD= string2coordlist(string_list_input[29])
COMP = string2coordlist(string_list_input[136])
NormP= string2coordlist(string_list_input[137])
COMQ = string2coordlist(string_list_input[145])
NormQ= string2coordlist(string_list_input[146])
IntAB = string2coordlist(string_list_input[4])
IntAC = string2coordlist(string_list_input[5])
IntAD = string2coordlist(string_list_input[6])
IntAP = string2coordlist(string_list_input[7])
IntAQ = string2coordlist(string_list_input[8])
IntB = string2coordlist(string_list_input[13])
IntC = string2coordlist(string_list_input[26])
IntD = string2coordlist(string_list_input[34])
IntP = string2coordlist(string_list_input[141])
IntQ = string2coordlist(string_list_input[149])
IntmBinding = NormA
IntrBinding = string2coordlist(string_list_input[3])

print("------------------------------------------------------------------------------")
inner_angle = angles(COMA,COMB,IntAB,IntB,NormA,NormB)
sigma = mag([IntAB[0] - IntB[0],IntAB[1] - IntB[1],IntAB[2] - IntB[2]])
print("Normal point for A is %.6f %.6f %.6f" %(NormA[0] - COMA[0], NormA[1] - COMA[1],NormA[2] - COMA[2]))
#print("Normal point for B is %.6f %.6f %.6f" %(NormB[0] - COMA[0], NormB[1] - COMA[1],NormB[2] - COMA[2]))
print("Angles for dimers:")
print("%.6f, %.6f, %.6f, %.6f, %.6f, %.6f" % (inner_angle[0], inner_angle[1], inner_angle[2], inner_angle[3], inner_angle[4],sigma))
inner_angle = angles(COMA,COMC,IntAC,IntC,NormA,NormC)
sigma = mag([IntAC[0] - IntC[0],IntAC[1] - IntC[1],IntAC[2] - IntC[2]])
#print("Normal point for A is %.6f %.6f %.6f" %(NormA[0] - COMA[0], NormA[1] - COMA[1],NormA[2] - COMA[2]))
#print("Normal point for C is %.6f %.6f %.6f" %(NormC[0] - COMA[0], NormC[1] - COMA[1],NormC[2] - COMA[2]))
print("Angles for trimer AC:")
print("%.6f, %.6f, %.6f, %.6f, %.6f, %.6f" % (inner_angle[0], inner_angle[1], inner_angle[2], inner_angle[3], inner_angle[4],sigma))
inner_angle = angles(COMA,COMD,IntAD,IntD,NormA,NormD)
sigma = mag([IntAD[0] - IntD[0],IntAD[1] - IntD[1],IntAD[2] - IntD[2]])
#print("Normal point for A is %.6f %.6f %.6f" %(NormA[0] - COMA[0], NormA[1] - COMA[1],NormA[2] - COMA[2]))
#print("Normal point for D is %.6f %.6f %.6f" %(NormD[0] - COMA[0], NormD[1] - COMA[1],NormD[2] - COMA[2]))
print("Angles for hexamer AD:")
print("%.6f, %.6f, %.6f, %.6f, %.6f, %.6f" % (inner_angle[0], inner_angle[1], inner_angle[2], inner_angle[3], inner_angle[4],sigma))
inner_angle = angles(COMA,COMP,IntAP,IntP,NormA,NormP)
sigma = mag([IntAP[0] - IntP[0],IntAP[1] - IntP[1],IntAP[2] - IntP[2]])
#print("Normal point for A is %.6f %.6f %.6f" %(NormA[0] - COMA[0], NormA[1] - COMA[1],NormA[2] - COMA[2]))
#print("Normal point for P is %.6f %.6f %.6f" %(NormP[0] - COMA[0], NormP[1] - COMA[1],NormP[2] - COMA[2]))
#print("Angles for hexamer AP:")
#print("%.6f, %.6f, %.6f, %.6f, %.6f, %.6f" % (inner_angle[0], inner_angle[1], inner_angle[2], inner_angle[3], inner_angle[4],sigma))
inner_angle = angles(COMA,COMQ,IntAQ,IntQ,NormA,NormQ)
sigma = mag([IntAQ[0] - IntQ[0],IntAQ[1] - IntQ[1],IntAQ[2] - IntQ[2]])
#print("Normal point for A is %.6f %.6f %.6f" %(NormA[0] - COMA[0], NormA[1] - COMA[1],NormA[2] - COMA[2]))
#print("Normal point for Q is %.6f %.6f %.6f" %(NormQ[0] - COMA[0], NormQ[1] - COMA[1],NormQ[2] - COMA[2]))
#print("Angles for trimer AQ:")
#print("%.6f, %.6f, %.6f, %.6f, %.6f, %.6f" % (inner_angle[0], inner_angle[1], inner_angle[2], inner_angle[3], inner_angle[4],sigma))
c_IntAB = [IntAB[0] - COMA[0],IntAB[1] - COMA[1],IntAB[2] - COMA[2]]
c_IntAC = [IntAC[0] - COMA[0],IntAC[1] - COMA[1],IntAC[2] - COMA[2]]
c_IntAD = [IntAD[0] - COMA[0],IntAD[1] - COMA[1],IntAD[2] - COMA[2]]
c_IntAP = [IntAP[0] - COMA[0],IntAP[1] - COMA[1],IntAP[2] - COMA[2]]
c_IntAQ = [IntAQ[0] - COMA[0],IntAQ[1] - COMA[1],IntAQ[2] - COMA[2]]
c_IntmBinding = [IntmBinding[0] - COMA[0],IntmBinding[1] - COMA[1],IntmBinding[2] - COMA[2]]
c_IntrBinding = [IntrBinding[0] - COMA[0],IntrBinding[1] - COMA[1],IntrBinding[2] - COMA[2]]
print("------------------------------------------------------------------------------")
print("COM %.6f, %.6f, %.6f" %(0.000000, 0.000000, 0.000000))
print("b %.6f %.6f %.6f" %(c_IntAB[0], c_IntAB[1],c_IntAB[2]))
print("c %.6f %.6f %.6f" %(c_IntAC[0], c_IntAC[1],c_IntAC[2]))
print("d %.6f %.6f %.6f" %(c_IntAD[0], c_IntAD[1],c_IntAD[2]))
print("p %.6f %.6f %.6f" %(c_IntAP[0], c_IntAP[1],c_IntAP[2]))
print("q %.6f %.6f %.6f" %(c_IntAQ[0], c_IntAQ[1],c_IntAQ[2]))
print("m %.6f %.6f %.6f" %(c_IntmBinding[0], c_IntmBinding[1],c_IntmBinding[2]))
print("r %.6f %.6f %.6f" %(c_IntrBinding[0], c_IntrBinding[1],c_IntrBinding[2]))
print("------------------------------------------------------------------------------")



