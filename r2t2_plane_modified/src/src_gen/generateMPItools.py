# -*- coding: utf-8 -*-
#Copyright (C) 2016 Timo Väisänen and University of Helsinki
#All rights reserved.
#The new BSD License is applied to this software, see LICENSE.txt

numList = []
macroList = []
macroList.append("!#MPITOOLS_SCRIPT001")
macroList.append("!#MPITOOLS_SCRIPT002")
macroList.append("!#MPITOOLS_SCRIPT101")
macroList.append("!#MPITOOLS_SCRIPT102")
macroList.append("!#MPITOOLS_SCRIPT201")
macroList.append("!#MPITOOLS_SCRIPT202")
macroList.append("!#MPITOOLS_SCRIPT301")
macroList.append("!#MPITOOLS_SCRIPT302")
macroList.append("!#MPITOOLS_SCRIPT401")
macroList.append("!#MPITOOLS_SCRIPT402")
macroList.append("!#MPITOOLS_SCRIPT501")
macroList.append("!#MPITOOLS_SCRIPT502")


proto_mpi_file = "proto_mpi_tools.f90"
mpi_output = "mpi_tools.f90"
definitions_file = "definitions.f90"



line1 = "!#MPITOOLS_SCRIPT1"
line2 = "!#MPITOOLS_SCRIPT2"



str_real="real(kind=rk)"
str_integer="integer"
str_complex="complex(kind=rk)"
str_real_all="real(kind=rk),allocatable"
str_logical="logical"


with open(proto_mpi_file) as myFile:
    for num, line in enumerate(myFile, 1):
        if line.strip() in macroList:
            numList.append(num)

def addType(typeList,UnpackList,AddressList,CopyArrayList,SubReadBufferList,line):
    valid = False
    alloc = False
    str1=""
    str2=""
    str3=""
    if str_real_all in line:
        valid = True
        alloc = True
        str1="SUBREADBUFFER(dB%"
        str2="COPYARRAY(dB%"
    elif str_real in line:
        valid = True
        str1="TYPELISTADDREAL"
        str2="MPIUNPACKREAL(dB%"
    elif str_complex in line:
        valid = True      
        str1="TYPELISTADDCOMPLEX"
        str2="MPIUNPACKCMPLX(dB%"
    elif str_logical in line:
        valid = True
        str1="TYPELISTADDLOGICAL"
        str2="MPIUNPACKLOGICAL(dB%"
    elif str_integer in line:
        valid = True
        str1="TYPELISTADDINTEGER"
        str2="MPIUNPACKINTEGER(dB%"
    if not valid:
        return
    str3="MPIADDRESS(dB%"
    b1, m, a = line.partition('::')
    b, m, a = a.partition('!')
    if not alloc:
        str1=str1
        str2=str2+b.replace(" ", "")+")"
        str3=str3+b.replace(" ", "")+")"
        typeList.append(str1)
        UnpackList.append(str2)
        AddressList.append(str3)
    else:
        b, m, a = b.partition('(:)')
        str1=str1+b.replace(" ", "")+")"
        str2=str2+b.replace(" ", "")+")"
        str3=str3+b.replace(" ", "")+")"      
        CopyArrayList.append(str1)
        SubReadBufferList.append(str2)


def removeLines(list0,start,end):
    for i in range(end-2,start-1,-1):
        #print lines[i]
        del list0[i]




record = False
typeList= []
UnpackList = []
AddressList = []
CopyArrayList = []
SubReadBufferList = []


ind=0
with open(definitions_file) as myFile:
    for line in myFile:
        line = line.rstrip() 
        if line2 in line:
            record = False 
        if record:
            addType(typeList,UnpackList,AddressList,CopyArrayList,SubReadBufferList,line)
        if line1 in line:
            record = True


lines = []


with open(proto_mpi_file, 'r+') as myFile:
    for line in myFile:
       lines.append(line.strip("\n"))


for i in range(len(macroList)-1,0,-2):
    print i,len(macroList)
    removeLines(lines,numList[i-1],numList[i])


ind1=lines.index("#define NUMOFARRAYSINDB 15")
lines[ind1]="#define NUMOFARRAYSINDB "+str(len(CopyArrayList))

ind1=lines.index("#define SIZEOFDATABLOCK 200")
lines[ind1]="#define SIZEOFDATABLOCK "+str(len(typeList)+2)


whitespace = "        "

def insertTo(whichList,lines,macroListInd,macroList,whitespace):
    ind1=lines.index(macroList[macroListInd])
    num=1
    for line in whichList:
        lines.insert(ind1+num,whitespace+line)
        num=num+1

insertTo(UnpackList,lines,0,macroList,whitespace)
insertTo(typeList,lines,2,macroList,whitespace)
ind3=lines.index(macroList[4])+1
str0="j1="+str(len(typeList))
#lines.insert(ind3,whitespace+str0)
#lines.insert(ind3,whitespace+"lengths(:)=1") 
ind4=lines.index(macroList[6])

insertTo(AddressList,lines,6,macroList,whitespace)
insertTo(CopyArrayList,lines,8,macroList,whitespace)
insertTo(SubReadBufferList,lines,10,macroList,whitespace)


with open(mpi_output, 'w') as myFile:
    for item in lines:
        myFile.write("%s\n" % item)

