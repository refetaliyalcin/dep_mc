# -*- coding: utf-8 -*-
#Copyright (C) 2016 Timo Väisänen and University of Helsinki
#All rights reserved.
#The new BSD License is applied to this software, see LICENSE.txt

macroList1 = []
macroList2 = []
macroList1.append("!#OUTPUT_GENERATOR_SCRIPT1")
macroList1.append("!#OUTPUT_GENERATOR_SCRIPT2")
rowList = []


class fill_output():

    lines = []

    def __init__(self):
        pass

    def make_output(self,proto_file,list1):
        #read everything from the proto_file
        i = 0
        with open(proto_file, 'r') as myFile:
            for line in myFile:
                i=i+1
                self.lines.append(line.strip("\n"))
                if line.strip() in macroList1:
                    rowList.append(i)
        #remove lines which are not needed
        self.removeLines(self.lines,rowList[0],rowList[1])

        i = rowList[0] 


        for j in list1:
            i=i
            self.lines.insert(i,j)
        i=i+1



    def removeLines(self,list0,start,end):
        for i in range(end-2,start-1,-1):
            del list0[i]

    def print_output(self,output_file):
        with open(output_file, 'w') as myFile:
            for item in self.lines:
                myFile.write("%s\n" % item)
        
