# -*- coding: utf-8 -*-
#Copyright (C) 2016 Timo Väisänen and University of Helsinki
#All rights reserved.
#The new BSD License is applied to this software, see LICENSE.txt

macroList1 = []
macroList2 = []
macroList1.append("!#INPUT_READER_SCRIPT1")
macroList1.append("!#INPUT_READER_SCRIPT2")
macroList2.append("!#INPUT_READER_SCRIPT3")
macroList2.append("!#INPUT_READER_SCRIPT4")
rowList1 = []
rowList2 = []
tab = "    "
class fill_reader():

    lines = []

    def __init__(self):
        pass

    #generate lines
    def make_reader(self,proto_file,list1,list2):
        #read everything from the proto_file
        i = 0
        with open(proto_file, 'r+') as myFile:
            for line in myFile:
                i=i+1
                self.lines.append(line.strip("\n"))
                if line.strip() in macroList1:
                    rowList1.append(i)
                if line.strip() in macroList2:
                    rowList2.append(i)
        #remove lines which are not needed
        self.removeLines(self.lines,rowList2[0],rowList2[1])
        self.removeLines(self.lines,rowList1[0],rowList1[1])
        
        i = rowList1[0]-1
        i=i+1
        str0 = tab+tab+"call init_input("+str(len(list1)+5)+")"
        self.lines.insert(i,str0)
        for j in list1:
            i=i+1
            self.lines.insert(i,j)
        
        ind = 0
        for i in self.lines:
            ind = ind+1
            if macroList2[0] in i:
                break
        i = ind-1
        for j in list2:
            i=i+1
            self.lines.insert(i,j)                    

    #Remove lines between script text
    def removeLines(self,list0,start,end):
        for i in range(end-2,start-1,-1):
            del list0[i]

    #save new input reader
    def print_definitions(self,output_file):
        with open(output_file, 'w') as myFile:
            for item in self.lines:
                myFile.write("%s\n" % item)
        
