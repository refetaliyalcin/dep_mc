# -*- coding: utf-8 -*-
#Copyright (C) 2016 Timo Väisänen and University of Helsinki
#All rights reserved.
#The new BSD License is applied to this software, see LICENSE.txt
import os.path
import fill_definitions as ifd
import fill_reader as fir
import fill_output as fore
reload(ifd)
reload(fir)
reload(fore)

output_format = "write(stream,'(A30,20X,"
output_format2 = "write(stream,'(A30,8X,"
output_format3 = "write(stream,'(A30,2X,"

STR_REAL4='real'
STR_REAL8='real(kind=rk)'
STR_INTEGER='integer'
STR_COMPLEX='complex(kind=rk)'
STR_LOGICAL='logical'
STR_CHARACTER='character'
STR_ALLOCATABLE = 'allocatable'

typeDict = {}
typeDict["REAL4"] = "real" 
typeDict["REAL8"] = "real(kind=rk)" 
typeDict["INTEGER"] = "integer" 
typeDict["COMPLEX"] = "complex(kind=rk)" 
typeDict["LOGICAL"] = "logical"  
typeDict["CHARACTER"] = "character" 

outputDict = {}
outputDict["REAL4"] = output_format+"E12.5)') " 
outputDict["REAL8"] = output_format+"E12.5)') " 
outputDict["INTEGER"] = output_format+"I12)') " 
outputDict["COMPLEX"] =  output_format2+"E12.5,E12.5)') " 
outputDict["LOGICAL"] = output_format+"L12)') "  
outputDict["CHARACTER"] = output_format3+"A30)') " 

outputDict2 = {}
outputDict2["REAL4"] = "" 
outputDict2["REAL8"] = "" 
outputDict2["INTEGER"] = "" 
outputDict2["COMPLEX"] =  "" 
outputDict2["LOGICAL"] = ""  
outputDict2["CHARACTER"] = "adjustl(trim(" 

outputDict3 = {}
outputDict3["REAL4"] = "" 
outputDict3["REAL8"] = "" 
outputDict3["INTEGER"] = "" 
outputDict3["COMPLEX"] =  "" 
outputDict3["LOGICAL"] = ""  
outputDict3["CHARACTER"] = "))" 


tab = "    "
space_def = 57
space_def2 = 65

#Hold the input parameters in a one place
class var_object():
    type0 = None
    allocatable = None
    read_name = None
    name = None
    value = None
    MPI_share = False
    read_from_file = None
    desc = None
    output_order = 0
    show_in_output = False

    def __init__(self):
        pass

    def set_type(self,str0):
        if str0.strip() not in typeDict:
            print "key: "+str0.strip()+" not in type dictionary"
        self.type0 = str0.strip();

    def set_allocatable(self,str0):     
        try:
            self.allocatable = int(float(str0));
            if self.allocatable < 0:
                print "ERROR: valid allocatable params are >=0"
            return True
        except ValueError:
            print str0
            print "ERROR: parameter allocatable has to be numeric"
            return False
   

    def set_output_order(self,str0):     
        try:
            self.output_order = int(float(str0));
            if self.output_order < 0:
                print "ERROR: valid output_order params are >=0"
            return True
        except ValueError:
            print str0
            print "ERROR: parameter output_order has to be numeric"
            return False    

    def set_read_name(self,str0):
        self.read_name = str0.strip();

    def set_name(self,str0):
        self.name = str0.strip();

    def set_show_in_output(self,str0):
        if(str0.lower().strip()=='true'):
            self.show_in_output = True
        else:
            self.show_in_output = False

    def set_value(self,str0):
        self.value = str0.strip();

    def set_MPI_share(self,str0):
        if(str0.lower().strip()=='true'):
            self.MPI_share = True
        else:
            self.MPI_share = False

    def set_read_from_file(self,str0):
        if(str0.lower().strip()=='true'):
            self.read_from_file = True
        else:
            self.read_from_file = False

    def set_desc(self,str0):
        self.desc = str0.strip();

    #Generate object
    def generate_var_object(self,list0):
        self.set_type(list0[0])
        self.set_allocatable(list0[1])
        self.set_read_name(list0[2])
        self.set_name(list0[3])
        self.set_value(list0[4])
        self.set_MPI_share(list0[5])
        self.set_read_from_file(list0[6])
        self.set_desc(list0[7])
        self.set_output_order(list0[8])
        self.set_show_in_output(list0[9])


    #Generate definitions for dataBlock
    def generate_definition_str(self):
        str0 = ""
        str_end = ""
        str0 = tab+tab+typeDict[self.type0]
        if self.type0 != "CHARACTER":
            if self.allocatable == 0:
                str0 = str0 +",allocatable"
                str_end = "(:)"
            elif self.allocatable >1:
                str_end = "("+str(self.allocatable)+")"
        else:
            if self.allocatable == 0:
                str0 = str0 +",allocatable"
                str_end = "(:)"
            elif self.allocatable >1:
                str0 = str0 +"(len="+str(self.allocatable)+")"

        str0 = str0 + " :: " + self.name + str_end

        for i in range(len(str0),space_def):
            str0 = str0+" "
        str0 = str0 + "!"+self.desc
        return str0


    #Generate required line input reader set1
    def generate_input_reader_line1(self):
        str0 = ""
        a = 0
        if self.type0 != "CHARACTER":
            str0 = tab+tab+"call setup_value(\""
            str0 = str0+self.read_name+"\","+self.value+")"
            a = len(str0)
        else:
            str0 = tab + tab+ "write(tmp,\"(A)\") \""
            str0 = str0 + self.value + "\"\n"
            str1 = tab + tab+"call setup_value(\""+self.read_name+ "\",tmp,64)"
            str0 = str0 + str1
            a = len(str1)
        
        for i in range(a,space_def2):
            str0 = str0+" "
        str0 = str0 + "!"+self.desc
        return str0
    #Generate required line input reader set2
    def generate_input_reader_line2(self):
        str0 = ""
        str_end = ""
        str0 = tab+tab+"call get_value(\""+ self.read_name +"\",dB%"
        if self.type0 == "CHARACTER":
            str0 = str0+self.name+",tmp)"    
        else:
            str0 = str0+self.name+")"
        
        
        for i in range(len(str0),space_def2):
            str0 = str0+" "
        str0 = str0 + "!"+self.desc
        return str0


    #Generate output line
    def generate_output_line(self):
        str0 = tab+tab+outputDict[self.type0]
        str0 = str0+"\""+self.read_name+"\","+outputDict2[self.type0]
        str0 = str0 +" dB%" + self.name + outputDict3[self.type0]
        return str0

#Tool which helps to create input for RTEngine
class Generate_variables(object):

    var_object_list = []
    def_object_list = []
    def_non_mpi_object_list = []
    inp_list1 = []
    inp_list2 = []


    def __init__(self):
        pass

    #Read parameters for input file default
    def generate_input_file0(self):
        proto_def   = "proto_definitions.f90"
        proto_mpi   = "proto_mpi_tools.f90"
        proto_reader = "proto_input_reader.f90"
        proto_writer = "proto_output_writer.f90"
        output_def = "definitions.f90"
        output_mpi = "output_mpi.f90"
        output_reader = "input_reader.f90"
        output_writer = "output_writer.f90"
        var_list = "var_list_n.csv"
        self.generate_input_file(proto_def,proto_mpi,proto_reader,proto_writer,output_def,output_mpi,output_reader,output_writer,var_list)  

    #Read parameters for input file
    def generate_input_file(self,proto_def,proto_mpi,proto_reader,proto_writer,output_def,output_mpi,output_reader,output_writer,var_list):
        errorCount = 0
        if(not os.path.isfile(proto_def)):
            errorCount=errorCount+1
            print "proto file for definitions does not exist"
        if(not os.path.isfile(proto_mpi)):
            errorCount=errorCount+1
            print "proto file for mpi_tools does not exist"
        if(not os.path.isfile(proto_reader)):
            errorCount=errorCount+1
            print "proto file for input_reader does not exist"
        if(not os.path.isfile(proto_writer)):
            errorCount=errorCount+1
            print "proto file for proto_writer does not exist"
        if(not os.path.isfile(var_list)):
            errorCount=errorCount+1
            print "var_list does not exist"      
        if(errorCount>0):
            print "abort"
            return
        self.read_var_list(var_list)
        self.read_inp_list()
        gen = ifd.fill_definitions()
        gen.make_definitions(proto_def,self.def_object_list, self.def_non_mpi_object_list)
        gen.print_definitions(output_def)
        gen = fir.fill_reader()
        gen.make_reader(proto_reader,self.inp_list1,self.inp_list2)
        gen.print_definitions(output_reader)
        gen = fore.fill_output()
        list2 = self.make_sorted_list(self.var_object_list)  
        list3 = self.read_output_list(list2)
        #print list3
        gen.make_output(proto_writer,list3)
        gen.print_output(output_writer)

    #make sorted list
    def make_sorted_list(self,list0):
        sorted_list = []
        ind = 0
        for i in self.var_object_list:
            
            if(i.show_in_output==True):
                found = False
                for j in range(0,ind):
                    if(sorted_list[j].output_order < i.output_order):
                        found = True
                        sorted_list.insert(j,i)
                        ind=ind+1
                        break
                if(not found):
                    sorted_list.append(i)
                    ind=ind+1
        return sorted_list

    

    def read_output_list(self,list0):
        list3 = []
        for i in list0:
            str0 = i.generate_output_line()
            #print str0
            list3.append(str0)
        return list3

    #Read and parse variables from the file "var_list"
    def read_var_list(self,var_list):
        with open(var_list, 'r') as f:
            content = f.readlines() 
            for j in range(0,len(content)-1):
                x = content[j].strip().split(",")
                c = var_object()
                c.generate_var_object(x)
                self.var_object_list.append(c)
        for i in self.var_object_list:
            var_str = i.generate_definition_str()
            if(i.MPI_share==True):
                self.def_object_list.append(var_str)
            else:
                self.def_non_mpi_object_list.append(var_str)
        pass
        
    #Read variables from the file "var_list"
    def read_inp_list(self):
        for i in self.var_object_list:
            if i.read_from_file==True :
                var_str1 = i.generate_input_reader_line1()
                var_str2 = i.generate_input_reader_line2()
                self.inp_list1.append(var_str1)
                self.inp_list2.append(var_str2)
        pass







