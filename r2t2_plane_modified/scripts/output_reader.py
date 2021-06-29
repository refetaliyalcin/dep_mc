import os
import numpy as np

class DataBlock():
    def __init__(self):
        self.rawCB = None
        self.rawRT = None
        self.rt = None
        self.the = None
        self.phi = None
        self.theb = None
        self.phib = None
        self.pol = None
        self.cb = None
        self.rtp = None
        self.cbp = None
        self.nthe = 0
        self.nphi = 0
        self.ntheb = 0
        self.nphib = 0
        self.nrays = 0
        self.nrays = None
        self.npol = 0

class OutputReader():
    #Init
    def __init__(self):
        pass

    #read file data to list
    def __read_file_to_mem(self,file_name):
        lines = []
        #Check whether the file exist
        if(not os.path.isfile(file_name)):
            print("File \"" + file_name + "\"does not exist.")
            return False, None
        else:
            print("File "+file_name+" found, proceed..")
        #Open and read
        with open(file_name, 'r') as f:
            for line in f:
                sline = line.strip("\n")
                list0 = sline.split()
                if('#' not in sline):
                    lines.append(list0)
        print("Reading succeeded...")
        return True,lines


    def read_data(self,file_name_rt,file_name_cb):
        dB = DataBlock()
        dB = self.__read_file_rt(dB, file_name_rt)
        dB = self.__read_file_cb(dB,file_name_cb)
        print("Reading succeeded, return to the caller")
        return dB

    #read rt file
    def __read_file_rt(self, dB, file_name):
        success, dB.rawRT = self.__read_file_to_mem(file_name)
        if(not success): pass
        dB.nrays,vals,idx = self.__extract_nrays_vals(dB.rawRT,"nrays")
        dB.rt,dB.the,dB.phi,dB.pol = self.__parse_data(vals,dB.rawRT,idx,3,7,dB)
        dB.rt = self.__normalize_data_rt(dB.nrays,dB.rt)
        print("rt data ready")
        return dB


    #read cb file
    def __read_file_cb(self,dB,file_name):
        success, dB.rawCB = self.__read_file_to_mem(file_name)
        if (not success): pass
        dB.ntheb,dB.nphib = self.__extract_info_from_cb(dB.rawCB)
        self.__parse_data_cb(dB,dB.rawCB)
        dB.cb = self.__normalize_data_cb(dB.nrays, dB.cb)
        print("cb data ready")
        return dB


    def __parse_data_cb(self,dB,data):
        dB.cb = np.zeros((len(dB.pol),3,dB.ntheb,dB.nphib,4))
        dB.theb = np.zeros(dB.ntheb)
        dB.phib = np.zeros(dB.nphib)
        idx = 1
        for i in range(0,len(dB.pol)):
            for j in range(0,3):
                for k in range(0,dB.ntheb):
                    dB.theb[k] = data[idx][0]
                    for m in range(0,dB.nphib):
                        dB.phib[m] = data[idx][1]
                        dB.cb[i,j,k,m,:] = data[idx][4:8]
                        idx=idx+1




    #Normalize data so that it takes account that the number of rays
    #may have varied
    def __normalize_data_rt(self,nrays,data):
        for i in range(0,len(nrays)):
            data[i,:,:,:]=data[i,:,:,:]
        return data

    def __normalize_data_cb(self,nrays,data):
        for i in range(0,len(nrays)):
            data[i,:,:,:,:]=data[i,:,:,:,:]
        return data


    def __parse_data(self,vals,data,sidx,fromcol,tocol,dB):
        dB.npol = vals[0]
        dB.nthe = vals[1]
        dB.nphi = vals[2]
        arr = np.zeros((dB.npol,dB.nthe,dB.nphi,4))
        pols = np.zeros(dB.npol)
        the = np.zeros(dB.nthe)
        phi = np.zeros(dB.nphi)
        idx = sidx+2
        for i in range(0,dB.npol):
            pols[i] = data[idx][2]
            for j in range(0,dB.nthe):
                the[j] = data[idx][0]
                for k in range(0,dB.nphi):
                    phi[k] = data[idx][1]
                    arr[i,j,k] = data[idx][fromcol:tocol]
                    idx=idx+1
        return arr,the,phi,pols

    def __extract_info_from_cb(self,data):
        npol,tmp,nthe,nphi  = np.array(data[0][:]).astype(int)
        return nthe,nphi


    def __extract_nrays_vals(self,data,str):
        idx = self.__find_line_num_containing(data,str)
        nrays = np.array(data[idx][1:]).astype(int)
        vals = np.array(data[idx+1][:]).astype(int)
        return nrays, vals,idx

    #Find the index of the line which contains nrays
    def __find_line_num_containing(self,list,str):
        for idx, line in enumerate(list):
            if str in line:
                return idx
        return idx
