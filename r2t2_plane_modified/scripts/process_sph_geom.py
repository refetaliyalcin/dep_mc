import numpy as np

class ProcessSphereDat():

    def process_data(self,dB):
        phi = np.arange(0,360,45)
        indx = [i for i, item in enumerate(dB.phi) if item in phi]
        indxb = [i for i, item in enumerate(dB.phib) if item in phi]
        if dB.phi[indx][0]!=0 or dB.phi[indx][1]!=45 or dB.phi[indx][2]!=90: print("ERROR in phi angles")
        if dB.phib[indxb][0] != 0 or dB.phib[indxb][1] != 45 or dB.phib[indxb][2] != 90: print("ERROR in phi angles")
        if(int(dB.pol[0])>3): print("ERROR, use Q+,Q-")
        MRT=np.zeros((1,dB.nthe,4,4))
        MBS = np.zeros((3, 1, dB.nthe, 4, 4))
        for i in range(0,dB.nthe):
            MRT[0,i,0,:] = MRT[0,i,0,:] + np.sum(dB.rt[0,i,indx[0::2],:],0)/4
            MRT[0, i, 1, :] = MRT[0, i, 1, :] + \
                              (-dB.rt[0, i, indx[0], :] - dB.rt[0, i, indx[4], :]
                               +dB.rt[0, i, indx[2], :]+dB.rt[0, i, indx[6], :])/ 4
            MRT[0, i, 2, :] = MRT[0, i, 2, :] + \
                              (-dB.rt[0, i, indx[3], :] - dB.rt[0, i, indx[7], :]
                               +dB.rt[0, i, indx[1], :]+dB.rt[0, i, indx[5], :])/ 4
            if(len(dB.pol)>=2):
                MRT[0, i, 3, 2] = MRT[0, i, 3, 2] + \
                              (dB.rt[1, i, indx[0], 2] + dB.rt[1, i, indx[4], 2])/2
                MRT[0, i, 3, 3] = np.sum(dB.rt[1, i, :, 3])/dB.nphi
            MRT[0, i, :, :]=MRT[0, i, :,:].transpose()

        MBS = np.zeros((1,3, dB.ntheb, 4, 4))
        for j in range(0,3):
            for i in range(0,dB.ntheb):
                MBS[0,j,i,0,:] = MBS[0,j,i,0,:] + np.sum(dB.cb[0,j,i,indxb[0::2],:],0)/4
                MBS[0,j, i, 1, :] = MBS[0,j, i, 1, :] + \
                                  (-dB.cb[0,j, i, indxb[0], :] - dB.cb[0,j, i, indxb[4], :]
                                   +dB.cb[0,j, i, indxb[2], :]+dB.cb[0,j, i, indxb[6], :])/ 4
                MBS[0,j, i, 2, :] = MBS[0,j, i, 2, :] + \
                                  (-dB.cb[0,j, i, indxb[3], :] - dB.cb[0,j, i, indxb[7], :]
                                   +dB.cb[0,j, i, indxb[1], :]+dB.cb[0,j, i, indxb[5], :])/ 4
                if(len(dB.pol)>=2):

                    MBS[0, j, i, 3, 2]=MBS[0, j, i, 3, 2]+ \
                                   (dB.cb[1,j,i, indxb[0], 2] + dB.cb[1,j, i, indxb[4], 2]) / 2
                    MBS[0, j, i, 3, 3] = np.sum(dB.cb[1, j, i, :, 3]) / dB.nphib
                MBS[0, j, i, :, :]=MBS[0, j, i, :,:].transpose()
        return MRT,MBS

    def combine_data(self,MRT,MBS,the,theb,limit):
        indx = [i for i, item in enumerate(the) if item > limit]
        indxb = [i for i, item in enumerate(theb) if item <= limit]
        array1 = np.concatenate((MRT[0,indx,:,:], MBS[0, 2, indxb[::-1], :, :]), axis=0)
        array2 = np.concatenate((the[indx], theb[indxb[::-1]]), axis=0)
        array2=180-array2
        data = np.zeros((len(array1),17))
        for i in range(0,len(array1)):
            data[i,0] = array2[i]
            data[i,1:] = np.resize(array1[i,:,:],16)
        return data
