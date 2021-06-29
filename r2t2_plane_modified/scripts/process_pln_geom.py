import numpy as np

class ProcessPlaneDat():

    def process_k(self,k):
        W=-1
        if(np.mod(k,2)==0): W=1
        k3 = int(np.ceil(k/2))
        return W,k3

    def process_data(self,dB):
        phi = np.arange(0,360,45)
        indx = [i for i, item in enumerate(dB.phi) if item in phi]
        indxb = [i for i, item in enumerate(dB.phib) if item in phi]
        if dB.phi[indx][0]!=0 or dB.phi[indx][1]!=45 or dB.phi[indx][2]!=90: print("ERROR in phi angles")
        if dB.phib[indxb][0] != 0 or dB.phib[indxb][1] != 45 or dB.phib[indxb][2] != 90: print("ERROR in phi angles")
        MRT=np.zeros((1,dB.nthe,dB.nphi,4,4))
        for i in range(0,dB.nthe):
            for j in range(0,dB.nphi):
                MRT[0,i,j,0,:] = MRT[0,i,j,0,:] + (dB.rt[0,i,j,:] + dB.rt[1,i,j,:])*0.5
                for k in dB.pol:
                    W,k3 = self.process_k(k)
                    MRT[0, i, j, k3, :] = MRT[0, i, j, k3, :] + W*dB.rt[int(k-1), i, j, :] * 0.5
                MRT[0, i,j, :, :] = MRT[0, i,j, :, :].transpose()

        MBS = np.zeros((1,3, dB.ntheb, dB.nphib,  4, 4))
        for i in range(0,dB.ntheb):
            for j in range(0,dB.nphib):
                for n in range(0,3):
                    MBS[0,n,i,j,0,:] = MBS[0,n,i,j,0,:] + (dB.cb[0,n,i,j,:] + dB.cb[1,n,i,j,:])*0.5
                    for k in dB.pol:
                        W,k3 = self.process_k(k)
                        MBS[0,n,i,j,k3,:] = MBS[0,n,i,j,k3,:] + W*dB.cb[int(k-1),n,i,j,:] * 0.5
                    MBS[0,n,i,j,:,:] = MBS[0,n, i,j,:,:].transpose()
        return MRT,MBS




    def combine_data(self,MRT,MBS,the,theb,phi,phib,limit,phiS):
        indx = [i for i, item in enumerate(the) if item > limit]
        indxb = [i for i, item in enumerate(theb) if item <= limit]
        indPhi= np.where(phi==phiS)
        indPhib= np.where(phib==phiS)

        array1a = np.concatenate((MRT[0,indx,indPhi,:,:], MBS[0, 2, indxb[::-1],indPhib, :, :]), axis=1)

        indPhi= np.where(phi==phiS+180)
        indPhib= np.where(phib==phiS+180)

        array1b = np.concatenate((MRT[0, indx, indPhi, :, :], MBS[0, 2, indxb[::-1], indPhib, :, :]), axis=1)

        array1b = array1b[0,::-1,:,:]
        array1a = array1a[0, :, :, :]


        array1 = np.concatenate((array1a,array1b),axis=0)

        array2a = np.concatenate((the[indx], theb[indxb[::-1]]), axis=0)
        array2a = 180 - array2a
        array2b = 180 + np.concatenate((the[indx], theb[indxb[::-1]]), axis=0)[::-1]
        array2 = np.concatenate((array2a,array2b),axis=0)





        data = np.zeros((len(array2),17))
        for i in range(0,len(array2)):
            data[i,0] = array2[i]
            data[i,1:] = np.resize(array1[i,:,:],16)
        return data

