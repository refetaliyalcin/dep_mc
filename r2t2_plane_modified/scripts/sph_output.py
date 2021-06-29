#Copyright (C) 2018 Timo Väisänen and University of Helsinki
#All rights reserved.
#The new BSD License is applied to this software, see LICENSE.txt
import numpy as np
import output_reader
import process_sph_geom

reader = output_reader.OutputReader()
dataprocessor = process_sph_geom.ProcessSphereDat()
dB = reader.read_data("rt.out","cb.out")
MRT,MBS = dataprocessor.process_data(dB)

data = dataprocessor.combine_data(MRT,MBS,dB.the,dB.theb,30.0)
np.savetxt("output.out",data)