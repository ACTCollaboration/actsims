import os,sys
import numpy as np

def mkdir(dirpath,comm=None):
    if comm is None:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
    exists = os.path.exists(dirpath)
    comm.Barrier()
    if comm.Get_rank()==0: 
        if not (exists):
            os.makedirs(dirpath)
    return exists
