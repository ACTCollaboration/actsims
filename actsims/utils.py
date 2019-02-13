import numpy as np
from pixell import enmap,enplot
from orphics import io


def plot(fname,imap,dg=8):
    img = enplot.plot(enmap.downgrade(imap,dg),grid=False)
    if fname is None: 
        enplot.show(img)
    else: 
        enplot.write(fname,img)
        print(io.bcolors.OKGREEN+"Saved high-res plot to", fname+io.bcolors.ENDC)
