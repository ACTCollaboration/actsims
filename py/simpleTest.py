
import simTools


#patch can be currently any of:




patch = 'deep8_s15_pa1_f150'



noiseSim = simTools.getActpolSim(iterationNum = 0, \
                                      simType = 'noise', \
                                      patch = patch)

cmbSim = simTools.getActpolSim(iterationNum = 0, \
                                      simType = 'cmb', \
                                      patch = patch)
