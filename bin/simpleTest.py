from __future__ import print_function
from actsims import simTools


#psa can be currently any of:
#                   [['deep1_s13_pa1_f150'] ,\
#                    ['deep5_s13_pa1_f150'],\
#                    ['deep6_s13_pa1_f150'],\
                   
#                    ['deep56_s14_pa1_f150'],\
#                    ['deep56_s14_pa2_f150'],\
#                    ['deep56_s15_pa1_f150'],\
#                    ['deep56_s15_pa2_f150'],\
#                    ['deep56_s15_pa3_f090'],\
#                    ['deep56_s15_pa3_f150'],\
#                    ['deep8_s15_pa1_f150'],\
#                    ['deep8_s15_pa2_f150'],\
#                    ['deep8_s15_pa3_f090'],\
#                    ['deep8_s15_pa3_f150']]



#Alternatively, you can pass the info separately, e.g: season = 's13', pa = 'pa1', freqGHz = 150, region = 'deep5'

patch = 'deep56'
season = 's15'
array = 'pa1'



noiseSim = simTools.getActpolSim(iterationNum = 0, \
                                      simType = 'noise', \
                                      patch = patch, \
                                 season = season, \
                                 array = array)

cmbSim = simTools.getActpolSim(iterationNum = 0, \
                                      simType = 'cmb', \
                                      patch = patch, \
                                 season = season, \
                                 array = array)

print("Noise sim shape is :",noiseSim.shape)
print("CMB sim shape is :",cmbSim.shape)
