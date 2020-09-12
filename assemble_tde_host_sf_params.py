'''
MIT License

Copyright (c) 2020 Nathaniel Jacob Roth

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE
'''

import numpy as np
import matplotlib.pyplot as plt
import glob
import json
from math import exp


filename_prefix = '/Users/nathanielroth/Dropbox/popsyn/'
#filename_suffix = '_parasfh_maxp.json'
filename_suffix = '_parasfh.json'

filename_list = glob.glob(filename_prefix + '*' + filename_suffix)

tde_names = []
tde_mstars = []
tde_uminusr = []
tde_tage = []
tde_tau = []
for fn in filename_list:

    print(fn)
    stripped_name = fn[len(filename_prefix):] # remove prefix
    stripped_name = stripped_name[:-len(filename_suffix)] # remove suffix
    tde_names.append(stripped_name)


    with open(fn) as json_file:
        data = json.load(json_file)
        tde_uminusr.append(data['u-r.0'][0]) # [0] is what you want, other elements are errors
        tde_mstars.append(data['mass'][0]) 
        tde_tage.append(data['tage'][0])
        tde_tau.append(data['tau'][0]) 
                

f = open('all_tde_host_data_popsyn.txt', 'w')

f.write("%-24s %-16s %-16s %-16s %-16s %-16s %-16s\n" % ('name', 'log10(mstar)', 'u-r color', 't_age (Gyr)', 'tau (Gyr)', 'SFR (msun yr^-1)', 'sSFR (yr^-1)'))

for i in range(len(tde_names)):

    #age_univ_gyr = 13.7
    #t_begin = age_univ_gyr - tde_tage[i]
    ##t_begin = tde_tage[i]

    #K = pow(10.,tde_mstars[i])/tde_tau[i] / (1. - exp(-1. * (age_univ_gyr- t_begin)/tde_tau[i])) *1.e-9 # solar mass / yr

    #SFR = K *  exp(-1. * (age_univ_gyr - t_begin)/tde_tau[i]) # solar mass / yr
    #sSFR = SFR / pow(10.,tde_mstars[i])

    #f.write('%-24s %-16.6f %-16.6f %-16.6f %-16.6f %-16.6e %-16.6e\n' % (tde_names[i],tde_mstars[i],tde_uminusr[i],tde_tage[i],tde_tau[i],SFR, sSFR))

    
    K = pow(10.,tde_mstars[i])/tde_tau[i] /(1. - exp(-tde_tage[i]/tde_tau[i])) *1.e-9 # solar mass / yr

    SFR = K *   exp(-tde_tage[i]/tde_tau[i]) 
    sSFR = SFR / pow(10.,tde_mstars[i])

    f.write('%-24s %-16.6f %-16.6f %-16.6f %-16.6f %-16.6e %-16.6e\n' % (tde_names[i],tde_mstars[i],tde_uminusr[i],tde_tage[i],tde_tau[i],SFR, sSFR))

    


tage = 6.05
tau = 0.29
mstar = pow(10.,10.46)

K = mstar/tau * 1. /(1. - exp(-tage/tau)) *1.e-9 # solar mass / yr

SFR = K *  exp(-tage/tau) 
sSFR = SFR / mstar


print(SFR,sSFR)

