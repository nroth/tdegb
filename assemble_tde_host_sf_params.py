import numpy as np
import matplotlib.pyplot as plt
import glob
import json


filename_prefix = '/Users/nathanielroth/Dropbox/popsyn/'
#filename_suffix = '_parasfh_maxp.json'
filename_suffix = '_parasfh.json'

filename_list = glob.glob(filename_prefix + '*' + filename_suffix)

tde_names = []
tde_mstars = []
tde_uminusr = []
for fn in filename_list:

    print(fn)
    stripped_name = fn[len(filename_prefix):] # remove prefix
    stripped_name = stripped_name[:-len(filename_suffix)] # remove suffix
    tde_names.append(stripped_name)


    with open(fn) as json_file:
        data = json.load(json_file)
        tde_uminusr.append(data['u-r.0'][0]) # [0] is mass, other elements are errors
        tde_mstars.append(data['mass'][0]) # [0] is mass, other elements are errors
                

f = open('all_tde_host_data_popsyn.txt', 'w')

f.write("%-24s %-16s %-16s\n" % ('name', 'log10(mstar)', 'u-r color'))

for i in range(len(tde_names)):

    f.write('%-24s %-16.6f %-16.6f\n' % (tde_names[i],tde_mstars[i],tde_uminusr[i]))
