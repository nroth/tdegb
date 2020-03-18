import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def plot_1d_histogram(prefix = './',hist_type = 'gals', axis = 'mstar_mendel'):
    filename = prefix + hist_type + '_' + axis + '_1d.hist'

    hist = np.loadtxt(filename)

    left_edges = hist[::,0]
    right_edges = hist[::,1]
    bin_centers = 0.5 * (left_edges + right_edges) # assumes uniform spacing
    widths = right_edges - left_edges
    counts = hist[::,2]

    print("Total counts is %d" %counts.sum())
    plt.bar(bin_centers,counts, width = widths, align='center',edgecolor='k',color='white')





def plot_2d_histogram(prefix ='./', hist_type = 'gals',axis1 = 'mstar_mendel', axis2 = 'UminusR', transpose_default = 0, cmap = "YlGnBu", cm_maxval = 1.0, cm_minval = 0., zlim = 0):

    switched_name = 0
    try:
        filename = "gals_" + axis1 + "_" + axis2 + "_2d.hist"
        hist_2d_counts = np.loadtxt(filename,skiprows=2,usecols=2) # fix the name, 2D
    except:
        try:
            filename = "gals_" + axis2 + "_" + axis1 + "_2d.hist"
            hist_2d_counts = np.loadtxt(filename,skiprows=2,usecols=2) # fix the name, 2D
            switched_name = 1
        except:
            print("Error: check names")
    
        

    h_grid = np.loadtxt(filename,max_rows=1) # edges
    v_grid = np.loadtxt(filename,skiprows=1, max_rows=1) #edges
    
    hist_2d_counts = hist_2d_counts.reshape((len(h_grid) - 1,len(v_grid) - 1)) # centers

    h_mesh,v_mesh = np.meshgrid(h_grid,v_grid,indexing='ij')
        
    cmap = plt.get_cmap(cmap)
    new_cmap = truncate_colormap(cmap,  maxval = cm_maxval, minval=cm_minval,)

    if (transpose_default == 0):
        plt.pcolormesh(h_mesh,v_mesh, hist_2d_counts,shading='gauraud',cmap=new_cmap)  
        plt.xlabel(axis2,fontsize=16)                                
        plt.ylabel(axis1,fontsize=16)   
    else:
        plt.pcolormesh(v_mesh,h_mesh, hist_2d_counts,shading='gauraud',cmap=new_cmap) 
        plt.xlabel(axis1,fontsize=16)                                
        plt.ylabel(axis2,fontsize=16)   

    
    if switched_name == 0:
        if transpose_default == 0:
            plt.xlabel(axis2,fontsize=16)                                
            plt.ylabel(axis1,fontsize=16)   
        else:
            plt.xlabel(axis1,fontsize=16)                                
            plt.ylabel(axis2,fontsize=16)  
    else:
        if transpose_default == 0:
            plt.xlabel(axis1,fontsize=16)                                
            plt.ylabel(axis2,fontsize=16)   
        else:
            plt.xlabel(axis2,fontsize=16)                                
            plt.ylabel(axis1,fontsize=16)  


    cbar = plt.colorbar()                                                                          
    max_count = np.amax(hist_2d_counts)                                                       
    min_count = np.amin(hist_2d_counts)  

    if (zlim == 0):
        zlim = max_count
    


    cbar.ax.tick_params(labelsize=18)                                                          
    cbar.set_label(r'log density (code units)',fontsize=16)  
    plt.clim(0,zlim)
                                                                                           

    print("Total counts is %d" % hist_2d_counts.sum())


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap 
