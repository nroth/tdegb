import numpy as np
import matplotlib
import matplotlib.pyplot as plt


class Binner:


    def __init__(self,data_path = './',hist_type = 'gals',save_path = './'):

        self.data_path = data_path
        self.hist_type = hist_type
        self.save_path = save_path

    def set_data_path(self,data_path):

        self.data_path = data_path
        
    def set_hist_type(self,hist_type):

        self.hist_type = hist_type

    def set_save_path(self,save_path):

        self.save_path = save_path


    def read_1d_histogram(self,axis):

        self.axis = axis
        self.working_filename = self.data_path + self.hist_type + '_' + self.axis + '_1d.hist'
        self.savefile_name = self.working_filename.replace('.hist','.png')

        self.hist_1d_counts = np.loadtxt(self.working_filename)


    def rescale_1d_counts(self,scalar):

        self.hist_1d_counts[:,2] *= scalar

    def plot_1d_histogram(self,axis,rescale_by_binwidth = False, value_only = False):


        left_edges = self.hist_1d_counts[::,0]
        right_edges = self.hist_1d_counts[::,1]
        bin_centers = 0.5 * (left_edges + right_edges) # assumes uniform spacing
        widths = right_edges - left_edges
        counts = self.hist_1d_counts[::,2]

        self.rescale_1d_by_binwidth = rescale_by_binwidth

        if(self.rescale_1d_by_binwidth):
            counts /= widths
            

        plt.xlabel(axis,fontsize=14)
        if (self.hist_type == 'flares'):
            plt.title("Host galaxies of detected flares",fontsize=14)
            if (self.rescale_1d_by_binwidth):
                plt.ylabel(r'Expected # of ZTF detections yr$^{-1}$ (bin interval)$^{-1}$',fontsize=14)
            else:
                plt.ylabel(r'Expected # of ZTF detections / yr / bin',fontsize=14)
        if (self.hist_type == 'vol_disrupt' or self.hist_type == 'gals'):
            plt.title("Catalogue galaxies",fontsize=14)
            plt.ylabel(r'Number of galaxies / bin',fontsize=14)

        print("Total counts is %f" % self.hist_1d_counts[:,2].sum())

        if (value_only):
            return self.hist_1d_counts, plt.errorbar(bin_centers,counts, xerr= widths/2., yerr = 0, elinewidth=1,color='k',fmt='o',markersize=3,label = 'model')
        else:
            return self.hist_1d_counts, plt.bar(bin_centers,counts, width = widths, align='center',edgecolor='k',color='white')
    

    def save_current_histogram(self,fn = ''):

        if (fn == 'pdf'):
            self.savefile_name = self.savefile_name.replace('.png','.pdf')
        else:
            if (fn != ''):
                self.savefile_name = self.data_path + fn
        print("saving to ",self.savefile_name)
        plt.savefig(self.savefile_name)

    def read_2d_histogram(self,axis1='log_mbh_sigma',axis2 = 'z',transpose_default = 0):
        try:
            self.working_filename = self.data_path + self.hist_type + "_" + axis1 + "_" + axis2 + "_2d.hist"
            self.savefile_name = self.working_filename.replace('.hist','.png')
            self.hist_2d_counts = np.loadtxt(self.working_filename,skiprows=2,usecols=2)
            self.switched_2d_name = 0
            print("Read in data from", self.working_filename)
        except:
            try:
                self.working_filename = self.data_path + self.hist_type + "_" + axis2 + "_" + axis1 + "_2d.hist"
                self.savefile_name = self.working_filename.replace('.hist','.png')
                self.hist_2d_counts = np.loadtxt(self.working_filename,skiprows=2,usecols=2)
                self.switched_2d_name = 1
                print("Read in data from ", self.working_filename)
            except:
                
                print("Error: check names", self.working_filename)
    
        self.transpose_default = transpose_default

        if self.switched_2d_name == 0:
            if self.transpose_default == 0:
                plt.xlabel(axis2,fontsize=16)                                
                plt.ylabel(axis1,fontsize=16)   
            else:
                plt.xlabel(axis1,fontsize=16)                                
                plt.ylabel(axis2,fontsize=16)  
        else:
            if self.transpose_default == 0:
                plt.xlabel(axis1,fontsize=16)                                
                plt.ylabel(axis2,fontsize=16)   
            else:
                plt.xlabel(axis2,fontsize=16)                                
                plt.ylabel(axis1,fontsize=16)  


    def rescale_2d_counts(self,scalar):

        self.hist_2d_counts *= scalar

    # trying to use gouraud shading will result in error: I think it needs a different grid shape
    def plot_2d_histogram(self, zlim = 0, cmap_name = 'jet',cm_maxval = 1.,cm_minval = 0.,bin_shading = 'flat'):


        h_grid = np.loadtxt(self.working_filename,max_rows=1) # edges
        v_grid = np.loadtxt(self.working_filename,skiprows=1, max_rows=1) #edges
    
        hist_2d_counts = self.hist_2d_counts.reshape((len(h_grid) - 1,len(v_grid) - 1)) # centers

        h_mesh,v_mesh = np.meshgrid(h_grid,v_grid,indexing='ij')

        cmap = plt.get_cmap(cmap_name)
        new_cmap = truncate_colormap(cmap,  maxval = cm_maxval, minval=cm_minval,)

        if (self.transpose_default == 0):
            hist2d = plt.pcolormesh(h_mesh,v_mesh, hist_2d_counts,shading=bin_shading,cmap=new_cmap)  
        else:
            hist2d = plt.pcolormesh(v_mesh,h_mesh, hist_2d_counts,shading=bin_shading,cmap=new_cmap)


        cbar = plt.colorbar()                                                                          
        max_count = np.amax(hist_2d_counts)                                                       
        min_count = np.amin(hist_2d_counts)  

        if (zlim == 0):
            zlim = max_count
        
        cbar.ax.tick_params(labelsize=18)                                                          
        cbar.set_label(r'rate',fontsize=16)  
        plt.clim(0,zlim)
        plt.title(self.hist_type)
                                                                                           

        print("Total counts is %f" % hist_2d_counts.sum())

        return hist_2d_counts, hist2d, cbar

    def plot_2d_contours(self, level_list = [0.1,0.3,0.5,0.7,0.9]):

        h_grid = np.loadtxt(self.working_filename,max_rows=1) # edges
        v_grid = np.loadtxt(self.working_filename,skiprows=1, max_rows=1) #edges


        # contours want bin-centers, not bin-edges like pcolormesh
        h_grid = np.array([0.5 * (h_grid[i] + h_grid[i+1]) for i in range(len(h_grid) - 1)])
        v_grid = np.array([0.5 * (v_grid[i] + v_grid[i+1]) for i in range(len(v_grid) - 1)])
    
        hist_2d_counts = self.hist_2d_counts.reshape((len(h_grid),len(v_grid))) # centers

        h_mesh,v_mesh = np.meshgrid(h_grid,v_grid,indexing='ij')

        max_count = np.amax(hist_2d_counts)                                                      
        min_count = np.amin(hist_2d_counts)
        #contour_levels = np.array([0.1 * max_count,0.3 * max_count,0.5 * max_count, 0.7 * max_count,  0.9 * max_count])
        contour_levels = max_count * np.array(level_list)

        if (self.transpose_default == 0):
            cont = plt.gca().contour(h_mesh,v_mesh,hist_2d_counts,contour_levels,colors='k',linewidths=2)
        else:
            cont = plt.gca().contour(v_mesh,h_mesh,hist_2d_counts,contour_levels,colors='k',linewidths=2)

        print("Total counts is %f" % hist_2d_counts.sum())

        return hist_2d_counts, cont
        
# static
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = matplotlib.colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),cmap(np.linspace(minval, maxval, n)))
    return new_cmap 


# old way
def plot_2d_histogram(prefix ='./', hist_type = 'gals',axis1 = 'log_mbh_sigma', axis2 = 'z', transpose_default = 0, cmap = "YlGnBu", cm_maxval = 1.0, cm_minval = 0., zlim = 0):

    switched_name = 0
    try:
        filename = prefix + hist_type + "_" + axis1 + "_" + axis2 + "_2d.hist"
        hist_2d_counts = np.loadtxt(filename,skiprows=2,usecols=2) 
    except:
        try:
            filename = prefix + hist_type + "_" + axis2 + "_" + axis1 + "_2d.hist"
            hist_2d_counts = np.loadtxt(filename,skiprows=2,usecols=2)
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
        plt.pcolormesh(h_mesh,v_mesh, hist_2d_counts,shading='flat',cmap=new_cmap)  
    else:
        plt.pcolormesh(v_mesh,h_mesh, hist_2d_counts,shading='flat',cmap=new_cmap) 

    # need to figure out where to do this when plotting multiple things
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
    cbar.set_label(r'rate',fontsize=16)  
    plt.clim(0,zlim)
                                                                                           

    print("Total counts is %d" % hist_2d_counts.sum())


# old way
def plot_2d_contours(prefix ='./', hist_type = 'gals',axis1 = 'mstar_mendel', axis2 = 'UminusR', transpose_default = 0, cmap = "YlGnBu", cm_maxval = 1.0, cm_minval = 0., zlim = 0):

    switched_name = 0
    try:
        filename = prefix + hist_type + "_" + axis1 + "_" + axis2 + "_2d.hist"
        hist_2d_counts = np.loadtxt(filename,skiprows=2,usecols=2) 
    except:
        try:
            filename = prefix + hist_type + "_" + axis2 + "_" + axis1 + "_2d.hist"
            hist_2d_counts = np.loadtxt(filename,skiprows=2,usecols=2) 
            switched_name = 1
        except:
            print("Error: check names")


    h_grid = np.loadtxt(filename,max_rows=1) # edges
    v_grid = np.loadtxt(filename,skiprows=1, max_rows=1) #edges

    # contours want bin-centers, not bin-edges like pcolormesh
    h_grid = np.array([0.5 * (h_grid[i] + h_grid[i+1]) for i in range(len(h_grid) - 1)])
    v_grid = np.array([0.5 * (v_grid[i] + v_grid[i+1]) for i in range(len(v_grid) - 1)])
    
    hist_2d_counts = hist_2d_counts.reshape((len(h_grid),len(v_grid))) # centers

    h_mesh,v_mesh = np.meshgrid(h_grid,v_grid,indexing='ij')

    max_count = np.amax(hist_2d_counts)                                                      
    min_count = np.amin(hist_2d_counts)
    contour_levels = np.array([0.1 * max_count,0.3 * max_count,0.5 * max_count, 0.7 * max_count,  0.9 * max_count])

    if (transpose_default == 0):
        plt.contour(h_mesh,v_mesh,hist_2d_counts,contour_levels,colors='k',linewidths=2)
    else:
        plt.contour(v_mesh,h_mesh,hist_2d_counts,contour_levels,colors='k',linewidths=2)

    print("Total counts is %f" % hist_2d_counts.sum())



