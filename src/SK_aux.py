class Skx_data:
    """Class containing the raw data from the skyrmion simulation.
    """
    def __init__(self,path='./',zipped=False,file_pref=[],file_suff=[]):
        import numpy as np
        import sys
        from SK_IO_control import file_handler
        self.pos=[]
        self.mom=[]
        self.ene=[]
        self.ang=[]
        self.times=[]
        self.Natoms=[]
        self.num_times=[]

        try:
            inp_dim=np.shape(file_pref)
        except:
            print('You have to give more than one file name')
            sys.exit()
        for ii in range(len(file_pref)):
            data=file_handler(path,zipped,file_pref[ii],file_suff[ii])
            if file_pref[ii].lower()=='moment'.lower(): 
                self.mom=data[:,1:]
                self.times=data[:,0]
            if file_pref[ii].lower()=='restart'.lower():
                self.mom=data
            if file_pref[ii].lower()=='coord':
                self.pos=data
                self.Natoms=len(self.pos[:,0])
            if file_pref[ii].lower()=='localenergy':
                self.ene=data
            del data
        try:
            self.num_times=int(len(self.times)/self.Natoms)
            if 'localenergy' in file_pref:
                ene_times=int(len(self.ene[:,0])/self.Natom)
                if ene_times !=self.num_times:
                    print('Energy and moments shapes do not conform')
                    sys.exit()
        except:
            pass
        return

    def find_bound(self):
        import numpy as np

        bounds = np.array([ np.amin(self.pos[:,0]),np.amax(self.pos[:,0]),     \
                            np.amin(self.pos[:,1]),np.amax(self.pos[:,1]),     \
                            np.amin(self.pos[:,2]),np.amax(self.pos[:,2])])
        return bounds

    def find_midpoint(self):
        import numpy as np

        mid = np.array([np.amax(self.pos[:,0])-np.amin(self.pos[:,0]),         \
                        np.amax(self.pos[:,1])-np.amin(self.pos[:,1]),         \
                        np.amax(self.pos[:,2])-np.amin(self.pos[:,2])])*0.5
        return mid

    def find_background(self):
        import numpy as np

        back = np.average(self.mom[:self.Natoms,2])

        if back>0.5:
            self.core=-1
        elif back<0.5:
            self.core=1
        else:
            print('Not able to properly determine background polarization, defaulting to +1')
            self.core=-1
        return self.core

    def angle_proj(self):
        import numpy as np
        self.ang=np.arctan2(self.mom[:,1]/np.sqrt(self.mom[:,1]**2+self.mom[:,0]**2),self.mom[:,0]/np.sqrt(self.mom[:,1]**2+self.mom[:,0]**2))
        return self.ang

def progress_bar(itr,max_itr,msg):
    """Progress bar to track the evolution of a certain process. Based on the implementation described in https://stackoverflow.com/questions/3160699/python-progress-bar.
    Args
    ----------
        - itr: (int) current iteration
        - max_itr: (int) maximum number of iterations of the process.
        - msg: (str) display message of the process.
    """
    import sys
    sys.stdout.write('\r')
    progress=(itr+1)/(max_itr)
    sys.stdout.write(msg+"[%-20s] %d%%" % ('='*int(20*progress),100*progress))
    sys.stdout.flush()
    if itr+1==max_itr:
        print(" Done!")
    return

def vector_to_HSV(mag_z,ang,coord,mesh):
    import numpy as np
    from colorsys import hsv_to_rgb
    from scipy.interpolate import griddata

    S=np.where(mag_z>0.0,mag_z, 0.0)
    S=np.ones_like(mag_z)-S
    V=np.where(mag_z<0.0,mag_z, 0.0)
    V=np.ones_like(mag_z)+V
    H=(ang-np.amin(ang))/(np.amax(ang)-np.amin(ang))
    RGB=np.asarray([hsv_to_rgb(H[j],S[j],V[j]) for j in range(np.shape(S)[0])])
    R_grid=griddata(coord,RGB[:,0],(mesh[0],mesh[1]),method='cubic')
    G_grid=griddata(coord,RGB[:,1],(mesh[0],mesh[1]),method='cubic')
    B_grid=griddata(coord,RGB[:,2],(mesh[0],mesh[1]),method='cubic')
    R_grid=(R_grid-np.amin(R_grid))/(np.amax(R_grid)-np.amin(R_grid))
    G_grid=(G_grid-np.amin(G_grid))/(np.amax(G_grid)-np.amin(G_grid))
    B_grid=(B_grid-np.amin(B_grid))/(np.amax(B_grid)-np.amin(B_grid))
    RGB_grid=np.dstack((R_grid,G_grid,B_grid))
    del H,S,V,RGB,R_grid,G_grid,B_grid
    return RGB_grid

def vector_to_HLS(mag_z,mag_x,mag_y,ang,coord,mesh):
    import numpy as np
    from colorsys import hls_to_rgb
    from scipy.interpolate import griddata

    H=(ang-np.amin(ang))/(np.amax(ang)-np.amin(ang))
    R=np.sqrt(mag_x**2+mag_x**2)
    R=(R-np.amin(R))/(np.amax(R)-np.amin(R))
    Z=(mag_z-np.amin(mag_z))/(np.amax(mag_z)-np.amin(mag_z))
    RGB=np.asarray([hls_to_rgb(H[j],Z[j],R[j]) for j in range(np.shape(R)[0])])
    R_grid=griddata(coord,RGB[:,0],(mesh[0],mesh[1]),method='cubic')
    G_grid=griddata(coord,RGB[:,1],(mesh[0],mesh[1]),method='cubic')
    B_grid=griddata(coord,RGB[:,2],(mesh[0],mesh[1]),method='cubic')
    RGB_grid=np.dstack((R_grid,G_grid,B_grid))
    del H,R,Z,RGB,R_grid,G_grid,B_grid
    return RGB_grid

def color_wheel(curr_ax,xticks,xlabels,fontsize):
    import numpy as np
    import matplotlib as mpl

    norm = mpl.colors.Normalize(0,2*np.pi) 
    npoints = 200
    theta = np.linspace(0,2*np.pi,npoints)
    rad   = np.linspace(0,1,2)
    rg,tg = np.meshgrid(rad,theta)
    curr_ax.pcolormesh(theta, rad, tg.T,norm=norm,cmap='hsv')
    curr_ax.set_xticks(xticks)
    curr_ax.set_yticklabels([])
    curr_ax.set_xticklabels(xlabels)
    curr_ax.tick_params(pad=0.0,labelsize=fontsize)
    curr_ax.axvline(0,lw=0.5,linestyle='--',color='black')
    curr_ax.axvline(np.pi*0.5,lw=0.5,linestyle='--',color='black')
    del norm,theta,rad,rg,tg
    pass

################################################################################
# Create the list of folders to loop over
################################################################################
def create_list(path,suffix='',prefix=''):

    import os
    items = sorted(os.listdir(path),key=natural_keys)
    curr_list=[]
    for item in items:
        if os.path.isdir(path+item):
            if len(suffix)>0:
                if item.endswith(suffix):
                    curr_list.append(item[:-len(suffix)])
            if len(prefix)>0:
                if item.startswith(prefix):
                    curr_list.append(item[len(prefix):])
            if len(suffix)==0 and len(prefix)==0:
                curr_list.append(item)
    return curr_list

def atof(text):
    try:
        retval = float(text)
    except ValueError:
        retval = text
    return retval

def natural_keys(text):
    import re
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    float regex comes from https://stackoverflow.com/a/12643073/190597
    '''
    return [ atof(c) for c in re.split(r'[+-]?([0-9]+(?:[.][0-9]*)?|[.][0-9]+)', text) ]