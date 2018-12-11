class Skx_geometry:
    """Class containing the methods to determine the geometrical information of the skyrmions.
    """
    def __init__(self):
        pass

    def triangulate_data(self,coord,ind_z=0):
        """Delaunay triangulation of the coordinates making use of 'qhull' which are oriented in counterclockwise fashion.
        Right now the triangulation is performed in slices over the z-coordinates of the system, i.e, in planes 
        with the chosen plane being selected by the value of the z-coordinate.
        Args
        ----------
            - coord: (float [3,N] array) coordinates of the atoms in the simulation box.
            - ind_z: (int) index of the z-coordinate over which one is slicing the sample. (default z=0)
        Returns
        ----------
            - del_val: (integer [N,3] array) vertices obtained from the Delaunay triangulation. 
            - ind_sub: (integer [N] array) list of the atoms which correspond to the selected plane.
        Author
        ----------
        Jonathan Chico
        """
        import numpy as np
        from scipy.spatial import Delaunay

        tol=1e-5

        ind_sub=np.where(abs(coord[:,2]-coord[ind_z,2])<tol)[0]

        del_vert=Delaunay(coord[ind_sub,0:2])

        del_val=del_vert.simplices

        return del_val,ind_sub


    def find_skyrmion_perimeter(self,coord,tol=0.0001):
        import numpy as np
        # Find the unique points in the y-axis
        unique_rij=np.unique(coord[:,1],axis=0)
        # Define the perimeter array
        perimeter=np.zeros([len(unique_rij)*2,3],dtype=np.float64)
        # Find the perimeter for the system
        for ii in range(len(unique_rij)):

            trim_indx=np.where(coord[:,1]==unique_rij[ii])
            # Find the maximum on the x-coordinate
            bound_x=[min(coord[trim_indx[0],0]),max(coord[trim_indx[0],0])]

            # Find the index for the minimum value
            current_index=np.where(np.sqrt((coord[:,0]-bound_x[0])**2+         \
                                            (coord[:,1]-unique_rij[ii])**2)<tol)
            # Assign the values to the perimeter array
            perimeter[ii,0]=bound_x[0]
            perimeter[ii,1]=unique_rij[ii]
            perimeter[ii,2]=coord[current_index[0],2]
            # Find the index for the maximum value
            current_index=np.where(np.sqrt((coord[:,0]-bound_x[1])**2+         \
                                            (coord[:,1]-unique_rij[ii])**2)<tol)
            # Assign the values to the perimeter array
            perimeter[ii+len(unique_rij),0]=bound_x[1]
            perimeter[ii+len(unique_rij),1]=unique_rij[ii]
            perimeter[ii+len(unique_rij),2]=coord[current_index[0],2]

        del trim_indx, unique_rij,current_index

        return perimeter

    def find_skyrmion_deviation(self,perimeter,Sk_center):
        import numpy as np
        bval=0.0
        aval=0.0
        ave_rad=0.0
        ave_dev=0.0
        Deviation=np.zeros(len(perimeter),dtype=np.float64)
        theta_angle=np.zeros(len(perimeter),dtype=np.float64)
        partial_radius=np.zeros([len(perimeter),4],dtype=np.float64)

        partial_radius[:,0]=perimeter[:,0]-Sk_center[0]
        partial_radius[:,1]=perimeter[:,1]-Sk_center[1]
        partial_radius[:,2]=perimeter[:,2]
        partial_radius[:,3]=np.sqrt(partial_radius[:,0]**2+\
                                    partial_radius[:,1]**2+\
                                    partial_radius[:,2]**2)
        theta_angle=np.arctan2(partial_radius[:,1],partial_radius[:,0])
        ave_rad=np.average(partial_radius[:,3])
        aval=2.0*np.average(partial_radius[:,3]*np.cos(theta_angle[:]))
        bval=2.0*np.average(partial_radius[:,3]*np.sin(theta_angle[:]))
        Deviation=partial_radius[:,3]-ave_rad-aval*np.cos(theta_angle[:])-     \
                    bval*np.sin(theta_angle[:])
        ave_dev=np.average(Deviation)

        del Deviation,theta_angle,partial_radius,aval,bval

        return ave_dev,ave_rad

    #---------------------------------------------------------------------------
    # Find the positions from grided data
    #---------------------------------------------------------------------------
    def find_smoothed(self,x_mesh,y_mesh,mom_grid,core,tol=-0.8,rad=50):
        """Find the position of the minimum from interpolated data, making use of the `minimum_filter` function or the `maximum_filter`. Based on similar function by Markus Hoffman.
        Args
        ----------
            - x_mesh: meshed coordinates along the x-direction
            - y_mesh: meshed coordinates along the y-direction
            - mom_grid: grided magnetization data along the z-direction
            - tol: (float) tolerance of the magnetization to find the minimum
            - rad: (int) number of points considered in the filter
        Results
        ----------
            - pos_grid: (float [3] array) position of the minimum along the z-direction
        Author
        ----------
        Jonathan Chico
        """
        import numpy as np
        from scipy.ndimage import minimum_filter,maximum_filter

        if core<0:
            min_val = (mom_grid == minimum_filter(mom_grid,rad))
            valx,valy=np.nonzero(min_val)
            #-------------------------------------------------------------------
            # Find where the minimum is below the tolerance threshold
            #-------------------------------------------------------------------
            indx=np.where(mom_grid[valx,valy]<tol)[0]
            if tol>0:
                tol=-1.0*core
        else:
            max_val = (mom_grid == maximum_filter(mom_grid,rad))
            valx,valy=np.nonzero(max_val)
            if tol<0:
                tol=-1.0*core
            #-------------------------------------------------------------------
            # Find where the minimum is below the tolerance threshold
            #-------------------------------------------------------------------
            indx=np.where(mom_grid[valx,valy]>tol)[0]
        pos_grid=np.zeros([3],dtype=np.float64)
        #-----------------------------------------------------------------------
        # Store the position of the minimum in an array
        #-----------------------------------------------------------------------
        pos_grid[0]=np.average(x_mesh[valx[indx[:]],valy[indx[:]]])
        pos_grid[1]=np.average(y_mesh[valx[indx[:]],valy[indx[:]]])
        #-----------------------------------------------------------------------
        # Cleanup
        #-----------------------------------------------------------------------
        try:
            del min_val
        except:
            pass
        try:
            del max_val
        except:
            pass
        del valx
        del valy
        return pos_grid
