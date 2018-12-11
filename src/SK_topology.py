def solid_angle(A,B,C):
    """ Calculation of the solid angle between 3 spins that form a triangle.
    Args
    ----------
        - A,B,C : (float [3] array) three 3D vectors
    Returns
    ----------
        - angle: solid angle between the 3 vectors in the vertices of a triangle
    Authors
    ----------
    Jonathan Chico
    """
    import numpy as np

    perp=np.cross(B,C)

    x =  A[:,0]*perp[:,0]+A[:,1]*perp[:,1]+A[:,2]*perp[:,2]
    y = 1.0 +   A[:,0]*B[:,0]+A[:,1]*B[:,1]+A[:,2]*B[:,2]+\
                A[:,0]*C[:,0]+A[:,1]*C[:,1]+A[:,2]*C[:,2]+\
                B[:,0]*C[:,0]+B[:,1]*C[:,1]+B[:,2]*C[:,2]
    angle = 2.0*np.arctan2(x,y)

    return angle

def topological_charge(del_vert,coord,ind_sub,mom):
    """ Calculation of the topological charge for ultra-narrow magnetic textures, via the
    calculation of the solid angle spanned by spins located in an oriented triangle. 
    Args
    ----------
        - del_vert: (int [N,3] array) array containing the ordered indices resulting from the delaunay triangulation
        - coord: (float [N,3] array) array of the coordinates of the spins in the simulation box
        - ind_sub: (integer [N] array) list of the atoms which correspond to the selected plane. 
        - mom: (float [N,3] array) magnetic moments of the spins in the simulation box.
    Returns
    ----------
        - charge: (float) topological charge of the magnetic configuration.
    Author
    ----------
    Jonathan Chico
    """
    import numpy as np

    charge=0.0
    temp_coord=coord[ind_sub,:]
    mag_mom=mom[ind_sub,:]

    ver_coord=np.zeros([3,len(del_vert),3],dtype=np.float64)

    ver_coord[0,:,:]=temp_coord[del_vert[:,0],:]
    ver_coord[1,:,:]=temp_coord[del_vert[:,1],:]
    ver_coord[2,:,:]=temp_coord[del_vert[:,2],:]

    AB=ver_coord[0,:,:]-ver_coord[1,:,:]
    AC=ver_coord[0,:,:]-ver_coord[2,:,:]

    # Now one needs to find the normal of the triangle
    triang_norm=np.cross(AB,AC)
    inv_norm_vec=1.0/np.sqrt(triang_norm[:,0]**2+triang_norm[:,1]**2+triang_norm[:,2]**2)
    # Then one calculates the cross product between them to find the normal
    triang_norm[:,0]=triang_norm[:,0]*inv_norm_vec
    triang_norm[:,1]=triang_norm[:,1]*inv_norm_vec
    triang_norm[:,2]=triang_norm[:,2]*inv_norm_vec
    sign_norm=np.sign(triang_norm[:,2])
    # Calculate the partial solid angle
    part_angle=solid_angle(mag_mom[del_vert[:,0],:],mag_mom[del_vert[:,1],:],mag_mom[del_vert[:,2],:])
    # Sum over the site contributions
    charge=np.sum(sign_norm*part_angle)
    charge = charge/(4.0*np.pi)
    return charge

def smooth_data(coord,mom,npoints,core):
    """Function cropping and interpolating data in the proximity of an skyrmion, such that 
    numerical derivatives can more accurately calculated.
    Args
    ----------
        - coord: (float [N,3] array) positions of the atoms in the lattice.
        - mom: (float [N,3] array) magnetic moments of the sample.
        - npoints: (int [2] array) number of points to be used in the data interpolation
    Returns
    ----------
        -x_mesh: (float [npoints[0],npoints[1]] array) interpolated positions of the atoms in the x-direction.
        -y_mesh: (float [npoints[0],npoints[1]] array) interpolated positions of the atoms in the y-direction.
        - mom_grid: (float [3,npoints[0],npoints[1]] array) interpolated magnetization.
    Author
    ----------
    Jonathan Chico
    """
    import numpy as np
    from scipy.interpolate import griddata
    import sys

    tol=1e-6

    indx=np.where(mom[:,2]< -1.0*core-tol)[0]
    x_grid=np.linspace(min(coord[indx,0]),max(coord[indx,0]),npoints[0])
    y_grid=np.linspace(min(coord[indx,1]),max(coord[indx,1]),npoints[1])
    # Create the mesh over which one will be interpolating
    x_mesh,y_mesh = np.meshgrid(x_grid,y_grid)
    find_bound=np.where(coord[:,0]==max(coord[indx,0]))[0]

    mom_grid=[]
    # Create the mesh over which one will be interpolating
    for ii in range(0,3):
        mom_grid.append(griddata(coord[indx,0:2],mom[indx,ii],(x_mesh, y_mesh),method='cubic',fill_value=max(mom[find_bound,ii])))

    return x_mesh,y_mesh,mom_grid


def dissipation_matrix(mag,coord,ind_z=0):
    """ Calculation of the dissipation matrix used in the Thiele equation, which is defined as
    .math::
        :nowrap:
        \begin{equation}
        D_{i,j}=\int\left(\frac{\partial \mathbf{m}}{\partial r_i}\cdot \frac{\partial \mathbf{m}}{\partial r_j} \right)
        \end{equation}
    Args
    ----------
        - mag: (float [N,3] array) magnetic moments of the spins in the simulation box.
        - coord: (float [N,3] array) array of the coordinates of the spins in the simulation box 
        - ind_z: (int) index of the z-coordinate over which one is slicing the sample. (default ind_z=0)
    Returns
    ----------
        - Disp_matrix: (float [3,3] array) Array of the dissipation matrix
    Authors
    ----------
        - Jonathan Chico
        - Markus Hoffmann
    """
    import numpy as np
    from scipy import integrate
    # Calculate the gradients of the magnetization
    grad=[]
    tol=1e-10
    for ii in range(0,3):
        grad.append(np.gradient(mag[ii],edge_order=2))

    Disp_full=np.zeros([3,3,np.shape(grad)[2],np.shape(grad)[3]],dtype=np.float64)
    for mu in range(np.shape(grad)[1]):
        for nu in range(np.shape(grad)[1]):
            Disp_full[mu,nu,ii,jj]= grad[0][mu][:][:]*grad[0][nu][:][:]+       \
                                    grad[1][mu][:][:]*grad[1][nu][:][:]+       \
                                    grad[2][mu][:][:]*grad[2][nu][:][:]
    Disp_matrix=np.zeros([3,3],dtype=np.float64)
    # Integrate the gradients to obtain the matrix
    Disp_matrix=integrate.simps(integrate.simps(Disp_full,axis=-1),axis=-1)
    return Disp_matrix

