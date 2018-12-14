def Skyrmion_vel_wrapper(Skx_control):
    """Wrapper function to handle the data layer of the calculation of time dependent
    properties of the skyrmions.
    Args
    ----------
        - Skx_control: (dict) dictionary containing the control parameters of the data handling and Figure plotting.
    Author
    ----------
    Jonathan Chico
    """
    import os
    import SK_aux 
    print("Studying skyrmion dynamics")
    #---------------------------------------------------------------------------
    # Create the output directory
    #---------------------------------------------------------------------------
    output_dir=Skx_control['Data_control']['out_path']+'ANALYSIS/'
    try:
        os.makedirs(output_dir)
    except:
        print("Output directory already exists")
    #---------------------------------------------------------------------------
    # Create an instance with the data of the skyrmion
    #---------------------------------------------------------------------------
    Skx_data=SK_aux.Skx_data(path=Skx_control['Data_control']['path'],         \
                        zipped=Skx_control['Data_control']['Misc']['zipped'],  \
                        file_pref=['coord','moment','localenergy'],            \
                        file_suff=['.out','.out','.out'])
    #---------------------------------------------------------------------------
    # Find the orientation of the background and core of the skyrmion
    #---------------------------------------------------------------------------
    Skx_data.find_background()
    #---------------------------------------------------------------------------
    # Make sure the positions are in nm
    #---------------------------------------------------------------------------
    Skx_data.pos=Skx_data.pos*Skx_control['Data_control']['Misc']['alat']
    #---------------------------------------------------------------------------
    # Make sure the times are in s
    #---------------------------------------------------------------------------
    Skx_data.times=\
    Skx_data.times*Skx_control['Data_control']['Skyrmion_velocity']['time_step']
    #---------------------------------------------------------------------------
    # Calculate the velocity, position and radius as a function of time
    #---------------------------------------------------------------------------
    Skyrmion_vel(output_dir,Skx_control,Skx_data)
    print("Velocity analysis done!")
    return

def get_bare_data(itr,Skx_data,Skx_plots,Skx_control,x_mesh,y_mesh,del_val,    \
    index_del):
    
    import numpy as np
    import SK_geometry
    from scipy.interpolate import griddata
    from SK_topology import topological_charge
    from SK_aux import progress_bar
    import SK_geometry

    start_pos=itr*Skx_data.Natoms
    final_pos=(itr+1)*Skx_data.Natoms
    Skx_geometry=SK_geometry.Skx_geometry()
    #---------------------------------------------------------------------------
    # Find the indexes that belong to the skyrmion core
    #---------------------------------------------------------------------------
    trim_ind=np.where(Skx_data.mom[start_pos:final_pos,2]<\
        Skx_control['Data_control']['Skx_threshold'])[0]
    temp_coord=np.zeros([len(trim_ind),3],dtype=np.float64)
    # Check if there is a skyrmion at all in this time step
    if len(trim_ind)>0:
        temp_coord=Skx_data.pos[trim_ind[:]%Skx_data.Natoms,:]
        #-----------------------------------------------------------------------
        # Store the data in temporary lists
        # This should allow one to avoid the possibility of having a point in
        # which there is no skyrmion and one is trying to calculate its velocity
        #-----------------------------------------------------------------------
        #-----------------------------------------------------------------------
        # Store the different times
        #-----------------------------------------------------------------------
        curr_time=Skx_data.times[start_pos]
        #-----------------------------------------------------------------------
        # Calculate the position of the skyrmion by averaging over the position
        # of the spins that belong to the core
        #-----------------------------------------------------------------------
        curr_pos_bare = np.average(temp_coord,axis=0)
        #-----------------------------------------------------------------------
        # Calculate the skyrmion position from interpolation of the
        # magnetization
        #-----------------------------------------------------------------------
        # Create the mesh over which one will be interpolating
        curr_mom_grid = griddata(Skx_data.pos[:,0:2],                          \
            Skx_data.mom[start_pos:final_pos,2],(x_mesh, y_mesh),method='cubic')
        # Now need to filter the magnetization to obtain only the skyrmion 
        # center
        curr_pos_grid = Skx_geometry.find_smoothed(x_mesh,y_mesh,curr_mom_grid,\
            Skx_data.core)
        #-----------------------------------------------------------------------
        # Calculate the perimeter of the skyrmion from the bare data
        #-----------------------------------------------------------------------
        perimeter=Skx_geometry.find_skyrmion_perimeter(temp_coord,tol=0.05)
        #-----------------------------------------------------------------------
        # Calculate the perimeter of the skyrmions from the interpolated data
        #-----------------------------------------------------------------------
        if Skx_data.core<0:
            if Skx_data.core<0 and \
                Skx_control['Data_control']['Skx_threshold']>0:
                t_sign=1.0
            if Skx_data.core<0 and \
            Skx_control['Data_control']['Skx_threshold']<0:
                t_sign=-1.0
            indx=np.where(np.less(curr_mom_grid,\
                t_sign*Skx_control['Data_control']['Skx_threshold']))
        else:
            if Skx_data.core>0 and \
                Skx_control['Data_control']['Skx_threshold']>0:
                t_sign=-1.0
            if Skx_data.core>0 and \
                Skx_control['Data_control']['Skx_threshold']<0:
                t_sign=1.0
            indx=np.where(np.greater(curr_mom_grid,\
                t_sign*Skx_control['Data_control']['Skx_threshold']))
        smooth_coord=np.zeros([len(indx[0]),3],dtype=np.float64)
        smooth_coord[:,0]=x_mesh[indx[0][:],indx[1][:]]
        smooth_coord[:,1]=y_mesh[indx[0][:],indx[1][:]]
        #-----------------------------------------------------------------------
        # Calculate the perimeter of the skyrmion from the interpolated data
        #-----------------------------------------------------------------------
        perimeter_smooth=Skx_geometry.find_skyrmion_perimeter(smooth_coord,    \
            tol=0.0005)
        #-----------------------------------------------------------------------
        # Determine the radius of the skyrmion from the bare data
        #-----------------------------------------------------------------------
        (ave_dev,curr_ave_rad)=\
        Skx_geometry.find_skyrmion_deviation(perimeter,curr_pos_bare[0:2])
        #------------------------------------------------------------------------
        # Determine the radius of the skyrmion
        #-----------------------------------------------------------------------
        (ave_grid_dev,curr_grid_rad)=\
        Skx_geometry.find_skyrmion_deviation(perimeter_smooth,\
            curr_pos_grid[0:2])
        #-----------------------------------------------------------------------
        # Calculate topological information
        #-----------------------------------------------------------------------
        if Skx_control['Data_control']['topology']['execute']:
            curr_charge=topological_charge(del_val,Skx_data.pos,index_del,     \
                Skx_data.mom[start_pos:final_pos,:])
        if Skx_control['Data_control']['Skyrmion_velocity']['plot_energy']:
            #-------------------------------------------------------------------
            # Masking the magnetization array so that only the skyrmion region 
            # in kept
            #-------------------------------------------------------------------
            if Skx_data.core<0:
                if Skx_data.core<0 and Skx_control['Data_control']['Skx_threshold']>0:
                    t_sign=1.0
                if Skx_data.core<0 and Skx_control['Data_control']['Skx_threshold']<0:
                    t_sign=-1.0
                masked_m   = np.ma.masked_where(np.greater(curr_mom_grid,\
                    t_sign*Skx_control['Data_control']['Skx_threshold']),\
                    curr_mom_grid)
            else:
                if Skx_data.core>0 and Skx_control['Data_control']['Skx_threshold']>0:
                    t_sign=-1.0
                if Skx_data.core>0 and Skx_control['Data_control']['Skx_threshold']<0:
                    t_sign=1.0
                masked_m   = np.ma.masked_where(np.less(curr_mom_grid,\
                    t_sign*Skx_control['Data_control']['Skx_threshold']),\
                    curr_mom_grid)
            #-------------------------------------------------------------------
            # Calculate the grided total energy
            #-------------------------------------------------------------------
            tot_curr_ene_grid=griddata(Skx_data.pos[:,0:2],                    \
                Skx_data.ene[start_pos:final_pos,0],(x_mesh, y_mesh),          \
                method='cubic')
            #-------------------------------------------------------------------
            # Calculate the grided decomposed energy
            #-------------------------------------------------------------------
            if Skx_control['Data_control']['Skyrmion_velocity']['comp_energy']:
                xc_curr_ene_grid=griddata(Skx_data.pos[:,0:2],                 \
                Skx_data.ene[start_pos:final_pos,1],(x_mesh, y_mesh),          \
                method='cubic')
                dm_curr_ene_grid=griddata(Skx_data.pos[:,0:2],                 \
                Skx_data.ene[start_pos:final_pos,2],(x_mesh, y_mesh),          \
                method='cubic')
            #-------------------------------------------------------------------
            # Calculating the force
            #-------------------------------------------------------------------
            tot_ene_grad   = np.gradient(tot_curr_ene_grid,edge_order=2)
            if Skx_control['Data_control']['Skyrmion_velocity']['comp_energy']:
                xc_ene_grad     = np.gradient(xc_curr_ene_grid,edge_order=2)
                dm_ene_grad     = np.gradient(dm_curr_ene_grid,edge_order=2)
            #-------------------------------------------------------------------
            # Using the mask to only display the forces in the skyrmion
            #-------------------------------------------------------------------
            ene_grad_x = -1.0*np.ma.masked_array(tot_ene_grad[0], masked_m.mask)
            ene_grad_y = -1.0*np.ma.masked_array(tot_ene_grad[1], masked_m.mask)
            if Skx_control['Data_control']['Skyrmion_velocity']['comp_energy']:
                #---------------------------------------------------------------
                # Exchange contribution
                #---------------------------------------------------------------
                xc_ene_grad_x = -1.0*np.ma.masked_array(xc_ene_grad[0],        \
                    masked_m.mask)
                xc_ene_grad_y = -1.0*np.ma.masked_array(xc_ene_grad[1],        \
                    masked_m.mask)
                #---------------------------------------------------------------
                # DM contribution
                #---------------------------------------------------------------
                dm_ene_grad_x = -1.0*np.ma.masked_array(dm_ene_grad[0],        \
                    masked_m.mask)
                dm_ene_grad_y = -1.0*np.ma.masked_array(dm_ene_grad[1],        \
                    masked_m.mask)
                #---------------------------------------------------------------
                # Save the forces
                #---------------------------------------------------------------
                curr_force=[ene_grad_x,ene_grad_y]
                curr_force_xc=[xc_ene_grad_x,xc_ene_grad_y]
                curr_force_dm=[dm_ene_grad_x,dm_ene_grad_y]
                #---------------------------------------------------------------
                # Cleanup
                #---------------------------------------------------------------
                del xc_curr_ene_grid,dm_curr_ene_grid,tot_curr_ene_grid,masked_m
            else:
                curr_force=[ene_grad_x,ene_grad_y]
                curr_force_xc=[]
                curr_force_dm=[]
                #---------------------------------------------------------------
                # Cleanup
                #---------------------------------------------------------------
                del tot_curr_ene_grid,masked_m
        else:
            curr_force=[]
            curr_force_xc=[]
            curr_force_dm=[]
        #-----------------------------------------------------------------------
        # Plot the contour and position of the center
        #-----------------------------------------------------------------------
        if Skx_control['Data_control']['Skyrmion_velocity']['plot']:
            Skx_plots.plot_time_profile(x_mesh,y_mesh,curr_mom_grid,itr,\
                Skx_data.num_times,curr_pos_grid)
    del perimeter,perimeter_smooth
    del temp_coord,smooth_coord,indx

    return curr_time,curr_pos_bare,curr_pos_grid,curr_ave_rad,curr_grid_rad,   \
    curr_charge,curr_mom_grid,curr_force,curr_force_xc,curr_force_dm

#-------------------------------------------------------------------------------
# Function to calculate the skyrmion Velocity
#-------------------------------------------------------------------------------
def Skyrmion_vel(out_folder,Skx_control,Skx_data):
    """Calculation of the time-dependent position, Velocity and geometric properties of the skyrmion
    The position is determined by calculating the average of the positions belonging to the skyrmion
    body. 
    The position is the smoothed making use of the sagvol filter included in 'scipy'
    The Velocity is calculated by taking both the raw and smoothed positions and calculating the 
    gradient using the time snapshots. The obtained Velocity from the smoothed data is then smoothed
    once more.
    The obtained results are the plotted and save into output files for further processing.
    Args
    ----------
        - diff_times: (int) number of different measured times.
        - nrAtoms: (int) number of atoms in the simulation box.
        - mom: (float [N,3] array) magnetization for each atom in the simulation box.
        - times: (float [N] array) array containing the measured times in seconds.
        - coord: (float [N,3] array) array containing the coordinates of each atom in the simulation box.
        - out_folder: (str) path where to store the output data from the analysis.
        - loc: (str) preffix indicating the direction of the applied external current.
        - imp: (str) name of the studied impurity.
        - bext: (str) magnitude of the applied external magnetic field.
        - stt_val: (str) magnitude of the applied external current.
        - m_z_tol: (float) magnitude defining the boundary of the magnetic moment which determines the spins belonging to the skyrmion.
        - Fig_control: (dict) dictionary containing control parameters for the plotting of figures.
        - npoints: (int [2] array) number of points used for the interpolation
        - topology: (boolean) flag to decide if the topological information should be calculated.
        - ind_sub: (int) index for a substrate atom
    Results
    ----------
        - pos_grid: (float [times,3] array) position of the center of the skyrmion in the interpolated grid
    Author
    ----------
    Jonathan Chico
    """
    import numpy as np
    from scipy.signal import savgol_filter
    from SK_aux import progress_bar
    import SK_plots
    import SK_geometry
    from SK_IO_control import SK_overview
    #---------------------------------------------------------------------------
    # Format for the different files
    #---------------------------------------------------------------------------
    pos_fmt = " {: 6d}   {: 4.8E}   {: 4.8E}   {: 4.8E}   {: 4.8E}   {: 4.8E}   {: 4.8E}\n"
    vel_fmt = " {: 6d}   {: 4.8E}   {: 4.8E}   {: 4.8E}   {: 4.8E}   {: 4.8E}\n"
    topo_fmt= " {: 6d}   {: 4.8E}   {: 4.8E}\n"
    head_pos_fmt   = " {:>6s}   {:>15s}   {:>15s}   {:>15s}   {:>15s}   {:>15s}   {:>15s}\n"
    head_vel_fmt = " {:>6s}   {:>15s}   {:>15s}   {:>15s}   {:>15s}   {:>15s}\n"
    head_topo_fmt= " {:>6s}   {:>15s}   {:>15s}\n"
    #---------------------------------------------------------------------------
    # Parameters for the plotting
    #---------------------------------------------------------------------------
    tol=0.1
    #---------------------------------------------------------------------------
    # Trimming the data to separate them in different times
    #---------------------------------------------------------------------------
    rad_bare=[]
    rad_grid=[]
    pos_bare=[]
    pos_grid=[]
    mag_data=[]
    force_data=[]
    time_array=[]
    topo_charge=[]
    pos_bare_mod=[]
    pos_grid_mod=[]
    force_data_xc=[]
    force_data_dm=[]
    Skx_geometry=SK_geometry.Skx_geometry()
    #---------------------------------------------------------------------------
    # Create an instance of the class for the plots
    #---------------------------------------------------------------------------
    Skx_plots=SK_plots.Skx_plots(Skx_control['Fig_control'])
    #---------------------------------------------------------------------------
    # Calculate the boundaries of the sample
    #---------------------------------------------------------------------------
    bounds=Skx_data.find_bound()
    #---------------------------------------------------------------------------
    # Check if the time dependent contour of the skyrmion will be plotted
    #---------------------------------------------------------------------------
    if Skx_control['Data_control']['Skyrmion_velocity']['plot']:
        uniq=np.unique(Skx_data.times)
        if (len(uniq)>1):
            step_size=(uniq[1]-uniq[0])
        else:
            step_size=1
        Skx_plots.create_time_plot(step_size,Skx_data.num_times,bounds[:4])
    #---------------------------------------------------------------------------
    # Interpolate the coordinates
    #---------------------------------------------------------------------------
    x_grid=np.linspace(bounds[0],bounds[1],                                    \
        Skx_control['Data_control']['Misc']['gridpoint'][0])
    y_grid=np.linspace(bounds[2],bounds[3],                                    \
        Skx_control['Data_control']['Misc']['gridpoint'][1])
    # Create the mesh over which one will be interpolating
    x_mesh,y_mesh = np.meshgrid(x_grid,y_grid)
    del x_grid, y_grid
    #---------------------------------------------------------------------------
    # Calculate the Delaunay tessellation of the coordinates
    #---------------------------------------------------------------------------
    if Skx_control['Data_control']['topology']['execute']:
        (del_val,index_del)=Skx_geometry.triangulate_data(Skx_data.pos,        \
            Skx_control['Data_control']['Misc']['sub_index'])
    else:
        del_val=[]
        index_del=[]
    #---------------------------------------------------------------------------
    # Loop along the different time steps
    #---------------------------------------------------------------------------
    for itr in range(0,Skx_data.num_times):
        #-----------------------------------------------------------------------
        # Function to gather the bare data 
        #-----------------------------------------------------------------------
        (curr_time,curr_pos_bare,curr_pos_grid,curr_ave_rad,curr_grid_rad,     \
        curr_charge,curr_mag_conf,curr_force,curr_force_xc,curr_force_dm)=     \
        get_bare_data(itr,Skx_data,Skx_plots,Skx_control,x_mesh,y_mesh,del_val,\
            index_del)
        #-----------------------------------------------------------------------
        # If the skyrmion is at the edge of the sample stop the calculation of 
        # the velocity to avoid repetitions
        #-----------------------------------------------------------------------
        if abs(curr_pos_bare[0]-bounds[1])<tol or \
            abs(curr_pos_bare[1]-bounds[3])<tol:
                break
        if abs(curr_pos_bare[0]-bounds[0])<tol or \
            abs(curr_pos_bare[1]-bounds[2])<tol:
                break
        time_array.append(curr_time)
        #-----------------------------------------------------------------------
        # Store the bare positions
        #-----------------------------------------------------------------------
        pos_bare.append(curr_pos_bare)
        pos_bare_mod.append(np.sqrt((pos_bare[-1][0]-pos_bare[0][0])**2+\
                                    (pos_bare[-1][1]-pos_bare[0][1])**2+\
                                    (pos_bare[-1][2]-pos_bare[0][2])**2))
        #-----------------------------------------------------------------------
        # Store the interpolated positions
        #-----------------------------------------------------------------------
        pos_grid.append(curr_pos_grid)
        pos_grid_mod.append(np.sqrt((pos_grid[-1][0]-pos_grid[0][0])**2+\
                                    (pos_grid[-1][1]-pos_grid[0][1])**2))
        #-----------------------------------------------------------------------
        # Store the radius
        #-----------------------------------------------------------------------
        rad_bare.append(curr_ave_rad)
        rad_grid.append(curr_grid_rad)
        #-----------------------------------------------------------------------
        # Store the topological information
        #-----------------------------------------------------------------------
        if Skx_control['Data_control']['topology']['execute']:
            topo_charge.append(curr_charge)
        #-----------------------------------------------------------------------
        # Store the information needed for the forces
        #-----------------------------------------------------------------------
        if Skx_control['Data_control']['Skyrmion_velocity']['plot_energy']:
            mag_data.append(curr_mag_conf)
            del curr_mag_conf
            force_data.append(curr_force)
            del curr_force
            #-------------------------------------------------------------------
            # If the decomposed forces should be plotted store their values
            #-------------------------------------------------------------------
            if Skx_control['Data_control']['Skyrmion_velocity']['comp_energy']:
                force_data_xc.append(curr_force_xc)
                del curr_force_xc
                force_data_dm.append(curr_force_dm)
                del curr_force_dm
        progress_bar(itr,Skx_data.num_times,'Time data ')
    #---------------------------------------------------------------------------
    # Find how many different times are found
    #---------------------------------------------------------------------------
    diff_times = len(time_array)
    #---------------------------------------------------------------------------
    # Transform to numpy arrays
    #---------------------------------------------------------------------------
    rad_bare        = np.asarray(rad_bare,dtype=np.float64)
    rad_grid        = np.asarray(rad_grid,dtype=np.float64)
    pos_bare        = np.asarray(pos_bare,dtype=np.float64)
    pos_grid        = np.asarray(pos_grid,dtype=np.float64)
    time_array      = np.asarray(time_array,dtype=np.float64)
    topo_charge     = np.asarray(topo_charge,dtype=np.float64)
    pos_bare_mod    = np.asarray(pos_bare_mod,dtype=np.float64)
    pos_grid_mod    = np.asarray(pos_grid_mod,dtype=np.float64)
    #---------------------------------------------------------------------------
    # Animate the forces plot
    #---------------------------------------------------------------------------
    if Skx_control['Data_control']['Skyrmion_velocity']['plot_energy']:
        # Create an instance of the visualization of the forces
        Force_anim=SK_plots.Force_Animation(                             \
            size=(Skx_control['Fig_control']['width'],\
                    Skx_control['Fig_control']['height']),                 \
            font_size=Skx_control['Fig_control']['fontsize'],                  \
            extent=(bounds[:4]),mesh=(x_mesh,y_mesh),             \
            comp_energy=Skx_control['Data_control']['Skyrmion_velocity']['comp_energy'],\
            imp_pos=Skx_control['Data_control']['Misc']['imp_pos'])
        #-----------------------------------------------------------------------
        # Definting the axis extensions at every time step for the visualization
        # of the forces. This ensures that the skyrmion is kept centered
        #-----------------------------------------------------------------------
        ax_data = [ pos_grid[:,0]-np.amax(rad_grid)*1.5,\
                    pos_grid[:,0]+np.amax(rad_grid)*1.5,\
                    pos_grid[:,1]-np.amax(rad_grid)*1.5,\
                    pos_grid[:,1]+np.amax(rad_grid)*1.5]
        ani_name=out_folder+\
        '/SK_force'+Skx_control['Data_control']['file_name']+'.mp4'
        Force_anim.save_force_plot(mag_data,diff_times,ax_data,ani_name,       \
            Skx_control['Data_control']['Skyrmion_velocity']['comp_energy'],   \
            force_data,force_data_xc,force_data_dm)
        del mag_data,ax_data,force_data,force_data_xc,force_data_dm
    #---------------------------------------------------------------------------
    # Save the contour plot of the skyrmion to file
    #---------------------------------------------------------------------------
    if Skx_control['Data_control']['Skyrmion_velocity']['plot']:
        fig_name=out_folder+\
        '/SK_prof'+Skx_control['Data_control']['file_name']+'.pdf'
        Skx_plots.save_time_plot(fig_name,\
            imp_pos=Skx_control['Data_control']['Misc']['imp_pos'])
    #---------------------------------------------------------------------------
    # Cleanup
    #---------------------------------------------------------------------------
    del x_mesh, y_mesh
    #---------------------------------------------------------------------------
    # Clean the topology data
    #---------------------------------------------------------------------------
    if Skx_control['Data_control']['topology']['execute']:
        del del_val,index_del
    #---------------------------------------------------------------------------
    # If there is more than one data point in the array continue
    #---------------------------------------------------------------------------
    if diff_times>3:
        #-----------------------------------------------------------------------
        # Opening file for the writting of data
        #-----------------------------------------------------------------------
        SK_pos_file_name=out_folder+\
        '/SK_pos'+Skx_control['Data_control']['file_name']+'.dat'
        SK_pos= open(SK_pos_file_name, 'w')
        #-----------------------------------------------------------------------
        # Write the data to file
        #-----------------------------------------------------------------------
        SK_pos.write(head_pos_fmt.\
        format("# Itr","Time [s]","R_x [nm]","R_y [nm]","R_z [nm]","\Delta R [nm]","rad [nm]"))
        for ii in range(0,diff_times):
            SK_pos.write(pos_fmt.format(ii,time_array[ii],                     \
                        pos_bare[ii,0],pos_bare[ii,1],pos_bare[ii,2],          \
                        pos_bare_mod[ii],rad_bare[ii]))
        SK_pos.close()
        #-----------------------------------------------------------------------
        # Opening file for the writting of data
        #-----------------------------------------------------------------------
        SK_pos_file_name=out_folder+\
        '/SK_pos_grid'+Skx_control['Data_control']['file_name']+'.dat'
        SK_pos= open(SK_pos_file_name, 'w')
        #-----------------------------------------------------------------------
        # Write the data to file
        #-----------------------------------------------------------------------
        SK_pos.write(head_pos_fmt.\
        format("# Itr","Time [s]","R_x [nm]","R_y [nm]","R_z [nm]","\Delta R [nm]","rad [nm]"))
        for ii in range(0,diff_times):
            SK_pos.write(pos_fmt.format(ii, time_array[ii],                    \
                        pos_grid[ii,0],pos_grid[ii,1],pos_grid[ii,2],          \
                        pos_grid_mod[ii],rad_grid[ii]))
        SK_pos.close()
        #-----------------------------------------------------------------------
        # Find the the value of the filter for the Velocity
        #-----------------------------------------------------------------------
        if int(diff_times*0.5)%2==0:
            vel_filter=int(diff_times*0.10)+1
        else:
            vel_filter=int(diff_times*0.10)
        if vel_filter<=5:
            vel_filter=7
        elif vel_filter%2==0:
            vel_filter=vel_filter+1
        #-----------------------------------------------------------------------
        # Calculation of the Velocity in different approaches
        #-----------------------------------------------------------------------
        vel_bare      = np.zeros([diff_times,3],dtype=np.float64)
        vel_grid      = np.zeros([diff_times,3],dtype=np.float64)
        vel_bare_mod  = np.zeros(diff_times,dtype=np.float64)
        vel_grid_mod  = np.zeros(diff_times,dtype=np.float64)
        #-----------------------------------------------------------------------
        # Calculating the Velocity from the bare skyrmion core position
        #-----------------------------------------------------------------------
        vel_bare=np.gradient(pos_bare[:,:],time_array,edge_order=2,axis=0)*1e-9
        vel_bare_mod=np.sqrt(vel_bare[:,0]**2+vel_bare[:,1]**2+vel_bare[:,2]**2)
        #-----------------------------------------------------------------------
        # Calculating the Velocity from the interpolated skyrmion core position
        #-----------------------------------------------------------------------
        vel_grid=np.gradient(pos_grid[:,:],time_array,edge_order=2,axis=0)*1e-9
        #-----------------------------------------------------------------------
        # File name to store the data for the bare position
        #----------------------------------------------------------------------
        file_name=out_folder+\
        '/SK_vel'+Skx_control['Data_control']['file_name']+'.dat'
        SK_vel = open(file_name, 'w')
        #-----------------------------------------------------------------------
        # Write data to file for the bare position
        #-----------------------------------------------------------------------
        SK_vel.write(head_vel_fmt.\
        format("# Itr","Time [s]","V_x [m/s]","V_y [m/s]","V_z [m/s]","V [m/s]"))
        for ii in range(0,diff_times):
            SK_vel.write(vel_fmt.format(ii,time_array[ii],                     \
            vel_bare[ii,0],vel_bare[ii,1],vel_bare[ii,2],vel_bare_mod[ii]))
        SK_vel.close()
        #-----------------------------------------------------------------------
        # File name to store the data for the interpolated position
        #-----------------------------------------------------------------------
        file_name=out_folder+\
        '/SK_vel_grid'+Skx_control['Data_control']['file_name']+'.dat'
        SK_vel = open(file_name, 'w')
        #-----------------------------------------------------------------------
        # Write data to file for the interpolated position
        #-----------------------------------------------------------------------
        SK_vel.write(head_vel_fmt.\
        format("# Itr","Time [s]","V_x [m/s]","V_y [m/s]","V_z [m/s]","V [m/s]"))
        for ii in range(0,diff_times):
            vel_grid_mod[ii]=np.sqrt(vel_grid[ii,:].dot(vel_grid[ii,:]))
            SK_vel.write(vel_fmt.format(ii,time_array[ii],                     \
                        vel_grid[ii,0],vel_grid[ii,1],vel_grid[ii,2],          \
                        vel_grid_mod[ii]))
        SK_vel.close()
        #-----------------------------------------------------------------------
        # Handling the smoothed data
        #-----------------------------------------------------------------------
        pos_filt        = np.zeros([diff_times,3],dtype=np.float64)
        vel_filt        = np.zeros([diff_times,3],dtype=np.float64)
        pos_filt_mod    = np.zeros(diff_times,dtype=np.float64)
        vel_filt_mod    = np.zeros(diff_times,dtype=np.float64)
        #-----------------------------------------------------------------------
        # Filtering the skyrmion position data
        #-----------------------------------------------------------------------
        pos_filt=savgol_filter(pos_bare[:,:],window_length=vel_filter,         \
            polyorder=3,axis=0)
        pos_filt_mod=np.sqrt((pos_filt[:,0]-pos_filt[0,0])**2+\
            (pos_filt[:,1]-pos_filt[0,1])**2+(pos_filt[:,2]-pos_filt[0,2])**2)
        #-----------------------------------------------------------------------
        # Filtering the skyrmion radius
        #-----------------------------------------------------------------------
        rad_filt=savgol_filter(rad_bare,window_length=vel_filter,polyorder=3,  \
            axis=0)
        #-----------------------------------------------------------------------
        # Calculating the Velocity of the filtered position data
        #-----------------------------------------------------------------------
        vel_filt=np.gradient(pos_filt[:,:],time_array,edge_order=2,axis=0)*1e-9
        vel_filt_mod=np.sqrt(vel_filt[:,0]**2+vel_filt[:,1]**2+vel_filt[:,2]**2)
        #-----------------------------------------------------------------------
        # File name to print out the data for the filtered position 
        #-----------------------------------------------------------------------
        file_name=out_folder+\
        '/SK_pos_ave'+Skx_control['Data_control']['file_name']+'.dat'
        #-----------------------------------------------------------------------
        # Write the filtered position data to file
        #-----------------------------------------------------------------------
        SK_pos= open(file_name, 'w')
        SK_pos.write(head_pos_fmt.\
        format("# Itr","Time [s]","R_x [nm]","R_y [nm]","R_z [nm]","R [nm]","rad [nm]"))
        for ii in range(0,diff_times):
            SK_pos.write(pos_fmt.format(ii,time_array[ii],                     \
                        pos_filt[ii,0],pos_filt[ii,1],pos_filt[ii,2],          \
                        pos_filt_mod[ii],rad_filt[ii]))
        SK_pos.close()
        #-----------------------------------------------------------------------
        # Write the filtered Velocity data to file
        #-----------------------------------------------------------------------
        file_name=out_folder+\
        '/SK_vel_ave'+Skx_control['Data_control']['file_name']+'.dat'
        SK_vel= open(file_name, 'w')
        SK_vel.write(head_vel_fmt.\
        format("# Itr","Time [s]","V_x [m/s]","V_y [m/s]","V_z [m/s]","V [m/s]"))
        for ii in range(0,diff_times):
            SK_vel.write(vel_fmt.format(ii,time_array[ii],                     \
                        vel_filt[ii,0],vel_filt[ii,1],vel_filt[ii,2],          \
                        vel_filt_mod[ii]))
        SK_vel.close()
        #-----------------------------------------------------------------------
        # Postprocessing the Velocity
        #-----------------------------------------------------------------------
        vel_smooth     = np.zeros([diff_times,3],dtype=np.float64)
        vel_smooth_mod = np.zeros(diff_times,dtype=np.float64)
        #-----------------------------------------------------------------------
        # Smoothing the Velocity data
        #-----------------------------------------------------------------------
        vel_smooth=savgol_filter(vel_filt[:,:],window_length=vel_filter,       \
            polyorder=3,axis=0)
        vel_smooth_mod=np.sqrt( vel_smooth[:,0]**2+vel_smooth[:,1]**2+         \
                                vel_smooth[:,2]**2)
        #-----------------------------------------------------------------------
        # File name to write the data
        #-----------------------------------------------------------------------
        file_name=out_folder+\
        '/SK_vel_smooth'+Skx_control['Data_control']['file_name']+'.dat'
        #-----------------------------------------------------------------------
        # Printing the smoothed Velocity data to file
        #-----------------------------------------------------------------------
        SK_vel = open(file_name, 'w')
        SK_vel.write(head_vel_fmt.format("# Itr","Time [s]","V_x [m/s]","V_y [m/s]","V_z [m/s]","V [m/s]"))
        for ii in range(0,diff_times):
            SK_vel.write(vel_fmt.format(ii,time_array[ii],                     \
                        vel_smooth[ii,0],vel_smooth[ii,1],vel_smooth[ii,2],    \
                        vel_smooth_mod[ii]))
        SK_vel.close()
        #-----------------------------------------------------------------------
        # Topological information
        #-----------------------------------------------------------------------
        if Skx_control['Data_control']['topology']['execute']:
            file_name=out_folder+\
            '/SK_topo'+Skx_control['Data_control']['file_name']+'.dat'
            SK_topo = open(file_name,'w')
            SK_topo.write(head_topo_fmt.format("# Itr","Time [s]","Q_{SK}"))
            for ii in range(0,diff_times):
                SK_topo.write(topo_fmt.format(ii,time_array[ii],topo_charge[ii]))
            SK_topo.close()
        #-----------------------------------------------------------------------
        # Collect the most important data into a file
        #-----------------------------------------------------------------------
        file_name=out_folder+\
        '/SK_overview'+Skx_control['Data_control']['file_name']+'.yml'
        SK_overview(file_name,pos_bare,pos_grid,pos_filt,vel_bare,            \
        vel_bare_mod,vel_grid,vel_grid_mod,vel_filt,vel_filt_mod,vel_smooth,   \
        vel_smooth_mod,rad_bare,rad_grid)
        #-----------------------------------------------------------------------
        # Plot the Velocity of the skyrmion and its components
        #-----------------------------------------------------------------------
        y_data=[[vel_smooth_mod,vel_filt_mod,vel_grid_mod,vel_bare_mod],       \
                [vel_smooth[:,0],vel_filt[:,0],vel_grid[:,0],vel_bare[:,0]],   \
                [vel_smooth[:,1],vel_filt[:,1],vel_grid[:,1],vel_bare[:,1]]]
        x_data=time_array*1e9
        labels=['Smoothed','Filtered','Interpolated','Bare']
        y_labels=\
                [r'$V_{Skx}$ [m/s]',r'$V_{Skx}^x$ [m/s]',r'$V_{Skx}^y$ [m/s]']
        #-----------------------------------------------------------------------
        # Plotting the actual figure
        #-----------------------------------------------------------------------
        fig_name=out_folder+\
        '/SK_vel'+Skx_control['Data_control']['file_name']+'.pdf'
        line_array=['-','--',':']
        Skx_plots.gridded_plot(x_data,y_data,line_array,labels,fig_name,       \
            y_labels)
        #-----------------------------------------------------------------------
        # Plot the positions of the skyrmion and its components
        #-----------------------------------------------------------------------
        y_data=[[pos_filt_mod,pos_grid_mod,pos_bare_mod],                      \
                [pos_filt[:,0],pos_grid[:,0],pos_bare[:,0]],                   \
                [pos_filt[:,1],pos_grid[:,1],pos_bare[:,1]]]
        labels=['Filtered','Interpolated','Bare']
        y_labels=\
            [r'$\Delta r_{Skx}$ [nm]',r'$r_{Skx}^x$ [nm]',r'$r_{Skx}^y$ [nm]']
        #-----------------------------------------------------------------------
        # Plotting the actual figure
        #-----------------------------------------------------------------------
        fig_name=out_folder+\
        '/SK_pos'+Skx_control['Data_control']['file_name']+'.pdf'
        line_array=['-','--',':']
        Skx_plots.gridded_plot(x_data,y_data,line_array,labels,fig_name,       \
            y_labels)
        #-----------------------------------------------------------------------
        # Make a plot of the modulus of the radius of the skyrmion
        #-----------------------------------------------------------------------
        y_data=[rad_filt,rad_grid,rad_bare]
        labels=['Filtered','Interpolated','Bare']
        line_array=['--',':','-']
        y_label=r'Radius$_\textrm{SK}$ [nm]'
        x_label=r'Time [ns]'
        fig_name=out_folder+\
        '/SK_rad'+Skx_control['Data_control']['file_name']+'.pdf'
        Skx_plots.line_plot(x_data,y_data,y_label,x_label,labels,line_array,   \
            fig_name)
        #-----------------------------------------------------------------------
        # Make a plot of the topological charge
        #-----------------------------------------------------------------------
        y_data=topo_charge
        labels=''
        line_array='-'
        y_label=r'$Q_{SK}$'
        x_label=r'Time [ns]'
        fig_name=out_folder+\
        '/SK_topo'+Skx_control['Data_control']['file_name']+'.pdf'
        Skx_plots.line_plot(x_data,y_data,y_label,x_label,labels,line_array,   \
            fig_name)
    #---------------------------------------------------------------------------
    # Cleanup
    #---------------------------------------------------------------------------
    del pos_filt,vel_filt,pos_filt_mod,vel_filt_mod
    del vel_bare,vel_grid,vel_bare_mod,vel_grid_mod
    del vel_smooth,vel_smooth_mod 
    del y_data,x_data
    del rad_bare,rad_grid,pos_bare,pos_grid,time_array,topo_charge
    del pos_bare_mod,pos_grid_mod
    return 
