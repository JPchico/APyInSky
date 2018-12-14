def profile_wrapper(Skx_control):
    """Wrapper function to plot the magnetization profile of the skyrmion for the 
    different impurities.
    Args
    ----------
        Skx_control: dictionary containing the control options for data and figure control

    Author
    ----------
    Jonathan Chico
    """
    import os
    import sys
    import numpy as np
    import SK_plots
    import SK_geometry
    import SK_aux
    from scipy.interpolate import griddata
    from SK_IO_control import SK_static_overview
    from SK_topology import topological_charge

    Skx_data=SK_aux.Skx_data(path=Skx_control['Data_control']['path'],         \
                        zipped=Skx_control['Data_control']['Misc']['zipped'],  \
                        file_pref=['coord','restart'],                         \
                        file_suff=['.out','.out'])
    Skx_data.pos=Skx_data.pos*Skx_control['Data_control']['Misc']['alat']
    Skx_data.find_background()

    tol=1e-2
    print("Ploting skyrmion profiles")
    #---------------------------------------------------------------------------
    # Create the output directory
    #---------------------------------------------------------------------------
    output_dir=Skx_control['Data_control']['out_path']+'ANALYSIS/'
    try:
        os.makedirs(output_dir)
    except:
        print("Output directory already exists")
    #---------------------------------------------------------------------------
    # Create instances of the needed classes
    #---------------------------------------------------------------------------
    Skx_plots=SK_plots.Skx_plots(Skx_control['Fig_control'])
    Skx_geometry=SK_geometry.Skx_geometry()
    #---------------------------------------------------------------------------
    # Calculate the topological charge
    #---------------------------------------------------------------------------
    if Skx_control['Data_control']['topology']['execute']:
        (del_val,index_del)=Skx_geometry.triangulate_data(Skx_data.pos,        \
            Skx_control['Data_control']['Misc']['sub_index'])
        topo_charge=topological_charge(del_val,Skx_data.pos,index_del,\
            Skx_data.mom)
        del del_val, index_del
    else:
        topo_charge=[]
    #---------------------------------------------------------------------------
    # First one needs to grid the coordinate data
    #---------------------------------------------------------------------------
    x_grid=np.linspace(np.amin(Skx_data.pos[:,0]),np.amax(Skx_data.pos[:,0]),  \
        Skx_control['Data_control']['Misc']['gridpoint'][0])
    y_grid=np.linspace(np.amin(Skx_data.pos[:,1]),np.amax(Skx_data.pos[:,1]),  \
        Skx_control['Data_control']['Misc']['gridpoint'][0])
    coord = np.meshgrid(x_grid,y_grid)
    del x_grid,y_grid
    #---------------------------------------------------------------------------
    # Grid the magnetization data
    #---------------------------------------------------------------------------
    mag_grid=[]
    index_list=[2,0,1]
    for ind in index_list:
        mag_grid.append(griddata(Skx_data.pos[:,0:2],Skx_data.mom[:,ind],      \
            (coord[0],coord[1]),method='cubic'))
    #---------------------------------------------------------------------------
    # Transform the vectors to RGB colors in the HSV colorspace
    #---------------------------------------------------------------------------
    RGB_grid=SK_aux.vector_to_HSV(Skx_data.mom[:,2],Skx_data.angle_proj(),     \
                                Skx_data.pos[:,0:2],coord)
    mag_grid.append(RGB_grid)
    del RGB_grid
    #---------------------------------------------------------------------------
    # Find the center of the skyrmion
    #---------------------------------------------------------------------------
    Skx_center=Skx_geometry.find_smoothed(coord[0],coord[1],mag_grid[0],       \
                                        Skx_data.core)
    #---------------------------------------------------------------------------
    # Correcting for possible user mistake of the skyrmion threshold
    #---------------------------------------------------------------------------
    if Skx_data.core<1:
        if Skx_control['Data_control']['Skx_threshold']<0:
            Skx_control['Data_control']['Skx_threshold']=                      \
                -1.0*Skx_control['Data_control']['Skx_threshold']
        indx=np.where(mag_grid[0]<Skx_control['Data_control']['Skx_threshold'])
    else:
        if Skx_control['Data_control']['Skx_threshold']>0:
            Skx_control['Data_control']['Skx_threshold']=                      \
                -1.0*Skx_control['Data_control']['Skx_threshold']
        indx=np.where(mag_grid[0]>Skx_control['Data_control']['Skx_threshold'])
    #---------------------------------------------------------------------------
    # Find the perimeter of the skyrmion
    #---------------------------------------------------------------------------
    smooth_coord=np.zeros([len(indx[0]),3],dtype=np.float64)
    smooth_coord[:,0]=coord[0][indx[0][:],indx[1][:]]
    smooth_coord[:,1]=coord[1][indx[0][:],indx[1][:]]
    perimeter=Skx_geometry.find_skyrmion_perimeter(smooth_coord,tol=0.0001)
    del smooth_coord,indx
    #---------------------------------------------------------------------------
    # Determine the radius of the skyrmion
    #---------------------------------------------------------------------------
    (ave_dev,Skx_rad)=\
    Skx_geometry.find_skyrmion_deviation(perimeter,Skx_center[0:2])
    del perimeter
    ind_x=np.where(np.abs(Skx_data.pos[:,0]-Skx_center[0])<tol)[0]
    ind_y=np.where(np.abs(Skx_data.pos[:,1]-Skx_center[1])<tol)[0]

    bounds=\
        [Skx_center[0]-Skx_rad*(1+Skx_control['Fig_control']['zoom_factor']),  \
         Skx_center[0]+Skx_rad*(1+Skx_control['Fig_control']['zoom_factor']),  \
         Skx_center[1]-Skx_rad*(1+Skx_control['Fig_control']['zoom_factor']),  \
         Skx_center[1]+Skx_rad*(1+Skx_control['Fig_control']['zoom_factor'])]
    #---------------------------------------------------------------------------
    # Plot a single profile
    #---------------------------------------------------------------------------
    output_file=output_dir+\
        'Skx_profile'+Skx_control['Data_control']['file_name']+'.pdf'
    Skx_plots.Skx_profile(coord,bounds,Skx_data.find_bound()[:4],mag_grid[0],  \
        output_file,Skx_data.pos[ind_x,1],Skx_data.pos[ind_y,0],               \
        Skx_data.mom[ind_x,2],Skx_data.mom[ind_y,2],Skx_center,                \
        imp_pos=Skx_control['Data_control']['Misc']['imp_pos'],                \
        labels=[r'$r_x$ [nm]',r'$r_y$ [nm]',r'$m_z$'])
    #---------------------------------------------------------------------------
    # Plot a 4x4 grid of magnetization components
    #---------------------------------------------------------------------------
    output_file=output_dir+\
        'Skx_profile_grid'+Skx_control['Data_control']['file_name']+'.pdf'
    Skx_plots.SKX_profile_4x4(coord,mag_grid,output_file,                      \
        imp_pos=Skx_control['Data_control']['Misc']['imp_pos'],                \
        extent=Skx_data.find_bound()[:4],bounds=bounds,                        \
        text=[r'$m_z$',r'$m_x$',r'$m_y$',r''],                                 \
        y_labels=[r'$r_y$ [nm]','',r'$r_y$ [nm]',''],                          \
        x_labels=['','',r'$r_x$ [nm]',r'$r_x$ [nm]'],                          \
        y_ticks=[True,False,True,False],x_ticks=[False,False,True,True])

    output_file=output_dir+\
        'SK_overview'+Skx_control['Data_control']['file_name']+'.yml'
    SK_static_overview(output_file,topo_charge,Skx_center[0:2],Skx_rad,ave_dev)

    print("Profiles done!")
    del ind_x,ind_y,mag_grid,Skx_data,bounds,ave_dev,Skx_rad,Skx_center,coord
    del Skx_geometry
    return
