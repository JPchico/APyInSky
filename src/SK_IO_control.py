def SK_init():
    """Initialization of variables to their defaults values. Reading of the inputfile to gather the necessary data for analysis.
    Returns
    ----------
        - Skx_control: (dict) dictionary containing all the needed control parameters
    Author
    ----------
    Jonathan Chico
    """
    Skx_control=set_control_defaults()
    Skx_control=input_gatherer(Skx_control)
    return Skx_control

def input_gatherer(Skx_control):
    """ Reads the input controller for the script from a yaml file which contains all the necessary input to run APyInSky. 
    The 'Skx_inp.yml' has the following information in it.
    - Data_control:
    - Fig_control:
        - dpi: (float) dpi of the figures printed.
        - height: (float) height of the figures printed.
        - width: (float) width of the figures printed.
        - fontsize: (float) size of the fonts used in the plots.
        - zoom_factor: (float) how zoomed to the skyrmion the profile is
    Returns
    ----------
        - Skx_control: (dict) dictionary containing all the needed control parameters
    Author
    ----------
    Jonathan Chico
    """
    import yaml
    import sys
    import numpy as np

    print("#"*80)
    print("# Welcome to APyInSky!")
    print("#"*80)
    try:
        with open("Skx_inp.yml", 'r') as stream:
            Skx_inp = yaml.load(stream)
    except:
        print("No 'Skx_inp.yml' file found. Shutting down")
        sys.exit()

    Skx_control=update(Skx_control,Skx_inp)
    del Skx_inp
    # Transform the list to a numpy array for easier manipulation
    try:
        Skx_control['Data_control']['Misc']['imp_pos']=\
        np.asarray(Skx_control['Data_control']['Misc']['imp_pos'])
        # Multiply by the lattice constant to be able to place it correctly in the plots
        Skx_control['Data_control']['Misc']['imp_pos']=\
        Skx_control['Data_control']['Misc']['alat']*Skx_control['Data_control']['Misc']['imp_pos']
    except:
        pass
    try:
        Skx_control['Data_control']['Skyrmion_velocity']['time_step']=\
        float(Skx_control['Data_control']['Skyrmion_velocity']['time_step'])
    except:
        Skx_control['Data_control']['Skyrmion_velocity']['time_step']=1e-16
        print('No time step given. Assuming 1e-16 s')
    try:
        if len(Skx_control['Data_control']['file_name'])>1:
            Skx_control['Data_control']['file_name']=\
            '_'+Skx_control['Data_control']['file_name']
    except:
        pass
    return Skx_control

def update(d, u):
    import collections
    """ Recursive function to update the keys of a dictionary such that if no values are given in the input the defaults are used. Small modificaions of the routine presented in
    https://stackoverflow.com/questions/3232943/update-value-of-a-nested-dictionary-of-varying-depth
    """
    for k, v in u.items():
        if isinstance(v, collections.Mapping):
            d[k] = update(d.get(k, {}), v)
        else:
            d[k] = v
    return d

def set_control_defaults():
    """
    Function to define the default values of the parameters, such that if the user does not define something the code can still be safely used.
    Author
    ----------
    Jonathan Chico
    """
    Skx_control=dict()
    Skx_control['Fig_control']=dict()
    Skx_control['Fig_control']['dpi']=800
    Skx_control['Fig_control']['height']=10
    Skx_control['Fig_control']['width']=16
    Skx_control['Fig_control']['fontsize']=28
    Skx_control['Fig_control']['zoom_factor']=0.15
    Skx_control['Fig_control']['tex_fonts']=False
    Skx_control['Data_control']=dict()
    Skx_control['Data_control']['path']='./'
    Skx_control['Data_control']['out_path']='./'
    Skx_control['Data_control']['Misc']=dict()
    Skx_control['Data_control']['Misc']['alat']=0.392
    Skx_control['Data_control']['Misc']['imp_pos']=[]
    Skx_control['Data_control']['Misc']['zipped']=False
    Skx_control['Data_control']['Misc']['sub_index']=0
    Skx_control['Data_control']['Misc']['gridpoint']=[500,500]
    Skx_control['Data_control']['Topology']=dict()
    Skx_control['Data_control']['Topology']['execute']=False
    Skx_control['Data_control']['Skyrmion_velocity']=dict()
    Skx_control['Data_control']['Skyrmion_velocity']['execute']=False
    Skx_control['Data_control']['Skyrmion_velocity']['time_step']=float(1e-16)
    Skx_control['Data_control']['Skyrmion_velocity']['plot']=False 
    Skx_control['Data_control']['Skyrmion_velocity']['plot_energy']=False
    Skx_control['Data_control']['Skyrmion_velocity']['comp_energy']=False
    Skx_control['Data_control']['profile']=dict()
    Skx_control['Data_control']['profile']['execute']=False
    Skx_control['Data_control']['Skx_threshold']=0.8
    Skx_control['Data_control']['file_name']=''
    return Skx_control

def SK_overview(file_name,pos_bare,pos_grid,pos_filt,vel_bare,vel_bare_mod,   \
    vel_grid,vel_grid_mod,vel_filt,vel_filt_mod,vel_smooth,vel_smooth_mod,     \
    rad_bare,rad_grid):
    """Function to write in a dictionary relevant quantities of the positions, velocities, and radii
    of the skyrmion.
    Args
    ----------
        - out_folder: (str) path for the output file
        - pos_bare: (float [N,3] array) positions of the skyrmion core from bare data
        - pos_grid: (float [N,3] array) positions of the skyrmion core from interpolated data
        - pos_filt: (float [N,3] array) positions of the skyrmion core from smoothed data
        - vel_bare: (float [N,3] array) velocity of the skyrmion from bare data
        - vel_bare_mod: (float [N] array) speed of the skyrmion from bare data
        - vel_grid: (float [N,3] array) velocity of the skyrmion from interpolated data
        - vel_grid_mod: (float [N] array) speed of the skyrmion from interpolated data
        - vel_filt: (float [N,3] array) velocity of the skyrmion from filtered data
        - vel_filt_mod: (float [N] array) speed of the skyrmion from filtered data
        - vel_smooth: (float [N,3] array) velocity of the skyrmion from smooth data
        - vel_smooth_mod: (float [N] array) speed of the skyrmion from smooth data
        - rad_bare: (float [N] array) radius of the skyrmion from bare data
        - rad_grid: (float [N] array) radius of the skyrmion from interpolated data
    Author
    ----------
    Jonathan Chico
    """
    import numpy as np
    import yaml
    from collections import OrderedDict

    ten_perc=int(len(vel_grid)*0.1)
    ninety_perc=int(len(vel_grid)*0.9)

    SK_summary=OrderedDict()
    #---------------------------------------------------------------------------
    # Position dictionary
    #---------------------------------------------------------------------------
    SK_summary['core_position']=OrderedDict()
    SK_summary['core_position']['method']=OrderedDict()
    #---------------------------------------------------------------------------
    # Bare positions
    #---------------------------------------------------------------------------
    SK_summary['core_position']['method']['bare']=OrderedDict()
    SK_summary['core_position']['method']['bare']['min']=[                     \
        float(np.amin(pos_bare[:,0])),\
        float(np.amin(pos_bare[:,1])),\
        float(np.amin(pos_bare[:,2]))\
    ]
    SK_summary['core_position']['method']['bare']['max']=[                     \
        float(np.amax(pos_bare[:,0])),\
        float(np.amax(pos_bare[:,1])),\
        float(np.amax(pos_bare[:,2]))\
    ]
    SK_summary['core_position']['method']['bare']['mean']=[                    \
        float(np.mean(pos_bare[:,0])),\
        float(np.mean(pos_bare[:,1])),\
        float(np.mean(pos_bare[:,2]))\
    ]
    #---------------------------------------------------------------------------
    # Grid positions
    #---------------------------------------------------------------------------
    SK_summary['core_position']['method']['grid']=OrderedDict()
    SK_summary['core_position']['method']['grid']['min']=[                     \
        float(np.amin(pos_grid[:,0])),\
        float(np.amin(pos_grid[:,1])),\
        float(np.amin(pos_grid[:,2]))\
    ]
    SK_summary['core_position']['method']['grid']['max']=[                     \
        float(np.amax(pos_grid[:,0])),\
        float(np.amax(pos_grid[:,1])),\
        float(np.amax(pos_grid[:,2]))\
    ]
    SK_summary['core_position']['method']['grid']['mean']=[                    \
        float(np.mean(pos_grid[:,0])),\
        float(np.mean(pos_grid[:,1])),\
        float(np.mean(pos_grid[:,2]))\
    ]
    #---------------------------------------------------------------------------
    # Filtered positions
    #---------------------------------------------------------------------------
    SK_summary['core_position']['method']['filt']=OrderedDict()
    SK_summary['core_position']['method']['filt']['min']=[                     \
        float(np.amin(pos_filt[:,0])),\
        float(np.amin(pos_filt[:,1])),\
        float(np.amin(pos_filt[:,2]))\
    ]
    SK_summary['core_position']['method']['filt']['max']=[                     \
        float(np.amax(pos_filt[:,0])),\
        float(np.amax(pos_filt[:,1])),\
        float(np.amax(pos_filt[:,2]))\
    ]
    SK_summary['core_position']['method']['filt']['mean']=[                    \
        float(np.mean(pos_filt[:,0])),\
        float(np.mean(pos_filt[:,1])),\
        float(np.mean(pos_filt[:,2]))\
    ]
    #---------------------------------------------------------------------------
    # Velocity dictionary
    #---------------------------------------------------------------------------
    SK_summary['vel']=OrderedDict()
    SK_summary['vel']['method']=OrderedDict()
    #---------------------------------------------------------------------------
    # Bare Velocity
    #---------------------------------------------------------------------------
    SK_summary['vel']['method']['bare']=OrderedDict()
    SK_summary['vel']['method']['bare']['min']=[                               \
        float(np.amin(vel_bare[:,0])),\
        float(np.amin(vel_bare[:,1])),\
        float(np.amin(vel_bare[:,2]))\
    ]
    SK_summary['vel']['method']['bare']['max']=[                               \
        float(np.amax(vel_bare[:,0])),\
        float(np.amax(vel_bare[:,1])),\
        float(np.amax(vel_bare[:,2]))\
    ]
    SK_summary['vel']['method']['bare']['mean']=[                              \
        float(np.mean(vel_bare[:,0])),\
        float(np.mean(vel_bare[:,1])),\
        float(np.mean(vel_bare[:,2]))\
    ]
    SK_summary['vel']['method']['bare']['t_ini']=[                             \
        float(np.mean(vel_bare[0:ten_perc,0])),\
        float(np.mean(vel_bare[0:ten_perc,1])),\
        float(np.mean(vel_bare[0:ten_perc,2]))\
    ]
    SK_summary['vel']['method']['bare']['t_fin']=[                             \
        float(np.mean(vel_bare[ninety_perc:,0])),\
        float(np.mean(vel_bare[ninety_perc:,1])),\
        float(np.mean(vel_bare[ninety_perc:,2]))\
    ]
    #---------------------------------------------------------------------------
    # Interpolated Velocity
    #---------------------------------------------------------------------------
    SK_summary['vel']['method']['grid']=OrderedDict()
    SK_summary['vel']['method']['grid']['min']=[                               \
        float(np.amin(vel_grid[:,0])),\
        float(np.amin(vel_grid[:,1])),\
        float(np.amin(vel_grid[:,2]))\
    ]
    SK_summary['vel']['method']['grid']['max']=[                               \
        float(np.amax(vel_grid[:,0])),\
        float(np.amax(vel_grid[:,1])),\
        float(np.amax(vel_grid[:,2]))\
    ]
    SK_summary['vel']['method']['grid']['mean']=[                              \
        float(np.mean(vel_grid[:,0])),\
        float(np.mean(vel_grid[:,1])),\
        float(np.mean(vel_grid[:,2]))\
    ]
    SK_summary['vel']['method']['grid']['t_ini']=[                             \
        float(np.mean(vel_grid[0:ten_perc,0])),\
        float(np.mean(vel_grid[0:ten_perc,1])),\
        float(np.mean(vel_grid[0:ten_perc,2]))\
    ]
    SK_summary['vel']['method']['grid']['t_fin']=[                             \
        float(np.mean(vel_grid[ninety_perc:,0])),\
        float(np.mean(vel_grid[ninety_perc:,1])),\
        float(np.mean(vel_grid[ninety_perc:,2]))\
    ]
    #---------------------------------------------------------------------------
    # Filtered Velocity
    #---------------------------------------------------------------------------
    SK_summary['vel']['method']['filt']=OrderedDict()
    SK_summary['vel']['method']['filt']['min']=[                               \
        float(np.amin(vel_filt[:,0])),\
        float(np.amin(vel_filt[:,1])),\
        float(np.amin(vel_filt[:,2]))\
    ]
    SK_summary['vel']['method']['filt']['max']=[                               \
        float(np.amax(vel_filt[:,0])),\
        float(np.amax(vel_filt[:,1])),\
        float(np.amax(vel_filt[:,2]))\
    ]
    SK_summary['vel']['method']['filt']['mean']=[                              \
        float(np.mean(vel_filt[:,0])),\
        float(np.mean(vel_filt[:,1])),\
        float(np.mean(vel_filt[:,2]))\
    ]
    SK_summary['vel']['method']['filt']['t_ini']=[                             \
        float(np.mean(vel_filt[0:ten_perc,0])),\
        float(np.mean(vel_filt[0:ten_perc,1])),\
        float(np.mean(vel_filt[0:ten_perc,2]))\
    ]
    SK_summary['vel']['method']['filt']['t_fin']=[                             \
        float(np.mean(vel_filt[ninety_perc:,0])),\
        float(np.mean(vel_filt[ninety_perc:,1])),\
        float(np.mean(vel_filt[ninety_perc:,2]))\
    ]
    #---------------------------------------------------------------------------
    # Smooth Velocity
    #---------------------------------------------------------------------------
    SK_summary['vel']['method']['smooth']=OrderedDict()
    SK_summary['vel']['method']['smooth']['min']=[                             \
        float(np.amin(vel_smooth[:,0])),\
        float(np.amin(vel_smooth[:,1])),\
        float(np.amin(vel_smooth[:,2]))\
    ]
    SK_summary['vel']['method']['smooth']['max']=[                             \
        float(np.amax(vel_smooth[:,0])),\
        float(np.amax(vel_smooth[:,1])),\
        float(np.amax(vel_smooth[:,2]))\
    ]
    SK_summary['vel']['method']['smooth']['mean']=[                            \
        float(np.mean(vel_smooth[:,0])),\
        float(np.mean(vel_smooth[:,1])),\
        float(np.mean(vel_smooth[:,2]))\
    ]
    SK_summary['vel']['method']['smooth']['t_ini']=[                           \
        float(np.mean(vel_smooth[0:ten_perc,0])),\
        float(np.mean(vel_smooth[0:ten_perc,1])),\
        float(np.mean(vel_smooth[0:ten_perc,2]))\
    ]
    SK_summary['vel']['method']['smooth']['t_fin']=[                           \
        float(np.mean(vel_smooth[ninety_perc:,0])),\
        float(np.mean(vel_smooth[ninety_perc:,1])),\
        float(np.mean(vel_smooth[ninety_perc:,2]))\
    ]
    #---------------------------------------------------------------------------
    # Speed dictionary
    #---------------------------------------------------------------------------
    SK_summary['speed']=OrderedDict()
    SK_summary['speed']['method']=OrderedDict()
    #---------------------------------------------------------------------------
    # Bare Speed
    #---------------------------------------------------------------------------
    SK_summary['speed']['method']['bare']=OrderedDict()
    SK_summary['speed']['method']['bare']['min']=float(np.amin(vel_bare_mod))
    SK_summary['speed']['method']['bare']['max']=float(np.amax(vel_bare_mod))
    SK_summary['speed']['method']['bare']['mean']=float(np.mean(vel_bare_mod))
    SK_summary['speed']['method']['bare']['t_ini']=                            \
        float(np.mean(vel_bare_mod[0:ten_perc]))
    SK_summary['speed']['method']['bare']['t_fin']=                            \
        float(np.mean(vel_bare_mod[ninety_perc:]))
    #---------------------------------------------------------------------------
    # Interpolated Speed
    #---------------------------------------------------------------------------
    SK_summary['speed']['method']['grid']=OrderedDict()
    SK_summary['speed']['method']['grid']['min']=float(np.amin(vel_grid_mod))
    SK_summary['speed']['method']['grid']['max']=float(np.amax(vel_grid_mod))
    SK_summary['speed']['method']['grid']['mean']=float(np.mean(vel_grid_mod))
    SK_summary['speed']['method']['grid']['t_ini']=                            \
        float(np.mean(vel_grid_mod[0:ten_perc]))
    SK_summary['speed']['method']['grid']['t_fin']=                            \
        float(np.mean(vel_grid_mod[ninety_perc:]))
    #---------------------------------------------------------------------------
    # Filtered Speed
    #---------------------------------------------------------------------------
    SK_summary['speed']['method']['filt']=OrderedDict()
    SK_summary['speed']['method']['filt']['min']=float(np.amin(vel_filt_mod))
    SK_summary['speed']['method']['filt']['max']=float(np.amax(vel_filt_mod))
    SK_summary['speed']['method']['filt']['mean']=float(np.mean(vel_filt_mod))
    SK_summary['speed']['method']['filt']['t_ini']=                            \
        float(np.mean(vel_filt_mod[0:ten_perc]))
    SK_summary['speed']['method']['filt']['t_fin']=                            \
        float(np.mean(vel_filt_mod[ninety_perc:]))
    #---------------------------------------------------------------------------
    # Smooth Speed
    #---------------------------------------------------------------------------
    SK_summary['speed']['method']['smooth']=OrderedDict()
    SK_summary['speed']['method']['smooth']['min']=                            \
        float(np.amin(vel_smooth_mod))
    SK_summary['speed']['method']['smooth']['max']=                            \
        float(np.amax(vel_smooth_mod))
    SK_summary['speed']['method']['smooth']['mean']=                           \
        float(np.mean(vel_smooth_mod))
    SK_summary['speed']['method']['smooth']['t_ini']=                          \
        float(np.mean(vel_smooth_mod[0:ten_perc]))
    SK_summary['speed']['method']['smooth']['t_fin']=                          \
        float(np.mean(vel_smooth_mod[ninety_perc:]))
    #---------------------------------------------------------------------------
    # Radius dictionary
    #---------------------------------------------------------------------------
    SK_summary['radius']=OrderedDict()
    SK_summary['radius']['method']=OrderedDict()
    #---------------------------------------------------------------------------
    # Bare Velocity
    #---------------------------------------------------------------------------
    SK_summary['radius']['method']['bare']=OrderedDict()
    SK_summary['radius']['method']['bare']['min']=float(np.amin(rad_bare))
    SK_summary['radius']['method']['bare']['max']=float(np.amax(rad_bare))
    SK_summary['radius']['method']['bare']['mean']=float(np.mean(rad_bare))
    SK_summary['radius']['method']['bare']['t_ini']=                           \
        float(np.mean(rad_bare[0:ten_perc]))
    SK_summary['radius']['method']['bare']['t_fin']=                           \
        float(np.mean(rad_bare[ninety_perc:]))
    #---------------------------------------------------------------------------
    # Interpolated Velocity
    #---------------------------------------------------------------------------
    SK_summary['radius']['method']['grid']=OrderedDict()
    SK_summary['radius']['method']['grid']['min']=float(np.amin(rad_grid))
    SK_summary['radius']['method']['grid']['max']=float(np.amax(rad_grid))
    SK_summary['radius']['method']['grid']['mean']=float(np.mean(rad_grid))
    SK_summary['radius']['method']['grid']['t_ini']=                           \
        float(np.mean(rad_grid[0:ten_perc]))
    SK_summary['radius']['method']['grid']['t_fin']=                           \
        float(np.mean(rad_grid[ninety_perc:]))

    yaml.add_representer(OrderedDict, lambda dumper, data: dumper.represent_mapping('tag:yaml.org,2002:map', data.items()))
    with open(file_name, 'w') as outfile:
        yaml.dump(SK_summary, outfile,default_flow_style=False)

    del SK_summary
    return

def SK_static_overview(file_name,topo_charge,pos,rad,deviation):
    
    import yaml
    from collections import OrderedDict

    SK_summary=OrderedDict()
    SK_summary['radius']=float(rad)
    SK_summary['topo_charge']=float(topo_charge)
    SK_summary['deviation']=float(deviation)
    SK_summary['core_pos']=[float(pos[0]),float(pos[1])]
    
    yaml.add_representer(OrderedDict, lambda dumper, data: dumper.represent_mapping('tag:yaml.org,2002:map', data.items()))
    with open(file_name, 'w') as outfile:
        yaml.dump(SK_summary, outfile,default_flow_style=False)

    del SK_summary
    return

def file_handler(path,zipped,file_prefix,file_suffix='.out'):
    """File wrapper to read, untar and extract the necessary data for the given function that is being analysed.
    Args
    ----------
        - path: (str) path to the desired folder
        - zipped: (boolean) indicating if the data is in a tar.gz
        - file_prefix: (str) prefix of the file that one wishes to extract data from.
            * coord: coordinates of all the atoms in the simulation box
            * restart: magnetic configuration of all the atoms in the simulation box.
            * moment: time evolution of the magnetic configuration for all the atoms in the simulation box.
            * posfile_clus: position of the atoms in the cluster.
            * totenergy: total energy per atom in the system.
            * rate_if: barriers and pre-factors from the GNEB calculation.
            * enfit_path: fitted energy landscape from GENB.
            * localenergy: site dependent energy
        - file_suffix: (str) suffix of the file to discriminate between input and output files that have the same prefix (default='out')
    Returns
    ----------
        - data: numpy array containing the data extracted from the file.
    Author
    ----------
    Jonathan Chico
    """
    import glob
    import tarfile
    import pandas as pd
    import numpy as np
    # If the file is in a tar.gz special care must be paid to it
    if zipped:
        # Find the name of the tar.gz file and pick the first one found, as there should only be one
        file_name=glob.glob(path+'/'+"*.tar.gz")[0]
        tar_data=tarfile.open(file_name)
        names=tar_data.getnames()
        members=tar_data.getmembers()
        # Find the first file which corresponds to the pattern that one is searching for
        tmp_name=pd.Series(names)
        target_file=tmp_name[tmp_name.str.match(file_prefix)].values
        # Make sure the suffix is included too
        target_file = [s for s in target_file if file_suffix in s]
        ind=names.index(target_file[0])
        data_file = tar_data.extractfile(members[ind])
    else:
       data_file=glob.glob(path+'/'+file_prefix+'*'+file_suffix)[0] 
    (skiprows,columns)=set_format(data_file,file_prefix,zipped)
    # Extract the needed data 
    data=pd.read_csv(data_file,header=None,delim_whitespace=True,             \
        usecols=columns,skiprows=skiprows).values
    return data

def set_format(data_file,file_prefix,zipped):
    """Function to determine if the UppASD data is in the new or legacy format. Determines the number of columns that are used and the number of rows that must be skipped.
    Args
    ----------
        - data_file: file identifier which is currently being read.
        - file_preffix: (str) name of file which is being read.
        - zipped: (boolean) indicating if the data is in a tar.gz
    Returns
    ----------
        - skiprows: (int) number of lines to be skipped when read
        - columns: (list) list with the number of the columns that are read for each type of file
    Author
    ----------
    Jonathan Chico
    """
    import numpy as np
    import pandas as pd

    line = pd.read_csv(data_file,header=None,delim_whitespace=True,nrows=1,    \
        usecols=[0]).values
    try:
        data_file.seek(0)
    except:
        pass
    #---------------------------------------------------------------------------
    # Check whether the file is in the new or old fmt
    #---------------------------------------------------------------------------
    try:
        comp=str(line[0])
    except:
        comp=str(line)

    if comp=='#':
        type_fmt='new'
    else:
        type_fmt='old'
    if file_prefix=='coord':
        columns=[1,2,3]
        skiprows=0
    if file_prefix=='restart':
        if type_fmt=='new':
            columns=[4,5,6]
            skiprows=7
        else:
            columns=[3,4,5]
            skiprows=1
    if file_prefix=='moment':
        if type_fmt=='new':
            columns=[0,4,5,6]
            skiprows=7
        else:
            columns=[0,2,3,4]
            skiprows=0
    if file_prefix=='posfile_clus':
        columns=[2,3,4]
        skiprows=0
    if file_prefix=='totenergy':
        columns=[1]
        skiprows=1
    if file_prefix=='rate_if':
        columns=[0,1]
        skiprows=1
    if file_prefix=='enfit_path':
        columns=[1,2]
        skiprows=1
    if file_prefix=='localenergy':
        skiprows=1
        columns=[3,4,6]
    return skiprows,columns