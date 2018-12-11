def main()
    """
    Example of how to generate a minimum input file for `APyInSky`
    """
    import yaml
    from collections import OrderedDict

    Skx_control=OrderedDict()
    Skx_control['Fig_control']=OrderedDict()
    Skx_control['Fig_control']['dpi']=800
    Skx_control['Fig_control']['height']=10
    Skx_control['Fig_control']['width']=16
    Skx_control['Fig_control']['fontsize']=28
    Skx_control['Fig_control']['zoom_factor']=0.0
    Skx_control['Fig_control']['tex_fonts']=False
    Skx_control['Data_control']=OrderedDict()
    Skx_control['Data_control']['path']='./'
    Skx_control['Data_control']['out_path']='./'
    Skx_control['Data_control']['Misc']=OrderedDict()
    Skx_control['Data_control']['Misc']['alat']=0.392
    Skx_control['Data_control']['Misc']['imp_pos']=[]
    Skx_control['Data_control']['Misc']['zipped']=False
    Skx_control['Data_control']['Misc']['sub_index']=0
    Skx_control['Data_control']['Misc']['gridpoint']=[500,500]
    Skx_control['Data_control']['Topology']=OrderedDict()
    Skx_control['Data_control']['Topology']['execute']=False
    Skx_control['Data_control']['Skyrmion_velocity']=OrderedDict()
    Skx_control['Data_control']['Skyrmion_velocity']['execute']=False
    Skx_control['Data_control']['Skyrmion_velocity']['time_step']=float(1e-16)
    Skx_control['Data_control']['Skyrmion_velocity']['plot']=False 
    Skx_control['Data_control']['Skyrmion_velocity']['plot_energy']=False
    Skx_control['Data_control']['Skyrmion_velocity']['comp_energy']=False
    Skx_control['Data_control']['profile']=OrderedDict()
    Skx_control['Data_control']['profile']['execute']=False
    Skx_control['Data_control']['Skx_threshold']=0.8
    Skx_control['Data_control']['file_name']=''

    yaml.add_representer(OrderedDict, lambda dumper, data: dumper.represent_mapping('tag:yaml.org,2002:map', data.items()))
    with open(file_name, 'w') as outfile:
        yaml.dump(SK_summary, outfile,default_flow_style=False)

if __name__ == '__main__':
    main()
