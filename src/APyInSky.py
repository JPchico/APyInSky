import SK_GNEB
import SK_profile
import SK_geometry
import SK_Speed_mom
import SK_Ene_Profile
from SK_IO_control import SK_init

def main():
    """Main wrapper handling the calls for the desired analysis. This will take care of
    reading the necessary control variables, and then calling the next layer of wrapper 
    functions for the desired type/types of analysis chosen in the control file.
    Author
    ----------
    Jonathan Chico
    """
    #---------------------------------------------------------------------------
    # Setup the input parameters
    #---------------------------------------------------------------------------
    Skx_control=SK_init()
    #---------------------------------------------------------------------------
    # Plot all the profiles
    #---------------------------------------------------------------------------
    if Skx_control['Data_control']['profile']['execute']:
        SK_profile.profile_wrapper(Skx_control)
    #---------------------------------------------------------------------------
    # Calculation of the skyrmion velocity under the influence of a current
    #---------------------------------------------------------------------------
    if Skx_control['Data_control']['Skyrmion_velocity']['execute']:
        SK_Speed_mom.Skyrmion_vel_wrapper(Skx_control)
    print("#"*80)
    print("# Skyrmion analysis done!")
    print("#"*80)
    return

if __name__ == '__main__':
    main()