import numpy as np
import matplotlib.pyplot as plt

class Skx_plots:
    """
    Class containing the necessary methods for the creation of plots of the analyzed skyrmions.
    It contains the capacity of making:
    - Stacked plots for quick visualization of vectorial data. 
    - Linear plots for visualization of individual data.
    - Profiles with marginals.
    - 4x4 Grids for different components of a vectorial quantity.
    Args
    ----------
    Fig_control: (dict) dictionary containing all the options for the plotting of figures.
    Author
    ----------
    Jonathan Chico
    """
    def __init__(self,Fig_control):

        self.symbol_size=max(Fig_control['width'],Fig_control['height'])*10
        self.label_font_size=Fig_control['fontsize']
        self.aux_font_size=Fig_control['fontsize']*0.75
        self.fig_size=(Fig_control['width'],Fig_control['height'])
        self.dpi=Fig_control['dpi']

        plt.rc('text', usetex=Fig_control['tex_fonts'])
        plt.rc('font', family='serif')
        return

    def gridded_plot(self,x_data,y_data,line_array,label_array,fig_name,       \
        y_labels):
        """
        Stacked figure to plot multiple subfigures which share a same x-axis.
        It is aimed at the plottin of vectorial quantities, as one can plot each
        component as as individual subplot.
        The method also handles the fact that several actors might be plotted in each of the subfigures.
        Args
        ----------
            - x_data: (1-D array) values of the x-axis. They are assumed to be the same for all plots.
            - y_data: (N-D array) values for the y-components for each of the subplots.
            - line_array: (1-D array) linestyles for each of the actors in the subplots.
            - label_array: (1-D str array) labels for the legends of each actor in the subplot
            - fig_name: (str) name to save the plot
            - y_labels: (1-D str array) labels for each of the actors in the subplots.
        Author
        ----------
        Jonathan Chico
        """

        from matplotlib import cm as cm
    
        colors = cm.Paired(np.linspace(0,1,len(y_data[0])))

        fig,grid = plt.subplots(len(y_data),sharex=True,figsize=self.fig_size)

        for ii in range(len(y_data)):
            for jj in range(len(y_data[ii])):
                if label_array[jj]!='Bare':
                    grid[ii].plot(x_data,y_data[ii][jj],line_array[ii],        \
                        color=colors[jj],alpha=0.75,lw=5,label=label_array[jj])
                else:
                    grid[ii].scatter(x_data,y_data[ii][jj],color=colors[jj],   \
                        alpha=0.75,s=self.symbol_size,lw=0.75,                 \
                        edgecolor='black',label=label_array[jj])
                grid[ii].set_ylabel(y_labels[ii],fontsize=self.label_font_size)
                grid[ii].tick_params(axis='x', colors='black',                 \
                    labelsize=self.aux_font_size,width=2)
                grid[ii].tick_params(axis='y', colors='black',                 \
                    labelsize=self.aux_font_size,width=2)
                #---------------------------------------------------------------
                # Plotting axis information
                #---------------------------------------------------------------
                for axis in ['top','bottom','left','right']:
                    grid[ii].spines[axis].set_linewidth(3)
        #-----------------------------------------------------------------------
        # Axis information of the y-direction of the Velocity
        #-----------------------------------------------------------------------
        grid[-1].set_xlabel('Time [ns]',fontsize=self.label_font_size)
        #-----------------------------------------------------------------------
        # Plotting the legends
        #-----------------------------------------------------------------------
        grid[1].legend(fontsize=self.aux_font_size,loc='center left',          \
            bbox_to_anchor=(1,0.5))
        plt.grid(False)
        fig.subplots_adjust(hspace=0)
        #-----------------------------------------------------------------------
        ## Figure parameters
        #-----------------------------------------------------------------------
        fig.savefig(fig_name,transparent=True,dpi=self.dpi,bbox_inches='tight')
        plt.gcf().clear()
        plt.cla()
        plt.close(fig)

        return

    def line_plot(self,x_data,y_data,y_label,x_label,label_array,line_array,   \
        fig_name):
        from matplotlib import cm as cm

        #-----------------------------------------------------------------------
        # Create the figure
        #-----------------------------------------------------------------------
        fig = plt.figure(figsize=self.fig_size)
        ax = fig.add_subplot(111)
        #-----------------------------------------------------------------------
        # Check the dimensions of the input data to see if it is a single line 
        # or a set of plots
        #-----------------------------------------------------------------------
        if len(np.shape(y_data))>1:
            colors = cm.Paired(np.linspace(0,1,len(y_data)))
            max_val=np.amax([np.amax(list_ent) for list_ent in y_data])
            min_val=np.amin([np.amin(list_ent) for list_ent in y_data])
            for ii in range(len(y_data)):
                if label_array[ii]!='Bare':
                    ax.plot(x_data,y_data[ii],line_array[ii],color=colors[ii], \
                        alpha=0.75,lw=5,label=label_array[ii])
                else:
                    ax.scatter(x_data,y_data[ii],color=colors[ii],alpha=0.75,  \
                        s=self.symbol_size,lw=0.75,edgecolor='black',          \
                        label=label_array[ii])
        else:
            max_val=np.amax(y_data)
            min_val=np.amin(y_data)
            colors = 'black' 
            if label_array!='Bare':
                ax.plot(x_data,y_data,line_array,color=colors,alpha=0.75,lw=5, \
                    label=label_array)
            else:
                ax.scatter(x_data,y_data,color=colors,alpha=0.75,              \
                    s=self.symbol_size,lw=0.75,edgecolor='black',              \
                    label=label_array)

        plt.legend(fontsize=self.aux_font_size)
        #-----------------------------------------------------------------------
        # Plotting axis information
        #-----------------------------------------------------------------------
        plt.xlabel(x_label,fontsize=self.label_font_size)
        plt.ylabel(y_label,fontsize=self.label_font_size)
        plt.ylim(min_val*0.9,max_val*1.1)
        ax.tick_params(axis='x',colors='black',labelsize=self.aux_font_size,   \
            width=2)
        ax.tick_params(axis='y',colors='black',labelsize=self.aux_font_size,   \
            width=2)
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(3)
        plt.grid(False)
        #-----------------------------------------------------------------------
        # Figure parameters
        #-----------------------------------------------------------------------
        fig.savefig(fig_name,transparent=True,dpi=self.dpi,bbox_inches='tight')
        plt.gcf().clear()
        plt.cla()
        plt.close(fig)
        return

    def create_time_plot(self,step_size,num_times,bounds):
        from matplotlib import colors, colorbar
        from matplotlib import cm as cm
        from matplotlib.gridspec import GridSpec

        cmap=cm.coolwarm
        self.fig_traj = plt.figure(figsize=self.fig_size)
        grid=GridSpec(1,2,figure=self.fig_traj,height_ratios=[1],              \
            width_ratios=[1,0.05])
        self.fig_traj.tight_layout(pad=0.0)
        grid.update(left=0.05, right=0.95, bottom=0.10, top=0.875, wspace=0.0, \
            hspace=0.0)
        self.ax_traj= plt.subplot(grid[0,0])
        self.ax_cbar= plt.subplot(grid[0,1])
        ########################################################################
        # Setting the x-axis parameters
        ########################################################################
        self.ax_traj.set_xlabel(r'$R_x$ [nm]',fontsize=self.label_font_size)
        self.ax_traj.tick_params(axis='x',colors='black',                      \
            labelsize=self.aux_font_size,width=2)
        ########################################################################
        # Setting the y-axis parameters
        ########################################################################
        self.ax_traj.set_ylabel(r'$R_y$ [nm]',fontsize=self.label_font_size)
        self.ax_traj.tick_params(axis='y',colors='black',                      \
            labelsize=self.aux_font_size,width=2)
        ########################################################################
        # Axis parameters
        ########################################################################
        self.ax_traj.axis(bounds)
        for axis in ['top','bottom','left','right']:
            self.ax_traj.spines[axis].set_linewidth(3)
        ########################################################################
        # Setting the cbar-axis parameters
        ########################################################################
        norm = colors.Normalize(vmin=0,vmax=step_size*num_times*1e9)
        cbar = colorbar.ColorbarBase(self.ax_cbar,cmap=cmap,norm=norm,         \
            spacing='uniform',orientation='vertical')
        cbar.set_label(r'Time [ns]',fontsize=self.label_font_size)
        cbar.ax.tick_params(labelsize=self.aux_font_size,width=2,colors='black')
        cbar.outline.set_linewidth(3)

        return

    ############################################################################
    # Function to plot the skyrmion profile as a function of time
    ############################################################################
    def plot_time_profile(self,x_grid,y_grid,mag_grid,itr,num_times,pos):
        """ Function to plot the time evolution of the skyrmion profile.
        Args
        ----------
            - Fig_control: (dict) dictionary controling the figure properties
            - pos: (float [N,3] array) position of the core of the skyrmion
        Author
        ----------
        Jonathan Chico
        """
        from matplotlib import cm as cm
        from matplotlib.colors import to_rgba_array

        color_set = cm.coolwarm(np.linspace(0,1,num_times))
        ########################################################################
        # Plot the skyrmion profile for the current time step
        ########################################################################
        self.ax_traj.contour(x_grid,y_grid,mag_grid,[0.0],linewidths=0.5,      \
            linestyles='solid',colors=to_rgba_array(color_set[itr]))
        if len(pos)>0:
            self.ax_traj.scatter(pos[0],pos[1],\
                color=to_rgba_array(color_set[itr]), \
                s=self.symbol_size*0.25,lw=0.01,edgecolor='black',alpha=0.5,   \
                zorder=9)
        self.ax_traj.set_aspect('equal')
        self.fig_traj.canvas.draw()

        return

    def save_time_plot(self,fig_name,imp_pos=[]):
        """Function to save the time evolution of the skyrmion shell and core positions.
        """
        if len(imp_pos)>0:
            self.ax_traj.scatter(imp_pos[0],imp_pos[1],color='#525252',\
                alpha=0.5,s=self.symbol_size,lw=1.00,edgecolor='black',\
                zorder=10)
        self.fig_traj.savefig(fig_name,transparent=True,dpi=self.dpi,          \
            bbox_inches='tight')
        self.fig_traj.clear()
        self.ax_traj.cla()
        plt.close(self.fig_traj)
        return

    def Skx_profile(self,coord,bounds,extent,mag,output_file,x_cut,y_cut,      \
        mag_x_cut,mag_y_cut,skx_cent=[],imp_pos=[],                            \
        labels=[r'$R_x$ [nm]',r'$R_y$ [nm]',r'$M_z$']):
        from matplotlib import cm as cm
        from matplotlib.ticker import NullFormatter
        from matplotlib.gridspec import GridSpec
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        #-----------------------------------------------------------------------
        # Create the figures and the axis 
        #-----------------------------------------------------------------------
        curr_size=max(self.fig_size)
        fig_prof = plt.figure(figsize=(curr_size,curr_size))
        grid = GridSpec(4,4,figure=fig_prof,hspace=0.0,wspace=0.0)
        ax_prof = fig_prof.add_subplot(grid[:-1, 1:])
        ax_profy = fig_prof.add_subplot(grid[:-1, 0], sharey=ax_prof)
        ax_profx = fig_prof.add_subplot(grid[-1, 1:], sharex=ax_prof)

        #-----------------------------------------------------------------------
        # no labels
        #-----------------------------------------------------------------------
        ax_profx.yaxis.set_major_formatter(NullFormatter())
        ax_profy.xaxis.set_major_formatter(NullFormatter())
        #-----------------------------------------------------------------------
        # Color map of the skyrmion profile
        #-----------------------------------------------------------------------
        prof = ax_prof.imshow(mag,interpolation='nearest',cmap=cm.coolwarm,    \
            extent=extent,aspect='equal',origin='lower')
        ax_profx.plot(x_cut,mag_x_cut,lw=5)
        ax_profy.plot(mag_y_cut,y_cut,lw=5)
        if len(skx_cent)>0:
            ax_profx.axvline(skx_cent[0],lw=2,linestyle='--',color='black')
            ax_profy.axhline(skx_cent[1],lw=2,linestyle='--',color='black')
        #-----------------------------------------------------------------------
        # Contour for the m=0 line
        #-----------------------------------------------------------------------
        ax_prof.contour(coord[0],coord[1],mag,levels=[0.0],linewidths=3,      \
            linestyles='dashed',colors='black')
        #-----------------------------------------------------------------------
        # indication of the impurity location
        #-----------------------------------------------------------------------
        if len(imp_pos)>0:
            ax_prof.scatter(imp_pos[0],imp_pos[1],color='#ffffff',alpha=0.75,  \
                s=self.symbol_size,lw=3.00,edgecolor='black')
        #-----------------------------------------------------------------------
        # Set the axis boundaries
        #-----------------------------------------------------------------------
        ax_prof.axis(bounds)
        ax_profx.axis([bounds[0],bounds[1],-1.05,1.05])
        ax_profy.axis([-1.05,1.05,bounds[2],bounds[3]])
        #-----------------------------------------------------------------------
        # Plot a vertical line from the center of the skyrmion
        #-----------------------------------------------------------------------
        ax_prof.axvline(skx_cent[0],ymin=0,ymax=0.5,zorder=20,lw=2,            \
            linestyle='--',color='black')
        #-----------------------------------------------------------------------
        # Plot a horizontal line from the center of the skyrmion
        #-----------------------------------------------------------------------
        ax_prof.axhline(skx_cent[1],xmin=0,xmax=0.5,zorder=20,lw=2,            \
            linestyle='--',color='black')
        #-----------------------------------------------------------------------
        # Set the axis parameters for the full profile
        #-----------------------------------------------------------------------
        ax_prof.tick_params(axis='y',which='both',left=False,right=False,      \
            labelleft=False,width=2)
        ax_prof.tick_params(axis='x',which='both',top=False,bottom=False,      \
            labelbottom=False,width=2)
        #-----------------------------------------------------------------------
        # Set the axis parameters for the x-cut of the profile
        #-----------------------------------------------------------------------
        ax_profx.set_xlabel(labels[0],fontsize=self.label_font_size)
        ax_profx.tick_params(axis='y',which='both',left=False,right=False,     \
            labelleft=False,width=2)
        ax_profx.tick_params(axis='x',colors='black',                          \
            labelsize=self.aux_font_size,width=2)
        #-----------------------------------------------------------------------
        # Set the axis parameters for the y-cut of the profile
        #-----------------------------------------------------------------------
        ax_profy.set_ylabel(labels[1],fontsize=self.label_font_size)
        ax_profy.tick_params(axis='x',which='both',bottom=False,top=False,     \
            labelbottom=False,width=2)
        ax_profy.tick_params(axis='y',colors='black',                          \
            labelsize=self.aux_font_size,width=2)
        #-----------------------------------------------------------------------
        # Set the width of the lines in the plots
        #-----------------------------------------------------------------------
        for axis in ['top','bottom','left','right']:
            ax_profx.spines[axis].set_linewidth(3)
        for axis in ['top','bottom','left','right']:
            ax_profy.spines[axis].set_linewidth(3)
        for axis in ['top','right','bottom','left']:
            ax_prof.spines[axis].set_linewidth(3)
        # We change the fontsize of minor ticks label
        cbaxes = inset_axes(ax_prof,width="30%",height="3%",loc=3,borderpad=1.0)
        cbar=plt.colorbar(prof,cax=cbaxes,orientation='horizontal',ticks=[-1,1])
        cbar.set_label(labels[2],fontsize=self.label_font_size,labelpad=-2)
        cbar.ax.xaxis.set_ticks_position('top')
        cbar.ax.xaxis.set_label_position('top')
        cbar.ax.tick_params(labelsize=self.aux_font_size,width=2,colors='black')
        cbar.outline.set_linewidth(3)
        #-----------------------------------------------------------------------
        # Save the plot to file
        #-----------------------------------------------------------------------
        fig_prof.savefig(output_file,transparent=True,dpi=self.dpi,            \
            bbox_inches='tight')
        plt.close(fig_prof)

        return

    def SKX_profile_4x4(self,coord,mag,output_file,imp_pos,extent,bounds,text, \
        y_labels,x_labels,y_ticks,x_ticks):
        from matplotlib import cm as cm
        from matplotlib.gridspec import GridSpec
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        from matplotlib.projections import get_projection_class
        from SK_aux import color_wheel
        #-----------------------------------------------------------------------
        # Create the figure and the axes
        #-----------------------------------------------------------------------
        curr_size=max(self.fig_size)
        fig_prof = plt.figure(figsize=(curr_size,curr_size))
        gs = GridSpec(4,4,figure=fig_prof,hspace=0.0,wspace=0.0)
        grid=[]
        grid.append(fig_prof.add_subplot(gs[0,0]))
        grid.append(fig_prof.add_subplot(gs[0,1],sharex=grid[0],sharey=grid[0]))
        grid.append(fig_prof.add_subplot(gs[1,0],sharey=grid[0],sharex=grid[0]))
        grid.append(fig_prof.add_subplot(gs[1,1],sharey=grid[1],sharex=grid[0]))
        #-----------------------------------------------------------------------
        # List of which axis should have thick lines
        #-----------------------------------------------------------------------
        ax_list=[['top','bottom','left','right'],['top','bottom','right'],['bottom','left','right'],['bottom','right']]
        #-----------------------------------------------------------------------
        # Loop over the plots
        #-----------------------------------------------------------------------
        for ii in range(0,4):
            #-------------------------------------------------------------------
            # Plot the x,y and z components of the magnetization
            #-------------------------------------------------------------------
            im=grid[ii].imshow(mag[ii],interpolation='lanczos',\
                    cmap=cm.coolwarm,extent=extent,aspect='equal',origin='lower')
            if ii<3:
                #---------------------------------------------------------------
                # Define the axis where the color bars are going to plotted
                #---------------------------------------------------------------
                cbaxes = inset_axes(grid[ii],width="30%",height="3%",loc=3,\
                    borderpad=1.0) 
                #---------------------------------------------------------------
                # Color bar properties
                #---------------------------------------------------------------
                cbar=plt.colorbar(im,cax=cbaxes,orientation='horizontal',      \
                    ticks=[-1,1])
                cbar.set_label(text[ii],fontsize=self.label_font_size*0.5,\
                    labelpad=-4)
                cbar.ax.xaxis.set_ticks_position('top')
                cbar.ax.xaxis.set_label_position('top')
                cbar.ax.tick_params(labelsize=self.aux_font_size*0.5,width=2,  \
                    colors='black')
                cbar.outline.set_linewidth(1.5)
            else:
                #---------------------------------------------------------------
                # Define the axis where the color bars are going to plotted
                #---------------------------------------------------------------
                cbaxes = inset_axes(grid[ii],width="20%",height="20%",loc=3,   \
                    axes_class=get_projection_class("polar"),borderpad=1.0) 
                #---------------------------------------------------------------
                # Set the colorbar as a color wheel
                #---------------------------------------------------------------
                color_wheel(cbaxes,[0,np.pi*0.5],[r'$m_x$',r'$m_y$'],          \
                            self.aux_font_size*0.5)
            #-------------------------------------------------------------------
            # Setting the parameters for the axis
            #-------------------------------------------------------------------
            grid[ii].tick_params(axis='y',which='both',left=y_ticks[ii],       \
                right=False,labelleft=y_ticks[ii],                             \
                labelsize=self.aux_font_size*0.5,width=2)
            grid[ii].tick_params(axis='x',which='both',top=False,              \
                bottom=x_ticks[ii],labelbottom=x_ticks[ii],                    \
                labelsize=self.aux_font_size*0.5,width=2)
            #-------------------------------------------------------------------
            # Plot an impurity if it is present
            #-------------------------------------------------------------------
            if len(imp_pos)>0:
                grid[ii].scatter(imp_pos[0],imp_pos[1],color='#ffffff',\
                    alpha=0.75,s=self.symbol_size,lw=3.00,edgecolor='black')
            #-------------------------------------------------------------------
            # Set the limits to the axis
            #-------------------------------------------------------------------
            grid[ii].axis(bounds)
            #-------------------------------------------------------------------
            # Set the linewidths for the axis
            #-------------------------------------------------------------------
            for axis in ['top','bottom','left','right']:
                grid[ii].spines[axis].set_linewidth(0)
            for axis in ax_list[ii]:
                grid[ii].spines[axis].set_linewidth(3)
            #-------------------------------------------------------------------
            # Set the axis labels
            #-------------------------------------------------------------------
            grid[ii].set_xlabel(x_labels[ii],fontsize=self.label_font_size*0.5)
            grid[ii].set_ylabel(y_labels[ii],fontsize=self.label_font_size*0.5)
        #-----------------------------------------------------------------------
        # Plot the contour at m=0
        #-----------------------------------------------------------------------
        grid[0].contour(coord[0],coord[1],mag[0],[0.0],linewidths=3,           \
            linestyles='dashed',colors='black')
        #-----------------------------------------------------------------------
        # Save the plot to file
        #-----------------------------------------------------------------------
        plt.savefig(output_file,transparent=True,dpi=self.dpi,\
            bbox_inches='tight')

        plt.close(fig_prof)

        return

class Force_Animation:
    """ Class dealing with the animation of the magnetic forces as a function of time. It displays the magnetic forces acting over a skyrmion as it traverses the sample. 
    It has the capacity to display both the total force only, as well as a decomposed visualization, where the exchange and DMI contributions are explicitly displayed.
    Calling this class will initialize several plots and axes which are defined to handle the visualization of these animations.
    Args:
    ----------
        - size: (int [2] array, optional) dimensions of the figure that will be animated.
        - font_size: (int, optional) fontsize used for the figures.
        - extent: (float [4] array, optional) size of the initial dummy imshow plot.
        - mesh: (float [2,N,N] array) interpolated positions of the atoms.
        - comp_energy: (boolean, optional) logical indicating if the decomposed energy will be plotted or not.
        - imp_pos: (float [2] array, optional) position of the impurity if there is one present.
    Author:
    ----------
    Jonathan Chico

    """
    def __init__(self,size=(16,10),font_size=28,extent=(0,1,0,1),mesh=([],[]), \
        comp_energy=False,imp_pos=[]):
        """ Class constructor for the animation of the forces.
        """
        from matplotlib.gridspec import GridSpec
        if comp_energy:
            #-------------------------------------------------------------------
            # Creating the plots and the axis
            #-------------------------------------------------------------------
            self.fig_force = plt.figure(figsize=size)
            gs=GridSpec(1,3,figure=self.fig_force,height_ratios=[1], width_ratios=[1,1,1])
            self.fig_force.tight_layout(pad=0.0)
            gs.update(left=0.075, right=0.925, bottom=0.00, top=1.00, wspace=0.0, hspace=0.0)
            # Defining the total force plot
            self.ax_force=plt.subplot(gs[0,0])
            # Setting a title for the total force plot
            self.ax_force.set_title("Total Force",fontsize=font_size)
            # Defining the total exchange plot
            self.ax_force_xc=plt.subplot(gs[0,1],sharey=self.ax_force)
            # Setting a title for the exchange force plot
            self.ax_force_xc.set_title("Exchange Force",fontsize=font_size)
            # Defining the DMI force plot
            self.ax_force_dm=plt.subplot(gs[0,2],sharey=self.ax_force)
            # Setting a title for the DMI force plot
            self.ax_force_dm.set_title("DMI Force",fontsize=font_size)
            #-------------------------------------------------------------------
            # Defining the axis properties
            #-------------------------------------------------------------------
            # Axis properties for the total force plot
            self.ax_force.set_xlabel(r'$R_x$ [nm]',fontsize=font_size)
            self.ax_force.set_ylabel(r'$R_y$ [nm]',fontsize=font_size)
            self.ax_force.tick_params(axis='x',colors='black',                 \
                labelsize=font_size*0.75,width=2)
            self.ax_force.tick_params(axis='y',colors='black',                 \
                labelsize=font_size*0.75,width=2)
            # Thickness of the axis lines
            for axis in ['top','bottom','left','right']:
                self.ax_force.spines[axis].set_linewidth(3)
            # Axis properties for the exchange force plot
            self.ax_force_xc.set_xlabel(r'$R_x$ [nm]',fontsize=font_size)
            self.ax_force_xc.tick_params(axis='x',colors='black',              \
                labelsize=font_size*0.75,width=2)
            # Removing the y-ticks for the internal plots
            self.ax_force_xc.yaxis.set_ticks_position('none')
            plt.setp(self.ax_force_xc.get_yticklabels(), visible=False)
            # Thickness of the axis lines
            for axis in ['top','bottom','left','right']:
                self.ax_force_xc.spines[axis].set_linewidth(3)
            # Axis properties for the DMI force plot
            self.ax_force_dm.set_xlabel(r'$R_x$ [nm]',fontsize=font_size)
            self.ax_force_dm.tick_params(axis='x',colors='black',              \
                labelsize=font_size*0.75,width=2)
            # Removing the y-ticks for the internal plots
            self.ax_force_dm.yaxis.set_ticks_position('none')
            plt.setp(self.ax_force_dm.get_yticklabels(), visible=False)
            # Thickness of the axis lines
            for axis in ['top','bottom','left','right']:
                self.ax_force_dm.spines[axis].set_linewidth(3)
            #-------------------------------------------------------------------
            # Creating the plots
            #-------------------------------------------------------------------
            # Skyrmion profile background for the total force
            self.mag_plot=self.ax_force.imshow([[],[]],cmap=plt.cm.coolwarm,   \
                aspect='equal',origin='lower',extent=extent)
            # Skyrmion profile background for the exchange force
            self.mag_plot_xc=self.ax_force_xc.imshow([[],[]],                  \
                cmap=plt.cm.coolwarm,aspect='equal',origin='lower',            \
                extent=extent)
            # Skyrmion profile background for the DMI force
            self.mag_plot_dm=self.ax_force_dm.imshow([[],[]],                  \
                cmap=plt.cm.coolwarm,aspect='equal',origin='lower',            \
                extent=extent)
            # Defining the stride for the quiver plot to have a less dense plot
            every_step=int(np.shape(mesh[0])[0]*0.005)
            self.skip=(slice(None,None,every_step),slice(None,None,every_step))
            # Plot for the total force
            self.force_plot=self.ax_force.quiver(mesh[0][self.skip],           \
                mesh[1][self.skip],[],[],pivot='mid',scale=0.275,width=0.005)
            # Plot for the exchange force
            self.force_plot_xc=self.ax_force_xc.quiver(mesh[0][self.skip],     \
                mesh[1][self.skip],[],[],pivot='mid',scale=0.275,width=0.005)
            # Plot for the DMI force
            self.force_plot_dm=self.ax_force_dm.quiver(mesh[0][self.skip],     \
                mesh[1][self.skip],[],[],pivot='mid',scale=0.275,width=0.005)
            if len(imp_pos)>0:
                symbol_size=max(size)*10
                self.ax_force.scatter(imp_pos[0],imp_pos[1],color='#525252',\
                    alpha=0.5,s=symbol_size,lw=1.00,edgecolor='black',zorder=10)
                self.ax_force_xc.scatter(imp_pos[0],imp_pos[1],color='#525252',\
                    alpha=0.5,s=symbol_size,lw=1.00,edgecolor='black',zorder=10)
                self.ax_force_dm.scatter(imp_pos[0],imp_pos[1],color='#525252',\
                    alpha=0.5,s=symbol_size,lw=1.00,edgecolor='black',zorder=10)
        else:
            #-------------------------------------------------------------------
            # Creating the plot and the axis
            #-------------------------------------------------------------------
            self.fig_force = plt.figure(figsize=size)
            gs=GridSpec(1,1,figure=self.fig_force,height_ratios=[1],           \
                width_ratios=[1])
            self.fig_force.tight_layout(pad=0.0)
            gs.update(left=0.05, right=0.95, bottom=0.10, top=0.875,           \
                wspace=0.0, hspace=0.0)
            self.ax_force=plt.subplot(gs[0,0])
            #-------------------------------------------------------------------
            # Defining the axis properties
            #-------------------------------------------------------------------
            self.ax_force.set_xlabel(r'$R_x$ [nm]',fontsize=font_size)
            self.ax_force.set_ylabel(r'$R_x$ [nm]',fontsize=font_size)
            self.ax_force.tick_params(axis='x',colors='black',                 \
                labelsize=font_size*0.75,width=2)
            self.ax_force.tick_params(axis='y',colors='black',                 \
                labelsize=font_size*0.75,width=2)
            for axis in ['top','bottom','left','right']:
                self.ax_force.spines[axis].set_linewidth(3)
            #-------------------------------------------------------------------
            # Creating the plots
            #-------------------------------------------------------------------
            symbol_size=max(size)*10
            self.mag_plot=self.ax_force.imshow([[],[]],cmap=plt.cm.coolwarm,   \
                aspect='equal',origin='lower',extent=extent)
            every_step=int(np.shape(mesh[0])[0]*0.005)
            self.skip=(slice(None,None,every_step),slice(None,None,every_step))
            self.force_plot=self.ax_force.quiver(mesh[0][self.skip],          \
                mesh[1][self.skip],[],[],pivot='mid',scale=0.3)
            if len(imp_pos)>0:
                self.ax_force.scatter(imp_pos[0],imp_pos[1],color='#525252',   \
                    alpha=0.5,s=symbol_size,lw=1.00,edgecolor='black',zorder=10)

    def update_force_plot(self,itr):
        """Function to update the data dealing the visualization of the magnetic configuration and the forces.
        The magnetic configuration is displayed via an `imshow` 2D map of the z-component of the magnetization. The local forces are shown via a `quiver` plot.
        Args:
        ----------
            - itr: (int) present iteration which is being visualized
        Returns
        ----------
            - self.force_plot: instance of the matplotlib plot. Notice it is necessary to provide a return even if the update will be done anyhow, as it is expected by the animation function.
        Author:
        ----------
        Jonathan Chico
        """
        from SK_aux import progress_bar
        #-----------------------------------------------------------------------
        # Updating the magnetic configuration
        #-----------------------------------------------------------------------
        self.mag_plot.set_array(self.mag_data[itr])
        #-----------------------------------------------------------------------
        # Updating the forces
        #-----------------------------------------------------------------------
        self.force_plot.set_UVC(self.force_data[itr][0][self.skip],            \
            self.force_data[itr][1][self.skip])
        #-----------------------------------------------------------------------
        # Update the axis
        #-----------------------------------------------------------------------
        self.ax_force.axis((self.curr_axis[0][itr],self.curr_axis[1][itr],\
                            self.curr_axis[2][itr],self.curr_axis[3][itr]))
        #-----------------------------------------------------------------------
        # Display the progress bar to time how long the rendering takes
        #-----------------------------------------------------------------------
        progress_bar(itr,self.max_iter,'Generating forces video ')
        return self.force_plot,

    def update_comp_force_plot(self,itr):
        """Function to update the data dealing the visualization of the magnetic configuration and the forces.
        The magnetic configuration is displayed via an `imshow` 2D map of the z-component of the magnetization. The local forces are shown via a `quiver` plot. In this case the forces will be decomposed and plotted as individual figures.
        Args:
        ----------
            - itr: (int) present iteration which is being visualized
        Returns
        ----------
            - self.force_plot: instance of the matplotlib plot. Notice it is necessary to provide a return even if the update will be done anyhow, as it is expected by the animation function.
        Author:
        ----------
        Jonathan Chico
        """
        from SK_aux import progress_bar
        #-----------------------------------------------------------------------
        # Updating the magnetic configuration
        #-----------------------------------------------------------------------
        self.mag_plot.set_array(self.mag_data[itr])
        self.mag_plot_xc.set_array(self.mag_data[itr])
        self.mag_plot_dm.set_array(self.mag_data[itr])
        #-----------------------------------------------------------------------
        # Updating the forces
        #-----------------------------------------------------------------------
        self.force_plot.set_UVC(self.force_data[itr][0][self.skip],            \
            self.force_data[itr][1][self.skip])
        self.force_plot_xc.set_UVC(self.force_data_xc[itr][0][self.skip],      \
            self.force_data_xc[itr][1][self.skip])
        self.force_plot_dm.set_UVC(self.force_data_dm[itr][0][self.skip],      \
            self.force_data_dm[itr][1][self.skip])
        #-----------------------------------------------------------------------
        # Update the axis
        #-----------------------------------------------------------------------
        self.ax_force.axis((self.curr_axis[0][itr],self.curr_axis[1][itr],\
                            self.curr_axis[2][itr],self.curr_axis[3][itr]))
        self.ax_force_xc.axis((self.curr_axis[0][itr],self.curr_axis[1][itr],\
                            self.curr_axis[2][itr],self.curr_axis[3][itr]))
        self.ax_force_dm.axis((self.curr_axis[0][itr],self.curr_axis[1][itr],\
                            self.curr_axis[2][itr],self.curr_axis[3][itr]))
        #-----------------------------------------------------------------------
        # Display the progress bar to time how long the rendering takes
        #-----------------------------------------------------------------------
        progress_bar(itr,self.max_iter,'Generating forces video ')
        return self.force_plot,

    def save_force_plot(self,mag_data,diff_times,axis_data,ani_name,comp_plot, \
        force_data,force_data_xc=[],force_data_dm=[]):
        """Function to produce the actual animation of the magnetization and the forces. 
        Args:
        ----------
            - mag_data: (float [times,N,N] array) interpolated data containing the z-direction of the magnetization for each time step
            - diff_times: (int) number of frames in the video
            - axis_data: (float [times,4] array) array containing the limits for the visualization of the data at each time step. It makes sure that the skyrmion is centered in the image.
            - ani_name: (str) name of the output video.
            - comp_plot: (boolean) logical indicating whether the decomposed forces will be animated or not.
            - force_data: (float [times,2,N,N] array) x-y components of the total magnetic force acting over the skyrmion.
            - force_data_xc: (float [times,2,N,N] array, optional) x-y components of the exchange magnetic force acting over the skyrmion.
            - force_data_dm: (float [times,2,N,N] array, optional) x-y components of the DMI magnetic force acting over the skyrmion.
        Author:
        ----------
        Jonathan Chico
        """
        import matplotlib.animation as anim
        import getpass
        #-----------------------------------------------------------------------
        # Save data to class instances so that it can be used by the update 
        # function
        #-----------------------------------------------------------------------
        self.mag_data=mag_data
        if comp_plot:
            self.force_data=force_data
            self.force_data_xc=force_data_xc
            self.force_data_dm=force_data_dm
        else:
            self.force_data=force_data
        self.curr_axis=axis_data
        self.max_iter=diff_times
        #-----------------------------------------------------------------------
        # Define the options of the writer
        #-----------------------------------------------------------------------
        Writer = anim.writers['ffmpeg']
        metadata=dict(title=ani_name, artist=getpass.getuser(),                \
            comment='Animation produced via APyInSky')
        writer = Writer(fps=15, metadata=metadata, bitrate=1800)
        #-----------------------------------------------------------------------
        # Call the animation creator
        #-----------------------------------------------------------------------
        if comp_plot:
            ani = anim.FuncAnimation(self.fig_force,                           \
                self.update_comp_force_plot,blit=True,frames=diff_times)
        else:
            ani = anim.FuncAnimation(self.fig_force,self.update_force_plot,    \
                blit=True,frames=diff_times)
        #-----------------------------------------------------------------------
        # Save the animation to file
        #-----------------------------------------------------------------------
        ani.save(ani_name, writer=writer)
        #-----------------------------------------------------------------------
        # Cleanup the unnecessary data
        #-----------------------------------------------------------------------
        del self.force_data, self.mag_data
        del self.ax_force,self.fig_force,ani
        if comp_plot:
            del self.force_data_xc,self.force_data_dm
        return
