#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 15:40:01 2022

@author: Luigi Sante Zampa
"""

# -----------------------------------------------------------------------------
# Import libraries
from nllgrid import NLLGrid
from . import utils as utl
from . import plot as lsz_plot
from . import raster_tools as rt

# -----------------------------------------------------------------------------
# Set the aliases for some libraries from the utils module
sys = utl.sys
np = utl.np
os = utl.os
shutil = utl.shutil
platform = utl.platform

# -----------------------------------------------------------------------------
# Set the path to the current directory
mdir = os.path.dirname(os.path.abspath(__file__))
if mdir not in sys.path:
    sys.path.append(mdir)

# -----------------------------------------------------------------------------
# Set the path separator
s = os.sep

# -----------------------------------------------------------------------------
# Global raster map of Earth (RGB)

raster_globe = mdir +s+ 'raster_globe.tiff'

# -----------------------------------------------------------------------------
# Global topography

raster_topo = mdir +s+ 'SRTMvf_250m_WGS84_CGIAR.tif'

# -----------------------------------------------------------------------------
# Deafault parameter string

# -----------------------------------------------------------------------------
#%% Parameters Headers

param_header = """
# =============================================================================
#  NonLinLoc programs control file
#
#  Anthony Lomax <anthony@alomax.net>
#
#  See "Control File" and "Running the Sample Location" pages
#     in the NonLicLoc on-line documentation:
#     http://www.alomax.net/nlloc
# =============================================================================
# 
# = comment
# 
# non-nested include files allowed, use:
# INCLUDE - Include
# optional, repeatable
# Syntax 1: INCLUDE includeFile
# Inserts text from another file at current positon in control file.
#
#    includeFile (string) path and name of file to include
# 
# =============================================================================
# =============================================================================
# Generic control file statements
# =============================================================================
#
# CONTROL - Control
# required, non-repeatable
# Syntax 1: CONTROL messageFlag randomNumberSeed
# Sets various general program control parameters.
#
#    messageFlag (integer, min:-1, default:1) sets the verbosity level for messages printed to the terminal ( -1 = completely silent, 0 = error messages only, 1 = 0 + higher-level warning and progress messages, 2 and higher = 1 + lower-level warning and progress messages + information messages, ...)
#    randomNumberSeed (integer) integer seed value for generating random number sequences (used by program NLLoc to generate Metropolis samples and by program Time2EQ to generate noisy time picks)
#
# CONTROL 1 100 
#  
# TRANS - Geographic Transformation
# required, non-repeatable
# Syntax 1: TRANS GLOBAL
# Syntax 2: TRANS SIMPLE latOrig longOrig rotAngle
# Syntax 3: TRANS NONE
# Syntax 4: TRANS SDC latOrig longOrig rotAngle
# Syntax 5: TRANS LAMBERT refEllipsoid latOrig longOrig firstStdParal secondStdParal rotAngle
# Syntax 6: TRANS TRANS_MERC refEllipsoid latOrig longOrig rotAngle
# Syntax 7: TRANS AZIMUTHAL_EQUIDIST refEllipsoid latOrig longOrig rotAngle
# Sets geographic to working coordinates transformation parameters. The GLOBAL option sets spherical regional/teleseismic mode, with no geographic transformation - most position values are input and used directly as latitude and longitude in degrees. The SIMPLE, SDC, LAMBERT and TRANS_MERC options make transformations of geographic coordinates into a Cartesian/rectangular system. The NONE transformation performs no geographic conversion.
#
#    latOrig (float, min:-90.0, max:90.0) latitude in decimal degrees of the rectangular coordinates origin
#    longOrig (float, min:-180.0, max:180.0) longitude in decimal degrees of the rectangular coordinates origin
#    rotAngle (float, min:-360.0, max:360.0) rotation angle of geographic north in degrees clockwise relative to the rectangular coordinates system Y-axis
#    refEllipsoid (choice: WGS-84 GRS-80 WGS-72 Australian Krasovsky International Hayford-1909 Clarke-1880 Clarke-1866 Airy Bessel Hayford-1830 Sphere) reference ellipsoid name
#    latOrig (float, min:-90.0, max:90.0) latitude in decimal degrees of the rectangular coordinates origin
#    longOrig (float, min:-180.0, max:180.0) longitude in decimal degrees of the rectangular coordinates origin
#    firstStdParal secondStdParal (float, min:-90.0, max:90.0) first and second standard parallels (meridians) in decimal degrees
#    rotAngle (float, min:-360.0, max:360.0) rotation angle of geographic north in degrees clockwise relative to the rectangular coordinates system Y-axis
#
# TRANS  NONE                                             
# 
# MAPLINE - Geographic Maplines
# optional, repeatable
# Syntax 1: MAPLINE formatType name red green blue lineStyle
# Specifies a file and drawing parameters for geographic line data.
#
#    formatType (choice: GMT_LATLON GMT_LONLAT XY_LONLAT GMT_LONLATDEPTH GMT_LONLATELEV_M GMT_GRD) line file format or GMT grd file format
#    name (string) full path and file name
#    red green blue (float, min:0.0, max:1.0) red, green and blue intensities (0.0-1.0) (not implemented)
#    lineStyle (choice: SOLID DASHED DOTTED DASHDOT) line style (not implemented)
#
# MAPLINE  GMT_LONLAT ./data_geog/alaska_coasts.xy  0.0 0.0 0.0  SOLID
# MAPLINE  GMT_LONLAT ./data_geog/alaska_rivers.xy  0.0 0.0 1.0  SOLID
# MAPLINE  GMT_LONLAT ./data_geog/alaska.xy  0.0 0.0 0.0  SOLID
# 
# =============================================================================
# END of Generic control file statements
# =============================================================================
# =============================================================================
"""
# -----------------------------------------------------------------------------
#%% Parameters Vel2Grid

param_vel2grid = """
# =============================================================================
# =============================================================================
# Vel2Grid control file statements
# =============================================================================
#
# VGOUT - Output File Root Name
# required, non-repeatable
# Syntax 1: VGOUT fileRoot
# Specifies the directory path and file root name (no extension) for the output velocity grid.
# fileRoot (string) full or relative path and file root name (no extension) for output
#
# VGOUT  ./model/tomo 
#
# VGTYPE - Wave Type
# required, repeatable
# Syntax 1: VGTYPE waveType
# Specifies the physical wave type for a velocity grid.
#
# WaveType (choice: P S) wave type
#
# VGTYPE P
#
# VGGRID - Grid Description
# required, non-repeatable
# Syntax 1: VGGRID xNum yNum zNum xOrig yOrig zOrig dx dy dz gridType
# Specifies the size and type of the 3D velocity grid.
#
# xNum yNum zNum (integer, min:2) number of grid nodes in the x, y and z directions
# xOrig yOrig zOrig (float) x, y and z location of the grid origin in km relative to the geographic origin.
# dx dy dz (float) grid node spacing in kilometers along the x, y and z axes
# gridType (choice: VELOCITY VELOCITY_METERS SLOWNESS VEL2 SLOW2 SLOW_2_METERS SLOW_LEN) 
# physical quantity to store on grid ( VELOCITY = km/s, 
# VELOCITY_METERS = m/s, SLOWNESS = s/km, 
# VEL2 = vel**2, SLOW2 = (s/km)**2, 
# SLOW_2_METERS = slow**2 ((s/m)**2), SLOW_LEN = slow*dx (sec)).
#
# Layer 2DGrid (NOTE: num_grid_x must be = 2 for 2D grids)
# VGGRID   175  114   52  -95.000000  -79.000000   -4.000000   1.00000   1.00000   1.00000 SLOW_LEN
#
# -----------------------------------------------------------------------------
# velocity model description
# -----------------------------------------------------------------------------
#
# VGINP path/model_file_in.mod SIMUL2K 
#
# LAYER - Velocity Model - Layer
# optional, repeatable
# Syntax 1: LAYER depth VpTop VpGrad VsTop VsGrad rhoTop rhoGrad
# Specifies a constant or gradient velocity layer.
#
#    depth (float) depth to top of layer (use negative values for layers above z=0)
#    VpTop VsTop rhoTop (float) P velocity, and S velocity in km/s and density in kg/m**3 at the top of the layer.
#    VpGrad VsGrad rhoGrad (float) Linear P velocity and S velocity gradients in km/s/km and density gradient in kg/m**3/km increasing directly downwards from the top of the layer.
#
# scak model used by the Alaska Earthquake Center (Silwal and Tape, 2016)
# https://agupubs.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2F2015JB012588&file=jgrb51546-sup-0001-supplementary.pdf
# LAYER depth Vp_top Vp_grad Vs_top Vs_grad den_top den_grad
# LAYER   -4.0     4.80    0    2.70    0    2.40    0
# LAYER    0.0     5.85    0    3.29    0    2.40    0
# LAYER    22.0    6.80    0    3.82    0    2.50    0
# LAYER    39.5    8.00    0    4.49    0    2.80    0
# =============================================================================
# END of Vel2Grid control file statements
# =============================================================================
# =============================================================================
"""
# -----------------------------------------------------------------------------
#%% Parameters Grid2Time

param_grd2time = """
# =============================================================================
# =============================================================================
# Grid2Time control file statements
# =============================================================================
#
# GTFILES - Input and Output File Root Name
# required, non-repeatable
# Syntax 1: GTFILES ttimeFileRoot outputFileRoot waveType iSwapBytesOnInput
# Specifies the directory path and file root name (no extension), and the wave type identifier for the input velocity grid and output time grids.
#
#    ttimeFileRoot (string) full or relative path and file root name (no extension) for input velocity grid (generated by program Vel2Grid)
#    outputFileRoot (string) full or relative path and file root name (no extension) for output travel-time and take-off angle grids
#    waveType (choice: P S) wave type
#    iSwapBytesOnInput (integer, min:0, max:1, default:0) flag to indicate if hi and low bytes of input velocity grid file should be swapped
#
# GTFILES  ./model/layer  ./time/layer P	  # uncomment to generate P travel times
#
# GTMODE - Program Modes
# required, non-repeatable
# Syntax 1: GTMODE gridMode angleMode
# Specifies several program run modes.
#
#    gridMode (choice: GRID3D GRID2D) grid type ( GRID3D for a 3D, Nx*Ny*Nz grid or GRID2D for a 2D, 2*Ny*Nz grid)
#    angleMode (choice: ANGLES_YES ANGLES_NO) sets if take-off angles are calculated and an angles grid is output ( ANGLES_YES for angles calulcation or ANGLES_NO for no angles calculation)
#
# GTMODE GRID3D ANGLES_YES
#
# -----------------------------------------------------------------------------
# description of source (e.g. seismic station) for calculating travel-time field
# -----------------------------------------------------------------------------
#
# GTSRCE - Source Description
# required, repeatable
# Syntax 1: GTSRCE label XYZ xSrce ySrce zSrce elev
# Syntax 2: GTSRCE label LATLON latSrce longSrce zSrce elev
# Syntax 3: GTSRCE label LATLONDM latDegSrce latMinSrce latDir longDegSrce longMinSrce longDir zSrce elev
# Syntax 4: GTSRCE label LATLONDS latDegSrce latMinSrce latSecSrce latDir longDegSrce longMinSrce longSecSrce longDir zSrce elev
# Specifies a source location. One time grid and one angles grid (if requested) will be generated for each source. Four formats are supported: XYZ (rectangular grid coordinates), LATLON (decimal degrees for latitude/longitude), LATLONDM (degrees + decimal minutes for latitude/longitude) and LATLONDS (degrees + minutes + decimal seconds for latitude/longitude).
#
#    label (string) source label ( i.e. a station code: ABC )
#    xSrce ySrce (float) x and y grid positions relative to geographic origin in kilometers for source
#    zSrce (float) z grid position (depth, positive DOWN) in kilometers for source
#    elev (float) elevation above z grid position (positive UP) in kilometers for source
#    latSrce (float, min:-90.0, max:90.0) latitude in decimal degrees for source (pos = North)
#    longSrce (float, min:-180.0, max:180.0) longitude in decimal degrees for source (pos = East)
#    latDegSrce latMinSrce latSecSrce (float) latitude degrees, minutes and seconds for source
#    longDegSrce longMinSrce longSecSrce (float) longitude degrees, minutes and seconds for source
#    latDir (choice: N S) geographic direction
#    longDir (choice: W E) geographic direction
#
# Examples:
#
# GTSRCE  STA   XYZ   27.25  -67.78  0.0  1.242
# GTSRCE  CALF  LATLON   43.753  6.922  0.0  1.242
# GTSRCE  JOU  LATLONDM  43 38.00 N  05 39.52 E   0.0   0.300
#
# Include file with STSRCE statements, this line can be placed anywhere in control file
# INCLUDE stations_0.stz 
# 
# GT_PLFD - Podvin and Lecomte Finite Difference
# required, non-repeatable, for Podvin and Lecomte finite difference, must not be present otherwise
# Syntax 1: GT_PLFD hs_eps_init message_flag
# Selects Podvin and Lecomte finite difference method and specifies method parameters.
#
#    hs_eps_init (float, min:0.0) fraction (typically 1.0E-3) defining the tolerated model inhomogeneity for exact initialization. A tolerance larger than 0.01 will potentially create errors larger than those involved by the F.D. scheme without any exact initialization.
#    message_flag (integer, min:0, max:2) Message flag (0:silent, 1:few messages, 2:verbose) A negative value inhibits "clever" initialization.
#
# GT_PLFD  1.0e-3  0
#
# =============================================================================
# END of Grid2Time control file statements
# =============================================================================
# =============================================================================
"""

# -----------------------------------------------------------------------------
#%% Parameters Loc

param_loc = """
# =============================================================================
# =============================================================================
# NLLoc control file statements
# =============================================================================
#
# LOCSIG NonLinLoc - ALomax Scientific
# LOCCOM Rilocalizzazione sintetici 
# LOCFILES obs/eq001.obs NLLOC_OBS           ./time/layer  ./loc/loc
# LOCHYPOUT SAVE_NLLOC_ALL  SAVE_HYPOINV_SUM 
# LOCSEARCH  OCT 10 10 4 0.01 20000 5000 0 1
# LOCGRID   175  114   48  -95.000000  -79.000000    0.000000   1.00000   1.00000   1.00000 PROB_DENSITY  SAVE
# LOCMETH EDT_OT_WT 9999.0 4 -1 -1  -1.730 -1 0.1
# LOCGAU 0.15 1.0 
# LOCGAU2 0.5 0.05 1.0 
# LOCPHASEID  P   P p G PN PG 
# LOCQUAL2ERR 0.1 0.5 1.0 2.0 99999.9 
# LOCPHSTAT 15 -1 9999.0 1.0 1.0 
# LOCANGLES ANGLES_NO 5 
# LOCMAG ML_HB 1.0 1.110 0.00189 
# LOCELEVCORR 0 4.0P 2.87S
# LOCSTAWT 0 -1
# LOCTOPO_SURFACE 
#
# =============================================================================
# END of NLLoc control file statements
# =============================================================================
# =============================================================================
"""
# -----------------------------------------------------------------------------
#%% Functions 

# -----------------------------------------------------------------------------
def create_vel2grid_file( control=[1, 54321], 
                          trans="NONE", 
                          vgout='layer',
                          path=mdir,
                          name='vel2grid.in',
                          vpvs=1.78,
                          vgtype=['P'],
                          vggrid=[],
                          vginp=None,
                          layer=False,
                          depth=[],
                          vp_top=[],
                          vp_grad=[],
                          vs_top=[],
                          vs_grad=[],
                          den_top=[],
                          den_grad=[],
                          reverse_depth=False,
                          run_file=False,
                          plot=False,
                          vmin=None,
                          vmax=None, 
                          copy_bin=True,
                          section=[ 0, None, None ],
                          printf=True,
                          print_prefix="\t" ) :

    os.makedirs( path, exist_ok=True )
    os.makedirs( path+s+'model', exist_ok=True )

    if ( type( vginp ) == str ) :
        in_model_name = vginp.split( s )[-1]
        if os.path.isfile( path +s+ in_model_name ) :
            if shutil._samefile( path+s+in_model_name, vginp ) is False :
                shutil.copy( vginp, path+s+in_model_name ) 
    if ( type( vginp ) in [list, tuple] ) and ( len(vginp)!=0 ) :
        in_model_name = vginp[0].split( s )[-1]
        if os.path.isfile( path +s+ in_model_name ) :
            if shutil._samefile( path+s+in_model_name, vginp[0] ) is False :
                shutil.copy( vginp[0], path+s+in_model_name )   

    if vginp is None :
        vginp = []
    li = ""

        # Print the parameter file if the printf flag is True
    if printf is True :
        print( f"\n\n{print_prefix}# -----------------------------------------" )
        print( f"{print_prefix}Vel2Grid parameters file :" +\
               f"\n{print_prefix + path +s+ name} " 
               f"\n{print_prefix}CONTENT:" )
    
    with open( path +s+ name, 'w' ) as f :
        
        for l in param_header.splitlines() :
            
            if ("CONTROL" in l) and ("Control" not in l) and ( "Syntax" not in l ):
                l = f"CONTROL {control[0]} {control[1]}"

            if ("TRANS " in l) and ("Syntax" not in l) and ("TRANS - Geographic" not in l):
                if type( trans ) == str :
                    l = f"TRANS {trans}"
                if type( trans ) in ( list, tuple ) :
                    l = " ".join(str(w) for w in trans)
                
            f.write( l + ' \n' )

            if l.startswith( '#' ) == False :
                print( print_prefix + l )            
            
        for l in param_vel2grid.splitlines() :
            
            if ( "VGOUT" in l ) and ( "- Output File Root Name" not in l ) and\
                ( "Syntax" not in l ):
                l = f"{print_prefix}VGOUT model{s+vgout}"
            if ( "VGTYPE" in l ) and ( "- Wave Type" not in l ) and ( "Syntax" not in l ) :
                l = ""
                for vgt in vgtype :
                    l = l + f"{print_prefix}VGTYPE {vgt} \n"
            if  ( vggrid != [] ) and ( "VGGRID" in l ) and ( "- Grid Description" not in l ) and\
                ( "Syntax" not in l ):
                if len( vggrid ) < 10 :
                    vggrid = vggrid + "SLOW_LEN"
                l = f"{print_prefix}VGGRID"
                for vgg in vggrid :
                    l = l + f" {vgg}"
            if  ( "VGINP" in l ) and ( vginp != [] ) and ( layer is False ):
                if type( vginp ) == str :
                    vginp = [ vginp ] + [ 'SIMUL2K' ] 
                l = f"{print_prefix}VGINP {in_model_name} {vginp[1]}"   
            
            if  ( layer is True ) and ( "LAYER" in l ) and \
                ( '- Velocity Model -' not in l ) and ( 'Syntax' not in l ) and \
                ( 'depth Vp_top' not in l ) and ( li == ""  ) :
                    if reverse_depth == True :
                        depth = [ i*-1 for i in depth ]
                    if vp_grad == [] :
                        vp_grad = [ 0 for i in vp_top ]
                    if vs_top == [] :
                        vs_top = [ i / vpvs for i in vp_top ]
                    if ( vs_grad == [] ) and ( vp_grad == [] ) :
                        vs_grad = [ 0 for i in vp_top ]
                    if ( vs_grad == [] ) and ( vp_grad != [] ) :
                        vs_grad = [ i / vpvs for i in vp_grad ]
                    if den_top == [] :
                        den_top = [ 0.31 * ((i*1e3)**0.25) for i in vp_top ]
                    if den_grad == [] :
                        den_grad = [ 0 for i in vp_top ]
                    for i, d in enumerate( depth ) :
                        li += f"{print_prefix}LAYER {depth[i]: >8.3f} {vp_top[i]: >8.4f} {vp_grad[i]: >8.4f} {vs_top[i]: >8.3f}"
                        li += f"{vs_grad[i]: >8.4f} {den_top[i]: >8.2f} {den_grad[i]: >8.2f} \n"
                    l = li
                    
            elif ( vginp == [] ) and ( layer == True ) and ( "LAYER" in l ) and \
                 ( '- Velocity Model -' not in l ) and ( 'Syntax' not in l ) and \
                 ( 'depth Vp_top' not in l ) and ( li != ""  ) :
                     continue
                
            f.write( l + ' \n' )

            if printf is True :
                if l.startswith( '#' ) == False :
                    print(  l )

    f.close()

    if run_file is True :

        if layer is False :
            run_nlloc( path=path,
            	       vel2grid=True,	 
                       vel2grid_file= path +s+ name,
                       gtmode="GRID3D ANGLES_YES", 
                       copy_bin=copy_bin )
        
        if layer is True :
            run_nlloc( path=path,
            	       vel2grid=True,	 
                       vel2grid_file=path +s+ name,
                       gtmode="GRID2D ANGLES_YES", 
                       copy_bin=copy_bin )

    if plot is True :
        _ = read_model_file( path +s+ 'model', plot=True, section=section,
                        vmin=vmin, vmax=vmax )

    return path +s+ 'model'

# -----------------------------------------------------------------------------
def create_model_file( model=None, 
                       step=None, 
                       float_type="FLOAT",
                       lim=None, 
                       plot=False, 
                       name='layer', 
                       method='nearest',
                       path=os.getcwd(), 
                       proj_name="SIMPLE", 
                       orig_lon=0, 
                       orig_lat=0, 
                       **kwargs ) :
    """
    Creates a model file for NLLocPy.

    Args:
        model (dict): A dictionary containing the model information. It should have the following keys:
            - 'x': Array of x-coordinates.
            - 'y': Array of y-coordinates.
            - 'z': Array of z-coordinates.
            - 'vp': Array of P-wave velocities.
            - 'vs' (optional): Array of S-wave velocities.
        step (float): The grid spacing for the model.
        float_type (str): The data type for the model grid. Default is "FLOAT".
        lim (tuple): The limits of the model grid in the form (xmin, xmax, ymin, ymax, zmin, zmax).
        plot (bool): Whether to plot the model after creating the file. Default is False.
        name (str): The name of the model file. Default is 'layer'.
        method (str): The interpolation method for gridding the model. Default is 'nearest'.
        path (str): The path where the model file will be saved. Default is the current working directory.
        proj_name (str): The name of the project. Default is "SIMPLE".
        orig_lon (float): The longitude of the origin. Default is 0.
        orig_lat (float): The latitude of the origin. Default is 0.
        section (list): The section to plot. It should be a list in the form [start, end, step]. Default is [0, None, None].

    Returns:
        tuple: A tuple containing the paths to the P-wave and S-wave model files, and the NLLGrid object.

    """
    os.makedirs( path, exist_ok=True )
    os.makedirs( path+s+'model', exist_ok=True )
    
    grd = NLLGrid()
    
    if model is not None :

        x = model['x']
        y = model['y']
        z = model['z']
        vp = model['vp']
        if 'vs' in model :
            vs = model['vs']
        else :
            vs = None
        
    if lim is not None :
        xmin, xmax, ymin, ymax, zmin, zmax = lim
    else :
        xmin, xmax = x.min(), x.max()
        ymin, ymax = y.min(), y.max()
        zmin, zmax = z.min(), z.max()
    
    if step == None :
        xstep = np.min( np.diff( np.unique( x ) ) )
        ystep = np.min( np.diff( np.unique( y ) ) )
        zstep = np.min( np.diff( np.unique( z ) ) )
        step = np.min( [ xstep, ystep, zstep ] )
    
    xu = np.arange( xmin-step, xmax+step+step/2, step )
    yu = np.arange( ymin-step, ymax+step+step/2, step )
    zu = np.arange( zmax+step, zmin-step+step/2, -step )
    
    Y, X, Z = np.meshgrid( yu, xu, zu )
    
    Vp = utl.sp.interpolate.griddata( (x,y,z), vp, (X,Y,Z), method=method )
    grd.array = ( 1/Vp[ :, :, : ] ) * step 
    # grd.array = Vp
    grd.dy = step #km
    grd.dx = step #km
    grd.dz = step #km
    grd.x_orig = xu.min()  #km
    grd.y_orig = yu.min()  #km
    grd.z_orig = -zu.max() #km
    grd.nx = len( xu )
    grd.ny = len( yu )
    grd.nz = len( zu )
    grd.type = 'SLOW_LEN'
    grd.float_type = float_type
    grd.orig_lat = orig_lat
    grd.orig_lon = orig_lon
    grd.proj_name = proj_name

    file_P_grid = path +s+'model' +s+ name + '.P.mod.grd'
    grd.basename = path +s+'model' +s+ name + '.P.mod'
    grd.write_hdr_file()
    grd.write_buf_file()

    if vs is not None :
        Vs = utl.sp.interpolate.griddata( (x,y,z), vs, (X,Y,Z), method=method )
        grd.array = ( 1/Vs[ :, :, : ] ) * step
        file_S_grid = path +s+'model' +s+ name + '.S.mod.grd'
        grd.basename = path +s+'model' +s+ name + '.S.mod'
        grd.write_hdr_file()
        grd.write_buf_file()
    else :
        file_S_grid = None 

    if plot == True :

        _ = read_model_file( path+s+'model', plot=True, **kwargs )

    return file_P_grid, file_S_grid, grd

# -----------------------------------------------------------------------------
def create_station_file( st_dict=None, 
                         st_name=None, 
                         x=None, 
                         y=None, 
                         lon=None,
                         lat=None,
                         z=0, 
                         elev=None, 
                         label="XYZ", 
                         path=os.getcwd(), 
                         file_name="stations.stz", 
                         elev_m2km=True, 
                         printf=False, 
                         **prj_kwargs ) :
    """
    Create a Nlloc stations file.

    Args:
        - st_dict (dict): A dictionary containing station information.

        - st_name (list): A list of station names.

        - x (array-like): An array-like object containing x-coordinates of the stations.

        - y (array-like): An array-like object containing y-coordinates of the stations.

        - z (float or int, optional): The z-coordinate of the stations. Defaults to 0.

        - elev (array-like, optional): An array-like object containing elevation values 
            of the stations. Defaults to None.

        - label (str, optional): The label for the coordinate system. Defaults to "XYZ".

        - path (str, optional): The path where the stations file will be created. 
            Defaults to the current working directory.
            
        - file_name (str, optional): The name of the stations file. 
            Defaults to "stations.stz".

        - elev_m2km (bool, optional): Whether to convert elevation values 
            from meters to kilometers. Defaults to True.

        - printf (bool, optional): Whether to print the contents of the stations file. 
            Defaults to False.

        - prj_kwargs (dict): Keyword arguments for the projection. 
            See the function cart2geo or geo2cart for more information 
            about the dictionary arguments. 

    Returns:
        str: The path of the created stations file.
    """

    os.makedirs( path, exist_ok=True )
    
    if st_dict is not None :

        for k in st_dict :
            if type( st_dict[ k ] ) in ( list, tuple ) :
                st_dict[ k ] = np.array( st_dict[ k ] )
            
        st_name = st_dict['st']
        
        if elev is None :
            elev = st_dict['z']
        if 'x' in st_dict :
            x = st_dict['x']
        else :
            x = None
        if 'y' in st_dict :
            y = st_dict['y']
        else :
            y = None
        if 'lon' in st_dict :
            lon = st_dict['lon']
        else :
            lon = None
        if 'lat' in st_dict :
            lat = st_dict['lat']
        else :
            lat = None    

    if type( label ) == str :
        label = np.full( len( st_name ), label )

    if lon is None :
        lon, lat = cart2geo( x, y, **prj_kwargs )

    if x is None :
        x, y = geo2cart( lon, lat, **prj_kwargs )

    if type( z ) in ( int, float ) :
        z = np.full( len( x ), z )
        
    if z is None :
        z = np.zeros( len( x ) ) 
        
    if type( elev ) in ( int, float ) :
        elev = np.full( len( x ), elev ) 
        
    if elev is None :
        elev = np.zeros( len( x ) ) 
        
    if elev_m2km is True :
        elev = elev / 1e3
    
    len_nm = []
    for nm in st_name :
        len_nm.append( len( nm ) )
    max_nm_len = np.max( len_nm )
    
    if max_nm_len < 4 :
        max_nm_len = 4
        
    with open( path +s+ file_name, 'w' ) as f :
        
        for i, nm in enumerate( st_name ) :
            
            if label[i] in [ "XYZ", "LONLAT" ] :
                l = f"GTSRCE {nm} {label[i]} {lon[i]:.6f} {lat[i]:.6f} " + \
                    f"{z[i]:.3f} {elev[i]:.3f}"
            if label[i] in [ "LATLON" ] :
                l = f"GTSRCE {nm} {label[i]} {lat[i]:.6f} {lon[i]:.6f} " + \
                    f"{z[i]:.3f} {elev[i]:.3f}"
            if label[i] in [ "XYZ" ] :
                l = f"GTSRCE {nm} {label[i]} {x[i]:.3f} {y[i]:.3f} " + \
                    f"{z[i]:.3f} {elev[i]:.3f}"
                
            f.write( l + ' \n')
        
        f.close()

    # Print the contents of the stations file if the printf flag is True
    if printf == True :
        print( f"\tNlloc stations file:\n{path +s+ file_name}")
        file_lines_list = utl.read_file( path +s+ file_name )
        print( '\t',*file_lines_list[:5], sep="\n" )
        if len( file_lines_list ) > 5 :
            print( f"...\n\tTotal lines: {len(file_lines_list)}")

    return path +s+ file_name

# -----------------------------------------------------------------------------
def create_grd2time_file( path, 
                          stFile, 
                          control=[1, 54321], 
                          trans=None, 
                          gtmode="GRID3D ANGLES_YES", 
                          ttimeFileRoot='model'+s+'layer', 
                          outputFileRoot='time'+s+'layer',
                          gt_plfd=[1.0e-3, 0], 
                          name="grd2time.in", 
                          wavetype="P",
                          all_files_in=True,
                          printf=True,
                          start_clean=True ) :
    
    """
    Creates a grd2time parameter file.
    For more information on the parameters, see the NonLinLoc documentation at:
    http://alomax.free.fr/nlloc/soft7.00/control.html#_Grid2Time_
    ( link to the "Grid2Time control file statements" section of the webpage http://alomax.free.fr/nlloc/ )

    Args:
        - path (str): The path where the parameter file will be created.
        
        - stFile (str or dict): The path to the station file or a dictionary 
            containing station information. Defaults to None.
        
        - control (list): A list containing two integers representing the Nlloc 
            control parameters. Defaults to [1, 54321].
        
        - trans (str): The transformation type. Defaults to "NONE".
        
        - gtmode (str): The GTMODE parameter value. Specifies the grid type and 
            whether take-off angles are calculated. Options are "GRID3D ANGLES_YES" 
            for a 3D grid with angles calculation, and "GRID2D ANGLES_NO" for a 2D 
            grid without angles calculation. Defaults to "GRID3D ANGLES_YES".
        
        - ttimeFileRoot (str): The root name of the ttime files. Specifies the 
            directory path and file root name (no extension), and the wave type 
            identifier for the input velocity grid. Defaults to 'model'+s+'layer'.
        
        - outputFileRoot (str): The root name of the output files. Specifies the 
            directory path and file root name (no extension) for output travel-time 
            and take-off angle grids. Defaults to 'time'+s+'layer'.
        
        - gt_plfd (list): A list containing two floats representing the GT_PLFD 
            parameters. Selects Podvin and Lecomte finite difference method and 
            specifies method parameters. The first float is the fraction defining 
            the tolerated model inhomogeneity for exact initialization. The second 
            integer is the message flag (0:silent, 1:few messages, 2:verbose). 
            Defaults to [1.0e-3, 0].
        
        - name (str): The name of the input parameter file 
            that is created with this function. Defaults to "grd2time.in".
        
        - wavetype (str): The wave type. Options are "P" for primary waves and 
            "S" for secondary waves. Defaults to "P".
        
        - all_files_in (bool): A flag indicating whether to copy all input files 
            to the specified path. Defaults to True.

        - printf (bool): A flag indicating whether to print the parameter file

        - start_clean (bool): A flag indicating whether to remove the directory

    Returns:
        str: The path to the created parameter file.
    """

    # Check if stFile is a string (path to the station file)
    if type( stFile ) == str :
        # Check if the path is absolute
        if os.path.isabs( stFile ) :
            # Check if the file exists
            if os.path.isfile( stFile ) is False :
                # Raise an error if the file does not exist
                raise ValueError( f"File {stFile} does not exist." )
            else :
                # Get the relative path of the file
                stFile_rel = os.path.relpath( stFile, start=path )
        else :
            # Check if the file exists in the given path
            if os.path.isfile( path +s+ stFile ) is False :
                # Raise an error if the file does not exist
                raise ValueError( f"File {path +s+ stFile} does not exist." )
            else :
                # If the file exists, set the relative path to the file name
                stFile_rel = stFile

    # Check if stFile is a dictionary (station information)
    if type( stFile ) == dict :
        # Set the station file name
        stfile_name = "stations.stz"
        # Create the station file
        stFile = create_station_file( st_dict=stFile, 
                                      path=path, 
                                      file_name=stfile_name )
        stFile_rel = os.path.relpath( stFile, start=path )

    # Get the directory and name of the ttime file
    ttimeFileDir = os.path.dirname( ttimeFileRoot )
    ttimeDirName = os.path.basename( ttimeFileDir )
    # Check if the ttime file directory is an absolute path
    if os.path.isabs( ttimeFileDir ) :
        # Check if the directory exists
        if os.path.isdir( ttimeFileDir ) :
            # Get the relative path of the ttime file
            ttimeFileRoot = os.path.relpath( ttimeFileDir, start=path ) +s+ ttimeFileRoot.split( s )[-1]
        else :
            # Raise an error if the directory does not exist
            raise ValueError( f"Directory {ttimeFileDir} does not exist." )
    else :
        # Check if the directory exists in the given path
        if os.path.isdir( path +s+ ttimeFileDir ) is False:
            # Raise an error if the directory does not exist
            raise ValueError( f"Directory {path +s+ ttimeFileDir} does not exist." )

    # Get the directory and name of the output file
    outputFileDir = os.path.dirname( outputFileRoot )
    outputDirName = os.path.basename( outputFileDir )
    # Check if the output file directory is an absolute path
    if os.path.isabs( outputFileDir ) :
        # Check if the directory exists
        if os.path.isdir( outputFileDir ) :
            # Get the relative path of the output file
            outputFileRoot = os.path.relpath( outputFileDir, start=path ) +s+ outputFileRoot.split( s )[-1]
        else :
            # Create the directory if it doesn't exist
            os.makedirs( path+s+outputDirName, exist_ok=True )
            # Set the output file root to the directory name and file name
            outputFileRoot = outputDirName +s+ outputFileRoot.split( s )[-1]
    else :
        # Check if the directory exists in the given path
        if os.path.isdir( path +s+ outputFileDir ) is False:
            # Create the directory if it doesn't exist
            os.makedirs( path+s+outputDirName, exist_ok=True )
            # Set the output file root to the directory name and file name
            outputFileRoot = outputDirName +s+ outputFileRoot.split( s )[-1]

    # Check if all files should be copied to the specified path
    if all_files_in is True :

        # Copy the station file to the specified path
        stfile_name = stFile.split( s )[-1]
        if utl.are_equal( stFile, path +s+ stfile_name ) is False :
            if os.path.isfile( path +s+ stfile_name  ) :
                os.remove( path +s+ stfile_name )
            shutil.copy( stFile, path +s+ stfile_name )
        stFile_rel = stfile_name

        # Copy the ttime file directory to the specified path
        if utl.are_equal( ttimeFileDir, path +s+ ttimeDirName ) is False :
            if os.path.isdir( path +s+ ttimeDirName ) :
                shutil.rmtree( path +s+ ttimeDirName )
            shutil.copytree( ttimeFileDir, path +s+ ttimeDirName )
        ttimeFileRoot =  ttimeDirName +s+ ttimeFileRoot.split( s )[-1]

        # Copy the output file directory to the specified path
        if utl.are_equal( outputFileDir, path +s+ outputDirName ) is False :
            if os.path.isdir( path +s+ outputDirName ) :
                shutil.rmtree( path +s+ outputDirName )
            shutil.copytree( outputFileDir, path +s+ outputDirName )
        outputFileRoot =  outputDirName +s+ outputFileRoot.split( s )[-1]

    # Check if the start_clean flag is True
    if start_clean is True :
        # Remove the content of the directory of output files
        utl.empty_directory( os.path.dirname( outputFileRoot ) )

    if not trans :
        print( ttimeFileDir)
        model = read_model_file( ttimeFileDir, plot=False )[0]
        trans = f"{model.proj_name} {model.orig_lat} {model.orig_lon} {model.map_rot}"

    # Print the parameter file if the printf flag is True
    if printf is True :
        print( "\n\n# -----------------------------------------" )
        print( f"Grid2time parameters file (V{wavetype}):\n{path +s+ name} \nCONTENT:" )

    # Open the parameter file for writing
    with open( path +s+ name, 'w' ) as f :

        # Write the parameter header to the file
        for l in param_header.splitlines() :
            
            if ("CONTROL" in l) and ("Control" not in l) and ( 'Syntax' not in l ) :
                l = f"CONTROL {control[0]} {control[1]}"

            if ("TRANS " in l) and ("Syntax" not in l) and ("TRANS - Geographic" not in l):
                
                if type( trans ) == str :
                    l = f"TRANS {trans}"
                if type( trans ) in ( list, tuple ) :
                    l = " ".join(str(w) for w in trans)
                
            f.write( l + ' \n' )
            
            if printf is True :
                if l.startswith( '#' ) == False :
                    print( l )

        # Write the grd2time parameters to the file
        for l in param_grd2time.splitlines() :
            
            if ("GTFILES" in l) and ("Input and Output File Root Name" not in l) and\
                ("Syntax" not in l) :
                l = f"GTFILES {ttimeFileRoot} {outputFileRoot} {wavetype}"
            
            if ("GTMODE" in l) and ("Program Modes" not in l) and ("Syntax" not in l) :
                l = f"GTMODE {gtmode}"    
                
            if ("GT_PLFD" in l) and ("Podvin" not in l) and ("Syntax" not in l) :
                l = f"GT_PLFD {gt_plfd[0]} {gt_plfd[1]}"                 

            if ("INCLUDE" in l) :
                l = f"INCLUDE {stFile_rel}"   
                
            f.write( l + ' \n' )

            if printf is True :
                if l.startswith( '#' ) == False :
                    print( l )

    f.close()
    
    return path +s+ name

# -----------------------------------------------------------------------------
def create_loc_file( path, 
                     obsFile, 
                     locgrid=[], 
                     obsFileType="NLLOC_OBS", 
                     control=[1, 54321], 
                     trans=None, 
                     name="loc.in", 
                     locsig=" NonLinLoc - ALomax Scientific", 
                     loccom=" ", 
                     lochypout="SAVE_NLLOC_ALL", 
                     locsearch=["OCT", 10, 10, 4, 0.01, 20000, 5000, 0, 1], 
                     ttimeFileRoot="time"+s+"layer", 
                     outputFileRoot="loc"+s+"loc", 
                     locmeth=["EDT_OT_WT", 9999.0, 4, -1, -1, -1.730, -1, 0.1], 
                     locgau=[0.15, 1.0, ''], 
                     locgau2=None, 
                     locqual2err=[0.1, 0.5, 1, 2, 99999.9], 
                     locphstat=[15, -1, 9999.0, 1.0, 1.0], 
                     locangles=["ANGLES_NO", 5], 
                     locmag=["ML_HB", 1.0, 1.11, 0.00189], 
                     locphaseid=["P", "p", "G", "PN", "PG", "S", "s", "G", "SN", "SG"], 
                     locelevcorr=[], 
                     loctopo_surface=[], 
                     flagDumpDecimation=1, 
                     all_files_in=True, 
                     printf=True, 
                     print_prefix="\t",
                     start_clean=True):
    
    """
    Create a location file for NonLinLoc.
    For more information on the parameters, see the NonLinLoc documentation at:
    http://alomax.free.fr/nlloc/soft7.00/control.html#_NLLoc_
    ( link to the "NLLoc control file statements" section of the webpage http://alomax.free.fr/nlloc/ )

    Parameters:
        - path (str): The path where the location file will be created.
        - obsFile (str): The path to the observation file or the name of the observation file if it exists in the given path.
        - locgrid (list, optional): The locgrid list specifying the grid parameters. Default is an empty list.
        - obsFileType (str, optional): The type of observation file. Default is "NLLOC_OBS".
        - control (list, optional): The control parameters for NonLinLoc. Default is [1, 54321].
        - trans (str, optional): The transformation type. Default is "NONE".
        - name (str, optional): The name of the location file. Default is "loc.in".
        - locsig (str, optional): The signature for the location file. Default is " NonLinLoc - ALomax Scientific".
        - loccom (str, optional): The comment for the location file. Default is an empty string.
        - lochypout (str, optional): The type of output for the hypocenter. Default is "SAVE_NLLOC_ALL".
        - locsearch (list, optional): The search parameters for NonLinLoc. Default is ["OCT", 10, 10, 4, 0.01, 20000, 5000, 0, 1].
        - ttimeFileRoot (str, optional): The root name of the travel time files. Default is "time"+s+"layer".
        - outputFileRoot (str, optional): The root name of the output files. Default is "loc"+s+"loc".
        - locmeth (list, optional): The method parameters for NonLinLoc. Default is ["EDT_OT_WT", 9999.0, 4, -1, -1, -1.730, -1, 0.1].
        - locgau (list, optional): The Gaussian parameters for NonLinLoc. Default is [0.15, 1.0, ''].
        - locgau2 (list, optional): The second Gaussian parameters for NonLinLoc. Default is [0.05, 0.5, ''].
        - locqual2err (list, optional): The quality to error conversion parameters for NonLinLoc. Default is [0.1, 0.5, 1, 2, 99999.9].
        - locphstat (list, optional): The travel-time statistics parameters for NonLinLoc. Default is [15, -1, 9999.0, 1.0, 1.0].
        - locangles (list, optional): The angle parameters for NonLinLoc. Default is ["ANGLES_NO", 5].
        - locmag (list, optional): The magnitude parameters for NonLinLoc. Default is ["ML_HB", 1.0, 1.11, 0.00189].
        - locphaseid (list, optional): The phase identification parameters for phases. Default is ["P", "p", "G", "PN", "PG", "S", "s", "G", "SN", "SG"].
        - locelevcorr (list, optional): The elevation correction parameters for NonLinLoc. Default is an empty list.
        - loctopo_surface (list, optional): The topographic surface parameters for NonLinLoc. Default is an empty list.
        - flagDumpDecimation (int, optional): The flag for dump decimation. Default is 1.
        - all_files_in (bool, optional): Flag to indicate if all files should be copied to the specified path. Default is True.
        - printf (bool, optional): Flag to indicate if the lines should be printed. Default is True.
        - start_clean (bool, optional): Flag to indicate if the output directory should be cleaned before creating the location file. Default is True.
        
        Returns:
            str: The path to the created location file.
    """
    os.makedirs( path, exist_ok=True )

    # Check if obsFile is a string (path to the station file)
    if type( obsFile ) == str :
        # Check if the path is absolute
        if os.path.isabs( obsFile ) :
            # Check if the file exists
            if os.path.isfile( obsFile ) is False :
                # Raise an error if the file does not exist
                raise ValueError( f"File {obsFile} does not exist." )
            else :
                # Get the relative path of the file
                obsFile_rel = os.path.relpath( obsFile, start=path )
        else :
            # Check if the file exists in the given path
            if os.path.isfile( path +s+ obsFile ) is False :
                # Raise an error if the file does not exist
                raise ValueError( f"File {path +s+ obsFile} does not exist." )
            else :
                # If the file exists, set the relative path to the file name
                obsFile_rel = obsFile

    # Get the directory and name of the ttime files
    ttimeFileDir = os.path.dirname( ttimeFileRoot )
    ttimeDirName = os.path.basename( ttimeFileDir )
    # Check if the ttime files directory is an absolute path
    if os.path.isabs( ttimeFileDir ) :
        # Check if the directory exists
        if os.path.isdir( ttimeFileDir ) :
            # Get the relative path of the ttime files
            ttimeFileRoot = os.path.relpath( ttimeFileDir, start=path ) +s+ ttimeFileRoot.split( s )[-1]
        else :
            # Raise an error if the directory does not exist
            raise ValueError( f"Directory {ttimeFileDir} does not exist." )
    else :
        # Check if the directory exists in the given path
        if os.path.isdir( path +s+ ttimeFileDir ) is False:
            # Raise an error if the directory does not exist
            raise ValueError( f"Directory {path +s+ ttimeFileDir} does not exist." )

    # Get the directory and name of the output files
    outputFileDir = os.path.dirname( outputFileRoot )
    outputDirName = os.path.basename( outputFileDir )
    # Check if the output files directory is an absolute path
    if os.path.isabs( outputFileDir ) :
        # Check if the directory exists
        if os.path.isdir( outputFileDir ) :
            # Get the relative path of the output files
            outputFileRoot = os.path.relpath( outputFileDir, start=path ) +s+ outputFileRoot.split( s )[-1]
        else :
            # Create the directory if it doesn't exist
            os.makedirs( path+s+outputDirName, exist_ok=True )
            # Set the output file root 
            outputFileRoot = outputDirName +s+ outputFileRoot.split( s )[-1]
    else :
        # Check if the directory exists in the given path
        if os.path.isdir( path +s+ outputFileDir ) is False:
            # Create the directory if it doesn't exist
            os.makedirs( path+s+outputDirName, exist_ok=True )
            # Set the output file root 
            outputFileRoot = outputDirName +s+ outputFileRoot.split( s )[-1]

    # Check if all files should be copied to the specified path
    if all_files_in is True :

        # Copy the observation file to the specified path
        stfile_name = obsFile.split( s )[-1]
        if utl.are_equal( obsFile, path +s+ stfile_name ) is False :
            if os.path.isfile( path +s+ stfile_name  ) :
                os.remove( path +s+ stfile_name )
            shutil.copy( obsFile, path +s+ stfile_name )
        obsFile_rel = stfile_name

        # Copy the ttime file directory to the specified path
        if utl.are_equal( ttimeFileDir, path +s+ ttimeDirName ) is False :
            if os.path.isdir( path +s+ ttimeDirName ) :
                shutil.rmtree( path +s+ ttimeDirName )
            shutil.copytree( ttimeFileDir, path +s+ ttimeDirName )
        ttimeFileRoot =  ttimeDirName +s+ ttimeFileRoot.split( s )[-1]

        # Copy the output file directory to the specified path
        if utl.are_equal( outputFileDir, path +s+ outputDirName ) is False :
            if os.path.isdir( path +s+ outputDirName ) :
                shutil.rmtree( path +s+ outputDirName )
            shutil.copytree( outputFileDir, path +s+ outputDirName )
        outputFileRoot =  outputDirName +s+ outputFileRoot.split( s )[-1]

    # Get the locgrid list
    # list : [ xNum, yNum, zNum, xOrig, yOrig, zOrig, dx, dy, dz, gridType, saveFlag ]
        
    # Create an error message for the locgrid parameter
    locgrid_err_string = f"locgrid parameter must have 11 elements.\n"+\
                         f"It can be a list or a tuple containing:\n"+\
                         f"[ xNum, yNum, zNum, xOrig, yOrig, zOrig, dx, dy, dz, gridType, saveFlag ]\n"+\
                         f"Or it can be a string with the path of the input vel. model\n" +\
                         f"Or it can be a string with the path and name of the input vel. model\n"+\
                         f"The locgrid parameter you provide, is instead:\n"+\
                         f"{locgrid}"
    
    if type( locgrid ) in ( list, tuple ) :
        if len( locgrid ) < 11 :
            raise ValueError( locgrid_err_string )

    elif type( locgrid ) is str :

        if os.path.isfile( locgrid ) :
            locgrid = create_locgrid( locgrid )

        elif os.path.isfile( path +s+ locgrid ) :
            locgrid = create_locgrid( path +s+ locgrid )

        elif len( locgrid.split() ) in [ 10, 11 ] :
            locgrid = [ i for i in locgrid.split() ]

        else :
            raise ValueError( locgrid_err_string )
    else :
        raise ValueError( locgrid_err_string )

    # Check if the start_clean flag is True
    if start_clean is True :
        # Remove the content of the directory of output files
        utl.empty_directory( os.path.dirname( outputFileRoot ) )

    if not trans :
        model = None
        # Check model in the ttime dir. if the ttime dir. is an absolute path
        if os.path.isdir( ttimeFileDir ) :
            model = read_model_file( ttimeFileDir )
        # Check if model is not empty
        if model :
            trans = f"{model[0].proj_name} {model[0].orig_lat}" +\
                    f"{model[0].orig_lon} {model[0].map_rot}"
        # Check model in the ttime dir. if the ttime dir. is a relative path 
        else :
            model = read_model_file( path +s+ ttimeFileDir )
        if model :
            trans = f"{model[0].proj_name} {model[0].orig_lat} " +\
                    f"{model[0].orig_lon} {model[0].map_rot}"
        else :
            raise ValueError( f"Transformation parameters not found." )

    if printf is True :
        print( f"\n\n{print_prefix}{print_prefix}# -----------------------------------------" )
        print( f"Loc. parameters file :\n{print_prefix}{path +s+ name}\nCONTENT:" )

    with open( path +s+ name, 'w' ) as f :
        
        for l in param_header.splitlines() :
            
            if ("CONTROL" in l) and ("Control" not in l) and ( 'Syntax' not in l )  :
                l = f"CONTROL {control[0]} {control[1]}"

            if ("TRANS " in l) and ("Syntax" not in l) and ("TRANS - Geographic" not in l):
                if type( trans ) == str :
                    l = f"TRANS {trans}"
                if type( trans ) in ( list, tuple ) :
                    l = " ".join(str(w) for w in trans)
                
            f.write( l + ' \n' )

            if printf is True :
                if l.startswith( '#' ) == False :
                    print( print_prefix + l )

        for l in param_loc.splitlines() :
            
            if ( locsig is not None ) and ( "LOCSIG" in l ) :
                l = f"LOCSIG {locsig}"
            if ( loccom is not None ) and ( "LOCCOM" in l ) :
                l = f"LOCCOM {loccom}"
            if "LOCFILES" in l :
                l = f"LOCFILES {obsFile_rel} {obsFileType} {ttimeFileRoot} {outputFileRoot}"                 
            if ( "LOCHYPOUT" in l ) and ( lochypout is not None ) :
                l = f"LOCHYPOUT {lochypout}"  
            if "LOCGRID" in l :
                if type( locgrid ) in ( list, tuple ) :
                    locgrid = [ str(i) for i in locgrid if 'LOCGRID' not in str(i) ]
                    l = f"LOCGRID {' '.join( locgrid ) }"
                if type( locgrid ) is str :
                    if "LOCGRID" in l :
                        l = f"{locgrid}"
                    else :
                        l = f"LOCGRID {locgrid}"
                l = f"LOCGRID {' '.join( locgrid ) }"
            if "LOCSEARCH" in l :
                if type( locsearch ) in ( list, tuple ) :
                    locsearch = [ str(i) for i in locsearch if 'LOCSEARCH' not in str(i) ]
                    l = f"LOCSEARCH {' '.join( locsearch )}"
                if type( locsearch ) is str :
                    if "LOCSEARCH" in l :
                        l = f"{locsearch}"
                    else :
                        l = f"LOCSEARCH {locsearch}"
            if "LOCMETH" in l :
                if type( locmeth ) in ( list, tuple ) :
                    locmeth = [ str(i) for i in locmeth if 'LOCMETH' not in str(i) ]
                    l = f"LOCMETH {' '.join( locmeth )}"
                if type( locmeth ) is str :
                    if "LOCMETH" in l :
                        l = f"{locmeth}"
                    else :
                        l = f"LOCMETH {locmeth}"
            if "LOCGAU " in l :
                if type( locgau ) in ( list, tuple ) :
                    locgau = [ str(i) for i in locgau if 'LOCGAU' not in str(i) ]
                    l = f"LOCGAU {' '.join( locgau ) }"
                if type( locgau ) is str :
                    if "LOCGAU" in l :
                        l = f"{locgau}"
                    else :
                        l = f"LOCGAU {locgau}"
            if "LOCGAU2" in l :
                if locgau2 :
                    if type( locgau2 ) in ( list, tuple ) :
                        locgau2 = [ str(i) for i in locgau2 if "LOCGAU2" not in str(i) ]
                        l = f"LOCGAU2 {' '.join( locgau2 )}"
                    if type( locgau2 ) is str :
                        if "LOCGAU2" in l :
                            l = f"{locgau2}"
                        else :
                            l = f"LOCGAU2 {locgau2}"
            if ('LOCPHASEID' in l ) :
                if type( locphaseid ) in ( list, tuple ) :
                    locphaseid = [ str(i) for i in locphaseid if 'LOCPHASEID' not in str(i) ]
                    l = f"LOCPHASEID {' '.join( locphaseid )}"
                if type( locphaseid ) is str :
                    if "LOCPHASEID" in l :
                        l = f"{locphaseid}"
                    else :
                        l = f"LOCPHASEID {locphaseid}"
            if "LOCQUAL2ERR" in l :
                if locqual2err is not None :
                    locqual2err = [ str(i) for i in locqual2err if 'LOCQUAL2ERR' not in str(i) ]
                    l = f"LOCQUAL2ERR {' '.join( locqual2err )}"
                if type( locqual2err ) is str :
                    if "LOCQUAL2ERR" in l :
                        l = f"{locqual2err}"
                    else :
                        l = f"LOCQUAL2ERR {locqual2err}"
            if "LOCPHSTAT" in l :
                if locphstat is not None :
                    locphstat = [ str(i) for i in locphstat if 'LOCPHSTAT' not in str(i)]
                    l = f"LOCPHSTAT {' '.join( locphstat )}"
                if type( locphstat ) is str :
                    if "LOCPHSTAT" in l :
                        l = f"{locphstat}"
                    else :
                        l = f"LOCPHSTAT {locphstat}"
            if "LOCANGLES" in l :
                if locangles is not None :
                    locangles = [ str(i) for i in locangles if 'LOCANGLES' not in str(i)]
                    l = f"LOCANGLES {' '.join( locangles )}"
                if type( locangles ) is str :
                    if "LOCANGLES" in l :
                        l = f"{locangles}"
                    else :
                        l = f"LOCANGLES {locangles}"
            if "LOCMAG" in l :
                if locmag is not None :
                    locmag = [ str(i) for i in locmag if 'LOCMAG' not in str(i)]
                    l = f"LOCMAG {' '.join( locmag )}"
                if type( locmag ) is str :
                    if "LOCMAG" in l :
                        l = f"{locmag}"
                    else :
                        l = f"LOCMAG {locmag}"
            if ( "LOCELEVCORR " in l ) and ( locelevcorr != [] ) :
                if type( locelevcorr ) in ( list, tuple ) :
                    locelevcorr = [ str(i) for i in locelevcorr if 'LOCELEVCORR' not in str(i)]
                    l = f"LOCELEVCORR {' '.join( locelevcorr )}"
                if type( locelevcorr ) is str :
                    if "LOCELEVCORR" in l :
                        l = f"{locelevcorr}"
                    else :
                        l = f"LOCELEVCORR {locelevcorr}"
            if ( "LOCTOPO_SURFACE "  in l ) and ( loctopo_surface != [] ) :
                if type( loctopo_surface ) in ( list, tuple ) :
                    loctopo_surface = [ str(i) for i in loctopo_surface if 'LOCTOPO_SURFACE' not in str(i)]
                    l = f"LOCTOPO_SURFACE {' '.join( loctopo_surface )}"
                if type( loctopo_surface ) is str :
                    if "LOCTOPO_SURFACE" in l :
                        l = f"{loctopo_surface}"
                    else :
                        l = f"LOCTOPO_SURFACE {loctopo_surface}"

            f.write( l + ' \n' )

            if printf is True :
                if l.startswith( '#' ) == False :
                    print( print_prefix + l )
    f.close() 
    
    return path +s+ name

# -----------------------------------------------------------------------------
def create_tt_file( file_path_name, 
                       tt_dict=None, 
                       eve=None,
                       st_name=None, 
                       instrument=None, 
                       component=None,
                       p_phase_onset=None, 
                       phase_descriptor=None, 
                       first_motion=None,
                       yy=None, mm=None, 
                       dd=None, 
                       h=None, 
                       m=None, 
                       s=None, 
                       err=None, 
                       errmag=None, 
                       coda_duration=None, 
                       amplitude=None, 
                       period=None,
                       priorwt=None ) :

    if tt_dict is None :
        tt_dict = []
    
    if eve is None :
        if 'eve' in tt_dict :
            eve = tt_dict[ 'eve' ]

    if st_name is None :
        if 'st' in tt_dict :
            st_name = tt_dict[ 'st' ]
        if 'name' in tt_dict :
            st_name = tt_dict[ 'name' ]
    
    if instrument is None :
        if 'instrument' in tt_dict  :
            instrument = tt_dict[ 'instrument' ]
        else :
            instrument = np.full( np.size(st_name), '?' )

    if component is None :
        if 'component' in tt_dict  :
            component = tt_dict[ 'component' ]
        else :
            component = np.full( np.size(st_name), '?' )

    if p_phase_onset is None :
        if 'p_phase_onset' in tt_dict  :
            p_phase_onset = tt_dict[ 'p_phase_onset' ]
        else :
            p_phase_onset = np.full( np.size(st_name), '?' )

    if phase_descriptor is None :
        if 'phase_descriptor' in tt_dict  :
            phase_descriptor = tt_dict[ 'phase_descriptor' ]
        else :
            phase_descriptor = np.full( np.size(st_name), 'P' )

    if first_motion is None :
        if 'first_motion' in tt_dict  :
            first_motion = tt_dict[ 'first_motion' ]
        else :
            first_motion = np.full( np.size(st_name), '?' )

    if yy is None :
        if 'yy' in tt_dict  :
            yy = tt_dict[ 'yy' ]
        else :
            print( 'ERROR: no year (yy) set' )
            return 

    if mm is None :
        if 'mm' in tt_dict  :
            mm = tt_dict[ 'mm' ]
        else :
            print( 'ERROR: no month (mm) set' )
            return 

    if dd is None :
        if 'dd' in tt_dict  :
            dd = tt_dict[ 'dd' ]
        else :
            print( 'ERROR: no day (dd) set' )
            return 

    if h is None :
        if 'h' in tt_dict  :
            h = tt_dict[ 'h' ]
        else :
            print( 'ERROR: no hour (h) set' )
            return 

    if m is None :
        if 'm' in tt_dict  :
            m = tt_dict[ 'm' ]
        else :
            print( 'ERROR: no minute (m) set' )
            return 

    if s is None :
        if 's' in tt_dict  :
            s = tt_dict[ 's' ]
        else :
            print( 'ERROR: no second (s) set' )
            return

    if err is None :
        if 'err' in tt_dict  :
            err = tt_dict[ 'err' ]
        else :
            err = np.full( np.size(st_name), 'GAU' )

    if errmag is None :
        if 'errmag' in tt_dict  :
            errmag = tt_dict[ 'errmag' ]
        else :
            errmag = np.full( np.size(st_name), 0.01 )
    else :
        if np.size( errmag ) == 1 :
            errmag = np.full( np.size(st_name), errmag )

    if coda_duration is None :
        if 'coda_duration' in tt_dict  :
            coda_duration = tt_dict[ 'coda_duration' ]
        else :
            coda_duration = np.full( np.size(st_name), -1.0 )

    if amplitude is None :
        if 'amplitude' in tt_dict  :
            amplitude = tt_dict[ 'amplitude' ]
        else :
            amplitude = np.full( np.size(st_name), -1.0)

    if period is None :
        if 'period' in tt_dict  :
            period = tt_dict[ 'period' ]
        else :
            period = np.full( np.size(st_name), -1.0 )

    if priorwt is None :
        if 'priorwt' in tt_dict  :
            priorwt = tt_dict[ 'priorwt' ]
        else :
            priorwt = np.full( np.size(st_name), 1.0 )

    yymmdd = np.full( np.size(yy), '        ' )
    hm = np.full( np.size(yy), '    ' )
    for i, _ in enumerate( yy ) :
        yymmdd[i] = str( int(yy[i]) ) + str( int(mm[i]) ).rjust(2, '0') + str( int(dd[i]) ).rjust(2, '0')
        hm[i] = str( int(h[i]) ).rjust(2, '0') + str( int(m[i]) ).rjust(2, '0')
    
    array = np.column_stack( ( eve, st_name, instrument, component, p_phase_onset, 
                               phase_descriptor, first_motion, yymmdd, hm, np.round(s,4), 
                               err, errmag, coda_duration, amplitude, period, priorwt ) )

    with open( file_path_name, 'w' ) as f :
        for i, l in enumerate( array ) :
            if (i!=0) and ( l[0] != array[i-1,0] ) :
                f.write('\n')
            f.write( f"{l[1]:<6} {l[2]:<4} {l[3]:<4} {l[4]:<1} {l[5]:<6} {l[6]:<1} "+\
                     f"{l[7]:>10} {l[8]:>4} {l[9]:>6} {l[10]:<3} {l[11]:>9} "+\
                     f"{l[12]:>9} {l[13]:>9} {l[14]:>9} {l[15]:>9} \n")
    
    f.close()


    return file_path_name

# -----------------------------------------------------------------------------
def read_station_file( stations_file ) :
    """
    Reads a station file and returns a dictionary containing the station information.

    Parameters:
    stations_file (str): The path to the station file.

    Returns:
    dict: A dictionary containing the station information with the following keys:
        - 'st': List of station names.
        - 'label': List of station formats, i.e., "XYZ", "LONLAT" or "LATLON".
        - 'x': List of x-coordinates.
        - 'y': List of y-coordinates.
        - 'z': List of z-coordinates.
        - 'elev': List of elevation values.
        """

    with open( stations_file, 'r' ) as f :

        file_lines = f.readlines()

        f.close()

    first_line = file_lines[0].split()

    label = first_line[2]

    st_dict = { "st":[], "label":[], "x":[], "y":[], "z":[], "elev":[] }

    if ( label == "XYZ" ) or ( label == "LONLAT" ) :
        xc = 3
        yc = 4

    if label == "LATLON" :
        xc = 4
        yc = 3

    for l in file_lines :

        st_dict[ 'st' ].append( l.split()[1] )
        st_dict[ 'label' ].append( l.split()[2] )
        st_dict[ 'x' ].append( float( l.split()[xc] ) )
        st_dict[ 'y' ].append( float( l.split()[yc] ) )
        st_dict[ 'z' ].append( float( l.split()[5] ) )
        st_dict[ 'elev' ].append( float( l.split()[6] ) )

    for k in st_dict :
        st_dict[ k ] = np.array( st_dict[ k ] )

    return st_dict

# -----------------------------------------------------------------------------
def read_hyp_tt( hyp_file ):
    """
    Read travel-time data from an hypo71 file and return a dictionary.

    Parameters:

        - hyp_file (str): The path to the file containing the travel time data.

    Returns:

        dict: A dictionary containing the travel time data, with the following keys:
            - 'eve': List of event indices
            - 'st': List of station names
            - 'instrument': List of instrument names
            - 'component': List of component names
            - 'p_phase_onset': List of P-phase onsets
            - 'phase_descriptor': List of phase descriptors
            - 'first_motion': List of first motion values
            - 'yy': List of year values
            - 'mm': List of month values
            - 'dd': List of day values
            - 'h': List of hour values
            - 'm': List of minute values
            - 's': List of second values
            - 'err': List of error values
            - 'errmag': List of error magnitude values
            - 'coda_duration': List of coda duration values
            - 'amplitude': List of amplitude values
            - 'period': List of period values
            - 'priorwt': List of prior weight values
    """

    with open( hyp_file, 'r' ) as f :
        lines = f.readlines()
    f.close()

    tt_dict = { 'eve':[], 'st':[],  'instrument':[], 
                'component':[], 'p_phase_onset':[],
                'phase_descriptor':[], 'first_motion':[], 
                'yy':[], 'mm':[], 'dd':[], 
                'h':[], 'm':[], 's':[], 'err':[], 'errmag':[], 
                'coda_duration':[], 'amplitude':[], 'period':[], 'priorwt':[] }

    e = 0

    for i, l in enumerate( lines ) :
        if not l.strip() :
            e += 1
            continue
        print( l )
        print( l.split() )
        print( l.split()[6] )
        tt_dict['eve'].append( e )    
        tt_dict['st'].append( l.split()[0] )
        tt_dict['instrument'].append( l.split()[1] )
        tt_dict['component'].append( l.split()[2] )
        tt_dict['p_phase_onset'].append( l.split()[3] )
        tt_dict['phase_descriptor'].append( l.split()[4].strip() )
        tt_dict['first_motion'].append( l.split()[5] )
        tt_dict['yy'].append( int( l.split()[6][:4] ) )
        tt_dict['mm'].append( int( l.split()[6][4:6] ) )
        tt_dict['dd'].append( int( l.split()[6][6:] ) ) 
        tt_dict['h'].append( int( l.split()[7][:2] ) )
        tt_dict['m'].append( int( l.split()[7][2:] )  )
        tt_dict['s'].append( float( l.split()[8] ) )
        tt_dict['err'].append( l.split()[9] )
        tt_dict['errmag'].append( float( l.split()[10] ) ) 
        tt_dict['coda_duration'].append( float( l.split()[11] ) ) 
        tt_dict['amplitude'].append( float( l.split()[12] )  )
        tt_dict['period'].append( float( l.split()[13] )  )
        tt_dict['priorwt'].append( float( l.split()[14] )  )

    for k in tt_dict :

        tt_dict[ k ] = np.array( tt_dict[ k ] )

    return tt_dict

# -----------------------------------------------------------------------------
def read_model_file( model_data, plot=False, 
                     section=[ 0, None, None ],
                     vmin=None, vmax=None ) :
    """
    Reads a model file or directory of model files and 
    optionally plots the velocity model.
    
    Parameters:
     - model_data (str): Path to the model file or directory containing model files.
     - plot (bool, optional): If True, plots the velocity model. Default is False.
     - section (list, optional): Specifies the section of the model to plot. 
        Default is [0, None, None].
     - vmin (float, optional): Minimum value for the plot color scale. Default is None.
     - vmax (float, optional): Maximum value for the plot color scale. Default is None.
    
    Returns:
     - list: A list of NLLGrid objects representing the model grids.
    
    Notes:
     - If `model_data` is a file, it reads the file and creates a single NLLGrid object.
     - If `model_data` is a directory, it reads all files 
        with a '.mod' extension in the directory and creates NLLGrid objects for each.
     - If `plot` is True, it plots the velocity model using 
        the specified section and color scale limits.
    """

    grids = []
    file_names = []

    if os.path.isfile( model_data ) :
        model_path = model_data.split('.mod')[0] + '.mod'
        grids.append( NLLGrid( model_data ) )

    if os.path.isdir( model_data ) :
        for f in os.listdir( model_data ) :
            if  '.mod' in f :
                model_path = model_data +s+ f.split('.mod')[0] + '.mod'
                if model_path not in file_names :
                    file_names.append( model_path )
                    grids.append( NLLGrid( model_path ) )

    if plot is True :

        n = len( grids )

        for i, A in enumerate( grids ) :

            # Convert from SLOW_LEN to VELOCITY
            AA = ( 1/A.array ) * A.dx 

            ex = A.get_extent()

            if (section[0]==None) and (section[1]==None) :
                AA = np.rot90(AA[:,:,section[2]], k=-1)
                lim = [ex[0], ex[1], ex[2], ex[3]]
            if (section[0]==None) and (section[2]==None) :
                AA = np.rot90(AA[:,section[1],:], k=-1)
                lim = [ex[2], ex[3], -ex[5], -ex[4]]
            if (section[1]==None) and (section[2]==None) :
                AA = np.rot90(AA[section[0],:,:], k=-1 )
                lim = [ex[0], ex[1], -ex[5], -ex[4]]
            ax = utl.plta( AA, axis=True, lim=lim, adjst_lim=False, new_fig=False,
                           sbplt=[n,1,i+1], print_stat=False, label='km / s',
                           vmin=vmin, vmax=vmax )
            ax.set_title( f"{file_names[i].split(os.sep)[-1]}" )

        
        lsz_plot.plt.tight_layout()
            
    return grids

# -----------------------------------------------------------------------------
def create_locgrid( model=None, 
                    xNum=None, 
                    yNum=None, 
                    zNum=None, 
                    xOrig=None, 
                    yOrig=None, 
                    zOrig=None, 
                    step=None, 
                    gridType='PROB_DENSITY', 
                    saveFlag='SAVE',
                    reverse_z=True, 
                    gridDim='3d', 
                    printf=False, 
                    print_prefix='\t') :

    """
    Description
    ----------
    This functions create a list with the parameters describing 
    the initial or nested 3D search grid.  \n
    The order of LOCGRID statements is critical (see Notes). \n
    LOCGRID  xNum yNum zNum xOrig yOrig zOrig dx dy dz gridType saveFlag. \n

    Parameters
    ----------
    model : python dictionary with keys 'x', 'y', and 'z', corresponding to objects 
        (list, tuple, array or pandas series) containing the coordinates values 
        of each element of the model. Deafault is None.

    xNum yNum zNum (integer, min:2) : number of grid nodes in the x, y and z directions. 
        If None these values are taken from the argument "model". Default is None.

    xOrig yOrig zOrig (float) : x, y and z location of the grid origin in km relative to the geographic origin. 
        Use a large, negative value ( i.e. -1.0e30 ) to indicate automatic positioning.

        of grid along corressponding direction (valid for nested grids only, may not be used for initial grid).
        If None these are taken from the argument "model". Default is None.

    dx dy dz (float) : grid node spacing along the x, y and z axes (kilometers units)
        If None these valus are taken from the argument "model". Default is None.

    gridType (choice: MISFIT PROB_DENSITY) : statistical quantity to calculate on grid
        Default is "PROB_DENSITY"

    saveFlag (choice: SAVE NO_SAVE) : specifies if the results of the search over this grid should be saved to disk

    gridDim (choice: '3d' or '2d') : specifies if the grid is 3D or 2D. Default is 3D. 
        It is used to calculate the grid dimensions if the argument "model" is a dictionary. 

    Returns
    -------
    list : [ xNum, yNum, zNum, xOrig, yOrig, zOrig, dx, dy, dz, gridType, saveFlag ]

    Notes
    -------
    1. The first LOCGRID is the initial search grid which may not have automatic positionig along any axes. 
        The succeeding LOCGRID statements may specify automatic positioning along one or more axes ( xOrig, yOrig, zOrig = -1.0e30 ), 
        but must all be sized ( i.e. dx*(xNum-1) , etc.) so that they can be fully contained within the preceeding grid. 
        The NLLoc program will attempt to translate a nested grid that intersects a boundary of the initial grid 
        so that it is contained inside of the initial grid; if this is not possible the location will be terminated prematurely.
    
    2. With automatic positioning ( xOrig, yOrig, zOrig = -1.0e30 ), a grid is shifted in x/y/z 
        so that it is centered on the minimum misfit hypocenter x/y/z of the preceeding grid.
    
    3. Each search over a grid with gridType = PROB_DENSITY is time consuming 
        and should generally only be used for a nested grid on which the full PDF is required and will be saved to disk. 
            Use gridType = MISFIT for the initial grid, for larger nested grids, 
            and for smaller nested grids in maximum-likelihood hypocenter searches ( i.e. where the PDF is not if interest).
    
    4. The 3D grid dimensions are in kilometers with Z positive down (left-handed coordinate system).
    
    5. The grid is dx*(xNum-1) km long in the x direction, and similarly for y and z.
    
    6. For 2D velocity and travel-time grids LOCGRID should be 3D and positioned absolutely in space, 
        thus xNum >> 2 and xOrig and zOrig are in general != 0.0
    """
    
    # Check if a model is provided
    if model is not None :

        # If the model is a dictionary
        if type( model ) == dict :

            # Get unique values of x, y, and z from the model
            xu = np.unique( model['x'] )
            yu = np.unique( model['y'] )
            zu = np.unique( model['z'] )

            # Calculate the minimum difference between the unique values
            step_test = np.abs( np.min( ( np.diff(xu).min(), np.diff(yu).min(), np.diff(zu).min() ) ) ) 

            # If the step size is not equal to the minimum difference
            if float( round(step,3 ) ) != float( round( step_test,3 ) ) : 
                # Copy the model
                model = model.copy()
                # If step is not defined, set it to the minimum difference
                if step is None :
                    step = step_test  
                # Generate new x, y, and z arrays with the defined step size
                xn = np.arange( model['x'].min(), model['x'].max(), step )
                yn = np.arange( model['y'].min(), model['y'].max(), step )
                zn = np.arange( model['z'].min(), model['z'].max(), step )   
                # Create a meshgrid and flatten it to 1D arrays
                X, Y, Z = np.meshgrid( np.round(xn), np.round(yn), np.round(zn) ) 
                model['x'] = X.ravel()
                model['y'] = Y.ravel()
                model['z'] = Z.ravel()
            
            # Get unique values of x, y, and z from the model
            xu = np.unique( model['x'] )
            yu = np.unique( model['y'] )
            zu = np.unique( model['z'] )
            
            # If reverse_z is True, multiply zu by -1
            if reverse_z is True :
                zu *= -1

            # If xNum, yNum, or zNum are not defined, set them to the length of the unique values
            if xNum is None :
                xNum = len( xu )
            if yNum is None :    
                yNum = len( yu )
            if zNum is None :
                zNum = len( zu )

            # If xOrig, yOrig, or zOrig are not defined, set them to the minimum of the unique values
            if xOrig is None :
                xOrig = xu.min()
            if yOrig is None :
                yOrig = yu.min()
            if zOrig is None :
                zOrig = zu.min()

            # If gridDim is '3d', calculate xNum, yNum, and zNum based on the maximum unique values, the origins, and the step size
            if gridDim=='3d' :

                xNum = int(np.floor( ( xu.max() - xOrig ) / step ))
                yNum = int(np.floor( ( yu.max() - yOrig ) / step ))
                zNum = int(np.floor( ( zu.max() - zOrig ) / step ))

            # If gridDim is '2d', calculate xNum, yNum, and zNum based on the maximum unique values, the origins, and the step size
            if gridDim=='2d' :
                xNum = 2
                yNum = np.max( ( int(np.floor( ( yu.max() - yOrig ) / step )), 
                                int(np.floor( ( xu.max() - xOrig ) / step )) ) ) 
                zNum = int(np.floor( ( zu.max() - zOrig ) / step ))

            # Print xNum, yNum, zNum, and step
            print( print_prefix, xNum, yNum, zNum, step )

        # If the model is a string
        elif type( model ) == str :

            # Read the model
            model = read_model_file( model )
            model = model[0]
            # Get xNum, yNum, zNum, xOrig, yOrig, zOrig, stepx, stepy, stepz from the model
            xNum = model.nx
            yNum = model.ny
            zNum = model.nz
            xOrig = model.x_orig
            yOrig = model.y_orig
            zOrig = model.z_orig
            stepx = model.dx
            stepy = model.dy
            stepz = model.dz
            # Calculate the mean step size
            step = np.mean( ( stepx, stepy, stepz ) )

        # If the model is a NLLGrid object
        elif str(type( model )).strip() == "<class 'nllgrid.NLLGrid.NLLGrid'>" :
            # Get xNum, yNum, zNum, xOrig, yOrig, zOrig, stepx, stepy, stepz from the model
            model = model
            xNum = model.nx
            yNum = model.ny
            zNum = model.nz
            xOrig = model.x_orig
            yOrig = model.y_orig
            zOrig = model.z_orig
            stepx = model.dx
            stepy = model.dy
            stepz = model.dz
            step = np.mean( ( stepx, stepy, stepz ) )

    # Convert xNum, yNum, and zNum to integers
    xNum = int( xNum )
    yNum = int( yNum )
    zNum = int( zNum )
    # Create a list with the grid parameters
    locgrid_list = [ xNum, yNum, zNum, xOrig, 
                     yOrig, zOrig, step, step, step, 
                     gridType, saveFlag ]

    # Create a list with the names of the grid parameters
    strpar = [ "xNum", "yNum", "zNum", "xOrig", "yOrig", 
            "zOrig", "dx", "dy", "dz", "gridType", "saveFlag" ] 

    # If printf is True, print the grid parameters
    if printf is True :
        print(f"\n\n{print_prefix}# -----------------------------------------")
        print( f"{print_prefix}LOCGRID:" )
        for i, n in enumerate( strpar ) :
            print( print_prefix, n,"=", locgrid_list[i] )

    # Return the list with the grid parameters
    return locgrid_list

# -----------------------------------------------------------------------------
def read_hyp_phs( file, fn=1e2, print_n=False ) :

    with open( file, 'r' ) as f :

        lines = f.readlines()

    n = 1
    c = np.empty( (0) )
    phs = {'eve':c, 'sts':c, 'yy':c, 'mm':c, 'dd':c, 'h':c, 'm':c, 'sp':c, 'ss':c, 'wp':c, 
           'ws':c, 'phsp':c, 'phss':c, 'rp':c, 'fmp':c, 'rs':c, 'fms':c }
    
    for i, l in enumerate( lines ) :

        if l[ 0 : 5 ].split() == [] :
            n += 1
            continue

        if l[ 0 : 5 ].split() != [] :
            phs['sts'] = np.append( phs['sts'], l[ 0 : 4 ].split()[0] )
        else : 
            phs['sts'] = np.append( phs['sts'], np.nan ) 

        if l[ 4 : 5  ] != ' ' :
            phs['rp'] = np.append( phs['rp'], l[ 4 : 5 ] )
        else : 
            phs['rp'] = np.append( phs['rp'], np.nan )

        if l[ 5 : 6  ] != ' ' :
            phs['phsp'] = np.append( phs['phsp'], l[ 5 : 6 ] )
        else : 
            phs['phsp'] = np.append( phs['phsp'], np.nan )

        if l[ 6 : 7  ] != ' ' :
            phs['fmp'] = np.append( phs['fmp'], l[ 6 : 7 ] )
        else : 
            phs['fmp'] = np.append( phs['fmp'], np.nan )

        if l[ 7 : 8 ].split() != [] :
            phs['wp'] = np.append( phs['wp'], int( l[ 7 : 8 ] ) )  
        else : 
            phs['wp'] = np.append(  phs['wp'], 0 ) 
        
        if l[ 9 : 11 ].split() != [] :
            phs['yy'] = np.append( phs['yy'] , float( l[ 9 : 11 ].split()[0] ) )  
        else : 
            phs['yy'] = np.append(  phs['yy'], np.nan )              
        
        if l[ 11 : 13 ].split() != [] :
            phs['mm'] = np.append( phs['mm'] , float( l[ 11 : 13 ].split()[0] ) ) 
        else : 
            phs['mm'] = np.append( phs['mm'], np.nan )             
        
        if l[ 13 : 15 ].split() != [] :
            phs['dd'] = np.append( phs['dd'] , float( l[ 13 : 15 ].split()[0] ) )  
        else : 
            phs['dd'] = np.append( phs['dd'], np.nan )             
        
        if l[ 15 : 17 ].split() != [] :
            phs['h'] = np.append( phs['h'] , float( l[ 15 : 17 ].split()[0] ) ) 
        else : 
            phs['h'] = np.append(  phs['h'], np.nan )           
        
        if l[ 17 : 19 ].split() != [] :
            phs['m'] = np.append( phs['m'] , float( l[ 17 : 19 ].split()[0] ) )   
        else : 
            phs['m'] = np.append(  phs['m'], np.nan )                 
        
        if l[ 19 : 24 ].split() != [] :
            phs['sp'] = np.append( phs['sp'] , float( l[ 19 : 24 ].split()[0] )/fn )  
        else : 
            phs['sp'] = np.append(  phs['sp'], np.nan )            
        
        if l[ 31 : 36 ].split() != [] :
            phs['ss'] = np.append( phs['ss'] , float( l[ 31 : 36 ].replace(' ', '0') )/fn )    
        else : 
            phs['ss'] = np.append(  phs['ss'], np.nan )   

        if l[ 36 : 37  ] != ' ' :
            phs['rs'] = np.append( phs['rs'], l[ 36 : 37 ] )
        else : 
            phs['rs'] = np.append( phs['rs'], np.nan )

        if l[ 37 : 38 ] != ' ' :
            phs['phss'] = np.append( phs['phss'] , l[ 37 : 38 ] )    
        else : 
            phs['phss'] = np.append(  phs['phss'], np.nan )               

        if l[ 38 : 39 ] != ' ' :
            phs['fms'] = np.append( phs['fms'] , l[ 38 : 39 ] )    
        else : 
            phs['fms'] = np.append(  phs['fms'], np.nan )

        if l[ 39 : 41].split() != [] :
            phs['ws'] = np.append( phs['ws'] , float( l[ 39 : 41].split()[0] ) )  
        else : 
            phs['ws'] = np.append(  phs['ws'], 0 )    

        phs['eve'] = np.append( phs['eve'] , n ) 

    for k in phs :
        phs[ k ] = np.array( phs[k] ) 
    
    if print_n is True :    
        print( f'N° of events : {int( phs["eve"][-1] )} ' )
        print( f'N° of records : {phs["eve"].size} ' )
        print( f'P waves phases : {np.sum(np.isfinite(phs["sp"]))} ' )
        print( f'S waves phases : {np.sum(np.isfinite(phs["ss"]))} ' )
    
    return phs

# -----------------------------------------------------------------------------
def read_hyp_out( file, plot=False, aspect='equal', size=1, reverse_z=True, reverse_x=True ) :

    f = open( file, 'r' )

    lines = f.readlines()

    hyp_out = { 'eve':[], 'yy':[], 'mm':[], 'dd':[], 'h':[], 'm':[], 's':[], 
                'lat':[], 'lon':[], 'z':[], 'lon_deg':[], 'lon_min':[], 'lat_deg':[], 
                'lat_min':[], 'mag':[], 'no':[], 'rmsres':[], 'x':[], 'y':[],
                'gap':[], 'dmin':[], 'rzdm':[], 'np':[], 'ns':[] } 
                
    for i, l in enumerate( lines ) :
        
        if i < 5 :
            continue
        hyp_out['eve'].append( int( l[:4] ) )
        hyp_out['yy'].append( int( l[4:7] ) )
        hyp_out['mm'].append( int( l[7:9] ) )
        hyp_out['dd'].append( int( l[9:11] ) )
        hyp_out['h'].append( int( l[11:14] ) )
        hyp_out['m'].append( int( l[14:16] ) )
        hyp_out['s'].append( float( l[16:22] )  )
        hyp_out['lat_deg'].append( int( l[22:25] ) ) 
        hyp_out['lat_min'].append( float( l[26:32] ) )
        hyp_out['lon_deg'].append( int( l[32:36] ) ) 
        hyp_out['lon_min'].append( float( l[37:43] ) )     
        hyp_out['z'].append( float( l[43:50] ) ) 
        hyp_out['mag'].append( float( l[50:57] ) )
        hyp_out['no'].append( int( l[57:61] ) )
        hyp_out['rmsres'].append( float( l[61:67] ) ) 
        hyp_out['x'].append( float( l[67:74] ) ) 
        hyp_out['y'].append( float( l[74:81] ) )     
   
    for k in hyp_out :
        hyp_out[ k ] = np.array( hyp_out[ k ] )
        
     
    hyp_out['lon'] = hyp_out['lon_deg'] + hyp_out['lon_min'] / 60
    hyp_out['lat'] = hyp_out['lat_deg'] + hyp_out['lat_min'] / 60
        
    if reverse_z is True :
        hyp_out['z'] = hyp_out['z'] * -1
        
    if reverse_x is True :
        hyp_out['x'] = hyp_out['x'] * -1
    
    f.close()
    
    if plot == True :
        
        f = lsz_plot.plt.figure()
        
        lsz_plot.plt.subplot( 311, aspect=aspect )
        lsz_plot.plt.scatter( hyp_out['x'], hyp_out['y'], s=size )
        lsz_plot.plt.title( 'X-Y axis' )
        
        lsz_plot.plt.subplot( 312, aspect=aspect )
        lsz_plot.plt.scatter( hyp_out['x'], hyp_out['z'], s=size ) 
        lsz_plot.plt.title( 'X-Z axis' )
        
        lsz_plot.plt.subplot( 313, aspect=aspect )
        lsz_plot.plt.scatter( hyp_out['y'], hyp_out['z'], s=size )
        lsz_plot.plt.title( 'Y-Z axis' )     
        
        f.tight_layout()
        
    return hyp_out

# -----------------------------------------------------------------------------
def read_tt_file( tt_file ) :
    """
    Read a travel-time file or observation file (generally '.obs') saved in Nlloc format, 
    and return a dictionary containing the parsed data.

    Parameters:
    tt_file (str): The path to the travel time file.

    Returns:
    dict: A dictionary containing the parsed travel time data. 
        The dictionary has the following keys:
        - 'eve': Numpy array of event indices.
        - 'st': Numpy array of station names.
        - 'instrument': Numpy array of instrument names.
        - 'component': Numpy array of component names.
        - 'p_phase_onset': Numpy array of P-phase onsets.
        - 'phase_descriptor': Numpy array of phase descriptors.
        - 'first_motion': Numpy array of first motion values.
        - 'yy': Numpy array of year values.
        - 'mm': Numpy array of month values.
        - 'dd': Numpy array of day values.
        - 'h': Numpy array of hour values.
        - 'm': Numpy array of minute values.
        - 's': Numpy array of second values.
        - 'err': Numpy array of error values.
        - 'errmag': Numpy array of error magnitude values.
        - 'coda_duration': Numpy array of coda duration values.
        - 'amplitude': Numpy array of amplitude values.
        - 'period': Numpy array of period values.
        - 'priorwt': Numpy array of prior weight values.
        - 'dt': Numpy array of combined date and time values.

        N.B. All arrays in the dictionary have one dimension and the same length.
    """

    with open( tt_file, 'r' ) as f :
        lines = f.readlines()
    f.close()

    tt_dict = { 'eve':[], 'st':[],  'instrument':[], 'component':[], 'p_phase_onset':[],
               'phase_descriptor':[], 'first_motion':[], 'yy':[], 'mm':[], 'dd':[], 
               'h':[], 'm':[], 's':[], 'err':[], 'errmag':[], 'coda_duration':[],
               'amplitude':[], 'period':[], 'priorwt':[] }

    e = 0

    for i, l in enumerate( lines ) :
        if not l.strip() :
            e += 1
            continue
        tt_dict['eve'].append( e )
        tt_dict['st'].append( l.split()[0] )
        tt_dict['instrument'].append( l.split()[1] )
        tt_dict['component'].append( l.split()[2] )
        tt_dict['p_phase_onset'].append( l.split()[3] )
        tt_dict['phase_descriptor'].append( l.split()[4].strip() )
        tt_dict['first_motion'].append( l.split()[5] )
        tt_dict['yy'].append( int( l.split()[6][:4] ) )
        tt_dict['mm'].append( int( l.split()[6][4:6] ) )
        tt_dict['dd'].append( int( l.split()[6][6:] ) ) 
        tt_dict['h'].append( int( l.split()[7][:2] ) )
        tt_dict['m'].append( int( l.split()[7][2:] )  )
        tt_dict['s'].append( float( l.split()[8] ) )
        tt_dict['err'].append( l.split()[9] )
        tt_dict['errmag'].append( float( l.split()[10] ) ) 
        tt_dict['coda_duration'].append( float( l.split()[11] ) ) 
        tt_dict['amplitude'].append( float( l.split()[12] )  )
        tt_dict['period'].append( float( l.split()[13] )  )
        tt_dict['priorwt'].append( float( l.split()[14] )  )

    for k in tt_dict :
        tt_dict[ k ] = np.array( tt_dict[ k ] )

    tt_dict['dt'] = utl.combine64( years = tt_dict['yy'], 
                                    months = tt_dict['mm'], 
                                    days = tt_dict['dd'], 
                                    hours = tt_dict['h'], 
                                    minutes = tt_dict['m'], 
                                    seconds = tt_dict['s'] )

    return tt_dict

# -----------------------------------------------------------------------------
def read_pun_file( file, plot=False, aspect='equal', size=1, reverse_z=True ) :

    f = open( file, 'r' )

    lines = f.readlines()

    hyp_out = { 'eve':[], 'yy':[], 'mm':[], 'dd':[], 
                'h':[], 'm':[], 's':[], 
                'lat':[], 'lon':[], 'z':[], 
                'lon_deg':[], 'lon_min':[], 'lat_deg':[], 
                'lat_min':[], 'mag':[], 'no':[], 'rms':[],
                'gap':[], 'dmin':[], 'erh':[], 'erz':[], 'qm':[] } 
                
    for i, l in enumerate( lines ) :
        
        if i == 0 :
            continue
        hyp_out['eve'].append( int( i ) )
        hyp_out['yy'].append( int( l[1:3] ) )
        hyp_out['mm'].append( int( l[3:5] ) )
        hyp_out['dd'].append( int( l[5:7] ) )
        hyp_out['h'].append( int( l[8:10] ) )
        hyp_out['m'].append( int( l[10:12] ) )
        hyp_out['s'].append( float( l[12:18].replace(' ', '0') ) )       
        hyp_out['lat_deg'].append( int( l[19:21] ) ) 
        hyp_out['lat_min'].append( float( l[22:27] ) )
        hyp_out['lon_deg'].append( int( l[29:31] ) )   
        hyp_out['lon_min'].append( float( l[32:37] ) )             
        hyp_out['z'].append( float( l[37:44] ) ) 
        hyp_out['mag'].append( float( l[44:51] ) )
        hyp_out['no'].append( int( l[51:54] ) )
        hyp_out['gap'].append( float( l[54:58] ) ) 
        hyp_out['dmin'].append( float( l[58:63] ) ) 
        hyp_out['rms'].append( float( l[63:68] ) ) 
        hyp_out['erh'].append( float( l[68:73] ) ) 
        hyp_out['erz'].append( float( l[73:78] ) ) 
        hyp_out['qm'].append( l[78:81] )  
   
    for k in hyp_out :
        hyp_out[ k ] = np.array( hyp_out[ k ] )
         
    hyp_out['lon'] = hyp_out['lon_deg'] + hyp_out['lon_min'] / 60
    hyp_out['lat'] = hyp_out['lat_deg'] + hyp_out['lat_min'] / 60         
        
    if reverse_z is True :
        hyp_out['z'] = hyp_out['z'] * -1     
    
    f.close()
    
    if plot == True :
        
        f = lsz_plot.plt.figure()
        
        lsz_plot.plt.subplot( 313, aspect=aspect )
        lsz_plot.plt.scatter( hyp_out['lon'], hyp_out['lat'], s=size )
        lsz_plot.plt.title( 'X-Y axis' )
        
        lsz_plot.plt.subplot( 311, aspect=aspect )
        lsz_plot.plt.scatter( hyp_out['lon'], hyp_out['z'], s=size ) 
        lsz_plot.plt.title( 'X-Z axis' )
        
        lsz_plot.plt.subplot( 312, aspect=aspect )
        lsz_plot.plt.scatter( hyp_out['lat'], hyp_out['z'], s=size )
        lsz_plot.plt.title( 'Y-Z axis' )     
        
        f.tight_layout()
        
    return hyp_out

# -----------------------------------------------------------------------------
def nlloc2hyp( tt_data, file, printf=False, locqual2err=[ 0.01, 0.02, 0.08, 0.2, 99999.9 ] ):
    """
    Convert travel-time data from Nlloc to Hypo71 observation file format.
    
    Args:
        - tt_data (str or dict): The travel time data. 
            It can be either a string representing the file path to the travel time data file, 
            or a dictionary containing the travel time data.
        
        - file (str): The file path to save the converted Hypo71 observation file.
        
        - printf (bool, optional): Whether to print the file path of the converted file. 
            Defaults to False.
        
        - locqual2err (list, optional): A list of error magnitude thresholds 
            for assigning weights to travel time data. 
            Defaults to [0.01, 0.02, 0.08, 0.2, 99999.9].
    
    Returns:
        str: The file path of the converted Hypo71 observation file.
    """

    path = os.path.dirname( file )
    file_name = os.path.basename( file )
    os.makedirs( path, exist_ok=True )

    if ( type( tt_data ) == str ) :

        tt_dict = read_tt_file( tt_data )[0]

    else :

        tt_dict = tt_data.copy()

    tt_dict[ 'w' ] = np.zeros( tt_dict[ 'errmag' ].shape ).astype( int )
    idx0 = tt_dict[ 'errmag' ] <= locqual2err[0]
    idx1 = ( tt_dict[ 'errmag' ] > locqual2err[0] ) & ( tt_dict[ 'errmag' ] <= locqual2err[1] )
    idx2 = ( tt_dict[ 'errmag' ] > locqual2err[1] ) & ( tt_dict[ 'errmag' ] <= locqual2err[2] )
    idx3 = ( tt_dict[ 'errmag' ] > locqual2err[2] ) & ( tt_dict[ 'errmag' ] <= locqual2err[3] )
    idx4 = ( tt_dict[ 'errmag' ] > locqual2err[3] )
    tt_dict[ 'w' ][idx0] = 0
    tt_dict[ 'w' ][idx1] = 1
    tt_dict[ 'w' ][idx2] = 2
    tt_dict[ 'w' ][idx3] = 3
    tt_dict[ 'w' ][idx4] = 4

    if 'mag' not in tt_dict :
        tt_dict['mag'] = np.full( tt_dict['eve'].shape, 3 )

    lines = []

    idxs = np.core.defchararray.find(tt_dict['phase_descriptor'], 'S') != -1

    for i, e in enumerate( tt_dict['eve'] ) :

        if ( i != 0 ) and ( tt_dict['eve'][i] != tt_dict['eve'][i-1] ) :
            line = '\n'
            lines.append(line)
            continue

        if "S" in tt_dict['phase_descriptor'][i]:
            continue

        else :

            eve_i = tt_dict['eve'][i]
            st_name_i = tt_dict['st'][i]
            si = (tt_dict['eve'] == eve_i) &\
                 (tt_dict['st'] == st_name_i) &\
                 (idxs)

            sn = "{:<4}".format(st_name_i)
            ppo = tt_dict['p_phase_onset'][i] if '?' not in tt_dict['p_phase_onset'][i] else " "
            pfm = tt_dict['first_motion'][i] if '?' not in tt_dict['first_motion'][i] else " "

            pw = f"{tt_dict['w'][i]}"

            pyy = str( tt_dict['yy'][i] )
            if len(pyy) < 2 :
                pyy = pyy.zfill(2) 
            if len(pyy) > 2 :
                pyy = pyy[-2:] 
            pmm = str( tt_dict['mm'][i] )
            if len(pmm) < 2 :
                pmm = pmm.zfill(2)
            pdd = str( tt_dict['dd'][i] )
            if len(pdd) < 2 :
                pdd = pdd.zfill(2)
            ph = str( tt_dict['h'][i] )
            if len(ph) < 2 :
                ph = ph.zfill(2)
            pm = str( tt_dict['m'][i] )
            if len(pm) < 2 :
                pm = pm.zfill(2)

            ps = f"{tt_dict['s'][i]: >5.2f}"

            line = f"{sn}{ppo}P{pfm}{pw} {pyy}{pmm}{pdd}{ph}{pm}{ps}"

            if np.sum(si) != 0:
                line = line + f"       {float(tt_dict['s'][si]): >5.2f} S {int(tt_dict['w'][si])}"

            lines.append(line)

    with open(path +s+ file_name, 'w') as f:
        f.write('\n'.join(lines))
    f.close()
    
    if printf == True :
        print( f"Hypo71 obs. file:\n{path +s+ file_name}" )

    return path +s+ file_name

# -----------------------------------------------------------------------------
def hyp2nlloc( hyp_data, 
               file_name=None, 
               path=os.getcwd(), 
               yy1=1900, 
               yy2=2000, 
               locqual2err=[0.01, 0.02, 0.08, 0.2, 99999.9], 
               fn=1 ) :
    """
    Convert Hypo71 format travel-time data to NLLoc format.

    Args:

        - hyp_data (str or dict): The input Hypoinverse data. 
            If it's a string, it should be the path to the Hypoinverse data file. 
            If it's a dictionary, it should be the Hypoinverse data in dictionary format.
        
        - file_name (str, optional): The name of the output file. 
            If not provided, the data will not be saved to a file. Defaults to None.
       
        - path (str, optional): The path where the output file will be saved. Defaults to ''.
        
        - yy1 (int, optional): The year offset for years greater than 50. Defaults to 1900.
        
        - yy2 (int, optional): The year offset for years less than or equal to 50. 
            Defaults to 2000.
        
        - locqual2err (list, optional): A list of error values corresponding 
            to different location qualities. Defaults to [0.01, 0.02, 0.08, 0.2, 99999.9].
        
        - fn (int, optional): The format number of the Hypoinverse data file. Defaults to 1.

    Returns:

        dict: The converted NLLoc format data in dictionary format.
    """

    if type( hyp_data ) == str : 
        hyp_dict = read_hyp_tt( hyp_data, fn=fn )
    else :
        hyp_dict = utl.copy.deepcopy( hyp_data )

    tt_dict = { 'eve':[], 'st':[],  'instrument':[], 
                'component':[], 'p_phase_onset':[],
                'phase_descriptor':[], 'first_motion':[], 
                'yy':[], 'mm':[], 'dd':[], 
                'h':[], 'm':[], 's':[], 'err':[], 
                'errmag':[], 'coda_duration':[],
                'amplitude':[], 'period':[], 'priorwt':[] }
    
    idx = hyp_dict[ 'yy' ] > 50 
    hyp_dict[ 'yy' ][idx] += yy1
    hyp_dict[ 'yy' ][~idx] += yy2

    datetim = utl.combine64( years=hyp_dict['yy'], months=hyp_dict['mm'], 
                             days=hyp_dict['dd'], 
                             hours=hyp_dict['h'], minutes=hyp_dict['m'], 
                             seconds=hyp_dict['ss'] )

    for i,_ in enumerate( hyp_dict['sts'] ) :

        tt_dict['eve'].append( int( hyp_dict['eve'][i] ) )
        tt_dict['st'].append( hyp_dict['sts'][i] )
        tt_dict['instrument'].append( '?' )
        tt_dict['component'].append( '?' )
        if str( hyp_dict['rp'][i] ) == 'nan' :  
            tt_dict['p_phase_onset'].append( '?' )
        else :
            tt_dict['p_phase_onset'].append( hyp_dict['rp'][i] )
        if str( hyp_dict['fmp'][i] ) == 'nan' :
            tt_dict['first_motion'].append( '?' )
        else :
            tt_dict['first_motion'].append( hyp_dict['fmp'][i] )
        tt_dict['yy'].append( int( hyp_dict['yy'][i] ) )
        tt_dict['mm'].append( int( hyp_dict['mm'][i] ) )
        tt_dict['dd'].append( int( hyp_dict['dd'][i] ) )
        tt_dict['h'].append( int( hyp_dict['h'][i] ) )
        tt_dict['m'].append( int( hyp_dict['m'][i] ) )
        tt_dict['s'].append( hyp_dict['sp'][i] )
        tt_dict['phase_descriptor'].append( 'P' )
        tt_dict['err'].append( 'GAU' )
        tt_dict['coda_duration'].append( '-1' )
        tt_dict['amplitude'].append( '-1' )
        tt_dict['period'].append( '-1' )
        if hyp_dict['wp'][i] == 0 :
            tt_dict[ 'errmag' ].append( locqual2err[0] )
        if hyp_dict['wp'][i] == 1 :
            tt_dict[ 'errmag' ].append( locqual2err[1] )
        if hyp_dict['wp'][i] == 2 :
            tt_dict[ 'errmag' ].append( locqual2err[2] )
        if hyp_dict['wp'][i] == 3 :
            tt_dict[ 'errmag' ].append( locqual2err[3] )
        if hyp_dict['wp'][i] == 4 :
            tt_dict[ 'errmag' ].append( locqual2err[4] )
            tt_dict['priorwt'].append( 0.0 ) 
        if hyp_dict['wp'][i] != 4 :
           tt_dict['priorwt'].append( 1.0 )  
        if hyp_dict['wp'][i] > 4 :
            tt_dict[ 'errmag' ].append( locqual2err[4] ) 

        if np.isfinite( hyp_dict['ss'][i] ) :

            datetim_ss = datetim[i].astype(str).replace('-',' ').replace('T',' ').replace(':',' ').split()
            datetim_yy = int( datetim_ss[0] )
            datetim_mm = int( datetim_ss[1] )
            datetim_dd = int( datetim_ss[2] )
            datetim_h = int( datetim_ss[3] )
            datetim_m = int( datetim_ss[4] )
            datetim_s = round( float( datetim_ss[5] ), 2 )

            tt_dict['eve'].append( int( hyp_dict['eve'][i] ) )
            tt_dict['st'].append( hyp_dict['sts'][i] )
            tt_dict['instrument'].append( '?' )
            tt_dict['component'].append( '?' )
            if str( hyp_dict['rs'][i] ) == 'nan' :
                tt_dict['p_phase_onset'].append( '?' )
            else :
                tt_dict['p_phase_onset'].append( hyp_dict['rs'][i] )
            if str( hyp_dict['fms'][i] )  == 'nan' :
                tt_dict['first_motion'].append( '?' )
            else :
                tt_dict['first_motion'].append( hyp_dict['fms'][i] )
            tt_dict['yy'].append( datetim_yy )
            tt_dict['mm'].append( datetim_mm )
            tt_dict['dd'].append( datetim_dd )
            tt_dict['h'].append( datetim_h )
            tt_dict['m'].append( datetim_m )
            tt_dict['s'].append( datetim_s ) 
            tt_dict['err'].append( 'GAU' )
            tt_dict['coda_duration'].append( '-1' )
            tt_dict['amplitude'].append( '-1' )
            tt_dict['period'].append( '-1' )
            tt_dict['phase_descriptor'].append( 'S' )

            if hyp_dict['ws'][i] == 0 :
                tt_dict[ 'errmag' ].append( locqual2err[0] )
            if hyp_dict['ws'][i] == 1 :
                tt_dict[ 'errmag' ].append( locqual2err[1] )
            if hyp_dict['ws'][i] == 2 :
                tt_dict[ 'errmag' ].append( locqual2err[2] )
            if hyp_dict['ws'][i] == 3 :
                tt_dict[ 'errmag' ].append( locqual2err[3] )
            if hyp_dict['ws'][i] == 4 :
                tt_dict[ 'errmag' ].append( locqual2err[4] )
                tt_dict['priorwt'].append( 0.0 ) 
            if hyp_dict['ws'][i] != 4 :
                tt_dict['priorwt'].append( 1.0 )
            if hyp_dict['ws'][i] > 4 :
                tt_dict[ 'errmag' ].append( locqual2err[4] )

    for k in list(tt_dict.keys()) :
        tt_dict[ k] = np.array( tt_dict[ k] )
        if np.size( tt_dict[k] ) == 0 :
            del tt_dict[k]

    if file_name is not None : 
        _ = create_tt_file( path +s+ file_name, tt_dict )

    return tt_dict 

# -----------------------------------------------------------------------------
def run_nlloc( path=os.getcwd(),
               vel2grid=False,
               grid2time=False, 
               nlloc=False,
               stFile=None, 
               obsFile=None, 
               modelFileRoot='model'+s+'layer',
               ttimeFileRoot='time'+s+'layer', 
               outputFileRoot='loc'+s+"loc",
               loc_file='loc.in', 
               locgrid=[], 
               grid2time_files=[], 
               vel2grid_file=None, 
               name='run_file.sh', 
               start_clean=True,
               all_files_in=True,
               run_file=True, 
               printf=True,
               print_prefix="\t",
               copy_bin=False, 
               control=[1, 54321], 
               trans=None, 
               gtmode="GRID3D ANGLES_YES", 
               gt_plfd=[1.0e-3, 0], 
               wavetype=["P","S"],
               step=None, 
               gridType='PROB_DENSITY', 
               saveFlag="SAVE", 
               zOrig=None, xNum=None, yNum=None, 
               zNum=None, xOrig=None, yOrig=None, 
               locsearch=["OCT", 10, 10, 4, 0.01, 100000, 5000, 0, 1], 
               locqual2err=[0.1, 0.5, 1, 2, 99999.9],
               obsFileType="NLLOC_OBS", 
               locsig=" NonLinLoc - ALomax Scientific",
               loccom=" ", 
               lochypout="SAVE_NLLOC_ALL",
               locmeth=["EDT_OT_WT", 9999.0, 4, -1, -1, -1.730, -1, 0.1], 
               locgau=[0.15, 1.0, ''],
               locgau2=None,
               locphstat=[15, -1, 9999.0, 1.0, 1.0],
               locangles=["ANGLES_NO", 5],
               locmag=["ML_HB", 1.0, 1.11, 0.00189],
               locphaseid=["P", "p", "G", "PN", "PG", "S", "s", "G", "SN", "SG"],
               locelevcorr=[],
               loctopo_surface=[],
               flagDumpDecimation=1 ):

    """
    Run NonLinLoc (NLLoc) earthquake location software.
    For more information on NLLoc parameters, see http://alomax.net/nlloc/

    Parameters:
        
        - path (str): The path where the NLLoc files will be created. 
            Default is the current working directory.

        - vel2grid (bool): Flag to indicate whether to perform Vel2Grid computation.
            Default is False.
        
        - grid2time (bool): Flag to indicate whether to perform Grid2Time computation. 
            Default is False.
        
        - nlloc (bool): Flag to indicate whether to perform earthquake location. 
            Default is False.
        
        - stFile (str): The path to the station file. Default is None.
        
        - obsFile (str or dict): The path to the observation file or a dictionary 
            containing travel-timeinformation. Default is None.
        
        - modelFileRoot (str): The root name of the model file. Default is 'model'+s+'layer'.
        
        - ttimeFileRoot (str): The root name of the travel time file. Default is 'time'+s+'layer'.
        
        - outputFileRoot (str): The root name of the output file. Default is 'loc'+s+'loc'.
        
        - loc_file (str): The path to the location file. Default is 'loc.in'.
        
        - grid2time_files (list): A list of paths to Grid2Time files. Default is an empty list.
        
        - vel2grid_file (str): The path to the Vel2Grid file. Default is None.
        
        - name (str): The name of the run file. Default is 'run_file.sh'.
        
        - start_clean (bool): Flag to indicate whether to start with a clean directory. 
            Default is True.
        
        - all_files_in (bool): Flag to indicate whether all files are in the same directory. 
            Default is True.
        
        - run_file (bool): Flag to indicate whether to create the run file. Default is True.
        
        - printf (bool): Flag to indicate whether to print the run file. Default is True.
        
        - print_prefix (str): The prefix for the print statement. Default is "\t".

        - copy_bin (bool): Flag to indicate whether to copy the NLLoc executables 
            into the input path. Default is True.
        
        - control (list): A list containing the control parameters for NLLoc. 
            Default is [1, 54321].
        
        - trans (str): The transformation type. Default is "NONE".
        
        - gtmode (str): The grid type and angle mode. Default is "GRID3D ANGLES_YES".
        
        - gt_plfd (list): A list containing the grid type and plane fit distance parameters. 
            Default is [1.0e-3, 0].
        
        - wavetype (list): A list of wave types to consider. Default is ["P", "S"].
        
        - step (float): The step size for the grid. Default is None.
        
        - gridType (str): The grid type. Default is "PROB_DENSITY".
        
        - saveFlag (str): The save flag. Default is "SAVE".
        
        - zOrig (float): The origin of the z-axis. Default is None.
        
        - xNum (int): The number of grid points in the x-direction. Default is None.
        
        - yNum (int): The number of grid points in the y-direction. Default is None.
        
        - zNum (int): The number of grid points in the z-direction. Default is None.
        
        - xOrig (float): The origin of the x-axis. Default is None.
        
        - yOrig (float): The origin of the y-axis. Default is None.
        
        - locsearch (list): A list containing the location search parameters. 
            Default is ["OCT", 10, 10, 4, 0.01, 100000, 5000, 0, 1].
        
        - locqual2err (list): A list containing the location quality to error parameters. 
            Default is [0.1, 0.5, 1, 2, 99999.9].
        
        - obsFileType (str): The observation file type. Default is "NLLOC_OBS".
        
        - locsig (str): The signature for the location. 
            Default is "NonLinLoc - ALomax Scientific".
        
        - loccom (str): The comment for the location. Default is an empty string.
        
        - lochypout (str): The output format for the hypocenter file. 
            Default is "SAVE_NLLOC_ALL".
        
        - locmeth (list): A list containing the location method parameters. 
            Default is ["EDT_OT_WT", 9999.0, 4, -1, -1, -1.730, -1, 0.1].
        
        - locgau (list): A list containing the Gaussian smoothing parameters 
            for the first iteration. Default is [0.15, 1.0, ''].
        
        - locgau2 (list): A list containing the Gaussian smoothing parameters 
            for the second iteration. Default is [0.05, 0.5, ''].
        
        - locphstat (list): A list containing the travel-timestatistics parameters. 
            Default is [15, -1, 9999.0, 1.0, 1.0].
        
        - locangles (list): A list containing the angle parameters. 
            Default is ["ANGLES_NO", 5].
        
        - locmag (list): A list containing the magnitude parameters. 
            Default is ["ML_HB", 1.0, 1.11, 0.00189].
        
        - locphaseid (list): A list containing the phase identification parameters for phases. 
            Default is ["P", "p", "G", "PN", "PG", "S", "s", "G", "SN", "SG"].
        
        - locelevcorr (list): A list containing the elevation correction parameters. 
            Default is an empty list.
        
        - loctopo_surface (list): A list containing the topographic surface parameters. 
            Default is an empty list.
        
        - flagDumpDecimation (int): The flag for dump decimation. Default is 1.

        - st_label (str): The station label. Default is "XYZ".
            This must be specified when the stFile is a dictionary
            that conteins station coordinates both as 'x', 'y' and 'lon', 'lat'.
            Possible st_labels are "XYZ" and "LATLON".

    Returns:
            
        - str: The path to the run file.
    """

    # --------------------------------------------------------------
    # Create the code directory
    if os.path.isdir( path ) is False :
        os.makedirs( path, exist_ok=True )

    # Create the run file f
    f = open( path +s+ name, 'w' )

    # If copy_bin is true, 
    # copy directory bin with all the NLLoc executables into the input path
    if copy_bin is True : 
        os.makedirs( path+s+'bin', exist_ok=True )
        for root, dirs, files in os.walk( mdir+s+'src'+s+'bin' ):
            for file in files:
                shutil.copyfile( mdir+s+'src'+s+'bin'+s+file, path+s+'bin'+s+file )
        bin = 'bin/'
    else :
        bin = ''

    # --------------------------------------------------------------
    # VEL2GRID 

    if vel2grid is True :

        # Create/Copy Vel2Grid file
        if vel2grid_file is not None :
            if "GRID3D" in gtmode :
                grid_dim = "3d"
            if "GRID2D" in gtmode :
                grid_dim = "2d"
            v2g_file_name = vel2grid_file.split( s )[-1]
            if os.path.isfile( path +s+ v2g_file_name ) :
                if shutil._samefile( path+s+v2g_file_name, vel2grid_file ) is False :
                    shutil.copy( vel2grid_file, path+s+v2g_file_name )

            # Write Vel2Grid coommand to run file
            l = "# Create input model grid"
            f.write( l + ' \n' )
            if ( grid_dim == "3d" ) or ( grid_dim == "3D" ) :
                f.write( bin+'Vel2Grid3D ' + v2g_file_name + ' \n' )  
            else :
                f.write( bin+'Vel2Grid ' + v2g_file_name + ' \n' )  

    # --------------------------------------------------------------
    # GRID2TIME

    # Create/Copy grid2time file 
    if grid2time is True :

        grid2time_file_names = []

        if grid2time_files != [] :
            for i in grid2time_files :
                g2t_file_name = i.split( s )[-1]
                if os.path.isfile( path +s+ g2t_file_name ) :
                    if shutil._samefile( path+s+g2t_file_name, i ) is False :
                        shutil.copy( i, path+s+g2t_file_name )  
                grid2time_file_names.append( g2t_file_name )

        else :
            for i, w in enumerate( wavetype ) :
                if w == "P" :
                    g2tfP = create_grd2time_file( stFile=stFile, path=path, name="grid2timeP.in", 
                                                  trans=trans, wavetype=w, gtmode=gtmode,
                                                  ttimeFileRoot=modelFileRoot, control=control,
                                                  outputFileRoot=ttimeFileRoot, gt_plfd=gt_plfd,
                                                  all_files_in=all_files_in,
                                                  start_clean=start_clean, printf=printf )
                if w == "S" :
                    g2tfS = create_grd2time_file( stFile=stFile, path=path, name="grid2timeS.in", 
                                                  trans=trans, wavetype=w, gtmode=gtmode,
                                                  ttimeFileRoot=modelFileRoot, control=control,
                                                  outputFileRoot=ttimeFileRoot, gt_plfd=gt_plfd,
                                                  all_files_in=all_files_in,
                                                  start_clean=start_clean, printf=printf )
            grid2time_file_names = [ g2tfP, g2tfS ]

        # Write Grid2Time coommand to run file
        for i, g2t in enumerate( grid2time_file_names ) :
            print( g2t )
            l = "# compute travel time from each stations to each grid point"
            f.write( l + ' \n' )
            f.write( bin+'Grid2Time ' + g2t + ' \n' )

    # --------------------------------------------------------------
    # NLLOC
    
    # Create/Copy loc files
    if ( nlloc is True ) and ( loc_file is not None ) :

        # Copy/Create Obs. file
        if obsFile is not None :

            if type( obsFile ) == str :
                obs_file_name = obsFile.split( s )[-1]
                if os.path.isfile( path +s+ obs_file_name ) :
                    if shutil._samefile( path+s+obs_file_name, obsFile ) is False :
                        shutil.copy( obsFile, path+s+obs_file_name )
                    else :
                        shutil.copy( obsFile, path+s+obs_file_name )
                    obsFile = path+s+obs_file_name

            if type( obsFile ) == dict :
                obs_file_name = "obs.in"
                obsFile = create_tt_file( path +s+ obs_file_name, tt_dict=obsFile )

        # Copy loc file 
        loc_file_name = loc_file.split( s )[-1]

        # Get the directory and name of the model file
        modelFileDir = os.path.dirname( modelFileRoot )
        ttimeFileDir = os.path.dirname( ttimeFileRoot ) 
        
        if os.path.isfile( path+s+loc_file_name ) \
           and not shutil._samefile( path+s+loc_file_name, loc_file ) \
           and os.path.isfile( loc_file ) :
            
            shutil.copy( loc_file, path+s+loc_file_name )

        # Create locgrid list
        else :

            # Initialize modgrid list empty
            modgrid = []
            # Initialize timgrid list empty
            timgrid = []
            # Initialize locgrid list empty
            locgrid = [ 0, 0, 0, 0, 0, 0, 1, 1, 1, gridType, saveFlag ]

            # Check if the model directory is an absolute path
            if os.path.isabs( modelFileDir ) and ( modgrid == [] ) :
                # Check if the directory exists
                if os.path.isdir( modelFileDir ) :
                    # Get the locgrid list from the model file
                    model_lst = read_model_file( modelFileDir )
                    if model_lst : 
                        modgrid = create_locgrid( model = model_lst[0] )

            # Check if the model directory is a relative path
            if ( os.path.isdir( path +s+ modelFileDir ) ) and ( modgrid == [] ) :
                # Get the locgrid list from the model file
                model_lst = read_model_file( path +s+ modelFileDir )
                if model_lst : 
                    modgrid = create_locgrid( model = model_lst[0] )

            # Check if the ttime directory is an absolute path
            if ( os.path.isabs( ttimeFileDir ) ) and ( timgrid == [] ) :
                # Check if the directory exists
                if os.path.isdir( ttimeFileDir ) :
                    # Get the locgrid list from the model file
                    tt_model_lst = read_model_file( ttimeFileDir )
                    if tt_model_lst : 
                        timgrid = create_locgrid( model = tt_model_lst[0] )

            # Check if the ttime directory is a relative path
            if ( os.path.isdir( path +s+ ttimeFileDir ) ) and ( timgrid == [] ) :
                # Get the locgrid list from the model file
                tt_model_lst = read_model_file( path +s+ ttimeFileDir )
                if tt_model_lst : 
                    timgrid = create_locgrid( model = tt_model_lst[0] )

            if timgrid == [] :
                if modgrid != [] :
                    timgrid = utl.copy.deepcopy( modgrid )
                else :
                    timgrid = utl.copy.deepcopy( locgrid )
    
            if timgrid != locgrid :
                locgrid = utl.copy.deepcopy( timgrid )

            if xNum is not None :
                locgrid[0] = xNum
            if yNum is not None :
                locgrid[1] = yNum
            if zNum is not None :
                locgrid[2] = zNum
            if xOrig is not None :
                locgrid[3] = xOrig
            if yOrig is not None :
                locgrid[4] = yOrig
            if zOrig is not None :
                locgrid[5] = zOrig
            if step is not None :
                locgrid[6] = step
                locgrid[7] = step
                locgrid[8] = step

            if ( timgrid != [] ) and\
                 ( locgrid[3] + locgrid[0] * locgrid[6] >\
                   timgrid[3] + timgrid[0] * timgrid[6] ) :
                 locgrid[0] =  int( np.floor( ( ( timgrid[3] + timgrid[0] * timgrid[6] ) -\
                            locgrid[3]  ) / locgrid[6] ) )
            if ( timgrid != [] ) and\
                    ( locgrid[4] + locgrid[1] * locgrid[7] >\
                        timgrid[4] + timgrid[1] * timgrid[7] ) :
                locgrid[1] =  int( np.floor( ( ( timgrid[4] + timgrid[1] * timgrid[7] ) -\
                            locgrid[4] ) / locgrid[7] ) )
            if ( timgrid != [] ) and\
                    ( locgrid[5] + locgrid[2] * locgrid[8] >\
                        timgrid[5] + timgrid[2] * timgrid[8] ) :
                locgrid[2] =  int( np.floor( ( ( timgrid[5] + timgrid[2] * timgrid[8] ) -\
                            locgrid[5] ) / locgrid[8] ) )
            
            # Create loc file
            loc_file = create_loc_file( path=path,
                                        obsFile=obsFile, 
                                        locgrid=locgrid, 
                                        name=loc_file_name,
                                        ttimeFileRoot=ttimeFileRoot,
                                        outputFileRoot=outputFileRoot,
                                        trans=trans, 
                                        locsearch=locsearch, 
                                        locqual2err=locqual2err,
                                        obsFileType=obsFileType, 
                                        locsig=locsig, 
                                        loccom=loccom,
                                        lochypout=lochypout, 
                                        locmeth=locmeth, 
                                        locgau=locgau,
                                        locgau2=locgau2, 
                                        locphstat=locphstat, 
                                        locangles=locangles,
                                        locmag=locmag, 
                                        locphaseid=locphaseid, 
                                        locelevcorr=locelevcorr, 
                                        loctopo_surface=loctopo_surface,
                                        flagDumpDecimation=flagDumpDecimation,
                                        all_files_in=all_files_in,
                                        start_clean=start_clean, printf=printf )

        # Write NLLoc coommand to run file
        l = "# locate earthquakes (s)"
        f.write( l + ' \n' )
        f.write( bin+'NLLoc ' + loc_file_name + ' \n' )
    
    # --------------------------------------------------------------
    # Close the run file
    f.close()

    # print run file
    if printf == True :
        print( f"\n{print_prefix}# --------------------------------------------------------------" )
        print( f"{print_prefix}Run file :\n{print_prefix}{path +s+ name}" )
        file_lines_list = utl.read_file(  path +s+ name )
        print(f"{print_prefix}CONTENT:\n")
        for lr in file_lines_list :
            print( f"{print_prefix}{lr}" )
        print( f"\n{print_prefix}# --------------------------------------------------------------\n" )

    # Make run file executable 
    os.system( 'chmod u+x ' + path +s+ name )

    # Run file call
    if run_file == True : 
        run_cmd( string2run=name, path=path, addbefore='./' )

    if printf == True :
        print( f"\n{print_prefix}# --------------------------------------------------------------" )

    return path +s+ name
    
# -----------------------------------------------------------------------------    
def read_loc( path, plot=False, aspect='auto', reverse_z=True, size=None ) :
    """
    Read earthquake location data from files in a given directory (usualy called loc).
    
    Args:
        path (str): The path to the directory containing the earthquake location files.
        plot (bool, optional): Whether to plot the earthquake locations. Defaults to False.
        aspect (str, optional): The aspect ratio of the plot. Defaults to 'auto'.
        reverse_z (bool, optional): Whether to reverse the sign of the z-coordinate. Defaults to True.
        size (None, optional): Not used in the function.

    Returns:
        
        tuple: A tuple containing two dictionaries: eve_dict and phs_dict.
            - eve_dict (dict): A dictionary containing earthquake event information.
            - phs_dict (dict): A dictionary containing phase arrival information.

    """

    file_name = [] 

    for root, dirs, files in os.walk( path ):
        for file in sorted(files):
            if ( file.endswith('.hyp') ) and ( '.sum.' not in file ) and ( 'last.' not in file ):
                file_name.append( file )

    eve_dict = { 'eve':[], 'x':[], 'y':[], 'z':[], 'yy':[], 'mm':[], 'dd':[], 'h':[], 
                 'm':[], 's':[], 'pmax':[], 'mfmin':[], 'mfmax':[], 'rms':[], 'nphs':[],
                 'gap':[], 'dist':[], 'mamp':[], 'mdur':[], 'lon':[], 'lat':[],
                 'ellaz1':[], 'dip1':[], 'len1':[], 'ellaz2':[], 'dip2':[],
                 'len2':[], 'len3':[], 'horunc':[], 'covxx':[], 'covyy':[], 'covzz':[], 
                 'errx':[], 'erry':[], 'errz':[] }

    phs_dict = { 'eve':[], 'st':[],  'instrument':[], 'component':[], 'p_phase_onset':[],
                 'phase_descriptor':[], 'first_motion':[], 'yy':[], 'mm':[], 'dd':[], 
                 'h':[], 'm':[], 's':[], 'err':[], 'errmag':[], 'coda_duration':[],
                 'amplitude':[], 'period':[], 'ttpred':[], 'res':[], 'lon':[], 'lat':[],
                 'weight':[], 'statx':[], 'staty':[], 'statz':[], 'sdist':[], 'sazim':[], 
                 'raz':[], 'rdip':[], 'rqual':[], 'tcorr':[], 'p0s1':[], 'x':[], 'y':[], 
                 'z':[], 'yy0':[], 'mm0':[], 'dd0':[], 'h0':[], 'm0':[], 's0':[], 
                 'errx':[], 'erry':[], 'errz':[], 'gap':[], 'tt':[] }

    n1 = 0

    for f in file_name :

        fi = open( path +s+ f, 'r' )
        linesi = fi.readlines()
        fi.close()
        n1 += 1

        for i, l in enumerate( linesi ) :
            if "HYPOCENTER" in l :
                
                eve_dict['eve'].append( n1 )
                lsplit = l.split() 
                for si, ls in enumerate( lsplit ) : 
                    if ls == 'x' :
                        eve_dict['x'].append( float( lsplit[ si + 1 ] ) )   
                    if ls == 'y' :
                        eve_dict['y'].append( float( lsplit[ si + 1 ] ) ) 
                    if ls == 'z' :
                        eve_dict['z'].append( float( lsplit[ si + 1 ] ) ) 
                        
            if "GEOGRAPHIC" in l :      
                lsplit = l.split() 
                for si, ls in enumerate( lsplit ) : 
                    if ls == 'Long' :
                        eve_dict['lon'].append( float( lsplit[ si + 1 ] ) )   
                    if ls == 'Lat' :
                        eve_dict['lat'].append( float( lsplit[ si + 1 ] ) )     
                    if ls == 'OT' :
                        eve_dict['yy'].append( int( lsplit[ si + 1 ] ) )  
                        eve_dict['mm'].append( int( lsplit[ si + 2 ] ) )  
                        eve_dict['dd'].append( int( lsplit[ si + 3 ] ) )
                        eve_dict['h'].append( int( lsplit[ si + 4 ] ) )
                        eve_dict['m'].append( int( lsplit[ si + 5 ] ) )
                        eve_dict['s'].append( float( lsplit[ si + 6 ] ) )
                        
            if "QUALITY" in l :      
                lsplit = l.split() 
                for si, ls in enumerate( lsplit ) : 
                    if ls == 'Pmax' :
                        eve_dict['pmax'].append( float( lsplit[ si + 1 ] ) )   
                    if ls == 'MFmin' :
                        eve_dict['mfmin'].append( float( lsplit[ si + 1 ] ) )     
                    if ls == 'MFmax' :
                        eve_dict['mfmax'].append( float( lsplit[ si + 1 ] ) )  
                    if ls == 'RMS' :
                        eve_dict['rms'].append( float( lsplit[ si + 1 ] ) )  
                    if ls == 'Nphs' :
                        eve_dict['nphs'].append( float( lsplit[ si + 1 ] ) ) 
                    if ls == 'Gap' :
                        eve_dict['gap'].append( float( lsplit[ si + 1 ] ) )
                    if ls == 'Dist' :
                        eve_dict['dist'].append( float( lsplit[ si + 1 ] ) )
                    if ls == 'Mamp' :
                        eve_dict['mamp'].append( float( lsplit[ si + 1 ] ) )
                    if ls == 'Mdur' :
                        eve_dict['mdur'].append( float( lsplit[ si + 1 ] ) )

            if "STATISTICS" in l :      
                lsplit = l.split() 
                for si, ls in enumerate( lsplit ) :
                    if ls == 'CovXX' :
                        eve_dict['covxx'].append( float( lsplit[ si + 1 ] ) ) 
                        eve_dict['errx'].append( np.sqrt( 3.53 * eve_dict['covxx'][-1] ) )
                    if ls == 'YY' :
                        eve_dict['covyy'].append( float( lsplit[ si + 1 ] ) )
                        eve_dict['erry'].append( np.sqrt( 3.53 * eve_dict['covyy'][-1] ) )
                    if ls == 'ZZ' :
                        eve_dict['covzz'].append( float( lsplit[ si + 1 ] ) ) 
                        eve_dict['errz'].append( np.sqrt( 3.53 * eve_dict['covzz'][-1] ) )
                    if ls == 'EllAz1' :
                        eve_dict['ellaz1'].append( float( lsplit[ si + 1 ] ) )   
                    if ls == 'Dip1' :
                        eve_dict['dip1'].append( float( lsplit[ si + 1 ] ) )     
                    if ls == 'Len1' :
                        eve_dict['len1'].append( float( lsplit[ si + 1 ] ) )  
                    if ls == 'Az2' :
                        eve_dict['ellaz2'].append( float( lsplit[ si + 1 ] ) )   
                    if ls == 'Dip2' :
                        eve_dict['dip2'].append( float( lsplit[ si + 1 ] ) )     
                    if ls == 'Len2' :
                        eve_dict['len2'].append( float( lsplit[ si + 1 ] ) ) 
                    if ls == 'Len3' :
                        eve_dict['len3'].append( float( lsplit[ si + 1 ] ) ) 

            if "QML_OriginUncertainty" in l :      
                lsplit = l.split() 
                for si, ls in enumerate( lsplit ) : 
                    if ls == 'horUnc' :
                        eve_dict['horunc'].append( float( lsplit[ si + 1 ] ) )

            if "PHASE ID Ins Cmp"  in l :
                lines_phs = linesi[ i+1: ]
                for iphs, lphs in enumerate( linesi[ i+1: ] ) :
                    if "END_PHASE" in lphs :
                        break
                    else :
                        phs_list = lphs.split()
                        phs_dict['eve'].append( eve_dict['eve'][-1] )
                        phs_dict['st'].append( phs_list[0] )
                        phs_dict['instrument'].append( phs_list[1] )
                        phs_dict['component'].append( phs_list[2] )
                        phs_dict['p_phase_onset'].append( phs_list[3] )
                        phs_dict['phase_descriptor'].append( phs_list[4] )
                        phs_dict['first_motion'].append( phs_list[5] )
                        phs_dict['yy'].append( int( phs_list[6][:4] ) )
                        phs_dict['mm'].append( int( phs_list[6][4:6] ) )
                        phs_dict['dd'].append( int( phs_list[6][6:] ) )
                        phs_dict['h'].append( int( phs_list[7][:2] ) )
                        phs_dict['m'].append( int( phs_list[7][2:] ) )
                        phs_dict['s'].append( float( phs_list[8] ) )
                        phs_dict['err'].append( phs_list[9] )
                        phs_dict['errmag'].append( float( phs_list[10] ) )
                        phs_dict['coda_duration'].append( float( phs_list[11] ) )
                        phs_dict['amplitude'].append( float( phs_list[12] ) )
                        phs_dict['period'].append( float( phs_list[13] ) )
                        phs_dict['ttpred'].append( float( phs_list[15] ) )
                        phs_dict['res'].append( float( phs_list[16] ) )
                        phs_dict['weight'].append( float( phs_list[17] ) )
                        phs_dict['statx'].append( float( phs_list[18] ) )
                        phs_dict['staty'].append( float( phs_list[19] ) )
                        phs_dict['statz'].append( float( phs_list[20] ) )
                        phs_dict['sdist'].append( float( phs_list[21] ) )
                        phs_dict['sazim'].append( float( phs_list[22] ) )
                        phs_dict['raz'].append( float( phs_list[23] ) )
                        phs_dict['rdip'].append( float( phs_list[24] ) )
                        phs_dict['rqual'].append( float( phs_list[25] ) )
                        phs_dict['tcorr'].append( float( phs_list[26] ) )
                        if 'S' in phs_list[4] :
                            phs_dict['p0s1'].append( 1 )
                        if 'P' in phs_list[4] :
                            phs_dict['p0s1'].append( 0 )
                        phs_dict['x'].append( float( eve_dict['x'][-1] ) )
                        phs_dict['y'].append( float( eve_dict['y'][-1] ) )
                        phs_dict['z'].append( float( eve_dict['z'][-1] ) )
                        phs_dict['yy0'].append( int( eve_dict['yy'][-1] ) )
                        phs_dict['mm0'].append( int( eve_dict['mm'][-1] ) )
                        phs_dict['dd0'].append( int( eve_dict['dd'][-1] ) )
                        phs_dict['h0'].append( int( eve_dict['h'][-1]) )
                        phs_dict['m0'].append( int( eve_dict['m'][-1] ) )
                        phs_dict['s0'].append( float( eve_dict['s'][-1] ) )
                        phs_dict['errx'].append( float( eve_dict['errx'][-1] ) )
                        phs_dict['erry'].append( float( eve_dict['erry'][-1] ) )
                        phs_dict['errz'].append( float( eve_dict['errz'][-1] ) ) 
                        phs_dict['gap'].append( float( eve_dict['gap'][-1] ) ) 
                        phs_dict['lon'].append( float( eve_dict['lon'][-1] ) ) 
                        phs_dict['lat'].append( float( eve_dict['lat'][-1] ) )
    
    for k in eve_dict :
        eve_dict[ k ] = np.asarray( eve_dict[ k ] )

    for k in phs_dict :
        phs_dict[ k ] = np.asarray( phs_dict[ k ] )
    
    if reverse_z == True :
        eve_dict['z'] = eve_dict['z'] * -1
        phs_dict['z'] = phs_dict['z'] * -1

    phs_dict['tp'] = np.full( phs_dict['eve'].size, np.nan )
    phs_dict['ts-tp'] = np.full( phs_dict['eve'].size, np.nan )
    phs_dict['tti'] = utl.combine64( years=phs_dict['yy'], months=phs_dict['mm'], 
                                    days=phs_dict['dd'], hours=phs_dict['h'], 
                                    minutes=phs_dict['m'], seconds=phs_dict['s'] )
    phs_dict['tt0'] = utl.combine64( years=phs_dict['yy0'], months=phs_dict['mm0'], 
                                    days=phs_dict['dd0'], hours=phs_dict['h0'], 
                                    minutes=phs_dict['m0'], seconds=phs_dict['s0'] )
    phs_dict['tt'] = ( phs_dict['tti'].astype(float) - phs_dict['tt0'].astype(float) )/1e9   

    if plot==True :
        
        utl.plt.subplot(131, aspect=aspect )
        utl.plt.scatter( eve_dict['x'], eve_dict['z'], s=eve_dict['errx'], c=eve_dict['errz'] ) 
        utl.plt.title( "Z/X [km]")
        
        utl.plt.subplot(132, aspect=aspect)
        utl.plt.scatter(  eve_dict['y'], eve_dict['z'], s=eve_dict['errx'], c=eve_dict['errz'] )
        utl.plt.title( "Z/Y [km]")
        
        utl.plt.subplot(133, aspect=aspect)
        utl.plt.scatter(  eve_dict['x'], eve_dict['y'], s=eve_dict['errx'], c=eve_dict['errz'] )
        utl.plt.colorbar()
        utl.plt.title( "Y/X [km]")
        utl.plt.tight_layout()
        
    return eve_dict, phs_dict

# -----------------------------------------------------------------------------
def plot_eve_st( xyz_events=None, 
                 xyz_stations=None, 
                 esize=None, ssize=None, 
                 ecolor='b', scolor='r', 
                  emarker='o', smarker='v', alpha=0.5,
                 aspect='auto', zst_scale=1, drop_duplicates=True, 
                 idxe=None, idxs=None, 
                 titles=['X / Z', 'Y / Z', 'Y / X'] ) :
    """
    Plot events and stations in 3D space.

    Parameters:
    - xyz_events (array-like or dict or None): The coordinates of the events in 3D space. 
      If array-like, it should have shape (N, 3) where N is the number of events and each row represents the (x, y, z) coordinates of an event.
      If dict, it should have keys 'x', 'y', and 'z' with corresponding values being arrays of the x, y, and z coordinates of the events.
      If None, no events will be plotted.
    - xyz_stations (array-like or dict or None): The coordinates of the stations in 3D space. 
      If array-like, it should have shape (M, 3) where M is the number of stations and each row represents the (x, y, z) coordinates of a station.
      If dict, it should have keys 'x', 'y', and 'z' with corresponding values being arrays of the x, y, and z coordinates of the stations.
      If None, no stations will be plotted.
    - size (float or None): Deprecated. Use esize and ssize instead.
    - esize (float or None): The size of the event markers.
    - ssize (float or None): The size of the station markers.
    - ecolor (str): The color of the event markers.
    - scolor (str): The color of the station markers.
    - emarker (str): The marker style for the event markers.
    - smarker (str): The marker style for the station markers.
    - alpha (float): The transparency of the markers.
    - aspect (str): The aspect ratio of the subplots.
    - zst_scale (float): The scaling factor for the z-coordinate of the stations.
    - drop_duplicates (bool): Whether to drop duplicate events and stations.
    - idxe (array-like or None): The indices of the events to plot.
    - idxs (array-like or None): The indices of the stations to plot.
    - titles (list): The titles of the subplots.

    Returns:
    - tuple: A tuple containing three subplots.
    """

    # Create three subplots with the same aspect ratio
    p1 = lsz_plot.plt.subplot( 131, aspect=aspect )
    p2 = lsz_plot.plt.subplot( 132, aspect=aspect )
    p3 = lsz_plot.plt.subplot( 133, aspect=aspect )

    # If events data is provided
    if xyz_events is not None :

        # If the data is a list or tuple, convert each element to a numpy array
        if type( xyz_events ) in (list, tuple) :
            xeve = np.asarray( xyz_events[0] )
            yeve = np.asarray( xyz_events[1] )
            zeve = np.asarray( xyz_events[2] )
        # If the data is a dictionary, convert each value to a numpy array
        if type( xyz_events ) == dict :
            xeve = np.asarray( xyz_events['x'] )
            yeve = np.asarray( xyz_events['y'] )
            zeve = np.asarray( xyz_events['z'] )
        # If the data is a numpy array, split it into x, y, and z coordinates
        if type( xyz_events ) == np.ndarray :
            xeve = xyz_events[:,0]
            yeve = xyz_events[:,1]
            zeve = xyz_events[:,2]

        # If idxe is provided, select the events with the given indices
        if idxe is not None :
            xeve = xeve[ idxe ]
            yeve = yeve[ idxe ]
            zeve = zeve[ idxe ]

        # If duplicates should be dropped, find the unique rows in the data
        if drop_duplicates == True :
            xyz_eu = np.unique( np.column_stack( ( xeve, yeve, zeve) ), axis=0 )
            xeve = xyz_eu[:,0]
            yeve = xyz_eu[:,1]
            zeve = xyz_eu[:,2] 

        # Plot the events on each subplot
        p1.scatter( xeve,  zeve, marker=emarker, c=ecolor, s=esize, label='events' )  
        p2.scatter( yeve,  zeve, marker=emarker, c=ecolor, s=esize, label='events' )
        p3.scatter( xeve,  yeve, marker=emarker, c=ecolor, s=esize, label='events' )  
    
    # If stations data is provided
    if xyz_stations is not None : 

        # If the data is a list or tuple, convert each element to a numpy array
        if type( xyz_stations ) in (list, tuple) :
            xst = np.asarray( xyz_stations[0] )
            yst = np.asarray( xyz_stations[1] ) 
            zst = np.asarray( xyz_stations[2] ) * zst_scale 
        # If the data is a dictionary, convert each value to a numpy array
        if type( xyz_stations ) == dict :
            xst = np.asarray( xyz_stations['x'] )
            yst = np.asarray( xyz_stations['y'] )
            zst = np.asarray( xyz_stations['z'] ) * zst_scale
        # If the data is a numpy array, split it into x, y, and z coordinates
        if type( xyz_stations ) == np.ndarray :
            xst = xyz_stations[:,0]
            yst = xyz_stations[:,1]
            zst = xyz_stations[:,2] * zst_scale

        # If idxs is provided, select the stations with the given indices
        if idxs is not None :
            xst = xst[ idxs ]
            yst = yst[ idxs ]
            zst = zst[ idxs ]

        # If duplicates should be dropped, find the unique rows in the data
        if drop_duplicates == True :
            xyz_su = np.unique( np.column_stack( ( xst, yst, zst) ), axis=0 )
            xst = xyz_su[:,0]
            yst = xyz_su[:,1]
            zst = xyz_su[:,2] 

        # Plot the stations on each subplot
        p1.scatter( xst,  zst, marker=smarker, c=scolor, s=ssize, label='stations', alpha=alpha )  
        p2.scatter( yst,  zst, marker=smarker, c=scolor, s=ssize, label='stations', alpha=alpha )
        p3.scatter( xst,  yst, marker=smarker, c=scolor, s=ssize, label='stations', alpha=alpha )      

    # Set the title for each subplot
    p1.set_title( titles[0] )
    p2.set_title( titles[1] )
    p3.set_title( titles[2] )

    lsz_plot.plt.tight_layout()

    return p1, p2, p3

# -----------------------------------------------------------------------------
def nlloc2simul( tt_data, 
                 reverse_z=True, 
                 path=os.getcwd(), 
                 file_name='tt_data.sim',
                 locqual2err=[0.01, 0.02, 0.08, 0.2, 99999.9], 
                 nzco=5, 
                 prj_cart='+proj=gnom', 
                 prj_geo=None, 
                 lon0=None, 
                 lat0=None,
                 change_lon0_sign=False, 
                 rota=0, 
                 prj2geo=True ) :

    """
    Converts NLLoc travel-time data to SIMUL format and writes it to a file.

    Parameters:
        - tt_data (str or DataFrame): Travel-time data or path to the file containing travel-time data.
        - reverse_z (bool, optional): If True, reverses the sign of the z-coordinate. Default is True.
        - path (str, optional): Directory path where the output file will be saved. Default is the current working directory.
        - file_name (str, optional): Name of the output file. Default is 'tt_data.sim'.
        - locqual2err (list, optional): List of error magnitude thresholds for location quality. Default is [0.01, 0.02, 0.08, 0.2, 99999.9].
        - nzco (int, optional): Number of zones for coordinate transformation. Default is 5.
        - prj_cart (str, optional): Cartographic projection string. Default is '+proj=gnom'.
        - prj_geo (str, optional): Geographic projection string. Default is None.
        - lon0 (float, optional): Longitude of the projection center. Default is None.
        - lat0 (float, optional): Latitude of the projection center. Default is None.
        - change_lon0_sign (bool, optional): If True, changes the sign of lon0. Default is False.
        - rota (float, optional): Rotation angle for the projection. Default is 0.
        - prj2geo (bool, optional): If True, converts projected coordinates to geographic coordinates. Default is True.

    Returns:
        - str: Path to the output file.
    """

    os.makedirs( path, exist_ok=True )

    file = path +s+ file_name

    if type( tt_data ) == str :
        tt_data = read_loc( path=tt_data, reverse_z=reverse_z )[1]
    else :
        tt_data = tt_data.copy()
        if reverse_z == True :
            tt_data['z'] = tt_data['z'] * -1 

    print( tt_data.keys() )

    if (prj2geo == True) and ('x' in tt_data) and ('y' in tt_data) :
        tt_data['lon'], tt_data['lat'] = cart2geo( tt_data['x'], tt_data['y'], 
                                                    prj_cart=prj_cart, prj_geo=prj_geo,
                                                    nzco=nzco, lon0=lon0, lat0=lat0,
                                                    change_lon0_sign=change_lon0_sign, rota=rota )

    tt_data[ 'lon_deg' ] = np.sign( tt_data[ 'lon' ] ) * np.divmod( np.abs( tt_data[ 'lon' ] ), 
                            np.ones(tt_data[ 'lon' ].shape) )[0]
    tt_data[ 'lon_deg' ] = tt_data[ 'lon_deg' ].astype(int)
    tt_data[ 'lon_min' ] = np.divmod( np.abs( tt_data[ 'lon' ] ), 
                                     np.ones(tt_data[ 'lon' ].shape) )[1] * 60  
    tt_data[ 'lat_deg' ] = np.sign( tt_data[ 'lat' ] ) * np.divmod( np.abs( tt_data[ 'lat' ] ), 
                                     np.ones(tt_data[ 'lat' ].shape) )[0]
    tt_data[ 'lat_deg' ] = tt_data[ 'lat_deg' ].astype(int)
    tt_data[ 'lat_min' ] = np.divmod( np.abs( tt_data[ 'lat' ] ), 
                                     np.ones(tt_data[ 'lat' ].shape) )[1] * 60   

    tt_data[ 'w' ] = np.zeros( tt_data[ 'errmag' ].shape ).astype( int )
    idx0 = tt_data[ 'errmag' ] <= locqual2err[0]
    idx1 = ( tt_data[ 'errmag' ] > locqual2err[0] ) & ( tt_data[ 'errmag' ] <= locqual2err[1] )
    idx2 = ( tt_data[ 'errmag' ] > locqual2err[1] ) & ( tt_data[ 'errmag' ] <= locqual2err[2] )
    idx3 = ( tt_data[ 'errmag' ] > locqual2err[2] ) & ( tt_data[ 'errmag' ] <= locqual2err[3] )
    idx4 = ( tt_data[ 'errmag' ] > locqual2err[3] )
    tt_data[ 'w' ][idx0] = 0
    tt_data[ 'w' ][idx1] = 1
    tt_data[ 'w' ][idx2] = 2
    tt_data[ 'w' ][idx3] = 3
    tt_data[ 'w' ][idx4] = 4   

    tt_data['dt'] = utl.combine64( years=tt_data['yy'], months=tt_data['mm'], days=tt_data['dd'],
                                hours=tt_data['h'], minutes=tt_data['m'], seconds=tt_data['s'] )

    tt_data['dt0'] = utl.combine64( years=tt_data['yy0'], months=tt_data['mm0'], days=tt_data['dd0'],
                                hours=tt_data['h0'], minutes=tt_data['m0'], seconds=tt_data['s0'] )

    tt_data['Ddt'] = ( tt_data['dt'] - tt_data['dt0'] ) / np.timedelta64(1, 's')    

    idx = tt_data['Ddt'] < 0 
    for k in tt_data :
        tt_data[k] = tt_data[k][~idx] 

    eveu = np.unique( tt_data['eve'] )

    if 'mag' not in tt_data :
        tt_data['mag'] = np.full( tt_data['eve'].shape, 3 )

    
    simul_dict = { 'eve':[], 'yy':[], 'mm':[], 'dd':[], 'h':[], 'm':[], 's':[], 'lat':[], 
                   'lon':[], 'z':[], 'mag':[], 'st':[], 'r':[], 'w':[], 'tt':[],
                   'lon_deg':[], 'lon_min':[], 'lat_deg':[], 'lat_min':[] }
    
    simul_dict['eve'] = tt_data['eve']
    simul_dict['yy'] = tt_data['yy0']
    simul_dict['mm'] = tt_data['mm0']
    simul_dict['dd'] = tt_data['dd0']
    simul_dict['h'] = tt_data['h0']
    simul_dict['m'] = tt_data['m0']
    simul_dict['s'] = tt_data['s0']
    simul_dict['lat'] = tt_data['lat']
    simul_dict['lon'] = tt_data['lon']
    simul_dict['z'] = tt_data['z']
    simul_dict['mag'] = tt_data['mag']
    simul_dict['st'] = tt_data['st']
    simul_dict['r'] = tt_data['p_phase_onset']
    simul_dict['w'] = tt_data['w']
    simul_dict['r'] = np.array( [ r.replace('?', 'e') for r in simul_dict['r'] ] )
    simul_dict['tt'] = tt_data['Ddt']

    idx = tt_data['p0s1'] == 1
    for i in np.where( idx )[0] :
        idx2 = ( tt_data['st'] == tt_data['st'][i] ) & ( tt_data['p0s1'] == 0 )
        if np.any( idx2 ) :
            simul_dict['tt'][i] = tt_data['Ddt'][i] - tt_data['Ddt'][idx2][0]

    idx = simul_dict['yy0'] < 2000
    simul_dict['yy0'][idx] = simul_dict['yy0'][idx] - 1900
    simul_dict['yy0'][~idx] = simul_dict['yy0'][~idx] - 2000

    with open( file, 'w' ) as f : 

        for e in eveu :

            i = tt_data['eve'] == e

            if tt_data['lat_deg'][i][0] > 0 :
                lat_dir = 'N'
            else :
                lat_dir = 'S'

            if tt_data['lon_deg'][i][0] > 0 :
                lon_dir = 'E'
            else :
                lon_dir = 'W'

            line_1 = f"{str(tt_data['yy0'][i][0])[-2:]: >2s}{tt_data['mm0'][i][0]: >2d}" + \
                     f"{tt_data['dd0'][i][0]: >2d}{tt_data['h0'][i][0]: >3d}" + \
                     f"{tt_data['m0'][i][0]: >2d} {tt_data['s0'][i][0]: >5.2f}" + \
                     f"{tt_data['lat_deg'][i][0]: >3d}{lat_dir}{tt_data['lat_min'][i][0]: >5.2f} " +\
                     f"{tt_data['lon_deg'][i][0]: >3d}{lon_dir}{tt_data['lon_min'][i][0]: >5.2f}" +\
                     f"{tt_data['z'][i][0]: >7.2f}"

            if 'mag' in  tt_data :

                line_1 += f"{tt_data['mag'][i][0]: >7.2f}"

            line_1 += " \n"

            f.write( line_1)

            line_i = ''
            for ei,_ in enumerate( tt_data['eve'][i] ) :

                if tt_data['p0s1'][i][ei] == 1 :
                    idx = ( tt_data['st'][i] == tt_data['st'][i][ei] ) &\
                          ( tt_data['p0s1'][i] == 0 )
                    if np.any( idx ) :
                        sim_tt = tt_data['Ddt'][i][ei] - tt_data['Ddt'][i][idx][0]
                else :
                    sim_tt = tt_data['Ddt'][i][ei]

                n = ei + 1

                if ( n % 6 != 0 ) :

                    end = ''

                if ( n % 6 == 0 ) or ( n == len( tt_data['eve'][i] ) ) :
                    
                    end = ' \n'

                nm = tt_data['st'][i][ei]
                while len( nm ) < 4 :
                    nm = nm + ' '

                pha = tt_data['p_phase_onset'][i][ei].replace('?', 'e') +\
                      tt_data['phase_descriptor'][i][ei] +\
                      tt_data['first_motion'][i][ei].replace('?', 'p')

                line_i = line_i + f"{nm}{pha}" +\
                                  f"{tt_data['w'][i][ei]: >1d}{sim_tt: >6.2f}"+end
                # line_i = line_i + f"{tt_dict['name'][i][ei]: 4}{tt_dict['r'][i][ei]: >3}" +\
                #                   f"{tt_dict['w'][i][ei]: >1d}{tt_dict['tt'][i][ei]: >6.2f}"+end

            f.write( line_i )
            f.write( ' \n' )

    f.close()

    return simul_dict, file

# -----------------------------------------------------------------------------
def cart2geo( xc, yc, prj_cart='+proj=gnom', prj_geo=None, 
              lon0=None, lat0=None, re=6378163, ell=298.26, nzco=None,
              rlt=0.99330647, rota=0, change_lon0_sign=False ) :
    """
    Converts Cartesian coordinates to geographic coordinates.

    Parameters:
        - xc (float): The x-coordinate in Cartesian system.
        - yc (float): The y-coordinate in Cartesian system.
        - prj_geo (str or int, optional): Projection string for the geographic coordinate system. Default is None. It is possible to use an integer to specify the epsg code insted.
        - prj_cart (str or int, optional): Projection string for the Cartesian coordinate system. Default is '+proj=gnom'. It is possible to use an integer to specify the epsg code insted.
        - lon0 (float, optional): The longitude of the origin point. Defaults to None.
        - lat0 (float, optional): The latitude of the origin point. Defaults to None.
        - re (float, optional): The equatorial radius of the Earth. Defaults to 6378163.
        - ell (float, optional): The eccentricity of the Earth. Defaults to 298.26.
        - nzco (int, optional): The flag for the conversion method. Defaults to 5.
        - rlt (float, optional): The ratio of the polar radius to the equatorial radius. Defaults to 0.99330647.
        - rota (float, optional): The rotation angle. Defaults to 0.
        - change_lon0_sign (bool, optional): Flag to change the sign of lon0. Defaults to False.

    Returns:
    - tuple: A tuple containing the longitude and latitude in the geographic system.
    """

    if change_lon0_sign == True :
        lon0 = np.copy( lon0 ) * -1
    
    if ( prj_cart is not None ) and ( nzco not in [ 5 ] ) :
        
        if prj_geo is None :
            prj_geo= f'+proj=longlat +a={re} +rf={ell}'
        
        if ( lon0 is not None ) and ( 'lon_0' not in prj_cart ):
            prj_cart = prj_cart + f' lon_0={lon0}'
    
        if ( lat0 is not None ) and ( 'lat_0' not in prj_cart ) :
            prj_cart = prj_cart + f' lat_0={lat0}'
            
        xg, yg = utl.prjxy( prj_cart, prj_geo, xc*1e3, yc*1e3 )    
        
    # -------------------------------------------------------------------------    
    if nzco == 5 :
        # ---------------------------------------------------------------------
        # The statement "nzco==5" uses short distance conversion 
        # as in Thurber original simulps program.
        # ---------------------------------------------------------------------
        
        rad = 2 * np.pi / 360         

        lat0m, lon0m, xlnkm, xltkm, rota, snr, csr = setorg( lat0, lon0, rota )
        
        fy = csr*yc - snr*xc
        fx = snr*yc + csr*xc
        
        fy = fy / xltkm
        plat = lat0m + fy

        xlt1 = np.arctan( rlt * np.tan( rad*( plat + lat0m ) / 120 ) )
        fx = fx / ( xlnkm * np.cos( xlt1 ) )
        plon = lon0m + fx

        yg = plat / 60
        xg = plon / 60
        
    return xg, yg

# -----------------------------------------------------------------------------
def geo2cart( xg, yg, prj_geo=None, prj_cart='+proj=gnom', rota=0,
              lon0=None, lat0=None, re=6378163.0, ell=298.26, nzco=None,
              rlt=0.99330647, change_lon0_sign=False ) :
    """
    Converts geographic coordinates to Cartesian coordinates.

    Parameters:
        - xg (float or array-like): X-coordinate(s) in geographic coordinate system.
        - yg (float or array-like): Y-coordinate(s) in geographic coordinate system.
        - prj_geo (str or int, optional): Projection string for the geographic coordinate system. Default is None. It is possible to use an integer to specify the epsg code insted.
        - prj_cart (str or int, optional): Projection string for the Cartesian coordinate system. Default is '+proj=gnom'. It is possible to use an integer to specify the epsg code insted.
        - rota (float, optional): Rotation angle in degrees. Default is 0.
        - lon0 (float, optional): Longitude of the origin point. Default is None.
        - lat0 (float, optional): Latitude of the origin point. Default is None.
        - re (float, optional): Equatorial radius of the Earth. Default is 6378163.0.
        - ell (float, optional): Ellipticity of the Earth. Default is 298.26.
        - nzco (int, optional): Coordinate system type. Default is None. Set nzco=5 to use short distance conversion as in Thurber original simulps program.
        - rlt (float, optional): Ratio of the polar radius to the equatorial radius. Default is 0.99330647.
        - change_lon0_sign (bool, optional): Flag to change the sign of lon0. Default is False.

    Returns:
    - xc (float or array-like): X-coordinate(s) in Cartesian coordinate system.
    - yc (float or array-like): Y-coordinate(s) in Cartesian coordinate system.
    """

    if change_lon0_sign == True :
        lon0 = np.copy( lon0 ) * -1
    
    if ( prj_cart is not None ) and ( nzco not in [ 5, 3 ] ) :
        
        prj_geo= f'+proj=longlat +a={re} +rf={ell}'
        
        if ( lon0 is not None ) and ( 'lon_0' not in prj_cart ):
            prj_cart = prj_cart + f' lon_0={lon0}'
    
        if ( lat0 is not None ) and ( 'lat_0' not in prj_cart ) :
            prj_cart = prj_cart + f' lat_0={lat0}'
            
        xc, yc = utl.prjxy( prj_geo, prj_cart, xg, yg )
        
        xc = xc/1e3 # m to km
        yc = yc/1e3 # m to km
     
    # -------------------------------------------------------------------------    
    if nzco == 5 :
        # ---------------------------------------------------------------------
        # The statement "nzco==5" uses short distance conversion 
        # as in Thurber original simulps program.
        # It calculates distance of station or event
        # from given coordinate origin in terms of (possibly rotated) 
        # cartesian coords x by using the short distance conversion factors 
        # from the function setorg
        # ---------------------------------------------------------------------
        
        if type(xg) in ( list, tuple ) :
            xg = np.array( xg )

        if type(yg) in ( list, tuple ) :
            yg = np.array( yg )            
        
        if prj_geo is not None :

            xg, yg = utl.prjxy( prj_geo, f'+proj=longlat +a={re} +rf={ell}', xg, yg )

        # deg to rad conversion factor
        rad = 2 * np.pi / 360 
        
        xgm = xg * 60
        ygm = yg * 60
        
        lat0m, lon0m, xlnkm, xltkm, rota, snr, csr = setorg( lat0, lon0, rota )
        
        # convert lat and lon differences to km
        x = xgm - lon0m
        y = ygm - lat0m
        xlt1 = np.arctan( rlt * np.tan( ( rad * ( ygm + lat0m ) / 120 ) ) )
        xc = x * xlnkm * np.cos( xlt1 ) 
        yc = y * xltkm
        
        # now do rotation
        if rota != 0 :
            ty = csr*yc + snr*xc
            xc = csr*xc - snr*yc
            yc = ty
       # ----------------------------------------------------------------------
       
    return xc, yc


# -----------------------------------------------------------------------------
def setorg( lat0, lon0, rota=0, rtl=0.99330647, re=6378163, ell=298.26 ) :
    """
    This function establishes the short distance conversion factors
    given the origin of coordinates
    the rotation angle is converted to radians also
    """
    
    # deg to rad conversion factor
    rad = 2 * np.pi / 360 
    
    # Convert lat0 and lon0 from decimal degrees to minuts
    lat0m = lat0 * 60
    lon0m = lon0 * 60
    
    # conversion factor for latitude
    dlt1 = np.arctan( rtl * np.tan( lat0m * rad / 60 ) )
    dlt2 = np.arctan( rtl * np.tan( ( lat0m + 1 ) * rad / 60 ) )
    dlt = dlt2 - dlt1
    r = ( re / 1e3 ) * ( 1 - ( ( np.sin( dlt1 )**2 ) / ell ) )
    xltkm = r * dlt 
    
    # conversion factor for longitude
    dlt = np.arccos( 1 - ( 1 - np.cos( rad/60 ) ) * np.cos( dlt1 )**2 )
    bc = r * dlt
    xlnkm = bc / np.cos( dlt1 )
    
    # convert coordinates with rotation cosines
    rota = rota * rad 
    snr = np.sin( rota )
    csr = np.cos( rota )
    
    return lat0m, lon0m, xlnkm, xltkm, rota, snr, csr

# -----------------------------------------------------------------------------
def transform_lim( lim, ttype='geo2cart', 
                   prj_geo=None, prj_cart='+proj=gnom', 
                   rota=0, lon0=None, lat0=None, 
                   re=6378163.0, ell=298.26, nzco=None,
                   rlt=0.99330647, change_lon0_sign=False ) :
    """
    Transforms the given limit coordinates based on the specified transformation type.

    Parameters:
    - lim (list): The limit coordinates [xmin, xmax, ymin, ymax].
    - ttype (str): The transformation type. Default is 'geo2cart'.
    - prj_geo (str): The projection string for the geographic coordinates. Default is None.
    - prj_cart (str): The projection string for the Cartesian coordinates. Default is '+proj=gnom'.
    - rota (float): The rotation angle in degrees. Default is 0.
    - lon0 (float): The longitude of the origin point. Default is None.
    - lat0 (float): The latitude of the origin point. Default is None.
    - re (float): The equatorial radius of the Earth. Default is 6378163.0.
    - ell (float): The eccentricity of the Earth. Default is 298.26.
    - nzco (float): The vertical coordinate of the north pole. Default is None.
    - rlt (float): The ratio of the length of the tangent to the length of the arc. Default is 0.99330647.
    - change_lon0_sign (bool): Whether to change the sign of lon0. Default is False.

    Returns:
    - list: The transformed limit coordinates [xmin, xmax, ymin, ymax].
    """
    
    lim = utl.copy.deepcopy( lim )

    xmin = lim[0]
    xmax = lim[1]
    ymin = lim[2]
    ymax = lim[3]

    if ttype == 'geo2cart' :

        xmin, ymin = geo2cart( xmin, ymin, prj_geo=prj_geo, prj_cart=prj_cart, 
                               rota=rota, 
                               lon0=lon0, lat0=lat0, re=re, ell=ell, nzco=nzco, 
                               rlt=rlt, change_lon0_sign=change_lon0_sign )

        xmax, ymax = geo2cart( xmax, ymax, prj_geo=prj_geo, prj_cart=prj_cart, 
                               rota=rota, 
                               lon0=lon0, lat0=lat0, re=re, ell=ell, nzco=nzco, 
                               rlt=rlt, change_lon0_sign=change_lon0_sign )
    if ttype == 'cart2geo' :

        xmin, ymin = cart2geo( xmin, ymin, prj_geo=prj_geo, prj_cart=prj_cart, 
                               rota=rota, 
                               lon0=lon0, lat0=lat0, re=re, ell=ell, nzco=nzco, 
                               rlt=rlt, change_lon0_sign=change_lon0_sign )
        
        xmax, ymax = cart2geo( xmax, ymax, prj_geo=prj_geo, prj_cart=prj_cart, 
                               rota=rota, 
                               lon0=lon0, lat0=lat0, re=re, ell=ell, nzco=nzco, 
                               rlt=rlt, change_lon0_sign=change_lon0_sign )
    lim[0] = xmin
    lim[1] = xmax
    lim[2] = ymin
    lim[3] = ymax

    return lim

# -----------------------------------------------------------------------------
def validate_tt_st( tt_data, 
                    station_data, 
                    path=None, 
                    file_name=None,
                    threshold_P=4, 
                    threshold_S=0,
                    new_tt_file=False, 
                    log_file=False, 
                    new_tt_dict=True, 
                    renum_eve=True, 
                    printf=False, 
                    keepS=True,
                    prefix='\t' ) :

    """
    Validates the travel-time data and station data.

    Parameters:

        - tt_data (str or dict): The travel-time data. 
            If a string is provided, it is assumed to be the path to a file containing the travel-time data. 
            If a dictionary is provided, it is assumed to be the travel-time data itself.
        
        - station_data (str or dict): The station data. 
            If a string is provided, it is assumed to be the path to a file containing the station data. 
            If a dictionary is provided, it is assumed to be the station data itself.
        
        - path (str, optional): The path where the new travel-time file will be created. 
            If not provided, the current working directory is used.
        
        - file_name (str, optional): The name of the new travel-time file. 
            If not provided, the default name 'tt_data.in' is used.
        
        - threshold_P (int, optional): The threshold for the number of P data required 
            for an event to be considered valid. Defaults to 4.

        - threshold_S (int, optional): The threshold for the number of S data required
            for an event to be considered valid. Defaults to 0.
        
        - new_tt_file (bool, optional): Whether to create a new travel-time file with the validated data. 
            Defaults to False.
        
        - log_file (bool, optional): Whether to create a log file with the validation results. 
            Defaults to False.
        
        - new_tt_dict (bool, optional): Whether to return the validated travel-time data as a dictionary. 
            Defaults to True.
        
        - change_lon0_sign (bool, optional): Whether to change the sign of the reference longitude. 
            Defaults to True.
        
        - renum_eve (bool, optional): Whether to renumber the events in the travel-time data. 
            Defaults to True.
        
        - printf (bool, optional): Whether to print the validation results. Defaults to False.
        
        - keepS (bool, optional): Whether to keep S data in the validated travel-time data. Defaults to True.

    Returns:
        
        - If new_tt_dict is True: A tuple containing the validated 
            travel-time data as a dictionary and a boolean array indicating the validity of each record.
        
        - If new_tt_dict is False: A boolean array indicating the validity of each record.
    """

    if type( tt_data ) == str :
        
        tt_dict = read_tt_file( tt_data )[0]

    else :

        tt_dict = tt_data.copy()
        
        if path is None :
            path = os.getcwd()
        
        if file_name is None :
            file_name = 'tt_data.in'
        
    if type( station_data ) == str :

        station_dict = read_station_file( station_data )

    else :
        
        station_dict = station_data.copy()

    if 'p0s1' not in tt_dict :
        tt_dict['p0s1'] = np.full( tt_dict['eve'].shape, 0 )
        for i, r in enumerate( tt_dict['phase_descriptor'] ) :
            if 'S' in r :
                tt_dict['p0s1'][i] = 1

    n_tt_dict = tt_dict.copy()
    un_st = np.unique( tt_dict['st'] )
    un_eve = np.unique( tt_dict['eve'] )
    
    log = ''
    idx = np.full( tt_dict['st'].shape, True )

    log = log + "VALIDATION LOG \n"
    log = log + "# ---------------------------------------------------\n"
    print( "\n\n\tVALIDATION LOG")
    print( f"{prefix}# ---------------------------------------------------\n" )
    
    # ---------------------------------------------------
    idxnst =utl.copy.copy( ~idx )
    for nm in un_st :
        
        idxtt = tt_dict['st'] == nm
        idxst = station_dict['st'] == nm
        
        if np.sum( idxst ) == 0 :
            string = f'Station {nm} is not in station file: {sum(idxtt)} tt data associated (S and/or P)'
            if printf :
                print( string )
            log = log + string + ' \n'
            
            for k in tt_dict :
                idxnst = idxnst | idxtt
                
    print( f'\tTotal travel-time data with station not found: {np.count_nonzero(idxnst)} / {np.size(idxnst)}' )
    log = log + f'Total travel-time data with station not found: {np.count_nonzero(idxnst)} / {np.size(idxnst)} \n' 

    log = log + "# --------------------------------------------------- \n"
    print( "\t# --------------------------------------------------- \n" )
    
    # Boolean array with True value indicating the correct records
    idx = idx & ~idxnst

    # ---------------------------------------------------
    idxs = np.full( tt_dict['st'].shape, False) 
    for i, r in enumerate( tt_dict['phase_descriptor'] ) :

        if 'S' in r :

            neve = tt_dict['st'][ tt_dict['eve'] == tt_dict['eve'][i] ]

            if sum( neve == tt_dict['st'][i] ) != 2 :

                idxs[ i ] = True
                string = f"S phase of event {tt_dict['dd'][i]}/{tt_dict['mm'][i]}/{tt_dict['yy'][i]},"+\
                         f" rec. by station {tt_dict['st'][i]}, does not have any associated P phase \n"
                log = log + string
                if printf :
                    print( '\t' + string )

    print( f"\tTotal S data with no matching P data: {np.count_nonzero( idxs )} / {np.size(idxs)}" )
    log = log + f"Total S data with no matching P data: {np.count_nonzero( idxs )} / {np.size(idxs)} \n"  

    log = log + "# --------------------------------------------------- \n"
    print( "\t# --------------------------------------------------- \n" )

    # Boolean array with True value indicating the correct records
    idx = idx & ~idxs

    # ---------------------------------------------------
    if keepS is False :

        idxsd = ( tt_dict[ 'p0s1' ] == 1 ) 

        print( f"\tTotal recorded S data: {np.count_nonzero(idxsd)} / {np.size(idxsd)}" )
        log = log + f"Total recorded S data: {np.count_nonzero(idxsd)} / {np.size(idxsd)} \n" 

        log = log + "# --------------------------------------------------- \n"
        print( "\t# --------------------------------------------------- \n" )

        # Boolean array with True value indicating the correct records
        idx = idx & ~idxsd

    # ---------------------------------------------------
    idxeP = np.full( tt_dict['eve'].shape, False)

    for i, e in enumerate( np.unique( tt_dict[ 'eve' ] ) ) :
        
        ne = np.count_nonzero( (tt_dict[ 'eve' ] == e) & (tt_dict[ 'p0s1' ] == 0) & idx )

        if ne < threshold_P :
            idxeP = idxeP | ( tt_dict[ 'eve' ] == e )

    un_eve_th = np.unique( tt_dict[ 'eve' ][ idxeP ] )
    print( f"\tTotal events with less than {threshold_P} P data: {np.size(un_eve_th)} / {np.size(un_eve)} " )
    print( f"\tTotal records from events with less than {threshold_P} P data: {np.count_nonzero(idxeP)} / {np.size(idxeP)} " )
    log = log + f"Total events with less than {threshold_P} P data: {np.size(un_eve_th)} / {np.size(un_eve)}\n"
    log = log + f"Total records with less than {threshold_P} P data: {np.count_nonzero(idxeP)} / {np.size(idxeP)}\n"

    log = log + "# --------------------------------------------------- \n"
    print( "\t# --------------------------------------------------- \n" )

    # Boolean array with True value indicating the correct records
    idx = idx & ~idxeP

    # ---------------------------------------------------
    idxeS = np.full( tt_dict['eve'].shape, False)

    for i, e in enumerate( np.unique( tt_dict[ 'eve' ] ) ) :
        
        ne = np.count_nonzero( (tt_dict[ 'eve' ] == e) & (tt_dict[ 'p0s1' ] == 1) & idx )

        if ne < threshold_S :
            idxeS = idxeS | ( tt_dict[ 'eve' ] == e )

    un_eve_th = np.unique( tt_dict[ 'eve' ][ idxeS ] )
    print( f"\tTotal events with less than {threshold_S} S data: {np.size(un_eve_th)} / {np.size(un_eve)} " )
    print( f"\tTotal records from events with less than {threshold_S} S data: {np.count_nonzero(idxeS)} / {np.size(idxeS)} " )
    log = log + f"Total events with less than {threshold_S} S data: {np.size(un_eve_th)} / {np.size(un_eve)}\n"
    log = log + f"Total records with less than {threshold_S} S data: {np.count_nonzero(idxeS)} / {np.size(idxeS)}\n"

    log = log + "# --------------------------------------------------- \n"
    print( "\t# --------------------------------------------------- \n" )

    # Boolean array with True value indicating the correct records
    idx = idx & ~idxeS

    # ---------------------------------------------------
    # Filtering the final dictionary
    for k in n_tt_dict :
        n_tt_dict[ k ] = n_tt_dict[ k ][ idx ]

    # ---------------------------------------------------
    # Print conclusions
    n_un_eve = np.unique( n_tt_dict[ 'eve' ] )
    print( f"\tTOTAL VALID DATA: {np.count_nonzero(idx)} / {np.size(idx)}" )
    log = log + f"\nTOTAL VALID DATA: {np.count_nonzero(idx)} / {np.size(idx)} \n"
    print( f"\tTOTAL VALID EVENTS: {np.size(n_un_eve)} / {np.size(un_eve)}" )
    log = log + f"TOTAL VALID EVENTS: {np.size(n_un_eve)} / {np.size(un_eve)} \n"
    log = log + "# --------------------------------------------------- \n"
    print( "\t# --------------------------------------------------- \n" )
 
    # ---------------------------------------------------
    # Create new travel-time file
    if new_tt_file is True :
        
        create_tt_file( path + utl.os.sep+ file_name, n_tt_dict )
          
    if log_file is True :

        with open( path +s+ 'log_tt_st', 'w' ) as f :
              f.write( log )
              f.close()

    # ---------------------------------------------------
    # Renumbering events
    if renum_eve is True :
        num = 1
        new_eve = np.zeros( n_tt_dict['eve'].size )
        
        for i, e in enumerate( n_tt_dict['eve'] ) :
            if i == 0 :
                new_eve[i] = num
            else :
                if n_tt_dict['eve'][i] != n_tt_dict['eve'][i-1] :
                    num += 1
                new_eve[i] = num
    
        n_tt_dict['eve'] = new_eve
    
    # ---------------------------------------------------
    # Return 
    if new_tt_dict is True : 
        return n_tt_dict, idx
    else :
        return idx

# -----------------------------------------------------------------------------
def wadati_plot( loc_data, s=1, plot=True, time_threshold=None ) :
    """
    Calculate and plot the Wadati diagram for seismic travel time data.

    Parameters:
    loc_data (str or dict): The seismic travel time data. If it is a string, it should be the file path to the data. If it is a dictionary, it should contain the following keys: 'eve', 'yy', 'mm', 'dd', 'h', 'm', 's', 'yy0', 'mm0', 'dd0', 'h0', 'm0', 's0', 'phase_descriptor', and 'st'.
    s (int, optional): The size of the scatter points in the plot. Default is 1.
    plot (bool, optional): Whether to plot the Wadati diagram. Default is True.
    time_threshold (float, optional): The time threshold for filtering the data points. Only data points within the threshold will be plotted. Default is None.

    Returns:
    vpvs (float): The Vp/Vs ratio calculated from the linear regression.
    iout (ndarray): A boolean array indicating which data points are within the time threshold.

    """

    if type( loc_data ) == str :

        tt_dict = loc_data( loc_data )[0]

    else :

        tt_dict = loc_data.copy()

    tt_dict['tp'] = np.full( tt_dict['eve'].size, np.nan )

    tt_dict['ts-tp'] = np.full( tt_dict['eve'].size, np.nan )

    tt_dict['tti'] = utl.combine64( years=tt_dict['yy'], months=tt_dict['mm'], 
                                    days=tt_dict['dd'], hours=tt_dict['h'], 
                                    minutes=tt_dict['m'], seconds=tt_dict['s'] )

    tt_dict['tt0'] = utl.combine64( years=tt_dict['yy0'], months=tt_dict['mm0'], 
                                    days=tt_dict['dd0'], hours=tt_dict['h0'], 
                                    minutes=tt_dict['m0'], seconds=tt_dict['s0'] )

    tt_dict['tt'] = ( tt_dict['tti'].astype(float) - tt_dict['tt0'].astype(float) )/1e9

    for i,e in enumerate( tt_dict['eve'] ) :
        if i == np.size( tt_dict['eve'] ) -2: 
            break

        if 'P' in tt_dict['phase_descriptor'][i] :
            if 'S' in tt_dict['phase_descriptor'][i+1] :
                if tt_dict['st'][i+1] == tt_dict['st'][i] :
                    tt_dict['tp'][i] = tt_dict['tt'][i]
                    tt_dict['ts-tp'][i] = tt_dict['tt'][i+1]

    # linear regression
    idx = np.isnan( tt_dict['ts-tp'] )
    X = tt_dict['tp'][~idx][:,np.newaxis]
    a = np.linalg.lstsq( X, tt_dict['ts-tp'][~idx], rcond=None )[0]
    vpvs = np.round( a[0], 2 )
    y = tt_dict[ 'tp' ] * a[ 0 ]

    iout = np.full( tt_dict['tp'].shape, True )
    if time_threshold is not None :
        My = y + time_threshold
        my = y - time_threshold
        iout = ( my < tt_dict['ts-tp'] ) & ( My > tt_dict['ts-tp'] ) 

    if plot == True :
        lsz_plot.plt.scatter( tt_dict['tp'][iout], tt_dict['ts-tp'][iout], s=s, c='b' )
        lsz_plot.plt.plot( tt_dict['tp'][~idx], y[~idx], c='r', label=f'Vp/Vs: {vpvs}' )

        if time_threshold is not None :

            lsz_plot.plt.plot( tt_dict['tp'][~idx], my[~idx], c='k', linestyle='dashed' )
            lsz_plot.plt.plot( tt_dict['tp'][~idx], My[~idx], c='k', linestyle='dashed' )
            lsz_plot.plt.scatter( tt_dict['tp'][~iout], tt_dict['ts-tp'][~iout], s=s, c='r' )
        lsz_plot.plt.xlabel( 'Tp [ sec ]' )
        lsz_plot.plt.ylabel( 'Ts - Tp [ sec ]' )
        lsz_plot.plt.legend()

    iout = np.isin( tt_dict['eve'], tt_dict['eve'][~iout & ~idx ] ) &\
           np.isin( tt_dict['st'], tt_dict['st'][~iout & ~idx ] )

    return vpvs, iout, tt_dict

# -------------------------------------------------------------------------
def map_rec_eve(xy_events=None, xy_stations=None, raster_map=raster_globe, 
                rcolor='r', rsize=None, rmarker='^',
                ecolor='b', esize=None, emarker='o', 
                lim=None, prjcode=4326, extend_lim=5, aspect='auto', rgb=True):
    """
    Plot the receivers and events on a map.

    Parameters:
    - xy_events (list, tuple, dict, np.ndarray): The coordinates of the events.
    - xy_stations (list, tuple, dict, np.ndarray): The coordinates of the stations.
    - raster_map: The raster map to be used as the background.
    - rcolor (str): The color of the receiver markers.
    - rsize: The size of the receiver markers.
    - rmarker (str): The marker style of the receiver markers.
    - ecolor (str): The color of the event markers.
    - esize: The size of the event markers.
    - emarker (str): The marker style of the event markers.
    - receiver: Not used.
    - lim (list): The limits of the map.
    - prjcode (int): The projection code.
    - extend_lim (int): The amount by which to extend the map limits.
    - aspect (str): The aspect ratio of the plot.
    - rgb (bool): Whether to use RGB color for the raster map.

    Returns:
    None
    """

    if type( xy_stations ) in ( list, tuple ) :
        latr = np.asarray( xy_stations[1] )
        lonr = np.asarray( xy_stations[0] )

    if type( xy_stations ) == dict :
        latr = np.asarray( xy_stations['lat'] )
        lonr = np.asarray( xy_stations['lon'] )

    if type( xy_stations ) == np.ndarray :
        latr = xy_stations[:,1]
        lonr = xy_stations[:,0]

    if xy_stations is None :
        latr = np.nan
        lonr = np.nan

    if type( xy_events ) in ( list, tuple ) :
        late = np.asarray( xy_events[1] )
        lone = np.asarray( xy_events[0] )

    if type( xy_events ) == dict :
        late = np.asarray( xy_events['lat'] )
        lone = np.asarray( xy_events['lon'] )

    if type( xy_events ) == np.ndarray :
        late = xy_events[:,1]
        lone = xy_events[:,0]

    if xy_events is None :
        latr = np.nan
        lonr = np.nan
        
    if lim is None :
        lim = [ np.min( ( np.min(lone), np.min(lonr) ) ),
                np.max( ( np.max(lone), np.max(lonr) ) ),
                np.min( ( np.min(late), np.min(latr) ) ),
                np.max( ( np.max(late), np.max(latr) ) ) ]
        
    if raster_map is not None : 
        map_lim = utl.extend_lim( lim, extend_lim )
        
        r_m_crop = lsz_plot.rt.raster_warp( raster_map, lim=map_lim, 
                                    lim_prjcode=prjcode, out_prjcode=prjcode, 
                                    close=False )
        ax = lsz_plot.rt.pltr( r_m_crop, axis=True, rgb=rgb )

    else :
        lsz_plot.plt.figure()
        
        ax = lsz_plot.plt.subplot( 111 )

    ax.scatter( lonr, latr, c=rcolor, s=rsize, marker=rmarker,
                            label='receivers', aspect=aspect )
    
    ax.scatter( lone, late, c=ecolor, s=esize, marker=emarker, 
                            label='events', aspect=aspect )
    
    ax.legend()

    lsz_plot.plt.show()

# -----------------------------------------------------------------------------
def split_tt_file( tt_data, path=None, name='tt_split', idx=None ) :
    """
    Splits a travel-timefile into multiple files based on unique event IDs.

    Args:
        
        - tt_data (str or dict): The travel-timedata. 
            If a string, it is assumed to be the path to the travel-timefile.
            If a dictionary, it is assumed to be the travel-timedata already loaded into memory.
        
        - path (str, optional): The path where the split travel-timefiles will be saved. 
            If not provided, the path of the
            original travel-timefile will be used.
        
        - name (str, optional): The base name for the split travel-timefiles. 
            Defaults to 'phs_split'.
        
        - idx (int or list, optional): The indices of the travel-timedata to be used for splitting. 
            If provided, only the specified indices will be used. 
            Defaults to None, which means all indices will be used.

    Returns:
        str: The path where the split travel-timefiles are saved.

    """

    # If the travel-timedata is a string, read the travel-timefile
    if type( tt_data ) == str :
        tt_dict = read_tt_file( tt_data )
    else :
        tt_dict = tt_data.copy()

    # If an index is provided, select the corresponding elements from each array in the travel-timedictionary
    if idx is not None :
        for k in tt_dict :
            tt_dict[k] = tt_dict[k][idx]

    # Find the unique events in the travel-timedictionary
    events = np.unique( tt_dict['eve'] )

    # If no path is provided, create a default path using the travel-timedata file name
    if path is None :
        path = os.path.split( tt_data )[0] +s+ 'tt_split'

    if os.path.exists( path ) == False :
        os.makedirs( path )

    # For each unique event
    for i in events :
        # Find the indices where the event matches the current unique event
        idx = tt_dict['eve'] == i

        tt_dict_e = {}
        # Create a travel-timedictionary for the current unique event
        for k in tt_dict :
            tt_dict_e[k] = tt_dict[k][idx]

        # Create a travel-time file for the current unique event
        create_tt_file( path +s+ name + f'_e{i}.obs', tt_dict_e )
        
    # Return the path where the travel-time files were created
    return path

# -----------------------------------------------------------------------------
def run_cmd( string2run, 
             path=None, 
             addbefore='',
             addafter=''):
    """
    Executes a given command in a specified directory.

    Parameters:
        
        - string2run (str): The command to be executed.
        
        - path (str, optional): The directory where the command should be run. 
            If not specified, the command will be run in the current working directory.
    
        - make_executable (bool, optional): If True, the command will be made executable 
            before running it. Defaults to False.

        - addbefore (str, optional): A string to be added before the command. Defaults to ''.

        - addafter (str, optional): A string to be added after the command. Defaults to ''.

    Returns:
        None

    N.B.
    This function changes the current working directory to the specified path, 
    makes the command executable if make_executable is True, and then runs the command. 
    If a path is specified, the function will change back to the original working directory 
    after the command is run.
    """

    if path is not None :
        # Go to the main code directory
        cwd = os.getcwd() # current working directory
        os.chdir( path )
    
    # Run job
    if platform.system() == 'Linux' :
        
        string2run = addbefore + string2run + addafter
        os.system( string2run )
    
    if path is not None :
        # Go back to the original working directory
        os.chdir( cwd )

# -----------------------------------------------------------------------------
def filt_dict( dictionary, idx, keys=None ):
    """
    Filter a dictionary by selecting specific indices and keys.

    Args:
        - dictionary (dict): The input dictionary to be filtered.
        - idx (list): The indices to select from the dictionary values.
        - keys (list, optional): The keys to filter the dictionary by. If None, all keys will be considered.

    Returns:
        - dict: The filtered dictionary.

    """

    dictionary = utl.copy.deepcopy( dictionary )

    for k in dictionary :

        if dictionary[ k ] is not None :

            if type( dictionary[ k ] ) in (list, tuple, set, int, float) :

                dictionary[ k ] = np.array( dictionary[ k ] )
            
            if np.size( dictionary[ k ] ) == np.size( idx ) :
                
                if keys :
                    if k in keys :
                        dictionary[ k ] = dictionary[ k ][idx]
                else :
                    dictionary[ k ] = dictionary[ k ][idx] 

    return dictionary

