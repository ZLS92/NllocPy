# -*- coding: utf-8 -*-lzampa
"""
Created on Wed Oct 14 18:58:52 2020

@author: lzampa
"""

# -----------------------------------------------------------------------------
import os
import pyproj as prj
import numpy as np
import scipy as sp
from matplotlib.widgets import LassoSelector
from matplotlib.colors import LightSource
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import itertools
import random
import time
import copy
from matplotlib.path import Path
import datetime
from osgeo import gdal, osr, ogr
import filecmp
import shutil
import platform
import tempfile
import sys
import html
from matplotlib import cm
import io 
import pdfkit

# Disable gdal exceptions
gdal.UseExceptions()

# -----------------------------------------------------------------------------
# Constants

G = 6.6742*1e-11 # [m3/(kg *s^2)]
M = 5.97*1e24 # [kg]
a_wgs84 = 6378137 # [m]
c_wgs84 = 6356758 # [m]
R_wgs84 = ((a_wgs84**2)*c_wgs84)**(1/3) # [m]
J2_wgs84 = 1.081874*1e-3
w_wgs84 = 7.292115*1e-5 # [rad/sec]

# -----------------------------------------------------------------------------
def cm2in( cm ) :
    
    inches = cm * 1/2.54 
    
    return inches 

# -----------------------------------------------------------------------------
def tmi(t1=None, t2=None):
    """
    Print time interval from time t1 to t2
    --- If "t2 = None" , then, "t2 = now"
    """
    
    if ( t1 is None ) and ( t2 is None ) :
        return time.time()
    
    if t2 == None :
        print( time.time() - t1 )
    else:
        print( t2 - t1 )
        
# -----------------------------------------------------------------------------
def stat( array, decimals=None, printf=True, out=tuple, multilines=False, sep='' ):
    """
    Print statistics of a numpy array
    """        

    Min, Max, Mean, Std = ( np.nanmin(array), np.nanmax(array),
                            np.nanmean(array), np.nanstd(array) )
    

    if decimals != None:
        Min, Max, Mean, Std = ( np.round(Min, decimals), np.round(Max, decimals),
                                np.round(Mean, decimals), np.round(Std, decimals) )   
    
    len_ = np.max( ( len(str(Min)), len(str(Max)), len(str(Mean)), len(str(Std)) ) )      
    
    string = f"Min:{Min}{sep} Max:{Max}{sep} Mean:{Mean}{sep} Std:{Std}"
   
    if multilines == True :
        string = f"Min  = {Min: >{len_}} \nMax  = {Max: >{len_}} \nMean = {Mean: >{len_}} \nStd  = {Std: >{len_}}"
    
    if printf == True :
        print( string ) 
        
    if out == str :
        return string
    
    if out == list :     
        return [ Min, Max, Mean, Std ]
    
    if out == tuple :     
        return Min, Max, Mean, Std

    if out == dict :     
        return { 'Min':Min, 'Max':Max, 'Mean':Mean, 'Std':Std } 
    
# -----------------------------------------------------------------------------
def dms2dd(degrees, minutes=0.0, seconds=0.0):
    """
    Convert degrees, minutes, and seconds to decimal degrees.

    Parameters:
    degrees (int, float, list, tuple, numpy.ndarray): Degrees.
    minutes (int, float, list, tuple, numpy.ndarray): Minutes.
    seconds (int, float, list, tuple, numpy.ndarray): Seconds.

    Returns:
    numpy.ndarray: Decimal degrees.
    """
    # Convert inputs to numpy arrays if they are not already
    if not isinstance(degrees, np.ndarray):
        degrees = np.array(degrees)
    if not isinstance(minutes, np.ndarray):
        minutes = np.array(minutes)
    if not isinstance(seconds, np.ndarray):
        seconds = np.array(seconds)

    return degrees + minutes / 60 + seconds / 3600

# -----------------------------------------------------------------------------
def dd2dms(decimal_degrees, precision=4 ):
    """
    Convert decimal degrees to degrees, minutes, and seconds.

    Parameters:
        - decimal_degrees (int, float, list, tuple, numpy.ndarray): Decimal degrees.
        - precision (int): Number of decimal places to round seconds to.

    Returns:
    tuple: Degrees, minutes, and seconds.
    """

    # Convert inputs to numpy arrays if they are not already
    if not isinstance(decimal_degrees, np.ndarray):
        decimal_degrees = np.array(decimal_degrees)

    degrees = np.floor(decimal_degrees)
    minutes = np.floor((decimal_degrees - degrees) * 60)
    seconds = (decimal_degrees - degrees - minutes / 60) * 3600
    seconds = np.round(seconds, precision)

    return degrees, minutes, seconds

# -----------------------------------------------------------------------------
def extend_lim( lim, d, method='percentage', sqr_area=False, plot=False, 
                c1='r', c2='b' ):
    """
    - lim = [xmin, xmax, ymin, ymax]
    - two methods aveilable:
        1) distance --> extend for a distance in the same units of lim
        2) percentage --> extend for a percentage of the total xy lengths of the area
    - sqr_area = True --> force the output boundaries to have length_X = length_Y
    """
    
    xd, yd = lim[1] - lim[0], lim[3] - lim[2]
    if sqr_area == True:
        if xd > yd: lim2 = [lim[0], lim[1], lim[2] - (xd - yd) / 2, lim[3] + (xd - yd) / 2]
        if yd > xd: lim2 = [lim[0] - (yd - xd) / 2, lim[1] + (yd - xd) / 2, lim[2], lim[3]]
        xd, yd = lim2[1] - lim2[0], lim2[3] - lim2[2]
    else: lim2 = lim

    if np.size(d) == 1: dx, dy = d, d
    if np.size(d) == 2: dx, dy = d[0], d[1]
    if method == 'percentage':
        xp, yp = xd * dx / 100, yd * dy / 100
        lim_new = [lim2[0] - xp, lim2[1] + xp, lim2[2] - yp, lim2[3] + yp]
    if method == 'distance': lim_new = [lim2[0] - dx, lim2[1] + dx, lim2[2] - dy, lim2[3] + dy]

    if plot == True:
        xplot1 = [lim[0], lim[0], lim[1], lim[1], lim[0]]
        yplot1 = [lim[2], lim[3], lim[3], lim[2], lim[2]]
        xplot2 = [lim_new[0], lim_new[0], lim_new[1], lim_new[1], lim_new[0]]
        yplot2 = [lim_new[2], lim_new[3], lim_new[3], lim_new[2], lim_new[2]]
        plt.plot(xplot1, yplot1, c=c1 )
        plt.plot(xplot2, yplot2, c=c2 )
        
    return lim_new

# -----------------------------------------------------------------------------
def prj_( prjcode ):
    """
    Convert epsg or proj4 string to proj4 code object
    """
    if type( prjcode )==int:
        prj4_obj = prj.CRS( 'epsg:'+str( prjcode ) )
        
    if type( prjcode )==str:
        prj4_obj = prj.CRS( prjcode )
        
    return prj4_obj

# -----------------------------------------------------------------------------
def prj_units( prjcode ):
    """
    Return EPGS or Proj4 code units ( i.e., meters or degree )
    """
    
    prj4_obj = prj_( prjcode )
    units = prj4_obj.axis_info[0].unit_name
        
    return units

# -----------------------------------------------------------------------------
def prjxy( prjcode_in, prjcode_out, x, y, z=None ):
    """
    Transform coordinates from a reference system to another
    """
    
    x = copy.copy( x )
    y = copy.copy( y )
    
    if prjcode_in == prjcode_out :
        prj_coord = x, y
        
    else :
        prj_in = prj_( prjcode_in )
        prj_out = prj_( prjcode_out )
        trans = prj.Transformer.from_crs( prj_in, prj_out, always_xy=True )
        if z != None :
            prj_coord = trans.transform( x, y, z )
        if z == None :
            prj_coord = trans.transform( x, y )

    return prj_coord

# -----------------------------------------------------------------------------
def prj2epsg( prjcode ) :
    
    crs = prj.CRS( prjcode )
    epsg = crs.to_epsg()
    
    return epsg

# -----------------------------------------------------------------------------
def lim_sort( lim, p=False ):

    lim_out = [lim[0], lim[2], lim[1], lim[3]]

    if p == True :
        print( lim_out )

    return lim_out

# -----------------------------------------------------------------------------
def lim2points( lim ) :

    x = np.array( ( lim[0], lim[1], lim[1], lim[0] ) )
    y = np.array( ( lim[3], lim[3], lim[2], lim[2] ) )

    return x, y

# -----------------------------------------------------------------------------
def prj_lim(lim, prjcode_in, prjcode_out, sort='xlyl'):
    """
    Transform limits of an area (i.e. [LonMin, LonMax, LatMin, LatMax]) from a
    reference system to another
    """
    if sort=='xlyl':
        x,y = prjxy(prjcode_in, prjcode_out, (lim[0],lim[1]), (lim[2],lim[3]))
    if sort=='xyl':
        x,y = prjxy(prjcode_in, prjcode_out, (lim[0],lim[2]), (lim[1],lim[3]))

    limf = [x[0],x[1],y[0],y[1]]

    return limf

# -----------------------------------------------------------------------------
def prj_centre(xy_centre, prj_type, ellps='WGS84'):
    """
    Creates a PROJ string code, with a reference system centered (i.e tangent)
    to a defined "cental point" of the earth surface: xy_centre,
    i.e. [Lon_centre, Lat_centre] (Lat&Lon in geographic coordinates)
    """
    x_centre, y_centre = xy_centre
    if type(prj_type)==str:
       prj_cent_str = f'''+proj={prj_type} +lat_0={x_centre} +lon_0={y_centre} +ellps={ellps}'''
    if type(prj_type)==int:
       prj_cent_str = f'''+epsg={prj_type} +lat_0={x_centre} +lon_0={y_centre} +ellps={ellps}'''

    return prj_cent_str

# -----------------------------------------------------------------------------
def xy_in_lim( x, y, lim, extend=0, method='percentage', plot=False, s=1, plot_lim=None ):
    """
    Return all xy points within the given limits [xmin, xmax, ymin, ymax]
    """
    
    if extend != 0 :
        lim = extend_lim( lim, d=extend, method=method, plot=plot_lim )
    
    idx = (x >= lim[0]) & (x <= lim[1]) & (y >= lim[2]) & (y <= lim[3])
    xl = x[ idx ]
    yl = y[ idx ]
    
    if plot is True :
        plt.scatter( x[~idx], y[~idx], c='b', s=s )
        plt.scatter( xl, yl, c='r', s=s )

    return xl, yl, idx

# -----------------------------------------------------------------------------
#@jit( forceobj=True )
def find_by_dist( xy1, xy2, dist, prj_dist=None, prj_in=4326, plot=False ):
    
    """
    Find all points of a set (xy1) within a distance (dist)
    from another set of points (xy2)
    """

    set1 = np.column_stack(xy1)
    set2 = np.column_stack(xy2)
    if prj_dist is not None:
       set1[:,0], set1[:,1] = prjxy(prj_in, prj_dist, set1[:,0], set1[:,1])
       set2[:,0], set2[:,1] = prjxy(prj_in, prj_dist, set2[:,0], set2[:,1])

    f_dist = lambda i: np.any(np.sqrt((i[0]-set2[:,0])**2 + (i[1]-set2[:,1])**2)<=dist)
    set1_idx = np.array((list(map(f_dist, set1))))
    set1_new = set1[set1_idx]

    if plot==True:
        plt.scatter(set1[:,0],set1[:,1], c='g', label='set1')
        plt.scatter(set2[:,0],set2[:,1], c='b', label='set2')
        plt.scatter(set1_new[:,0],set1_new[:,1], marker='+', c='r',
                    label=f'''set1 points, {dist} distant from set2 points''')
        plt.legend()

    if prj_dist is not None:
        set1_new[:,0], set1_new[:,1] = prjxy(prj_dist, prj_in, set1_new[:,0], set1_new[:,1])

    return set1_idx, set1_new

# -----------------------------------------------------------------------------
def min_dist( x, y, prjcode_in=4326, prjcode_out=4326, 
              data_type='scattered_points', unique=True ) :
    """
    Find the minimum average distance among a set of scattered points.
    The result is a dictionary contatnjing also other statistics.
    """
    
    if prjcode_in != prjcode_out :
        x, y = prjxy( prjcode_in, prjcode_out, x, y )

    if data_type == 'grids' :
        dx, dy, mdxy = stepxy( x, y )
        min_dist_dict = {'dx':dx, 'dy':dy, 'mean':mdxy }        
        
    if  data_type == 'scattered_points' :   
        points_all = np.column_stack( ( x.ravel(), y.ravel() ) )
        if unique is True:
            points = np.unique( points_all, axis=0 )
        else:
            points = points_all
        N = points.shape[0]
        md = np.zeros( N )

        # for i in range( N ):
        #     apoints = np.delete( points, i, 0 )
        #     d = np.sqrt( ( points[(i, 0)] - apoints[:, 0] ) ** 2 + \
        #                  ( points[(i, 1)] - apoints[:, 1] ) ** 2 )
        #     md[ i ] = np.min( d )

        # --------------------------------------
        kdtree = sp.spatial.cKDTree( points )
        distances, _ = kdtree.query( points, k=2 )
        md = np.min( distances[:, 1:], axis=1 )
        # --------------------------------------

        meand = np.mean(md)
        mind = np.min(md)
        maxd = np.max(md)
        stdd = np.std(md)
        moded = sp.stats.mode(md, nan_policy='omit', keepdims=True).mode[0]
        
        # moded = sp_stat.mode(md, nan_policy='omit')[0][0]
        min_dist_dict = { 'mean':meand,'val':md, 'min':mind, 
                           'max':maxd, 'std':stdd, 'mode':moded,
                           'dist':md }
    
    return min_dist_dict

# -----------------------------------------------------------------------------
def xy2XY( x, y, step=None, lim=None, dist='mean', method='distance', extend=False,
           plot=False ) :
    """
    Convert scattered points (x, y) into a regular grid (X, Y) using interpolation.
    
    Parameters:
    - x, y: arrays of scattered points coordinates.
    - step: grid step size. If None, it is set to the mean distance between points.
    - lim: limits of the grid area. If None, it is determined based on the points.
    - dist: type of distance used to determine the step size. Default is 'mean'.
    - method: method used to extend the grid area. Default is 'distance'.
    - extend: flag to indicate whether to extend the grid area. Default is False.
    - plot: flag to indicate whether to plot the grid. Default is False.
    
    Returns:
    - X, Y: arrays of the regular grid coordinates.
    """
    
    if step is None :
        step = min_dist( x, y )['mean']
        
    if lim is None :
        lim = xy2lim( x, y, extend=extend, method=method )
        
    xm = np.arange( lim[0], lim[1]+step, step )  
    ym = np.arange( lim[3], lim[2]-step, -step ) 
    
    X, Y = np.meshgrid( xm, ym )
    
    if plot is True :
        plt.figure()
        Z = np.zeros( X.shape )
        plt.pcolor( X, Y, Z, edgecolor='k', vmin=0, vmax=1, cmap='Greys' )
        plt.scatter( x, y, c='r', s=0.5 )
        
# -----------------------------------------------------------------------------
def xyz2xy( xyz, xy, method='nearest', algebra='diff', rescale=False, fillnan=True,
            lim=None, plot=False ) :

    x, y, z = xyz
    
    nan = np.isnan( z )
    
    x1 = x[ ~nan ]
    y1 = y[ ~nan ]
    z1 = z[ ~nan ]
    
    if len( xy ) == 2 :
        x2, y2 = xy
    if len( xy ) > 2 : 
        x2, y2, z2 = xy
    
    if np.size( x1 ) < 4 :
        method = 'nearest'
    
    try:    
        zn = sp.interpolate.griddata( ( x1, y1 ), z1, ( x2, y2 ), method=method, 
                                      rescale=rescale )
    except : 
        zn = sp.interpolate.griddata( ( x1, y1 ), z1, ( x2, y2 ), method='nearest', 
                                      rescale=rescale )        
            
    if fillnan is True :
        
        if ( np.size( zn ) == 1 ) and ( np.any( np.isnan( zn ) ) ) :
            zn = sp.interpolate.griddata( ( x1, y1 ), z1, ( x2, y2 ), method='nearest' ) 
            
        if ( np.size( zn ) > 1 ) and ( np.any( np.isnan( zn ) ) ) :   
            idx = np.isnan( zn )
            zn[ idx ] = sp.interpolate.griddata( ( x1, y1 ), z1, ( x2[idx], y2[idx] ), 
                                                  method='nearest' )
            
            
    if len( xy ) > 2 : 
        if algebra == 'diff':
            za = z2 - zn
    else:
        za = None         

    if plot is True :
        plt.figure()
        sts = stat( zn )
        plt.scatter( xy[0], xy[1], c=zn, vmin=sts[2]-sts[3]*2, vmax=sts[2]+sts[3]*2, cmap='rainbow' )
        plt.colorbar()

    return zn, za

# -----------------------------------------------------------------------------
def find_by_att( dat1, dat2, col_xy, delta_xy, col_att, delta_att, condition=True,
                 prj_in=4326, prj_dist=None, unique=True, plot=False ) :

    dat1c = np.copy(dat1)
    dat2c = np.copy(dat2)
    
    if prj_dist is not None:
        dat1c[:,col_xy[0]], dat1c[:,col_xy[1]] = prjxy(prj_in, prj_dist, dat1[:,col_xy[0]], dat1[:,col_xy[1]])
        dat2c[:,col_xy[0]], dat2c[:,col_xy[1]] = prjxy(prj_in, prj_dist, dat2[:,col_xy[0]], dat2[:,col_xy[1]])
        
    mask1 = np.full(dat1c.shape[0], False)
    idx21 = []
    uniq = False
    col = col_xy + col_att
    delta = delta_xy + delta_att
    
    for n, i in enumerate(dat2c):
        mask0 = np.full(dat1c.shape[0], True)
        
        for c, d in zip(col, delta):
            msk = (dat1c[:,c] >= i[c] - d ) & ( dat1c[:,c] <= i[c] + d )
            mask0 = mask0 & msk
            
        if np.any(mask0):
            idx21.append([n, mask0])
            if ( np.sum(idx21[-1][1]) > 1 ) and ( unique is True ):
               uniq= uniq | True
               dist = np.sqrt( ( dat1c[mask0,col_xy[0]]-i[1] )**2  + 
                               ( dat1c[mask0,col_xy[1]]-i[0] )**2  )
               idx210N = np.nonzero( mask0*1 )[0]
               min_dist = np.nanmin( dist )
               unid = idx210N[ dist == min_dist ]                   
               idx21[-1][1] = unid[0] 
            if ( np.sum(idx21[-1][1] ) == 1 ): 
                idx21[-1][1] = np.where(mask0)[0].item()
            if ( np.sum(idx21[-1][1]) > 1 ) and ( unique is False ):
                uniq = uniq | True
                idx21[-1][1] = np.where(mask0).tolist()
        mask1 =  mask1 | mask0


    if unique is True:
        idx21 = np.asarray(idx21)

    if condition==False: mask1 = ~mask1
    dat1_msk = dat1[mask1]

    print( 'All data: ', len(dat1) )
    print( 'True: ', len(dat1_msk) )
    print( 'False: ', len(dat1)-len(dat1_msk) )

    if plot==True:
        plt.scatter(dat1[:,col_xy[0]], dat1[:,col_xy[1]], c='g', label='set1')
        plt.scatter(dat2[:,col_xy[0]], dat2[:,col_xy[1]], c='b', label='set2')
        plt.scatter(dat1_msk[:,col_xy[0]], dat1_msk[:,col_xy[1]], marker='+', c='r',
                    label='set1 points, defined by condition')
        plt.legend()

    return mask1, dat1_msk, idx21

# -----------------------------------------------------------------------------
def ell_radius( lat, radians=False ) :
    
    if radians is False :
        lat = np.copy( np.deg2rad( lat ) )
        
    num = (( a_wgs84**2 * np.cos(lat))**2) + (( c_wgs84**2 * np.sin(lat))**2 )   
    den = (( a_wgs84 * np.cos(lat))**2) + (( c_wgs84 * np.sin(lat))**2 )  
    
    R = np.sqrt( num / den )
    
    return R

# -----------------------------------------------------------------------------
def local_sph_raduis(lat):
    """
    Radius of the local sphere (based on latitude)

    Ref:
    - http://www2.ing.unipi.it/~a009220/lezioni/LI_ING_EDILE/AA1011/MATERIALE_DIDATTICO/APPUNTI/Geodesia.pdf

    """

    lat_r = np.radians( lat )
    
    e2 = ( a_wgs84**2 - c_wgs84**2 ) / a_wgs84**2
    
    N = a_wgs84 / np.sqrt( 1 - e2 * np.sin( lat_r )**2 )
    
    rho = a_wgs84 * ( 1 - e2 ) / np.sqrt( ( 1 - e2 * np.sin( lat_r )**2 )**3 )
    
    R = np.sqrt( rho * N )

    return R

# -----------------------------------------------------------------------------
def m2deg(md, n_digits=9, R=R_wgs84, lat=None):
    """
    Convert metric distances to decimal degrees (spherical approximation)

    lat = if not None, it calculate the distance using the Local Sphere Approximation,
            tangent to the wgs84 ellipsoid at the given latitude
    """

    if lat is None:   
        radi = R
    else:   
        radi = local_sph_raduis(lat)

    dd = np.round( md*360 / ( 2 * np.pi * radi ), n_digits)

    return dd

# -----------------------------------------------------------------------------
def deg2m(dd, lat=None, n_digits=3, R=R_wgs84):
    """
    Convert deg distances to meters (spherical approximation)

    lat = if not None, it calculate the distance using the Local Sphere Approximation,
            tangent to the wgs84 ellipsoid at the given latitude
    """

    if lat is None:   
        radi = R
    else:   
        radi = local_sph_raduis(lat)

    md = np.round( dd * 2 * np.pi * radi / 360, n_digits )

    return md

# -----------------------------------------------------------------------------
def pad_array( array, padw=25, mode='edge', alpha=None, constant_values=np.nan,
               plot=False, vmin=None, vmax=None, method='gdal', iter=1, 
               ptype='percentage', sqr_area=False ):
    
    ny, nx = array.shape
    
    if type( padw ) in ( int , float ) :
        padw = ( padw, padw ) 
    
    if ptype == 'percentage' :
        padw = [ np.int( ny * padw[0] / 100 ), np.int( nx * padw[1] / 100 ) ]
        
    if sqr_area == True:
        if ny > nx : padw[1] = padw[1] + ny - nx
        if ny < nx : padw[0] = padw[0] + nx - ny
            
    
    if mode in ('surface', 'gdal'): 
        pad_array = np.pad( array, pad_width=( ( padw[0], padw[0] ), ( padw[1], padw[1] ) ), 
                            mode='constant', constant_values=np.nan )
        pad_array = fillnan( pad_array, method=mode, iter=iter )      
    elif mode == 'constant' : 
        pad_array = np.pad( array, pad_width=((padw[0], padw[0]), (padw[1], padw[1])), 
                            mode=mode, constant_values=constant_values )        
    else: 
        pad_array = np.pad( array, pad_width=((padw[0], padw[0]), (padw[1], padw[1])), 
                            mode=mode )
    
    pnx, pny = pad_array.shape
    y0, y1, x0, x1 = (padw[0], pnx - padw[0], padw[1], pny - padw[1])
    original_shape_indx = (y0, y1, x0, x1)
        
    if alpha is not None:
        pad_array = taper( pad_array, alpha=alpha )[0]
        
    if plot == True:
        plta( array, sbplt=[1,2,1], vmin=vmin, vmax=vmax )
        plta( pad_array, sbplt=[1,2,2], vmin=vmin, vmax=vmax )
        
    return pad_array, original_shape_indx

# -----------------------------------------------------------------------------
def pad_xx_yy( pad_arr, xx, yy, plot=False ):
    
    dx = abs( np.mean(np.diff( xx, axis=1 ) ) )
    dy = abs( np.mean(np.diff( yy, axis=0 ) ) )
    
    a, c = pad_arr[0], pad_arr[1]
    xmin1, xmax1 = np.min(xx), np.max(xx)
    ymin1, ymax1 = np.min(yy), np.max(yy)
    xmin2, xmax2 = xmin1 - dx * c[2], xmax1 + dx * (a.shape[1] - c[3])
    ymax2, ymin2 = ymax1 + dy * c[0], ymin1 - dy * (a.shape[0] - c[1])
    xv = np.linspace(xmin2, xmax2, a.shape[1])
    yv = np.linspace(ymax2, ymin2, a.shape[0])
    xxn, yyn = np.meshgrid( xv, yv )
    
    if plot == True:
        plt.figure()
        plt.scatter((xxn.flatten()), (yyn.flatten()), c='b')
        plt.scatter((xx.flatten()), (yy.flatten()), c='r')
        
    return xxn, yyn
    
# -----------------------------------------------------------------------------
def taper( array, alpha=0.5, plot=False, vmin=None, vmax=None ):

    from scipy import signal
    
    nx, ny = array.shape[1], array.shape[0]
    tape_x = signal.tukey(nx, alpha=alpha)
    tape_y = signal.tukey(ny, alpha=alpha)
    
    t_xx, t_yy = np.meshgrid( tape_x, tape_y )
    taper_filt = t_xx * t_yy
    tarray = taper_filt * ( array - np.nanmean( array ) )
    tarray = tarray + np.nanmean( array )
          
    
    if plot == True:
        
        plta( array, sbplt=[1,3,1], vmin=vmin, vmax=vmax )
        plta( taper_filt, sbplt=[1,3,2], vmin=vmin, vmax=vmax )
        plta( tarray, sbplt=[1,3,3], vmin=vmin, vmax=vmax )
    
    return tarray, taper_filt        

# -----------------------------------------------------------------------------
def crop_pad(array, original_shape_indx, plot=False, vmin=None, vmax=None):
    
    array_crop = array[original_shape_indx[0]:original_shape_indx[1],
    original_shape_indx[2]:original_shape_indx[3]]
    
    if plot == True:
        
        plta( array, sbplt=[1,2,1], vmin=vmin, vmax=vmax )
        plta( array_crop, sbplt=[1,2,2], vmin=vmin, vmax=vmax )
        
    return array_crop

# -----------------------------------------------------------------------------
def xy2lim( x, y, prjcode_in=4326, prjcode_out=4326, extend=False, d=0,
            method='distance', sqr_area='False', plot=False ):

    if prjcode_in != prjcode_out:
        x,y = prjxy(prjcode_in, prjcode_out, x, y)

    lim = [ np.min(x), np.max(x), np.min(y), np.max(y) ]

    if extend is True:
        lim = extend_lim(lim, d, method, sqr_area )

    if plot is True:
        xplot = [lim[0], lim[0], lim[1], lim[1], lim[0]]
        yplot = [lim[2], lim[3], lim[3], lim[2], lim[2]]
        plt.scatter(x, y, marker='+')
        plt.plot(xplot, yplot, c='r')

    return lim

# -----------------------------------------------------------------------------
def absolute_file_paths( directory ):

    path = os.path.abspath(directory)
    
    return [entry.path for entry in os.scandir(path) if entry.is_file()]

# -----------------------------------------------------------------------------
def statistics(array, decimals=2):
    """
    Calculate the minimum, maximum, mean, and standard deviation of an array.
    Parameters:
    - array: numpy.ndarray
        The input array.
    - decimals: int, optional
        The number of decimal places to round the statistics to. Default is 2.
    Returns:
    - Min: float
        The minimum value of the array.
    - Max: float
        The maximum value of the array.
    - Mean: float
        The mean value of the array.
    - Std: float
        The standard deviation of the array.
    """
    Min = np.round( np.nanmin(array), decimals )
    Max = np.round( np.nanmax(array), decimals )
    Mean = np.round( np.nanmean(array), decimals )
    Std = np.round( np.nanstd(array), decimals )

    return Min, Max, Mean, Std

# -----------------------------------------------------------------------------
def plta( array, vmin=None, vmax=None, tit=None, lim=None, stat=6, sbplt=[],
          cmap='rainbow', axis=False, new_fig=True, contours=[], adjst_lim=True,
          flipud=False, hillshade=False, ve=2, aspect='auto', blend_mode='overlay',
          mask=None, points=None, pc='k', ps=1, label=None, label_size='large',
          xlabel=None, ylabel=None, x_ax=True, y_ax=True, letter=None, 
          xlett=0, ylett=0, colorbar=True, print_stat=True, alpha=1, lines=None,
          lc='k', lett_size='large', lett_colour='k', cc=None, out_lim=None,
          cl=True ) :
    
    """
    Plot array and print statistics
    """
    
    array = np.copy( array )
    
    if mask is not None :
        array[ mask ] = np.nan
        
    Min, Max, Mean, Std = np.nanmin(array), np.nanmax(array), \
                          np.nanmean(array), np.nanstd(array) 
    
    cmpa = plt.cm.get_cmap( cmap )
    
    if stat != None:
        if type(tit) == int : stat = tit
        Min, Max, Mean, Std = np.round(Min, stat), np.round(Max, stat), \
                              np.round(Mean, stat), np.round(Std, stat) 
                              
        if print_stat == True :
        
            print(f"Min:{Min} Max:{Max} Mean:{Mean} Std:{Std}")

    if sbplt != [] :
        r, c, n = sbplt
    else :
        r, c, n = 1, 1, 1
        
    if ( new_fig is True ) and ( n == 1 ) :
        plt.figure()  
        
    if sbplt != []:
        plt.subplot( r, c, n )
        
    if vmin == None:
        vmin = Mean - 2 * Std
    if vmax == None:
        vmax = Mean + 2 * Std
    
    if lim != None :
        if len( lim ) == 2 :
            lim = xy2lim( x=lim[0], y=lim[1] )
        dx = ( lim[1] - lim[0] ) / array.shape[0] 
        dy = ( lim[3] - lim[2] ) / array.shape[1] 
        
        if adjst_lim == True :
            lim[0], lim[2] = lim[0] - dx/2, lim[2] - dy/2 
            lim[1], lim[3] = lim[1] + dx/2, lim[3] + dy/2
    else :
        dx, dy = 1, 1
        
    if flipud == True :
        array = np.flipud( array )
    
    ax = plt.imshow( array, vmin=vmin, 
                     vmax=vmax, aspect=aspect, 
                     cmap=cmap, extent=lim,
                     alpha=alpha )
    plt.xlabel( xlabel, fontsize=label_size )
    plt.ylabel( ylabel, fontsize=label_size )
    
    if letter is not None :
        plt.annotate(letter, xy=(xlett,ylett), xycoords='axes fraction', size=lett_size,
                     c=lett_colour ) 
    
    ax = plt.gca()
    if colorbar is True :  
        cax = make_axes_locatable(ax).append_axes( 'right', size='2.5%', pad=0.05 ) 
        clb = plt.colorbar( cax=cax )
        
        if label is not None :
            clb.set_label( label, size=label_size )
                         
    if hillshade == True :
        ls = LightSource(azdeg=315, altdeg=45)
        rgb = ls.shade( array, cmap=cmpa, blend_mode=blend_mode, vert_exag=ve,
                        vmin=vmin, vmax=vmax, dx=1, dy=1 )
        alpha2 = np.isnan( array )
        rgb[ alpha2 ] = np.nan
        ax.imshow( rgb, extent=lim, alpha=alpha, aspect=aspect  )
    
    if points is not None :
        ax.scatter( points[0], points[1], c=pc, s=ps )

    if lines is not None :
        for line in lines :  
            ax.plot( line[0], line[1], c=lc)
        ax.set_xlim( lim[0], lim[1] )
        ax.set_ylim( lim[2], lim[3] )
    
    if contours != [] :
        cont = ax.contour( np.flipud(array), levels=contours, extent=lim, colors=cc,
                    linestyles='solid')  
        
        if cl is True :
            plt.clabel( cont )
        
    if axis is False:
        ax.axes.xaxis.set_visible(False)
        ax.axes.yaxis.set_visible(False)   

    if x_ax is False:
        ax.axes.xaxis.set_visible(False)
    if y_ax is False:
        ax.axes.yaxis.set_visible(False)        
        
    if type(tit) == str:
        ax.set_title(tit)
    if type(tit) == int:
        Min, Max, Mean, Std = np.round(Min, tit), np.round(Max, tit), \
                              np.round(Mean, tit), np.round(Std, tit)
        tit = f"[ Min={Min}  Max={Max}  Mean={Mean}  Std={Std} ]"
        ax.set_title(tit)
        
    if ( out_lim is not None ) and ( out_lim != [] ) :
        ax.set_xlim( [ out_lim[0], out_lim[1] ] )
        ax.set_ylim( [ out_lim[2], out_lim[3] ] )
        
    return ax

# -----------------------------------------------------------------------------
def map_profile( xy_start, xy_end, X, Y, Z=[], step=None, 
                 yaxis_sep=50, method='nearest',
                 colors=None,
                 legend_labels=None,
                 legend_label_size=None,
                 markers=None,
                 linestyle=None,
                 y_labels=None,
                 legend_ncol=1,
                 legend_loc='best',
                 x_label=None,
                 ax_label_size=None,
                 font='Tahoma',
                 y_ax_lim=None,
                 x_ax_lim=None,
                 plot=True,
                 subplot=(1, 1, 1),
                 text=None, text_size=14,
                 text_boxstyle='round',
                 text_facecolor='wheat',
                 text_alpha=0.5,
                 text_position=(0.95, 0.05),
                 xaxis_visible=True,
                 prjcode_in=4326,
                 smooth=1,
                 prjcode_out=None,
                 m2km=False,
                 deg2km=False,
                 prjcode_limits=None,
                 x_axis = True,
                 plot_map=False,
                 vmin=None,
                 vmax=None,
                 y_axis_type='multiple',
                 new_fig=True):

    from matplotlib import rcParams

    rcParams['font.family'] = font
    if prjcode_out is None:
        prjcode_out = prjcode_in

    if prjcode_limits is None:
        prjcode_limits = prjcode_in

    if prjcode_out != prjcode_in:
        X, Y = prjxy(prjcode_in, prjcode_out, X, Y)

    if prjcode_limits != prjcode_out:
        xy_start = prjxy(prjcode_limits, prjcode_out, xy_start[0], xy_start[1])
        xy_end = prjxy(prjcode_limits, prjcode_out, xy_end[0], xy_end[1])

    if step is None:
        step = min_dist(X, Y)['mean']
        print(step)
    x0, x1 = xy_start[0], xy_end[0]
    y0, y1 = xy_start[1], xy_end[1]
    theta = np.arctan((y1 - y0) / ((x1 - x0)+1e-24) )
    length = np.sqrt((y1 - y0) ** 2 + (x1 - x0) ** 2)
    lp = np.arange(0, length, step/smooth)
    xp = x0 + lp * np.cos(theta)
    yp = y0 + lp * np.sin(theta)
    
    if m2km is True :
        lp /= 1000
        
    if deg2km is True :
        lp = deg2m(lp)/1000  
        
    profiles = []
    for i, arr in enumerate(Z):
#        arr[np.isnan(arr)]=0
        lim = extend_lim( [ np.min((x0,x1)), np.max((x0,x1)), 
                            np.min((y0,y1)), np.max((y0,y1)) ], 2, sqr_area=True )
        xl, yl, idx = xy_in_lim( X, Y, lim )
        al = arr[idx]
        profiles.append(sp.interpolate.griddata((xl.ravel(), yl.ravel()), (al.ravel()),
                                                 (xp, yp), method=method, fill_value=np.nan))
        
    if plot is True:
        if new_fig is True :
            plt.figure()
        if colors is None:
            colors = ['#' + ''.join([random.choice('0123456789ABCDEF') for j in range(6)]) for i in range(len(profiles))]
        if legend_labels is None:
            legend_labels = ['function_' + str(i) for i in range(len(profiles))]
        if markers == 'random':
            markers_rand = itertools.cycle( ( ',', '+', '.', '*', '<', '>', 'h', 'p', 's',
                                              '.', 'H', 'D' ) )
        if linestyle is None :
            linestyle = ['solid' for i in range( len( profiles) ) ]
            
        lns = []
        ax_orig = plt.subplot(subplot[0], subplot[1], subplot[2])
        
        for i, (arr, color, lnsty) in enumerate(zip(profiles, colors, linestyle)):         
            if ( i == 0 ) or ( y_axis_type == 'unique') :
                ax = ax_orig
            else:
                ax = ax_orig.twinx()
                ax.spines['right'].set_position(('outward', yaxis_sep * (i - 1)))
            if markers == 'random' :
                ln = ax.plot(lp, arr, color=color, marker=(next(markers_rand)), linestyle=lnsty)
            else:
                ln = ax.plot(lp, arr, color=color, linestyle=lnsty )
            if y_axis_type == 'unique' :    
                ax.tick_params(axis='y', colors='k')
            else :
                ax.tick_params(axis='y', colors=color )
            if y_labels is not None:
                plt.ylabel(y_labels, fontsize=ax_label_size)
            if x_label is not None:
                plt.xlabel(x_label, fontsize=ax_label_size)
            lns += ln
            if y_ax_lim is not None :
                ax.set_ylim( y_ax_lim )      
            if x_ax_lim is not None :
                ax.set_xlim( x_ax_lim )                          

        if text is not None:
            props = dict(boxstyle=text_boxstyle, facecolor=text_facecolor, alpha=text_alpha)
            ax_orig.annotate(text, xy=text_position, xycoords='axes fraction',
              fontsize=text_size,
              bbox=props)
            plt.tight_layout()

        if legend_labels is not None:
            ax_orig.legend(lns, legend_labels, fontsize=legend_label_size, 
                                 ncol=legend_ncol, loc=legend_loc, framealpha=1)
            ax_orig.axes.get_xaxis().set_visible(xaxis_visible)
    
        if x_axis is False :
            plt.gca().axes.xaxis.set_ticklabels([])
            plt.gca().xaxis.set_major_locator(plt.NullLocator())
            
    if plot_map == True:
        plt.figure()
        if vmin == None:
            vmin = np.nanmean(Z[0]) - 2 * np.nanstd(Z[0])
        if vmax == None:
            vmax = np.nanmean(Z[0]) + 2 * np.nanstd(Z[0])
        plt.imshow(Z[0], vmin=vmin, vmax=vmax, cmap='rainbow',
                   extent=(X.min(), X.max(), Y.min(), Y.max()))
        plt.plot(xp, yp, color='k')            
    
    pr_list = []
    for i in profiles:
        pr_list.append(np.column_stack((xp, yp, lp, i)))

#    plt.tight_layout()

    return pr_list

# -----------------------------------------------------------------------------
def XYZ_crop2lim( X, Y, Z, lim ) :

    if len( Z.shape ) == 2:

        xi = np.where( ( X[0, :] >= lim[0]) & (X[0, :] <= lim[1] ) )
        yi = np.where( ( Y[:, 0] >= lim[2]) & (Y[:, 0] <= lim[3] ) )

        Xc = X[ np.min(yi):np.max(yi), np.min(xi):np.max(xi) ]
        Yc = Y[ np.min(yi):np.max(yi), np.min(xi):np.max(xi) ]
        Zc = Z[ np.min(yi):np.max(yi), np.min(xi):np.max(xi) ]
        

        return [ Xc, Yc, Zc ]

    if len( Z.shape ) == 1:

        xi = np.where( ( X >= lim[0] ) & ( X <= lim[1] ) )
        yi = np.where( ( Y >= lim[2] ) & ( Y <= lim[3] ) )

        xc = X[ ( xi, yi ) ]
        yc = Y[ ( xi, yi ) ]
        zc = Z[ ( xi, yi ) ]

        return [ xc, yc, zc ]

# -----------------------------------------------------------------------------
def isiterable(p_object):
    """
    Check if an object is iterable
    """

    try:
        it = iter(p_object)
    except TypeError:
        return False
    return True

# -----------------------------------------------------------------------------
def stepxy( xarray, yarray ) :
    
    " Distance between grid data points "
    
    cx = int( xarray.shape[1]/2 )
    cy = int( yarray.shape[0]/2 )
    
    arx = xarray[ cy:cy+2, cx:cx+2 ]
    ary = yarray[ cy:cy+2, cx:cx+2 ]
    
    dx = np.sqrt( ( arx[0,1] - arx[0,0] )**2 + ( ary[0,1] - ary[0,0] )**2 )
    dy = np.sqrt( ( arx[0,0] - arx[1,0] )**2 + ( ary[0,0] - ary[1,0] )**2 )
    
    dmxy = np.mean( ( dx, dy ) )
        
    return dx, dy, dmxy

# -----------------------------------------------------------------------------
def grid_fill( data, invalid=None ):
    """
    Replace the value of invalid 'data' cells (indicated by 'invalid') 
    by the value of the nearest valid data cell

    Input:
        data:    numpy array of any dimension
        invalid: a binary array of same shape as 'data'. True cells set where data
                 value should be replaced.
                 If None (default), use: invalid  = np.isnan(data)

    Output: 
        Return a filled array. 
    """

    if invalid is None: 
        invalid = np.isnan(data)

    ind = sp.ndimage.distance_transform_edt(invalid, return_distances=False, return_indices=True)
    
    return data[ tuple(ind) ]

# -----------------------------------------------------------------------------
def mask_edges( mask, plot=False, c='k', s=None, plt_array=None, vmin=None, vmax=None ) :
    
    """
    Returne x-y coordinates of the edges of a "mask-array" ( boolean or 0/1 array )
    
    0 0 0 0 0 0 0     0 0 0 0 0 0 0     ( (1,1)
    0 1 1 1 1 1 0     0 1 1 1 1 1 0       (1,2)
    0 1 1 1 1 1 0     0 1 0 0 0 1 0       (1,3) 
    0 1 1 1 1 1 0 ==> 0 1 0 0 0 1 0 ==>   (1,4)
    0 1 1 1 1 1 0     0 1 0 0 0 1 0       (1,5)
    0 1 1 1 1 1 0     0 1 1 1 1 1 0       (2,1)
    0 0 0 0 0 0 0     0 0 0 0 0 0 0       (2,5)... )
    """
    
    if type( mask ) == bool :
        mask = mask*1
    
    fil = [[-1,-1,-1],
           [-1, 8,-1],
           [-1,-1,-1]]
    
    edges = np.where( sp.ndimage.convolve(mask*1, fil, mode='constant') > 1)
    
    x_ind, y_ind = edges[1], edges[0] 
        
    if plot is True :
        if plt_array is not None :
            plt.imshow( plt_array, vmin=vmin, vmax=vmax, cmap='rainbow')
            plt.colorbar()
        plt.scatter( x_ind, y_ind, c=c, s=s )
        
    return x_ind, y_ind

# -----------------------------------------------------------------------------
def normalize( array, b=1, a=0 ) :
    
    Min = np.nanmin(  array )
    Max = np.nanmax( array )
    
    norm = ( b - a ) * ( ( array - Min ) / ( Max - Min ) ) + a
        
    return norm  

# -----------------------------------------------------------------------------
def fillnan( array, xy=None, method='nearest', size=3, iter=1, maxSearchDist=None, 
             plot=False, vmin=None, vmax=None, edges=False, smooth=0, tension=0.2 ) :
    
    zz = np.copy( array )
    
    if np.all( np.isfinite( zz ) ):
        zzfn = zz
    
    if ( array.shape[1] > 1 ) and ( xy is None ):
        xi, yi = np.arange(0, zz.shape[1], 1), np.arange(0, zz.shape[0], 1)
        xx, yy = np.meshgrid( xi, yi )
        
    if xy is not None:
        xx, yy = xy[0], xy[1]
        
    zz_nan = np.isnan( zz )  
    
    x, y, z = xx.flatten(), yy.flatten(), zz.flatten()
    nan = np.isnan(z)
    notnan = np.invert(nan)
    xn, yn = x[nan], y[nan]
    xi, yi, zi = x[notnan], y[notnan], z[notnan]
    
    if method == 'nearest':
        zfn = sp.interpolate.griddata( ( xi, yi), zi, (xn, yn), method='nearest')
        zn = np.copy(z)
        zn[nan] = zfn
        zzfn = zn.reshape( zz.shape )
        
    if method == 'nearest2D' :
        zzfn = grid_fill( array )
        
    if method == 'mean':
        zfn = np.nanmean( z )
        z[ np.isnan( z ) ] = zfn
        zzfn = z.reshape( zz.shape )
        
    if type( method ) in ('float', 'int'):
        zfn = method
        z[ np.isnan( z ) ] = zfn
        zzfn = z.reshape( np.shape( zz ) )
        
    if method == 'gdal':
        zz[ np.isnan( zz ) ] = 1e-24
        rx, ry, _ = stepxy( xx, yy )
        driver = gdal.GetDriverByName('GTiff')
        raster = driver.Create( '/vsimem/raster_filled.vrt', zz.shape[1], zz.shape[0],
                                 bands = 1, eType=gdal.GDT_Float32 )
    
        raster.SetGeoTransform( [ xx.min(), rx, 0.0, yy.max(), 0.0, -ry ] )
        raster.GetRasterBand( 1 ).SetNoDataValue( 1e-24 )
        
        out_band = raster.GetRasterBand( 1 )
        out_band.WriteArray( zz )        
        
        if maxSearchDist is None: 
            maxSearchDist = int( np.max( zz.shape ) )
            
        gdal.FillNodata( targetBand = raster.GetRasterBand(1), maskBand=None, 
                         maxSearchDist=maxSearchDist,
                         smoothingIterations=iter )
        
        zzfn = raster.GetRasterBand( 1 ).ReadAsArray()
        zzfn[ zzfn == 1e-24 ] = np.nan
        raster = None
        
    if method in ( 'gauss', 'uniform' ) :
        zn = np.copy( zz ) - np.nanmean( zz )
        zn = grid_fill( zn, invalid=None )
        if method == 'gauss' :
            func = sp.ndimage.gaussian_filter
        if method == 'uniform' :
            func = sp.ndimage.uniform_filter
        
        for i in range( iter ) :
            FG = func( zn, size )
            WF = func( ~zz_nan*1.0, size )
            DFG =  zn - FG
            zn = FG + DFG * WF 
        
        zzfn = zn + np.nanmean( zz )

    if plot == True:
        plta( array, sbplt=[1, 3, 1], tit='original', vmin=vmin, vmax=vmax )
        plta( zzfn, sbplt=[1, 3, 2], tit='fill', vmin=vmin, vmax=vmax )
        plta( array - zzfn, sbplt=[1, 3, 3], tit='differences' )
        
    return zzfn

# -----------------------------------------------------------------------------
def block_m( x, y, z, wind_size, method='mean', data_type='vector', lim=None,  
             prjcode_in=4326, prjcode_out=4326, nan=False, adjst_lim=True, 
             plot=False, s1=None, s2=None ) :         
    """
    Block Mean or Block Median
    """

    if prjcode_in != prjcode_out:
        x, y = prjxy(prjcode_in, prjcode_out, x, y)
    
    if lim == None:
        lim = [np.min(x), np.max(x), np.min(y), np.max(y)]
        
#    x_blk = np.linspace(lim[0], lim[1], int( abs( lim[0]-lim[1] ) / wind_size ) )
#    y_blk = np.linspace(lim[2], lim[3], int( abs( lim[2]-lim[3] ) / wind_size ) )          
         
    x_blk = np.arange(lim[0], lim[1], wind_size)
    y_blk = np.arange(lim[2], lim[3], wind_size) 

    Xb, Yb = np.meshgrid(x_blk, y_blk)
    xb, yb, zb = np.hsplit( np.zeros( ( (Xb.shape[0]-1) * (Xb.shape[1]-1), 3 ) ), 3 )
    
    n = 0
    for idx, _ in np.ndenumerate( Xb ) :
     
        i, j = idx
        if ( i == Xb.shape[0]-1 ) or ( j == Xb.shape[1]-1 ) :
            continue
        
        win_idx = ( ( x > Xb[ ( i, j ) ] ) & ( x < Xb[ ( i, j + 1 ) ]) & \
                    ( y > Yb[ ( i, j ) ] ) & ( y < Yb[ ( i + 1, j ) ] ) )
        
        if method == 'mean':              
            xb[n] = np.mean( x[ win_idx ] )
            yb[n] = np.mean( y[ win_idx ] )
            zb[n] = np.mean( z[ win_idx ] )
            
        if method == 'median':              
            xb[n] = np.median( x[ win_idx ] )
            yb[n] = np.median( y[ win_idx ] )
            zb[n] = np.median( z[ win_idx ] )
            
        n += 1    

    if data_type == 'grid' :
        if adjst_lim == True :
            x_grid = np.arange(np.min(xb), np.max(xb)+wind_size, wind_size)
            y_grid = np.arange(np.min(yb), np.max(yb)+wind_size, wind_size)
        if adjst_lim == False :
            x_grid = np.linspace(np.min(xb), np.max(xb), int( ( np.max(xb) - np.min(xb) ) / wind_size ) )
            y_grid = np.linspace(np.min(yb), np.max(yb), int( ( np.max(yb) - np.min(yb) ) / wind_size ) )
        Xg, Yg = np.meshgrid( x_grid, y_grid )
        Zg = xyz2xy( ( xb.ravel(), yb.ravel(), zb.ravel() ), ( Xg, Yg ), method='linear' )[0]
        xb, yb, zb = Xg, Yg, Zg

        if ( y[0,0] - y[-1,0] ) * ( yb[0,0] - yb[-1,0] ) < 0 :
            yb = np.flipud( yb )
            zb = np.flipud( zb )
        
    if data_type == 'vector' :
        if nan == False :
            xb = xb[ ( ~np.isnan( zb ) ) ]
            yb = yb[ ( ~np.isnan( zb ) ) ]
            zb = zb[ ( ~np.isnan( zb ) ) ]
    
    if plot == True:
        plt.figure()
        plt.scatter(x, y, c='b', s=s1)
        plt.scatter(xb, yb, c='r', s=s2)  
        
    return xb, yb, zb

# -----------------------------------------------------------------------------
def array2csv( array, headers=[], sep=',', fmt='% 15.6f', 
               path_name='new_dataset', nan=None, nan_c=2 ) :

    path = os.path.dirname( path_name )
    if path != '' :
        os.makedirs( path, exist_ok=True )
    
    if nan != None :
        inan = np.isnan( array[ :, nan_c] )
        if ( type( nan ) == bool ) and ( nan == False ) :
            array = array[ ~inan ]
        if type( nan ) in ( float, int ) :
            array[ :, nan_c ] = nan
            
    hd = ''
    fmt_str = ''
    n_col = array.shape[1]
    
    if type( fmt ) == str :
        fmt = [fmt]
        fmt = [ fmt[0] for i in range( n_col ) ]
        
    for i, _ in enumerate( fmt ):
        fmt_str = fmt_str + fmt[i] + sep
    fmt_str = fmt_str[:-1]
        
    if headers != []:        
        for i, h in enumerate( headers ):
            space = fmt[i]
            for ch in ['% ', 'f', 'd', 'i', 'e', '> ' ] :
                space = space.replace(ch, '.')    
            space = space.split('.')[1]
            hd = hd + f'{h:>{space}}' + sep
        hd = hd[:-1]   
     
    np.savetxt( path_name, array, header=hd, fmt=fmt_str, 
                comments='', delimiter=sep )
    
    abs_path = os.path.abspath(path_name)

    return abs_path

# -----------------------------------------------------------------------------
def read_csv( csv_file, sep=' ', header=None, skiprows=[], skipcols=[], 
              encoding="utf8", n_col=None, adjust_last_col_fmt=True, comments=[],
              col_size=[], new_sep=',', printf=False, decimals=2, space=12,
              rows=[], cols=[], show=False, nan=np.nan ) :
    
    dictionary = {}
    
    f = open( csv_file, "r", encoding=encoding, errors='ignore' )
    
    lines = f.readlines()
    
    if col_size != [] :

        for i, l in enumerate( lines ) :
            
            n = 0
            c = 0
            for sz in col_size :
                
                c += sz 
                
                l = l[:c+n] + new_sep + l[c+n:]
                n += 1
                lines[ i ] = l
                
        n_col = n
        sep = new_sep
    
    if printf is True :
        
        if rows == []:
            
            rows = range( len(lines) )
            
        for i in rows :
            
            print( lines[i], end="" )
            
    data = None
    hd = None
    
    if comments != [] :
        
        if skiprows == [] :
            
            skiprows = comments  
            
        comments = [ lines[i] for i in comments ]
        
    if (type(header) == int) and (header > 0) :
        
        skiprows = skiprows + [n for n in range(header)]
        
    if type( skiprows ) == int :
        
        skiprows = [n for n in range(skiprows)]
        
    for i, l in enumerate( lines ) :   
        
        if ( i in skiprows ) :
            
            continue

        if sep not in l :
            
            continue
        
        if skipcols != [] :
            
            for i, li in enumerate(l) :
                
                if i in skipcols :
                    
                    l = l[:i] + sep + l[(i+1):]
                    
        if i == header :
            
            hd = [ h.strip() for h in l.split( sep ) if ( (h!='') and (h.strip()!='') ) ]
            
            continue
        
        if hd != None :
            
            n_col = len( hd )
            
        if ( l != "\n" ) and ( data == None ) and ( n_col == None ) :
            
            n_col = len( [ n.strip() for n in l.split( sep ) if ( (n!='') and (n.strip()!='') ) ] )
            
        if ( n_col != None ) and ( data == None ) :
            
            data = [ [] for i in range( n_col ) ]
        
        if ( l != "\n" ) and ( data != None ) :
            
            h_fmt = l.split( sep )
            
            if col_size != [] :
                
                line_list = [ n.strip() for n in l.split( sep ) ]
            
            else :
                
                line_list = [ n.strip() for n in l.split( sep ) if ( (n!='') and (n.strip()!='') ) ]
                    
            for c in range( n_col ) :
                
                val = line_list[c]
                
                if ( val == '' ) or ( val.strip() == '' ) :
                    val = nan     
                data[c].append( val )  
    f.close() 
    
    
    if hd == None :
        
        hd = [ 'c' + str(n+1) for n in range( n_col) ]
        
    for i, h in enumerate( hd ) :
        
        try:
            
            dictionary[h] = np.asarray( data[i], dtype=float )
            
        except ValueError:
            
            dictionary[h] = np.array( data[i] )    
            
    array = dict2array( dictionary )[0]   
    
    fmt = []
    hfmt = []
    beg = 0
    for i, j in enumerate(h_fmt) :
        
        if i != 0 :
            
            if ( i!= len(h_fmt)-1 ) and ( j == '' ) and ( h_fmt[i-1] != '' ) :
                
                hfmt.append( ' '.join(h_fmt[ beg:i ]) )
                
                beg=i
                
            if ( i == len(h_fmt)-1 ) and ( beg!= 0 ) :
                
                hfmt.append( ' '.join(h_fmt[ beg: ] ) )
                
    if hfmt == [] :
        
        hfmt = h_fmt            
                    
    for i, k in enumerate( hfmt ) :  
        
        if ( i == len( hfmt ) - 1 ) and ( adjust_last_col_fmt == True ) : 
            len_l = str( len( k ) -1 )
        
        else :    
            len_l = str( len( k ) )
        string_str = k.strip()
        
        try:
            _ = int( string_str )
            fmt.append('% '+len_l+'d')
        except:
            
            try :
                _ = float( string_str )
                ndec = str( len( string_str.split( '.' )[1] ) )
                fmt.append('% '+len_l+'.'+ndec+'f')
            
            except :
                fmt.append('% '+len_l+'s')  

    if header == None:
        hd = None 
        
    if show is True :
        _ = print_table( dictionary, 
                         decimals=decimals, 
                         space=space, 
                         rows=rows,
                         cols=cols )
        
    return dictionary, array, hd, fmt, comments, data     

# -----------------------------------------------------------------------------
def csv2dict( csv_file, sep=',', header=None, skiprows=[], encoding="utf8", 
              n_col=None, adjust_last_col_fmt=True ) :
    
    d = read_csv( csv_file, sep=sep, header=header, skiprows=skiprows, 
                  encoding=encoding, n_col=n_col, 
                  adjust_last_col_fmt=adjust_last_col_fmt )[0]

    return d

# -----------------------------------------------------------------------------
def csv2array( csv_file, sep=',', header=None, skiprows=[], encoding="utf8", 
               n_col=None, adjust_last_col_fmt=True ) :
    
    a = read_csv( csv_file, sep=sep, header=header, skiprows=skiprows, 
                  encoding=encoding, n_col=n_col, 
                  adjust_last_col_fmt=adjust_last_col_fmt )[1]

    return a
   
# ----------------------------------------------------------------------------- 
def dict2array( dictionary, exclude=[] ):
    """
    Convert a dictionary of arrays into a 2D numpy array.

    Args:
        dictionary (dict): The dictionary containing arrays.
        show (bool, optional): Whether to display the intermediate steps. Defaults to False.
        exclude (list, optional): List of keys to exclude from the conversion. Defaults to [].

    Returns:
        tuple: A tuple containing the converted 2D numpy array and the headers.

    """
    headers = []

    # Find the maximum length among all arrays in the dictionary
    max_length = max(len(dictionary[k]) for k in dictionary if k not in exclude)

    for i, k in enumerate( dictionary ):

        if k in exclude:
            continue

        # If the current array is shorter than the maximum length, extend it with None values
        if len(dictionary[k]) < max_length:
            dictionary[k].extend([np.nan] * (max_length - len(dictionary[k])))

        # Try to convert to float, if not possible, fall back to object dtype
        try:
            dictionary[k] = np.array(dictionary[k], dtype=float)
        except ValueError:
            dictionary[k] = np.array(dictionary[k], dtype=object)

        if i == 0:
            array = np.copy(dictionary[k])
        else:
            array = np.column_stack((array, dictionary[k]))

        headers.append(k)

    return array, headers

# -----------------------------------------------------------------------------            
def dict2csv( dictionary, sep=',', 
              fmt='% 15.6f', 
              path_name='new_dataset.csv' ) :
    
    array, headers = dict2array( dictionary )

    abs_path = array2csv( array, headers=headers, 
                          sep=sep, fmt=fmt, 
                          path_name=path_name )

    return abs_path

# -----------------------------------------------------------------------------            
def addc2csv( csv_file, new_c=[], new_hd=[], new_fmt=[], sep=',', rep_fmt=None,
              path_name=None, rep_headers=None, header=None, skiprows=[], 
              encoding="utf8", n_col=None, adjust_last_col_fmt=True ) :

    if path_name == None :
        path_name = csv_file

    csvf = read_csv( csv_file, sep=sep, header=header, skiprows=skiprows, 
                     encoding=encoding, n_col=n_col, 
                     adjust_last_col_fmt=adjust_last_col_fmt )

    array = csvf[1]
    headers = csvf[2]
    fmt = csvf[3]
    n_c = array.shape[1]
    new_array = np.copy( array )

    if new_hd == [] :
       new_hd = [ 'c_'+str(n_c+i) for i in range( len( new_c ) ) ] 
    if new_fmt == [] :
       new_fmt = [ '% 15.6f' for i in range( len( new_c ) ) ]   
    if type( new_hd ) == str :
        new_hd = [ new_hd ] * len( new_c )    
    if type( new_fmt ) == str :
        new_fmt = [ new_fmt ] * len( new_c )     
   
    for i, col in enumerate( new_c ) :
        new_array = np.column_stack( ( new_array, col ) )
        if headers != None :
            headers.append( new_hd[i] )
        fmt.append( new_fmt[i] )

    if rep_fmt != None :
        fmt = rep_fmt
        
    if rep_headers != None :
        headers = rep_headers

    abs_path = array2csv( new_array, headers=headers, sep=sep, fmt=fmt, 
                          path_name=path_name )
    
    return abs_path
   
# ----------------------------------------------------------------------------- 
class SelectFromCollection(object):
    # Select indices from a matplotlib collection using `LassoSelector`.
    #
    # Selected indices are saved in the `ind` attribute. This tool fades out the
    # points that are not part of the selection (i.e., reduces their alpha
    # values). If your collection has alpha < 1, this tool will permanently
    # alter the alpha values.
    #
    # Note that this tool selects collection objects based on their *origins*
    # (i.e., `offsets`).
    #
    # Parameters
    # ----------
    # ax : :class:`~matplotlib.axes.Axes`
    #     Axes to interact with.
    #
    # collection : :class:`matplotlib.collections.Collection` subclass
    #     Collection you want to select from.
    #
    # alpha_other : 0 <= float <= 1
    #     To highlight a selection, this tool sets all selected points to an
    #     alpha value of 1 and non-selected points to `alpha_other`.
    #

    def __init__(self, ax, collection, alpha_other=0.3, facecolors=None):
        self.canvas = ax.figure.canvas
        self.collection = collection
        self.alpha_other = alpha_other

        self.xys = collection.get_offsets()
        self.Npts = len(self.xys)

        # Ensure that we have separate colors for each object
        self.fc = collection.get_facecolors()
        if len(self.fc) == 0:
            raise ValueError('Collection must have a facecolor')
        elif len(self.fc) == 1:
            self.fc = np.tile(self.fc, (self.Npts, 1))

        if facecolors is not None: self.fc = facecolors

        self.lasso = LassoSelector(ax, onselect=self.onselect)
        self.ind = []

    def onselect(self, verts):
        path = Path(verts)
        self.ind = np.nonzero(path.contains_points(self.xys))[0]
        self.fc[:, -1] = self.alpha_other
        self.fc[self.ind, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()

    def disconnect(self):
        self.lasso.disconnect_events()
        self.fc[:, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()

# -----------------------------------------------------------------------------
def rolling_win_2d( arr, win_shape ) :

    if ( win_shape[0] % 2 == 0 ) or ( win_shape[1] % 2 == 0 ) :
        raise NameError('Both values in win_shape must be odd integers !')

    r_extra = np.floor(win_shape[0] / 2).astype(int)
    c_extra = np.floor(win_shape[1] / 2).astype(int)
    a = np.empty((arr.shape[0] + 2 * r_extra, arr.shape[1] + 2 * c_extra))
    a[:] = np.nan
    a[r_extra:-r_extra, c_extra:-c_extra] = arr    

    s = (a.shape[0] - win_shape[0] + 1,) + (a.shape[1] - win_shape[1] + 1,) + win_shape
    strides = a.strides + a.strides

    windows = np.lib.stride_tricks.as_strided(a, shape=s, strides=strides)
    windows = windows.reshape( arr.size, win_shape[0] * win_shape[1] )

    return windows

# -----------------------------------------------------------------------------
def median_filt( array, radius=1, padw=0, pmode='linear_ramp', 
                 plot=False, vmin=None, vmax=None, iter=1 ):

    ar_pad, original_shape_indx = pad_array(array, padw, pmode)
    fw = np.ones( ( radius * 2 + 1, radius * 2 + 1 ) )

    for i in range( iter ) :
        ar_pad = sp.ndimage.median_filter( ar_pad, footprint=fw,
                                               mode='nearest' )

    ar_filt = crop_pad( ar_pad, original_shape_indx )

    if plot == True:

        plta( array, sbplt=[1, 3, 1], tit='Original')
        plta( ar_filt, vmin, vmax, sbplt=[1, 3, 2], tit='Filtered')
        plta( array - ar_filt, vmin, vmax, sbplt=[1, 3, 3], tit='Differences')
        plt.tight_layout() 

    return ar_filt

# -----------------------------------------------------------------------------
def mean_filt( array, radius=1, padw=0, pmode='linear_ramp', 
               plot=False, vmin=None, vmax=None, iter=1 ) :
    
    from scipy import signal

    ar_pad, original_shape_indx = pad_array(array, padw, pmode)
    fw = np.ones( ( radius*2 + 1, radius*2 + 1) ) / ( radius * 2 + 1 ) ** 1

    for i in range( iter ) : 
        ar_pad = signal.convolve2d( ar_pad, fw, mode='same' )

    ar_filt = crop_pad( ar_pad, original_shape_indx )

    if plot == True:

        plta( array, sbplt=[1, 3, 1], tit='Original')
        plta( ar_filt, vmin, vmax, sbplt=[1, 3, 2], tit='Filtered')
        plta( array - ar_filt, vmin, vmax, sbplt=[1, 3, 3], tit='Differences')
        plt.tight_layout() 

    return ar_filt

#------------------------------------------------------------------------------
def hanning_filt( array, padw=0, pmode='linear_ramp', 
                  plot=False, vmin=None, vmax=None, iter=1 ) :

    from scipy import signal

    ar_pad, original_shape_indx = pad_array(array, padw, pmode)
    fw = np.array( [ [ 0, 1, 0 ], [ 1, 2, 1 ], [ 0, 1, 0 ] ] ) / 6

    for i in range( iter ) : 
        ar_pad = signal.convolve2d( ar_pad, fw, mode='same' )
    
    ar_filt = crop_pad( ar_pad, original_shape_indx )
    
    if plot == True:
        
        plta( array, sbplt=[1, 3, 1], tit='Original')
        plta( ar_filt, vmin, vmax, sbplt=[1, 3, 2], tit='Filtered')
        plta( array - ar_filt, vmin, vmax, sbplt=[1, 3, 3], tit='Differences')
        plt.tight_layout() 
        
    return ar_filt


# -----------------------------------------------------------------------------
def gauss_filt( array, radius=1, sigma=1, padw=0, pmode='linear_ramp', 
                alpha=None, iter=1, plot=False, vmin=None, vmax=None ):

    ar_pad, original_shape_indx = pad_array(array, padw=padw, mode=pmode, alpha=alpha)

    x, y = np.mgrid[-radius:radius + 1, -radius:radius + 1]
#    normal = 1 / (2.0 * np.pi * sigma ** 2)
#    fw = np.exp(-((x ** 2 + y ** 2) / (2.0 * sigma ** 2))) * normal
    if radius is not None :
        w = radius * 2 + 1
        truncate = ( ( ( w - 1 ) / 2 ) - 0.5 ) / sigma
    else :
        truncate = 4

    for i in range( iter ) : 
#        ar_pad = signal.convolve2d(ar_pad, fw, mode='same')
        ar_pad = sp.ndimage.gaussian_filter(ar_pad, sigma, truncate=truncate)

    ar_filt = crop_pad( ar_pad, original_shape_indx )

    if plot == True:

        plta( array, sbplt=[1, 3, 1], tit='Original')
        plta( ar_filt, vmin, vmax, sbplt=[1, 3, 2], tit='Filtered')
        plta( array - ar_filt, vmin, vmax, sbplt=[1, 3, 3], tit='Differences')
        plt.tight_layout() 

    return ar_filt

# -----------------------------------------------------------------------------
def resampling_filt( array, factor, padw=0, pmode='gdal', 
                     alpha=None, spl_order=2, plot=False, vmin=None, vmax=None ) :

    if padw == 0 :
        padw = [ 0, 0 ]
    if array.shape[0] % 2 == 0  :
        padw[0] += 1 
    if array.shape[1] % 2 == 0  :  
        padw[1] += 1 

    ar_pad, original_shape_indx = pad_array( array, padw=padw, mode=pmode, alpha=alpha, ptype='' )
    
    # Down_Sampling
    ar_pad_dw = resampling( ar_pad, 1/factor, spl_order=spl_order )
    
    # Upsampling
    ar_pad_up = resampling( ar_pad_dw, factor, spl_order=spl_order )   
    
    ar_filt = crop_pad( ar_pad_up, original_shape_indx )
    
    if plot == True:
        
        plta( array, sbplt=[1, 3, 1], tit='Original')
        plta( ar_filt, vmin, vmax, sbplt=[1, 3, 2], tit='Filtered')
        plta( array - ar_filt, vmin, vmax, sbplt=[1, 3, 3], tit='Differences')
        plt.tight_layout() 
        
    return ar_filt    

# -----------------------------------------------------------------------------
def std_filt( array, padw=0, pmode='gdal', radius=1, n=1, 
              alpha=None, plot=False, vmin=None, vmax=None ) :
    
    ar_pad, original_shape_indx = pad_array(array, padw=padw, mode=pmode, alpha=alpha)
    
    win_shape = radius * 2 + 1, radius * 2 + 1 
    windows = rolling_win_2d( ar_pad, win_shape )
    w_mean = np.nanmean( windows, axis=1 )
    w_std = np.nanstd( windows, axis=1 )
    w_filt = w_mean + np.sign( w_mean ) * w_std
    ar_pad_filt = w_filt.reshape( ar_pad.shape )
    
    ar_filt = crop_pad( ar_pad_filt, original_shape_indx )    
    
    if plot == True:
        
        plta( array, sbplt=[1, 3, 1], tit='Original')
        plta( ar_filt, vmin, vmax, sbplt=[1, 3, 2], tit='Filtered')
        plta( array - ar_filt, vmin, vmax, sbplt=[1, 3, 3], tit='Differences')
        plt.tight_layout() 
        
    return ar_filt        

#------------------------------------------------------------------------------
def filt2d( array, radius=1, padw=0, pmode='linear_ramp', plot=False, vmin=None, 
            vmax=None, iter=1, ftype='mean', sigma=1, factor=2, mask=None,
            fill=None ):

    if type( ftype ) is str :
        ftype = [ ftype ]
    
    if  type( iter ) is int :
        iter = [ iter for ft in ftype ]
        
    if fill != None :
        nan = np.isnan( array )
        array = fillnan( array, method=fill )
    
    ar_filt = np.copy( array )
    for i, f in enumerate( ftype ) :
        
        if ftype[i] == 'mean' :
            ar_filt = mean_filt( ar_filt, radius=radius, padw=padw, pmode=pmode, 
                                 iter=iter[i] )  
        if ftype[i] == 'hanning' :
            ar_filt = hanning_filt( ar_filt, padw=padw, pmode=pmode, iter=iter[i] )
            
        if ftype[i] == 'median' :
            ar_filt = median_filt( ar_filt, radius=radius, padw=padw, pmode=pmode, 
                                   iter=iter[i] )      
        if ftype[i] == 'gauss' :
            ar_filt = gauss_filt( ar_filt, radius=radius, padw=padw, pmode=pmode, 
                                   iter=iter[i], sigma=sigma ) 
        if ftype[i] == 'resamplig' :
            ar_filt = resampling_filt( ar_filt, factor=factor, padw=padw, pmode=pmode ) 
            
    if fill != None :
          ar_filt[ nan ] = np.nan   
            
    if mask is not None:
        ar_filt[ mask ] = np.nan
        array = np.copy(array)
        array[ mask ] = np.nan            

    if plot == True:
        
        plta( array, sbplt=[1, 3, 1], tit='Original')
        plta( ar_filt, vmin, vmax, sbplt=[1, 3, 2], tit='Filtered')
        plta( array - ar_filt, vmin, vmax, sbplt=[1, 3, 3], tit='Differences')
        plt.tight_layout() 
        
    return ar_filt 

# -----------------------------------------------------------------------------
def lim2grid( xyz, lim, step=None, plot=False, vmin=None, vmax=None, prjcode_in=4326,
              prjcode_out=4326, method='linear', blkm=False, filt=False, radius=None,
              nan=True, padw=0 ) :
    """
    Convert a set of XYZ coordinates to a regular grid within specified limits.

    Parameters:
    xyz (tuple or list): The XYZ coordinates. If len(xyz) <= 2, it should be a tuple (x, y).
                         If len(xyz) > 2, it should be a tuple (x, y, z).
    lim (tuple or list): The limits of the grid in the form (xmin, xmax, ymin, ymax).
    step (float): The grid spacing. If not provided, it will be calculated as the mean distance between points.
    plot (bool): Whether to plot the resulting grid. Default is False.
    vmin (float): The minimum value for the plot color scale. Default is None.
    vmax (float): The maximum value for the plot color scale. Default is None.
    prjcode_in (int): The input projection code. Default is 4326 (WGS84).
    prjcode_out (int): The output projection code. Default is 4326 (WGS84).
    method (str): The method used for interpolation. Default is 'linear'.
    blkm (bool): Whether to perform block mean interpolation. Default is False.
    filt (bool): Whether to apply a 2D filter to the interpolated grid. Default is False.
    radius (int): The radius of the filter. If not provided, it will be calculated based on the grid spacing.
    nan (bool): Whether to remove NaN values from the interpolated grid. Default is True.
    padw (int): The padding width for the filter. Default is 0.

    Returns:
    If len(xyz) <= 2:
    - X (ndarray): The X coordinates of the grid.
    - Y (ndarray): The Y coordinates of the grid.

    If len(xyz) > 2:
    - X (ndarray): The X coordinates of the grid.
    - Y (ndarray): The Y coordinates of the grid.
    - Z (ndarray): The interpolated values of the grid.
    """

    if len( xyz ) <= 2 :
        x, y = xyz
    else :
        x, y, z = xyz
    
    if prjcode_in != prjcode_out :
        x, y = prjxy( prjcode_in, prjcode_out,  x, y )
        
    x, y, idx = xy_in_lim( x, y, lim )
    if len( xyz ) > 2 :
        z = z[idx]

    if step == None :
        step = min_dist( x, y )['mean']
        
    if ( blkm == True ) and  ( len( xyz ) > 2 ) :
        x, y, z = block_m( x, y, z, step*2, method='mean', data_type='vector', lim=lim )  
        
    xg = np.linspace( lim[0], lim[1], int( (lim[1]-lim[0])/step ) ) 
    yg = np.linspace( lim[3], lim[2], int( (lim[3]-lim[2])/step ) )
    X, Y = np.meshgrid( xg, yg )
    
    if len( xyz ) > 2 :
        Z = xyz2xy( ( x, y, z ), ( X, Y ), method=method, fillnan=False )[0]
        
        if filt == True :
            if radius == None :
                radius = int( min_dist( x, y )['mean'] / step )
            Z = filt2d( Z, radius=radius, ftype='mean', padw=padw  )    
        
    if plot == True :
        plta( Z, vmin=vmin, vmax=vmax )
        
    if len( xyz ) <= 2 :
        
        return X, Y
    
    if len( xyz ) > 2 :
        
        if nan == False :
            X = X[ ~np.isnan( Z ) ]
            Y = Y[ ~np.isnan( Z ) ]
            Z = Z[ ~np.isnan( Z ) ]
            
        return X, Y, Z

# -----------------------------------------------------------------------------
def mask2D( xyr=None, xyzgrid=None, array=None, array_and_lim=None, mask=None,
            prjcode=4326, plot=False, vmin=None, vmax=None, s=0,
            convexhull=False, cut2edges=False, in_out='in' ):

    if xyr is not None:
        x, y, r = xyr

    if xyzgrid is not None:
        xx, yy, zz = xyzgrid[0], xyzgrid[1], xyzgrid[2]
        xa, ya, za = xx.ravel(), yy.ravel(), zz.ravel()

    if array is not None :
        zz = np.copy( array )
        x, y = np.arange( zz.shape[1] ), np.arange( zz.shape[0] )
        xx, yy = np.meshgrid( x, y )
        xa, ya, za = xx.ravel(), yy.ravel(), zz.ravel()

    if array_and_lim is not None: 
        zz = array_and_lim[0]
        lim = array_and_lim[1]
        nx, ny = zz.shape[1], zz.shape[0]
        xi, yi = np.linspace(lim[0], lim[1], nx), np.linspace(lim[3], lim[2], ny)
        xx, yy = np.meshgrid(xi, yi)
        xa, ya, za = xx.ravel(), yy.ravel(), zz.ravel()

    if ( xyr is not None ) and ( convexhull is False ) :
        zm = np.empty( np.shape( za ) )
        zm[:] = np.nan
        for p in zip( x, y ):
            d = np.sqrt((p[0] - xa[:]) ** 2 + (p[1] - ya[:]) ** 2)
            zm[d <= r] = za[(d <= r)]

        ZM = zm.reshape( np.shape( zz ) )

    if ( xyr is not None ) and ( convexhull is True ) :
       mhull = xy_in_hull( xa, ya, ( x, y ), buffer=r, plot=False )[0]  
       Mhull = mhull.reshape( np.shape( zz ) ) 
       ZM = np.copy( zz )
       ZM[ ~Mhull ] = np.nan

    if xyr is None:
        ZM = zz

    if mask is not None:
        ZM[ np.isnan( mask ) ] = np.nan

    isfinit = np.isfinite( ZM )

    if plot == True:
        plta(zz, vmin, vmax, tit='original', sbplt=[1, 2, 1])
        plta(ZM, vmin, vmax, tit='mask', sbplt=[1, 2, 2])

    return ZM, isfinit

#------------------------------------------------------------------------------
def resampling( array, factor, spl_order=1, dimention='2D', mode='nearest', 
                plot=False, vmin=None, vmax=None, cval=0.0, nan=True ) :

    if type( array ) in ( list, tuple ) :

        if len( array ) == 3 :
            ax, ay, az = array[0], array[1], array[2]
        if len( array ) == 1 :
            az = np.copy( array )
        IsXy = True
    else :
        az = np.copy( array )
        ax, ay = None, None
        IsXy = False

    if type( factor ) in ( list, tuple ) :
        factor = factor[0]/array.shape[0], factor[1]/array.shape[1]

    if ( np.isnan( az ).any() ) and ( spl_order > 1 ) :
        inan = np.isnan( az )
        az = fillnan( az, method='nearest' )
        IsNan = True
    else :
        IsNan = False

    azr = sp.ndimage.zoom( az, factor, order=spl_order, mode=mode, prefilter=True ) 

    if IsNan == True :
        az [ inan ] = np.nan
        inanr = sp.ndimage.zoom( inan*1, factor, order=0, mode=mode, prefilter=True ).astype(bool)
        azr[ inanr ] = np.nan

    if plot == True :

        if dimention == '2D' :
            plta( az, sbplt=[1,2,1], tit='Original', vmin=vmin, vmax=vmax )
            plta( azr, sbplt=[1,2,2], tit='Resampled', vmin=vmin, vmax=vmax )

        if dimention == '1D' :
            plt.figure()
            plt.plot( az, c='k', label='Original' )
            plt.plot( azr, c='b', label='Resampled' )

    if IsXy == True :

        axr = sp.ndimage.zoom( ax, factor, order=1 )  
        ayr = sp.ndimage.zoom( ay, factor, order=1 ) 

        if nan == False :
            axr = axr[ ~inanr ]
            ayr = axr[ ~inanr ]
            azr = axr[ ~inanr ]

        return [ axr, ayr, azr ]

    else :

        return azr

# -----------------------------------------------------------------------------
def neighbors_points( xy1, xy2, radius, method='circle', plot=False, s=2 ) :

    x1_i, y1_i = xy1 
    shape_i = x1_i.shape
    x2, y2 = xy2
    lim2 = [ np.min(x2)-radius*2, np.max(x2)+radius*2,
             np.min(y2)-radius*2, np.max(y2)+radius*2 ]

    x1, y1, idx_i = xy_in_lim( x1_i, y1_i, lim2 )
    s1, s2 = x1.size, x2.size
    x1, y1 = x1.ravel(), y1.ravel()
    x2, y2 = x2.ravel(), y2.ravel()

    if method == 'circle' :

        x1, y1 = x1.reshape( (s1,1) ), y1.reshape( (s1,1) )
        x2, y2 = x2.reshape( (1,s2) ), y2.reshape( (1,s2) )
        X1 = np.repeat( x1, s2, 1 )
        Y1 = np.repeat( y1, s2, 1 )
        X2 = np.repeat( x2, s1, 0 )
        Y2 = np.repeat( y2, s1, 0 )

        R2 = (X2 - X1)**2 + (Y2 - Y1)**2

        idx_f = np.sum( R2 < radius**2, axis=1 ) > 0
        idx_i[ idx_i ] = idx_f
        idx = idx_i.reshape( shape_i )

    if method == 'square' :

            idx = np.zeros( shape_i, dtype=bool )
            for i in range( s1 ) :
                idx = idx | ( ( x2 - x1[i] )**2 + ( y2 - y1[i] )**2 < radius**2 )

    x_new, y_new = x1_i[idx], y1_i[idx]

    if plot is True :

        plt.scatter( x1, y1, c='k', s=s, label='original points' )
        plt.scatter( x_new, y_new, c='r', s=s, label='selected points' )

    return x_new, y_new, idx

# -----------------------------------------------------------------------------
def xy2convexhull( x, y, plot=False, close=True, size=0.5, color='b' ) :

    points = np.column_stack( ( x.ravel(), y.ravel() ) ) 

    hull = sp.spatial.ConvexHull( points )

    xh, yh = points[hull.vertices,0], points[hull.vertices,1]

    if close is True :
        xh, yh = np.append( xh, xh[0] ), np.append( yh, yh[0] ) 

    if plot is True :

        plt.scatter( x, y, s=size, c=color  )
        plt.plot(xh, yh, 'r--', lw=2)

    return [ hull, xh, yh ]

# -----------------------------------------------------------------------------
def xy_in_hull( x, y, hull, buffer=0, plot=False ) :

    original_shape = np.copy( x.shape )

    x, y = x.ravel(), y.ravel()

    if type( hull ) in ( list, tuple ) :
        xp, yp = hull[0], hull[1]
        hull,_,_ = xy2convexhull( xp, yp )

    xy = np.column_stack( ( x, y ) )

    in_hull = np.all( np.add(np.dot(xy, hull.equations[:,:-1].T ),
                      hull.equations[:,-1] ) <= buffer, axis=1 )

    xh, yh = xp[hull.vertices], yp[hull.vertices]
    xh, yh = np.append( xh, xh[0] ), np.append( yh, yh[0] ) 

    x, y = x.reshape( original_shape ), y.reshape( original_shape )
    in_hull = in_hull.reshape( original_shape ) 

    if plot == True :
        plt.plot( xh, yh, 'r--', lw=2 )
        plt.scatter( x[in_hull], y[in_hull] )

    return in_hull, x[in_hull], y[in_hull] 


# -----------------------------------------------------------------------------
def mask2hull( mask, nan=None, plot=False, close=True, size=0.5, color='b' ) :
    
    if type( mask ) in ( list, tuple ) :
        x, y, mask = mask
    else:
        x = np.arange( mask.shape[1] )
        y = np.arange( mask.shape[0] )
        x, y = np.meshgrid( x, y )
        
    if mask.dtype != 'bool' :
        if nan is not None :
            mask = mask != nan
        else :
            mask = np.isfinite( mask )

    x = x[ mask ]
    y = np.flipud(y)[ mask ]  
    
    hull_list = xy2convexhull( x, y, plot=plot, close=close, size=size, color=color )
    
    return hull_list
    
# -----------------------------------------------------------------------------
def filt_voids( xyzgrid, xyr, step=None, plot=False, vmin=None, vmax=None, 
                method='mean' ) :  

    _, isf = mask2D( xyr, xyzgrid, convexhull=True ) 
#    _, isf2 = mask2D( xyr, xyzgrid ) 
    array = np.copy( xyzgrid[2] )
    m = np.argwhere( isf )
#    plta( m*1.0, lim=(xyzgrid[0],xyzgrid[1]), points=(xyr[0],xyr[1]) )

    xx, yy, zz = xyzgrid
    xp, yp, r = np.copy( xyr )
    x, y = np.concatenate( ( xp, xx[~isf] ) ), np.concatenate( ( yp, yy[~isf] ) )

    if method == 'mean' :
        func = np.nanmean
    if method == 'median' :
        func = np.nanmedian

    if step is None : 
        step = min_dist( x, y )['mean']

    for i,j in m :
        isin = np.full( x.shape, False )
        inc = 0
        while  np.sum( isin ) < 3 :
            radius = r + step * inc    
            isin = ( x > xx[i,j] - radius ) & ( x < xx[i,j] + radius ) & \
                   ( y > yy[i,j] - radius ) & ( y < yy[i,j] + radius )
            inc += 1

        if inc == 0 : 
            continue

        win = ( xx > xx[i,j] - radius ) & ( xx < xx[i,j] + radius ) & \
              ( yy > yy[i,j] - radius ) & ( yy < yy[i,j] + radius )

        array[ i,j ] = func( array[ win ] )

    if plot == True :
        plta( xyzgrid[2], lim=(xyzgrid[0],xyzgrid[1]), points=(xyr[0],xyr[1]), 
              vmin=vmin, vmax=vmax, tit='Original', sbplt=[1,3,1] ) 
        plta( array, lim=(xyzgrid[0],xyzgrid[1]), points=(xyr[0],xyr[1]), 
              vmin=vmin, vmax=vmax, tit='Filtered', sbplt=[1,3,2] )
        plta( xyzgrid[2]-array, lim=(xyzgrid[0],xyzgrid[1]), points=(xyr[0],xyr[1]), 
              vmin=None, vmax=None, tit='Differences', sbplt=[1,3,3] )        

    return array

# -----------------------------------------------------------------------------
def mosaic_array( array_list, plot=False, vmin=None, vmax=None, ref_shape=0,
                 spl_order=1) :

    shape = array_list[ ref_shape ].shape

    array_mos = np.copy( array_list[ ref_shape ] )

    new_list = []

    for i, a in enumerate( array_list ) :
        shape_ratio = shape[0]/a.shape[0], shape[1]/a.shape[1]
        ar = sp.ndimage.zoom( a, shape_ratio, order=spl_order )
        new_list.append( ar )
        isn = np.isnan( array_mos ) 
        array_mos[isn] = ar[isn]

    if plot == True :
        plta( array_mos, vmin=vmin, vmax=vmax)    

    return array_mos, new_list

# -----------------------------------------------------------------------------
def del_in_lim(x, y, array, lim, remove='in'):
    """
    Delete points within or outside a specified limit.

    Args:
        x (array-like): X-coordinates of the points.
        y (array-like): Y-coordinates of the points.
        array (array-like): Array of values associated with the points.
        lim (float): Limit value [xmin, xmax, ymin, ymax].
        remove (str, optional): Specifies whether to remove points within ('in') or outside ('out') the limit. Defaults to 'in'.

    Returns:
        list: A list containing the filtered x-coordinates, y-coordinates, and array values.
    """

    ar = np.copy( array )
    x = np.asarray( x )
    y = np.asarray( y )
    idx = xy_in_lim( x, y, lim )[2]

    if remove == 'in' :
        ar = ar[ ~idx ] 
        x = x[ ~idx ]
        y = y[ ~idx ]

    if remove == 'out' :
        ar = ar[ idx ]
        x = x[ idx ]
        y = y[ idx ]

    return [ x, y, ar ]
    
# -----------------------------------------------------------------------------
def XY2bbox( X, Y, sx, sy=None ) :
    
    if sy is None : 
        sy = sx
    
    minbbX = np.nanmin( X ) - sx/2
    maxbbX = np.nanmax( X ) + sx/2
    
    minbbY = np.nanmin( Y ) - sy/2
    maxbbY = np.nanmax( Y ) + sy/2
    
    bbox = [ minbbX, maxbbX, minbbY, maxbbY ]
    
    return bbox

# -----------------------------------------------------------------------------
def julian2date( jday, year, array=True, array_type=dict ) :
    
    if ( type(jday) is int ) and ( type(year) is int ) :
        date = datetime.datetime(int(year), 1, 1)+datetime.timedelta(days=int(jday) -1) 
        
    else :
        date = []
        for y, d in zip( year, jday ) :
            date.append( datetime.datetime(int(round(y)), 1, 1)+\
                         datetime.timedelta(days=int(round(d)) -1) )
            
        if array is True :
            if array_type == 'datetime64[D]' :
                date = np.array( date, dtype='datetime64[D]' )
            if array_type == str :
                date = np.array( date, dtype='datetime64[D]' )
                date = date.astype( str )
            if array_type == dict :
                date = np.array( [ [ x.year, x.month, x.day ] for x in date ] )
                date = { 'year':date[:,0], 'month':date[:,1], 'day':date[:,2] }
            if array_type == 'array' :
                date = np.array( [ [ x.year, x.month, x.day ] for x in date ] )

    return date

# -----------------------------------------------------------------------------
def date2julian( yy, mm, dd ) :
    
    if type( yy ) in ( float, int ) :
        yy = np.array( [ yy ] )
        mm = np.array( [ mm ] )
        dd = np.array( [ dd ] )
        
    jday = np.zeros( yy.shape )    
    for i,_ in enumerate( yy ):
        
        day = datetime.date( yy[i], mm[i], dd[i] ) 
        day = day.toordinal() + 1721424.5
        jday[i] = day
    
    return jday

# -----------------------------------------------------------------------------
def date2datetime( date, time, fmt='%d-%m-%y %H:%M:%S', array=True,
                   array_type='datetime64[ms]' ) :
    
    if ( type( date ) is str ) and ( type( time ) is str ) :
        date_time_str = date + ' ' + time 
        date_time = datetime.strptime( date_time_str, fmt )
        
    else :   
        date_time = []
        for d, t in zip( date, time ) :
            date_time_str = d + ' ' + t 
            date_time.append( datetime.datetime.strptime( date_time_str, fmt ) ) 
            
        if array is True :
            # a = pd.to_datetime( date_time )
            np_date_time = np.array( date_time, dtype=array_type )   
            
        if array is True :
            if array_type == 'datetime64[D]' :
                date = np.array( date, dtype='datetime64[D]' )
            if array_type == str :
                date = np.array( date, dtype='datetime64[D]' )
                date = date.astype( str )
            if array_type == dict :
                date = np.array( [ [ x.year, x.month, x.day ] for x in date ] )
                date = { 'year':date[:,0], 'month':date[:,1], 'day':date[:,2] }
            if array_type == 'array' :
                date = np.array( [ [ x.year, x.month, x.day ] for x in date ] )
    
    return np_date_time

# -----------------------------------------------------------------------------
def combine64(years, months=1, days=1, hours=None, minutes=None,
              seconds=0, milliseconds=None, microseconds=None, 
              nanoseconds=None, weeks=None ):
    """
    Combines the given time components into a single numpy datetime type.
    
    Args:
        years (int or array-like): The number of years.
        months (int or array-like, optional): The number of months. Defaults to 1.
        days (int or array-like, optional): The number of days. Defaults to 1.
        hours (int or array-like, optional): The number of hours.
        minutes (int or array-like, optional): The number of minutes.
        seconds (int or array-like, optional): The number of seconds. Defaults to 0.
        milliseconds (int or array-like, optional): The number of milliseconds.
        microseconds (int or array-like, optional): The number of microseconds.
        nanoseconds (int or array-like, optional): The number of nanoseconds.
        weeks (int or array-like, optional): The number of weeks.
    
    Returns:
        numpy.datetime64: The combined datetime type.
    """
    
    years = np.asarray(years) - 1970
    months = np.asarray(months) - 1
    days = np.asarray(days) - 1
    
    types = ( '<M8[Y]', '<m8[M]', '<m8[D]', '<m8[W]', '<m8[h]',
              '<m8[m]', '<m8[s]', '<m8[ms]', '<m8[us]', '<m8[ns]' )
    
    if np.any( np.mod( seconds, 1 ) != 0 ) :
        nanoseconds = seconds.copy() * 1e9
        seconds=None
    
    vals = ( years, months, days, weeks, hours, minutes, seconds,
             milliseconds, microseconds, nanoseconds )
    
    datetime_type = sum(np.asarray(v, dtype=t) for t, v in zip(types, vals)
                    if v is not None )
    
    return datetime_type

# -----------------------------------------------------------------------------
def read_file(file, rows=None):
    """
    Read lines from a file.

    Args:
        file (str): The path to the file.
        rows (int or list, optional): The rows to read from the file. 
        If None, all rows are read. 
        If an integer, read that number of rows starting from the first row. 
        If a list, read the rows with the specified indices. Defaults to None.

    Returns:
        list: The lines read from the file.
    """
    f = open( file, 'r' )
    
    lines = [line.rstrip() for line in f]
    
    if type( rows ) == int :
        rows = list( range( rows ) )
        
    if ( rows == [] ) or ( rows is None ) :
        rows = list( range( len( lines ) ) )   
        
    lines = [ lines[index] for index in rows ]
    
    f.close()
    
    return lines

# -----------------------------------------------------------------------------
def open_dir( path ) :
    
    if os.name == "nt" :
        os.system( 'start ' + path )
        
    if os.name == 'posix' :
        os.system( 'nautilus ' + path + '&' )

# -----------------------------------------------------------------------------
def del_files(path, string):
    """
    Delete files and directories that contain the specified string in the given path.

    Args:
        path (str): The path to search for files and directories.
        string (str): The string to match in the file or directory names.

    Returns:
        None
    """

    for file in os.listdir( path ):

        full_path = os.path.join( path, file )
        
        if string in file:

            if os.path.isfile( full_path ):

                os.remove( full_path )

            elif os.path.isdir( full_path ):
        
                empty_directory( full_path )
                shutil.rmtree( full_path )

# -----------------------------------------------------------------------------
def ls(path):
    """
    List all files in the specified directory.

    Args:
        path (str): The path to the directory.

    Returns:
        None
    """
    for file in os.listdir( path ) :
        
        print( file )
        
# -----------------------------------------------------------------------------
def xyz2grid( x, y, z, lim=None, extend=None, extend_method='percentage', 
              sqr_area=False, gstep=None, blkm=None, method='linear', prjcode_in=4326,
              prjcode_out=4326, filt=None, filt_radius=1, msk_radius=None, msk_shp=None,
              in_out='in', fill_nan=None, plot=False, vmin=None, vmax=None,
              pltxy=False, s=None, msk_arr=None, adjst_lim=False, iter=1, tension=0.35,
              filt_factor=2, filt_sigma=1, padw=0, pmode='gdal' ):
              
   
    if prjcode_in != prjcode_out:
        x, y = prjxy( prjcode_in, prjcode_out, x, y )  
    
    if lim is None:
        lim = [np.min(x), np.max(x), np.min(y), np.max(y)]
        xl, yl, zl = x, y, z
    else:
        xl, yl, indx = xy_in_lim( x, y, lim, extend=33 )
        zl = z[indx]
        
    if extend is not None:
        lim = extend_lim(lim, extend, extend_method, sqr_area )

    if gstep == None:
        gstep = min_dist( xl, yl )['mean']
    
    if blkm != None :    
        if blkm == True:
            xl, yl, zl = block_m( xl, yl, zl, gstep, lim=lim )
        if type( blkm ) in ( int, float ) :
            xl, yl, zl = block_m( xl, yl, zl, blkm, lim=lim )
    
    if adjst_lim is False :    
        xg = np.arange( lim[0], lim[1], gstep )
        yg = np.arange( lim[3], lim[2], -gstep )
    if adjst_lim is True :   
        xg = np.linspace( lim[0], lim[1], int( ( lim[1] - lim[0] ) / gstep ) )
        yg = np.linspace( lim[3], lim[2], int( ( lim[3] - lim[2] ) / gstep ) )
    xx, yy = np.meshgrid(xg, yg)
       

#       zz = scy.interpolate.griddata(points, zl, (xx, yy), method=method)
    zz = xyz2xy( ( xl, yl, zl ), (xx, yy), method=method, fillnan=False )[0]

    if fill_nan is not None :
        if ( type(fill_nan) == bool ) and ( fill_nan is True ) :
            fill_nan = 'gdal'
        zz = fillnan( zz, xy=( xx, yy ), method=fill_nan, iter=iter ) 
        
    if filt is not None :
        zz = filt2d( zz, ftype=filt, iter=iter, radius=filt_radius, 
                         factor=filt_factor, sigma=filt_sigma, padw=padw, pmode=pmode )        
    
    if msk_radius is not None:
        zz = mask2D( ( xl, yl, msk_radius ), (xx, yy, zz) )[0]
        
    # if msk_shp is not None:
    #     zz = rt.mask_array( xx, yy, zz, msk_shp, prjcode=prjcode_out )

    if msk_arr is not None:
        zz[np.isnan( msk_arr ) ] = np.nan         

    if plot == True:
        if pltxy == False:
            plta( zz, vmin=vmin, vmax=vmax, cmap='rainbow' )
        if pltxy == True:
            limp = xy2lim( xx, yy ) 
            plta( zz, vmin=vmin, vmax=vmax, cmap='rainbow', lim=limp, points=[ xl, yl ] )
            
    return [ xx, yy, zz ], [ x, y, z ]
    
# -----------------------------------------------------------------------------
def mesh3Dmodel( x=None, y=None, z=None, lim=None, step=None, plot=False, 
                 sbplt=111, out_type='dictionary', size=None, centre_xy=False, 
                 xm=None, ym=None, crop2lim=True, pad=None, padw=None ) :
    """
    Generate a 3D mesh model based on the given parameters.

    Parameters:
    - x (ndarray, optional): X-coordinates of the mesh points.
    - y (ndarray, optional): Y-coordinates of the mesh points.
    - z (ndarray, optional): Z-coordinates of the mesh points.
    - lim (ndarray, optional): Limits of the mesh in the form [xmin, xmax, ymin, ymax, zmin, zmax].
    - step (float or tuple, optional): Step size for the mesh in each dimension. If a single value is provided, it is used for all dimensions.
    - plot (bool, optional): Whether to plot the mesh.
    - sbplt (int, optional): Subplot configuration for the plot.
    - out_type (str, optional): Output type of the mesh data. Possible values are '3darray', '1darray', and 'dictionary'.
    - size (float, optional): Size of the points in the plot.
    - centre_xy (bool, optional): Whether to center the mesh in the X and Y dimensions based on the mean values of the limits.
    - xm (float, optional): X-coordinate of the center point. If not provided, the mean of the X limits is used.
    - ym (float, optional): Y-coordinate of the center point. If not provided, the mean of the Y limits is used.
    - crop2lim (bool, optional): Whether to crop the mesh points to fit within the limits.
    - pad (float or tuple, optional): Padding size for the mesh in each dimension. If a single value is provided, it is used for all dimensions.
    - padw (float or tuple, optional): Padding step size for the mesh in each dimension. If a single value is provided, it is used for all dimensions.

    Returns:
    - If out_type is '3darray', returns X, Y, Z as 3D arrays representing the mesh points.
    - If out_type is '1darray', returns X, Y, Z as 1D arrays representing the mesh points, along with the shape of the mesh.
    - If out_type is 'dictionary', returns a dictionary containing 'x', 'y', and 'z' as 1D arrays representing the mesh points, along with the shape of the mesh.
    """

    # If step is a single number, convert it to a tuple for all dimensions
    if type(step) in (int, float) :
        step = step, step, step   

    # Make a copy of the limits
    lim = np.copy( lim )          

    # If limits and step are provided and centering is requested
    if ( lim is not None ) and ( step is not None ) and ( centre_xy is True ) :
        # If no center x-coordinate is provided, calculate it as the mean of x limits
        if xm is None :
            xm = np.round( np.mean( ( lim[1], lim[0] ) ) )
        # Calculate the range of x-coordinates
        xr = np.max( (np.ceil((xm-lim[0])/step[0]), np.ceil((lim[1]-xm)/step[0]) ) )
        # Update the x limits based on the center and range
        lim[0] = xm - step[0] * xr
        lim[1] = xm + step[0] * xr
        # Do the same for y-coordinates
        if ym is None :
            ym = np.round( np.mean( ( lim[2], lim[3] ) ) )
        yr = np.max( (np.ceil((ym-lim[2])/step[1]), np.ceil((lim[3]-ym)/step[1]) ) )
        lim[2] = ym - step[1] * yr
        lim[3] = ym + step[1] * yr

    # If step is a single number, convert it to a tuple for all dimensions
    if type(step) in (int, float) :
        step = step, step, step

    # If x-coordinates are provided, use them, otherwise generate them based on limits and step
    if x is not None :
        xu = np.unique( x )
    else :
        xu = np.arange( lim[0], lim[1]+step[0], step[0] )

    # Do the same for y and z coordinates
    if y is not None :
        yu = np.unique( y )
    else :
        yu = np.arange( lim[2], lim[3]+step[1], step[1] )

    if z is not None :
        zu = np.unique( z )
    else :
        zu = np.arange( lim[5], lim[4]-step[2], -step[2] )

    # If cropping is requested, crop the coordinates to fit within the limits
    if ( lim is not None ) and ( crop2lim is True ) :
        xu = xu[ (xu>=lim[0]) & (xu<=lim[1]) ]
        yu = yu[ (yu>=lim[2]) & (yu<=lim[3]) ]
        zu = zu[ (zu>=lim[4]) & (zu<=lim[5]) ]

    # If padding is requested
    if pad is not None :
        # If pad is a single number, convert it to a tuple for all dimensions
        if type( pad ) in ( int, float ) :
            padx, pady, padz = int(pad), int(pad), int(pad)
        else :
            padx, pady, padz = int(pad[0]), int(pad[1]), int(pad[2])
        # If no padding step is provided, use the same as the mesh step
        if padw is None :
            padwx, padwy, padwz = step[0], step[1], step[2]
        # If padw is a single number, convert it to a tuple for all dimensions
        if type( padw ) in ( int, float ) :
            padwx, padwy, padwz = padw, padw, padw 
        # If padw is a tuple, unpack it
        if type( padw ) in ( list, tuple ) :
            padwx, padwy, padwz = padw[0], padw[1], padw[2]

        # Add padding to the coordinates
        xu = np.concatenate( ( np.arange(xu.min(),xu.min()-padx*padwx*1.5, -padwx)[1:], xu ) )
        yu = np.concatenate( ( np.arange(yu.min(),yu.min()-pady*padwy*1.5, -padwy)[1:], yu ) )
        zu = np.concatenate( ( np.arange(zu.max(),zu.max()+padz*padwz*1.5, padwz)[1:][::-1] , zu ) )

        xu = np.concatenate( ( xu, np.arange(xu.max(),xu.max()+padx*padwx*1.5, padwx)[1:] ) )
        yu = np.concatenate( ( yu, np.arange(yu.max(),yu.max()+pady*padwy*1.5, padwy)[1:] ) )
        zu = np.concatenate( ( zu, np.arange(zu.min(),zu.min()-padz*padwz*1.5, -padwz)[1:]) )

    # Generate the 3D mesh
    X, Y, Z = np.meshgrid( xu, yu, zu )
    # Get the shape of the mesh
    shape = X.shape

    # If plotting is requested
    if plot is True :
        # If no size is provided, use 1
        if size is None :
            size = 1
        # Create a 3D subplot
        ax = plt.subplot(sbplt, projection='3d')
        # Plot the mesh points
        ax.scatter( X[:,int(shape[1]/2),:].ravel(), 
                    Y[:,int(shape[1]/2),:].ravel(), 
                    Z[:,int(shape[1]/2),:].ravel(), s=size, c='b' )
        ax.scatter( X[int(shape[0]/2),:,:].ravel(), 
                    Y[int(shape[0]/2),:,:].ravel(), 
                    Z[int(shape[0]/2),:,:].ravel(), s=size, c='b' )

    # Return the mesh data in the requested format
    if out_type == '3darray' :
        return X, Y, Z
    if out_type == '1darray' :
        return X.ravel(), Y.ravel(), Z.ravel(), shape
    if out_type == 'dictionary' :
        return { 'x':X.ravel(), 'y':Y.ravel(), 'z':Z.ravel() }, shape
    
# -----------------------------------------------------------------------------
def datetime2time(  dates, unit='s', ref_time='1970-01-01T00:00:00' ) :
    """
    Convert a numpy array of datetime64 objects to time since the ref_time.

    Parameters:
    dates (numpy.ndarray): Array of datetime64 objects.
    unit (str): Unit of time. Can be 's' (seconds), 'm' (minutes), 'h' (hours), 
                'D' (days), 'M' (months), or 'Y' (years).
    ref_time (str): starting reference time from wich to count; 
                    default Unix epoch (i.e., '1970-01-01T00:00:00').

    Returns:
    numpy.ndarray: Array of time since the Unix epoch in the specified units.
    """

    if unit == 'Y':
        # Convert datetime64 objects to days since the Unix epoch, then to years
        time = (dates - np.datetime64(ref_time)) / np.timedelta64(1, 'D')
        time /= 365.25 # Average length of a year in the Gregorian calendar

    elif unit == 'M':
        # Convert datetime64 objects to days since the Unix epoch, then to months
        time = (dates - np.datetime64(ref_time)) / np.timedelta64(1, 'D')
        time /= 30.44 # Average length of a month in the Gregorian calendar
    
    else:
        # Convert datetime64 objects to time since the Unix epoch in the specified units
        time = (dates - np.datetime64(ref_time)) / np.timedelta64(1, unit)

    return time


# -----------------------------------------------------------------------------
def are_equal(dir_file_1, dir_file_2):
    """
    Check if two files or directories are equal.

    Args:
        dir_file_1 (str): The path to the first file or directory.
        dir_file_2 (str): The path to the second file or directory.

    Returns:
        bool: True if the files or directories are equal, False otherwise.
    """

    # If dir_file_1 does not exist, raise a FileNotFoundError
    if not os.path.exists(dir_file_1):
        raise FileNotFoundError(f"{dir_file_1} does not exist")

    # If dir_file_2 does not exist, return False
    if not os.path.exists(dir_file_2):
        return False

    # If the absolute paths of dir_file_1 and dir_file_2 are not the same, return False
    if os.path.abspath(dir_file_1) != os.path.abspath(dir_file_2):
        return False

    # Check if both dir_file_1 and dir_file_2 are files
    if os.path.isfile(dir_file_1) and os.path.isfile(dir_file_2):
        # If they are files, compare them directly and store the result in 'value'
        value = filecmp.cmp(dir_file_1, dir_file_2, shallow=False)
    else:
        # If they are not files (i.e., they are directories), create a set of all file paths in dir_file_1
        set1 = set(os.path.join(root, file) 
                for root, dirs, files in os.walk(dir_file_1) 
                for file in files)
        
        # Create a set of all file paths in dir_file_2
        set2 = set(os.path.join(root, file) 
                for root, dirs, files in os.walk(dir_file_2) 
                for file in files)
        
        # Compare the contents of the files in the directories. If the paths and the contents of the files are the same,
        # 'value' will be True. Otherwise, it will be False.
        value = all(filecmp.cmp(file1, file2, shallow=False) 
                    for file1, file2 in zip(sorted(set1), sorted(set2)))
        
    # Return the result of the comparison
    return value

# -----------------------------------------------------------------------------
def empty_directory( dir_path ):
    """
    Deletes all files and subdirectories within a given directory.

    Args:
        dir_path (str): The path to the directory.

    Raises:
        OSError: If there is an error while deleting the files or directories.

    """

    if os.path.isdir( dir_path ):
        # Iterate over all files and directories in the given directory
        for filename in os.listdir(dir_path):
            # Construct the full path of the current file or directory
            file_path = os.path.join(dir_path, filename)
            try:
                # If the current item is a file or a symbolic link, delete it
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                # If the current item is a directory, delete it and all its contents
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                # If an error occurred while trying to delete the item, print the error message
                print('Failed to delete %s. Reason: %s' % (file_path, e))

# -----------------------------------------------------------------------------
def copy_dir_content(src_dir, dst_dir, include_strings=None, exclude_strings=None):
    """
    Copy the contents of a directory to another directory, 
    with options to include or exclude specific strings.

    Args:
        src_dir (str): The source directory path.
        dst_dir (str): The destination directory path.
        include_strings (list, optional): List of strings to include. Defaults to None.
        exclude_strings (list, optional): List of strings to exclude. Defaults to None.

    Returns:
        None
    """
    # Create the destination directory if it does not exist
    os.makedirs(dst_dir, exist_ok=True)

    # List all files and directories in the source directory
    for item in os.listdir(src_dir):
        s = os.path.join(src_dir, item)
        d = os.path.join(dst_dir, item)

        # Check if the item should be included or excluded
        if include_strings and not any(inc in item for inc in include_strings):
            continue
        if exclude_strings and any(exc in item for exc in exclude_strings):
            continue

        # If item is a file, copy it
        if os.path.isfile(s):
            shutil.copy2(s, d)
        
        # Else if item is a directory, copy the directory
        elif os.path.isdir(s):
            if os.path.isdir(d):  # If directory already exists, remove it
                shutil.rmtree(d)
            shutil.copytree(s, d)
            
# -----------------------------------------------------------------------------
def print_table( table, headers=None, 
                 space=12, decimals=2, 
                 rows=20, cols=[],
                 return_array=False, 
                 idx=None, 
                 title=None,
                 row_index=True,
                 col_index=True, 
                 center_title=False,
                 return_str=False,
                 colsxrow=6,
                 path_name=None, 
                 printf=True ):
    """
    Prints a table of values with optional headers, formatting, and row/column selection.

    Args:
        - table (list, tuple, dict, numpy.ndarray): The table of values to be printed.
        - headers (list, optional): The headers for each column. Defaults to None.
        - space (int, optional): The width of each column. Defaults to 12.
        - decimals (int, float, list, optional): The number of decimal places to display for each column.
            If a single value is provided, it will be applied to all columns. If a list is provided,
            each value will be applied to the corresponding column. Defaults to 2.
        - rows (int, list, optional): The indices of the rows to be included in the table. Defaults to [].
        - cols (int, list, optional): The indices of the columns to be included in the table. Defaults to [].
        - idx (boolean array, optional): The indices of the rows to be included in the table. Defaults to None.
        - title (str, optional): The title to be printed above the table. Defaults to None.
        - row_index (bool, optional): Whether to include a column of row numbers. Defaults to True.
        - col_index (bool, optional): Whether to include a row of column numbers. Defaults to True.
        - center_title (bool, optional): Whether to center the title above the table. Defaults to False.
        - return_str (bool, optional): Whether to return the table as a string. Defaults to False.
        - return_array (bool, optional): Whether to return the table as a numpy array. Defaults to False.
        - colsxrow (int, optional): The maximum number of columns to be printed per row. Defaults to 6.
        - path_name (str, optional): The path and name of the file to save the table. Defaults to None.
        - printf (bool, optional): Whether to print the table. Defaults to True.

    Returns:
        - None: If return_str is False and return_array is False.
        - str: If return_str is True.
        - numpy.ndarray: If return_array is True.
    """

    # Initialize an empty string to hold the output
    output = ""

    # Check if the input table is a numpy array
    if type(table) == np.ndarray:
        # Create a copy of the table
        table = table.copy()

    # Check if the input table is a dictionary
    if type(table) == dict:
        table = copy.copy(table)
        # If there are keys to exclude, create a copy of the table and remove those keys
        if cols != []:
            for i, c in enumerate(cols):
                if type(c) == str:
                    cols[i] = list(table.keys()).index(c) 
        for k in table.keys():
            if len( k ) > space:
                space = len( k ) + 2
        # Convert the dictionary to an array
        table, headers = dict2array(table)

    # Check if the input table is a list or a tuple
    if type(table) in (list, tuple):
        # Convert each element of the table to a numpy array and stack them column-wise
        for i, t in enumerate(table):
            if type(t) in ( list, tuple, int, float ):
                ti = np.array(t)
            if i == 0:
                array = ti
            else:
                array = np.column_stack((array, ti))
        # Replace the original table with the newly created 2D array
        table = array.copy()

    # If the table is a 1D array and headers are not provided, reshape it into a 2D array
    if ( len(table.shape) == 1 or (len(table.shape) == 2 and table.shape[0] == 1) ) and ( not headers ):
        table = table.reshape(-1, 1)

    # If the table is a 1D array and headers are provided, reshape it into a 2D array 
    # with one row and as many columns as the length of the headers
    if ( len(table.shape) == 1 or (len(table.shape) == 2 and table.shape[0] == 1) ) and headers :
        table = table.reshape(1, len(headers) )

    # If decimals is a single number, create a list of that number repeated for each column of the table
    if type(decimals) in (int, float):
        decimals = [decimals] * len(table[0])

    # Create a boolean array indicating whether each element of the table is an integer
    is_int = np.array([[isinstance(val, (int, float)) and val % 1 == 0 for val in row] for row in table])

    # If idx is not None, select the specified rows from the table
    if idx is not None :
        table = table[idx]

    if ( rows is not None ) and ( rows != [] ) and ( rows is not False ) :
        # If rows or cols is a single number, create a list of numbers up to that number
        if isinstance(rows, int):
            rows = list(range(rows))
        if isinstance(rows, np.ndarray):
            if rows.dtype == bool:
                rows = np.where(rows)[0].tolist()

    if ( cols is not None ) and ( cols != [] ) :

        if headers is not None:
            headers = [ headers[i] for i in cols ]
            
        table = table[:, cols]

    # If rows or cols is not specified, 
    # create a list of numbers up to the number of rows or columns in the table    
    if not rows:
        if isinstance(table[0], str):
            rows = list(range(len(table))) 
        else:
            rows = list(range(len(table[:,0])))
        
    if not cols:
        if isinstance(table[0], str):
            cols = list(range(len(table))) 
        else:
            cols = list(range(len(table[0])))  

    # If row_index is True, add a column of row numbers to the table
    len_row_index = len( str( len( rows ) ) ) 

    if row_index:
        starting_space_len = len_row_index + 1
        starting_space_str = " " * ( starting_space_len )
    else:
        starting_space_len = 0
        starting_space_str = ""

    # Print the table headers
    if headers is None : 
        if col_index:
            output += starting_space_str
            for j in cols:
                output += f"% {space}d" % j
            output += "\n"
        output += starting_space_str
        for j in range(len(table[0])):
            output += "-"*space
        output += "\n"

    else :
        output += starting_space_str
        for j in headers:
            output += f"% {space}s" % j
        output += "\n"
        if col_index:
            output += starting_space_str
            for j in cols:
                output += f"% {space}d" % j
            output += "\n"
        output += starting_space_str
        for j in range( len( table[0] ) ):
            output += "-"*space
        output += "\n"     

    # Print the contents of the table
    for i in range( table.shape[0] ):

        if i not in rows:
            continue

        if row_index:
            output += f"%{len_row_index}d|" % (i) # Row nums
        
        for j in range( table.shape[1] ):

            if isinstance(table[i][j], (int, float)):
            
                if is_int[i][j]:
                    ft = 'd'
                else:
                    ft = f".{decimals[j]}f"
                output += f"% {space}{ft}" % (table[i][j])
            
            else:
                # Truncate the string if it's longer than `space`
                str_val = str(table[i][j])
                if len(str_val) > space:
                    str_val = str_val[:space-1]

                output += f"% {space}s" % (str_val)

        output += "\n"

    if not colsxrow:
        colsxrow = table.shape[1]

    if table.shape[1] > colsxrow:

        sections = range(0, table.shape[1], colsxrow)
        out_section = ""

        # Split the string into lines
        lines = output.split('\n')

        for i, section in enumerate(sections):

            if i == 0 :
                for line in lines:
                    chunck = line[ : space *colsxrow + starting_space_len ] + '\n'
                    out_section += chunck

            else:
                for line in lines:
                    chunck = line[ section * space + starting_space_len : section * space + space * colsxrow+starting_space_len ] +'\n'
                    if row_index:
                        chunck = line[ : starting_space_len ] + chunck
                    out_section += chunck

        output = out_section


    # If a title is provided, print it centered above the table
    if title is not None:
        if center_title:
            title_str = "{:^{}}\n".format(title, space * len(table[0]))
        else:
            title_str = title + "\n"
        output = title_str + output

    if printf is True :
        print(output)

    print( "\nTotal rows: " + str(table.shape[0]) + '\n' )

    if ( path_name is not None ) and ( path_name is not False ) :
        with open( path_name, 'w' ) as f:
            f.write( output )
        f.close()

    # If return_str is True, return the output string
    if return_str is True :
        return output

    # If return_str is True, return the table
    if return_array is True :
        return table

# -----------------------------------------------------------------------------
def is_last( iterable ):
    """
    Returns an iterator that yields a tuple of the current element 
    and a boolean value indicating whether 
    it is the last element in the iterable.
    
    Args:
        - iterable: An iterable object.
    
    Yields:
        - tuple: A tuple containing the current element 
            and a boolean value indicating whether it is the last element.
    """

    it = iter( iterable )
    
    last = next( it )

    for val in it:
        
        yield last, False
        
        last = val

    yield last, True

# -----------------------------------------------------------------------------
def save_output_to_html():

    # Temporary files to hold output
    stdout_temp = sys.stdout
    stderr_temp = sys.stderr
    sys.stdout = open('stdout.txt', 'w')
    sys.stderr = open('stderr.txt', 'w')

    # Run the script
    sys.stdout.seek(0)
    sys.stderr.seek(0)
    stdout_content = sys.stdout.read()
    stderr_content = sys.stderr.read()
    
    # Print the captured content to the console
    stdout_temp.write(stdout_content)
    stderr_temp.write(stderr_content)

    # Close and restore stdout and stderr
    sys.stdout.close()
    sys.stderr.close()
    sys.stdout, sys.stderr = stdout_temp, stderr_temp

    # HTML escape only to a necessary extent
    stdout_html = html.escape(stdout_content, quote=False)
    stderr_html = html.escape(stderr_content, quote=False)

    # Check for images
    image_files = [f for f in os.listdir() if f.endswith(('.png', '.jpg'))]

    html_content = f"""
    <html>
    <head><title>Script Output</title></head>
    <body>
    <h1>Output</h1>
    <h2>Standard Output</h2>
    <pre>{stdout_html}</pre>
    <h2>Standard Error</h2>
    <pre>{stderr_html}</pre>
    <h2>Figures</h2>
    """
    
    for image in image_files:
        html_content += f'<img src="{image}" alt="{image}"><br>'
    
    html_content += "</body></html>"

    # Save HTML file
    with open('output.html', 'w') as f_html:
        f_html.write(html_content)
    
    return f_html, html_content

# ----------------------------------------------------------------------------- 
def fit_segmented_line( x, y, 
                        num_segments=None, 
                        threshold=1e-2, 
                        null_intercept=False,
                        keep_nodes=[], 
                        move_nodes=[],
                        place_nodes=[], 
                        plot=False, 
                        learning_rate=0.01 ):
    """
    Fits a segmented line to the given data points.

    Parameters:

        - x (array-like): The x-coordinates of the data points.
      
        - y (array-like): The y-coordinates of the data points.
      
        - num_segments (int, optional): The number of line segments to fit. 
            If not specified, the algorithm will automatically determine the number of segments.
        
        - threshold (float, optional): The convergence threshold for the algorithm. 
            The algorithm stops when the difference in standard deviation between 
                iterations is below this threshold.
        
        - null_intercept (bool, optional): Whether to include a null intercept in the line fitting. 
            If True, the line will pass through the origin (0, 0).
        
        - plot (bool, optional): Whether to plot the fitted line segments.

    Returns:

        - array-like: An array of shape (num_segments, 2) 
            containing the x and y coordinates 
            of the fitted line segments.
    """

    if null_intercept :
        x = np.vstack([0] + x.tolist()).ravel()
        y = np.vstack([0] + y.tolist()).ravel()

    coefficients = np.polyfit( x, y, 1 ) 
    std = np.std( y )
    diff = std.copy()
    x_pol = np.linspace( 0, x.max(), 100 )
    y_pol = np.polyval( coefficients, x_pol )
    y_fit = np.polyval( coefficients, x )

    i = 0
    while diff > threshold :
        i += 1
        Mx = np.ones((x.shape[0], i+1 ))
        for j in range( i+1 ):
            Mx[:, j] = x.ravel()**j
        coefficients, _, _, _ = np.linalg.lstsq( Mx, y, rcond=None )
        y_fit = np.polyval( np.flip( coefficients ), x )
        res = y - y_fit
        stdi = np.std( res )
        diff = std - stdi
        std = stdi

    y_pol = np.polyval( np.flip( coefficients ), x_pol )

    # Calculate the first and second derivatives of y_pol
    y_pol_fd = np.diff( y_pol )[:-1]
    x_pol_fd = x_pol[1:-1] - ( x_pol[1] - x_pol[0] ) / 2
    y_pol_sd = np.diff( y_pol_fd )[:-1]
    x_pol_sd = x_pol_fd[1:-1] - ( x_pol_fd[1] - x_pol_fd[0] ) / 2

    y_pol_int = np.interp( x_pol_sd, x_pol_fd, y_pol_fd )

    # Calculate the curvature radius
    curvature = np.abs( y_pol_sd ) / ( 1 + y_pol_int**2 )**1.5
    # curvature = np.concatenate(([0,0], curvature, [0,0]))

    # Calculate the differences between consecutive elements
    diff = np.diff(curvature)

    # Find where the differences change sign
    sign_changes = np.diff(np.sign(diff))

    # The relative maxima are where the sign changes from positive to negative
    maxima_indices = np.where(sign_changes == -2)[0] + 1

    # Get the x_pol_sd values at the maxima
    curv_rel_max_x = x_pol_sd[maxima_indices]
    curv_rel_max = curvature[maxima_indices]

    if num_segments is None :
        num_segments = len( curv_rel_max ) + 1

    # Select changing points
    isort = np.flip( np.argsort( curv_rel_max ) )
    xes = np.sort( curv_rel_max_x[isort][: num_segments-1] )
    yei = np.interp( xes, x_pol, y_pol )
    xe = np.concatenate( ( [x_pol[0]], xes, [x_pol[-1]] ) ) 
    ye = np.concatenate( ( [y_pol[0]], yei, [y_pol[-1]] ) ) 
    segments = np.column_stack( (xe, ye) )

    if move_nodes :
        for i, mn in enumerate( move_nodes ) :
            ni = mn[0]
            shift = mn[1]
            print( segments[ni,0] )
            segments[ni,0] = segments[ni,0] + shift 
            print( segments[ni,0] )

    if keep_nodes :
        if 0 not in keep_nodes :
            keep_nodes.insert(0, 0)
        if np.size(xe)-1 not in keep_nodes :
            keep_nodes.append( np.size(xe)-1 )
        segments = segments[ keep_nodes ]

    if place_nodes :
        for i, pn in enumerate( np.column_stack( place_nodes ) ) :
            ni = int( pn[0] )
            xn = pn[1]
            yn = np.polyval( np.flip( coefficients ), xn )
            if i >= segments.shape[0]-1 :
                segments = np.vstack( (segments, [xn, yn]) )
            else :
                segments[ni,0], segments[ni,1] = xn, yn

    print( segments )

    # Objetive function
    objf = lambda x, y, segments_x, segments_y :\
        np.sum( ( y - np.interp( x, segments_x, segments_y ) )**2 )

    # Gradient descent optimization
    def gradient_descent( objf, 
                          x, y, segments_x, segments_y, 
                          learning_rate=learning_rate, 
                          max_iterations=10000, 
                          tolerance=1e-6 ):
        
        """Gradient descent optimization."""
        
        residuals = []
        i = 0
        while True:
            gradients = np.zeros_like(segments_y)

            for j in range(len(segments_y)):
                segments_y[j] += learning_rate
                loss_plus = objf(x, y, segments_x, segments_y)
                segments_y[j] -= 2 * learning_rate
                loss_minus = objf(x, y, segments_x, segments_y)
                gradients[j] = (loss_plus - loss_minus) / (2 * learning_rate)
                segments_y[j] += learning_rate
            
            segments_y -= learning_rate * gradients
            
            residuals.append(objf(x, y, segments_x, segments_y))
            
            if i > 0 and abs(residuals[-1] - residuals[-2]) < tolerance:
                print(f'Converged after {i} iterations.')
                break
            i += 1
            if i >= max_iterations:
                print('Maximum number of iterations reached.')
                break

        return segments_y

    # Initial guess for the segmented line
    initial_guess = segments[:,1]

    # Perform gradient descent optimization
    optimized_segments_y = gradient_descent( objf, x, y, 
                                             segments[:,0], 
                                             initial_guess )

    # Combine segments_x with optimized_segments_y
    fitted_segments = np.column_stack((segments[:,0], optimized_segments_y))

    if plot == True :

        plt.figure()
        
        # Create a color map
        colors = cm.rainbow(np.linspace(0, 1, len(fitted_segments)-1))

        ymin = np.min( [ np.min(y), np.min(fitted_segments[:,1] ) ] )
        ymin = ymin - 0.1 * np.abs(ymin)
        for i in range(len(fitted_segments)-1):
            seg = fitted_segments[i:i+2]

            # Calculate the slope of the line
            slope = (seg[1,1] - seg[0,1]) / (seg[1,0] - seg[0,0])

            # Plot the line segment with a unique color
            plt.plot(seg[:,0], seg[:,1], color=colors[i], label=f'Slope: {slope:.2f}')

            # Add a vertical line from the end of the segment to the x-axis, except for the last segment
            if i < len(fitted_segments)-2:
                plt.vlines(seg[1,0], ymin, seg[1,1], linestyles='dashed')

                plt.text( seg[1,0], ymin, f'{seg[1,0]:.2f}', ha='left', va='bottom' )

        # Add a legenid
        plt.legend()

        # Plot the scatter plot with the grayscale color map
        plt.scatter(x, y, color='k', marker='o', facecolors='none')

        ymax = np.max( plt.gca().get_ylim() )
        plt.gca().set_ylim([ymin, ymax])

    return fitted_segments

# -----------------------------------------------------------------------------
def create_log_file( main_function, 
                     main_args=(),
                     file_name='log_file', 
                     add2name='_log',
                     convert_to_pdf=True ):
    """"
    Executes a given function and captures its stdout output, 
    saving it to an HTML log file.
    This function redirects the standard output to capture all print statements and 
    matplotlib figures generated during the execution of the provided main_function. 
    The captured output is then formatted into an HTML file, which includes the 
    execution date and time, and saved with the specified file name.

    Args:
        - main_function (function): The main function to execute and capture output from.

        - file_name (str, optional): The base name of the log file to create. 
            Defaults to 'log_file'. If the provided file name 
            does not have an '.html' extension, it will be replaced with '.html'.
        
        - add2name (str, optional): The string to append to the file name.
        
        - convert_to_pdf (bool, optional): 
            Whether to convert the HTML log file to a PDF file.

    Returns:
        - str: The name of the created HTML log file.
    """

    base_name, ext = os.path.splitext( file_name )
    base_name = base_name + add2name

    # Replace the extension if it is different from .html
    if ext != '.html':
        out_file_name = base_name + '.html'
    else:
        out_file_name = base_name + ext

    class DualOutput:
        def __init__(self):
            self.terminal = sys.stdout
            self.buffer = io.StringIO()

        def write(self, message):
            self.terminal.write(message)
            self.buffer.write(message)

        def flush(self):
            self.terminal.flush()

        def add_figure(self, filename):
            # Control image width via inline CSS
            self.buffer.write(f'<img src="{filename}" alt="{filename}" style="max-width:640px; width:100%;">\n')

    dual_output = DualOutput()
    original_stdout = sys.stdout  # Keep track of the original stdout
    sys.stdout = dual_output

    original_savefig = plt.savefig
    def savefig_wrapper(filename, *args, **kwargs):
        dual_output.add_figure(filename)
        original_savefig(filename, *args, **kwargs)

    plt.savefig = savefig_wrapper

    try:
        main_function(*main_args)  # Execute the main function
    finally:
        output = dual_output.buffer.getvalue()
        sys.stdout = original_stdout  # Reset stdout
        plt.savefig = original_savefig
        current_datetime = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        html_content = f"""
        <html>
        <head><title>{file_name} Output</title></head>
        <body>
        <h1>{file_name} Output</h1>
        <p>Generated on: {current_datetime}</p>
        <pre>{output}</pre>
        </body>
        </html>
        """
        with open(out_file_name, "w") as file:
            file.write(html_content)

    if convert_to_pdf:
        pdf_file_name = base_name + '.pdf'

        try:
            pdfkit.from_file(out_file_name, pdf_file_name)
            out_file_name = pdf_file_name
        
        except Exception as e:
            print(f"Error converting HTML to PDF: {e}")

    return out_file_name