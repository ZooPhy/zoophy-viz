'''
Utility function for pgmt
'''
from __future__ import print_function

import datetime
import imp
import json  # # used for importing JSONs
import sys
import time
import io
import pandas as pd
import matplotlib.path as mpltPath
from os import listdir, makedirs
from os.path import join, exists

import bezier  # # custom arbitrary order Bezier curves
import matplotlib as mpl
import matplotlib.patheffects as path_effects  # # for elegant text
import numpy as np
import unidecode  # # for removing diacritics from example geoJSON
from matplotlib import pyplot as plt
from matplotlib import gridspec  # # for composite figures
from matplotlib.collections import PatchCollection  # # for polygons too
from matplotlib.patches import Polygon  # # for polygons
from scipy.interpolate import \
    UnivariateSpline  # # used to smooth counts of lineages in each location at any given time
from scipy.interpolate import \
    interp1d  # # used to linearly interpolate between data points used in colouring polygons
from sklearn.decomposition import \
    IncrementalPCA  # # used to identify PCA1 when automatically producing a colour map

bt = imp.load_source('baltic', 'baltic.py') ## point to where baltic repo was cloned

def calendarTimeline(start_date, end_date, infmt='%Y-%m-%d', outfmt='%Y-%b',
                     optfmt=None, month_step=1):
    """
    Given two calendar dates returns a list of calendar dates at monthly (by default) intervals.
    """
    current_date = datetime.datetime.strptime(start_date, infmt)
    ending_date = datetime.datetime.strptime(end_date, infmt)
    timeline = []
    while current_date <= ending_date:
        # strftime is buggy with dates before 1900
        if current_date.year >= 1900:
            if optfmt and current_date.month != 1:
                d = datetime.datetime.strftime(current_date, optfmt)
            else:
                d = datetime.datetime.strftime(current_date, outfmt)
            timeline.append(d)
        carry, new_month = divmod(current_date.month - 1 + month_step, 12)
        new_month += 1
        current_date = current_date.replace(year=current_date.year + carry, month=new_month)
    # print(timeline)
    return timeline

def desaturate(rgb,sat,lum):
    if isinstance(rgb,unicode):
        rgb=mpl.colors.hex2color(rgb)
    elif np.array(rgb).shape[-1]==4:
        rgb=rgb[:-1]

    assert 0.0<=sat<=1.0,'Desaturation outside allowed interval'
    assert np.array(rgb).shape[-1]==3,'Incorrect colour array shape'
    hsv=mpl.colors.rgb_to_hsv(rgb)
    hsv[1]=sat
    hsv[2]=lum
    
    return mpl.colors.hsv_to_rgb(hsv)

def get_colours(aggregation, pop_centres, active_centres):
    '''
    Group stuff
    '''
    group_coordinates={}
    normalised_coordinates={}
    for high_order in set(aggregation.values()): ## iterate over high order group
        locs=[loc for loc in pop_centres.keys() if aggregation[loc]==high_order] ## get locations in group
        sorted_locations=sorted(locs,key=lambda loc:pop_centres[loc][0]) ## sort locations in group
        group_coordinates[high_order]=np.mean([pop_centres[loc] for loc in sorted_locations],axis=0)
        # print("WHA", high_order, len(locs), [pop_centres[loc][0] for loc in locs])
        
        if len(locs)>2:
            loc_pca = IncrementalPCA(n_components=2) ## create PCA object
            loc_pca.fit([pop_centres[loc] for loc in locs]) ## fit population centroid coordinates
            loc_transformed=map(tuple,loc_pca.transform([pop_centres[loc] for loc in locs])) ## get PCA coordinates for population centroids
            sorted_transform=sorted(map(tuple,loc_transformed),key=lambda x:x[0])
            
            for l,loc in enumerate(locs):
                normalised_coordinates[loc]=sorted_transform.index(loc_transformed[l])/float(len(locs)-1)
        else:
            for l,loc in enumerate(locs):
                normalised_coordinates[loc]=(0.5+float(l))/float(len(locs))
            
    sorted_group=sorted(list(set(aggregation.values())),key=lambda x:group_coordinates[x][0]) ## sort high-order groups according to longitude

    ipca = IncrementalPCA(n_components=2) ## create PCA object
    ipca.fit([group_coordinates[g] for g in sorted_group]) ## fit population centroid coordinates
    transformed=ipca.transform([group_coordinates[g] for g in sorted_group]) ## get PCA coordinates for population centroids
    XS=transformed[:,0] ## get PCA component 1
    xnorm=mpl.colors.Normalize(0,len(XS)-1) ## create normalization

    # cmap=mpl.cm.viridis ## base colour map
    # import colormaps as cmaps
    # plt.register_cmap(name='viridis', cmap=cmaps.viridis)
    # plt.set_cmap(cmaps.viridis)
    cmap=mpl.cm.get_cmap(name='plasma')
    # cmap=mpl.colo #cmaps.viridis

    group_colours={high_order:None for high_order in sorted_group}
    for group,pos in zip(sorted_group,np.linspace(0.0,1.0,len(sorted_group))): ## iterate over transformed coordinates
        if group in active_centres:
            clist=[]
            clist.append([0.0, desaturate(cmap(pos), 0.3, 0.95)])
            clist.append([0.5, desaturate(cmap(pos), 0.6, 0.8)])
            clist.append([1.0, desaturate(cmap(pos), 0.95, 0.5)])
        else:
            clist=[]
            clist.append([0.0, desaturate(mpl.cm.get_cmap(name='gray')(1), 0.0, 1.0)])
            clist.append([0.5, desaturate(mpl.cm.get_cmap(name='gray')(1), 0.0, 1.0)])
            clist.append([1.0, desaturate(mpl.cm.get_cmap(name='gray')(1), 0.0, 1.0)])
        # print(group, clist)
        group_colours[group]=mpl.colors.LinearSegmentedColormap.from_list(group,clist) ## create colour map
    
    return group_colours, normalised_coordinates, sorted_group

def gen_fig_1_2(args, aggregation, pop_centres, group_colours, normalised_coordinates, polygons,
                international_border, sorted_group):
    '''
    Generate the figures 1 and 2
    '''
    fig,ax = plt.subplots(figsize=(15,14),facecolor='w') ## start figure
    for loc in pop_centres.keys():
        # print("QWA:", loc, "pop:", pop_centres[loc])
        x,y=pop_centres[loc]
        group=aggregation[loc]
        cmap=group_colours[group]
        c=cmap(normalised_coordinates[loc])
        
        ax.scatter(x,y,s=200,facecolor=c,edgecolor='k',zorder=20) ## plot population centre

        effects=[path_effects.Stroke(linewidth=4, foreground='white'),
                    path_effects.Stroke(linewidth=0.5, foreground='k')] ## black text, white outline
        ax.text(x,y,loc,size=10,va='bottom',ha='left',zorder=100,path_effects=effects) ## add its name
        
    for i,loc in enumerate(polygons.keys()): ## iterate over locations
        ax.add_collection(PatchCollection(polygons[loc],facecolor='slategrey',edgecolor='w',lw=1,zorder=1)) ## add polygon

    for stretch in international_border: ## plot borders
        xs,ys=zip(*stretch)
        ax.plot(xs,ys,color='k',zorder=11,lw=2)
        ax.plot(xs,ys,color='w',zorder=10,lw=5)

    ax.spines['top'].set_visible(False) ## suppress axes and their labels
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.tick_params(size=0)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
        
    ax.set_aspect(1)
    # plt.show()
    plt.savefig(args.output+'/1.png')   # save the figure to file

    fig,ax = plt.subplots(figsize=(len(sorted_group),8),facecolor='w') ## start figure
    print(sorted_group)
    for g,group in enumerate(sorted_group):
        for i in np.linspace(0.0,1.0,15):
            ax.scatter(g,i,s=700,color=group_colours[group](i),zorder=20) ## plot a column of points with colours of the colour map

    ax.set_xticks(range(len(sorted_group)))
    ax.set_xticklabels(sorted_group,rotation=90)

    # plt.show()
    plt.savefig(args.output+'/2.png')   # save the figure to file

def gen_fig_3_4(args, polygons, aggregation, pop_centres, bt_tree, group_colours, normalised_coordinates):
    '''
    Generate figure 3 and 4
    '''
    tipSize=20
    branchWidth=2
    locTrait='states' ## trait name for location in tree
    # locTrait='region' ## trait name for location in tree

    plt.figure(figsize=(20,10),facecolor='w') ## start figure
    gs = gridspec.GridSpec(2, 1,height_ratios=[3,1],wspace=0.0,hspace=0.0) ## two row figure
    ax1 = plt.subplot(gs[0]) ## ax1 is tree
    ax2 = plt.subplot(gs[1],sharex=ax1) ## ax2 is a function for lineage presence circle radius
    for k in bt_tree.Objects: ## iterate over branches in tree
        loc=k.traits[locTrait] ## get branch location
        group=aggregation[loc] ## get location's group
        cmap=group_colours[group] ## get group's colour map
        c=cmap(normalised_coordinates[loc]) ## get colour for location based on index of location along PCA1
        
        y=k.y ## y coordinate of branch
        yp=k.parent.y
        
        x=k.absoluteTime ## x coordinate of branch
        xp=k.parent.absoluteTime ## x coordinate of parent branch
        
        if k.branchType=='leaf': ## if leaf
            ax1.scatter(x,y,s=tipSize,facecolor=c,edgecolor='none',zorder=101) ## plot colour tip circle
            ax1.scatter(x,y,s=tipSize*2,facecolor='k',edgecolor='none',zorder=100) ## black outline underneath every tip
            
        elif k.branchType=='node': ## if node
            yl=k.children[0].y ## first child y coordinate
            yr=k.children[-1].y ## last child y coordinate
            
            if xp==0.0:
                xp=x

            ax1.plot([x,x],[yl,yr],color=c,lw=branchWidth,zorder=98) ## plot vertical bar
            
        ax1.plot([xp,x],[y,y],color=c,lw=branchWidth,zorder=98) ## plot branch
        
    ax1.tick_params(size=0,labelsize=0) ## no ticks or tick labels on all axes

    circle_function={loc:None for loc in pop_centres.keys()}

    root_date=bt_tree.Objects[0].absoluteTime ## date of tree root
    last_tip=root_date+bt_tree.treeHeight ## last tip date in tree
    x_grid=np.linspace(root_date,last_tip,100) ## 100 equally spaced points between root and last tip
        
    for loc in pop_centres.keys(): ## for each location
        group=aggregation[loc] ## get group
        cmap=group_colours[group] ## get colour map
        c=cmap(normalised_coordinates[loc]) ## get location's colour
        
        ys=[bt_tree.countLineages(t,condition=lambda x:x.traits[locTrait]==loc) for t in x_grid] ## calculate number of lineages in the location at each time point
        ys=np.array(ys)/20.0 ## scale down
        spl=UnivariateSpline(x_grid,ys) ## use univariate spline...
        spl.set_smoothing_factor(1.1) ## with some smoothing
        
        ys=[sum([k.length for k in bt_tree.Objects if k.traits[locTrait]==loc and k.absoluteTime<t]) for t in x_grid] ## calculate cumulative time spent in location
        ys=np.array(ys)/25.0 ## scale down
        spl=UnivariateSpline(x_grid,ys) ## use univariate spline...    
        
        circle_function[loc]=spl ## assign smoothed function to location
        ax2.plot(x_grid,circle_function[loc](x_grid),color=c) ## plot circle radii over time

    # plt.show()
    plt.savefig(args.output+'/3.png')   # save the figure to file


    plt.figure(figsize=(10,10),facecolor='w') ## start figure
    gs = gridspec.GridSpec(1, 1,wspace=0.0,hspace=0.0)
    ax3 = plt.subplot(gs[0]) ## ax3 is a map showing circle sizes at their maximum extent for each location

    for loc in polygons.keys():
        ax3.add_collection(PatchCollection(polygons[loc],facecolor='w',edgecolor='grey')) ## plot polygons

    for loc in pop_centres.keys(): ## iterate over locations
        circle_radius=circle_function[loc] ## get circle radius function for location
        time_radii=[circle_radius(x) for x in x_grid if circle_radius(x)>0.0] ## calculate the radii for location over time, if they're not 0.0
        
        if len(time_radii)>0: ## if lineages are present
            radius_min,radius_max=min(time_radii),max(time_radii) ## get minimum and maximum radius values
        
            circle_min=plt.Circle(pop_centres[loc],radius=radius_min,alpha=0.5,facecolor='steelblue',edgecolor='k',lw=1,zorder=101)
            circle_max=plt.Circle(pop_centres[loc],radius=radius_max,alpha=0.5,facecolor='indianred',edgecolor='k',lw=1,zorder=100)
            ax3.add_patch(circle_min) ## plot circles at their minimum and maximum
            ax3.add_patch(circle_max)
    ax3.plot() ## make sure polygons will be plotted
    ax3.set_aspect(1) ## equal aspect

    # plt.show()
    plt.savefig(args.output+'/4.png')   # save the figure to file
    return circle_function

def get_circle_function(args, polygons, aggregation, pop_centres, bt_tree, group_colours, normalised_coordinates):
    '''
    Get circle function
    '''
    locTrait='states' ## trait name for location in tree
    # locTrait='region' ## trait name for location in tree

    circle_function={loc:None for loc in pop_centres.keys()}

    root_date=bt_tree.Objects[0].absoluteTime ## date of tree root
    last_tip=root_date+bt_tree.treeHeight ## last tip date in tree
    x_grid=np.linspace(root_date,last_tip,100) ## 100 equally spaced points between root and last tip
        
    for loc in pop_centres.keys(): ## for each location
        group=aggregation[loc] ## get group
        cmap=group_colours[group] ## get colour map
        c=cmap(normalised_coordinates[loc]) ## get location's colour
        
        ys=[bt_tree.countLineages(t,condition=lambda x:x.traits[locTrait]==loc) for t in x_grid] ## calculate number of lineages in the location at each time point
        ys=np.array(ys)/20.0 ## scale down
        spl=UnivariateSpline(x_grid,ys) ## use univariate spline...
        spl.set_smoothing_factor(1.1) ## with some smoothing
        
        ys=[sum([k.length for k in bt_tree.Objects if k.traits[locTrait]==loc and k.absoluteTime<t]) for t in x_grid] ## calculate cumulative time spent in location
        ys=np.array(ys)/25.0 ## scale down
        spl=UnivariateSpline(x_grid,ys) ## use univariate spline...    
        
        circle_function[loc]=spl ## assign smoothed function to location

    return circle_function

def Bezier_control(pointA,pointB,height,frac):
    """ 
    Given a line defined by 2 points A & B, 
    find a third point at a given distance (height) that defines a line perpendicular to line AB which intercepts AB at fraction (frac) along AB.
    Equation derived by Luiz Max Fagundes de Carvalho (University of Edinburgh).
    """
    x1,y1=pointA
    x2,y2=pointB

    sign=1
    if x1>x2:
        sign=-1

    slope = (y2-y1) / (x2-x1)
    d=np.sqrt((y2-y1)**2 + (x2-x1)**2) ## distance between points
    
    h=np.sqrt(height**2+(d*frac)**2) ## distance between desired height and point along line

    n1=x1+h*np.cos(np.arctan(height/float(d)/frac)+np.arctan(slope))*sign ## magic
    n2=y1+h*np.sin(np.arctan(height/float(d)/frac)+np.arctan(slope))*sign

    return (n1,n2) ## return third point's coordinate

def gen_fig_5_6(args, polygons, aggregation, pop_centres, bt_tree, group_colours,
                normalised_coordinates, international_border):
    '''
    Generate figure 5 and 6
    '''
    plt.figure(figsize=(15,18),facecolor='w') ## start figure
    gs = gridspec.GridSpec(2, 1,height_ratios=[4,1]) ## define subplots
    ax1 = plt.subplot(gs[0]) ## map here
    ax2 = plt.subplot(gs[1]) ## plot to show relationship between distances between migration centres and the distance of the control point to the straight line connecting migration centres

    migration_function={i:{j:None for j in pop_centres.keys() if i!=j} for i in pop_centres.keys()} ## matrix of pairs of locations that will contain the Bezier function
    bez_points_dict={i:{j:None for j in pop_centres.keys() if i!=j} for i in pop_centres.keys()}
    np_bez_points_dict={i:{j:None for j in pop_centres.keys() if i!=j} for i in pop_centres.keys()}

    all_distances=[] ## keep track of all distances
    control=lambda d:30*(1/d)**0.5 ## function that will convert the distance between migration points into a distance to be used for finding the control point for a Bezier function

    for i in pop_centres.keys(): ## iterate over locations
        A=pop_centres[i] ## fetch origin
        lon,lat=A
        groupA=aggregation[i] ## get origin's high order group
        c=group_colours[groupA](normalised_coordinates[i]) ## fetch colour map of the group, colour of location is determined by index of location's coordinates along PCA1
        ax1.scatter(lon,lat,s=300,facecolor=c,edgecolor='w',lw=2,zorder=100) ## plot migration centre
        
        for j in pop_centres.keys(): ## iterate over locations again
            if i!=j: ## if not self
                B=pop_centres[j] ## fetch destination
                d=np.sqrt((B[1]-A[1])**2 + (B[0]-A[0])**2) ## calculate distance between location A and location B
                all_distances.append(d) ## remember actual distince
                bez_points=[A] ## Bezier curve will start at point A
                bez_points.append(Bezier_control(A,B,control(d),0.01)) ## head towards the first control point, which is perpendicular to line AB at distance control(d), and 0.01 of the way along line AB
                bez_points.append(Bezier_control(A,B,0.0,0.1)) ## head towards second control point, directly on the line AB, 0.1 of the way along line AB
                bez_points.append(B) ## Bezier curve will finish at point B
                np_bez_points = np.transpose(np.asfortranarray(bez_points))
                bez_points_dict[i][j] = bez_points
                np_bez_points_dict[i][j] = np_bez_points
                curve = bezier.Curve(np_bez_points,degree=2) ## Bezier curve object
                migration_function[i][j]=curve.evaluate_multi ## only interested in calling the evaluate_multi function of the curve, which will return a list of coordinates, given a numpy array of fractions along the line


    for loc in polygons.keys():
        ax1.add_collection(PatchCollection(polygons[loc],facecolor='lightgrey',edgecolor='w',zorder=1)) ## plot location polygons

    Bezier_smooth=35 ## number of coordinates at which to compute the Bezier curve (more = smoother curve)
    eval_Bezier=np.linspace(0.0,1.0,Bezier_smooth) ## numpy array of 35 values going from 0.0 to 1.0

    for locA in pop_centres.keys():
        for locB in pop_centres.keys()[:pop_centres.keys().index(locA)]: ## iterate over pairs of locations
            A=pop_centres[locA]
            B=pop_centres[locB]
            d=np.sqrt((B[1]-A[1])**2 + (B[0]-A[0])**2) ## calculate distance between locations A and B
            sD=sorted(all_distances)
            if locA!=locB and (d in sD[:7] or d in sD[-7:]): ## if locations aren't the same and are some of the longest  or shortest distances - plot away
                migration=migration_function[locA][locB](eval_Bezier) ## compute Bezier curve coordinates along the path
                # Added to fix the bezier error
                migration=[tuple(row) for row in np.transpose(migration)]
                for q in range(len(migration)-1): ## iterate through Bezier line segments with fading alpha and reducing width
                    x1,y1=migration[q] ## coordinates of current point
                    x2,y2=migration[q+1] ## coordinates of next point
                    segL=(q+1)/float(len(migration)) ## fraction along length of Bezier line
                    if aggregation[locA]!=aggregation[locB]: ## locations in different high-order groups
                        fc=group_colours[aggregation[locA]](normalised_coordinates[locA]) ## colour by origin colour
                    else: ## locations in same high-order group
                        fc='k' ## colour black
                        
                    ax1.plot([x1,x2],[y1,y2],lw=7*segL,alpha=1,color=fc,zorder=99,solid_capstyle='round') ## plot actual lineage with width proportional to position along Bezier curve
                    ax1.plot([x1,x2],[y1,y2],lw=10*segL,alpha=1,color='w',zorder=98,solid_capstyle='round') ## plot white outline underneath

    for stretch in international_border: ## plot international border
        xs,ys=zip(*stretch)
        ax1.plot(xs,ys,color='k',zorder=11,lw=2)
        ax1.plot(xs,ys,color='w',zorder=10,lw=5)
        
    ax1.set_aspect(1) ## equal aspect ratio

    ax2.plot(sorted(all_distances),map(control,sorted(all_distances)),ls='--',color='k') ## plot distances between points against their control point distances
    ax2.set_ylabel('control point distance')
    ax2.set_xlabel('distance between migration centres')

    # plt.show()
    plt.savefig(args.output+'/5.png')   # save the figure to file

    # Constants
    smooth=10 ## how many segments Bezier lines will have

    tracking_length=30.0 ## number of days over which to plot the lineage
    depth=tracking_length/365.0 ## convert to decimal years
    departure_condition = lambda f:f-0.5 ## determines how far away (in Bezier fraction) the tail of the migrating lineage is
    transition_point=0.5 ## determines the fixed time point along a branch at which migration happens

    plt.figure(figsize=(15,15),facecolor='w') ## start figure
    gs = gridspec.GridSpec(1, 1) ## define subplots
    ax1 = plt.subplot(gs[0]) ## map here

    for loc in polygons.keys():
        ax1.add_collection(PatchCollection(polygons[loc],facecolor='lightgrey',edgecolor='w',zorder=1)) ## plot location polygons

    time_point=bt.decimalDate('2016-04-01') ## arbitrary time point at which to render the frame
    # locTrait='region' ## trait name for location in tree
    locTrait='states' ## trait name for location in tree
    for k in bt_tree.Objects: ## iterate over branches
        locA=''
        if k.parent.traits.has_key(locTrait):
            locA=k.parent.traits[locTrait] ## get location of parent
        locB=k.traits[locTrait] ## get location of branch
        ## if parent location is present, locations don't match, and branch is intersected
        if locA!='' and locA!=locB and k.parent.absoluteTime<time_point<=k.absoluteTime:
            transition=k.parent.absoluteTime+k.length*transition_point ## migration time
            
            begin_bezier=1-(transition-time_point)/float(depth) ## migrating lineage's head
            end_bezier=departure_condition(begin_bezier) ## migrating lineage's tail
            begin_bezier,end_bezier = np.clip([begin_bezier,end_bezier],0.0,1.0) ## restrict to interval [0,1]
            
            if end_bezier<=1.0 and begin_bezier>=0.0: ## if lineage should still be visible (head departed and tail hasn't arrived)
                points=np.linspace(end_bezier,begin_bezier,smooth) ## get a bunch of points between the start and end of migrating lineage
                bezier_line=migration_function[locA][locB](points) ## get coordinates of the migrating lineage
                # Added to fix the bezier error
                bezier_line=[tuple(row) for row in np.transpose(bezier_line)]
                for q in range(len(bezier_line)-1): ## iterate through Bezier line segments with fading alpha and reducing width
                    x1,y1=bezier_line[q] ## get coordinates for current segment's start
                    x2,y2=bezier_line[q+1] ## get coordinates for current segment's end
                    segL=(q+1)/float(len(bezier_line)) ## fraction along length of Bezier line
                    
                    if aggregation[locA]!=aggregation[locB]: ## locations in different high-order groups
                        fc=group_colours[aggregation[locA]](normalised_coordinates[locA]) ## colour by origin colour
                    else: ## locations in same high-order group
                        fc='k' ## colour black
                    
                    ax1.plot([x1,x2],[y1,y2],lw=7*segL,alpha=1,color=fc,zorder=99,solid_capstyle='round') ## plot actual lineage
                    ax1.plot([x1,x2],[y1,y2],lw=10*segL,alpha=1,color='w',zorder=98,solid_capstyle='round') ## plot underlying white background to help lineages stand out
        
    for stretch in international_border: ## plot international border
        xs,ys=zip(*stretch)
        ax1.plot(xs,ys,color='k',zorder=11,lw=2)
        ax1.plot(xs,ys,color='w',zorder=10,lw=5)
        
    ax1.set_aspect(1) ## equal aspect ratio

    # plt.show()
    plt.savefig(args.output+'/6.png')   # save the figure to file
    return migration_function


def get_migration_function(args, polygons, aggregation, pop_centres, bt_tree, group_colours,
                normalised_coordinates, international_border):
    '''
    Get migration function
    '''
    migration_function={i:{j:None for j in pop_centres.keys() if i!=j} for i in pop_centres.keys()} ## matrix of pairs of locations that will contain the Bezier function
    bez_points_dict={i:{j:None for j in pop_centres.keys() if i!=j} for i in pop_centres.keys()}
    np_bez_points_dict={i:{j:None for j in pop_centres.keys() if i!=j} for i in pop_centres.keys()}

    ## function that will convert the distance between migration points into a distance to be used for finding the control point for a Bezier function
    # control=lambda d:30*(1/d)**0.5 
    # The position of the control point should be a function of the size of the map as well as the two points
    # So, we'll scale as per map size
    scale_factor = -1
    for i in range(2):
        m=min([min([p.xy[:,i].min() for p in polygons[loc]]) for loc in pop_centres.keys()])
        M=max([max([p.xy[:,i].max() for p in polygons[loc]]) for loc in pop_centres.keys()])
        scale_factor = abs(M-m) if scale_factor < abs(M-m) else scale_factor
    print("scale_factor:", scale_factor)
    # TODO: Control this function if the bezier lines looks weird
    control=lambda d:scale_factor*(1/(2*d))**0.5

    for i in pop_centres.keys(): ## iterate over locations
        A=pop_centres[i] ## fetch origin
        for j in pop_centres.keys(): ## iterate over locations again
            if i!=j: ## if not self
                B=pop_centres[j] ## fetch destination
                d=np.sqrt((B[1]-A[1])**2 + (B[0]-A[0])**2) ## calculate distance between location A and location B
                bez_points=[A] ## Bezier curve will start at point A
                bez_points.append(Bezier_control(A,B,control(d),0.01)) ## head towards the first control point, which is perpendicular to line AB at distance control(d), and 0.01 of the way along line AB
                bez_points.append(Bezier_control(A,B,0.0,0.1)) ## head towards second control point, directly on the line AB, 0.1 of the way along line AB
                bez_points.append(B) ## Bezier curve will finish at point B
                np_bez_points = np.transpose(np.asfortranarray(bez_points))
                bez_points_dict[i][j] = bez_points
                np_bez_points_dict[i][j] = np_bez_points
                curve = bezier.Curve(np_bez_points,degree=2) ## Bezier curve object
                migration_function[i][j]=curve.evaluate_multi ## only interested in calling the evaluate_multi function of the curve, which will return a list of coordinates, given a numpy array of fractions along the line
    # sys.exit(0)
    return migration_function

def dd_to_human(date):
    month = int((date-int(date))*12)+1
    month = '0'+str(month) if month<10 else str(month)
    date_human = str(int(date)) + '-' + month + '-01'
    return date_human

def gen_animation(args, polygons, aggregation, pop_centres, group_colours, normalised_coordinates,
                  bt_tree, circle_function, polygon_data, polygon_max,
                  migration_function, international_border):
    t0 = time.time() ## time how long animation takes
    # ----- Global settings -----
    dpi = 90 ## dots per inch for each .png (90 used in the final version)
    every = 6 ## put months labels every number of months
    # locTrait='region' ## trait name for location in tree
    locTrait = 'states' ## trait name for location in tree
    tracking_length = 30.0 ## number of days over which to plot the lineage
    depth = tracking_length/365.0 ## convert to decimal years
    animation_duration = 60.0 ## seconds
    num_frames = 200 # number of frames to generate
    # disable secondary plot until we figure out what to put in there
    disable_chart = True # disable the chart below tree

    # ----- Map settings -----
    edges = 0.03 ## fraction of space around the map
    aspect_ratio = 1 # Aspect ration of the map
    ## determines how far away (in Bezier fraction) the tail of the migrating lineage is
    departure_condition = lambda f:f-0.5
    ## how many segments Bezier lines will have
    smooth = 10

    # ----- Tree settings -----
    tipSize = 20 ## size of tip circles in tree
    branchWidth = 2 ## line width of tree branches
    transition_point=0.5 ## determines the fixed time point along a branch at which migration happens

    # ----- Chart settings -----
    root_date=bt_tree.Objects[0].absoluteTime ## date of tree root
    last_tip=root_date+bt_tree.treeHeight ## last tip date in tree

    # Find start and end date for animation
    start_dd = min([x.absoluteTime for x in bt_tree.Objects])
    end_dd = max([x.absoluteTime for x in bt_tree.Objects])
    # Add 4% time for padding
    end_dd = end_dd + (end_dd-start_dd)*.05
    start = dd_to_human(start_dd)
    end = dd_to_human(end_dd)

    # set animation 
    # set frames based on number of events per time period
    animation_grid = []
    events = [x.absoluteTime for x in bt_tree.Objects]
    events.sort()
    print(len(events), min(events), max(events), max(events)-min(events), "years")
    time_frames = {x:0 for x in set([int(y) for y in events])}
    for event in events:
        time_frames[int(event)] += 1
    if int(start_dd) == int(end_dd):
        animation_grid += list(np.linspace(int(start_dd), int(start_dd)+1, num_frames))
    else:
        for time_period in range(int(start_dd), int(end_dd)):
            frames_in_period = 1
            if time_frames.has_key(time_period):
                frames_in_period += int(time_frames[time_period]*num_frames/len(events))+1
            start_period = bt.decimalDate(dd_to_human(time_period))
            end_period = bt.decimalDate(dd_to_human(time_period+1))
            animation_grid += list(np.linspace(start_period, end_period, frames_in_period))
            # print(time_period, frames_in_period)

    # animation_grid = list(np.linspace(bt.decimalDate(start), bt.decimalDate(end), num_frames))
    print("Animation Grid", animation_grid)
    print('Start of animation: %.2f\nEnd: %.2f'%(min(animation_grid), max(animation_grid)))
    print("", int(start_dd), int(end_dd))
    print('Number of frames to animate: %d'%(len(animation_grid)))
    # print(time_frames)
    # sys.exit()

    ## identify extremes of map based on edges of polygons
    lims=[]
    for i in range(2):
        m=min([min([p.xy[:,i].min() for p in polygons[loc]]) for loc in pop_centres.keys()])
        M=max([max([p.xy[:,i].max() for p in polygons[loc]]) for loc in pop_centres.keys()])
        lims.append([m,M])

    ## width and height of figure is difference in x and y limits
    max_width=35 ## map figure width
    max_height=20 ## map figure width
    map_width, map_height = [int(abs(np.diff(l))) for l in lims]
    ratio = map_width/float(map_height) ## aspect ratio of map
    print("map_width", map_width, "map_height", map_height, "ratio", ratio)
    # print("figwidth", figWidth*2, "figheight", figWidth*ratio)
    
    # Choose layout based on the map layout
    if ratio > 1.5:
        # Map is wide. Use landscape layout
        width = max_width
        height = width*9/16
        plt.figure(figsize=(width, height), facecolor='w')
        gs = gridspec.GridSpec(2, 2, width_ratios=[2, 1], height_ratios=[4,1],
                               wspace=0.0, hspace=0.0)
    elif ratio < 0.75:
        # Map is tall. Use Portrait layout
        height = max_height
        width = height*4/3
        plt.figure(figsize=(width, height), facecolor='w')
        gs = gridspec.GridSpec(2, 2, width_ratios=[2, 3], height_ratios=[4,1],
                               wspace=0.0, hspace=0.0)
    else:  
        # Map is somewhat square. Use square layout
        height = max_height
        width = height*4/3
        plt.figure(figsize=(width, height), facecolor='w')
        gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[4,1],
                               wspace=0.0, hspace=0.0)
    # sys.exit(0)
    ax1 = plt.subplot(gs[0:, 0]) ## ax1 is map
    if disable_chart:
        ax2 = plt.subplot(gs[0:, 1]) ## ax2 is tree
    else:
        ax2 = plt.subplot(gs[0, 1]) ## ax2 is tree
        ax3 = plt.subplot(gs[1, 1]) ## ax3 is tertiary data

    global travelers ## the animation will need to have information to traveling lineages
    travelers=[x for x in bt_tree.Objects if x.parent!=bt_tree.root and x.traits[locTrait]!=x.parent.traits[locTrait]] ## find lineages that have travelled - they're what's going to be animated
    print('\nNumber of travelling lineages: %d (%.3f%% of all lineages)'%(len(travelers),len(travelers)/float(len(bt_tree.Objects))*100))

    for loc in pop_centres.keys():
        circle=plt.Circle(pop_centres[loc],radius=0,label='%s_circle'%(loc),alpha=0.2,facecolor='indianred',edgecolor='k',lw=1,zorder=100) ## circle that tracks the number of lineages at location with radius 0 to begin with
        ax1.add_patch(circle) ## add circle to plot

    for loc in polygons.keys(): ## iterate over locations, plot the initial setup
        poly=PatchCollection(polygons[loc],facecolor='w',edgecolor='grey',lw=1,label='%s_polygon'%(loc),zorder=95) ## location polygon
        ax1.add_collection(poly) ## add location polygon to plot

    for stretch in international_border: ## plot international border
        xs,ys=zip(*stretch)
        ax1.plot(xs,ys,color='k',zorder=97,label='border_bg',lw=2)
        ax1.plot(xs,ys,color='w',zorder=96,label='border',lw=5)

    ax1.spines['top'].set_visible(False) ## remove borders and axis labels
    ax1.spines['right'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.tick_params(size=0) ## no ticks
    ax1.set_xticklabels([]) ## no tick labels
    ax1.set_yticklabels([])
    ax1.set_aspect(aspect_ratio) ## aspect ratio == 1

    ## identify extremes of map based on edges of polygons
    xlims, ylims = lims[0], lims[1]
    xlims = [xlims[0]-map_width*edges,xlims[1]+map_width*edges] ## x axis limits for map
    ylims = [ylims[0]-map_height*edges,ylims[1]+map_height*edges] ## y axis limits for map
    ## set plot limits
    ax1.set_ylim(ylims)
    ax1.set_xlim(xlims)

    xtimeline=calendarTimeline(start,end,'%Y-%m-%d','%Y-%m-%d') ## create timeline from start to end of animation delimited with months
    xpos=[bt.decimalDate(b) for b in xtimeline] ## convert calendar dates to decimal time
    xlabels=[(bt.convertDate(b,'%Y-%m-%d','%b\n%Y') if '-01-' in b else bt.convertDate(b,'%Y-%m-%d','%b')) if (int(b.split('-')[1])+every-1)%every==0 else '' for b in xtimeline] ## month or year-month (if January) tick labels every given number of months

    ################
    ## Secondary plot begins - tertiary data
    ################
    if not disable_chart:
        try:
            for i,group in enumerate(sorted_group): ## iterate through locations
                # import colormaps as cmaps
                # plt.register_cmap(name='viridis', cmap=cmaps.viridis)
                # plt.set_cmap(cmaps.viridis)
                # greyColour=mpl.cm.Greys ## get cmap
                greyColour=mpl.cm.get_cmap(name='Greys') ## get cmap
                c=group_colours[group](0.5) ## fixed colour
                locs=[loc for loc in aggregation.keys() if aggregation[loc]==group] ## get locations present in group
                animation_data=[sum([polygon_data[loc](day) for loc in locs if np.isnan(polygon_data[loc](day))==False]) for day in animation_grid] ## get sum of polygon data for each high order group, at every frame
                
                grey_colour=greyColour((i+1)/float(len(sorted_group)+2)) ## custom grey colour

                ax3.plot(animation_grid,animation_data,lw=3.3,color=grey_colour,zorder=2,label='BW') ## plot tertiary data for group in grey
                ax3.plot(animation_grid,animation_data,lw=3,color=c,zorder=1,label='Colour') ## plot tertiary data for group in colour, but underneath grey lines

            ax3.axvline(xpos[0],color='k',lw=3,label='time',zorder=100) ## add time arrow to indicate current time

            ax3.set_xticks([x+1/24.0 for x in xpos]) ## add ticks, tick labels and month markers
            ax3.set_xticklabels(xlabels) ## set x axis tick labels
            [ax3.axvspan(xpos[x],xpos[x]+1/12.,facecolor='k',edgecolor='none',alpha=0.04) for x in range(0,len(xpos),2)] ## add grey vertical bars for x axis

            ax3.xaxis.tick_bottom() ## make tertiary data plot pretty
            ax3.yaxis.tick_left()
            ax3.spines['top'].set_visible(False)
            ax3.spines['right'].set_visible(False)
            
            ax3.set_xlim(min(xpos),max(xpos)) ## set plot x axis limits
        except NameError:
            pass

    ################
    ## Secondary plot ends - tertiary data
    ################

    ax2.axvline(xpos[0],color='k',lw=3,label='time',zorder=200) ## add time arrow to indicate current time

    ax2.set_xticks([x+1/24.0 for x in xpos]) ## add ticks, tick labels and month markers
    ax2.set_xticklabels(xlabels) ## set x axis labels
    [ax2.axvspan(xpos[x],xpos[x]+1/12.,facecolor='k',edgecolor='none',alpha=0.04) for x in range(0,len(xpos),2)] ## set grey vertical bars for x axis

    ax2.xaxis.tick_bottom() ## make tree plot pretty
    ax2.yaxis.tick_left()
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)

    ax2.tick_params(axis='x',size=0) ## no ticks
    ax2.tick_params(axis='y',size=0)
    ax2.set_xticklabels([]) ## no tick labels
    ax2.set_yticklabels([])

    ax2.set_xlim(min(xpos),max(xpos)) ## axis limits
    ax2.set_ylim(-bt_tree.ySpan*0.01,bt_tree.ySpan*1.01)
    ################
    ## Tertiary plot ends - TREE
    ################

    def animate(frame):
        #### Primary plotting (map)
        ## reset lines (except borders) and texts in the plot
        ax1.lines=[line for line in ax1.lines if 'border' in line.get_label()]
        ax1.texts=[]
        
        if len(animation_grid)-1==frame: ## if at last frame
            next_time=animation_grid[frame] ## current frame is next frame
        else:
            next_time=animation_grid[frame+1] ## next frame
            
        current_time=animation_grid[frame]
        
        effects=[path_effects.Stroke(linewidth=4, foreground='white'),
                    path_effects.Stroke(linewidth=0.5, foreground='k')] ## black text, white outline

        ax1.text(0.05,0.1,'Time: %.3f'%(current_time),size=30,transform=ax1.transAxes,zorder=1000,path_effects=effects) ## add text to indicate current time point
        ax1.text(0.05,0.0,'Created by @evogytis. Modified for ZooPhy by @amagge.',size=18,ha='left',va='bottom',transform=ax1.transAxes,zorder=1000,path_effects=effects) ## "watermark"
        
        exists=[k for k in bt_tree.Objects if k.parent!=bt_tree.root and k.parent.absoluteTime<=current_time<=k.absoluteTime] ## identify lineages that exist at current timeslice

        lineage_locations=[c.traits[locTrait] for c in exists if c.traits[locTrait]!='Not Available'] ## identify locations where lineages are present
        presence=set(lineage_locations) ## all locations where lineages currently are in the tree

        circle=[c.set_radius(0) for c in ax1.patches if '_circle' in c.get_label()] ## reset circle sizes

        for region in presence: ## iterate through every region where a lineage exists
            circle=[c for c in ax1.patches if c.get_label()=='%s_circle'%(region)][0] ## fetch circle of location
            circle.set_radius(circle_function[region](current_time)) ## adjust radius

        for k in travelers: ## iterate through travelling lineages
            transition=k.parent.absoluteTime+k.length*transition_point ## fixed point at which lineage reaches destination
            
            if current_time-depth<transition<current_time+depth: ## make sure transition is within period of animation
                frac=1-(transition-current_time)/float(depth) ## frac will go from 0.0 to 2.0
                ori=k.parent.traits[locTrait] ## fetch locations
                des=k.traits[locTrait]  
                
                if aggregation[ori]!=aggregation[des]:
                    fc=group_colours[aggregation[ori]](normalised_coordinates[ori])
                else:
                    fc='k' ## colour line black unless it's travelling between countries

                lineage_tail=departure_condition(frac) ## Bezier line begins at half a fraction along the path
                lineage_head=frac
                
                lineage_tail,lineage_head=np.clip([lineage_tail,lineage_head],0.0,1.0) ## clip fractions to be within interval (0.0,1.0)
                
                eval_Bezier=np.linspace(lineage_tail,lineage_head,smooth)
                line=migration_function[ori][des](eval_Bezier)
                # Added to fix the bezier error
                line=[tuple(row) for row in np.transpose(line)]
                if lineage_head>0.0 and lineage_tail<1.0: ## only plot if line begins before destination
                    for q in range(len(line)-1): ## iterate through Bezier line segments with fading alpha and reducing width
                        x1,y1=line[q] ## get coordinates for current segment's start
                        x2,y2=line[q+1] ## get coordinates for current segment's end
                        segL=(q+1)/float(len(line)) ## fraction along length of Bezier line
                        ax1.plot([x1,x2],[y1,y2],lw=7*segL,alpha=1,color=fc,zorder=99,solid_capstyle='round') ## plot actual lineage
                        ax1.plot([x1,x2],[y1,y2],lw=10*segL,alpha=1,color='w',zorder=98,solid_capstyle='round') ## plot underlying white background to help lineages stand out

        for i,loc in enumerate(polygons.keys()): ## iterate over locations
            group=aggregation[loc]
            cmap=group_colours[group] ## get cmap
            ## get interpolated tertiary data as fraction of maximum
            frac=polygon_data[loc](current_time)/polygon_max[group]
            c=cmap(frac) ## get colour
            loc_polygons=[p for p in ax1.collections if p.get_label()=='%s_polygon'%(loc)] ## get location polygon
            for polygon in loc_polygons: ## for each component polygon
                polygon.set_facecolor(c) ## change the colour of locations based on tertiary value
        
        update=10 ## update progress bar every X frames
        
        #### Secondary plotting (tree)
        Ls2=[x for x in ax2.lines if 'Colour' not in str(x.get_label())] ## fetch all lines in plot
        partials=[x for x in ax2.lines if 'partial' in str(x.get_label())] ## fetch all tree branches in progress
        finished_lines=[x for x in ax2.lines if 'finished' in str(x.get_label())] ## fetch all tree branches that are finished
        finished_points=[x for x in ax2.collections if 'finished' in str(x.get_label())] ## fetch all tip circles that are finished
        
        finished_labels=[str(x.get_label()) for x in finished_lines]+[str(x.get_label()) for x in finished_points] ## combine everything that's finished (branches + tip circles)
        partial_labels=[str(x.get_label()) for x in partials] ## partially plotted branches
        
        if frame>0 and frame%update==0: ## progress bar
            # clear_output()
            timeElapsed=(time.time() - t0)/60.0
            progress=int((frame*(50/float(len(animation_grid)))))
            percentage=frame/float(len(animation_grid))*100
            rate=timeElapsed/float(frame)
            ETA=rate*(len(animation_grid)-frame)
            sys.stdout.write("[%-50s] %6.2f%%  frame: %5d %10s  time: %5.2f min  ETA: %5.2f min (%6.5f s/operation) %s %s %s" % ('='*progress,percentage,frame,animation_grid[frame],timeElapsed,ETA,rate,len(partials),len(finished_lines),len(finished_points)))
            sys.stdout.flush()

        ####
        ## COMMENT this bit out if you don't want the tree to appear out of the time arrow
        ####
        for ap in bt_tree.Objects: ## iterate over branches
            idx='%s'%(ap.index) ## get unique id of branch
            xp=ap.parent.absoluteTime ## get parent's time

            x=ap.absoluteTime ## get branch's time
            y=ap.y ## get branch's y coordinate

            location=ap.traits[locTrait] ## get branch's location
            group=aggregation[location] ## get location's group
            cmap=group_colours[group] ## get group colour map
            c=cmap(normalised_coordinates[location]) ## get colour
            if xp<=current_time<x: ## branch is intersected by time arrow
                if 'partial_%s'%(idx) in partial_labels: ## if branch was partially drawn before
                    l=[w for w in partials if 'partial_%s'%(idx)==str(w.get_label())][-1] ## get branch line
                    l.set_data([xp,current_time],[y,y]) ## adjust its end coordinate to be time arrow
                else: ## branch is intersected, but not drawn before
                    ax2.plot([xp,current_time],[y,y],lw=branchWidth,color=c,zorder=99,label='partial_%s'%(ap.index)) ## draw branch ending at time arrow, label as partially drawn
                    
            if x<=current_time: ## time arrow passed branch - add it to finished class
                if 'partial_%s'%(idx) in partial_labels: ## if branch has been partially drawn before
                    l=[w for w in partials if 'partial_%s'%(idx)==str(w.get_label())][-1] ## get branch
                    l.set_data([xp,x],[y,y]) ## set end coordinate to be actual end coordinate
                    l.set_label('finished_%s'%(idx)) ## set its label to finished
                    
                if 'finished_%s'%(idx) not in finished_labels: ## branch has not been drawn before at all
                    ax2.plot([xp,x],[y,y],lw=branchWidth,color=c,zorder=99,label='finished_%s'%(ap.index)) ## draw branch, add to finished class
                    
                if 'partial_%s'%(idx) in partial_labels or 'finished_%s'%(idx) not in finished_labels: 
                    if ap.branchType=='leaf': ## if leaf
                        ax2.scatter(x,y,s=tipSize,facecolor=c,edgecolor='none',zorder=102,label='finished_%s'%(ap.index)) ## add tip circle
                        ax2.scatter(x,y,s=tipSize*2,facecolor='k',edgecolor='none',zorder=101,label='finished_%s'%(ap.index)) ## add tip circle outline underneath
                    elif ap.branchType=='node': ## if node
                        yl=ap.children[0].y ## get y coordinates of first and last child
                        yr=ap.children[-1].y
                        ax2.plot([x,x],[yl,yr],lw=branchWidth,color=c,zorder=99,label='finished_%s'%(ap.index)) ## plot vertical bar for node
                    
        for l in Ls2: ## iterate over lines in tree
            if 'time' in l.get_label(): ## if line is time arrow
                l.set_data([current_time,current_time],[0,1]) ## adjust time arrow
        if not disable_chart:
            try:
                Ls3=[x for x in ax3.lines if 'Colour' not in str(x.get_label())] ## fetch all the lines with labels in cases plot
                for l in Ls3:
                    if 'time' in l.get_label(): ## if line is time arrow
                        l.set_data([current_time,current_time],[0,1]) ## adjust time arrow
                    else: ## not time arrow - tertiary plot lines
                        d=l.get_xydata() ## fetch all line data

                        for e in range(len(d)-1): ## iterate over points
                            x_now=d[:,0][e] ## get coordinates of current positions
                            y_now=d[:,1][e]

                            x_nex=d[:,0][e+1] ## cordinates of next position
                            y_nex=d[:,1][e+1]
                            if x_now<current_time: ## time arrow passed x coordinate
                                d[:,0][e]=current_time # adjust x coordinate so it's sitting on top of time arrow
                                d[:,1][e]=y_now+((y_nex-y_now)/(x_nex-x_now))*(current_time-x_now) ## interpolate where new y position is
                        l.set_data(d.transpose()) ## set new data
            except NameError:
                pass
        ## save individual frames for stitching up using 3rd party software (e.g. FFMpeg)
        plt.savefig(args.output+'/frames/ani_frame_%05d.png'%(frame), format='png',bbox_inches='tight',dpi=dpi)

    # regular execution
    for frame in range(len(animation_grid)): ## iterate through each frame
        animate(frame) ## animate will modify the map, tree and cases
        # plt.savefig(args.output+'/frames/ani_frame_%05d.png'%(frame), format='png',bbox_inches='tight',dpi=dpi) ## save individual frames for stitching up using 3rd party software (e.g. FFMpeg)

    # tried multithreading : it doesn't seem to work for now because apparently it has to be launched by the main thread

    print('\n\nDONE!')

    print('\nTime taken: %.2f minutes'%((time.time() - t0)/60.0))

    fps=int(len(animation_grid)/animation_duration) ## frame rate

    # plt.show()
    # plt.savefig(args.output+'/7.png')   # save the figure to file
