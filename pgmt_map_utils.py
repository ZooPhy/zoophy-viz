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

WORLD_COUNTRIES = "resources/maps/world.geojson"
WORLD_COUNTRIES_DIR = "resources/maps/countries"
WORLD_CONTINENTS_DIR = "resources/maps/continents"
WORLD_SUBREGIONS_DIR = "resources/maps/subregions"
WORLD_REGIONS_DIR = "resources/maps/regions"

def removeDiacritics(string):
    """
    Removes diacritic marks from unicode.
    """
    output = None
    if isinstance(string, str):
        output = string
    elif isinstance(string, unicode):
        output = unidecode.unidecode(string)
    return output

def get_global_mpl_settings():
    """
    Get global MatPlotLib settings
    """
    typeface = 'Helvetica Neue' ## set default matplotlib font and font size
    mpl.rcParams['axes.labelweight'] = 300
    mpl.rcParams['font.family'] = typeface
    mpl.rcParams['font.weight'] = 100
    # mpl.rcParams['font.weight'] = 'ultralight'
    mpl.rcParams['font.size'] = 22

    frame = '<iframe style="border: 0; width: 400px; height: 472px;" src="https://bandcamp.com/EmbeddedPlayer/album=29809561/size=large/bgcol=333333/linkcol=e99708/artwork=small/transparent=true/" seamless><a href="http://romowerikoito.bandcamp.com/album/nawam-r">NAWAMAR by Romowe Rikoito</a></iframe>'
    return mpl, frame

def get_point_coordinates(args):
    '''
    Loads the coordinates files
    '''
    # Load location of infected hosts
    loihs = {} ## dictionary with point coordinates
    if args.locs:
        print("Loading coordinates from", args.locs)
        pop_file = open(args.locs, 'r')
        for line in pop_file: ## iterate over coordinates
            loc, lat, lon = line.strip().split('\t')
            loihs[loc] = (float(lon), float(lat))
        pop_file.close()
    print("Loihs", loihs)
    return loihs

def load_polygons(map_loc):
    """
    Load polygons for a given file 
    """
    with io.open(map_loc, mode='r', encoding="ISO-8859-1") as data_file:
        json_map = json.load(data_file)

    features = json_map['features']
    location_points = {}  ## location points will be stored here
    polygons = {}         ## polygons will be stored here
    ## key name for each feature
    if 'name' in features[0]['properties']:
        locName = 'name'
    elif 'Name' in features[0]['properties']:
        locName = 'Name'
    elif 'NAME' in features[0]['properties']:
        locName = 'NAME'
    else:
        print("Name property not found in GeoJSON.")
        sys.exit(0)

    for loc in features: ## iterate through features (locations)
        if loc['geometry']:
            # print(loc['geometry'])
            poly = np.asarray(loc['geometry']['coordinates']) ## get coordinates
            ## standardised location name (remove diacritics)
            # print("Loading", loc['properties'][locName])
            location = removeDiacritics(loc['properties'][locName])
            polygons[location] = []
            location_points[location] = []
            if loc['geometry']['type'] == 'MultiPolygon': ## multiple parts detected
                for part in np.asarray(poly): ## iterate over each component polygon
                    for coords in np.asarray(part): ## iterate over coordinates
                        coords = np.array(coords)
                        xs = coords[:, 0] ## longitudes
                        ys = coords[:, 1] ## latitudes
                        # TODO: There are six countries that are both in Eastern and Western
                        # Hemisphere's across the 180th meridian. Plotting those may be a pain.
                        # We should be normalizing them based on the country's polygons and 
                        # other polygons in the spread, but for now ignore such polygons.
                        # Ignore polygons in Alaska that in the eastern hemisphere
                        if location == "Alaska":
                            if any(i > 0 for i in xs):
                                break
                        # Ignore polygons in Russia that in the western hemisphere
                        if location in ["Russia", "Russian Federation"]:
                            if any(i < 0 for i in xs):
                                # print("Skipping", len(xs), len(ys), xs)
                                break
                        ## append coordinates to location's list of coordinates
                        location_points[location].append(np.vstack(zip(xs, ys)))
            if loc['geometry']['type'] == 'Polygon': ## location is single part
                for coords in np.asarray(poly): ## iterate over coordinates
                    coords = np.array(coords)
                    xs = coords[:, 0] ## longitudes
                    ys = coords[:, 1] ## latitudes
                    ## append coordinates to location's list of coordinates
                    location_points[location].append(np.vstack(zip(xs, ys)))
            complete_location = []
            ## iterate and create a polygon for each component of a location
            for part in location_points[location]:
                complete_location.append(Polygon(part, True))
            ## assign list of polygons to a location
            polygons[location] = complete_location
    # print('%d polygons loaded:\n%s'%(len(polygons.keys()), polygons.keys()))
    # sys.exit(0)
    return polygons

def get_loc_dist(map_loc, loih_list, points_list):
    """
    Get the polygons that have the location of infected hosts's coordinates
    """
    polygons = load_polygons(map_loc)
    # for loc_name, polygon in polygons.items():
    #     print(loc_name, len(polygon))
    loihs_map = {}
    found = np.array([False for _ in range(len(loih_list))])
    for location, multi_polygon in polygons.items():
        inside = np.array([False for _ in range(len(loih_list))])
        for polygon in multi_polygon:
            points_inside = polygon.get_path().contains_points(points_list)
            inside = inside | np.array(points_inside)
        found = found | inside
        locs_found = [loih_list[i] for i,x in enumerate(inside) if x]
        if len(locs_found) > 0:
            perc_found = len(locs_found)*100/len(loih_list)
            print(str(perc_found)+"%", "of locations found in", location, "in", map_loc)
            loihs_map[location] = locs_found
    return loihs_map

def get_best_map(args, loihs):
    """
    Choose appropriate map based on location of infected hosts
    """
    loih_list, points_list = [], []
    for loih, points in loihs.items():
        loih_list.append(loih)
        points_list.append(points)
    # print(loih_list, points_list)
    print("Choosing best map based on coordinates")
    map_loc = None

    # First scan world by countries and check if almost all of the
    # coordinates are in a single country
    # Get the distribution of loihs by country polygons
    loihs_map = get_loc_dist(WORLD_COUNTRIES, loih_list, points_list)
    if len(loihs_map) ==  1:
        # if there is just one country try to find the country's map
        loihs_polygon = loihs_map.keys()[0]
        print("All locations were found in", loihs_polygon)
        country_map = join(WORLD_COUNTRIES_DIR, loihs_polygon+".geojson")
        if exists(country_map):
            map_loc = country_map

    # if not found try to find a subregion that fits all loihs
    if not map_loc:
        print("Searching in Subregions")
        for geomap in listdir(WORLD_SUBREGIONS_DIR):
            loihs_map = get_loc_dist(join(WORLD_SUBREGIONS_DIR, geomap), loih_list, points_list)
            loihs_set = set([y for _, x in loihs_map.items() for y in x])
            if len(loihs_set) == len(loih_list):
                print("All locations were found in", geomap)
                map_loc = join(WORLD_SUBREGIONS_DIR, geomap)

    # if not found try to find a continent that fits all loihs
    if not map_loc:
        print("Searching in Continents")
        for geomap in listdir(WORLD_CONTINENTS_DIR):
            loihs_map = get_loc_dist(join(WORLD_CONTINENTS_DIR, geomap), loih_list, points_list)
            loihs_set = set([y for _, x in loihs_map.items() for y in x])
            if len(loihs_set) == len(loih_list):
                print("All locations were found in", geomap)
                map_loc = join(WORLD_CONTINENTS_DIR, geomap)

    # if not found try to find a region that fits all loihs
    if not map_loc:
        print("Searching in Regions")
        for geomap in listdir(WORLD_REGIONS_DIR):
            loihs_map = get_loc_dist(join(WORLD_REGIONS_DIR, geomap), loih_list, points_list)
            loihs_set = set([y for _, x in loihs_map.items() for y in x])
            if len(loihs_set) == len(loih_list):
                print("All locations were found in", geomap)
                map_loc = join(WORLD_REGIONS_DIR, geomap)

    trim_maps = False
    if not map_loc:
        map_loc = WORLD_COUNTRIES
        trim_maps = True
    print("Loading polygons using map:", map_loc)

    loihs_map = get_loc_dist(WORLD_COUNTRIES, loih_list, points_list)
    polygons = load_polygons(map_loc)

    # Trim maps to focus well and avoid distraction
    if trim_maps:
        # For now just exclude the polygons if they are not part of the migrations
        # TODO: if world map, exclude polygons outside min and max limits
        for loc, _ in polygons.items():
            if not loc in loihs_map.keys():
                print("removing", loc)
                del polygons[loc]
    # sys.exit(0)
    return polygons


def get_polygon_and_loihs(args):
    """
    Get the polygon from the specified geoJSON map file
    """
    loihs = get_point_coordinates(args)
    if not args.map or (args.map and args.map == "auto"):
        polygons = get_best_map(args, loihs)
    else:
        print("Loading map", args.map)
        polygons = load_polygons(args.map)
    
    # remove blacklisted polygons
    if args.blst:
        blacklist_file = open(args.blst, 'r')
        for line in blacklist_file: ## iterate over spreadsheet of coordinates
            if line[0] != "#": ## if country matches a location for which polygon is available
                bl_loc = line.strip()
                if polygons.has_key(bl_loc): ## ignore blacklisted locations
                    print("Ignoring '", bl_loc, "' specified in Blacklist")
                    del polygons[bl_loc]
    # sys.exit(0)
    return polygons, loihs

def get_aggregation_and_members(args, polygons):
    '''
    Gets the aggregation object and members
    '''
    pop_centres = {}
    # pop_centres = get_point_coordinates(args, polygons)

    for loc in polygons.keys(): ## iterate through locations
        ## if capital coordinates were not found during the previous pass
        if pop_centres.has_key(loc) == False:
            ## assign pop centroid as mean of all polygon coordinates
            pop_centres[loc] = np.mean([np.mean(p.get_verts(), axis=0) for p in polygons[loc]],
                                       axis=0)
            #locations.append(loc)
    print('migration centres:\n%s\n'%(pop_centres.keys()))
    # print('\ndictionary format for migration centre coordinates:\n%s\n'%(pop_centres))

    # TODO: Load aggregation from file
    ## no higher-order aggregation (each location maps to itself)
    aggregation = {centre:centre for centre in pop_centres.keys()}
    assert len(aggregation.keys())>=len(pop_centres.keys()),'Not all pop_centres locations were assigned to a higher order group in aggregation dictionary'

    members={group:[loc for loc in aggregation.keys() if aggregation[loc]==group] for group in set(aggregation.values())}

    # print('dictionary format for aggregation of locations:\n%s'%(aggregation))

    # minValue, maxValue = 0, 1000 ## minimum and maximum observed values
    # minObs, maxObs = 2, 20 ## minimum and maximum observations per location

    # Load tertiary data if any
    loc_data={loc:{} for loc in pop_centres.keys()} ## no tertiary data
    ## create random tertiary data for each polygon
    # loc_data={loc:{'%d-%02d-%02d'%(np.random.randint(2011,2018), np.random.randint(1,12), np.random.randint(1,28)):np.random.randint(minValue,maxValue) for t in range(np.random.randint(minObs,maxObs))} for loc in pop_centres.keys()}

    # print('polygon data input dictionary:\n%s\n'%(loc_data))

    data_timelines = {} ## will contain a sorted list of time points for which data are available in a location
    polygon_data = {} ## will contain tertiary data available for each location, in their temporal order

    for loc in loc_data.keys():
        ## for each location remember the sorted timeline for which data are available
        data_timelines[loc] = sorted(loc_data[loc].keys(), key=lambda w:bt.decimalDate(w))
        
        if len(data_timelines[loc]) > 1: ## if there are entries
            ## for each location generate an interpolation function that will be able to return a number given a time point
            x = map(bt.decimalDate,data_timelines[loc])
            y = [loc_data[loc][day] for day in data_timelines[loc]]
            polygon_data[loc] = interp1d(x, y, bounds_error=False)
        else: ## not enough entries for interpolation - function will return 0.0 whenever called
            polygon_data[loc]=lambda d:0.0

    polygon_max={} ## will contain the maximum value observed in a given high order group at any point
    for group in set(aggregation.values()): ## iterate over higher order group designation
        if max([len(loc_data[loc].values()) for loc in members[group]])>0: ## if group has sufficient data
            polygon_max[group]=float(max([max(loc_data[loc].values()) for loc in members[group]])) ## for each higher-order group identify the maximum existing value (used to produce numbers in range [0,1] for colour maps later on)
        else:
            polygon_max[group]=1.0 ## insufficient data - maximum is set to 1.0

    print('maximum polygon data values:\n%s\n'%(polygon_max))
    return aggregation, members, pop_centres, polygon_data, polygon_max

def load_tree(args, pop_centres, polygons, loihs):
    '''
    Load the tree file i.e. the labeled phylogeny
    '''
    if args.tree:
        tree_file = args.tree
        print("Tree param", tree_file[-4:])
        if tree_file[-4:] == "tree":
            print("Loading tree file", tree_file)
            # treeFile=open('baltic_example.json','r')
            # bt_tree = bt.loadNexus(tree_file,tip_regex='\_([0-9\-]+)$',verbose=True) ## import tree from nexus
            # bt_tree = bt.loadNexus(tree_file,verbose=True)
            bt_tree = bt.loadNexus(tree_file,tip_regex=r'_([0-9]{4}\.[0-9]{0,10})_', verbose=False)
        elif tree_file[-4:] == "json":
            print("Loading json file", tree_file)
            translation={'name':'strain','height':'tvalue'} ## define translations between JSON and baltic
            bt_tree=bt.loadJSON(tree_file,translation) ## import tree from JSON
            highestTip=max([bt.decimalDate(k.traits['date']) for k in bt_tree.Objects if k.branchType=='leaf']) ## find highest tip
            bt_tree.setAbsoluteTime(highestTip) ## set absolute time for each branch in tree
        else:
            print("Unknown format for the tree file. Exiting.")
            sys.exit(0)
    else:
        print("Unknown tree  file. Exiting.")
        tree_file = 'baltic_example.json'
    # print("\n\n", [k.traits.keys() for k in bt_tree.Objects], "\n\n")
    
    traits_set = set()
    for k in bt_tree.Objects:
        for key in k.traits.keys():
            traits_set.add(key)
    print("\n\nTraits found:", traits_set, "\n")

    # locTrait='region' ## trait name for location in tree
    locTrait='states' ## trait name for location in tree
    tree_locs=list(set([k.traits[locTrait] for k in bt_tree.Objects])) ## fetch all trait states in tree
    print("\ntree_locs", tree_locs)
    
    # np.random.shuffle(tree_locs) ## randomise order of trait states from the tree
    # pretend={t:m for t,m in zip(tree_locs,pop_centres.keys())} ## assign each tree state to a polygon name

    # loihs = get_point_coordinates(args)
    loih_list, points_list = [], []
    for loih, points in loihs.items():
        loih_list.append(loih)
        points_list.append(points)
    # print(loih_list, points_list)
    loihs_map = {}
    found = np.array([False for _ in range(len(loihs))])
    for location, multi_polygon in polygons.items():
        inside = np.array([False for _ in range(len(loih_list))])
        for polygon in multi_polygon:
            points_inside = polygon.get_path().contains_points(points_list)
            inside = inside | np.array(points_inside)
        found = found | inside
        locs_found = [loih_list[i] for i,x in enumerate(inside) if x]
        for loc in locs_found:
            print("Location", loc, "found in area", location)
            # TODO: implement cases where multiple points are allowed in a polygon
            loihs_map[loc] = location
    # Ideally all locations should be found, exceptions are usually islands
    # When not found, assign to the polygon with the nearest centroid
    locs_not_found = [loih_list[i] for i,x in enumerate(found) if not x]
    for loc in locs_not_found:
        print("Location", loc, "was not found inside any area loaded")
        min_dist = float('Inf')
        for location, centroid in pop_centres.items():
            dist = plt.mlab.dist(loihs[loc], centroid)
            if dist < min_dist:
                min_dist = dist
                loihs_map[loc] = location
        print("Location", loc, "was found close to", loihs_map[loc])
    ## iterate over every branch in tree
    ## recode tree trait states into those from the map
    for k in bt_tree.Objects:
        k.traits[locTrait]=loihs_map[k.traits[locTrait]]

    tree_locs=set([k.traits[locTrait] for k in bt_tree.Objects]) ## get current set of trait states in tree
    print('\n%d tree locations: %s'%(len(tree_locs),', '.join(tree_locs)))
    print('\n%d map locations (defined migration centres): %s'%(len(pop_centres.keys()),', '.join(pop_centres.keys())))
    assert tree_locs.issubset(set(pop_centres.keys())),'Coordinates for state(s) %s present in the tree not provided'%(', '.join(tree_locs.difference(set(pop_centres.keys())))) ## every state present in the tree has to correspond to a set of coordinates
    # sys.exit(0)
    return bt_tree, loihs_map

def load_borders(polygons, aggregation):
    '''
    Load the borders
    '''
    print("\nLoading international borders from ", len(polygons), "polygons and", len(aggregation), "aggregations")
    border_sharing={locA:{locB:False for locB in polygons} for locA in polygons}
    international_border=[]
    for i,locA in enumerate(polygons): ## iterate over polygons
        for j,locB in enumerate(polygons): ## iterate over polygons again
            if i>j and aggregation[locA]!=aggregation[locB]: ## only do one pairwise comparison
                # print(i,j,locA,locB)
                for partA in polygons[locA]: ## iterate over constituent polygons of location A
                    setA=set(map(tuple,partA.get_verts())) ## convert coordinates of polygon into a set of tuples
                    for partB in polygons[locB]:  ## iterate over constituent polygons of location B
                        setB=set(map(tuple,partB.get_verts())) ## convert coordinates of polygon into a set of tuples
                        temp=[] ## temporary storage
                        if setA.isdisjoint(setB)==False: ## if polygon A and polygon B share coordinates
                            for point in map(tuple,partA.get_verts()): ## for each point in polygon A
                                if point in setA.intersection(setB) and point not in temp: ## if it is also present in polygon B and not seen yet - add to temporary list
                                    temp.append(point)
                                    
                        if len(temp)>0: ## if overlap between polygon A and B is not zero
                            international_border.append(temp) ## add to list of international border coordinates
    distance=lambda p:abs(p[0][0]-p[1][0])+abs(p[0][1]-p[1][1]) ## function to get Euclidean distance
    for s,stretch in enumerate(international_border): ## iterate over every segment of international border
        distances=[distance((stretch[i],stretch[i-1])) for i in range(len(stretch))] ## get distances between consecutive points
        M=max(distances) ## get largest distance
        idx=distances.index(M) ## identify where the largest distance is
        if idx!=0: ## if the largest distance is not between the first and last point of the border stretch... 
            international_border[s]=international_border[s][idx:]+international_border[s][:idx] ## reorganise the stretch such that the largest distance occurs between first and last points of the stretch
    print("loaded %s international borders"%(len(international_border)))
    return international_border
