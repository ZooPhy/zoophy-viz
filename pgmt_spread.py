'''
Spread object for pgmt
'''
from __future__ import print_function

from pgmt_utils import *

class Location(object):
    '''Class for storing location object in the coordinates file'''
    def __init__(self):
        self.name = None
        self.lon = None
        self.lat = None
        self.polygon = None

class Spread(object):
    '''Class for storing model hyperparameters'''
    def __init__(self):
        self.mplm = None
        self.frame = None
        self.polygons = None
        self.aggregation = None
        self.members = None
        self.pop_centres = None
        self.polygon_data = None
        self.polygon_max = None
        self.bt_tree = None
        self.group_colours = None
        self.normalised_coordinates = None
        self.sorted_group = None
        self.international_border = None
        self.circle_function = None
        self.migration_function = None

    def get_global_mpl_settings(self):
        '''Desc'''
        self.mplm, self.frame = get_global_mpl_settings()

    def get_polygons(self, args):
        '''Desc'''
        self.polygons = get_polygons(args)

    def get_aggregation_and_members(self, args):
        '''Desc'''
        response = get_aggregation_and_members(args, self.polygons)
        self.aggregation, self.members, self.pop_centres, self.polygon_data, self.polygon_max = response

    def load_tree(self, args):
        '''Desc'''
        print("Loading tree")
        self.bt_tree = load_tree(args, self.pop_centres, self.polygons)

    def get_colours(self):
        '''Desc'''
        print("Getting colours for population centres")
        response = get_colours(self.aggregation, self.pop_centres, self.active_centres)
        self.group_colours, self.normalised_coordinates, self.sorted_group = response

    def load_borders(self):
        '''Desc'''
        print("Loading borders. This can take a while.")
        self.international_border = load_borders(self.polygons, self.aggregation)

    def get_circle_function(self, args):
        '''Desc'''
        # Now generate figures 3 and 4
        self.circle_function = get_circle_function(args, self.polygons, self.aggregation,
                                           self.pop_centres, self.bt_tree,
                                           self.group_colours, self.normalised_coordinates)

    def get_migration_function(self, args):
        '''Desc'''
        # Now generate figures 5 and 6
        self.migration_function = get_migration_function(args, self.polygons, self.aggregation,
                                              self.pop_centres, self.bt_tree, self.group_colours,
                                              self.normalised_coordinates,
                                              self.international_border)

    def gen_animation(self):
        '''Desc'''
        # Now generate animation
        gen_animation(self.polygons, self.aggregation, self.pop_centres,
                      self.group_colours, self.normalised_coordinates, self.bt_tree,
                      self.circle_function, self.polygon_data, self.polygon_max,
                      self.migration_function, self.international_border)
