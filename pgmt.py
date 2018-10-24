'''
Program to generate a visualization from BEAST files
Created by @evogytis
Modified for ZooPhy by @amagge
'''

from __future__ import print_function
import argparse

from pgmt_utils import *
from pgmt_map_utils import *
from pgmt_spread import Spread

def plot_pgmt(args):
    '''Main method for plotting the map and tree'''
    # Create frames directory
    if not exists(args.output):
        makedirs(args.output)
    output_dirs = join(args.output, "frames")
    if not exists(output_dirs):
        makedirs(output_dirs)
    mplm, frame = get_global_mpl_settings()
    polygons, loihs = get_polygon_and_loihs(args)
    # get data objects
    response = get_aggregation_and_members(args, polygons)
    aggregation, members, pop_centres, polygon_data, polygon_max = response
    print("Loading tree")
    bt_tree, loihs_map = load_tree(args, pop_centres, polygons, loihs)

    print("Getting colours for population centres")
    active_centres = set([x for _, x in loihs_map.items()])
    response = get_colours(aggregation, pop_centres, active_centres)
    group_colours, normalised_coordinates, sorted_group = response

    # Load borders - takes time for many regions
    print("Loading borders. This can take a while.")
    # international_border = load_borders(polygons, aggregation)
    international_border = []

    # Now generate figures 1 and 2
    print("Generating animation")
    # gen_fig_1_2(args, aggregation, pop_centres, group_colours, normalised_coordinates, polygons,
    #             international_border, sorted_group)

    # Now generate figures 3 and 4
    # circle_function = gen_fig_3_4(args, polygons, aggregation, pop_centres, bt_tree,
    #                               group_colours, normalised_coordinates)
    circle_function = get_circle_function(args, polygons, aggregation, pop_centres,
                                          bt_tree, group_colours, normalised_coordinates)

    # Now generate figures 5 and 6
    # migration_function = gen_fig_5_6(args, polygons, aggregation, pop_centres, bt_tree,
    #                                  group_colours, normalised_coordinates, international_border)
    migration_function = get_migration_function(args, polygons, aggregation, pop_centres, bt_tree,
                                                group_colours, normalised_coordinates, international_border)

    # Now generate animation
    gen_animation(args, polygons, aggregation, pop_centres, group_colours, normalised_coordinates,
                  bt_tree, circle_function, polygon_data, polygon_max,
                  migration_function, international_border)

# TODO: Restructure and move implementation to class based loading
# def plot_pmt(args):
#     '''Main method for plotting the map and tree'''
#     # Create output directories if absent
#     if not exists(args.output):
#         makedirs(args.output)
#     output_dirs = join(args.output, "frames")
#     if not exists(output_dirs):
#         makedirs(output_dirs)

#     # First load the coordinates file and determine the best map to load polygons
#     polygons, loihs = get_polygon_and_loihs(args)

#     # Load the spread object
#     # spread = Spread()
#     print("Loading tree")
#     bt_tree, loihs_map = load_tree(args, pop_centres, polygons, loihs)

#     print("Getting colours for population centres")
#     active_centres = set([x for _, x in loihs_map.items()])

def main():
    '''Main method : parse arguments'''
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--tree', type=str, required=True,
                        help="Path to the BEAST tree file")
    parser.add_argument('-m', '--map', type=str, default='auto',
                        help="GeoJSON file to be used as Map")
    parser.add_argument('-l', '--locs', type=str, required=True,
                        help="Location coordinates file (See ReadMe for format)")
    parser.add_argument('-b', '--blst', type=str,
                        default='blacklist.txt',
                        help="Blacklisted locations file")
    parser.add_argument('-o', '--output', type=str, default='output',
                        help="Location where output files should be stored (Default: PWD/output/)")
    args = parser.parse_args()
    print(args)
    plot_pgmt(args)

if __name__ == '__main__':
    main()
