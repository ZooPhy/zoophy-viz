#!/bin/sh

# Path to the tree file generated by BEAST
tree_file=$1
# Path to the coordinates file create by zoophy
coords_file=$2
# Job ID created by zoophy. A directory will be created
# with this name and will be used for generating frames
# and creating the video file
work_dir=$3

# Path to the maps file
maps_file="auto"

if [ "$4" ]
then
    maps_file=$4
fi


# output_dir=output/$work_dir
output_dir=$work_dir

# Generate the frames for the spread
python -u pgmt.py -t $tree_file -l $coords_file -o $output_dir -m $maps_file > $output_dir-viz.log 2>&1


# Stitch the frames into a video
fps=15
frames_dir="$output_dir/frames"
output_file="$output_dir/spread.mp4"
scale="2160:trunc(ow/a/2)*2"

ffmpeg -y -framerate $fps -start_number 0 -i $frames_dir/ani_frame_%05d.png -pix_fmt yuv420p -b:a 64k -vf scale=$scale $output_file
