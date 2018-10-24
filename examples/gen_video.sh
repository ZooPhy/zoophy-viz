#!/bin/sh

fps=15
output_dir=$1
frames_dir="$output_dir/frames"
output_file="$output_dir/spread.mp4"
scale="2160:trunc(ow/a/2)*2"

ffmpeg -y -framerate $fps -start_number 0 -i $frames_dir/ani_frame_%05d.png -pix_fmt yuv420p -b:a 64k -vf scale=$scale $output_file
