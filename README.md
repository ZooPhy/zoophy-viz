# Generate a video for phylogenetic spread from a BEAST .tree and .coords file
This is a standalone executable version of Curonia notebook created by [@evogytis](https://github.com/evogytis). This has been mainly modified for integration into the ZooPhy pipeline. Check out his [original notebook](https://github.com/blab/baltic/blob/master/curonia.ipynb) if you like having more control over the visualization. Curonia's resulting [video of the West African Ebola virus epidemic](https://www.youtube.com/watch?v=86LCDRRHaJY) won the first place in John Hunter Excellence in Plotting Contest - 2018.

Note:
The original plot contains data from tertiary sources that is currently unavailable in ZooPhy. Hence, we only plot the spread and tree for now.

<span style="color:red;">The repository is currently under development with constant refactoring for better and extendeble code structure. Check back in a few days for a stable version.</span>

# Currently Implemented Enhancements
The original visualization requires a few inputs from the user such as the map, period of animation that needs to be animated for a longer duration, and so on. This project intends to eliminate some of these manual inputs for purposes of using it in a pipeline. For now, we've implemented the following:

1) Map selection: Given the coordinates file, the system chooses the relevant map automatically and trims them if necessary.

2) Frame animations: Not all events are equally interesting, so we automatically load frames based on the concentration of events in the timelines. 

# Future Enhancements
This work is still under development and will require some more features to be implemented for better utility in the future. Essentially, the following are TODOs.

1) More maps at country level: Currently only around 30 maps are made available at the country level, we plan to add the rest from OpenStreetMap.

2) Better coordinate normalization: For the 6 countries across the 180 degree meridian, its ideal to modify the coordinates in the other hemisphere for better display.

3) Better coordinate projection: For a spread containing USA and Japan, its probably more appropriate to show USA on the left and Japan on the right as opposed to the typical mercator projection.

3) Service available as a webpage: It would be good to have the service available as a webpage where users can upload their files and have the results emailed back to them once ready.

# Installation
Firstly, clone this repository.
```
git clone https://github.com/ZooPhy/zoophy-viz.git
```

Following are software requirements for this implmentation.
1) ```python``` (tested with version 2.7)

The system depends on a few python packages that can be installed using the command:
```
cd zoophy-viz
pip install --upgrade -r requirements.txt
```

2) ```FFmpeg``` to stitch images 

On Ubuntu you can install FFmpeg using:
```
sudo apt-get install ffmpeg
```
On mac you can install FFmpeg using:
```
brew install ffmpeg
```

For other OS, check [FFmpeg Documentation](https://www.ffmpeg.org/download.html).

# Running the script
What you need to use this script. 
1) An output file i.e. ```.tree``` file from BEAST. The .tree file needs to have discrete location ```states```.
2) A co-ordinates ```.txt``` file where the location identifier is mapped to the location trait.

```
./gen_pgmt_spread.sh [beast-tree-file] [coordinates-file] [working-directory]
```


Usage for the example included in the repository:
```
./gen_pgmt_spread.sh examples/ebola-liberia/ebola-liberia.tree examples/ebola-liberia/ebola-liberia-coords.txt output/ebola-liberia-plot
```
This should create a video file called ```spread.mp4``` under the ```output/ebola-liberia-plot``` directory.

# Credits

1) Code sources: [Curonia Notebook](https://github.com/blab/baltic/blob/master/curonia.ipynb), [Baltic Trees](https://github.com/blab/baltic/blob/master/baltic.py)

2) Map sources: [OpenStreetMap](https://www.openstreetmap.org/), [SpreaD3](https://github.com/phylogeography/SpreaD3), [Grokbase GeoJSON Maps](https://grokbase.com/t/gg/d3-js/1372gq18j9/geojson-maps), [click_that_hood](https://github.com/codeforamerica/click_that_hood), [Vizzuality](https://github.com/Vizzuality/growasia_calculator)
