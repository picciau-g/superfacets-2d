Executable name: Superfacets

It produces an oversegmentation of a triangle mesh (considered as its surface) into small patches, according to a distance metric which combines an approximation of the geodesic distance on the surface and an angular term as approximation of the surface curvature. It also allows for the visualization of the output segmentation applied to the object, and to save a screenshot of this result.

Compile:
In the folder there are source files and headers, and a .pro file. Under Linux-based systems it can be compiled by typing 'make' into the directory where the files are stored.

Execute:
The program is triggered from command line by typing './Superfacets': if no parameter is passed, then it just shows instructions about how to pass the parameters which are necessary to run the segmentation, and then exits.

There are several parameters which can be passed to the program, here we mention the main ones:
	'-m' [path/to/meshfile] loads the mesh stored in path/to/meshfile. Notice that it can be one between a '.off' or a '.tri' file.
	
	'-a' [alpha] states the weight of the angular term over the angular one (higher alpha means higher importance to the curvature). In our experiments, we have established a good value for this parameter in the range [0...500]
	
	'-nseg' [number of segments] allows to specify how many regions we want in the output segmentation
	
	'-r' [radius] is mutually exclulsive with the '-nseg' option, since it allows to specify the approximate size of each segment (expressed as a percentage of the bounding box diagonal)
	
	'-flood' [initmode] if the approximate size of a segment is specified, this parameter allows to choose between a Dijkstra-based initialization (initmode=1) and a regular subdivision of the volume enclosed by the mesh bounding box (initmode=0)
	
	'-out' [path/to/output/file] saves the output segmentation (the index of the cluster index of each face per line) to file path/to/output/file
	
	'-vis' opens a window which display a visual output of the segmentation once it is done. If we press the key 's' when this window is open, the image currently showed is saved to a file
	
	'-ov' [path/to/segmentation/file], where 'path/to/output/file' is an already computed segmentation file, only opens a window to display the segmentation without computing a new one
	
	
Examples of usage

./Superfacets -m "mesh.off" -nseg 150 -a 200 -out "output.seg" -vis
	subdivides model "mesh.off" in 150 surface patches with weigth of angular term = 200, then saves the result to file "output.seg" and visualizes it
	
./Superfacets -m "mesh.off" -ov "output.seg"
	loads model "mesh.off" and applies segmentation in "output.seg" to it (it must have been computed on the same model)
	
	
