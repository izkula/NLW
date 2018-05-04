// ImageJ macro to apply moco plugin to stabilize input video.
// Assumes 512x512 resolution input images.
// Args:
// path: string. path to tif containing image stack
// Output:
//   Registered image stack at <path>_registered.tif

//Parse command line argument
args = split(getArgument(), ",");

image_path = args[0]
out_path = args[1]
number_images = args[2]
print(out_path)

//here we go!


//Load the image stack
run("Image Sequence...", "open="+image_path+ " sort");
run("8-bit");

//Rename
rename("z")

//Crop
makeRectangle(98, 20, 580, 500);
run("Crop");

//Reduce Images Size
//run("Size...", "width=300 height=250 depth=5666 frames constrain average interpolation=Bilinear");

//Generate Template Image
run("Z Project...", "stop=700 projection=[Average Intensity]");
rename("AVG")


//Run MOCO
run("moco ", "value=43 downsample_value=1 template=AVG stack=z log=None plot=[No plot]");

//Save
saveAs("Tiff", out_path);

//Close
close('*')

run("Quit");
