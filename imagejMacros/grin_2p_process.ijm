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
print(out_path)

//here we go!

//Load the image stack
run("Image Sequence...", "open="+image_path+ " sort");

//Crop
makeRectangle(96, 5, 796, 505);
run("Crop");

//Reduce Images Size
Stack.getDimensions(width, height, channels, slices, frames); 
print(slices)
run("Size...", "width=300 height=250 depth=slices frames constrain average interpolation=Bilinear");

//Generate Template Image
run("Z Project...", "stop=300 projection=[Average Intensity]");
run("16-bit");

//Run MOCO
run("moco ", "value=60 downsample_value=1 template=AVG_z0 stack=z0 log=None plot=[No plot]");

//Save
saveAs("Tiff", out_path);

//Close
close('*')

run("Quit");
