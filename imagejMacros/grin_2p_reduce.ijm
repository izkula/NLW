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
open(image_path);

//Reduce Images Sampling
run("Reduce...", "reduction=3");

//Save
saveAs("Tiff", out_path);

//Close
close('*')

run("Quit");
