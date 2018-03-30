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

//Reduce SampleRate
Stack.getDimensions(width, height, channels, slices, frames); 
run("Grouped Z Project...", "projection=[Average Intensity] group=6");

//Save
saveAs("Tiff", out_path);

//Close
close('*')

run("Quit");
