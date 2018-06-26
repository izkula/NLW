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
print(number_images)
print(out_path)

//here we go!


//Load the image stack
run("Image Sequence...", "open="+image_path+ " sort");
run("8-bit");

//Rename
//rename("z")

getDimensions(width, height, channelCount, sliceCount, frameCount);

//Crop
makeRectangle(98, 20, 580, 480);
run("Crop");


//Reduce Images Size
run("Size...", "width=200 height=166 depth=sliceCount constrain average interpolation=Bilinear");
rename("z")

//Generate Template Image
run("Z Project...", "stop=250 projection=[Average Intensity]");
rename("AVG")

//Run MOCO
run("moco ", "value=30 downsample_value=1 template=AVG stack=z log=None plot=[No plot]");


//Grouped Z Project you must make sure it is an even number stack
run("Grouped Z Project...", "projection=[Average Intensity] group=2");

//Save
saveAs("Tiff", out_path);

//Close
close('*')

run("Quit");
