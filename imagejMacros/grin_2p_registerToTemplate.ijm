// ImageJ macro to apply moco plugin to stabilize input video.
// Assumes 512x512 resolution input images.
// Args:
// path: string. path to tif containing image stack
// Output:
//   Registered image stack at <path>_registered.tif

//Parse command line argument
args = split(getArgument(), ",");

image_path = args[0]
template_path = args[1]
out_path = args[2]
print(out_path)

//here we go!

//Load the image stack
run("Image Sequence...", "open="+image_path+ " sort");

//Load the template image
open(template_path)
print(out_path)


//Run MOCO
run("moco ", "value=43 downsample_value=1 template=template.tif stack=z0_processed.tif log=None plot=[No plot]");

//Remane Registered Stack
rename("z0_processed_reg.tif")

//Save
saveAs("Tiff", out_path);

//Close
close('*')

run("Quit");
