// ImageJ macro to apply moco plugin to stabilize input video.
// Assumes 512x512 resolution input images.
// Args:
// path: string. path to tif containing image stack
// Output:
//    A template image to register all sessions

//Parse command line argument
args = split(getArgument(), ",");

image_path = args[0]
out_path = args[1]
print(out_path)

//here we go!

//Load the image stack
run("Image Sequence...", "open="+image_path+ " sort");


//Generate Template Image
run("Z Project...", "stop=300 projection=[Average Intensity]");

rename("template.tif");
//Save
saveAs("Tiff", out_path);

//Close
close('*')

run("Quit");
