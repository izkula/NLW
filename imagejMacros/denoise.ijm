// ImageJ macro to apply denoise input video and save as multipage tiff
// Assumes 512x512 resolution input images.
// Args:
//   path_denoised: string. path to tif containing denoised image stack
// Output:
//   Denoised image stack at <path_denoised>_denoise.tif
// Example usage (mac):
//   /Applications/Fiji.app/Contents/MacOS/ImageJ-macosx  -macro /Users/samueljyang/Dropbox/dump/fiji_macros/turboreg.ijm "/Users/samueljyang/Dropbox/dump/oeg/denoised_small.tif"

// Parse this command line argument
args = split(getArgument()," "); 
path_denoised =args[0];


// Set this path manually for debugging in ImageJ or FIJI
//  path_denoised="/Users/samueljyang/Dropbox/dump/oeg/denoised_small.tif";


// Load the denoised image stack
//open(path_denoised);

out_path = args[1];
print(path_denoised);
print(out_path);


//run("Image Sequence...", "open="+path_denoised+" sort");
//////run("Image Sequence...", "open="+path_denoised+" file='.tif' sort");
run("Image Sequence...", "open="+path_denoised+" file='.tif' sort");


run("PureDenoise ", "parameters='3 1' estimation='Auto Global' ");
rename("denoised");

saveAs("Tiff", out_path+"/denoised.tif")


run("Quit");