// ImageJ macro to apply moco plugin to stabilize input video.
// Assumes videos need to be cropped to 512 x 512. 
// Args:
//   path_denoised: string. path to tif containing denoised image stack
// Output:
//   Registered image stack at <path_denoised>_registered.tif
// Example usage (mac):
//   /Applications/Fiji.app/Contents/MacOS/ImageJ-macosx  -macro /Users/samueljyang/Dropbox/dump/fiji_macros/turboreg.ijm "/Users/samueljyang/Dropbox/dump/oeg/denoised_small.tif"

// Parse this command line argument
args = split(getArgument(),","); 
path_denoised =args[0];
print(path_denoised);

// Set this path manually for debugging in ImageJ or FIJI 
//  path_denoised="/Users/samueljyang/Dropbox/dump/oeg/denoised_small.tif";
// Load the denoised image stack
// open(path_denoised);

out_path = args[1];
print(out_path);

template_path = args[2]
print(template_path)

first_frame = args[3]
num_frames = args[4]
downsample_val = args[5]; // 1/downsampleFactor
print(parseInt(first_frame));
print(parseInt(num_frames));
print(parseFloat(downsample_val));


run("Image Sequence...", "open="+path_denoised+ " number=" + parseInt(num_frames) + " starting=" + parseInt(first_frame) + " sort");
//run("Image Sequence...", "open="+path_denoised+" sort");
//run("Image Sequence...", "open="+path_denoised+" file='.tif' sort");

makeRectangle(100, 0, 512, 512);
run("Crop");

run("Scale...", "x=.5 y=.5 z=1.0 width=256 height=256 depth=" + parseInt(num_frames) + " interpolation=Bilinear average process create title=downsamp");
selectWindow("downsamp");

//run("PureDenoise ", "parameters='3 1' estimation='Auto Global' ");
//rename("denoised");
//saveAs("Tiff", out_path+"/denoised.tif")

// batch mode is faster by not displaying images
/////setBatchMode(true); //It appears as if this has to come *after* the denoising...

////open(out_path+"/denoised.tif");

for(i=1; i<=nSlices;i++) {
	setSlice(i);
	setMetadata("Label",i-1);
}

rename("denoised");

//Stack.getDimensions(width, height, channels, slices, frames); 


// Use the median image as the reference (target) image for alignment
/*
path_median= path_denoised +"_median.tif";
run("Z Project...", "start=1 stop=" + parseInt(num_frames) + " projection=Median");
*/

print(template_path);
//run("Image Sequence...", "open="+template_path+" file='.tif' sort");
//run("Z Project...", "start=1 stop=1 projection=Median");
open(template_path)


run("8-bit");
rename("median");

selectWindow("denoised");
run("8-bit");

run("moco ", "value=102 downsample_value=1 template=median stack=denoised log=None plot=[No plot]");

//run("Scale...", "x=1 y=1 z=0.333333 width=256 height=256 depth=" + parseInt(num_frames)/3 + " interpolation=Bilinear average process create title=downsamp");
run("Scale...", "x=1 y=1 z=" + downsample_val + " width=256 height=256 depth=" + parseInt(num_frames)*downsample_val + " interpolation=Bilinear average process create title=downsamp");

path_output = out_path+"/registered_noDenoise.tif";
saveAs("Tiff", path_output);
rename("registered");

path_median= path_denoised +"_median.tif";
run("Z Project...", "start=1 stop=110 projection=Median");
rename("median");
path_output = out_path+"/median.tif";
saveAs("Tiff", path_output);

selectWindow("registered");
run("Z Project...", "projection=[Average Intensity]");
rename("mean");
path_output = out_path+"/mean.tif";
saveAs("Tiff", path_output);

run("Quit");