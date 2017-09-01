// ImageJ macro to apply turboreg plugin to stabilize input video.
// Assumes 512x512 resolution input images.
// Args:
//   path_denoised: string. path to tif containing denoised image stack
// Output:
//   Registered image stack at <path_denoised>_registered.tif
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
print(out_path);


//run("Image Sequence...", "open="+path_denoised+" sort");
run("Image Sequence...", "open="+path_denoised+" file='.tif' sort");


run("PureDenoise ", "parameters='3 1' estimation='Auto Global' ");
rename("denoised");

saveAs("Tiff", out_path+"/denoised.tif")

// batch mode is faster by not displaying images
setBatchMode(true); //It appears as if this has to come *after* the denoising...

//while (nImages>0) { 
//	selectImage(nImages); 
//    close(); 
//} 


////open(out_path+"/denoised.tif");

for(i=1; i<=nSlices;i++) {
	setSlice(i);
	setMetadata("Label",i-1);
}

rename("denoised");

Stack.getDimensions(width, height, channels, slices, frames); 


// Use the median image as the reference (target) image for alignment
path_median= path_denoised +"_median.tif";
run("Z Project...", "start=1 stop=110 projection=Median");
rename("median");


// Convert denoised image stack into individual frames
selectWindow("denoised");

run("Stack to Images");


// Run TurboReg to align each of the denoised frames to the reference frame
for(i=0;i < slices; i++) {
  print(i);
  window_median="median";
  window_denoised=i;
  run("TurboReg ", "-align -window "+window_denoised+" 0 0 511 511 -window "+window_median+" 0 0 511 511 -translation 256 256 256 256 -showOutput");
  run("Make Substack...", "  slices=1");
  close("Output*");
  rename(window_denoised+"_registered");
}
//setBatchMode(false); // TODO: figure out how to use these properly
//updateDisplay();

// Combine registered frames into an image stack
run("Images to Stack", "name=Stack title=[] use");
// The first ~half of the images are the input images; select only the registered images
run("Make Substack...", "  slices="+1+slices+1+"-"+2*slices+1);

// Save registered image stack as a single .tif file
rename("registered");
path_output = out_path+"/registered.tif";
saveAs("Tiff", path_output);

run("Quit");

