arg = getArgument()
print("Running batch analysis with arguments:")
print(arg)

run("Image Sequence...", "open=/media/Data/data/BpodImageData/vglut1m15/Image2PShapeOlfGNG/Dec10_2015/Session1/Trial00001/vglut1m15_Image2PShapeOlfGNG_Dec10_2015_Session1-001/vglut1m15_Image2PShapeOlfGNG_Dec10_2015_Session1-001_Cycle00001_Ch2_000001.ome.tif sort");

run("PureDenoise ", "parameters='3 1' estimation='Auto Global' ");
rename("Denoised");

//saveAs("Tiff", dir)

//while (nImages>0) { 
//	selectImage(nImages); 
//        close(); 
//} 