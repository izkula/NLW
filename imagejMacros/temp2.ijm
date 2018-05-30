//Run MOCO
run("moco ", "value=43 downsample_value=1 template=template stack=z0_processed log=None plot=[No plot]");

//Remane Registered Stack
rename("z0_processed_reg.tif")
