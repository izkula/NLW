name="2zz"

rename("original")


selectWindow("original");
run("Z Project...", "projection=[Min Intensity]");
rename("min")


imageCalculator("Subtract create 32-bit stack", "original","min");
rename(name);

selectWindow("original");
close()
selectWindow("min");
close()

