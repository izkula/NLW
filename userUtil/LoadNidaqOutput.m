function [t, ch] = LoadNidaqOutput(outputFname_daq, nChannels)

fid = fopen(outputFname_daq,'r');
[data,count] = fread(fid,[nChannels+1,inf],'double');
fclose(fid);

t = data(1,1:end);
ch = data(2:end,1:end);