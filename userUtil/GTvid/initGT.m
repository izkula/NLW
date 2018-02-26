function [ vid ] = initGT()
%INITGT Initialize GT2750 camera

devid = 1; %%% Can check that this is the correct one with imaqhwinfo('gige')
vid = videoinput('gige', devid); 
disp('Initialized camera')

end

