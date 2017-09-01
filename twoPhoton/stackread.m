function result = stackread(filename)
  info = imfinfo(filename);
  num_images = numel(info);
  
  for k = 1:num_images
    result(:,:,k) = double(imread(filename, k)) / 65535;
  end
end
