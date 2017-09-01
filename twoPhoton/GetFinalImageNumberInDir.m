function finalNumber = GetFinalImageNumberInDir(targetDir)
%%% Extracts the largest number of a named tif file in the targetDir
%%% directory. (This should only be used on folders that contain
%%% a sequence of tif files). 



imgs = dir(targetDir);

imgNums = [];
for i=1:numel(imgs)
    if ~strcmp(imgs(i).name, '.') && ~strcmp(imgs(i).name, '..')
        imgIdx = strsplit(imgs(i).name,'.tif');
        if numel(imgIdx) > 1
            imgIdx = str2num(imgIdx{1});
            imgNums(i) = imgIdx;
        end
    end
end

if numel(imgNums) > 0
    finalNumber = max(imgNums);
else
    finalNumber = 0;
end