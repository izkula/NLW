function [subjects] = get_unique_subjects(f)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

subjects = {};
j = 1;
for i = 1:numel(f)
    C = strsplit(f(i).name,{'_','.'});
    try
    if strmatch(C{1}(1),'m')
    subjects{j} = lower(C{1});
    j = j+1;
    end
    end
end
subjects(strcmp('',subjects)) = [];
subjects = unique(subjects); %unique subjects in the whole project
       

end

