%% select intial points
movie = stackread('sample_movie.tif');

im_maxproj = max(movie,[],3);

points = select_cells(im_maxproj, 2*movie);
%% select more points
points = select_cells(im_maxproj, 2*movie, points);
