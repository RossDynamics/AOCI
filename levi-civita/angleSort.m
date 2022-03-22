function [sortedPts, indices, sortedthetas] = angleSort(points, mu)
%ANGLESORT Converts a set of Cartesian, barycentric points to secondary-
%centered polar coordinates and then sorts them based on their theta value.
%The second output argument is the sorting indices, and the third is the
%sorted theta values.

polarpts = cart2m2polar(points, mu);
[sortedthetas, indices] = sort(polarpts(2,:));
sortedPts = points(:,indices);

end

