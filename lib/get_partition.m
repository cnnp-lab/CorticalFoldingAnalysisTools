function this_partition = get_partition(point, bb, ratio)
%gets the index of a partition in which lies a (x, y, z) point. bb is the
%3D bounding-box of the entire 3D spectrum within which our analysis is
%being made (in this case, a bounding box that envelops the entire brain
%mesh). ratio is cubic root of the total amount of partitions in this space.
rangeX = bb([4]) - bb([1]);
rangeY = bb([5]) - bb([2]);
rangeZ = bb([6]) - bb([3]);

incX = rangeX / ratio;
incY = rangeY / ratio;
incZ = rangeZ / ratio;

nx = point(1) - bb([1]);
ny = point(2) - bb([2]);
nz = point(3) - bb([3]);

currX = ceil(nx / incX);
currY = ceil(ny / incY);
currZ = ceil(nz / incZ);

%fprintf('This vertex is located at %f, %f, %f\n', point(1), point(2), point(3));
%fprintf('This vertexs partition starts at %d, %d, %d\n', currX, currY, currZ);

ratioY = rangeY / rangeX;
ratioZ = rangeZ / rangeX;

this_partition = currX * (ratio)^2 + currY * (ratio) + currZ;
%this_partition = currX * (ceil(ratioY * ratio) * ceil(ratioZ * ratio)) + currY * ceil(ratioZ * ratio) + currZ;
