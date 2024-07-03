function [parts, partsCount] = give_parts_to_vertices(vertices, bb, ratio)

%get the partition of each vertex, storing it on the parts vertex
for k =  1:size(vertices)
    k_current = vertices(k,:);
    currPartInd = get_partition(k_current, bb, ratio);
    parts(k,:) = [currPartInd, k];
end
%sort the rows of the partition array so that we can be sure of what we'll
%do next
parts = sortrows(parts);
ratioY = (bb(5) - bb(2))/ (bb(4) - bb(1));
ratioZ = (bb(6) - bb(3))/ (bb(4) - bb(1));
parts_prev = parts(1);      %look at first partition
%for y = 1:((ratioY * (ratio + 1) * ratioZ * (1 + ratio) * (1 + ratio)))
for y = 1:((ratio + 1) ^ 3)
    partsCount(y,:) = [0, 1];
end

%each element of partsCount signifies [wh,hm], where
% --- wh = which index, on the sorted parts array, does the vertices tagged
% -------- by this partition start appearing
% --- hm = how many vertices belong to this partition. we need this to know
% -------- how many steps we take on the parts array during our labeling 

for y = 1:size(parts)
    y_current = parts(y);
    if y_current == parts_prev %if same, increment count of elements in partition
        partsCount(y_current,[1]) = partsCount(y_current,[1]) + 1;
    else
        whereStarts = partsCount(parts_prev,[1]) + partsCount(parts_prev, [2]);
        partsCount(y_current,:) = [1, whereStarts];
        parts_prev = y_current;
    end
end