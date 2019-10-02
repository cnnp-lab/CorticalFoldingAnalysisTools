function labels = find_smooth_labels_parts(toLabel, labeled, label, parts, partsCount, bb, ratio)


%sizenow is just to get around matlab's obnoxiousness
sizenow = size(toLabel);
%so we can preallocate labels
labels = zeros(sizenow(1), 1);
%set the closestID to 1, so that in the unlikely event that NO vertex lies
%in any partition, they'll all think the closest vertex is the first one.
closestID = 1;
for t = 1:size(toLabel)
    %get_partition function gives partition index.  it's self explanatory
    partOfThisVertex = get_partition(toLabel(t,:), bb, ratio);
    %partsCount is an auxiliary array that quickens the search for a part
    whereThisPartIs = partsCount(partOfThisVertex,[1 2]);
    %tentative distance to the closest vertex is infinity, of course.
    closest = Inf;       
    %now for each vertex of the full mesh in the partition of this one,
   
    for m = 1:whereThisPartIs(1);
        %use data from partsCount (whereThisPartIs) to know which vertices
        %to compare this one to
        lookinAt = parts(m + whereThisPartIs(2) - 1,:);
        %find the distance
        distance = dist(labeled(lookinAt(2),:), toLabel(t,:));
        %decide new tentative distance and vertexID of the distance
        if  distance < closest
            closest = distance;
            closestID = lookinAt(2);
        end
    end
    
    %after we've looked at every vertex in the partition, set the label of
    %this one to be the same as the closest one we've found.
    labels(t) = label(closestID);
end  