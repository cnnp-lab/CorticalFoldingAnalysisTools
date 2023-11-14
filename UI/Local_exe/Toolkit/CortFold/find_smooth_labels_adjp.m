function labels = find_smooth_labels_adjp(toLabel, labeled, label, parts, partsCount, bb, ratio)
sizenow = size(toLabel);
labels = zeros(sizenow(1), 1);
closestID = 1;
%Essentially the same as find_smooth_labels_parts. The only differente is
%that a third for-end loop is added below, and it accounts for the
%different partitions. Instead of having one partition index to look at, we
%have an array of 8, 18 or 27 parts, depending on where the current one is.
for t = 1:size(toLabel)
    [allParts, allPartsLength] = get_adj_parts(toLabel(t,:), bb, ratio);
    closest = Inf;      
    for n = 1:(allPartsLength - 1) %extra for loop
        whereThisPartIs = partsCount(allParts(n),[1 2]);
        for m = 1:whereThisPartIs(1);
            lookinAt = parts(m + whereThisPartIs(2) - 1,:);
            distance = mydist(labeled(lookinAt(2),:), toLabel(t,:));

            if  distance < closest
                closest = distance;
                closestID = lookinAt(2);
            end
        end
        
    labels(t) = label(closestID);
    end
end  