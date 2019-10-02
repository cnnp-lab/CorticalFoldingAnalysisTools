function labels = find_smooth_labels_brute(toLabel, labeled, label)
%so we can preallocate labels
labels = zeros(length(toLabel), 1);
for t = 1:size(toLabel)
    if mod(t, 10) == 0
        time = toc;
        fprintf('Time for last 10 vertices = %f\n', time);
        tic;
    end
    %tentative distance to the closest vertex is infinity, of course.
    closest = Inf;       
    %for each vertex,
    for m = 1:size(labeled)
        %find the distance
        distance = dist(labeled(m,:), toLabel(t,:));
        %decide new tentative distance and vertexID of the distance
        if  distance < closest
            closest = distance;
            closestID = m;
        end
    end
    %after we've looked at every vertex, set the label of
    %this one to be the same as the closest one we've found.
    labels(t) = label(closestID);
    fprintf('%d. distance weve found = %f\n', t, closest);
end  