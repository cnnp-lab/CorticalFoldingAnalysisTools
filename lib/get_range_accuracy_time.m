function timeVsRatios = get_range_accuracy_time(fm, fml, sm, gt, method, ratio_i, ratio_f)

bb = get_bounding_box([fm;sm]);
%partitions the 3D space by tagging vertices, and gets info on where each
%part lies on the parts array, storing that data in partsCount.

timeVsRatios = zeros(27,3);

for i = ratio_i:ratio_f
    tic;
    [parts, partsCount] = give_parts_to_vertices(fm, bb, i);
    labels_parts = method(sm, fm, fml, parts, partsCount, bb, i);
    time = toc
    accuracy = compute_accuracy(gt, labels_parts);
    %jarilson = [i, accuracy, time]
    timeVsRatios(i - ratio_i + 1,:) = [i, time, accuracy];
end

%timeVsRatios = jarilson;