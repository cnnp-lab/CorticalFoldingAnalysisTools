function boundingboxes = find_bb(fullSurf, label)

lowerleftin = [-Inf, -Inf, -Inf];
upperrightout = [Inf, Inf, Inf];

%label the labels
counter = 1;
for j = 1:size(label)
    if any(label(j) ~= labelOfLabels(counter))
        counter = counter + 1;
        labelOfLabels(counter) = label(j);
    end
end
        

for i = 1:size(fullSurf)
    
end



