function bb = get_bounding_box(vertices)

bb = zeros(6, 1); %preallocate for efficiency
greatestX = -Inf;
smallestX = Inf;
greatestY = -Inf;
smallestY = Inf;
greatestZ = -Inf;
smallestZ = Inf; % set initial values for greater bounding box determination

for k = 1:size(vertices)
    k_current = vertices(k,:);
    if k_current(1) > greatestX
        greatestX = k_current(1);
    end
    if k_current(1) < smallestX
        smallestX = k_current(1);
    end
    if k_current(2) > greatestY
        greatestY = k_current(2);
    end
    if k_current(2) < smallestY
        smallestY = k_current(2);
    end
    if k_current(3) > greatestZ
        greatestZ = k_current(3);
    end
    if k_current(3) < smallestZ
        smallestZ = k_current(3);
    end
end
bb(1) = smallestX - 0.1; % a little 
bb(2) = smallestY - 0.1; % 0.1 margin
bb(3) = smallestZ - 0.1; % just to 
bb(4) = greatestX + 0.1; % account for
bb(5) = greatestY + 0.1; % floating
bb(6) = greatestZ + 0.1; % point error