function [allParts, counter] = get_adj_parts(point, bb, ratio)

incX = (bb(4) - bb(1)) / ratio;
incY = (bb(5) - bb(2)) / ratio;
incZ = (bb(6) - bb(3)) / ratio;

counter = 1;

allParts = zeros(27, 1);

xValue = point(1) - (2 * incX);
for x = 1:3
    xValue = xValue + incX;
    yValue = point(2) - (2 * incY);
    for y = 1:3
        yValue = yValue + incY;
        zValue = point(3) - (2 * incZ);
        for z = 1:3
            zValue = zValue + incZ;
            if xValue > bb(1) && xValue < bb(4) && yValue > bb(2) && yValue < bb(5) && zValue > bb(3) && zValue < bb(6)
                allParts(counter) = get_partition([xValue, yValue, zValue], bb, ratio);
                counter = counter + 1;
            end
        end
    end
end

%shiftedLeftX = point(t,1) - incX;
%shiftedRightX = point(t,1) + incX;
%shiftedDownY = point(t,2) - incY;
%shiftedUpY = point(t,2) + incY;
%shiftedInZ = point(t,3) - incZ;
%shiftedOutZ = point(t,3) + incZ;


%pointInQuestion = point;
%if shiftedLeftX > bb(1)
%    pointInQuestion(1) = shiftedLeftX;
%    allParts(counter) = get_partition(pointInQuestion, bb, ratio);
%    counter = counter + 1;
%    if shiftedDownY > bb(2)
%        pointInQuestion(2) = shiftedDownY;
%        allParts(counter) = get_partition(pointInQuestion);
%        counter = counter + 1;
%    end
%    
%    if shiftedUpY < bb(5)
%        pointInQuestion(2) = shiftedUpY;
%        allParts(counter) = get_partition(pointInQuestion);
%        counter = counter + 1;
%    end

% if shifted(
