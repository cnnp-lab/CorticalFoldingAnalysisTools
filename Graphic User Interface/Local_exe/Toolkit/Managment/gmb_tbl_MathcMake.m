function idx = gmb_tbl_MathcMake (tbl,vars,vals)

% Initialize the values
L = size(tbl,1);
I = ones(L,1);

% Matchmake the values of the variables of interes
for i=1:length(vars)
    I = I & strcmp(string(tbl.(vars{i})),vals{i});
    if sum(I)==0
        continue;
    end
end

idx = find(I);