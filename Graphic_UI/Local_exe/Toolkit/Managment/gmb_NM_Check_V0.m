function NM = gmb_NM_Check_V0 (NM,ext,path)

flag = 0;


if isfile(fullfile(path,[NM,ext])) || isfolder(fullfile(path,NM))
    cnt = 1;
    flag = 1;
end

while flag
    if isfile(fullfile(path,[NM,'_',num2str(cnt),ext])) ||...
           isfolder(fullfile(path,[NM,'_',num2str(cnt)])) 
        cnt = cnt+1;
    else
        NM = [NM,'_',num2str(cnt)];
        flag = 0;
    end
end