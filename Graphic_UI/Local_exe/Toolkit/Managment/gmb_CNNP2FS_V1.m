function gmb_CNNP2FS_V1 (tbl,f_nm,path)

% Detect the number of datasets avaialable
DtSt = unique(tbl.dataset);

% Detect if full hemisphere
hm_dx = tbl.Atlas=="hemisphere";

% Create Root repository for all outout
if ~isfolder(fullfile(path,f_nm))
    mkdir(fullfile(path,f_nm));
end

for i=1:length(DtSt)
    % Values for this dataset
    idx = tbl.dataset == DtSt(i);
    I = find(idx);

    % Generate the repository for this dataset
    if ~isfolder(fullfile(path,f_nm, DtSt(i)))
        mkdir(fullfile(path,f_nm, DtSt(i)));
    end

    % Detect the number of specific Atlas
    Atl = unique(tbl.Atlas(idx & ~hm_dx));

    % Parameters of the study. Divide in this many files
    Param = string(tbl.Properties.VariableNames(9:end));

    if isempty(Atl) %Only Hemisphere data present
        % From Hemis & Region to ROIs
        roi = replace(tbl.Hemisphere(idx),"left","lh");
        roi = replace(roi,"right","rh");

        % All the regions
        Locs = unique(roi);

        % Subjects for the dataset and Atlas pair
        Sbj = tbl.SubjectID(I(roi==Locs(1)));
        if ~(tbl.Session(I(roi==Locs(1))) == "-")
            Sbj = strcat(Sbj,".",tbl.Session(I(roi==Locs(1))));
        end

        % For each Scale generate a separate file
        Scales = unique(tbl.Scale);

        % A table for each parameter and Scale
        for j=1:length(Param)
            for s=1:length(Scales)
                % Output table for the parameters
                out_tbl =table();

                % Add the subjects
                out_tbl.Subjects = Sbj;

                % For each ROI add a colum to the table
                for k=1:length(Locs)
                    k_dx = I(roi==Locs(k) & tbl.Scale==Scales(s));
                    out_tbl.(Locs(k)) = tbl.(Param(j))(k_dx);
                end

                % Remove nans and write the table
                fin_nm = fullfile(path, ...
                    f_nm, ...
                    DtSt(i), ...
                    [char(Param(j)),'_sc=',num2str(Scales(s)),'.csv']);
                if ~isfile(fin_nm)
                    if size(out_tbl,1)>1
                        T = convertvars(out_tbl, @isnumeric, @gmb_tbl_nanblank_V0);
                        writetable(T,fin_nm)
                    else
                        T = convertvars(out_tbl([1,1],:), @isnumeric, @gmb_tbl_nanblank_V0);
                        writetable(T(1,:),fin_nm)
                    end
                else
                    if size(out_tbl,1)>1
                        T = convertvars(out_tbl, @isnumeric, @gmb_tbl_nanblank_V0);
                        writetable(T,fin_nm,'WriteMode','Append','WriteVariableNames',false,'WriteRowNames',true);
                    else
                        T = convertvars(out_tbl([1,1],:), @isnumeric, @gmb_tbl_nanblank_V0);
                        writetable(T(1,:),fin_nm,'WriteMode','Append','WriteVariableNames',false,'WriteRowNames',true);
                    end
                end
            end
        end

    else
        for l=1:length(Atl)
            % Generate the Atlas folder
            mkdir(fullfile(path,f_nm, DtSt(i),Atl(l)))

            % The subjects of interest: th Specific Atlas and the Hemispheres
            ldx = idx & (hm_dx | tbl.Atlas == Atl(l));
            L = find(ldx);

            % Include the lobe on the region of interest
            roi = replace(tbl.Hemisphere(ldx),"left","lh");
            roi = replace(roi,"right","rh");

            % Include the Region on ROI
            reg = tbl.Region(ldx);
            delim = repmat("-",length(reg),1);
            delim(reg == "hemisphere") = "";
            reg(reg == "hemisphere") = "";

            roi = strcat(roi,delim,reg);

            % All the regions
            Locs = unique(roi);

            % Subjects for the dataset and Atlas pair
            Sbj = tbl.SubjectID(I(roi==Locs(1)));
            if ~(tbl.Session(I(roi==Locs(1))) == "-")
                Sbj = strcat(Sbj,".",tbl.Session(I(roi==Locs(1))));
            end

            % For each Scale generate a separate file
            scl = tbl.Scale(ldx);
            Scales = unique(scl);

            % A table for each parameter
            for j=1:length(Param)
                for s=1:length(Scales)
                    % Output table for the parameters
                    out_tbl =table();

                    % Add the subjects
                    out_tbl.Subjects = Sbj;

                    % For each ROI add a colum to the table
                    for k=1:length(Locs)
                        k_dx = I(roi==Locs(k) & scl==Scales(s));
                        out_tbl.(Locs(k)) = tbl.(Param(j))(L(roi==Locs(k)));
                    end

                    % Remove nans and write the table
                    fin_nm = fullfile(path, ...
                        f_nm, ...
                        DtSt(i), ...
                        Atl(l),...
                        [char(Param(j)),'_sc=',num2str(Scales(s)),'.csv']);
                    if ~isfile(fin_nm)
                        if size(out_tbl,1)>1
                            T = convertvars(out_tbl, @isnumeric, @gmb_tbl_nanblank_V0);
                            writetable(T,fin_nm)
                        else
                            T = convertvars(out_tbl([1,1],:), @isnumeric, @gmb_tbl_nanblank_V0);
                            writetable(T(1,:),fin_nm)
                        end
                    else
                        if size(out_tbl,1)>1
                            T = convertvars(out_tbl, @isnumeric, @gmb_tbl_nanblank_V0);
                            writetable(T,fin_nm,'WriteMode','Append','WriteVariableNames',false,'WriteRowNames',true)
                        else
                            T = convertvars(out_tbl([1,1],:), @isnumeric, @gmb_tbl_nanblank_V0);
                            writetable(T(1,:),fin_nm,'WriteMode','Append','WriteVariableNames',false,'WriteRowNames',true)
                        end
                    end
                end
            end
        end
    end
end

