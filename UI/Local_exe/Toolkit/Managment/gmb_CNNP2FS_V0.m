function gmb_CNNP2FS_V0 (tbl,f_nm,path)

% Detect the number of datasets avaialable
DtSt = unique(tbl.dSet);

% Detect if full hemisphere
hm_dx = tbl.Atlas=="hemisphere";

% Create Root repository for all outout
mkdir(fullfile(path,f_nm));

for i=1:length(DtSt)
    % Values for this dataset
    idx = tbl.dSet == DtSt(i);
    I = find(idx);

    % Generate the repository for this dataset
    mkdir(fullfile(path,f_nm, DtSt(i)));

    % Detect the number of specific Atlas
    Atl = unique(tbl.Atlas(idx & ~hm_dx));

    % Parameters of the study. Divide in this many files
    Param = string(tbl.Properties.VariableNames(8:end));

    if isempty(Atl) %Only Hemisphere dta present
        % From Hemis & Region to ROIs
        roi = replace(tbl.Hemisphere(idx),"left","lh");
        roi = replace(roi,"right","rh");

        % All the regions
        Locs = unique(roi);

        % Subjects for the dataset and Atlas pair
        Sbj = tbl.Subjects(I(roi==Locs(1)));
        
        % A table for each parameter
        for j=1:length(Param)
            % Output table for the parameters
            out_tbl =table();

            % Add the subjects
            out_tbl.Subjects = Sbj;
            
            % For each ROI add a colum to the table
            for k=1:length(Locs)
                out_tbl.(Locs(k)) = tbl.(Param(j))(I(roi==Locs(k)));
            end
            
            % Remove nans and write the table
            if size(out_tbl,1)>1
                T = convertvars(out_tbl, @isnumeric, @gmb_tbl_nanblank_V0);
                writetable(T,fullfile(path,f_nm, DtSt(i),[char(Param(j)),'.csv']))
            else
                T = convertvars(out_tbl([1,1],:), @isnumeric, @gmb_tbl_nanblank_V0);
                writetable(T(1,:),fullfile(path,f_nm, DtSt(i),[char(Param(j)),'.csv']))
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
            Sbj = tbl.Subjects(L(roi==Locs(1)));

            % A table for each parameter
            for j=1:length(Param)
                % Output table for the parameters
                out_tbl =table();

                % Add the subjects
                out_tbl.Subjects = Sbj;

                % For each ROI add a colum to the table
                for k=1:length(Locs)
                    out_tbl.(Locs(k)) = tbl.(Param(j))(L(roi==Locs(k)));
                end

                % Remove nans and write the table
                if size(out_tbl,1)>1
                    T = convertvars(out_tbl, @isnumeric, @gmb_tbl_nanblank_V0);
                    writetable(T,fullfile(path,f_nm, DtSt(i),Atl(l),[char(Param(j)),'.csv']))
                else
                    T = convertvars(out_tbl([1,1],:), @isnumeric, @gmb_tbl_nanblank_V0);
                    writetable(T(1,:),fullfile(path,f_nm, DtSt(i),Atl(l),[char(Param(j)),'.csv']))
                end
            end
        end
    end
end

