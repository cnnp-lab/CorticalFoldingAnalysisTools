function Scr_Rep = gmb_CortFold_Master_V1_1(path_0,out_nm_str)

T0 = datetime("now");

if nargin<1
    path_0 = cd;
end

% Path with necessary Toolkits
addpath(genpath(fullfile(path_0,'Toolkit')))

% Path to the report file for erros and other issues
rep_nm = gmb_NM_Check_V0 ('Report','.txt',fullfile(path_0,'OUTPUT'));
if ~isfolder(fullfile(path_0,'OUTPUT'))
    mkdir(fullfile(path_0,'OUTPUT'));
end
rep_path = fullfile(path_0,'OUTPUT',[rep_nm,'.txt']);
report = [];

%% Configuration file reading and setting up the data reading
% Path to the configuration file generated with the UI
conf_path = fullfile(path_0,'Config_Fold.csv');

% Check presence of file
if ~isfile (conf_path)
    report = [report;"AVORTED => Configuration file not found"];
    writematrix(report,rep_path)
    return;
end

% Structure of the expected configuration files
Conf_Vars = {'Root','FS','Subj','Ses','Mode'};

% Delimiter of subject and Session fields
Conf_Sep = ' | ';

% Available Hemisphere mode
Hemi_MD = ["left","right","both","avg","sum"];

% Avaialable Atlas currently
ATL_dx = "FSDK";

% Identifier to determine that there is no Session subfolder on the subject
nSes_id = '_*_';

%% Load the configuration data
Conf = readtable(conf_path,"Delimiter","comma");

% Check the quality of the configuratoin file
flag = 0;
if length(Conf.Properties.VariableNames)~=length(Conf_Vars) % Number of variables
    flag = 1;
else
    % Check variables naming matching dpected
    flbg = zeros(size(Conf_Vars));
    for i=1:length(Conf_Vars)
        flbg = flbg | strcmp(Conf.Properties.VariableNames,Conf_Vars{i});
    end
    if sum(flbg) ~= length(flbg)
        flag = 1;
    end
end

if flag
    % Not valid configuration file
    report = [report;"AVORTED => Configuration file not valid"];
    writematrix(report,rep_path)
    return;
end

% Determine which of the atlases are used for Lobe estimation
ATL_md = ["","hemisphere"]; % Hemisphere estimations specific naming
if Conf.Mode(1)>0
    ATL_md(1) = ATL_dx(Conf.Mode(1));
end

% Ovewrite mode for the geneerated dataset
% Avaialble options: 0->Nothing 1->Append 2->Overwrite
OW_md = Conf.Mode(end-2);

% Extraction mode
% Avaialble options: 0->Not Extract 1->FS format 2->CNNP format
EX_md = Conf.Mode(end-1);

% Determine the scales to expect based on the scaling mode
targetScale= 0.225:0.025:.375;
nIter=[5 5 4 4 4 4 4];
scale_FullSet = [];
for i=1:length(targetScale)
    scale_FullSet = [scale_FullSet,targetScale(i).*(2.^(0:nIter(i)))];
end

% Request the name of the output to the user if not provided
if nargin<2 && EX_md>0
    prompt = {'Enter output name'};
    dlgtitle = 'Input';
    fieldsize = [1 45];
    definput = {'CFpar'};
    answer = inputdlg(prompt,dlgtitle,fieldsize,definput);
    out_nm_str{1} = answer{1};
    out_nm_str{2} = fullfile(path_0,'OUTPUT');
end

% Definr the output file name
switch EX_md
    case 1
        rmdir(fullfile(out_nm_str{2},out_nm_str{1}),'s')
        out_nm_str{1} = [out_nm_str{1},'_FS'];

        % Check if there are other output files to not overwrite
        out_nm_str{1} = gmb_NM_Check_V0 (out_nm_str{1},'.csv',out_nm_str{2});
    case 2
        out_nm_str{1} = [out_nm_str{1},'_CNNP'];

        % Check if there are other output files to not overwrite
        out_nm_str{1} = gmb_NM_Check_V0 (out_nm_str{1},'.csv',out_nm_str{2});
        
        out_path = fullfile(out_nm_str{2},[out_nm_str{1},'.csv']);
end

%% Run the paramters estimator thrugh all the included datasets and files

% Detect the field with proper info
I = zeros(size(Conf(:,1)));
for i=1:(length(Conf_Vars)-1)
    I = I | ~cellfun(@isempty,Conf.(Conf_Vars{i}));
end
idx = find(I)';

% When both OW_md & EX_md are 0 => Sreening and report mode activated
Scr_Md = ~OW_md & ~EX_md;
Scr_Rep = table(repmat("",length(idx),1),zeros(length(idx),1),...
    zeros(length(idx),1),'VariableNames',{'dSet','N','cnt'});
cnt = 0;

% Show Progress
H = waitbar(0,'','Name','PROGRESS');
H.Children.Title.Interpreter = 'none';

% Estimate the total amount of time required
aux = Conf.Subj(I);% Subject Set
a = cellfun(@split, aux, repmat({Conf_Sep},size(aux)), 'UniformOutput', false);

% Number of subjects
s = cellfun(@length,cellfun(@unique,a, 'UniformOutput', false));
t = [];
tot_T = sum(s*5);

% Flag indicating that No data has been stored yet
store_flag = 0;

% Go though all datasets included
for i=idx
    % Naming the dataset
    [~,DtSt_nm,~] = fileparts(Conf.Root{i});

    % Update waitbar
    H.Name = string([DtSt_nm,'(',num2str(i),'/',num2str(sum(I)),') ',num2str(round(tot_T/60)),' mins.']);

    %% Check points

    % Check if dataset available
    FS_path = fullfile(Conf.Root{i},Conf.FS{i});
    if ~isfolder(FS_path)
        report = [report; string(['WARNING => Dataset ',DtSt_nm,' not found';' '])];
        continue;
    end

    % Extract all subject & Sessions pair
    Subj= split(Conf.Subj{i}, Conf_Sep);
    Ses = split(Conf.Ses{i} , Conf_Sep);

    % Check inconsistent structure
    if length(Subj)~=length(Ses)
        report = [report; string(['WARNING => Dataset ',DtSt_nm,' skiped due to inconsistent file location';' '])];
        continue;
    end

    %% Analysis of the dataset

    % Indicate that the analysis of the dataset has started
    report = [report; string(['Dataset ',DtSt_nm,' initiated'])];

    % All subjects available on the dataset
    IDs = unique(Subj);

    % Initialize the Dataset field of the Screening Output
    if Scr_Md
        cnt = cnt+1;
        Scr_Rep.dSet(cnt) = string(DtSt_nm);
        Scr_Rep.cnt (cnt) = 0;
        Scr_Rep.N (cnt)   = length(IDs);
    end

    for j=1:length(IDs)
        % Meassure the time of a single itteration
        tic;

        % Update Waitbar
        H.Name = string([DtSt_nm,'(',num2str(i),'/',num2str(sum(I)),') ',num2str(round(tot_T,2)),' mins.']);
        waitbar(j/length(IDs),H,['Progress: ',num2str(j),'/',num2str(length(IDs))]);

        % Repository for subject parameters to extract and estimate
        SbSs_tbl = [];

        % Directory with the files per subject
        Sub_path = fullfile(FS_path,IDs{j});

        % Location of for the Repository csv
        SubParam_file = fullfile(Sub_path,[IDs{j},'_CFpar.csv']);

        % Read the existing data if needed & avaialable
        SbOr_tbl = [];
        if ~isfile(SubParam_file)
            if OW_md == 1
                OW_md = 2;
            end
        elseif OW_md ~= 2
            % Load previouslly existing data
            opts = detectImportOptions(SubParam_file);
            opts.VariableTypes = repmat({'double'},1,length(opts.VariableTypes));
            opts.VariableTypes(1:7) = repmat({'string'},1,7);
            SbOr_tbl = readtable(SubParam_file,opts);
        end

        % Indicate that the analysis of the dataset has started
        report = [report; string(['   Subject ',IDs{j},' initiated'])];


        % Instances of this subect (Number of Sessions)
        jdx = find(strcmp(Subj,IDs{j}));

        % Indicator if data extracted for he desired configuration
        CnfMd_flg = 0;
        for k=1:length(jdx)
            % Paramtetres extraction mode: Lobe | Hemisphere
            for p=1:2 % Current 2 parameter modes, Scales not implemented yet
                if Conf.Mode(p) ~= 0
                    %% Determine the number of instances for each case
                    switch p
                        case 1 % For Lobes
                            scale_N = 0; % For LOBES multiscale not implemented
                            
                        case 2 % For Hemispheres
                            scl_cnt = dec2bin(Conf.Mode(3),3);
                            scale_N = [];
                            
                            if scl_cnt(3)=='1'
                                scale_N = [scale_N,0];
                            end
                            
                            if scl_cnt(2)=='1'
                                scale_N = [scale_N,scale_FullSet(1:22)];
                            end

                            if scl_cnt(1)=='1'
                                scale_N = [scale_N,scale_FullSet(23:end)];
                            end
                            
                    end
                    if ~isempty(SbOr_tbl)
                        switch p
                            case 1 % For Lobes
                                aux = load(fullfile(path_0,'Toolkit','CortFold','Atlas',[char(ATL_md(p)),'.mat']));
                                P_n = length(unique(aux.Map.Lobe_cd));

                            case 2 % For Hemispheres
                                P_n = 1;
                        end

                        % Fields of the matrix to compare
                        vars = SbOr_tbl.Properties.VariableNames([3:6,8]);

                        % Values to compare to
                        vals = ["",ATL_md(p)];

                        % Add Session
                        if ~strcmp(Ses{jdx(k)},'_*_')
                            vals = [string(Ses{jdx(k)}),vals];
                        else
                            vals = ["-",vals];
                        end

                        % Add Finess
                        switch p
                            case 1 % Lobe
                                vals = [vals,"lobe"];
                            case 2 % Hemisph
                                vals = [vals,ATL_md(p)];
                        end

                        p_dx = [];
                        flag = [];
                        vals = [vals,"0"];
                        for sc = 1:length(scale_N)
                            vals(end) = string(scale_N(sc));
                            % Special case for hemisphere parameter
                            % "Both" generates 2x the data (right and left)
                            if strcmp(Hemi_MD(Conf.Mode(end)),"both")
                                % left hemisphere
                                vals(2) = "left";
                                p_dx_1 = gmb_tbl_MathcMake (SbOr_tbl,vars,vals);

                                % left hemisphere
                                vals(2) = "right";
                                p_dx_2 = gmb_tbl_MathcMake (SbOr_tbl,vars,vals);

                                % Both hemispheres
                                aux = sort([p_dx_2;p_dx_1],'ascend');
                                p_dx = [p_dx,aux];
                                flag(sc) = length(aux) == P_n*2;
                            else
                                vals(2) = Hemi_MD(Conf.Mode(end));
                                aux = gmb_tbl_MathcMake (SbOr_tbl,vars,vals);
                                p_dx = [p_dx,aux];
                                flag(sc) = length(p_dx) == P_n;
                            end
                        end
                    else
                        p_dx = [];
                        flag = 0;
                    end



                    %% Include data was not found for the configuration
                    % Count the configurations that have not been estimated
                    CnfMd_flg = CnfMd_flg+~(length(flag)==sum(flag));

                    %% Data extraction
                    switch OW_md
                        case 0
                            if length(flag)==sum(flag)
                                % Extract them from existing file
                                tbl = SbOr_tbl(p_dx,[4,7:end]);

                                aux = '      EXTRACTED => ';
                                if ~strcmp(Ses{jdx(k)},'_*_')
                                    aux = [aux,Ses{jdx(k)},'/'];
                                end
                                aux = [aux,char(Hemi_MD(Conf.Mode(end))),' data obtained'];
                                report = [report; string(aux)];
                            else
                                aux = '      WARNIN => ';
                                if ~strcmp(Ses{jdx(k)},'_*_')
                                    aux = [aux,Ses{jdx(k)},'/'];
                                end
                                aux = [aux,char(Hemi_MD(Conf.Mode(end))),' data not found'];
                                report = [report; string(aux)];
                                tbl = [];
                            end
                        case 1 % Check data presence and then decide if procede

                            % Estimate parameters if needed
                            if length(flag)==sum(flag)
                                % Extract them from existing file
                                tbl = SbOr_tbl(p_dx,[4,7:end]);

                                aux = '      EXTRACTED => ';
                                if ~strcmp(Ses{jdx(k)},'_*_')
                                    aux = [aux,Ses{jdx(k)},'/'];
                                end
                                aux = [aux,char(Hemi_MD(Conf.Mode(end))),' data obtained'];
                                report = [report; string(aux)];

                            else
                                % Estimate parameters
                                [tbl,report] = gmb_CF_2table_V1(Sub_path,...% Path to subject data
                                    Ses{jdx(k)},...% Name/Code of the Session
                                    report,...% Text for the final report
                                    p,...% Code indicating Lobe(1) Hemisphere(2)
                                    Hemi_MD(Conf.Mode(end)),...% Treatmetn of the hemisphere data
                                    ATL_md(p),...% For Lobe estimation, the atlas to be used
                                    Conf.Mode(3));% Scaling Mode: None, High &| resol.
                            end

                        case 2
                            % Estimate parameters
                            [tbl,report] = gmb_CF_2table_V1(Sub_path,...% Path to subject data
                                Ses{jdx(k)},...% Name/Code of the Session
                                report,...% Text for the final report
                                p,...% Code indicating Lobe(1) Hemisphere(2)
                                Hemi_MD(Conf.Mode(end)),...% Treatmetn of the hemisphere data
                                ATL_md(p),...% For Lobe estimation, the atlas to be used
                                Conf.Mode(3));% Scaling Mode: None, High &| resol.
                    end

                    % Complete the table and store it
                    if ~isempty(tbl)
                        n_row = size(tbl,1);

                        % Add the atlas identification
                        switch p
                            case 1 % Lobe
                                tbl.Finess = repmat("lobe",n_row,1);
                            case 2 % Hemisph
                                tbl.Finess = repmat(ATL_md(p),n_row,1);
                        end
                        tbl = movevars(tbl,"Finess",'After',1);

                        % Add the atlas identification
                        tbl.Atlas = repmat(ATL_md(p),n_row,1);
                        tbl = movevars(tbl,"Atlas",'After',1);

                        % Add the Session [FOR NOW JUST EMPTY]
                        if ~strcmp(Ses{jdx(k)},'_*_')
                            aux = string(Ses{jdx(k)});
                        else
                            aux = "-";
                        end

                        tbl.Session = repmat(aux,n_row,1);
                        tbl = movevars(tbl,"Session",'Before',1);

                        % Add the Subject IDs
                        tbl.SubjectID = repmat(string(IDs{j}),n_row,1);
                        tbl = movevars(tbl,"SubjectID",'Before',1);

                        % Add the Dataset
                        tbl.dataset = repmat(string(DtSt_nm),n_row,1);
                        tbl = movevars(tbl,"dataset",'Before',1);

                        % Add this data to the subject file
                        SbSs_tbl = [SbSs_tbl;tbl];
                    end
                end
            end
        end

        % Add results to the Sreeninng report table
        % Increase by 1 the number of complete subjects
        if Scr_Md && CnfMd_flg==0
            Scr_Rep.cnt (cnt) = Scr_Rep.cnt (cnt)+1;
        end

        % Store subject parameteres on its folder
        % NaN to blanck
        switch  OW_md
            case 1
                TAB = SbOr_tbl;
                vars = SbSs_tbl.Properties.VariableNames(1:8);
                for k =1:size(SbSs_tbl,1)
                    vals = SbSs_tbl(k,1:8);
                    aux = gmb_tbl_MathcMake (TAB,vars,vals.Variables);
                    if any(aux)
                        TAB(aux(1),:) = SbSs_tbl(k,:);
                    else
                        TAB = [TAB;SbSs_tbl(k,:)];
                    end

                end
                
            case 2
                TAB = SbSs_tbl;
        end

        if OW_md>0
            if size(TAB,1)>1
                T = convertvars(TAB, @isnumeric, @gmb_tbl_nanblank_V0);
                writetable(T,fullfile(Sub_path,[IDs{j},'_CFpar.csv']))
            elseif size(TAB,1)==1
                T = convertvars(TAB([1,1],:), @isnumeric, @gmb_tbl_nanblank_V0);
                writetable(T(1,:),fullfile(Sub_path,[IDs{j},'_CFpar.csv']))
            end
        end

        % Store the subject estimated data on the extraction table
        if ~isempty(SbSs_tbl)
            switch EX_md
                case 1 % Store data with FreeSurfer file format
                    % Convert the current format
                    gmb_CNNP2FS_V1 (SbSs_tbl,out_nm_str{1},out_nm_str{2});

                case 2 % Store data with CNNP propietary file format

                    if ~store_flag % First individual, generate the report file
                        % Store the parameters
                        % NaN to blanck
                        if size(SbSs_tbl,1)>1
                            T = convertvars(SbSs_tbl, @isnumeric, @gmb_tbl_nanblank_V0);
                            writetable(T,out_path)
                        else
                            T = convertvars(SbSs_tbl([1,1],:), @isnumeric, @gmb_tbl_nanblank_V0);
                            writetable(T(1,:),out_path)
                        end

                        % Set the flag to append
                        store_flag = 1;
                    else
                        % Store the parameters
                        % NaN to blanck
                        if size(SbSs_tbl,1)>1
                            T = convertvars(SbSs_tbl, @isnumeric, @gmb_tbl_nanblank_V0);
                            writetable(T,out_path,'WriteMode','Append','WriteVariableNames',false,'WriteRowNames',true);
                        else
                            T = convertvars(SbSs_tbl([1,1],:), @isnumeric, @gmb_tbl_nanblank_V0);
                            writetable(T(1,:),out_path,'WriteMode','Append','WriteVariableNames',false,'WriteRowNames',true);
                        end
                    end
            end
        end

        % Estimate the time it took to run a single subject
        t(end+1) = toc;
        s(i) = s(i)-1;
        tot_T = sum(s*mean(t)/60);
    end
end

%% Write the final report and save the generated values
switch EX_md
    case 0
        if OW_md~=0
            report = [report; " ";" ";"Parameters estimated"];
        end

    case 1 % Store data with FreeSurfer file format
        report = [report; " ";" ";"Parameters stored on folder:";string(fullfile(out_nm_str{2},out_nm_str{1}))];

    case 2 % Store data with CNNP propietary file format
        report = [report; " ";" ";"Parameters stored on file:";string(fullfile(out_nm_str{2},[out_nm_str{1},'.csv']))];
end

T1 = datetime("now");
report = [report; " ";" ";string(['Time that the analysis took: ',char(string(T1-T0))])];

% Store the report
if ~Scr_Md
    writematrix(report,rep_path)
    msgbox("Analysis Done","Done!!","help");
end

close(H);






