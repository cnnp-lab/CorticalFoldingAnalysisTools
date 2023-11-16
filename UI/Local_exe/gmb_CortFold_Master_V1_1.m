function gmb_CortFold_Master_V1_1(path_0,out_nm_str)

if nargin<1
    path_0 = cd;
end

% Path with necessary Toolkits
addpath(genpath(fullfile(path_0,'ToolKit')))

% Path to the report file for erros and other issues
rep_nm = gmb_NM_Check_V0 ('Report','.txt',fullfile(path_0,'OUTPUT'));
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

% Delimiter of subject and sesion fields
Conf_Sep = ' | ';

% Available Hemisphere mode
Hemi_MD = ["left","right","both","avg","sum"];

% Avaialable Atlas currently
ATL_dx = "LUT";

% Identifier to determine that there is no sesion subfolder on the subject
nSes_id = '_*_';

%% Load the ocnfiguration data
Conf = readtable(conf_path);

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
OW_md = Conf.Mode(3);

% Extraction mode
% Avaialble options: 0->Not Extract 1->FS format 2->CNNP format
EX_md = Conf.Mode(4);

% When both OW_md & EX_md are 0 => Sreening and report mode activated
Scr_Md = ~OW_md & EX_md;
Scr_Md_cnt = 0;

% Request the name of the output to the user if not provided
if nargin<2 && EX_md>0
    prompt = {'Enter output name'};
    dlgtitle = 'Input';
    fieldsize = [1 45];
    definput = {'CFpar'};
    answer = inputdlg(prompt,dlgtitle,fieldsize,definput);
    out_nm_str = answer{1};
end

%% Run the paramters estimator thrugh all the included datasets and files

% Detect the field with proper info
I = zeros(size(Conf(:,1)));
for i=1:(length(Conf_Vars)-1)
    I = I | ~cellfun(@isempty,Conf.(Conf_Vars{i}));
end

% Table storing all the cortical folding parameters desired
CFextr_tbl = [];

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

% Go though all datasets included
for i=find(I)
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

    % Extract all subject & sesions pair
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


        % Instances of this subect (Number of sesions)
        jdx = find(strcmp(Subj,IDs{j}));
        for k=1:length(jdx)

            % Paramtetres extraction mode: Lobe | Hemisphere
            for p=1:2 % Current 2 parameter modes, Scales not implemented yet
                if Conf.Mode(p) ~= 0
                    %% Determine the number of instances for each case

                    if ~isempty(SbOr_tbl)
                        switch p
                            case 1 % For Lobes
                                aux = load(fullfile(path_0,'ToolKit','CortFold','Atlas',[char(ATL_md(p)),'.mat']));
                                P_n = length(unique(aux.Map(:,2)));
                            case 2 % For Hemispheres
                                P_n = 1;
                        end

                        % Fields of the matrix to compare
                        vars = SbOr_tbl.Properties.VariableNames(3:6);

                        % Values to compare to
                        vals = ["",ATL_md(p)];

                        % Add sesion
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
                            p_dx = sort([p_dx_2;p_dx_1],'ascend');
                            flag = length(p_dx) == P_n*2;
                        else
                            vals(2) = Hemi_MD(Conf.Mode(end));
                            p_dx = gmb_tbl_MathcMake (SbOr_tbl,vars,vals);
                            flag = length(p_dx) == P_n;
                        end
                    else
                        p_dx = [];
                        flag = 0;
                    end

                    
                    %% Data extraction
                    switch OW_md
                        case 0
                            if flag
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
                            if flag
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
                                    Ses{jdx(k)},...% Name/Code of the sesion
                                    report,...% Text for the final report
                                    p,...% Code indicating Lobe(1) Hemisphere(2)
                                    Hemi_MD(Conf.Mode(end)),...% Treatmetn of the hemisphere data
                                    ATL_md(p));% For Lobe estimation, the atlas to be used
                            end

                        case 2
                            % Estimate parameters
                            [tbl,report] = gmb_CF_2table_V1(Sub_path,...% Path to subject data
                                Ses{jdx(k)},...% Name/Code of the sesion
                                report,...% Text for the final report
                                p,...% Code indicating Lobe(1) Hemisphere(2)
                                Hemi_MD(Conf.Mode(end)),...% Treatmetn of the hemisphere data
                                ATL_md(p));% For Lobe estimation, the atlas to be used
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

                        % Add the Sesion [FOR NOW JUST EMPTY]
                        if ~strcmp(Ses{jdx(k)},'_*_')
                            aux = string(Ses{jdx(k)});
                        else
                            aux = "-";
                        end

                        tbl.Sesion = repmat(aux,n_row,1);
                        tbl = movevars(tbl,"Sesion",'Before',1);

                        % Add the Subject IDs
                        tbl.Subjects = repmat(string(IDs{j}),n_row,1);
                        tbl = movevars(tbl,"Subjects",'Before',1);

                        % Add the Dataset
                        tbl.dSet = repmat(string(DtSt_nm),n_row,1);
                        tbl = movevars(tbl,"dSet",'Before',1);

                        % Add this data to the subject file
                        SbSs_tbl = [SbSs_tbl;tbl];
                    end
                end
            end
        end

        % Store subject parameteres on its folder
        % NaN to blanck
        if OW_md>0
            if size(SbSs_tbl,1)>1
                T = convertvars(SbSs_tbl, @isnumeric, @gmb_tbl_nanblank_V0);
                writetable(T,fullfile(Sub_path,[IDs{j},'_CFpar.csv']))
            elseif size(SbSs_tbl,1)==1
                T = convertvars(SbSs_tbl([1,1],:), @isnumeric, @gmb_tbl_nanblank_V0);
                writetable(T(1,:),fullfile(Sub_path,[IDs{j},'_CFpar.csv']))
            end
        end

        % Store the subject estimated data on the extraction table
        CFextr_tbl =[CFextr_tbl;SbSs_tbl];

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
        out_nm_str = [out_nm_str,'_FS'];

        % Check if there are other output files to not overwrite
        out_nm_str = gmb_NM_Check_V0 (out_nm_str,'.csv',fullfile(path_0,'OUTPUT'));

        % Convert the current format
        gmb_CNNP2FS_V0 (CFextr_tbl,out_nm_str,fullfile(path_0,'OUTPUT'));

        % Store the report
        report = [report; " ";" ";"Parameters stored on folder:";string(fullfile(cd,'OUTPUT',out_nm_str))];

    case 2 % Store data with CNNP propietary file format
        out_nm_str = [out_nm_str,'_CNNP'];

        % Check if there are other output files to not overwrite
        out_nm_str = gmb_NM_Check_V0 (out_nm_str,'.csv',fullfile(path_0,'OUTPUT'));

        % Store the parameters
        % NaN to blanck
        if size(CFextr_tbl,1)>1
            T = convertvars(CFextr_tbl, @isnumeric, @gmb_tbl_nanblank_V0);
            writetable(T,fullfile(path_0,'OUTPUT',[out_nm_str,'.csv']))
        else
            T = convertvars(CFextr_tbl([1,1],:), @isnumeric, @gmb_tbl_nanblank_V0);
            writetable(T(1,:),fullfile(path_0,'OUTPUT',[out_nm_str,'.csv']))
        end

        report = [report; " ";" ";"Parameters stored on file:";string(fullfile(cd,'OUTPUT',[out_nm_str,'.csv']))];
end

% Store the report
if ~Scr_Md
    writematrix(report,rep_path)
end

close(H);

msgbox("Analysis Done","Done!!","help");


