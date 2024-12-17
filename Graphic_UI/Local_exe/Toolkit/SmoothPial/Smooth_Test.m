function [Smth_cnt,path] = Smooth_Test(path_0)

% Path with necessary Toolkits
addpath(genpath(fullfile(path_0,'ToolKit')))

% Path to the configuration file generated with the UI
conf_path = fullfile(path_0,'Config_Fold.csv');

% Load the configuration
Conf = readtable(conf_path,"Delimiter","comma");

% Structure of the expected configuration files
Conf_Vars = {'Root','FS','Subj','Ses','Mode'};

% Delimiter of subject and Session fields
Conf_Sep = ' | ';

% Detect the field with proper info
I = zeros(size(Conf(:,1)));
for i=1:(length(Conf_Vars)-1)
    I = I | ~cellfun(@isempty,Conf.(Conf_Vars{i}));
end
idx = find(I)';

% Counting the number of appropriate files
Smth_cnt =zeros(sum(I),1);
cnt = 0;

% Paths to the subjects/sesions without the smooth surface file
path = cell(0);
p_nt = 0;

switch Conf.Mode(end) 
    case 1
        Hemi_MD = "lh";
    case 2
        Hemi_MD = "rh";
    otherwise
        Hemi_MD = ["lh","rh"];
end

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

% Pathc the memory leaking to the RAM
clear functions;

for i=idx
    cnt = cnt+1;

    % Naming the dataset
    [~,DtSt_nm,~] = fileparts(Conf.Root{i});

    % Update waitbar
    H.Name = string([DtSt_nm,'(',num2str(i),'/',num2str(sum(I)),') ',num2str(round(tot_T/60)),' mins.']);
    
    % Check if dataset available
    FS_path = fullfile(Conf.Root{i},Conf.FS{i});
    
    % Extract all subject & Sessions pair
    Subj= split(Conf.Subj{i}, Conf_Sep);
    Ses = split(Conf.Ses{i} , Conf_Sep);

    % All subjects available on the dataset
    IDs = unique(Subj);
    for j=1:length(IDs)
        % Meassure the time of a single itteration
        tic;

        % Update Waitbar
        H.Name = string([DtSt_nm,'(',num2str(i),'/',num2str(sum(I)),') ',num2str(round(tot_T,2)),' mins.']);
        waitbar(j/length(IDs),H,['Progress: ',num2str(j),'/',num2str(length(IDs))]);

        % Path to subject data
        Sub_path = fullfile(FS_path,IDs{j});

        % Instances of this subect (Number of Sessions)
        jdx = find(strcmp(Subj,IDs{j}));
        flag = 0;
        for k=1:length(jdx)
            % Add the subject data path
            if ~strcmp(Ses{jdx(k)},'_*_')
                DT_path = fullfile(Sub_path,Ses{jdx(k)});
            else
                DT_path = Sub_path;
            end

            % Go for each selected hemisphere
            for l=1:length(Hemi_MD)
                % Path to the pial-smooth file
                Smt_path = fullfile(DT_path,'surf',[char(Hemi_MD(l)),'.pial-outer-smoothed']);

                % Evaluate if the 'smooth' file is present adn it works
                try
                    evalc("[opialv,opialf] = freesurfer_read_surf(Smt_path)");
                catch me
                %if ~isfile(Smt_path)
                    % Trigger flag that at least 1 of the smooth surfaces
                    % not ok/missing
                    flag = flag+1;

                    % Add the path of the subject to the set
                    p_nt = p_nt+1;
                    path{p_nt} = fullfile(DT_path,'surf',char(Hemi_MD(l)));
                end
                % Pathc the memory leaking to the RAM
                clear functions;
            end
        end
        
        % If all the required files for the subject present, +1
        if ~flag
            Smth_cnt(cnt) = Smth_cnt(cnt)+1;
        end

        % Estimate the time it took to run a single subject
        t(end+1) = toc;
        s(i) = s(i)-1;
        tot_T = sum(s*mean(t)/60);

    end
end
path = unique(path);

close(H);

%clear functions