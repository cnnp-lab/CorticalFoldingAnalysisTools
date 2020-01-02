
function [tbl, corrupt_ids] ...
= extract_FreeSurferLobes_features(subjdir, ids, oldtbl, varargin)
% This function extracts the numeric data from FreeSurfer subjects into
% one table on a lobe-basis, which is returned and can optionally be saved.
% The table format is chosen so that it can easily be ID-joined with a
% metadata table and merged with further subjects later.
%
% The following columns are included:
% - SubjID, Hemisphere and Lobe function as an unique index
% - Lobe is a numeric value: 0 1 2 3 4 5 are CC F P T O Insula respectively
% the remaining columns are the extracted features:
% - AvgCortThickness: average pial thickness as from ?h.thickness,
%   corrected for areas towards the Corpus callosum (CC; vertices where T ~ 0) and
%   based on a area-weighted estimate. See Wang 2016 PNAS for details. This
%   measure is different from the FS estimated average thickness, and tends
%   to be systematically higher. Does not impact the scaling behaviour
%   though (as it is a systematic difference).
% - PialArea: corrected for CC areas; from ?h.pial; note, again different
%   from FS output, as FS includes CC area
% - SmoothPialArea: corrected for CC areas; as from ?h.pial-outer-smoothed
% - WhiteArea: corrected for CC areas; as from ?h.white
% - ConvexHullArea: the convex hull of the pial is calculated in matlab
%   by the convhull-function; corrected for CC areas
% - PialFullArea: not corrected for CC areas
% - WhiteFullArea: not corrected for CC areas
% - SmoothPialFullArea: not corrected for CC areas
% - ConvexHullFullArea: not corrected for CC areas
% - PialFullVol: Volume of the closed ?h.pial mesh
% - WhiteFullVol: Volume of the closed ?h.white mesh
% - GreymatterVol: The difference of the former two
% - SmoothPialFullVol: Volume of the closed ?h.pial-outer-smoothed mesh
%
% Returns:
% The table as described above and a list of IDs of corrupt subjects, i.e.
% such that have a missing FreeSurfer file, see above.
%
% Arguments:
% - SUBJDIR: char or string; path to the folder that contains the
%   FreeSurfer subjects (one folder per subject)
%   (the absolute path with a trailing slash has to be given!)
% - IDS: string array; IDs of subjects to include in the table
%   (or to update from oldtbl, see below)
% - OLDTBL: table; the old table will be merged with the new table
%   (replacing existing IDs; assuming the same columns to be present);
%   if no merging is wanted, give an empty table; see also option 'replace'
%
% Options:
% - 'format': char or string (default ''); format of folder names per subject,
%   e.g. 'Patient_*' or 'C_*.gz', where '*' will be replaced by the IDs
% - 'replace': if true (default), replace rows of existing IDs in 'oldtbl'
%   with the new information; if false, ignore those IDs while loading
% - 'hemi': how to treat hemispheres; 'avg' for averaging, 'sum' for summing
%   (caution: might make sense for area/volume but not for e.g. thickness),
%   'left' or 'right' for including only one of them respectively,
%   or 'both' for having separate rows for left and right hemispheres.
%   The column 'Hemisphere' indicates 'avg', 'sum' or 'left'/'right'.
%   'both' is the default option.
% - 'verbose': if true (default), output from the lib scripts is printed
% - 'saveas': filename to store the table as, <saveas>.mat
% - 'libdir': string; path to the lib-folder including the scripts for
%   extraction; by default it is assumed to be a subfolder of this
%
% Dependencies:
% The lib-folder should be in the same directory as this file,
% otherwise the option 'libdir' has to point to it.
%
% Licence: CC-BY
% 
% Yujiang Wang, September 2016 (extractMaster_Hemisphere.m)
% Tobias Ludwig & Yujiang Wang, September 2019
% Newcastle University, School of Computing, CNNP Lab (www.cnnp-lab.com)

%% version
VERSION_ID = 1; % for furture compatibility checks; not implemented yet
N_FEATURES = 13; % number of features = number of table columns - 2


%% parse options
p = inputParser;
p.addParameter('format', '',      @(x) isstring(x) | ischar(x));
p.addParameter('hemi',   'both',  @(x) isstring(x) | ischar(x));
p.addParameter('replace', true,   @(x) islogical(x));
p.addParameter('verbose', true,   @(x) islogical(x));
p.addParameter('saveas', '',      @(x) isstring(x) | ischar(x));
p.addParameter('libdir', './lib', @(x) isstring(x) | ischar(x));
p.parse(varargin{:});
param = p.Results;

% type security (for input string/char shouldn't matter, internally it does)
param.format = string(param.format);
param.hemi   = string(param.hemi);
param.saveas = char(param.saveas);
param.libdir = char(param.libdir);

% check argument validity
if ~isstring(ids(1))
    error(['The IDs have to be given as string array! Got ' ids(1)]);
end
if size(oldtbl, 2)>0 && size(oldtbl, 2) ~= N_FEATURES + 2
    error(['Merging with old table impossible, number of columns mismatch.' ...
    ' Check with the documentation which features are currently supported.'])
end

addpath(param.libdir)
addpath([param.libdir '/FSmatlab/'])
addpath([param.libdir '/iso2mesh/'])

% load look-up-table for FS lobe labels vs the labels we use here (0-5):
load(['LUT_lobes.mat'])

%% make table
tbl = table();

% hemisphere handling
switch(param.hemi)
    case "left"
        lrstr = 'l';
    case "right"
        lrstr = 'r';
    otherwise % avg, sum, both
        lrstr = 'lr';
end
leftright = ["left" "right"];

% if oldtbl given and replace == false, only load the new ids
if ~isempty(oldtbl) && ~param.replace
    ids = setdiff(unique(ids), unique(oldtbl.SubjectID));
end

corrupt_ids = strings(); % empty string array = [""], first will be deleted
corrupt = false;

%% Load subject data
tic % start timer
for iter = 1:length(ids)
    % make foldername
    id = ids(iter);
    if param.format == ""
        fn = char(id);
    else
        fn = char(strrep(param.format, '*', id));
    end
    
    % init each feature for all 6 lobes and two hemispheres
    z = zeros(6,2);
    AvgThickness = z;
    TotalArea = z;
    SmoothArea = z;
    WhiteArea = z;
    PialVol = z;
    WhiteVol = z;
    K_CHwoB = z; % gaussian curvature
    
    for lr = 1:length(lrstr)
        % load all relevant files
        pathpre = [subjdir, fn, '/surf/', lrstr(lr)];
        try
            if param.verbose
                [thickness, ~]  = read_curv([pathpre, 'h.thickness']);
                [pialv,pialf]   = freesurfer_read_surf([pathpre, 'h.pial']);
                [whitev,whitef] = freesurfer_read_surf([pathpre, 'h.white']);
                [opialv,opialf] = freesurfer_read_surf([pathpre, 'h.pial-outer-smoothed']);
                % original version: [opialv,opialf]= freesurfer_read_surf([pathpre, POSname]); %% TODO what is POSname???
                [~,label,~]     = read_annotation([subjdir, fn, '/label/', lrstr(lr), 'h.aparc.annot']);
            else % suppress all output from lib scripts (except errors/warnings)
                evalc("[thickness, ~]  = read_curv([pathpre, 'h.thickness'])");
                evalc("[pialv,pialf]   = freesurfer_read_surf([pathpre, 'h.pial'])");
                evalc("[whitev,whitef] = freesurfer_read_surf([pathpre, 'h.white'])");
                evalc("[opialv,opialf] = freesurfer_read_surf([pathpre, 'h.pial-outer-smoothed'])");
                evalc("[~,label,~]     = read_annotation([subjdir, fn, '/label/', lrstr(lr), 'h.aparc.annot'])");
            end
        catch me
            % if one subject is corrupt, we don't want to fail and stop,
            % but instead we print a warning and skip this subject
            warning(['A file of subject ' char(id) ' was not found: ' me.message ...
                '\nThe whole subject is skipped and excluded from the table.']);
            corrupt_ids(end+1) = id; % store the corrupt id
            corrupt = true; % flag to advance outer subject for loop
            continue
        end
        
        %% match labels to smooth surface
        label_smooth = matchSurfLabel(label,pialv,opialv);
        
        %% relabel for lobes
        nlabel_total_lobe  = label;
        nlabel_smooth_lobe = label_smooth;
        for k = 0:5
            idtoget = LUT_lobes(LUT_lobes(:,2)==k, 1);

            for l = 1:length(idtoget)
                nlabel_total_lobe(nlabel_total_lobe==idtoget(l))   = k;
                nlabel_smooth_lobe(nlabel_smooth_lobe==idtoget(l)) = k;
            end
        end
        
        %% extract measures for lobes
        [AvgThicknessLR, TotalAreaLR, SmoothAreaLR, WhiteAreaLR, PialVolLR, WhiteVolLR, ~]= ...
            getBasicMeasuresForParts(nlabel_total_lobe, ...
                                     nlabel_smooth_lobe, ...
                                     pialf,pialv,opialf,opialv, ...
                                     thickness,whitef,whitev);
        
        if corrupt
            corrupt = false; % toggle flag
            continue; % skip this subject
        end
        
        AvgThickness(:,lr)= AvgThicknessLR;
        TotalArea(:,lr)   = TotalAreaLR;
        SmoothArea(:,lr)  = SmoothAreaLR;
        WhiteArea(:,lr)   = WhiteAreaLR;
        PialVol(:,lr)     = PialVolLR;
        WhiteVol(:,lr)    = WhiteVolLR;
        
        %% calculate gaussian curvature (K) from convex hull (CV) of downsampled opial
        % downsample Ae
        chr = 0.2; % keepratio TODO rename? why so low?
        [opialv_ds,opialf_ds] = meshresample(opialv,opialf,chr);
        
        % reassign labels
        nlabel_smooth_lobe_ds = matchSurfLabel(nlabel_smooth_lobe,opialv,opialv_ds);
        
        KLR_from_CH_of_ds_opial=zeros(6,1);
        Area_from_CH_of_ds_opial=zeros(6,1);
        
        for l = 0:5
            ov_ids = find(nlabel_smooth_lobe_ds==l);
            Lobepoints = opialv_ds(ov_ids,:);
            if length(Lobepoints) > 4
                CHSf_ds = convhull(Lobepoints);
                Gru = getGaussianCurvPart(CHSf_ds,Lobepoints,1:length(ov_ids));

                isBoundary = zeros(length(ov_ids),1);
                for kl = 1:length(ov_ids)
                    %where does idFL pop up in opialf?
                    fid = opialf_ds(:,1)==ov_ids(kl) | ...
                          opialf_ds(:,2)==ov_ids(kl) | ...
                          opialf_ds(:,3)==ov_ids(kl);

                    neighbourIDs = unique(opialf_ds(fid,:));
                    ncolors = numel(unique(nlabel_smooth_lobe_ds(neighbourIDs)));
                    if ncolors > 1
                        isBoundary(kl) = 1;
                    end
                end

                B_ids = find(isBoundary==1);
                [TotalAreaCapi,~] = calcPartAreai(CHSf_ds,Lobepoints,1:size(Lobepoints,1));
                Btri = ismember(CHSf_ds(:,1),B_ids) & ...
                       ismember(CHSf_ds(:,2),B_ids) & ...
                       ismember(CHSf_ds(:,3),B_ids);
                
                KLR_from_CH_of_ds_opial(l+1)  = sum(Gru(isBoundary==0));
                Area_from_CH_of_ds_opial(l+1) = sum(TotalAreaCapi(~Btri));
            end
        end
        
        K_CHwoB(:,lr) = KLR_from_CH_of_ds_opial;
    end
    
    
    %% add data to table
    if param.hemi == "avg" || param.hemi == "sum"
        if param.hemi == "avg"
            fun = mean;
        else
            fun = sum;
        end
        
        for lobe = 1:6
            row = table();
            row.SubjectID = id;
            row.Hemisphere = param.hemi;
            row.Lobe = lobe - 1;
            
            row.AvgCortThickness  = fun(AvgThickness(lobe,:));
            row.PialArea          = fun(TotalArea(lobe,:));
            row.SmoothPialArea    = fun(SmoothArea(lobe,:));
            row.WhiteArea         = fun(WhiteArea(lobe,:));
            row.GaussianCurvature = fun(K_CHwoB(lobe,:));
            
            row.PialVol           = fun(PialVol(lobe,:));
            row.WhiteVol          = fun(WhiteVol(lobe,:));
            row.GreymatterVol     = fun(PialVol(lobe,:) - WhiteVol(lobe,:));
            tbl = [tbl; row];
        end
    
    else % for "both", "right" and "left"
        for lobe = 1:6
            for lr = 1:length(lrstr)
                row = table();
                row.SubjectID = id;
                row.Hemisphere = leftright(lr);
                row.Lobe = lobe - 1;

                row.AvgCortThickness  = AvgThickness(lobe,lr);
                row.PialArea          = TotalArea(lobe,lr);
                row.SmoothPialArea    = SmoothArea(lobe,lr);
                row.WhiteArea         = WhiteArea(lobe,lr);
                row.GaussianCurvature = K_CHwoB(lobe,lr);

                row.PialVol           = PialVol(lobe,lr);
                row.WhiteVol          = WhiteVol(lobe,lr);
                row.GreymatterVol     = PialVol(lobe,lr) - WhiteVol(lobe,lr);
                
                tbl = [tbl; row];
            end
        end
    end
end

if param.verbose
    toc % print stopwatch
end

corrupt_ids = corrupt_ids(2:end); % cut off the first ""
if ~isempty(corrupt_ids)
    disp("The following subjects were excluded from the table as corrupt:")
    disp(corrupt_ids(2:end))
end

%% merging old table and saving (if wished)

% save before merging as well in case of merge conflicts!
if ~isempty(param.saveas) % TODO only merge security if saveas given?
    save([param.saveas '.mat'], 'tbl');
    if param.verbose
        disp(['Wrote table to ' param.saveas ' before merging.']),
    end
end
    
if ~isempty(oldtbl)
    tbl = [tbl; oldtbl]; % stack the new table on top of the old
    [~,ia,~] = unique(tbl(:, 1:3), 'rows'); % find (SubjID,Hemisphere,Lobe)-rows
    tbl = tbl(ia, :);    % of these, take the first occurence (the new one)
end
if param.verbose
    disp('Merge with old table successful.');
end
    
if ~isempty(param.saveas)
    save([param.saveas '.mat'], 'tbl');
    if param.verbose
        disp(['Wrote table to ' param.saveas ' after merging.']),
    end
end

end