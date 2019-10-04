%% TODO
% - prevent merge conflicts by checking the table dimensions
% - make format an option and oldtbl an argument (?)

function [tbl, corrupt_ids] ...
= extract_FreeSurferHemi_features(subjdir, ids, oldtbl, varargin)
% This function extracts the numeric data from FreeSurfer subjects into
% one table, which is returned and can optionally be saved.
% The table format is chosen so that it can easily be ID-joined with a
% metadata table and merged with further subjects later.
%
% The following columns are included:
% - SubjID and Hemisphere function as an unique index
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

corrupt_ids = [];
corrupt = false;

%% Load subject data
tic % start timer
iter = 1;
while iter <= length(ids) % weird matlab behaviour: cannot for-iterate over string array, hence using while
    % make foldername
    id = ids(iter);
    iter = iter + 1;
    if param.format == ""
        fn = char(id);
    else
        fn = char(strrep(param.format, '*', id));
    end
    
    for lr = 1:length(lrstr)
        % load all relevant files
        try
            if param.verbose
                [thickness, ~]  = read_curv([subjdir, fn, '/surf/', lrstr(lr), 'h.thickness']);
                [pialv,pialf]   = freesurfer_read_surf([subjdir, fn, '/surf/', lrstr(lr), 'h.pial']);
                [whitev,whitef] = freesurfer_read_surf([subjdir, fn, '/surf/', lrstr(lr), 'h.white']);
                [opialv,opialf] = freesurfer_read_surf([subjdir, fn, '/surf/', lrstr(lr), 'h.pial-outer-smoothed']);
            else % suppress all output from lib scripts (except errors/warnings)
                evalc("[thickness, ~]  = read_curv([subjdir, fn, '/surf/', lrstr(lr), 'h.thickness'])");
                evalc("[pialv,pialf]   = freesurfer_read_surf([subjdir, fn, '/surf/', lrstr(lr), 'h.pial'])");
                evalc("[whitev,whitef] = freesurfer_read_surf([subjdir, fn, '/surf/', lrstr(lr), 'h.white'])");
                evalc("[opialv,opialf] = freesurfer_read_surf([subjdir, fn, '/surf/', lrstr(lr), 'h.pial-outer-smoothed'])");
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
        
        % calculate full areas --------------------------------------------
        tic
        % pial
        pial_area=calcTriangleArea(pialf,pialv);
        pialFullArea(lr)=sum(pial_area);

        % opial (outer pial)
        opial_area=calcTriangleArea(opialf,opialv);
        opialFullArea(lr)=sum(opial_area);

        % convex hull
        convHull = convhull(pialv(:,1),pialv(:,2),pialv(:,3));
        CH_area=calcTriangleArea(convHull,pialv);
        CHFullArea(lr)=sum(CH_area);
        
        % white
        white_area=calcTriangleArea(whitef,whitev);
        whiteFullArea(lr)=sum(white_area);

        % calculate CC areas-----------------------------------------------
        cc=find(thickness==0); % corpus callosum + brain stem + ...
        
        % for pial
        pialCCArea(lr)=calcPartArea(pialf,pialv,cc);

        % for white
        whiteCCArea(lr)=calcPartArea(whitef,whitev,cc);
        
        % for opial
        lblcc=zeros(length(pialv),1);
        lblcc(cc)=1;

        bb = get_bounding_box([pialv;opialv]);
        ratio = 32;
        [parts, partsCount] = give_parts_to_vertices(pialv, bb, ratio);

        labels_adjParts = find_smooth_labels_adjp(opialv, pialv, lblcc, parts, partsCount, bb, ratio);
        opialCCArea(lr) = calcPartArea(opialf,opialv,find(labels_adjParts==1));

        
        % calculate thickness----------------------------------------------
        
        %find out which faces are relevant
        [liaa]=ismember(pialf,cc);
        fid=sum(liaa,2);
        fid_notCC=fid==0;
        
        % make thickness facebased
        thicknessf=makeFacebased(thickness,pialf);
        
        % get weight
        w=pial_area(fid_notCC)./sum(pial_area(fid_notCC)) + white_area(fid_notCC)./sum(white_area(fid_notCC));
        w=w/2;
        
        % get thickness
        thicknessPial(lr)=sum(w.*thicknessf(fid_notCC));
        
        % calculate Volumes------------------------------------------------
        pialFullVol(lr)=calcMeshVol(pialf,pialv);
        whiteFullVol(lr)=calcMeshVol(whitef,whitev);
        GMVol(lr)=pialFullVol(lr)-whiteFullVol(lr); % TODO read from FS (?)
        opialFullVol(lr)=calcMeshVol(opialf,opialv);
        toc
    end
    
    if corrupt
        corrupt = false; % toggle flag
        continue; % skip this subject
    end
    
    %% add data to table
    if param.hemi == "avg" || param.hemi == "sum"
        if param.hemi == "avg"
            fun = mean;
        else
            fun = sum;
        end
        
        row = table();
        row.SubjectID = id;
        row.Hemisphere = param.hemi;
            
        row.AvgCortThickness   = fun(thicknessPial);
        row.PialArea           = fun(pialFullArea-pialCCArea);
        row.SmoothPialArea     = fun(opialFullArea-opialCCArea);
        row.WhiteArea          = fun(whiteFullArea-whiteCCArea);
        row.ConvexHullArea     = fun(CHFullArea-opialCCArea);
        
        row.PialFullArea       = fun(pialFullArea);
        row.WhiteFullArea      = fun(whiteFullArea);
        row.SmoothPialFullArea = fun(opialFullArea);
        row.ConvexHullFullArea = fun(CHFullArea);
            
        row.PialFullVol        = fun(pialFullVol);
        row.WhiteFullVol       = fun(whiteFullVol);
        row.GreymatterVol      = fun(GMVol);
        row.SmoothPialFullVol  = fun(opialFullVol);
            
        tbl = [tbl; row];
    else % for "both", "right" and "left"
        for lr = 1:length(lrstr)
            row = table();
            row.SubjectID = id;
            row.Hemisphere = leftright(lr);
            
            row.AvgCortThickness   = thicknessPial(lr);
            row.PialArea           = pialFullArea(lr)-pialCCArea(lr);
            row.SmoothPialArea     = opialFullArea(lr)-opialCCArea(lr);
            row.WhiteArea          = whiteFullArea(lr)-whiteCCArea(lr);
            row.ConvexHullArea     = CHFullArea(lr)-opialCCArea(lr);
            
            row.PialFullArea       = pialFullArea(lr);
            row.WhiteFullArea      = whiteFullArea(lr);
            row.SmoothPialFullArea = opialFullArea(lr);
            row.ConvexHullFullArea = CHFullArea(lr);
            
            row.PialFullVol        = pialFullVol(lr);
            row.WhiteFullVol       = whiteFullVol(lr);
            row.GreymatterVol      = GMVol(lr);
            row.SmoothPialFullVol  = opialFullVol(lr);
            
            tbl = [tbl; row];
        end
    end
end
if param.verbose
    toc % print stopwatch
end

if ~isempty(corrupt)
    disp("The following subjects were excluded from the table as corrupt:")
    disp(corrupt_ids)
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
    [~,ia,~] = unique(tbl(:, 1:2), 'rows'); % find (SubjID,Hemisphere)-rows
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