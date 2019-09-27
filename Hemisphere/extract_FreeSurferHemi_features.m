%% TODO
% - error handling (file not found, ...)
% - option to not overwrite IDs from oldtbl but to exclude them from loading
% - options for output: verbose loading, progress bar?, none
% - save before merging to prevent data loss in merge conflicts
% - prevent merge conflicts by checking the table dimensions


function tbl = extract_FreeSurferHemi_features(subjdir, ids, format, varargin)
% This function extracts the numeric data from FreeSurfer subjects into
% one table, which is returned and can optionally be saved.
% The table format is chosen so that it can easily be ID-joined with a
% metadata table and merged with further subjects later.
%
% The following columns are included:
% - SubjID and Hemisphere function as an index
% the features are:
% - TODO explain features
%
%
%
% Arguments:
% - 'subjdir': string; path to the folder that contains the
%   FreeSurfer subjects (one folder per subject)
% - 'ids': string array; IDs of subjects to include in the table (or to
%   update, if oldtbl is given, see below)
% - 'format': string; format of the folder names per subject, e.g.
%   'Patient_*' or 'C_*.gz', where '*' will be replaced by the IDs

% Options:
% - 'libdir': string; path to the lib-folder including the scripts for
%   extraction; by default it is assumed to be a subfolder of this
% - 'oldtbl': if given, the old table will be merged with the new table
%   (replacing existing IDs; assuming the same columns to be present) % TODO what if more/different cols?
% - 'saveas': filename to store the table as, <saveas>.mat
% - 'hemi': how to handle hemispheres; 'avg' for averaging, 'sum' for summing
%   (caution: might make sense for area/volume but not for e.g. thickness),
%   'left' or 'right' for including only one of them respectively,
%   or 'both' for having separate rows for left and right hemispheres.
%   The column 'Hemisphere' indicates 'avg', 'sum' or 'left'/'right'.
%   'both' is the default option.
%
% Dependencies:
% The lib-folder should included as a subdirectory of the directory of this
% file.
%
% Licence: CC-BY
% 
% Yujiang Wang, September 2016 (extractMaster_Hemisphere.m)
% Tobias Ludwig, September 2019
% Newcastle University, School of Computing, CNNP Lab (www.cnnp-lab.com)


%% parse options TODO
p = inputParser;
p.addParameter('libdir', './lib', @(x) isstring(x) | ischar(x));
p.addParameter('oldtbl', table(), @(x) istable(x));
p.addParameter('saveas', '',      @(x) isstring(x) | ischar(x));
p.addParameter('hemi',   'both',  @(x) isstring(x));
p.parse(varargin{:});
param = p.Results;

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
leftright = ["left" "right"]; % XXX

%% Load subject data
tic
iter = 1;
while iter <= length(ids) % weird matlab behaviour: cannot for-iterate over string array, hence using while
    % make foldername
    id = ids(iter);
    iter = iter + 1;
    fn = char(strrep(format, '*', id));
    
    % check if 'SmoothPialArea' exists % TODO necessary???
    if exist([subjdir, fn, '/surf/rh.pial-outer-smoothed'], 'file') ~= 2
        warning(['File not found: ' subjdir, fn, '/surf/rh.pial-outer-smoothed. ' ...
            'Subject ' char(id) ' will be excluded.']); % TODO the whole subject?
        continue;
    end
    
    for lr = 1:length(lrstr)
        % load all relevant files
        [thickness, ~]  = read_curv([subjdir, fn, '/surf/', lrstr(lr), 'h.thickness']);
        [pialv,pialf]   = freesurfer_read_surf([subjdir, fn, '/surf/', lrstr(lr), 'h.pial']);
        [whitev,whitef] = freesurfer_read_surf([subjdir, fn, '/surf/', lrstr(lr), 'h.white']);
        [opialv,opialf] = freesurfer_read_surf([subjdir, fn, '/surf/', lrstr(lr), 'h.pial-outer-smoothed']);
        
        % calculate full areas --------------------------------------------
        
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
        GMVol(lr)=pialFullVol(lr)-whiteFullVol(lr); % TODO read from FS
        opialFullVol(lr)=calcMeshVol(opialf,opialv);
    end
    
    
    %% add data to table
    if param.hemi == "avg" | param.hemi == "sum"
        if param.hemi == "avg"
            fun = mean;
        else
            fun = sum;
        end
        
        row = table();
        row.SubjectID = id;
        row.Hemisphere = param.hemi;
            
        row.AvgCortThickness  = fun(thicknessPial);
        row.PialArea          = fun(pialFullArea-pialCCArea);
        row.SmoothPialArea    = fun(opialFullArea-opialCCArea);
        row.WhiteArea         = fun(whiteFullArea-whiteCCArea);
           
        row.PialFullArea      = fun(pialFullArea);
        row.WhiteFullArea     = fun(whiteFullArea);
        row.SmoothPialFullArea= fun(opialFullArea);
        row.ConvexHullFullArea= fun(CHFullArea);
        row.ConvexHullArea    = fun(CHFullArea-opialCCArea);
            
        row.PialFullVol       = fun(pialFullVol);
        row.WhiteFullVol      = fun(whiteFullVol);
        row.GreymatterVol     = fun(GMVol);
        row.SmoothPialFullVol = fun(opialFullVol);
            
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
            
            row.PialFullArea       = pialFullArea(lr);
            row.WhiteFullArea      = whiteFullArea(lr);
            row.SmoothPialFullArea = opialFullArea(lr);
            row.ConvexHullFullArea = CHFullArea(lr);
            row.ConvexHullArea     = CHFullArea(lr)-opialCCArea(lr);
            
            row.PialFullVol        = pialFullVol(lr);
            row.WhiteFullVol       = whiteFullVol(lr);
            row.GreymatterVol      = GMVol(lr);
            row.SmoothPialFullVol  = opialFullVol(lr);
            
            tbl = [tbl; row];
        end
    end
end
toc

%% merging old table (if existent) and saving (if wished)

% TODO save before merging as well? in case of merge conflicts!
if ~isempty(param.saveas)
    save([param.saveas '.mat'], 'tbl');
end

if ~isempty(param.oldtbl)
    tbl = [tbl; param.oldtbl]; % stack the new table on top of the old
    [~,ia,~] = unique(tbl(:, 1:2), 'rows'); % find unique (SubjID,Hemisphere)-rows
    tbl = tbl(ia, :);    % of these, take the first occurence (the new one)
end

if ~isempty(param.saveas)
    save([param.saveas '.mat'], 'tbl');
end
end