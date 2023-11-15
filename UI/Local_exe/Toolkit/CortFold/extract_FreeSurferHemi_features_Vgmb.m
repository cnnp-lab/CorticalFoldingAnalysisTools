
function [tbl, corrupt] ...
    = extract_FreeSurferHemi_features_Vgmb(subjdir,varargin)
% This function extracts the numeric data from FreeSurfer subjects into
% one table on a hemisphere-basis, which is returned and can optionally be saved.
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

%% parse options
p = inputParser;
p.addParameter('hemi',   'both',  @(x) isstring(x) | ischar(x));
p.addParameter('verbose', true,   @(x) islogical(x));
p.parse(varargin{:});
param = p.Results;

% type security (for input string/char shouldn't matter, internally it does)
param.hemi   = string(param.hemi);

%% make table
tbl = table();

% hemisphere handling
switch(param.hemi)
    case "left"
        lrstr = 'l';
        leftright = ["left"];
    case "right"
        lrstr = 'r';
        leftright = ["right"];
    otherwise % avg, sum, both
        lrstr = 'lr';
        leftright = ["left" "right"];
end

corrupt = zeros(1,length(lrstr)); % Inidicated if the ID is corrupted or not

%% Load subject data

for lr = 1:length(lrstr)
    % load all relevant files
    pathpre = [subjdir, '/surf/', lrstr(lr)];
    try
        if param.verbose
            [thickness, ~]  = read_curv([pathpre, 'h.thickness']);
            [pialv,pialf]   = freesurfer_read_surf([pathpre, 'h.pial']);
            [whitev,whitef] = freesurfer_read_surf([pathpre, 'h.white']);
            [opialv,opialf] = freesurfer_read_surf([pathpre, 'h.pial-outer-smoothed']);
        else % suppress all output from lib scripts (except errors/warnings)
            evalc("[thickness, ~]  = read_curv([pathpre, 'h.thickness'])");
            evalc("[pialv,pialf]   = freesurfer_read_surf([pathpre, 'h.pial'])");
            evalc("[whitev,whitef] = freesurfer_read_surf([pathpre, 'h.white'])");
            evalc("[opialv,opialf] = freesurfer_read_surf([pathpre, 'h.pial-outer-smoothed'])");
        end
    catch me
        corrupt(lr) = 1; % store the corrupt Hemisphere files
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
end

%% add data to table

row_template = table('Size',[1 19],'VariableTypes',[repmat({'string'},1,2),repmat({'double'},1,17)],...
    'VariableNames',{'Hemisphere','Region',...
    'AvgCortThickness','PialArea','SmoothPialArea','WhiteArea','GaussianCurvature','ConvexHullArea',...
    'PialFullArea','WhiteFullArea','SmoothPialFullArea','ConvexHullFullArea',...
    'PialVol','WhiteVol','GreymatterVol',...
    'PialFullVol','WhiteFullVol','GreymatterFullVol','SmoothPialFullVol'});

if param.hemi == "avg" || param.hemi == "sum"
        if param.hemi == "avg"
            fun = mean;
        else
            fun = sum;
        end

        row = row_template;
        row.Hemisphere = param.hemi;
        row.Region = 'hemisphere';

        if sum(corrupt)==0
            row.AvgCortThickness   = fun(thicknessPial);
            row.PialArea           = fun(pialFullArea-pialCCArea);
            row.SmoothPialArea     = fun(opialFullArea-opialCCArea);
            row.WhiteArea          = fun(whiteFullArea-whiteCCArea);
            row.GaussianCurvature  = nan;
            row.ConvexHullArea     = fun(CHFullArea-opialCCArea);

            row.PialFullArea       = fun(pialFullArea);
            row.WhiteFullArea      = fun(whiteFullArea);
            row.SmoothPialFullArea = fun(opialFullArea);
            row.ConvexHullFullArea = fun(CHFullArea);

            row.PialVol           = nan;
            row.WhiteVol          = nan;
            row.GreymatterVol     = nan;

            row.PialFullVol       = fun(pialFullVol);
            row.WhiteFullVol      = fun(whiteFullVol);
            row.GreymatterFullVol = fun(GMVol);
            row.SmoothPialFullVol = fun(opialFullVol);
        end

        tbl = [tbl; row];

else % for "both", "right" and "left"
    for lr = 1:length(lrstr)
        
        row = row_template;
        row.Hemisphere = leftright(lr);
        row.Region = 'hemisphere';

        if corrupt(lr) == 0
            row.AvgCortThickness   = thicknessPial(lr);
            row.PialArea           = pialFullArea(lr)-pialCCArea(lr);
            row.SmoothPialArea     = opialFullArea(lr)-opialCCArea(lr);
            row.WhiteArea          = whiteFullArea(lr)-whiteCCArea(lr);
            row.GaussianCurvature  = nan;
            row.ConvexHullArea     = CHFullArea(lr)-opialCCArea(lr);

            row.PialFullArea       = pialFullArea(lr);
            row.WhiteFullArea      = whiteFullArea(lr);
            row.SmoothPialFullArea = opialFullArea(lr);
            row.ConvexHullFullArea = CHFullArea(lr);

            row.PialVol           = nan;
            row.WhiteVol          = nan;
            row.GreymatterVol     = nan;

            row.PialFullVol       = pialFullVol(lr);
            row.WhiteFullVol      = whiteFullVol(lr);
            row.GreymatterFullVol = GMVol(lr);
            row.SmoothPialFullVol = opialFullVol(lr);
        end

        tbl = [tbl; row];
    end
end

end





