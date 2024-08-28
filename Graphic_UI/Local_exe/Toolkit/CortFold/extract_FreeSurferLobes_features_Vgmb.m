
function [tbl, corrupt] ...
    = extract_FreeSurferLobes_features_Vgmb(subjdir, varargin)

% THE ENTIRE DESCRIPTION NEED A REVISION ONCE IT IS ADAPTED TO THE NEW PIPELINE

% This function extracts the numeric data from FreeSurfer subjects into
% one table on a lobe-basis, which is returned and can optionally be saved.
% The table format is chosen so that it can easily be ID-joined with a
% metadata table and merged with further subjects later.
%
% The following columns are included:
% - SubjID, Hemisphere and Lobe function as an unique index
% - Lobe is a numeric value: 0 1 2 3 4 5 are for Corpus Callosum, Frontal lobe,
%   Parietal lobe, Temporal lobe, Occipital lobe and Insula respectively
% the remaining columns are the extracted features:
% - AvgCortThickness: average pial thickness as from ?h.thickness,
%   corrected for areas towards the Corpus callosum (CC; vertices where T ~ 0) and
%   based on a area-weighted estimate. See Wang 2016 PNAS for details. This
%   measure is different from the FS estimated average thickness, and tends
%   to be systematically higher. Does not impact the scaling behaviour
%   though (as it is a systematic difference).
% - PialArea: from ?h.pial
% - SmoothPialArea: from ?h.pial-outer-smoothed
% - WhiteArea: from ?h.white
% - GaussianCurvature: the Gaussian curvature of the smoothed (downsampled)
%   convex hull (without boundaries), i.e. a cap around the lobe; it is
%   crucial to correct the other measures for curvature when comparing relative lobe sizes
%   since curvature captures "how much" of a sphere is captured in a lobe
%   (Gauss-Bonnet). See Wang 2019 Comms Biol for details. Note also that
%   the total curvature of all the lobes will sum up to more than 4*pi as
%   we integrate every lobe's exposed surface separately without their boundaries.
%   However the correction term is still 4*pi for everyone, as we correct
%   every lobe back to a sphere. I.e the sphere is the reference point.
% - PialVol: Volume between the ?h.pial mesh and the origin (like a cone)
% - WhiteVol: Volume between the ?h.white mesh and the origin (like a cone)
% - GreymatterVol: The difference of the former two
%
% Returns:
% The table as described above and a list of IDs of corrupt subjects, i.e.
% such that have a missing FreeSurfer file, see above.
% Note that "lobe" 0 (corpus callosum) and 5 (insula) are generally not
% used. You can choose to include the insula PialArea in the other lobes
% (see Wang 2019 Comms Biol), but it doesn't really make a huge difference.
% Ignore any thickness & curvature estimates for the corpus callosum.
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
% ISO2MESH libaray needs to be downloaded&added to lib!!!: http://iso2mesh.sourceforge.net/
% Licence: CC-BY
%
% Yujiang Wang, September 2016 (extractMaster_Hemisphere.m)
% Tobias Ludwig & Yujiang Wang, January 2020
% Newcastle University, School of Computing, CNNP Lab (www.cnnp-lab.com)

%% parse options
p = inputParser;
p.addParameter('hemi',   'both',  @(x) isstring(x) | ischar(x));
p.addParameter('verbose', true,   @(x) islogical(x));
p.addParameter('atlas',   'FSDK',  @(x) isstring(x) | ischar(x));
p.parse(varargin{:});
param = p.Results;

% type security (for input string/char shouldn't matter, internally it does)
param.hemi   = string(param.hemi);
param.atlas  = string(param.atlas);

% load look-up-table for FS lobe labels vs the labels we use here (0-5):
switch param.atlas
    case "FSDK"
        aux = load('FSDK.mat');
        Atlas = aux.Map;

        % Lobes to considere
        [aux,idx,~] = unique(Atlas.Lobe_cd);
        LB_dx = aux';

        % Lobes namings for LB_dx
        %LB_nm = string(LB_dx);
        LB_nm = Atlas.Lobe_nm(idx);
    otherwise

end

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
tic % start timer

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
    pathpre = [subjdir, '/surf/', lrstr(lr)];
    try
        if param.verbose
            [thickness, ~]  = read_curv([pathpre, 'h.thickness']);
            [pialv,pialf]   = freesurfer_read_surf([pathpre, 'h.pial']);
            [whitev,whitef] = freesurfer_read_surf([pathpre, 'h.white']);
            [opialv,opialf] = freesurfer_read_surf([pathpre, 'h.pial-outer-smoothed']);
            [~,label,~]     = read_annotation([subjdir, '/label/', lrstr(lr), 'h.aparc.annot']);
        else % suppress all output from lib scripts (except errors/warnings)
            evalc("[thickness, ~]  = read_curv([pathpre, 'h.thickness'])");
            evalc("[pialv,pialf]   = freesurfer_read_surf([pathpre, 'h.pial'])");
            evalc("[whitev,whitef] = freesurfer_read_surf([pathpre, 'h.white'])");
            evalc("[opialv,opialf] = freesurfer_read_surf([pathpre, 'h.pial-outer-smoothed'])");
            evalc("[~,label,~]     = read_annotation([subjdir, '/label/', lrstr(lr), 'h.aparc.annot'])");
        end
    catch me
        corrupt(lr) = 1; % store the corrupt Hemisphere files
        continue
    end

    %% match labels to smooth surface
    label_smooth = matchSurfLabel(label,pialv,opialv);

    %% relabel for lobes
    nlabel_total_lobe  = label;
    nlabel_smooth_lobe = label_smooth;
    for k = LB_dx
        idtoget = Atlas.ROI_cd(Atlas.Lobe_cd==k);

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


    AvgThickness(:,lr)= AvgThicknessLR;
    TotalArea(:,lr)   = TotalAreaLR;
    SmoothArea(:,lr)  = SmoothAreaLR;
    WhiteArea(:,lr)   = WhiteAreaLR;
    PialVol(:,lr)     = PialVolLR;
    WhiteVol(:,lr)    = WhiteVolLR;

    %% calculate gaussian curvature (K) from convex hull (CV) of downsampled opial
    % downsample Ae
    chr = 0.2; % keepratio - downsample to 20% of original mesh to speed up computation later.
    [opialv_ds,opialf_ds] = meshresample(opialv,opialf,chr);

    % reassign labels
    nlabel_smooth_lobe_ds = matchSurfLabel(nlabel_smooth_lobe,opialv,opialv_ds);

    KLR_from_CH_of_ds_opial=zeros(6,1);
    Area_from_CH_of_ds_opial=zeros(6,1);

    for l = LB_dx
        ov_ids = find(nlabel_smooth_lobe_ds==l);
        Lobepoints = opialv_ds(ov_ids,:);
        if length(Lobepoints) > 4
            CHSf_ds = convhull(Lobepoints);
            Gru = getGaussianCurvPart(CHSf_ds,Lobepoints,1:length(ov_ids));

            isBoundary = zeros(length(ov_ids),1);
            for kl = 1:length(ov_ids)
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

    for lobe = (LB_dx+1)
        row = table();
        row.Hemisphere = param.hemi;
        row.Region = LB_nm(lobe);
        row.Scale = 0;

        row.AvgCortThickness  = fun(AvgThickness(lobe,:));
        row.PialArea          = fun(TotalArea(lobe,:));
        row.SmoothPialArea    = fun(SmoothArea(lobe,:));
        row.WhiteArea         = fun(WhiteArea(lobe,:));
        row.GaussianCurvature = fun(K_CHwoB(lobe,:));
        row.ConvexHullArea    = nan;

        row.PialFullArea      = nan;
        row.WhiteFullArea     = nan;
        row.SmoothPialFullArea= nan;
        row.ConvexHullFullArea= nan;

        row.PialVol           = fun(PialVol(lobe,:));
        row.WhiteVol          = fun(WhiteVol(lobe,:));
        row.GreymatterVol     = fun(PialVol(lobe,:) - WhiteVol(lobe,:));

        row.PialFullVol       = nan;
        row.WhiteFullVol      = nan;
        row.GreymatterFullVol = nan;
        row.SmoothPialFullVol = nan;

        tbl = [tbl; row];
    end

else % for "both", "right" and "left"
    for lobe = (LB_dx+1)
        for lr = 1:length(lrstr)
            row = table();
            row.Hemisphere = leftright(lr);
            row.Region = LB_nm(lobe);
            row.Scale = 0;

            row.AvgCortThickness  = AvgThickness(lobe,lr);
            row.PialArea          = TotalArea(lobe,lr);
            row.SmoothPialArea    = SmoothArea(lobe,lr);
            row.WhiteArea         = WhiteArea(lobe,lr);
            row.GaussianCurvature = K_CHwoB(lobe,lr);
            row.ConvexHullArea    = nan;

            row.PialFullArea      = nan;
            row.WhiteFullArea     = nan;
            row.SmoothPialFullArea= nan;
            row.ConvexHullFullArea= nan;

            row.PialVol           = PialVol(lobe,lr);
            row.WhiteVol          = WhiteVol(lobe,lr);
            row.GreymatterVol     = PialVol(lobe,lr) - WhiteVol(lobe,lr);


            row.PialFullVol      = nan;
            row.WhiteFullVol     = nan;
            row.GreymatterFullVol= nan;
            row.SmoothPialFullVol= nan;


            tbl = [tbl; row];
        end
    end
end

end
