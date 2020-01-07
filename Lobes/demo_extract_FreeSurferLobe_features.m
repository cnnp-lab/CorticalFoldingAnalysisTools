%% DEMO FILE
% We show how to load controls and patients into separate tables and then
% merge both tables and join them with a corresponding metadata file.
% Note that in general there is no need to do it separately, but here we do
% to demonstrate the merging process.

%The lib folder and subfolders need to be added to path to work. Note that
%we included the freesurfer matlab functions again in the lib for
%completeness. Please refer to the freesurfer manual for further
%information on this.

%You can see some example outputs in the subfolder demo (as .mat files 
%containing tables).

addpath(genpath('../lib'))

%% load controls
% subjdir contains the FreeSurfer files, one directory per subject
% (the absolute path with a trailing slash has to be given!)
subjdir='/home/tobi/Downloads/Newcastle/CorticalFoldingAnalysisTools/';

% IDs are in the form 'C001', 'C002', ...
% you have to select which subjects you want to calculate by string IDs!
ids = num2str([1:30]', '%03u'); % formating ID (leading zeros)
ids = string([repmat('C', [length(ids),1]), ids]);

% pattern for foldername of each subject
format='*_T1.nii.gz'; % '*' will be replaced by the ID - can leave empty if you directly provide the directory names

% How to handle hemispheres? 'avg' for averaging, 'sum' for summing
% (caution: might make sense for area/volume but not for e.g. thickness),
% 'left' or 'right' for including only one of them respectively,
% or 'both' for having separate rows for left and right hemispheres.
% The column 'Hemisphere' indicates 'avg', 'sum' or 'left'/'right'.
hemi = "both"; %recommended

% run extraction - this is the main function doing the data extraction!
% as we have no old table that should be used as a basis, we give an empty
% table as a new argument
%tbl = extract_FreeSurferLobes_features(subjdir, ids, table(), 'format', format, 'hemi', hemi);

%% load patients and merge with controls
%ids = ["930"];
ids = ["100307"];
% merging is done by specifying the old table to add to/replacing old by new ids,
% saving is not needed but prevents data loss in case of merge errors
% since before the newly extracted data is backuped before merging
[tbl, corrupt] = extract_FreeSurferLobes_features(subjdir, ids, table(), ...
                                            'format', format, 'hemi', hemi, ...
                                            'replace', true, ... % default
                                            'verbose', true, ... % default
                                            'libdir', '../lib' ...
                                            ...'saveas', 'tbl'
                                            );
%% join with metadata

% load metadata table
load('../demo/metadata.mat'); % will load variable 'metadata'

% join table with metadata table by SubjectID (common column in the tables)
tbl_meta = join(tbl, metadata, 'Keys', 'SubjectID');

%% how to plot an example in 100307 in the left hemisphere

lH=tbl(tbl.Hemisphere=="left",:);

Ae_uncorrected=lH.SmoothPialArea(2:5);%the order of the lobes is: Frontal lobe, Parietal lobe, Temporal lobe, Occipital lobe
At_uncorrected=lH.PialArea(2:5);
T=lH.AvgCortThickness(2:5);
K=lH.GaussianCurvature(2:5);

correction_term=4*pi./K;

Ae_corrected=Ae_uncorrected.*correction_term;
At_corrected=At_uncorrected.*correction_term;

x=log10(Ae_corrected);
y=log10(At_corrected.*sqrt(T));
figure(1)
scatter(x,y)
lsline
[b,bint]=regress(y,[x ones(size(x))]);
disp(['Slope is between ' num2str(bint(1,1)) ' and ' num2str(bint(1,2))])

%for more examples of how the plots in the paper here 
%https://www.nature.com/articles/s42003-019-0421-7 were produced, see:
%https://doi.org/10.5281/zenodo.2595060
