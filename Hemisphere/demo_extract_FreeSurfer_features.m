%% DEMO FILE
% We show how to load controls and patients into separate tables and then
% merge both tables and join them with a corresponding metadata file.
% Note that in general there is no need to do it separately, but here we do
% to demonstrate the merging process.

%The lib folder and subfolders need to be added to path to work. Note that
%we included the freesurfer matlab functions again in the lib for
%completeness. Please refer to the freesurfer manual for further
%information on this.

addpath(genpath('../lib'))

%% load controls
% subjdir contains the FreeSurfer files, one directory per subject
subjdir='/bigstorage/Beth/TLE/';

% IDs are in the form 'C001', 'C002', ...
% you have to select which subjects you want to calculate by string IDs!
ids = num2str([1:30]', '%03u'); % formating ID (leading zeros)
ids = string([repmat('C', [length(ids),1]), ids]);

% pattern for foldername of each subject
format='*_T1.nii.gz'; % '*' will be replaced by the ID

% How to handle hemispheres? 'avg' for averaging, 'sum' for summing
% (caution: might make sense for area/volume but not for e.g. thickness),
% 'left' or 'right' for including only one of them respectively,
% or 'both' for having separate rows for left and right hemispheres.
% The column 'Hemisphere' indicates 'avg', 'sum' or 'left'/'right'.
hemi = "both";

% run extraction - this is the main function doing the data extraction!
controls = extract_FreeSurferHemi_features(subjdir, ids, format, 'hemi', hemi);

%% load patients and merge with controls
ids = ["816","844","871","873","877","878","887","902","909","946","968", ...
    "977","983","985","986","996","997","1001","1005","1012","1015","1022", ...
    "1035","1038","1041","1052","1055","1059","1065","1068","1082","1083", ...
    "1102","1118","1120","1131","1136","1144","1147","1148","1151","1152", ...
    "1179","1200","1210","1236","1252","1258","1260","1285","1295","1347","1359"];

% merging is done by specifying the old table to add to/replacing old by new ids
tle_controls = extract_FreeSurferHemi_features(subjdir, ids, format, ...
                                           'hemi', hemi, ...
                                           'saveas', 'tle_controls', ...
                                           'oldtbl', controls); % merge!
                                       
%% join with metadata

% load metadata table
load('./demo/metadata.mat'); % will load variable 'metadata'

% join table with metadata table by SubjectID (common column in the tables)
tle_controls = join(tle_controls, metadata, 'Keys', 'SubjectID');