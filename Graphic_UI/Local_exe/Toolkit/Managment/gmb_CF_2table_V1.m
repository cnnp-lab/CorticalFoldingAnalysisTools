function [tbl,report] = gmb_CF_2table_V1(DT_path,Ses,report,CF_MD,Hemi_MD,atlas,SCL)

% Check if there is a Session
ses_cnt = strcmp(Ses,'_*_');
if ~ses_cnt
    DT_path = fullfile(DT_path,Ses);
end

tbl =[];
corrupt = [];

% Check if no-Scaling is requested
SCL_MD = dec2bin(SCL,3);
if SCL_MD(3) == '1'
    switch CF_MD
        case 1 % Lobewise analysis
            [a, b] = extract_FreeSurferLobes_features_Vgmb(DT_path,...
                'hemi', Hemi_MD,'verbose',false,'atlas',atlas);

            tbl =[tbl;a];
            corrupt = [corrupt,b];
            CF = 'Lobe';

        case 2 % Hemispherewise analysis
            [a, b] = extract_FreeSurferHemi_features_Vgmb(DT_path,...
                'hemi', Hemi_MD,'verbose',false);

            tbl =[tbl;a];
            corrupt = [corrupt,b];
            CF = 'Hemipshere';
    end
end

% Cjhcc is scaling is requested
switch SCL_MD([1,2])
    case '10'%High Scales
        targetScale=0.3250:0.025:.375;
        nIter=[4 4 4];
    case '01'%Low Scales
        targetScale=0.225:0.025:0.3000;
        nIter=[5 5 4 4 ];
    case '11'%All Sacles
        targetScale= 0.225:0.025:.375;
        nIter=[5 5 4 4 4 4 4];
    case '00'% No Scaling
        targetScale = [];
        nIter = [];
end

if ~isempty(nIter)
    switch CF_MD
        case 1 % Lobewise analysis
            % Not implemented yet
            
            CF = 'Lobe';

        case 2 % Hemispherewise analysis
            % Not implemented yet
            [a, b] = fastEstimateScale(targetScale,nIter,DT_path,Hemi_MD,0);
            
            tbl =[tbl;a];
            corrupt = [corrupt,b];
            CF = 'Hemipshere';
    end
end

if sum(corrupt)==0
    % Message if all goes well
    aux = '      ESTIMATED => ';
    if ~ses_cnt
        aux = [aux,Ses,'/'];
    end
    aux = [aux,CF,' complete'];
    report = [report; string(aux)];

else
    % For the report message
    switch Hemi_MD
        case {"left","right"}
            Lobe_cod = Hemi_MD;
        otherwise
            Lobe_cod ={"left","right"};
    end

    % Error mesages if corrupted folders & Remove the nonvalid data
    for l=1:length(corrupt)
        if corrupt(l)
            % Report the issue
            aux = '      WARNING => ';
            if ~ses_cnt
                aux = [aux,Ses,'/'];
            end
            aux = [aux,char(Lobe_cod{l}),' corrupted during ',CF,' estimation'];
            report = [report; string(aux)];

            switch Hemi_MD
                case {"avg","sum"}
                    ldx = strcmp(tbl.Hemisphere,Hemi_MD);

                otherwise
                    ldx = strcmp(tbl.Hemisphere,Lobe_cod{l});
            end
            tbl(ldx,:) = [];
        end
    end
end

