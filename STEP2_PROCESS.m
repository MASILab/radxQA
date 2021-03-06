%% Phantom QA
% Monitor results of QA analysis - log SNR of QA scans at VUMC
% The purpose of this report is to aid in monitoring the results of QA
% scans performed at the Vanderbilt University Medical Center.  
%
% (1) Search in QA_dir/NIFTIS/ for files (.nii + .txt) 
% (2) If "ACR_T1"/"T1-axial" or "EPI_STABILITY", process
% (3) outputs saved with unique name in QA_dir/PROCESSED/
%
% Inputs:
%   none: hard-coded to QA directory on server
%
% Outputs:
%   (1) csv file for upload to REDCap 
%   (2) PNGs/JPGs of slice, with ROIs, plots
%   (3) .sh with python code to execute
%
%
% T1 measures:
%
%
% EPI Stability measures:
%
%
% Packages required: 
% Functions required: 
% Files required: 
%
% See also: 
%
% Author:  Allen Newton (allen.t.newton@vanderbilt.edu)
% Date:    XXX
% Version: 1.0
% Changelog:
% v2: 18-Apr-2018 - schillkg: new functions, new measures
% 31-Jan-2018  - schillkg
%       modified from VCH_QA_analysis.m
%
%------------- BEGIN CODE --------------

%% Preliminaries

QA_dir = ['/blazepool' filesep 'radx' filesep 'QA'];
NIFTI_dir = [QA_dir filesep 'NIFTIS'];
OUT_dir = [QA_dir filesep 'PROCESSED'];
addpath(genpath([QA_dir filesep 'radxQA']));

%%%%%%%%%%%

if exist([OUT_dir filesep 'running.txt'], 'file') == 2
    disp('another instance running')
    exit
else
    disp('program will run...')
    try
        % create text file running.txt
        fileID = fopen([OUT_dir filesep 'running.txt'],'w');
        fprintf(fileID,'this is running');
        fclose(fileID);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % code goes here



if exist([OUT_dir filesep 'Results_ALL.mat'], 'file') == 2
    load([OUT_dir filesep 'Results_ALL.mat']);
    s_all_index = size(s_all,2)+1;
else
    s_all = struct('name',[],'type',[],'scanner',[],'scan_date',[],'series_description',[], ...
        'reso',[],'GeometricAccuracy',[],'GeometricDistortion',[],'HCSpatialRes',[],'SliceThick',[],'SlicePositionS1',[], ...
        'SlicePositionS11',[],'PIU',[],'PSG',[],'LowContrastDet',[],'SNR',[],'SFNR',[], ...
        'PF',[],'Drift',[],'RDC',[]);
    s_all_index = 1;
end

%% ACR_T1
% search through text files, with 'ACR_T1' or 'T1-axial' in name
% *** do not want ACR T1 CLEAR ***
% loop through all files

files1 = dir(fullfile(NIFTI_dir,'*ACR_T1*.nii'));
fileNames1 = {files1.name};
files2 = dir(fullfile(NIFTI_dir,'*T1-axial.nii'));
fileNames2 = {files2.name};
fileNames = [fileNames1 fileNames2];

for ii = 1:length(fileNames)
    
    try
        
    filename = fileNames{ii};
    [~,name,ext] = fileparts(filename);
    disp(name)
    
     % if name exists in Results, continue to next
    if sum(strcmp({s_all.name}, name))~=0
        disp('RESULTS ALL Already Contains this info..........')
        continue
    end
    
%     % if QA_dir/name.csv exists, continue to next file
%     if exist([OUT_dir filesep name '.csv'], 'file') == 2
%         disp('CSV output already exists...')
%         continue
%     end
    
    outfolder = [OUT_dir filesep name];
    mkdir(outfolder);
    
    % load nii
    nii = load_nii([NIFTI_dir filesep name '.nii']);
    % load text file, grab scanner, scan_date, series_description
    filetext = fileread([NIFTI_dir filesep name '_info.txt']);
    
    expr = '[^\n]*Subject:[^\n]*';
    matches = regexp(filetext,expr,'match');
    scanner = strrep(matches{1},'Subject: ','');
    
    expr = '[^\n]*Series date:[^\n]*';
    datematches = regexp(filetext,expr,'match');
    scan_date = strrep(datematches{1},'Series date: ','');
    
    expr = '[^\n]*Series description:[^\n]*';
    series_matches = regexp(filetext,expr,'match');
    series_description = strrep(series_matches{1},'Series description: ','');
    
    % get image volume
    img=double(squeeze(nii.img(:,:,:,1,1)));
    % get pixel size
    reso = nii.hdr.dime.pixdim(2:4);
    
    % (1) Geometric Accuracy
    % Slice 1 diameter measurement (*2)
    % Slice 5 diameter measurement (*4)
    [s1d1,s1d2] = GeometricAccuracy_S1(img,reso); pause(.7);
    export_fig([OUT_dir filesep name filesep 'GeometricAccuracy', '.png'],'-png'); close all;  pause(.7);

    % 2 measures + image
    % 188-192 acceptance
    %[s5d1,s5d2,s5d3,s5d4] = GeometricAccuracy_S5(img,reso);
    
    % (1b) global geometric distortion
    [A,B,C,D,Distortion] = GlobalDistortion_S7(img,reso);  pause(.7);
    % 4 diameters and image
    export_fig([OUT_dir filesep name filesep 'GeometricDistortion', '.png'],'-png'); close all;  pause(.7);
    
    % (2) High Contrast Spatial Resolution
    [HCSR,pf_ggd] = HighContrastSpatialRes(img,reso);  
    % HCSR is .9, 1, 1.1 and pass fail for row/col respectively
    
    % (3) Slice Thickness
    [slc_thk,pf_slthick] = SliceThicknessAccuracy(img,reso);  
    % slc_thk is thickness in mm, pf_slthick between 4.3 and 5.7
    
    % (4) Slice Position Accuracy - Slice 1
    [l_diff,sl_l_diff,pf_hdl] = SlicePositionAccuracyS1(img,reso);  
    % length diff, slice displacement, pass/fail
    
    [l_diff11,sl_l_diff11,pf_hdl11] = SlicePositionAccuracyS11(img,reso);  
    % length diff, slice displacement, pass/fail
    
    % (5) Image Intensity Uniformity
    [PIU,PIU2]=PercentIntegralUniformity(img,reso);
    % >87.5 for 1.5, >82 for 3T
    
    % (6) Percent signal ghosting
    [PSG,pf_psg]=PercentSignalGhosting(img,reso);  pause(.7);
    % signal ghost, pass/fail (.025), and image
    export_fig([OUT_dir filesep name filesep 'Ghosting', '.png'],'-png'); close all;  pause(.7);

    % (7) low contrast detectability
    % slice 8-11
    slice = 11;
    [truth11,spokes11] = LowContrastDetectability(img,reso,slice);  pause(.7);
    export_fig([OUT_dir filesep name filesep 'ContrastBest', '.png'],'-png'); close all;  pause(.7);
     slice = 10;
    [truth10,spokes10] = LowContrastDetectability(img,reso,slice); close all;
    slice = 9;
    [truth9,spokes9] = LowContrastDetectability(img,reso,slice); close all;
     slice = 8;
    [truth8,spokes8] = LowContrastDetectability(img,reso,slice);  pause(.7);
    export_fig([OUT_dir filesep name filesep 'ContrastWorse', '.png'],'-png'); close all;  pause(.7);
    spokes_total = spokes8+spokes9+spokes10+spokes11;
    
    [SNR] = SignalToNoiseRatioPhantom(img,reso); pause(.7);
    export_fig([OUT_dir filesep name filesep 'SNR', '.png'],'-png'); close all; pause(.7);

    % save outputs in own matrix (.mat)
    s = struct;
    s.name = name;
    s.type = 'T1';
    s.scanner = scanner;
    s.scan_date = scan_date;
    s.series_description = series_description;
    s.reso = reso;
    s.GeometricAccuracy = [s1d1,s1d2];
    s.GeometricDistortion = [A,B,C,D,Distortion];
    s.HCSpatialRes = HCSR;
    s.SliceThick = [slc_thk,pf_slthick];
    s.SlicePositionS1 = [l_diff,sl_l_diff,pf_hdl];
    s.SlicePositionS11 = [l_diff11,sl_l_diff11,pf_hdl11];
    s.PIU = [PIU,PIU2];
    s.PSG = [PSG,pf_psg];
    s.LowContrastDet = spokes_total;
    s.SNR = SNR;
    s.SFNR = [];
    s.PF = [];
    s.Drift = [];
    s.RDC = [];
    
    save([outfolder filesep 'Results.mat'],'s')
    pause(.7);
    
    % save outputs in overall matrix (.mat)
    s_all(s_all_index) = s;
    s_all_index = s_all_index + 1;
    
    % save outputs in text file
    % create sh with python code    
    
%     % save outputs in QA_dir/PROCESSED (OUT_dir)
%     matrix1 = {'record_id';'scan_date';'scanner';'series_description';'series_snr';'phantom_data_complete'};
%     matrix2 = {'1';[num2str(scan_date(1:4)) '-' num2str(scan_date(5:6)) '-' num2str(scan_date(7:8))];scanner;series_description;num2str(SNR);'1'};
%     fid = fopen([OUT_dir filesep name '.csv'], 'w' );
%     fprintf(fid,'%s,%s,%s,%s,%s,%s\n%s,%s,%s,%s,%s,%s',matrix1{:},matrix2{:});
%     fclose(fid);
%     
%     % create sh with python code to execute, place in OUT_dir
%     % example: python phantom2redcap.py csv_path image_path
%     unixtr = sprintf('#!/bin/bash\npython %s/phantom2redcap.py %s %s',[QA_dir filesep 'EXTRA'],[OUT_dir filesep name '.csv'],[OUT_dir filesep name, '.png']);
%     fid = fopen([OUT_dir filesep name '.sh'],'wt');
%     fprintf(fid,'%s',unixtr);
%     fclose(fid);
       
    clear filename name ext nii filetext img data s filename name ext outfolder
    
    save([OUT_dir filesep 'Results_ALL.mat'],'s_all');    
    pause(.7);
    
    catch
        disp('failed for some reason...')
    end
    
end

%% 

files = dir(fullfile(NIFTI_dir,'*EPI_STABILITY*.nii'));
fileNames = {files.name};

for ii = 1:length(fileNames)
    
    try
    filename = fileNames{ii};
    [~,name,ext] = fileparts(filename);
    disp(name)
    
    % if name exists in Results, continue to next
    %load([OUT_dir filesep 'Results_ALL.mat']);
    if sum(strcmp({s_all.name}, name))~=0
        disp('RESULTS ALL Already Contains this info..........')
        continue
    end
    
%     % if QA_dir/name.csv exists, continue to next file
%     if exist([OUT_dir filesep name '.csv'], 'file') == 2
%         disp('CSV output already exists...')
%         continue
%     end
    
    outfolder = [OUT_dir filesep name];
    mkdir(outfolder);
    
    % load nii
    nii = load_nii([NIFTI_dir filesep name '.nii']);
    % load text file, grab scanner, scan_date, series_description
    filetext = fileread([NIFTI_dir filesep name '_info.txt']);
    
    expr = '[^\n]*Subject:[^\n]*';
    matches = regexp(filetext,expr,'match');
    scanner = strrep(matches{1},'Subject: ','');
    
    expr = '[^\n]*Series date:[^\n]*';
    datematches = regexp(filetext,expr,'match');
    scan_date = strrep(datematches{1},'Series date: ','');
    
    expr = '[^\n]*Series description:[^\n]*';
    series_matches = regexp(filetext,expr,'match');
    series_description = strrep(series_matches{1},'Series description: ','');
    
     % get image volume
    img=double(squeeze(nii.img(:,:,:,:,1)));
    % get pixel size
    reso = nii.hdr.dime.pixdim(2:5);
    
    % (1) geometric accuracy
    [A,B,C,D,Distortion,slice] = GlobalDistortion_EPI(img,reso); pause(.7);
    export_fig([OUT_dir filesep name filesep 'GeometricDistortion', '.png'],'-png'); close all; pause(.7);
    
    % (2) Image Intensity Uniformity
    [PIU,PIU2]=PercentIntegralUniformity_EPI(img,reso,slice);
    
    % (3) Percent Signal Ghosting
    [PSG,pf_psg]=PercentSignalGhosting_EPI(img,reso,slice); pause(.7);
    export_fig([OUT_dir filesep name filesep 'Ghosting', '.png'],'-png'); close all; pause(.7);
    
    % 4)Signal to Fluctionation Noise Ratio Summary Value
    [SFNR,SNR,PF,Drift,RDC]=FriedmanMeasures(img,reso,slice);
    
    % save outputs in own matrix (.mat)
    s = struct;
    s.name = name;
    s.type = 'EPI STABILITY';
    s.scanner = scanner;
    s.scan_date = scan_date;
    s.series_description = series_description;
    s.reso = reso;
    s.GeometricAccuracy = [];
    s.GeometricDistortion = [A,B,C,D,Distortion,slice];
    s.HCSpatialRes = [];
    s.SliceThick = [];
    s.SlicePositionS1 = [];
    s.SlicePositionS11 = [];
    s.PIU = [PIU,PIU2];
    s.PSG = [PSG,pf_psg];
    s.LowContrastDet = [];
    s.SNR = SNR;
    s.SFNR = SFNR;
    s.PF = PF;
    s.Drift = Drift;
    s.RDC = RDC;
    
    save([outfolder filesep 'Results.mat'],'s')
    pause(.7);
    
    s_all(s_all_index) = s;
    s_all_index = s_all_index + 1;
    
%     % save outputs in QA_dir/PROCESSED (OUT_dir)
%     matrix1 = {'record_id';'scan_date';'scanner';'series_description';'series_snr';'phantom_data_complete'};
%     matrix2 = {'1';[num2str(scan_date(1:4)) '-' num2str(scan_date(5:6)) '-' num2str(scan_date(7:8))];scanner;series_description;num2str(SNR);'1'};
%     fid = fopen([OUT_dir filesep name '.csv'], 'w' );
%     fprintf(fid,'%s,%s,%s,%s,%s,%s\n%s,%s,%s,%s,%s,%s',matrix1{:},matrix2{:});
%     fclose(fid);
%     
%     % create sh with python code to execute, place in OUT_dir
%     % example: python phantom2redcap.py csv_path image_path
%     unixtr = sprintf('#!/bin/bash\npython %s/phantom2redcap.py %s %s',[QA_dir filesep 'EXTRA'],[OUT_dir filesep name '.csv'],[OUT_dir filesep name, '.png']);
%     fid = fopen([OUT_dir filesep name '.sh'],'wt');
%     fprintf(fid,'%s',unixtr);
%     fclose(fid);
    
    clear filename name ext nii filetext img data s filename name ext outfolder

    pause(.7);
    save([OUT_dir filesep 'Results_ALL.mat'],'s_all');
    catch
        disp('failed for some reason....')
    end
    
end


 %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % delete running.txt file
        delete([OUT_dir filesep 'running.txt'])
    catch
        % if failed, still delete file
        % delete running.txt file
        delete([OUT_dir filesep 'running.txt'])
    end
end % if running.txt exists
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





