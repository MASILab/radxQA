%% Phantom QA
% Monitor results of QA analysis - log SNR of QA scans at VUMC
% The purpose of this report is to aid in monitoring the results of QA
% scans performed at the Vanderbilt University Medical Center.  
%
% (1) Search in QA_dir/NIFTIS/ for files (.nii + .txt) 
% (2) If "ACR_T1" or "EPI_STABILITY", process
% (3) outputs saved with unique name in QA_dir/PROCESSED/
%
% Inputs:
%   none: hard-coded to QA directory on server
%
% Outputs:
%   (1) csv file for upload to REDCap (SNRs, scan info)
%   (2) PNG/JPG of slice, with ROIs drawn
%   (3) .sh with python code to execute
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
% 31-Jan-2018  - schillkg
%       modified from VCH_QA_analysis.m
%
%------------- BEGIN CODE --------------

%% Preliminaries

QA_dir = ['/Volumes' filesep 'GRAID' filesep 'BlazeFake' filesep 'QA'];
NIFTI_dir = [QA_dir filesep 'NIFTIS'];
OUT_dir = [QA_dir filesep 'PROCESSED'];
addpath(genpath([QA_dir filesep 'EXTRA']));

%% ACR_T1
% search through text files, with 'ACR_T1' or 'T1-axial' in name
% *** do not want ACR T1 CLEAR ***
% loop through all files
% process each, calculate spatial SNR
% output date of scan, scanner, scan_name, SNR, image

files = dir(fullfile(NIFTI_dir,'*ACR_T1*.nii'));
fileNames = {files.name};

for ii = 1:length(fileNames)
    filename = fileNames{ii};
    [~,name,ext] = fileparts(filename);
    disp(name)
    % if QA_dir/name.csv exists, continue to next file
    if exist([OUT_dir filesep name '.csv'], 'file') == 2
        disp('CSV output already exists...')
        continue
    end
    
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
    
    % calculations
    img=double(squeeze(nii.img(:,:,:,1,1)));
    
    kernelsize=5;
    edges=img-convn(img,ones(kernelsize,kernelsize,1),'same')./(kernelsize^2);
    [~,index]=sort(squeeze(sum(sum(abs(edges),1),2)),'ascend');
    SliceToAnalyze=index(1);
    
    data=squeeze(img(:,:,SliceToAnalyze,1,1));
    wbmask=mywbmasker2fcn(data);
    obj=zeros(size(wbmask));
    centerpt=round(findcentroid(data));
    obj=IMcircle(size(obj),centerpt,40);
    bgd=~bwmorph(wbmask,'dilate',20);
    obj=obj>0;
    bgd=bgd>0;
    SNR=mean(data(obj))/std(data(bgd));
    
    % make figure
    h=figure;
    data(bgd)=mean(data(:)); % this sets the background to a non-zeros number for display purposes;
    mymontage2(squeeze(data(:,:,round(size(data,3)/2),1,1)),'io',squeeze(obj(:,:,round(size(data,3)/2),1,1)),'c',squeeze(bgd(:,:,round(size(data,3)/2),1,1)),'r');
    text(10,10,'std.dev.','Color',[1 0.5 0.5]);
    text(round(mean(xlim)),round(mean(ylim)),'mean','Color','c');
    set(gca,'XColor','k');
    export_fig([OUT_dir filesep name, '.png'],'-png');

    % save outputs in QA_dir/PROCESSED (OUT_dir)
    matrix1 = {'record_id';'scan_date';'scanner';'series_description';'series_snr';'phantom_data_complete'};
    matrix2 = {'1';[num2str(scan_date(1:4)) '-' num2str(scan_date(5:6)) '-' num2str(scan_date(7:8))];scanner;series_description;num2str(SNR);'1'};
    fid = fopen([OUT_dir filesep name '.csv'], 'w' );
    fprintf(fid,'%s,%s,%s,%s,%s,%s\n%s,%s,%s,%s,%s,%s',matrix1{:},matrix2{:});
    fclose(fid);
    
    % create sh with python code to execute, place in OUT_dir
    % example: python phantom2redcap.py csv_path image_path
    unixtr = sprintf('#!/bin/bash\npython %s/phantom2redcap.py %s %s',[QA_dir filesep 'EXTRA'],[OUT_dir filesep name '.csv'],[OUT_dir filesep name, '.png']);
    fid = fopen([OUT_dir filesep name '.sh'],'wt');
    fprintf(fid,'%s',unixtr);
    fclose(fid);
       
    clear filename name ext nii filetext img data
end

%% EPI_STABILITY
% search through text files, with 'ACR_T1' or 'T1-axial' in name
% *** do not want ACR T1 CLEAR ***
% loop through all files
% process each, calculate spatial SNR
% output date of scan, scanner, scan_name, SNR, image

files = dir(fullfile(NIFTI_dir,'*EPI_STABILITY*.nii'));
fileNames = {files.name};

for ii = 1:length(fileNames)
    filename = fileNames{ii};
    [~,name,ext] = fileparts(filename);
    disp(name)
    % if QA_dir/name.csv exists, continue to next file
    if exist([OUT_dir filesep name '.csv'], 'file') == 2
        disp('CSV output already exists...')
        continue
    end
    
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
    
    % calculations
    img = double(squeeze(nii.img));
    tSNR=mean(img,4)./std(img,[],4);
    tSNR(tSNR==Inf) = NaN;
    
    IM=mean(img,4);
    kernelsize=5;
    edges=IM-convn(IM,ones(kernelsize,kernelsize,1),'same')./(kernelsize^2);
    [~,index]=sort(squeeze(sum(sum(abs(edges),1),2)),'ascend');
    SliceToAnalyze=index(1);
    
    data=squeeze(img(:,:,SliceToAnalyze,1));
    tSNR=squeeze(tSNR(:,:,SliceToAnalyze));
    wbmask=mywbmasker2fcn(data);
    obj=zeros(size(wbmask));
    centerpt=round(findcentroid(data));
    obj=IMcircle(size(data),centerpt,20);
    stability=mean(tSNR(obj));
    
    % make figure
    h=figure;
    mymontage2(data,'io',obj,'c');
    text(10,10,'std.dev.','Color',[1 0.5 0.5]);
    text(round(mean(xlim)),round(mean(ylim)),'mean','Color','c');
    axis on
    set(gca,'XColor','k');
    set(gcf,'position',[1 1 600 600]);
    export_fig([OUT_dir filesep name, '.png'],'-png');
    
    % save outputs in QA_dir/PROCESSED (OUT_dir)
    matrix1 = {'record_id';'scan_date';'scanner';'series_description';'series_snr';'phantom_data_complete'};
    matrix2 = {'1';[num2str(scan_date(1:4)) '-' num2str(scan_date(5:6)) '-' num2str(scan_date(7:8))];scanner;series_description;num2str(SNR);'1'};
    fid = fopen([OUT_dir filesep name '.csv'], 'w' );
    fprintf(fid,'%s,%s,%s,%s,%s,%s\n%s,%s,%s,%s,%s,%s',matrix1{:},matrix2{:});
    fclose(fid);
    
    % create sh with python code to execute, place in OUT_dir
    % example: python phantom2redcap.py csv_path image_path
    unixtr = sprintf('#!/bin/bash\npython %s/phantom2redcap.py %s %s',[QA_dir filesep 'EXTRA'],[OUT_dir filesep name '.csv'],[OUT_dir filesep name, '.png']);
    fid = fopen([OUT_dir filesep name '.sh'],'wt');
    fprintf(fid,'%s',unixtr);
    fclose(fid);
    
    clear filename name ext nii filetext img
end

    
    
    
    
    
    
    


