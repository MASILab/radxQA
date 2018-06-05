%% Step5: Post-process Results
% grab processed data, make plots for
% (1) most recent, single time point - T1, EPI Stability
% (2) last week results - T1, EPI Stability
% (3) last month results - T1, EPI Stability
% (4) last 100 days results - T1, EPI Stability
% (5) last forever results - T1, EPI Stability
% make plots and text for error values
% place in Var/www/html/QA/

%% Preliminaries

QA_dir = ['/Volumes' filesep 'GRAID' filesep 'BlazeFake' filesep 'QA'];
NIFTI_dir = [QA_dir filesep 'NIFTIS'];
OUT_dir = [QA_dir filesep 'PROCESSED'];
RESULTS_dir = [QA_dir filesep 'Var' filesep 'www' filesep 'html' filesep 'QA'];
addpath(genpath([QA_dir filesep 'EXTRA']));

% load s_all
load([OUT_dir filesep 'Results_ALL.mat']);

% get all dates
Dates = zeros(size(s_all,2),1);
for i=1:size(s_all,2)
    Dates(i) = datenum(datetime(s_all(i).scan_date,'InputFormat','yyyyMMdd'));
    %Dates(i) = datenum(s_all(i).scan_date,'yyyyMMdd'); % this doesn't handle year change well?
end

%% T1 axial
%T1_Axial_100days(s_all,Dates,RESULTS_dir);
%T1_Axial_AllTime(s_all,RESULTS_dir);
%T1_Axial_30days(s_all,Dates,RESULTS_dir);
%T1_Axial_7days(s_all,Dates,RESULTS_dir);
T1_Axial_Eval_v2(s_all,Dates,RESULTS_dir);
T1_Axial_1days_v2(s_all,Dates,RESULTS_dir,OUT_dir);

%% EPI

EPI_Eval_v2(s_all,Dates,RESULTS_dir);
EPI_1days_v2(s_all,Dates,RESULTS_dir,OUT_dir);





