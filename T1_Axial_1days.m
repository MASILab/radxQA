function []=T1_Axial_1days(s_all,Dates,RESULTS_dir,OUT_dir)

% grab all with s_all.type = 'T1'
T1_binary = strcmp({s_all.type}, 'T1')

% T1: get most recent date, 100 days, all time dates
T1MostRecent = max(Dates(T1_binary==1));
T1MostRecentIndex = find(Dates==T1MostRecent & T1_binary');

clear GA Dist t1 t2 HCSR SliceThick SlicePos PIU PSG LCOD SNR

name = s_all(T1MostRecentIndex).name;

load([OUT_dir filesep name filesep 'Results.mat']) % load S

A = imread([OUT_dir filesep name filesep 'SNR.png']);
B = imread([OUT_dir filesep name filesep 'ContrastWorse.png']);
C = imread([OUT_dir filesep name filesep 'ContrastBest.png']);
D = imread([OUT_dir filesep name filesep 'Ghosting.png']);
E = imread([OUT_dir filesep name filesep 'GeometricDistortion.png']);
F = imread([OUT_dir filesep name filesep 'GeometricAccuracy.png']);

if (sum(s.GeometricAccuracy<188) + sum(s.GeometricAccuracy>192)) > 0; pfg = 'FAIL'; else pfg = 'PASS'; end   
if (max(s.HCSpatialRes)) > 1; pfhc = 'FAIL'; else pfhc = 'PASS'; end   
if s.SliceThick(2) ==0; pfst = 'FAIL'; else pfst = 'PASS'; end      
if s.SlicePositionS1(3) ==0; pfsp = 'FAIL'; else pfsp = 'PASS'; end      
if s.PIU(1) < .82; pfpiu = 'FAIL'; else pfpiu = 'PASS'; end      
if s.PSG(2) ==0; pfpg = 'FAIL'; else pfpg = 'PASS'; end      
if s.LowContrastDet < 37; pflc = 'FAIL'; else pflc = 'PASS'; end      
[aa bb] =max(abs(s.GeometricAccuracy-190)); 

data = {'Scan Name',s.name,'N/A'; ...
    'Scan Type',s.type ,'N/A'; ...
    'Scanner',s.scanner, 'N/A'; ...
    'Date',s.scan_date, 'N/A'; ...
    'Description',s.series_description, 'N/A'; ...
    'Resolution (xy)',num2str(s.reso(1)), 'N/A'; ...
    'Geometric Accuracy',num2str(s.GeometricAccuracy(bb)), pfg; ...
    'Geometric Distrotion',num2str(s.GeometricDistortion(5)), 'N/A'; ...
    'Spatial Res',num2str(max(s.HCSpatialRes)), pfhc; ...
    'Slice Thickness',num2str(s.SliceThick(1)), pfst; ...
    'Slice Position',num2str(s.SlicePositionS1(2)), pfsp; ...
    '% Uniformity',num2str(s.PIU(1)), pfpiu; ...
    '% Signal Ghosting',num2str(s.PSG(1)), pfpg; ...
    '# spokes detected',num2str(s.LowContrastDet), pflc; ...
    'SNR',num2str(s.SNR), 'N/A'}   ;   

close all;
f = figure(5); set(gcf, 'units','normalized','outerposition',[0.05 0.05 .9 .9]); 
t = uitable(f,'Data',data,'ColumnWidth',{100,300,70})
subplot(1,2,2),plot(3)
pos = get(subplot(1,2,2),'position');
delete(subplot(2,1,2))
set(t,'units','normalized')
set(t,'position',pos)
axis off; 
subplot(3,4,1); imagesc(A); axis equal; axis off; title('SNR ROIs');
subplot(3,4,2); imagesc(D); axis equal; axis off; title('Ghosting ROIs');
subplot(3,4,5); imagesc(B); axis equal; axis off; colormap gray; title('Worst Case Contrast');
subplot(3,4,6); imagesc(C); axis equal; axis off; colormap gray; title('Best Case Contrast');
subplot(3,4,9); imagesc(E); axis equal; axis off; title('Distortion Measures');
subplot(3,4,10); imagesc(F); axis equal; axis off; title('Geometric Accuracy');

saveas(gcf,[RESULTS_dir filesep 'T1_SingleDay.fig'])
print_current_figure(400,[RESULTS_dir filesep 'T1_SingleDay.png'])
pause(2); close all;

%%

fileID = fopen([RESULTS_dir filesep 'T1_SingleDay.txt'],'w');
fprintf(fileID,['Scanner Type: ' s.type ' \n']);
fprintf(fileID,['Scanner: ' s.scanner ' \n']);
fprintf(fileID,['Date: ' s.scan_date ' \n']);
fprintf(fileID,['Description: ' s.series_description ' \n']);
fprintf(fileID,['Resolution (xy): ' num2str(s.reso(1)) ' \n']);
fprintf(fileID,['Geometric Accuracy: ' pfg ' : ' num2str(s.GeometricAccuracy(bb)) ' \n']);
fprintf(fileID,['Geometric Distortion: ' num2str(s.GeometricDistortion(5)) ' \n']);
fprintf(fileID,['Spatial Res: ' pfhc ' : ' num2str(max(s.HCSpatialRes)) ' \n']);
fprintf(fileID,['Slice Thickness: ' pfst ' : ' num2str(s.SliceThick(1)) ' \n']);
fprintf(fileID,['Slice Position: ' pfsp ' : ' num2str(s.SlicePositionS1(2)) ' \n']);
fprintf(fileID,['Percent Uniformity: ' pfpiu ' : ' num2str(s.PIU(1)) ' \n']);
fprintf(fileID,['Percent Signal Ghosting: ' pfpg ' : ' num2str(s.PSG(1)) ' \n']);
fprintf(fileID,['Number spokes detected: ' pflc ' : ' num2str(s.LowContrastDet) ' \n']);
fprintf(fileID,['SNR: ' num2str(s.SNR) ' \n']);



fclose(fileID);





