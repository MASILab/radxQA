function [] = EPI_1days(s_all,Dates,RESULTS_dir,OUT_dir)


EPI_binary = strcmp({s_all.type}, 'EPI STABILITY');
EPIIndex = find(EPI_binary);

clear GA Dist t1 t2 HCSR SliceThick SlicePos PIU PSG LCOD SNR
clear GA Dist t1 t2 PIU PSG SNR SFNR PF Drift RDC

ind=1;
for i=1:length(EPIIndex)
    index = EPIIndex(i);
    %if sum(strfind(s_all(index).name, 'ACR'))==0;
    GA(ind) = mean(s_all(index).GeometricDistortion(1:4));
    Dist(ind) = s_all(index).GeometricDistortion(5);
    t1(ind) = datetime(s_all(index).scan_date,'InputFormat','yyyyMMdd');
    t2(ind) = datenum(datetime(s_all(index).scan_date,'InputFormat','yyyyMMdd'));
    PIU(ind) = s_all(index).PIU(1);
    PSG(ind) = s_all(index).PSG(1);
    SNR(ind) = s_all(index).SNR;
    SFNR(ind) = s_all(index).SFNR;
    PF(ind) = s_all(index).PF;
    Drift(ind) = s_all(index).Drift;
    RDC(ind) = s_all(index).RDC;
    ind=ind+1;
    %end
end

[t1_sorted,t1_order]=sort(t1);

% get most recent date, 100 days, all time dates
EPIMostRecent = max(Dates(EPI_binary==1));
EPIMostRecentIndex = find(Dates==EPIMostRecent & EPI_binary');

name = s_all(EPIMostRecentIndex).name;
clear s
load([OUT_dir filesep name filesep 'Results.mat']) % load S

D = imread([OUT_dir filesep name filesep 'Ghosting.png']);
E = imread([OUT_dir filesep name filesep 'GeometricDistortion.png']);

clear data;
data = {'Scan Name',s.name,[],''; ...
    'Scan Type',s.type ,[],''; ...
    'Scanner',s.scanner,[], ''; ...
    'Date',s.scan_date, [],''; ...
    'Description',s.series_description,[], ''; ...
    'Resolution (xy)',num2str(s.reso(1)),[], '';}   ;   

Param = GA; rownum = 7; name = 'Geometric Accuracy'
Val_sorted = Param(t1_order);
[Shewart,message] = ShewartEval(Val_sorted);
data{rownum,1} = name; data{rownum,2}=Param(end); data{rownum,3}=Shewart; data{rownum,4}=message;

Param = Dist; rownum = 8; name = 'Distortion'
Val_sorted = Param(t1_order);
[Shewart,message] = ShewartEval(Val_sorted);
data{rownum,1} = name; data{rownum,2}=Param(end); data{rownum,3}=Shewart; data{rownum,4}=message;

Param = PIU; rownum = 9; name = 'Image Intensity Uniformity'
Val_sorted = Param(t1_order);
[Shewart,message] = ShewartEval(Val_sorted);
data{rownum,1} = name; data{rownum,2}=Param(end); data{rownum,3}=Shewart; data{rownum,4}=message;

Param = PSG; rownum = 10; name = 'Ghosting Ratio'
Val_sorted = Param(t1_order);
[Shewart,message] = ShewartEval(Val_sorted);
data{rownum,1} = name; data{rownum,2}=Param(end); data{rownum,3}=Shewart; data{rownum,4}=message;

Param = SNR; rownum = 11; name = 'Signal To Noise Ratio'
Val_sorted = Param(t1_order);
[Shewart,message] = ShewartEval(Val_sorted);
data{rownum,1} = name; data{rownum,2}=Param(end); data{rownum,3}=Shewart; data{rownum,4}=message;

Param = SFNR; rownum = 12; name = 'Signal To Fluctuation Noise Ratio'
Val_sorted = Param(t1_order);
[Shewart,message] = ShewartEval(Val_sorted);
data{rownum,1} = name; data{rownum,2}=Param(end); data{rownum,3}=Shewart; data{rownum,4}=message;

Param = PF; rownum = 13; name = 'Percent Fluctuation'
Val_sorted = Param(t1_order);
[Shewart,message] = ShewartEval(Val_sorted);
data{rownum,1} = name; data{rownum,2}=Param(end); data{rownum,3}=Shewart; data{rownum,4}=message;

Param = Drift; rownum = 14; name = 'Drift'
Val_sorted = Param(t1_order);
[Shewart,message] = ShewartEval(Val_sorted);
data{rownum,1} = name; data{rownum,2}=Param(end); data{rownum,3}=Shewart; data{rownum,4}=message;

Param = RDC; rownum = 15; name = 'Radius of Decorrelation'
Val_sorted = Param(t1_order);
[Shewart,message] = ShewartEval(Val_sorted);
data{rownum,1} = name; data{rownum,2}=Param(end); data{rownum,3}=Shewart; data{rownum,4}=message;

close all;
f = figure(5); set(gcf, 'units','normalized','outerposition',[0.05 0.05 .9 .9]); 
t = uitable(f,'Data',data,'ColumnWidth',{170,110,50,280})
t.ColumnName = {'Measurement','Value','Shewart','Actionable?'};
subplot(1,2,2),plot(3)
pos = get(subplot(1,2,2),'position');
delete(subplot(2,1,2))
set(t,'units','normalized')
set(t,'position',pos)
axis off; 

subplot(2,2,1); imagesc(D); axis equal; axis off; title('Ghosting ROIs');
subplot(2,2,3); imagesc(E); axis equal; axis off; title('Distortion Measures');

saveas(gcf,[RESULTS_dir filesep 'EPI_SingleDay.fig'])
print_current_figure(400,[RESULTS_dir filesep 'EPI_SingleDay.png'])
pause(2); close all;


