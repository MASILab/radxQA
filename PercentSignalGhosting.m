function [PSG,pf_hdl]=PercentSignalGhosting(img,reso)
%ROI: 1=water ROI, 2=top ROI, 3=bottom ROI, 4=left ROI,5=right ROI

% get slice 7
Slice = squeeze(img(:,:,7));
pxl_sz = reso(:);
I = imrotate(flip(Slice,1),90);
per = .1;
I_max=double(max(max(I)));
[hist_cnt,hist_int]=hist(I(:),0:I_max);
hist_sample_start=round(I_max*per);%percentage of max as min to exclude air intenisty
%2.find max within specified histogram range
[int_cnt,int_pk]=max(hist_cnt(hist_sample_start:size(hist_cnt,2)));%v2
int_pk=int_pk+hist_sample_start-1;%v2

water_mean=int_pk;
I_bin = I > water_mean/2;

% find center of phantom
band_per=[.3 .6];
% sum each colum position across central rows, get low (L) and high (R) indices
row_band=zeros(1,size(I_bin,2));
for i=round(band_per(1,1)*size(I_bin,1)): round(band_per(1,2)*size(I_bin,1))
    row_band=row_band+double(I_bin(i,:));
end
ind_row_low=find(row_band>0,1);%v3
ind_row_high=find(row_band>0,1,'last');%v3

% do same across columns
col_band=zeros(size(I_bin,1),1);
for i=round(band_per(1,1)*size(I_bin,2)):...
        round(band_per(1,2)*size(I_bin,2))
    col_band=col_band+double(I_bin(:,i));
end
ind_col_low=find(col_band>0,1);%v3
ind_col_high=find(col_band>0,1,'last');%v3

ind_centre=[round((ind_col_high-ind_col_low)/2+ind_col_low) round((ind_row_high-ind_row_low)/2+ind_row_low)];%phantom centre (y,x)

%6.create circular ROI around phantom centre with user specified radius
 area_ROI_water=200;
 
 radius_cm=sqrt(area_ROI_water/pi);
radius_mm=radius_cm*10;
radius_pxl=radius_mm/pxl_sz(1,1);%from this line, use ellipse ROI if image
theta=0:0.01:2*pi;               %is anisotropic (see test.m for working)
x=radius_pxl*cos(theta)+ind_centre(1,2);
y=radius_pxl*sin(theta)+ind_centre(1,1)+10;%HW:centre of ROI is 10 pxls%below phantom centre

BW=roipoly(I,x,y);

I_ROI_water = I.*BW;
I_ROI_water_sum=sum(I(BW));%this finds sum of water within mask
I_ROI_water_mean=I_ROI_water_sum/size(I(BW),1);%mean of water ROI
I_ROI_water_sigma=std(double(I(BW)));%std of water ROI
mu_S7=I_ROI_water_mean;

%2.use mean and std of water to mask image
I_bin = I > mu_S7/2;

% find center of phantom
band_per=[.3 .6];
% sum each colum position across central rows, get low (L) and high (R) indices
row_band=zeros(1,size(I_bin,2));
for i=round(band_per(1,1)*size(I_bin,1)): round(band_per(1,2)*size(I_bin,1))
    row_band=row_band+double(I_bin(i,:));
end
ind_row_low=find(row_band>0,1);%v3
ind_row_high=find(row_band>0,1,'last');%v3

% do same across columns
col_band=zeros(size(I_bin,1),1);
for i=round(band_per(1,1)*size(I_bin,2)):...
        round(band_per(1,2)*size(I_bin,2))
    col_band=col_band+double(I_bin(:,i));
end
ind_col_low=find(col_band>0,1);%v3
ind_col_high=find(col_band>0,1,'last');%v3

ind_centre=[round((ind_col_high-ind_col_low)/2+ind_col_low) round((ind_row_high-ind_row_low)/2+ind_row_low)];%phantom centre (y,x)

FOV_row = size(I,1);
FOV_col = size(I,2);

%6.calculate normal ellipse ROI pixel length & width
w_cm=sqrt(10/(pi));%HW:ellipse are=10cm2%v4
l_cm=4*w_cm;%HW:length:width=4:1
w_mm=w_cm*10;
l_mm=l_cm*10;
w_pxl=w_mm/pxl_sz(1,1);
l_pxl=l_mm/pxl_sz(2,1);

%%

%7.create ROI mask image
[centre_coord_2,contour_xy_2,I_bin_2]=MakeEllipseROI(I,ind_centre,ind_col_low,0,[w_pxl,l_pxl],pxl_sz,2);
[centre_coord_3,contour_xy_3,I_bin_3]=MakeEllipseROI(I,ind_centre,ind_col_high,FOV_row,[w_pxl,l_pxl],pxl_sz,3);
[centre_coord_4,contour_xy_4,I_bin_4]=MakeEllipseROI(I,ind_centre,ind_row_low,0,[w_pxl,l_pxl],pxl_sz,4);
[centre_coord_5,contour_xy_5,I_bin_5]=MakeEllipseROI(I,ind_centre,ind_row_high,FOV_col,[w_pxl,l_pxl],pxl_sz,5);
%8.apply mask & get the mean intensity value & std
I_ROI_2 = I.*I_bin_2;
I_ROI_2_sum=sum(I(I_bin_2));
I_ROI_2_mean=I_ROI_2_sum/size(I(I_bin_2),1);
I_ROI_2_sigma=std(double(I(I_bin_2)));
I_ROI_3=I.*I_bin_3;
I_ROI_3_sum=sum(I(I_bin_3));
I_ROI_3_mean=I_ROI_3_sum/size(I(I_bin_3),1);
I_ROI_3_sigma=std(double(I(I_bin_3)));
I_ROI_4=I.*I_bin_4;
I_ROI_4_sum=sum(I(I_bin_4));
I_ROI_4_mean=I_ROI_4_sum/size(I(I_bin_4),1);
I_ROI_4_sigma=std(double(I(I_bin_4)));
I_ROI_5=I.*I_bin_5;
I_ROI_5_sum=sum(I(I_bin_5));
I_ROI_5_mean=I_ROI_5_sum/size(I(I_bin_5),1);
I_ROI_5_sigma=std(double(I(I_bin_5)));
%9.find PSG
PSG=abs(((I_ROI_2_mean+I_ROI_3_mean)-(I_ROI_4_mean+I_ROI_5_mean))/(2*mu_S7));
%10.check if passed test
if PSG<=0.025
    pf_hdl=1;
else
    pf_hdl=0;
end
%11.display the ROI on image
figure;
imshow(I,[]);
hold on
plot(contour_xy_2(1,:),contour_xy_2(2,:),'Color','r','LineWidth',1);
plot(contour_xy_3(1,:),contour_xy_3(2,:),'Color','r','LineWidth',1);
plot(contour_xy_4(1,:),contour_xy_4(2,:),'Color','r','LineWidth',1);
plot(contour_xy_5(1,:),contour_xy_5(2,:),'Color','r','LineWidth',1);
hold off
end