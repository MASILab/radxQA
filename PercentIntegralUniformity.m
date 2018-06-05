function [PIU,PIU2]=PercentIntegralUniformity(img,reso)

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

C_x=ind_centre(1,2);
C_y=ind_centre(1,1)+10;%+10 to get ROI centre y cood
R=radius_pxl;
radius_cm_small=sqrt(1/pi);%1cm2 circle
radius_mm_small=radius_cm_small*10;
radius_pxl_small= radius_mm_small/pxl_sz(1,1); %this can be reused for max mean intensity
r=radius_pxl_small;
center=[C_x C_y];
radius=r;
n=1000;

THETA = linspace(0, 2 * pi, n);
RHO = ones(1, n) * radius;
[X Y] = pol2cart(THETA, RHO);
xcc = X + center(1);
ycc = Y + center(2);

BW_small=roipoly(I,xcc,ycc);
        
I_ROI_small = I.*BW_small;

I_ROI_small_sum=sum(sum(I_ROI_small));%this finds sum of intensity within mask
I_ROI_small_mean(1,1)=I_ROI_small_sum/size(find(I_ROI_small),1);%mean of min ROI

numlapis=((2*R)-(R+r))/(2*r);

% figure;
%             imshow(I_ROI_water,[]);hold on;
%             plot(xcc,ycc,'-','linewidth',2,'color',0.5.*rand(1,3));
%             hold on;
            
cnt=2;
for cnt1=1:numlapis
    lapis(cnt1)=cnt1*6;
    
    center=[C_x C_y];
    radius=cnt1*2*r;
    n=lapis(cnt1)+1;
    
    THETA = linspace(0, 2 * pi, n);
    RHO = ones(1, n) * radius;
    [X Y] = pol2cart(THETA, RHO);
    xcoor = X + center(1);
    ycoor = Y + center(2);
    for cnt2=1:lapis(cnt1)
          center=[xcoor(cnt2) ycoor(cnt2)];
    radius=r;
    n=1000;
     
     THETA = linspace(0, 2 * pi, n);
    RHO = ones(1, n) * radius;
    [X Y] = pol2cart(THETA, RHO);
    xc = X + center(1);
    yc = Y + center(2);
        
        BW_small=roipoly(I,xc,yc);
        
        I_ROI_small = I.*BW_small;
        
        I_ROI_small_sum=sum(sum(I_ROI_small));%this finds sum of intensity within mask
        I_ROI_small_mean(cnt,1)=I_ROI_small_sum/size(find(I_ROI_small),1);%mean of min ROI
        cnt=cnt+1;
        
%             plot(xc,yc,'-','linewidth',2,'color',0.5.*rand(1,3));
        
    end
end
PIU=1-(max(I_ROI_small_mean)-min(I_ROI_small_mean))/(max(I_ROI_small_mean)+min(I_ROI_small_mean));

% alternative way
% smooth size of radius (pixels), then find max and min in 
B = imgaussfilt(I,radius_pxl_small,'FilterSize',11);
B = imgaussfilt(I,radius_pxl_small);

B_masked = B.*BW;
high = max(B_masked(BW>0));
low = min(B_masked(BW>0));
PIU2 = 1 - ((high-low)/(high+low));







end