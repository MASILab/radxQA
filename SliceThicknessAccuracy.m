function [slc_thk,pf_slthick] = SliceThicknessAccuracy(img,reso)

Slice = squeeze(img(:,:,1));
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

ind_centre=[round((ind_col_high-ind_col_low)/2+ind_col_low) round((ind_row_high-ind_row_low)/2+ind_row_low)];%centre of phantom
centre=ind_centre; tol = water_mean/4;
% find horizontal boundary of ramp
%1.read in up direction and stops at big intensity change
mu=I(centre(1,1),centre(1,2));
for i=1:size(I,1)/2
    %disp(i)
    mu=(I(centre(1,1)-i,centre(1,2))+mu)/2;
    diff_mu=abs(mu-I(centre(1,1),centre(1,2)));
    if diff_mu>=tol
        break%remember to +1 when defining bndry, see below
    end
end
bndry_low=centre(1,1)-i+1;%top bndry
%2.read in down direction and stops at big intensity change
mu=I(centre(1,1),centre(1,2));
for i=1:size(I,1)/2
    %disp(i)
    mu=(I(centre(1,1)+i,centre(1,2))+mu)/2;
    diff_mu=abs(mu-I(centre(1,1),centre(1,2)));
    if diff_mu>=tol
        break%remember to -1 when defining bndry, see below
    end
end
bndry_high=centre(1,1)+i-1;%bottom bndry

bndry_low=bndry_low+1;%NB
bndry_high=bndry_high-1;%NB

integertest=~mod((bndry_high-bndry_low)/2,1);

if integertest%if even rows
    ROI_top_ind=[bndry_low,(bndry_high-bndry_low)/2+bndry_low];
    ROI_bottom_ind=[(bndry_high-bndry_low)/2+bndry_low,bndry_high];
else%if odd rows
    ROI_top_ind=[bndry_low,round((bndry_high-bndry_low)/2+bndry_low)-1];
    ROI_bottom_ind=[ceil((bndry_high-bndry_low)/2+bndry_low),bndry_high];
end

% average value of a rectangular ROI, that is specified by the input boundary. 
bndry_row = ROI_top_ind+[1 0];
bndry_col = [ind_centre(1,2)-10,ind_centre(1,2)+10];
%1.take sum by rows
sum_row=sum(I(bndry_row(1,1):bndry_row(1,2),...
    bndry_col(1,1):bndry_col(1,2)),1);
%2.take sum by col to find total summation
sum_total=sum(sum_row,2);
%3.find mean
ROI_top=sum_total/((bndry_row(1,2)-bndry_row(1,1)+1)*...
    (bndry_col(1,2)-bndry_col(1,1)+1));

bndry_row = ROI_bottom_ind-[1 0];
bndry_col = [ind_centre(1,2)-10,ind_centre(1,2)+10];
sum_row=sum(I(bndry_row(1,1):bndry_row(1,2),...
    bndry_col(1,1):bndry_col(1,2)),1);
%2.take sum by col to find total summation
sum_total=sum(sum_row,2);
%3.find mean
ROI_bottom=sum_total/((bndry_row(1,2)-bndry_row(1,1)+1)*...
    (bndry_col(1,2)-bndry_col(1,1)+1));

% fprintf('Top boundary of top ROI is %i\n',ROI_top_ind(1,1));%display bndry
% fprintf('Bottom boundary of top ROI is %i\n',ROI_top_ind(1,2));%on screen
% fprintf('Top boundary of Bottom ROI is %i\n',ROI_bottom_ind(1,1));
% fprintf('Bottom boundary of Bottom ROI is %i\n',ROI_bottom_ind(1,2));


ROI_check=abs(ROI_top-ROI_bottom)/mean([ROI_top,ROI_bottom]);

if ROI_check>=1 %if >1 then water inside ROI, reduce boundary by 1 pixel
    ROI_top=fun_ACR_FindAveRectROI(I,...
        [ROI_top_ind(1,1)+1 ROI_top_ind(1,2)],...
        [ind_centre(1,2)-20,ind_centre(1,2)+20]);%HW:hori l of ROI=40 pxls
    ROI_bottom=fun_ACR_FindAveRectROI(I,...
        [ROI_bottom_ind(1,1) ROI_bottom_ind(1,2)-1],...
        [ind_centre(1,2)-20,ind_centre(1,2)+20]);%HW:hori l of ROI=40 pxls
    ROI_check=abs(ROI_top-ROI_bottom)/mean([ROI_top,ROI_bottom]);
end

ramp_mean=(ROI_top+ROI_bottom)/2;

%4.threshold image based on ramp mean intensity
I_bin_l = I > ramp_mean/2;

%I_bin_l=add_threshold(I,ramp_mean/2);%mask image to measure ramp length

%5.define middle of ramp
half_of_bndry=round((bndry_high-bndry_low)/2+bndry_low);
%++++++++++++++++++v3 start+++++++++++++++++++++++
%6.define 2 n-by-1 vectors to store the pxl length of 2 ROI
dummy1=zeros();
dummy2=zeros();
cnt=1;
for i=bndry_low:half_of_bndry-1
    centre_hori=ind_centre(1,2);
    ind_row_ROI=i;
    %1.find peaks of row intensity vector, result is left side of peak
    row_ROI=I_bin_l(ind_row_ROI,:);
    
    pks=0;
    cnt_pks=1;
    %2.search 1 within row vector
    for ii=2:size(row_ROI,2)
        if row_ROI(1,ii-1)==0 && row_ROI(1,ii)==1
            pks(1,cnt_pks)=ii;
            cnt_pks=cnt_pks+1;
        end
    end
    
    %2.find the peak closest to ramp centre, this is one isde of length
    distance=abs(pks-centre_hori);
    [~,min_distance_ind]=min(distance);
    l_a=pks(min_distance_ind);
    %3.find the other side of length based on the choice
    for ii=l_a:size(row_ROI,2)
        if row_ROI(1,ii)==0;
            break
        end
        l_b=ii;
    end
    %4.find the pixel distance
    l=abs(l_a-l_b);
    dummy1(cnt,1) = l;
    cnt=cnt+1;
end

cnt=1;
for i=half_of_bndry-1:bndry_high
    centre_hori=ind_centre(1,2);
    ind_row_ROI=i;
    %1.find peaks of row intensity vector, result is left side of peak
    row_ROI=I_bin_l(ind_row_ROI,:);
    % [~,pks]=findpeaks(double(row_ROI));%double is required
    pks=0;
    cnt_pks=1;
    %2.search 1 within row vector
    for ii=2:size(row_ROI,2)
        if row_ROI(1,ii-1)==0 && row_ROI(1,ii)==1
            pks(1,cnt_pks)=ii;
            cnt_pks=cnt_pks+1;
        end
    end
    %2.find the peak closest to ramp centre, this is one isde of length
    distance=abs(pks-centre_hori);
    [~,min_distance_ind]=min(distance);
    l_a=pks(min_distance_ind);
    %3.find the other side of length based on the choice
    for ii=l_a:size(row_ROI,2)
        if row_ROI(1,ii)==0;
            break
        end
        l_b=ii;
    end
    %4.find the pixel distance
    l=abs(l_a-l_b);
    dummy2(cnt,1)=l;
    
    cnt=cnt+1;
end

%7.convert pxl to mm
dummy1=dummy1*pxl_sz(1,1);
dummy2=dummy2*pxl_sz(1,1);

%8.find the length most close to 5mm
dummy11=abs(dummy1-50);
dummy22=abs(dummy2-50);

[~,ind]=min(dummy11);
l_top=dummy1(ind);
[~,ind]=min(dummy22);
I_bottom=dummy2(ind);
slc_thk=0.2*((l_top*I_bottom)/(l_top+I_bottom));

%9.tell user if test pass/fail and create pass/fail handle
if slc_thk>=4.3&&slc_thk<=5.7

    pf_slthick=1;
else

    pf_slthick=0;
end








end