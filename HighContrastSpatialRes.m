function [HCSR,pf_hdl] = HighContrastSpatialRes(img,reso)
%function [HCSR,pf_hdl]=fun_ACR_2_S1(dir_name,file_name,visual,manual,imag_check,myContrast)
Slice = squeeze(img(:,:,1));
pxl_sz = reso(:);
I = imrotate(flip(Slice,1),90);
myContrast = .15;
hole_coord=round([-23.44 29.30;-6.84 36.13;0 29.30;15.63 36.13;...
    23.44 30.27;38.09 36.13]/pxl_sz(1));%input vec num in mm

%3.use Otsu's method to threshold image
HCSR_vec=zeros(2,3);%predefine result vector

per = .1;
visual = 1;
% Find the water intensity peak on the image intensity histogram.
% This intensity peak will be used to threshold image.

%1.specify intensity range
I_max=double(max(max(I)));
[hist_cnt,hist_int]=hist(I(:),0:I_max);
hist_sample_start=round(I_max*per);%percentage of max as min to exclude air intenisty
%2.find max within specified histogram range
[int_cnt,int_pk]=max(hist_cnt(hist_sample_start:size(hist_cnt,2)));%v2
int_pk=int_pk+hist_sample_start-1;%v2
mu = int_pk;
% threshold
I_bin = I > mu/2;

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
% phantom center
ind_centre=[round((ind_col_high-ind_col_low)/2+ind_col_low) round((ind_row_high-ind_row_low)/2+ind_row_low)];

%5.find the pixel location of UL & LR corner of 3 pairs
hole_loc=zeros(6,2);
for i=1:6
    hole_loc(i,:)=hole_coord(i,:)+fliplr(ind_centre);
end
sample_wdth=round(10.74/pxl_sz(1));%HW:sample width based on pxlsz=0.9766mm
sample_wdth_09=round(9.766/pxl_sz(1));%shorter for 0.9mm one

%6.plot all rows of UL holes in all 3 pairs into matrix
hole_prof_row_1=zeros(sample_wdth,sample_wdth);
hole_prof_row_2=zeros(sample_wdth,sample_wdth);
hole_prof_row_3=zeros(sample_wdth_09,sample_wdth_09);

for i=1:sample_wdth%HW:the row number is fixed to suit order in hole_loc
    hole_prof_row_1(i,:)=...
        I(hole_loc(1,2)+i-1,hole_loc(1,1):hole_loc(1,1)+sample_wdth-1);
    hole_prof_row_2(i,:)=...
        I(hole_loc(3,2)+i-1,hole_loc(3,1):hole_loc(3,1)+sample_wdth-1);
    %             hole_prof_row_3(i,:)=...
    %                 I(hole_loc(5,2)+i-1,hole_loc(5,1):hole_loc(5,1)+sample_wdth-1);
end

for i=1:sample_wdth_09
    hole_prof_row_3(i,:)=...
        I(hole_loc(5,2)+i-1,hole_loc(5,1):hole_loc(5,1)+sample_wdth_09-1);
end


%7.plot all cols of LR holes in all 3 pairs into matrix
hole_prof_col_1=zeros(sample_wdth,sample_wdth);
hole_prof_col_2=zeros(sample_wdth,sample_wdth);
hole_prof_col_3=zeros(sample_wdth_09,sample_wdth_09);
for i=1:sample_wdth%HW:the row number is fixed to suit order in hole_loc
    hole_prof_col_1(i,:)=I(hole_loc(2,2):hole_loc(2,2)+sample_wdth-1,...
        hole_loc(2,1)-sample_wdth+i);
    hole_prof_col_2(i,:)=I(hole_loc(4,2):hole_loc(4,2)+sample_wdth-1,...
        hole_loc(4,1)-sample_wdth+i);
    %             hole_prof_col_3(i,:)=I(hole_loc(6,2):hole_loc(6,2)+sample_wdth-1,...
    %                 hole_loc(6,1)-sample_wdth+i);
end
for i=1:sample_wdth_09
    hole_prof_col_3(i,:)=I(hole_loc(6,2):hole_loc(6,2)+sample_wdth_09-1,...
        hole_loc(6,1)-sample_wdth+i);
end

%====================v5 start====================
%8.find peaks & compare to user's contrast
for i=1:size(hole_prof_row_1,1)
    %disp(i)
    dummy=hole_prof_row_1(i,:);
    %             [pks_r_1 locs_r_1]=findpeaks(dummy);
    [locs_r_1, pks_r_1]=peakfinder(dummy,0);%v6
    if size(pks_r_1,2)~=4
        continue
    elseif size(pks_r_1,2)==4
        vly_r_1=[min(dummy(locs_r_1(1,1)+1:locs_r_1(1,2)-1)) ...
            min(dummy(locs_r_1(1,2)+1:locs_r_1(1,3)-1)) ...
            min(dummy(locs_r_1(1,3)+1:locs_r_1(1,4)-1))];
        dummy_vly=[vly_r_1(1,1),vly_r_1(1,1),vly_r_1(1,2),...
            vly_r_1(1,2),vly_r_1(1,3),vly_r_1(1,3)];
        dummy_pks=[pks_r_1(1,1),pks_r_1(1,2),pks_r_1(1,2),...
            pks_r_1(1,3),pks_r_1(1,3),pks_r_1(1,4)];
        for l=1:6
            dummy_diff(1,l)=abs(dummy_vly(1,l)-dummy_pks(1,l))/...
                (dummy_vly(1,l)+dummy_pks(1,l));
        end
        if sum(dummy_diff<myContrast*ones(1,6))==0
            HCSR_vec(1,1)=1;
            break
        end
    end
end


for i=1:size(hole_prof_row_2,1)
    %disp(i)
    dummy=hole_prof_row_2(i,:);
    %             [pks_r_1 locs_r_1]=findpeaks(dummy);
    [locs_r_1 pks_r_1]=peakfinder(dummy,0);%v6
    if size(pks_r_1,2)~=4
        continue
    elseif size(pks_r_1,2)==4
        vly_r_1=[min(dummy(locs_r_1(1,1)+1:locs_r_1(1,2)-1)) ...
            min(dummy(locs_r_1(1,2)+1:locs_r_1(1,3)-1)) ...
            min(dummy(locs_r_1(1,3)+1:locs_r_1(1,4)-1))];
        dummy_vly=[vly_r_1(1,1),vly_r_1(1,1),vly_r_1(1,2),...
            vly_r_1(1,2),vly_r_1(1,3),vly_r_1(1,3)];
        dummy_pks=[pks_r_1(1,1),pks_r_1(1,2),pks_r_1(1,2),...
            pks_r_1(1,3),pks_r_1(1,3),pks_r_1(1,4)];
        for l=1:6
            dummy_diff(1,l)=abs(dummy_vly(1,l)-dummy_pks(1,l))/...
                (dummy_vly(1,l)+dummy_pks(1,l));
        end
        if sum(dummy_diff<myContrast*ones(1,6))==0
            HCSR_vec(1,2)=1;

            break
        end
    end
end

for i=1:size(hole_prof_row_3,1)
    %disp(i)
    dummy=hole_prof_row_3(i,:);
    %             [pks_r_1 locs_r_1]=findpeaks(dummy);
    [locs_r_1 pks_r_1]=peakfinder(dummy,0);%v6
    if size(pks_r_1,2)~=4
        continue
    elseif size(pks_r_1,2)==4
        vly_r_1=[min(dummy(locs_r_1(1,1)+1:locs_r_1(1,2)-1)) ...
            min(dummy(locs_r_1(1,2)+1:locs_r_1(1,3)-1)) ...
            min(dummy(locs_r_1(1,3)+1:locs_r_1(1,4)-1))];
        dummy_vly=[vly_r_1(1,1),vly_r_1(1,1),vly_r_1(1,2),...
            vly_r_1(1,2),vly_r_1(1,3),vly_r_1(1,3)];
        dummy_pks=[pks_r_1(1,1),pks_r_1(1,2),pks_r_1(1,2),...
            pks_r_1(1,3),pks_r_1(1,3),pks_r_1(1,4)];
        for l=1:6
            dummy_diff(1,l)=abs(dummy_vly(1,l)-dummy_pks(1,l))/...
                (dummy_vly(1,l)+dummy_pks(1,l));
        end
        if sum(dummy_diff<myContrast*ones(1,6))==0
            HCSR_vec(1,3)=1;

            break
        end
    end
end

for i=1:size(hole_prof_col_1,1)
    %disp(i)
    dummy=hole_prof_col_1(i,:);
    %             [pks_r_1 locs_r_1]=findpeaks(dummy);
    [locs_r_1 pks_r_1]=peakfinder(dummy,0);%v6
    if size(pks_r_1,2)~=4
        continue
    elseif size(pks_r_1,2)==4
        vly_r_1=[min(dummy(locs_r_1(1,1)+1:locs_r_1(1,2)-1)) ...
            min(dummy(locs_r_1(1,2)+1:locs_r_1(1,3)-1)) ...
            min(dummy(locs_r_1(1,3)+1:locs_r_1(1,4)-1))];
        dummy_vly=[vly_r_1(1,1),vly_r_1(1,1),vly_r_1(1,2),...
            vly_r_1(1,2),vly_r_1(1,3),vly_r_1(1,3)];
        dummy_pks=[pks_r_1(1,1),pks_r_1(1,2),pks_r_1(1,2),...
            pks_r_1(1,3),pks_r_1(1,3),pks_r_1(1,4)];
        for l=1:6
            dummy_diff(1,l)=abs(dummy_vly(1,l)-dummy_pks(1,l))/...
                (dummy_vly(1,l)+dummy_pks(1,l));
        end
        if sum(dummy_diff<myContrast*ones(1,6))==0
            HCSR_vec(2,1)=1;
           
            break
        end
    end
end

for i=1:size(hole_prof_col_2,1)
    %disp(i)
    dummy=hole_prof_col_2(i,:);
    %             [pks_r_1 locs_r_1]=findpeaks(dummy);
    [locs_r_1 pks_r_1]=peakfinder(dummy,0);%v6
    if size(pks_r_1,2)~=4
        continue
    elseif size(pks_r_1,2)==4
        vly_r_1=[min(dummy(locs_r_1(1,1)+1:locs_r_1(1,2)-1)) ...
            min(dummy(locs_r_1(1,2)+1:locs_r_1(1,3)-1)) ...
            min(dummy(locs_r_1(1,3)+1:locs_r_1(1,4)-1))];
        dummy_vly=[vly_r_1(1,1),vly_r_1(1,1),vly_r_1(1,2),...
            vly_r_1(1,2),vly_r_1(1,3),vly_r_1(1,3)];
        dummy_pks=[pks_r_1(1,1),pks_r_1(1,2),pks_r_1(1,2),...
            pks_r_1(1,3),pks_r_1(1,3),pks_r_1(1,4)];
        for l=1:6
            dummy_diff(1,l)=abs(dummy_vly(1,l)-dummy_pks(1,l))/...
                (dummy_vly(1,l)+dummy_pks(1,l));
        end
        if sum(dummy_diff<myContrast*ones(1,6))==0
            HCSR_vec(2,2)=1;

            break
        end
    end
end


for i=1:size(hole_prof_col_3,1)
    dummy=hole_prof_col_3(i,:);
    %             [pks_r_1 locs_r_1]=findpeaks(dummy);
    [locs_r_1 pks_r_1]=peakfinder(dummy,0);%v6
    if size(pks_r_1,2)~=4
        continue
    elseif size(pks_r_1,2)==4
        vly_r_1=[min(dummy(locs_r_1(1,1)+1:locs_r_1(1,2)-1)) ...
            min(dummy(locs_r_1(1,2)+1:locs_r_1(1,3)-1)) ...
            min(dummy(locs_r_1(1,3)+1:locs_r_1(1,4)-1))];
        dummy_vly=[vly_r_1(1,1),vly_r_1(1,1),vly_r_1(1,2),...
            vly_r_1(1,2),vly_r_1(1,3),vly_r_1(1,3)];
        dummy_pks=[pks_r_1(1,1),pks_r_1(1,2),pks_r_1(1,2),...
            pks_r_1(1,3),pks_r_1(1,3),pks_r_1(1,4)];
        for l=1:6
            dummy_diff(1,l)=abs(dummy_vly(1,l)-dummy_pks(1,l))/...
                (dummy_vly(1,l)+dummy_pks(1,l));
        end
        if sum(dummy_diff<myContrast*ones(1,6))==0
            HCSR_vec(2,3)=1;

            break
        end
    end
end


dummy=sum(HCSR_vec,2);
if dummy(1,1)==3
    HCSR(1,1)=0.9;
elseif dummy(1,1)==2
    HCSR(1,1)=1.0;
elseif dummy(1,1)==1
    HCSR(1,1)=1.1;
end
if dummy(2,1)==3
    HCSR(2,1)=0.9;
elseif dummy(2,1)==2
    HCSR(2,1)=1.0;
elseif dummy(2,1)==1
    HCSR(2,1)=1.1;
end

%13.pass/fail handle
if HCSR(1,1)<=1%v4
    pf_hdl(1,1)=1;
else
    pf_hdl(1,1)=0;
end
if HCSR(2,1)<=1%v4
    pf_hdl(1,2)=1;
else
    pf_hdl(1,2)=0;
end

end
