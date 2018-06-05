function [l_diff,sl_l_diff,pf_hdl]= SlicePositionAccuracyS11(img,reso)

% get slice 1
Slice = squeeze(img(:,:,11));
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

ind_centre=[round((ind_col_high-ind_col_low)/2+ind_col_low) round((ind_row_high-ind_row_low)/2+ind_row_low)];%phantom centre (y,x)

%5.go down by 10mm & sample the wedge left&right bndry
l_pxl=round(10/pxl_sz(1,1));%HW:look at 10mm below top bndry to get width

sample_row=I_bin(ind_col_low+l_pxl,ind_centre(1,2)-2*l_pxl:ind_centre(1,2)+2*l_pxl);
bndry_l=find(1-sample_row,1,'first')+ind_centre(1,2)-2*l_pxl-1;
bndry_r=find(1-sample_row,1,'last')+ind_centre(1,2)-2*l_pxl-1;

%8.find bndry of each wedge
integertest=~mod((bndry_r-bndry_l)/2,1);
if integertest%if even rows
    bndry_1=(bndry_r-bndry_l)/2+bndry_l-1;%HW
    bndry_2=bndry_r-(bndry_r-bndry_l)/2+1;%HW
else%if odd rows
    bndry_1=floor((bndry_r-bndry_l)/2+bndry_l)-1;%HW
    bndry_2=round((bndry_r-bndry_l)/2+bndry_l)+1;%HW
end
l_wedge_bndry=[bndry_l bndry_1];
r_wedge_bndry=[bndry_2 bndry_r];

%9.find length inside each wedge and put into vectors
dummy1=zeros();
dummy2=zeros();
row_start=ind_col_low+l_pxl;
row_end=ind_col_low+round(44/pxl_sz(1,1));%v7
cnt_ind=1;
for j=l_wedge_bndry(1,1):l_wedge_bndry(1,2)
    cnt=1;
    for i=row_start:row_end
        if I_bin(i,j)<1
            cnt=cnt+1;
        else
            break;
        end
    end
    dummy1(1,cnt_ind)=cnt;
    cnt_ind=cnt_ind+1;
end
cnt_ind=1;
for j=r_wedge_bndry(1,1):r_wedge_bndry(1,2)
    cnt=1;
    for i=row_start:row_end
        if I_bin(i,j)<1
            cnt=cnt+1;
        else
            break;
        end
    end
    dummy2(1,cnt_ind)=cnt;
    cnt_ind=cnt_ind+1;
end

%10.search length result & delete the one outside std (false result)
dummy1_mu=mean(dummy1);
dummy1_std=std(dummy1);
dummy2_mu=mean(dummy2);
dummy2_std=std(dummy2);
dummy1_diff=abs(dummy1-dummy1_mu);
dummy2_diff=abs(dummy2-dummy2_mu);
dummy1_true=find(dummy1_diff<=dummy1_std);
dummy1=mean(dummy1(dummy1_true));%v8:take mean of true value as result
dummy2_true=find(dummy2_diff<=dummy2_std);
dummy2=mean(dummy2(dummy2_true));%v8:take mean of true value as result

%11.convert pxl to mm
dummy1=dummy1*pxl_sz(2,1);
dummy2=dummy2*pxl_sz(2,1);

%12.find the length difference using difference between mean lengths
length_l=mean(dummy1);
length_r=mean(dummy2);
l_diff=length_r-length_l;%use right minus left
sl_l_diff = l_diff/2;

%12.show result and create pass/fail handle
if l_diff<5||l_diff>-5
    pf_hdl=1;
else
    pf_hdl=0;
end

end