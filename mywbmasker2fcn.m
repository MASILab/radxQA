function wbmask=mywbmasker2fcn(IM,conn)
% Returns the largest or second largest cluster of connected voxels in IM exceeding a certain cutoff value.
% The cutoff used based on analysis of the histogram of image values.
% Idea: alter this to minimize the ratio of sum(wbperim(@#$))/# of voxels in mask
switch ndims(IM)
    case 2
        IM=IM(:,:,1);
    case 4
        IM=squeeze(IM(:,:,:,1));
end
if ~exist('conn','var')
    conn=18;
end
% begin by thresholding the image
[count,xloc]=hist(IM(:),100);
count=movingaveragefcn(count,5);
count(1:4)=count(5);
changes=zeros(size(count));
for ii=2:length(count)
    if count(ii)>count(ii-1)
        changes(ii)=1;
    elseif count(ii)<count(ii-1)
        changes(ii)=-1;
    end
end
changes2=abs(changes(2:end)-changes(1:end-1));
index=find(changes2>0);
wbmask=IM>xloc(index(2));
origwbmask=wbmask;
% keep only the 2nd largest cluster of voxels 
L=bwlabeln(wbmask,conn);
clusters=unique(L(:));
for ii=1:length(clusters);
    clustersize(ii)=length(find(L==clusters(ii)));
end
sorted=sort(clustersize,'descend');
wbmask=(L==clusters(find(clustersize==sorted(2)))); %this finds the 2nd largest cluster because the background is likely the largest cluster.
% fill in isolated pixels 
if ndims(IM)==3
    sizes=size(wbmask);
    for ii=1:sizes(3)
        tmpwb=squeeze(wbmask(:,:,ii));
        wbmask(:,:,ii)=imfill(wbmask(:,:,ii),'holes');
    end
elseif ndims(IM==2)
    wbmask=imfill(wbmask,'holes');
end
% test to see if more than half of the slice corners are included in the
% wbmask.  If so, I will make an assumption that the background was smaller
% than the brain, and that I need to keep the largest cluster of connected
% components, not the second largest.
expression='wbmask(';
for ii=1:ndims(wbmask)
    if ii<3
        expression=cat(2,expression,'[1,end],');
    else
        expression=cat(2,expression,':,');
    end
end
expression=[expression(1:end-1),');'];
tmp=(eval(expression));
if round(mean(tmp(:)))==1
    wbmask=origwbmask;
    % keep only the largest cluster of voxels 
    L=bwlabeln(wbmask,conn);
    clusters=unique(L(:));
    for ii=1:length(clusters);
        clustersize(ii)=length(find(L==clusters(ii)));
    end
    sorted=sort(clustersize,'descend');
    wbmask=(L==clusters(find(clustersize==sorted(1)))); % this finds the largest cluster because the background has been removed.
    % fill in isolated pixels 
    if ndims(IM)==3
        sizes=size(wbmask);
        for ii=1:sizes(3)
            tmpwb=squeeze(wbmask(:,:,ii));
            wbmask(:,:,ii)=imfill(wbmask(:,:,ii),'holes');
        end
    elseif ndims(IM==2)
        wbmask=imfill(wbmask,'holes');
    end
end
end
