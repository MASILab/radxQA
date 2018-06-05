function centerpt=findcentroid(IM)
tmp=sum(IM,3);
centerpt=[0,0];
tmpp=squeeze(sum(tmp,2));
tmpp=tmpp(:);
for ii=1:size(tmpp,1);
    centerpt(1)=centerpt(1)+(ii*tmpp(ii));
end
centerpt(1)=centerpt(1)./sum(IM(:));
tmpp=squeeze(sum(tmp,1));
tmpp=tmpp(:);
for ii=1:size(tmpp,1);
    centerpt(2)=centerpt(2)+(ii*tmpp(ii));
end
centerpt(2)=centerpt(2)./sum(IM(:));
end
