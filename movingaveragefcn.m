function [count]=movingaveragefcn(count,numpts)
a=1;
for ii=1:numpts
    b(ii)=1/numpts;
end
count=filter(b,a,count);
end
