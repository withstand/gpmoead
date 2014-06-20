function [pf, ps] = findparetosolution(val,pop)
pf = [];
ps = [];

nval = size(val,1);
for iPf=1:nval
    do = val([1:iPf-1 iPf+1:nval],:)-...
        repmat(val(iPf,:), nval-1,1)>=0;
    if sum(sum(do,2)==0) == 0
        pf = [pf; val(iPf,:)]; %#ok<*AGROW>
        ps = [ps; pop(iPf,:)];
    end
end
end