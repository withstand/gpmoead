function [x, actualLength] = wsdesign(m,n)

% 
nn=2;
while nchoosek(nn+n-1,n-1) < m
    nn = nn + 1;
end

x = dist_m_n(nn,n) / nn;

actualLength = nchoosek(nn+n-1,n-1);
if actualLength ~= m
    fprintf('Request weights length %d, actual weights length %d.\n', m, nchoosek(nn+n-1,n-1));
end


end

%% 
function dis = dist_m_n(m,n)
if n==1
    dis = m;
elseif m==0
    dis = zeros(1,n);
else
    dis = [];
    for i = 0:m
        disx = dist_m_n(m-i,n-1);
        kk = size(disx,1);
        dist = [repmat(i, kk,1) disx];
        dis = [dis;dist];
    end
end

end