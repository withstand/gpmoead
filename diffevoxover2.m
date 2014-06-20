function child = diffevoxover2(zero,one,two, domain)
%% Produce a child with 3 individuals
% Input
%   zero, one, two  -   1 * nvar
%   domain          -   lx, ux,  2 * nvar
% Output
%   child           -   new child


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameter to control the evaluation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rate = 0.5;

child = zero + rate * (two - one);
idx = child < domain(1,:);
remidy = domain(1,:) + rand * (zero - domain(1,:));
child(idx) = remidy(idx);
idx = child > domain(2,:);
remidy = domain(2,:) - rand * (domain(2,:) - zero);
child(idx) = remidy(idx);

end