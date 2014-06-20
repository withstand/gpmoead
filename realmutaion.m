function ind = realmutaion(x, nvar,domain)
rate = 1.0 / nvar;

etam = 20;

rnvar = rand(1,nvar);
rnd = rand(1,nvar);
mut_idx = rnvar < rate;
left_idx = rnd <= 0.5;
right_idx = rnd > 0.5;

lm_idx = mut_idx & left_idx;
rm_idx = mut_idx & right_idx;
% yl = zeros(size(x));
% yu = ones(size(x));
yl = domain(1,:);
yu = domain(2,:);

d1 = (x - yl) ./ (yu-yl);
d2 = (yu - x) ./ (yu-yl);
mut_pow = 1/(etam + 1);

x1 = x + (yu-yl) .* ((2*rnd+ (1-2*rnd).* (1-d1).^(etam+1)  ) .^ mut_pow - 1);
x2 = x + (yu-yl) .* (1 - (2*(1-rnd) + 2*(rnd-0.5) .* (1-d2) .^(etam+1)) .^ mut_pow);

ind = x;
ind(lm_idx) = x1(lm_idx);
ind(rm_idx) = x2(rm_idx);

ind(ind<yl)=yl(ind<yl);
ind(ind>yu)=yu(ind>yu);
end