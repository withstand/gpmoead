function [p,P] = matingselection(npop, nbr, size,type)
%% Randomly select size individuals amone neighbour of n or the whole population
% Output
%   p   -   size * 1
% Input
%   nbr -   neighbours
%   npop-   popsize
%   size-   select size individuals
%   type-   selection type
%               1 among neighbours
%               2 among the whole population


if type == 1   %neighbourhood
    assert(size <= numel(nbr), ...
        'Can not select %d individuals in %d neighbours.', size, numel(nbr));
    pc = randperm(numel(nbr));
    p = nbr(pc(1:size))';
    P = nbr;
else
    pc = randperm(npop);
    p = pc(1:size)';
    P = 1:npop;
end

end