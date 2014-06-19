function mopRet = moeadde(mop)
mopRet = mop;
%% Input
%   1) npop:    population size npop
%   2) nnbr:    neighbour size nnbr
%   3) nvar:    parameter size nvar
%   4) nobj:    objective size nobj
%   5) domain:  [lower_bound; upper_bound]
%   6) func:    objective function - return [val, std] pair for individual
%   7) fitness: structure to select fitness function with different arguments


%%
% Retrieve data from mop and store them locally
%
npop = mop.npop;
nobj = mop.nobj;
nvar = mop.nvar;
nnbr = mop.nnbr;
domain = mop.domain;
L = domain(1,:);
S = domain(2,:) - L;
func = mop.func;

% data fields that might not exist
% fitnessStruct = getfieldwithdefault(mop,'fitnessStruct',[]);
constraintFunc = getfieldwithdefault(mop,'constraint',[]);
nbrs = getfieldwithdefault(mop,'nbrs',[]);
pop = getfieldwithdefault(mop,'pop',[]);
weightArray = getfieldwithdefault(mop, 'weightArray',[]);
val = getfieldwithdefault(mop,'val',[]);


%% Initialize
ticID = tic;

%   1) weightArray:         the weight vectors
if isempty(weightArray)
    [weightArray, actualLength] = wsdesign(npop,nobj);
else
    actualLength = size(weightArray,1);
end

% Could change npop slightly according to the way of designing weightArray
if actualLength ~= npop
    npop = actualLength;
end


%   2) nbrs:           the neighours for each individual
if isempty(nbrs)
    nbrs = zeros(npop, nnbr);
    for wi = 1:npop
        D = pdist2(weightArray(wi,:), weightArray, 'euclidean');
        [~,I] = sort(D);
        nbrs(wi,:) = I(1:nnbr);
    end
end


%   3) pop:         initial population
if isempty(pop)
    pop01 = lhsdesign(npop, nvar);
    %     pop01 = rand(npop, nvar);
    pop = pop01 .* repmat(S, npop,1) + repmat(L, npop,1);
    
    % regenerate pop for it
    if ~isempty(constraintFunc)
        for iPop = 1:npop
            while 1
                if constraintFunc(pop(iPop,:))
                    break;
                end
                newPop = rand(1, nvar) .* S + L;
                pop(iPop,:) = newPop;
            end
        end
    end
end

%   4) [val,std]:   initial estimations for each individual
if isempty(val)
    val = func(pop);
end

%   5) vr:          ideal objectives
vr = min(val,[],1);

%   6) fitness:     fitness of
fitness = zeros(npop,1);
% fitness = fit_func(val,std, vr);

% fprintf('Initialize problem in %f seconds.\n', toc(ticID));
fprintf('%6s\t%12s\t%12s\t%12s\t%12s\n',...
    '   ', '    ', '   ','Fitness','    ');
fprintf('%6s\t%12s\t%12s\t%12s\t%12s\n',...
    'Iter', 'Time', 'max','min','mean');

for iGen = 1:mop.ngen
    ticIteration = tic;
    %% Update
    perm = randperm(npop);
    %
    for i = 1:npop
        
        %how is this work for MOEA/D-DE, not really needed for it, remove it
        %if isempty(mop.fitness_struct.y_ref)
        %    fitnessStruct.y_ref = val;
        %end
        
        n = perm(i);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % parameter to control the evaluation %
        realb = 0.8;                          %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        realb = 1 - exp(-iGen / (mop.ngen/10));  
        if (rand < realb)
            type = 1;  % neighbourhood
        else
            type = 2;  % whole population
        end
        
        %% select the indexes of mating parents
        [p,P] = matingselection(npop,nbrs(n,:),2,type);
        
        while 1
            %% produce a child solution
            child = diffevoxover2(pop(i,:),pop(p(1),:),...
                pop(p(2),:), domain);
            %% apply polynomial mutation
            child = realmutaion(child, nvar,domain);
            if isempty(constraintFunc) || constraintFunc(child)
                break
            end
        end
        %% evaluate the child solution
        val_c = func(child);
        
        %% update the reference points and other solution in the neighbourhood
        % or the whole population
        vr = min([vr; val_c]);
        fitness = evaluatefitness(weightArray,val,vr);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % parameter to control the evaluation          %
        nupdate = min(nnbr,max(5, round(nnbr / 5)));   %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sel = randperm(numel(P), nupdate);
        selected = P(sel);
        for si = selected
            fitness_c =evaluatefitness(weightArray(si,:),val_c, vr);
            if fitness(si) < fitness_c  %maximize fitness
                pop(si,:) = child;
                val(si,:) = val_c;
                fitness(si) = fitness_c;
            end
        end
    end
    fprintf('%6d\t%12.4f\t%12.4f\t%12.4f\t%12.4f\n',...
        iGen, toc(ticIteration), max(fitness),min(fitness), mean(fitness));
end
mopRet.npop = npop; % could be changed by the weightArray design
mopRet.pop = pop;
mopRet.val = val;
mopRet.fitness = fitness;
mopRet.weightArray = weightArray;
mopRet.nbrs = nbrs;

% It's arguable whethere to obtain the solutions or produced a Pareto solution from evaluated points
% PF = [];
% PS = [];
% for i=1:npop
%     do = val([1:i-1 i+1:npop],:)-repmat(val(i,:), npop-1,1)>=-1e-3;
%     if sum(sum(do,2)==0) == 0
%         PF = [PF; val(i,:)];
%         PS = [PS; pop(i,:)];
%     end
% end
PF = val;
PS = pop;

mopRet.PF = PF;
mopRet.PS = PS;

end


%% Evaluate the fitness function with a fitnessStructure
% For MOEA/D-DE, no need for ei, just negte for negative Tchebycheff decomposition
function fit = evaluatefitness( ws, y_est, idealpoint)
ns = size(y_est,1);
obj_dist = max(ws .* (y_est - repmat(idealpoint,ns,1)),[],2);
fit = -obj_dist;
end


% % %%
% % function ei_ret = ei(weight, y_est, y_std, idealpoint, g_ref)
% % obj_dist = te(weight, y_est, y_std, idealpoint);
% % u = (g_ref - obj_dist(:,1)) ./ obj_dist(:,2);
% % ei_ret = (g_ref - obj_dist(:,1)) .* normcdf(u) + obj_dist(:,2) .* normpdf(u);
% % end
% 
% 
% function obj_dist = te(weight, y_est, idealpoint)
% [sampleSize,~]=size(y_est);
% % if isempty(y_std)
% %     y_std = zeros(size(y_est));
% % end
% [ns,~] = size(weight);
% 
% if ns~=sampleSize
%     %     if ns == 1
%     %         weight = repmat(weight,sampleSize,1);
%     %         ns = sampleSize;
%     %     elseif sampleSize ==1
%     %         y_est = repmat(y_est, ns, 1);
%     % %         y_std = repmat(y_std, ns ,1);
%     %         %         sampleSize = ns;
%     %     else
%     error('Currently not working for multiple weight vs. multiple y.');
%     %     end
% end
% 
% 
% % obj_dist = zeros(ns,1);
% obj_dist = max(weight .* (y_est - repmat(idealpoint,ns,1)),[],2);
% 
% 
% % for i = 1:ns
% %     lambda = weight(i, :);
% %     yihat = lambda .* (y_est(i,:) - idealpoint);
% %
% % % %     yistd = lambda .* y_std(i,:);
% % %     od = yihat(1);
% % % %     od(2) = yistd(1);
% % %     for jj = 2:objDim
% % %         od = max(od(1), yihat(jj));
% % %     end
% %     obj_dist(i) = max(yihat);
% % end
% 
% 
% %     function dist = gvmax(mu1, std1, mu2, std2)
% %         tau = sqrt(std1.^2 + std2.^2);
% %         if sum(tau) == 0
% %             dist=[max(mu1, mu2) 0];
% %         else
% %             alpha = (mu1 - mu2) ./ tau;
% %             dist = [mu1.*normcdf(alpha)+mu2.*normcdf(-alpha)+...
% %                 tau.*normpdf(alpha)  sqrt(...
% %                 (mu1.^2+std1.^2).*normcdf(alpha)+...
% %                 (mu2.^2+std2.^2).*normcdf(-alpha)+...
% %                 (mu1+mu2).*tau.*normpdf(alpha))];
% %         end
% %     end
% end




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