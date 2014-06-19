function mopRet = moeaddeforgp(mop)
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
fitnessStruct = getfieldwithdefault(mop,'fitnessStruct',[]);
constraintFunc = getfieldwithdefault(mop,'constraint',[]);
nbrs = getfieldwithdefault(mop,'nbrs',[]);
pop = getfieldwithdefault(mop,'pop',[]);
weightArray = getfieldwithdefault(mop, 'weightArray',[]);
val = getfieldwithdefault(mop,'val',[]);
std = getfieldwithdefault(mop,'std',[]);


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
    %     pop01 = lhsdesign(npop, nvar);
    pop01 = rand(npop, nvar);
    pop = pop01 .* repmat(S, npop,1) + repmat(L, npop,1);
end
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

% for iEvaluated = 1:npop
%     while 1
%         if ~mop.isclose(pop(iEvaluated,:))
%             break;
%         end
%         newDesign = rand(1, nvar) .* S + L;
%         pop(iEvaluated,:) = newDesign;
%     end
% end


%   4) [val,std]:   initial estimations for each individual
if isempty(val)
    [val, std] = func(pop);
end

%   5) vr:          ideal objectives
% vr = fitnessStruct.vr;

%   6) fitness:     fitness of
fitness = zeros(npop,1);
% fitness = fit_func(val,std, vr);

% fprintf('Initialize problem in %f seconds.\n', toc(ticID));
fprintf('%6s\t%12s\t%12s\t%12s\t%12s\n',...
    '   ', '    ', '   ','Fitness','    ');
fprintf('%6s\t%12s\t%12s\t%12s\t%12s\t%6s\n',...
    'Iter', 'Time', 'max','min','mean','Stall');
stallGeneration = 0;
for iGen = 1:mop.ngen
    ticIteration = tic;
    %% Update
    perm = randperm(npop);
    %
    isStall = 1;
    for i = 1:npop
        
        %how is this work for MOEA/D-DE, not really needed for it, remove it
        %if isempty(mop.fitness_struct.y_ref)
        %    fitnessStruct.y_ref = val;
        %end
        
        n = perm(i);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % parameter to control the evaluation %
        realb = exp(-iGen / (mop.ngen/2));                          %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (rand < realb)
            type = 2;   % whole population
        else
            type = 1;   % neighbourhood
        end
        
        %% select the indexes of mating parents
        [p,P] = matingselection(npop,nbrs(n,:),2,type);
        
        while 1
            %% produce a child solution
            child = diffevoxover2(pop(i,:),pop(p(1),:),...
                pop(p(2),:), domain);
            %% apply polynomial mutation
            child = realmutaion(child, nvar,domain);
%             if ~mop.isclose(child) && ...
            if isempty(constraintFunc) || constraintFunc(child)
                break
            end
        end
        %% evaluate the child solution
        [val_c, std_c] = func(child);
        
        %% update the reference points and other solution in the neighbourhood
        % or the whole population
        %         vr = min([vr; val_c-3*std_c]);
        fitness = evaluatefitness(fitnessStruct.name,weightArray,val,std,...
            fitnessStruct.idealObjective, fitnessStruct.gRef);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % parameter to control the evaluation          %
        nupdate = min(nnbr,max(4, round(nnbr / 5)));   %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sel = randperm(numel(P), nupdate);
        selected = P(sel);
        for si = selected
            fitness_c =evaluatefitness(fitnessStruct.name, ...
                weightArray(si,:),val_c, std_c,...
                fitnessStruct.idealObjective, fitnessStruct.gRef(si));
            if fitness(si) < fitness_c  %maximize fitness
                isStall = 0;
                stallGeneration = 0;
                pop(si,:) = child;
                val(si,:) = val_c;
                std(si,:) = std_c;
                fitness(si) = fitness_c;
            end
        end
    end
    
    if isStall
        stallGeneration = stallGeneration + 1;
    end
    

    fprintf('%6d\t%12.4f\t%12.4f\t%12.4f\t%12.4f\t%6d\n',...        
        iGen, toc(ticIteration), max(fitness),min(fitness), ...
        mean(fitness), stallGeneration);
    
    
    if stallGeneration > 20
        break;
    end
end
mopRet.npop = npop; % could be changed by the weightArray design
mopRet.pop = pop;
mopRet.val = val;
mopRet.std = std;
mopRet.fitness = fitness;
% mopRet.weightArray = weightArray;
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
function fit = evaluatefitness(name,ws, y_est, y_std, idealPoint, gRef)
[sampleSize,objDim] = size(y_est);

switch name
    case 'eite'
        objRet = te(ws, y_est, y_std, idealPoint);
        fit = ei(objRet, gRef);
    case 'negte'        
%         objRet = te(ws, y_est, y_std, idealPoint);
%         fit = -objRet(:,1);         
        fit = -max(ws .* (y_est - repmat(idealPoint, sampleSize, 1)),[],2); 
    case 'wsei'    
        eiY = zeros(sampleSize, objDim);
        for i=1:objDim
            eiY(:,i) = ei([y_est(:,i) y_std(:,i)], idealPoint(i));
        end
        fit = sum(ws .* eiY, 2);     
        % eiY reference zero, no more expected improvement   
    otherwise
end
end


%%
function eiRet = ei(objDistribution, gRef)
u = (gRef - objDistribution(:,1)) ./ objDistribution(:,2);
eiRet = (gRef -objDistribution(:,1)) .* normcdf(u) +...
    objDistribution(:,2) .* normpdf(u);
end


function objDist = te(weight, y_est, y_std, idealpoint)
[sampleSize,objDim]=size(y_est);
% if isempty(y_std)
%     y_std = zeros(size(y_est));
% end
[ns,~] = size(weight);

if ns~=sampleSize
    error('Currently not working for multiple weight vs. multiple y.');
end
% obj_dist = zeros(ns,2);

yhat = weight .* (y_est - repmat(idealpoint, ns, 1));
ystd = weight .* y_std;
objDist = [yhat(:,1) ystd(:,1)];
for jj=2:objDim
    objDist = gvmax(objDist(:,1), objDist(:,2), yhat(:,jj), ystd(:,jj));
end

% for i = 1:ns
%     lambda = weight(i, :);
%     yihat = lambda .* (y_est(i,:) - idealpoint);
%     yistd = lambda .* y_std(i,:);
%     od(1) = yihat(1);
%     od(2) = yistd(1);
%     for jj = 2:objDim
%         od = gvmax(od(1),od(2), yihat(jj),yistd(jj));
%     end
%     obj_dist(i,:) = od;
% end

    function dist = gvmax(mu1, std1, mu2, std2)
        tau = sqrt(std1.^2 + std2.^2);
        if sum(tau) == 0
            dist=[max(mu1, mu2) 0];
        else
            alpha = (mu1 - mu2) ./ tau;
            normcdfAlpha = normcdf(alpha);
            normpdfAlpha = normpdf(alpha);
            dist = [mu1.*normcdfAlpha+mu2.*(1-normcdfAlpha)+...
                tau.*normpdfAlpha  sqrt(...
                (mu1.^2+std1.^2).*normcdfAlpha+...
                (mu2.^2+std2.^2).*(1-normcdfAlpha)+...
                (mu1+mu2).*tau.*normpdfAlpha)];
        end
    end
end




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