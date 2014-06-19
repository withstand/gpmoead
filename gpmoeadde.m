function mopForMoeadde = gpmoeadde(mop)
%% Gaussian Process Model based MOEA/D-DE
closeThreshold = 1e-6;

%% handle the structure of mop
ticID = tic;

mopForMoeadde = mop;

%% Input
%   1) npop:    population size npop
%   2) nnbr:    neighbour size nnbr
%   3) nvar:    parameter size nvar
%   4) nobj:    objective size nobj
%   5) domain:  [lower_bound; upper_bound]
%   6) func:    objective function - return [val, std] pair for individual
%   7) fitnessStruct: fitness function with different arguments

%%
% Retrieve data from mop and store them locally
%

%%
% parameter to solve the problem
npop = mop.npop;
% ngen = mop.ngen;
nnbr = mop.nnbr;


%%
% parameter for the problem
nobj = mop.nobj;
nvar = mop.nvar;
domain = mop.domain;
L = domain(1,:);
S = domain(2,:) - L;
func = mop.func;

% data fields that might not exist
constraintFunc = getfieldwithdefault(mop,'constraint',[]);
nbrs = getfieldwithdefault(mop,'nbrs',[]);
weightArray = getfieldwithdefault(mop, 'weightArray',[]);

initialSize = getfieldwithdefault(mop, 'initialSize', 11 * nvar ...
    -1);
evaluatedDesignArray = getfieldwithdefault(mop,'evaluatedDesignArray',[]);
evaluatedObjective = getfieldwithdefault(mop, 'evaluatedObjective', []);
%% Initialize evaluatedSet
if isempty(evaluatedDesignArray)
    evaluatedDesignArray = lhsdesign(initialSize, nvar) .* ...
        repmat(S, initialSize, 1) + repmat(L, initialSize, 1);
end

if ~isempty(constraintFunc)
    for iEvaluated = 1:initialSize
        while 1
            if constraintFunc(evaluatedDesignArray(iEvaluated,:))
                break;
            end
            newDesign = rand(1, nvar) .* S + L;
            evaluatedDesignArray(iEvaluated,:) = newDesign;
        end
    end
end


%% Evaluate evaluatedSet
if isempty(evaluatedObjective)
    evaluatedObjective = func(evaluatedDesignArray);
end


if isempty(weightArray)
    [weightArray, actualLength] = wsdesign(npop, nobj);
else
    actualLength = size(weightArray,1);
end
if actualLength ~= npop
    npop = actualLength;
    mop.npop = npop;
end

if isempty(nbrs)
    nbrs = zeros(npop, nnbr);
    for wi = 1:npop
        D = pdist2(weightArray(wi,:), weightArray, 'euclidean');
        [~,I] = sort(D);
        nbrs(wi,:) = I(1:nnbr);
    end
end

%% Stop criteria

mopForMoeadde.isclose = @isclose;

mopForMoeadde.weightArray = weightArray;



mopForMoeadde.nbrs = nbrs;
mopForMoeadde.fitness = inf;

% while max(mopForMoeadde.fitness) > 0.1
for iGeneration = 1:20
    eeID = tic;
    % buildCGPmodel
    cgp = cell(1,nobj);
    funcs = cell(1, nobj);
    for iObj = 1:nobj
        cgp{iObj} = clusteredgptrain(evaluatedDesignArray, ...
            evaluatedObjective(:,iObj));
        funcs{iObj} = cgp{iObj}.gp_func;
    end
    % buildsurrogatedmop
    
    
    %     mopForMoeadde.pop = rand(npop, nvar) .* repmat(S, npop,1) + ...
    %         repmat(L,npop, 1);
    %
    if rand < 0.0
        mopForMoeadde.fitnessStruct.name = 'negte';
    else
        mopForMoeadde.fitnessStruct.name = 'ei';
    end
    mopForMoeadde.pop = [];
    mopForMoeadde.func = mofunction(funcs);
    %     mopForMoeadde.fitnessStruct.yRef = evaluatedObjective;
    idealObjective = min(evaluatedObjective,[], 1);
    gRef = zeros(npop,1);
    for iPop = 1:npop
        gRefArray = zeros(size(evaluatedObjective,1),1);
        for jEval = 1: size(evaluatedObjective,1)
            gRefArray(jEval) = max(weightArray(iPop,:) .* ...
                (evaluatedObjective(jEval,:) - idealObjective));
        end
        gRef(iPop) = min(gRefArray);
    end
    mopForMoeadde.fitnessStruct.idealObjective = idealObjective;
    mopForMoeadde.fitnessStruct.gRef = gRef;
    
    %% iterat moeadde
    
    % init population if needed
    % evaluated popution if needed
    mopForMoeadde = moeaddeforgp(mopForMoeadde);
    % locate candidates
    
    candidateArray = selectCandidate(mopForMoeadde.pop, ...
        10,...%         mopForMoeadde.npop,...%10,...
        mopForMoeadde.fitness);
    % evaluate and added to evaluatedSet
    candidateObjective = func(candidateArray);
    if size(candidateArray,1)>0
        figure(1)
        plot(evaluatedObjective(:,1), evaluatedObjective(:,2),'gx', ...
            candidateObjective(:,1), candidateObjective(:,2),'ro',...
            idealObjective(1), idealObjective(2),'ks')
        figure(2)
        plot(evaluatedDesignArray(:,1), evaluatedDesignArray(:,2),'gx',...
            candidateArray(:,1), candidateArray(:,2), 'ro');
        drawnow
    end
    
    
    
    
    evaluatedDesignArray = [evaluatedDesignArray; ...
        candidateArray]; %#ok<*AGROW>
    evaluatedObjective = [evaluatedObjective; ...
        candidateObjective];
    
    fprintf('Iteration of exploration and exploit, time %d, allocated candidates:%d.\n',...
        toc(eeID), size(candidateArray,1));
    %     cmdStr = ['save result_', datestr(now,30), '.mat', ...
    %         ' mopForMoeadde evaluatedDesignArray evaluatedObjective'];
    %     eval(cmdStr);
    %
    
end

mopForMoeadde.time = toc(ticID);




    function candidate = selectCandidate(mopPop, candidateSize, fitnessArray)
        
        candidate = [];
        n = size(mopPop,1);
        %         [~, index] = sort(fitnessArray,1,'descend');
        pd = pdist2(mopPop, [mopPop; evaluatedDesignArray]);
        for i=1:n
            pd(i,i) = inf;
        end
        %         pd(1:n, 1:n) = inf;
        distance = min(pd,[],2);
        [~, index] = sort(distance,1, 'descend');
        for i = 1:n
            sel = index(i);
            %             dist1 = pdist2(mopPop(i,:), evaluatedDesignArray);
            %             if isempty(candidate)
            %                 dist2 = inf;
            %             else
            %                 dist2 = pdist2(mopPop(i,:), candidate);
            %             end
            %             if min(dist2) > closeThreshold
            if distance(sel) > closeThreshold
                candidate = [candidate; mopPop(sel,:)];
            end
            %             end
            if size(candidate,1) >= candidateSize
                return;
            end
        end
        
        %         function dist = distance(vec, array)
        %             dist = min(pdist2(vec, array, 'euclidean'));
        %         end
    end



    function [close, closestDistance] = isclose(designPoint)
        closestDistance = min(pdist2(designPoint, evaluatedDesignArray, 'euclidean'));
        close = closestDistance < closeThreshold;
    end


end

