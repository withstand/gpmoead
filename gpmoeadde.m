function mopForMoeadde = gpmoeadde(mop)
%% Gaussian Process Model based MOEA/D-DE
xThreshold = 1e-2;

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
% fprintf('%10s\t%15s\t%15s\t%15s\t%15s\n',...
%     'Iter', 'Time', 'Evaluated','Candidate','max fitness');
% while max(mopForMoeadde.fitness) > 0.1
for iGeneration = 1:10
    
    if size(evaluatedDesignArray,1) > 50 * nobj;
        break;
    end
    
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
    
%     mopForMoeadde.fitnessStruct.name = 'negte';
        mopForMoeadde.fitnessStruct.name = 'eite';
%         mopForMoeadde.fitnessStruct.name = 'wsei';
    disp(['Call MOEA/D-DE with fitness:' mopForMoeadde.fitnessStruct.name]);
    mopForMoeadde.pop = [];
    mopForMoeadde.func = mofunction(funcs);
    %     mopForMoeadde.fitnessStruct.yRef = evaluatedObjective;
    idealObjective = min(evaluatedObjective,[], 1);
    gRef = zeros(npop,1);
    for iPop = 1:npop
        evaluatedSize = size(evaluatedObjective,1);
        gRef(iPop) = min(max(repmat(weightArray(iPop,:),evaluatedSize ,1) .* ...
            (evaluatedObjective - repmat(idealObjective, evaluatedSize,1)),...
            [],2));
    end
    mopForMoeadde.fitnessStruct.idealObjective = idealObjective;
    mopForMoeadde.fitnessStruct.gRef = gRef;
    

    %% iterat moeadde
    
    % init population if needed
    % evaluated popution if needed
    mopForMoeadde = moeaddeforgp(mopForMoeadde);
    % locate candidates
    
    [candidateArray, candidateFitness] = selectCandidate(mopForMoeadde.pop, ...
        10,...%         mopForMoeadde.npop,...%10,...
        mopForMoeadde.fitness);
    % evaluate and added to evaluatedSet
    candidateObjective = func(candidateArray);
    [candidateObjectiveEst, candidateStdEst] = mopForMoeadde.func(...
        candidateArray);
    
    if size(candidateArray,1)>0
        titleStr = ['evaluated=' num2str(size(evaluatedDesignArray,1)),...
            '; candidate=', num2str(size(candidateArray,1))];
        figure(1)
        plot(evaluatedObjective(:,1), evaluatedObjective(:,2),'gx', ...
            candidateObjective(:,1), candidateObjective(:,2),'ro',...
            idealObjective(1), idealObjective(2),'ks','MarkerSize',10)
        [PF,~] = findparetosolution([evaluatedObjective; candidateObjective],...
            [evaluatedDesignArray; candidateArray]);
        
        hold on

        
        load kno1pf
        plot(kno1PF(:,1), kno1PF(:,2),'rv','LineWidth',2,'MarkerSize',6);
        
        PF = sortrows(PF);
        plot(PF(:,1),PF(:,2),'k','LineWidth',2);        
        
        hold off
        
        title(titleStr);
%         figure(2)
%         plot(evaluatedDesignArray(:,1), evaluatedDesignArray(:,2),'gx',...
%             candidateArray(:,1), candidateArray(:,2), 'ro','MarkerSize',10);
%         title(titleStr);
        drawnow
    end
    
    
    
    
    evaluatedDesignArray = [evaluatedDesignArray; ...
        candidateArray]; %#ok<*AGROW>
    evaluatedObjective = [evaluatedObjective; ...
        candidateObjective];
    
    fprintf('-----------------------------------------------------------------------------\n');
    fprintf('%10s\t%15s\t%15s\t%15s\t%15s\n',...
        'Iter', 'Time', 'Evaluated','Candidate','max fitness');
    fprintf('%10d\t%15.4f\t%15.4f\t%15.4f\t%15.4f\n',...
        iGeneration, toc(eeID), size(evaluatedDesignArray,1),...
        size(candidateArray,1), max(mopForMoeadde.fitness));
    
    disp([candidateArray   candidateObjective candidateObjectiveEst ...
        candidateStdEst candidateFitness]);
    
    %     fprintf('Iteration of exploration and exploit, time %d, allocated candidates:%d.\n',...
    %         toc(eeID), size(candidateArray,1));
    %     cmdStr = ['save result_', datestr(now,30), '.mat', ...
    %         ' mopForMoeadde evaluatedDesignArray evaluatedObjective'];
    %     eval(cmdStr);
    %
    
end
[PF,PS] = findparetosolution(evaluatedObjective,...
            evaluatedDesignArray);
mopForMoeadde.PF = PF;
mopForMoeadde.PS = PS;
mopForMoeadde.evaluated = evaluatedDesignArray;
mopForMoeadde.evaluatedObjective = evaluatedObjective;
mopForMoeadde.time = toc(ticID);




    function [candidate, cf] = selectCandidate(mopPop,...
            candidateSize, fitnessArray)
        
        candidate = [];
        cf = [];
        n = size(mopPop,1);
        
        pd = pdist2(mopPop, [mopPop; evaluatedDesignArray]);
        for i=1:n
            pd(i,i) = inf;
        end
        distance = min(pd,[],2);
        
        %% sort strategy
        [~, index] = sort(fitnessArray,1,'descend');
        
        %% max distance strategy
%         [~, index] = sort(distance,1, 'descend');
        
        %% random strategy
%         index = randperm(n);
        for i = 1:n
            sel = index(i);
            %             dist1 = pdist2(mopPop(i,:), evaluatedDesignArray);
            %             if isempty(candidate)
            %                 dist2 = inf;
            %             else
            %                 dist2 = pdist2(mopPop(i,:), candidate);
            %             end
            %             if min(dist2) > closeThreshold
            if distance(sel) > xThreshold 
                candidate = [candidate; mopPop(sel,:)];
                cf = [cf; fitnessArray(sel)];
            end
            %             end
            if size(candidate,1) >= candidateSize
                return;
            end
        end
    end



    function [close, closestDistance] = isclose(designPoint)
        closestDistance = min(pdist2(designPoint, evaluatedDesignArray, 'euclidean'));
        close = closestDistance <= closeThreshold;
    end


    function [pf, ps] = findparetosolution(val,pop)
        pf = [];
        ps = [];

        nval = size(val,1);
        for iPf=1:nval
            do = val([1:iPf-1 iPf+1:nval],:)-...
                repmat(val(iPf,:), nval-1,1)>=-1e-3;
            if sum(sum(do,2)==0) == 0
                pf = [pf; val(iPf,:)];
                ps = [ps; pop(iPf,:)];
            end
        end
    end

end

