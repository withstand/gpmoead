function mopForMoeadde = gpmoeadde(mop)
%% Gaussian Process Model based MOEA/D-DE
xThreshold = 1e-3;

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

% mopForMoeadde.isclose = @isclose;

mopForMoeadde.weightArray = weightArray;
mopForMoeadde.nbrs = nbrs;
mopForMoeadde.fitness = inf;

isFirstBuild = 1;
cgp = cell(1,nobj);
funcs = cell(1, nobj);
% fprintf('%10s\t%15s\t%15s\t%15s\t%15s\n',...
%     'Iter', 'Time', 'Evaluated','Candidate','max fitness');
% while max(mopForMoeadde.fitness) > 0.1
for iGeneration = 1:1000
    
    if size(evaluatedDesignArray,1) > 50 * nvar;
        break;
    end
    
    eeID = tic;
    % buildCGPmodel
    disp('Building meta-model.....');
    if isFirstBuild
        for iObj = 1:nobj
            cgp{iObj} = clusteredgptrain(evaluatedDesignArray, ...
                evaluatedObjective(:,iObj),[],...
                struct('algorithm','de','refresh',0));
            funcs{iObj} = cgp{iObj}.gp_func;
        end
        isFirstBuild = 0;
    else
        for iObj = 1:nobj
            cgp{iObj} = clusteredgptrain(evaluatedDesignArray, ...
                evaluatedObjective(:,iObj),[],...
                struct('algorithm','pd','refresh',0), cgp{iObj});
            funcs{iObj} = cgp{iObj}.gp_func;
        end
    end
    % buildsurrogatedmop
    
    mopForMoeadde.fitnessStruct.name = 'negte';
    %     mopForMoeadde.fitnessStruct.name = 'eite';
    %             mopForMoeadde.fitnessStruct.name = 'wsei';
    mopForMoeadde.pop = [];
    mopForMoeadde.func = mofunction(funcs);
    %     mopForMoeadde.evaluated = evaluatedDesignArray;
    %     mopForMoeadde.evaluatedObjective = evaluatedObjective;
    
%     showmop(newfigure2ax,mopForMoeadde);
    
    %     pause
    
    %     comparetwomop(mop, mopForMoeadde,evaluatedDesignArray);
    %     mopForMoeadde.fitnessStruct.yRef = evaluatedObjective;
%     idealObjective = min(evaluatedObjective,[], 1);
%     gRef = zeros(npop,1);
%     for iPop = 1:npop
%         evaluatedSize = size(evaluatedObjective,1);
%         gRef(iPop) = min(max(repmat(weightArray(iPop,:),evaluatedSize ,1) .* ...
%             (evaluatedObjective - repmat(idealObjective, evaluatedSize,1)),...
%             [],2));
%     end
%     mopForMoeadde.fitnessStruct.idealObjective = idealObjective;
%     mopForMoeadde.fitnessStruct.gRef = gRef;
    
    
    %% iterat moeadde
    disp(['Call MOEA/D-DE with fitness:' mopForMoeadde.fitnessStruct.name]);
    
    % init population if needed
    % evaluated popution if needed
    %     mopFake = moeadde(mopForMoeadde);
    %     mopForMoeadde = moeaddeforgp(mopForMoeadde);
    mopForMoeadde = moeadde(mopForMoeadde);
    
    % locate candidates
    
    candidateArray = selectCandidate(mopForMoeadde, 10, evaluatedDesignArray, xThreshold);
    %,...%         mopForMoeadde.npop,...%10,...
    %         mopForMoeadde.fitness);
    % evaluate and added to evaluatedSet
    candidateObjective = func(candidateArray);
    %     [candidateObjectiveEst, candidateStdEst] = mopForMoeadde.func(...
    %         candidateArray);
    
    %     if size(candidateArray,1)>0
    %         titleStr = ['evaluated=' num2str(size(evaluatedDesignArray,1)),...
    %             '; candidate=', num2str(size(candidateArray,1))];
    %         figure(1)
    %         clf
    %         plot(candidateObjective(:,1), candidateObjective(:,2),'rp',...
    %             'MarkerSize',15);
    %
    %         idealObjective = min([evaluatedObjective; candidateObjective],[], 1);
    %         hold on;
    %         plot(idealObjective(1), idealObjective(2),'ks','MarkerSize',10)
    %         hold off
    %
    %
    % %         [PF,~] = findparetosolution([evaluatedObjective; candidateObjective],...
    % %             [evaluatedDesignArray; candidateArray]);
    %
    %
    %
    %         if exist([lower(mop.name), 'pf.mat'],'file')
    %             ret = [];
    %             name = lower(mop.name);
    %             commandStr = ['ret = load(''', name, 'pf'');'];
    %             %         load kno1pf
    %             eval(commandStr);
    %
    %             commandStr = ['plot(ret.',...
    %                 name, 'PF(:,1), ret.',...
    %                 name, 'PF(:,2),''rv'',''LineWidth'',2,''MarkerSize'',2);'];
    %             hold on
    %             eval(commandStr);
    %             hold off
    %         end
    %
    %
    % %         hold on
    % %         PF = sortrows(PF);
    % %         plot(PF(:,1),PF(:,2),'k*','MarkerSize',12);
    % %         hold off
    %
    %         %         figure(2)
    %         %         plot(evaluatedDesignArray(:,1), evaluatedDesignArray(:,2),'gx',...
    %         %             candidateArray(:,1), candidateArray(:,2), 'ro','MarkerSize',10);
    %         %         title(titleStr);
    %
    %         hold on
    %         trueValueOfPs = mop.func(mopForMoeadde.PS);
    %         plot(mopForMoeadde.PF(:,1),mopForMoeadde.PF(:,2),'go',...
    %             trueValueOfPs(:,1), trueValueOfPs(:,2), 'g+')
    %         hold off
    %
    % %         hold on
    % %         trueValueOfPs = mop.func(mopFake.PS);
    % %         plot(mopFake.PF(:,1),mopFake.PF(:,2),'co',...
    % %             trueValueOfPs(:,1), trueValueOfPs(:,2), 'c+')
    % %         hold off
    %         title(titleStr);
    %
    %         drawnow
    %     end
    
    
    
    
    evaluatedDesignArray = [evaluatedDesignArray; ...
        candidateArray]; %#ok<*AGROW>
    evaluatedObjective = [evaluatedObjective; ...
        candidateObjective];
    
    fprintf('-----------------------------------------------------------------------------\n');
    fprintf('%10s\t%15s\t%15s\t%15s\t%15s\n',...
        'Iter', 'Time', 'Evaluated','Candidate','max fitness');
    fprintf('%10d\t%15.4f\t%15.4f\t%15d\t%15.4f\n',...
        iGeneration, toc(eeID), size(evaluatedDesignArray,1),...
        size(candidateArray,1), max(mopForMoeadde.fitness));
    
    %     disp([candidateArray   candidateObjective candidateObjectiveEst ...
    %         candidateStdEst candidateFitness]);
    
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







%     function [close, closestDistance] = isclose(designPoint)
%         closestDistance = min(pdist2(designPoint, evaluatedDesignArray, 'euclidean'));
%         close = closestDistance <= closeThreshold;
%     end




% end

