function mopRet = gpmoeadde(mop)
%% Gaussian Process Model based MOEA/D-DE
    
    
%% handle the structure of mop
    ticID = tic;
    
    mopRet = mop;
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
    initialSize = getfieldwithdefault(mop, 'initialSize', 11 * nvar ...
                                           -1);
    evaluatedSet = getfieldwithdefault(mop,'evaluatedSet',[]);
    
    %% Initialize evaluatedSet
    evaluatedDesignArray = lhsdesign(initialSize, nvar) .* ...
        repmat(S, initialSize, 1) + repmat(L, initialSize, 1);

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
    evaluatedObjective = func(pop);
    
    
    if isempty(weightArray)
        [weightArray, actualLenght] = wsdesign(npop, nobj);            
    else
        actualLength = size(weightArray,1);
    end
    if actualLength ~= nop
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
    mopForMoeadde = mop;
    mopForMoeadde.fitnessStruct.name = 'ei';
    mopForMoeadde.weightArray = weightArray;
    mopForMoeadde.pop = lhsdesign(npop, nvar) .* repmat(S, npop,1) + ...
        repmat(L,npop, 1);
    mopForMoeadde.nbrs = nbrs;
    while 1
        
        % buildCGPmodel
        cgp = clusteredgptrain(evaluatedDesignArray, evaluatedObjective);
        
        % buildsurrogatedmop
        mopForMoeadde.func = cgp.gp_func;
        mopForMoeadde.fitnessSturct.yRef = evaluatedObjective;        

        
        
        %% iterat moeadde

        % init population if needed
        % evaluated popution if needed
        mopForMoeadde = moeadde(mopForMoeadde);        
        % locate candidates
        
        candidateArray = selectCandidate(mopForMoeadde.pop, 10);
        % evaluate and added to evaluatedSet        
        candidateObjective = func(candidateArray);
        evaluatedDesignArray = [evaluatedDesignArray; ...
                            candidateArray];
        evaluatedObjective = [evaluatedObjective; ...
                            candidateObjective];        

    end



end
