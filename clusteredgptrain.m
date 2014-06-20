function clusteredGPModel = clusteredgptrain(data, y, oldCgp)
% Train a clustered Gaussian Process model with multi-dimensional data-output
%
%   cgp = clusteredgptrain(data, y, oldCgp)
%   Input:
%       data:   design points, n * nvar
%       y:      objective values, n * nobj
%       oldCgp: old cgp, used to transfer options or update an old model
%               with new data
%   Output:
%       cgp:    trained model
%
if nargin < 3 % all options set to default, using de alogrithm to start
    oldCgp = [];
end
% if nargin<4
%     gpOpts = [];
% end
% if nargin < 3
%     clusteringOpts = [];
% end
clusteredGPModel = oldCgp;

%% data properties
[nx, nvar] = size(data);
[ny, nobj] = size(y);
assert(nx == ny);
n = nx;



%% clustering options
L1 = getfieldwithdefault(clusteredGPModel,'L1',80);
L2 = getfieldwithdefault(clusteredGPModel,'L2',20);
if n <= L1
    clusteringSize = 1;
else
    clusteringSize = 1 + ceil((n-L1)/L2);
end
% clusteringSize = max(clusteringSize, 1);
% build a cgm for each objective

cgpm = cell(1, nobj);
%% Build Clustering
fprintf('Building C-fuzzy clustered Gaussian Process models.\n');
fprintf('%10s\t%10s\t%10s\n',...
    'Objective', 'Model No.', 'Time(s)');
for iObj = 1:nobj
    currentCGPM.models = cell(1,clusteringSize);
    if clusteringSize > 1
        [centers,U,~] = fcm(data,clusteringSize,[2.0,100,1e-5,0]);
        % obtain the indexs show the relation to different cluster
        [~,idx] = sort(U,2,'descend');
        currentCGPM.centers = centers;
        currentCGPM.idx = idx(:,1:L1);
    else
        currentCGPM.idx = 1:n;
    end
    
    for i = 1:clusteringSize
        idxi = currentCGPM.idx(i,:);
        tic;
        if isempty(oldCgp) || ...
                 isempty(getfieldwithdefault(oldCgp,'cgpm',[])) || ...
                 numel(oldCgp.cgpm) ~= nobj || ...
                 isempty(getfieldwithdefault(oldCgp.cgpm{iObj},'models',[])) || ...
                 numel(oldCgp.cgpm{iObj}.models) < i
            currentCGPM.models{i} = gptrain(data(idxi,:),y(idxi, iObj),...
                amendstruct(clusteredGPModel,...
                struct('algorithm','de')));
        else
            currentCGPM.models{i} = gptrain(data(idxi,:),y(idxi, iObj),...
                amendstruct(clusteredGPModel.cgpm{iObj}.models{i}, ...
                struct('algorithm','pd')));
        end
%         fprintf('Building model %d, time used %f.\n', i, toc);
        fprintf('%10d\t%10d\t%10.2f\n',...
             iObj, i,  toc);
    end
    cgpm{iObj} = currentCGPM;
end  % iObj for each objective

clusteredGPModel.cgpm = cgpm;
clusteredGPModel.data = data;
clusteredGPModel.y = y;
clusteredGPModel.clusteringSize = clusteringSize;
clusteredGPModel.L1 = L1;
clusteredGPModel.L2 = L2;
clusteredGPModel.datai = @getcluster;
clusteredGPModel.yi = @gety;
clusteredGPModel.func = @func;
clusteredGPModel.nobj = nobj;
clusteredGPModel.nvar = nvar;
clusteredGPModel.domain = getfieldwithdefault(clusteredGPModel,...
    'domain', [min(data); max(data)]);

    function yi = gety(ci)
        ci_index  = clusteredGPModel.cgpm{1}.idx(ci,:);
        yi = clusteredGPModel.y(ci_index,:);
    end


    function datai = getcluster(ci)
        ci_index  = clusteredGPModel.cgpm{1}.idx(ci,1:clusteredGPModel.L1);
        datai = clusteredGPModel.data(ci_index,:);
    end

    function [hatY, hatStd] = func(xStar)
        n = size(xStar,1);
        hatY = zeros(n, nobj);
        hatStd = zeros(n, nobj);
        for iobj = 1:nobj
            yStar = clusteredGPEst(clusteredGPModel.cgpm{iobj}, xStar);
            hatY(:,iobj) = yStar(:,1);
            hatStd(:, iobj) = yStar(:,2);
        end
    end



    function yStar = clusteredGPEst(cgp, xStar)
        [m,~] = size(xStar);
        yStar = zeros(m,2);
        parClusteringSize = clusteredGPModel.clusteringSize;
        parCenter = getfieldwithdefault(cgp,'centers',[]);
        closest = zeros(m,1);
        for ii=1:m  % m points to estimate
            % find the nearest centers to test
            if  parClusteringSize > 1
                d = parCenter - repmat(xStar(ii,:), parClusteringSize,1);
                d = diag(d * d');
                closestArray = find(d == min(d));
                % In case there are more cloest center points, just pick the first one.
                closest(ii) = closestArray(1);
            else
                closest(ii) = 1;
            end
            yStar(ii,:) = cgp.models{closest(ii)}.gp_func(xStar(ii,:));
        end
        
    end
end


