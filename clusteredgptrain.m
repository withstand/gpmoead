function cgp = clusteredgptrain(data, y, clusteringOpts, gpOpts, oldCgp)
if nargin<5
    oldCgp = [];
end
if nargin<4
    gpOpts = [];
end
if nargin < 3
    clusteringOpts = [];
end


%% data properties
[sampleSize, ~] = size(data);
assert(sampleSize==numel(y));

%% clustering options
L1 = getfieldwithdefault(clusteringOpts,'L1',80);
L2 = getfieldwithdefault(clusteringOpts,'L2',20);
cgp.clusteringOpts.L1 = L1;
cgp.clusteringOpts.L2 = L2;
if sampleSize <= L1
    clusteringSize = 1;
else
    clusteringSize = 1 + ceil((sampleSize-L1)/L2);
end
% clusteringSize = max(clusteringSize, 1);
cgp.clusteringOpts.clusteringSize = clusteringSize;

%% gp parameters
cgp.gpOpts = gpOpts;

%% Build Clustering
cgp.models = cell(clusteringSize,1);
if clusteringSize > 1
    [centers,U,~] = fcm(data,clusteringSize,[2.0,100,1e-5,0]);
    % obtain the indexs show the relation to different cluster
    [~,idx] = sort(U,2,'descend');
    cgp.centers = centers;
    cgp.idx = idx(:,1:L1);
else
    cgp.idx = 1:sampleSize;
end
for i = 1:clusteringSize
    idxi = cgp.idx(i,:);
    tic;
    if isempty(oldCgp)
        cgp.models{i} = gptrain(data(idxi,:),y(idxi), cgp.gpOpts);
    else
        if i > numel(oldCgp.models)
             cgp.models{i} = gptrain(data(idxi,:),y(idxi),...
                amendstruct(gpOpts,...
                struct('algorithm','de')));           
        else            
            cgp.models{i} = gptrain(data(idxi,:),y(idxi),...
                amendstruct(oldCgp.models{i}.gpOpts, gpOpts));
        end
    end
    fprintf('Building model %d, time used %f.\n', i, toc);
end


cgp.gp_func = @clusteredGPEst;
cgp.datai = @get_cluster;
cgp.yi = @get_y;
cgp.func = @func; 

    function yi = get_y(ci)
        ci_index  = cgp.idx(ci,1:cgp.clusteringOpts.L1);
        yi = y(ci_index);
    end


    function datai = get_cluster(ci)
        ci_index  = cgp.idx(ci,1:cgp.clusteringOpts.L1);
        datai = data(ci_index,:);
    end

    function [hat_y, hat_std] = func(xStar)
        yStar = clusteredGPEst(xStar);
        hat_y = yStar(:,1);
        hat_std = yStar(:,2);
    end
    function yStar = clusteredGPEst(xStar)
        
        [m,~] = size(xStar);
        
        yStar = zeros(m,2);
        
        parClusteringSize = cgp.clusteringOpts.clusteringSize;
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
%        
% %        parModels = [cgp.models{closest}];
% 
%        for ii = 1:m
%             % gpt = gp_predict(cgp.models{closest}, xStar(ii,:));
%         end
        
    end
end


