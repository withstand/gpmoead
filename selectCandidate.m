function [candidate] = selectcandidate(mop,...
    candidateSize, evaluatedDesignArray,xThreshold)

candidate = [];
%         cf = [];
n = mop.npop;

pd = pdist2(mop.pop, [mop.pop; evaluatedDesignArray]);
for i=1:n
    pd(i,i) = inf;
end
distance = min(pd,[],2);

%% sort strategy
%         [~, index] = sort(fitnessArray,1,'descend');
kk = round(n/ candidateSize);

%% max distance strategy
%         [~, index] = sort(distance,1, 'descend');

%% random strategy
%         index = randperm(n);

for i = 1:candidateSize
    if i==candidateSize
        toSelect = (i-1)*kk+1:n;
    else
        toSelect = (i-1)*kk+1:1:i*kk;
    end
    
    %             dist1 = pdist2(mopPop(i,:), evaluatedDesignArray);
    %             if isempty(candidate)
    %                 dist2 = inf;
    %             else
    %                 dist2 = pdist2(mopPop(i,:), candidate);
    %             end
    %             if min(dist2) > closeThreshold
    for j = randperm(numel(toSelect))
        sel = toSelect(j);
        if distance(sel) >= xThreshold
            candidate = [candidate; mop.pop(sel,:)]; %#ok<*AGROW>
            break;
            %                     cf = [cf; mop.fitness(sel)];
        end
        %             end
        %                 if size(candidate,1) >= candidateSize
        %                     return;
        %                 end
    end
end
end
