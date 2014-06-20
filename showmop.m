function showmop(axesHandles, mop)

domain =  mop.domain;
L = domain(1,:);
U = domain(2,:);
S = U - L;
dS = S / 100;
xgv = L(1):dS(1):U(1);
ygv = L(2):dS(2):U(2);


[x,y] = meshgrid(xgv, ygv);
nn = numel(x);
snn = sqrt(nn);


z = mop.func([reshape(x,nn,1) reshape(y,nn,1)]);

samplePoint = getfieldwithdefault(mop,'evaluated',[]);
% sampleZ = mop.func(samplePoint);
sampleZ = getfieldwithdefault(mop,'evaluatedObjective',[]);


for iObj = 1:mop.nobj
    %     surf(x,y,reshape(z1(:,iObj), snn,snn));
    %     hold on
    %     surf(x,y,reshape(z2(:,iObj), snn,snn));
    %     hold off
    surf(axesHandles(iObj),x,y,reshape(z(:,iObj), snn,snn));
    if ~isempty(samplePoint) && ~isempty(sampleZ)
        hold(axesHandles(iObj),'on');
        plot3(axesHandles(iObj), samplePoint(:,1), ...
            samplePoint(:,2), ...
            sampleZ(:,iObj),'ko','MarkerFaceColor','k');
        hold(axesHandles(iObj),'off');
    end
    %     legend({[' mop1'],...
    %         [' mop2']});
    title(axesHandles(iObj),['objective ', num2str(iObj)]);
end
drawnow

% disp('Compare the mop.....');
% pause

end