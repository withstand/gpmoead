function comparetwomop(mop1, mop2, samplePoint)
assert(mop1.nvar==2 && mop2.nvar==2);

domain =  mop1.domain;
L = domain(1,:);
U = domain(2,:);
S = U - L;
dS = S / 100;
xgv = L(1):dS(1):U(1);
ygv = L(2):dS(2):U(2);


[x,y] = meshgrid(xgv, ygv);
nn = numel(x);
snn = sqrt(nn);


z1 = mop1.func([reshape(x,nn,1) reshape(y,nn,1)]);
z2 = mop2.func([reshape(x,nn,1) reshape(y,nn,1)]);
for iObj = 1:mop1.nobj
    figure(iObj+32767);
%     surf(x,y,reshape(z1(:,iObj), snn,snn));
%     hold on
%     surf(x,y,reshape(z2(:,iObj), snn,snn));              
%     hold off
    surf(x,y,reshape(z1(:,iObj)-z2(:,iObj), snn,snn));
    zmin = min(z1(:,iObj)-z2(:,iObj));
    if nargin>=3
        hold on
        plot3(samplePoint(:,1), samplePoint(:,2), 1.5*zmin*ones(size(samplePoint,1),1),'ko');
        hold off
    end
%     legend({[' mop1'],...
%         [' mop2']});
    title(['objective ', num2str(iObj)]);
    
    
    
end
drawnow

% disp('Compare the mop.....');
% pause

end