% mop = testmop('zdt2',2);
clear all
mops = {testmop('kno1'), testmop('kno3'), testmop('sch2'),...
    testmop('zdt1',10),testmop('zdt2',10)};
for i=1:numel(mops)
    % mop = testmop('kno3',2);
    mops{i}.ngen = 500;
    mops{i}.npop = 400;
    mops{i}.nnbr = 2;
    mops{i} = gpmoeadde(mops{i});
    figure;plot(mops{i}.PF(:,1),mopsdd810412{i}.PF(:,2),'go')
end