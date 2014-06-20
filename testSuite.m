% mop = testmop('zdt2',2);
clear all
mops = {testmop('kno1'), testmop('kno3'), testmop('sch2'),...
    testmop('zdt1',2),testmop('zdt2',2)};
for i=[1, 4 5]
    % mop = testmop('kno3',2);
    mops{i}.ngen = 50;
    mops{i}.npop = 200;
    mops{i}.nnbr = 2;
    mops{i} = gpmoeadde(mops{i});
    figure;plot(mops{i}.PF(:,1),mops{i}.PF(:,2),'go');
    drawnow;
end