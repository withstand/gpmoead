mop = testmop('kno1');
mop.ngen=100;mop.npop = 300;mop.nnbr=10;
mop = moeadde(mop);
figure
subplot 211
plot(mop.PF(:,1),mop.PF(:,2),'x')
subplot 212
plot(mop.PF(:,1),mop.PF(:,2),'x')


mop = testmop('sch2');
mop.ngen=100;mop.npop = 300;mop.nnbr=10;
mop = moeadde(mop);
figure
plot(mop.PF(:,1),mop.PF(:,2),'x')



mop = testmop('zdt1',10);
mop.ngen=100;mop.npop = 300;mop.nnbr=10;
mop = moeadde(mop);
figure
plot(mop.PF(:,1),mop.PF(:,2),'x')



mop = testmop('zdt2',10);
mop.ngen=100;mop.npop = 300;mop.nnbr=10;
mop = moeadde(mop);
figure
plot(mop.PF(:,1),mop.PF(:,2),'x')


mop = testmop('kno3');
mop.ngen=100;mop.npop = 300;mop.nnbr=10;
mop = moeadde(mop);
figure
plot3(mop.PF(:,1),mop.PF(:,2),mop.PF(:,3),'x')
