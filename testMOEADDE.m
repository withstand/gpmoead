




mop = testmop('kno1');
mop.ngen=400;mop.npop = 300;mop.nnbr=10;
mop = moeadde(mop);
figure
subplot 211
plot(mop.PF(:,1),mop.PF(:,2),'x')
subplot 212
plot(mop.PS(:,1),mop.PS(:,2),'x')


mop = testmop('sch2');
mop.ngen=400;mop.npop = 300;mop.nnbr=10;
mop = moeadde(mop);
figure
plot(mop.PF(:,1),mop.PF(:,2),'x')



mop = testmop('zdt1',2);
mop.ngen=400;mop.npop = 300;mop.nnbr=2;
mop = moeadde(mop);
figure
plot(mop.PF(:,1),mop.PF(:,2),'x')

zdt1PF = mop.PF;
save zdt1pf zdt1PF

mop = testmop('zdt2',2);
mop.ngen=400;mop.npop = 300;mop.nnbr=10;
mop = moeadde(mop);
figure
plot(mop.PF(:,1),mop.PF(:,2),'x')
zdt2PF = mop.PF;
save zdt2pf zdt2PF

mop = testmop('kno3');
mop.ngen=100;mop.npop = 300;mop.nnbr=10;
mop = moeadde(mop);
figure
plot3(mop.PF(:,1),mop.PF(:,2),mop.PF(:,3),'x')
