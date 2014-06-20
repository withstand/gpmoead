function ah = newfigure2ax
fig = figure;
pos = get(fig,'Position');
pos(1) = max(50,pos(1)/4);
pos(3) = pos(3) * 2;
set(fig,'Position',pos);
ah(1)=subplot(1,2,1);
ah(2)=subplot(1,2,2);