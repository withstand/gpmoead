function s = amendstruct(sOld, sNew)
s = sOld;
names = fieldnames(sNew);
for i = 1:numel(names)
    s.(names{i}) = sNew.(names{i});
end
