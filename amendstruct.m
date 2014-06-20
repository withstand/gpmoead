function s = amendstruct(sOld, sNew)
% Fix fields in the first struct with fields of the second struct.
% s = amendstruct(sOld, sNew)
% Add all new fields in sNew to sOld
%   change all exist fields in sOld with sNew value
% Return a new struct.
%
% For example
%   s1 = struct('a',1,'b',2);
%   s2 = struct('a',12,'c',22);
%   s3 = amendstruct(s1,s2);
%   % s3 == struct('a',12,'b',2,'c',22);

s = sOld;
names = fieldnames(sNew);
for i = 1:numel(names)
    s.(names{i}) = sNew.(names{i});
end
