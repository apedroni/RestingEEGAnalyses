% gets fields out of a structure
function X = getStruct(Structure,structchain) 

strch = ['x.' structchain];
X = cell2mat(arrayfun(@(x) eval(strch), Structure,'UniformOutput',0));

end

m.m = [1 2]
m.b = [2 3]

cell2mat(arrayfun(@(x) m.m(2), m,'UniformOutput',0));