function inspect_epsilon_mat(path)
% INSPECT_EPSILON_MAT  Print numeric arrays (vectors & matrices) in a .mat

raw = load(path);
fnTop = fieldnames(raw);
if numel(fnTop)==1 && isstruct(raw.(fnTop{1})), raw = raw.(fnTop{1}); end

disp('Numeric arrays found (name : size):');
list_numeric(raw,"");

end

function list_numeric(S, prefix)
fn = fieldnames(S);
for k = 1:numel(fn)
    name = fn{k};
    if strlength(prefix)>0, full = prefix+"."+name; else, full = string(name); end
    v = S.(name);
    if isstruct(v)
        list_numeric(v, full);
    elseif isnumeric(v)
        sz = size(v);
        fprintf('  %-30s : %s\n', full, mat2str(sz));
    end
end
end