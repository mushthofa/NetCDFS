function list = trimlist(list, maxn)
    idxempty = cellfun(@(x) isempty(x), list);
    list(idxempty) = [];
    if(length(list)>maxn)
        list = list(1:maxn);
    end
end