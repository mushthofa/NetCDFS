function [comp, num, subgraph] = subconncomp(g, nodes)

n = sum(nodes);
subgraph = g(nodes, nodes);
subgraph = (subgraph+subgraph') > 0;

comp = makeSet(n);
visited = zeros(n, 1);
for nn = 1:n
    if(~visited(nn))
        q = nn;
        while ~isempty(q)
            node = q(1);
            visited(node) = 1;
            adj = find(subgraph(node, :));
            for i=1:length(adj)
                comp = unionSet(comp, node, adj(i));
                if(~visited(adj(i)))
                    q = [q; adj(i)];
                end
            end
            q = q(2:end);
        end
    end
end

comp = flattenSet(comp);
num = length(unique(comp));

end



function par = makeSet(n)
    par = 1:n;
end

function [p, par] = findSet(par, i)
    if(par(i) ~= i)
        par(i) = findSet(par, par(i));
    end
    p = par(i);
end

function par = unionSet(par, i, j)
    [ir, par] = findSet(par, i);
    [jr, par] = findSet(par, j);
    par(ir) = jr;
end

function par = flattenSet(par)
    for i=1:length(par)    
        [~, par] = findSet(par, i);
    end
end
    

  