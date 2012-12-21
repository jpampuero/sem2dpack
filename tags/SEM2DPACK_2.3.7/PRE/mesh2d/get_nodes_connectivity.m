function conn = get_nodes_connectivity(m)

[Q,nel] = size(m.enod);
np = size(m.coor,2);
for k=1:np
  conn{k} = [];
end
for e=1:nel,
  for k=1:Q,
    conn{m.enod(k,e)}(end+1) = e;
  end
end
