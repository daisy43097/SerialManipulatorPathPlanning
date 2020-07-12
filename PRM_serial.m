function [path, V, E] = PRM_serial(qI, qG, n, K, B, max, L)
qI = [round(qI(1),4);round(qI(2),4);round(qI(3),4);round(qI(4),4)];
qG = [round(qG(1),4);round(qG(2),4);round(qG(3),4);round(qG(4),4)];
V = [qI];
E = [];
EI = [];
weights = [];
[~,c] = size(B);
for i = 1:c % ordering obstacles into ccw 
    D = B{1,i};
    o = ordering(D);
    B{1,i} = o;
end
while size(V,2) < n+1
    q = [max*rand(1); max*rand(1); max*rand(1); max*rand(1)];
    if isconfig_collision_free(q, B, L)
        V = [V q];
    end 
end 
d = 1e10;
for j = 2:n+1 % find the closest vertex from qI
    if config_dist(qI, V(:,j)) < d
        qnear = V(:,j);
        d = config_dist(qI, V(:,j));
        E = [qI;qnear];
        EI = [find_config_index(qI,V); find_config_index(qnear,V)];
    end 
end 
weights = [config_dist(E(1:4,1), E(5:8,1))];
for k = 1:n
    S = V;
    S(:,k) = [];
    Nq = k_closest_config(V(:,k), S, K);
    for m = 1:K
        if ismember([V(:,i);Nq(:,m)]', E', 'rows') == 0 && isintersect_linkpolygon(V(:,k), Nq(:,m), B, L) == 0
            E = [E [V(:,k);Nq(:,m)]];
            EI = [EI [find_config_index(V(:,k),V); find_config_index(Nq(:,m),V)]];
            weights = [weights [config_dist(V(:,k), Nq(:,m))]];
        end 
    end 
end 
d = 1e10;
for x = 1:n+1 % find the closest vertex from qG
    if config_dist(qG, V(:,x))<d
        qgnear = V(:,x);
        d = config_dist(qG, V(:,x));
    end 
end 
if isintersect_linkpolygon(qgnear, qG, B, L) == 0
    V = [V qG];
    E = [E [qgnear; qG]];
    EI = [EI [find_config_index(qgnear,V);find_config_index(qG,V)]];
    weights = [weights [config_dist(qgnear, qG)]];
else 
    error('Fail! Increase sample size!')
end
G = ones(n+2)*Inf;
for i = 1 : n+2
    G(i,i) = 0;
end
for i = 1:size(EI,2)
    G(EI(1,i),EI(2,i)) = weights(i);
    G(EI(2,i),EI(1,i)) = weights(i);
end 
try
    route = dijkstra(G, 1, n+2);
    path = V(:,route);
catch
    error('Please increase sample size!')
end
