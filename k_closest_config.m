function Nq = k_closest_config(q, S, K)
D = [];
for i = 1:size(S,2)
    s = S(:,i);
    D = [D config_dist(q,s)];
end 
[~,I] = mink(D, K);
Nq = S(:,I);
