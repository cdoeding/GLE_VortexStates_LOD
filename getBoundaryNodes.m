function B = getBoundaryNodes(p,x_a,x_b,condition)
%GETBOUNDARYNODES Summary of this function goes here
%   Detailed explanation goes here

N = size(p,1);
B = ones(N,1);

if strcmp(condition,'Dirichlet')
    for i = 1:N
       if p(i,1) == x_a || p(i,1) == x_b || p(i,2) == x_a || p(i,2) == x_b
           B(i) = 0;
       end
    end
end

B = spdiags(B,0,N,N);

end

