
Dist = abs(repmat(beta(aa:bb), 1, N-5) - repmat(pp(1:(N-5)), N-5, 1));

[n n1] = size(Dist);

NN = sparse(ones(1,n));
LL = speye(n);
       
ub = ones(n^2,1);
lb = zeros(n^2,1);

b = ones(2*n,1);

A = full(get_BlockDiagonal_L(NN,LL,n, false));
ctype = repmat('C',1,n^2);
       
c = reshape(Dist, [n^2 1]);
%initiate cplex problem
cplex = Cplex('MILP_struct');
%objective
cplex.Model.sense = 'minimize';
cplex.Model.obj = c;
cplex.Model.lb = lb;
cplex.Model.ub = ub;
cplex.Model.A = A;
cplex.Model.lhs = b;
cplex.Model.rhs = b;
cplex.Model.ctype = ctype;

cplex.solve();
x1 = cplex.Solution.x;

Beta_out = zeros(n, 1);

for i = 1:n
    
    init = (i-1)*n + 1;
    and = (i-1)*n + n;
    
    index = find(x1(init:and) > 0);
    Beta_out(index) = pp(i);
       
end
