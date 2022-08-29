

Dist_LB = abs(repmat(beta_LB(aa:bb), 1, N-5) - repmat(pp(1:(N-5)), N-5, 1));
Dist_UB = abs(repmat(beta_UB(aa:bb), 1, N-5) - repmat(pp(1:(N-5)), N-5, 1));

[n n1] = size(Dist_LB);

NN = sparse(ones(1,n));
LL = speye(n);
       
ub = ones(n^2,1);
lb = zeros(n^2,1);

b = ones(2*n,1);

A = full(get_BlockDiagonal_L(NN,LL,n, false));
ctype = repmat('C',1,n^2);
       
c = reshape(Dist_LB, [n^2 1]);
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

Beta_out_LB = zeros(n, 1);

for i = 1:n
    
    init = (i-1)*n + 1;
    and = (i-1)*n + n;
    
    index = find(x1(init:and) > 0);
    Beta_out_LB(index) = pp(i);
       
end

%[Beta_in_LB Beta_out_LB pp([Sigma_Opt(1) Sigma_Opt(N:-1:2)])']

%CTV_p(Beta_out_LB)
%CTV(pp, [Sigma_Opt(1) Sigma_Opt(N:-1:2)])

                       
c = reshape(Dist_UB, [n^2 1]);
%initiate cplex problem
cplex_UB = Cplex('MILP_struct');
%objective
cplex_UB.Model.sense = 'minimize';
cplex_UB.Model.obj = c;
cplex_UB.Model.lb = lb;
cplex_UB.Model.ub = ub;
cplex_UB.Model.A = A;
cplex_UB.Model.lhs = b;
cplex_UB.Model.rhs = b;
cplex_UB.Model.ctype = ctype;

cplex_UB.solve();
x1 = cplex_UB.Solution.x;

Beta_out_UB = zeros(n, 1);

for i = 1:n
    
    init = (i-1)*n + 1;
    and = (i-1)*n + n;
    
    index = find(x1(init:and) > 0);
    Beta_out_UB(index) = pp(i);
    
end


     