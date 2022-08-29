function [LB_improved] = Solve_LB_improved(F_tilde_LB, F_tilde_UB, LB_NN, LB, UB, theta, N)

    c = [F_tilde_LB, F_tilde_UB, LB_NN];
    cc = [LB, UB, theta];
    d = [(6/(N-6))*(LB^2), (6/(N-6))*(UB^2), (6/(N-6))*(theta^2)];
    dd = c + d;
    
    cvx_clear
    cvx_begin
    variable x(3,1); 
    maximize ( dd*x - (6/(N-6))*((cc*x)*(cc*x)) ) ;
    subject to% 
        sum(x) == 1;
        x >= 0;
    cvx_end

    %[F_tilde_LB F_tilde_UB LB_NN LB UB theta]
    
    LB_improved = dd*x - (6/(N-6))*(cc*x)^2;
    
    %x
end