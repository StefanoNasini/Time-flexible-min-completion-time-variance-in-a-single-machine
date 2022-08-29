%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve all the collection of instances
%--------------------------------------------------------------------------
% To run the experiment execute from the matlab terminal the command MAIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fclose('all');
close all;
clear all;
clc;

format longG;

delete('.\output\*.*');

addpath('C:\Users\s.nasini\Dropbox\cvx')
addpath('C:\Users\s.nasini\Dropbox\Rabia+Stefano\Matching')

addpath ('C:\Program Files\IBM\ILOG\CPLEX_Studio129\cplex\matlab\x64_win64')
addpath ('C:\Program Files\IBM\ILOG\CPLEX_Studio129\cplex\examples\src\matlab')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Arguments:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Solve_bound:      If 0, the theta^LB is solved. If 1, the theta^UB is solved;
% VI_type           An integer from 1 to 6, denoting the type of VI to used
%                   (see list of VI types below);
% Mean_bound:       If 1, the mean_bound is solved and appended as a VI to 
%                   the problem. If 0, the problem is solved directly;
% N:                Number of jobs;
%
%--------------------------------------------------------------------------
% List of VI types
%--------------------------------------------------------------------------
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p_bar_case = 1;

%--------------------------------------------------------------------------
% Many jobs experiment
%--------------------------------------------------------------------------

N_set = [1000];
L_set = [100];
S_set = [2];
coef = [5];
GenType = [0 1];
Replication = 1;

NBS = length(coef);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build data files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid1 = fopen('.\output\SummaryTable.log', 'a+');
fprintf(fid1,'Instance N GenType L S coef replication P(N-1) p_bar_gap CVT_AS CTV_opt LB_vs LB_nn LB_impr_nonconvex F_C_hat_nop1 F_C_hat PasteBeta_C_hat k_heu k_opt_C_hat L_p1 U_p1 L_Mann U_Mann mean(dist_C_hat) max(dist_C_hat) min(dist_C_hat) CPU_QP_C_hat CPU_heur \n');

fid2 = fopen('.\output\Beta_sol_nop1.log', 'a+');
fid3 = fopen('.\output\Beta_sol.log', 'a+');

Instance = 1;

%coef = [0.6 0.7 0.8 0.9 1 1.25 1.5 1.75 2];

for N = N_set,

    outputstr = [];
    for i = 1:N
        outputstr = [outputstr ['%-8.2f']];
    end

    for type = GenType
        
        for ii = 1:Replication,
               
            PP = PP_generator(N, coef, ii, type);

            for jj = 1:NBS

                pp = PP{jj};

                %pp = [10 13 18 22 63 104 145 186 227 268 309 390 391 432 486 3514 3634 6980];
                
                % pp = [1 10 17 18 22 63 104 145 186 227 268 309 390 391 432 473 3000 3514 5831 5832];
                %
                % for k = 1:(N-1)
                %     pp(k) = min(pp(k+1), round(pp(k) + k*rand()) );
                % end  
                
                % Case 0;
                % pp = [10 13 18 22 63 104 145 186 227 268 309 390 391 432 473 3514 4831 6980];

                % Case 1
                % pp = [10 13 18 22 63 104 145 186 227 268 309 390 391 432 473 3514 6000 6980];
                
                %N = 12;
                %pp = [7 24 47 58 69 76 80 104 106 120 121 566]; 
                % Case 0 (Either seq1 = (n, n-2, n-3, ..., n-4, n-1) or seq2 = (n, n-2, n-4, ..., n-3, n-1))
                %pp = [7 24 47 58 69 76 80 104 106 120 269 566]; % Case 2 (seq1 = (n, n-2, n-3, ..., n-4, n-1))
                %pp = [7 24 47 58 69 76 80 104 106 120 566 566]; % Case 1 (seq3 = (n, n-2, n-3, n-4, ... , n-1))
                
                p_bar = pp(N-2) + (4*pp(N-3) + sum(pp(1:(N-4))))/2;
                
                beta = sum(pp(1:(N-7)));
                
                G = min(pp(2:N) - pp(1:(N-1)));
                
                P_3_4 = pp(N-3) - pp(N-4);
                P_5_6 = pp(N-5) - pp(N-6);
                
                f1 = ((7*pp(N-7)+beta)/(7*(3*pp(N-7)+beta)))*(14*((3*pp(N-7)+ beta)/(7*pp(N-7)+beta) - ((N-4)/N))*P_3_4  - 2*P_5_6 + 2*((2*N-7)/N)*beta -4*beta*((3*pp(N-5)+4*pp(N-6)+beta)/(7*pp(N-7)+ beta)) + 14*((N-7)/N)*pp(N-7)) ;
                f2 = ((7*pp(N-7)+beta)/(7*(3*pp(N-7)+beta)))*(14*((3*pp(N-7)+ beta)/(7*pp(N-7)+beta) - ((N-4)/N))*P_3_4  - 14*((pp(N-7)+beta)/(7*pp(N-7)+ beta))*P_5_6 - (14*beta/N) + 14*((N-7)/N)*pp(N-7)) ;
                f3 = ((7*pp(N-6)+beta)/(7*(3*pp(N-6)+beta)))*(8*((beta)/(7*pp(N-6)+beta) - ((N-7)/N))*P_3_4  + 7*((pp(N-6)+7*beta)/(7*pp(N-6)+ beta))*pp(N-5) + 2*((6*N-49)/N)*pp(N-6) - (68*beta/7) + 2*(((beta^2)*(34*N - 49))/(7*N*(7*pp(N-6) + beta))) ) ;
                f4 =  2*((N-6)/N)*(3*P_3_4 - 2*P_5_6);
                
                p_bar_bar = min([pp(N-2) + max([f1, f2, f3, f4]), p_bar]);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Select case to overwrite p(N-1)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if p_bar_case == 0
                    
                    rand_p_n_1 = rand();
                    pp(N-1) = floor(rand_p_n_1*pp(N-2) + (1-rand_p_n_1)*p_bar_bar);
                    
                    pp(N) = pp(N-1) + 1;
                    
                end
                if p_bar_case == 2
                    
                    rand_p_n_1 = rand();
                    pp(N-1) = ceil(rand_p_n_1*p_bar_bar + (1-rand_p_n_1)*p_bar);
                    
                    pp(N) = pp(N-1) + 1;
                    
                end
                if p_bar_case == 1

                    rand_p_n_1 = rand();
                    
                    pp(N-1) = 1 + ceil(p_bar);
                    pp(N) = pp(N-1) + 1;
                    
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Common parameters
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                MS = sum(pp);

                if abs(floor(N/2) - N/2) > 0
                    Sigma_AS = [N:(-2):1 2:(2):(N-1)];
                else
                    Sigma_AS = [N:(-2):1 1:(2):(N-1)];
                end

                tic
                Sigma_Opt = H_NC(pp, Sigma_AS);
                CPU_time_heur = toc;
                
                pp_opt = pp([Sigma_Opt(1) Sigma_Opt(N:-1:2)]);

                CTV_AS = CTV(pp, Sigma_AS);
                CTV_opt = CTV(pp, Sigma_Opt);
                LB_nn = LB_basic(pp, 0);
                LB_vs = LB_VS(pp);
                
                m = (N-6)/2;

                if abs(floor(N/2) - N/2) > 0
                    
                    LB0 = 0;
                    for i = 1:m
                        LB0 = LB0 + (sum(pp(1:2*i)))^2;
                    end
 
                else
                    
                    LB0 = 0;
                    for i = 1:m
                        LB0 = LB0 + (sum(pp(1:2*i-1)))^2;
                    end
                    
                end

                LB0 = LB0/2;
                
                k_heu = find(Sigma_Opt == 1);
                
                u_Mann = zeros(N-5, 1);
                v_Mann = zeros(N-5, 1);
                pp_decr = pp(N:-1:1);
                
                for kk = 4:1:(N-2)
                      
                    kkk = kk-3;
                    
                    u_Mann(kkk) = pp_decr(3) + sum((2:1:(kk-2)).*pp_decr(4:kk)) + (kk-1)*pp_decr(N) - pp_decr(2) -  sum((2:1:(N-kk)).*pp_decr((kk+1):1:(N-1))) + ((N-1)/2)*(pp_decr(N-kk+2) - pp_decr(N));
   
                end
                for kk = 5:1:(N-1)
                      
                    kkk = kk-4;
                    
                    v_Mann(kkk) = -pp_decr(2) - sum((2:1:(N-kk)).*pp_decr(4:N-kk+2)) - (N-kk+1)*pp_decr(N) + pp_decr(3) +  sum((2:1:(kk-2)).*pp_decr((N-kk+3):1:(N-1))) + ((N-1)/2)*(pp_decr(kk) - pp_decr(N));
                       
                end
                
                
                LL_Mann = 3 + max(find(u_Mann <= 0));
                UU_Mann = 4 + min(find(v_Mann >= 0));
                    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Solve the beta_problems
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if p_bar_case == 0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    AA = zeros(N-5, 1);
                    BB = zeros(N-5, 1);
                    uu_1 = zeros(N-5, 1);
                    uu_2 = zeros(N-5, 1);
                    vv_1 = zeros(N-5, 1);
                    vv_2 = zeros(N-5, 1);
                    
                    for k = 4:(N-2)
                        
                        kk = k - 3;
                        
                        AA(kk) = pp(N-2) - pp(N-1) + (k-1)*pp(1) + ((N-1)/2)*(pp(N-k) - pp(1));
                        BB(kk) = pp(N-1) - pp(N-2) + ((N-2)/2)*(pp(N-k) - pp(1));
                        
                        uu_1(kk) = AA(kk) + 2*pp(N-3) - 2*pp(N-4) + sum(((1:(k-4))+2).*pp(N-(1:(k-4))-4)) - sum((N - k - (1:(N-k-2)) - 1).*pp((1:(N-k-2))+1));
                        uu_2(kk) = AA(kk) + 2*pp(N-4) - 2*pp(N-3) + sum(((1:(k-4))+2).*pp(N-(1:(k-4))-4)) - sum((N - k - (1:(N-k-2)) - 1).*pp((1:(N-k-2))+1));
                        
                        vv_1(kk) = BB(kk) - 2*pp(N-3) + 2*pp(N-4) + (N-k+1)*pp(1) + sum( (N - k - (1:(N-k-2)) + 1).*pp(N - (1:(N-k-2)) - 4) ) - sum(((1:(k-4))+2).*pp(N-(1:(k-4))-4));
                        vv_2(kk) = BB(kk) - 2*pp(N-4) + 2*pp(N-3) + (N-k+1)*pp(1) + sum( (N - k - (1:(N-k-2)) + 1).*pp(N - (1:(N-k-2)) - 4) ) - sum(((1:(k-4))+2).*pp(N-(1:(k-4))-4));
                        
                    end
                    
                    LL_1 = 5 + max(find(uu_1 <= 0));
                    LL_2 = 5 + max(find(uu_2 <= 0));
                    
                    UU_1 = 3 + max(find(vv_1 >= 0));
                    UU_2 = 3 + max(find(vv_2 >= 0));

                    %LL_1 = 4;
                    %UU_1 = 10;
                        
                    %LL_2 = 4;
                    %UU_2 = 10;
                    
                    %[uu_1; vv_1]
                    
                    B = 5*((N-5)/N)*(pp(N-5))^2 + 2*pp(N-5)*(5*pp(N) + 4*pp(N-2) + 3*pp(N-4)+ 2*pp(N-3) + pp(N-1)+ ((3*N-5)/N)*sum(pp(1:(N-5))) )  + sum(pp(1:(N-5)))*(2*pp(N) + 2*pp(N-2) + 2*pp(N-4) + ((N-1)/N)*sum(pp(1:(N-5)) )) ;  
                    A = (1/(2*N))*(2*N*pp(N) + 2*pp(N-1) + 2*(N-1)*pp(N-2) + 2*(N-2)*pp(N-3) + (N-1)*pp(N-4) + (N-1)*sum(pp(1:(N-5))));
                    D = (4*pp(N) + pp(N-1) + 3*pp(N-2) + 8*pp(N-3)/N + 4*(N-2)*pp(N-4)/N + 2*sum(pp(1:(N-5)))  )/4;

                    LB_1 = (MS+pp(N))/2 - ((N-2)/(2*N))*(pp(N-1) - pp(N-2)) ;
                    LB_2 = LB_1;                        

                    UB_1 = A ;
                    UB_2 = B/(2*(sum(pp(1:(N-5))) + 5*pp(N-5)));

                    theta_1 = (3*pp(N) + 2*pp(N-2) + pp(N-3) + 3*MS - 2*pp(N-1) - pp(N-4))/6;
                    theta_2 = (3*pp(N) + 2*pp(N-2) + pp(N-4) + 3*MS - 2*pp(N-1) - pp(N-3))/6;

                    mean_p = mean(pp(1:(N-5)));   
                   
                    for S = S_set,
                          
                        Non_convex_beta_bound = 0;
                        
                        %----------------------------------------------
                        % Solve the beta-problem using the quadratic formulation
                        %----------------------------------------------
                        
                        G_s = zeros(S,1);
                        for s = 1:S
                            G_s(s) = min(pp(2:N) - pp(1:(N-1)));
                        end
                
                        L = 2*(N-6);

                        %----------------------------------------------
                        % CALL AMPL
                        %----------------------------------------------

                        Call_AMPL_one;
                          
                        %----------------------------------------------
                                                
                        Incumbent_1 = 1.0e+20;
                        Incumbent_2 = 1.0e+20;
                        
                        k_opt_1 = 0;
                        k_opt_2 = 0;
                        
                        aa = 1;
                        bb = (N-5);
                        
                        if min(F_1_nop1) < min(F_1_nop1)
                            
                            k_opt_beta_F = find(F_1_nop1 == min(F_1_nop1)) + LL_1-1;

                            k_opt_beta = find(beta_1_nop1 == min(beta_1_nop1)) + 3;
                            
                            F_C_hat_nop1 = min(F_1_nop1);
                            
                            k = k_opt_beta;
                            
                            beta = beta_1;
                                                            
                            Beta_round_C_hat;
                            
                            Beta_out = [pp(N); pp(N-2); pp(N-3); Beta_out; pp(N-4); pp(N-1)];
                            beta_print = [pp(N); pp(N-2); pp(N-3); beta; pp(N-4); pp(N-1)];
                            beta_1_nop1_print = [pp(N); pp(N-2); pp(N-3); beta_1_nop1; pp(N-4); pp(N-1)];
                            
                            CTV_k = CTV_p(Beta_out);
                                       
                            F_C_hat_1 = F_1;
                            Incumbent_1 = CTV_k;
                            beta_incumbent = Beta_out;
                                
                            Dist_1 = abs(beta(aa:bb)' - pp_opt(4:(N-2)))./pp_opt(4:(N-2));     
                            
                            k_opt_1 = k;

                        else
                            
                            k_opt_beta_F = find(F_2_nop1 == min(F_2_nop1)) + LL_2-1;

                            k_opt_beta = find(beta_2_nop1 == min(beta_2_nop1)) + 3;
                            
                            F_C_hat_nop1 = min(F_2_nop1);
                                              
                            beta = beta_2;
                            
                            k = k_opt_beta;
                                                        
                            Beta_round_C_hat;

                            Beta_out = [pp(N); pp(N-2); pp(N-4); Beta_out; pp(N-3); pp(N-1)];
                            beta_print = [pp(N); pp(N-2); pp(N-4); beta; pp(N-3); pp(N-1)];
                            beta_1_nop1_print = [pp(N); pp(N-2); pp(N-4); beta_1_nop1; pp(N-3); pp(N-1)];
                            
                            CTV_k = CTV_p(Beta_out);
                            
                            F_C_hat_2 = F_2;
                            Incumbent_2 = CTV_k;
                            beta_incumbent = Beta_out;
                                
                            Dist_2 = abs(beta(aa:bb)' - pp_opt(4:(N-2)))./pp_opt(4:(N-2));                        

                            k_opt_2 = k;
                                
                           
                        end
                        
                        if Incumbent_1 < Incumbent_2

                            F_C_hat = F_C_hat_1;
                            Dist_C_hat = Dist_1;
                            k_opt_C_hat = k_opt_1;
                            PasteBeta_C_hat = Incumbent_1;

                        else
                            
                            F_C_hat = F_C_hat_2;
                            Dist_C_hat = Dist_2;
                            k_opt_C_hat = k_opt_2;
                            PasteBeta_C_hat = Incumbent_2;
                        end
                        
                        fprintf(fid2, [outputstr ' \r\n'], beta_print);                        
                        fprintf(fid3, [outputstr ' \r\n'], beta_1_nop1_print);
                        %fprintf(fid1,'\r\n');
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Print the CPLEX output
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        fprintf(fid1,'%-8d %-8d %-8d %-8d %-8d %-8.2f %-8d %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f \r\n', Instance, N, type, L, S, coef(jj), ii, pp(N-1), p_bar, CTV_AS, CTV_opt, LB_vs, LB_nn, Non_convex_beta_bound, F_C_hat_nop1, F_C_hat, PasteBeta_C_hat, k_heu, k_opt_C_hat, LL_1, UU_1, LL_Mann, UU_Mann, mean(Dist_C_hat), max(Dist_C_hat), min(Dist_C_hat), CPU_time_QP_C_hat, CPU_time_heur);
                        %fprintf(fid1,'\r\n');
                            
                        Instance = Instance + 1;  
                    
                    end % S
                    
                end %case

                if p_bar_case == 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    AA = zeros(N-5, 1);
                    BB = zeros(N-5, 1);
                    uu_1 = zeros(N-5, 1);
                    vv_1 = zeros(N-5, 1);

                    for k = 4:(N-2)
                        
                        kk = k - 3;
                        
                        AA(kk) = pp(N-2) - pp(N-1) + (k-1)*pp(1) + ((N-1)/2)*(pp(N-k) - pp(1));
                        BB(kk) = pp(N-1) - pp(N-2) + ((N-2)/2)*(pp(N-k) - pp(1));
                        
                        uu_1(kk) = AA(kk) + 2*pp(N-3) - 2*pp(N-4) + sum(((1:(k-4))+2).*pp(N-(1:(k-4))-4)) - sum((N - k - (1:(N-k-2)) - 1).*pp((1:(N-k-2))+1));
                        vv_1(kk) = BB(kk) - 2*pp(N-3) + 2*pp(N-4) + (N-k+1)*pp(1) + sum( (N - k - (1:(N-k-2)) + 1).*pp(N - (1:(N-k-2)) - 4) ) - sum(((1:(k-4))+2).*pp(N-(1:(k-4))-4));

                    end
                    
                    LL_1 = 5 + max(find(uu_1 <= 0));
                    LL_2 = LL_1;
                    
                    UU_1 = 3 + max(find(vv_1 >= 0));
                    UU_2 = UU_1;
                    
                    %LL_1 = 4;
                    %UU_1 = 10;
                        
                    %LL_2 = 4;
                    %UU_2 = 10;
                    
                    %[uu_1; vv_1]
                    
                    B = 5*((N-5)/N)*(pp(N-5))^2 + 2*pp(N-5)*(5*pp(N) + 4*pp(N-2) + 3*pp(N-4)+ 2*pp(N-3) + pp(N-1)+ ((3*N-5)/N)*sum(pp(1:(N-5))) )  + sum(pp(1:(N-5)))*(2*pp(N) + 2*pp(N-2) + 2*pp(N-4) + ((N-1)/N)*sum(pp(1:(N-5)) )) ;  
                    A = (1/(2*N))*(2*N*pp(N) + 2*pp(N-1) + 2*(N-1)*pp(N-2) + 2*(N-2)*pp(N-3) + (N-1)*pp(N-4) + (N-1)*sum(pp(1:(N-5))));
                    D = (4*pp(N) + pp(N-1) + 3*pp(N-2) + 8*pp(N-3)/N + 4*(N-2)*pp(N-4)/N + 2*sum(pp(1:(N-5)))  )/4;

                    LB_1 = (MS+pp(N))/2 - ((N-2)/(2*N))*(pp(N-1) - pp(N-2)) ;                      
                    UB_1 = A ;
                    theta_1 = (3*pp(N) + 2*pp(N-2) + pp(N-3) + 3*MS - 2*pp(N-1) - pp(N-4))/6;

                    mean_p = mean(pp(1:(N-5)));   
                   
                    for S = S_set,
                         
                        Non_convex_beta_bound = 0;
                        
                        %----------------------------------------------
                        % Solve the beta-problem using the quadratic formulation
                        %---------------------------h-------------------
                        
                        L = 2*(N-6);

                        %----------------------------------------------
                        % CALL AMPL
                        %----------------------------------------------

                        Call_AMPL_one;
                          
                        %----------------------------------------------
                        
                        aa = 1;
                        bb = (N-5);

                        k_opt_beta = find(beta_nop1 == min(beta_nop1)) + 3;
                        
                        F_C_hat_nop1 = min(F_nop1);
                        k = k_opt_beta;
                        Beta_round_C_hat;
                        
                        Beta_out = [pp(N); pp(N-2); pp(N-3); Beta_out; pp(N-4); pp(N-1)];
                        beta_print = [pp(N); pp(N-2); pp(N-3); beta; pp(N-4); pp(N-1)];
                        beta_1_nop1_print = [pp(N); pp(N-2); pp(N-3); beta_nop1; pp(N-4); pp(N-1)];
                            
                        CTV_k = CTV_p(Beta_out);
                        
                        F_C_hat_1 = F;
                        beta_incumbent = Beta_out;
                                
                        Dist_1 = abs(beta(aa:bb)' - pp_opt(4:(N-2)))./pp_opt(4:(N-2));     
                            
                        k_opt_1 = k;

                        F_C_hat = F_C_hat_1;
                        Dist_C_hat = Dist_1;
                        k_opt_C_hat = k_opt_1;
                        PasteBeta_C_hat = CTV_k;

                        fprintf(fid2, [outputstr ' \r\n'], beta_print);                        
                        fprintf(fid3, [outputstr ' \r\n'], beta_1_nop1_print);
                        %fprintf(fid1,'\r\n');
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Print the CPLEX output
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        fprintf(fid1,'%-8d %-8d %-8d %-8d %-8d %-8.2f %-8d %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f \r\n', Instance, N, type, L, S, coef(jj), ii, pp(N-1), p_bar, CTV_AS, CTV_opt, LB_vs, LB_nn, Non_convex_beta_bound, F_C_hat_nop1, F_C_hat, PasteBeta_C_hat, k_heu, k_opt_C_hat, LL_1, UU_1, LL_Mann, UU_Mann, mean(Dist_C_hat), max(Dist_C_hat), min(Dist_C_hat), CPU_time_QP_C_hat, CPU_time_heur);
                        %fprintf(fid1,'\r\n');
                            
                        Instance = Instance + 1;  
                        
                    end % S
                                        
                end %case
                
                if p_bar_case == 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    AA = zeros(N-5, 1);
                    BB = zeros(N-5, 1);
                    uu_1 = zeros(N-5, 1);
                    vv_1 = zeros(N-5, 1);
                    
                    for k = 4:(N-2)
                        
                        kk = k - 3;
                        
                        AA(kk) = pp(N-2) - pp(N-1) + (k-1)*pp(1) + ((N-1)/2)*(pp(N-k) - pp(1));
                        BB(kk) = pp(N-1) - pp(N-2) + ((N-2)/2)*(pp(N-k) - pp(1));
                        
                        uu_1(kk) = AA(kk) + 2*pp(N-3) - 3*pp(N-4) + sum(((1:(k-5))+3).*pp(N-(1:(k-5))-4)) - sum((N - k - (1:(N-k-1)) + 1).*pp((1:(N-k-1))+1));
                        vv_1(kk) = BB(kk) - 2*pp(N-3) - 3*pp(N-4) + (N-k)*pp(1) + sum( (N - k - (1:(N-k-2)) ).*pp(N - (1:(N-k-2)) - 4) ) - sum(((1:(k-4))+3).*pp(N-(1:(k-4))-4));

                    end
                    
                    LL_1 = 4 + max(find(uu_1 <= 0));
                    UU_1 = 5 + max(find(vv_1 >= 0));
                                        
                    B = 5*((N-5)/N)*(pp(1))^2 + 2*pp(1)*(5*pp(N) + 4*pp(N-2) + 3*pp(N-3)+ 2*pp(N-4) + pp(N-1)+ ((2*N-5)/N)*sum(pp(1:(N-5))) ) +  sum(pp(1:(N-5)))*(2*pp(N) + 2*pp(N-2) + 2*pp(N-3) + 2*pp(N-4) + ((N-1)/N)*sum(pp(1:(N-5)) )) ;
                    A = 5*((N-5)/N)*(pp(N-5))^2 + 2*pp(N-5)*(5*pp(N) + 4*pp(N-2) + 3*pp(N-3)+ 2*pp(N-4) + pp(N-1)+ ((2*N-5)/N)*sum(pp(1:(N-5))) ) +  sum(pp(1:(N-5)))*(2*pp(N) + 2*pp(N-2) + 2*pp(N-3) + 2*pp(N-4) + ((N-1)/N)*sum(pp(1:(N-5)) )) ;  
                    D1 = -5*((N-5)/N)*(pp(1))^2 + 2*pp(1)*(5*pp(N) + 4*pp(N-2) + 3*pp(N-3)+ 2*pp(N-4) + pp(N-1)+ ((N+5)/N)*sum(pp(1:(N-5))) ) +  sum(pp(1:(N-5)))*(2*pp(N) + 2*pp(N-2) + 2*pp(N-3) + 2*pp(N-4) + ((N+1)/N)*sum(pp(1:(N-5)) )) ; 
                    D2 = -5*((N-5)/N)*(pp(N-5))^2 + 2*pp(N-5)*(5*pp(N) + 4*pp(N-2) + 3*pp(N-3)+ 2*pp(N-4) + pp(N-1)+ ((N+5)/N)*sum(pp(1:(N-5))) ) +  sum(pp(1:(N-5)))*(2*pp(N) + 2*pp(N-2) + 2*pp(N-3) + 2*pp(N-4) + ((N+1)/N)*sum(pp(1:(N-5)) )) ; 

                    T1 = (1/(2*N))*( 2*N*pp(N) + 2*pp(N-1) + 2*(N-1)*pp(N-2) + 2*(N-2)*pp(N-3) + (N-1)*sum(pp(1:(N-4)))  );
                    T2 = min(D1/(2*(sum(pp(1:(N-5))) + 5*pp(N-5))),D2/(2*(sum(pp(1:(N-5))) + 5*pp(1))));

                    %LB = (MS+pp(N))/2 - ((N-2)/(2*N))*(pp(N-1) - pp(N-2)) ;
                    LB = max(T1, T2);

                    UB = max(A/(2*(sum(pp(1:(N-5))) + 5*pp(N-5))), B/(2*(sum(pp(1:(N-5))) + 5*pp(1)))) ;

                    theta_6 = (4*pp(N) + 3*pp(N-2) + 2*pp(N-3) + pp(N-4)+ 2*MS - pp(N-1))/6;

                    mean_p = mean(pp(1:(N-5)));  
                    
                    for S = S_set,
                        
                        Non_convex_beta_bound = 0;
                        
                        %----------------------------------------------
                        % Solve the beta-problem using the quadratic formulation
                        %----------------------------------------------
                        
                        L = 2*(N-6);
                        
                        %----------------------------------------------
                        % CALL AMPL
                        %----------------------------------------------

                        Call_AMPL_one;
                          
                        %----------------------------------------------
                        
                        aa = 1;
                        bb = (N-5);

                        k_opt_beta = find(beta_nop1 == min(beta_nop1)) + 3;
                        
                        F_C_hat_nop1 = min(F_nop1);
                        k = k_opt_beta;
                        Beta_round_C_hat;
                        
                        Beta_out = [pp(N); pp(N-2); pp(N-3); pp(N-4); Beta_out; pp(N-1)];
                        beta_print = [pp(N); pp(N-2); pp(N-3); pp(N-4); beta; pp(N-1)];
                        beta_1_nop1_print = [pp(N); pp(N-2); pp(N-3); pp(N-4); beta_nop1; pp(N-1)];
                            
                        CTV_k = CTV_p(Beta_out);
                        
                        F_C_hat_1 = F;
                        beta_incumbent = Beta_out;
                                
                        Dist_1 = abs(beta(aa:bb)' - pp_opt(4:(N-2)))./pp_opt(4:(N-2));     
                            
                        k_opt_1 = k;

                        F_C_hat = F_C_hat_1;
                        Dist_C_hat = Dist_1;
                        k_opt_C_hat = k_opt_1;
                        PasteBeta_C_hat = CTV_k;

                        fprintf(fid2, [outputstr ' \r\n'], beta_print);                        
                        fprintf(fid3, [outputstr ' \r\n'], beta_1_nop1_print);
                        %fprintf(fid1,'\r\n');
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Print the CPLEX output
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        fprintf(fid1,'%-8d %-8d %-8d %-8d %-8d %-8.2f %-8d %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f \r\n', Instance, N, type, L, S, coef(jj), ii, pp(N-1), p_bar, CTV_AS, CTV_opt, LB_vs, LB_nn, Non_convex_beta_bound, F_C_hat_nop1, F_C_hat, PasteBeta_C_hat, k_heu, k_opt_C_hat, LL_1, UU_1, LL_Mann, UU_Mann, mean(Dist_C_hat), max(Dist_C_hat), min(Dist_C_hat), CPU_time_QP_C_hat, CPU_time_heur);
                        %fprintf(fid1,'\r\n');
                            
                        Instance = Instance + 1;
                    
                    end %S
                    
                end %case
            end
        end
    end
end

fclose(fid1);  
fclose(fid2);  
fclose(fid3);  
