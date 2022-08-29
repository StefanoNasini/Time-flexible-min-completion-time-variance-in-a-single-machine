function [Sigma] = H_NC(p, sigma)

N = length(p);

Incumbent_0 = CTV(p, sigma);
GAP = 1;
epsilon = 1.0e-5;

Iteration = 1;

pp = p(sigma);

while GAP > epsilon

    for m = 1:(N-1)
        for s = (m+1):N

                %pp = p(sigma);

                %C_bar = 0;
                %for j = 1:N
                %    C_j = sum(pp(1:j));
                %    C_bar = C_bar + C_j;
                %end
                %C_bar = C_bar/N;
                
                C_bar = [N:-1:1]*pp'/N;

                P_sm = pp(s) - pp(m);
                f_sm = P_sm*(s-m);
                cost_sm = (f_sm^2)/N + P_sm*(P_sm - 2*f_sm/N)*(s-m);

                %SUM = 0;
                %for j = (m+1):s
                %    C_j = sum(pp(1:j));
                %    SUM = SUM + C_j - C_bar;
                %end
                
                SUM = [repmat((s-m),1,m+1) (s-m-1):-1:1]*pp(1:s)' - (s-m)*C_bar;
                
                P = sum(pp((m+1):s));

                Delta_sm = cost_sm + 2*P_sm*SUM - P*(2*P_sm);

                if Delta_sm < 0
                    x = sigma(s);
                    y = pp(s);
                    
                    sigma(s) = sigma(m);
                    sigma(m) = x;
                    
                    pp(s) = pp(m);
                    pp(m) = y;  
                    
                end

        end
    end
    
    Incumbent_1 = CTV(p, sigma);
    
    GAP = (Incumbent_0 - Incumbent_1)/Incumbent_0;
    
    Incumbent_0 = Incumbent_1;
    
    sprintf('Iteration %d, Gap = %-8.10f', Iteration, GAP)
    
    Iteration = Iteration + 1;
        
end

Sigma = sigma;

end