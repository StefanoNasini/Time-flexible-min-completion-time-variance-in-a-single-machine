function [CTV_LB_basic] = LB_basic(p, a)

    N = length(p);

    m = (N-6)/2;
    MS = sum(p);
    
    P_bar = p(N-2) + (4*p(N-3) + sum(p(1:(N-4))))/2;
    
    if abs(floor(N/2) - N/2) > 0
            LB0 = 0;
            for i = 1:m
                LB0 = LB0 + (sum(p(1:2*i)))^2;
            end
            LB0 = LB0/2;
    else
        LB0 = 0;
        for i = 1:m
            LB0 = LB0 + (sum(p(1:2*i-1)))^2;
        end
        LB0 = LB0/2;
    end
        
    if a == 0 
        
        LL = 6*p(N)^2 + 6*(p(N)+p(N-2))^2 + 6*(p(N)+p(N-2) + p(N-3))^2 + 6*MS^2 + 6*(MS - p(N-1))^2;

        if p(N-1) <= P_bar
            L0 = LL + 6*(MS - p(N-1) - p(N-4))^2 - (3*MS + 3*p(N) - 2*p(N-1) + 2*p(N-2) + p(N-3) - p(N-4))^2;
        else
            L0 = LL + 6*(p(N)+p(N-2) + p(N-3) + p(N-4))^2 - (2*MS + 4*p(N) - p(N-1) + 3*p(N-2) + 2*p(N-3) + p(N-4))^2;
        end
        
    else
        
        LL = 6*p(N)^2 + 6*(p(N)+p(N-2))^2 + 6*(p(N)+p(N-2) + p(N-4))^2 + 6*MS^2 + 6*(MS - p(N-1))^2;

        if p(N-1) <= P_bar
            L0 = LL + 6*(MS - p(N-1) - p(N-3))^2 - (3*MS + 3*p(N) - 2*p(N-1) + 2*p(N-2) + p(N-4) - p(N-3))^2;
        else
            L0 = LL + 6*(p(N)+p(N-2) + p(N-3) + p(N-4))^2 - (2*MS + 4*p(N) - p(N-1) + 3*p(N-2) + 2*p(N-3) + p(N-4))^2;
        end
        
    end
    
    CTV_LB_basic = (L0/6 + LB0)/N;
    
end
