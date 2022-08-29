function [CTV_LB_VS] = LB_VS(p)

    N = length(p);

    m = (N-4)/2;
    MS = sum(p);

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
    
    T = (2*p(N) + p(N-2) + 2*MS - p(N-1))/4;
    L0 = ((p(N) - T)^2 + (p(N) + p(N-2)- T)^2 + (MS - T)^2 + (MS - p(N-1) - T)^2); 
      
    CTV_LB_VS = (L0 + LB0)/N;
    
end