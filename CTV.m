function [CTV_val] = CTV(p, sigma)

    N = length(p);
    C_bar = 0;
    C = zeros(1,N);
    pp = p(sigma);
    
    for j = 1:N
        C(j) = sum(pp(1:j));
        C_bar = C_bar + C(j);
    end
    C_bar = C_bar/N;

    CTV_val = 0;
    for j = 1:N
        CTV_val = CTV_val + (C(j) - C_bar)^2;
    end
    CTV_val = CTV_val/N;
         
end