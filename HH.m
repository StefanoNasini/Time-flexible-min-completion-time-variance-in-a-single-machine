function [h Yk ggg] = HH(k,gg)

global N;
global Beta;
global Gamma;
global gamma_star;

if k == N
    h = 0;
    Yk = 0;
    ggg = 0;
    return
end

[HHH_former y_former gg_former] = HH(k + 1, gg);
[HHH_latter y_latter gg_latter] = HH(k + 1, gg + Gamma(k));

Former = HHH_former + Beta(k)*gg;
Latter = HHH_latter + Beta(k)*(gamma_star(k) - gg);

if Former >= Latter
    Yk = [y_former 0];
    ggg = [gg_former gg];
else
    Yk = [y_former 1];
    ggg = [gg_former (gg + Gamma(k))];
end

h = max( Former, Latter); 

end