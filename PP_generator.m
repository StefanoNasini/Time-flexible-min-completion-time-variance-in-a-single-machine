function [PP] = PP_generator(N, coef, seed, type)

    rng(seed);
    
    NBS = length(coef);
    PP = cell(1,NBS);  
    
    if type == 0

        for i = 1:NBS

            r = floor(N*coef(i));
            pp = zeros(1,N-3);

            for k = 1:N
                pp(k) = round(1 + (r-1)*rand());
            end    

            pp = sort(pp);   

            %pp(N-2) = round(sum(pp(1:(N-3)))/3);
            %pp(N-1) = round(sum(pp(1:(N-3)))/2);
            %pp(N) = pp(N-1) + 1;
            
            PP{i} = pp;

        end

    else
        
        for i = 1:NBS
            
            pp = zeros(1,N);
            pp(1) = round(coef(i)*rand(1,1) + 1);
            
            for ii = 2:N
                pp(ii) = pp(ii-1) + round(coef(i)*rand(1,1) + 1);
            end
            
            pp(N-2) = round(sum(pp(1:(N-3)))/3);
            pp(N-1) = round(sum(pp(1:(N-3)))/2);
            pp(N) = pp(N-1) + 1;
                    
            PP{i} = pp;
            
        end
        
    end
end
