
if p_bar_case == 0
    
    fid = fopen('.\cases.dat', 'w');          
    fprintf(fid,'\nparam Case:= %g;\nparam LB0:= %-8.10f; \nparam LLB_1:= %-8.10f; \nparam LLB_2:= %-8.10f; \nparam UUB_1:= %-8.10f; \nparam UUB_2:= %-8.10f; \nparam Ttheta_1:= %-8.10f; \nparam Ttheta_2:= %-8.10f; \nparam N:= %g; \nparam L:= %d;\nparam S:= %d;\nparam LL_1:= %d;\nparam UU_1:= %d;\nparam LL_2:= %d;\nparam UU_2:= %d;\n', 0, LB0, LB_1, LB_2, UB_1, UB_2, theta_1, theta_2, N, L, S, LL_1, UU_1, LL_2, UU_2);

    fprintf(fid,'\n\nparam q:=');

    for i=1:N
        fprintf(fid,'\n%d %g', i, pp(i));
    end
    fprintf(fid,';\n');

    fprintf(fid,'\n\nparam G:=');
    for s=1:S
        fprintf(fid,'\n%d %g', s, G_s(s));
    end
    fprintf(fid,';\n');

    fclose(fid);   

    %----------------------------------------------

    if C_bound == 0
        
        tic
        [status, result] = system('ampl Main_beta_hat.run');             
        CPU_time_QP = toc;
        
        result
        
        AMPL_F_LB_1 = fopen('.\output\AMPL_F_LB_1.txt','r');
        AMPL_F_UB_1 = fopen('.\output\AMPL_F_UB_1.txt','r');
        AMPL_F_LB_2 = fopen('.\output\AMPL_F_LB_2.txt','r');
        AMPL_F_UB_2 = fopen('.\output\AMPL_F_UB_2.txt','r'); 

        AMPL_F_LB_1_nop1 = fopen('.\output\AMPL_F_LB_1_nop1.txt','r');
        AMPL_F_UB_1_nop1 = fopen('.\output\AMPL_F_UB_1_nop1.txt','r');
        AMPL_F_LB_2_nop1 = fopen('.\output\AMPL_F_LB_2_nop1.txt','r');
        AMPL_F_UB_2_nop1 = fopen('.\output\AMPL_F_UB_2_nop1.txt','r'); 
        
        AMPL_beta_LB_1_nop1 = fopen('.\output\AMPL_F_LB_1_beta_nop1.txt','r');
        AMPL_beta_UB_1_nop1 = fopen('.\output\AMPL_F_UB_1_beta_nop1.txt','r');
        AMPL_beta_LB_2_nop1 = fopen('.\output\AMPL_F_LB_2_beta_nop1.txt','r');
        AMPL_beta_UB_2_nop1 = fopen('.\output\AMPL_F_UB_2_beta_nop1.txt','r');                         
        
        AMPL_beta_LB_1 = fopen('.\output\AMPL_F_LB_1_beta.txt','r');
        AMPL_beta_UB_1 = fopen('.\output\AMPL_F_UB_1_beta.txt','r');
        AMPL_beta_LB_2 = fopen('.\output\AMPL_F_LB_2_beta.txt','r');
        AMPL_beta_UB_2 = fopen('.\output\AMPL_F_UB_2_beta.txt','r'); 
        
        F_LB_1 = fscanf(AMPL_F_LB_1,'%f');
        F_UB_1 = fscanf(AMPL_F_UB_1,'%f');
        F_LB_2 = fscanf(AMPL_F_LB_2,'%f');
        F_UB_2 = fscanf(AMPL_F_UB_2,'%f');  
        
        F_LB_1_nop1 = fscanf(AMPL_F_LB_1_nop1,'%f');
        F_UB_1_nop1 = fscanf(AMPL_F_UB_1_nop1,'%f');
        F_LB_2_nop1 = fscanf(AMPL_F_LB_2_nop1,'%f');
        F_UB_2_nop1 = fscanf(AMPL_F_UB_2_nop1,'%f');
        
        beta_LB_1_nop1 = fscanf(AMPL_beta_LB_1_nop1,'%f');
        beta_UB_1_nop1 = fscanf(AMPL_beta_UB_1_nop1,'%f');
        beta_LB_2_nop1 = fscanf(AMPL_beta_LB_2_nop1,'%f');
        beta_UB_2_nop1 = fscanf(AMPL_beta_UB_2_nop1,'%f');
        
        beta_LB_1 = fscanf(AMPL_beta_LB_1,'%f');
        beta_UB_1 = fscanf(AMPL_beta_UB_1,'%f');
        beta_LB_2 = fscanf(AMPL_beta_LB_2,'%f');
        beta_UB_2 = fscanf(AMPL_beta_UB_2,'%f');
        
        fclose(AMPL_F_LB_1);
        fclose(AMPL_F_UB_1);
        fclose(AMPL_F_LB_2);
        fclose(AMPL_F_UB_2);
        
        fclose(AMPL_F_LB_1_nop1);
        fclose(AMPL_F_UB_1_nop1);
        fclose(AMPL_F_LB_2_nop1);
        fclose(AMPL_F_UB_2_nop1);
        
        fclose(AMPL_beta_LB_1_nop1);
        fclose(AMPL_beta_UB_1_nop1);
        fclose(AMPL_beta_LB_2_nop1);
        fclose(AMPL_beta_UB_2_nop1); 
        
        fclose(AMPL_beta_LB_1);
        fclose(AMPL_beta_UB_1);
        fclose(AMPL_beta_LB_2);
        fclose(AMPL_beta_UB_2);  
                        
        delete('.\output\AMPL_F_LB_1.txt');
        delete('.\output\AMPL_F_UB_1.txt');
        delete('.\output\AMPL_F_LB_2.txt');
        delete('.\output\AMPL_F_UB_2.txt');

        delete('.\output\AMPL_F_LB_1_nop1.txt');
        delete('.\output\AMPL_F_UB_1_nop1.txt');
        delete('.\output\AMPL_F_LB_2_nop1.txt');
        delete('.\output\AMPL_F_UB_2_nop1.txt');
        
        delete('.\output\AMPL_F_LB_1_beta.txt');
        delete('.\output\AMPL_F_UB_1_beta.txt');
        delete('.\output\AMPL_F_LB_2_beta.txt');
        delete('.\output\AMPL_F_UB_2_beta.txt');

        delete('.\output\AMPL_F_LB_1_beta_nop1.txt');
        delete('.\output\AMPL_F_UB_1_beta_nop1.txt');
        delete('.\output\AMPL_F_LB_2_beta_nop1.txt');
        delete('.\output\AMPL_F_UB_2_beta_nop1.txt');
                      
    else
        
        tic
        [status, result] = system('ampl Main_beta_hat_C_hat.run');             
        CPU_time_QP_C_hat = toc;

        result
        
        AMPL_F_1 = fopen('.\output\AMPL_F_1.txt','r');
        AMPL_F_2 = fopen('.\output\AMPL_F_2.txt','r');
                        
        AMPL_F_1_nop1 = fopen('.\output\AMPL_F_1_nop1.txt','r');
        AMPL_F_2_nop1 = fopen('.\output\AMPL_F_2_nop1.txt','r');
        
        AMPL_beta_1_nop1 = fopen('.\output\AMPL_F_1_beta_nop1.txt','r');
        AMPL_beta_2_nop1 = fopen('.\output\AMPL_F_2_beta_nop1.txt','r');                     
        
        AMPL_beta_1 = fopen('.\output\AMPL_F_1_beta.txt','r');
        AMPL_beta_2 = fopen('.\output\AMPL_F_2_beta.txt','r');
        
        F_1 = fscanf(AMPL_F_1,'%f');
        F_2 = fscanf(AMPL_F_2,'%f');  
        
        F_1_nop1 = fscanf(AMPL_F_1_nop1,'%f');
        F_2_nop1 = fscanf(AMPL_F_2_nop1,'%f'); 
        
        beta_1_nop1 = fscanf(AMPL_beta_1_nop1,'%f');
        beta_2_nop1 = fscanf(AMPL_beta_2_nop1,'%f');
        
        beta_1 = fscanf(AMPL_beta_1,'%f');
        beta_2 = fscanf(AMPL_beta_2,'%f');
        
        fclose(AMPL_F_1);
        fclose(AMPL_F_2); 
        fclose(AMPL_F_1_nop1);
        fclose(AMPL_F_2_nop1);         
        fclose(AMPL_beta_1_nop1);
        fclose(AMPL_beta_2_nop1);          
        fclose(AMPL_beta_1);
        fclose(AMPL_beta_2);  
        
        delete('.\output\AMPL_F_1.txt');
        delete('.\output\AMPL_F_2.txt');
        delete('.\output\AMPL_F_1_nop1.txt');
        delete('.\output\AMPL_F_2_nop1.txt');
        delete('.\output\AMPL_F_1_beta.txt');
        delete('.\output\AMPL_F_2_beta.txt');
        delete('.\output\AMPL_F_1_beta_nop1.txt');
        delete('.\output\AMPL_F_2_beta_nop1.txt');
    
    end
end


if p_bar_case == 2
    
    G_s = zeros(S,1);
    for s = 1:S
        G_s(s) = min(pp(2:N) - pp(1:(N-1)));
    end
    
    L = 2*(N-6);
    
    fid = fopen('.\cases.dat', 'w');          
    fprintf(fid,'\nparam Case:= %g; \nparam LB0:= %-8.10f; \nparam LLB_1:= %-8.10f; \nparam UUB_1:= %-8.10f; \nparam Ttheta_1:= %-8.10f; \nparam N:= %g; \nparam L:= %d;\nparam S:= %d;\nparam LL_1:= %d;\nparam UU_1:= %d;\n', 2, LB0, LB_1, UB_1, theta_1, N, L, S, LL_1, UU_1);
    
    fprintf(fid,'\n\nparam q:=');
    for i=1:N
        fprintf(fid,'\n%d %g', i, pp(i));
    end
    
    fprintf(fid,';\n');
    fprintf(fid,'\n\nparam G:=');
    for s=1:S
        fprintf(fid,'\n%d %g', s, G_s(s));
    end
    fprintf(fid,';\n');
    fclose(fid); 
    
    %----------------------------------------------

    if C_bound == 0                                                     
        
        tic
        [status, result] = system('ampl Main_beta_hat.run');             
        CPU_time_QP = toc;
        
        result
        
        AMPL_beta_LB = fopen('.\output\AMPL_F_LB_beta.txt','r');
        AMPL_beta_UB = fopen('.\output\AMPL_F_UB_beta.txt','r');
        AMPL_beta_LB_nop1 = fopen('.\output\AMPL_F_LB_beta_nop1.txt','r');
        AMPL_beta_UB_nop1 = fopen('.\output\AMPL_F_UB_beta_nop1.txt','r');              
        AMPL_F_LB = fopen('.\output\AMPL_F_LB.txt','r');
        AMPL_F_UB = fopen('.\output\AMPL_F_UB.txt','r');
        AMPL_F_LB_nop1 = fopen('.\output\AMPL_F_LB_nop1.txt','r');
        AMPL_F_UB_nop1 = fopen('.\output\AMPL_F_UB_nop1.txt','r');
        
        beta_LB = fscanf(AMPL_beta_LB,'%f');
        beta_UB = fscanf(AMPL_beta_UB,'%f');
        beta_LB_nop1 = fscanf(AMPL_beta_LB_nop1,'%f');
        beta_UB_nop1 = fscanf(AMPL_beta_UB_nop1,'%f');
        F_LB = fscanf(AMPL_F_LB,'%f');
        F_UB = fscanf(AMPL_F_UB,'%f');
        F_LB_nop1 = fscanf(AMPL_F_LB_nop1,'%f');
        F_UB_nop1 = fscanf(AMPL_F_UB_nop1,'%f');
        
        fclose(AMPL_beta_LB);
        fclose(AMPL_beta_UB);                        
        fclose(AMPL_beta_LB_nop1);
        fclose(AMPL_beta_UB_nop1);
        fclose(AMPL_F_LB);
        fclose(AMPL_F_UB);      
        fclose(AMPL_F_LB_nop1);
        fclose(AMPL_F_UB_nop1);  
        
        delete('.\output\AMPL_F_LB.txt');
        delete('.\output\AMPL_F_UB.txt');
        delete('.\output\AMPL_F_LB_nop1.txt');
        delete('.\output\AMPL_F_UB_nop1.txt');
        delete('.\output\AMPL_F_LB_beta.txt');
        delete('.\output\AMPL_F_UB_beta.txt');
        delete('.\output\AMPL_F_LB_beta_nop1.txt');
        delete('.\output\AMPL_F_UB_beta_nop1.txt');
                        
    else

        tic
        [status, result] = system('ampl Main_beta_hat_C_hat.run');             
        CPU_time_QP_C_hat = toc;

        result
        
        AMPL_F = fopen('.\output\AMPL_F.txt','r');  
        AMPL_F_nop1 = fopen('.\output\AMPL_F_nop1.txt','r');  
        AMPL_beta_nop1 = fopen('.\output\AMPL_F_beta_nop1.txt','r');                 
        AMPL_beta = fopen('.\output\AMPL_F_beta.txt','r');
        
        F = fscanf(AMPL_F,'%f');   
        F_nop1 = fscanf(AMPL_F_nop1,'%f');   
        beta_nop1 = fscanf(AMPL_beta_nop1,'%f');       
        beta = fscanf(AMPL_beta,'%f');    
        
        fclose(AMPL_F);  
        fclose(AMPL_F_nop1);  
        fclose(AMPL_beta_nop1);        
        fclose(AMPL_beta);
        
        delete('.\output\AMPL_F.txt');
        delete('.\output\AMPL_F_nop1.txt');
        delete('.\output\AMPL_F_beta.txt');
        delete('.\output\AMPL_F_beta_nop1.txt');
            
    end
end


if p_bar_case == 1
    
    G_s = zeros(S,1);
    for s = 1:S
        G_s(s) = min(pp(2:N) - pp(1:(N-1)));
    end
    
    L = 2*(N-6);
    fid = fopen('.\cases.dat', 'w');          
    fprintf(fid,'\nparam Case:= %g; \nparam LB0:= %-8.10f; \nparam LLB:= %-8.10f; \nparam UUB:= %-8.10f; \nparam Ttheta:= %-8.10f; \nparam N:= %g; \nparam L:= %d;\nparam S:= %d;\nparam LL_1:= %d;\nparam UU_1:= %d;\n', 1, LB0, LB, UB, theta_6, N, L, S, LL_1, UU_1);

    fprintf(fid,'\n\nparam q:=');
    for i=1:N
        fprintf(fid,'\n%d %g', i, pp(i));
    end
    fprintf(fid,';\n');
                        
    fprintf(fid,'\n\nparam G:=');
    for s=1:S
        fprintf(fid,'\n%d %g', s, G_s(s));
    end
    fprintf(fid,';\n');
    fclose(fid); 
    
    
    if C_bound == 0     
        
        %----------------------------------------------
        % Quadratic formulation
        %----------------------------------------------
        
        tic
        [status, result] = system('ampl Main_beta_hat.run');             
        CPU_time_QP = toc;
        
        result
        
        AMPL_beta_LB = fopen('.\output\AMPL_F_LB_beta.txt','r');
        AMPL_beta_UB = fopen('.\output\AMPL_F_UB_beta.txt','r');
        AMPL_beta_LB_nop1 = fopen('.\output\AMPL_F_LB_beta_nop1.txt','r');
        AMPL_beta_UB_nop1 = fopen('.\output\AMPL_F_UB_beta_nop1.txt','r');              
        AMPL_F_LB = fopen('.\output\AMPL_F_LB.txt','r');
        AMPL_F_UB = fopen('.\output\AMPL_F_UB.txt','r');
        AMPL_F_LB_nop1 = fopen('.\output\AMPL_F_LB_nop1.txt','r');
        AMPL_F_UB_nop1 = fopen('.\output\AMPL_F_UB_nop1.txt','r');
        
        beta_LB = fscanf(AMPL_beta_LB,'%f');
        beta_UB = fscanf(AMPL_beta_UB,'%f');
        beta_LB_nop1 = fscanf(AMPL_beta_LB_nop1,'%f');
        beta_UB_nop1 = fscanf(AMPL_beta_UB_nop1,'%f');
        F_LB = fscanf(AMPL_F_LB,'%f');
        F_UB = fscanf(AMPL_F_UB,'%f');
        F_LB_nop1 = fscanf(AMPL_F_LB_nop1,'%f');
        F_UB_nop1 = fscanf(AMPL_F_UB_nop1,'%f');
        
        fclose(AMPL_beta_LB);
        fclose(AMPL_beta_UB);                        
        fclose(AMPL_beta_LB_nop1);
        fclose(AMPL_beta_UB_nop1);
        fclose(AMPL_F_LB);
        fclose(AMPL_F_UB);      
        fclose(AMPL_F_LB_nop1);
        fclose(AMPL_F_UB_nop1);  
        
        delete('.\output\AMPL_F_LB.txt');
        delete('.\output\AMPL_F_UB.txt');
        delete('.\output\AMPL_F_LB_nop1.txt');
        delete('.\output\AMPL_F_UB_nop1.txt');
        delete('.\output\AMPL_F_LB_beta.txt');
        delete('.\output\AMPL_F_UB_beta.txt');
        delete('.\output\AMPL_F_LB_beta_nop1.txt');
        delete('.\output\AMPL_F_UB_beta_nop1.txt');
                        
    else
        
        tic
        [status, result] = system('ampl Main_beta_hat_C_hat.run');             
        CPU_time_QP_C_hat = toc;

        result
        
        AMPL_F = fopen('.\output\AMPL_F.txt','r');  
        AMPL_F_nop1 = fopen('.\output\AMPL_F_nop1.txt','r');  
        AMPL_beta_nop1 = fopen('.\output\AMPL_F_beta_nop1.txt','r');                 
        AMPL_beta = fopen('.\output\AMPL_F_beta.txt','r');
        
        F = fscanf(AMPL_F,'%f');   
        F_nop1 = fscanf(AMPL_F_nop1,'%f');   
        beta_nop1 = fscanf(AMPL_beta_nop1,'%f');       
        beta = fscanf(AMPL_beta,'%f');    
        
        fclose(AMPL_F);  
        fclose(AMPL_F_nop1);  
        fclose(AMPL_beta_nop1);        
        fclose(AMPL_beta);
        
        delete('.\output\AMPL_F.txt');
        delete('.\output\AMPL_F_nop1.txt');
        delete('.\output\AMPL_F_beta.txt');
        delete('.\output\AMPL_F_beta_nop1.txt');
        
    end
    
end

    
