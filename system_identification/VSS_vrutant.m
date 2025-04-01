clc;
clear all;
close all;

sys_desired = [86 -294 -287 -262 -120 140 438 641 613 276 -325 -1009 -1487 -1451 -680 856 2954 5206 7106 8192 8192 7106 5206 2954 856 -680 -1451 -1487 -1009 -325 276 613 641 438 140 -120 -262 -287 -294 86] * 2^(-15);

for itr=1:100
    
    x=randn(1,60000);
    % Model for the LMS Algorithm
    model_coeff = zeros(1,length(sys_desired));
    % Model for the VSS-LMS Algorithm
    model_coeff_vss = zeros(1,length(sys_desired));
    model_tap = zeros(1,length(sys_desired));

    noise_floor = 40;
    sys_opt = filter(sys_desired,1,x)+awgn(x,noise_floor)-x;
    
    
    input_var = var(x);
    
    
    mu_max = 1/(input_var*length(sys_desired));
    
    
    mu_LMS = 0.0004;
    
    
    mu_min = mu_LMS;

   
    mu_VSS(1)=1; %initial value of mu for VSS
    alpha  = 0.97;
    gamma = 4.8e-4;
    
    for i=1:length(x)        
        %% LMS Algorithm
        model_tap=[x(i) model_tap(1:end-1)];
        model_out(i) = model_tap * model_coeff';
        e_LMS(i)=sys_opt(i)-model_out(i);
        model_coeff = model_coeff + mu_LMS * e_LMS(i) * model_tap;

        %% Variable stepsize
        model_out_vss(i) = model_tap * model_coeff_vss';
        e_vss(i) = sys_opt(i) - model_out_vss(i);
        model_coeff_vss = model_coeff_vss + mu_VSS(i) * e_vss(i) * model_tap;
        mu_VSS(i+1) = alpha * mu_VSS(i) + gamma * e_vss(i) * e_vss(i) ;

        if (mu_VSS(i+1)>mu_max)
            mu_VSS(i+1)=mu_max;%max
        elseif(mu_VSS(i+1)<mu_min)
            mu_VSS(i+1)= mu_min;
        else
            mu_VSS(i+1) = mu_VSS(i+1) ;
        end
        
    end
    
    
    err_LMS(itr,:) = e_LMS.^2;
    err_VSS(itr,:) = e_vss.^2;
    
    
    clc
    disp(char(strcat('iteration no : ',{' '}, num2str(itr) )))
end


figure;
plot(10*log10(mean(err_LMS)),'-b');
hold on;
plot(10*log10(mean(err_VSS)), '-r');
title('Comparison of LMS and VSS - Q - Newton Algorithms'); xlabel('iterations');ylabel('MSE(dB)');legend('LMS Algorithm','VSS - Q - Newton Algorithms')
grid on;