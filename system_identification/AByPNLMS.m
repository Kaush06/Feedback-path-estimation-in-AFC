%AByPNLMS
% clc;
% clear all;
% close all;

for itr = 1:5
    n = 1*10^(4);
    c = rand(1,n)-0.5;
    x = filter(1,[1 -0.8],c);
    noise = awgn(x,30)-x;
    delta=var(x);
    w = [zeros(1,35) -0.005 -0.006 -0.008 -0.05 -0.1 -0.07 0.01  0.15 0.25 -0.16 -0.23 0.15 0.15 -0.05 -0.14 0.1 0.075 -0.04 -0.08 0 0.05 0.04 0.025 0.01 -0.02 -0.01 0.03 0.01 -0.01 -0.025 -0.01 0.025 0.02 -0.02 -0.01 0.01 0.02 0.01 -0.01 -0.01 0.01 0.015 -0.008 -0.0075 0.005 -0.02 0 0.015 0.015 0.0003 -0.0003 0.0003 -0.003 0.003 0.003 0 0 -0.00003 0.00003 0 0.008 0.005 0.0008 -0.01 -0.01];
    % a_1 = load('C:\Users\skaus\OneDrive\Desktop\BTP\HA_Code\acoustic_paths\g1.mat');
    % w = a_1.g1;
    L = length(w);
    kk = 0.275;%best result


    mu_ABxPNLMS = 0.36;
    mu_IPNLMS = 0.36;
    mu_NLMS = 0.36;
    mu_AByPNLMS = 0.36;

    beta = 0.95;%fix
    gamma_AByPNLMS = 0.225;
    gamma_ABxPNLMS = 0.225;

    sigma2_v = 0.001;%fix
    sigma2_e = 0;

    sys_tap = zeros(1,length(w));
w_est_AByPNLMS = 0.00001*ones(1,length(w));   %Estimated Weight Vector AByPNLMS
        
    model_tap = zeros(1,length(w));
    sys_out = zeros(1,length(w));

        for i=1:n

        sys_tap = [x(i) sys_tap(1:end-1)];
        sys_out(i) = sys_tap*w' + noise(i);
        model_tap = [x(i) model_tap(1:end-1)];

        model_out_AByPNLMS(i) = model_tap*w_est_AByPNLMS';
        e_AByPNLMS(i) = sys_out(i)-model_out_AByPNLMS(i);
        sigma2_e = beta*sigma2_e + (1-beta)*(e_AByPNLMS(i))*(e_AByPNLMS(i));
        lambda = (sigma2_e - sigma2_v)/(sigma2_e + 0.001);
        zeta = kk*lambda*(norm(w_est_AByPNLMS,1))/L;

        a1 = 0;
        N=0; 
        for l = 1:L
             if (abs(w_est_AByPNLMS(l))<=zeta)
                 a1=a1+abs(w_est_AByPNLMS(l));
                 N=N+1;
             end
         end
         g1 = (1-gamma_AByPNLMS)/(2*L) + (1+gamma_AByPNLMS)*a1/((N)*(2*norm(w_est_AByPNLMS,1) + 0.001));
         g21 = (1-g1*N)/(L-N);
         for l = 1:L
             if (N==0)
                 p_AByPNLMS(l) = 1/L;
             elseif (abs(w_est_AByPNLMS(l))<=zeta)
                 p_AByPNLMS(l) = g1;
             else
                 p_AByPNLMS(l) = g21;
             end
         end
         P_AByPNLMS = diag(p_AByPNLMS);

         %%Weight Update step
         w_est_AByPNLMS = w_est_AByPNLMS + mu_AByPNLMS*((model_tap)*(P_AByPNLMS)/((model_tap)*P_AByPNLMS*(model_tap') + 0.001)).*e_AByPNLMS(i);
          % [w_sorted,index_sorted] = sort(w_est_AByPNLMS,"descend","ComparisonMethod","abs");
         mdl_AByPNLMS(i,:) = w_est_AByPNLMS;

        end

        for j=1:n
        E_AByPNLMS(itr,j)=mean((w-mdl_AByPNLMS(j,:)).^2);
        end
    error_AByPNLMS(itr,:)=e_AByPNLMS.^2;
    clc
    iteration = itr
end

E1_AByPNLMS=mean(error_AByPNLMS);

figure;
hold on
title('MSE Plot');xlabel('Iterations');ylabel('MSE (dB)');grid on;
plot(10*log10(E1_AByPNLMS),'y');
hold off
legend('AByPNLMS');
figure;
hold on
grid on;title('MSD Plot');xlabel('Iterations');ylabel('MSD (dB)')

plot(10*log10(mean(E_AByPNLMS)),'y');
legend('AByPNLMS');
hold off