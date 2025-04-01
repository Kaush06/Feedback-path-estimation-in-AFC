 %IPNLMS minus sparsity aware
% clc;
% clear all;
% close all;

for itr = 1:5
    n = 2*10^(4);
    c = rand(1,n)-0.5;
    x = filter(1,[1 -0.8],c);
    noise = awgn(x,30)-x;
    delta=var(x);
    w = [zeros(1,35) -0.005 -0.006 -0.008 -0.05 -0.1 -0.07 0.01  0.15 0.25 -0.16 -0.23 0.15 0.15 -0.05 -0.14 0.1 0.075 -0.04 -0.08 0 0.05 0.04 0.025 0.01 -0.02 -0.01 0.03 0.01 -0.01 -0.025 -0.01 0.025 0.02 -0.02 -0.01 0.01 0.02 0.01 -0.01 -0.01 0.01 0.015 -0.008 -0.0075 0.005 -0.02 0 0.015 0.015 0.0003 -0.0003 0.0003 -0.003 0.003 0.003 0 0 -0.00003 0.00003 0 0.008 0.005 0.0008 -0.01 -0.01];
    % a_1 = load('C:\Users\skaus\OneDrive\Desktop\BTP\HA_Code\acoustic_paths\g1.mat');
    % w = a_1.g1;
    L = length(w);
    
    mu_ABxPNLMS = 0.36;
    mu_IPNLMS = 0.36;
    mu_NLMS = 0.36;
    mu_AByPNLMS = 0.36;

    sys_tap = zeros(1,length(w));

w_est_IPNLMS = 0.00001*ones(1,length(w));      %Estimated Weight Vector IPNLMS
        
    model_tap = zeros(1,length(w));
    sys_out = zeros(1,length(w));

    alpha=0.5;

    for i=1:n

        sys_tap = [x(i) sys_tap(1:end-1)];
        sys_out(i) = sys_tap*w' + noise(i);
        model_tap = [x(i) model_tap(1:end-1)];

        model_out_IPNLMS(i) = model_tap*w_est_IPNLMS';

        e_IPNLMS(i) = sys_out(i) - model_out_IPNLMS(i);
    
      W1=0;
      for k= 1:L
           W1=W1+abs(w_est_IPNLMS(k));
      end
       temp = W1;
       P_temp = zeros(1,L);
      for j=1:L
         P_temp(j)=((1-alpha)/(2*L))+((1+alpha)*abs(w_est_IPNLMS(j))/(2*W1+0.001));
      end
      phi_IPNLMS=diag(P_temp);

      w_est_IPNLMS  =  w_est_IPNLMS  +  (mu_IPNLMS/(model_tap*phi_IPNLMS*model_tap'+delta)) *(phi_IPNLMS*model_tap')'.* e_IPNLMS(i);  
    mdl_IPNLMS(i,:) = w_est_IPNLMS;
    
    end
        for j=1:n
        E_IPNLMS(itr,j)=mean((w-mdl_IPNLMS(j,:)).^2);
    end

    error_IPNLMS(itr,:)=e_IPNLMS.^2;
    clc
    iteration = itr
end

E1_IPNLMS=mean(error_IPNLMS);

figure;
hold on
title('MSE Plot');xlabel('Iterations');ylabel('MSE (dB)');grid on;

plot(10*log10(E1_IPNLMS),'g');
hold off
legend('IPNLMS');
figure;
hold on
grid on;title('MSD Plot');xlabel('Iterations');ylabel('MSD (dB)')

plot(10*log10(mean(E_IPNLMS)),'g');
legend('IPNLMS');
hold off