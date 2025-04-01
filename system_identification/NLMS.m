%NLMS
% clc;
% clear all;
% close all;

for itr = 1:5
    n = 1*10^(4);
    c = rand(1,n)-0.5;
    d = filter(1,[1 -0.8],c);
    noise = awgn(d,30)-d;
    delta=var(d);
    w = [zeros(1,35) -0.005 -0.006 -0.008 -0.05 -0.1 -0.07 0.01  0.15 0.25 -0.16 -0.23 0.15 0.15 -0.05 -0.14 0.1 0.075 -0.04 -0.08 0 0.05 0.04 0.025 0.01 -0.02 -0.01 0.03 0.01 -0.01 -0.025 -0.01 0.025 0.02 -0.02 -0.01 0.01 0.02 0.01 -0.01 -0.01 0.01 0.015 -0.008 -0.0075 0.005 -0.02 0 0.015 0.015 0.0003 -0.0003 0.0003 -0.003 0.003 0.003 0 0 -0.00003 0.00003 0 0.008 0.005 0.0008 -0.01 -0.01];
    L = length(w);

    mu_ABxPNLMS = 0.36;
    mu_IPNLMS = 0.36;
    mu_NLMS = 0.36;
    mu_AByPNLMS = 0.36;

    sys_tap = zeros(1,length(w));

    w_est_NLMS = 0.00001*ones(1,length(w));      %Estimated Weight Vector NLMS
    
    model_tap = zeros(1,length(w));
    D = zeros(1,length(w));

    for i=1:n

        sys_tap = [d(i) sys_tap(1:end-1)];
        D(i) = sys_tap*w' + noise(i);
        model_tap = [d(i) model_tap(1:end-1)];

        model_out_NLMS(i) = model_tap*w_est_NLMS';
        
        e_NLMS(i)=D(i)-model_out_NLMS(i);
        w_est_NLMS = w_est_NLMS + (mu_NLMS/(model_tap*model_tap' + 0.001))*model_tap*e_NLMS(i);
        mdl_NLMS(i,:) = w_est_NLMS;
    end

    for j=1:n
        E_NLMS(itr,j)=mean((w-mdl_NLMS(j,:)).^2);
    end
    error_NLMS(itr,:)=e_NLMS.^2;
    clc
    iteration = itr
end

E1_NLMS=mean(error_NLMS);

figure;
hold on
title('MSE Plot');xlabel('Iterations');ylabel('MSE (dB)');grid on;

plot(10*log10(E1_NLMS),'b');
hold off
legend('NLMS');
figure;
hold on
grid on;title('MSD Plot');xlabel('Iterations');ylabel('MSD (dB)')

plot(10*log10(mean(E_NLMS)),'b');
legend('NLMS');
hold off