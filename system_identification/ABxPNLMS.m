%ABxPNLMS

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
    P = 1;    %Number of samples per partitions
    M = L/P;   %Number of Partitions
    kk = 0.275;%best result
    LL = 150;
    k_flag = 0;

    mu_ABxPNLMS = 0.36;
    mu_IPNLMS = 0.36;
    mu_NLMS = 0.36;
    mu_AByPNLMS = 0.36;

    gamma_ABxPNLMS = 0.5;

    sigma2_v = 0.001;%fix
    sigma2_e = 0;

    sys_tap = zeros(1,length(w));
    w_est_ABxPNLMS = 0.00001*ones(1,length(w));      %Estimated Weight Vector ABxPNLMS
        
    model_tap = zeros(1,length(w));
    sys_out = zeros(1,length(w));

        for i=1:n

        sys_tap = [x(i) sys_tap(1:end-1)];
        sys_out(i) = sys_tap*w' + noise(i);
        model_tap = [x(i) model_tap(1:end-1)];
            
        model_out_ABxPNLMS(i) = model_tap*w_est_ABxPNLMS';

        [w_sorted,index_sorted] = sort(w_est_ABxPNLMS,"descend","ComparisonMethod","abs");
         g = zeros(1,L);
         s = 0;
         for m = 1:M
             for j = (m-1)*P+1:m*P
                temp(j - (m-1)*P) = w_sorted(j);
             end
            del = norm(temp,2);
            for xx = (m-1)*P+1:m*P
                g(index_sorted(xx))= del;
            end
            s =s+del;
         end
         
         pp1 = (1-gamma_ABxPNLMS)/(2*L);
         pp2 = ((1+gamma_ABxPNLMS)*g)/(2*s*P + 0.001);
         P_final = pp1+pp2;
         phi_ABxPNLMS=diag(P_final);

         %%Weight Update step

         e_ABxPNLMS(i) = sys_out(i) - model_out_ABxPNLMS(i);
         w_est_ABxPNLMS = w_est_ABxPNLMS + mu_ABxPNLMS*((model_tap)*(phi_ABxPNLMS)/((model_tap)*phi_ABxPNLMS*(model_tap') + 0.001)).*e_ABxPNLMS(i);


         mdl_ABxPNLMS(i,:) = w_est_ABxPNLMS;
        end

    for j=1:n
    E_ABxPNLMS(itr,j)=mean((w-mdl_ABxPNLMS(j,:)).^2);
    end
    error_ABxPNLMS(itr,:)=e_ABxPNLMS.^2;
    clc
    iteration = itr
end

E1_ABxPNLMS=mean(error_ABxPNLMS);

% figure;
% hold on
% title('MSE Plot');xlabel('Iterations');ylabel('MSE (dB)');grid on;
% plot(10*log10(E1_ABxPNLMS),'b');
% hold off
% legend('ABxPNLMS');
figure;
hold on
grid on;title('MSD Plot');xlabel('Iterations');ylabel('MSD (dB)')
plot(10*log10(mean(E_ABxPNLMS)),'b');
legendStr = sprintf('ABxPNLMS_P=%.1f', P);
legend(legendStr);