% Feedback cancellation in  Digital Hearing Aid
% September 2017
% Somanath Pradhan
% Use of this code with out prior permission is a crime
%%  *************************************************************************
% clc;
% close all;
% clear all;
disp('Simulation Started');

a_1=load('Normal_Impulse.mat');% Impulse response
aa_1=(a_1.Normal_Impulse)';

a_2=load('Varying_Impulse.mat');% Impulse response
aa_2=(a_2.Varying_Impulse)';

s1=audioread('D:\BTP\code\HA_Code\Wideband_Speech\sp01.wav');
s2=audioread('D:\BTP\code\HA_Code\Wideband_Speech\sp02.wav');
s3=audioread('D:\BTP\code\HA_Code\Wideband_Speech\sp03.wav');
s4=audioread('D:\BTP\code\HA_Code\Wideband_Speech\sp04.wav');
s5=audioread('D:\BTP\code\HA_Code\Wideband_Speech\sp05.wav');
s6=audioread('D:\BTP\code\HA_Code\Wideband_Speech\sp06.wav');
s7=audioread('D:\BTP\code\HA_Code\Wideband_Speech\sp07.wav');
s8=audioread('D:\BTP\code\HA_Code\Wideband_Speech\sp08.wav');
s9=audioread('D:\BTP\code\HA_Code\Wideband_Speech\sp09.wav');
s10=audioread('D:\BTP\code\HA_Code\Wideband_Speech\sp10.wav');
s11=audioread('D:\BTP\code\HA_Code\Wideband_Speech\sp11.wav');
s12=audioread('D:\BTP\code\HA_Code\Wideband_Speech\sp12.wav');
s13=audioread('D:\BTP\code\HA_Code\Wideband_Speech\sp13.wav');
s14=audioread('D:\BTP\code\HA_Code\Wideband_Speech\sp14.wav');
s15=audioread('D:\BTP\code\HA_Code\Wideband_Speech\sp15.wav');
s16=audioread('D:\BTP\code\HA_Code\Wideband_Speech\sp16.wav');
s17=audioread('D:\BTP\code\HA_Code\Wideband_Speech\sp17.wav');
s18=audioread('D:\BTP\code\HA_Code\Wideband_Speech\sp18.wav');
s19=audioread('D:\BTP\code\HA_Code\Wideband_Speech\sp19.wav');
s20=audioread('D:\BTP\code\HA_Code\Wideband_Speech\sp20.wav');
s21=audioread('D:\BTP\code\HA_Code\Wideband_Speech\sp21.wav');
s22=audioread('D:\BTP\code\HA_Code\Wideband_Speech\sp22.wav');
s23=audioread('D:\BTP\code\HA_Code\Wideband_Speech\sp23.wav');
s24=audioread('D:\BTP\code\HA_Code\Wideband_Speech\sp24.wav');
s25=audioread('D:\BTP\code\HA_Code\Wideband_Speech\sp25.wav');
s26=audioread('D:\BTP\code\HA_Code\Wideband_Speech\sp26.wav');
s27=audioread('D:\BTP\code\HA_Code\Wideband_Speech\sp27.wav');
s28=audioread('D:\BTP\code\HA_Code\Wideband_Speech\sp28.wav');
s29=audioread('D:\BTP\code\HA_Code\Wideband_Speech\sp29.wav');
s30=audioread('D:\BTP\code\HA_Code\Wideband_Speech\sp30.wav');

% Input Signal
Input=0.05*[s1;s2;s3;s4;s5;s6;s7;s8;s9;s10];%s11;s12;s13;s14;s15;s16;s17;s18;s19;s20;s21;s22;s23;s24;s25;s26;s27;s28;s29;s30];
x=Input;

N=length(Input(1:200000));%Input sample length

P = 5;
L=length(aa_1);%Length of Impulse response of acoustic paths
M = L/P;

F_dB=40;%forward path gain in dB scale
F_amp = 8;%10.^(F_dB./20);%forward path gain in normal scale
delay_forw=48;
delta = 1e-3;
mu_h=0.36;

Nfreq=L;
H_1=fft(aa_1,Nfreq);
H_2=fft(aa_2,Nfreq);

for itr=1:1
e_delay = zeros(N+delay_forw,1); 
feedback_tap= zeros(L,1);
feedback_cancel_tap= zeros(L,1);
h_cap= zeros(L,1);

Ma=21;% Order of AR model
Frame_size = 160;
tap=zeros(Frame_size,1); % frame size for Levinson-Durbin recurtion-1

gamma_ABxPNLMS = 0.225;
% delta = 0.001;
error_tap=zeros(Ma,1);
loudspeaker_tap=zeros(Ma,1);
cancel_tap=zeros(L,1);
cancel_update_tap=zeros(L,1);

h_hat= zeros(L,1);
% H_hilbert=[4.54486433487e-05, -0.03191360456117,-3.31277815639e-05, -0.02604783588745,-1.594271701847e-06, -0.03732840213574,3.469425875731e-06, -0.05311074599375,-1.033785626931e-06, -0.07669953750135,-5.966834480336e-06,  -0.1168538267669,1.666246290771e-05,  -0.2058058100572,-1.893028029711e-05,  -0.6344726293921, 0,   0.6344726293921,1.893028029711e-05,   0.2058058100572,-1.666246290771e-05,   0.1168538267669,5.966834480336e-06,  0.07669953750135, 1.033785626931e-06,  0.05311074599375,-3.469425875731e-06,  0.03732840213574, 1.594271701847e-06,  0.02604783588745, 3.31277815639e-05,  0.03191360456117, -4.54486433487e-05]';
% Hilbert_length=length(H_hilbert);
% Hilbert_Delay=16;
% u_tap=zeros(Hilbert_length,1); 

for n=1:N
    
    h=aa_1;
    u(n)=F_amp*e_delay(n);%loudspeaker output

    %% Frequency Shifting
    
%     uu(n)=F_amp*e_delay(n);%loudspeaker output
%     
%     uu1(n+Hilbert_Delay)=uu(n);
%     uu1_u(n)=uu1(n)cos(2*pi*n-3);% Multiply Cos function with shift 3 Hz left side
%   
%     u_tap=[uu(n); u_tap(1:end-1)];
%     uu2(n)=u_tap'*H_hilbert;
%     uu2_u(n)=uu2(n)sin(2*pi*n-3);% Multiply Sin function with shift 3 Hz left side
%     u(n)=uu1_u(n)+uu2_u(n);% Frequency Shifting Operation
%% 
    feedback_tap=[u(n); feedback_tap(1:end-1)];
    feedback_cancel_tap=[u(n); feedback_cancel_tap(1:end-1)];
    f(n)=feedback_tap'*h;%output of feedback path-1
    
    m(n) = x(n) + f(n); %output of microphone-1
    f_cap(n)=feedback_cancel_tap'*h_cap;%output of feedback canceller-1
    e(n)=m(n)-f_cap(n); % subtract feedback canceller-1 output from microphone-1 output
    e_delay(n+delay_forw)=e(n); %insert delay in the forward path
    
    cancel_tap=[u(n); cancel_tap(1:end-1)]; %Input vector to  adaptive shadow filter-1
    f_hat(n)=cancel_tap'*h_hat; % output of adaptive shadow filter
    error(n)=m(n)-f_hat(n);%error used to update shadow filter-1
    
    tap=[e(n);tap(1:end-1,1)];
    [r,lg] = xcorr(tap,'biased');
    r(lg<0) = [];
    L1 = levinson(r,Ma-1); %Levinson-Durbin algorithm-1
    L_PEF=L1';
    
    error_tap=[error(n); error_tap(1:end-1)];
    error_prefilter(n)=error_tap'*L_PEF;%error signal prefiltered
    
    loudspeaker_tap=[u(n); loudspeaker_tap(1:end-1)];
    loudspeaker_prefiltered(n)=loudspeaker_tap'*L_PEF;%prefiltering of loudspeaker signal
    cancel_update_tap=[loudspeaker_prefiltered(n); cancel_update_tap(1:end-1)]; %Input vector to  adaptive shadow filter-1
       
    %% h_hat = h_hat + (mu_h  / (norm(cancel_update_tap)^2 + delta)) * cancel_update_tap .* error_prefilter(n);%update shadow filter-1
    
    %% ABxPNLMS for system update
    [w_sorted,index_sorted] = sort(h_hat,"descend","ComparisonMethod","abs");
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
         pp2 = ((1+gamma_ABxPNLMS)*g)/(2*s*P + delta);
         P_final = pp1+pp2;
         phi_ABxPNLMS=diag(P_final);

         %%Weight Update step

         % error_prefilter(i) = sys_out(i) - model_out_ABxPNLMS(i);
         h_hat = h_hat + mu_h*((phi_ABxPNLMS)*(cancel_update_tap)/((cancel_update_tap')*phi_ABxPNLMS*(cancel_update_tap) + delta)).*error_prefilter(n);
    %%
         h_cap=h_hat;
   
    
    HH_1_cap = fft( h_cap,Nfreq);
    
    diff = H_1(1:Nfreq/2+1) - HH_1_cap(1:Nfreq/2+1);% fft difference-1
    
    MSG(n)=20*log10 ( min (1./  abs(diff)   )  );% maximum stable gain
    
    ASG(n)=MSG(n)-( 20*log10( min( 1./ (  abs(H_1(1:Nfreq/2+1)))  )  )  );% Added stable gain
     
    MIS(n)=10*log10( ( (norm(diff)^2)  ) / ( (norm(H_1(1:Nfreq/2))^2)  ));% Misalignment

end
clc
iteration=itr

end
%% Performance Measures
fs=8000;  
t = linspace(0,(N-1)/fs,N);

 enhanced=u(end-80000:end);
 audiowrite('Loudspeaker.wav',enhanced,8000)
 clean=Input(end-80000:end);
 audiowrite('Mic.wav',clean,8000)
 PESQ_Score=pesq('Mic.wav','Loudspeaker.wav');
 
 

figure;
plot(t,MIS,'r','linewidth',2);grid on;
xlabel('Time (sec)');
ylabel('Misalignment (dB)');

figure;
plot(t,MSG,'r','linewidth',2);grid on;
xlabel('Time (sec)');
ylabel('MSG (dB)');

figure;
plot(t,ASG,'r','linewidth',2);grid on;
xlabel('Time (sec)');
ylabel('ASG (dB)');


figure;
plot(h);
hold on;
plot(h_cap,'r')

display('Simulation Completed')