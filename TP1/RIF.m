% Script for computing the BER for BPSK/QPSK modulation in ISI Channels
% 
close all;
clear all;

%% Simulation parameters
% On décrit ci-après les paramètres généraux de la simulation

%Frame length
M=4; %2:BPSK, 4: QPSK
N  = 1000000; % Number of transmitted bits or symbols
Es_N0_dB = [0:100]; % Eb/N0 values
%Multipath channel parameters
hc=[1 0.8*exp(1i*pi/3) 0.3*exp(1i*pi/6) ];%0.1*exp(1i*pi/12)];%ISI channel
%hc=[0.04, -0.05, 0.07, -0.21, -0.5, 0.72, 0.36, 0, 0.21, 0.03, 0.07];
%a=1.2;
%hc=[1 -a];
Lc=length(hc);%Channel length
ChannelDelay=0; %delay is equal to number of non causal taps

%Preallocations
nErr_zfinf=zeros(1,length(Es_N0_dB));
for ii = 1:length(Es_N0_dB)

   % BPSK symbol generations
%    bits = rand(1,N)>0.5; % generating 0,1 with equal probability
%    s = 1-2*bits; % BPSK modulation following: {0 -> +1; 1 -> -1} 
   
    % QPSK symbol generations
   bits = rand(2,N)>0.5; % generating 0,1 with equal probability
   s = 1/sqrt(2)*((1-2*bits(1,:))+1j*(1-2*bits(2,:))); % QPSK modulation following the BPSK rule for each quadatrure component: {0 -> +1; 1 -> -1} 
   sigs2=var(s);
   
   % Channel convolution: equivalent symbol based representation
   z = conv(hc,s);  
   
   %Generating noise
   sig2b=10^(-Es_N0_dB(ii)/10);
   %n = sqrt(sig2b)*randn(1,N+Lc-1); % white gaussian noise, BPSK Case
    n = sqrt(sig2b/2)*randn(1,N+Lc-1)+1j*sqrt(sig2b/2)*randn(1,N+Lc-1); % white gaussian noise, QPSK case
   
   % Adding Noise
   y = z + n; % additive white gaussian noise

   %% zero forcing equalization
   % We now study ZF equalization
   
   %Unconstrained ZF equalization, only if stable inverse filtering
   
   
   %%
   % 
   %  The unconstrained ZF equalizer, when existing is given by 
   % 
   % $w_{,\infty,zf}=\frac{1}{h(z)}$
   % 
   %%
    s_zf=filter(1,hc,y);%if stable causal filter is existing
    bhat_zf = zeros(2,length(bits));
    bhat_zf(1,:)= real(s_zf(1:N)) < 0;
    bhat_zf(2,:)= imag(s_zf(1:N)) < 0;
    nErr_zfinfdirectimp(1,ii) = size(find([bits(:)- bhat_zf(:)]),1);
    %Otherwise, to handle the non causal case
    Nzf=200; % Indice ou on s'arrete d'approximer ??
    [r, p, k]=residuez(1, hc);
    [w_zfinf]=ComputeRI( Nzf, r, p, k );
    s_zf=conv(w_zfinf,y);
    bhat_zf = zeros(2,length(bits));
    bhat_zf(1,:)= real(s_zf(Nzf:N+Nzf-1)) < 0;
    bhat_zf(2,:)= imag(s_zf(Nzf:N+Nzf-1)) < 0;
    
    nErr_zfinf(1,ii) = size(find([bits(:)- bhat_zf(:)]),1);



   %%MMSE
    deltac= zeros(1,2*Lc-1) ;
    deltac(Lc) = 1 ;
    Nmmse = 200 ;%causal part
    [r,p,k]=residuez(fliplr(conj(hc)),(conv(hc,fliplr(conj(hc)))+(sig2b/sigs2)*deltac));
    [w_mmseinf]=ComputeRI(Nmmse,r,p,k) ;
    s_zf=conv(w_mmseinf,y);
    bhat_mmse = zeros(2,length(bits));
    bhat_mmse(1,:)= real(s_zf(Nmmse:N+Nmmse-1)) < 0;
    bhat_mmse(2,:)= imag(s_zf(Nmmse:N+Nmmse-1)) < 0;
    nErr_zfmmse(1,ii) = size(find([bits(:)- bhat_mmse(:)]),1);

    %%RIFZF

    Nw = 10;
    H=toeplitz([hc(1) zeros(1,Nw-1)].',[hc,zeros(1,Nw-1)]).';
    Hd = inv(H'*H)*H';
    P = H*Hd;
    [~,D] = max(diag(abs(P)));
    w_RIP = Hd(:,D);
    s_rip=conv(w_RIP,y);
    s_rip = s_rip(D:end); %astuce du prof
    bhat_rip = zeros(2,length(bits));
    bhat_rip(1,:)= real(s_rip(1:end-(Nw+1))) < 0; %indices ?
    bhat_rip(2,:)= imag(s_rip(1:end-(Nw+1))) < 0; %indices ?
    nErr_rip(1,ii) = size(find([bits(:)- bhat_rip(:)]),1);
  
    %%RIFMMSE

    %deltac= zeros(1,2*Lc-1) ;
    %deltac(Lc) = 1 ;
    %[r,p,k]=residuez(fliplr(conj(hc)),(conv(hc,fliplr(conj(hc)))+(sig2b/sigs2)*deltac));
    %w_mmserif= w_RIP ;
    %s_mmserif=conv(w_mmserif,y);
    %bhat_mmserif = zeros(2,length(bits));
    %bhat_mmserif(1,:)= real(s_mmserif(1:end-(Nw+1))) < 0;
    %bhat_mmserif(2,:)= imag(s_mmserif(1:end-(Nw+1))) < 0;
    %nErr_mmserif(1,ii) = size(find([bits(:)- bhat_mmserif(:)]),1);

end
simBer_rip = nErr_rip/N/log2(M); % simulated ber
%simBer_mmserif = nErr_mmserif/N/log2(M); % simulated ber
simBer_zfinfmmse = nErr_zfmmse/N/log2(M); % simulated ber
simBer_zfinfdirectimp = nErr_zfinfdirectimp/N/log2(M); % simulated ber
simBer_zfinf = nErr_zfinf/N/log2(M); % simulated ber
% plot

figure
semilogy(Es_N0_dB,simBer_zfinfdirectimp(1,:),'bs-','Linewidth',2);
hold on
semilogy(Es_N0_dB,simBer_zfinf(1,:),'rs-','Linewidth',2);
hold on
semilogy(Es_N0_dB,simBer_zfinfmmse(1,:),'ks-','Linewidth',2);
hold on
semilogy(Es_N0_dB,simBer_rip(1,:),'ys-','Linewidth',2);
%hold on
%semilogy(Es_N0_dB,simBer_mmserif(1,:)-0.0001,'gs-','Linewidth',2);
axis([0 50 10^-6 0.5])
grid on
legend('sim-zf-direct','sim-zf-inf','sim-mmse-inf','sim-rif');
xlabel('E_s/N_0, dB');
ylabel('Bit Error Rate');
title('Bit error probability curve for QPSK in ISI with ZF MMSE and RIF equalizers')

figure
title('Impulse response')
stem(real(w_mmseinf))
hold on
stem(real(w_mmseinf),'r-')
ylabel('Amplitude');
xlabel('time index')



