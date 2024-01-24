clear all;
close all;
%Simulation parameters
 %
 %modulation parameters
 Es_N0_dB = [0:200]; % Eb/N0 values
 M = 4; %Modulation order
 Ns = 50000;
 N= log2(M)*Ns;
 %Viterbi decoding parameters
 const = qammod((0:M-1)',M) ; %reference Gray QPSK constellation
 tblen = 16; %Traceback depth
 nsamp = 1; %Oversampling rate
 preamble = [];
 postamble =[];
 %Channel Parameters
 hc = [0.623; 0.489+0.234i; 0.398i; 0.21];
 %Message generation
 %
 bits= randi([0 1] ,N,1) ;
 s = qammod(bits, M,'InputType','bit','UnitAveragePower',true) ;
 %Channel output
 z = filter(hc,1 ,s) ;
 nErr_mlinf=zeros(1,length(Es_N0_dB));
 for ii = 1:length(Es_N0_dB)
 %Generating noise
 sig2b=10^(-Es_N0_dB(ii)/10);
 %n = sqrt(sig2b)*randn(1,N+Lc-1); % white gaussian noise, BPSK Case
 n = sqrt(sig2b/2)*randn(1,N/2)+1j*sqrt(sig2b/2)*randn(1,N/2); % white gaussian noise, QPSK case
   
 % Adding Noise
 y = z + n'; % additive white gaussian noise %!!!!!!!

 s_ml = mlseeq(y,hc,const,tblen, 'rst' ,nsamp,[] ,[])' ;
 bhat_ml = zeros(2,N/2);
 bhat_ml(1,:)= -real(s_ml(1:N/2)) < 0;
 bhat_ml(2,:)= -imag(s_ml(1:N/2)) < 0;
 bhat_ml_vec = cat(1,bhat_ml(1,:),bhat_ml(2,:));
 nErr_mlinf(1,ii) = size(find([bits(:)- bhat_ml_vec(:)]),1);
 end
 simBer_mlinf = nErr_mlinf/N/log2(M); % simulated ber

% plot

figure
semilogy(Es_N0_dB,simBer_mlinf(1,:),'bs-','Linewidth',2);
axis([0 50 10^-6 0.5])
grid on
legend('sim-ml');
xlabel('E_s/N_0, dB');
ylabel('Bit Error Rate');
title('Bit error probability curve for QPSK in ISI with ML equalizer')
