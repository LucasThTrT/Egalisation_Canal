% Script for computing the BER for QAM modulation in ISI Channels with FDE
% 
close all;
clear all;

%% Simulation parameters
% On décrit ci-aprés les parametres généraux de la simulation

%modulation parameters
M = 4; %Modulation order 
Nframe = 10000;
Nfft=1024;
Ncp=8;
Ns=Nframe*(Nfft+Ncp);
N= log2(M)*Nframe*Nfft;

%Channel Parameters
Eb_N0_dB = [0:20]; % Eb/N0 

%Multipath channel parameters
hc=[1 -0.9];
Lc=length(hc);%Channel length
H=fft(hc,Nfft);

%Preallocations
nErr_zffde=zeros(1,length(Eb_N0_dB));
nErr_mmsefde=zeros(1,length(Eb_N0_dB));

for ii = 1:length(Eb_N0_dB)

   %Message generation
   bits= randi([0 1],N,1);
   s = qammod(bits,M,'InputType','bit');
   sigs2=var(s);
   
   % Add CP
   % On rajoute un préfixe cyclique (CP) à chaque bloc de données à
   % l'émetteur (voir rapportTP2.txt pour plus de précision)

   % ici on rajoute Ncp bits derrière chaque symbole OFDM
   % Afin d'éviter les intérférence inter-symboles
   smat=reshape(s,Nfft,Nframe);
   smatcp=[smat(end-Ncp+1:end,:);smat];
   scp=reshape(smatcp,1,(Nfft+Ncp)*Nframe);
   
    % Channel convolution: equivalent symbol based representation
   z = filter(hc,1,scp);  
   
   %Generating noise
   sig2b=10^(-Eb_N0_dB(ii)/10);
   
   %n = sqrt(sig2b)*randn(1,N+Lc-1); % white gaussian noise, 
   n = sqrt(sig2b/2)*randn(1,Ns)+1j*sqrt(sig2b/2)*randn(1,Ns); % white gaussian noise, 
   
    % Noise addition
   ycp = z + n; % additive white gaussian noise  

   %remove CP
   % On a ajouté dans CP Ncp symboles après chaque symbole OFDM
   % On va donc créer un matrice prenant uniquement les symboles sans leur
   % prolongement
   ycp = reshape(ycp, Nfft+Ncp,Nframe);
   ycp_reduit = ycp(Ncp + 1:end,:);
   
   % FDE
   % EGALISATION PAR FILTRAGE DANS LE DOMAINE FRÉQUENTIEL
   % On va compenser les efftes du canal de transmision qui peuvnt causer
   % des distorsioins et des interférences sur les données transmises.

   % On passe d'abord dans le domaine fréquentiel :
   % N-point DFT
   % Transormée de Fourier Discrète
   y = fft(ycp_reduit, Nfft);  % Transformée rapide de Fourier

   % Egalisateur RIF car fft donne une fenetre de taille fini


   % ZF
   % FDE par critères ZF
   % Calcul de W :
   w_zf = conj(H)./(abs(H).*abs(H));

   % On obtient alors notre signal après reception dans le filtre linéaire
   y_recu = y.*w_zf.';

   % On Repasse dans le domaine temporel pour pouvoir faire la prise de
   % décision 

   y_temporel = reshape(ifft(y_recu),1,Nfft*Nframe).';
   
   % Prise de décision
   bhat_zf = qamdemod(y_temporel,4,'OutputType','bit');


   % On rentre les informations d'erreurs
   nErr_zffde(1,ii) = size(find([bits(:)- bhat_zf(:)]),1);



   % MMSE
   % FDE par critères MMSE
   % Calcul de W :
   w_mmse = conj(H)./(abs(H).*abs(H)+sig2b/sigs2);

   % On obtient alors notre signal après reception dans le filtre linéaire
   y_recu = y.*w_mmse.';

   % On Repasse dans le domaine temporel pour pouvoir faire la prise de
   % décision 

   y_temporel = reshape(ifft(y_recu),1,Nfft*Nframe).';
   
   % Prise de décision
   bhat_mmse = qamdemod(y_temporel,4,'OutputType','bit');

   nErr_mmsefde(1,ii) = size(find([bits(:)- bhat_mmse(:)]),1);

end

simBer_zf = nErr_zffde/N; % simulated ber
simBer_mmse = nErr_mmsefde/N;

% plot

figure
semilogy(Eb_N0_dB,simBer_zf(1,:),'bs-','Linewidth',2);
hold on
semilogy(Eb_N0_dB,simBer_mmse(1,:),'rd-','Linewidth',2);
axis([0 70 10^-6 0.5])
grid on
legend('sim-zf-fde','sim-mmse-fde');
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');
title('Bit error probability curve for QPSK in ISI with ZF and MMSE equalizers')
