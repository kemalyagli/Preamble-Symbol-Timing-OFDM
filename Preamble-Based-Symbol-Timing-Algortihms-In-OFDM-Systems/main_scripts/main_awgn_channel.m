clear variables; close all; clc
addpath("functions/metrics/");
addpath("functions/preambles/");

% Variables
nFFT = 256; % FFT length
nDSC = 208; % number of subcarriers
nBitPerSym = 208; % bits per OFDM symbol
nSym = 80; % number of symbols
cpLength = 32; % cyclic prefix length
u = 1; % sequence index for CAZAC preambles

% SNR
ebN0dB = 0:2:20; % bit to noise ratio
esN0dB = ebN0dB + 10*log10(nDSC/nFFT) + 10*log10(nFFT/(nFFT+cpLength)); % symbol to noise ratio conversion

% Predefined arrays for storing results
schmidlMeanAwgn = zeros(1,length(esN0dB));
schmidlMeanSquareAwgn = zeros(1, length(ebN0dB));
minnMeanAwgn = zeros(1,length(esN0dB));
minnMeanSquareAwgn = zeros(1, length(ebN0dB));
parkMeanAwgn = zeros(1,length(esN0dB));
parkMeanSquareAwgn = zeros(1, length(ebN0dB));
renMeanAwgn = zeros(1,length(esN0dB));
renMeanSquareAwgn = zeros(1, length(ebN0dB));
kimMeanAwgn = zeros(1,length(esN0dB));
kimMeanSquareAwgn = zeros(1, length(ebN0dB));
proposedMeanAwgn  = zeros(1,length(esN0dB));
proposedMeanSquareAwgn = zeros(1, length(ebN0dB));


% Preamble energy calculations
es = mean(abs(schmidl_preamble(nFFT)).^2);
em = mean(abs(minn_preamble(nFFT)).^2);
ep = mean(abs(park_preamble(nFFT)).^2);
er = mean(abs(ren_preamble(nFFT)).^2);
ek = mean(abs(kim_preamble(nFFT)).^2);
epro = mean(abs(proposed_preamble(u,nFFT)).^2);

% Normalize preambles
schmidlPreamble = sqrt(1/es)*schmidl_preamble(nFFT);
schmidlCp = schmidlPreamble(length(schmidlPreamble)-cpLength+1:end);
minnPreamble = sqrt(1/em)*minn_preamble(nFFT);
minnCP = minnPreamble(length(minnPreamble)-cpLength+1:end);
parkPreamble = sqrt(1/ep)*park_preamble(nFFT);
parkCP = parkPreamble(length(parkPreamble)-cpLength+1:end);
[renPreamble, pn] = ren_preamble(nFFT);
renPreamble = sqrt(1/er)*renPreamble;
renCP = renPreamble(length(renPreamble)-cpLength+1:end);
kimPreamble = sqrt(1/ek)*kim_preamble(nFFT);
kimCP = kimPreamble(length(kimPreamble)-cpLength+1:end);
proposedPreamble = sqrt(1/epro)*proposed_preamble(u,nFFT);
proposedCP = proposedPreamble(length(proposedPreamble)-cpLength+1:end);

len = length([schmidlCp schmidlPreamble]);

% Simulation loop for different SNR values
for ii = 1:length(ebN0dB) 
   
   sigma = 10^(-esN0dB(ii)/20);

   for iteration=1:10000  % Transmitter
   ipBit = rand(1,nBitPerSym*nSym) > 0.5; % random 1's and 0's
   ipMod = 2*ipBit-1; % BPSK modulation 0 --> -1, 1 --> +1
   ipMod = reshape(ipMod,nBitPerSym,nSym).'; % grouping into multiple symbols

   % Assigning modulated symbols to subcarriers from [-26 to -1, +1 to +26]
   xF = [zeros(nSym,(nFFT-nDSC)/2) ipMod(:,[1:nBitPerSym/2]) zeros(nSym,1) ipMod(:,[nBitPerSym/2+1:nBitPerSym]) zeros(nSym,(nFFT-nDSC)/2-1)] ;
   
   % Taking FFT, the term (nFFT/sqrt(nDSC)) is for normalizing the power of transmit symbol to 1 
   xt = (nFFT/sqrt(nDSC))*ifft(fftshift(xF.')).';

   % Appending cylic prefix
   xt = [xt(:,nFFT-cpLength+1:nFFT) xt];

   % Concatenating multiple symbols to form a long vector
   xt = reshape(xt.',1,nSym*(nFFT+cpLength));%CP eklenmiÅŸ 280 sembol vektÃ¶r formunda(1 frame)
   
   % Preamble starting point
   pbn=len+1;
   % Symbol starting point
   sbn=pbn+len;

   % AWGN Channel Simulation
   xtPreambleSchmidlAwgn = [xt(1:len) schmidlCp schmidlPreamble xt(len+1:end)];
   xtPreambleMinnAwgn = [xt(1:len) minnCP minnPreamble xt(len+1:end)];
   xtPreambleParkAwgn = [xt(1:len) parkCP parkPreamble xt(len+1:end)];
   xtPreambleRenAwgn = [xt(1:len) renCP renPreamble xt(len+1:end)];
   xtPreambleKimAwgn = [xt(1:len) kimCP kimPreamble xt(len+1:end)];
   xtPreambleProposedAwgn = [xt(1:len) proposedCP proposedPreamble xt(len+1:end)];

   % Generate Gaussian noise with unit variance and zero mean
   nt = 1/sqrt(2)*(randn(1,length(xtPreambleSchmidlAwgn)) + 1i*randn(1,length(xtPreambleSchmidlAwgn)));

   % Add noise to the transmitted signal
   ytSchmidlAwgn = xtPreambleSchmidlAwgn + sigma*nt;
   ytMinnAwgn = xtPreambleMinnAwgn + sigma*nt;
   ytParkAwgn = xtPreambleParkAwgn + sigma*nt;
   ytRenAwgn = xtPreambleRenAwgn + sigma*nt;
   ytKimAwgn = xtPreambleKimAwgn + sigma*nt;
   ytProposedAwgn = xtPreambleProposedAwgn + sigma*nt;

   % Synchronization metrics for each preamble
   msAwgn=schmidl_sync_metric(ytSchmidlAwgn,nFFT);
   [metricSchmidlAwgn, sbnEstSchmidlAwgn] = max(msAwgn(sbn-nFFT:sbn+nFFT));
   sbnEstSchmidlAwgn=sbnEstSchmidlAwgn+nFFT+2*cpLength;
   pbnEstSchmidlAwgn = sbnEstSchmidlAwgn-len;

   mmAwgn=minn_sync_metric(ytMinnAwgn,nFFT);
   [metricMinnAwgn, sbnEstMinnAwgn] = max(mmAwgn(sbn-nFFT:sbn+nFFT));
   sbnEstMinnAwgn=sbnEstMinnAwgn+nFFT+2*cpLength;
   pbnEstMinnAwgn = sbnEstMinnAwgn-len;
   
   mpAwgn=park_sync_metric(ytParkAwgn,nFFT);
   [metricParkAwgn, sbnEstParkAwgn] = max(mpAwgn(sbn-nFFT:sbn+nFFT));
   sbnEstParkAwgn=sbnEstParkAwgn+nFFT+2*cpLength;
   pbnEstParkAwgn = sbnEstParkAwgn-len;
    

   mrAwgn=ren_sync_metric(ytRenAwgn,nFFT, pn);
   [metricRenAwgn, sbnEstRenAwgn] = max(mrAwgn(sbn-nFFT:sbn+nFFT));
   sbnEstRenAwgn=sbnEstRenAwgn+nFFT+2*cpLength;
   pbnEstRenAwgn = sbnEstRenAwgn-len; 

   mkAwgn=kim_sync_metric(ytKimAwgn,nFFT);
   [metricKimAwgn, sbnEstKimAwgn] = max(mkAwgn(sbn-nFFT:sbn+nFFT));
   sbnEstKimAwgn=sbnEstKimAwgn+nFFT+2*cpLength;
   pbnEstKimAwgn = sbnEstKimAwgn-len;

   mproAwgn=minn_sync_metric(ytProposedAwgn,nFFT);
   [metricProAwgn, sbnEstProAwgn] = max(mproAwgn(sbn-nFFT:sbn+nFFT));
   sbnEstProAwgn=sbnEstProAwgn+nFFT+2*cpLength;
   pbnEstProAwgn = sbnEstProAwgn-len;

   %MSE
   schmidlMeanAwgn(ii) = schmidlMeanAwgn(ii) + (sbn - sbnEstSchmidlAwgn);
   schmidlMeanSquareAwgn(ii) = schmidlMeanSquareAwgn(ii) + power((sbn - sbnEstSchmidlAwgn),2);
   minnMeanAwgn(ii) = minnMeanAwgn(ii) + (sbn - sbnEstMinnAwgn);
   minnMeanSquareAwgn(ii) = minnMeanSquareAwgn(ii) + power((sbn - sbnEstMinnAwgn),2); 
   parkMeanAwgn(ii) = parkMeanAwgn(ii) + (sbn - sbnEstParkAwgn);
   parkMeanSquareAwgn(ii) = parkMeanSquareAwgn(ii) + power((sbn - sbnEstParkAwgn),2);  
   renMeanAwgn(ii) = renMeanAwgn(ii) + (sbn - sbnEstRenAwgn);
   renMeanSquareAwgn(ii) = renMeanSquareAwgn(ii) + power((sbn - sbnEstRenAwgn),2);
   kimMeanAwgn(ii) = kimMeanAwgn(ii) + (sbn - sbnEstKimAwgn);
   kimMeanSquareAwgn(ii) = kimMeanSquareAwgn(ii) + power((sbn - sbnEstKimAwgn),2);
   proposedMeanAwgn(ii) = proposedMeanAwgn(ii) + (sbn - sbnEstProAwgn);
   proposedMeanSquareAwgn(ii) = proposedMeanSquareAwgn(ii) + power((sbn - sbnEstProAwgn),2);

   end
   
end


schmidlMeanSquareAwgn = schmidlMeanSquareAwgn/iteration;
schmidlMeanAwgn = schmidlMeanAwgn/iteration;
minnMeanSquareAwgn = minnMeanSquareAwgn/iteration;
minnMeanAwgn = minnMeanAwgn/iteration;
parkMeanSquareAwgn = parkMeanSquareAwgn/iteration;
parkMeanAwgn = parkMeanAwgn/iteration;
renMeanSquareAwgn = renMeanSquareAwgn/iteration;
renMeanAwgn = renMeanAwgn/iteration;
kimMeanSquareAwgn = kimMeanSquareAwgn/iteration;
kimMeanAwgn = kimMeanAwgn/iteration;
proposedMeanSquareAwgn = proposedMeanSquareAwgn/iteration;
proposedMeanAwgn = proposedMeanAwgn/iteration;







