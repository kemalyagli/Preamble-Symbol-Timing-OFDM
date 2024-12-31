clear variables; close all; clc
addpath("functions/metrics/");
addpath("functions/preambles/");

% Variables
nFFT = 256; % FFT length
nDSC = 208; % number of subcarriers
nBitPerSym = 208; % bits per OFDM symbol
nSym = 80; % number of symbols 
cpLength = 16;
u=1; % sequence index for CAZAC preambles

%SNR
ebN0db = 0:5:20; % bit to noise ratio
esN0db = ebN0db + 10*log10(nDSC/nFFT) + 10*log10(nFFT/(nFFT+cpLength)); % converting to symbol to noise ratio

% Predefined arrays for storing results
schmidlMeanAwgn  = zeros(1,length(esN0db));
schmidlMeanSquareAwgn=zeros(1, length(ebN0db));
minnMeanAwgn  = zeros(1,length(esN0db));
minnMeanSquareAwgn=zeros(1, length(ebN0db));
parkMeanAwgn  = zeros(1,length(esN0db));
parkMeanSquareAwgn=zeros(1, length(ebN0db));
renMeanAwgn  = zeros(1,length(esN0db));
renMeanSquareAwgn=zeros(1, length(ebN0db));
kimMeanAwgn  = zeros(1,length(esN0db));
kimMeanSquareAwgn=zeros(1, length(ebN0db));
proposedMeanAwgn  = zeros(1,length(esN0db));
proposedMeanSquareAwgn=zeros(1, length(ebN0db));

schmidlMeanRayleigh  = zeros(1,length(esN0db));
schmidlMeanSquareRayleigh=zeros(1, length(ebN0db));
minnMeanRayleigh  = zeros(1,length(esN0db));
minnMeanSquareRayleigh=zeros(1, length(ebN0db));
parkMeanRayleigh  = zeros(1,length(esN0db));
parkMeanSquareRayleigh=zeros(1, length(ebN0db));
renMeanRayleigh  = zeros(1,length(esN0db));
renMeanSquareRayleigh=zeros(1, length(ebN0db));
kimMeanRayleigh  = zeros(1,length(esN0db));
kimMeanSquareRayleigh=zeros(1, length(ebN0db));
proposedMeanRayleigh  = zeros(1,length(esN0db));
proposedMeanSquareRayleigh=zeros(1, length(ebN0db));

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
for ii = 1:length(ebN0db)
    
sigma = 10^(-esN0db(ii)/20);

   for iteration=1:10000   % Transmitter
   ipBit = rand(1,nBitPerSym*nSym) > 0.5; % random 1's and 0's
   ipMod = 2*ipBit-1; % BPSK modulation 0 --> -1, 1 --> +1
   ipMod = reshape(ipMod,nBitPerSym,nSym).'; % grouping into multiple symbols

   % Assigning modulated symbols to subcarriers from [-26 to -1, +1 to +26]
   xF = [zeros(nSym,(nFFT-nDSC)/2) ipMod(:,[1:nBitPerSym/2]) zeros(nSym,1) ipMod(:,[nBitPerSym/2+1:nBitPerSym]) zeros(nSym,(nFFT-nDSC)/2-1)] ;
   
   % Taking FFT, the term (nFFT/sqrt(nDSC)) is for normalizing the power of transmit symbol to 1 
   xt = (nFFT/sqrt(nDSC))*ifft(fftshift(xF.')).';

   % Appending cylic prefix
   xt = [xt(:,[nFFT-cpLength+1:nFFT]) xt];

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

   % Multipath channel simulation
   nTap = 10;
   T=0:nTap-1;
%  ht_normal = 1/sqrt(2)*1/sqrt(nTap).*(randn(nSym+1,nTap) +
%  1i*randn(nSym+1,nTap)); % constant channel coefficient
   ht_exp = exp(-T/nTap)/sqrt(2)*1/sqrt(nTap).*(randn(nSym+1,nTap) + 1i*randn(nSym+1,nTap)); % exponantially delaying channel coefficient
   ht_exp = ht_exp*sqrt(1/sum(mean(abs(ht_exp).^2)));
    
   % cp-ofdm symbol
   xtPreambleSchmidlRayleigh = [schmidlCp schmidlPreamble xt];
   xtPreambleMinnRayleigh = [minnCP minnPreamble xt];
   xtPreambleParkRayleigh = [parkCP parkPreamble xt];
   xtPreambleRenRayleigh = [renCP renPreamble xt];
   xtPreambleKimRayleigh = [kimCP kimPreamble xt];
   xtPreambleProposedRayleigh = [proposedCP proposedPreamble xt];
   
   % formatting the received vector into symbols
   xtPreambleSchmidlRayleighSym = reshape(xtPreambleSchmidlRayleigh.',nFFT+cpLength,nSym+1).'; 
   xtPreambleMinnRayleighSym = reshape(xtPreambleMinnRayleigh.',nFFT+cpLength,nSym+1).'; 
   xtPreambleParkRayleighSym = reshape(xtPreambleParkRayleigh.',nFFT+cpLength,nSym+1).'; 
   xtPreambleRenRayleighSym = reshape(xtPreambleRenRayleigh.',nFFT+cpLength,nSym+1).'; 
   xtPreambleKimRayleighSym = reshape(xtPreambleKimRayleigh.',nFFT+cpLength,nSym+1).';
   xtPreambleProposedRayleighSym = reshape(xtPreambleProposedRayleigh.',nFFT+cpLength,nSym+1).';
  
   % convolution of each symbol with the random channel
   for jj = 1:nSym+1
      xhtPreambleSchmidlRayleighSym(jj,:) = conv(ht_exp(jj,:),xtPreambleSchmidlRayleighSym(jj,:));
      xhtPreambleMinnRayleighSym(jj,:) = conv(ht_exp(jj,:),xtPreambleMinnRayleighSym(jj,:));
      xhtPreambleParkRayleighSym(jj,:) = conv(ht_exp(jj,:),xtPreambleParkRayleighSym(jj,:));
      xhtPreambleRenRayleighSym(jj,:) = conv(ht_exp(jj,:),xtPreambleRenRayleighSym(jj,:));
      xhtPreambleKimRayleighSym(jj,:) = conv(ht_exp(jj,:),xtPreambleKimRayleighSym(jj,:));
      xhtPreambleProposedRayleighSym(jj,:) = conv(ht_exp(jj,:),xtPreambleProposedRayleighSym(jj,:));
   end

   % Concatenating multiple symbols to form a long vector
   xhtPreambleSchmidl = reshape(xhtPreambleSchmidlRayleighSym.',1,(nSym+1)*(nFFT+cpLength+nTap-1));
   xhtPreambleMinn = reshape(xhtPreambleMinnRayleighSym.',1,(nSym+1)*(nFFT+cpLength+nTap-1));
   xhtPreamblePark = reshape(xhtPreambleParkRayleighSym.',1,(nSym+1)*(nFFT+cpLength+nTap-1));
   xhtPreambleRen = reshape(xhtPreambleRenRayleighSym.',1,(nSym+1)*(nFFT+cpLength+nTap-1));
   xhtPreambleKim = reshape(xhtPreambleKimRayleighSym.',1,(nSym+1)*(nFFT+cpLength+nTap-1));
   xhtPreambleProposed = reshape(xhtPreambleProposedRayleighSym.',1,(nSym+1)*(nFFT+cpLength+nTap-1));

   nt_r = 1/sqrt(2)*(randn(1,length(xhtPreambleSchmidl)) + 1i*randn(1,length(xhtPreambleSchmidl)));

   ytPreambleSchmidl = xhtPreambleSchmidl + 10^(-esN0db(ii)/20)*nt_r;
   ytPreambleMinn = xhtPreambleMinn + 10^(-esN0db(ii)/20)*nt_r;
   ytPreamblePark = xhtPreamblePark + 10^(-esN0db(ii)/20)*nt_r;
   ytPreambleRen = xhtPreambleRen + 10^(-esN0db(ii)/20)*nt_r;
   ytPreambleKim = xhtPreambleKim + 10^(-esN0db(ii)/20)*nt_r;
   ytPreambleProposed = xhtPreambleProposed + 10^(-esN0db(ii)/20)*nt_r;

   ytPreambleSchmidl2 = [ytPreambleSchmidl(end-len+1:end) ytPreambleSchmidl];% xt cp-ofdm sembolü
   ytPreambleMinn2 = [ytPreambleSchmidl(end-len+1:end) ytPreambleMinn];
   ytPreamblePark2 = [ytPreambleSchmidl(end-len+1:end) ytPreamblePark];
   ytPreambleRen2 = [ytPreambleSchmidl(end-len+1:end) ytPreambleRen];
   ytPreambleKim2 = [ytPreambleSchmidl(end-len+1:end) ytPreambleKim];
   ytPreambleProposed2 = [ytPreambleSchmidl(end-len+1:end) ytPreambleProposed];



   msRayleigh=schmidl_sync_metric(ytPreambleSchmidl2,nFFT);
   [metricSchmidlRayleigh, sbnEstSchmidlRayleigh] = max(msRayleigh(sbn-nFFT:sbn+nFFT));
   sbnEstSchmidlRayleigh=sbnEstSchmidlRayleigh+nFFT+2*cpLength;
   pbnEstSchmidlRayleigh = sbnEstSchmidlRayleigh-len;  

   mmRayleigh=minn_sync_metric(ytPreambleMinn2,nFFT);
   [metricMinnRayleigh, sbnEstMinnRayleigh] = max(mmRayleigh(sbn-nFFT:sbn+nFFT));
   sbnEstMinnRayleigh=sbnEstMinnRayleigh+nFFT+2*cpLength;
   pbnEstMinnRayleigh = sbnEstMinnRayleigh-len; 

   mpRayleigh=park_sync_metric(ytPreamblePark2,nFFT);
   [metricParkRayleigh, sbnEstParkRayleigh] = max(mpRayleigh(sbn-nFFT:sbn+nFFT));
   sbnEstParkRayleigh=sbnEstParkRayleigh+nFFT+2*cpLength;
   pbnEstParkRayleigh = sbnEstParkRayleigh-len;  

   mrRayleigh=ren_sync_metric(ytPreambleRen2,nFFT, pn);
   [metricRenRayleigh, sbnEstRenRayleigh] = max(mrRayleigh(sbn-nFFT:sbn+nFFT));
   sbnEstRenRayleigh=sbnEstRenRayleigh+nFFT+2*cpLength;
   pbnEstRenRayleigh = sbnEstRenRayleigh-len; 

   mkRayleigh=kim_sync_metric(ytPreambleKim2,nFFT);
   [metricKimRayleigh, sbnEstKimRayleigh] = max(mkRayleigh(sbn-nFFT:sbn+nFFT));
   sbnEstKimRayleigh=sbnEstKimRayleigh+nFFT+2*cpLength;
   pbnEstKimRayleigh = sbnEstKimRayleigh-len; 

   mproRayleigh=minn_sync_metric(ytPreambleProposed2,nFFT);
   [metricProRayleigh, sbnEstProRayleigh] = max(mproRayleigh(sbn-nFFT:sbn+nFFT));
   sbnEstProRayleigh=sbnEstProRayleigh+nFFT+2*cpLength;
   pbnEstProRayleigh = sbnEstProRayleigh-len; 

   %MSE
   schmidlMeanRayleigh(ii) = schmidlMeanRayleigh(ii) + (sbn - sbnEstSchmidlRayleigh);
   schmidlMeanSquareRayleigh(ii) = schmidlMeanSquareRayleigh(ii) + power((sbn - sbnEstSchmidlRayleigh),2);
   minnMeanRayleigh(ii) = minnMeanRayleigh(ii) + (sbn - sbnEstMinnRayleigh);
   minnMeanSquareRayleigh(ii) = minnMeanSquareRayleigh(ii) + power((sbn - sbnEstMinnRayleigh),2); 
   parkMeanRayleigh(ii) = parkMeanRayleigh(ii) + (sbn - sbnEstParkRayleigh);
   parkMeanSquareRayleigh(ii) = parkMeanSquareRayleigh(ii) + power((sbn - sbnEstParkRayleigh),2); 
   renMeanRayleigh(ii) = renMeanRayleigh(ii) + (sbn - sbnEstRenRayleigh);
   renMeanSquareRayleigh(ii) = renMeanSquareRayleigh(ii) + power((sbn - sbnEstRenRayleigh),2);
   kimMeanRayleigh(ii) = kimMeanRayleigh(ii) + (sbn - sbnEstKimRayleigh);
   kimMeanSquareRayleigh(ii) = kimMeanSquareRayleigh(ii) + power((sbn - sbnEstKimRayleigh),2);
   proposedMeanRayleigh(ii) = proposedMeanRayleigh(ii) + (sbn - sbnEstProRayleigh);
   proposedMeanSquareRayleigh(ii) = proposedMeanSquareRayleigh(ii) + power((sbn - sbnEstProRayleigh),2);

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

schmidlMeanSquareRayleigh = schmidlMeanSquareRayleigh/iteration;
schmidlMeanRayleigh = schmidlMeanRayleigh/iteration;
minnMeanSquareRayleigh = minnMeanSquareRayleigh/iteration;
minnMeanRayleigh = minnMeanRayleigh/iteration;
parkMeanSquareRayleigh = parkMeanSquareRayleigh/iteration;
parkMeanRayleigh = parkMeanRayleigh/iteration;
renMeanSquareRayleigh = renMeanSquareRayleigh/iteration;
renMeanRayleigh = renMeanRayleigh/iteration;
kimMeanSquareRayleigh = kimMeanSquareRayleigh/iteration;
kimMeanRayleigh = kimMeanRayleigh/iteration;
proposedMeanSquareRayleigh = proposedMeanSquareRayleigh/iteration;
proposedMeanRayleigh = proposedMeanRayleigh/iteration;



