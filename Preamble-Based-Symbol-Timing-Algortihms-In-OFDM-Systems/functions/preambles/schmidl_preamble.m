% Generates Schmidl Preamble in Time Domain

function a=schmidl_preamble(nFFTSize) 
%IEEE 802.11 standard
% for each symbol bits a1 to a52 are assigned to subcarrier 
% index [-26 to -1 1 to 26] 
% A short OFDM training symbol consists of 12 subcarriers
inputFFTShortPreamble = [zeros(1,8) 1+1i 0 0 0  -1-1i 0 0 0 ... % [-32:-17]
 1+1i 0 0 0  -1-1i 0 0 0 -1-1i 0 0 0   1+1i 0 0 0 ...             % [-16:-1]
 0   0 0 0  -1-1i 0 0 0 -1-1i 0 0 0   1+1i 0 0 0 ...             % [0:15]
 1+1i 0 0 0   1+1i 0 0 0  1+1i 0 0 0   0   0 0 0 ];              % [16:31]
inputiFFTShortPreamble = sqrt(13/6)*fftshift(inputFFTShortPreamble); 
% taking ifft
outputiFFTShortPreamble = ifft(inputiFFTShortPreamble,nFFTSize); % generate 64 sample sequence
a=outputiFFTShortPreamble;

BPSKTable=[1,-1];
    buf=BPSKTable(randi([0,1],nFFTSize/2,1)+1); 
    x=zeros(1,nFFTSize); 
    index = 1; 
    for n=1:2:nFFTSize 
         x(n)=buf(index); 
         index=index+1; 
    end
    a = ifft(x);   %[A A]
end
