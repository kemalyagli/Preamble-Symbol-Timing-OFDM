%Zaman Domaininde Ren Preamble Üretir.

function [r,pn]=ren_SPreamble(nFFTSize) 

u=1;
a = zeros(1,nFFTSize/2);
    for n=0:nFFTSize/2-1
       a(n+1) =exp(2*sqrt(-1)*pi*u*n*n/nFFTSize);
    end
A = ifft(a);

pn=(-1).^round(rand(1,round(nFFTSize))); 


r=[A A].*pn;

% BPSKTable=[1,-1];
%     buf=BPSKTable(randi([0,1],nFFTSize/2,1)+1); 
%     x=zeros(1,nFFTSize); 
%     index = 1; 
%     for n=1:2:nFFTSize 
%          x(n)=buf(index); 
%          index=index+1; 
%     end
%     a = ifft(x);   %[A A]
end
