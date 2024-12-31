% Generates Park Preamble in Time Domain

function d=park_preamble(nFFTSize) 
    BPSKTable=[1,-1];
    buf=BPSKTable(randi([0,1],nFFTSize/2,1)+1); 
    x=zeros(1,nFFTSize/2); 
    index = 1; 
    for n=1:2:nFFTSize/2 
         x(n)=buf(index); 
         index=index+1; 
    end
    a = ifft(x); 
    d= [a(1:nFFTSize/4) fliplr(a(1:nFFTSize/4)) conj(a(1:nFFTSize/4)) conj(fliplr(a(1:nFFTSize/4)))];%[A B A* B*]
end