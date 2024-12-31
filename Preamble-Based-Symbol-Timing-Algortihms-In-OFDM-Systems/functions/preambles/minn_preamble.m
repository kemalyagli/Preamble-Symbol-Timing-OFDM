% Generates Minn Preamble in Time Domain

function c=minn_preamble(nFFTSize) 
    BPSKTable=[1,-1];
    buf=BPSKTable(randi([0,1],nFFTSize/2,1)+1); 
    x=zeros(1,nFFTSize/2); 
    index = 1; 
    for n=1:2:nFFTSize/2 
         x(n)=buf(index); 
         index=index+1; 
    end
    a = ifft(x); 
    c= [a -a]; %[A A -A -A]
end