% Generates Proposed Preamble in Time Domain

function S=proposed_preamble(u,N) 
    a = zeros(1,N/4);
    for n=0:N/4-1
       a(n+1) =exp(4*sqrt(-1)*pi*u*(n)*(n)/N);
    end
       A = a;%ifft(a);
       S = [A A circshift(A,N/8) circshift(A,N/8)];
end
