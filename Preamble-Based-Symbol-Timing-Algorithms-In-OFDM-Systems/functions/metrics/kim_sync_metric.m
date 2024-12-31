% Calculating the Kim Sync Metric

function M=kim_sync_metric(rm, nFFTSize)

for k=nFFTSize/2+2:length(rm)/40
        P3=0;
        R3=0;
    for j=1:nFFTSize/2
        P3(j)=rm(k-j-1+nFFTSize/2).*rm(k+j+nFFTSize/2);
        R3(j)=abs(rm(j+k+nFFTSize/2)).^2;
    end
        Pk(k+nFFTSize)=sum(P3);
        Rk(k+nFFTSize)=sum(R3);
end
        %Computing timing metric
        M=(((abs(Pk)).^2))./(Rk.^2);
end
