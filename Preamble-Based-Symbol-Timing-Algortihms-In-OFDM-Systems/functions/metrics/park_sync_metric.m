% Calculating the Park Sync Metric

function M=park_sync_metric(rm, nFFTSize)
    for k=nFFTSize/2+3:length(rm)/32
            P3=0;
            R3=0;
        for j=1:nFFTSize/2+1
            P3(j)=rm(k-j-1).*rm(k+j);
            R3(j)=abs(rm(j+k)).^2;
        end
            Pk(k+nFFTSize/2)=sum(P3);
            Rk(k+nFFTSize/2)=sum(R3);
    end
%Computing timing metric
M=(((abs(Pk)).^2))./(Rk.^2);
end
