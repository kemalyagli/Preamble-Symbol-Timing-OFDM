% Calculating the Schmidl Sync Metric

function [M,Pk]=schmidl_sync_metric(rm, nFFTSize) 
for k=1:length(rm)/40
         Pm=0; 
         Rm=0; 
     for j=1:nFFTSize/2
         Pm(j)=conj(rm(j+k)).*rm(j+k+nFFTSize/2); 
         Rm(j)=conj(rm(j+k+nFFTSize/2)).*rm(j+k+nFFTSize/2); 
     end 
         Pk(k+nFFTSize+1)=sum(Pm); 
         Rk(k+nFFTSize+1)=sum(Rm); 
end 
        %Computing the metric 
        M=(((abs(Pk)).^2))./(Rk.^2);
end


