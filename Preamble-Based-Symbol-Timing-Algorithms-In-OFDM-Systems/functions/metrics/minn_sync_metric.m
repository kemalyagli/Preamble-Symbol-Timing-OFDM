% Calculating the Minn Sync Metric

function M=minn_sync_metric(rm, nFFTSize) 
    for k=1:length(rm)/40
     Pm1=0; 
     Rm1=0; 
     Rm2=0;
     Pm2=0; 
         for j=1:nFFTSize/4
         Pm1(j)=conj(rm(j+k)).*rm(k+j+nFFTSize/4); 
         Pm2(j)=conj(rm(j+k+nFFTSize/2)).*rm(k+j+3*nFFTSize/4);
         Rm1(j)=abs(rm(k+j+nFFTSize/4)).^2;
         Rm2(j)=abs(rm(k+j+3*nFFTSize/4)).^2;
         end 
     Pk(k+nFFTSize+1)=sum([Pm1 Pm2]); 
     Rk(k+nFFTSize+1)=sum([Rm1 Rm2]); 
    end 
%Computing the metric 
M=((abs(Pk).^2))./(Rk.^2);
end

