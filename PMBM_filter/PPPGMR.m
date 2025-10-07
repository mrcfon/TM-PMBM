function [w_b,meanB,covB]=PPPGMR(windowSize,p_s,Nx,F,Q,weights_b,mean_b,cov_b)
% Poisson Point Process (PPP) Gaussian Mixture (GM) reduction for birth model
% Author: Marco Fontana

weightsList=weights_b;
meansList=mean_b;
covsList(:,:,1)=cov_b;
%we go through prediction for k+1,...,k+windowSize
for iUpds=1:windowSize-1
    %prediction
    weightsList=p_s*weightsList;
    meansList=F*meansList;
    for iC=1:length(weightsList)
        covsList(:,:,iC)=F*covsList(:,:,iC)*F'+Q;
    end
    %add a new component
    weightsList=[weightsList,weights_b];
    meansList=[meansList,mean_b];
    covsList=cat(3,covsList,cov_b);
end
%GM reduction via moment matching
weight_fused=0;
mean_fused=zeros(size(mean_b));
cov_fused=zeros(size(cov_b));
for iP=1:length(weightsList)
    weight_fused=weight_fused+weightsList(iP);
    mean_fused=mean_fused+weightsList(iP)*meansList(:,iP);
    cov_fused=cov_fused+weightsList(iP)*(covsList(:,:,iP)+(meansList(:,iP)*meansList(:,iP)'));
end
mean_fused=mean_fused/weight_fused;
cov_fused=cov_fused/weight_fused-(mean_fused*mean_fused');

%we add the PPP component at the end of the time window, of type 2
w_b=weight_fused;
meanB=mean_fused;
covB=cov_fused;
end