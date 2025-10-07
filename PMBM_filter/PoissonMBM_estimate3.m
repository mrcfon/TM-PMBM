function [X_estimate,metadata]=PoissonMBM_estimate3(filter_upd,cum_n_meas)

%Author: Angel F. Garcia-Fernandez | Modifications: Marco Fontana

%We select the hypothesis (with fixed cardinality, i.e. considering MBM01 expansion) with highest
%weight 

if(~isempty(filter_upd.globHyp))
    
    globHyp=filter_upd.globHyp;
    globHypWeight=filter_upd.globHypWeight;
    Nhyp=length(globHypWeight);
    N_tracks=length(filter_upd.tracks);

    
    weight_tot=globHypWeight;
    for j=1:Nhyp
        for i=1:N_tracks
            hyp=globHyp(j,i);
            if(hyp~=0)
                eB=filter_upd.tracks{i}.eB(hyp);
                if(eB>0.5)
                    weight_tot(j)=weight_tot(j)*eB;
                else
                    weight_tot(j)=weight_tot(j)*(1-eB);
                end
            end           
        end
    end
    
    [m,index_hyp]=max(weight_tot);    
    %We construct the estimate
    X_estimate=[];
    metadata=[];
    for i=1:N_tracks
        hyp=globHyp(index_hyp,i);
        if(hyp~=0)
            eB=filter_upd.tracks{i}.eB(hyp);
            if(eB>0.5)
                X_estimate=[X_estimate;filter_upd.tracks{i}.meanB{hyp}];
                meta_tmp=filter_upd.tracks{i}.t_ini(:,2:3);%[init_window,meas_idx]
                track_idx=meta_tmp(2)+cum_n_meas(meta_tmp(1));
                metadata=[metadata;track_idx,i,hyp];
            end
        end     
    end
   
else
    X_estimate=[];
end