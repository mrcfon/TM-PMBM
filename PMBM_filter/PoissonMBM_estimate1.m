function [X_estimate,metadata]=PoissonMBM_estimate1(filter_upd,existence_estimation_threshold...
    ,cum_n_meas)

%Author: Angel F. Garcia Fernandez and Marco Fontana


%Option 1: one picks global hypothesis with highest weight and then takes the
%estimates above a certain threshold  
X_estimate=[];
metadata=[];
if(~isempty(filter_upd.globHyp))
    globHypWeight=filter_upd.globHypWeight;
    [~,index]=max(globHypWeight);

    HypMax=filter_upd.globHyp(index,:);

    for i=1:length(HypMax)
        hyp_i=HypMax(i);
        if(hyp_i>0)
            Existence=filter_upd.tracks{i}.eB(hyp_i);
            if(Existence>existence_estimation_threshold)
                
                %we estimate the entire state, position and velocity, at
                %the end of the time window
                X_estimate_i=filter_upd.tracks{i}.meanB{hyp_i};
                X_estimate=[X_estimate;X_estimate_i];
                meta_tmp=filter_upd.tracks{i}.t_ini(:,2:3);%[init_window,meas_idx]
                track_idx=meta_tmp(2)+cum_n_meas(meta_tmp(1));
                metadata=[metadata;track_idx,i,hyp_i];
            end
        end
    end
end
