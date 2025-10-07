function   [squared_gospa,gospa_loc,gospa_mis,gospa_fal]=ComputeGOSPAerror...
    (X_estimate,X_truth,t_birth,t_death,c_gospa,nIter,windowSize)
%Author: Angel F. Garcia-Fernandez | Modifications: Marco Fontana
%Computes the squared GOSPA error (alpha=2 and its decomposition, which is only possible for this choice of alpha)

%We need to convert the ground truth and estimated set into a different
%format to apply GOSPA code.

if nIter==1
    k=1;
else
    k=(nIter-1)*windowSize+1;
end
alive_targets_index=and(k>=t_birth,k<t_death);
Nalive_targets=sum(alive_targets_index);
%we select the locations of the alive targets
alive_targets_pos_index=[alive_targets_index;false(1,length(t_birth));alive_targets_index;false(1,length(t_birth))];
X_pos_truth_k=reshape(X_truth(alive_targets_pos_index(:),nIter),2,Nalive_targets);

%we consider just the locations
X_estimate_pos=X_estimate(1:2:end);
N_est=length(X_estimate_pos)/2; %two dimensions x and y
X_estimate_pos_reshape=reshape(X_estimate_pos,2,N_est);

[d_gospa, ~, decomp_cost] = GOSPA(X_pos_truth_k, X_estimate_pos_reshape, 2, c_gospa, 2);

squared_gospa=d_gospa^2;
gospa_loc=decomp_cost.localisation;
gospa_mis=decomp_cost.missed;
gospa_fal=decomp_cost.false;
