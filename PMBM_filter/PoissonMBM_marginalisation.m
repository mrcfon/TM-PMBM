function filter_upd_marg=PoissonMBM_marginalisation(filter_updTrj,Nx)
%Author: Marco Fontana
filter_upd_marg=filter_updTrj;

%Marginalisation Poisson
filter_upd_marg.meanPois=filter_upd_marg.meanPois(Nx+1:end,:);
filter_upd_marg.covPois=filter_upd_marg.covPois(Nx+1:end,Nx+1:end,:);
filter_upd_marg.typeP = [];

Ntracks=length(filter_upd_marg.tracks);
for i=1:Ntracks
    idx = filter_upd_marg.tracks{i}.eB ~= 0;
    beta_i = filter_upd_marg.tracks{i}.beta;
    eB_i = filter_upd_marg.tracks{i}.eB;
    meanB_i = filter_upd_marg.tracks{i}.meanB;
    covB_i = filter_upd_marg.tracks{i}.covB;
    sel_beta = beta_i(idx);
    sel_eB = eB_i(idx);
    sel_mean = meanB_i(idx);
    sel_cov = covB_i(idx);
    for k = 1:numel(sel_mean)
        sel_eB(k) = sel_eB(k) * sel_beta{k}(end);
        sel_mean{k} = sel_mean{k}(end-Nx+1:end);
        sel_cov{k} = sel_cov{k}(end-Nx+1:end,end-Nx+1:end);
    end
    filter_upd_marg.tracks{i}.eB(idx) = sel_eB;
    filter_upd_marg.tracks{i}.meanB(idx) = sel_mean;
    filter_upd_marg.tracks{i}.covB(idx) = sel_cov;
    filter_upd_marg.tracks{i}.beta=[];
end