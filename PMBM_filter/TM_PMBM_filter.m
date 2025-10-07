% Demo with the implementation of the 
% Trajectory Measurement Poisson multi-Bernoulli mixture (TM-PMBM) filter 
% described in
% M. Fontana, Á. F. García-Fenández, S. Maskell, 
% "Poisson multi-Bernoulli mixture filter for trajectory measurements". 
% Available: https://arxiv.org/abs/2504.08421

% Based on the original implementation of the PMBM filter descibed in 
% Á. F. García-Fenández, J. L. Williams, K. Granstrom and L. Svensson, 
% "Poisson Multi-Bernoulli Mixture Filter: Direct Derivation and Implementation," 
% in IEEE Transactions on Aerospace and Electronic Systems, 
% vol. 54, no. 4, pp. 1883-1901, Aug. 2018.

% Performance is measured with the GOSPA metric (alpha=2) and its
% decomposition into localisation errors for properly detected targets,
% and costs for missed targets and false targets (only possible for alpha=2).

% A. S. Rahmathullah, Á. F. García-Fenández and L. Svensson, 
% "Generalized optimal sub-pattern assignment metric," 
% 2017 20th International Conference on Information Fusion (Fusion),
% Xi'an, 2017, pp. 1-8.
% Short video on GOSPA: https://www.youtube.com/watch?v=M79GTTytvCM

% Copyright (c) 2023, Marco Fontana, Angel F. Garcia-Fernandez and Simon
% Maskell
% All rights reserved.

clear
addpath('..\GOSPA code')
addpath('..\Assignment')

rand('seed',9)
randn('seed',9)

%Filter parameters
T_pruning=0.0001; %Threshold for pruning multi-Bernoulli mixtures weights
T_pruningPois=10^(-5); %Threshold for pruning PHD of the Poisson component
Nhyp_max=200;  %Maximum number of hypotheses (MBM components)
gating_threshold=9; %Threshold for gating
existence_threshold=1e-5; %Existence threshold: Bernoulli components with existence below this threshold are removed
distance_threshold=0.25; %Distance threshold for merge, pairs of Bernoullis below this threshold are merged

type_estimator=1; %Choose Estimator 1, 2 or 3 as defined in the paper
existence_estimation_threshold1=0.1; %Only for esimator 1

windowSize=5;
Nmc=30; % Number of MC runs

t_th=zeros(1,1);
rms_gospa_tot_th=zeros(2,1);
rms_gospa_loc_th=zeros(2,1);
rms_gospa_false_th=zeros(2,1);
rms_gospa_mis_th=zeros(2,1);

ScenarioWilliams15tracklets_GM;
%ScenarioReal_tracklets_GM;

NtimeWindows=floor((Nsteps-1)/windowSize)+1; %k=0 is not is computed separately
time_timestep=zeros(Nmc,NtimeWindows);

%GOSPA errors for the estimations at the end of each time window
squared_gospa_t_tot=zeros(1,NtimeWindows);
squared_gospa_loc_t_tot=zeros(1,NtimeWindows); %Localisation error
squared_gospa_false_t_tot=zeros(1,NtimeWindows); %False target error
squared_gospa_mis_t_tot=zeros(1,NtimeWindows); %Misdetection error

%We go through all Monte Carlo runs
for i=1:Nmc
    tic
    %initlialisation on the target state space
    filter_predTrj.weightPois=weightsTrj_b0;
    filter_predTrj.meanPois=meansTrj_b0;
    filter_predTrj.covPois=covsTrj_b0;
    filter_predTrj.typeP=types_b0;
    filter_predTrj.tracks=cell(0,1);
    filter_predTrj.globHyp=[];
    filter_predTrj.globHypWeight=[];

    z_trj_t=cell(2,NtimeWindows); %measurements in trajectory space
    tracklets_t=cell(1,NtimeWindows);
    n_measTrj_kW=zeros(1,NtimeWindows);
    
    %time window ending at for k=1 (window0, indexed at kW=1)
    [z,zType,X_timeWindow]=CreateTrajectoryMeasurement(...
        X_truth(:,1:windowSize+1),...
        t_birth,t_death,windowSize,p_d,l_clutter,sensorP,1,H,chol_R,Nx,Nz);
    z0_meas_idx=find(zType~=2);
    n_measTrj_kW(1)=length(z0_meas_idx);
    %tracklet measurements
    z_trj_t{1,1}=[zeros(Nz,n_measTrj_kW(1));z(1:Nz,z0_meas_idx)];
    z_trj_t{2,1}=2*ones(1,n_measTrj_kW(1));
    tracklets_t{1}=[zeros(Nx,size(X_timeWindow,2));X_timeWindow(1:Nx,:)];
    X_truth_sampled=X_truth(:,1);

    %Simulate measurements
    for k=1:windowSize:Nsteps-windowSize
        [z,zType,X_timeWindow]=CreateTrajectoryMeasurement(...
            X_truth(:,k:k+windowSize),...
            t_birth,t_death,windowSize,p_d,l_clutter,sensorP,k,H,chol_R,Nx,Nz);   
        nIter=ceil(k/windowSize)+1;
        %tracklet measurements
        z_trj_t{1,nIter}=z;
        z_trj_t{2,nIter}=zType;
        tracklets_t{nIter}=X_timeWindow;
        %we save just the ground-truth at the end of the time window
        X_truth_sampled=[X_truth_sampled,X_truth(:,k+windowSize)];
        n_measTrj_kW(nIter)=size(z,2);
    end
    cum_n_Trj=[0,cumsum(n_measTrj_kW)];

    %Perform filtering
    for nIter=1:NtimeWindows
        if nIter==1 %the first iteration consider the meas ending at k=1
            iniK=0;
            finK=1;
        else
            iniK=(nIter-2)*windowSize+1;
            finK=(nIter-1)*windowSize+1;
        end

        %Update in the trajectory space
        z=z_trj_t{1,nIter};
        zType_t=z_trj_t{2,nIter};
        filter_updTrj=PoissonMBM_update(filter_predTrj,z,zType_t,windowSize,...
            H,R,p_d,iniK,finK,gating_threshold,intensity_clutter,Nhyp_max);
        %Marginalisation for tracklet states
        filter_updTrj=PoissonMBM_marginalisation(filter_updTrj,Nx);
       
        %State estimation         
        switch type_estimator
            case 1
                [X_estimateTrj,metadataTrj]=PoissonMBM_estimate1(filter_updTrj,existence_estimation_threshold1,cum_n_Trj);
            case 2
                [X_estimateTrj,metadataTrj]=PoissonMBM_estimate2(filter_updTrj,cum_n_Trj);
            case 3
                [X_estimateTrj,metadataTrj]=PoissonMBM_estimate3(filter_updTrj,cum_n_Trj);
        end

        %Computation of squared GOSPA error and its decomposition
        %Obtain ground truth state
        [squared_gospaTrj,gospa_locTrj,gospa_misTrj,gospa_falTrj]=ComputeGOSPAerror(X_estimateTrj,X_truth_sampled,t_birth,t_death,c_gospa,nIter,windowSize);

        %We sum the squared errors
        squared_gospa_t_tot(nIter)=squared_gospa_t_tot(nIter) + squared_gospaTrj;
        squared_gospa_loc_t_tot(nIter)=squared_gospa_loc_t_tot(nIter) + gospa_locTrj;
        squared_gospa_false_t_tot(nIter)=squared_gospa_false_t_tot(nIter) + gospa_falTrj;
        squared_gospa_mis_t_tot(nIter)=squared_gospa_mis_t_tot(nIter) + gospa_misTrj;
        
        if Nmc==1
            visMode=0;
            if exist('scenarioF','var') == 0
                scenarioF=figure;
            end
            DrawEstimates(X_truth,t_birth,t_death,X_estimateTrj,...
                tracklets_t,XLim,YLim,z_trj_t,Nx,Nz,nIter,metadataTrj,0);
        end

        %We project the PMBM updated density into a PMB density
        %filter_updTrj=PMB_projection(filter_updTrj);

        %Hypothesis reduction, pruning, normalisation
        filter_updTrj=PoissonMBM_pruning(filter_updTrj,T_pruning,T_pruningPois,Nhyp_max,existence_threshold);
        
        %Prediction on trajectories space
        filter_predTrj=PoissonMBM_pred(filter_updTrj,F,Q,p_s,weightsTrj_b,meansTrj_b,covsTrj_b,types_b);
    end
    t=toc;
    disp(['Completed iteration number ', num2str(i),' time ', num2str(t), ' sec'])
end
%Root mean square GOSPA errors at each time step
rms_gospa_t=sqrt(squared_gospa_t_tot/Nmc);
rms_gospa_loc_t=sqrt(squared_gospa_loc_t_tot/Nmc);
rms_gospa_false_t=sqrt(squared_gospa_false_t_tot/Nmc);
rms_gospa_mis_t=sqrt(squared_gospa_mis_t_tot/Nmc);

%Root mean square GOSPA errors across all time steps
rms_gospa_tot=sqrt(sum(squared_gospa_t_tot)/(Nmc*NtimeWindows))
rms_gospa_loc_tot=sqrt(sum(squared_gospa_loc_t_tot)/(Nmc*NtimeWindows))
rms_gospa_false_tot=sqrt(sum(squared_gospa_false_t_tot)/(Nmc*NtimeWindows))
rms_gospa_mis_tot=sqrt(sum(squared_gospa_mis_t_tot)/(Nmc*NtimeWindows))

gospaTotF=figure;
subplot(221)
plot(1:NtimeWindows,rms_gospa_t,'Linewidth',1.3)
grid on
xlabel('Time step')
ylabel('RMS GOSPA error')

subplot(222)
plot(1:NtimeWindows,rms_gospa_loc_t,'Linewidth',1.3)
grid on
xlabel('Time step')
ylabel('RMS GOSPA localisation error')

subplot(223)
plot(1:NtimeWindows,rms_gospa_false_t,'Linewidth',1.3)
grid on
xlabel('Time step')
ylabel('RMS GOSPA false target error')

subplot(224)
plot(1:NtimeWindows,rms_gospa_mis_t,'Linewidth',1.3)
grid on
xlabel('Time step')
ylabel('RMS GOSPA missed target error')