%% SIMULATION PARAMETERS
Scenario_number=1;
numtruth=4; % number of targets

Nx=4; %Single target state dimension
Nz=2; %Single measurement dimension
Nsteps=251; %Considered number of time steps in the simulation
sensorP.Area=[100 100];
c_gospa=10; %Parameter c of the GOSPA metric. We also consider p=2 and alpha=2
%for visualisation
XLim=[0,sensorP.Area(1)];
YLim=[0,sensorP.Area(2)];

%% MODEL PARAMETERS
sensorP.q=0.01;
sensorP.T = 0.2;

%ground-truth model
F=kron(eye(2),[1 sensorP.T;0 1]);
Q=sensorP.q*kron(eye(2),[sensorP.T^3/3 sensorP.T^2/2; sensorP.T^2/2 sensorP.T]);
p_s=0.99;

H.partial=kron(eye(Nz),[1,0]);
H.full=kron(eye(2),H.partial);
H.p2=kron([0,1],H.partial);
H.p1=kron([1,0],H.partial);

R.partial=0.01*eye(Nz);
R.full=kron(eye(2),R.partial);
chol_R.partial=chol(R.partial)';
chol_R.full=chol(R.full)';

%Probability of detection
p_d.full=0.9; %prob. of detection of a tracklet measurement
p_d.partial=0.7; %prob. of full detection given a full measurement
p_d.alpha=0.5*(1-p_d.partial);
p_d.partialNC=p_d.full*(1-p_d.alpha); %eq. prob. of detection for target implementation
p_d.decomp = 1-(p_d.full*p_d.alpha+(1-p_d.full))^2;

%clutter parameters
intensity_clutter=1e-6;
l_clutter.dist=ones(1,3)/3; %clutter distribution among the tracklet types
l_clutter.intensity=intensity_clutter*sensorP.Area(1)*sensorP.Area(2);

%% SIMULATION
[X_truth,t_birth,t_death]=TrajectoryWilliams15(Scenario_number,Nsteps,F,numtruth,Q,sensorP.Area);
N_targets_alive_t=zeros(1,Nsteps);
for ik=1:Nsteps
    N_targets_alive_t(ik)=nnz(ik>=t_birth)-nnz(ik>=t_death);
end

%% PPP for birth process
%PPP component(s) for one time step
nPHD_comp=1; %n. components at each time step
mean_b=[sensorP.Area(1)/2;0;sensorP.Area(2)/2;0];
cov_b(1:4,1:4,1)=diag([sensorP.Area(1)/2 1 sensorP.Area(2)/2 1].^2);
weights_b0=3; %weight of the PPP component at k=1
weights_b=0.005; %weight of the PPP component at k>1

%first time step (k=0)
weightsTrj_b0=weights_b0;
meansTrj_b0=kron([0;1],mean_b);
covsTrj_b0=kron([0,0;0,1],cov_b);
types_b0=2*ones(1,nPHD_comp);

weightsTar_b0=weights_b0;
meansTar_b0=mean_b;
covsTar_b0=cov_b;

%computing the PPP components on the time window via moment matching
[weightsTrj_b_merged,mean_b_merged,cov_b_merged]=PPPGMR(windowSize,p_s,Nx,...
    F,Q,weights_b,mean_b,cov_b);
weightsTrj_b=sum(weightsTrj_b_merged);
meansTrj_b=kron([0;1],mean_b_merged);
covsTrj_b=kron([0,0;0,1],cov_b_merged);
types_b=2*ones(1,nPHD_comp);
%adaptng for target state implemenation
weightsTar_b=sum(weightsTrj_b);
meansTar_b=mean_b_merged;
covsTar_b=cov_b_merged;

%% TRACKING MODEL
sensorP.T=sensorP.T*windowSize; %we consider a step=window in the simulation
F=kron(eye(2),[1 sensorP.T;0 1]);
Q=sensorP.q*kron(eye(2),[sensorP.T^3/3 sensorP.T^2/2; sensorP.T^2/2 sensorP.T]);
p_s=p_s^(windowSize); %probability of survival of in the trajectory space

%% Visualisation
% figure(5)
% clf
% plot(X_truth(1,t_birth(1):t_death(1)-1),X_truth(3,t_birth(1):t_death(1)-1),'b','Linewidth',1.3)
% hold on
% plot(X_truth(1,t_birth(1)),X_truth(3,t_birth(1)),'xb','Linewidth',1.3)
% plot(X_truth(1,t_birth(1):5:(t_death(1)-1)),X_truth(3,t_birth(1):5:(t_death(1)-1)),'ob','Linewidth',1.3)
% 
% plot(X_truth(5,t_birth(2):t_death(2)-1),X_truth(7,t_birth(2):t_death(2)-1),'r','Linewidth',1.3)
% plot(X_truth(5,t_birth(2)),X_truth(7,t_birth(2)),'xr','Linewidth',1.3)
% plot(X_truth(5,t_birth(2):5:t_death(2)-1),X_truth(7,t_birth(2):5:t_death(2)),'or','Linewidth',1.3)
% 
% 
% plot(X_truth(9,t_birth(3):t_death(3)-1),X_truth(11,t_birth(3):t_death(3)-1),'g','Linewidth',1.3)
% plot(X_truth(9,t_birth(3)),X_truth(11,t_birth(3)),'xg','Linewidth',1.3)
% plot(X_truth(9,t_birth(3):5:t_death(3)-1),X_truth(11,t_birth(3):5:t_death(3)-1),'og','Linewidth',1.3)
% 
% plot(X_truth(13,t_birth(4):t_death(4)-1),X_truth(15,t_birth(4):t_death(4)-1),'black','Linewidth',1.3)
% plot(X_truth(13,t_birth(4)),X_truth(15,t_birth(4)),'xblack','Linewidth',1.3)
% plot(X_truth(13,t_birth(4):5:t_death(4)-1),X_truth(15,t_birth(4):5:t_death(4)-1),'oblack','Linewidth',1.3)
% 
% %axis([0 Area(1) 0 Area(2)])
% hold off
% xlabel('x position (m)')
% ylabel('y position (m)')
% axis equal
% %axis([110 180 110 180])
% grid on