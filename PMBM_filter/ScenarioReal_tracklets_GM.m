%% SIMULATION PARAMETERS
Nx=4; %Single target state dimension
Nz=2; %Single measurement dimension
Nsteps=251; %Considered number of time steps in the simulation
sensorP.Area=[600 400];
c_gospa=10; %Parameter c of the GOSPA metric. We also consider p=2 and alpha=2
%for visualisation
XLim=[0,sensorP.Area(1)];
YLim=[0,sensorP.Area(2)];

%% MODEL PARAMETERS
%Measurement parameters
sensorP.q=0.01;
sensorP.T=0.2;

%ground-truth model
F=kron(eye(2),[1 sensorP.T;0 1]);
Q=sensorP.q*kron(eye(2),[sensorP.T^3/3 sensorP.T^2/2; sensorP.T^2/2 sensorP.T]);
mu=0.001;
lambda=0.8;
p_s=exp(-mu*sensorP.T); %probability of survival in the ground-truth

H.partial=kron(eye(Nz),[1,0]);
H.full=kron(eye(2),H.partial);
H.p2=kron([0,1],H.partial);
H.p1=kron([1,0],H.partial);

R.partial=0.1*eye(Nz);
R.full=kron(eye(2),R.partial);
chol_R.partial=chol(R.partial)';
chol_R.full=chol(R.full)';

%Probability ofPccc detection
p_d.full=0.9; %prob. of detection of a tracklet measurement
p_d.partial=0.7; %prob. of full detection given a full measurement
alpha=0.5*(1-p_d.partial);
p_d.partialNC=p_d.full*(1-alpha); %eq. prob. of detection for target implementation
p_d.decomp = 1-(p_d.full*alpha+(1-p_d.full))^2;

%clutter parameters
intensity_clutter=1e-6;
l_clutter.dist=ones(1,3)/3; %clutter distribution among the tracklet types
l_clutter.intensity=intensity_clutter*sensorP.Area(1)*sensorP.Area(2);

%% PPP components for simulation and prior
mean_b=[sensorP.Area(1)/2;0;sensorP.Area(2)/2;0];
cov_b(1:4,1:4,1)=diag([sensorP.Area(1)/2 1 sensorP.Area(2)/2 1].^2);

%% BIRTH PARAMETERS for simulation
%We generate the time intervals
delta_tk_t=ones(Nsteps,1)*sensorP.T;
tk=cumsum(delta_tk_t); %Time at which is measurement is obtained

%Poisson birth
lambda_new_born_t=lambda/mu*(1-exp(-mu*delta_tk_t));
N_new_born_t=poissrnd(lambda_new_born_t);

%% SIMULATION
N_targets_tot=sum(N_new_born_t);
%Life span for each target is exponential distrubuted
life_span_continuous=exprnd(1/mu,N_targets_tot,1);

%X_truth matrix with the ground truth set of trajectories/targets
X_truth=zeros(Nx*N_targets_tot,Nsteps);
t_birth=zeros(1,N_targets_tot); %t_birth is the discrete time step at which the target is born
t_death=zeros(1,N_targets_tot); %t_death is the discrete time step the target is already dead

N_new_born_t_prov=N_new_born_t; 
time_aux=1; %Auxiliary variable
N_targets_alive_t=zeros(1,Nsteps);

P_a=diag([30 1 30 1].^2); %to spread the location fo appearence of the targets
for i=1:N_targets_tot
    while(N_new_born_t_prov(time_aux)==0)
        time_aux=time_aux+1;
    end
    N_new_born_t_prov(time_aux)=N_new_born_t_prov(time_aux)-1;
    t_birth(i)=time_aux;
    t_death_continuous=life_span_continuous(i)+tk(time_aux);
    t_death_discrete=find(tk<t_death_continuous);
    t_death_discrete=t_death_discrete(end)+1; %We sum one
    t_death(i)=t_death_discrete;
    
    N_targets_alive_t(t_birth(i):t_death(i)-1)=N_targets_alive_t(t_birth(i):t_death(i)-1)+1;

    %Sample of target state at appearing time
    X_a=mean_b+chol(P_a)'*randn(Nx,1);
    X_truth(4*i-3:4*i,time_aux)=X_a;
    
    %Sample of target states at the following time steps
    for k=time_aux+1:t_death(i)-1
        X_truth(4*i-3:4*i,k)=F*X_truth(4*i-3:4*i,k-1)+chol(Q)'*randn(Nx,1);
    end
end
%% PPP for birth process
%PPP component(s) for one time step
nPHD_comp=1; %n. components at each time step
mean_b=[sensorP.Area(1)/2;0;sensorP.Area(2)/2;0];
cov_b(1:4,1:4,1)=diag([sensorP.Area(1)/2 1 sensorP.Area(1)/2 1].^2);
weights_b0=lambda_new_born_t(1); %weight of the PPP component at k=1
weights_b=lambda_new_born_t(1); %weight of the PPP component at k>1

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
sensorP.T=sensorP.T*windowSize;
F=kron(eye(2),[1 sensorP.T;0 1]); 
Q=sensorP.q*kron(eye(2),[sensorP.T^3/3 sensorP.T^2/2; sensorP.T^2/2 sensorP.T]);
p_s=p_s^(windowSize); %probability of survival of in the trajectory space

%% VISUALISATION
%uncomment to plot the scenario
% scenarioVis=figure;
% set(gcf, 'DefaultTextInterpreter', 'latex');
% set(gcf, 'DefaultLegendInterpreter', 'latex');
% col_set=jet(Nsteps);
% for k=1:Nsteps
%     for i=1:N_targets_tot
%         if (k>=t_birth(i))&&(k<t_death(i))
%             plot(X_truth(4*i-3,k),X_truth(4*i-1,k),'Color',col_set(k,:),'Marker','.','MarkerSize',12)
%             hold on
%             if t_birth(i)==k
%                 text(X_truth(4*i-3,t_birth(i))-20,X_truth(4*i-1,t_birth(i))+20,num2str(i))
%             end
%         end
%     end
%     title(['time=',num2str(k),', n tar alive=',num2str(N_targets_alive_t(k)),' n tar=',num2str(N_targets_tot)])
%     xlim([0,sensorP.Area(1)])
%     ylim([0,sensorP.Area(2)])
% 
%     xlabel('x position (m)')
%     ylabel('y position (m)')
%     axis equal
%     axis([200 400 100 350])
%     grid on
% 
%     pause(0.1);
% end

