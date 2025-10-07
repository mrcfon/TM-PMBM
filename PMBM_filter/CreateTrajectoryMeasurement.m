function [z,zType,X_timeWindow]=CreateTrajectoryMeasurement(X_multi_k,...
    t_birth,t_death,windowSize,p_d,l_clutter,sensorP,k,H,chol_R,Nx,Nz)
%Author: Angel F. Garcia-Fernandez | Modifications: Marco Fontana
%W=whatever, S=start, E=end, F=full, P=partial

kNext=k+windowSize;
NzFull=2*Nz; %2-state measurements

%existence
existence_targetsS=and(t_birth<=k,t_death>k); %does the target exist at t1?
existence_targetsE=and(t_birth<=kNext,t_death>kNext); %does the target exist at t2?

existence_targetsF=(existence_targetsS&existence_targetsE); %is it a full state?
existence_targetsP1=(existence_targetsS&~existence_targetsE); %is it a type 1 partial state?
existence_targetsP2=(~existence_targetsS&existence_targetsE); %is it a type 2 partial state?

existence_targetsW=(existence_targetsS|existence_targetsE); %does it exist in the time window, in any way?
zType=zeros(1,length(existence_targetsW));

%detection
detected_targetsW=rand(1,length(existence_targetsW))<p_d.full; %for all kinds of measurements, it corresponds to p^D
detected_targetsP1=(detected_targetsW&existence_targetsP1); %detected measurements of type 1 partial state
detected_targetsP2=(detected_targetsW&existence_targetsP2); %detected measurements of type 2 partial state
detected_targetsF=(detected_targetsW&existence_targetsF); %detected full measurements
detected_targetsW=(detected_targetsW&existence_targetsW); %detected measurements (generic)

%3 cases for full measurements
fullTrackletDist=rand(1,length(existence_targetsW));
probPartialMisDet=(1-p_d.partial)/2;
misdetected_targetsS=fullTrackletDist<probPartialMisDet;
misdetected_targetsE=fullTrackletDist>=p_d.partial+probPartialMisDet;

misdetected_targetsS=detected_targetsF&misdetected_targetsS;
misdetected_targetsE=detected_targetsF&misdetected_targetsE;

%adding partially misdetected measurements to the sets
detected_targetsP1=detected_targetsP1|misdetected_targetsE;
detected_targetsP2=detected_targetsP2|misdetected_targetsS;
detected_targetsF(misdetected_targetsS|misdetected_targetsE)=0;

zType(detected_targetsP1)=1; %detected just at the beginning
zType(detected_targetsP2)=2; %detected just at the end
zType(detected_targetsF)=3;

detected_targetsW_idx=find(detected_targetsW);
X_multi_k_r=reshape(X_multi_k(:,1),Nx,[]);
X_multi_kNext_r=reshape(X_multi_k(:,end),Nx,[]);
X_timeWindow=[X_multi_k_r;X_multi_kNext_r];
z_timeWindow=zeros(NzFull,size(X_timeWindow,2));
for iX=detected_targetsW_idx
    %get the start point to compute the initial speed
    if existence_targetsP2(iX)==1
        iX_idx=(iX-1)*Nx+1;
        timeIni=t_birth(iX)-k+1;
        X_multi_k_r(:,iX)=X_multi_k(iX_idx:iX_idx+Nx-1,timeIni);
    end
    %get the end point to compute the final speed
    if existence_targetsP1(iX)==1
        iX_idx=(iX-1)*Nx+1;
        timeEnd=t_death(iX)-k;
        X_multi_kNext_r(:,iX)=X_multi_k(iX_idx:iX_idx+Nx-1,timeEnd);
    end
    z_timeWindow(:,iX)=H.full*[X_multi_k_r(:,iX);X_multi_kNext_r(:,iX)]+chol_R.full*randn(size(H.full,1),1);
end

z_timeWindow(:,~detected_targetsW)=[];
zType(~detected_targetsW)=[];
%apply partial misdetections
z_timeWindow(1:Nz,zType==2)=0;
z_timeWindow(Nz+1:end,zType==1)=0;

%clutter
N_clutter_measurementsS=poissrnd(l_clutter.intensity*l_clutter.dist(1));
N_clutter_measurementsE=poissrnd(l_clutter.intensity*l_clutter.dist(2));
N_clutter_measurementsF=poissrnd(l_clutter.intensity*l_clutter.dist(3));
N_clutter=N_clutter_measurementsS+N_clutter_measurementsE+N_clutter_measurementsF;
z_clutter_timeWindow=zeros(NzFull,N_clutter);
zType=[zType,2*ones(1,N_clutter_measurementsE),1*ones(1,N_clutter_measurementsS),...
    3*ones(1,N_clutter_measurementsF)];

z_clutter=[sensorP.Area(1)*rand(1,N_clutter);rand(1,N_clutter);...
    sensorP.Area(2)*rand(1,N_clutter);rand(1,N_clutter)];
z_clutter_timeWindow(:,1:N_clutter_measurementsE)=...
    [zeros(Nz,N_clutter_measurementsE);z_clutter(1:2:end,1:N_clutter_measurementsE)];
z_clutter_timeWindow(:,N_clutter_measurementsE+1:end)=...
    [z_clutter(1:2:end,N_clutter_measurementsE+1:end);zeros(Nz,N_clutter_measurementsS+N_clutter_measurementsF)];

F=kron(eye(2),[1 sensorP.T;0 1]);
Q=sensorP.q*kron(eye(2),[sensorP.T^3/3 sensorP.T^2/2; sensorP.T^2/2 sensorP.T]);
for iClutter=1:N_clutter_measurementsF
    idx_iC=N_clutter_measurementsS+N_clutter_measurementsE+iClutter;
    z_iC=z_clutter(1:Nx,idx_iC); %clutter target state at time k 
    z_iC_kNext=F*z_iC+chol(Q)'*randn(Nx,1);
    z_clutter_timeWindow(Nz+1:end,N_clutter_measurementsS+N_clutter_measurementsE+iClutter)=...
        z_iC_kNext(1:2:end);
end
z=[z_timeWindow,z_clutter_timeWindow];
