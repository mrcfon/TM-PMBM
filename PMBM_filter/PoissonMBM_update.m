function filter_updTrj=PoissonMBM_update(filter_pred,z,zType,windowSize,...
    H,R,p_d,iniK,finK,gating_threshold,intensity_clutter,Nhyp_max)
%Author: Marco Fontana and Angel F. Garcia-Fernandez
%mexMurty by David F. Crouse

[Nz,Nx]=size(H.partial);
NxF=Nx*2;
NzF=Nz*2;
selV=cell(1,3);
selV{1}=(1:Nz);
selV{2}=(Nz+1:NzF);
selV{3}=(1:NzF);
kW=ceil(finK/windowSize);

Nprev_tracks=length(filter_pred.tracks);

%Poisson update
filter_updTrj.weightPois=(1-p_d.full)*filter_pred.weightPois;
filter_updTrj.meanPois=filter_pred.meanPois;
filter_updTrj.covPois=filter_pred.covPois;
filter_updTrj.typeP=filter_pred.typeP;

%New track generation
Nnew_tracks=size(z,2); %Number of new tracks (some may have existence probability equal to zero)

for m=1:size(z,2)
    %We go through the Poisson components
    switch zType(m)
        case 1
            comp2Check=find(filter_updTrj.typeP==3);
            startZ_m=iniK;
            startNewState=iniK;
            Hbar=H.p1;
            Rbar=R.partial;
        case 2
            comp2Check=(1:length(filter_pred.weightPois));
            startZ_m=finK;
            startNewState=finK;
            Hbar=[]; %depends on the trajectory state type
            Rbar=R.partial;
        case 3
            comp2Check=find(filter_updTrj.typeP==3);
            startZ_m=iniK;
            startNewState=iniK;
            Hbar=H.full;
            Rbar=R.full;
    end
    %select the detected part of the measurement
    z_m=z(selV{zType(m)},m);
    v_q=zeros(1,length(filter_updTrj.weightPois));

    for i=comp2Check
        [Hbar,selStates_i,alpha]=PPPparam(filter_pred.typeP(i),zType(m),Nx,NxF,H,Hbar,p_d);
        %select the state in which the windowed trajectory exists
        mean_i=filter_pred.meanPois(selStates_i,i);
        cov_i=filter_pred.covPois(selStates_i,selStates_i,i);

        if zType(m)==3
            z_maha=z_m(selV{2});
            Hmaha=H.p2;
            Rmaha=R.partial;
        else
            z_maha=z_m;
            Hmaha=Hbar;
            Rmaha=Rbar;
        end
        
        S_pred_i_maha=Hmaha*cov_i*Hmaha'+Rmaha;
        z_pred_i_maha=Hmaha*mean_i;
        maha=(z_maha-z_pred_i_maha)'/S_pred_i_maha*(z_maha-z_pred_i_maha);


        S_pred_i=Hbar*cov_i*Hbar'+Rbar;
        z_pred_i=Hbar*mean_i;
        if (maha<gating_threshold)
            v_q(i)=p_d.full*alpha*mvnpdf(z_m,z_pred_i,S_pred_i)*filter_pred.weightPois(i);
        end
        
    end

    eB=sum(v_q);
    if eB>0
        %A new track is created for this measurement
        [~,q_star]=max(v_q);
        
        weightB=eB+intensity_clutter;
        %compute mean and covariance taking the component q_star
        [Hbar,selStates_i,~]=PPPparam(filter_pred.typeP(q_star),zType(m),Nx,NxF,H,Hbar,p_d);
        mean_i=filter_pred.meanPois(selStates_i,q_star);
        cov_i=filter_pred.covPois(selStates_i,selStates_i,q_star);
        
        S_pred_i=Hbar*cov_i*Hbar'+Rbar;
        z_pred_i=Hbar*mean_i;

        K_pred_i=cov_i*Hbar'/S_pred_i;
        P_u_i=(eye(size(Hbar,2))-K_pred_i*Hbar)*cov_i;
        P_u_i=(P_u_i+P_u_i')/2;
        x_u_i=mean_i+K_pred_i*(z_m-z_pred_i);

        filter_updTrj.tracks{Nprev_tracks+m}.meanB{1}=x_u_i;
        filter_updTrj.tracks{Nprev_tracks+m}.covB{1}=P_u_i;
        filter_updTrj.tracks{Nprev_tracks+m}.eB=eB/weightB;
        filter_updTrj.tracks{Nprev_tracks+m}.t_ini=[startNewState,kW,m];%[init_time,init_window,meas_idx]
        filter_updTrj.tracks{Nprev_tracks+m}.aHis{1}=m;
        filter_updTrj.tracks{Nprev_tracks+m}.weightBLog_k=log(weightB);
        filter_updTrj.tracks{Nprev_tracks+m}.beta{1}=[0,1];
    else
        %The Bernoulli component has existence probability zero (it is
        %clutter). It will be removed by pruning
        weightB=intensity_clutter;
        
        filter_updTrj.tracks{Nprev_tracks+m}.eB=0;
        filter_updTrj.tracks{Nprev_tracks+m}.weightBLog_k=log(weightB);
        
        filter_updTrj.tracks{Nprev_tracks+m}.meanB{1}=zeros(length(selV{zType(m)}),1);
        filter_updTrj.tracks{Nprev_tracks+m}.covB{1}=zeros(length(selV{zType(m)}),length(selV{zType(m)}));
        filter_updTrj.tracks{Nprev_tracks+m}.t_ini=[startZ_m,m];
        filter_updTrj.tracks{Nprev_tracks+m}.aHis{1}=m;
        filter_updTrj.tracks{Nprev_tracks+m}.beta{1}=0;
    end 
end



%We update all previous tracks with the measurements
for i=1:Nprev_tracks
    Nhyp_i=length(filter_pred.tracks{i}.eB);
    filter_updTrj.tracks{i}.t_ini=filter_pred.tracks{i}.t_ini;
    
    for j=1:Nhyp_i %We go through all hypotheses
        
        mean_j=filter_pred.tracks{i}.meanB{j};
        cov_j=filter_pred.tracks{i}.covB{j};
        eB_j=filter_pred.tracks{i}.eB(j);
        aHis_j=filter_pred.tracks{i}.aHis{j};
        beta_j=filter_pred.tracks{i}.beta{j};
        
        %Misdetection hypotheses
        filter_updTrj.tracks{i}.meanB{j}=mean_j;
        filter_updTrj.tracks{i}.covB{j}=cov_j;
        filter_updTrj.tracks{i}.eB(j)=eB_j*(1-p_d.full)/(1-eB_j+eB_j*(1-p_d.full));
        filter_updTrj.tracks{i}.aHis{j}=[aHis_j,0];
        filter_updTrj.tracks{i}.beta{j}=beta_j;
        filter_updTrj.tracks{i}.weightBLog_k(j)=log(1-eB_j*p_d.full); %Weight only at time step k (beta=1)

        %KF moments for all the possible cases
        [z_pred_j,inv_S_pred_j,K_pred_j,S_pred_j,Hlist]=MBparam(mean_j,cov_j,H,R);  
        
        %We go through all measurements
        for m=1:size(z,2)
            index_hyp=j+Nhyp_i*m;
            zTm=zType(m);
            switch zTm
                case 1
                    %we compute only mean and cov of the full state
                    beta=[beta_j(1),beta_j(2)*(1-p_d.partial)/2];
                    alpha=sum(beta);
                    beta=beta/sum(beta);
                case 2
                    beta=[0,1];
                    alpha=beta_j(2)*(1-p_d.partial)/2;
                case 3
                    beta=[0,1];
                    alpha=beta_j(2)*p_d.partial;
            end

            %gating
            z_m=z(selV{zTm},m);
            mahaV=ones(1,2)*(gating_threshold+1);
            if zTm==3
                mahaV(1)=mahaDist(z_m(selV{1}),z_pred_j{1},inv_S_pred_j{1});
                mahaV(2)=mahaDist(z_m(selV{2}),z_pred_j{2},inv_S_pred_j{2});
            else
                %[i,m,zTm]
                mahaV(zTm)=mahaDist(z_m,z_pred_j{zTm},inv_S_pred_j{zTm});
            end
            
            if (mahaV(2)<gating_threshold)||(mahaV(1)<gating_threshold)
                %We create the component
                P_u=(eye(size(Hlist{zTm},2))-K_pred_j{zTm}*Hlist{zTm})*cov_j;
                P_u=(P_u+P_u')/2;
                x_u=mean_j+K_pred_j{zTm}*(z_m-z_pred_j{zTm});
                
                filter_updTrj.tracks{i}.meanB{index_hyp}=x_u;
                filter_updTrj.tracks{i}.covB{index_hyp}=P_u;
                filter_updTrj.tracks{i}.eB(index_hyp)=1;
                filter_updTrj.tracks{i}.aHis{index_hyp}=[aHis_j,m];
                filter_updTrj.tracks{i}.beta{index_hyp}=beta;
                
                %Update of the weights
                filter_updTrj.tracks{i}.weightBLog_k(index_hyp)=log(eB_j*p_d.full*alpha...
                    *mvnpdf(z_pred_j{zTm},z_m,S_pred_j{zTm})); %could be -Inf
            else
                filter_updTrj.tracks{i}.weightBLog_k(index_hyp)=-Inf;
            end
        end
    end
    
end

%Update of global hypotheses

if(Nprev_tracks==0)
    %No measurement hypothesised to be received from a target yet (No
    %previous tracks)
    %We create one global hypothesis
    filter_updTrj.globHypWeight=1;
    filter_updTrj.globHyp=ones(1,Nnew_tracks);
    
    if(Nnew_tracks==0)
        filter_updTrj.tracks=cell(0,1);
    end
        
else
    
    globWeightLog=[];
    globHyp=[];
    
    for p=1:length(filter_pred.globHypWeight)
        %We create cost matrix for each global hypothesis.
        cost_matrix_log=-Inf(Nprev_tracks+Nnew_tracks,size(z,2));
        cost_misdetection=zeros(Nprev_tracks,1);
        
        for i=1:Nprev_tracks
            index_hyp=filter_pred.globHyp(p,i); %Hypothesis for track i in p global hypothesis
            Nhyp_i=length(filter_pred.tracks{i}.eB);
            %We generate the cost matrix for measurements
            
            if(index_hyp~=0)
                              
                index_max=length(filter_updTrj.tracks{i}.weightBLog_k);
                indices=index_hyp+Nhyp_i*(1:size(z,2));
                indices_c=indices<=index_max;
                indices_c=indices(indices_c);
                
                weights_log=filter_updTrj.tracks{i}.weightBLog_k(indices_c);
                
                %We remove the weight of the misdetection hypothesis to use Murty (Later this weight is added).
                cost_matrix_log(i,1:length(indices_c))=weights_log-filter_updTrj.tracks{i}.weightBLog_k(index_hyp);
                cost_misdetection(i)=filter_updTrj.tracks{i}.weightBLog_k(index_hyp);
                              
                
            end
        end
        
        %New targets
        for i=Nprev_tracks+1:Nnew_tracks+Nprev_tracks
            weights_log=filter_updTrj.tracks{i}.weightBLog_k;
            index=filter_updTrj.tracks{i}.aHis{1};
            cost_matrix_log(i,index)=weights_log;
        end
        
                
        %We remove -Inf rows and columns for performing optimal assignment. We take them
        %into account for indexing later.
        %Columns that have only one value different from Inf are not fed
        %into Murty either.
               
        cos_matrix_inf=~isinf(cost_matrix_log);    
        indices_stay_column=find(sum(cos_matrix_inf,1)>1); %There is one after inequality

        indices_stay_row=find(sum(cos_matrix_inf(:,indices_stay_column),2)>0);    %There is a zero after inequality
        
        fixed_indices_stay_column=find(sum(cos_matrix_inf,1)==1); %These belong to the final hypotheses but are not fed into Murty's algorithms
        [fixed_indices_stay_row,~]=find(cos_matrix_inf(:,fixed_indices_stay_column));
        
      
        %cost_fixed is the overall cost of the hypotheses that always belong to the output 
        cost_fixed=sum(cost_matrix_log(fixed_indices_stay_row+size(cost_matrix_log,1)*(fixed_indices_stay_column'-1)));      
        cost_matrix_log_trimmed=cost_matrix_log(indices_stay_row,indices_stay_column);       
        
       
        
        globWeightLog_pred=log(filter_pred.globHypWeight(p));
        
        if(isempty(cost_matrix_log_trimmed))
            %One hypothesis: All targets are misdetected, no need for Murty                 
            opt_indices_trans=zeros(1,size(z,2));
            opt_indices_trans(fixed_indices_stay_column)=fixed_indices_stay_row; 
            globWeightLog=[globWeightLog,sum(cost_misdetection)+cost_fixed+globWeightLog_pred];
      
        else
            %Number of new global hypotheses from this global hypothesis
            kbest=ceil(Nhyp_max*filter_pred.globHypWeight(p)); 
               
            
            %We run MURTY algorithm (making use of Hungarian algorithm to
            %solve the assginment problem). We call the function by using transpose and negative value
            
            [opt_indices,~,nlcost]=kBest2DAssign(-cost_matrix_log_trimmed',kbest);
            opt_indices=opt_indices';
            nlcost=nlcost';
            
            %Optimal indices without removing Inf rows
            opt_indices_trans=Inf(size(opt_indices,1),size(z,2));           
            for i=1:size(opt_indices,1)               
                opt_indices_trans(i,indices_stay_column)=indices_stay_row(opt_indices(i,:));
                
              %We add the single trajectory hypotheses that belong to all k
               %max
               opt_indices_trans(i,fixed_indices_stay_column)=fixed_indices_stay_row; 
            end
 
            globWeightLog=[globWeightLog,-nlcost+sum(cost_misdetection)+cost_fixed+globWeightLog_pred];
        end
        
        %We write the corresponding global hypothesis
        globHypProv=zeros(size(opt_indices_trans,1),Nprev_tracks+Nnew_tracks);
        
        for i=1:Nprev_tracks
            index_hyp=filter_pred.globHyp(p,i); %Hypothesis for track i in p global hypothesis
            Nhyp_i=length(filter_pred.tracks{i}.eB);            
            for j=1:size(opt_indices_trans,1)
                index_track= find(opt_indices_trans(j,:)==i);             
                if(isempty(index_track))
                    index=index_hyp;
                else
                    index=index_hyp+Nhyp_i*index_track;                   
                end
                globHypProv(j,i)=index;
            end
        end
        
        
        for i=Nprev_tracks+1:Nnew_tracks+Nprev_tracks
            
            for j=1:size(opt_indices_trans,1)
                index_track=find(opt_indices_trans(j,:)==i);
                if(isempty(index_track))
                    index=0;
                else
                    index=1;
                end       
                globHypProv(j,i)=index; %It is either 0 or 1
            end
        end
        
        globHyp=[globHyp;globHypProv];
        
   
    end
    filter_updTrj.globHyp=globHyp;
    %Normalisation of weights of global hypotheses
    globWeight=exp(globWeightLog-max(max(globWeightLog)));
    globWeight=globWeight/sum(globWeight);
    filter_updTrj.globHypWeight=globWeight;

    [~,sortUpd]=sort(filter_updTrj.globHypWeight(),'descend');
    filter_updTrj.globHyp=filter_updTrj.globHyp(sortUpd,:);
    filter_updTrj.globHypWeight=filter_updTrj.globHypWeight(sortUpd);
end
%-------------------------------------------------------------------------
function [Hbar,selStates_i,alpha]=PPPparam(typeP,zType,Nx,NxF,H,Hbar_old,p_d)
alpha=1;
Hbar=Hbar_old;
if typeP==3
    selStates_i=(1:NxF);
    switch zType
        case 1
            alpha=(1-p_d.partial)/2;
        case 2
            Hbar=H.p2;
            alpha=(1-p_d.partial)/2;
        case 3
            alpha=p_d.partial;
    end
else
    selStates_i=(Nx+1:NxF);
    if zType==2
        Hbar=H.partial;
    end
end
if alpha<0
    warning('Negative propbability. Check the condition p_d.partial<0.5')
end

function [z_pred_j,inv_S_pred_j,K_pred_j,S_pred_j,Hlist]=MBparam(mean_j,cov_j,H,R)
inv_S_pred_j=cell(1,3);
K_pred_j=cell(1,3);
z_pred_j=cell(1,3);
S_pred_j=cell(1,3);

Hlist=cell(1,3);
Hlist{1}=H.p1;
Hlist{3}=H.full;
Hlist{2}=H.p2;

for iCase=1:3
    Hbar=Hlist{iCase};
    if iCase==3
        Rbar=R.full;
    else
        Rbar=R.partial;
    end

    z_pred_j{iCase}=Hbar*mean_j;
    S_pred_j{iCase}=Hbar*cov_j*Hbar'+Rbar;
    
    %Inverse S_pred_j
    Vs= chol(S_pred_j{iCase}); 
    inv_sqrt_S= inv(Vs); 
    inv_S_pred_j{iCase}= inv_sqrt_S*inv_sqrt_S';
    K_pred_j{iCase}=cov_j*Hbar'*inv_S_pred_j{iCase};
end
%----------------------------------------
function mDist=mahaDist(z,z_pred,inv_S_pred)
    mDist=(z-z_pred)'*inv_S_pred*(z-z_pred);