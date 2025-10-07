function DrawEstimates(X_truth,t_birth,t_death,X_estimateTrj,...
    X_timeWindow,XLim,YLim,z,Nx,Nz,kVal,metadataTrj,visMode)
%Author Marco Fontana and Ángel F. García-Fernández
%This function plots the output of the filter at each time step
%XLim and YLim set the limits for the representation

noAttchedLabel=1;
metadataTrj=[metadataTrj,ones(size(metadataTrj,1),1)*kVal];%[trackID,trackID_inList,BerIdx,time]
if kVal==1
    clf
    axis([XLim(1) XLim(2),YLim(1) YLim(2)])
    xlabel('x position (m)')
    ylabel('y position (m)')
    grid on
    hold on
    for i=1:size(X_truth,1)/Nx
            plot(X_truth((i-1)*Nx+1,t_birth(i)),X_truth((i-1)*Nx+3,t_birth(i)),'Color','b','Marker','*','MarkerSize',1)
            plot(X_truth((i-1)*Nx+1,t_birth(i):t_death(i)-1),X_truth((i-1)*Nx+3,t_birth(i):t_death(i)-1),'LineStyle','-','Color','b','LineWidth',1)
            text(X_truth((i-1)*Nx+1,t_birth(i))+1,X_truth((i-1)*Nx+3,t_birth(i))+1,num2str(i),'color','b')    
            plot(X_truth((i-1)*Nx+1,t_death(i)-1),X_truth((i-1)*Nx+3,t_death(i)-1),'Color','b','Marker','+','MarkerSize',1)
    end
    
    if (~isempty(X_timeWindow))
        for iZ=1:length(X_timeWindow)
            X_t=X_timeWindow{iZ};
            for iT=1:size(X_t,2)
                if X_t(Nx+1,iT)~=0
                    plot(X_t(Nx+1,iT),X_t(Nx+3,iT),'Color','r','Marker','.','MarkerSize',12,'LineWidth',1)
                end
            end
        end
    end

    %We plot the measurements
    if ~isempty(z)
        for iZ=1:size(z,2)
            z_t=z{1,iZ};
            for iT=1:size(z_t,2)
                if z_t(Nz+1,iT)==0
                    plot(z_t(1,iT),z_t(2,iT),'Color','k','Marker','x','MarkerSize',12,'LineWidth',1)
                    if visMode==2
                        text(z_t(1,iT)+1,z_t(2,iT)+1,num2str(iT),'color','k')
                    end
                elseif z_t(1,iT)==0
                    plot(z_t(Nz+1,iT),z_t(Nz+2,iT),'Color','k','Marker','*','MarkerSize',12,'LineWidth',1)
                    if visMode==2
                        text(z_t(Nz+1,iT)+1,z_t(Nz+2,iT)+1,num2str(iT),'color','k')
                    end
                else
                    plot(z_t([1,Nz+1],iT),z_t([2,Nz+2],iT),'Color','k','LineWidth',2)
                    if visMode==2
                        text((z_t(1,iT)+z_t(Nz+1,iT))*0.5,(z_t(2,iT)+z_t(Nz+2,iT))*0.5,num2str(iT),'color','k')
                    end
                end
            end
        end
    end
end

if ~isempty(X_estimateTrj)
    X_binary=X_estimateTrj~=0;
    colLines=lines(max(metadataTrj(:,1)));
    for k=1:size(X_binary,2)
        X_binario_k=X_binary(:,k);
        X_binario_m=reshape(X_binario_k,4,length(X_binario_k)/Nx);
        indices=find(sum(X_binario_m,1)>0);
        X_filtrado_k=X_estimateTrj(:,k);
       
        for i=1:length(indices)     
            plot(X_filtrado_k(4*indices(i)-3),X_filtrado_k(4*indices(i)-1),'o','Color',colLines(metadataTrj(i,1),:),'LineWidth',2)
            if visMode>0
                text(X_filtrado_k(4*indices(i)-3)+1,X_filtrado_k(4*indices(i)-1)+1,mat2str(metadataTrj(i,1)),'Color',colLines(metadataTrj(i,1),:))
            end
        end
    end
end