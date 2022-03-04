% The environmental selection of MaOEA/AS
function [Zminn,Population,Flagg,T1,T2] = EnvironmentalSelection(Population,N,M,AveDcP_P,Zmin,Flagg,T1,T2)
    %% Normalization
    Zmax   = max(Population.objs,[],1);
    PopObj = (Population.objs-repmat(Zmin,size(Population.objs,1),1))./repmat(Zmax-Zmin,size(Population.objs,1),1);
    %% Convergence calculation
    Zmin1   = min(PopObj,[],1);
    PopObj = PopObj - repmat(Zmin1,length(Population),1);
    Con    = sqrt(sum(PopObj.^2,2));
    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(Population.objs,Population.cons,N);    
    %% Select the solutions within and before the last front
    Last = find(FrontNo<=MaxFNo);
    nLL=size(Last,2);
    %Corner_Selection---------------------------------------------------
    PopObj1=PopObj(Last,:);
    AngleCM=acos(1-pdist2(PopObj1,eye(M),'cosine'));
    [~,CM]=min(AngleCM,[],1);
    Con(Last(CM))=0; 
	%% Calculate the angle between each two solutions
    Angle = acos(1-pdist2(PopObj1,PopObj1,'cosine'));
    Angle(logical(eye(nLL))) = inf;
    maxDAPT=max(max(Angle));
    
    %% Deletion strategy   
    Del = false(1,size(PopObj1,1)); 
    dLL=nLL-N;
    while sum(Del) < dLL
        Remain   = find(~Del);
        [DAPT1,Di]=sort(Angle(Remain,Remain),2);
        [~,index]=min(DAPT1(:,1)); 
        Di2=Di(index,1);
        if Flagg==1
            if Con(Last(Remain(index))) < Con(Last(Remain(Di2))) 
                Del(Remain(Di2))=true;
            elseif Con(Last(Remain(index))) > Con(Last(Remain(Di2)))
                Del(Remain(index))=true;
            else
                if DAPT1(index,2) < DAPT1(Di2,2)
                    Del(Remain(index))=true;
                elseif DAPT1(index,2) > DAPT1(Di2,2)
                    Del(Remain(Di2))=true;
                else
                    ppp=[Di2 index];
                    Del(Remain(ppp(randperm(numel(ppp),1))))=true;
                end
            end      
        else
            if (Con(Last(Remain(index))))*DAPT1(Di2,2)/maxDAPT<=(Con(Last(Remain(Di2))))*DAPT1(index,2)/maxDAPT
                Del(Remain(Di2))=true;
            else
                Del(Remain(index))=true;
            end
        end 
    end 
    Population = Population(Last(~Del));
    
    Zminn = min(Population.objs,[],1);
    PopObj = Population.objs;
    [N,M]  = size(PopObj);
    W      = max(1e-6,PopObj./repmat(sum(PopObj,2),1,M));
    PopObj = PopObj - repmat(Zminn,N,1);
    ASF    = max(PopObj./W,[],2);
    AveDcP_Con=mean(ASF); % Calculate the average convergence
    
    %% Adaptive switching strategy 
    if Flagg==1
        if abs(AveDcP_Con-AveDcP_P) < T1
            Flagg=0;
            T1=1.4*T1;
            if T1 > 1
                T1=0.5;
            end
        end
    else
        if abs(AveDcP_Con-AveDcP_P) >= T2
            Flagg=1;
            T2=0.95*T2;
            if T2 < 0.05
                T2=0.1;
            end
        end
    end
end