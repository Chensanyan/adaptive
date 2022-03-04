% The mating selection of I_MaOEA_CSS
function [AveDcP_P,MatingPool] = MatingSelection(PopObj,Zmin)
    %% Calculate the ASF value of each solution
    [N,M] = size(PopObj);
    W     = max(1e-6,PopObj./repmat(sum(PopObj,2),1,M));
    PopObj = PopObj - repmat(Zmin,N,1);
    ASF    = max(PopObj./W,[],2);
    %% Calculate the minimum angle of each solution to others
    Angle = acos(1-pdist2(PopObj,PopObj,'cosine'));
    Angle(logical(eye(N))) = inf;
    Amin  = min(Angle,[],2); 
   
    AveDcP_P=mean(ASF); % Calculate the average convergence
    
    %% Mating selection
     MatingPool = zeros(1,N);
      for i = 1 : N
        p = randperm(N,2);
        if ASF(p(1)) < ASF(p(2)) 
            MatingPool(i) = p(1);
        elseif ASF(p(1)) > ASF(p(2))
            MatingPool(i) = p(2);
        else
            if Amin(p(1)) > Amin(p(2))
                MatingPool(i) = p(1);
            elseif Amin(p(1)) < Amin(p(2))
                MatingPool(i) = p(2);
            else
                pp=[p(1) p(2)];
                MatingPool(i) =pp(randperm(numel(pp),1));
            end
        end
      end
end