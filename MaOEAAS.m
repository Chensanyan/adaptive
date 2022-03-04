classdef MaOEAAS < ALGORITHM
% <multi/many> <real/binary/permutation>
%  MaOEA/AS  My own algorithm 
  methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population = Problem.Initialization();
            Zmin       = min(Population.objs,[],1);
            Flag       = 1;
            T1         =0.008;
            T2         =0.2;
            %% Optimization
            while Algorithm.NotTerminated(Population)
                [AveDcP_P,MatingPool] = MatingSelection(Population.objs,Zmin);
                Offspring  = OperatorGA(Population(MatingPool));    
                Zmin       = min([Zmin;Offspring.objs],[],1);
                [Zmin,Population,Flag,T1,T2] = EnvironmentalSelection([Population,Offspring],Problem.N,Problem.M,AveDcP_P,Zmin,Flag,T1,T2);
            end
        end
    end
end