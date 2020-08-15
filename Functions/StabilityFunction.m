function [S1, S01, S001,SOutput] = StabilityFunction(Pmatrix, BIOMASS, ASS, PROD, DNAT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes stability of a soil food web. It is a recreation
% of the methods by De Ruiter et al. (1995). It uses output and the same
% input as the AllFlux.m function. 
% Lisa Ceelen, 2020, Student at University of Amsterdam.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input 
% Indices of parameters must correspond with indices of species
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Parameter data
%       Assimilation Efficiency       ASS - (kgC kgC^-1)
%       Production Efficiency         PROD - (kgC kgC^-1)
%       Natural Death Rate            DNAT - (KgC year^-1)
%   Biomass data                      BIOMASS - (KgC)
%   PredationMatrix (Pmatrix from AllFlux.m function)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3*1 array of output with S values 1, 0.1 and 0.01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization
size = length(BIOMASS);
S = [1,0.1,0.01];       
SOutput = [0,0,0];

% Make matrix without diagonal
Top = (Pmatrix*-1)./BIOMASS;
Bottom = ((transpose(PROD).*transpose(ASS)).*transpose(Pmatrix))./BIOMASS;

NoDiagOrig = Top+Bottom;
% Make diagonals
Diag1 = (-1*DNAT).*eye(size);
Diag01 = (-0.1*DNAT).*eye(size);
Diag001 = (-0.01*DNAT).*eye(size);

% Make random generator for 3 x 1000 matrices, S1, S01, S001 
for s = 1:3
    for matrix = 1:1000
    % Randomize a matrix
    NoDiag = NoDiagOrig;
        for i = 1:size
            for j = i+1:size
                Randomizer = rand*2;
                NoDiag(i,j) = NoDiag(i,j)*Randomizer;
                NoDiag(j,i) = NoDiag(j,i)*Randomizer;
            end
        end
        % Add diagonal
        Diag = (-1*S(s)*DNAT).*eye(size);

        StabilityTest = NoDiag+ Diag;

        % Compute eigenvalues
        Eigenvalues = eig(StabilityTest);
            % get real part 
            RealEig = real(Eigenvalues);
        % Check if negative and Add to counter stable/unstable
        if sum(RealEig>0)== 0 % if true it is unstable
            SOutput(s) = SOutput(s)+1;
        else
        end
    end
end

S1 = SOutput(1);
S01 = SOutput(2);
S001 = SOutput(3);

    