function [Species,CassSP,CexcrSP,CprodSP,CminSP,NassSP,NexcrSP,NprodSP,NminSP,OutputSP,Pmatrix] = ...
    AllFlux(ASS,PROD,DNAT,CN,BIOMASS,W,Species)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script calculates the C and N fluxes between species in the soil 
% food web, it is a recreation of De Ruiter et al., 1993.
% Lisa Ceelen, 2019, Student at University of Amsterdam. 
% Can take in any foodweb and compute the flux parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input 
% Indices of parameters must correspond with indices of species
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Parameter data
%       Assimilation Efficiency       ASS - (kgC kgC^-1)
%       Production Efficiency         PROD - (kgC kgC^-1)
%       Natural Death Rate            DNAT - (KgC year^-1)
%       C:N Ratio                     CN - (KgC KgN^-1)
%   Biomass data                      BIOMASS - (KgC)      
%   Feeding preferences matrix        W - No unit  
%   Species list                      Species - (Character list)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 9 X*1 arrays of totals of:      (X is number of species)
% Species, Carbon Assimilation, Carbon Excretion, Carbon Production, 
% Carbon Mineralisation, Nitrogen Assimilation, Nitrogen Excretion, 
% Nitrogen Production and Nitrogen Mineralisation.
% 
% Table with: 
% Species,CassSP,CexcrSP,CprodSP,CminSP,NassSP,NexcrSP,NprodSP,NminSP
% 
% Matrix with all feeding rates between functional groups
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%% Initializing %%%%%
NumSpecies = height(Species);
% Feeding Rates (kgC ha^-1 year^-1)% starts as zero matrix 
% like preferences matrix predator to prey %
F = zeros(1,NumSpecies);       % feeding rate pred i 
Fr = zeros(NumSpecies);      % feeding rate prey j to i

% Predeation Death rates = sum of all Feeding rates on prey(kgC ha^-1 year^-1)% 
% is a zero vector so the first calculation uses no predation, after first 
% calculation second index is used. 
P = zeros(1,NumSpecies); 
Pmatrix = Fr;
% Carbon fluxes % calculated per Fr, starts as a zero MATRIX
Cass = Fr;
Cexcr = Fr;
Cprod = Fr;
Cmin = Fr;

% Nitrogen fluxes % calculated per Fr, starts as a zero MATRIX
Nass = Fr;
Nexcr = Fr;
Nprod = Fr;
Nmin = Fr;

%%%%% Calculations %%%%%
% First calculate feeding rates % Predation
% Run as cascade from top species: 
% For each i (Predator) calculate for each j (prey) the
% biomass loss because of feeding by i. index j will
% correspond to index i for the species, so solution will be added to
% value in P[i]. 
for i = 1:NumSpecies
    
    % calculate feeding for predator i % (kgC ha^-1 year^-1)
    F(i) = (DNAT(i)*BIOMASS(i)+P(i))/(ASS(i)*PROD(i));
    
    for j = i+1:NumSpecies
        if W(i,j) > 0
        % calculate feeding on prey j by predator i % (kgC ha^-1 year^-1)
        Fr(i,j) = ((W(i,j)*BIOMASS(j))/ sum(W(i,:).*transpose(BIOMASS)))*F(i);
        Fr(isnan(Fr))=0;
        % add feeding rate to Predation rate and Predation matrix
        P(j) = P(j)+Fr(i,j);            %(kgC ha^-1 year^-1)
        Pmatrix(i,j)= Fr(i,j);          %For function output
        
        % Calculate Carbon mineralization
        Cass(i,j) = Fr(i,j)*ASS(i);         % (kgC ha^-1 year^-1)
        Cexcr(i,j) = Fr(i,j)*(1-ASS(i));    % (kgC ha^-1 year^-1) 
        
        Cprod(i,j) = Fr(i,j)*ASS(i)*PROD(i);    % (kgC ha^-1 year^-1)
        Cmin(i,j)= Fr(i,j)*ASS(i)*(1-PROD(i));  % (kgC ha^-1 year^-1)
        
        % Calculate Nitrogen mineralization
        Nass(i,j)= Fr(i,j)*(ASS(i)/CN(j));          % (kgN ha^-1 year^-1)
        Nexcr(i,j) = Fr(i,j)*((1-ASS(i))/CN(j));    % (kgN ha^-1 year^-1)
        
        Nprod(i,j) = Fr(i,j)*ASS(i)*(PROD(i)/CN(i));% (kgN ha^-1 year^-1)
        Nmin(i,j) = Fr(i,j)*ASS(i)*((1/CN(j))-(PROD(i)/CN(i)));  % (kgN ha^-1 year^-1)
        
        else
        end
    end
end
%%% OUTPUT ARRANGE %%%
% C and N per species (predator)
CassSP = sum(Cass,2);
CexcrSP = sum(Cexcr,2);
CprodSP = sum(Cprod,2);
CminSP = sum(Cmin,2);

NassSP = sum(Nass,2);
NexcrSP = sum(Nexcr,2);
NprodSP = sum(Nprod,2);
NminSP = sum(Nmin,2);

CassSP1 = array2table(CassSP);
CexcrSP1 = array2table(CexcrSP);
CprodSP1 = array2table(CprodSP);
CminSP1 = array2table(CminSP);

NassSP1 = array2table(NassSP);
NexcrSP1 = array2table(NexcrSP);
NprodSP1 = array2table(NprodSP);
NminSP1 = array2table(NminSP);

OutputSP = [Species,CassSP1,CexcrSP1,CprodSP1,CminSP1,NassSP1,NexcrSP1,NprodSP1,NminSP1];