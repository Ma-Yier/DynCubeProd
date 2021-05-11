function [minFlux,maxFlux,B,dynam,Brange,Trange,Frange,TMY,currentP,stoplabel,rxnList]...
    = DynCubeProd(varargin)
%DynCubeProd is a function of DynCubeProd that identifies a set of unused
%reactions for production of target metabolites.
%
%function [minFlux,maxFlux,B,dynam,Brange,Trange,Frange,TMY,currentP,stoplabel]...
%    = DynCubeProd(model, targetMet,  glucoseRxn, oxygenRxn, biomassRxn, options)
%
%INPUTS
% model     COBRA model structure containing the following required fields to perform DynCubeProd.
%   rxns                    Rxns in the model
%   mets                    Metabolites in the model
%   S                       Stoichiometric matrix (sparse)
%   b                       RHS of Sv = b (usually zeros)
%   c                       Objective coefficients
%   lb                      Lower bounds for fluxes
%   ub                      Upper bounds for fluxes
%   rev                     Reversibility of fluxes
%
% targetMet   target metabolites
%             (e.g., {'glu__D_c','gly_c'} )
% glucoseRxn  Reaction representing glucose uptake (e.g., EX_glc__D_e)
% oxygenRxn   Reaction representing oxygen uptake (e.g., EX_o2_e)
% biomassRxn  Reaction representing biomass objective function
%                       (e.g., BIOMASS_Ec_iJO1366_core_53p95M)
%
%OPTIONAL INPUTS
% GUR    Glucose uptake ratio (Default: 10)
% OUR    Oxygen uptake ratio (Default: 10)
% minGrowth   The minimum value of biomass objective function that the
%             designed strain must achieve. (Default: 0.01)
% P           paramters for the cube size (Default:16)
%
%OUTPUTS
% minFlux   The minimum values of the target metabolite production
%           obtained by FVA.
% maxFlux   The maximum values of the target metabolite production
%           obtained by FVA.
% B   Biomass growth rate after applying the candidate unused reactions
% dynam    A matrix with PXPXP size, in which the value of elements represent 
%         whether it is feasible to find a solution in this sub-cube. 
%         0 represents it is infeasible otherwise it is feasible when the 
%         element is 1.
% b    The position of found knockout strategy on biomass axis, the b-th
%     interval on biomass axis.
% t    The position of found knockout strategy on target metabolite axis,
%     the t-th interval on target metabolite axis.
% f    The position of found knockout strategy on sum of absolute values of
%     toal fluxes axis, the f-th interval on this axis.
% TMY    The theoretical maximum yield of target metabolite
% currentP    The value of P when candidate knockout strategy is found
% stoplabel   1 or 0. 1 means funciton ends before all cubes are computed,
%             0 means all cubes are computed.
% rxnList    List of found unused reactions
%
%
% May 12, 2021   Tamura Takeyuki, Ma Yier
%

s=size(varargin,2);
if size(varargin,2)<5
    error('''model'',''target'',''glucoseRxn'',''oxygenRxn''.''biomassRxn'' must be specified.')
end
model=varargin{1};
targetMet=varargin{2};
glucoseRxnID=findRxnIDs(model,varargin{3});
if glucoseRxnID==0
    error('invalid glucoseRxn name')
end
oxygenRxnID=findRxnIDs(model,varargin{4});
if oxygenRxnID==0
    error('invalid oxygenRxn name')
end
biomassRxnID=findRxnIDs(model,varargin{5});
if biomassRxnID==0
    error('invalid biomassRxn name')
end
GUR=10;
OUR=10;
minGrowth=0.01;
P=16;
for i=3:floor(s/2)
    if strcmp(varargin{2*i},'GUR')==1
        GUR=varargin{2*i+1};
    elseif strcmp(varargin{2*i},'OUR')==1
        OUR=varargin{2*i+1};
    elseif strcmp(varargin{2*i},'minGrowth')==1
        minGrowth=varargin{2*i+1};
    elseif strcmp(varargin{2*i},'P')==1
        P=varargin{2*i+1};
    else
        error('Options must be a subset of {type, GUR, OUR, minGrowth, P}')
    end
end

rxnNum=CMet(model,glucoseRxnID,oxygenRxnID,targetMet,GUR,OUR);



for i=1:size(targetMet,2)
       [minFlux(i),maxFlux(i),B(i),Brange(i),Trange(i),Frange(i),TMY(i),dynam{i},currentP(i),stoplabel(i),rxnList{i}]=...
           BTconstraintSearch(model,targetMet{i},...
           GUR,OUR,minGrowth,glucoseRxnID,oxygenRxnID,biomassRxnID,P,rxnNum); 
end


%save('DynCubeProd.mat');
end

