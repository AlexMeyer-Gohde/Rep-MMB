function [ys,params,check] = US_FV10_rep_steadystate(ys,exo,M_,options_)
% function [ys,params,check] = NK_baseline_steadystate(ys,exo,M_,options_)
% computes the steady state for the NK_baseline.mod and uses a numerical
% solver to do so
% Inputs:
%   - ys        [vector] vector of initial values for the steady state of
%                   the endogenous variables
%   - exo       [vector] vector of values for the exogenous variables
%   - M_        [structure] Dynare model structure
%   - options   [structure] Dynare options structure
%
% Output:
%   - ys        [vector] vector of steady state values for the the endogenous variables
%   - params    [vector] vector of parameter values
%   - check     [scalar] set to 0 if steady state computation worked and to
%                    1 of not (allows to impose restrictions on parameters)

% Copyright (C) 2013-2020 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

% read out parameters to access them with their name
NumberOfParameters = M_.param_nbr;
for ii = 1:NumberOfParameters
    paramname = M_.param_names{ii};
    eval([ paramname ' = M_.params(' int2str(ii) ');']);
end
% initialize indicator
check = 0;

%% Enter model equations here

options=optimset(); % set options for numerical solver

% the steady state computation follows FV (2006), section 4.1
PI=PIbar;
u=1;
q=1;
d=1;
phi=1;
m=0;
zeta=1;
LambdaYd= (LambdaA+alppha*Lambdamu)/(1-alppha);
mu_z=exp(LambdaYd);
mu_I=exp(Lambdamu);
mu_A=exp(LambdaA);

%set the parameter Lambdax
Lambdax=mu_z;

%set the parameter gammma1
gammma1=mu_z*mu_I/betta-(1-delta);
if gammma1<0 % parameter violates restriction; Preventing this cannot be implemented via prior restriction as it is a composite of different parameters and the valid prior region has unknown form
    check=1; %set failure indicator
    return; %return without updating steady states
end


r=1*gammma1;
R=1+(PI*mu_z/betta-1);

%set Rbar
Rbar=R;

PIstar=((1-thetap*PI^(-(1-epsilon)*(1-chi)))/(1-thetap))^(1/(1-epsilon));
PIstarw=((1-thetaw*PI^(-(1-chiw)*(1-eta))*mu_z^(-(1-eta)))/(1-thetaw))^(1/(1-eta));

mc=(epsilon-1)/epsilon*(1-betta*thetap*PI^((1-chi)*epsilon))/(1-betta*thetap*PI^(-(1-epsilon)*(1-chi)))*PIstar;
w=(1-alppha)*(mc*(alppha/r)^alppha)^(1/(1-alppha));
wstar=w*PIstarw;
vp=(1-thetap)/(1-thetap*PI^((1-chi)*epsilon))*PIstar^(-epsilon);
vw=(1-thetaw)/(1-thetaw*PI^((1-chiw)*eta)*mu_z^eta)*PIstarw^(-eta);
tempvaromega=alppha/(1-alppha)*w/r*mu_z*mu_I;

[ld,fval,exitflag]=fzero(@(ld)(1-betta*thetaw*mu_z^(eta-1)*PI^(-(1-chiw)*(1-eta)))/(1-betta*thetaw*mu_z^(eta*(1+gammma))*PI^(eta*(1-chiw)*(1+gammma)))...
-(eta-1)/eta*wstar/(varpsi*PIstarw^(-eta*gammma)*ld^gammma)*((1-h*mu_z^(-1))^(-1)-betta*h*(mu_z-h)^(-1))*...
((mu_A*mu_z^(-1)*vp^(-1)*tempvaromega^alppha-tempvaromega*(1-(1-delta)*(mu_z*mu_I)^(-1)))*ld-vp^(-1)*Phi)^(-1),0.25,options);
if exitflag <1
    %indicate the SS computation was not sucessful; this would also be detected by Dynare
    %setting the indicator here shows how to use this functionality to
    %filter out parameter draws
    check=1; %set failure indicator
    return; %return without updating steady states
end


l=vw*ld;
k=tempvaromega*ld;
x=(1-(1-delta)*(mu_z*mu_I)^(-1))*k;
yd=(mu_A/mu_z*k^alppha*ld^(1-alppha)-Phi)/vp;
c=(mu_A*mu_z^(-1)*vp^(-1)*tempvaromega^alppha-tempvaromega*(1-(1-delta)*(mu_z*mu_I)^(-1)))*ld-vp^(-1)*Phi;
lambda=(1-h*betta*mu_z^(-1))*(1-h/mu_z)^(-1)/c;
F=yd-1/(1-alppha)*w*ld;
f=(eta-1)/eta*wstar*PIstarw^(-eta)*lambda*ld/(1-betta*thetaw*mu_z^(eta-1)*PI^(-(1-chiw)*(1-eta)));
f2=varpsi*d*phi*PIstarw^(-eta*(1+gammma))*ld^(1+gammma)/(1-betta*thetaw*(PI^chiw/PI)^(-eta*(1+gammma))*(wstar/wstar*mu_z)^(eta*(1+gammma)));

g1=lambda*mc*yd/(1-betta*thetap*PI^((1-chi)*epsilon));
g2=epsilon/(epsilon-1)*g1;

yg = 1;

%% end own model equations

params=NaN(NumberOfParameters,1);
for iter = 1:length(M_.params) %update parameters set in the file
    eval([ 'params(' num2str(iter) ') = ' M_.param_names{iter} ';' ])
end

NumberOfEndogenousVariables = M_.orig_endo_nbr; %auxiliary variables are set automatically
for ii = 1:NumberOfEndogenousVariables
    varname = M_.endo_names{ii};
    eval(['ys(' int2str(ii) ') = ' varname ';']);
end
