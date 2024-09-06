var  c consumption_obs epsil Fp Fw g gdp_obs h hours_obs i inflation_obs investment_obs kbar lambdaf lambdaz
     muup muzstar phi pi pinvest_obs pitarget pstar q Re Re_obs RealRe_obs rk Rk s u wage_obs wtilde wstar 
     zetac zetai y
     cf hf i_f kbarf lambdazf phif qf Ref  rkf Rkf  sf uf Fwf wtildef Fpf;




varexo e_epsil e_lambdaf e_muup e_muzstar e_pitarget e_zetac e_zetai e_g;     // e_xp



parameters 


       actil_p adptil_p adytil_p alpha_p aptil_p aytil_p b_p beta_p c_p delta_p epsil_p etag_p g_p
           i_p iota_p iotamu_p iotaw_p lambdaf_p lambdaw_p mu_p muup_p muzstar_p pi_p pibar_p pitarget_p
           psik_p psil_p psiL_p Re_p rhoepsil_p rhog_p rholambdaf_p rhomuup_p rhomuzstar_p rhopitarget_p
           rhoterm_p rhozetac_p rhozetai_p rhotil_p Sdoupr_p sigmaL_p sigmaa_p stdepsil_p stdg_p stdlambdaf_p
           stdmuup_p stdmuzstar_p stdpitarget_p stdterm_p stdxp_p stdzetac_p stdzetai_p tauc_p taud_p tauk_p
           taul_p tauo_p term_p upsil_p we_p xip_p xiw_p zeta_p zetac_p zetai_p;

% 2. Set parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% when stopshock = 1, then non-risk shocks are all turned off
@# define stopshock = 0

% when stopsignal = 1, then signals on risk are turned off
@# define stopsignal = 0

% when stopunant = 1, then unanticipated risk shock turned off
@# define stopunant = 0

% when signal_corr_nonzero = 1, sig_corr_p can be non zero.
@# define signal_corr_nonzero = 1
@#define possible_signals = ["8"]
 %["0", "1", "2", "3", "4", "5", "6", "7", "8"]

stdlambdaf_p   =   0.010955960700000  ;  
stdmuup_p      =   0.003985791800000  ;	 
stdg_p         =   0.022812723300000  ;	 
stdmuzstar_p   =   0.007148784800000  ;	 
stdgamma_p     =   0.008103540300000  ;	 
stdepsil_p     =   0.004633811900000  ;	 
stdxp_p        =   0.489344620900000  ;	 
stdzetac_p     =   0.023325355100000  ;	 
stdzetai_p     =   0.054964824400000  ;	 
stdterm_p      =   0.001603753000000  ;	 
% Place holder for net worth 
xiw_p          =   0.812796311400000  ;	 
b_p            =   0.735843822700000  ;	 
Fomegabar_p    =   0.005588569300000  ;	 
mu_p           =   0.214894511100000  ;	 
sigmaa_p       =   2.535553419500000  ;	 
Sdoupr_p       =  10.780000003400000  ; 
xip_p          =   0.741218603400000  ;	 
aptil_p        =   2.396495942700000  ;	 
rhotil_p       =   0.850296450300000  ;	 
iota_p         =   0.897367052100000  ;	 
iotaw_p        =   0.489073535900000  ;	 
iotamu_p       =   0.936565280700000  ;	 
adytil_p       =   0.364943654300000  ;
 
rholambdaf_p   = 0.9108528528580380 ;	 
rhomuup_p      = 0.9870257396836700 ;	 
rhog_p         = 0.9427215849959780 ;	 
rhomuzstar_p   = 0.1459051086113400 ;	 
rhoepsil_p     = 0.8089285617540170 ;	 
rhosigma_p     = 0.9706370265612010 ;	 
rhozetac_p     = 0.8968400853887450 ;	 
rhozetai_p     = 0.9086616567125290 ;	 
rhoterm_p      = 0.9743991813961140 ;	 
stdsigma1_p    = 0.0700061676650730 ;	 
	
		 

// Calibrated parameters.
actil_p           = 0;
adptil_p          = 0;
alpha_p           = 0.4;
aytil_p           = 0;
beta_p            = 0.998704208591811;
bigtheta_p        = 0.005;
c_p               = 1.545858551297361;
delta_p           = 0.025;
epsil_p           = 1;
etag_p            = 0.2043;
g_p               = 0.586751768198739;
gamma_p           = 0.985;
i_p               = 0.739400293322006;
lambdaf_p         = 1.2;
lambdaw_p         = 1.05;
muup_p            = 1;
muzstar_p         = 1.004124413586981;
pi_p              = 1.006010795406775;
pibar_p           = 1.006010795406775;
pitarget_p        = 1.006010795406775;
psik_p            = 0;
psil_p            = 0;
psiL_p            = 0.7705;
Re_p              = 0.011470654984364;
rhogamma_p        = 0;
rhopitarget_p     = 0.975;
sigmaL_p          = 1;
sigma_p           = 0.327545843119697;
stdpitarget_p     = 0.0001;
tauc_p            = 0.047;
taud_p            = 0;
tauk_p            = 0.32;
taul_p            = 0.241;
tauo_p            = 1;
term_p            = 1;
upsil_p           = 1.004223171829000;
we_p              = 0.005;
zeta_p            = 1;
zetac_p           = 1;
zetai_p           = 1;



//---------------------------------------------------------------------
// 3. Model declaration
//---------------------------------------------------------------------

model;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary expressions.  These simplify the equations without adding
% additional variables.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # pitilde    = (pitarget^iota_p * pi(-1)^(1-iota_p) * pibar_p^(1-iota_p-(1-iota_p)));
  # Kp         = (Fp * ((1 - xip_p * (pitilde / pi)^(1/(1-lambdaf))) / (1 - xip_p))^(1-lambdaf));

 # Kpf        = Fpf ;
 
 # rk_ss      = (steady_state(rk)); # rkf_ss     = (steady_state(rkf));

 # a          = (rk_ss * (exp(sigmaa_p * (u - 1)) - 1) / sigmaa_p);

 # af         = (rk_ss * (exp(sigmaa_p * (uf - 1)) - 1) / sigmaa_p);

 # pitildep1  = (pitarget(+1)^iota_p * pi^(1-iota_p) * pibar_p^(1-iota_p-(1-iota_p)));
  
 # yz         = (pstar^(lambdaf/(lambdaf-1)) * (epsil * (u * kbar(-1) / (muzstar * upsil_p))^alpha_p 
                 * (h * wstar^(lambdaw_p/(lambdaw_p-1)))^(1-alpha_p) - phi));
 # yzf         = ((epsil * (uf * kbarf(-1) / (muzstar * upsil_p))^alpha_p* hf^(1-alpha_p) - phif));  

 # Kpp1       = (Fp(+1) * ((1 - xip_p * (pitildep1 / pi(+1))^(1/(1-lambdaf(+1)))) / (1 - xip_p))^(1-lambdaf(+1)));
 # pitildewp1 = (pitarget(+1)^iotaw_p * pi^(1-iotaw_p) * pibar_p^(1-iotaw_p-(1-iotaw_p)));
 # piwp1      = (pi(+1) * muzstar(+1) * wtilde(+1) / wtilde);
 # piw        = (pi * muzstar * wtilde / wtilde(-1));
 # pitildew   = (pitarget^iotaw_p * pi(-1)^(1-iotaw_p) * pibar_p^(1-iotaw_p-(1-iotaw_p)));
 # Kwp1       = (((1 - xiw_p * (pitildewp1 / piwp1 * muzstar_p^(1-iotamu_p) * muzstar(+1)^iotamu_p)^(1/(1-lambdaw_p))) 
                 / (1-xiw_p))^(1-lambdaw_p*(1+sigmaL_p)) * wtilde(+1) * Fw(+1) / psiL_p);
 
# Kwp1f       = (wtildef(+1) * Fwf(+1) / psiL_p);

# Kw         = (((1 - xiw_p * (pitildew / piw * muzstar_p^(1-iotamu_p) * muzstar^iotamu_p)^(1/(1-lambdaw_p))) 
                 / (1 - xiw_p))^(1-lambdaw_p*(1+sigmaL_p)) * wtilde * Fw / psiL_p);
  
# Kwf         = (wtildef * Fwf / psiL_p);

# S          = (exp(sqrt(Sdoupr_p / 2)*(zetai * muzstar * upsil_p * i / i(-1) - muzstar_p * upsil_p))
                 + exp(-sqrt(Sdoupr_p / 2) * (zetai * muzstar * upsil_p * i/i(-1) - muzstar_p * upsil_p)) - 2);

# Sf         = (exp(sqrt(Sdoupr_p / 2)*(zetai * muzstar * upsil_p * i_f / i_f(-1) - muzstar_p * upsil_p))
               + exp(-sqrt(Sdoupr_p / 2) * (zetai * muzstar * upsil_p * i_f/i_f(-1) - muzstar_p * upsil_p)) - 2);

# Spr        = (sqrt(Sdoupr_p / 2) * (exp(sqrt(Sdoupr_p / 2) * (zetai * muzstar * upsil_p * i / i(-1) - muzstar_p * upsil_p)) 
                 - exp(-sqrt(Sdoupr_p / 2) * (zetai * muzstar * upsil_p * i / i(-1) - muzstar_p * upsil_p))));
 
# Sprf       = (sqrt(Sdoupr_p / 2) * (exp(sqrt(Sdoupr_p / 2) * (zetai * muzstar * upsil_p * i_f / i_f(-1) - muzstar_p * upsil_p)) 
                 - exp(-sqrt(Sdoupr_p / 2) * (zetai * muzstar * upsil_p * i_f / i_f(-1) - muzstar_p * upsil_p))));

# Sprp1      = (sqrt(Sdoupr_p / 2) * (exp(sqrt(Sdoupr_p / 2) * (zetai(+1) * muzstar(+1) * upsil_p * i(+1) / i - muzstar_p * upsil_p)) 
                 - exp(-sqrt(Sdoupr_p / 2) * (zetai(+1) * muzstar(+1) * upsil_p * i(+1) / i - muzstar_p * upsil_p))));
  
# Sprp1f     = (sqrt(Sdoupr_p / 2) * (exp(sqrt(Sdoupr_p / 2) * (zetai(+1) * muzstar(+1) * upsil_p * i_f(+1) / i_f - muzstar_p * upsil_p)) 
                 - exp(-sqrt(Sdoupr_p / 2) * (zetai(+1) * muzstar(+1) * upsil_p * i_f(+1) / i_f - muzstar_p * upsil_p))));
  
# Rk_ss      = (steady_state(Rk)); 
# Rkf_ss     = (steady_state(Rkf));
  
# kbar_ss    = (steady_state(kbar));
# kbarf_ss   = (steady_state(kbarf));

# h_ss       = (steady_state(h));
# hf_ss      = (steady_state(hf));

# g_ss       = (etag_p * (steady_state(c) + steady_state(i)) / (1 - etag_p));
# gf_ss       = (etag_p * (steady_state(cf) + steady_state(i_f)) / (1 - etag_p));

# yf = g + cf + i_f / muup;
y = g + c + i / muup;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Equations characterizing equilibrium.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eqn 1: Law of motion for \latex{p^*}

pstar = ((1 - xip_p) * (Kp / Fp)^(lambdaf/(1-lambdaf)) 
            + xip_p * ((pitilde / pi) * pstar(-1))^(lambdaf/(1-lambdaf)))^((1-lambdaf)/lambdaf);
  
% Eqn 2: Law of motion for \latex{F_p}.
    
 Fp = zetac * lambdaz * yz + (pitildep1 / pi(+1))^(1/(1-lambdaf(+1))) * beta_p * xip_p * Fp(+1);
     
 Fpf = zetac *lambdazf * yzf;
  
% Eqn 3: Law of motion for \latex{K_p}
  % This error is ignored.
  %  Kp = zetac * lambdaf * lambdaz * pstar^(lambdaf/(lambdaf-1)) * yz * s 
  %       + beta_p * xip_p * (pitildep1 / pi(+1))^(lambdaf(+1)/(1-lambdaf(+1))) * Kpp1;
    
Kp  = zetac * lambdaf * lambdaz * yz * s + beta_p * xip_p * (pitildep1 / pi(+1))^(lambdaf(+1)/(1-lambdaf(+1))) * Kpp1;

Kpf = zetac * lambdaf * lambdazf * yzf * sf;
 
 % Eqn 4: Relationship between \latex{K_p} and \latex{F_p}.
  % This equation is contained in the definitions of Kp and Kpp1 in the
  % auxiliary equations.
  
  % Eqn 5: Law of motion for \latex{F_w}.
 
Fw = zetac * lambdaz * wstar^(lambdaw_p/(lambdaw_p-1)) * h * (1 - taul_p) / lambdaw_p 
         + beta_p * xiw_p * muzstar_p^((1-iotamu_p)/(1-lambdaw_p)) * (muzstar(+1)^(iotamu_p/(1-lambdaw_p)-1))
         * pitildewp1^(1/(1-lambdaw_p)) / pi(+1) * (1 / piwp1)^(lambdaw_p/(1-lambdaw_p))  *  Fw(+1);
   
Fwf = zetac * lambdazf * hf * (1 - taul_p) / lambdaw_p;
 
  % Eqn 6: Law of motion for \latex{K_w}.
 Kw = zetac * (wstar^(lambdaw_p/(lambdaw_p-1)) * h)^(1+sigmaL_p) * zeta_p + beta_p * xiw_p
         * (pitildewp1 / piwp1 * muzstar_p^(1-iotamu_p) * muzstar(+1)^iotamu_p)^(lambdaw_p*(1+sigmaL_p)/(1-lambdaw_p)) * Kwp1;
Kwf = zetac * (hf)^(1+sigmaL_p) * zeta_p;

  % Eqn 7: Relationship between \latex{F_w} and \latex{K_w}.
  % This equation is contained in the definitions of Kw and Kwp1 in the
  % auxiliary equations.
  
  % Eqn 8: Law of motion of \latex{w^*}
    wstar = ((1 - xiw_p) * ( ((1 - xiw_p * (pitildew / piw * muzstar_p^(1-iotamu_p) * muzstar^iotamu_p)^(1/(1-lambdaw_p))) 
            / (1 - xiw_p))^lambdaw_p ) + xiw_p * (pitildew / piw * muzstar_p^(1-iotamu_p) * muzstar^iotamu_p 
            * wstar(-1))^(lambdaw_p/(1-lambdaw_p)))^(1/(lambdaw_p/(1-lambdaw_p)));
  
  % Eqn 9: Efficiency condition for setting captial utilization
    rk  = tauo_p * rk_ss * exp(sigmaa_p * (u  - 1));
    rkf = tauo_p * rk_ss * exp(sigmaa_p * (uf - 1));

  % Eqn 10: Rental rate on capital 
    rk  = alpha_p * epsil * ((upsil_p * muzstar * h * wstar^(lambdaw_p/(lambdaw_p-1)) /(u  * kbar(-1)) )^(1 - alpha_p)) * s;
    rkf = alpha_p * epsil * ((upsil_p * muzstar * hf                                  /(uf * kbarf(-1)))^(1 - alpha_p)) * sf;

  % Eqn 11: Marginal Cost 
    s = (rk  / alpha_p)^alpha_p * (wtilde  / (1 - alpha_p))^(1-alpha_p) / epsil;
   sf = (rkf / alpha_p)^alpha_p * (wtildef / (1 - alpha_p))^(1-alpha_p) / epsil;


  % Eqn 12: Resource constraint
   // @# if cee == 0
   // yz = g + c + i / muup + tauo_p * a * kbar(-1) / (muzstar * upsil_p) + d + bigtheta_p * (1 - gamma) * (n - we_p) / gamma;
   // @# else
    yz  = g + c  + i   / muup + tauo_p * a  * kbar(-1)  / (muzstar * upsil_p) ;
    yzf = g + cf + i_f / muup + tauo_p * af * kbarf(-1) / (muzstar * upsil_p);

   // @# endif  

  % Eqn 13: Law of motion for capital
    kbar  = (1 - delta_p) * kbar(-1)  / (muzstar * upsil_p) + (1 - S ) * i;
    kbarf = (1 - delta_p) * kbarf(-1) / (muzstar * upsil_p) + (1 - Sf) * i_f;

  % Eqn 14: Household FOC w.r.t. risk-free bonds
    0 = beta_p * zetac(+1) * lambdaz(+1)  / (muzstar(+1) * pi(+1))  *(1 + (1 - taud_p) * Re ) - zetac * lambdaz;
    0 = beta_p * zetac(+1) * lambdazf(+1) / (muzstar(+1) * pibar_p) *(1 + (1 - taud_p) * Ref) - zetac * lambdazf;

  % Eqn 15: Household FOC w.r.t. consumption
    (1 + tauc_p) * zetac * lambdaz = muzstar * zetac / (c  * muzstar - b_p * c(-1))   - b_p * beta_p * zetac(+1) / (c(+1)  * muzstar(+1) - b_p * c);
   (1 + tauc_p) * zetac * lambdazf = muzstar * zetac / (cf * muzstar - b_p * cf(-1))  - b_p * beta_p * zetac(+1) / (cf(+1) * muzstar(+1) - b_p * cf);

  % Eqn 16: FOC for capital
   // @# if cee == 0
    %0 = (1 - Gp1 - omegabar(+1) * (1 - Fp1)) * (1 + Rk(+1)) / (1 + Re) + (1 - Fp1) / (1 - Fp1 - mu_p * omegabar(+1) 
    %    * normpdf((log(omegabar(+1)) + sigma^2 / 2) / sigma) / omegabar(+1) / sigma) * ((1 + Rk(+1)) / (1 + Re) * ((1 - mu_p) * Gp1 
    %    + omegabar(+1) * (1 - Fp1)) - 1);

  //  0 = (1 - Gammap1) * (1 + Rk(+1)) / (1 + Re) + Gammaprp1 / (Gammaprp1 - mu_p * Gprp1) * ((1 + Rk(+1)) / (1 + Re) * (Gammap1 - mu_p * Gp1) - 1);

 //   @# else
    0 = beta_p * zetac(+1) * lambdaz (+1) / (muzstar(+1) * pi(+1))  * (1 + Rk(+1))  - zetac * lambdaz;
    0 = beta_p * zetac(+1) * lambdazf(+1) / (muzstar(+1) * pibar_p) * (1 + Rkf(+1)) - zetac * lambdazf;
  //  @# endif
  
  % Eqn 17: Definition of return of entrepreneurs, Rk
    1 + Rk  = ((1 - tauk_p) * (u  * rk  - tauo_p * a ) + (1 - delta_p) * q ) * pi     / (upsil_p * q(-1) ) + tauk_p * delta_p;
    1 + Rkf = ((1 - tauk_p) * (uf * rkf - tauo_p * af) + (1 - delta_p) * qf) *pibar_p / (upsil_p * qf(-1)) + tauk_p * delta_p;
  % Eqn 18: Household FOC w.r.t. investment
    0 = - zetac * lambdaz  / muup + lambdaz  * zetac * q   * (-Spr * zetai * i * muzstar * upsil_p / i(-1) + 1 - S)
        + beta_p * zetac(+1) * lambdaz(+1) * q(+1) * Sprp1 * (zetai(+1) * i(+1) * muzstar(+1) * upsil_p / i)^2 / (muzstar(+1) * upsil_p);
  
    0 = - zetac * lambdazf / muup + lambdazf * zetac * qf * (-Sprf * zetai * i_f * muzstar * upsil_p / i_f(-1) + 1 - Sf)
        + beta_p * zetac(+1) * lambdazf(+1) * qf(+1) * Sprp1f * (zetai(+1) * i_f(+1) * muzstar(+1) * upsil_p / i_f)^2 / (muzstar(+1) * upsil_p);
  
 % Eqn 19: Definition of yz.  
  % This equation is represented in the definition of yz in the definition
  % of the auxiliary equations.
  
   
  % Eqn 20: Monetary Policy Rule
    log(Re / Re_p) = rhotil_p * log(Re(-1) / Re_p) + (1 / Re_p) * (1 - rhotil_p) 
                     * (aptil_p * pi_p * log(pi(+1) / pitarget) + (1 / 4) * adytil_p * muzstar_p * log(gdp_obs));
  % monetary policy rule with short term interest rate:
  //  log(Re / Re_p) = rhotil_p * log(Re(-1) / Re_p) + (1 / Re_p) * (1 - rhotil_p) * pi_p * log(pitarget / pi_p) 
  //                   + (1 / Re_p) * (1 - rhotil_p) * aptil_p * pi_p * (log(pi(+1)) - log(pitarget)) 
  //                   + (1 / (4 * Re_p)) * (1 - rhotil_p) * adytil_p * muzstar_p * ((c_p * log(c / c(-1)) 
  //                   + i_p * log(i / i(-1)) - i_p * log(muup / muup(-1)) + g_p * log(g / g(-1)) ) / ( (c_p+i_p)/(1-etag_p) ) 
  //                   + log(muzstar / muzstar_p)) + (1 / Re_p) * (1 - rhotil_p) * adptil_p * log(pi / pi(-1)) 
  //                   - (1 / (4 * Re_p)) * (1 - rhotil_p) * aytil_p * (c_p * log(c / c_p) + i_p * log(i / i_p) 
  //                   - i_p * log(muup) + g_p * log(g / g_p)) / ((c_p+i_p)/(1-etag_p))
  //                   + (1 - @{stopshock}) * (1 / (400 * Re_p)) * e_xp;
  
  % Eqn 21: GDP
  % This is not used.  It is only used in the manuscript in the monetary
  % policy rule.
  % # y = g + c + i / muup;
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Observation equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  consumption_obs = c / c(-1) * muzstar / muzstar_p;
  gdp_obs         = (c + i / muup + g) / (c(-1) + i(-1) / muup(-1) + g(-1)) * muzstar / muzstar_p;
  hours_obs       = h / h_ss;
  inflation_obs   =  pi / pi_p;
  investment_obs  = i / i(-1) * muzstar / muzstar_p; 
  pinvest_obs     = muup(-1) / muup;
  Re_obs          = exp(Re - Re_p);
  RealRe_obs      = ((1 + Re) / pi(+1))/((1 + Re_p) / pi_p);
  wage_obs        = wtilde / wtilde(-1) * muzstar / muzstar_p;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shock equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  log(epsil / epsil_p)       = rhoepsil_p    * log(epsil(-1) / epsil_p)       + (1 - @{stopshock}) * e_epsil;
  log(g / g_ss)              = rhog_p        * log(g(-1) / g_ss)              + (1 - @{stopshock}) * e_g;
  log(lambdaf / lambdaf_p)   = rholambdaf_p  * log(lambdaf(-1) / lambdaf_p)   + (1 - @{stopshock}) * e_lambdaf;
  log(muup / muup_p)         = rhomuup_p     * log(muup(-1) / muup_p)         + (1 - @{stopshock}) * e_muup;
  log(muzstar / muzstar_p)   = rhomuzstar_p  * log(muzstar(-1) / muzstar_p)   + (1 - @{stopshock}) * e_muzstar;
  log(pitarget / pitarget_p) = rhopitarget_p * log(pitarget(-1) / pitarget_p) + (1 - @{stopshock}) * e_pitarget;
  log(zetac / zetac_p)       = rhozetac_p    * log(zetac(-1) / zetac_p)       + (1 - @{stopshock}) * e_zetac;
  log(zetai / zetai_p)       = rhozetai_p    * log(zetai(-1) / zetai_p)       + (1 - @{stopshock}) * e_zetai;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Signal equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This ensures zero profits.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  phi = steady_state(phi);
  phif = steady_state(phif);

end;
initval;
c  =               		 1.72001;
consumption_obs=   		 1;
epsil          =   		 1;
Fp             =   		 8.09892;
Fw             =   		 2.32845;
g              =   		 0.771211;
gdp_obs        =   		 1;
h              =   		 1.08818;
hours_obs      =   		 1;
i              =   		 1.28368;
inflation_obs  =   		 1;
investment_obs =   		 1;
kbar           =   		 38.7956;
lambdaf        =   		 1.2;
lambdaz        =   		 0.557268;
muup           =   		 1;
muzstar        =   		 1.00412;
phi            =   		 0.754979;
pi             =   		 1.00601;
pinvest_obs    =   		 1;
pitarget       =   		 1.00601;
pstar          =   		 1;
q              =   		 1;
Re             =   		 0.0114707;
Re_obs         =   		 1;
RealRe_obs     =   		 1;
rk             =   		 0.0392464;
Rk             =   		 0.0114707;
s              =   		 0.833333;
u              =   		 1;
wage_obs       =   		 1;
wtilde         =   		 2.0814;
wstar          =   		 1;
zetac          =   		 1;
zetai          =   		 1;
y              =   		 3.7749;
cf             =   		 1.72001;
hf             =   		 1.08818;
i_f            =   		 1.28368;
kbarf          =   		 38.7956;
lambdazf       =   		 0.557268;
phif           =   		 0.754979;
qf             =   		 1;
Ref            =   		 0.0114707;
rkf            =   		 0.0392464;
Rkf            =   		 0.0114707;
sf             =   		 0.833333;
uf             =   		 1;
Fwf            =   		 0.438346;
wtildef        =   		 2.0814;
Fpf            =   		 2.10363;
e_g            =   		 0;

end;
steady(maxit=1000,solve_algo=1);
resid;
%check;

shocks; 
//var e_xp;       stderr stdxp_p;
var e_lambdaf;  stderr stdlambdaf_p;
var e_pitarget; stderr stdpitarget_p;
var e_muup;     stderr stdmuup_p;
var e_g;        stderr stdg_p;
var e_muzstar;  stderr stdmuzstar_p;
var e_epsil;    stderr stdepsil_p;
var e_zetac;    stderr stdzetac_p;
var e_zetai;    stderr stdzetai_p;
%var interest_; stderr 0;
%var fiscal_;  stderr (stdg_p/coffispol);
end;
stoch_simul (AR=100,IRF=0, noprint,nograph);
