//**************************************************************************
// This file has been transformed from its MMB version to be a replication 
// .mod file in June 2022.
//
//**************************************************************************
	// A New Comparative Approach to Macroeconomic Modeling and Policy Analysis
	//
	// Volker Wieland, Tobias Cwik, Gernot J. Mueller, Sebastian Schmidt and 
	// Maik Wolters
	//
	// Working Paper, 2009
	//**************************************************************************
	
	// Model: CA_ToTEM - The Bank of Canada New Quarterly Projection and Policy Analysis Model
	
	// Further references:
	// Murchison, S. and A. Rennison. 2006. ToTEM: The Bank of Canada’s New Canadian Projection
	// Model. Bank of Canada Technical Report
	
	
	
	// This file simulates the dynamic response of Bank of Canada's ToTEM model
	// to specific shocks. 
	
	
	var cbal,cbal_cor,defn_yn,defn_yn_ss,den_w,deprk,dr1rown,gbn_cor,gbn_yn,gbn_yn_ss,
	    gn_yn_ss,inff,infprow,infq,infq2,infq_pg,infq_row,infrow,infrow2,infrow_ss,gn_yn_cor
	    la,la0_cor,la5_cor,la9_cor,ladotrow,la_inv,la_tot,lc,lc0_cor,lc995_cor,lcfl,lchm,lckl,
	    lcklc,lcklcm,lcm,lcom,lcomc,lcomd,lcomdc,lcomdg,lcomdinv,lcomdx,lcomx,lcomy,lcshare,
	    lcshare_add,lcshare_sreq,lcx8,lc_cor,lc_cor_tot,leec,leeg,leeinv,leex,lforex,lforexn0_cor,
	    lforexn_cor,lforexn_cor_tot,lg,lgbn_r,lgkl,lgklc,lgklcm,lgm,lhaw,lhawc,lhawcom,lhawcom_sreq,
	    lhawc_sreq,lhawg,lhawg_sreq,lhawinv,lhawinv_sreq,lhawx,lhawx_sreq,lhaw_sreq,lindex_w,linv,
	    linvc,linvcom,linvg,linvinv,linvkl,linvklc,linvklcm,linvm,linvx,linv_add,linv_cor,linv_cor1,
	    lk,lkcom,lkg,lkinv,lktotal,lkx,lk_cor,ll,llabshare,llabshare_add,llabshare_sreq,llagrlpc,llc,
	    llcom,llg,llinv,llx,lm,lm0_cor,lmanx,lmcg_r,lmcinv_r,lmcm_r,lmcsum_r,lmcx_r,lmc_r,lmpk,lmpkcom,
	    lmpkg,lmpkinv,lmpkx,lmpl,lm_cor,lm_cor_tot,lpcomdsum_r,lpcomd_r,lpcomrow_prow,lpcom_r,lpcrow_r,
	    lpc_cor,lpc_r,lpfx_r,lpg_r,lpinvsum_r,lpinv_r,lpmanx_r,lpm_fitsum_r,lpm_fit_r,lpm_r,lprow2_cor,
	    lprow2_r,lprow_cor,lpxrow_prow,lpxrow_r,lpx_r,lp_r,lqcom_r,lqg_r,lqinv_r,lqx_r,lq_r,ltdnr_r,ltincr_r,
	    ltransf_r,lw,lwn_r,lwsum_r,lx,lxdc0_cor,lxdc_cor,lxdc_cor_tot,lxdg_cor,lxdinv_cor,lxdx_cor,lxf_cor,lxkl,
	    lxklc,lxklcm,lxm,lxw0_cor,lxw_cor,lxw_cor_tot,lxw_r,lx_cor,ly,lykl,lyrow,lyrow_cor,lyrow_gap,lyrow_sreq,
	    ly_gap,nfa,nfa_ss,num_w,pcomrow_prow_ss,pertarget,pertran,profit,r1,r1n,r1n_ss,r1n_tran,r1row,r1rown,
	    r1rown_des,r1row_ss,rfn,rf_prem,rgb,rgbn,rgbn_cor,rgbn_ss,rgb_ss,rhn,rh_prem,risk,risk_cor,risk_dyn,risk_ss,
	    rkn,target,tdn,tdn_cor,tdn_ss,tinc,tincr_yn_ss,tinc_cor,tinc_ss,transf_r_cor,transf_yn,transf_yn_ss,u,ucom,
	    ug,uinv,usact,usum,utotal,ux,w_inf;
	
	

	
	
	varexo lyrow_shk,lpcomrow_shk,lxdc_shk,lforexn_shk,lc_shk,la_shk,gn_yn_shk;%
	
	
  
	
	     
	parameters 
	           alpha_c,alpha_c2,alpha_com,alpha_g,alpha_hawc,alpha_hawcom,alpha_hawg,alpha_hawinv,alpha_hawx,
	           alpha_inv,alpha_x,a_cbal,a_gbn,a_k,a_ladotus,a_lc,a_lc2,a_lforexn,a_lg,a_linv,a_lm,a_lpcomrow,a_lprow,
	           a_lyrow,a_pc,a_r1n,a_rf_prem,a_rgbn,a_rh_prem,a_tdn,a_tinc,a_transf_r,a_xdc,a_xdg,a_xdinv,a_xdx,a_xf,a_xw,
	           b,beta,beta_hawc,beta_hawcom,beta_hawg,beta_hawinv,beta_hawx,bqy,calvo_pcomd,calvo_pcpix,calvo_pcrow,
	           calvo_pfx,calvo_pg,calvo_pinv,calvo_pm,calvo_px,calvo_w,chi,chi2,chi3,chi_k,delta,delta_c,deprk_shk,
	           dyn_pcomd,dyn_pcpix,dyn_pcrow,dyn_pfx,dyn_pg,dyn_pinv,dyn_pm,dyn_px,dyn_w,exog_lpcomrow_prow1_shk,e_w,
	           fiscal1,fiscal2,fiscal3,fiscal4,gbn_shk,gbn_yn_ss_shk,gf,gf2,gn_yn_ss_shk,habit,habitrow,hmshar,h_alpha,
	           infrow_shk,la0_shk,la5_shk,la9_shk,ladot,land,la_c,la_g,la_inv_shk,la_x,lc0_shk,lc995_shk,lforexn0_shk,
	           lm0_shk,lpcomrow_prow_shk,lprow2_shk,lprow_shk,lpxrow_prow_shk,lrelc_a,lrelc_gbn_yn_ss,lrelc_lk,lrelc_lyres,
	           lrelc_nfa,lrelc_nfa_ss,lrelc_pcom,lrelc_yrow,lrellab_a,lrellab_gbn_yn_ss,lrellab_lk,lrellab_lyres,lrellab_nfa,
	           lrellab_nfa_ss,lrellab_pcom,lrellab_yrow,lxcom_shk,lxdc0_shk,lxw0_shk,lyrow_sreq_shk,nfa_ss_shk,omega,
	           pcomrowshare,psi,psicom,psig,psiinv,psix,r1rown_shk,r1_ss,r1_ss_shk,rho,risk_shk,rk_prem,sigma,sigma_com,smooth1,
	           sreqswitch,sub,sub2,sub3,subrow,tau,tau2,theta,theta_c2,theta_com,tinc_ss_shk,totvarpi,transf_r_ss_shk,usact_shk,
	           xi,z_lag, r1n_ss_ss, infq_ss ly_gap_ss,ly_ss,pertran_shk,pertarget_shk,target_shk,lp_res_shk,lpc_shk,lxdinv_shk,lxdg_shk,
	           lxdx_shk,lxf_shk,lxw_shk,ly_res_shk,linv_shk,lx_shk,lm_shk,r1n_tran_shk,rh_prem_shk,
	           rf_prem_shk,lk_shk,transf_r_shk,tdn_shk,rgbn_shk,tinc_shk,cbal_shk,nfa_shk;
	
	
	///////
	//  Parameters for monetary policy rule
	r1n_ss_ss = 1.25599099e-02;
	infq_ss = 4.95065682e-03;
	ly_gap_ss = 0.00000000;
	ly_ss = 4.71663730e-01;
	//////
	
	alpha_c = 7.76262151e-01;
	alpha_c2 = 9.56185291e-01;
	alpha_com = 2.60000000e-01;
	alpha_g = 9.09840353e-01;
	alpha_hawc = 8.00000000e-01;
	alpha_hawcom = 8.00000000e-01;
	alpha_hawg = 8.00000000e-01;
	alpha_hawinv = 8.00000000e-01;
	alpha_hawx = 8.00000000e-01;
	alpha_inv = 5.80161502e-01;
	alpha_x = 6.08246440e-01;
	a_cbal = 9.00000000e-01;
	a_gbn = 9.00000000e-01;
	a_k = 0.00000000e+00;
	a_ladotus = 0.00000000e+00;
	a_lc = 6.70000000e-01;
	a_lc2 = -2.50000000e-01;
	a_lforexn = 6.57033000e-01;
	a_lg = 0.00000000e+00;
	a_linv = 9.00000000e-01;
	a_lm = 8.49927000e-01;
	a_lpcomrow = 0.00000000e+00;
	a_lprow = -1.00000000e-01;
	a_lyrow = 7.00000000e-01;
	a_pc = 8.50000000e-01;
	a_r1n = 0.00000000e+00;
	a_rf_prem = 0.00000000e+00;
	a_rgbn = 6.00000000e-01;
	a_rh_prem = 0.00000000e+00;
	a_tdn = 8.00000000e-01;
	a_tinc = 9.50000000e-01;
	a_transf_r = 9.50000000e-01;
	a_xdc = 7.29816000e-01;
	a_xdg = 5.00000000e-01;
	a_xdinv = 2.75003000e-01;
	a_xdx = -1.30854000e-01;
	a_xf = 8.00000000e-01;
	a_xw = 0.00000000e+00;
	b = 2.50000000e+00;
	beta = 9.92497256e-01;
	beta_hawc = 1.00000000e-01;
	beta_hawcom = 1.00000000e-01;
	beta_hawg = 1.00000000e-01;
	beta_hawinv = 1.00000000e-01;
	beta_hawx = 1.00000000e-01;
	bqy = 0.00000000e+00;
	calvo_pcomd = 7.89856000e-01;
	calvo_pcpix = 7.00000000e-01;
	calvo_pcrow = 7.00000000e-01;
	calvo_pfx = 8.50000000e-01;
	calvo_pg = 7.00000000e-01;
	calvo_pinv = 7.00000000e-01;
	calvo_pm = 9.20000000e-01;
	calvo_px = 7.00000000e-01;
	calvo_w = 8.50000000e-01;
	chi = 3.00000000e+01;
	chi2 = 5.00373300e+00;
	chi3 = 0.00000000e+00;
	chi_k = 1.00000000e+01;
	delta = 5.00000000e-01;
	delta_c = 9.18332451e-01;
	deprk_shk = 0.00000000e+00;
	dyn_pcomd = 2.03979000e-01;
	dyn_pcpix = 1.80000000e-01;
	dyn_pcrow = 6.00000000e-01;
	dyn_pfx = 9.00000000e-01;
	dyn_pg = 1.80000000e-01;
	dyn_pinv = 1.80000000e-01;
	dyn_pm = 9.00000000e-01;
	dyn_px = 1.80000000e-01;
	dyn_w = 1.80000000e-01;
	exog_lpcomrow_prow1_shk = 0.00000000e+00;
	e_w = 4.00000000e+01;
	fiscal1 = 6.99999990e-01;
	fiscal2 = 1.00000000e+00;
	fiscal3 = 0.00000000e+00;
	fiscal4 = 9.40000000e-01;
	gbn_shk = 0.00000000e+00;
	gbn_yn_ss_shk = 0.00000000e+00;
	gf = 1.20000000e+00;
	gf2 = 5.00000000e-01;
	gn_yn_ss_shk = 0.00000000e+00;
	habit = 7.50000000e-01;
	habitrow = 7.00000000e-01;
	hmshar = 2.00000000e-01;
	h_alpha = 5.00000000e-01;
	infrow_shk = 0.00000000e+00;
	la0_shk = 0.00000000e+00;
	la5_shk = 0.00000000e+00;
	la9_shk = 0.00000000e+00;
	ladot = 0.00000000e+00;
	land = 1.90589716e-01;
	la_c = -1.60155516e-01;
	la_g = -1.66835398e-01;
	la_inv_shk = 0.00000000e+00;
	la_x = -2.36117826e-01;
	lc0_shk = 0.00000000e+00;
	lc995_shk = 0.00000000e+00;
	lforexn0_shk = 0.00000000e+00;
	lm0_shk = 0.00000000e+00;
	lpcomrow_prow_shk = 0.00000000e+00;
	lprow2_shk = 0.00000000e+00;
	lprow_shk = 0.00000000e+00;
	%lpxrow_prow_shk = 0.00000000e+00;
	lpxrow_prow_shk = 0.035521;
	lrelc_a = 2.06587440e-03;
	lrelc_gbn_yn_ss = 1.56537880e-03;
	lrelc_lk = 5.27127110e-08;
	lrelc_lyres = -0.00000000e+00;
	lrelc_nfa = 3.94573480e-08;
	lrelc_nfa_ss = 3.27733850e-01;
	lrelc_pcom = 3.31821890e-02;
	lrelc_yrow = -2.58232400e-03;
	lrellab_a = -1.13717560e-01;
	lrellab_gbn_yn_ss = 1.26247010e-02;
	lrellab_lk = 4.84368800e-09;
	lrellab_lyres = 0.00000000e+00;
	lrellab_nfa = 1.04659200e-09;
	lrellab_nfa_ss = 2.27663420e-01;
	lrellab_pcom = -1.95730590e-02;
	lrellab_yrow = 1.42146930e-01;
	lxcom_shk = 0.00000000e+00;
	lxdc0_shk = 0.00000000e+00;
	lxw0_shk = 0.00000000e+00;
	lyrow_sreq_shk = 0.00000000e+00;
	nfa_ss_shk = 0.00000000e+00;
	omega = 4.70368217e-01;
	pertran_shk=0.00000000e+00;
	pertarget_shk=0.00000000e+00;
	pcomrowshare = 1.00000000e-02;
	psi = -5.44250407e-02;
	psicom = -5.44250407e-02;
	psig = -5.44250407e-02;
	psiinv = -5.44250407e-02;
	psix = -5.44250407e-02;
	r1rown_shk = 0.00000000e+00;
	r1_ss = 7.55946120e-03;
	r1_ss_shk = 0.00000000e+00;
	rho = 5.00000000e+00;
	risk_shk = 0.00000000e+00;
	rk_prem = 2.31324936e-01;
	sigma = 5.00000000e-01;
	sigma_com = 4.00000000e-01;
	smooth1 = 8.00000000e-01;
	sreqswitch = 1.00000000e+00;
	sub = 9.00000000e-01;
	sub2 = 6.00000000e-01;
	sub3 = 6.41258730e+00;
	subrow = 1.80000000e+00;
	tau = 8.26484240e-03;
	tau2 = 8.26484240e-03;
	theta = 5.00000000e-01;
	theta_c2 = 4.00000000e-01;
	theta_com = 2.04262000e-01;
	tinc_ss_shk = 0.00000000e+00;
	totvarpi = 7.50000000e-01;
	transf_r_ss_shk = 0.00000000e+00;
	usact_shk = 0.00000000e+00;
	xi = 2.30000000e-01;
	z_lag = 4.75000000e-01;
	target_shk=0.00000000e+00;
	lp_res_shk=0.00000000e+00;
	lpc_shk=0.00000000e+00;
	lxdinv_shk=0.00000000e+00;
	lxdg_shk=0.00000000e+00;
	lxdx_shk=0.00000000e+00;
	lxf_shk=0.00000000e+00;
	lxw_shk=0.00000000e+00;
	ly_res_shk=0.00000000e+00;
	linv_shk=0.00000000e+00;
	lx_shk=0.00000000e+00;
	lm_shk=0.00000000e+00;
	r1n_tran_shk=0.00000000e+00;
	rh_prem_shk=0.00000000e+00;
	rf_prem_shk=0.00000000e+00;
	lk_shk=0.00000000e+00;
	transf_r_shk=0.00000000e+00;
	tdn_shk=0.00000000e+00;
	rgbn_shk=0.00000000e+00;
	tinc_shk=0.00000000e+00;
	cbal_shk=0.00000000e+00;
	nfa_shk=0.00000000e+00;
	

	
	
	model;
	

	
	
	// Original Model Code:
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//% Consumption goods production sector %
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	 exp(lckl) = (delta^(1/sigma)*(exp(la_tot)*exp(llc))^((sigma-1)/sigma)+(1-delta)^(1/sigma)*(u*exp(lk))^((sigma-1)/sigma))^(sigma/(sigma-1)) ;
	 exp(lcklc) = (delta_c^(1/sigma_com)*exp(lckl)^((sigma_com-1)/sigma_com)+(1-delta_c)^(1/sigma_com)*exp(lcomdc)^((sigma_com-1)/sigma_com))^(sigma_com/(sigma_com-1)) ;
	 exp(lcklcm) = (alpha_c^(1/theta)*exp(la_c+lcklc)^((theta-1)/theta)+(1-alpha_c)^(1/theta)*exp(lcm-lm_cor_tot)^((theta-1)/theta))^(theta/(theta-1)) ;
	 exp(lcx8) = exp(lcklcm)-chi/2*(exp(linvc-linvc(-1))-exp(ladot))^2*exp(linvc)-psi*(1-exp(rho*(u-1)))*exp(lk)*exp(lpinv_r)-totvarpi/2*(exp(lhawc)/exp(lhawc(-1))-1)^2*exp(la_tot)-chi_k/2*(exp(linvc-lk)-4*(exp(ladot)+deprk-1))^2*exp(lk) ;
	 lpcomd_r = lmc_r+1/sigma_com*(log(1-delta_c)+lcklc-lcomdc)+1/theta*(log(alpha_c)+lcklcm-lcklc)+(theta-1)/theta*la_c ;
	 lpm_fit_r = lmc_r+1/theta*(log(1-alpha_c)+lcklcm-lcm+lm_cor_tot) ;
	 lmpl = 1/theta*(log(alpha_c)+lcklcm-lcklc)+(theta-1)/theta*la_c+1/sigma*(log(delta)+lckl-llc)+(sigma-1)/sigma*la_tot+1/sigma_com*(log(delta_c)+lcklc-lckl) ;
	 lmc_r = lwn_r-lmpl ;
	 llc = leec+lhawc ;
	 exp(lmc_r)*(exp(log(alpha_c)+lcklcm-lcklc)^(1/theta)*exp(la_c)^((theta-1)/theta)*(delta_c*exp(lcklc-lckl))^(1/sigma_com)*(exp(lckl+log(delta))/(exp(la_tot)*exp(llc)))^(1/sigma)*exp(la_tot+leec)-totvarpi/exp(lhawc)*exp(la_tot)*(exp(h_alpha*(lhawc-lhawc(-1))+(1-h_alpha)*(lhawc(-1)-lhawc(-2)))-1)*exp(h_alpha*(lhawc-lhawc(-1))+(1-h_alpha)*(lhawc(-1)-lhawc(-2)))*h_alpha)-totvarpi/exp(lhawc)*exp(lmc_r(1)+infq(1))/(1+rfn(1))*exp(la_tot(1))*(exp(h_alpha*(lhawc(1)-lhawc)+(1-h_alpha)*(lhawc-lhawc(-1)))-1)*exp(h_alpha*(lhawc(1)-lhawc)+(1-h_alpha)*(lhawc-lhawc(-1)))*(1-2*h_alpha)-totvarpi/exp(lhawc)*exp(lmc_r(2)+infq(2)+infq(1))/((1+rfn(1))*(1+rfn))*exp(la_tot(2))*(exp(h_alpha*(lhawc(2)-lhawc(1))+(1-h_alpha)*(lhawc(1)-lhawc))-1)*exp(h_alpha*(lhawc(2)-lhawc(1))+(1-h_alpha)*(lhawc(1)-lhawc))*(h_alpha-1) = exp(lwn_r) ;
	 lmpk = log(-psi*rho)+rho*(u-1)+lpinv_r ;
	 exp(lk) = ((1-deprk(-1))*exp(lk(-1))+exp(linvc(-1))/4)*(1+lk_cor) ;
	 exp(lq_r) = (exp(lmc_r(1)+infq(1))*(u(1)*exp(lmpk(1))+chi_k/2*(exp(linvc-lk)^2-(4*(exp(ladot)+deprk-1))^2)-psi*(1-exp(rho*(u(1)-1)))*exp(lpinv_r(1)))+(1-deprk(1))*exp(lq_r(1)+infq(1)))/(1+rkn) ;
	 lmpk = 1/theta*(log(alpha_c)+lcklcm-lcklc)+(theta-1)/theta*la_c+1/sigma*(log(1-delta)+lckl-lk-log(u))+1/sigma_com*(log(delta_c)+lcklc-lckl) ;
	 exp(lq_r)-exp(lpinv_r) = exp(lmc_r)*(chi_k*(exp(linvc-lk)-4*(exp(ladot)+deprk-1))+chi/2*(exp(linvc-linv_cor-linvc(-1))-exp(ladot))^2+exp(linvc-linv_cor-linvc(-1))*chi*(exp(linvc-linv_cor-linvc(-1))-exp(ladot)))-1/(1+rkn)*exp(lmc_r(1)+infq(1))*exp(linvc(1)-linvc-linv_cor)^2*chi*(exp(linvc(1)-linvc-linv_cor)-exp(ladot)) ;
	 infq = dyn_pcpix/(1+beta*dyn_pcpix)*infq(-1)+beta/(1+beta*dyn_pcpix)*infq(1)+(1-dyn_pcpix)*(1-beta)/(1+beta*dyn_pcpix)*pertarget+xi*(1-calvo_pcpix)*(1-beta*calvo_pcpix)/((1+beta*dyn_pcpix)*calvo_pcpix)*lmc_r+lxdc_cor_tot ;
	 
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//% Investment goods production sector  %
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	 exp(linvkl) = (delta^(1/sigma)*(exp(la_tot)*exp(llinv))^((sigma-1)/sigma)+(1-delta)^(1/sigma)*(uinv*exp(lkinv))^((sigma-1)/sigma))^(sigma/(sigma-1)) ;
	 exp(linvklc) = (delta_c^(1/sigma_com)*exp(linvkl)^((sigma_com-1)/sigma_com)+(1-delta_c)^(1/sigma_com)*exp(lcomdinv)^((sigma_com-1)/sigma_com))^(sigma_com/(sigma_com-1)) ;
	 exp(linvklcm) = (alpha_inv^(1/theta)*exp(la_inv+linvklc)^((theta-1)/theta)+(1-alpha_inv)^(1/theta)*exp(linvm-lm_cor_tot)^((theta-1)/theta))^(theta/(theta-1)) ;
	 exp(linv) = exp(linvklcm)-chi/2*(exp(linvinv-linvinv(-1))-exp(ladot))^2*exp(linvinv)-psiinv*(1-exp(rho*(uinv-1)))*exp(lkinv)-totvarpi/2*(exp(lhawinv)/exp(lhawinv(-1))-1)^2*exp(la_tot)-chi_k/2*(exp(linvinv-lkinv)-4*(exp(ladot)+deprk-1))^2*exp(lkinv) ;
	 lpcomd_r = lmcinv_r+1/sigma_com*(log(1-delta_c)+linvklc-lcomdinv)+1/theta*(log(alpha_inv)+linvklcm-linvklc)+(theta-1)/theta*la_inv ;
	 lpm_fit_r = lmcinv_r+1/theta*(log(1-alpha_inv)+linvklcm-linvm+lm_cor_tot) ;
	 lmcinv_r = lwn_r-(1/theta*(log(alpha_inv)+linvklcm-linvklc)+(theta-1)/theta*la_inv+1/sigma*(log(delta)+linvkl-llinv)+(sigma-1)/sigma*la_tot+1/sigma_com*(log(delta_c)+linvklc-linvkl)) ;
	 llinv = leeinv+lhawinv ;
	 exp(lmcinv_r)*(exp(log(alpha_inv)+linvklcm-linvklc)^(1/theta)*exp(la_inv)^((theta-1)/theta)*(delta_c*exp(linvklc-linvkl))^(1/sigma_com)*(exp(linvkl+log(delta))/(exp(la_tot)*exp(llinv)))^(1/sigma)*exp(la_tot+leeinv)-totvarpi/exp(lhawinv)*exp(la_tot)*(exp(h_alpha*(lhawinv-lhawinv(-1))+(1-h_alpha)*(lhawinv(-1)-lhawinv(-2)))-1)*exp(h_alpha*(lhawinv-lhawinv(-1))+(1-h_alpha)*(lhawinv(-1)-lhawinv(-2)))*h_alpha)-totvarpi/exp(lhawinv)*exp(lmcinv_r(1)+infq(1))/(1+rfn(1))*exp(la_tot(1))*(exp(h_alpha*(lhawinv(1)-lhawinv)+(1-h_alpha)*(lhawinv-lhawinv(-1)))-1)*exp(h_alpha*(lhawinv(1)-lhawinv)+(1-h_alpha)*(lhawinv-lhawinv(-1)))*(1-2*h_alpha)-totvarpi/exp(lhawinv)*exp(lmcinv_r(2)+infq(2)+infq(1))/((1+rfn(1))*(1+rfn))*exp(la_tot(2))*(exp(h_alpha*(lhawinv(2)-lhawinv(1))+(1-h_alpha)*(lhawinv(1)-lhawinv))-1)*exp(h_alpha*(lhawinv(2)-lhawinv(1))+(1-h_alpha)*(lhawinv(1)-lhawinv))*(h_alpha-1) = exp(lwn_r) ;
	 lmpkinv = log(-psiinv*rho)+rho*(uinv-1) ;
	 exp(lkinv) = ((1-deprk(-1))*exp(lkinv(-1))+exp(linvinv(-1))/4)*(1+lk_cor) ;
	 exp(lqinv_r) = (exp(lmcinv_r(1)+infq(1))*(uinv(1)*exp(lmpkinv(1))+chi_k/2*(exp(linvinv-lkinv)^2-(4*(exp(ladot)+deprk-1))^2)-psiinv*(1-exp(rho*(uinv(1)-1))))+(1-deprk(1))*exp(lqinv_r(1)+infq(1)))/(1+rkn) ;
	 lmpkinv = 1/theta*(log(alpha_inv)+linvklcm-linvklc)+(theta-1)/theta*la_inv+1/sigma*(log(1-delta)+linvkl-lkinv-log(uinv))+1/sigma_com*(log(delta_c)+linvklc-linvkl) ;
	 exp(lqinv_r)-exp(lpinv_r) = exp(lmcinv_r)*(chi_k*(exp(linvinv-lkinv)-4*(exp(ladot)+deprk-1))+chi/2*(exp(linvinv-linv_cor-linvinv(-1))-exp(ladot))^2+exp(linvinv-linv_cor-linvinv(-1))*chi*(exp(linvinv-linv_cor-linvinv(-1))-exp(ladot)))-1/(1+rkn)*exp(lmcinv_r(1)+infq(1))*exp(linvinv(1)-linvinv-linv_cor)^2*chi*(exp(linvinv(1)-linvinv-linv_cor)-exp(ladot)) ;
	 lpinv_r-lpinv_r(-1)+infq = dyn_pinv/(1+beta*dyn_pinv)*(lpinv_r(-1)+infq(-1)-lpinv_r(-2))+beta/(1+beta*dyn_pinv)*(lpinv_r(1)+infq(1)-lpinv_r)+(1-dyn_pinv)*(1-beta)/(1+beta*dyn_pinv)*pertarget+xi*(1-calvo_pinv)*(1-beta*calvo_pinv)/((1+beta*dyn_pinv)*calvo_pinv)*(lmcinv_r-lpinv_r)+lxdinv_cor ;
	 
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//% Government goods production sector  %
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	 exp(lgkl) = (delta^(1/sigma)*(exp(la_tot)*exp(llg))^((sigma-1)/sigma)+(1-delta)^(1/sigma)*(ug*exp(lkg))^((sigma-1)/sigma))^(sigma/(sigma-1)) ;
	 exp(lgklc) = (delta_c^(1/sigma_com)*exp(lgkl)^((sigma_com-1)/sigma_com)+(1-delta_c)^(1/sigma_com)*exp(lcomdg)^((sigma_com-1)/sigma_com))^(sigma_com/(sigma_com-1)) ;
	 exp(lgklcm) = (alpha_g^(1/theta)*exp(la_g+lgklc)^((theta-1)/theta)+(1-alpha_g)^(1/theta)*exp(lgm-lm_cor_tot)^((theta-1)/theta))^(theta/(theta-1)) ;
	 exp(lg) = exp(lgklcm)-chi/2*(exp(linvg-linvg(-1))-exp(ladot))^2*exp(linvg)-psig*(1-exp(rho*(ug-1)))*exp(lkg)*exp(lpinv_r-lpg_r)-totvarpi/2*(exp(lhawg)/exp(lhawg(-1))-1)^2*exp(la_tot)-chi_k/2*(exp(linvg-lkg)-4*(exp(ladot)+deprk-1))^2*exp(lkg) ;
	 lpcomd_r = lmcg_r+1/sigma_com*(log(1-delta_c)+lgklc-lcomdg)+1/theta*(log(alpha_g)+lgklcm-lgklc)+(theta-1)/theta*la_g ;
	 lpm_fit_r = lmcg_r+1/theta*(log(1-alpha_g)+lgklcm-lgm+lm_cor_tot) ;
	 lmcg_r = lwn_r-(1/theta*(log(alpha_g)+lgklcm-lgklc)+(theta-1)/theta*la_g+1/sigma*(log(delta)+lgkl-llg)+(sigma-1)/sigma*la_tot+1/sigma_com*(log(delta_c)+lgklc-lgkl)) ;
	 llg = leeg+lhawg ;
	 exp(lmcg_r)*(exp(log(alpha_g)+lgklcm-lgklc)^(1/theta)*exp(la_g)^((theta-1)/theta)*(delta_c*exp(lgklc-lgkl))^(1/sigma_com)*(exp(lgkl+log(delta))/(exp(la_tot)*exp(llg)))^(1/sigma)*exp(la_tot+leeg)-totvarpi/exp(lhawg)*exp(la_tot)*(exp(h_alpha*(lhawg-lhawg(-1))+(1-h_alpha)*(lhawg(-1)-lhawg(-2)))-1)*exp(h_alpha*(lhawg-lhawg(-1))+(1-h_alpha)*(lhawg(-1)-lhawg(-2)))*h_alpha)-totvarpi/exp(lhawg)*exp(lmcg_r(1)+infq(1))/(1+rfn(1))*exp(la_tot(1))*(exp(h_alpha*(lhawg(1)-lhawg)+(1-h_alpha)*(lhawg-lhawg(-1)))-1)*exp(h_alpha*(lhawg(1)-lhawg)+(1-h_alpha)*(lhawg-lhawg(-1)))*(1-2*h_alpha)-totvarpi/exp(lhawg)*exp(lmcg_r(2)+infq(2)+infq(1))/((1+rfn(1))*(1+rfn))*exp(la_tot(2))*(exp(h_alpha*(lhawg(2)-lhawg(1))+(1-h_alpha)*(lhawg(1)-lhawg))-1)*exp(h_alpha*(lhawg(2)-lhawg(1))+(1-h_alpha)*(lhawg(1)-lhawg))*(h_alpha-1) = exp(lwn_r) ;
	 lmpkg = log(-psig*rho)+rho*(ug-1)-lpg_r+lpinv_r ;
	 exp(lkg) = ((1-deprk(-1))*exp(lkg(-1))+exp(linvg(-1))/4)*(1+lk_cor) ;
	 exp(lqg_r) = (exp(lmcg_r(1)+infq(1))*(ug(1)*exp(lmpkg(1))+chi_k/2*(exp(linvg-lkg)^2-(4*(exp(ladot)+deprk-1))^2)-psig*(1-exp(rho*(ug(1)-1)))*exp(lpinv_r(1)-lpg_r(1)))+(1-deprk(1))*exp(lqg_r(1)+infq(1)))/(1+rkn) ;
	 lmpkg = 1/theta*(log(alpha_g)+lgklcm-lgklc)+(theta-1)/theta*la_g+1/sigma*(log(1-delta)+lgkl-lkg-log(ug))+1/sigma_com*(log(delta_c)+lgklc-lgkl) ;
	 exp(lqg_r)-exp(lpinv_r) = exp(lmcg_r)*(chi_k*(exp(linvg-lkg)-4*(exp(ladot)+deprk-1))+chi/2*(exp(linvg-linv_cor-linvg(-1))-exp(ladot))^2+exp(linvg-linv_cor-linvg(-1))*chi*(exp(linvg-linv_cor-linvg(-1))-exp(ladot)))-1/(1+rkn)*exp(lmcg_r(1)+infq(1))*exp(linvg(1)-linvg-linv_cor)^2*chi*(exp(linvg(1)-linvg-linv_cor)-exp(ladot)) ;
	 infq_pg = dyn_pg/(1+beta*dyn_pg)*infq_pg(-1)+beta/(1+beta*dyn_pg)*infq_pg(1)+(1-dyn_pg)*(1-beta)/(1+beta*dyn_pg)*pertarget+xi*(1-calvo_pg)*(1-beta*calvo_pg)/((1+beta*dyn_pg)*calvo_pg)*(lmcg_r-lpg_r)+lxdg_cor ;
	 
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%  Export goods production sector %
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	 exp(lxkl) = (delta^(1/sigma)*(exp(la_tot)*exp(llx))^((sigma-1)/sigma)+(1-delta)^(1/sigma)*(ux*exp(lkx))^((sigma-1)/sigma))^(sigma/(sigma-1)) ;
	 exp(lxklc) = (delta_c^(1/sigma_com)*exp(lxkl)^((sigma_com-1)/sigma_com)+(1-delta_c)^(1/sigma_com)*exp(lcomdx)^((sigma_com-1)/sigma_com))^(sigma_com/(sigma_com-1)) ;
	 exp(lxklcm) = (alpha_x^(1/theta)*exp(la_x+lxklc)^((theta-1)/theta)+(1-alpha_x)^(1/theta)*exp(lxm-lm_cor_tot)^((theta-1)/theta))^(theta/(theta-1)) ;
	 exp(lmanx) = exp(lxklcm)-chi/2*(exp(linvx-linvx(-1))-exp(ladot))^2*exp(linvx)-psix*(1-exp(rho*(ux-1)))*exp(lkx)*exp(lpinv_r-lpmanx_r)-totvarpi/2*(exp(lhawx)/exp(lhawx(-1))-1)^2*exp(la_tot)-chi_k/2*(exp(linvx-lkx)-4*(exp(ladot)+deprk-1))^2*exp(lkx) ;
	 lpcomd_r = lmcx_r+1/sigma_com*(log(1-delta_c)+lxklc-lcomdx)+1/theta*(log(alpha_x)+lxklcm-lxklc)+(theta-1)/theta*la_x ;
	 lpm_fit_r = lmcx_r+1/theta*(log(1-alpha_x)+lxklcm-lxm+lm_cor_tot) ;
	 lmcx_r = lwn_r-(1/theta*(log(alpha_x)+lxklcm-lxklc)+(theta-1)/theta*la_x+1/sigma*(log(delta)+lxkl-llx)+(sigma-1)/sigma*la_tot+1/sigma_com*(log(delta_c)+lxklc-lxkl)) ;
	 llx = leex+lhawx ;
	 exp(lmcx_r)*(exp(log(alpha_x)+lxklcm-lxklc)^(1/theta)*exp(la_x)^((theta-1)/theta)*(delta_c*exp(lxklc-lxkl))^(1/sigma_com)*(exp(lxkl+log(delta))/(exp(la_tot)*exp(llx)))^(1/sigma)*exp(la_tot+leex)-totvarpi/exp(lhawx)*exp(la_tot)*(exp(h_alpha*(lhawx-lhawx(-1))+(1-h_alpha)*(lhawx(-1)-lhawx(-2)))-1)*exp(h_alpha*(lhawx-lhawx(-1))+(1-h_alpha)*(lhawx(-1)-lhawx(-2)))*h_alpha)-totvarpi/exp(lhawx)*exp(lmcx_r(1)+infq(1))/(1+rfn(1))*exp(la_tot(1))*(exp(h_alpha*(lhawx(1)-lhawx)+(1-h_alpha)*(lhawx-lhawx(-1)))-1)*exp(h_alpha*(lhawx(1)-lhawx)+(1-h_alpha)*(lhawx-lhawx(-1)))*(1-2*h_alpha)-totvarpi/exp(lhawx)*exp(lmcx_r(2)+infq(2)+infq(1))/((1+rfn(1))*(1+rfn))*exp(la_tot(2))*(exp(h_alpha*(lhawx(2)-lhawx(1))+(1-h_alpha)*(lhawx(1)-lhawx))-1)*exp(h_alpha*(lhawx(2)-lhawx(1))+(1-h_alpha)*(lhawx(1)-lhawx))*(h_alpha-1) = exp(lwn_r) ;
	 lmpkx = log(-psix*rho)+rho*(ux-1)-lpmanx_r+lpinv_r ;
	 exp(lkx) = ((1-deprk(-1))*exp(lkx(-1))+exp(linvx(-1))/4)*(1+lk_cor) ;
	 exp(lqx_r) = (exp(lmcx_r(1)+infq(1))*(ux(1)*exp(lmpkx(1))+chi_k/2*(exp(linvx-lkx)^2-(4*(exp(ladot)+deprk-1))^2)-psix*(1-exp(rho*(ux(1)-1)))*exp(lpinv_r(1)-lpmanx_r(1)))+(1-deprk(1))*exp(lqx_r(1)+infq(1)))/(1+rkn) ;
	 lmpkx = 1/theta*(log(alpha_x)+lxklcm-lxklc)+(theta-1)/theta*la_x+1/sigma*(log(1-delta)+lxkl-lkx-log(ux))+1/sigma_com*(log(delta_c)+lxklc-lxkl) ;
	 exp(lqx_r)-exp(lpinv_r) = exp(lmcx_r)*(chi_k*(exp(linvx-lkx)-4*(exp(ladot)+deprk-1))+chi/2*(exp(linvx-linv_cor-linvx(-1))-exp(ladot))^2+exp(linvx-linv_cor-linvx(-1))*chi*(exp(linvx-linv_cor-linvx(-1))-exp(ladot)))-1/(1+rkn)*exp(lmcx_r(1)+infq(1))*exp(linvx(1)-linvx-linv_cor)^2*chi*(exp(linvx(1)-linvx-linv_cor)-exp(ladot)) ;
	 lpmanx_r+infq-lpmanx_r(-1) = dyn_px/(1+beta*dyn_px)*(lpmanx_r(-1)+infq(-1)-lpmanx_r(-2))+beta/(1+beta*dyn_px)*(lpmanx_r(1)+infq(1)-lpmanx_r)+(1-dyn_px)*(1-beta)/(1+beta*dyn_px)*pertarget+xi*(1-calvo_px)*(1-beta*calvo_px)/((1+beta*dyn_px)*calvo_px)*(lmcx_r-lpmanx_r)+lxdx_cor ;
	 lmanx = log(omega)+gf*(lprow2_r-lpfx_r)+lyrow+gf2*lyrow_gap+lx_cor ;
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//% Commodity goods production sector %
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	 exp(lykl) = (delta^(1/sigma)*(exp(la_tot)*exp(llcom))^((sigma-1)/sigma)+(1-delta)^(1/sigma)*exp(log(ucom)+lkcom)^((sigma-1)/sigma))^(sigma/(sigma-1)) ;
	 exp(lcomy)^((theta_com-1)/theta_com) = alpha_com^(1/theta_com)*exp(lykl)^((theta_com-1)/theta_com)+(1-alpha_com)^(1/theta_com)*(exp(la_tot)*land)^((theta_com-1)/theta_com) ;
	 exp(lcom) = exp(lcomy)-chi/2*(exp(linvcom-linvcom(-1))-exp(ladot))^2*exp(linvcom)-chi2/2*(exp(llcom-llcom(-1))-1)^2*exp(la_tot)-psicom*(1-exp(rho*(ucom-1)))*exp(lkcom)*exp(lpinv_r-lpcom_r)-chi_k/2*(exp(linvcom-lkcom)-4*(exp(ladot)+deprk-1))^2*exp(lkcom) ;
	 lpcom_r = lp_r+lpcomrow_prow+lforex ;
	 exp(lwn_r) = exp(lpcom_r)*exp(log(alpha_com)+lcomy-lykl)^(1/theta_com)*exp(log(delta)+lykl-llcom)^(1/sigma)*exp(la_tot)^(1-1/sigma)-exp(lpcom_r)*chi2*(exp(llcom-llcom(-1))-1)/exp(la_tot+llcom(-1))+exp(lpcom_r(1)+infq(1))/(1+rfn)*chi2/2*(exp(llcom(1)-llcom)-1)*exp(llcom(1)-llcom)/exp(llcom)*exp(la_tot(1)) ;
	 exp(llcom) = exp(lhawcom) ;
	 lmpkcom = log(-psicom*rho)+rho*(ucom-1)-lpcom_r+lpinv_r ;
	 exp(lkcom) = ((1-deprk(-1))*exp(lkcom(-1))+exp(linvcom(-1))/4)*(1+lk_cor) ;
	 exp(lqcom_r) = (exp(lpcom_r(1)+infq(1))*(exp(log(ucom(1))+lmpkcom(1))+chi_k/2*(exp(linvcom-lkcom)^2-(4*(exp(ladot)+deprk-1))^2)-exp(lpinv_r(1)-lpcom_r(1))*psicom*(1-exp(rho*(ucom(1)-1))))+(1-deprk(1))*exp(lqcom_r(1)+infq(1)))/(1+rkn) ;
	 exp(lmpkcom) = exp(log(alpha_com)+lcomy-lykl)^(1/theta_com)*exp(log(1-delta)+lykl-lkcom-log(ucom))^(1/sigma) ;
	 exp(lqcom_r)-exp(lpinv_r) = exp(lpcom_r)*(chi_k*(exp(linvcom-lkcom)-4*(exp(ladot)+deprk-1))+chi/2*(exp(linvcom-linv_cor-linvcom(-1))-exp(ladot))^2+exp(linvcom-linv_cor-linvcom(-1))*chi*(exp(linvcom-linv_cor-linvcom(-1))-exp(ladot)))-1/(1+rkn)*exp(lpcom_r(1)+infq(1))*exp(linvcom(1)-linvcom-linv_cor)*chi*(exp(linvcom(1)-linvcom-linv_cor)-exp(ladot)) ;
	 lpcomd_r+infq-lpcomd_r(-1) = dyn_pcomd/(1+beta*dyn_pcomd)*(lpcomd_r(-1)+infq(-1)-lpcomd_r(-2))+beta/(1+beta*dyn_pcomd)*(lpcomd_r(1)+infq(1)-lpcomd_r)+(1-dyn_pcomd)*(1-beta)/(1+beta*dyn_pcomd)*pertarget+(1-calvo_pcomd)*(1-beta*calvo_pcomd)/((1+beta*dyn_pcomd)*calvo_pcomd)*(lpcom_r-lpcomd_r) ;
	 
	//%%%%%%%%%%%%%%%%%
	//% Import sector %
	//%%%%%%%%%%%%%%%%%
	 lpm_fit_r+infq-lpm_fit_r(-1) = dyn_pm/(1+beta*dyn_pm)*(lpm_fit_r(-1)+infq(-1)-lpm_fit_r(-2))+beta/(1+beta*dyn_pm)*(lpm_fit_r(1)-lpm_fit_r+infq(1))+(1-dyn_pm)*(1-beta)/(1+beta*dyn_pm)*pertarget+(1-calvo_pm)*(1-beta*calvo_pm)/((1+beta*dyn_pm)*calvo_pm)*(lmcm_r-lpm_fit_r) ;
	 exp(lmcm_r) = exp(lpm_r) ;
	 lpm_r = lpxrow_r+lp_r+lforex+lxf_cor ;
	
	//%%%%%%%%%%%%%%%%%%%%%%%%
	//% Foreign import price %
	//%%%%%%%%%%%%%%%%%%%%%%%%
	 lpfx_r+infq_row-lpfx_r(-1) = dyn_pfx/(1+beta*dyn_pfx)*(lpfx_r(-1)+infq_row(-1)-lpfx_r(-2))+beta/(1+beta*dyn_pfx)*(lpfx_r(1)+infq_row(1)-lpfx_r)+(1-calvo_pfx)*(1-beta*calvo_pfx)/((1+beta*dyn_pfx)*calvo_pfx)*(lpmanx_r-lforex-lp_r-lpfx_r) ;
	
	//%%%%%%%
	//% UIP %
	//%%%%%%%
	 lforex = z_lag*(lforex(-1)+infq_row+lp_r(-1)-lp_r-infq+2*(pertarget-infrow_ss))+(1-z_lag)*(lforex(1)-infq_row(1)+lp_r(1)-lp_r+infq(1))+log(1+r1rown)+log(1+risk)-log(1+r1n)-log(1+r1row_ss)-log(1+risk_ss)+log(1+r1_ss) ;
	 
	//%%%%%%%%%%%%%%
	//% Households %
	//%%%%%%%%%%%%%%
	 lxw_r-lwn_r = sub2/(e_w+sub2)*(log(sub3*e_w/(e_w-1))+num_w-den_w) ;
	 exp(lwn_r) = ((1-calvo_w)*exp(lxw_r)^(1-e_w)+calvo_w*exp(lwn_r(-1)-infq+lindex_w)^(1-e_w))^(1/(1-e_w))*exp(lxw_cor_tot) ;
	 lindex_w = dyn_w*(lwn_r(-1)+infq(-1)-lwn_r(-2))+(1-dyn_w)*(ladot+pertarget) ;
	 exp(num_w) = exp(ll)^(1+1/sub2)*exp(lpc_r)*(exp(lcfl)-habit*exp(lcfl(-1)))*(1+tinc)/exp(lwn_r)+1/(beta*exp(lc_cor_tot))*calvo_w*exp(lwn_r(1)+infq(1)-lwn_r-lindex_w)^(e_w*(1+1/sub2))*exp(lwn_r(1)+infq(1)-lwn_r)*exp(llagrlpc(1)-lpc_r(1)-infq(1)-llagrlpc+lpc_r)*exp(num_w(1)) ;
	 exp(den_w) = exp(lhaw)*(1-tdn)*(1-hmshar)+1/(beta*exp(lc_cor_tot))*calvo_w*exp(lwn_r(1)+infq(1)-lwn_r-lindex_w)*exp(lindex_w)*exp(llagrlpc(1)-llagrlpc-infq(1)-lpc_r(1)+lpc_r)*exp(den_w(1)) ;
	 llagrlpc+log(1+tinc) = -1/sub*log(exp(lcfl)-habit*exp(lcfl(-1)))+(1/sub-1)/(1+1/sub2)*sub3*exp(ll)^(1+1/sub2) ;
	 log(beta*exp(lc_cor_tot)*(1+rhn)) = llagrlpc-lpc_r-llagrlpc(1)+lpc_r(1)+infq(1) ;
	 lcx8 = log(alpha_c2)-theta_c2*(-lpc_r)+lc+(1-theta_c2)*lpc_cor ;
	 lcomc = log(1-alpha_c2)-theta_c2*(lpcomd_r-lpc_r)+lc ;
	 exp(lpc_r) = (alpha_c2*exp(lpc_cor)^(1-theta_c2)+(1-alpha_c2)*exp(lpcomd_r)^(1-theta_c2))^(1/(1-theta_c2)) ;
	 
	//%%%%%%%%%%%%%%%%%%%
	//% Current Account %
	//%%%%%%%%%%%%%%%%%%%
	 cbal = exp(lforex(1)-infq_row(1)+lp_r(1)+infq(1)-lforex-lp_r)*(1+risk)*(r1rown*nfa(-1)*4*exp(lp_r(-1)-infq+ly(-1))/exp(lp_r+ly)+(1+r1rown)*(exp(lpx_r+lx)-exp(lpm_r+lm))/exp(lp_r+ly))+cbal_cor ;
	 nfa = exp(lforex(1)-infq_row(1)+infq(1)+lp_r(1)-lforex-lp_r)*(1+risk)*nfa(-1)*exp(lp_r(-1)-infq+ly(-1))/exp(lp_r+ly)+0.25*cbal+nfa_shk ;
	 risk = risk_ss+risk_dyn ;
	 risk_dyn = tau2*(exp(-(nfa-nfa_ss))-1)+lforexn_cor_tot+risk_cor ;
	 risk_ss = tau*(exp(-nfa_ss)-1) ;
	 
	//%%%%%%%%%%%%%%
	//% Government %
	//%%%%%%%%%%%%%%
	 lg = fiscal4*(ladot+lg(-1))+(1-fiscal4)*(log(gn_yn_ss)-lpg_r+lp_r+ly)+gn_yn_cor ;%
	 lgbn_r = log((1+rgbn)^0.25*exp(lgbn_r(-1)-infq)+(exp(lg+lpg_r)+exp(ltransf_r)-exp(ltdnr_r)-exp(ltincr_r))/4)+gbn_cor ;
	 gbn_yn = exp(lgbn_r)/exp(ly+lp_r) ;
	 defn_yn = (exp(lg+lpg_r)+exp(ltransf_r)-exp(ltdnr_r)-exp(ltincr_r))/exp(ly+lp_r) ;
	 defn_yn_ss = 4*(1-(1+rgbn_ss)^0.25/exp(ladot+pertarget))*gbn_yn_ss ;
	 transf_yn = transf_yn_ss+transf_r_cor+fiscal3*ly_gap ;
	 exp(ltransf_r) = transf_yn*exp(ly+lp_r) ;
	 tdn = fiscal1*tdn(-1)+(1-fiscal1)*(tdn_ss+fiscal2*(gbn_yn(1)-gbn_yn_ss(1)))+tdn_cor ;
	 tdn_ss = (defn_yn_ss-gn_yn_ss-transf_yn_ss+tincr_yn_ss)*(-1)/llabshare_sreq ;
	 exp(ltdnr_r) = tdn*exp(lwn_r+lhaw) ;
	 tinc = tinc_ss+tinc_cor ;
	 tinc_ss = tincr_yn_ss/lcshare_sreq ;
	 exp(ltincr_r) = tinc*exp(lpc_r+lc) ;
	
	//%%%%%%%%%%%%%%%
	//% Definitions %
	//%%%%%%%%%%%%%%%
	 lcshare = exp(lc+lpc_r)/exp(ly+lp_r) ;
	 lcshare_sreq = lcshare_sreq(-1)+lcshare_add ;
	 lcshare_add = lrelc_pcom*lpcomrow_prow_shk+lrelc_a*la_shk+lrelc_yrow*lyrow_sreq_shk+lrelc_lyres*ly_res_shk+lrelc_lk*lk_shk+lrelc_nfa*nfa_shk+lrelc_nfa_ss*nfa_ss_shk+lrelc_gbn_yn_ss*gbn_yn_ss_shk ;
	 llabshare = exp(lhaw+lwn_r)/exp(ly+lp_r) ;
	 llabshare_sreq = llabshare_sreq(-1)+llabshare_add ;
	 llabshare_add = lrellab_pcom*lpcomrow_prow_shk+lrellab_a*la_shk+lrellab_yrow*lyrow_sreq_shk+lrellab_lyres*ly_res_shk+lrellab_lk*lk_shk+lrellab_nfa*nfa_shk+lrellab_nfa_ss*nfa_ss_shk+lrellab_gbn_yn_ss*gbn_yn_ss_shk ;
	 exp(lchm) = hmshar*(exp(lwn_r+lhaw)-exp(ltdnr_r)+exp(ltransf_r))/((1+tinc)*exp(lpc_r)) ;
	 r1n = smooth1*r1n(-1)+(1-smooth1)*r1n_ss+(1-smooth1)*b*(infq2-pertarget-pertran)+(1-smooth1)*bqy*ly_gap ;
	 infq2 = infq(2) ;
	 ly_gap = 0.126*log(u)+0.017*log(uinv)+0.06*log(ug)+0.04*log(ux)+0.3*(llc-lhawc_sreq)+0.04*(llinv-lhawinv_sreq)+0.14*(llg-lhawg_sreq)+0.09*(llx-lhawx_sreq)+0.01*(llcom-lhawcom_sreq) ;
	 lhawc_sreq = alpha_hawc*lhawc_sreq(-1)+beta_hawc*lhawc_sreq(+1)+(1-alpha_hawc-beta_hawc)*lhawc ;
	 lhawinv_sreq = alpha_hawinv*lhawinv_sreq(-1)+beta_hawinv*lhawinv_sreq(+1)+(1-alpha_hawinv-beta_hawinv)*lhawinv ;
	 lhawg_sreq = alpha_hawg*lhawg_sreq(-1)+beta_hawg*lhawg_sreq(+1)+(1-alpha_hawg-beta_hawg)*lhawg ;
	 lhawx_sreq = alpha_hawx*lhawx_sreq(-1)+beta_hawx*lhawx_sreq(+1)+(1-alpha_hawx-beta_hawx)*lhawx ;
	 lhawcom_sreq = alpha_hawcom*lhawcom_sreq(-1)+beta_hawcom*lhawcom_sreq(+1)+(1-alpha_hawcom-beta_hawcom)*lhawcom ;
	 lmcsum_r = lmc_r+beta*lmcsum_r(+1) ;
	 lpm_fitsum_r = lpm_fit_r+beta*lpm_fitsum_r(+1) ;
	 lpinvsum_r = lpinv_r+beta*lpinvsum_r(+1) ;
	 lpcomdsum_r = lpcomd_r+beta*lpcomdsum_r(+1) ;
	 lwsum_r = lwn_r-la+beta*lwsum_r(+1) ;
	 usum = u+beta*usum(+1) ;
	 log(1+r1n) = log(1+r1)+lpc_r(1)+infq(1)-lpc_r ;
	 1+r1n_ss = (1+r1_ss)*exp(pertarget) ;
	 1+rkn = (1+rfn)*(1+rk_prem) ;
	 1+rhn = (1+r1n)*(1+rh_prem) ;
	 1+rfn = (1+r1n)*(1+rf_prem) ;
	 1+rgb = (1+rgbn)/exp(pertarget)^4 ;
	 1+rgbn = (1+r1n+0.01)^4+rgbn_cor ;
	 1+rgb_ss = (1+rgbn_ss)/exp(pertarget)^4 ;
	 1+rgbn_ss = (1+r1n_ss+0.01)^4 ;
	 inff = 4*infq ;
	 infq_pg = lpg_r+infq-lpg_r(-1) ;
	 w_inf = 4*(lwn_r+infq-lwn_r(-1)) ;
	 lw = lwn_r-lp_r ;
	 profit = 1-exp(lwn_r+lhaw-lp_r-ly) ;
	
	//%%%%%%%%%%%%%%%%%%%
	//% Market clearing %
	//%%%%%%%%%%%%%%%%%%%
	 exp(ly+lp_r) = exp(lpc_r+lc)+exp(lpinv_r+linv)+exp(lpx_r+lx)-exp(lpm_r+lm)+exp(lpg_r+lg)-lp_res_shk ;
	 exp(ly) = exp(ly(-1))*((exp(lpc_r(-1))*exp(lc)+exp(lpinv_r(-1))*exp(linv)+exp(lpx_r(-1))*exp(lx)-exp(lpm_r(-1))*exp(lm)+exp(lpg_r(-1))*exp(lg))/(exp(lpc_r(-1))*exp(lc(-1))+exp(lpinv_r(-1))*exp(linv(-1))+exp(lpx_r(-1))*exp(lx(-1))-exp(lpm_r(-1))*exp(lm(-1))+exp(lpg_r(-1))*exp(lg(-1)))*(exp(lpc_r)*exp(lc)+exp(lpinv_r)*exp(linv)+exp(lpx_r)*exp(lx)-exp(lpm_r)*exp(lm)+exp(lpg_r)*exp(lg))/(exp(lpc_r)*exp(lc(-1))+exp(lpinv_r)*exp(linv(-1))+exp(lpx_r)*exp(lx(-1))-exp(lpm_r)*exp(lm(-1))+exp(lpg_r)*exp(lg(-1))))^0.5*exp(ly_res_shk) ;
	 exp(lc) = exp(lcfl)+exp(lchm) ;
	 exp(linv) = exp(linvc)+exp(linvcom)+exp(linvx)+exp(linvg)+exp(linvinv) ;
	 exp(lpx_r+lx) = exp(lpmanx_r+lmanx)+exp(lcomx+lpcom_r) ;
	 exp(lx) = exp(lcomx)+exp(lmanx) ;
	 exp(lm) = exp(lcm)+exp(lxm)+exp(lgm)+exp(linvm) ;
	 exp(lcom) = exp(lcomd)+exp(lcomx) ;
	 exp(lcomd) = exp(lcomdc)+exp(lcomc)+exp(lcomdx)+exp(lcomdg)+exp(lcomdinv) ;
	 exp(ll) = exp(llcom)+exp(llc)+exp(llx)+exp(llg)+exp(llinv) ;
	 exp(lhaw) = exp(lhawc)+exp(lhawcom)+exp(lhawx)+exp(lhawg)+exp(lhawinv) ;
	 exp(lhaw_sreq) = exp(lhawc_sreq)+exp(lhawcom_sreq)+exp(lhawx_sreq)+exp(lhawg_sreq)+exp(lhawinv_sreq) ;
	 exp(lktotal+log(utotal)) = exp(lk+log(u))+exp(lkinv+log(uinv))+exp(lkg+log(ug))+exp(lkx+log(ux))+exp(lkcom+log(ucom)) ;
	 exp(lktotal) = exp(lk)+exp(lkinv)+exp(lkg)+exp(lkx)+exp(lkcom) ;
	 
	 
	//%%%%%%%%%%%%%%%%%%%%%%%
	//% Exogenous Processes %
	//%%%%%%%%%%%%%%%%%%%%%%%
	 la = la(-1)+ladot+la_shk+a_ladotus*lyrow_sreq_shk ;
	 la_tot = la+la0_cor+la5_cor+la9_cor ;
	 la0_cor = la0_shk ;
	 la5_cor = 0.5*la5_cor(-1)+la5_shk ;
	 la9_cor = 0.9*la9_cor(-1)+la9_shk ;
	 la_inv = la_inv(-1)+la_inv_shk ;
	 lpc_cor = a_pc*lpc_cor(-1)+lpc_shk ;
	 lxdc_cor_tot = lxdc_cor+lxdc0_cor ;
	 lxdc_cor = a_xdc*lxdc_cor(-1)+lxdc_shk ;
	 lxdc0_cor = lxdc0_shk ;
	 lxdinv_cor = a_xdinv*lxdinv_cor(-1)+lxdinv_shk ;
	 lxdx_cor = a_xdx*lxdx_cor(-1)+lxdx_shk ;
	 lxdg_cor = a_xdg*lxdg_cor(-1)+lxdg_shk ;
	 lxf_cor = a_xf*lxf_cor(-1)+lxf_shk ;
	 lxw_cor_tot = lxw_cor+lxw0_cor ;
	 lxw_cor = a_xw*lxw_cor(-1)+lxw_shk ;
	 lxw0_cor = lxw0_shk ;
	 pertarget = pertarget(-1) - pertarget_shk ;
	 target = target(-1)+target_shk ;
	 pertran = a_r1n*pertran(-1)+pertran_shk ;
	 r1n_tran = a_r1n*r1n_tran(-1)+r1n_tran_shk ;
	 rh_prem = a_rh_prem*rh_prem(-1)+rh_prem_shk ;
	 rf_prem = a_rf_prem*rf_prem(-1)+rf_prem_shk ;
	 lc_cor_tot = lc_cor+lc0_cor+lc995_cor ;
	 lc_cor = a_lc*lc_cor(-1)-lc_shk+a_lc2*(lyrow_gap-a_lc*lyrow_gap(-1)) ;
	 lc0_cor = lc0_shk ;
	 lc995_cor = 0.995*lc995_cor(-1)+lc995_shk ;
	 linv_cor = linv_cor1+linv_add ;
	 linv_add = chi3*lmc_r ;
	 linv_cor1 = a_linv*linv_cor1(-1)+linv_shk ;
	 gn_yn_cor = a_lg*gn_yn_cor(-1)+gn_yn_shk ;
	 gn_yn_ss = 1*gn_yn_ss(-1)+gn_yn_ss_shk ;
	 lx_cor = 0.89*lx_cor(-1)+lx_shk ;
	 lm_cor_tot = lm_cor+lm0_cor ;
	 lm_cor = a_lm*lm_cor(-1)+lm_shk ;
	 lm0_cor = lm0_cor(-1)+lm0_shk ;
	 lforexn_cor_tot = lforexn_cor+lforexn0_cor ;
	 lforexn_cor = a_lforexn*lforexn_cor(-1)+lforexn_shk ;
	 lforexn0_cor = lforexn0_shk ;
	 cbal_cor = a_cbal*cbal_cor(-1)+cbal_shk ;
	 nfa_ss = 1.8*nfa_ss(-1)-0.8*nfa_ss(-2)+nfa_ss_shk ;
	 risk_cor = 1.925*risk_cor(-1)-0.925*risk_cor(-2)+risk_shk ;
	 gbn_cor = a_gbn*gbn_cor(-1)+gbn_shk ;
	 gbn_yn_ss = 1*gbn_yn_ss(-1)+gbn_yn_ss_shk ;
	 rgbn_cor = a_rgbn*rgbn_cor(-1)+rgbn_shk ;
	 transf_r_cor = a_transf_r*transf_r_cor(-1)+transf_r_shk ;
	 transf_yn_ss = transf_yn_ss(-1)+transf_r_ss_shk ;
	 tdn_cor = a_tdn*tdn_cor(-1)+tdn_shk ;
	 tinc_cor = a_tinc*tinc_cor(-1)+tinc_shk ;
	 tincr_yn_ss = tincr_yn_ss(-1)+tinc_ss_shk ;
	 lk_cor = a_k*lk_cor(-1)+lk_shk ;
	 deprk = deprk(-1)+deprk_shk ;
	 usact = usact(-1)+usact_shk ;
	
	//%%%%%%%%%%%%%%%%%
	//% Foreign Block %
	//%%%%%%%%%%%%%%%%%
	lpcomrow_prow-log(pcomrow_prow_ss) = a_lpcomrow*(lpcomrow_prow(-1)-log(pcomrow_prow_ss))+lpcomrow_shk+exog_lpcomrow_prow1_shk ;
	log(pcomrow_prow_ss) = 1.29832*log(pcomrow_prow_ss(-1))-0.29832*log(pcomrow_prow_ss(-2))+lpcomrow_prow_shk ;
	lyrow = lyrow_sreq+lyrow_gap ;
	lyrow_sreq = ladotrow+1.2*lyrow_sreq(-1)-0.2*lyrow_sreq(-2)+lyrow_sreq_shk-0.006*sreqswitch*(lpcomrow_prow-lpcomrow_prow(-1)) ;
	ladotrow = ladot ;
	lyrow_gap = habitrow/(1+habitrow)*(lyrow_gap(-1)-(lyrow_sreq-lyrow_sreq(-1)))+1/(1+habitrow)*lyrow_gap(1)-(1-habitrow)/((1+habitrow)*subrow)*(r1row-r1row_ss)+lyrow_cor ;
	r1rown_des = (1+r1row_ss)^4*exp(infrow_ss)^4-1+0.7*lyrow_gap+4*(4*infrow(4)-4*infrow_ss(4)) ;
	dr1rown = (1-0.66)*(r1rown_des-r1rown_des(-1))+(1-0.66)*(1-0.67)*(r1rown_des(-1)-((1+r1rown(-1))^4-1))+0.66*0.67*dr1rown(-1)+r1rown_shk ;
	(1+r1rown)^4-1 = (1+r1rown(-1))^4-1+dr1rown ;
	r1rown = (1+r1row)*exp(infrow(1))-1 ;
	infq_row = infprow ;
	4*infprow = (1-0.644)*4*infprow(-1)+0.644*4*infrow+0.114*lyrow_gap+0.021*(lpcomrow_prow-lpcomrow_prow(-1)) ;
	infrow = dyn_pcrow/(1+beta*dyn_pcrow)*infrow(-1)+beta/(1+beta*dyn_pcrow)*infrow(1)+(1-dyn_pcrow)*(1-beta)/(1+beta*dyn_pcrow)*infrow_ss+(1-calvo_pcrow)*(1-beta*calvo_pcrow)/((1+beta*dyn_pcrow)*calvo_pcrow)*0.06*(lyrow_gap+0*(pcomrowshare*lpcomrow_prow+(1-pcomrowshare)*0-lpcrow_r))+lprow_cor ;
	infrow2 = dyn_pcrow/(1+beta*dyn_pcrow)*infrow2(-1)+beta/(1+beta*dyn_pcrow)*infrow2(1)+(1-dyn_pcrow)*(1-beta)/(1+beta*dyn_pcrow)*infrow_ss+(1-calvo_pcrow)*(1-beta*calvo_pcrow)/((1+beta*dyn_pcrow)*calvo_pcrow)*0.06*lyrow_gap+lprow2_cor ;
	lpcrow_r = lpcrow_r(-1)-infq_row+infrow ;
	lprow2_r-lprow2_r(-1) = infrow2-infq_row ;
	lpxrow_prow = 0.85*lpxrow_prow(-1)+lpxrow_prow_shk ;
	lpxrow_r = lpxrow_prow ;
	lyrow_cor = a_lyrow*lyrow_cor(-1)+lyrow_shk ;
	lprow_cor = a_lprow*lprow_cor(-1)+lprow_shk ;
	lprow2_cor = a_lprow*lprow2_cor(-1)+lprow2_shk ;
	r1row_ss = r1row_ss(-1)+r1_ss_shk ;
	infrow_ss = infrow_ss(-1)+infrow_shk ;
	end;

	%resid;
	options_.solve_tolf=1e-5;
	steady(maxit=1000);
	shocks;
	var lc_shk; stderr 0.0125*100;
	var lforexn_shk; stderr 0.012*100;
	var lpcomrow_shk; stderr 0.6*100;
	var la_shk; stderr 0.01*100;
	var lyrow_shk; stderr 0.002*100;
	var lxdc_shk; stderr 0.002*100;
	var gn_yn_shk; stderr 0.02*100;
	end;
	options_.Schur_vec_tol = 1e-6;
	
	stoch_simul (irf = 0,AR=100, noprint, nograph);