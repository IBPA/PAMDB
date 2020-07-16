function real_benchmark_fitting
	%sim_x_opt = simultaneous_fitting(1,1000,50000);
	%print_plotting_sim_data(sim_x_opt);
	sim_LOOCV(10,100,5000);
end
function output = simf_doi_10_1128_AEM_00791_07_figure3a(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pBAD = x(1);
	ALPHA_RBS_pC = x(2);
	K_arabinose = x(3);
	K_pBAD = x(4);
	alpha_pBAD = x(5);
	alpha_pC = x(6);
	beta_pBAD = x(7);
	n_arabinose = x(8);
	n_pBAD = x(9);
	scale_GFPuv_doi_10_1128_AEM_00791_07 = x(10);
	for i = 1:length(input)
		pC = 20*alpha_pC;
		arabinose = input(i);
		araC = ALPHA_RBS_pC*pC*K_arabinose^n_arabinose/(K_arabinose^n_arabinose + arabinose^n_arabinose);
		pBAD = 20*(beta_pBAD + (alpha_pBAD - beta_pBAD)*K_pBAD^n_pBAD/(K_pBAD^n_pBAD + araC^n_pBAD));
		GFPuv = scale_GFPuv_doi_10_1128_AEM_00791_07*ALPHA_RBS_pBAD*pBAD;
		output(i) = GFPuv;
	end
end

function output = simf_doi_10_1038_nature12148_figure16a(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0030 = x(11);
	K_arabinose = x(3);
	K_pBAD = x(4);
	alpha_pBAD = x(5);
	alpha_pLlacO_1 = x(12);
	beta_pBAD = x(7);
	n_arabinose = x(8);
	n_pBAD = x(9);
	scale_mCherry_doi_10_1038_nature12148 = x(13);
	for i = 1:length(input)
		pLlacO_1 = 5*alpha_pLlacO_1;
		arabinose = input(i);
		araC = ALPHA_BBa_B0030*pLlacO_1*K_arabinose^n_arabinose/(K_arabinose^n_arabinose + arabinose^n_arabinose);
		pBAD = 15*(beta_pBAD + (alpha_pBAD - beta_pBAD)*K_pBAD^n_pBAD/(K_pBAD^n_pBAD + araC^n_pBAD));
		mCherry = scale_mCherry_doi_10_1038_nature12148*ALPHA_BBa_B0030*pBAD;
		output(i) = mCherry;
	end
end

function output = simf_doi__10_1093_nar_25_6_1203_figure4a(x, input)
	output = zeros(1,length(input));
	ALPHA_RBSII = x(14);
	ALPHA_RBS_pN25 = x(15);
	KE1_pLtetO_1 = x(16);
	KE2_pLtetO_1 = x(17);
	K_aTc = x(18);
	alpha_pLtetO_1 = x(19);
	alpha_pN25 = x(20);
	beta_pLtetO_1 = x(21);
	n_aTc = x(22);
	scale_luc_doi__10_1093_nar_25_6_1203 = x(23);
	for i = 1:length(input)
		pN25 = 1*alpha_pN25;
		aTc = input(i);
		tetR = ALPHA_RBS_pN25*pN25*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pLtetO_1 = 15*(beta_pLtetO_1 + (alpha_pLtetO_1 - beta_pLtetO_1)/(1 + tetR^4/KE1_pLtetO_1 + tetR^2/KE2_pLtetO_1));
		luc = scale_luc_doi__10_1093_nar_25_6_1203*ALPHA_RBSII*pLtetO_1;
		output(i) = luc;
	end
end

function output = simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6(x, input)
	output = zeros(1,length(input));
	ALPHA_RBSII = x(14);
	ALPHA_RBS_placI = x(24);
	KE1_pLAC = x(25);
	KE2_pLAC = x(26);
	K_IPTG = x(27);
	alpha_pLAC = x(28);
	alpha_placIq = x(29);
	beta_pLAC = x(30);
	n_IPTG = x(31);
	scale_eYFP_http___dspace_mit_edu_handle_1721_1_8228 = x(32);
	for i = 1:length(input)
		placIq = 10*alpha_placIq;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placIq*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		lacI_free = ALPHA_RBS_placI*placIq;
		pLAC = 10*(beta_pLAC + (alpha_pLAC - beta_pLAC)/(1 + lacI^4/KE1_pLAC + (lacI_free^4 - lacI^4)/KE2_pLAC));
		eYFP = scale_eYFP_http___dspace_mit_edu_handle_1721_1_8228*ALPHA_RBSII*pLAC;
		output(i) = eYFP;
	end
end

function output = simf_doi_10_1073pnas_0408507102_figure2A(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_Unknown_Hooshangi = x(33);
	ALPHA_RBS_placI = x(24);
	KE1_pLtetO_1 = x(16);
	KE2_pLtetO_1 = x(17);
	K_aTc = x(18);
	alpha_pLtetO_1 = x(19);
	alpha_placIq = x(29);
	beta_pLtetO_1 = x(21);
	n_aTc = x(22);
	scale_eYFP_doi_10_1073pnas_0408507102 = x(34);
	for i = 1:length(input)
		placIq = 10*alpha_placIq;
		aTc = input(i);
		tetR = ALPHA_RBS_placI*placIq*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pLtetO_1 = 10*(beta_pLtetO_1 + (alpha_pLtetO_1 - beta_pLtetO_1)/(1 + tetR^4/KE1_pLtetO_1 + tetR^2/KE2_pLtetO_1));
		eYFP = scale_eYFP_doi_10_1073pnas_0408507102*ALPHA_RBS_Unknown_Hooshangi*pLtetO_1;
		output(i) = eYFP;
	end
end

function output = simf_doi_10_1073_pnas_0606717104_figure4_6(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pLAC = x(35);
	ALPHA_RBS_placI = x(24);
	KE1_pLAC = x(25);
	KE2_pLAC = x(26);
	K_IPTG = x(27);
	alpha_pLAC = x(28);
	alpha_placI = x(36);
	beta_pLAC = x(30);
	n_IPTG = x(31);
	scale_lacZ_doi_10_1073_pnas_0606717104 = x(37);
	for i = 1:length(input)
		placI = 1*alpha_placI;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placI*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		lacI_free = ALPHA_RBS_placI*placI;
		pLAC = 1*(beta_pLAC + (alpha_pLAC - beta_pLAC)/(1 + lacI^4/KE1_pLAC + (lacI_free^4 - lacI^4)/KE2_pLAC));
		lacZ = scale_lacZ_doi_10_1073_pnas_0606717104*ALPHA_RBS_pLAC*pLAC;
		output(i) = lacZ;
	end
end

function output = simf_doi_10_1186_1754_1611_3_4_figure3_1(x)
	ALPHA_BBa_B0032 = x(38);
	alpha_BBa_J23101 = x(39);
	scale_GFPmut3b_doi_10_1186_1754_1611_3_4 = x(40);
	BBa_J23101 = 10*alpha_BBa_J23101;
	GFPmut3b = scale_GFPmut3b_doi_10_1186_1754_1611_3_4*ALPHA_BBa_B0032*BBa_J23101;
	output = GFPmut3b;
end

function output = simf_doi_10_1186_1754_1611_3_4_figure3_2(x)
	ALPHA_BBa_B0032 = x(38);
	alpha_BBa_J23116 = x(41);
	scale_GFPmut3b_doi_10_1186_1754_1611_3_4 = x(40);
	BBa_J23116 = 10*alpha_BBa_J23116;
	GFPmut3b = scale_GFPmut3b_doi_10_1186_1754_1611_3_4*ALPHA_BBa_B0032*BBa_J23116;
	output = GFPmut3b;
end

function output = simf_doi_10_1186_1754_1611_3_4_figure3_3(x)
	ALPHA_BBa_B0032 = x(38);
	alpha_BBa_J23150 = x(42);
	scale_GFPmut3b_doi_10_1186_1754_1611_3_4 = x(40);
	BBa_J23150 = 10*alpha_BBa_J23150;
	GFPmut3b = scale_GFPmut3b_doi_10_1186_1754_1611_3_4*ALPHA_BBa_B0032*BBa_J23150;
	output = GFPmut3b;
end

function output = simf_doi_10_1186_1754_1611_3_4_figure3_4(x)
	ALPHA_BBa_B0032 = x(38);
	alpha_BBa_J23151 = x(43);
	scale_GFPmut3b_doi_10_1186_1754_1611_3_4 = x(40);
	BBa_J23151 = 10*alpha_BBa_J23151;
	GFPmut3b = scale_GFPmut3b_doi_10_1186_1754_1611_3_4*ALPHA_BBa_B0032*BBa_J23151;
	output = GFPmut3b;
end

function output = simf_doi_10_1186_1754_1611_3_4_figure3_5(x)
	ALPHA_BBa_B0032 = x(38);
	alpha_BBa_J23102 = x(44);
	scale_GFPmut3b_doi_10_1186_1754_1611_3_4 = x(40);
	BBa_J23102 = 10*alpha_BBa_J23102;
	GFPmut3b = scale_GFPmut3b_doi_10_1186_1754_1611_3_4*ALPHA_BBa_B0032*BBa_J23102;
	output = GFPmut3b;
end

function output = simf_doi_10_1186_1754_1611_3_4_figure3_6(x)
	ALPHA_BBa_B0032 = x(38);
	alpha_pLtetO_1 = x(19);
	scale_GFPmut3b_doi_10_1186_1754_1611_3_4 = x(40);
	pLtetO_1 = 10*alpha_pLtetO_1;
	GFPmut3b = scale_GFPmut3b_doi_10_1186_1754_1611_3_4*ALPHA_BBa_B0032*pLtetO_1;
	output = GFPmut3b;
end

function output = simf_doi_10_1186_1754_1611_3_4_figure3_7(x)
	ALPHA_BBa_B0032 = x(38);
	alpha_pLlacO_1 = x(12);
	scale_GFPmut3b_doi_10_1186_1754_1611_3_4 = x(40);
	pLlacO_1 = 10*alpha_pLlacO_1;
	GFPmut3b = scale_GFPmut3b_doi_10_1186_1754_1611_3_4*ALPHA_BBa_B0032*pLlacO_1;
	output = GFPmut3b;
end

function output = simf_doi_10_1093_nar_gku593_figure3D_1(x)
	ALPHA_BBa_B0030 = x(11);
	alpha_BBa_J23109 = x(45);
	scale_GFPmut3b_doi_10_1093_nar_gku593 = x(46);
	BBa_J23109 = 10*alpha_BBa_J23109;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku593*ALPHA_BBa_B0030*BBa_J23109;
	output = GFPmut3b;
end

function output = simf_doi_10_1093_nar_gku593_figure3D_2(x)
	ALPHA_BBa_B0030 = x(11);
	alpha_BBa_J23114 = x(47);
	scale_GFPmut3b_doi_10_1093_nar_gku593 = x(46);
	BBa_J23114 = 10*alpha_BBa_J23114;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku593*ALPHA_BBa_B0030*BBa_J23114;
	output = GFPmut3b;
end

function output = simf_doi_10_1093_nar_gku593_figure3D_3(x)
	ALPHA_BBa_B0030 = x(11);
	alpha_BBa_J23116 = x(41);
	scale_GFPmut3b_doi_10_1093_nar_gku593 = x(46);
	BBa_J23116 = 10*alpha_BBa_J23116;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku593*ALPHA_BBa_B0030*BBa_J23116;
	output = GFPmut3b;
end

function output = simf_doi_10_1093_nar_gku593_figure3D_4(x)
	ALPHA_BBa_B0030 = x(11);
	alpha_BBa_J23115 = x(48);
	scale_GFPmut3b_doi_10_1093_nar_gku593 = x(46);
	BBa_J23115 = 10*alpha_BBa_J23115;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku593*ALPHA_BBa_B0030*BBa_J23115;
	output = GFPmut3b;
end

function output = simf_doi_10_1093_nar_gku593_figure3D_5(x)
	ALPHA_BBa_B0030 = x(11);
	alpha_BBa_J23105 = x(49);
	scale_GFPmut3b_doi_10_1093_nar_gku593 = x(46);
	BBa_J23105 = 10*alpha_BBa_J23105;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku593*ALPHA_BBa_B0030*BBa_J23105;
	output = GFPmut3b;
end

function output = simf_doi_10_1093_nar_gku593_figure3D_6(x)
	ALPHA_BBa_B0030 = x(11);
	alpha_BBa_J23106 = x(50);
	scale_GFPmut3b_doi_10_1093_nar_gku593 = x(46);
	BBa_J23106 = 10*alpha_BBa_J23106;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku593*ALPHA_BBa_B0030*BBa_J23106;
	output = GFPmut3b;
end

function output = simf_doi_10_1093_nar_gku593_figureS6(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0030 = x(11);
	ALPHA_RBS_pC = x(2);
	K_arabinose = x(3);
	K_pBAD = x(4);
	alpha_pBAD = x(5);
	alpha_pC = x(6);
	beta_pBAD = x(7);
	n_arabinose = x(8);
	n_pBAD = x(9);
	scale_GFPmut3b_doi_10_1093_nar_gku593 = x(46);
	for i = 1:length(input)
		pC = 10*alpha_pC;
		arabinose = input(i);
		araC = ALPHA_RBS_pC*pC*K_arabinose^n_arabinose/(K_arabinose^n_arabinose + arabinose^n_arabinose);
		pBAD = 10*(beta_pBAD + (alpha_pBAD - beta_pBAD)*K_pBAD^n_pBAD/(K_pBAD^n_pBAD + araC^n_pBAD));
		GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku593*ALPHA_BBa_B0030*pBAD;
		output(i) = GFPmut3b;
	end
end

function output = simf_doi_10_1093_nar_gku1388_figure2A_1(x)
	ALPHA_BBa_B0030 = x(11);
	alpha_BBa_J23101 = x(39);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(51);
	BBa_J23101 = 10*alpha_BBa_J23101;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*BBa_J23101;
	output = GFPmut3b;
end

function output = simf_doi_10_1093_nar_gku1388_figure2A_2(x)
	ALPHA_BBa_B0030 = x(11);
	alpha_BBa_J23106 = x(50);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(51);
	BBa_J23106 = 10*alpha_BBa_J23106;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*BBa_J23106;
	output = GFPmut3b;
end

function output = simf_doi_10_1093_nar_gku1388_figure2A_3(x)
	ALPHA_BBa_B0030 = x(11);
	alpha_BBa_J23105 = x(49);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(51);
	BBa_J23105 = 10*alpha_BBa_J23105;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*BBa_J23105;
	output = GFPmut3b;
end

function output = simf_doi_10_1093_nar_gku1388_figure2A_4(x)
	ALPHA_BBa_B0030 = x(11);
	alpha_BBa_J23115 = x(48);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(51);
	BBa_J23115 = 10*alpha_BBa_J23115;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*BBa_J23115;
	output = GFPmut3b;
end

function output = simf_doi_10_1093_nar_gku1388_figure2A_5(x)
	ALPHA_BBa_B0030 = x(11);
	alpha_BBa_J23114 = x(47);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(51);
	BBa_J23114 = 10*alpha_BBa_J23114;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*BBa_J23114;
	output = GFPmut3b;
end

function output = simf_doi_10_1093_nar_gku1388_figure2A_6(x)
	ALPHA_BBa_B0030 = x(11);
	alpha_BBa_J23117 = x(52);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(51);
	BBa_J23117 = 10*alpha_BBa_J23117;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*BBa_J23117;
	output = GFPmut3b;
end

function output = simf_doi_10_1093_nar_gku1388_figure2B_1(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0030 = x(11);
	KE1_pTET_star_ = x(53);
	KE2_pTET_star_ = x(54);
	K_aTc = x(18);
	alpha_BBa_J23117 = x(52);
	alpha_pTET_star_ = x(55);
	beta_pTET_star_ = x(56);
	n_aTc = x(22);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(51);
	for i = 1:length(input)
		BBa_J23117 = 10*alpha_BBa_J23117;
		aTc = input(i);
		tetR = ALPHA_BBa_B0030*BBa_J23117*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pTET_star_ = 10*(beta_pTET_star_ + (alpha_pTET_star_ - beta_pTET_star_)/(1 + tetR^4/KE1_pTET_star_ + tetR^2/KE2_pTET_star_));
		GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*pTET_star_;
		output(i) = GFPmut3b;
	end
end

function output = simf_doi_10_1093_nar_gku1388_figure2B_2(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0030 = x(11);
	KE1_pTET_star_ = x(53);
	KE2_pTET_star_ = x(54);
	K_aTc = x(18);
	alpha_BBa_J23114 = x(47);
	alpha_pTET_star_ = x(55);
	beta_pTET_star_ = x(56);
	n_aTc = x(22);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(51);
	for i = 1:length(input)
		BBa_J23114 = 10*alpha_BBa_J23114;
		aTc = input(i);
		tetR = ALPHA_BBa_B0030*BBa_J23114*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pTET_star_ = 10*(beta_pTET_star_ + (alpha_pTET_star_ - beta_pTET_star_)/(1 + tetR^4/KE1_pTET_star_ + tetR^2/KE2_pTET_star_));
		GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*pTET_star_;
		output(i) = GFPmut3b;
	end
end

function output = simf_doi_10_1038_nature09565_figure1C(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0033 = x(57);
	ALPHA_BBa_B0034 = x(58);
	K_arabinose = x(3);
	K_pBAD = x(4);
	alpha_BBa_J23117 = x(52);
	alpha_pBAD = x(5);
	beta_pBAD = x(7);
	n_arabinose = x(8);
	n_pBAD = x(9);
	scale_eYFP_doi_10_1038_nature09565 = x(59);
	for i = 1:length(input)
		BBa_J23117 = 15*alpha_BBa_J23117;
		arabinose = input(i);
		araC = ALPHA_BBa_B0034*BBa_J23117*K_arabinose^n_arabinose/(K_arabinose^n_arabinose + arabinose^n_arabinose);
		pBAD = 15*(beta_pBAD + (alpha_pBAD - beta_pBAD)*K_pBAD^n_pBAD/(K_pBAD^n_pBAD + araC^n_pBAD));
		eYFP = scale_eYFP_doi_10_1038_nature09565*ALPHA_BBa_B0033*pBAD;
		output(i) = eYFP;
	end
end

function output = simf_doi_10_1038_nature09565_figureS1C(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0032 = x(38);
	ALPHA_BBa_B0033 = x(57);
	KE1_pLtetO_1 = x(16);
	KE2_pLtetO_1 = x(17);
	K_aTc = x(18);
	alpha_BBa_J23117 = x(52);
	alpha_pLtetO_1 = x(19);
	beta_pLtetO_1 = x(21);
	n_aTc = x(22);
	scale_eYFP_doi_10_1038_nature09565 = x(59);
	for i = 1:length(input)
		BBa_J23117 = 15*alpha_BBa_J23117;
		aTc = input(i);
		tetR = ALPHA_BBa_B0032*BBa_J23117*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pLtetO_1 = 15*(beta_pLtetO_1 + (alpha_pLtetO_1 - beta_pLtetO_1)/(1 + tetR^4/KE1_pLtetO_1 + tetR^2/KE2_pLtetO_1));
		eYFP = scale_eYFP_doi_10_1038_nature09565*ALPHA_BBa_B0033*pLtetO_1;
		output(i) = eYFP;
	end
end

function output = simf_doi_10_1038_nbt_2401_figureS11(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_A = x(60);
	KE1_pTAC = x(61);
	KE2_pTAC = x(62);
	K_IPTG = x(27);
	alpha_BBa_J23101 = x(39);
	alpha_pTAC = x(63);
	beta_pTAC = x(64);
	n_IPTG = x(31);
	scale_sfGFP_doi_10_1038_nbt_2401 = x(65);
	for i = 1:length(input)
		BBa_J23101 = 5*alpha_BBa_J23101;
		IPTG = input(i);
		lacI = ALPHA_RBS_A*BBa_J23101*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		lacI_free = ALPHA_RBS_A*BBa_J23101;
		pTAC = 5*(beta_pTAC + (alpha_pTAC - beta_pTAC)/(1 + lacI^4/KE1_pTAC + (lacI_free^4 - lacI^4)/KE2_pTAC));
		sfGFP = scale_sfGFP_doi_10_1038_nbt_2401*ALPHA_RBS_A*pTAC;
		output(i) = sfGFP;
	end
end

function output = simf_doi_10_1038_nbt_2401_figureS13(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_A = x(60);
	K_arabinose = x(3);
	K_pBAD = x(4);
	alpha_BBa_J23101 = x(39);
	alpha_pBAD = x(5);
	beta_pBAD = x(7);
	n_arabinose = x(8);
	n_pBAD = x(9);
	scale_sfGFP_doi_10_1038_nbt_2401 = x(65);
	for i = 1:length(input)
		BBa_J23101 = 5*alpha_BBa_J23101;
		arabinose = input(i);
		araC = ALPHA_RBS_A*BBa_J23101*K_arabinose^n_arabinose/(K_arabinose^n_arabinose + arabinose^n_arabinose);
		pBAD = 5*(beta_pBAD + (alpha_pBAD - beta_pBAD)*K_pBAD^n_pBAD/(K_pBAD^n_pBAD + araC^n_pBAD));
		sfGFP = scale_sfGFP_doi_10_1038_nbt_2401*ALPHA_RBS_A*pBAD;
		output(i) = sfGFP;
	end
end

function output = simf_doi_10_1038_nature11516_figureS8A(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pC = x(2);
	ALPHA_RBS_psicA = x(66);
	K_arabinose = x(3);
	K_pBAD = x(4);
	alpha_pBAD = x(5);
	alpha_pC = x(6);
	beta_pBAD = x(7);
	n_arabinose = x(8);
	n_pBAD = x(9);
	scale_mRFP1_doi_10_1038_nature11516 = x(67);
	for i = 1:length(input)
		pC = 15*alpha_pC;
		arabinose = input(i);
		araC = ALPHA_RBS_pC*pC*K_arabinose^n_arabinose/(K_arabinose^n_arabinose + arabinose^n_arabinose);
		pBAD = 15*(beta_pBAD + (alpha_pBAD - beta_pBAD)*K_pBAD^n_pBAD/(K_pBAD^n_pBAD + araC^n_pBAD));
		mRFP1 = scale_mRFP1_doi_10_1038_nature11516*ALPHA_RBS_psicA*pBAD;
		output(i) = mRFP1;
	end
end

function output = simf_doi_10_1038_nature11516_figureS8B(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_placI = x(24);
	ALPHA_RBS_psicA = x(66);
	KE1_pTAC = x(61);
	KE2_pTAC = x(62);
	K_IPTG = x(27);
	alpha_pTAC = x(63);
	alpha_placIq = x(29);
	beta_pTAC = x(64);
	n_IPTG = x(31);
	scale_mRFP1_doi_10_1038_nature11516 = x(67);
	for i = 1:length(input)
		placIq = 15*alpha_placIq;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placIq*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		lacI_free = ALPHA_RBS_placI*placIq;
		pTAC = 15*(beta_pTAC + (alpha_pTAC - beta_pTAC)/(1 + lacI^4/KE1_pTAC + (lacI_free^4 - lacI^4)/KE2_pTAC));
		mRFP1 = scale_mRFP1_doi_10_1038_nature11516*ALPHA_RBS_psicA*pTAC;
		output(i) = mRFP1;
	end
end

function output = simf_doi_10_1038_nbt_2149_Supplement_page17(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pET_29b = x(68);
	ALPHA_RBS_placI = x(24);
	KE1_placUV5 = x(69);
	KE2_placUV5 = x(70);
	K_IPTG = x(27);
	alpha_placI = x(36);
	alpha_placUV5 = x(71);
	beta_placUV5 = x(72);
	n_IPTG = x(31);
	scale_mRFP1_doi_10_1038_nbt_2149 = x(73);
	for i = 1:length(input)
		placI = 10*alpha_placI;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placI*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		lacI_free = ALPHA_RBS_placI*placI;
		placUV5 = 10*(beta_placUV5 + (alpha_placUV5 - beta_placUV5)/(1 + lacI^4/KE1_placUV5 + (lacI_free^4 - lacI^4)/KE2_placUV5));
		mRFP1 = scale_mRFP1_doi_10_1038_nbt_2149*ALPHA_RBS_pET_29b*placUV5;
		output(i) = mRFP1;
	end
end

function x_opt = simultaneous_fitting(ensemble_size, ite_num, func_ite_num)
	parameter_name = {'ALPHA_RBS_pBAD'; 'ALPHA_RBS_pC'; 'K_arabinose'; 'K_pBAD'; 'alpha_pBAD'; 'alpha_pC'; 'beta_pBAD'; 'n_arabinose'; 'n_pBAD'; 'scale_GFPuv_doi_10_1128_'; 'ALPHA_BBa_B0030'; 'alpha_pLlacO_1'; 'scale_mCherry_doi_10_103'; 'ALPHA_RBSII'; 'ALPHA_RBS_pN25'; 'KE1_pLtetO_1'; 'KE2_pLtetO_1'; 'K_aTc'; 'alpha_pLtetO_1'; 'alpha_pN25'; 'beta_pLtetO_1'; 'n_aTc'; 'scale_luc_doi__10_1093_n'; 'ALPHA_RBS_placI'; 'KE1_pLAC'; 'KE2_pLAC'; 'K_IPTG'; 'alpha_pLAC'; 'alpha_placIq'; 'beta_pLAC'; 'n_IPTG'; 'scale_eYFP_http___dspace'; 'ALPHA_RBS_Unknown_Hoosha'; 'scale_eYFP_doi_10_1073pn'; 'ALPHA_RBS_pLAC'; 'alpha_placI'; 'scale_lacZ_doi_10_1073_p'; 'ALPHA_BBa_B0032'; 'alpha_BBa_J23101'; 'scale_GFPmut3b_doi_10_11'; 'alpha_BBa_J23116'; 'alpha_BBa_J23150'; 'alpha_BBa_J23151'; 'alpha_BBa_J23102'; 'alpha_BBa_J23109'; 'scale_GFPmut3b_doi_10_10'; 'alpha_BBa_J23114'; 'alpha_BBa_J23115'; 'alpha_BBa_J23105'; 'alpha_BBa_J23106'; 'scale_GFPmut3b_doi_10_10'; 'alpha_BBa_J23117'; 'KE1_pTET_star_'; 'KE2_pTET_star_'; 'alpha_pTET_star_'; 'beta_pTET_star_'; 'ALPHA_BBa_B0033'; 'ALPHA_BBa_B0034'; 'scale_eYFP_doi_10_1038_n'; 'ALPHA_RBS_A'; 'KE1_pTAC'; 'KE2_pTAC'; 'alpha_pTAC'; 'beta_pTAC'; 'scale_sfGFP_doi_10_1038_'; 'ALPHA_RBS_psicA'; 'scale_mRFP1_doi_10_1038_'; 'ALPHA_RBS_pET_29b'; 'KE1_placUV5'; 'KE2_placUV5'; 'alpha_placUV5'; 'beta_placUV5'; 'scale_mRFP1_doi_10_1038_'};
	lb = [0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01 0.03 0.002 0.01 0.03 0.03 0.001 0.001 0.001 0.002 0.002 0 1 0.01 0.03 0.001 0.001 0.001 0.002 0.002 0 1 0.01 0.03 0.01 0.03 0.002 0.01 0.03 0.002 0.01 0.002 0.002 0.002 0.002 0.002 0.01 0.002 0.002 0.002 0.002 0.01 0.002 0.001 0.001 0.002 0 0.03 0.03 0.01 0.03 0.001 0.001 0.002 0 0.01 0.03 0.01 0.03 0.001 0.001 0.002 0 0.01];
	ub = [5 5 50 50 50 50 1 4 4 10 5 50 10 5 5 50 50 50 50 50 1 4 10 5 50 50 50 50 50 1 4 10 5 10 5 50 10 5 50 10 50 50 50 50 50 10 50 50 50 50 10 50 50 50 50 1 5 5 10 5 50 50 50 1 10 5 10 5 50 50 50 1 10];
	[CI_lb, CI_ub, x_opt, solution_ensemble] = fit_a_model(ensemble_size, ite_num, func_ite_num, @simultaneous_fitting_error, lb, ub);
	parameter_fileID = fopen('Plot_data/sim_CI_Ellis.dat','w');
	for i = 1:length(parameter_name)
		fprintf(parameter_fileID, '%s\t', parameter_name{i});
		fprintf(parameter_fileID, '%f\t', lb(i));
		fprintf(parameter_fileID, '%f\t', ub(i));
		fprintf(parameter_fileID, '%f\t', CI_lb(i));
		fprintf(parameter_fileID, '%f\t', CI_ub(i));
		fprintf(parameter_fileID, '%f\n', x_opt(i));
	end
	fclose(parameter_fileID);
	print_parameter_to_json('JSON/inferred_parameter_value.json', lb, ub, CI_lb, CI_ub, x_opt, solution_ensemble, parameter_name);
	SA_fileID = fopen('Plot_data/sensitivity_analysis_Ellis.dat','w');
	opt_error = norm(simultaneous_fitting_error(x_opt));
	pertubation = [0.01 0.05 0.1 0.5 2 10 50 100];
	error_diff = pertubation;
	for i = 1:length(parameter_name)
		x = x_opt;
		for j = 1:length(pertubation)
			x(i) = pertubation(j)*x_opt(i);
			error = norm(simultaneous_fitting_error(x));
			error_diff(j) = abs(error - opt_error)/opt_error;
		end
		fprintf(SA_fileID, '%s	', parameter_name{i});
		fprintf(parameter_fileID, '%f\t', mean(error_diff));
		fprintf(parameter_fileID, '%f\n', std(error_diff)/sqrt(length(error_diff)));
	end
	fclose(SA_fileID);
end
function error = simultaneous_fitting_error(x)
	error = zeros(169,1);
	current_base_index = 0;
	%==============
	arabinose = [0.003 0.005 0.05 0.5 1];
	GFPuv = [0.044 0.089 0.556 1 1];
	output = 1*(simf_doi_10_1128_AEM_00791_07_figure3a(x, arabinose) - GFPuv);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.003 0.008 0.025 0.075 0.2 0.7 2 6 20 50];
	mCherry = [0.015 0.015 0.015 0.02 0.06 0.48 0.84 0.93 0.96 0.98 0.98 1];
	output = 1*(simf_doi_10_1038_nature12148_figure16a(x, arabinose) - mCherry);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 1 2 3 4 5 6 7 8 9 10 12 13 17 18 20 25 35 60];
	luc = [0 0.001 0.001 0.001 0.002 0.004 0.007 0.011 0.014 0.032 0.036 0.054 0.089 0.643 0.893 1 1 1 1];
	output = 1*(simf_doi__10_1093_nar_25_6_1203_figure4a(x, aTc) - luc);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0.001 0.01 0.03 0.05 0.075 0.1 0.15 0.2 0.3 0.5 1];
	eYFP = [0.015 0.014 0.023 0.062 0.098 0.187 0.239 0.374 0.572 0.869 1];
	output = 1*(simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6(x, IPTG) - eYFP);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [9.258 13.887 20.831 32.403 46.29 83.322 162.015 231.45 370.32 648.06 925.8 1388.7];
	eYFP = [0.003 0.003 0.003 0.003 0.005 0.01 0.05 0.167 0.333 0.667 0.833 1];
	output = 1*(simf_doi_10_1073pnas_0408507102_figure2A(x, aTc) - eYFP);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.005 0.01 0.02 0.05 0.1 0.2 1];
	lacZ = [0.001 0.003 0.043 0.5 0.984 0.999 1 1];
	output = 1*(simf_doi_10_1073_pnas_0606717104_figure4_6(x, IPTG) - lacZ);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) = 1.45*(simf_doi_10_1186_1754_1611_3_4_figure3_1(x) - 0.69);
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) = 55.769*(simf_doi_10_1186_1754_1611_3_4_figure3_2(x) - 0.018);
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) = 5.179*(simf_doi_10_1186_1754_1611_3_4_figure3_3(x) - 0.193);
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) = 2.5*(simf_doi_10_1186_1754_1611_3_4_figure3_4(x) - 0.4);
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) = 1.526*(simf_doi_10_1186_1754_1611_3_4_figure3_5(x) - 0.655);
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) = 1.115*(simf_doi_10_1186_1754_1611_3_4_figure3_6(x) - 0.897);
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) = 1*(simf_doi_10_1186_1754_1611_3_4_figure3_7(x) - 1);
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) = 560.2*(simf_doi_10_1093_nar_gku593_figure3D_1(x) - 0.002);
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) = 56.02*(simf_doi_10_1093_nar_gku593_figure3D_2(x) - 0.018);
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) = 28.01*(simf_doi_10_1093_nar_gku593_figure3D_3(x) - 0.036);
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) = 20.007*(simf_doi_10_1093_nar_gku593_figure3D_4(x) - 0.05);
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) = 11.671*(simf_doi_10_1093_nar_gku593_figure3D_5(x) - 0.086);
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) = 4.668*(simf_doi_10_1093_nar_gku593_figure3D_6(x) - 0.214);
	current_base_index = current_base_index + 1;
	%==============
	arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2];
	GFPmut3b = [0 0.002 0.041 0.479 0.949 0.997 1 1 1];
	output = 1*(simf_doi_10_1093_nar_gku593_figureS6(x, arabinose) - GFPmut3b);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) = 1*(simf_doi_10_1093_nar_gku1388_figure2A_1(x) - 1);
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) = 4.643*(simf_doi_10_1093_nar_gku1388_figure2A_2(x) - 0.215);
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) = 10.484*(simf_doi_10_1093_nar_gku1388_figure2A_3(x) - 0.095);
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) = 21.667*(simf_doi_10_1093_nar_gku1388_figure2A_4(x) - 0.046);
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) = 54.167*(simf_doi_10_1093_nar_gku1388_figure2A_5(x) - 0.018);
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) = 325*(simf_doi_10_1093_nar_gku1388_figure2A_6(x) - 0.003);
	current_base_index = current_base_index + 1;
	%==============
	aTc = [2 4 10 25 50 100 200];
	GFPmut3b = [0.002 0.014 0.167 0.666 0.817 0.841 0.831];
	output = 1.189*(simf_doi_10_1093_nar_gku1388_figure2B_1(x, aTc) - GFPmut3b);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [10 25 50 100 200];
	GFPmut3b = [0.001 0.022 0.176 0.411 0.461];
	output = 2.171*(simf_doi_10_1093_nar_gku1388_figure2B_2(x, aTc) - GFPmut3b);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.015 0.05 0.5 5 10];
	eYFP = [0.002 0.002 0.057 0.628 0.999 1 1];
	output = 1*(simf_doi_10_1038_nature09565_figure1C(x, arabinose) - eYFP);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 0.025 0.25 2.5 25 100 250];
	eYFP = [0.005 0.006 0.007 0.03 0.297 0.388 0.398];
	output = 2.515*(simf_doi_10_1038_nature09565_figureS1C(x, aTc) - eYFP);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1];
	sfGFP = [0.025 0.025 0.035 0.055 0.115 0.2 0.3 0.4 0.65 1];
	output = 1*(simf_doi_10_1038_nbt_2401_figureS11(x, IPTG) - sfGFP);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.1 1 2 5 7 10 12.5 25 37.5 50 62.5];
	sfGFP = [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425];
	output = 2.353*(simf_doi_10_1038_nbt_2401_figureS13(x, arabinose) - sfGFP);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0 0.002 0.01 0.04 0.16 1 4 25];
	mRFP1 = [0.003 0.003 0.003 0.017 0.23 0.82 0.939 0.941 0.941];
	output = 1.062*(simf_doi_10_1038_nature11516_figureS8A(x, arabinose) - mRFP1);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0 0 0.001 0.004 0.025 0.1 0.5 2.5];
	mRFP1 = [0.027 0.027 0.027 0.03 0.067 0.46 0.884 0.991 1];
	output = 1*(simf_doi_10_1038_nature11516_figureS8B(x, IPTG) - mRFP1);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.025 0.05 0.1 0.25 0.5];
	mRFP1 = [0.02 0.031 0.092 0.306 0.745 0.898 0.969 1 1];
	output = 1*(simf_doi_10_1038_nbt_2149_Supplement_page17(x, IPTG) - mRFP1);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
end

function print_plotting_sim_data(sim_x_opt)
	sim_pair_fileID = fopen('Plot_data/sim_pair_Ellis.dat','w');
	json_fileID = fopen('JSON/simulation_data_Ellis.json','w');
	%%%%%%%%%%%%%
	arabinose = [0.003 0.005 0.05 0.5 1];
	GFPuv = [0.044 0.089 0.556 1 1];
	output_sim = simf_doi_10_1128_AEM_00791_07_figure3a(sim_x_opt,arabinose);
	sim_pair_matrix = [output_sim; GFPuv];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1128/AEM.00791-07', 'figure3a', arabinose, output_sim, {'pC = 20*{alpha_pC}'; 'araC = {ALPHA_RBS_pC}*pC*{K_arabinose}^{n_arabinose}/({K_arabinose}^{n_arabinose} + arabinose^{n_arabinose})'; 'pBAD = 20*({beta_pBAD} + ({alpha_pBAD} - {beta_pBAD})*{K_pBAD}^{n_pBAD}/({K_pBAD}^{n_pBAD} + araC^{n_pBAD}))'; 'GFPuv = {scale_GFPuv_doi_10_1128_AEM_00791_07}*{ALPHA_RBS_pBAD}*pBAD'});
	%%%%%%%%%%%%%
	arabinose = [0 0.001 0.003 0.008 0.025 0.075 0.2 0.7 2 6 20 50];
	mCherry = [0.015 0.015 0.015 0.02 0.06 0.48 0.84 0.93 0.96 0.98 0.98 1];
	output_sim = simf_doi_10_1038_nature12148_figure16a(sim_x_opt,arabinose);
	sim_pair_matrix = [output_sim; mCherry];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1038/nature12148', 'figure16a', arabinose, output_sim, {'pLlacO_1 = 5*{alpha_pLlacO_1}'; 'araC = {ALPHA_BBa_B0030}*pLlacO_1*{K_arabinose}^{n_arabinose}/({K_arabinose}^{n_arabinose} + arabinose^{n_arabinose})'; 'pBAD = 15*({beta_pBAD} + ({alpha_pBAD} - {beta_pBAD})*{K_pBAD}^{n_pBAD}/({K_pBAD}^{n_pBAD} + araC^{n_pBAD}))'; 'mCherry = {scale_mCherry_doi_10_1038_nature12148}*{ALPHA_BBa_B0030}*pBAD'});
	%%%%%%%%%%%%%
	aTc = [0 1 2 3 4 5 6 7 8 9 10 12 13 17 18 20 25 35 60];
	luc = [0 0.001 0.001 0.001 0.002 0.004 0.007 0.011 0.014 0.032 0.036 0.054 0.089 0.643 0.893 1 1 1 1];
	output_sim = simf_doi__10_1093_nar_25_6_1203_figure4a(sim_x_opt,aTc);
	sim_pair_matrix = [output_sim; luc];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi: 10.1093/nar/25.6.1203', 'figure4a', aTc, output_sim, {'pN25 = 1*{alpha_pN25}'; 'tetR = {ALPHA_RBS_pN25}*pN25*{K_aTc}^{n_aTc}/({K_aTc}^{n_aTc} + aTc^{n_aTc})'; 'pLtetO_1 = 15*({beta_pLtetO_1} + ({alpha_pLtetO_1} - {beta_pLtetO_1})/(1 + tetR^4/{KE1_pLtetO_1} + tetR^2/{KE2_pLtetO_1}))'; 'luc = {scale_luc_doi__10_1093_nar_25_6_1203}*{ALPHA_RBSII}*pLtetO_1'});
	%%%%%%%%%%%%%
	IPTG = [0.001 0.01 0.03 0.05 0.075 0.1 0.15 0.2 0.3 0.5 1];
	eYFP = [0.015 0.014 0.023 0.062 0.098 0.187 0.239 0.374 0.572 0.869 1];
	output_sim = simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; eYFP];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'http://dspace.mit.edu/handle/1721.1/8228', 'figure4-6', IPTG, output_sim, {'placIq = 10*{alpha_placIq}'; 'lacI = {ALPHA_RBS_placI}*placIq*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'lacI_free = {ALPHA_RBS_placI}*placIq'; 'pLAC = 10*({beta_pLAC} + ({alpha_pLAC} - {beta_pLAC})/(1 + lacI^4/{KE1_pLAC} + (lacI_free^4 - lacI^4)/{KE2_pLAC}))'; 'eYFP = {scale_eYFP_http___dspace_mit_edu_handle_1721_1_8228}*{ALPHA_RBSII}*pLAC'});
	%%%%%%%%%%%%%
	aTc = [9.258 13.887 20.831 32.403 46.29 83.322 162.015 231.45 370.32 648.06 925.8 1388.7];
	eYFP = [0.003 0.003 0.003 0.003 0.005 0.01 0.05 0.167 0.333 0.667 0.833 1];
	output_sim = simf_doi_10_1073pnas_0408507102_figure2A(sim_x_opt,aTc);
	sim_pair_matrix = [output_sim; eYFP];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1073pnas.0408507102', 'figure2A', aTc, output_sim, {'placIq = 10*{alpha_placIq}'; 'tetR = {ALPHA_RBS_placI}*placIq*{K_aTc}^{n_aTc}/({K_aTc}^{n_aTc} + aTc^{n_aTc})'; 'pLtetO_1 = 10*({beta_pLtetO_1} + ({alpha_pLtetO_1} - {beta_pLtetO_1})/(1 + tetR^4/{KE1_pLtetO_1} + tetR^2/{KE2_pLtetO_1}))'; 'eYFP = {scale_eYFP_doi_10_1073pnas_0408507102}*{ALPHA_RBS_Unknown_Hooshangi}*pLtetO_1'});
	%%%%%%%%%%%%%
	IPTG = [0 0.005 0.01 0.02 0.05 0.1 0.2 1];
	lacZ = [0.001 0.003 0.043 0.5 0.984 0.999 1 1];
	output_sim = simf_doi_10_1073_pnas_0606717104_figure4_6(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; lacZ];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1073/pnas.0606717104', 'figure4-6', IPTG, output_sim, {'placI = 1*{alpha_placI}'; 'lacI = {ALPHA_RBS_placI}*placI*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'lacI_free = {ALPHA_RBS_placI}*placI'; 'pLAC = 1*({beta_pLAC} + ({alpha_pLAC} - {beta_pLAC})/(1 + lacI^4/{KE1_pLAC} + (lacI_free^4 - lacI^4)/{KE2_pLAC}))'; 'lacZ = {scale_lacZ_doi_10_1073_pnas_0606717104}*{ALPHA_RBS_pLAC}*pLAC'});
	%%%%%%%%%%%%%
	output_sim = simf_doi_10_1186_1754_1611_3_4_figure3_1(sim_x_opt);
	point_pair = [output_sim 	0.69];
	fprintf(sim_pair_fileID, '%f \t %f\n', point_pair);
	print_simulated_data_to_json(json_fileID, 'doi:10.1186/1754-1611-3-4', 'figure3-1', [], output_sim, {'BBa_J23101 = 10*{alpha_BBa_J23101}'; 'GFPmut3b = {scale_GFPmut3b_doi_10_1186_1754_1611_3_4}*{ALPHA_BBa_B0032}*BBa_J23101'});
	%%%%%%%%%%%%%
	output_sim = simf_doi_10_1186_1754_1611_3_4_figure3_2(sim_x_opt);
	point_pair = [output_sim 	0.018];
	fprintf(sim_pair_fileID, '%f \t %f\n', point_pair);
	print_simulated_data_to_json(json_fileID, 'doi:10.1186/1754-1611-3-4', 'figure3-2', [], output_sim, {'BBa_J23116 = 10*{alpha_BBa_J23116}'; 'GFPmut3b = {scale_GFPmut3b_doi_10_1186_1754_1611_3_4}*{ALPHA_BBa_B0032}*BBa_J23116'});
	%%%%%%%%%%%%%
	output_sim = simf_doi_10_1186_1754_1611_3_4_figure3_3(sim_x_opt);
	point_pair = [output_sim 	0.193];
	fprintf(sim_pair_fileID, '%f \t %f\n', point_pair);
	print_simulated_data_to_json(json_fileID, 'doi:10.1186/1754-1611-3-4', 'figure3-3', [], output_sim, {'BBa_J23150 = 10*{alpha_BBa_J23150}'; 'GFPmut3b = {scale_GFPmut3b_doi_10_1186_1754_1611_3_4}*{ALPHA_BBa_B0032}*BBa_J23150'});
	%%%%%%%%%%%%%
	output_sim = simf_doi_10_1186_1754_1611_3_4_figure3_4(sim_x_opt);
	point_pair = [output_sim 	0.4];
	fprintf(sim_pair_fileID, '%f \t %f\n', point_pair);
	print_simulated_data_to_json(json_fileID, 'doi:10.1186/1754-1611-3-4', 'figure3-4', [], output_sim, {'BBa_J23151 = 10*{alpha_BBa_J23151}'; 'GFPmut3b = {scale_GFPmut3b_doi_10_1186_1754_1611_3_4}*{ALPHA_BBa_B0032}*BBa_J23151'});
	%%%%%%%%%%%%%
	output_sim = simf_doi_10_1186_1754_1611_3_4_figure3_5(sim_x_opt);
	point_pair = [output_sim 	0.655];
	fprintf(sim_pair_fileID, '%f \t %f\n', point_pair);
	print_simulated_data_to_json(json_fileID, 'doi:10.1186/1754-1611-3-4', 'figure3-5', [], output_sim, {'BBa_J23102 = 10*{alpha_BBa_J23102}'; 'GFPmut3b = {scale_GFPmut3b_doi_10_1186_1754_1611_3_4}*{ALPHA_BBa_B0032}*BBa_J23102'});
	%%%%%%%%%%%%%
	output_sim = simf_doi_10_1186_1754_1611_3_4_figure3_6(sim_x_opt);
	point_pair = [output_sim 	0.897];
	fprintf(sim_pair_fileID, '%f \t %f\n', point_pair);
	print_simulated_data_to_json(json_fileID, 'doi:10.1186/1754-1611-3-4', 'figure3-6', [], output_sim, {'pLtetO_1 = 10*{alpha_pLtetO_1}'; 'GFPmut3b = {scale_GFPmut3b_doi_10_1186_1754_1611_3_4}*{ALPHA_BBa_B0032}*pLtetO_1'});
	%%%%%%%%%%%%%
	output_sim = simf_doi_10_1186_1754_1611_3_4_figure3_7(sim_x_opt);
	point_pair = [output_sim 	1];
	fprintf(sim_pair_fileID, '%f \t %f\n', point_pair);
	print_simulated_data_to_json(json_fileID, 'doi:10.1186/1754-1611-3-4', 'figure3-7', [], output_sim, {'pLlacO_1 = 10*{alpha_pLlacO_1}'; 'GFPmut3b = {scale_GFPmut3b_doi_10_1186_1754_1611_3_4}*{ALPHA_BBa_B0032}*pLlacO_1'});
	%%%%%%%%%%%%%
	output_sim = simf_doi_10_1093_nar_gku593_figure3D_1(sim_x_opt);
	point_pair = [output_sim 	0.002];
	fprintf(sim_pair_fileID, '%f \t %f\n', point_pair);
	print_simulated_data_to_json(json_fileID, 'doi:10.1093/nar/gku593', 'figure3D-1', [], output_sim, {'BBa_J23109 = 10*{alpha_BBa_J23109}'; 'GFPmut3b = {scale_GFPmut3b_doi_10_1093_nar_gku593}*{ALPHA_BBa_B0030}*BBa_J23109'});
	%%%%%%%%%%%%%
	output_sim = simf_doi_10_1093_nar_gku593_figure3D_2(sim_x_opt);
	point_pair = [output_sim 	0.018];
	fprintf(sim_pair_fileID, '%f \t %f\n', point_pair);
	print_simulated_data_to_json(json_fileID, 'doi:10.1093/nar/gku593', 'figure3D-2', [], output_sim, {'BBa_J23114 = 10*{alpha_BBa_J23114}'; 'GFPmut3b = {scale_GFPmut3b_doi_10_1093_nar_gku593}*{ALPHA_BBa_B0030}*BBa_J23114'});
	%%%%%%%%%%%%%
	output_sim = simf_doi_10_1093_nar_gku593_figure3D_3(sim_x_opt);
	point_pair = [output_sim 	0.036];
	fprintf(sim_pair_fileID, '%f \t %f\n', point_pair);
	print_simulated_data_to_json(json_fileID, 'doi:10.1093/nar/gku593', 'figure3D-3', [], output_sim, {'BBa_J23116 = 10*{alpha_BBa_J23116}'; 'GFPmut3b = {scale_GFPmut3b_doi_10_1093_nar_gku593}*{ALPHA_BBa_B0030}*BBa_J23116'});
	%%%%%%%%%%%%%
	output_sim = simf_doi_10_1093_nar_gku593_figure3D_4(sim_x_opt);
	point_pair = [output_sim 	0.05];
	fprintf(sim_pair_fileID, '%f \t %f\n', point_pair);
	print_simulated_data_to_json(json_fileID, 'doi:10.1093/nar/gku593', 'figure3D-4', [], output_sim, {'BBa_J23115 = 10*{alpha_BBa_J23115}'; 'GFPmut3b = {scale_GFPmut3b_doi_10_1093_nar_gku593}*{ALPHA_BBa_B0030}*BBa_J23115'});
	%%%%%%%%%%%%%
	output_sim = simf_doi_10_1093_nar_gku593_figure3D_5(sim_x_opt);
	point_pair = [output_sim 	0.086];
	fprintf(sim_pair_fileID, '%f \t %f\n', point_pair);
	print_simulated_data_to_json(json_fileID, 'doi:10.1093/nar/gku593', 'figure3D-5', [], output_sim, {'BBa_J23105 = 10*{alpha_BBa_J23105}'; 'GFPmut3b = {scale_GFPmut3b_doi_10_1093_nar_gku593}*{ALPHA_BBa_B0030}*BBa_J23105'});
	%%%%%%%%%%%%%
	output_sim = simf_doi_10_1093_nar_gku593_figure3D_6(sim_x_opt);
	point_pair = [output_sim 	0.214];
	fprintf(sim_pair_fileID, '%f \t %f\n', point_pair);
	print_simulated_data_to_json(json_fileID, 'doi:10.1093/nar/gku593', 'figure3D-6', [], output_sim, {'BBa_J23106 = 10*{alpha_BBa_J23106}'; 'GFPmut3b = {scale_GFPmut3b_doi_10_1093_nar_gku593}*{ALPHA_BBa_B0030}*BBa_J23106'});
	%%%%%%%%%%%%%
	arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2];
	GFPmut3b = [0 0.002 0.041 0.479 0.949 0.997 1 1 1];
	output_sim = simf_doi_10_1093_nar_gku593_figureS6(sim_x_opt,arabinose);
	sim_pair_matrix = [output_sim; GFPmut3b];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1093/nar/gku593', 'figureS6', arabinose, output_sim, {'pC = 10*{alpha_pC}'; 'araC = {ALPHA_RBS_pC}*pC*{K_arabinose}^{n_arabinose}/({K_arabinose}^{n_arabinose} + arabinose^{n_arabinose})'; 'pBAD = 10*({beta_pBAD} + ({alpha_pBAD} - {beta_pBAD})*{K_pBAD}^{n_pBAD}/({K_pBAD}^{n_pBAD} + araC^{n_pBAD}))'; 'GFPmut3b = {scale_GFPmut3b_doi_10_1093_nar_gku593}*{ALPHA_BBa_B0030}*pBAD'});
	%%%%%%%%%%%%%
	output_sim = simf_doi_10_1093_nar_gku1388_figure2A_1(sim_x_opt);
	point_pair = [output_sim 	1];
	fprintf(sim_pair_fileID, '%f \t %f\n', point_pair);
	print_simulated_data_to_json(json_fileID, 'doi:10.1093/nar/gku1388', 'figure2A-1', [], output_sim, {'BBa_J23101 = 10*{alpha_BBa_J23101}'; 'GFPmut3b = {scale_GFPmut3b_doi_10_1093_nar_gku1388}*{ALPHA_BBa_B0030}*BBa_J23101'});
	%%%%%%%%%%%%%
	output_sim = simf_doi_10_1093_nar_gku1388_figure2A_2(sim_x_opt);
	point_pair = [output_sim 	0.215];
	fprintf(sim_pair_fileID, '%f \t %f\n', point_pair);
	print_simulated_data_to_json(json_fileID, 'doi:10.1093/nar/gku1388', 'figure2A-2', [], output_sim, {'BBa_J23106 = 10*{alpha_BBa_J23106}'; 'GFPmut3b = {scale_GFPmut3b_doi_10_1093_nar_gku1388}*{ALPHA_BBa_B0030}*BBa_J23106'});
	%%%%%%%%%%%%%
	output_sim = simf_doi_10_1093_nar_gku1388_figure2A_3(sim_x_opt);
	point_pair = [output_sim 	0.095];
	fprintf(sim_pair_fileID, '%f \t %f\n', point_pair);
	print_simulated_data_to_json(json_fileID, 'doi:10.1093/nar/gku1388', 'figure2A-3', [], output_sim, {'BBa_J23105 = 10*{alpha_BBa_J23105}'; 'GFPmut3b = {scale_GFPmut3b_doi_10_1093_nar_gku1388}*{ALPHA_BBa_B0030}*BBa_J23105'});
	%%%%%%%%%%%%%
	output_sim = simf_doi_10_1093_nar_gku1388_figure2A_4(sim_x_opt);
	point_pair = [output_sim 	0.046];
	fprintf(sim_pair_fileID, '%f \t %f\n', point_pair);
	print_simulated_data_to_json(json_fileID, 'doi:10.1093/nar/gku1388', 'figure2A-4', [], output_sim, {'BBa_J23115 = 10*{alpha_BBa_J23115}'; 'GFPmut3b = {scale_GFPmut3b_doi_10_1093_nar_gku1388}*{ALPHA_BBa_B0030}*BBa_J23115'});
	%%%%%%%%%%%%%
	output_sim = simf_doi_10_1093_nar_gku1388_figure2A_5(sim_x_opt);
	point_pair = [output_sim 	0.018];
	fprintf(sim_pair_fileID, '%f \t %f\n', point_pair);
	print_simulated_data_to_json(json_fileID, 'doi:10.1093/nar/gku1388', 'figure2A-5', [], output_sim, {'BBa_J23114 = 10*{alpha_BBa_J23114}'; 'GFPmut3b = {scale_GFPmut3b_doi_10_1093_nar_gku1388}*{ALPHA_BBa_B0030}*BBa_J23114'});
	%%%%%%%%%%%%%
	output_sim = simf_doi_10_1093_nar_gku1388_figure2A_6(sim_x_opt);
	point_pair = [output_sim 	0.003];
	fprintf(sim_pair_fileID, '%f \t %f\n', point_pair);
	print_simulated_data_to_json(json_fileID, 'doi:10.1093/nar/gku1388', 'figure2A-6', [], output_sim, {'BBa_J23117 = 10*{alpha_BBa_J23117}'; 'GFPmut3b = {scale_GFPmut3b_doi_10_1093_nar_gku1388}*{ALPHA_BBa_B0030}*BBa_J23117'});
	%%%%%%%%%%%%%
	aTc = [2 4 10 25 50 100 200];
	GFPmut3b = [0.002 0.014 0.167 0.666 0.817 0.841 0.831];
	output_sim = simf_doi_10_1093_nar_gku1388_figure2B_1(sim_x_opt,aTc);
	sim_pair_matrix = [output_sim; GFPmut3b];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1093/nar/gku1388', 'figure2B-1', aTc, output_sim, {'BBa_J23117 = 10*{alpha_BBa_J23117}'; 'tetR = {ALPHA_BBa_B0030}*BBa_J23117*{K_aTc}^{n_aTc}/({K_aTc}^{n_aTc} + aTc^{n_aTc})'; 'pTET_star_ = 10*({beta_pTET_star_} + ({alpha_pTET_star_} - {beta_pTET_star_})/(1 + tetR^4/{KE1_pTET_star_} + tetR^2/{KE2_pTET_star_}))'; 'GFPmut3b = {scale_GFPmut3b_doi_10_1093_nar_gku1388}*{ALPHA_BBa_B0030}*pTET_star_'});
	%%%%%%%%%%%%%
	aTc = [10 25 50 100 200];
	GFPmut3b = [0.001 0.022 0.176 0.411 0.461];
	output_sim = simf_doi_10_1093_nar_gku1388_figure2B_2(sim_x_opt,aTc);
	sim_pair_matrix = [output_sim; GFPmut3b];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1093/nar/gku1388', 'figure2B-2', aTc, output_sim, {'BBa_J23114 = 10*{alpha_BBa_J23114}'; 'tetR = {ALPHA_BBa_B0030}*BBa_J23114*{K_aTc}^{n_aTc}/({K_aTc}^{n_aTc} + aTc^{n_aTc})'; 'pTET_star_ = 10*({beta_pTET_star_} + ({alpha_pTET_star_} - {beta_pTET_star_})/(1 + tetR^4/{KE1_pTET_star_} + tetR^2/{KE2_pTET_star_}))'; 'GFPmut3b = {scale_GFPmut3b_doi_10_1093_nar_gku1388}*{ALPHA_BBa_B0030}*pTET_star_'});
	%%%%%%%%%%%%%
	arabinose = [0 0.001 0.015 0.05 0.5 5 10];
	eYFP = [0.002 0.002 0.057 0.628 0.999 1 1];
	output_sim = simf_doi_10_1038_nature09565_figure1C(sim_x_opt,arabinose);
	sim_pair_matrix = [output_sim; eYFP];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1038/nature09565', 'figure1C', arabinose, output_sim, {'BBa_J23117 = 15*{alpha_BBa_J23117}'; 'araC = {ALPHA_BBa_B0034}*BBa_J23117*{K_arabinose}^{n_arabinose}/({K_arabinose}^{n_arabinose} + arabinose^{n_arabinose})'; 'pBAD = 15*({beta_pBAD} + ({alpha_pBAD} - {beta_pBAD})*{K_pBAD}^{n_pBAD}/({K_pBAD}^{n_pBAD} + araC^{n_pBAD}))'; 'eYFP = {scale_eYFP_doi_10_1038_nature09565}*{ALPHA_BBa_B0033}*pBAD'});
	%%%%%%%%%%%%%
	aTc = [0 0.025 0.25 2.5 25 100 250];
	eYFP = [0.005 0.006 0.007 0.03 0.297 0.388 0.398];
	output_sim = simf_doi_10_1038_nature09565_figureS1C(sim_x_opt,aTc);
	sim_pair_matrix = [output_sim; eYFP];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1038/nature09565', 'figureS1C', aTc, output_sim, {'BBa_J23117 = 15*{alpha_BBa_J23117}'; 'tetR = {ALPHA_BBa_B0032}*BBa_J23117*{K_aTc}^{n_aTc}/({K_aTc}^{n_aTc} + aTc^{n_aTc})'; 'pLtetO_1 = 15*({beta_pLtetO_1} + ({alpha_pLtetO_1} - {beta_pLtetO_1})/(1 + tetR^4/{KE1_pLtetO_1} + tetR^2/{KE2_pLtetO_1}))'; 'eYFP = {scale_eYFP_doi_10_1038_nature09565}*{ALPHA_BBa_B0033}*pLtetO_1'});
	%%%%%%%%%%%%%
	IPTG = [0 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1];
	sfGFP = [0.025 0.025 0.035 0.055 0.115 0.2 0.3 0.4 0.65 1];
	output_sim = simf_doi_10_1038_nbt_2401_figureS11(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; sfGFP];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1038/nbt.2401', 'figureS11', IPTG, output_sim, {'BBa_J23101 = 5*{alpha_BBa_J23101}'; 'lacI = {ALPHA_RBS_A}*BBa_J23101*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'lacI_free = {ALPHA_RBS_A}*BBa_J23101'; 'pTAC = 5*({beta_pTAC} + ({alpha_pTAC} - {beta_pTAC})/(1 + lacI^4/{KE1_pTAC} + (lacI_free^4 - lacI^4)/{KE2_pTAC}))'; 'sfGFP = {scale_sfGFP_doi_10_1038_nbt_2401}*{ALPHA_RBS_A}*pTAC'});
	%%%%%%%%%%%%%
	arabinose = [0.1 1 2 5 7 10 12.5 25 37.5 50 62.5];
	sfGFP = [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425];
	output_sim = simf_doi_10_1038_nbt_2401_figureS13(sim_x_opt,arabinose);
	sim_pair_matrix = [output_sim; sfGFP];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1038/nbt.2401', 'figureS13', arabinose, output_sim, {'BBa_J23101 = 5*{alpha_BBa_J23101}'; 'araC = {ALPHA_RBS_A}*BBa_J23101*{K_arabinose}^{n_arabinose}/({K_arabinose}^{n_arabinose} + arabinose^{n_arabinose})'; 'pBAD = 5*({beta_pBAD} + ({alpha_pBAD} - {beta_pBAD})*{K_pBAD}^{n_pBAD}/({K_pBAD}^{n_pBAD} + araC^{n_pBAD}))'; 'sfGFP = {scale_sfGFP_doi_10_1038_nbt_2401}*{ALPHA_RBS_A}*pBAD'});
	%%%%%%%%%%%%%
	arabinose = [0 0 0.002 0.01 0.04 0.16 1 4 25];
	mRFP1 = [0.003 0.003 0.003 0.017 0.23 0.82 0.939 0.941 0.941];
	output_sim = simf_doi_10_1038_nature11516_figureS8A(sim_x_opt,arabinose);
	sim_pair_matrix = [output_sim; mRFP1];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1038/nature11516', 'figureS8A', arabinose, output_sim, {'pC = 15*{alpha_pC}'; 'araC = {ALPHA_RBS_pC}*pC*{K_arabinose}^{n_arabinose}/({K_arabinose}^{n_arabinose} + arabinose^{n_arabinose})'; 'pBAD = 15*({beta_pBAD} + ({alpha_pBAD} - {beta_pBAD})*{K_pBAD}^{n_pBAD}/({K_pBAD}^{n_pBAD} + araC^{n_pBAD}))'; 'mRFP1 = {scale_mRFP1_doi_10_1038_nature11516}*{ALPHA_RBS_psicA}*pBAD'});
	%%%%%%%%%%%%%
	IPTG = [0 0 0 0.001 0.004 0.025 0.1 0.5 2.5];
	mRFP1 = [0.027 0.027 0.027 0.03 0.067 0.46 0.884 0.991 1];
	output_sim = simf_doi_10_1038_nature11516_figureS8B(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; mRFP1];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1038/nature11516', 'figureS8B', IPTG, output_sim, {'placIq = 15*{alpha_placIq}'; 'lacI = {ALPHA_RBS_placI}*placIq*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'lacI_free = {ALPHA_RBS_placI}*placIq'; 'pTAC = 15*({beta_pTAC} + ({alpha_pTAC} - {beta_pTAC})/(1 + lacI^4/{KE1_pTAC} + (lacI_free^4 - lacI^4)/{KE2_pTAC}))'; 'mRFP1 = {scale_mRFP1_doi_10_1038_nature11516}*{ALPHA_RBS_psicA}*pTAC'});
	%%%%%%%%%%%%%
	IPTG = [0 0.001 0.005 0.01 0.025 0.05 0.1 0.25 0.5];
	mRFP1 = [0.02 0.031 0.092 0.306 0.745 0.898 0.969 1 1];
	output_sim = simf_doi_10_1038_nbt_2149_Supplement_page17(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; mRFP1];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1038/nbt.2149', 'Supplement-page17', IPTG, output_sim, {'placI = 10*{alpha_placI}'; 'lacI = {ALPHA_RBS_placI}*placI*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'lacI_free = {ALPHA_RBS_placI}*placI'; 'placUV5 = 10*({beta_placUV5} + ({alpha_placUV5} - {beta_placUV5})/(1 + lacI^4/{KE1_placUV5} + (lacI_free^4 - lacI^4)/{KE2_placUV5}))'; 'mRFP1 = {scale_mRFP1_doi_10_1038_nbt_2149}*{ALPHA_RBS_pET_29b}*placUV5'});
	fclose(sim_pair_fileID);
	fclose(json_fileID);
end

function sim_LOOCV(ensemble_size, ite_num, func_ite_num)
	lb = [0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01 0.03 0.002 0.01 0.03 0.03 0.001 0.001 0.001 0.002 0.002 0 1 0.01 0.03 0.001 0.001 0.001 0.002 0.002 0 1 0.01 0.03 0.01 0.03 0.002 0.01 0.03 0.002 0.01 0.002 0.002 0.002 0.002 0.002 0.01 0.002 0.002 0.002 0.002 0.01 0.002 0.001 0.001 0.002 0 0.03 0.03 0.01 0.03 0.001 0.001 0.002 0 0.01 0.03 0.01 0.03 0.001 0.001 0.002 0 0.01];
	ub = [5 5 50 50 50 50 1 4 4 10 5 50 10 5 5 50 50 50 50 50 1 4 10 5 50 50 50 50 50 1 4 10 5 10 5 50 10 5 50 10 50 50 50 50 50 10 50 50 50 50 10 50 50 50 50 1 5 5 10 5 50 50 50 1 10 5 10 5 50 50 50 1 10];
	LOOCV_sim_data_fileID = fopen('Plot_data/LOOCV_sim_data_Ellis.dat','w');
	LOO_model_name_list = {'doi_10_1186_1754_1611_3_4_figure3_1', 'doi_10_1186_1754_1611_3_4_figure3_2', 'doi_10_1186_1754_1611_3_4_figure3_6', 'doi_10_1186_1754_1611_3_4_figure3_7', 'doi_10_1093_nar_gku593_figure3D_2', 'doi_10_1093_nar_gku593_figure3D_3', 'doi_10_1093_nar_gku593_figure3D_4', 'doi_10_1093_nar_gku593_figure3D_5', 'doi_10_1093_nar_gku593_figure3D_6', 'doi_10_1093_nar_gku593_figureS6', 'doi_10_1093_nar_gku1388_figure2A_1', 'doi_10_1093_nar_gku1388_figure2A_2', 'doi_10_1093_nar_gku1388_figure2A_3', 'doi_10_1093_nar_gku1388_figure2A_4', 'doi_10_1093_nar_gku1388_figure2A_5', 'doi_10_1093_nar_gku1388_figure2A_6', 'doi_10_1093_nar_gku1388_figure2B_1', 'doi_10_1093_nar_gku1388_figure2B_2', 'doi_10_1038_nature09565_figureS1C', 'doi_10_1038_nbt_2401_figureS11', 'doi_10_1038_nbt_2401_figureS13', 'doi_10_1038_nature11516_figureS8A', 'doi_10_1038_nature11516_figureS8B'};
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1186_1754_1611_3_4_figure3_1, lb, ub);
	desired_output = 0.69;
	output = simf_doi_10_1186_1754_1611_3_4_figure3_1(x_opt);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1186_1754_1611_3_4_figure3_1\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1186_1754_1611_3_4_figure3_2, lb, ub);
	desired_output = 0.018;
	output = simf_doi_10_1186_1754_1611_3_4_figure3_2(x_opt);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1186_1754_1611_3_4_figure3_2\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1186_1754_1611_3_4_figure3_6, lb, ub);
	desired_output = 0.897;
	output = simf_doi_10_1186_1754_1611_3_4_figure3_6(x_opt);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1186_1754_1611_3_4_figure3_6\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1186_1754_1611_3_4_figure3_7, lb, ub);
	desired_output = 1;
	output = simf_doi_10_1186_1754_1611_3_4_figure3_7(x_opt);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1186_1754_1611_3_4_figure3_7\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1093_nar_gku593_figure3D_2, lb, ub);
	desired_output = 0.018;
	output = simf_doi_10_1093_nar_gku593_figure3D_2(x_opt);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1093_nar_gku593_figure3D_2\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1093_nar_gku593_figure3D_3, lb, ub);
	desired_output = 0.036;
	output = simf_doi_10_1093_nar_gku593_figure3D_3(x_opt);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1093_nar_gku593_figure3D_3\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1093_nar_gku593_figure3D_4, lb, ub);
	desired_output = 0.05;
	output = simf_doi_10_1093_nar_gku593_figure3D_4(x_opt);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1093_nar_gku593_figure3D_4\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1093_nar_gku593_figure3D_5, lb, ub);
	desired_output = 0.086;
	output = simf_doi_10_1093_nar_gku593_figure3D_5(x_opt);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1093_nar_gku593_figure3D_5\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1093_nar_gku593_figure3D_6, lb, ub);
	desired_output = 0.214;
	output = simf_doi_10_1093_nar_gku593_figure3D_6(x_opt);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1093_nar_gku593_figure3D_6\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1093_nar_gku593_figureS6, lb, ub);
	arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2];
	desired_output = [0 0.002 0.041 0.479 0.949 0.997 1 1 1];
	output = simf_doi_10_1093_nar_gku593_figureS6(x_opt, arabinose);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1093_nar_gku593_figureS6\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2A_1, lb, ub);
	desired_output = 1;
	output = simf_doi_10_1093_nar_gku1388_figure2A_1(x_opt);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1093_nar_gku1388_figure2A_1\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2A_2, lb, ub);
	desired_output = 0.215;
	output = simf_doi_10_1093_nar_gku1388_figure2A_2(x_opt);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1093_nar_gku1388_figure2A_2\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2A_3, lb, ub);
	desired_output = 0.095;
	output = simf_doi_10_1093_nar_gku1388_figure2A_3(x_opt);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1093_nar_gku1388_figure2A_3\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2A_4, lb, ub);
	desired_output = 0.046;
	output = simf_doi_10_1093_nar_gku1388_figure2A_4(x_opt);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1093_nar_gku1388_figure2A_4\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2A_5, lb, ub);
	desired_output = 0.018;
	output = simf_doi_10_1093_nar_gku1388_figure2A_5(x_opt);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1093_nar_gku1388_figure2A_5\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2A_6, lb, ub);
	desired_output = 0.003;
	output = simf_doi_10_1093_nar_gku1388_figure2A_6(x_opt);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1093_nar_gku1388_figure2A_6\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2B_1, lb, ub);
	aTc = [2 4 10 25 50 100 200];
	desired_output = [0.002 0.014 0.167 0.666 0.817 0.841 0.831];
	output = simf_doi_10_1093_nar_gku1388_figure2B_1(x_opt, aTc);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1093_nar_gku1388_figure2B_1\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2B_2, lb, ub);
	aTc = [10 25 50 100 200];
	desired_output = [0.001 0.022 0.176 0.411 0.461];
	output = simf_doi_10_1093_nar_gku1388_figure2B_2(x_opt, aTc);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1093_nar_gku1388_figure2B_2\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1038_nature09565_figureS1C, lb, ub);
	aTc = [0 0.025 0.25 2.5 25 100 250];
	desired_output = [0.005 0.006 0.007 0.03 0.297 0.388 0.398];
	output = simf_doi_10_1038_nature09565_figureS1C(x_opt, aTc);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1038_nature09565_figureS1C\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1038_nbt_2401_figureS11, lb, ub);
	IPTG = [0 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1];
	desired_output = [0.025 0.025 0.035 0.055 0.115 0.2 0.3 0.4 0.65 1];
	output = simf_doi_10_1038_nbt_2401_figureS11(x_opt, IPTG);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1038_nbt_2401_figureS11\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1038_nbt_2401_figureS13, lb, ub);
	arabinose = [0.1 1 2 5 7 10 12.5 25 37.5 50 62.5];
	desired_output = [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425];
	output = simf_doi_10_1038_nbt_2401_figureS13(x_opt, arabinose);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1038_nbt_2401_figureS13\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1038_nature11516_figureS8A, lb, ub);
	arabinose = [0 0 0.002 0.01 0.04 0.16 1 4 25];
	desired_output = [0.003 0.003 0.003 0.017 0.23 0.82 0.939 0.941 0.941];
	output = simf_doi_10_1038_nature11516_figureS8A(x_opt, arabinose);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1038_nature11516_figureS8A\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1038_nature11516_figureS8B, lb, ub);
	IPTG = [0 0 0 0.001 0.004 0.025 0.1 0.5 2.5];
	desired_output = [0.027 0.027 0.027 0.03 0.067 0.46 0.884 0.991 1];
	output = simf_doi_10_1038_nature11516_figureS8B(x_opt, IPTG);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1038_nature11516_figureS8B\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	fclose(LOOCV_sim_data_fileID);
end
function error = LOOCV_sim_error_doi_10_1186_1754_1611_3_4_figure3_1(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 0, 0);
end
function error = LOOCV_sim_error_doi_10_1186_1754_1611_3_4_figure3_2(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 1, 1);
end
function error = LOOCV_sim_error_doi_10_1186_1754_1611_3_4_figure3_6(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 2, 2);
end
function error = LOOCV_sim_error_doi_10_1186_1754_1611_3_4_figure3_7(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 3, 3);
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku593_figure3D_2(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 4, 4);
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku593_figure3D_3(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 5, 5);
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku593_figure3D_4(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 6, 6);
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku593_figure3D_5(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 7, 7);
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku593_figure3D_6(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 8, 8);
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku593_figureS6(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 9, 17);
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2A_1(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 18, 18);
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2A_2(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 19, 19);
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2A_3(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 20, 20);
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2A_4(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 21, 21);
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2A_5(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 22, 22);
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2A_6(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 23, 23);
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2B_1(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 24, 30);
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2B_2(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 31, 35);
end
function error = LOOCV_sim_error_doi_10_1038_nature09565_figureS1C(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 36, 42);
end
function error = LOOCV_sim_error_doi_10_1038_nbt_2401_figureS11(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 43, 52);
end
function error = LOOCV_sim_error_doi_10_1038_nbt_2401_figureS13(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 53, 63);
end
function error = LOOCV_sim_error_doi_10_1038_nature11516_figureS8A(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 64, 72);
end
function error = LOOCV_sim_error_doi_10_1038_nature11516_figureS8B(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 73, 81);
end

