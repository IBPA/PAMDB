function real_benchmark_fitting
	sim_x_opt = simultaneous_fitting(1,1000,50000);
	print_plotting_sim_data(sim_x_opt);
	%sim_LOOCV(2,100,1000);
end
function output = simf_doi_10_1128_AEM_00791_07_figure3a(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pBAD = x(1);
	ALPHA_RBS_pC = x(2);
	K_L_arabinose = x(3);
	K_pBAD = x(4);
	alpha_pBAD = x(5);
	alpha_pC = x(6);
	beta_pBAD = x(7);
	n_L_arabinose = x(8);
	n_pBAD = x(9);
	scale_GFPuv_doi_10_1128_AEM_00791_07 = x(10);
	for i = 1:length(input)
		pC = 20*alpha_pC;
		L_arabinose = input(i);
		araC = ALPHA_RBS_pC*pC*K_L_arabinose^n_L_arabinose/(K_L_arabinose^n_L_arabinose + L_arabinose^n_L_arabinose);
		pBAD = 20*(beta_pBAD + (alpha_pBAD - beta_pBAD)*K_pBAD^n_pBAD/(K_pBAD^n_pBAD + araC^n_pBAD));
		GFPuv = scale_GFPuv_doi_10_1128_AEM_00791_07*ALPHA_RBS_pBAD*pBAD;
		output(i) = GFPuv;
	end
end

function output = simf_doi_10_1038_nature12148_figure16a(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0030 = x(11);
	K_L_arabinose = x(3);
	K_pBAD = x(4);
	alpha_pBAD = x(5);
	alpha_pLlacO_1 = x(12);
	beta_pBAD = x(7);
	n_L_arabinose = x(8);
	n_pBAD = x(9);
	scale_mCherry_doi_10_1038_nature12148 = x(13);
	for i = 1:length(input)
		pLlacO_1 = 5*alpha_pLlacO_1;
		L_arabinose = input(i);
		araC = ALPHA_BBa_B0030*pLlacO_1*K_L_arabinose^n_L_arabinose/(K_L_arabinose^n_L_arabinose + L_arabinose^n_L_arabinose);
		pBAD = 15*(beta_pBAD + (alpha_pBAD - beta_pBAD)*K_pBAD^n_pBAD/(K_pBAD^n_pBAD + araC^n_pBAD));
		mCherry = scale_mCherry_doi_10_1038_nature12148*ALPHA_BBa_B0030*pBAD;
		output(i) = mCherry;
	end
end

function output = simf_http___jb_asm_org_content_177_14_4121_short_figure4a(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pBAD = x(1);
	ALPHA_RBS_pC = x(2);
	K_L_arabinose = x(3);
	K_pBAD = x(4);
	alpha_pBAD = x(5);
	alpha_pC = x(6);
	beta_pBAD = x(7);
	n_L_arabinose = x(8);
	n_pBAD = x(9);
	scale_AP_http___jb_asm_org_content_177_14_4121_short = x(14);
	for i = 1:length(input)
		pC = 20*alpha_pC;
		L_arabinose = input(i);
		araC = ALPHA_RBS_pC*pC*K_L_arabinose^n_L_arabinose/(K_L_arabinose^n_L_arabinose + L_arabinose^n_L_arabinose);
		pBAD = 20*(beta_pBAD + (alpha_pBAD - beta_pBAD)*K_pBAD^n_pBAD/(K_pBAD^n_pBAD + araC^n_pBAD));
		AP = scale_AP_http___jb_asm_org_content_177_14_4121_short*ALPHA_RBS_pBAD*pBAD;
		output(i) = AP;
	end
end

function output = simf_doi_10_1093_nar_25_6_1203_figure4a(x, input)
	output = zeros(1,length(input));
	ALPHA_RBSII = x(15);
	ALPHA_RBS_pN25 = x(16);
	K_aTc = x(17);
	K_pLtetO_1 = x(18);
	alpha_pLtetO_1 = x(19);
	alpha_pN25 = x(20);
	beta_pLtetO_1 = x(21);
	n_aTc = x(22);
	n_pLtetO_1 = x(23);
	scale_luc_doi_10_1093_nar_25_6_1203 = x(24);
	for i = 1:length(input)
		pN25 = 1*alpha_pN25;
		aTc = input(i);
		tetR = ALPHA_RBS_pN25*pN25*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pLtetO_1 = 15*(beta_pLtetO_1 + (alpha_pLtetO_1 - beta_pLtetO_1)*K_pLtetO_1^n_pLtetO_1/(K_pLtetO_1^n_pLtetO_1 + tetR^n_pLtetO_1));
		luc = scale_luc_doi_10_1093_nar_25_6_1203*ALPHA_RBSII*pLtetO_1;
		output(i) = luc;
	end
end

function output = simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6(x, input)
	output = zeros(1,length(input));
	ALPHA_RBSII = x(15);
	ALPHA_RBS_placI = x(25);
	K_IPTG = x(26);
	K_pLAC = x(27);
	alpha_pLAC = x(28);
	alpha_placIq = x(29);
	beta_pLAC = x(30);
	n_IPTG = x(31);
	n_pLAC = x(32);
	scale_eYFP_http___dspace_mit_edu_handle_1721_1_8228 = x(33);
	for i = 1:length(input)
		placIq = 10*alpha_placIq;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placIq*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pLAC = 10*(beta_pLAC + (alpha_pLAC - beta_pLAC)*K_pLAC^n_pLAC/(K_pLAC^n_pLAC + lacI^n_pLAC));
		eYFP = scale_eYFP_http___dspace_mit_edu_handle_1721_1_8228*ALPHA_RBSII*pLAC;
		output(i) = eYFP;
	end
end

function output = simf_doi_10_1073pnas_0408507102_figure2A(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_Unknown_Hooshangi = x(34);
	ALPHA_RBS_placI = x(25);
	K_aTc = x(17);
	K_pLtetO_1 = x(18);
	alpha_pLtetO_1 = x(19);
	alpha_placIq = x(29);
	beta_pLtetO_1 = x(21);
	n_aTc = x(22);
	n_pLtetO_1 = x(23);
	scale_eYFP_doi_10_1073pnas_0408507102 = x(35);
	for i = 1:length(input)
		placIq = 10*alpha_placIq;
		aTc = input(i);
		tetR = ALPHA_RBS_placI*placIq*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pLtetO_1 = 10*(beta_pLtetO_1 + (alpha_pLtetO_1 - beta_pLtetO_1)*K_pLtetO_1^n_pLtetO_1/(K_pLtetO_1^n_pLtetO_1 + tetR^n_pLtetO_1));
		eYFP = scale_eYFP_doi_10_1073pnas_0408507102*ALPHA_RBS_Unknown_Hooshangi*pLtetO_1;
		output(i) = eYFP;
	end
end

function output = simf_doi_10_1073_pnas_1316298111_figure2B(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pLAC = x(36);
	ALPHA_RBS_placI = x(25);
	K_IPTG = x(26);
	K_pLAC = x(27);
	alpha_pLAC = x(28);
	alpha_placI = x(37);
	beta_pLAC = x(30);
	n_IPTG = x(31);
	n_pLAC = x(32);
	scale_lacZ_doi_10_1073_pnas_1316298111 = x(38);
	for i = 1:length(input)
		placI = 1*alpha_placI;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placI*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pLAC = 1*(beta_pLAC + (alpha_pLAC - beta_pLAC)*K_pLAC^n_pLAC/(K_pLAC^n_pLAC + lacI^n_pLAC));
		lacZ = scale_lacZ_doi_10_1073_pnas_1316298111*ALPHA_RBS_pLAC*pLAC;
		output(i) = lacZ;
	end
end

function output = simf_doi_10_1073_pnas_0606717104_figure4_6(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pLAC = x(36);
	ALPHA_RBS_placI = x(25);
	K_IPTG = x(26);
	K_pLAC = x(27);
	alpha_pLAC = x(28);
	alpha_placI = x(37);
	beta_pLAC = x(30);
	n_IPTG = x(31);
	n_pLAC = x(32);
	scale_lacZ_doi_10_1073_pnas_0606717104 = x(39);
	for i = 1:length(input)
		placI = 1*alpha_placI;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placI*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pLAC = 1*(beta_pLAC + (alpha_pLAC - beta_pLAC)*K_pLAC^n_pLAC/(K_pLAC^n_pLAC + lacI^n_pLAC));
		lacZ = scale_lacZ_doi_10_1073_pnas_0606717104*ALPHA_RBS_pLAC*pLAC;
		output(i) = lacZ;
	end
end

function output = simf_doi_10_1186_1754_1611_3_4_figure3_1(x)
	ALPHA_BBa_B0032 = x(40);
	alpha_BBa_J23101 = x(41);
	scale_GFPmut3b_doi_10_1186_1754_1611_3_4 = x(42);
	BBa_J23101 = 10*alpha_BBa_J23101;
	GFPmut3b = scale_GFPmut3b_doi_10_1186_1754_1611_3_4*ALPHA_BBa_B0032*BBa_J23101;
	output = GFPmut3b;
end

function output = simf_doi_10_1186_1754_1611_3_4_figure3_2(x)
	ALPHA_BBa_B0032 = x(40);
	alpha_BBa_J23116 = x(43);
	scale_GFPmut3b_doi_10_1186_1754_1611_3_4 = x(42);
	BBa_J23116 = 10*alpha_BBa_J23116;
	GFPmut3b = scale_GFPmut3b_doi_10_1186_1754_1611_3_4*ALPHA_BBa_B0032*BBa_J23116;
	output = GFPmut3b;
end

function output = simf_doi_10_1186_1754_1611_3_4_figure3_3(x)
	ALPHA_BBa_B0032 = x(40);
	alpha_BBa_J23150 = x(44);
	scale_GFPmut3b_doi_10_1186_1754_1611_3_4 = x(42);
	BBa_J23150 = 10*alpha_BBa_J23150;
	GFPmut3b = scale_GFPmut3b_doi_10_1186_1754_1611_3_4*ALPHA_BBa_B0032*BBa_J23150;
	output = GFPmut3b;
end

function output = simf_doi_10_1186_1754_1611_3_4_figure3_4(x)
	ALPHA_BBa_B0032 = x(40);
	alpha_BBa_J23151 = x(45);
	scale_GFPmut3b_doi_10_1186_1754_1611_3_4 = x(42);
	BBa_J23151 = 10*alpha_BBa_J23151;
	GFPmut3b = scale_GFPmut3b_doi_10_1186_1754_1611_3_4*ALPHA_BBa_B0032*BBa_J23151;
	output = GFPmut3b;
end

function output = simf_doi_10_1186_1754_1611_3_4_figure3_5(x)
	ALPHA_BBa_B0032 = x(40);
	alpha_BBa_J23102 = x(46);
	scale_GFPmut3b_doi_10_1186_1754_1611_3_4 = x(42);
	BBa_J23102 = 10*alpha_BBa_J23102;
	GFPmut3b = scale_GFPmut3b_doi_10_1186_1754_1611_3_4*ALPHA_BBa_B0032*BBa_J23102;
	output = GFPmut3b;
end

function output = simf_doi_10_1186_1754_1611_3_4_figure3_6(x)
	ALPHA_BBa_B0032 = x(40);
	alpha_pLtetO_1 = x(19);
	scale_GFPmut3b_doi_10_1186_1754_1611_3_4 = x(42);
	pLtetO_1 = 10*alpha_pLtetO_1;
	GFPmut3b = scale_GFPmut3b_doi_10_1186_1754_1611_3_4*ALPHA_BBa_B0032*pLtetO_1;
	output = GFPmut3b;
end

function output = simf_doi_10_1186_1754_1611_3_4_figure3_7(x)
	ALPHA_BBa_B0032 = x(40);
	alpha_pLlacO_1 = x(12);
	scale_GFPmut3b_doi_10_1186_1754_1611_3_4 = x(42);
	pLlacO_1 = 10*alpha_pLlacO_1;
	GFPmut3b = scale_GFPmut3b_doi_10_1186_1754_1611_3_4*ALPHA_BBa_B0032*pLlacO_1;
	output = GFPmut3b;
end

function output = simf_doi_10_1093_nar_gku593_figure3D_1(x)
	ALPHA_BBa_B0030 = x(11);
	alpha_BBa_J23109 = x(47);
	scale_GFPmut3b_doi_10_1093_nar_gku593 = x(48);
	BBa_J23109 = 10*alpha_BBa_J23109;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku593*ALPHA_BBa_B0030*BBa_J23109;
	output = GFPmut3b;
end

function output = simf_doi_10_1093_nar_gku593_figure3D_2(x)
	ALPHA_BBa_B0030 = x(11);
	alpha_BBa_J23114 = x(49);
	scale_GFPmut3b_doi_10_1093_nar_gku593 = x(48);
	BBa_J23114 = 10*alpha_BBa_J23114;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku593*ALPHA_BBa_B0030*BBa_J23114;
	output = GFPmut3b;
end

function output = simf_doi_10_1093_nar_gku593_figure3D_3(x)
	ALPHA_BBa_B0030 = x(11);
	alpha_BBa_J23116 = x(43);
	scale_GFPmut3b_doi_10_1093_nar_gku593 = x(48);
	BBa_J23116 = 10*alpha_BBa_J23116;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku593*ALPHA_BBa_B0030*BBa_J23116;
	output = GFPmut3b;
end

function output = simf_doi_10_1093_nar_gku593_figure3D_4(x)
	ALPHA_BBa_B0030 = x(11);
	alpha_BBa_J23115 = x(50);
	scale_GFPmut3b_doi_10_1093_nar_gku593 = x(48);
	BBa_J23115 = 10*alpha_BBa_J23115;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku593*ALPHA_BBa_B0030*BBa_J23115;
	output = GFPmut3b;
end

function output = simf_doi_10_1093_nar_gku593_figure3D_5(x)
	ALPHA_BBa_B0030 = x(11);
	alpha_BBa_J23105 = x(51);
	scale_GFPmut3b_doi_10_1093_nar_gku593 = x(48);
	BBa_J23105 = 10*alpha_BBa_J23105;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku593*ALPHA_BBa_B0030*BBa_J23105;
	output = GFPmut3b;
end

function output = simf_doi_10_1093_nar_gku593_figure3D_6(x)
	ALPHA_BBa_B0030 = x(11);
	alpha_BBa_J23106 = x(52);
	scale_GFPmut3b_doi_10_1093_nar_gku593 = x(48);
	BBa_J23106 = 10*alpha_BBa_J23106;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku593*ALPHA_BBa_B0030*BBa_J23106;
	output = GFPmut3b;
end

function output = simf_doi_10_1093_nar_gku593_figureS6(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0030 = x(11);
	ALPHA_RBS_pC = x(2);
	K_L_arabinose = x(3);
	K_pBAD = x(4);
	alpha_pBAD = x(5);
	alpha_pC = x(6);
	beta_pBAD = x(7);
	n_L_arabinose = x(8);
	n_pBAD = x(9);
	scale_GFPmut3b_doi_10_1093_nar_gku593 = x(48);
	for i = 1:length(input)
		pC = 10*alpha_pC;
		L_arabinose = input(i);
		araC = ALPHA_RBS_pC*pC*K_L_arabinose^n_L_arabinose/(K_L_arabinose^n_L_arabinose + L_arabinose^n_L_arabinose);
		pBAD = 10*(beta_pBAD + (alpha_pBAD - beta_pBAD)*K_pBAD^n_pBAD/(K_pBAD^n_pBAD + araC^n_pBAD));
		GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku593*ALPHA_BBa_B0030*pBAD;
		output(i) = GFPmut3b;
	end
end

function output = simf_doi_10_1093_nar_gku1388_figure2A_1(x)
	ALPHA_BBa_B0030 = x(11);
	alpha_BBa_J23101 = x(41);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(53);
	BBa_J23101 = 10*alpha_BBa_J23101;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*BBa_J23101;
	output = GFPmut3b;
end

function output = simf_doi_10_1093_nar_gku1388_figure2A_2(x)
	ALPHA_BBa_B0030 = x(11);
	alpha_BBa_J23106 = x(52);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(53);
	BBa_J23106 = 10*alpha_BBa_J23106;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*BBa_J23106;
	output = GFPmut3b;
end

function output = simf_doi_10_1093_nar_gku1388_figure2A_3(x)
	ALPHA_BBa_B0030 = x(11);
	alpha_BBa_J23105 = x(51);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(53);
	BBa_J23105 = 10*alpha_BBa_J23105;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*BBa_J23105;
	output = GFPmut3b;
end

function output = simf_doi_10_1093_nar_gku1388_figure2A_4(x)
	ALPHA_BBa_B0030 = x(11);
	alpha_BBa_J23115 = x(50);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(53);
	BBa_J23115 = 10*alpha_BBa_J23115;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*BBa_J23115;
	output = GFPmut3b;
end

function output = simf_doi_10_1093_nar_gku1388_figure2A_5(x)
	ALPHA_BBa_B0030 = x(11);
	alpha_BBa_J23114 = x(49);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(53);
	BBa_J23114 = 10*alpha_BBa_J23114;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*BBa_J23114;
	output = GFPmut3b;
end

function output = simf_doi_10_1093_nar_gku1388_figure2A_6(x)
	ALPHA_BBa_B0030 = x(11);
	alpha_BBa_J23117 = x(54);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(53);
	BBa_J23117 = 10*alpha_BBa_J23117;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*BBa_J23117;
	output = GFPmut3b;
end

function output = simf_doi_10_1093_nar_gku1388_figure2B_1(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0030 = x(11);
	K_aTc = x(17);
	K_pTET_star_ = x(55);
	alpha_BBa_J23117 = x(54);
	alpha_pTET_star_ = x(56);
	beta_pTET_star_ = x(57);
	n_aTc = x(22);
	n_pTET_star_ = x(58);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(53);
	for i = 1:length(input)
		BBa_J23117 = 10*alpha_BBa_J23117;
		aTc = input(i);
		tetR = ALPHA_BBa_B0030*BBa_J23117*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pTET_star_ = 10*(beta_pTET_star_ + (alpha_pTET_star_ - beta_pTET_star_)*K_pTET_star_^n_pTET_star_/(K_pTET_star_^n_pTET_star_ + tetR^n_pTET_star_));
		GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*pTET_star_;
		output(i) = GFPmut3b;
	end
end

function output = simf_doi_10_1093_nar_gku1388_figure2B_2(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0030 = x(11);
	K_aTc = x(17);
	K_pTET_star_ = x(55);
	alpha_BBa_J23114 = x(49);
	alpha_pTET_star_ = x(56);
	beta_pTET_star_ = x(57);
	n_aTc = x(22);
	n_pTET_star_ = x(58);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(53);
	for i = 1:length(input)
		BBa_J23114 = 10*alpha_BBa_J23114;
		aTc = input(i);
		tetR = ALPHA_BBa_B0030*BBa_J23114*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pTET_star_ = 10*(beta_pTET_star_ + (alpha_pTET_star_ - beta_pTET_star_)*K_pTET_star_^n_pTET_star_/(K_pTET_star_^n_pTET_star_ + tetR^n_pTET_star_));
		GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*pTET_star_;
		output(i) = GFPmut3b;
	end
end

function output = simf_doi_10_1038_msb4100173_figure3_1(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_Unknown = x(59);
	ALPHA_SAL_RBS = x(60);
	K_Sal = x(61);
	K_pSAL = x(62);
	alpha_pCONST_UNKNOWN = x(63);
	alpha_pSAL = x(64);
	beta_pSAL = x(65);
	n_Sal = x(66);
	n_pSAL = x(67);
	scale_gfp_doi_10_1038_msb4100173 = x(68);
	for i = 1:length(input)
		pCONST_UNKNOWN = 15*alpha_pCONST_UNKNOWN;
		Sal = input(i);
		nahR = ALPHA_RBS_Unknown*pCONST_UNKNOWN;
		Sal_nahR = nahR*K_Sal^n_Sal/(K_Sal^n_Sal + Sal^n_Sal);
		pSAL = 15*(beta_pSAL + (alpha_pSAL - beta_pSAL)*Sal_nahR^n_pSAL/(K_pSAL^n_pSAL + Sal_nahR^n_pSAL));
		gfp = scale_gfp_doi_10_1038_msb4100173*ALPHA_SAL_RBS*pSAL;
		output(i) = gfp;
	end
end

function output = simf_doi_10_1038_msb4100173_figure3_2(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pC = x(2);
	ALPHA_parent_RBS = x(69);
	K_L_arabinose = x(3);
	K_pBAD = x(4);
	alpha_pBAD = x(5);
	alpha_pC = x(6);
	beta_pBAD = x(7);
	n_L_arabinose = x(8);
	n_pBAD = x(9);
	scale_gfp_doi_10_1038_msb4100173 = x(68);
	for i = 1:length(input)
		pC = 15*alpha_pC;
		L_arabinose = input(i);
		araC = ALPHA_RBS_pC*pC*K_L_arabinose^n_L_arabinose/(K_L_arabinose^n_L_arabinose + L_arabinose^n_L_arabinose);
		pBAD = 15*(beta_pBAD + (alpha_pBAD - beta_pBAD)*K_pBAD^n_pBAD/(K_pBAD^n_pBAD + araC^n_pBAD));
		gfp = scale_gfp_doi_10_1038_msb4100173*ALPHA_parent_RBS*pBAD;
		output(i) = gfp;
	end
end

function output = simf_doi_10_1534_genetics_112_147199_figure3(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pLAC = x(36);
	ALPHA_RBS_placI = x(25);
	ALPHA_T710 = x(70);
	K_IPTG = x(26);
	K_pLAC = x(27);
	K_pLlacO_1 = x(71);
	alpha_pLAC = x(28);
	alpha_pLlacO_1 = x(12);
	alpha_placI = x(37);
	beta_pLAC = x(30);
	beta_pLlacO_1 = x(72);
	n_IPTG = x(31);
	n_pLAC = x(32);
	n_pLlacO_1 = x(73);
	scale_cfp_doi_10_1534_genetics_112_147199 = x(74);
	for i = 1:length(input)
		placI = 1*alpha_placI;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placI*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pLAC = 1*(beta_pLAC + (alpha_pLAC - beta_pLAC)*K_pLAC^n_pLAC/(K_pLAC^n_pLAC + lacI^n_pLAC));
		lacZ = ALPHA_RBS_pLAC*pLAC;
		pLlacO_1 = 1*(beta_pLlacO_1 + (alpha_pLlacO_1 - beta_pLlacO_1)*K_pLlacO_1^n_pLlacO_1/(K_pLlacO_1^n_pLlacO_1 + lacI^n_pLlacO_1));
		cfp = scale_cfp_doi_10_1534_genetics_112_147199*ALPHA_T710*pLlacO_1;
		output(i) = cfp;
	end
end

function output = simf_doi_10_1016_j_ymben_2015_03_009_figure4_b(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_placI = x(25);
	ALPHA_g10 = x(75);
	K_IPTG = x(26);
	K_pTrc = x(76);
	alpha_pTrc = x(77);
	alpha_placIq = x(29);
	beta_pTrc = x(78);
	n_IPTG = x(31);
	n_pTrc = x(79);
	scale_gfp_LAA_doi_10_1016_j_ymben_2015_03_009 = x(80);
	for i = 1:length(input)
		placIq = 15*alpha_placIq;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placIq*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pTrc = 15*(beta_pTrc + (alpha_pTrc - beta_pTrc)*K_pTrc^n_pTrc/(K_pTrc^n_pTrc + lacI^n_pTrc));
		gfp_LAA = scale_gfp_LAA_doi_10_1016_j_ymben_2015_03_009*ALPHA_g10*pTrc;
		output(i) = gfp_LAA;
	end
end

function output = simf_doi_10_1016_j_ymben_2015_03_009_figure4_c(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pLUX = x(81);
	ALPHA_RBS_pluxR = x(82);
	K_AHL = x(83);
	K_pLUX = x(84);
	alpha_pLUX = x(85);
	alpha_pluxR = x(86);
	beta_pLUX = x(87);
	n_AHL = x(88);
	n_pLUX = x(89);
	scale_gfp_LAA_doi_10_1016_j_ymben_2015_03_009 = x(80);
	for i = 1:length(input)
		pluxR = 15*alpha_pluxR;
		AHL = input(i);
		luxR = ALPHA_RBS_pluxR*pluxR;
		AHL_luxR = luxR*K_AHL^n_AHL/(K_AHL^n_AHL + AHL^n_AHL);
		pLUX = 15*(beta_pLUX + (alpha_pLUX - beta_pLUX)*AHL_luxR^n_pLUX/(K_pLUX^n_pLUX + AHL_luxR^n_pLUX));
		gfp_LAA = scale_gfp_LAA_doi_10_1016_j_ymben_2015_03_009*ALPHA_RBS_pLUX*pLUX;
		output(i) = gfp_LAA;
	end
end

function output = simf_doi_10_1016_j_ces_2012_12_016_figure4(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0034 = x(90);
	K_C_3OC12HSL = x(91);
	K_pLUX = x(84);
	alpha_pLUX = x(85);
	alpha_pLtetO_1 = x(19);
	beta_pLUX = x(87);
	n_C_3OC12HSL = x(92);
	n_pLUX = x(89);
	scale_GFPmut3b_doi_10_1016_j_ces_2012_12_016 = x(93);
	for i = 1:length(input)
		pLtetO_1 = 200*alpha_pLtetO_1;
		lasR = ALPHA_BBa_B0034*pLtetO_1;
		C_3OC12HSL = input(i);
		C_3OC12HSL_lasR = lasR*K_C_3OC12HSL^n_C_3OC12HSL/(K_C_3OC12HSL^n_C_3OC12HSL + C_3OC12HSL^n_C_3OC12HSL);
		pLUX = 200*(beta_pLUX + (alpha_pLUX - beta_pLUX)*C_3OC12HSL_lasR^n_pLUX/(K_pLUX^n_pLUX + C_3OC12HSL_lasR^n_pLUX));
		GFPmut3b = scale_GFPmut3b_doi_10_1016_j_ces_2012_12_016*ALPHA_BBa_B0034*pLUX;
		output(i) = GFPmut3b;
	end
end

function output = simf_doi_10_1038_msb_2009_75_figure3_1(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_UNKNOWN = x(94);
	ALPHA_RBS_pLAC = x(36);
	ALPHA_RBS_placI = x(25);
	K_IPTG = x(26);
	K_pLAC = x(27);
	alpha_pLAC = x(28);
	alpha_placI = x(37);
	beta_pLAC = x(30);
	n_IPTG = x(31);
	n_pLAC = x(32);
	scale_gfp_LVA_doi_10_1038_msb_2009_75 = x(95);
	for i = 1:length(input)
		placI = 1*alpha_placI;
		pLAC = 200*alpha_pLAC;
		gfp_LVA = scale_gfp_LVA_doi_10_1038_msb_2009_75*ALPHA_RBS_UNKNOWN*pLAC;
		output(i) = gfp_LVA;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placI*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pLAC = 1*(beta_pLAC + (alpha_pLAC - beta_pLAC)*K_pLAC^n_pLAC/(K_pLAC^n_pLAC + lacI^n_pLAC));
		lacZ = ALPHA_RBS_pLAC*pLAC;
	end
end

function output = simf_doi_10_1016_j_bios_2012_08_011_figure2d(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0032 = x(40);
	ALPHA_BBa_B0034 = x(90);
	K_AHL = x(83);
	K_pLUX = x(84);
	alpha_pLUX = x(85);
	alpha_pTET = x(96);
	beta_pLUX = x(87);
	n_AHL = x(88);
	n_pLUX = x(89);
	scale_GFPmut3b_doi_10_1016_j_bios_2012_08_011 = x(97);
	for i = 1:length(input)
		pTET = 10*alpha_pTET;
		luxR = ALPHA_BBa_B0034*pTET;
		AHL = input(i);
		AHL_luxR = luxR*K_AHL^n_AHL/(K_AHL^n_AHL + AHL^n_AHL);
		pLUX = 10*(beta_pLUX + (alpha_pLUX - beta_pLUX)*AHL_luxR^n_pLUX/(K_pLUX^n_pLUX + AHL_luxR^n_pLUX));
		GFPmut3b = scale_GFPmut3b_doi_10_1016_j_bios_2012_08_011*ALPHA_BBa_B0032*pLUX;
		output(i) = GFPmut3b;
	end
end

function output = simf_doi_10_1093_nar_gku964_figure3A_1(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0032 = x(40);
	ALPHA_BBa_B0034 = x(90);
	K_AHL = x(83);
	K_pLUX = x(84);
	alpha_pLUX = x(85);
	alpha_pLtetO_1 = x(19);
	beta_pLUX = x(87);
	n_AHL = x(88);
	n_pLUX = x(89);
	scale_mRFP1_doi_10_1093_nar_gku964 = x(98);
	for i = 1:length(input)
		pLtetO_1 = 200*alpha_pLtetO_1;
		luxR = ALPHA_BBa_B0034*pLtetO_1;
		GFPmut3b = ALPHA_BBa_B0034*pLtetO_1;
		AHL = input(i);
		AHL_luxR = luxR*K_AHL^n_AHL/(K_AHL^n_AHL + AHL^n_AHL);
		pLUX = 200*(beta_pLUX + (alpha_pLUX - beta_pLUX)*AHL_luxR^n_pLUX/(K_pLUX^n_pLUX + AHL_luxR^n_pLUX));
		mRFP1 = scale_mRFP1_doi_10_1093_nar_gku964*ALPHA_BBa_B0032*pLUX;
		output(i) = mRFP1;
	end
end

function output = simf_doi_10_1093_nar_gku964_figure3A_2(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0032 = x(40);
	K_AHL = x(83);
	K_pLUX = x(84);
	alpha_pLUX = x(85);
	alpha_pLtetO_1 = x(19);
	beta_pLUX = x(87);
	n_AHL = x(88);
	n_pLUX = x(89);
	scale_mRFP1_doi_10_1093_nar_gku964 = x(98);
	for i = 1:length(input)
		pLtetO_1 = 200*alpha_pLtetO_1;
		luxR = ALPHA_BBa_B0032*pLtetO_1;
		GFPmut3b = ALPHA_BBa_B0032*pLtetO_1;
		AHL = input(i);
		AHL_luxR = luxR*K_AHL^n_AHL/(K_AHL^n_AHL + AHL^n_AHL);
		pLUX = 200*(beta_pLUX + (alpha_pLUX - beta_pLUX)*AHL_luxR^n_pLUX/(K_pLUX^n_pLUX + AHL_luxR^n_pLUX));
		mRFP1 = scale_mRFP1_doi_10_1093_nar_gku964*ALPHA_BBa_B0032*pLUX;
		output(i) = mRFP1;
	end
end

function output = simf_doi_10_1093_nar_gku964_figure3A_3(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0032 = x(40);
	ALPHA_BBa_B0033 = x(99);
	K_AHL = x(83);
	K_pLUX = x(84);
	alpha_pLUX = x(85);
	alpha_pLtetO_1 = x(19);
	beta_pLUX = x(87);
	n_AHL = x(88);
	n_pLUX = x(89);
	scale_mRFP1_doi_10_1093_nar_gku964 = x(98);
	for i = 1:length(input)
		pLtetO_1 = 200*alpha_pLtetO_1;
		luxR = ALPHA_BBa_B0033*pLtetO_1;
		GFPmut3b = ALPHA_BBa_B0033*pLtetO_1;
		AHL = input(i);
		AHL_luxR = luxR*K_AHL^n_AHL/(K_AHL^n_AHL + AHL^n_AHL);
		pLUX = 200*(beta_pLUX + (alpha_pLUX - beta_pLUX)*AHL_luxR^n_pLUX/(K_pLUX^n_pLUX + AHL_luxR^n_pLUX));
		mRFP1 = scale_mRFP1_doi_10_1093_nar_gku964*ALPHA_BBa_B0032*pLUX;
		output(i) = mRFP1;
	end
end

function output = simf_doi_10_1021_bp010141k_figure3_1(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_Unknown = x(59);
	ALPHA_RBS_placI = x(25);
	K_IPTG = x(26);
	K_pTrc = x(76);
	alpha_pTrc = x(77);
	alpha_placIq = x(29);
	beta_pTrc = x(78);
	n_IPTG = x(31);
	n_pTrc = x(79);
	scale_gfp_uv_doi_10_1021_bp010141k = x(100);
	for i = 1:length(input)
		placIq = 20*alpha_placIq;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placIq*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pTrc = 20*(beta_pTrc + (alpha_pTrc - beta_pTrc)*K_pTrc^n_pTrc/(K_pTrc^n_pTrc + lacI^n_pTrc));
		gfp_uv = scale_gfp_uv_doi_10_1021_bp010141k*ALPHA_RBS_Unknown*pTrc;
		output(i) = gfp_uv;
	end
end

function output = simf_doi_10_1021_sb300055e_figure4a_1(x)
	ALPHA_RBS_Glk = x(101);
	alpha_BBa_J23114 = x(49);
	scale_Glk_doi_10_1021_sb300055e = x(102);
	BBa_J23114 = 1*alpha_BBa_J23114;
	Glk = scale_Glk_doi_10_1021_sb300055e*ALPHA_RBS_Glk*BBa_J23114;
	output = Glk;
end

function output = simf_doi_10_1021_sb300055e_figure4a_2(x)
	ALPHA_RBS_Glk = x(101);
	alpha_BBa_J23117 = x(54);
	scale_Glk_doi_10_1021_sb300055e = x(102);
	BBa_J23117 = 1*alpha_BBa_J23117;
	Glk = scale_Glk_doi_10_1021_sb300055e*ALPHA_RBS_Glk*BBa_J23117;
	output = Glk;
end

function output = simf_doi_10_1021_sb300055e_figure4a_3(x)
	ALPHA_RBS_Glk = x(101);
	alpha_BBa_J23109 = x(47);
	scale_Glk_doi_10_1021_sb300055e = x(102);
	BBa_J23109 = 1*alpha_BBa_J23109;
	Glk = scale_Glk_doi_10_1021_sb300055e*ALPHA_RBS_Glk*BBa_J23109;
	output = Glk;
end

function output = simf_doi_10_1038_msb4100081_figure3b(x, input)
	output = zeros(1,length(input));
	ALPHA_RBSII = x(15);
	K_aTc = x(17);
	K_pLtetO_1 = x(18);
	alpha_pLtetO_1 = x(19);
	alpha_plac_ara_1 = x(103);
	beta_pLtetO_1 = x(21);
	n_aTc = x(22);
	n_pLtetO_1 = x(23);
	scale_gfpmut3_1_doi_10_1038_msb4100081 = x(104);
	for i = 1:length(input)
		plac_ara_1 = 5*alpha_plac_ara_1;
		aTc = input(i);
		tetR = ALPHA_RBSII*plac_ara_1*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pLtetO_1 = 15*(beta_pLtetO_1 + (alpha_pLtetO_1 - beta_pLtetO_1)*K_pLtetO_1^n_pLtetO_1/(K_pLtetO_1^n_pLtetO_1 + tetR^n_pLtetO_1));
		gfpmut3_1 = scale_gfpmut3_1_doi_10_1038_msb4100081*ALPHA_RBSII*pLtetO_1;
		output(i) = gfpmut3_1;
	end
end

function output = simf_doi_10_1186_1754_1611_6_9_figure2_b(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_UTR = x(105);
	ALPHA_RBS_pLAMBDA = x(106);
	ALPHA_RBS_placI = x(25);
	K_IPTG = x(26);
	K_pLAMBDA = x(107);
	K_pLlacO_1 = x(71);
	alpha_pLAMBDA = x(108);
	alpha_pLlacO_1 = x(12);
	alpha_placI = x(37);
	beta_pLAMBDA = x(109);
	beta_pLlacO_1 = x(72);
	n_IPTG = x(31);
	n_pLAMBDA = x(110);
	n_pLlacO_1 = x(73);
	scale_gfp_doi_10_1186_1754_1611_6_9 = x(111);
	for i = 1:length(input)
		placI = 15*alpha_placI;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placI*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pLlacO_1 = 15*(beta_pLlacO_1 + (alpha_pLlacO_1 - beta_pLlacO_1)*K_pLlacO_1^n_pLlacO_1/(K_pLlacO_1^n_pLlacO_1 + lacI^n_pLlacO_1));
		cI = ALPHA_RBS_UTR*pLlacO_1;
		pLAMBDA = 10*(beta_pLAMBDA + (alpha_pLAMBDA - beta_pLAMBDA)*K_pLAMBDA^n_pLAMBDA/(K_pLAMBDA^n_pLAMBDA + cI^n_pLAMBDA));
		gfp = scale_gfp_doi_10_1186_1754_1611_6_9*ALPHA_RBS_pLAMBDA*pLAMBDA;
		output(i) = gfp;
	end
end

function output = simf_doi_10_1186_1752_0509_5_111_figure3_a(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pBAD = x(1);
	ALPHA_RBS_pC = x(2);
	K_L_arabinose = x(3);
	K_pBAD = x(4);
	alpha_pBAD = x(5);
	alpha_pC = x(6);
	beta_pBAD = x(7);
	n_L_arabinose = x(8);
	n_pBAD = x(9);
	scale_gfpmut2_doi_10_1186_1752_0509_5_111 = x(112);
	for i = 1:length(input)
		pC = 1*alpha_pC;
		L_arabinose = input(i);
		araC = ALPHA_RBS_pC*pC*K_L_arabinose^n_L_arabinose/(K_L_arabinose^n_L_arabinose + L_arabinose^n_L_arabinose);
		pBAD = 5*(beta_pBAD + (alpha_pBAD - beta_pBAD)*K_pBAD^n_pBAD/(K_pBAD^n_pBAD + araC^n_pBAD));
		gfpmut2 = scale_gfpmut2_doi_10_1186_1752_0509_5_111*ALPHA_RBS_pBAD*pBAD;
		output(i) = gfpmut2;
	end
end

function output = simf_doi_10_1186_1752_0509_5_111_figure3_b(x, input)
	output = zeros(1,length(input));
	ALPHA_RBSII = x(15);
	ALPHA_RBS_pBAD = x(1);
	ALPHA_RBS_pN25 = x(16);
	K_L_arabinose = x(3);
	K_pBAD = x(4);
	K_pLtetO_1 = x(18);
	alpha_pBAD = x(5);
	alpha_pLtetO_1 = x(19);
	alpha_pN25 = x(20);
	beta_pBAD = x(7);
	beta_pLtetO_1 = x(21);
	n_L_arabinose = x(8);
	n_pBAD = x(9);
	n_pLtetO_1 = x(23);
	scale_gfpmut2_doi_10_1186_1752_0509_5_111 = x(112);
	for i = 1:length(input)
		pN25 = 1*alpha_pN25;
		tetR = ALPHA_RBS_pN25*pN25;
		L_arabinose = input(i);
		pLtetO_1 = 15*(beta_pLtetO_1 + (alpha_pLtetO_1 - beta_pLtetO_1)*K_pLtetO_1^n_pLtetO_1/(K_pLtetO_1^n_pLtetO_1 + tetR^n_pLtetO_1));
		araC = ALPHA_RBSII*pLtetO_1*K_L_arabinose^n_L_arabinose/(K_L_arabinose^n_L_arabinose + L_arabinose^n_L_arabinose);
		pBAD = 5*(beta_pBAD + (alpha_pBAD - beta_pBAD)*K_pBAD^n_pBAD/(K_pBAD^n_pBAD + araC^n_pBAD));
		gfpmut2 = scale_gfpmut2_doi_10_1186_1752_0509_5_111*ALPHA_RBS_pBAD*pBAD;
		output(i) = gfpmut2;
	end
end

function output = simf_dx_doi_org_10_1021_sb400152n_figure2(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0034 = x(90);
	K_AHL = x(83);
	K_pLUX = x(84);
	alpha_pCAT = x(113);
	alpha_pLUX = x(85);
	beta_pLUX = x(87);
	n_AHL = x(88);
	n_pLUX = x(89);
	scale_eYFP_dx_doi_org_10_1021_sb400152n = x(114);
	for i = 1:length(input)
		pCAT = 10*alpha_pCAT;
		luxR = ALPHA_BBa_B0034*pCAT;
		AHL = input(i);
		AHL_luxR = luxR*K_AHL^n_AHL/(K_AHL^n_AHL + AHL^n_AHL);
		pLUX = 10*(beta_pLUX + (alpha_pLUX - beta_pLUX)*AHL_luxR^n_pLUX/(K_pLUX^n_pLUX + AHL_luxR^n_pLUX));
		eYFP = scale_eYFP_dx_doi_org_10_1021_sb400152n*ALPHA_BBa_B0034*pLUX;
		output(i) = eYFP;
	end
end

function output = simf_doi_10_1073_pnas_252535999_figure3(x, input)
	output = zeros(1,length(input));
	ALPHA_RBSII = x(15);
	ALPHA_RBS_placI = x(25);
	K_IPTG = x(26);
	K_pLAC = x(27);
	K_pLAMBDA = x(107);
	alpha_pLAC = x(28);
	alpha_pLAMBDA = x(108);
	alpha_placIq = x(29);
	beta_pLAC = x(30);
	beta_pLAMBDA = x(109);
	n_IPTG = x(31);
	n_pLAC = x(32);
	n_pLAMBDA = x(110);
	scale_eYFP_doi_10_1073_pnas_252535999 = x(115);
	for i = 1:length(input)
		placIq = 10*alpha_placIq;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placIq*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pLAC = 10*(beta_pLAC + (alpha_pLAC - beta_pLAC)*K_pLAC^n_pLAC/(K_pLAC^n_pLAC + lacI^n_pLAC));
		cI = ALPHA_RBSII*pLAC;
		eCFP = ALPHA_RBSII*pLAC;
		pLAMBDA = 15*(beta_pLAMBDA + (alpha_pLAMBDA - beta_pLAMBDA)*K_pLAMBDA^n_pLAMBDA/(K_pLAMBDA^n_pLAMBDA + cI^n_pLAMBDA));
		eYFP = scale_eYFP_doi_10_1073_pnas_252535999*ALPHA_RBSII*pLAMBDA;
		output(i) = eYFP;
	end
end

function output = simf_doi_10_1038_35002131_figure5a(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_D = x(116);
	ALPHA_RBS_E = x(117);
	K_IPTG = x(26);
	K_pTrc_2 = x(118);
	alpha_pLs1con = x(119);
	alpha_pTrc_2 = x(120);
	beta_pTrc_2 = x(121);
	n_IPTG = x(31);
	n_pTrc_2 = x(122);
	scale_GFPmut3b_doi_10_1038_35002131 = x(123);
	for i = 1:length(input)
		pLs1con = 15*alpha_pLs1con;
		IPTG = input(i);
		lacI = ALPHA_RBS_D*pLs1con*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pTrc_2 = 15*(beta_pTrc_2 + (alpha_pTrc_2 - beta_pTrc_2)*K_pTrc_2^n_pTrc_2/(K_pTrc_2^n_pTrc_2 + lacI^n_pTrc_2));
		GFPmut3b = scale_GFPmut3b_doi_10_1038_35002131*ALPHA_RBS_E*pTrc_2;
		output(i) = GFPmut3b;
	end
end

function output = simf_doi_10_1016_j_ymben_2012_08_006_figure7(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0034 = x(90);
	ALPHA_RBS_native = x(124);
	ALPHA_RBS_tetR = x(125);
	K_aTc = x(17);
	K_pJ23114_lacO = x(126);
	K_pTET = x(127);
	alpha_pC_tetR_UNKNOWN = x(128);
	alpha_pJ23114_lacO = x(129);
	alpha_pTET = x(96);
	beta_pJ23114_lacO = x(130);
	beta_pTET = x(131);
	n_aTc = x(22);
	n_pJ23114_lacO = x(132);
	n_pTET = x(133);
	scale_glk_doi_10_1016_j_ymben_2012_08_006 = x(134);
	for i = 1:length(input)
		pC_tetR_UNKNOWN = 1*alpha_pC_tetR_UNKNOWN;
		aTc = input(i);
		tetR = ALPHA_RBS_tetR*pC_tetR_UNKNOWN*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pTET = 1*(beta_pTET + (alpha_pTET - beta_pTET)*K_pTET^n_pTET/(K_pTET^n_pTET + tetR^n_pTET));
		lacI = ALPHA_BBa_B0034*pTET;
		pJ23114_lacO = 1*(beta_pJ23114_lacO + (alpha_pJ23114_lacO - beta_pJ23114_lacO)*K_pJ23114_lacO^n_pJ23114_lacO/(K_pJ23114_lacO^n_pJ23114_lacO + lacI^n_pJ23114_lacO));
		glk = scale_glk_doi_10_1016_j_ymben_2012_08_006*ALPHA_RBS_native*pJ23114_lacO;
		output(i) = glk;
	end
end

function output = simf_doi_10_1006_plas_2000_1477_figure5(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pLAC = x(36);
	ALPHA_RBS_placI = x(25);
	K_IPTG = x(26);
	K_pLAC = x(27);
	alpha_pLAC = x(28);
	alpha_placI = x(37);
	beta_pLAC = x(30);
	n_IPTG = x(31);
	n_pLAC = x(32);
	scale_lacZ_doi_10_1006_plas_2000_1477 = x(135);
	for i = 1:length(input)
		placI = 1*alpha_placI;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placI*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pLAC = 15*(beta_pLAC + (alpha_pLAC - beta_pLAC)*K_pLAC^n_pLAC/(K_pLAC^n_pLAC + lacI^n_pLAC));
		lacZ = scale_lacZ_doi_10_1006_plas_2000_1477*ALPHA_RBS_pLAC*pLAC;
		output(i) = lacZ;
	end
end

function output = simf_doi_10_1038_nbt1413_figure3(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0032 = x(40);
	ALPHA_BBa_B0034 = x(90);
	K_AHL = x(83);
	K_pLUX = x(84);
	alpha_pLUX = x(85);
	alpha_pLtetO_1 = x(19);
	beta_pLUX = x(87);
	n_AHL = x(88);
	n_pLUX = x(89);
	scale_GFPmut3b_doi_10_1038_nbt1413 = x(136);
	for i = 1:length(input)
		pLtetO_1 = 10*alpha_pLtetO_1;
		luxR = ALPHA_BBa_B0034*pLtetO_1;
		AHL = input(i);
		AHL_luxR = luxR*K_AHL^n_AHL/(K_AHL^n_AHL + AHL^n_AHL);
		pLUX = 10*(beta_pLUX + (alpha_pLUX - beta_pLUX)*AHL_luxR^n_pLUX/(K_pLUX^n_pLUX + AHL_luxR^n_pLUX));
		GFPmut3b = scale_GFPmut3b_doi_10_1038_nbt1413*ALPHA_BBa_B0032*pLUX;
		output(i) = GFPmut3b;
	end
end

function output = simf_doi_10_1093_nar_gkt052_figure3b(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0034 = x(90);
	K_C_3OC12HSL = x(91);
	K_pLasRV = x(137);
	alpha_BBa_J23100 = x(138);
	alpha_pLasRV = x(139);
	beta_pLasRV = x(140);
	n_C_3OC12HSL = x(92);
	n_pLasRV = x(141);
	scale_GFPmut3b_doi_10_1093_nar_gkt052 = x(142);
	for i = 1:length(input)
		BBa_J23100 = 200*alpha_BBa_J23100;
		lasR = ALPHA_BBa_B0034*BBa_J23100;
		C_3OC12HSL = input(i);
		C_3OC12HSL_lasR = lasR*K_C_3OC12HSL^n_C_3OC12HSL/(K_C_3OC12HSL^n_C_3OC12HSL + C_3OC12HSL^n_C_3OC12HSL);
		pLasRV = 200*(beta_pLasRV + (alpha_pLasRV - beta_pLasRV)*C_3OC12HSL_lasR^n_pLasRV/(K_pLasRV^n_pLasRV + C_3OC12HSL_lasR^n_pLasRV));
		GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gkt052*ALPHA_BBa_B0034*pLasRV;
		output(i) = GFPmut3b;
	end
end

function output = simf_doi_10_1093_nar_gkt052_figure3c(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0034 = x(90);
	K_C_3OC12HSL = x(91);
	K_pLasR3 = x(143);
	alpha_BBa_J23100 = x(138);
	alpha_pLasR3 = x(144);
	beta_pLasR3 = x(145);
	n_C_3OC12HSL = x(92);
	n_pLasR3 = x(146);
	scale_GFPmut3b_doi_10_1093_nar_gkt052 = x(142);
	for i = 1:length(input)
		BBa_J23100 = 200*alpha_BBa_J23100;
		lasR = ALPHA_BBa_B0034*BBa_J23100;
		C_3OC12HSL = input(i);
		C_3OC12HSL_lasR = lasR*K_C_3OC12HSL^n_C_3OC12HSL/(K_C_3OC12HSL^n_C_3OC12HSL + C_3OC12HSL^n_C_3OC12HSL);
		pLasR3 = 200*(beta_pLasR3 + (alpha_pLasR3 - beta_pLasR3)*C_3OC12HSL_lasR^n_pLasR3/(K_pLasR3^n_pLasR3 + C_3OC12HSL_lasR^n_pLasR3));
		GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gkt052*ALPHA_BBa_B0034*pLasR3;
		output(i) = GFPmut3b;
	end
end

function output = simf_doi_10_1093_nar_gkt052_figure3d(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0034 = x(90);
	K_C_3OC12HSL = x(91);
	K_pLasRL = x(147);
	alpha_BBa_J23100 = x(138);
	alpha_pLasRL = x(148);
	beta_pLasRL = x(149);
	n_C_3OC12HSL = x(92);
	n_pLasRL = x(150);
	scale_GFPmut3b_doi_10_1093_nar_gkt052 = x(142);
	for i = 1:length(input)
		BBa_J23100 = 200*alpha_BBa_J23100;
		lasR = ALPHA_BBa_B0034*BBa_J23100;
		C_3OC12HSL = input(i);
		C_3OC12HSL_lasR = lasR*K_C_3OC12HSL^n_C_3OC12HSL/(K_C_3OC12HSL^n_C_3OC12HSL + C_3OC12HSL^n_C_3OC12HSL);
		pLasRL = 200*(beta_pLasRL + (alpha_pLasRL - beta_pLasRL)*C_3OC12HSL_lasR^n_pLasRL/(K_pLasRL^n_pLasRL + C_3OC12HSL_lasR^n_pLasRL));
		GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gkt052*ALPHA_BBa_B0034*pLasRL;
		output(i) = GFPmut3b;
	end
end

function output = simf_doi_10_1038_nmeth_2205_figureS1_1(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_KDL = x(151);
	K_Mg2_plus__plus_ = x(152);
	K_pMgrB = x(153);
	alpha_pMgrB = x(154);
	beta_pMgrB = x(155);
	n_Mg2_plus__plus_ = x(156);
	n_pMgrB = x(157);
	scale_gfp_doi_10_1038_nmeth_2205 = x(158);
	for i = 1:length(input)
		Mg2_plus__plus_ = input(i);
		pMgrB = 15*(beta_pMgrB + (alpha_pMgrB - beta_pMgrB)*K_pMgrB^n_pMgrB/(K_pMgrB^n_pMgrB + Mg2_plus__plus_^n_pMgrB));
		gfp = scale_gfp_doi_10_1038_nmeth_2205*ALPHA_RBS_KDL*pMgrB;
		output(i) = gfp;
	end
end

function output = simf_doi_10_1038_nmeth_2205_figureS1_2(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_KDL = x(151);
	ALPHA_RBS_pC = x(2);
	K_L_arabinose = x(3);
	K_pBAD = x(4);
	alpha_pBAD = x(5);
	alpha_pC = x(6);
	beta_pBAD = x(7);
	n_L_arabinose = x(8);
	n_pBAD = x(9);
	scale_gfp_doi_10_1038_nmeth_2205 = x(158);
	for i = 1:length(input)
		pC = 1*alpha_pC;
		L_arabinose = input(i);
		araC = ALPHA_RBS_pC*pC*K_L_arabinose^n_L_arabinose/(K_L_arabinose^n_L_arabinose + L_arabinose^n_L_arabinose);
		pBAD = 15*(beta_pBAD + (alpha_pBAD - beta_pBAD)*K_pBAD^n_pBAD/(K_pBAD^n_pBAD + araC^n_pBAD));
		gfp = scale_gfp_doi_10_1038_nmeth_2205*ALPHA_RBS_KDL*pBAD;
		output(i) = gfp;
	end
end

function output = simf_doi_10_1038_nmeth_2205_figureS1_3(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_KDL = x(151);
	ALPHA_RBS_pN25 = x(16);
	K_aTc = x(17);
	K_pLtetO_1 = x(18);
	alpha_pLtetO_1 = x(19);
	alpha_pN25 = x(20);
	beta_pLtetO_1 = x(21);
	n_aTc = x(22);
	n_pLtetO_1 = x(23);
	scale_gfp_doi_10_1038_nmeth_2205 = x(158);
	for i = 1:length(input)
		pN25 = 1*alpha_pN25;
		aTc = input(i);
		tetR = ALPHA_RBS_pN25*pN25*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pLtetO_1 = 15*(beta_pLtetO_1 + (alpha_pLtetO_1 - beta_pLtetO_1)*K_pLtetO_1^n_pLtetO_1/(K_pLtetO_1^n_pLtetO_1 + tetR^n_pLtetO_1));
		gfp = scale_gfp_doi_10_1038_nmeth_2205*ALPHA_RBS_KDL*pLtetO_1;
		output(i) = gfp;
	end
end

function output = simf_doi_10_1038_nmeth_2205_figureS1_41(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_KDL = x(151);
	ALPHA_RBS_placI = x(25);
	K_IPTG = x(26);
	K_pLlacO_1 = x(71);
	alpha_pLlacO_1 = x(12);
	alpha_placIq = x(29);
	beta_pLlacO_1 = x(72);
	n_IPTG = x(31);
	n_pLlacO_1 = x(73);
	scale_gfp_doi_10_1038_nmeth_2205 = x(158);
	for i = 1:length(input)
		placIq = 1*alpha_placIq;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placIq*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pLlacO_1 = 15*(beta_pLlacO_1 + (alpha_pLlacO_1 - beta_pLlacO_1)*K_pLlacO_1^n_pLlacO_1/(K_pLlacO_1^n_pLlacO_1 + lacI^n_pLlacO_1));
		gfp = scale_gfp_doi_10_1038_nmeth_2205*ALPHA_RBS_KDL*pLlacO_1;
		output(i) = gfp;
	end
end

function output = simf_doi_10_1038_nmeth_2205_figureS1_42(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_KDL = x(151);
	ALPHA_RBS_placI = x(25);
	K_IPTG = x(26);
	K_pTrc_2 = x(118);
	alpha_pTrc_2 = x(120);
	alpha_placIq = x(29);
	beta_pTrc_2 = x(121);
	n_IPTG = x(31);
	n_pTrc_2 = x(122);
	scale_gfp_doi_10_1038_nmeth_2205 = x(158);
	for i = 1:length(input)
		placIq = 1*alpha_placIq;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placIq*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pTrc_2 = 15*(beta_pTrc_2 + (alpha_pTrc_2 - beta_pTrc_2)*K_pTrc_2^n_pTrc_2/(K_pTrc_2^n_pTrc_2 + lacI^n_pTrc_2));
		gfp = scale_gfp_doi_10_1038_nmeth_2205*ALPHA_RBS_KDL*pTrc_2;
		output(i) = gfp;
	end
end

function output = simf_doi_10_1038_nchembio_1411_figure4_1(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_placI = x(25);
	ALPHA_rbs1 = x(159);
	K_IPTG = x(26);
	K_pTAC = x(160);
	alpha_pTAC = x(161);
	alpha_placI = x(37);
	beta_pTAC = x(162);
	n_IPTG = x(31);
	n_pTAC = x(163);
	scale_eYFP_doi_10_1038_nchembio_1411 = x(164);
	for i = 1:length(input)
		placI = 10*alpha_placI;
		luxR = ALPHA_rbs1*placI;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placI*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pTAC = 10*(beta_pTAC + (alpha_pTAC - beta_pTAC)*K_pTAC^n_pTAC/(K_pTAC^n_pTAC + lacI^n_pTAC));
		eYFP = scale_eYFP_doi_10_1038_nchembio_1411*ALPHA_rbs1*pTAC;
		output(i) = eYFP;
	end
end

function output = simf_doi_10_1038_nchembio_1411_figure4_2(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0034 = x(90);
	ALPHA_RBS_placI = x(25);
	ALPHA_rbs1 = x(159);
	K_IPTG = x(26);
	K_pTAC = x(160);
	alpha_J23119_tetO = x(165);
	alpha_pTAC = x(161);
	alpha_placI = x(37);
	beta_pTAC = x(162);
	n_IPTG = x(31);
	n_pTAC = x(163);
	scale_eYFP_doi_10_1038_nchembio_1411 = x(164);
	for i = 1:length(input)
		J23119_tetO = 10*alpha_J23119_tetO;
		placI = 10*alpha_placI;
		luxR = ALPHA_rbs1*placI;
		IPTG = input(i);
		eYFP = scale_eYFP_doi_10_1038_nchembio_1411*ALPHA_BBa_B0034*J23119_tetO;
		output(i) = eYFP;
		lacI = ALPHA_RBS_placI*placI*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pTAC = 10*(beta_pTAC + (alpha_pTAC - beta_pTAC)*K_pTAC^n_pTAC/(K_pTAC^n_pTAC + lacI^n_pTAC));
		tetR = ALPHA_rbs1*pTAC;
	end
end

function output = simf_doi_10_1093_nar_gks583_figure3(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_PtetA = x(166);
	ALPHA_RBS_pN25 = x(16);
	K_PtetA = x(167);
	K_aTc = x(17);
	alpha_PtetA = x(168);
	alpha_pN25 = x(20);
	beta_PtetA = x(169);
	n_PtetA = x(170);
	n_aTc = x(22);
	scale_mRFP1_doi_10_1093_nar_gks583 = x(171);
	for i = 1:length(input)
		pN25 = 1*alpha_pN25;
		aTc = input(i);
		tetR = ALPHA_RBS_pN25*pN25*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		PtetA = 1*(beta_PtetA + (alpha_PtetA - beta_PtetA)*K_PtetA^n_PtetA/(K_PtetA^n_PtetA + tetR^n_PtetA));
		mRFP1 = scale_mRFP1_doi_10_1093_nar_gks583*ALPHA_RBS_PtetA*PtetA;
		output(i) = mRFP1;
	end
end

function output = simf_doi_10_1002_bit_20371_figure3(x, input)
	output = zeros(1,length(input));
	ALPHA_RBSII = x(15);
	ALPHA_RBS_placI = x(25);
	K_C4HSL = x(172);
	K_pLAMBDA_R_O12 = x(173);
	K_qsc = x(174);
	alpha_pLAMBDA_R_O12 = x(175);
	alpha_placIq = x(29);
	alpha_qsc = x(176);
	beta_pLAMBDA_R_O12 = x(177);
	beta_qsc = x(178);
	n_C4HSL = x(179);
	n_pLAMBDA_R_O12 = x(180);
	n_qsc = x(181);
	scale_eYFP_doi_10_1002_bit_20371 = x(182);
	for i = 1:length(input)
		placIq = 10*alpha_placIq;
		C4HSL = input(i);
		rhlR = ALPHA_RBS_placI*placIq;
		C4HSL_rhlR = rhlR*K_C4HSL^n_C4HSL/(K_C4HSL^n_C4HSL + C4HSL^n_C4HSL);
		qsc = 10*(beta_qsc + (alpha_qsc - beta_qsc)*C4HSL_rhlR^n_qsc/(K_qsc^n_qsc + C4HSL_rhlR^n_qsc));
		cI_LVA = ALPHA_RBSII*qsc;
		pLAMBDA_R_O12 = 15*(beta_pLAMBDA_R_O12 + (alpha_pLAMBDA_R_O12 - beta_pLAMBDA_R_O12)*K_pLAMBDA_R_O12^n_pLAMBDA_R_O12/(K_pLAMBDA_R_O12^n_pLAMBDA_R_O12 + cI_LVA^n_pLAMBDA_R_O12));
		eYFP = scale_eYFP_doi_10_1002_bit_20371*ALPHA_RBSII*pLAMBDA_R_O12;
		output(i) = eYFP;
	end
end

function output = simf_doi_10_1073_pnas_1314114111_figure5_b(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pTAC = x(183);
	ALPHA_RBS_placI = x(25);
	K_IPTG = x(26);
	K_pTAC = x(160);
	alpha_pTAC = x(161);
	alpha_placIq = x(29);
	beta_pTAC = x(162);
	n_IPTG = x(31);
	n_pTAC = x(163);
	scale_GFPmut3b_doi_10_1073_pnas_1314114111 = x(184);
	for i = 1:length(input)
		placIq = 20*alpha_placIq;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placIq*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pTAC = 20*(beta_pTAC + (alpha_pTAC - beta_pTAC)*K_pTAC^n_pTAC/(K_pTAC^n_pTAC + lacI^n_pTAC));
		GFPmut3b = scale_GFPmut3b_doi_10_1073_pnas_1314114111*ALPHA_RBS_pTAC*pTAC;
		output(i) = GFPmut3b;
	end
end

function output = simf_doi_10_1073_pnas_1314114111_figure5_d(x, input)
	output = zeros(1,length(input));
	ALPHA_RBSII = x(15);
	ALPHA_RBS_pN25 = x(16);
	K_aTc = x(17);
	K_pLtetO_1 = x(18);
	alpha_pLtetO_1 = x(19);
	alpha_pN25 = x(20);
	beta_pLtetO_1 = x(21);
	n_aTc = x(22);
	n_pLtetO_1 = x(23);
	scale_mCherry_doi_10_1073_pnas_1314114111 = x(185);
	for i = 1:length(input)
		pN25 = 1*alpha_pN25;
		aTc = input(i);
		tetR = ALPHA_RBS_pN25*pN25*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pLtetO_1 = 10*(beta_pLtetO_1 + (alpha_pLtetO_1 - beta_pLtetO_1)*K_pLtetO_1^n_pLtetO_1/(K_pLtetO_1^n_pLtetO_1 + tetR^n_pLtetO_1));
		mCherry = scale_mCherry_doi_10_1073_pnas_1314114111*ALPHA_RBSII*pLtetO_1;
		output(i) = mCherry;
	end
end

function output = simf_doi_10_1038_msb_2013_58_figureS14(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_Unknown = x(59);
	K_AHL = x(83);
	K_pLUX = x(84);
	alpha_pLUX = x(85);
	alpha_pluxR = x(86);
	beta_pLUX = x(87);
	n_AHL = x(88);
	n_pLUX = x(89);
	scale_sfGFP_doi_10_1038_msb_2013_58 = x(186);
	for i = 1:length(input)
		pluxR = 10*alpha_pluxR;
		AHL = input(i);
		luxR = ALPHA_RBS_Unknown*pluxR;
		AHL_luxR = luxR*K_AHL^n_AHL/(K_AHL^n_AHL + AHL^n_AHL);
		pLUX = 10*(beta_pLUX + (alpha_pLUX - beta_pLUX)*AHL_luxR^n_pLUX/(K_pLUX^n_pLUX + AHL_luxR^n_pLUX));
		sfGFP = scale_sfGFP_doi_10_1038_msb_2013_58*ALPHA_RBS_Unknown*pLUX;
		output(i) = sfGFP;
	end
end

function output = simf_doi_10_1021_sb400131a_figure4c_1(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pN25 = x(16);
	ALPHA_UTR1 = x(187);
	K_aTc = x(17);
	K_pLtetO_1 = x(18);
	alpha_pLtetO_1 = x(19);
	alpha_pN25 = x(20);
	beta_pLtetO_1 = x(21);
	n_aTc = x(22);
	n_pLtetO_1 = x(23);
	scale_deGFP_doi_10_1021_sb400131a = x(188);
	for i = 1:length(input)
		pN25 = 1*alpha_pN25;
		aTc = input(i);
		tetR = ALPHA_RBS_pN25*pN25*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pLtetO_1 = 10*(beta_pLtetO_1 + (alpha_pLtetO_1 - beta_pLtetO_1)*K_pLtetO_1^n_pLtetO_1/(K_pLtetO_1^n_pLtetO_1 + tetR^n_pLtetO_1));
		deGFP = scale_deGFP_doi_10_1021_sb400131a*ALPHA_UTR1*pLtetO_1;
		output(i) = deGFP;
	end
end

function output = simf_doi_10_1021_sb400131a_figure4c_2(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_placI = x(25);
	ALPHA_UTR1 = x(187);
	K_IPTG = x(26);
	K_pLlacO_1 = x(71);
	alpha_pLlacO_1 = x(12);
	alpha_placIq = x(29);
	beta_pLlacO_1 = x(72);
	n_IPTG = x(31);
	n_pLlacO_1 = x(73);
	scale_deGFP_doi_10_1021_sb400131a = x(188);
	for i = 1:length(input)
		placIq = 1*alpha_placIq;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placIq*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pLlacO_1 = 10*(beta_pLlacO_1 + (alpha_pLlacO_1 - beta_pLlacO_1)*K_pLlacO_1^n_pLlacO_1/(K_pLlacO_1^n_pLlacO_1 + lacI^n_pLlacO_1));
		deGFP = scale_deGFP_doi_10_1021_sb400131a*ALPHA_UTR1*pLlacO_1;
		output(i) = deGFP;
	end
end

function output = simf_doi_10_1186_1471_2105_13_S4_S11_figure2_1(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0034 = x(90);
	K_AHL = x(83);
	K_pLUX = x(84);
	alpha_pLUX = x(85);
	alpha_pLtetO_1 = x(19);
	beta_pLUX = x(87);
	n_AHL = x(88);
	n_pLUX = x(89);
	scale_mRFP1_doi_10_1186_1471_2105_13_S4_S11 = x(189);
	for i = 1:length(input)
		pLtetO_1 = 10*alpha_pLtetO_1;
		luxR = ALPHA_BBa_B0034*pLtetO_1;
		AHL = input(i);
		AHL_luxR = luxR*K_AHL^n_AHL/(K_AHL^n_AHL + AHL^n_AHL);
		pLUX = 10*(beta_pLUX + (alpha_pLUX - beta_pLUX)*AHL_luxR^n_pLUX/(K_pLUX^n_pLUX + AHL_luxR^n_pLUX));
		mRFP1 = scale_mRFP1_doi_10_1186_1471_2105_13_S4_S11*ALPHA_BBa_B0034*pLUX;
		output(i) = mRFP1;
	end
end

function output = simf_doi_10_1186_1471_2105_13_S4_S11_figure2_2(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0034 = x(90);
	K_AHL = x(83);
	K_pLUX = x(84);
	alpha_pLUX = x(85);
	alpha_pLtetO_1 = x(19);
	beta_pLUX = x(87);
	n_AHL = x(88);
	n_pLUX = x(89);
	scale_mRFP1_doi_10_1186_1471_2105_13_S4_S11 = x(189);
	for i = 1:length(input)
		pLtetO_1 = 5*alpha_pLtetO_1;
		luxR = ALPHA_BBa_B0034*pLtetO_1;
		AHL = input(i);
		AHL_luxR = luxR*K_AHL^n_AHL/(K_AHL^n_AHL + AHL^n_AHL);
		pLUX = 5*(beta_pLUX + (alpha_pLUX - beta_pLUX)*AHL_luxR^n_pLUX/(K_pLUX^n_pLUX + AHL_luxR^n_pLUX));
		mRFP1 = scale_mRFP1_doi_10_1186_1471_2105_13_S4_S11*ALPHA_BBa_B0034*pLUX;
		output(i) = mRFP1;
	end
end

function output = simf_doi_10_1186_1471_2105_13_S4_S11_figure2_3(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0034 = x(90);
	K_AHL = x(83);
	K_pLUX = x(84);
	alpha_pLUX = x(85);
	alpha_pLtetO_1 = x(19);
	beta_pLUX = x(87);
	n_AHL = x(88);
	n_pLUX = x(89);
	scale_mRFP1_doi_10_1186_1471_2105_13_S4_S11 = x(189);
	for i = 1:length(input)
		pLtetO_1 = 1*alpha_pLtetO_1;
		luxR = ALPHA_BBa_B0034*pLtetO_1;
		AHL = input(i);
		AHL_luxR = luxR*K_AHL^n_AHL/(K_AHL^n_AHL + AHL^n_AHL);
		pLUX = 1*(beta_pLUX + (alpha_pLUX - beta_pLUX)*AHL_luxR^n_pLUX/(K_pLUX^n_pLUX + AHL_luxR^n_pLUX));
		mRFP1 = scale_mRFP1_doi_10_1186_1471_2105_13_S4_S11*ALPHA_BBa_B0034*pLUX;
		output(i) = mRFP1;
	end
end

function output = simf_doi_10_1038_nature03461_figure2b_1(x, input)
	output = zeros(1,length(input));
	ALPHA_RBSII = x(15);
	ALPHA_pLUX_RBS = x(190);
	K_AHL = x(83);
	K_pLAC = x(27);
	K_pLUX = x(84);
	alpha_pLAC = x(28);
	alpha_pLUX = x(85);
	alpha_pLUX_L = x(191);
	beta_pLAC = x(30);
	beta_pLUX = x(87);
	n_AHL = x(88);
	n_pLAC = x(32);
	n_pLUX = x(89);
	scale_gfp_LVA_doi_10_1038_nature03461 = x(192);
	for i = 1:length(input)
		pLUX_L = 15*alpha_pLUX_L;
		luxR_G2F = ALPHA_pLUX_RBS*pLUX_L;
		AHL = input(i);
		AHL_luxR = luxR_G2F*K_AHL^n_AHL/(K_AHL^n_AHL + AHL^n_AHL);
		pLUX = 15*(beta_pLUX + (alpha_pLUX - beta_pLUX)*AHL_luxR^n_pLUX/(K_pLUX^n_pLUX + AHL_luxR^n_pLUX));
		lacI_M1 = ALPHA_RBSII*pLUX;
		pLAC = 15*(beta_pLAC + (alpha_pLAC - beta_pLAC)*K_pLAC^n_pLAC/(K_pLAC^n_pLAC + lacI_M1^n_pLAC));
		gfp_LVA = scale_gfp_LVA_doi_10_1038_nature03461*ALPHA_RBSII*pLAC;
		output(i) = gfp_LVA;
	end
end

function output = simf_doi_10_1038_nature03461_figure2b_2(x, input)
	output = zeros(1,length(input));
	ALPHA_RBSII = x(15);
	ALPHA_pLUX_RBS = x(190);
	K_AHL = x(83);
	K_pLAC = x(27);
	K_pLUX = x(84);
	alpha_pLAC = x(28);
	alpha_pLUX = x(85);
	alpha_pLUX_L = x(191);
	beta_pLAC = x(30);
	beta_pLUX = x(87);
	n_AHL = x(88);
	n_pLAC = x(32);
	n_pLUX = x(89);
	scale_gfp_LVA_doi_10_1038_nature03461 = x(192);
	for i = 1:length(input)
		pLUX_L = 15*alpha_pLUX_L;
		luxR = ALPHA_pLUX_RBS*pLUX_L;
		AHL = input(i);
		AHL_luxR = luxR*K_AHL^n_AHL/(K_AHL^n_AHL + AHL^n_AHL);
		pLUX = 15*(beta_pLUX + (alpha_pLUX - beta_pLUX)*AHL_luxR^n_pLUX/(K_pLUX^n_pLUX + AHL_luxR^n_pLUX));
		lacI_M1 = ALPHA_RBSII*pLUX;
		pLAC = 15*(beta_pLAC + (alpha_pLAC - beta_pLAC)*K_pLAC^n_pLAC/(K_pLAC^n_pLAC + lacI_M1^n_pLAC));
		gfp_LVA = scale_gfp_LVA_doi_10_1038_nature03461*ALPHA_RBSII*pLAC;
		output(i) = gfp_LVA;
	end
end

function output = simf_doi_10_1038_nature03461_figure2b_3(x, input)
	output = zeros(1,length(input));
	ALPHA_RBSII = x(15);
	ALPHA_pLUX_RBS = x(190);
	K_AHL = x(83);
	K_pLAC = x(27);
	K_pLUX = x(84);
	alpha_pLAC = x(28);
	alpha_pLUX = x(85);
	alpha_pLUX_L = x(191);
	beta_pLAC = x(30);
	beta_pLUX = x(87);
	n_AHL = x(88);
	n_pLAC = x(32);
	n_pLUX = x(89);
	scale_gfp_LVA_doi_10_1038_nature03461 = x(192);
	for i = 1:length(input)
		pLUX_L = 15*alpha_pLUX_L;
		luxR = ALPHA_pLUX_RBS*pLUX_L;
		AHL = input(i);
		AHL_luxR = luxR*K_AHL^n_AHL/(K_AHL^n_AHL + AHL^n_AHL);
		pLUX = 15*(beta_pLUX + (alpha_pLUX - beta_pLUX)*AHL_luxR^n_pLUX/(K_pLUX^n_pLUX + AHL_luxR^n_pLUX));
		lacI_M1 = ALPHA_RBSII*pLUX;
		pLAC = 15*(beta_pLAC + (alpha_pLAC - beta_pLAC)*K_pLAC^n_pLAC/(K_pLAC^n_pLAC + lacI_M1^n_pLAC));
		gfp_LVA = scale_gfp_LVA_doi_10_1038_nature03461*ALPHA_RBSII*pLAC;
		output(i) = gfp_LVA;
	end
end

function output = simf_doi_10_1073_pnas_0704256104_figureS1b_1(x, input)
	output = zeros(1,length(input));
	ALPHA_RBSII = x(15);
	K_C4HSL = x(172);
	K_pRHL = x(193);
	alpha_pRHL = x(194);
	alpha_placIq = x(29);
	beta_pRHL = x(195);
	n_C4HSL = x(179);
	n_pRHL = x(196);
	scale_gfp_LVA_doi_10_1073_pnas_0704256104 = x(197);
	for i = 1:length(input)
		placIq = 10*alpha_placIq;
		C4HSL = input(i);
		rhlR = ALPHA_RBSII*placIq;
		C4HSL_rhlR = rhlR*K_C4HSL^n_C4HSL/(K_C4HSL^n_C4HSL + C4HSL^n_C4HSL);
		pRHL = 10*(beta_pRHL + (alpha_pRHL - beta_pRHL)*C4HSL_rhlR^n_pRHL/(K_pRHL^n_pRHL + C4HSL_rhlR^n_pRHL));
		gfp_LVA = scale_gfp_LVA_doi_10_1073_pnas_0704256104*ALPHA_RBSII*pRHL;
		output(i) = gfp_LVA;
	end
end

function output = simf_doi_10_1073_pnas_0704256104_figureS1b_2(x, input)
	output = zeros(1,length(input));
	ALPHA_RBSII = x(15);
	K_C_3OC12HSL = x(91);
	K_pLAS = x(198);
	alpha_pLAS = x(199);
	alpha_placIq = x(29);
	beta_pLAS = x(200);
	n_C_3OC12HSL = x(92);
	n_pLAS = x(201);
	scale_gfp_LVA_doi_10_1073_pnas_0704256104 = x(197);
	for i = 1:length(input)
		placIq = 10*alpha_placIq;
		C_3OC12HSL = input(i);
		lasR = ALPHA_RBSII*placIq;
		C_3OC12HSL_lasR = lasR*K_C_3OC12HSL^n_C_3OC12HSL/(K_C_3OC12HSL^n_C_3OC12HSL + C_3OC12HSL^n_C_3OC12HSL);
		pLAS = 10*(beta_pLAS + (alpha_pLAS - beta_pLAS)*C_3OC12HSL_lasR^n_pLAS/(K_pLAS^n_pLAS + C_3OC12HSL_lasR^n_pLAS));
		gfp_LVA = scale_gfp_LVA_doi_10_1073_pnas_0704256104*ALPHA_RBSII*pLAS;
		output(i) = gfp_LVA;
	end
end

function output = simf_doi_10_1038_nature09565_figure1C(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0033 = x(99);
	ALPHA_BBa_B0034 = x(90);
	K_L_arabinose = x(3);
	K_pBAD = x(4);
	alpha_BBa_J23117 = x(54);
	alpha_pBAD = x(5);
	beta_pBAD = x(7);
	n_L_arabinose = x(8);
	n_pBAD = x(9);
	scale_eYFP_doi_10_1038_nature09565 = x(202);
	for i = 1:length(input)
		BBa_J23117 = 15*alpha_BBa_J23117;
		L_arabinose = input(i);
		araC = ALPHA_BBa_B0034*BBa_J23117*K_L_arabinose^n_L_arabinose/(K_L_arabinose^n_L_arabinose + L_arabinose^n_L_arabinose);
		pBAD = 15*(beta_pBAD + (alpha_pBAD - beta_pBAD)*K_pBAD^n_pBAD/(K_pBAD^n_pBAD + araC^n_pBAD));
		eYFP = scale_eYFP_doi_10_1038_nature09565*ALPHA_BBa_B0033*pBAD;
		output(i) = eYFP;
	end
end

function output = simf_doi_10_1038_nature09565_figureS1C(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0032 = x(40);
	ALPHA_BBa_B0033 = x(99);
	K_aTc = x(17);
	K_pLtetO_1 = x(18);
	alpha_BBa_J23117 = x(54);
	alpha_pLtetO_1 = x(19);
	beta_pLtetO_1 = x(21);
	n_aTc = x(22);
	n_pLtetO_1 = x(23);
	scale_eYFP_doi_10_1038_nature09565 = x(202);
	for i = 1:length(input)
		BBa_J23117 = 15*alpha_BBa_J23117;
		aTc = input(i);
		tetR = ALPHA_BBa_B0032*BBa_J23117*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pLtetO_1 = 15*(beta_pLtetO_1 + (alpha_pLtetO_1 - beta_pLtetO_1)*K_pLtetO_1^n_pLtetO_1/(K_pLtetO_1^n_pLtetO_1 + tetR^n_pLtetO_1));
		eYFP = scale_eYFP_doi_10_1038_nature09565*ALPHA_BBa_B0033*pLtetO_1;
		output(i) = eYFP;
	end
end

function output = simf_doi_10_1093_nar_gks597_figureS3_left(x, input)
	output = zeros(1,length(input));
	ALPHA_D103 = x(203);
	ALPHA_RBS_placI = x(25);
	K_IPTG = x(26);
	K_pTAC = x(160);
	alpha_pTAC = x(161);
	alpha_placIq = x(29);
	beta_pTAC = x(162);
	n_IPTG = x(31);
	n_pTAC = x(163);
	scale_mRFP1_doi_10_1093_nar_gks597 = x(204);
	for i = 1:length(input)
		placIq = 5*alpha_placIq;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placIq*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pTAC = 5*(beta_pTAC + (alpha_pTAC - beta_pTAC)*K_pTAC^n_pTAC/(K_pTAC^n_pTAC + lacI^n_pTAC));
		mRFP1 = scale_mRFP1_doi_10_1093_nar_gks597*ALPHA_D103*pTAC;
		output(i) = mRFP1;
	end
end

function output = simf_doi_10_1093_nar_gks597_figureS3_right(x, input)
	output = zeros(1,length(input));
	ALPHA_D103 = x(203);
	ALPHA_RBS_pC_tetR_cassette = x(205);
	K_aTc = x(17);
	K_pTET_tetR_cassette = x(206);
	alpha_pC_tetR_cassette = x(207);
	alpha_pTET_tetR_cassette = x(208);
	beta_pTET_tetR_cassette = x(209);
	n_aTc = x(22);
	n_pTET_tetR_cassette = x(210);
	scale_mRFP1_doi_10_1093_nar_gks597 = x(204);
	for i = 1:length(input)
		pC_tetR_cassette = 5*alpha_pC_tetR_cassette;
		aTc = input(i);
		tetR = ALPHA_RBS_pC_tetR_cassette*pC_tetR_cassette*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pTET_tetR_cassette = 5*(beta_pTET_tetR_cassette + (alpha_pTET_tetR_cassette - beta_pTET_tetR_cassette)*K_pTET_tetR_cassette^n_pTET_tetR_cassette/(K_pTET_tetR_cassette^n_pTET_tetR_cassette + tetR^n_pTET_tetR_cassette));
		mRFP1 = scale_mRFP1_doi_10_1093_nar_gks597*ALPHA_D103*pTET_tetR_cassette;
		output(i) = mRFP1;
	end
end

function output = simf_doi_10_1038_nbt_2401_figureS11(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_A = x(211);
	K_IPTG = x(26);
	K_pTAC = x(160);
	alpha_BBa_J23101 = x(41);
	alpha_pTAC = x(161);
	beta_pTAC = x(162);
	n_IPTG = x(31);
	n_pTAC = x(163);
	scale_sfGFP_doi_10_1038_nbt_2401 = x(212);
	for i = 1:length(input)
		BBa_J23101 = 5*alpha_BBa_J23101;
		IPTG = input(i);
		lacI = ALPHA_RBS_A*BBa_J23101*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pTAC = 5*(beta_pTAC + (alpha_pTAC - beta_pTAC)*K_pTAC^n_pTAC/(K_pTAC^n_pTAC + lacI^n_pTAC));
		sfGFP = scale_sfGFP_doi_10_1038_nbt_2401*ALPHA_RBS_A*pTAC;
		output(i) = sfGFP;
	end
end

function output = simf_doi_10_1038_nbt_2401_figureS13(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_A = x(211);
	K_L_arabinose = x(3);
	K_pBAD = x(4);
	alpha_BBa_J23101 = x(41);
	alpha_pBAD = x(5);
	beta_pBAD = x(7);
	n_L_arabinose = x(8);
	n_pBAD = x(9);
	scale_sfGFP_doi_10_1038_nbt_2401 = x(212);
	for i = 1:length(input)
		BBa_J23101 = 5*alpha_BBa_J23101;
		L_arabinose = input(i);
		araC = ALPHA_RBS_A*BBa_J23101*K_L_arabinose^n_L_arabinose/(K_L_arabinose^n_L_arabinose + L_arabinose^n_L_arabinose);
		pBAD = 5*(beta_pBAD + (alpha_pBAD - beta_pBAD)*K_pBAD^n_pBAD/(K_pBAD^n_pBAD + araC^n_pBAD));
		sfGFP = scale_sfGFP_doi_10_1038_nbt_2401*ALPHA_RBS_A*pBAD;
		output(i) = sfGFP;
	end
end

function output = simf_doi_10_1038_nature11516_figureS8A(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pC = x(2);
	ALPHA_RBS_psicA = x(213);
	K_L_arabinose = x(3);
	K_pBAD = x(4);
	alpha_pBAD = x(5);
	alpha_pC = x(6);
	beta_pBAD = x(7);
	n_L_arabinose = x(8);
	n_pBAD = x(9);
	scale_mRFP1_doi_10_1038_nature11516 = x(214);
	for i = 1:length(input)
		pC = 15*alpha_pC;
		L_arabinose = input(i);
		araC = ALPHA_RBS_pC*pC*K_L_arabinose^n_L_arabinose/(K_L_arabinose^n_L_arabinose + L_arabinose^n_L_arabinose);
		pBAD = 15*(beta_pBAD + (alpha_pBAD - beta_pBAD)*K_pBAD^n_pBAD/(K_pBAD^n_pBAD + araC^n_pBAD));
		mRFP1 = scale_mRFP1_doi_10_1038_nature11516*ALPHA_RBS_psicA*pBAD;
		output(i) = mRFP1;
	end
end

function output = simf_doi_10_1038_nature11516_figureS8B(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_placI = x(25);
	ALPHA_RBS_psicA = x(213);
	K_IPTG = x(26);
	K_pTAC = x(160);
	alpha_pTAC = x(161);
	alpha_placIq = x(29);
	beta_pTAC = x(162);
	n_IPTG = x(31);
	n_pTAC = x(163);
	scale_mRFP1_doi_10_1038_nature11516 = x(214);
	for i = 1:length(input)
		placIq = 15*alpha_placIq;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placIq*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pTAC = 15*(beta_pTAC + (alpha_pTAC - beta_pTAC)*K_pTAC^n_pTAC/(K_pTAC^n_pTAC + lacI^n_pTAC));
		mRFP1 = scale_mRFP1_doi_10_1038_nature11516*ALPHA_RBS_psicA*pTAC;
		output(i) = mRFP1;
	end
end

function output = simf_doi_10_1038_nature11516_figureS8C(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS2 = x(215);
	ALPHA_RBS_psicA = x(213);
	K_AHL = x(83);
	K_pLUX = x(84);
	alpha_BBa_J23100 = x(138);
	alpha_pLUX = x(85);
	beta_pLUX = x(87);
	n_AHL = x(88);
	n_pLUX = x(89);
	scale_mRFP1_doi_10_1038_nature11516 = x(214);
	for i = 1:length(input)
		BBa_J23100 = 15*alpha_BBa_J23100;
		AHL = input(i);
		luxR = ALPHA_RBS2*BBa_J23100;
		AHL_luxR = luxR*K_AHL^n_AHL/(K_AHL^n_AHL + AHL^n_AHL);
		pLUX = 15*(beta_pLUX + (alpha_pLUX - beta_pLUX)*AHL_luxR^n_pLUX/(K_pLUX^n_pLUX + AHL_luxR^n_pLUX));
		mRFP1 = scale_mRFP1_doi_10_1038_nature11516*ALPHA_RBS_psicA*pLUX;
		output(i) = mRFP1;
	end
end

function output = simf_doi_10_1038_nature11516_figureS8D(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS2 = x(215);
	ALPHA_RBS_psicA = x(213);
	K_AHL = x(83);
	K_pLUX_star_ = x(216);
	alpha_BBa_J23100 = x(138);
	alpha_pLUX_star_ = x(217);
	beta_pLUX_star_ = x(218);
	n_AHL = x(88);
	n_pLUX_star_ = x(219);
	scale_mRFP1_doi_10_1038_nature11516 = x(214);
	for i = 1:length(input)
		BBa_J23100 = 15*alpha_BBa_J23100;
		AHL = input(i);
		luxR = ALPHA_RBS2*BBa_J23100;
		AHL_luxR = luxR*K_AHL^n_AHL/(K_AHL^n_AHL + AHL^n_AHL);
		pLUX_star_ = 15*(beta_pLUX_star_ + (alpha_pLUX_star_ - beta_pLUX_star_)*AHL_luxR^n_pLUX_star_/(K_pLUX_star_^n_pLUX_star_ + AHL_luxR^n_pLUX_star_));
		mRFP1 = scale_mRFP1_doi_10_1038_nature11516*ALPHA_RBS_psicA*pLUX_star_;
		output(i) = mRFP1;
	end
end

function output = simf_doi_10_1038_nbt_2149_Supplement_page17(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pET_29b = x(220);
	ALPHA_RBS_placI = x(25);
	K_IPTG = x(26);
	K_placUV5 = x(221);
	alpha_placI = x(37);
	alpha_placUV5 = x(222);
	beta_placUV5 = x(223);
	n_IPTG = x(31);
	n_placUV5 = x(224);
	scale_mRFP1_doi_10_1038_nbt_2149 = x(225);
	for i = 1:length(input)
		placI = 10*alpha_placI;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placI*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		placUV5 = 10*(beta_placUV5 + (alpha_placUV5 - beta_placUV5)*K_placUV5^n_placUV5/(K_placUV5^n_placUV5 + lacI^n_placUV5));
		mRFP1 = scale_mRFP1_doi_10_1038_nbt_2149*ALPHA_RBS_pET_29b*placUV5;
		output(i) = mRFP1;
	end
end

function output = simf_doi_10_1186_1754_1611_5_12_Supplement_page1(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pET_29b = x(220);
	ALPHA_RBS_placI = x(25);
	K_IPTG = x(26);
	K_pTrc = x(76);
	alpha_pTrc = x(77);
	alpha_placIq = x(29);
	beta_pTrc = x(78);
	n_IPTG = x(31);
	n_pTrc = x(79);
	scale_mRFP1_doi_10_1186_1754_1611_5_12 = x(226);
	for i = 1:length(input)
		placIq = 10*alpha_placIq;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placIq*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pTrc = 10*(beta_pTrc + (alpha_pTrc - beta_pTrc)*K_pTrc^n_pTrc/(K_pTrc^n_pTrc + lacI^n_pTrc));
		mRFP1 = scale_mRFP1_doi_10_1186_1754_1611_5_12*ALPHA_RBS_pET_29b*pTrc;
		output(i) = mRFP1;
	end
end

function output = simf_doi_10_1186_1754_1611_5_12_Supplement_page2(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pET_29b = x(220);
	ALPHA_RBS_placI = x(25);
	K_IPTG = x(26);
	K_pTrc = x(76);
	alpha_pTrc = x(77);
	alpha_placIq = x(29);
	beta_pTrc = x(78);
	n_IPTG = x(31);
	n_pTrc = x(79);
	scale_GFP_uv_doi_10_1186_1754_1611_5_12 = x(227);
	for i = 1:length(input)
		placIq = 20*alpha_placIq;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placIq*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pTrc = 20*(beta_pTrc + (alpha_pTrc - beta_pTrc)*K_pTrc^n_pTrc/(K_pTrc^n_pTrc + lacI^n_pTrc));
		GFP_uv = scale_GFP_uv_doi_10_1186_1754_1611_5_12*ALPHA_RBS_pET_29b*pTrc;
		output(i) = GFP_uv;
	end
end

function output = simf_doi_10_1186_1754_1611_5_12_Supplement_page3(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pET_29b = x(220);
	ALPHA_RBS_placI = x(25);
	K_IPTG = x(26);
	K_pTrc = x(76);
	alpha_pTrc = x(77);
	alpha_placIq = x(29);
	beta_pTrc = x(78);
	n_IPTG = x(31);
	n_pTrc = x(79);
	scale_mRFP1_doi_10_1186_1754_1611_5_12 = x(226);
	for i = 1:length(input)
		placIq = 15*alpha_placIq;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placIq*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pTrc = 15*(beta_pTrc + (alpha_pTrc - beta_pTrc)*K_pTrc^n_pTrc/(K_pTrc^n_pTrc + lacI^n_pTrc));
		mRFP1 = scale_mRFP1_doi_10_1186_1754_1611_5_12*ALPHA_RBS_pET_29b*pTrc;
		output(i) = mRFP1;
	end
end

function output = simf_doi_10_1186_1754_1611_5_12_Supplement_page4(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pET_29b = x(220);
	ALPHA_RBS_placI = x(25);
	K_IPTG = x(26);
	K_pTrc = x(76);
	alpha_pTrc = x(77);
	alpha_placIq = x(29);
	beta_pTrc = x(78);
	n_IPTG = x(31);
	n_pTrc = x(79);
	scale_mRFP1_doi_10_1186_1754_1611_5_12 = x(226);
	for i = 1:length(input)
		placIq = 5*alpha_placIq;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placIq*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pTrc = 5*(beta_pTrc + (alpha_pTrc - beta_pTrc)*K_pTrc^n_pTrc/(K_pTrc^n_pTrc + lacI^n_pTrc));
		mRFP1 = scale_mRFP1_doi_10_1186_1754_1611_5_12*ALPHA_RBS_pET_29b*pTrc;
		output(i) = mRFP1;
	end
end

function output = simf_doi_10_1186_1754_1611_5_12_Supplement_page5(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pC_tetR = x(228);
	ALPHA_RBS_pET_29b = x(220);
	K_aTc = x(17);
	K_pTET_operon = x(229);
	alpha_pC_tetR = x(230);
	alpha_pTET_operon = x(231);
	beta_pTET_operon = x(232);
	n_aTc = x(22);
	n_pTET_operon = x(233);
	scale_mRFP1_doi_10_1186_1754_1611_5_12 = x(226);
	for i = 1:length(input)
		pC_tetR = 10*alpha_pC_tetR;
		aTc = input(i);
		tetR = ALPHA_RBS_pC_tetR*pC_tetR*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pTET_operon = 10*(beta_pTET_operon + (alpha_pTET_operon - beta_pTET_operon)*K_pTET_operon^n_pTET_operon/(K_pTET_operon^n_pTET_operon + tetR^n_pTET_operon));
		mRFP1 = scale_mRFP1_doi_10_1186_1754_1611_5_12*ALPHA_RBS_pET_29b*pTET_operon;
		output(i) = mRFP1;
	end
end

function output = simf_doi_10_1186_1754_1611_5_12_Supplement_page6(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pC_tetR = x(228);
	ALPHA_RBS_pET_29b = x(220);
	K_aTc = x(17);
	K_pTET_operon = x(229);
	alpha_pC_tetR = x(230);
	alpha_pTET_operon = x(231);
	beta_pTET_operon = x(232);
	n_aTc = x(22);
	n_pTET_operon = x(233);
	scale_GFP_uv_doi_10_1186_1754_1611_5_12 = x(227);
	for i = 1:length(input)
		pC_tetR = 20*alpha_pC_tetR;
		aTc = input(i);
		tetR = ALPHA_RBS_pC_tetR*pC_tetR*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pTET_operon = 20*(beta_pTET_operon + (alpha_pTET_operon - beta_pTET_operon)*K_pTET_operon^n_pTET_operon/(K_pTET_operon^n_pTET_operon + tetR^n_pTET_operon));
		GFP_uv = scale_GFP_uv_doi_10_1186_1754_1611_5_12*ALPHA_RBS_pET_29b*pTET_operon;
		output(i) = GFP_uv;
	end
end

function output = simf_doi_10_1186_1754_1611_5_12_Supplement_page7(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pC_tetR = x(228);
	ALPHA_RBS_pET_29b = x(220);
	K_aTc = x(17);
	K_pTET_operon = x(229);
	alpha_pC_tetR = x(230);
	alpha_pTET_operon = x(231);
	beta_pTET_operon = x(232);
	n_aTc = x(22);
	n_pTET_operon = x(233);
	scale_mRFP1_doi_10_1186_1754_1611_5_12 = x(226);
	for i = 1:length(input)
		pC_tetR = 15*alpha_pC_tetR;
		aTc = input(i);
		tetR = ALPHA_RBS_pC_tetR*pC_tetR*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pTET_operon = 15*(beta_pTET_operon + (alpha_pTET_operon - beta_pTET_operon)*K_pTET_operon^n_pTET_operon/(K_pTET_operon^n_pTET_operon + tetR^n_pTET_operon));
		mRFP1 = scale_mRFP1_doi_10_1186_1754_1611_5_12*ALPHA_RBS_pET_29b*pTET_operon;
		output(i) = mRFP1;
	end
end

function output = simf_doi_10_1186_1754_1611_5_12_Supplement_page8(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pC_tetR = x(228);
	ALPHA_RBS_pET_29b = x(220);
	K_aTc = x(17);
	K_pTET_operon = x(229);
	alpha_pC_tetR = x(230);
	alpha_pTET_operon = x(231);
	beta_pTET_operon = x(232);
	n_aTc = x(22);
	n_pTET_operon = x(233);
	scale_mRFP1_doi_10_1186_1754_1611_5_12 = x(226);
	for i = 1:length(input)
		pC_tetR = 5*alpha_pC_tetR;
		aTc = input(i);
		tetR = ALPHA_RBS_pC_tetR*pC_tetR*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pTET_operon = 5*(beta_pTET_operon + (alpha_pTET_operon - beta_pTET_operon)*K_pTET_operon^n_pTET_operon/(K_pTET_operon^n_pTET_operon + tetR^n_pTET_operon));
		mRFP1 = scale_mRFP1_doi_10_1186_1754_1611_5_12*ALPHA_RBS_pET_29b*pTET_operon;
		output(i) = mRFP1;
	end
end

function output = simf_doi_10_1186_1754_1611_5_12_Supplement_page17(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pET_29b = x(220);
	ALPHA_RBS_placI = x(25);
	K_IPTG = x(26);
	K_placUV5 = x(221);
	alpha_placIq = x(29);
	alpha_placUV5 = x(222);
	beta_placUV5 = x(223);
	n_IPTG = x(31);
	n_placUV5 = x(224);
	scale_mRFP1_doi_10_1186_1754_1611_5_12 = x(226);
	for i = 1:length(input)
		placIq = 10*alpha_placIq;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placIq*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		placUV5 = 10*(beta_placUV5 + (alpha_placUV5 - beta_placUV5)*K_placUV5^n_placUV5/(K_placUV5^n_placUV5 + lacI^n_placUV5));
		mRFP1 = scale_mRFP1_doi_10_1186_1754_1611_5_12*ALPHA_RBS_pET_29b*placUV5;
		output(i) = mRFP1;
	end
end

function output = simf_doi_10_1186_1754_1611_5_12_Supplement_page18(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pET_29b = x(220);
	ALPHA_RBS_placI = x(25);
	K_IPTG = x(26);
	K_placUV5 = x(221);
	alpha_placIq = x(29);
	alpha_placUV5 = x(222);
	beta_placUV5 = x(223);
	n_IPTG = x(31);
	n_placUV5 = x(224);
	scale_GFP_uv_doi_10_1186_1754_1611_5_12 = x(227);
	for i = 1:length(input)
		placIq = 20*alpha_placIq;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placIq*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		placUV5 = 20*(beta_placUV5 + (alpha_placUV5 - beta_placUV5)*K_placUV5^n_placUV5/(K_placUV5^n_placUV5 + lacI^n_placUV5));
		GFP_uv = scale_GFP_uv_doi_10_1186_1754_1611_5_12*ALPHA_RBS_pET_29b*placUV5;
		output(i) = GFP_uv;
	end
end

function output = simf_doi_10_1186_1754_1611_5_12_Supplement_page19(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pET_29b = x(220);
	ALPHA_RBS_placI = x(25);
	K_IPTG = x(26);
	K_placUV5 = x(221);
	alpha_placIq = x(29);
	alpha_placUV5 = x(222);
	beta_placUV5 = x(223);
	n_IPTG = x(31);
	n_placUV5 = x(224);
	scale_mRFP1_doi_10_1186_1754_1611_5_12 = x(226);
	for i = 1:length(input)
		placIq = 15*alpha_placIq;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placIq*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		placUV5 = 15*(beta_placUV5 + (alpha_placUV5 - beta_placUV5)*K_placUV5^n_placUV5/(K_placUV5^n_placUV5 + lacI^n_placUV5));
		mRFP1 = scale_mRFP1_doi_10_1186_1754_1611_5_12*ALPHA_RBS_pET_29b*placUV5;
		output(i) = mRFP1;
	end
end

function output = simf_doi_10_1186_1754_1611_5_12_Supplement_page20(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pET_29b = x(220);
	ALPHA_RBS_placI = x(25);
	K_IPTG = x(26);
	K_placUV5 = x(221);
	alpha_placIq = x(29);
	alpha_placUV5 = x(222);
	beta_placUV5 = x(223);
	n_IPTG = x(31);
	n_placUV5 = x(224);
	scale_mRFP1_doi_10_1186_1754_1611_5_12 = x(226);
	for i = 1:length(input)
		placIq = 5*alpha_placIq;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placIq*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		placUV5 = 5*(beta_placUV5 + (alpha_placUV5 - beta_placUV5)*K_placUV5^n_placUV5/(K_placUV5^n_placUV5 + lacI^n_placUV5));
		mRFP1 = scale_mRFP1_doi_10_1186_1754_1611_5_12*ALPHA_RBS_pET_29b*placUV5;
		output(i) = mRFP1;
	end
end

function output = simf_doi_10_1186_1754_1611_5_12_Supplement_page21(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pET_29b = x(220);
	ALPHA_RBS_placI = x(25);
	K_IPTG = x(26);
	K_pLlacO_1 = x(71);
	alpha_pLlacO_1 = x(12);
	alpha_placIq = x(29);
	beta_pLlacO_1 = x(72);
	n_IPTG = x(31);
	n_pLlacO_1 = x(73);
	scale_mRFP1_doi_10_1186_1754_1611_5_12 = x(226);
	for i = 1:length(input)
		placIq = 10*alpha_placIq;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placIq*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pLlacO_1 = 10*(beta_pLlacO_1 + (alpha_pLlacO_1 - beta_pLlacO_1)*K_pLlacO_1^n_pLlacO_1/(K_pLlacO_1^n_pLlacO_1 + lacI^n_pLlacO_1));
		mRFP1 = scale_mRFP1_doi_10_1186_1754_1611_5_12*ALPHA_RBS_pET_29b*pLlacO_1;
		output(i) = mRFP1;
	end
end

function output = simf_doi_10_1186_1754_1611_5_12_Supplement_page22(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pET_29b = x(220);
	ALPHA_RBS_placI = x(25);
	K_IPTG = x(26);
	K_pLlacO_1 = x(71);
	alpha_pLlacO_1 = x(12);
	alpha_placIq = x(29);
	beta_pLlacO_1 = x(72);
	n_IPTG = x(31);
	n_pLlacO_1 = x(73);
	scale_GFP_uv_doi_10_1186_1754_1611_5_12 = x(227);
	for i = 1:length(input)
		placIq = 20*alpha_placIq;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placIq*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pLlacO_1 = 20*(beta_pLlacO_1 + (alpha_pLlacO_1 - beta_pLlacO_1)*K_pLlacO_1^n_pLlacO_1/(K_pLlacO_1^n_pLlacO_1 + lacI^n_pLlacO_1));
		GFP_uv = scale_GFP_uv_doi_10_1186_1754_1611_5_12*ALPHA_RBS_pET_29b*pLlacO_1;
		output(i) = GFP_uv;
	end
end

function output = simf_doi_10_1186_1754_1611_5_12_Supplement_page23(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pET_29b = x(220);
	ALPHA_RBS_placI = x(25);
	K_IPTG = x(26);
	K_pLlacO_1 = x(71);
	alpha_pLlacO_1 = x(12);
	alpha_placIq = x(29);
	beta_pLlacO_1 = x(72);
	n_IPTG = x(31);
	n_pLlacO_1 = x(73);
	scale_mRFP1_doi_10_1186_1754_1611_5_12 = x(226);
	for i = 1:length(input)
		placIq = 15*alpha_placIq;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placIq*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pLlacO_1 = 15*(beta_pLlacO_1 + (alpha_pLlacO_1 - beta_pLlacO_1)*K_pLlacO_1^n_pLlacO_1/(K_pLlacO_1^n_pLlacO_1 + lacI^n_pLlacO_1));
		mRFP1 = scale_mRFP1_doi_10_1186_1754_1611_5_12*ALPHA_RBS_pET_29b*pLlacO_1;
		output(i) = mRFP1;
	end
end

function output = simf_doi_10_1186_1754_1611_5_12_Supplement_page24(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pET_29b = x(220);
	ALPHA_RBS_placI = x(25);
	K_IPTG = x(26);
	K_pLlacO_1 = x(71);
	alpha_pLlacO_1 = x(12);
	alpha_placIq = x(29);
	beta_pLlacO_1 = x(72);
	n_IPTG = x(31);
	n_pLlacO_1 = x(73);
	scale_mRFP1_doi_10_1186_1754_1611_5_12 = x(226);
	for i = 1:length(input)
		placIq = 5*alpha_placIq;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placIq*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pLlacO_1 = 5*(beta_pLlacO_1 + (alpha_pLlacO_1 - beta_pLlacO_1)*K_pLlacO_1^n_pLlacO_1/(K_pLlacO_1^n_pLlacO_1 + lacI^n_pLlacO_1));
		mRFP1 = scale_mRFP1_doi_10_1186_1754_1611_5_12*ALPHA_RBS_pET_29b*pLlacO_1;
		output(i) = mRFP1;
	end
end

function output = simf_doi_10_1186_1754_1611_5_12_Supplement_page25(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pET_29b = x(220);
	ALPHA_RBS_placI = x(25);
	K_IPTG = x(26);
	K_pT7lacO = x(234);
	alpha_pT7lacO = x(235);
	alpha_placIq = x(29);
	beta_pT7lacO = x(236);
	n_IPTG = x(31);
	n_pT7lacO = x(237);
	scale_mRFP1_doi_10_1186_1754_1611_5_12 = x(226);
	for i = 1:length(input)
		placIq = 10*alpha_placIq;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placIq*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pT7lacO = 10*(beta_pT7lacO + (alpha_pT7lacO - beta_pT7lacO)*K_pT7lacO^n_pT7lacO/(K_pT7lacO^n_pT7lacO + lacI^n_pT7lacO));
		mRFP1 = scale_mRFP1_doi_10_1186_1754_1611_5_12*ALPHA_RBS_pET_29b*pT7lacO;
		output(i) = mRFP1;
	end
end

function output = simf_doi_10_1186_1754_1611_5_12_Supplement_page26(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pET_29b = x(220);
	ALPHA_RBS_placI = x(25);
	K_IPTG = x(26);
	K_pT7lacO = x(234);
	alpha_pT7lacO = x(235);
	alpha_placIq = x(29);
	beta_pT7lacO = x(236);
	n_IPTG = x(31);
	n_pT7lacO = x(237);
	scale_GFP_uv_doi_10_1186_1754_1611_5_12 = x(227);
	for i = 1:length(input)
		placIq = 20*alpha_placIq;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placIq*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pT7lacO = 20*(beta_pT7lacO + (alpha_pT7lacO - beta_pT7lacO)*K_pT7lacO^n_pT7lacO/(K_pT7lacO^n_pT7lacO + lacI^n_pT7lacO));
		GFP_uv = scale_GFP_uv_doi_10_1186_1754_1611_5_12*ALPHA_RBS_pET_29b*pT7lacO;
		output(i) = GFP_uv;
	end
end

function output = simf_doi_10_1186_1754_1611_5_12_Supplement_page27(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pET_29b = x(220);
	ALPHA_RBS_placI = x(25);
	K_IPTG = x(26);
	K_pT7lacO = x(234);
	alpha_pT7lacO = x(235);
	alpha_placIq = x(29);
	beta_pT7lacO = x(236);
	n_IPTG = x(31);
	n_pT7lacO = x(237);
	scale_mRFP1_doi_10_1186_1754_1611_5_12 = x(226);
	for i = 1:length(input)
		placIq = 15*alpha_placIq;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placIq*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pT7lacO = 15*(beta_pT7lacO + (alpha_pT7lacO - beta_pT7lacO)*K_pT7lacO^n_pT7lacO/(K_pT7lacO^n_pT7lacO + lacI^n_pT7lacO));
		mRFP1 = scale_mRFP1_doi_10_1186_1754_1611_5_12*ALPHA_RBS_pET_29b*pT7lacO;
		output(i) = mRFP1;
	end
end

function output = simf_doi_10_1186_1754_1611_5_12_Supplement_page28(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pET_29b = x(220);
	ALPHA_RBS_placI = x(25);
	K_IPTG = x(26);
	K_pT7lacO = x(234);
	alpha_pT7lacO = x(235);
	alpha_placIq = x(29);
	beta_pT7lacO = x(236);
	n_IPTG = x(31);
	n_pT7lacO = x(237);
	scale_mRFP1_doi_10_1186_1754_1611_5_12 = x(226);
	for i = 1:length(input)
		placIq = 5*alpha_placIq;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placIq*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pT7lacO = 5*(beta_pT7lacO + (alpha_pT7lacO - beta_pT7lacO)*K_pT7lacO^n_pT7lacO/(K_pT7lacO^n_pT7lacO + lacI^n_pT7lacO));
		mRFP1 = scale_mRFP1_doi_10_1186_1754_1611_5_12*ALPHA_RBS_pET_29b*pT7lacO;
		output(i) = mRFP1;
	end
end

function output = simf_doi_10_1186_1754_1611_5_12_Supplement_page29(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pC = x(2);
	ALPHA_RBS_pET_29b = x(220);
	K_L_arabinose = x(3);
	K_pBAD = x(4);
	alpha_pBAD = x(5);
	alpha_pC = x(6);
	beta_pBAD = x(7);
	n_L_arabinose = x(8);
	n_pBAD = x(9);
	scale_mRFP1_doi_10_1186_1754_1611_5_12 = x(226);
	for i = 1:length(input)
		pC = 1*alpha_pC;
		L_arabinose = input(i);
		araC = ALPHA_RBS_pC*pC*K_L_arabinose^n_L_arabinose/(K_L_arabinose^n_L_arabinose + L_arabinose^n_L_arabinose);
		pBAD = 10*(beta_pBAD + (alpha_pBAD - beta_pBAD)*K_pBAD^n_pBAD/(K_pBAD^n_pBAD + araC^n_pBAD));
		mRFP1 = scale_mRFP1_doi_10_1186_1754_1611_5_12*ALPHA_RBS_pET_29b*pBAD;
		output(i) = mRFP1;
	end
end

function output = simf_doi_10_1186_1754_1611_5_12_Supplement_page30(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pC = x(2);
	ALPHA_RBS_pET_29b = x(220);
	K_L_arabinose = x(3);
	K_pBAD = x(4);
	alpha_pBAD = x(5);
	alpha_pC = x(6);
	beta_pBAD = x(7);
	n_L_arabinose = x(8);
	n_pBAD = x(9);
	scale_GFP_uv_doi_10_1186_1754_1611_5_12 = x(227);
	for i = 1:length(input)
		pC = 1*alpha_pC;
		L_arabinose = input(i);
		araC = ALPHA_RBS_pC*pC*K_L_arabinose^n_L_arabinose/(K_L_arabinose^n_L_arabinose + L_arabinose^n_L_arabinose);
		pBAD = 20*(beta_pBAD + (alpha_pBAD - beta_pBAD)*K_pBAD^n_pBAD/(K_pBAD^n_pBAD + araC^n_pBAD));
		GFP_uv = scale_GFP_uv_doi_10_1186_1754_1611_5_12*ALPHA_RBS_pET_29b*pBAD;
		output(i) = GFP_uv;
	end
end

function output = simf_doi_10_1186_1754_1611_5_12_Supplement_page31(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pC = x(2);
	ALPHA_RBS_pET_29b = x(220);
	K_L_arabinose = x(3);
	K_pBAD = x(4);
	alpha_pBAD = x(5);
	alpha_pC = x(6);
	beta_pBAD = x(7);
	n_L_arabinose = x(8);
	n_pBAD = x(9);
	scale_mRFP1_doi_10_1186_1754_1611_5_12 = x(226);
	for i = 1:length(input)
		pC = 1*alpha_pC;
		L_arabinose = input(i);
		araC = ALPHA_RBS_pC*pC*K_L_arabinose^n_L_arabinose/(K_L_arabinose^n_L_arabinose + L_arabinose^n_L_arabinose);
		pBAD = 15*(beta_pBAD + (alpha_pBAD - beta_pBAD)*K_pBAD^n_pBAD/(K_pBAD^n_pBAD + araC^n_pBAD));
		mRFP1 = scale_mRFP1_doi_10_1186_1754_1611_5_12*ALPHA_RBS_pET_29b*pBAD;
		output(i) = mRFP1;
	end
end

function output = simf_doi_10_1186_1754_1611_5_12_Supplement_page32(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pC = x(2);
	ALPHA_RBS_pET_29b = x(220);
	K_L_arabinose = x(3);
	K_pBAD = x(4);
	alpha_pBAD = x(5);
	alpha_pC = x(6);
	beta_pBAD = x(7);
	n_L_arabinose = x(8);
	n_pBAD = x(9);
	scale_mRFP1_doi_10_1186_1754_1611_5_12 = x(226);
	for i = 1:length(input)
		pC = 1*alpha_pC;
		L_arabinose = input(i);
		araC = ALPHA_RBS_pC*pC*K_L_arabinose^n_L_arabinose/(K_L_arabinose^n_L_arabinose + L_arabinose^n_L_arabinose);
		pBAD = 5*(beta_pBAD + (alpha_pBAD - beta_pBAD)*K_pBAD^n_pBAD/(K_pBAD^n_pBAD + araC^n_pBAD));
		mRFP1 = scale_mRFP1_doi_10_1186_1754_1611_5_12*ALPHA_RBS_pET_29b*pBAD;
		output(i) = mRFP1;
	end
end

function output = simf_doi_10_1371_journal_pone_0039407_figure6B(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0034 = x(90);
	ALPHA_RBS_pLAC = x(36);
	ALPHA_RBS_placI = x(25);
	K_IPTG = x(26);
	K_pLAC = x(27);
	K_pLlacO_1 = x(71);
	alpha_pLAC = x(28);
	alpha_pLlacO_1 = x(12);
	alpha_placI = x(37);
	beta_pLAC = x(30);
	beta_pLlacO_1 = x(72);
	n_IPTG = x(31);
	n_pLAC = x(32);
	n_pLlacO_1 = x(73);
	scale_mRFP1_doi_10_1371_journal_pone_0039407 = x(238);
	for i = 1:length(input)
		placI = 1*alpha_placI;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placI*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pLAC = 1*(beta_pLAC + (alpha_pLAC - beta_pLAC)*K_pLAC^n_pLAC/(K_pLAC^n_pLAC + lacI^n_pLAC));
		lacZ = ALPHA_RBS_pLAC*pLAC;
		pLlacO_1 = 5*(beta_pLlacO_1 + (alpha_pLlacO_1 - beta_pLlacO_1)*K_pLlacO_1^n_pLlacO_1/(K_pLlacO_1^n_pLlacO_1 + lacI^n_pLlacO_1));
		mRFP1 = scale_mRFP1_doi_10_1371_journal_pone_0039407*ALPHA_BBa_B0034*pLlacO_1;
		output(i) = mRFP1;
	end
end

function output = simf_doi_10_1073_pnas_1301301110_sd03_xls_1(x)
	ALPHA_BBa_B0032 = x(40);
	alpha_BBa_J23117 = x(54);
	scale_sfGFP_doi_10_1073_pnas_1301301110 = x(239);
	BBa_J23117 = 10*alpha_BBa_J23117;
	sfGFP = scale_sfGFP_doi_10_1073_pnas_1301301110*ALPHA_BBa_B0032*BBa_J23117;
	output = sfGFP;
end

function output = simf_doi_10_1073_pnas_1301301110_sd03_xls_2(x)
	ALPHA_BBa_B0032 = x(40);
	alpha_BBa_J23101 = x(41);
	scale_sfGFP_doi_10_1073_pnas_1301301110 = x(239);
	BBa_J23101 = 10*alpha_BBa_J23101;
	sfGFP = scale_sfGFP_doi_10_1073_pnas_1301301110*ALPHA_BBa_B0032*BBa_J23101;
	output = sfGFP;
end

function x_opt = simultaneous_fitting(ensemble_size, ite_num, func_ite_num)
	parameter_name = {'ALPHA_RBS_pBAD'; 'ALPHA_RBS_pC'; 'K_L_arabinose'; 'K_pBAD'; 'alpha_pBAD'; 'alpha_pC'; 'beta_pBAD'; 'n_L_arabinose'; 'n_pBAD'; 'scale_GFPuv_doi_10_1128_'; 'ALPHA_BBa_B0030'; 'alpha_pLlacO_1'; 'scale_mCherry_doi_10_103'; 'scale_AP_http___jb_asm_o'; 'ALPHA_RBSII'; 'ALPHA_RBS_pN25'; 'K_aTc'; 'K_pLtetO_1'; 'alpha_pLtetO_1'; 'alpha_pN25'; 'beta_pLtetO_1'; 'n_aTc'; 'n_pLtetO_1'; 'scale_luc_doi_10_1093_na'; 'ALPHA_RBS_placI'; 'K_IPTG'; 'K_pLAC'; 'alpha_pLAC'; 'alpha_placIq'; 'beta_pLAC'; 'n_IPTG'; 'n_pLAC'; 'scale_eYFP_http___dspace'; 'ALPHA_RBS_Unknown_Hoosha'; 'scale_eYFP_doi_10_1073pn'; 'ALPHA_RBS_pLAC'; 'alpha_placI'; 'scale_lacZ_doi_10_1073_p'; 'scale_lacZ_doi_10_1073_p'; 'ALPHA_BBa_B0032'; 'alpha_BBa_J23101'; 'scale_GFPmut3b_doi_10_11'; 'alpha_BBa_J23116'; 'alpha_BBa_J23150'; 'alpha_BBa_J23151'; 'alpha_BBa_J23102'; 'alpha_BBa_J23109'; 'scale_GFPmut3b_doi_10_10'; 'alpha_BBa_J23114'; 'alpha_BBa_J23115'; 'alpha_BBa_J23105'; 'alpha_BBa_J23106'; 'scale_GFPmut3b_doi_10_10'; 'alpha_BBa_J23117'; 'K_pTET_star_'; 'alpha_pTET_star_'; 'beta_pTET_star_'; 'n_pTET_star_'; 'ALPHA_RBS_Unknown'; 'ALPHA_SAL_RBS'; 'K_Sal'; 'K_pSAL'; 'alpha_pCONST_UNKNOWN'; 'alpha_pSAL'; 'beta_pSAL'; 'n_Sal'; 'n_pSAL'; 'scale_gfp_doi_10_1038_ms'; 'ALPHA_parent_RBS'; 'ALPHA_T710'; 'K_pLlacO_1'; 'beta_pLlacO_1'; 'n_pLlacO_1'; 'scale_cfp_doi_10_1534_ge'; 'ALPHA_g10'; 'K_pTrc'; 'alpha_pTrc'; 'beta_pTrc'; 'n_pTrc'; 'scale_gfp_LAA_doi_10_101'; 'ALPHA_RBS_pLUX'; 'ALPHA_RBS_pluxR'; 'K_AHL'; 'K_pLUX'; 'alpha_pLUX'; 'alpha_pluxR'; 'beta_pLUX'; 'n_AHL'; 'n_pLUX'; 'ALPHA_BBa_B0034'; 'K_C_3OC12HSL'; 'n_C_3OC12HSL'; 'scale_GFPmut3b_doi_10_10'; 'ALPHA_RBS_UNKNOWN'; 'scale_gfp_LVA_doi_10_103'; 'alpha_pTET'; 'scale_GFPmut3b_doi_10_10'; 'scale_mRFP1_doi_10_1093_'; 'ALPHA_BBa_B0033'; 'scale_gfp_uv_doi_10_1021'; 'ALPHA_RBS_Glk'; 'scale_Glk_doi_10_1021_sb'; 'alpha_plac_ara_1'; 'scale_gfpmut3_1_doi_10_1'; 'ALPHA_RBS_UTR'; 'ALPHA_RBS_pLAMBDA'; 'K_pLAMBDA'; 'alpha_pLAMBDA'; 'beta_pLAMBDA'; 'n_pLAMBDA'; 'scale_gfp_doi_10_1186_17'; 'scale_gfpmut2_doi_10_118'; 'alpha_pCAT'; 'scale_eYFP_dx_doi_org_10'; 'scale_eYFP_doi_10_1073_p'; 'ALPHA_RBS_D'; 'ALPHA_RBS_E'; 'K_pTrc_2'; 'alpha_pLs1con'; 'alpha_pTrc_2'; 'beta_pTrc_2'; 'n_pTrc_2'; 'scale_GFPmut3b_doi_10_10'; 'ALPHA_RBS_native'; 'ALPHA_RBS_tetR'; 'K_pJ23114_lacO'; 'K_pTET'; 'alpha_pC_tetR_UNKNOWN'; 'alpha_pJ23114_lacO'; 'beta_pJ23114_lacO'; 'beta_pTET'; 'n_pJ23114_lacO'; 'n_pTET'; 'scale_glk_doi_10_1016_j_'; 'scale_lacZ_doi_10_1006_p'; 'scale_GFPmut3b_doi_10_10'; 'K_pLasRV'; 'alpha_BBa_J23100'; 'alpha_pLasRV'; 'beta_pLasRV'; 'n_pLasRV'; 'scale_GFPmut3b_doi_10_10'; 'K_pLasR3'; 'alpha_pLasR3'; 'beta_pLasR3'; 'n_pLasR3'; 'K_pLasRL'; 'alpha_pLasRL'; 'beta_pLasRL'; 'n_pLasRL'; 'ALPHA_RBS_KDL'; 'K_Mg2_plus__plus_'; 'K_pMgrB'; 'alpha_pMgrB'; 'beta_pMgrB'; 'n_Mg2_plus__plus_'; 'n_pMgrB'; 'scale_gfp_doi_10_1038_nm'; 'ALPHA_rbs1'; 'K_pTAC'; 'alpha_pTAC'; 'beta_pTAC'; 'n_pTAC'; 'scale_eYFP_doi_10_1038_n'; 'alpha_J23119_tetO'; 'ALPHA_RBS_PtetA'; 'K_PtetA'; 'alpha_PtetA'; 'beta_PtetA'; 'n_PtetA'; 'scale_mRFP1_doi_10_1093_'; 'K_C4HSL'; 'K_pLAMBDA_R_O12'; 'K_qsc'; 'alpha_pLAMBDA_R_O12'; 'alpha_qsc'; 'beta_pLAMBDA_R_O12'; 'beta_qsc'; 'n_C4HSL'; 'n_pLAMBDA_R_O12'; 'n_qsc'; 'scale_eYFP_doi_10_1002_b'; 'ALPHA_RBS_pTAC'; 'scale_GFPmut3b_doi_10_10'; 'scale_mCherry_doi_10_107'; 'scale_sfGFP_doi_10_1038_'; 'ALPHA_UTR1'; 'scale_deGFP_doi_10_1021_'; 'scale_mRFP1_doi_10_1186_'; 'ALPHA_pLUX_RBS'; 'alpha_pLUX_L'; 'scale_gfp_LVA_doi_10_103'; 'K_pRHL'; 'alpha_pRHL'; 'beta_pRHL'; 'n_pRHL'; 'scale_gfp_LVA_doi_10_107'; 'K_pLAS'; 'alpha_pLAS'; 'beta_pLAS'; 'n_pLAS'; 'scale_eYFP_doi_10_1038_n'; 'ALPHA_D103'; 'scale_mRFP1_doi_10_1093_'; 'ALPHA_RBS_pC_tetR_casset'; 'K_pTET_tetR_cassette'; 'alpha_pC_tetR_cassette'; 'alpha_pTET_tetR_cassette'; 'beta_pTET_tetR_cassette'; 'n_pTET_tetR_cassette'; 'ALPHA_RBS_A'; 'scale_sfGFP_doi_10_1038_'; 'ALPHA_RBS_psicA'; 'scale_mRFP1_doi_10_1038_'; 'ALPHA_RBS2'; 'K_pLUX_star_'; 'alpha_pLUX_star_'; 'beta_pLUX_star_'; 'n_pLUX_star_'; 'ALPHA_RBS_pET_29b'; 'K_placUV5'; 'alpha_placUV5'; 'beta_placUV5'; 'n_placUV5'; 'scale_mRFP1_doi_10_1038_'; 'scale_mRFP1_doi_10_1186_'; 'scale_GFP_uv_doi_10_1186'; 'ALPHA_RBS_pC_tetR'; 'K_pTET_operon'; 'alpha_pC_tetR'; 'alpha_pTET_operon'; 'beta_pTET_operon'; 'n_pTET_operon'; 'K_pT7lacO'; 'alpha_pT7lacO'; 'beta_pT7lacO'; 'n_pT7lacO'; 'scale_mRFP1_doi_10_1371_'; 'scale_sfGFP_doi_10_1073_'};
	lb = [0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01 0.03 0.002 0.01 0.01 0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01 0.03 0.01 0.03 0.002 0.01 0.01 0.03 0.002 0.01 0.002 0.002 0.002 0.002 0.002 0.01 0.002 0.002 0.002 0.002 0.01 0.002 0.001 0.002 0 1 0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01 0.03 0.03 0.001 0 1 0.01 0.03 0.001 0.002 0 1 0.01 0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.03 0.001 1 0.01 0.03 0.01 0.002 0.01 0.01 0.03 0.01 0.03 0.01 0.002 0.01 0.03 0.03 0.001 0.002 0 1 0.01 0.01 0.002 0.01 0.01 0.03 0.03 0.001 0.002 0.002 0 1 0.01 0.03 0.03 0.001 0.001 0.002 0.002 0 0 1 1 0.01 0.01 0.01 0.001 0.002 0.002 0 1 0.01 0.001 0.002 0 1 0.001 0.002 0 1 0.03 0.001 0.001 0.002 0 1 1 0.01 0.03 0.001 0.002 0 1 0.01 0.002 0.03 0.001 0.002 0 1 0.01 0.001 0.001 0.001 0.002 0.002 0 0 1 1 1 0.01 0.03 0.01 0.01 0.01 0.03 0.01 0.01 0.03 0.002 0.01 0.001 0.002 0 1 0.01 0.001 0.002 0 1 0.01 0.03 0.01 0.03 0.001 0.002 0.002 0 1 0.03 0.01 0.03 0.01 0.03 0.001 0.002 0 1 0.03 0.001 0.002 0 1 0.01 0.01 0.01 0.03 0.001 0.002 0.002 0 1 0.001 0.002 0 1 0.01 0.01];
	ub = [5 5 50 50 50 50 1 4 4 10 5 50 10 10 5 5 50 50 50 50 1 4 4 10 5 50 50 50 50 1 4 4 10 5 10 5 50 10 10 5 50 10 50 50 50 50 50 10 50 50 50 50 10 50 50 50 1 4 5 5 50 50 50 50 1 4 4 10 5 5 50 1 4 10 5 50 50 1 4 10 5 5 50 50 50 50 1 4 4 5 50 4 10 5 10 50 10 10 5 10 5 10 50 10 5 5 50 50 1 4 10 10 50 10 10 5 5 50 50 50 1 4 10 5 5 50 50 50 50 1 1 4 4 10 10 10 50 50 50 1 4 10 50 50 1 4 50 50 1 4 5 50 50 50 1 4 4 10 5 50 50 1 4 10 50 5 50 50 1 4 10 50 50 50 50 50 1 1 4 4 4 10 5 10 10 10 5 10 10 5 50 10 50 50 1 4 10 50 50 1 4 10 5 10 5 50 50 50 1 4 5 10 5 10 5 50 50 1 4 5 50 50 1 4 10 10 10 5 50 50 50 1 4 50 50 1 4 10 10];
	[CI_lb, CI_ub, x_opt, solution_ensemble] = fit_a_model(ensemble_size, ite_num, func_ite_num, @simultaneous_fitting_error, lb, ub);
	parameter_fileID = fopen('Plot_data/sim_CI.dat','w');
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
	SA_fileID = fopen('Plot_data/sensitivity_analysis.dat','w');
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
	error = zeros(765,1);
	current_base_index = 0;
	%==============
	L_arabinose = [0.003 0.005 0.05 0.5 1];
	GFPuv = [0.044 0.089 0.556 1 1];
	output = 1*(simf_doi_10_1128_AEM_00791_07_figure3a(x, L_arabinose) - GFPuv);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	L_arabinose = [0 0.001 0.003 0.008 0.025 0.075 0.2 0.7 2 6 20 50];
	mCherry = [0.015 0.015 0.015 0.02 0.06 0.48 0.84 0.93 0.96 0.98 0.98 1];
	output = 1*(simf_doi_10_1038_nature12148_figure16a(x, L_arabinose) - mCherry);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	L_arabinose = [0 0.001 0.013 0.133 1.33 13.3];
	AP = [0.004 0.004 0.021 0.412 0.823 1];
	output = 1*(simf_http___jb_asm_org_content_177_14_4121_short_figure4a(x, L_arabinose) - AP);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 1 2 3 4 5 6 7 8 9 10 12 13 17 18 20 25 35 60];
	luc = [0 0.001 0.001 0.001 0.002 0.004 0.007 0.011 0.014 0.032 0.036 0.054 0.089 0.643 0.893 1 1 1 1];
	output = 1*(simf_doi_10_1093_nar_25_6_1203_figure4a(x, aTc) - luc);
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
	output = 1*(simf_doi_10_1073_pnas_1316298111_figure2B(x, IPTG) - lacZ);
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
	L_arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2];
	GFPmut3b = [0 0.002 0.041 0.479 0.949 0.997 1 1 1];
	output = 1*(simf_doi_10_1093_nar_gku593_figureS6(x, L_arabinose) - GFPmut3b);
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
	Sal = [0 0.001 0.002 0.01 0.039 0.16 0.62];
	gfp = [0.005 0.023 0.115 0.308 0.462 0.615 0.769];
	output = 1.3*(simf_doi_10_1038_msb4100173_figure3_1(x, Sal) - gfp);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	L_arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.3];
	gfp = [0 0.001 0.002 0.019 0.308 0.769 1];
	output = 1*(simf_doi_10_1038_msb4100173_figure3_2(x, L_arabinose) - gfp);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0.001 0.005 0.01 0.02 0.03 0.05 0.1 0.2 1];
	cfp = [0.091 0.182 0.327 0.4 0.727 0.909 1 1 1];
	output = 1*(simf_doi_10_1534_genetics_112_147199_figure3(x, IPTG) - cfp);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.01 0.05 0.1 1 8];
	gfp_LAA = [0.125 0.125 0.417 0.625 0.75 0.75 0.683];
	output = 1.333*(simf_doi_10_1016_j_ymben_2015_03_009_figure4_b(x, IPTG) - gfp_LAA);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	AHL = [0 0 0 0 0 0 0 0 0];
	gfp_LAA = [0.183 0.183 0.208 0.167 0.25 0.5 0.833 1 1];
	output = 1*(simf_doi_10_1016_j_ymben_2015_03_009_figure4_c(x, AHL) - gfp_LAA);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	C_3OC12HSL = [0 0 0 0.001 0.001 0.005];
	GFPmut3b = [0.091 0.136 0.182 0.545 0.818 1];
	output = 1*(simf_doi_10_1016_j_ces_2012_12_016_figure4(x, C_3OC12HSL) - GFPmut3b);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0.005 0.01 0.02 0.04 0.1];
	gfp_LVA = [0.004 0.004 0.004 0.2 1];
	output = 1*(simf_doi_10_1038_msb_2009_75_figure3_1(x, IPTG) - gfp_LVA);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	AHL = [0 0 0 0 0 0 0.002];
	GFPmut3b = [0.018 0.018 0.018 0.212 0.939 1 0.97];
	output = 1*(simf_doi_10_1016_j_bios_2012_08_011_figure2d(x, AHL) - GFPmut3b);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	AHL = [0 0 0 0 0.001 0.01 0.1];
	mRFP1 = [0.024 0.081 0.372 0.763 0.874 0.889 0.891];
	output = 1.123*(simf_doi_10_1093_nar_gku964_figure3A_1(x, AHL) - mRFP1);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	AHL = [0 0 0 0 0.001 0.01 0.1];
	mRFP1 = [0.015 0.052 0.295 0.755 0.911 0.931 0.934];
	output = 1.071*(simf_doi_10_1093_nar_gku964_figure3A_2(x, AHL) - mRFP1);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	AHL = [0 0 0 0 0.001 0.01 0.1];
	mRFP1 = [0.006 0.01 0.05 0.319 0.82 0.98 1];
	output = 1*(simf_doi_10_1093_nar_gku964_figure3A_3(x, AHL) - mRFP1);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.01 0.1 1 10];
	gfp_uv = [0.14 0.14 0.76 0.96 1];
	output = 1*(simf_doi_10_1021_bp010141k_figure3_1(x, IPTG) - gfp_uv);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) = 1*(simf_doi_10_1021_sb300055e_figure4a_1(x) - 1);
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) = 2.333*(simf_doi_10_1021_sb300055e_figure4a_2(x) - 0.429);
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) = 8.4*(simf_doi_10_1021_sb300055e_figure4a_3(x) - 0.119);
	current_base_index = current_base_index + 1;
	%==============
	aTc = [0.3 0.7 1 2 3 4 5 10 25 100];
	gfpmut3_1 = [0.039 0.196 0.294 0.51 0.667 0.843 0.902 0.902 0.98 1];
	output = 1*(simf_doi_10_1038_msb4100081_figure3b(x, aTc) - gfpmut3_1);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0 0 0 0.001 0.003 0.01 0.04 0.2 1];
	gfp = [0.211 0.211 0.237 0.237 0.263 0.447 1 0.921 0.868 0.842];
	output = 1*(simf_doi_10_1186_1754_1611_6_9_figure2_b(x, IPTG) - gfp);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	L_arabinose = [0.03 0.08 0.2 0.5 1.2 3 8 20 50 125 300];
	gfpmut2 = [0.014 0.083 0.181 0.361 0.5 0.694 0.778 0.889 1 0.917 0.778];
	output = 1*(simf_doi_10_1186_1752_0509_5_111_figure3_a(x, L_arabinose) - gfpmut2);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	L_arabinose = [0.03 0.08 0.2 0.5 1.2 3 8 20 50 125 300];
	gfpmut2 = [0.014 0.014 0.014 0.017 0.021 0.028 0.061 0.167 0.389 0.597 0.556];
	output = 1.674*(simf_doi_10_1186_1752_0509_5_111_figure3_b(x, L_arabinose) - gfpmut2);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	AHL = [0 0 0 0 0 0.001 0.004 0.01 0.04 0.1 0.4 1];
	eYFP = [0.001 0.029 0.114 0.286 0.429 0.643 0.743 0.829 0.857 0.9 0.929 1];
	output = 1*(simf_dx_doi_org_10_1021_sb400152n_figure2(x, AHL) - eYFP);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.002 0.01 0.02 0.1 1];
	eYFP = [1 0.867 0.733 0.667 0.333 0.2 0.133];
	output = 1*(simf_doi_10_1073_pnas_252535999_figure3(x, IPTG) - eYFP);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0.001 0.01 0.02 0.03 0.04 0.06 0.1 0.3 0.6 1 3 6 10];
	GFPmut3b = [0.011 0.011 0.011 0.011 0.022 0.056 0.111 0.444 0.722 0.889 1 1 1];
	output = 1*(simf_doi_10_1038_35002131_figure5a(x, IPTG) - GFPmut3b);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0.1 1 10 100 1000];
	glk = [1 0.786 0.357 0.071 0.036];
	output = 1*(simf_doi_10_1016_j_ymben_2012_08_006_figure7(x, aTc) - glk);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0.001 0.01 0.02 0.03 0.05 0.1 0.2 0.5 1];
	lacZ = [0.017 0.042 0.1 0.15 0.3 0.6 0.85 0.975 1];
	output = 1*(simf_doi_10_1006_plas_2000_1477_figure5(x, IPTG) - lacZ);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	AHL = [0 0 0 0.001 0.005 0.01 0.1 1 10 100];
	GFPmut3b = [0 0.013 0.108 0.343 0.873 0.954 0.999 1 1 1];
	output = 1*(simf_doi_10_1038_nbt1413_figure3(x, AHL) - GFPmut3b);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	C_3OC12HSL = [0 0.001 0.01 0.1 1 10];
	GFPmut3b = [0 0.225 0.693 0.874 0.898 0.9];
	output = 1.111*(simf_doi_10_1093_nar_gkt052_figure3b(x, C_3OC12HSL) - GFPmut3b);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	C_3OC12HSL = [0 0.001 0.01 0.1 1 10];
	GFPmut3b = [0 0.048 0.446 0.927 0.995 1];
	output = 1*(simf_doi_10_1093_nar_gkt052_figure3c(x, C_3OC12HSL) - GFPmut3b);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	C_3OC12HSL = [0 0.001 0.01 0.1 1 10];
	GFPmut3b = [0 0.001 0.15 0.299 0.3 0.3];
	output = 3.332*(simf_doi_10_1093_nar_gkt052_figure3d(x, C_3OC12HSL) - GFPmut3b);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	Mg2_plus__plus_ = [0 0.001 0.005 0.05 0.5 5 50];
	gfp = [0.054 0.05 0.047 0.043 0.028 0.017 0.01];
	output = 18.413*(simf_doi_10_1038_nmeth_2205_figureS1_1(x, Mg2_plus__plus_) - gfp);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	L_arabinose = [0.017 0.027 0.11 0.27 0.67 1.3 6.7];
	gfp = [0.007 0.009 0.017 0.06 0.121 0.115 0.103];
	output = 8.286*(simf_doi_10_1038_nmeth_2205_figureS1_2(x, L_arabinose) - gfp);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0.25 1 4 16 64 240];
	gfp = [0.002 0.002 0.003 0.009 0.036 0.043];
	output = 23.2*(simf_doi_10_1038_nmeth_2205_figureS1_3(x, aTc) - gfp);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0 0.001 0.01 0.1 1];
	gfp = [0.009 0.034 0.069 0.086 0.095 0.121];
	output = 8.286*(simf_doi_10_1038_nmeth_2205_figureS1_41(x, IPTG) - gfp);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0 0.001 0.01 0.1 1];
	gfp = [0.129 0.129 0.129 0.345 0.905 1];
	output = 1*(simf_doi_10_1038_nmeth_2205_figureS1_42(x, IPTG) - gfp);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1 0.15 0.2 1];
	eYFP = [0.001 0.002 0.003 0.005 0.008 0.016 0.024 0.032 0.06 0.1 0.14 0.2];
	output = 5*(simf_doi_10_1038_nchembio_1411_figure4_1(x, IPTG) - eYFP);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1 0.15 0.2 1];
	eYFP = [1 1 0.8 0.2 0.04 0.028 0.02 0.012 0.01 0.008 0.008 0.008];
	output = 1*(simf_doi_10_1038_nchembio_1411_figure4_2(x, IPTG) - eYFP);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 15];
	mRFP1 = [0.08 1];
	output = 1*(simf_doi_10_1093_nar_gks583_figure3(x, aTc) - mRFP1);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	C4HSL = [0 0 0 0 0.001 0.004 0.01 0.04];
	eYFP = [1 1 1 0.667 0.433 0.167 0.033 0.013];
	output = 1*(simf_doi_10_1002_bit_20371_figure3(x, C4HSL) - eYFP);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0.001 0.005 0.01 0.015 0.025 0.04 0.05 0.1 0.2 0.4];
	GFPmut3b = [0.01 0.1 0.27 0.35 0.5 0.7 0.9 0.9 1 1];
	output = 1*(simf_doi_10_1073_pnas_1314114111_figure5_b(x, IPTG) - GFPmut3b);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [9.2 13.8 18.4 23 27.6 32.2 46 92 138];
	mCherry = [0.111 0.333 0.533 0.667 0.756 0.867 0.944 1 1];
	output = 1*(simf_doi_10_1073_pnas_1314114111_figure5_d(x, aTc) - mCherry);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	AHL = [0 0 0 0 0];
	sfGFP = [0.06 0.06 0.09 0.15 1];
	output = 1*(simf_doi_10_1038_msb_2013_58_figureS14(x, AHL) - sfGFP);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0.46 4.6 46 460 4600];
	deGFP = [0.025 0.025 0.075 0.25 0.8];
	output = 1.25*(simf_doi_10_1021_sb400131a_figure4c_1(x, aTc) - deGFP);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0 0.001 0.01 0.1 1];
	deGFP = [0.075 0.075 0.075 0.25 0.875 1];
	output = 1*(simf_doi_10_1021_sb400131a_figure4c_2(x, IPTG) - deGFP);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	AHL = [0 0 0 0 0 0 0 0 0 0 0];
	mRFP1 = [0.006 0.086 0.126 0.343 0.571 0.657 0.857 0.914 0.943 0.971 1];
	output = 1*(simf_doi_10_1186_1471_2105_13_S4_S11_figure2_1(x, AHL) - mRFP1);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	AHL = [0 0 0 0 0 0 0 0 0 0 0];
	mRFP1 = [0.001 0.011 0.017 0.034 0.057 0.069 0.086 0.097 0.103 0.109 0.114];
	output = 8.75*(simf_doi_10_1186_1471_2105_13_S4_S11_figure2_2(x, AHL) - mRFP1);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	AHL = [0 0 0 0 0 0 0 0 0 0 0];
	mRFP1 = [0.001 0.001 0.001 0.002 0.007 0.014 0.026 0.03 0.033 0.034 0.035];
	output = 28.226*(simf_doi_10_1186_1471_2105_13_S4_S11_figure2_3(x, AHL) - mRFP1);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	AHL = [0.1 0.3 0.5 1 2 4 10 100 1000];
	gfp_LVA = [1 0.947 0.716 0.579 0.158 0.158 0.105 0.053 0.011];
	output = 1*(simf_doi_10_1038_nature03461_figure2b_1(x, AHL) - gfp_LVA);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	AHL = [0.1 0.3 0.7 1.2 2 5 10 20 30 50 100 500 1000];
	gfp_LVA = [0.947 0.916 0.937 0.916 0.905 0.842 0.684 0.526 0.368 0.284 0.242 0.137 0.126];
	output = 1.056*(simf_doi_10_1038_nature03461_figure2b_2(x, AHL) - gfp_LVA);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	AHL = [0.2 0.5 1 2 5 10 50 100 200 300 400 600 1000 2000 3000];
	gfp_LVA = [0.842 0.832 0.8 0.789 0.737 0.695 0.632 0.611 0.526 0.421 0.316 0.263 0.232 0.179 0.126];
	output = 1.188*(simf_doi_10_1038_nature03461_figure2b_3(x, AHL) - gfp_LVA);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	C4HSL = [1 10 100 200 500 1000 5000 10000 50000 100000];
	gfp_LVA = [0.001 0.001 0.002 0.003 0.004 0.025 0.15 0.4 0.5 0.6];
	output = 1.667*(simf_doi_10_1073_pnas_0704256104_figureS1b_1(x, C4HSL) - gfp_LVA);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	C_3OC12HSL = [0.1 1 5 8 10 20 50 80 100 1000];
	gfp_LVA = [0.001 0.001 0.001 0.004 0.025 0.045 0.25 0.45 1 1];
	output = 1*(simf_doi_10_1073_pnas_0704256104_figureS1b_2(x, C_3OC12HSL) - gfp_LVA);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	L_arabinose = [0 0.001 0.015 0.05 0.5 5 10];
	eYFP = [0.002 0.002 0.057 0.628 0.999 1 1];
	output = 1*(simf_doi_10_1038_nature09565_figure1C(x, L_arabinose) - eYFP);
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
	IPTG = [0 0.002 0.01 0.05 0.1 0.25 1];
	mRFP1 = [0.054 0.107 0.179 0.571 0.693 0.893 1];
	output = 1*(simf_doi_10_1093_nar_gks597_figureS3_left(x, IPTG) - mRFP1);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 0.05 0.1 0.25 0.5 1 10 50];
	mRFP1 = [0.032 0.079 0.129 0.193 0.214 0.25 0.279 0.357];
	output = 2.8*(simf_doi_10_1093_nar_gks597_figureS3_right(x, aTc) - mRFP1);
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
	L_arabinose = [0.1 1 2 5 7 10 12.5 25 37.5 50 62.5];
	sfGFP = [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425];
	output = 2.353*(simf_doi_10_1038_nbt_2401_figureS13(x, L_arabinose) - sfGFP);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	L_arabinose = [0 0 0.002 0.01 0.04 0.16 1 4 25];
	mRFP1 = [0.002 0.002 0.002 0.01 0.137 0.487 0.557 0.558 0.558];
	output = 1.791*(simf_doi_10_1038_nature11516_figureS8A(x, L_arabinose) - mRFP1);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0 0 0.001 0.004 0.025 0.1 0.5 2.5];
	mRFP1 = [0.016 0.016 0.016 0.018 0.04 0.273 0.525 0.588 0.593];
	output = 1.686*(simf_doi_10_1038_nature11516_figureS8B(x, IPTG) - mRFP1);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	AHL = [0 0 0 0 0 0 0 0 0];
	mRFP1 = [0.102 0.102 0.128 0.329 0.89 0.988 0.999 1 1];
	output = 1*(simf_doi_10_1038_nature11516_figureS8C(x, AHL) - mRFP1);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	AHL = [0 0 0 0 0 0 0 0 0];
	mRFP1 = [0.02 0.02 0.026 0.072 0.241 0.278 0.282 0.283 0.283];
	output = 3.536*(simf_doi_10_1038_nature11516_figureS8D(x, AHL) - mRFP1);
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
	%==============
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	mRFP1 = [0.214 0.679 0.764 0.821 0.821 0.836];
	output = 1.197*(simf_doi_10_1186_1754_1611_5_12_Supplement_page1(x, IPTG) - mRFP1);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	GFP_uv = [0.051 0.492 0.61 0.678 0.746 0.847];
	output = 1.18*(simf_doi_10_1186_1754_1611_5_12_Supplement_page2(x, IPTG) - GFP_uv);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	mRFP1 = [0.296 0.914 0.964 1 1 1];
	output = 1*(simf_doi_10_1186_1754_1611_5_12_Supplement_page3(x, IPTG) - mRFP1);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	mRFP1 = [0.129 0.407 0.443 0.5 0.5 0.5];
	output = 2*(simf_doi_10_1186_1754_1611_5_12_Supplement_page4(x, IPTG) - mRFP1);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 0.013 0.025 0.05 0.1 0.2];
	mRFP1 = [0.007 0.429 0.714 0.929 0.986 0.964];
	output = 1.014*(simf_doi_10_1186_1754_1611_5_12_Supplement_page5(x, aTc) - mRFP1);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 0.013 0.025 0.05 0.1 0.2];
	GFP_uv = [0.002 0.288 0.542 0.797 0.932 0.966];
	output = 1.035*(simf_doi_10_1186_1754_1611_5_12_Supplement_page6(x, aTc) - GFP_uv);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 0.013 0.025 0.05 0.1 0.2];
	mRFP1 = [0.007 0.271 0.643 0.807 0.929 0.893];
	output = 1.077*(simf_doi_10_1186_1754_1611_5_12_Supplement_page7(x, aTc) - mRFP1);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 0.013 0.025 0.05 0.1 0.2];
	mRFP1 = [0.007 0.429 0.557 0.571 0.557 0.543];
	output = 1.75*(simf_doi_10_1186_1754_1611_5_12_Supplement_page8(x, aTc) - mRFP1);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	mRFP1 = [0.021 0.157 0.2 0.286 0.3 0.279];
	output = 3.333*(simf_doi_10_1186_1754_1611_5_12_Supplement_page17(x, IPTG) - mRFP1);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	GFP_uv = [0.017 0.09 0.115 0.153 0.161 0.183];
	output = 5.463*(simf_doi_10_1186_1754_1611_5_12_Supplement_page18(x, IPTG) - GFP_uv);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	mRFP1 = [0.036 0.25 0.343 0.486 0.5 0.536];
	output = 1.867*(simf_doi_10_1186_1754_1611_5_12_Supplement_page19(x, IPTG) - mRFP1);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	mRFP1 = [0.014 0.086 0.107 0.143 0.143 0.143];
	output = 7*(simf_doi_10_1186_1754_1611_5_12_Supplement_page20(x, IPTG) - mRFP1);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	mRFP1 = [0.014 0.286 0.421 0.557 0.614 0.571];
	output = 1.628*(simf_doi_10_1186_1754_1611_5_12_Supplement_page21(x, IPTG) - mRFP1);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	GFP_uv = [0.008 0.064 0.105 0.203 0.242 0.339];
	output = 2.95*(simf_doi_10_1186_1754_1611_5_12_Supplement_page22(x, IPTG) - GFP_uv);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	mRFP1 = [0.014 0.25 0.4 0.593 0.65 0.664];
	output = 1.505*(simf_doi_10_1186_1754_1611_5_12_Supplement_page23(x, IPTG) - mRFP1);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	mRFP1 = [0.014 0.229 0.307 0.393 0.393 0.414];
	output = 2.414*(simf_doi_10_1186_1754_1611_5_12_Supplement_page24(x, IPTG) - mRFP1);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	mRFP1 = [0.129 0.357 0.729 0.843 0.821 0.714];
	output = 1.186*(simf_doi_10_1186_1754_1611_5_12_Supplement_page25(x, IPTG) - mRFP1);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	GFP_uv = [0.034 0.619 0.805 1 0.983 0.881];
	output = 1*(simf_doi_10_1186_1754_1611_5_12_Supplement_page26(x, IPTG) - GFP_uv);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.004 0.02 0.1 0.5];
	mRFP1 = [0.093 0.321 0.714 0.943 0.943 0.957];
	output = 1.045*(simf_doi_10_1186_1754_1611_5_12_Supplement_page27(x, IPTG) - mRFP1);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.004 0.02 0.1 0.5];
	mRFP1 = [0.093 0.321 0.629 0.679 0.586 0.486];
	output = 1.474*(simf_doi_10_1186_1754_1611_5_12_Supplement_page28(x, IPTG) - mRFP1);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	L_arabinose = [0 0.001 0.01 0.1 1 20];
	mRFP1 = [0 0.004 0.029 0.086 0.143 0.179];
	output = 5.6*(simf_doi_10_1186_1754_1611_5_12_Supplement_page29(x, L_arabinose) - mRFP1);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	L_arabinose = [0 0.006 0.013 0.05 0.1 0.5];
	GFP_uv = [0.008 0.017 0.068 0.153 0.305 0.661];
	output = 1.513*(simf_doi_10_1186_1754_1611_5_12_Supplement_page30(x, L_arabinose) - GFP_uv);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	L_arabinose = [0 0.001 0.01 0.1 1 20];
	mRFP1 = [0 0.007 0.05 0.15 0.329 0.464];
	output = 2.154*(simf_doi_10_1186_1754_1611_5_12_Supplement_page31(x, L_arabinose) - mRFP1);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	L_arabinose = [0 0.001 0.01 0.1 1 20];
	mRFP1 = [0.004 0.004 0.014 0.032 0.064 0.079];
	output = 12.727*(simf_doi_10_1186_1754_1611_5_12_Supplement_page32(x, L_arabinose) - mRFP1);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0.005 0.01 0.02 0.03 0.05 0.1 0.2 0.3 0.8 1];
	mRFP1 = [0.004 0.038 0.192 0.423 0.75 0.904 0.904 0.908 0.981 1];
	output = 1*(simf_doi_10_1371_journal_pone_0039407_figure6B(x, IPTG) - mRFP1);
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) = 2.372*(simf_doi_10_1073_pnas_1301301110_sd03_xls_1(x) - 0.422);
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) = 1*(simf_doi_10_1073_pnas_1301301110_sd03_xls_2(x) - 1);
	current_base_index = current_base_index + 1;
end

function print_plotting_sim_data(sim_x_opt)
	sim_pair_fileID = fopen('Plot_data/sim_pair.dat','w');
	json_fileID = fopen('JSON/simulation_data.json','w');
	%%%%%%%%%%%%%
	L_arabinose = [0.003 0.005 0.05 0.5 1];
	GFPuv = [0.044 0.089 0.556 1 1];
	output_sim = simf_doi_10_1128_AEM_00791_07_figure3a(sim_x_opt,L_arabinose);
	sim_pair_matrix = [output_sim; GFPuv];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1128/AEM.00791-07', 'figure3a', L_arabinose, output_sim, {'pC = 20*{alpha_pC}'; 'araC = {ALPHA_RBS_pC}*pC*{K_L_arabinose}^{n_L_arabinose}/({K_L_arabinose}^{n_L_arabinose} + L_arabinose^{n_L_arabinose})'; 'pBAD = 20*({beta_pBAD} + ({alpha_pBAD} - {beta_pBAD})*{K_pBAD}^{n_pBAD}/({K_pBAD}^{n_pBAD} + araC^{n_pBAD}))'; 'GFPuv = {scale_GFPuv_doi_10_1128_AEM_00791_07}*{ALPHA_RBS_pBAD}*pBAD'});
	%%%%%%%%%%%%%
	L_arabinose = [0 0.001 0.003 0.008 0.025 0.075 0.2 0.7 2 6 20 50];
	mCherry = [0.015 0.015 0.015 0.02 0.06 0.48 0.84 0.93 0.96 0.98 0.98 1];
	output_sim = simf_doi_10_1038_nature12148_figure16a(sim_x_opt,L_arabinose);
	sim_pair_matrix = [output_sim; mCherry];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1038/nature12148', 'figure16a', L_arabinose, output_sim, {'pLlacO_1 = 5*{alpha_pLlacO_1}'; 'araC = {ALPHA_BBa_B0030}*pLlacO_1*{K_L_arabinose}^{n_L_arabinose}/({K_L_arabinose}^{n_L_arabinose} + L_arabinose^{n_L_arabinose})'; 'pBAD = 15*({beta_pBAD} + ({alpha_pBAD} - {beta_pBAD})*{K_pBAD}^{n_pBAD}/({K_pBAD}^{n_pBAD} + araC^{n_pBAD}))'; 'mCherry = {scale_mCherry_doi_10_1038_nature12148}*{ALPHA_BBa_B0030}*pBAD'});
	%%%%%%%%%%%%%
	L_arabinose = [0 0.001 0.013 0.133 1.33 13.3];
	AP = [0.004 0.004 0.021 0.412 0.823 1];
	output_sim = simf_http___jb_asm_org_content_177_14_4121_short_figure4a(sim_x_opt,L_arabinose);
	sim_pair_matrix = [output_sim; AP];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'http://jb.asm.org/content/177/14/4121.short', 'figure4a', L_arabinose, output_sim, {'pC = 20*{alpha_pC}'; 'araC = {ALPHA_RBS_pC}*pC*{K_L_arabinose}^{n_L_arabinose}/({K_L_arabinose}^{n_L_arabinose} + L_arabinose^{n_L_arabinose})'; 'pBAD = 20*({beta_pBAD} + ({alpha_pBAD} - {beta_pBAD})*{K_pBAD}^{n_pBAD}/({K_pBAD}^{n_pBAD} + araC^{n_pBAD}))'; 'AP = {scale_AP_http___jb_asm_org_content_177_14_4121_short}*{ALPHA_RBS_pBAD}*pBAD'});
	%%%%%%%%%%%%%
	aTc = [0 1 2 3 4 5 6 7 8 9 10 12 13 17 18 20 25 35 60];
	luc = [0 0.001 0.001 0.001 0.002 0.004 0.007 0.011 0.014 0.032 0.036 0.054 0.089 0.643 0.893 1 1 1 1];
	output_sim = simf_doi_10_1093_nar_25_6_1203_figure4a(sim_x_opt,aTc);
	sim_pair_matrix = [output_sim; luc];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1093/nar/25.6.1203', 'figure4a', aTc, output_sim, {'pN25 = 1*{alpha_pN25}'; 'tetR = {ALPHA_RBS_pN25}*pN25*{K_aTc}^{n_aTc}/({K_aTc}^{n_aTc} + aTc^{n_aTc})'; 'pLtetO_1 = 15*({beta_pLtetO_1} + ({alpha_pLtetO_1} - {beta_pLtetO_1})*{K_pLtetO_1}^{n_pLtetO_1}/({K_pLtetO_1}^{n_pLtetO_1} + tetR^{n_pLtetO_1}))'; 'luc = {scale_luc_doi_10_1093_nar_25_6_1203}*{ALPHA_RBSII}*pLtetO_1'});
	%%%%%%%%%%%%%
	IPTG = [0.001 0.01 0.03 0.05 0.075 0.1 0.15 0.2 0.3 0.5 1];
	eYFP = [0.015 0.014 0.023 0.062 0.098 0.187 0.239 0.374 0.572 0.869 1];
	output_sim = simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; eYFP];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'http://dspace.mit.edu/handle/1721.1/8228', 'figure4-6', IPTG, output_sim, {'placIq = 10*{alpha_placIq}'; 'lacI = {ALPHA_RBS_placI}*placIq*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'pLAC = 10*({beta_pLAC} + ({alpha_pLAC} - {beta_pLAC})*{K_pLAC}^{n_pLAC}/({K_pLAC}^{n_pLAC} + lacI^{n_pLAC}))'; 'eYFP = {scale_eYFP_http___dspace_mit_edu_handle_1721_1_8228}*{ALPHA_RBSII}*pLAC'});
	%%%%%%%%%%%%%
	aTc = [9.258 13.887 20.831 32.403 46.29 83.322 162.015 231.45 370.32 648.06 925.8 1388.7];
	eYFP = [0.003 0.003 0.003 0.003 0.005 0.01 0.05 0.167 0.333 0.667 0.833 1];
	output_sim = simf_doi_10_1073pnas_0408507102_figure2A(sim_x_opt,aTc);
	sim_pair_matrix = [output_sim; eYFP];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1073pnas.0408507102', 'figure2A', aTc, output_sim, {'placIq = 10*{alpha_placIq}'; 'tetR = {ALPHA_RBS_placI}*placIq*{K_aTc}^{n_aTc}/({K_aTc}^{n_aTc} + aTc^{n_aTc})'; 'pLtetO_1 = 10*({beta_pLtetO_1} + ({alpha_pLtetO_1} - {beta_pLtetO_1})*{K_pLtetO_1}^{n_pLtetO_1}/({K_pLtetO_1}^{n_pLtetO_1} + tetR^{n_pLtetO_1}))'; 'eYFP = {scale_eYFP_doi_10_1073pnas_0408507102}*{ALPHA_RBS_Unknown_Hooshangi}*pLtetO_1'});
	%%%%%%%%%%%%%
	IPTG = [0 0.005 0.01 0.02 0.05 0.1 0.2 1];
	lacZ = [0.001 0.003 0.043 0.5 0.984 0.999 1 1];
	output_sim = simf_doi_10_1073_pnas_1316298111_figure2B(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; lacZ];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1073/pnas.1316298111', 'figure2B', IPTG, output_sim, {'placI = 1*{alpha_placI}'; 'lacI = {ALPHA_RBS_placI}*placI*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'pLAC = 1*({beta_pLAC} + ({alpha_pLAC} - {beta_pLAC})*{K_pLAC}^{n_pLAC}/({K_pLAC}^{n_pLAC} + lacI^{n_pLAC}))'; 'lacZ = {scale_lacZ_doi_10_1073_pnas_1316298111}*{ALPHA_RBS_pLAC}*pLAC'});
	%%%%%%%%%%%%%
	IPTG = [0 0.005 0.01 0.02 0.05 0.1 0.2 1];
	lacZ = [0.001 0.003 0.043 0.5 0.984 0.999 1 1];
	output_sim = simf_doi_10_1073_pnas_0606717104_figure4_6(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; lacZ];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1073/pnas.0606717104', 'figure4-6', IPTG, output_sim, {'placI = 1*{alpha_placI}'; 'lacI = {ALPHA_RBS_placI}*placI*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'pLAC = 1*({beta_pLAC} + ({alpha_pLAC} - {beta_pLAC})*{K_pLAC}^{n_pLAC}/({K_pLAC}^{n_pLAC} + lacI^{n_pLAC}))'; 'lacZ = {scale_lacZ_doi_10_1073_pnas_0606717104}*{ALPHA_RBS_pLAC}*pLAC'});
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
	L_arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2];
	GFPmut3b = [0 0.002 0.041 0.479 0.949 0.997 1 1 1];
	output_sim = simf_doi_10_1093_nar_gku593_figureS6(sim_x_opt,L_arabinose);
	sim_pair_matrix = [output_sim; GFPmut3b];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1093/nar/gku593', 'figureS6', L_arabinose, output_sim, {'pC = 10*{alpha_pC}'; 'araC = {ALPHA_RBS_pC}*pC*{K_L_arabinose}^{n_L_arabinose}/({K_L_arabinose}^{n_L_arabinose} + L_arabinose^{n_L_arabinose})'; 'pBAD = 10*({beta_pBAD} + ({alpha_pBAD} - {beta_pBAD})*{K_pBAD}^{n_pBAD}/({K_pBAD}^{n_pBAD} + araC^{n_pBAD}))'; 'GFPmut3b = {scale_GFPmut3b_doi_10_1093_nar_gku593}*{ALPHA_BBa_B0030}*pBAD'});
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
	print_simulated_data_to_json(json_fileID, 'doi:10.1093/nar/gku1388', 'figure2B-1', aTc, output_sim, {'BBa_J23117 = 10*{alpha_BBa_J23117}'; 'tetR = {ALPHA_BBa_B0030}*BBa_J23117*{K_aTc}^{n_aTc}/({K_aTc}^{n_aTc} + aTc^{n_aTc})'; 'pTET_star_ = 10*({beta_pTET_star_} + ({alpha_pTET_star_} - {beta_pTET_star_})*{K_pTET_star_}^{n_pTET_star_}/({K_pTET_star_}^{n_pTET_star_} + tetR^{n_pTET_star_}))'; 'GFPmut3b = {scale_GFPmut3b_doi_10_1093_nar_gku1388}*{ALPHA_BBa_B0030}*pTET_star_'});
	%%%%%%%%%%%%%
	aTc = [10 25 50 100 200];
	GFPmut3b = [0.001 0.022 0.176 0.411 0.461];
	output_sim = simf_doi_10_1093_nar_gku1388_figure2B_2(sim_x_opt,aTc);
	sim_pair_matrix = [output_sim; GFPmut3b];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1093/nar/gku1388', 'figure2B-2', aTc, output_sim, {'BBa_J23114 = 10*{alpha_BBa_J23114}'; 'tetR = {ALPHA_BBa_B0030}*BBa_J23114*{K_aTc}^{n_aTc}/({K_aTc}^{n_aTc} + aTc^{n_aTc})'; 'pTET_star_ = 10*({beta_pTET_star_} + ({alpha_pTET_star_} - {beta_pTET_star_})*{K_pTET_star_}^{n_pTET_star_}/({K_pTET_star_}^{n_pTET_star_} + tetR^{n_pTET_star_}))'; 'GFPmut3b = {scale_GFPmut3b_doi_10_1093_nar_gku1388}*{ALPHA_BBa_B0030}*pTET_star_'});
	%%%%%%%%%%%%%
	Sal = [0 0.001 0.002 0.01 0.039 0.16 0.62];
	gfp = [0.005 0.023 0.115 0.308 0.462 0.615 0.769];
	output_sim = simf_doi_10_1038_msb4100173_figure3_1(sim_x_opt,Sal);
	sim_pair_matrix = [output_sim; gfp];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1038/msb4100173', 'figure3-1', Sal, output_sim, {'pCONST_UNKNOWN = 15*{alpha_pCONST_UNKNOWN}'; 'nahR = {ALPHA_RBS_Unknown}*pCONST_UNKNOWN'; 'Sal_nahR = nahR*{K_Sal}^{n_Sal}/({K_Sal}^{n_Sal} + Sal^{n_Sal})'; 'pSAL = 15*({beta_pSAL} + ({alpha_pSAL} - {beta_pSAL})*Sal_nahR^{n_pSAL}/({K_pSAL}^{n_pSAL} + Sal_nahR^{n_pSAL}))'; 'gfp = {scale_gfp_doi_10_1038_msb4100173}*{ALPHA_SAL_RBS}*pSAL'});
	%%%%%%%%%%%%%
	L_arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.3];
	gfp = [0 0.001 0.002 0.019 0.308 0.769 1];
	output_sim = simf_doi_10_1038_msb4100173_figure3_2(sim_x_opt,L_arabinose);
	sim_pair_matrix = [output_sim; gfp];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1038/msb4100173', 'figure3-2', L_arabinose, output_sim, {'pC = 15*{alpha_pC}'; 'araC = {ALPHA_RBS_pC}*pC*{K_L_arabinose}^{n_L_arabinose}/({K_L_arabinose}^{n_L_arabinose} + L_arabinose^{n_L_arabinose})'; 'pBAD = 15*({beta_pBAD} + ({alpha_pBAD} - {beta_pBAD})*{K_pBAD}^{n_pBAD}/({K_pBAD}^{n_pBAD} + araC^{n_pBAD}))'; 'gfp = {scale_gfp_doi_10_1038_msb4100173}*{ALPHA_parent_RBS}*pBAD'});
	%%%%%%%%%%%%%
	IPTG = [0.001 0.005 0.01 0.02 0.03 0.05 0.1 0.2 1];
	cfp = [0.091 0.182 0.327 0.4 0.727 0.909 1 1 1];
	output_sim = simf_doi_10_1534_genetics_112_147199_figure3(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; cfp];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1534/genetics.112.147199', 'figure3', IPTG, output_sim, {'placI = 1*{alpha_placI}'; 'lacI = {ALPHA_RBS_placI}*placI*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'pLAC = 1*({beta_pLAC} + ({alpha_pLAC} - {beta_pLAC})*{K_pLAC}^{n_pLAC}/({K_pLAC}^{n_pLAC} + lacI^{n_pLAC}))'; 'lacZ = {ALPHA_RBS_pLAC}*pLAC'; 'pLlacO_1 = 1*({beta_pLlacO_1} + ({alpha_pLlacO_1} - {beta_pLlacO_1})*{K_pLlacO_1}^{n_pLlacO_1}/({K_pLlacO_1}^{n_pLlacO_1} + lacI^{n_pLlacO_1}))'; 'cfp = {scale_cfp_doi_10_1534_genetics_112_147199}*{ALPHA_T710}*pLlacO_1'});
	%%%%%%%%%%%%%
	IPTG = [0 0.001 0.01 0.05 0.1 1 8];
	gfp_LAA = [0.125 0.125 0.417 0.625 0.75 0.75 0.683];
	output_sim = simf_doi_10_1016_j_ymben_2015_03_009_figure4_b(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; gfp_LAA];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1016/j.ymben.2015.03.009', 'figure4-b', IPTG, output_sim, {'placIq = 15*{alpha_placIq}'; 'lacI = {ALPHA_RBS_placI}*placIq*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'pTrc = 15*({beta_pTrc} + ({alpha_pTrc} - {beta_pTrc})*{K_pTrc}^{n_pTrc}/({K_pTrc}^{n_pTrc} + lacI^{n_pTrc}))'; 'gfp_LAA = {scale_gfp_LAA_doi_10_1016_j_ymben_2015_03_009}*{ALPHA_g10}*pTrc'});
	%%%%%%%%%%%%%
	AHL = [0 0 0 0 0 0 0 0 0];
	gfp_LAA = [0.183 0.183 0.208 0.167 0.25 0.5 0.833 1 1];
	output_sim = simf_doi_10_1016_j_ymben_2015_03_009_figure4_c(sim_x_opt,AHL);
	sim_pair_matrix = [output_sim; gfp_LAA];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1016/j.ymben.2015.03.009', 'figure4-c', AHL, output_sim, {'pluxR = 15*{alpha_pluxR}'; 'luxR = {ALPHA_RBS_pluxR}*pluxR'; 'AHL_luxR = luxR*{K_AHL}^{n_AHL}/({K_AHL}^{n_AHL} + AHL^{n_AHL})'; 'pLUX = 15*({beta_pLUX} + ({alpha_pLUX} - {beta_pLUX})*AHL_luxR^{n_pLUX}/({K_pLUX}^{n_pLUX} + AHL_luxR^{n_pLUX}))'; 'gfp_LAA = {scale_gfp_LAA_doi_10_1016_j_ymben_2015_03_009}*{ALPHA_RBS_pLUX}*pLUX'});
	%%%%%%%%%%%%%
	C_3OC12HSL = [0 0 0 0.001 0.001 0.005];
	GFPmut3b = [0.091 0.136 0.182 0.545 0.818 1];
	output_sim = simf_doi_10_1016_j_ces_2012_12_016_figure4(sim_x_opt,C_3OC12HSL);
	sim_pair_matrix = [output_sim; GFPmut3b];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1016/j.ces.2012.12.016', 'figure4', C_3OC12HSL, output_sim, {'pLtetO_1 = 200*{alpha_pLtetO_1}'; 'lasR = {ALPHA_BBa_B0034}*pLtetO_1'; 'C_3OC12HSL_lasR = lasR*{K_C_3OC12HSL}^{n_C_3OC12HSL}/({K_C_3OC12HSL}^{n_C_3OC12HSL} + C_3OC12HSL^{n_C_3OC12HSL})'; 'pLUX = 200*({beta_pLUX} + ({alpha_pLUX} - {beta_pLUX})*C_3OC12HSL_lasR^{n_pLUX}/({K_pLUX}^{n_pLUX} + C_3OC12HSL_lasR^{n_pLUX}))'; 'GFPmut3b = {scale_GFPmut3b_doi_10_1016_j_ces_2012_12_016}*{ALPHA_BBa_B0034}*pLUX'});
	%%%%%%%%%%%%%
	IPTG = [0.005 0.01 0.02 0.04 0.1];
	gfp_LVA = [0.004 0.004 0.004 0.2 1];
	output_sim = simf_doi_10_1038_msb_2009_75_figure3_1(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; gfp_LVA];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1038/msb.2009.75', 'figure3-1', IPTG, output_sim, {'placI = 1*{alpha_placI}'; 'pLAC = 200*{alpha_pLAC}'; 'gfp_LVA = {scale_gfp_LVA_doi_10_1038_msb_2009_75}*{ALPHA_RBS_UNKNOWN}*pLAC'; 'lacI = {ALPHA_RBS_placI}*placI*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'pLAC = 1*({beta_pLAC} + ({alpha_pLAC} - {beta_pLAC})*{K_pLAC}^{n_pLAC}/({K_pLAC}^{n_pLAC} + lacI^{n_pLAC}))'; 'lacZ = {ALPHA_RBS_pLAC}*pLAC'});
	%%%%%%%%%%%%%
	AHL = [0 0 0 0 0 0 0.002];
	GFPmut3b = [0.018 0.018 0.018 0.212 0.939 1 0.97];
	output_sim = simf_doi_10_1016_j_bios_2012_08_011_figure2d(sim_x_opt,AHL);
	sim_pair_matrix = [output_sim; GFPmut3b];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1016/j.bios.2012.08.011', 'figure2d', AHL, output_sim, {'pTET = 10*{alpha_pTET}'; 'luxR = {ALPHA_BBa_B0034}*pTET'; 'AHL_luxR = luxR*{K_AHL}^{n_AHL}/({K_AHL}^{n_AHL} + AHL^{n_AHL})'; 'pLUX = 10*({beta_pLUX} + ({alpha_pLUX} - {beta_pLUX})*AHL_luxR^{n_pLUX}/({K_pLUX}^{n_pLUX} + AHL_luxR^{n_pLUX}))'; 'GFPmut3b = {scale_GFPmut3b_doi_10_1016_j_bios_2012_08_011}*{ALPHA_BBa_B0032}*pLUX'});
	%%%%%%%%%%%%%
	AHL = [0 0 0 0 0.001 0.01 0.1];
	mRFP1 = [0.024 0.081 0.372 0.763 0.874 0.889 0.891];
	output_sim = simf_doi_10_1093_nar_gku964_figure3A_1(sim_x_opt,AHL);
	sim_pair_matrix = [output_sim; mRFP1];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1093/nar/gku964', 'figure3A-1', AHL, output_sim, {'pLtetO_1 = 200*{alpha_pLtetO_1}'; 'luxR = {ALPHA_BBa_B0034}*pLtetO_1'; 'GFPmut3b = {ALPHA_BBa_B0034}*pLtetO_1'; 'AHL_luxR = luxR*{K_AHL}^{n_AHL}/({K_AHL}^{n_AHL} + AHL^{n_AHL})'; 'pLUX = 200*({beta_pLUX} + ({alpha_pLUX} - {beta_pLUX})*AHL_luxR^{n_pLUX}/({K_pLUX}^{n_pLUX} + AHL_luxR^{n_pLUX}))'; 'mRFP1 = {scale_mRFP1_doi_10_1093_nar_gku964}*{ALPHA_BBa_B0032}*pLUX'});
	%%%%%%%%%%%%%
	AHL = [0 0 0 0 0.001 0.01 0.1];
	mRFP1 = [0.015 0.052 0.295 0.755 0.911 0.931 0.934];
	output_sim = simf_doi_10_1093_nar_gku964_figure3A_2(sim_x_opt,AHL);
	sim_pair_matrix = [output_sim; mRFP1];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1093/nar/gku964', 'figure3A-2', AHL, output_sim, {'pLtetO_1 = 200*{alpha_pLtetO_1}'; 'luxR = {ALPHA_BBa_B0032}*pLtetO_1'; 'GFPmut3b = {ALPHA_BBa_B0032}*pLtetO_1'; 'AHL_luxR = luxR*{K_AHL}^{n_AHL}/({K_AHL}^{n_AHL} + AHL^{n_AHL})'; 'pLUX = 200*({beta_pLUX} + ({alpha_pLUX} - {beta_pLUX})*AHL_luxR^{n_pLUX}/({K_pLUX}^{n_pLUX} + AHL_luxR^{n_pLUX}))'; 'mRFP1 = {scale_mRFP1_doi_10_1093_nar_gku964}*{ALPHA_BBa_B0032}*pLUX'});
	%%%%%%%%%%%%%
	AHL = [0 0 0 0 0.001 0.01 0.1];
	mRFP1 = [0.006 0.01 0.05 0.319 0.82 0.98 1];
	output_sim = simf_doi_10_1093_nar_gku964_figure3A_3(sim_x_opt,AHL);
	sim_pair_matrix = [output_sim; mRFP1];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1093/nar/gku964', 'figure3A-3', AHL, output_sim, {'pLtetO_1 = 200*{alpha_pLtetO_1}'; 'luxR = {ALPHA_BBa_B0033}*pLtetO_1'; 'GFPmut3b = {ALPHA_BBa_B0033}*pLtetO_1'; 'AHL_luxR = luxR*{K_AHL}^{n_AHL}/({K_AHL}^{n_AHL} + AHL^{n_AHL})'; 'pLUX = 200*({beta_pLUX} + ({alpha_pLUX} - {beta_pLUX})*AHL_luxR^{n_pLUX}/({K_pLUX}^{n_pLUX} + AHL_luxR^{n_pLUX}))'; 'mRFP1 = {scale_mRFP1_doi_10_1093_nar_gku964}*{ALPHA_BBa_B0032}*pLUX'});
	%%%%%%%%%%%%%
	IPTG = [0 0.01 0.1 1 10];
	gfp_uv = [0.14 0.14 0.76 0.96 1];
	output_sim = simf_doi_10_1021_bp010141k_figure3_1(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; gfp_uv];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1021/bp010141k', 'figure3-1', IPTG, output_sim, {'placIq = 20*{alpha_placIq}'; 'lacI = {ALPHA_RBS_placI}*placIq*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'pTrc = 20*({beta_pTrc} + ({alpha_pTrc} - {beta_pTrc})*{K_pTrc}^{n_pTrc}/({K_pTrc}^{n_pTrc} + lacI^{n_pTrc}))'; 'gfp_uv = {scale_gfp_uv_doi_10_1021_bp010141k}*{ALPHA_RBS_Unknown}*pTrc'});
	%%%%%%%%%%%%%
	output_sim = simf_doi_10_1021_sb300055e_figure4a_1(sim_x_opt);
	point_pair = [output_sim 	1];
	fprintf(sim_pair_fileID, '%f \t %f\n', point_pair);
	print_simulated_data_to_json(json_fileID, 'doi:10.1021/sb300055e', 'figure4a-1', [], output_sim, {'BBa_J23114 = 1*{alpha_BBa_J23114}'; 'Glk = {scale_Glk_doi_10_1021_sb300055e}*{ALPHA_RBS_Glk}*BBa_J23114'});
	%%%%%%%%%%%%%
	output_sim = simf_doi_10_1021_sb300055e_figure4a_2(sim_x_opt);
	point_pair = [output_sim 	0.429];
	fprintf(sim_pair_fileID, '%f \t %f\n', point_pair);
	print_simulated_data_to_json(json_fileID, 'doi:10.1021/sb300055e', 'figure4a-2', [], output_sim, {'BBa_J23117 = 1*{alpha_BBa_J23117}'; 'Glk = {scale_Glk_doi_10_1021_sb300055e}*{ALPHA_RBS_Glk}*BBa_J23117'});
	%%%%%%%%%%%%%
	output_sim = simf_doi_10_1021_sb300055e_figure4a_3(sim_x_opt);
	point_pair = [output_sim 	0.119];
	fprintf(sim_pair_fileID, '%f \t %f\n', point_pair);
	print_simulated_data_to_json(json_fileID, 'doi:10.1021/sb300055e', 'figure4a-3', [], output_sim, {'BBa_J23109 = 1*{alpha_BBa_J23109}'; 'Glk = {scale_Glk_doi_10_1021_sb300055e}*{ALPHA_RBS_Glk}*BBa_J23109'});
	%%%%%%%%%%%%%
	aTc = [0.3 0.7 1 2 3 4 5 10 25 100];
	gfpmut3_1 = [0.039 0.196 0.294 0.51 0.667 0.843 0.902 0.902 0.98 1];
	output_sim = simf_doi_10_1038_msb4100081_figure3b(sim_x_opt,aTc);
	sim_pair_matrix = [output_sim; gfpmut3_1];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1038/msb4100081', 'figure3b', aTc, output_sim, {'plac_ara_1 = 5*{alpha_plac_ara_1}'; 'tetR = {ALPHA_RBSII}*plac_ara_1*{K_aTc}^{n_aTc}/({K_aTc}^{n_aTc} + aTc^{n_aTc})'; 'pLtetO_1 = 15*({beta_pLtetO_1} + ({alpha_pLtetO_1} - {beta_pLtetO_1})*{K_pLtetO_1}^{n_pLtetO_1}/({K_pLtetO_1}^{n_pLtetO_1} + tetR^{n_pLtetO_1}))'; 'gfpmut3_1 = {scale_gfpmut3_1_doi_10_1038_msb4100081}*{ALPHA_RBSII}*pLtetO_1'});
	%%%%%%%%%%%%%
	IPTG = [0 0 0 0 0.001 0.003 0.01 0.04 0.2 1];
	gfp = [0.211 0.211 0.237 0.237 0.263 0.447 1 0.921 0.868 0.842];
	output_sim = simf_doi_10_1186_1754_1611_6_9_figure2_b(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; gfp];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1186/1754-1611-6-9', 'figure2-b', IPTG, output_sim, {'placI = 15*{alpha_placI}'; 'lacI = {ALPHA_RBS_placI}*placI*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'pLlacO_1 = 15*({beta_pLlacO_1} + ({alpha_pLlacO_1} - {beta_pLlacO_1})*{K_pLlacO_1}^{n_pLlacO_1}/({K_pLlacO_1}^{n_pLlacO_1} + lacI^{n_pLlacO_1}))'; 'cI = {ALPHA_RBS_UTR}*pLlacO_1'; 'pLAMBDA = 10*({beta_pLAMBDA} + ({alpha_pLAMBDA} - {beta_pLAMBDA})*{K_pLAMBDA}^{n_pLAMBDA}/({K_pLAMBDA}^{n_pLAMBDA} + cI^{n_pLAMBDA}))'; 'gfp = {scale_gfp_doi_10_1186_1754_1611_6_9}*{ALPHA_RBS_pLAMBDA}*pLAMBDA'});
	%%%%%%%%%%%%%
	L_arabinose = [0.03 0.08 0.2 0.5 1.2 3 8 20 50 125 300];
	gfpmut2 = [0.014 0.083 0.181 0.361 0.5 0.694 0.778 0.889 1 0.917 0.778];
	output_sim = simf_doi_10_1186_1752_0509_5_111_figure3_a(sim_x_opt,L_arabinose);
	sim_pair_matrix = [output_sim; gfpmut2];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1186/1752-0509-5-111', 'figure3-a', L_arabinose, output_sim, {'pC = 1*{alpha_pC}'; 'araC = {ALPHA_RBS_pC}*pC*{K_L_arabinose}^{n_L_arabinose}/({K_L_arabinose}^{n_L_arabinose} + L_arabinose^{n_L_arabinose})'; 'pBAD = 5*({beta_pBAD} + ({alpha_pBAD} - {beta_pBAD})*{K_pBAD}^{n_pBAD}/({K_pBAD}^{n_pBAD} + araC^{n_pBAD}))'; 'gfpmut2 = {scale_gfpmut2_doi_10_1186_1752_0509_5_111}*{ALPHA_RBS_pBAD}*pBAD'});
	%%%%%%%%%%%%%
	L_arabinose = [0.03 0.08 0.2 0.5 1.2 3 8 20 50 125 300];
	gfpmut2 = [0.014 0.014 0.014 0.017 0.021 0.028 0.061 0.167 0.389 0.597 0.556];
	output_sim = simf_doi_10_1186_1752_0509_5_111_figure3_b(sim_x_opt,L_arabinose);
	sim_pair_matrix = [output_sim; gfpmut2];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1186/1752-0509-5-111', 'figure3-b', L_arabinose, output_sim, {'pN25 = 1*{alpha_pN25}'; 'tetR = {ALPHA_RBS_pN25}*pN25'; 'pLtetO_1 = 15*({beta_pLtetO_1} + ({alpha_pLtetO_1} - {beta_pLtetO_1})*{K_pLtetO_1}^{n_pLtetO_1}/({K_pLtetO_1}^{n_pLtetO_1} + tetR^{n_pLtetO_1}))'; 'araC = {ALPHA_RBSII}*pLtetO_1*{K_L_arabinose}^{n_L_arabinose}/({K_L_arabinose}^{n_L_arabinose} + L_arabinose^{n_L_arabinose})'; 'pBAD = 5*({beta_pBAD} + ({alpha_pBAD} - {beta_pBAD})*{K_pBAD}^{n_pBAD}/({K_pBAD}^{n_pBAD} + araC^{n_pBAD}))'; 'gfpmut2 = {scale_gfpmut2_doi_10_1186_1752_0509_5_111}*{ALPHA_RBS_pBAD}*pBAD'});
	%%%%%%%%%%%%%
	AHL = [0 0 0 0 0 0.001 0.004 0.01 0.04 0.1 0.4 1];
	eYFP = [0.001 0.029 0.114 0.286 0.429 0.643 0.743 0.829 0.857 0.9 0.929 1];
	output_sim = simf_dx_doi_org_10_1021_sb400152n_figure2(sim_x_opt,AHL);
	sim_pair_matrix = [output_sim; eYFP];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'dx.doi.org/10.1021/sb400152n', 'figure2', AHL, output_sim, {'pCAT = 10*{alpha_pCAT}'; 'luxR = {ALPHA_BBa_B0034}*pCAT'; 'AHL_luxR = luxR*{K_AHL}^{n_AHL}/({K_AHL}^{n_AHL} + AHL^{n_AHL})'; 'pLUX = 10*({beta_pLUX} + ({alpha_pLUX} - {beta_pLUX})*AHL_luxR^{n_pLUX}/({K_pLUX}^{n_pLUX} + AHL_luxR^{n_pLUX}))'; 'eYFP = {scale_eYFP_dx_doi_org_10_1021_sb400152n}*{ALPHA_BBa_B0034}*pLUX'});
	%%%%%%%%%%%%%
	IPTG = [0 0.001 0.002 0.01 0.02 0.1 1];
	eYFP = [1 0.867 0.733 0.667 0.333 0.2 0.133];
	output_sim = simf_doi_10_1073_pnas_252535999_figure3(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; eYFP];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1073/pnas.252535999', 'figure3', IPTG, output_sim, {'placIq = 10*{alpha_placIq}'; 'lacI = {ALPHA_RBS_placI}*placIq*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'pLAC = 10*({beta_pLAC} + ({alpha_pLAC} - {beta_pLAC})*{K_pLAC}^{n_pLAC}/({K_pLAC}^{n_pLAC} + lacI^{n_pLAC}))'; 'cI = {ALPHA_RBSII}*pLAC'; 'eCFP = {ALPHA_RBSII}*pLAC'; 'pLAMBDA = 15*({beta_pLAMBDA} + ({alpha_pLAMBDA} - {beta_pLAMBDA})*{K_pLAMBDA}^{n_pLAMBDA}/({K_pLAMBDA}^{n_pLAMBDA} + cI^{n_pLAMBDA}))'; 'eYFP = {scale_eYFP_doi_10_1073_pnas_252535999}*{ALPHA_RBSII}*pLAMBDA'});
	%%%%%%%%%%%%%
	IPTG = [0.001 0.01 0.02 0.03 0.04 0.06 0.1 0.3 0.6 1 3 6 10];
	GFPmut3b = [0.011 0.011 0.011 0.011 0.022 0.056 0.111 0.444 0.722 0.889 1 1 1];
	output_sim = simf_doi_10_1038_35002131_figure5a(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; GFPmut3b];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1038/35002131', 'figure5a', IPTG, output_sim, {'pLs1con = 15*{alpha_pLs1con}'; 'lacI = {ALPHA_RBS_D}*pLs1con*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'pTrc_2 = 15*({beta_pTrc_2} + ({alpha_pTrc_2} - {beta_pTrc_2})*{K_pTrc_2}^{n_pTrc_2}/({K_pTrc_2}^{n_pTrc_2} + lacI^{n_pTrc_2}))'; 'GFPmut3b = {scale_GFPmut3b_doi_10_1038_35002131}*{ALPHA_RBS_E}*pTrc_2'});
	%%%%%%%%%%%%%
	aTc = [0.1 1 10 100 1000];
	glk = [1 0.786 0.357 0.071 0.036];
	output_sim = simf_doi_10_1016_j_ymben_2012_08_006_figure7(sim_x_opt,aTc);
	sim_pair_matrix = [output_sim; glk];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1016/j.ymben.2012.08.006', 'figure7', aTc, output_sim, {'pC_tetR_UNKNOWN = 1*{alpha_pC_tetR_UNKNOWN}'; 'tetR = {ALPHA_RBS_tetR}*pC_tetR_UNKNOWN*{K_aTc}^{n_aTc}/({K_aTc}^{n_aTc} + aTc^{n_aTc})'; 'pTET = 1*({beta_pTET} + ({alpha_pTET} - {beta_pTET})*{K_pTET}^{n_pTET}/({K_pTET}^{n_pTET} + tetR^{n_pTET}))'; 'lacI = {ALPHA_BBa_B0034}*pTET'; 'pJ23114_lacO = 1*({beta_pJ23114_lacO} + ({alpha_pJ23114_lacO} - {beta_pJ23114_lacO})*{K_pJ23114_lacO}^{n_pJ23114_lacO}/({K_pJ23114_lacO}^{n_pJ23114_lacO} + lacI^{n_pJ23114_lacO}))'; 'glk = {scale_glk_doi_10_1016_j_ymben_2012_08_006}*{ALPHA_RBS_native}*pJ23114_lacO'});
	%%%%%%%%%%%%%
	IPTG = [0.001 0.01 0.02 0.03 0.05 0.1 0.2 0.5 1];
	lacZ = [0.017 0.042 0.1 0.15 0.3 0.6 0.85 0.975 1];
	output_sim = simf_doi_10_1006_plas_2000_1477_figure5(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; lacZ];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1006/plas.2000.1477', 'figure5', IPTG, output_sim, {'placI = 1*{alpha_placI}'; 'lacI = {ALPHA_RBS_placI}*placI*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'pLAC = 15*({beta_pLAC} + ({alpha_pLAC} - {beta_pLAC})*{K_pLAC}^{n_pLAC}/({K_pLAC}^{n_pLAC} + lacI^{n_pLAC}))'; 'lacZ = {scale_lacZ_doi_10_1006_plas_2000_1477}*{ALPHA_RBS_pLAC}*pLAC'});
	%%%%%%%%%%%%%
	AHL = [0 0 0 0.001 0.005 0.01 0.1 1 10 100];
	GFPmut3b = [0 0.013 0.108 0.343 0.873 0.954 0.999 1 1 1];
	output_sim = simf_doi_10_1038_nbt1413_figure3(sim_x_opt,AHL);
	sim_pair_matrix = [output_sim; GFPmut3b];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1038/nbt1413', 'figure3', AHL, output_sim, {'pLtetO_1 = 10*{alpha_pLtetO_1}'; 'luxR = {ALPHA_BBa_B0034}*pLtetO_1'; 'AHL_luxR = luxR*{K_AHL}^{n_AHL}/({K_AHL}^{n_AHL} + AHL^{n_AHL})'; 'pLUX = 10*({beta_pLUX} + ({alpha_pLUX} - {beta_pLUX})*AHL_luxR^{n_pLUX}/({K_pLUX}^{n_pLUX} + AHL_luxR^{n_pLUX}))'; 'GFPmut3b = {scale_GFPmut3b_doi_10_1038_nbt1413}*{ALPHA_BBa_B0032}*pLUX'});
	%%%%%%%%%%%%%
	C_3OC12HSL = [0 0.001 0.01 0.1 1 10];
	GFPmut3b = [0 0.225 0.693 0.874 0.898 0.9];
	output_sim = simf_doi_10_1093_nar_gkt052_figure3b(sim_x_opt,C_3OC12HSL);
	sim_pair_matrix = [output_sim; GFPmut3b];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1093/nar/gkt052', 'figure3b', C_3OC12HSL, output_sim, {'BBa_J23100 = 200*{alpha_BBa_J23100}'; 'lasR = {ALPHA_BBa_B0034}*BBa_J23100'; 'C_3OC12HSL_lasR = lasR*{K_C_3OC12HSL}^{n_C_3OC12HSL}/({K_C_3OC12HSL}^{n_C_3OC12HSL} + C_3OC12HSL^{n_C_3OC12HSL})'; 'pLasRV = 200*({beta_pLasRV} + ({alpha_pLasRV} - {beta_pLasRV})*C_3OC12HSL_lasR^{n_pLasRV}/({K_pLasRV}^{n_pLasRV} + C_3OC12HSL_lasR^{n_pLasRV}))'; 'GFPmut3b = {scale_GFPmut3b_doi_10_1093_nar_gkt052}*{ALPHA_BBa_B0034}*pLasRV'});
	%%%%%%%%%%%%%
	C_3OC12HSL = [0 0.001 0.01 0.1 1 10];
	GFPmut3b = [0 0.048 0.446 0.927 0.995 1];
	output_sim = simf_doi_10_1093_nar_gkt052_figure3c(sim_x_opt,C_3OC12HSL);
	sim_pair_matrix = [output_sim; GFPmut3b];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1093/nar/gkt052', 'figure3c', C_3OC12HSL, output_sim, {'BBa_J23100 = 200*{alpha_BBa_J23100}'; 'lasR = {ALPHA_BBa_B0034}*BBa_J23100'; 'C_3OC12HSL_lasR = lasR*{K_C_3OC12HSL}^{n_C_3OC12HSL}/({K_C_3OC12HSL}^{n_C_3OC12HSL} + C_3OC12HSL^{n_C_3OC12HSL})'; 'pLasR3 = 200*({beta_pLasR3} + ({alpha_pLasR3} - {beta_pLasR3})*C_3OC12HSL_lasR^{n_pLasR3}/({K_pLasR3}^{n_pLasR3} + C_3OC12HSL_lasR^{n_pLasR3}))'; 'GFPmut3b = {scale_GFPmut3b_doi_10_1093_nar_gkt052}*{ALPHA_BBa_B0034}*pLasR3'});
	%%%%%%%%%%%%%
	C_3OC12HSL = [0 0.001 0.01 0.1 1 10];
	GFPmut3b = [0 0.001 0.15 0.299 0.3 0.3];
	output_sim = simf_doi_10_1093_nar_gkt052_figure3d(sim_x_opt,C_3OC12HSL);
	sim_pair_matrix = [output_sim; GFPmut3b];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1093/nar/gkt052', 'figure3d', C_3OC12HSL, output_sim, {'BBa_J23100 = 200*{alpha_BBa_J23100}'; 'lasR = {ALPHA_BBa_B0034}*BBa_J23100'; 'C_3OC12HSL_lasR = lasR*{K_C_3OC12HSL}^{n_C_3OC12HSL}/({K_C_3OC12HSL}^{n_C_3OC12HSL} + C_3OC12HSL^{n_C_3OC12HSL})'; 'pLasRL = 200*({beta_pLasRL} + ({alpha_pLasRL} - {beta_pLasRL})*C_3OC12HSL_lasR^{n_pLasRL}/({K_pLasRL}^{n_pLasRL} + C_3OC12HSL_lasR^{n_pLasRL}))'; 'GFPmut3b = {scale_GFPmut3b_doi_10_1093_nar_gkt052}*{ALPHA_BBa_B0034}*pLasRL'});
	%%%%%%%%%%%%%
	Mg2_plus__plus_ = [0 0.001 0.005 0.05 0.5 5 50];
	gfp = [0.054 0.05 0.047 0.043 0.028 0.017 0.01];
	output_sim = simf_doi_10_1038_nmeth_2205_figureS1_1(sim_x_opt,Mg2_plus__plus_);
	sim_pair_matrix = [output_sim; gfp];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1038/nmeth.2205', 'figureS1-1', Mg2_plus__plus_, output_sim, {'pMgrB = 15*({beta_pMgrB} + ({alpha_pMgrB} - {beta_pMgrB})*{K_pMgrB}^{n_pMgrB}/({K_pMgrB}^{n_pMgrB} + Mg2_plus__plus_^{n_pMgrB}))'; 'gfp = {scale_gfp_doi_10_1038_nmeth_2205}*{ALPHA_RBS_KDL}*pMgrB'});
	%%%%%%%%%%%%%
	L_arabinose = [0.017 0.027 0.11 0.27 0.67 1.3 6.7];
	gfp = [0.007 0.009 0.017 0.06 0.121 0.115 0.103];
	output_sim = simf_doi_10_1038_nmeth_2205_figureS1_2(sim_x_opt,L_arabinose);
	sim_pair_matrix = [output_sim; gfp];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1038/nmeth.2205', 'figureS1-2', L_arabinose, output_sim, {'pC = 1*{alpha_pC}'; 'araC = {ALPHA_RBS_pC}*pC*{K_L_arabinose}^{n_L_arabinose}/({K_L_arabinose}^{n_L_arabinose} + L_arabinose^{n_L_arabinose})'; 'pBAD = 15*({beta_pBAD} + ({alpha_pBAD} - {beta_pBAD})*{K_pBAD}^{n_pBAD}/({K_pBAD}^{n_pBAD} + araC^{n_pBAD}))'; 'gfp = {scale_gfp_doi_10_1038_nmeth_2205}*{ALPHA_RBS_KDL}*pBAD'});
	%%%%%%%%%%%%%
	aTc = [0.25 1 4 16 64 240];
	gfp = [0.002 0.002 0.003 0.009 0.036 0.043];
	output_sim = simf_doi_10_1038_nmeth_2205_figureS1_3(sim_x_opt,aTc);
	sim_pair_matrix = [output_sim; gfp];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1038/nmeth.2205', 'figureS1-3', aTc, output_sim, {'pN25 = 1*{alpha_pN25}'; 'tetR = {ALPHA_RBS_pN25}*pN25*{K_aTc}^{n_aTc}/({K_aTc}^{n_aTc} + aTc^{n_aTc})'; 'pLtetO_1 = 15*({beta_pLtetO_1} + ({alpha_pLtetO_1} - {beta_pLtetO_1})*{K_pLtetO_1}^{n_pLtetO_1}/({K_pLtetO_1}^{n_pLtetO_1} + tetR^{n_pLtetO_1}))'; 'gfp = {scale_gfp_doi_10_1038_nmeth_2205}*{ALPHA_RBS_KDL}*pLtetO_1'});
	%%%%%%%%%%%%%
	IPTG = [0 0 0.001 0.01 0.1 1];
	gfp = [0.009 0.034 0.069 0.086 0.095 0.121];
	output_sim = simf_doi_10_1038_nmeth_2205_figureS1_41(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; gfp];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1038/nmeth.2205', 'figureS1-41', IPTG, output_sim, {'placIq = 1*{alpha_placIq}'; 'lacI = {ALPHA_RBS_placI}*placIq*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'pLlacO_1 = 15*({beta_pLlacO_1} + ({alpha_pLlacO_1} - {beta_pLlacO_1})*{K_pLlacO_1}^{n_pLlacO_1}/({K_pLlacO_1}^{n_pLlacO_1} + lacI^{n_pLlacO_1}))'; 'gfp = {scale_gfp_doi_10_1038_nmeth_2205}*{ALPHA_RBS_KDL}*pLlacO_1'});
	%%%%%%%%%%%%%
	IPTG = [0 0 0.001 0.01 0.1 1];
	gfp = [0.129 0.129 0.129 0.345 0.905 1];
	output_sim = simf_doi_10_1038_nmeth_2205_figureS1_42(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; gfp];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1038/nmeth.2205', 'figureS1-42', IPTG, output_sim, {'placIq = 1*{alpha_placIq}'; 'lacI = {ALPHA_RBS_placI}*placIq*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'pTrc_2 = 15*({beta_pTrc_2} + ({alpha_pTrc_2} - {beta_pTrc_2})*{K_pTrc_2}^{n_pTrc_2}/({K_pTrc_2}^{n_pTrc_2} + lacI^{n_pTrc_2}))'; 'gfp = {scale_gfp_doi_10_1038_nmeth_2205}*{ALPHA_RBS_KDL}*pTrc_2'});
	%%%%%%%%%%%%%
	IPTG = [0 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1 0.15 0.2 1];
	eYFP = [0.001 0.002 0.003 0.005 0.008 0.016 0.024 0.032 0.06 0.1 0.14 0.2];
	output_sim = simf_doi_10_1038_nchembio_1411_figure4_1(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; eYFP];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1038/nchembio.1411', 'figure4-1', IPTG, output_sim, {'placI = 10*{alpha_placI}'; 'luxR = {ALPHA_rbs1}*placI'; 'lacI = {ALPHA_RBS_placI}*placI*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'pTAC = 10*({beta_pTAC} + ({alpha_pTAC} - {beta_pTAC})*{K_pTAC}^{n_pTAC}/({K_pTAC}^{n_pTAC} + lacI^{n_pTAC}))'; 'eYFP = {scale_eYFP_doi_10_1038_nchembio_1411}*{ALPHA_rbs1}*pTAC'});
	%%%%%%%%%%%%%
	IPTG = [0 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1 0.15 0.2 1];
	eYFP = [1 1 0.8 0.2 0.04 0.028 0.02 0.012 0.01 0.008 0.008 0.008];
	output_sim = simf_doi_10_1038_nchembio_1411_figure4_2(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; eYFP];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1038/nchembio.1411', 'figure4-2', IPTG, output_sim, {'J23119_tetO = 10*{alpha_J23119_tetO}'; 'placI = 10*{alpha_placI}'; 'luxR = {ALPHA_rbs1}*placI'; 'eYFP = {scale_eYFP_doi_10_1038_nchembio_1411}*{ALPHA_BBa_B0034}*J23119_tetO'; 'lacI = {ALPHA_RBS_placI}*placI*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'pTAC = 10*({beta_pTAC} + ({alpha_pTAC} - {beta_pTAC})*{K_pTAC}^{n_pTAC}/({K_pTAC}^{n_pTAC} + lacI^{n_pTAC}))'; 'tetR = {ALPHA_rbs1}*pTAC'});
	%%%%%%%%%%%%%
	aTc = [0 15];
	mRFP1 = [0.08 1];
	output_sim = simf_doi_10_1093_nar_gks583_figure3(sim_x_opt,aTc);
	sim_pair_matrix = [output_sim; mRFP1];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1093/nar/gks583', 'figure3', aTc, output_sim, {'pN25 = 1*{alpha_pN25}'; 'tetR = {ALPHA_RBS_pN25}*pN25*{K_aTc}^{n_aTc}/({K_aTc}^{n_aTc} + aTc^{n_aTc})'; 'PtetA = 1*({beta_PtetA} + ({alpha_PtetA} - {beta_PtetA})*{K_PtetA}^{n_PtetA}/({K_PtetA}^{n_PtetA} + tetR^{n_PtetA}))'; 'mRFP1 = {scale_mRFP1_doi_10_1093_nar_gks583}*{ALPHA_RBS_PtetA}*PtetA'});
	%%%%%%%%%%%%%
	C4HSL = [0 0 0 0 0.001 0.004 0.01 0.04];
	eYFP = [1 1 1 0.667 0.433 0.167 0.033 0.013];
	output_sim = simf_doi_10_1002_bit_20371_figure3(sim_x_opt,C4HSL);
	sim_pair_matrix = [output_sim; eYFP];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1002/bit.20371', 'figure3', C4HSL, output_sim, {'placIq = 10*{alpha_placIq}'; 'rhlR = {ALPHA_RBS_placI}*placIq'; 'C4HSL_rhlR = rhlR*{K_C4HSL}^{n_C4HSL}/({K_C4HSL}^{n_C4HSL} + C4HSL^{n_C4HSL})'; 'qsc = 10*({beta_qsc} + ({alpha_qsc} - {beta_qsc})*C4HSL_rhlR^{n_qsc}/({K_qsc}^{n_qsc} + C4HSL_rhlR^{n_qsc}))'; 'cI_LVA = {ALPHA_RBSII}*qsc'; 'pLAMBDA_R_O12 = 15*({beta_pLAMBDA_R_O12} + ({alpha_pLAMBDA_R_O12} - {beta_pLAMBDA_R_O12})*{K_pLAMBDA_R_O12}^{n_pLAMBDA_R_O12}/({K_pLAMBDA_R_O12}^{n_pLAMBDA_R_O12} + cI_LVA^{n_pLAMBDA_R_O12}))'; 'eYFP = {scale_eYFP_doi_10_1002_bit_20371}*{ALPHA_RBSII}*pLAMBDA_R_O12'});
	%%%%%%%%%%%%%
	IPTG = [0.001 0.005 0.01 0.015 0.025 0.04 0.05 0.1 0.2 0.4];
	GFPmut3b = [0.01 0.1 0.27 0.35 0.5 0.7 0.9 0.9 1 1];
	output_sim = simf_doi_10_1073_pnas_1314114111_figure5_b(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; GFPmut3b];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1073/pnas.1314114111', 'figure5-b', IPTG, output_sim, {'placIq = 20*{alpha_placIq}'; 'lacI = {ALPHA_RBS_placI}*placIq*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'pTAC = 20*({beta_pTAC} + ({alpha_pTAC} - {beta_pTAC})*{K_pTAC}^{n_pTAC}/({K_pTAC}^{n_pTAC} + lacI^{n_pTAC}))'; 'GFPmut3b = {scale_GFPmut3b_doi_10_1073_pnas_1314114111}*{ALPHA_RBS_pTAC}*pTAC'});
	%%%%%%%%%%%%%
	aTc = [9.2 13.8 18.4 23 27.6 32.2 46 92 138];
	mCherry = [0.111 0.333 0.533 0.667 0.756 0.867 0.944 1 1];
	output_sim = simf_doi_10_1073_pnas_1314114111_figure5_d(sim_x_opt,aTc);
	sim_pair_matrix = [output_sim; mCherry];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1073/pnas.1314114111', 'figure5-d', aTc, output_sim, {'pN25 = 1*{alpha_pN25}'; 'tetR = {ALPHA_RBS_pN25}*pN25*{K_aTc}^{n_aTc}/({K_aTc}^{n_aTc} + aTc^{n_aTc})'; 'pLtetO_1 = 10*({beta_pLtetO_1} + ({alpha_pLtetO_1} - {beta_pLtetO_1})*{K_pLtetO_1}^{n_pLtetO_1}/({K_pLtetO_1}^{n_pLtetO_1} + tetR^{n_pLtetO_1}))'; 'mCherry = {scale_mCherry_doi_10_1073_pnas_1314114111}*{ALPHA_RBSII}*pLtetO_1'});
	%%%%%%%%%%%%%
	AHL = [0 0 0 0 0];
	sfGFP = [0.06 0.06 0.09 0.15 1];
	output_sim = simf_doi_10_1038_msb_2013_58_figureS14(sim_x_opt,AHL);
	sim_pair_matrix = [output_sim; sfGFP];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1038/msb.2013.58', 'figureS14', AHL, output_sim, {'pluxR = 10*{alpha_pluxR}'; 'luxR = {ALPHA_RBS_Unknown}*pluxR'; 'AHL_luxR = luxR*{K_AHL}^{n_AHL}/({K_AHL}^{n_AHL} + AHL^{n_AHL})'; 'pLUX = 10*({beta_pLUX} + ({alpha_pLUX} - {beta_pLUX})*AHL_luxR^{n_pLUX}/({K_pLUX}^{n_pLUX} + AHL_luxR^{n_pLUX}))'; 'sfGFP = {scale_sfGFP_doi_10_1038_msb_2013_58}*{ALPHA_RBS_Unknown}*pLUX'});
	%%%%%%%%%%%%%
	aTc = [0.46 4.6 46 460 4600];
	deGFP = [0.025 0.025 0.075 0.25 0.8];
	output_sim = simf_doi_10_1021_sb400131a_figure4c_1(sim_x_opt,aTc);
	sim_pair_matrix = [output_sim; deGFP];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1021/sb400131a', 'figure4c-1', aTc, output_sim, {'pN25 = 1*{alpha_pN25}'; 'tetR = {ALPHA_RBS_pN25}*pN25*{K_aTc}^{n_aTc}/({K_aTc}^{n_aTc} + aTc^{n_aTc})'; 'pLtetO_1 = 10*({beta_pLtetO_1} + ({alpha_pLtetO_1} - {beta_pLtetO_1})*{K_pLtetO_1}^{n_pLtetO_1}/({K_pLtetO_1}^{n_pLtetO_1} + tetR^{n_pLtetO_1}))'; 'deGFP = {scale_deGFP_doi_10_1021_sb400131a}*{ALPHA_UTR1}*pLtetO_1'});
	%%%%%%%%%%%%%
	IPTG = [0 0 0.001 0.01 0.1 1];
	deGFP = [0.075 0.075 0.075 0.25 0.875 1];
	output_sim = simf_doi_10_1021_sb400131a_figure4c_2(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; deGFP];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1021/sb400131a', 'figure4c-2', IPTG, output_sim, {'placIq = 1*{alpha_placIq}'; 'lacI = {ALPHA_RBS_placI}*placIq*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'pLlacO_1 = 10*({beta_pLlacO_1} + ({alpha_pLlacO_1} - {beta_pLlacO_1})*{K_pLlacO_1}^{n_pLlacO_1}/({K_pLlacO_1}^{n_pLlacO_1} + lacI^{n_pLlacO_1}))'; 'deGFP = {scale_deGFP_doi_10_1021_sb400131a}*{ALPHA_UTR1}*pLlacO_1'});
	%%%%%%%%%%%%%
	AHL = [0 0 0 0 0 0 0 0 0 0 0];
	mRFP1 = [0.006 0.086 0.126 0.343 0.571 0.657 0.857 0.914 0.943 0.971 1];
	output_sim = simf_doi_10_1186_1471_2105_13_S4_S11_figure2_1(sim_x_opt,AHL);
	sim_pair_matrix = [output_sim; mRFP1];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1186/1471-2105-13-S4-S11', 'figure2-1', AHL, output_sim, {'pLtetO_1 = 10*{alpha_pLtetO_1}'; 'luxR = {ALPHA_BBa_B0034}*pLtetO_1'; 'AHL_luxR = luxR*{K_AHL}^{n_AHL}/({K_AHL}^{n_AHL} + AHL^{n_AHL})'; 'pLUX = 10*({beta_pLUX} + ({alpha_pLUX} - {beta_pLUX})*AHL_luxR^{n_pLUX}/({K_pLUX}^{n_pLUX} + AHL_luxR^{n_pLUX}))'; 'mRFP1 = {scale_mRFP1_doi_10_1186_1471_2105_13_S4_S11}*{ALPHA_BBa_B0034}*pLUX'});
	%%%%%%%%%%%%%
	AHL = [0 0 0 0 0 0 0 0 0 0 0];
	mRFP1 = [0.001 0.011 0.017 0.034 0.057 0.069 0.086 0.097 0.103 0.109 0.114];
	output_sim = simf_doi_10_1186_1471_2105_13_S4_S11_figure2_2(sim_x_opt,AHL);
	sim_pair_matrix = [output_sim; mRFP1];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1186/1471-2105-13-S4-S11', 'figure2-2', AHL, output_sim, {'pLtetO_1 = 5*{alpha_pLtetO_1}'; 'luxR = {ALPHA_BBa_B0034}*pLtetO_1'; 'AHL_luxR = luxR*{K_AHL}^{n_AHL}/({K_AHL}^{n_AHL} + AHL^{n_AHL})'; 'pLUX = 5*({beta_pLUX} + ({alpha_pLUX} - {beta_pLUX})*AHL_luxR^{n_pLUX}/({K_pLUX}^{n_pLUX} + AHL_luxR^{n_pLUX}))'; 'mRFP1 = {scale_mRFP1_doi_10_1186_1471_2105_13_S4_S11}*{ALPHA_BBa_B0034}*pLUX'});
	%%%%%%%%%%%%%
	AHL = [0 0 0 0 0 0 0 0 0 0 0];
	mRFP1 = [0.001 0.001 0.001 0.002 0.007 0.014 0.026 0.03 0.033 0.034 0.035];
	output_sim = simf_doi_10_1186_1471_2105_13_S4_S11_figure2_3(sim_x_opt,AHL);
	sim_pair_matrix = [output_sim; mRFP1];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1186/1471-2105-13-S4-S11', 'figure2-3', AHL, output_sim, {'pLtetO_1 = 1*{alpha_pLtetO_1}'; 'luxR = {ALPHA_BBa_B0034}*pLtetO_1'; 'AHL_luxR = luxR*{K_AHL}^{n_AHL}/({K_AHL}^{n_AHL} + AHL^{n_AHL})'; 'pLUX = 1*({beta_pLUX} + ({alpha_pLUX} - {beta_pLUX})*AHL_luxR^{n_pLUX}/({K_pLUX}^{n_pLUX} + AHL_luxR^{n_pLUX}))'; 'mRFP1 = {scale_mRFP1_doi_10_1186_1471_2105_13_S4_S11}*{ALPHA_BBa_B0034}*pLUX'});
	%%%%%%%%%%%%%
	AHL = [0.1 0.3 0.5 1 2 4 10 100 1000];
	gfp_LVA = [1 0.947 0.716 0.579 0.158 0.158 0.105 0.053 0.011];
	output_sim = simf_doi_10_1038_nature03461_figure2b_1(sim_x_opt,AHL);
	sim_pair_matrix = [output_sim; gfp_LVA];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1038/nature03461', 'figure2b-1', AHL, output_sim, {'pLUX_L = 15*{alpha_pLUX_L}'; 'luxR_G2F = {ALPHA_pLUX_RBS}*pLUX_L'; 'AHL_luxR = luxR_G2F*{K_AHL}^{n_AHL}/({K_AHL}^{n_AHL} + AHL^{n_AHL})'; 'pLUX = 15*({beta_pLUX} + ({alpha_pLUX} - {beta_pLUX})*AHL_luxR^{n_pLUX}/({K_pLUX}^{n_pLUX} + AHL_luxR^{n_pLUX}))'; 'lacI_M1 = {ALPHA_RBSII}*pLUX'; 'pLAC = 15*({beta_pLAC} + ({alpha_pLAC} - {beta_pLAC})*{K_pLAC}^{n_pLAC}/({K_pLAC}^{n_pLAC} + lacI_M1^{n_pLAC}))'; 'gfp_LVA = {scale_gfp_LVA_doi_10_1038_nature03461}*{ALPHA_RBSII}*pLAC'});
	%%%%%%%%%%%%%
	AHL = [0.1 0.3 0.7 1.2 2 5 10 20 30 50 100 500 1000];
	gfp_LVA = [0.947 0.916 0.937 0.916 0.905 0.842 0.684 0.526 0.368 0.284 0.242 0.137 0.126];
	output_sim = simf_doi_10_1038_nature03461_figure2b_2(sim_x_opt,AHL);
	sim_pair_matrix = [output_sim; gfp_LVA];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1038/nature03461', 'figure2b-2', AHL, output_sim, {'pLUX_L = 15*{alpha_pLUX_L}'; 'luxR = {ALPHA_pLUX_RBS}*pLUX_L'; 'AHL_luxR = luxR*{K_AHL}^{n_AHL}/({K_AHL}^{n_AHL} + AHL^{n_AHL})'; 'pLUX = 15*({beta_pLUX} + ({alpha_pLUX} - {beta_pLUX})*AHL_luxR^{n_pLUX}/({K_pLUX}^{n_pLUX} + AHL_luxR^{n_pLUX}))'; 'lacI_M1 = {ALPHA_RBSII}*pLUX'; 'pLAC = 15*({beta_pLAC} + ({alpha_pLAC} - {beta_pLAC})*{K_pLAC}^{n_pLAC}/({K_pLAC}^{n_pLAC} + lacI_M1^{n_pLAC}))'; 'gfp_LVA = {scale_gfp_LVA_doi_10_1038_nature03461}*{ALPHA_RBSII}*pLAC'});
	%%%%%%%%%%%%%
	AHL = [0.2 0.5 1 2 5 10 50 100 200 300 400 600 1000 2000 3000];
	gfp_LVA = [0.842 0.832 0.8 0.789 0.737 0.695 0.632 0.611 0.526 0.421 0.316 0.263 0.232 0.179 0.126];
	output_sim = simf_doi_10_1038_nature03461_figure2b_3(sim_x_opt,AHL);
	sim_pair_matrix = [output_sim; gfp_LVA];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1038/nature03461', 'figure2b-3', AHL, output_sim, {'pLUX_L = 15*{alpha_pLUX_L}'; 'luxR = {ALPHA_pLUX_RBS}*pLUX_L'; 'AHL_luxR = luxR*{K_AHL}^{n_AHL}/({K_AHL}^{n_AHL} + AHL^{n_AHL})'; 'pLUX = 15*({beta_pLUX} + ({alpha_pLUX} - {beta_pLUX})*AHL_luxR^{n_pLUX}/({K_pLUX}^{n_pLUX} + AHL_luxR^{n_pLUX}))'; 'lacI_M1 = {ALPHA_RBSII}*pLUX'; 'pLAC = 15*({beta_pLAC} + ({alpha_pLAC} - {beta_pLAC})*{K_pLAC}^{n_pLAC}/({K_pLAC}^{n_pLAC} + lacI_M1^{n_pLAC}))'; 'gfp_LVA = {scale_gfp_LVA_doi_10_1038_nature03461}*{ALPHA_RBSII}*pLAC'});
	%%%%%%%%%%%%%
	C4HSL = [1 10 100 200 500 1000 5000 10000 50000 100000];
	gfp_LVA = [0.001 0.001 0.002 0.003 0.004 0.025 0.15 0.4 0.5 0.6];
	output_sim = simf_doi_10_1073_pnas_0704256104_figureS1b_1(sim_x_opt,C4HSL);
	sim_pair_matrix = [output_sim; gfp_LVA];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1073/pnas.0704256104', 'figureS1b-1', C4HSL, output_sim, {'placIq = 10*{alpha_placIq}'; 'rhlR = {ALPHA_RBSII}*placIq'; 'C4HSL_rhlR = rhlR*{K_C4HSL}^{n_C4HSL}/({K_C4HSL}^{n_C4HSL} + C4HSL^{n_C4HSL})'; 'pRHL = 10*({beta_pRHL} + ({alpha_pRHL} - {beta_pRHL})*C4HSL_rhlR^{n_pRHL}/({K_pRHL}^{n_pRHL} + C4HSL_rhlR^{n_pRHL}))'; 'gfp_LVA = {scale_gfp_LVA_doi_10_1073_pnas_0704256104}*{ALPHA_RBSII}*pRHL'});
	%%%%%%%%%%%%%
	C_3OC12HSL = [0.1 1 5 8 10 20 50 80 100 1000];
	gfp_LVA = [0.001 0.001 0.001 0.004 0.025 0.045 0.25 0.45 1 1];
	output_sim = simf_doi_10_1073_pnas_0704256104_figureS1b_2(sim_x_opt,C_3OC12HSL);
	sim_pair_matrix = [output_sim; gfp_LVA];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1073/pnas.0704256104', 'figureS1b-2', C_3OC12HSL, output_sim, {'placIq = 10*{alpha_placIq}'; 'lasR = {ALPHA_RBSII}*placIq'; 'C_3OC12HSL_lasR = lasR*{K_C_3OC12HSL}^{n_C_3OC12HSL}/({K_C_3OC12HSL}^{n_C_3OC12HSL} + C_3OC12HSL^{n_C_3OC12HSL})'; 'pLAS = 10*({beta_pLAS} + ({alpha_pLAS} - {beta_pLAS})*C_3OC12HSL_lasR^{n_pLAS}/({K_pLAS}^{n_pLAS} + C_3OC12HSL_lasR^{n_pLAS}))'; 'gfp_LVA = {scale_gfp_LVA_doi_10_1073_pnas_0704256104}*{ALPHA_RBSII}*pLAS'});
	%%%%%%%%%%%%%
	L_arabinose = [0 0.001 0.015 0.05 0.5 5 10];
	eYFP = [0.002 0.002 0.057 0.628 0.999 1 1];
	output_sim = simf_doi_10_1038_nature09565_figure1C(sim_x_opt,L_arabinose);
	sim_pair_matrix = [output_sim; eYFP];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1038/nature09565', 'figure1C', L_arabinose, output_sim, {'BBa_J23117 = 15*{alpha_BBa_J23117}'; 'araC = {ALPHA_BBa_B0034}*BBa_J23117*{K_L_arabinose}^{n_L_arabinose}/({K_L_arabinose}^{n_L_arabinose} + L_arabinose^{n_L_arabinose})'; 'pBAD = 15*({beta_pBAD} + ({alpha_pBAD} - {beta_pBAD})*{K_pBAD}^{n_pBAD}/({K_pBAD}^{n_pBAD} + araC^{n_pBAD}))'; 'eYFP = {scale_eYFP_doi_10_1038_nature09565}*{ALPHA_BBa_B0033}*pBAD'});
	%%%%%%%%%%%%%
	aTc = [0 0.025 0.25 2.5 25 100 250];
	eYFP = [0.005 0.006 0.007 0.03 0.297 0.388 0.398];
	output_sim = simf_doi_10_1038_nature09565_figureS1C(sim_x_opt,aTc);
	sim_pair_matrix = [output_sim; eYFP];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1038/nature09565', 'figureS1C', aTc, output_sim, {'BBa_J23117 = 15*{alpha_BBa_J23117}'; 'tetR = {ALPHA_BBa_B0032}*BBa_J23117*{K_aTc}^{n_aTc}/({K_aTc}^{n_aTc} + aTc^{n_aTc})'; 'pLtetO_1 = 15*({beta_pLtetO_1} + ({alpha_pLtetO_1} - {beta_pLtetO_1})*{K_pLtetO_1}^{n_pLtetO_1}/({K_pLtetO_1}^{n_pLtetO_1} + tetR^{n_pLtetO_1}))'; 'eYFP = {scale_eYFP_doi_10_1038_nature09565}*{ALPHA_BBa_B0033}*pLtetO_1'});
	%%%%%%%%%%%%%
	IPTG = [0 0.002 0.01 0.05 0.1 0.25 1];
	mRFP1 = [0.054 0.107 0.179 0.571 0.693 0.893 1];
	output_sim = simf_doi_10_1093_nar_gks597_figureS3_left(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; mRFP1];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1093/nar/gks597', 'figureS3-left', IPTG, output_sim, {'placIq = 5*{alpha_placIq}'; 'lacI = {ALPHA_RBS_placI}*placIq*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'pTAC = 5*({beta_pTAC} + ({alpha_pTAC} - {beta_pTAC})*{K_pTAC}^{n_pTAC}/({K_pTAC}^{n_pTAC} + lacI^{n_pTAC}))'; 'mRFP1 = {scale_mRFP1_doi_10_1093_nar_gks597}*{ALPHA_D103}*pTAC'});
	%%%%%%%%%%%%%
	aTc = [0 0.05 0.1 0.25 0.5 1 10 50];
	mRFP1 = [0.032 0.079 0.129 0.193 0.214 0.25 0.279 0.357];
	output_sim = simf_doi_10_1093_nar_gks597_figureS3_right(sim_x_opt,aTc);
	sim_pair_matrix = [output_sim; mRFP1];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1093/nar/gks597', 'figureS3-right', aTc, output_sim, {'pC_tetR_cassette = 5*{alpha_pC_tetR_cassette}'; 'tetR = {ALPHA_RBS_pC_tetR_cassette}*pC_tetR_cassette*{K_aTc}^{n_aTc}/({K_aTc}^{n_aTc} + aTc^{n_aTc})'; 'pTET_tetR_cassette = 5*({beta_pTET_tetR_cassette} + ({alpha_pTET_tetR_cassette} - {beta_pTET_tetR_cassette})*{K_pTET_tetR_cassette}^{n_pTET_tetR_cassette}/({K_pTET_tetR_cassette}^{n_pTET_tetR_cassette} + tetR^{n_pTET_tetR_cassette}))'; 'mRFP1 = {scale_mRFP1_doi_10_1093_nar_gks597}*{ALPHA_D103}*pTET_tetR_cassette'});
	%%%%%%%%%%%%%
	IPTG = [0 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1];
	sfGFP = [0.025 0.025 0.035 0.055 0.115 0.2 0.3 0.4 0.65 1];
	output_sim = simf_doi_10_1038_nbt_2401_figureS11(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; sfGFP];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1038/nbt.2401', 'figureS11', IPTG, output_sim, {'BBa_J23101 = 5*{alpha_BBa_J23101}'; 'lacI = {ALPHA_RBS_A}*BBa_J23101*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'pTAC = 5*({beta_pTAC} + ({alpha_pTAC} - {beta_pTAC})*{K_pTAC}^{n_pTAC}/({K_pTAC}^{n_pTAC} + lacI^{n_pTAC}))'; 'sfGFP = {scale_sfGFP_doi_10_1038_nbt_2401}*{ALPHA_RBS_A}*pTAC'});
	%%%%%%%%%%%%%
	L_arabinose = [0.1 1 2 5 7 10 12.5 25 37.5 50 62.5];
	sfGFP = [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425];
	output_sim = simf_doi_10_1038_nbt_2401_figureS13(sim_x_opt,L_arabinose);
	sim_pair_matrix = [output_sim; sfGFP];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1038/nbt.2401', 'figureS13', L_arabinose, output_sim, {'BBa_J23101 = 5*{alpha_BBa_J23101}'; 'araC = {ALPHA_RBS_A}*BBa_J23101*{K_L_arabinose}^{n_L_arabinose}/({K_L_arabinose}^{n_L_arabinose} + L_arabinose^{n_L_arabinose})'; 'pBAD = 5*({beta_pBAD} + ({alpha_pBAD} - {beta_pBAD})*{K_pBAD}^{n_pBAD}/({K_pBAD}^{n_pBAD} + araC^{n_pBAD}))'; 'sfGFP = {scale_sfGFP_doi_10_1038_nbt_2401}*{ALPHA_RBS_A}*pBAD'});
	%%%%%%%%%%%%%
	L_arabinose = [0 0 0.002 0.01 0.04 0.16 1 4 25];
	mRFP1 = [0.002 0.002 0.002 0.01 0.137 0.487 0.557 0.558 0.558];
	output_sim = simf_doi_10_1038_nature11516_figureS8A(sim_x_opt,L_arabinose);
	sim_pair_matrix = [output_sim; mRFP1];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1038/nature11516', 'figureS8A', L_arabinose, output_sim, {'pC = 15*{alpha_pC}'; 'araC = {ALPHA_RBS_pC}*pC*{K_L_arabinose}^{n_L_arabinose}/({K_L_arabinose}^{n_L_arabinose} + L_arabinose^{n_L_arabinose})'; 'pBAD = 15*({beta_pBAD} + ({alpha_pBAD} - {beta_pBAD})*{K_pBAD}^{n_pBAD}/({K_pBAD}^{n_pBAD} + araC^{n_pBAD}))'; 'mRFP1 = {scale_mRFP1_doi_10_1038_nature11516}*{ALPHA_RBS_psicA}*pBAD'});
	%%%%%%%%%%%%%
	IPTG = [0 0 0 0.001 0.004 0.025 0.1 0.5 2.5];
	mRFP1 = [0.016 0.016 0.016 0.018 0.04 0.273 0.525 0.588 0.593];
	output_sim = simf_doi_10_1038_nature11516_figureS8B(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; mRFP1];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1038/nature11516', 'figureS8B', IPTG, output_sim, {'placIq = 15*{alpha_placIq}'; 'lacI = {ALPHA_RBS_placI}*placIq*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'pTAC = 15*({beta_pTAC} + ({alpha_pTAC} - {beta_pTAC})*{K_pTAC}^{n_pTAC}/({K_pTAC}^{n_pTAC} + lacI^{n_pTAC}))'; 'mRFP1 = {scale_mRFP1_doi_10_1038_nature11516}*{ALPHA_RBS_psicA}*pTAC'});
	%%%%%%%%%%%%%
	AHL = [0 0 0 0 0 0 0 0 0];
	mRFP1 = [0.102 0.102 0.128 0.329 0.89 0.988 0.999 1 1];
	output_sim = simf_doi_10_1038_nature11516_figureS8C(sim_x_opt,AHL);
	sim_pair_matrix = [output_sim; mRFP1];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1038/nature11516', 'figureS8C', AHL, output_sim, {'BBa_J23100 = 15*{alpha_BBa_J23100}'; 'luxR = {ALPHA_RBS2}*BBa_J23100'; 'AHL_luxR = luxR*{K_AHL}^{n_AHL}/({K_AHL}^{n_AHL} + AHL^{n_AHL})'; 'pLUX = 15*({beta_pLUX} + ({alpha_pLUX} - {beta_pLUX})*AHL_luxR^{n_pLUX}/({K_pLUX}^{n_pLUX} + AHL_luxR^{n_pLUX}))'; 'mRFP1 = {scale_mRFP1_doi_10_1038_nature11516}*{ALPHA_RBS_psicA}*pLUX'});
	%%%%%%%%%%%%%
	AHL = [0 0 0 0 0 0 0 0 0];
	mRFP1 = [0.02 0.02 0.026 0.072 0.241 0.278 0.282 0.283 0.283];
	output_sim = simf_doi_10_1038_nature11516_figureS8D(sim_x_opt,AHL);
	sim_pair_matrix = [output_sim; mRFP1];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1038/nature11516', 'figureS8D', AHL, output_sim, {'BBa_J23100 = 15*{alpha_BBa_J23100}'; 'luxR = {ALPHA_RBS2}*BBa_J23100'; 'AHL_luxR = luxR*{K_AHL}^{n_AHL}/({K_AHL}^{n_AHL} + AHL^{n_AHL})'; 'pLUX_star_ = 15*({beta_pLUX_star_} + ({alpha_pLUX_star_} - {beta_pLUX_star_})*AHL_luxR^{n_pLUX_star_}/({K_pLUX_star_}^{n_pLUX_star_} + AHL_luxR^{n_pLUX_star_}))'; 'mRFP1 = {scale_mRFP1_doi_10_1038_nature11516}*{ALPHA_RBS_psicA}*pLUX_star_'});
	%%%%%%%%%%%%%
	IPTG = [0 0.001 0.005 0.01 0.025 0.05 0.1 0.25 0.5];
	mRFP1 = [0.02 0.031 0.092 0.306 0.745 0.898 0.969 1 1];
	output_sim = simf_doi_10_1038_nbt_2149_Supplement_page17(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; mRFP1];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1038/nbt.2149', 'Supplement-page17', IPTG, output_sim, {'placI = 10*{alpha_placI}'; 'lacI = {ALPHA_RBS_placI}*placI*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'placUV5 = 10*({beta_placUV5} + ({alpha_placUV5} - {beta_placUV5})*{K_placUV5}^{n_placUV5}/({K_placUV5}^{n_placUV5} + lacI^{n_placUV5}))'; 'mRFP1 = {scale_mRFP1_doi_10_1038_nbt_2149}*{ALPHA_RBS_pET_29b}*placUV5'});
	%%%%%%%%%%%%%
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	mRFP1 = [0.214 0.679 0.764 0.821 0.821 0.836];
	output_sim = simf_doi_10_1186_1754_1611_5_12_Supplement_page1(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; mRFP1];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1186/1754-1611-5-12', 'Supplement-page1', IPTG, output_sim, {'placIq = 10*{alpha_placIq}'; 'lacI = {ALPHA_RBS_placI}*placIq*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'pTrc = 10*({beta_pTrc} + ({alpha_pTrc} - {beta_pTrc})*{K_pTrc}^{n_pTrc}/({K_pTrc}^{n_pTrc} + lacI^{n_pTrc}))'; 'mRFP1 = {scale_mRFP1_doi_10_1186_1754_1611_5_12}*{ALPHA_RBS_pET_29b}*pTrc'});
	%%%%%%%%%%%%%
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	GFP_uv = [0.051 0.492 0.61 0.678 0.746 0.847];
	output_sim = simf_doi_10_1186_1754_1611_5_12_Supplement_page2(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; GFP_uv];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1186/1754-1611-5-12', 'Supplement-page2', IPTG, output_sim, {'placIq = 20*{alpha_placIq}'; 'lacI = {ALPHA_RBS_placI}*placIq*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'pTrc = 20*({beta_pTrc} + ({alpha_pTrc} - {beta_pTrc})*{K_pTrc}^{n_pTrc}/({K_pTrc}^{n_pTrc} + lacI^{n_pTrc}))'; 'GFP_uv = {scale_GFP_uv_doi_10_1186_1754_1611_5_12}*{ALPHA_RBS_pET_29b}*pTrc'});
	%%%%%%%%%%%%%
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	mRFP1 = [0.296 0.914 0.964 1 1 1];
	output_sim = simf_doi_10_1186_1754_1611_5_12_Supplement_page3(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; mRFP1];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1186/1754-1611-5-12', 'Supplement-page3', IPTG, output_sim, {'placIq = 15*{alpha_placIq}'; 'lacI = {ALPHA_RBS_placI}*placIq*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'pTrc = 15*({beta_pTrc} + ({alpha_pTrc} - {beta_pTrc})*{K_pTrc}^{n_pTrc}/({K_pTrc}^{n_pTrc} + lacI^{n_pTrc}))'; 'mRFP1 = {scale_mRFP1_doi_10_1186_1754_1611_5_12}*{ALPHA_RBS_pET_29b}*pTrc'});
	%%%%%%%%%%%%%
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	mRFP1 = [0.129 0.407 0.443 0.5 0.5 0.5];
	output_sim = simf_doi_10_1186_1754_1611_5_12_Supplement_page4(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; mRFP1];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1186/1754-1611-5-12', 'Supplement-page4', IPTG, output_sim, {'placIq = 5*{alpha_placIq}'; 'lacI = {ALPHA_RBS_placI}*placIq*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'pTrc = 5*({beta_pTrc} + ({alpha_pTrc} - {beta_pTrc})*{K_pTrc}^{n_pTrc}/({K_pTrc}^{n_pTrc} + lacI^{n_pTrc}))'; 'mRFP1 = {scale_mRFP1_doi_10_1186_1754_1611_5_12}*{ALPHA_RBS_pET_29b}*pTrc'});
	%%%%%%%%%%%%%
	aTc = [0 0.013 0.025 0.05 0.1 0.2];
	mRFP1 = [0.007 0.429 0.714 0.929 0.986 0.964];
	output_sim = simf_doi_10_1186_1754_1611_5_12_Supplement_page5(sim_x_opt,aTc);
	sim_pair_matrix = [output_sim; mRFP1];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1186/1754-1611-5-12', 'Supplement-page5', aTc, output_sim, {'pC_tetR = 10*{alpha_pC_tetR}'; 'tetR = {ALPHA_RBS_pC_tetR}*pC_tetR*{K_aTc}^{n_aTc}/({K_aTc}^{n_aTc} + aTc^{n_aTc})'; 'pTET_operon = 10*({beta_pTET_operon} + ({alpha_pTET_operon} - {beta_pTET_operon})*{K_pTET_operon}^{n_pTET_operon}/({K_pTET_operon}^{n_pTET_operon} + tetR^{n_pTET_operon}))'; 'mRFP1 = {scale_mRFP1_doi_10_1186_1754_1611_5_12}*{ALPHA_RBS_pET_29b}*pTET_operon'});
	%%%%%%%%%%%%%
	aTc = [0 0.013 0.025 0.05 0.1 0.2];
	GFP_uv = [0.002 0.288 0.542 0.797 0.932 0.966];
	output_sim = simf_doi_10_1186_1754_1611_5_12_Supplement_page6(sim_x_opt,aTc);
	sim_pair_matrix = [output_sim; GFP_uv];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1186/1754-1611-5-12', 'Supplement-page6', aTc, output_sim, {'pC_tetR = 20*{alpha_pC_tetR}'; 'tetR = {ALPHA_RBS_pC_tetR}*pC_tetR*{K_aTc}^{n_aTc}/({K_aTc}^{n_aTc} + aTc^{n_aTc})'; 'pTET_operon = 20*({beta_pTET_operon} + ({alpha_pTET_operon} - {beta_pTET_operon})*{K_pTET_operon}^{n_pTET_operon}/({K_pTET_operon}^{n_pTET_operon} + tetR^{n_pTET_operon}))'; 'GFP_uv = {scale_GFP_uv_doi_10_1186_1754_1611_5_12}*{ALPHA_RBS_pET_29b}*pTET_operon'});
	%%%%%%%%%%%%%
	aTc = [0 0.013 0.025 0.05 0.1 0.2];
	mRFP1 = [0.007 0.271 0.643 0.807 0.929 0.893];
	output_sim = simf_doi_10_1186_1754_1611_5_12_Supplement_page7(sim_x_opt,aTc);
	sim_pair_matrix = [output_sim; mRFP1];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1186/1754-1611-5-12', 'Supplement-page7', aTc, output_sim, {'pC_tetR = 15*{alpha_pC_tetR}'; 'tetR = {ALPHA_RBS_pC_tetR}*pC_tetR*{K_aTc}^{n_aTc}/({K_aTc}^{n_aTc} + aTc^{n_aTc})'; 'pTET_operon = 15*({beta_pTET_operon} + ({alpha_pTET_operon} - {beta_pTET_operon})*{K_pTET_operon}^{n_pTET_operon}/({K_pTET_operon}^{n_pTET_operon} + tetR^{n_pTET_operon}))'; 'mRFP1 = {scale_mRFP1_doi_10_1186_1754_1611_5_12}*{ALPHA_RBS_pET_29b}*pTET_operon'});
	%%%%%%%%%%%%%
	aTc = [0 0.013 0.025 0.05 0.1 0.2];
	mRFP1 = [0.007 0.429 0.557 0.571 0.557 0.543];
	output_sim = simf_doi_10_1186_1754_1611_5_12_Supplement_page8(sim_x_opt,aTc);
	sim_pair_matrix = [output_sim; mRFP1];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1186/1754-1611-5-12', 'Supplement-page8', aTc, output_sim, {'pC_tetR = 5*{alpha_pC_tetR}'; 'tetR = {ALPHA_RBS_pC_tetR}*pC_tetR*{K_aTc}^{n_aTc}/({K_aTc}^{n_aTc} + aTc^{n_aTc})'; 'pTET_operon = 5*({beta_pTET_operon} + ({alpha_pTET_operon} - {beta_pTET_operon})*{K_pTET_operon}^{n_pTET_operon}/({K_pTET_operon}^{n_pTET_operon} + tetR^{n_pTET_operon}))'; 'mRFP1 = {scale_mRFP1_doi_10_1186_1754_1611_5_12}*{ALPHA_RBS_pET_29b}*pTET_operon'});
	%%%%%%%%%%%%%
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	mRFP1 = [0.021 0.157 0.2 0.286 0.3 0.279];
	output_sim = simf_doi_10_1186_1754_1611_5_12_Supplement_page17(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; mRFP1];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1186/1754-1611-5-12', 'Supplement-page17', IPTG, output_sim, {'placIq = 10*{alpha_placIq}'; 'lacI = {ALPHA_RBS_placI}*placIq*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'placUV5 = 10*({beta_placUV5} + ({alpha_placUV5} - {beta_placUV5})*{K_placUV5}^{n_placUV5}/({K_placUV5}^{n_placUV5} + lacI^{n_placUV5}))'; 'mRFP1 = {scale_mRFP1_doi_10_1186_1754_1611_5_12}*{ALPHA_RBS_pET_29b}*placUV5'});
	%%%%%%%%%%%%%
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	GFP_uv = [0.017 0.09 0.115 0.153 0.161 0.183];
	output_sim = simf_doi_10_1186_1754_1611_5_12_Supplement_page18(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; GFP_uv];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1186/1754-1611-5-12', 'Supplement-page18', IPTG, output_sim, {'placIq = 20*{alpha_placIq}'; 'lacI = {ALPHA_RBS_placI}*placIq*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'placUV5 = 20*({beta_placUV5} + ({alpha_placUV5} - {beta_placUV5})*{K_placUV5}^{n_placUV5}/({K_placUV5}^{n_placUV5} + lacI^{n_placUV5}))'; 'GFP_uv = {scale_GFP_uv_doi_10_1186_1754_1611_5_12}*{ALPHA_RBS_pET_29b}*placUV5'});
	%%%%%%%%%%%%%
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	mRFP1 = [0.036 0.25 0.343 0.486 0.5 0.536];
	output_sim = simf_doi_10_1186_1754_1611_5_12_Supplement_page19(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; mRFP1];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1186/1754-1611-5-12', 'Supplement-page19', IPTG, output_sim, {'placIq = 15*{alpha_placIq}'; 'lacI = {ALPHA_RBS_placI}*placIq*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'placUV5 = 15*({beta_placUV5} + ({alpha_placUV5} - {beta_placUV5})*{K_placUV5}^{n_placUV5}/({K_placUV5}^{n_placUV5} + lacI^{n_placUV5}))'; 'mRFP1 = {scale_mRFP1_doi_10_1186_1754_1611_5_12}*{ALPHA_RBS_pET_29b}*placUV5'});
	%%%%%%%%%%%%%
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	mRFP1 = [0.014 0.086 0.107 0.143 0.143 0.143];
	output_sim = simf_doi_10_1186_1754_1611_5_12_Supplement_page20(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; mRFP1];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1186/1754-1611-5-12', 'Supplement-page20', IPTG, output_sim, {'placIq = 5*{alpha_placIq}'; 'lacI = {ALPHA_RBS_placI}*placIq*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'placUV5 = 5*({beta_placUV5} + ({alpha_placUV5} - {beta_placUV5})*{K_placUV5}^{n_placUV5}/({K_placUV5}^{n_placUV5} + lacI^{n_placUV5}))'; 'mRFP1 = {scale_mRFP1_doi_10_1186_1754_1611_5_12}*{ALPHA_RBS_pET_29b}*placUV5'});
	%%%%%%%%%%%%%
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	mRFP1 = [0.014 0.286 0.421 0.557 0.614 0.571];
	output_sim = simf_doi_10_1186_1754_1611_5_12_Supplement_page21(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; mRFP1];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1186/1754-1611-5-12', 'Supplement-page21', IPTG, output_sim, {'placIq = 10*{alpha_placIq}'; 'lacI = {ALPHA_RBS_placI}*placIq*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'pLlacO_1 = 10*({beta_pLlacO_1} + ({alpha_pLlacO_1} - {beta_pLlacO_1})*{K_pLlacO_1}^{n_pLlacO_1}/({K_pLlacO_1}^{n_pLlacO_1} + lacI^{n_pLlacO_1}))'; 'mRFP1 = {scale_mRFP1_doi_10_1186_1754_1611_5_12}*{ALPHA_RBS_pET_29b}*pLlacO_1'});
	%%%%%%%%%%%%%
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	GFP_uv = [0.008 0.064 0.105 0.203 0.242 0.339];
	output_sim = simf_doi_10_1186_1754_1611_5_12_Supplement_page22(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; GFP_uv];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1186/1754-1611-5-12', 'Supplement-page22', IPTG, output_sim, {'placIq = 20*{alpha_placIq}'; 'lacI = {ALPHA_RBS_placI}*placIq*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'pLlacO_1 = 20*({beta_pLlacO_1} + ({alpha_pLlacO_1} - {beta_pLlacO_1})*{K_pLlacO_1}^{n_pLlacO_1}/({K_pLlacO_1}^{n_pLlacO_1} + lacI^{n_pLlacO_1}))'; 'GFP_uv = {scale_GFP_uv_doi_10_1186_1754_1611_5_12}*{ALPHA_RBS_pET_29b}*pLlacO_1'});
	%%%%%%%%%%%%%
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	mRFP1 = [0.014 0.25 0.4 0.593 0.65 0.664];
	output_sim = simf_doi_10_1186_1754_1611_5_12_Supplement_page23(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; mRFP1];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1186/1754-1611-5-12', 'Supplement-page23', IPTG, output_sim, {'placIq = 15*{alpha_placIq}'; 'lacI = {ALPHA_RBS_placI}*placIq*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'pLlacO_1 = 15*({beta_pLlacO_1} + ({alpha_pLlacO_1} - {beta_pLlacO_1})*{K_pLlacO_1}^{n_pLlacO_1}/({K_pLlacO_1}^{n_pLlacO_1} + lacI^{n_pLlacO_1}))'; 'mRFP1 = {scale_mRFP1_doi_10_1186_1754_1611_5_12}*{ALPHA_RBS_pET_29b}*pLlacO_1'});
	%%%%%%%%%%%%%
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	mRFP1 = [0.014 0.229 0.307 0.393 0.393 0.414];
	output_sim = simf_doi_10_1186_1754_1611_5_12_Supplement_page24(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; mRFP1];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1186/1754-1611-5-12', 'Supplement-page24', IPTG, output_sim, {'placIq = 5*{alpha_placIq}'; 'lacI = {ALPHA_RBS_placI}*placIq*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'pLlacO_1 = 5*({beta_pLlacO_1} + ({alpha_pLlacO_1} - {beta_pLlacO_1})*{K_pLlacO_1}^{n_pLlacO_1}/({K_pLlacO_1}^{n_pLlacO_1} + lacI^{n_pLlacO_1}))'; 'mRFP1 = {scale_mRFP1_doi_10_1186_1754_1611_5_12}*{ALPHA_RBS_pET_29b}*pLlacO_1'});
	%%%%%%%%%%%%%
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	mRFP1 = [0.129 0.357 0.729 0.843 0.821 0.714];
	output_sim = simf_doi_10_1186_1754_1611_5_12_Supplement_page25(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; mRFP1];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1186/1754-1611-5-12', 'Supplement-page25', IPTG, output_sim, {'placIq = 10*{alpha_placIq}'; 'lacI = {ALPHA_RBS_placI}*placIq*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'pT7lacO = 10*({beta_pT7lacO} + ({alpha_pT7lacO} - {beta_pT7lacO})*{K_pT7lacO}^{n_pT7lacO}/({K_pT7lacO}^{n_pT7lacO} + lacI^{n_pT7lacO}))'; 'mRFP1 = {scale_mRFP1_doi_10_1186_1754_1611_5_12}*{ALPHA_RBS_pET_29b}*pT7lacO'});
	%%%%%%%%%%%%%
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	GFP_uv = [0.034 0.619 0.805 1 0.983 0.881];
	output_sim = simf_doi_10_1186_1754_1611_5_12_Supplement_page26(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; GFP_uv];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1186/1754-1611-5-12', 'Supplement-page26', IPTG, output_sim, {'placIq = 20*{alpha_placIq}'; 'lacI = {ALPHA_RBS_placI}*placIq*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'pT7lacO = 20*({beta_pT7lacO} + ({alpha_pT7lacO} - {beta_pT7lacO})*{K_pT7lacO}^{n_pT7lacO}/({K_pT7lacO}^{n_pT7lacO} + lacI^{n_pT7lacO}))'; 'GFP_uv = {scale_GFP_uv_doi_10_1186_1754_1611_5_12}*{ALPHA_RBS_pET_29b}*pT7lacO'});
	%%%%%%%%%%%%%
	IPTG = [0 0.001 0.004 0.02 0.1 0.5];
	mRFP1 = [0.093 0.321 0.714 0.943 0.943 0.957];
	output_sim = simf_doi_10_1186_1754_1611_5_12_Supplement_page27(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; mRFP1];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1186/1754-1611-5-12', 'Supplement-page27', IPTG, output_sim, {'placIq = 15*{alpha_placIq}'; 'lacI = {ALPHA_RBS_placI}*placIq*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'pT7lacO = 15*({beta_pT7lacO} + ({alpha_pT7lacO} - {beta_pT7lacO})*{K_pT7lacO}^{n_pT7lacO}/({K_pT7lacO}^{n_pT7lacO} + lacI^{n_pT7lacO}))'; 'mRFP1 = {scale_mRFP1_doi_10_1186_1754_1611_5_12}*{ALPHA_RBS_pET_29b}*pT7lacO'});
	%%%%%%%%%%%%%
	IPTG = [0 0.001 0.004 0.02 0.1 0.5];
	mRFP1 = [0.093 0.321 0.629 0.679 0.586 0.486];
	output_sim = simf_doi_10_1186_1754_1611_5_12_Supplement_page28(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; mRFP1];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1186/1754-1611-5-12', 'Supplement-page28', IPTG, output_sim, {'placIq = 5*{alpha_placIq}'; 'lacI = {ALPHA_RBS_placI}*placIq*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'pT7lacO = 5*({beta_pT7lacO} + ({alpha_pT7lacO} - {beta_pT7lacO})*{K_pT7lacO}^{n_pT7lacO}/({K_pT7lacO}^{n_pT7lacO} + lacI^{n_pT7lacO}))'; 'mRFP1 = {scale_mRFP1_doi_10_1186_1754_1611_5_12}*{ALPHA_RBS_pET_29b}*pT7lacO'});
	%%%%%%%%%%%%%
	L_arabinose = [0 0.001 0.01 0.1 1 20];
	mRFP1 = [0 0.004 0.029 0.086 0.143 0.179];
	output_sim = simf_doi_10_1186_1754_1611_5_12_Supplement_page29(sim_x_opt,L_arabinose);
	sim_pair_matrix = [output_sim; mRFP1];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1186/1754-1611-5-12', 'Supplement-page29', L_arabinose, output_sim, {'pC = 1*{alpha_pC}'; 'araC = {ALPHA_RBS_pC}*pC*{K_L_arabinose}^{n_L_arabinose}/({K_L_arabinose}^{n_L_arabinose} + L_arabinose^{n_L_arabinose})'; 'pBAD = 10*({beta_pBAD} + ({alpha_pBAD} - {beta_pBAD})*{K_pBAD}^{n_pBAD}/({K_pBAD}^{n_pBAD} + araC^{n_pBAD}))'; 'mRFP1 = {scale_mRFP1_doi_10_1186_1754_1611_5_12}*{ALPHA_RBS_pET_29b}*pBAD'});
	%%%%%%%%%%%%%
	L_arabinose = [0 0.006 0.013 0.05 0.1 0.5];
	GFP_uv = [0.008 0.017 0.068 0.153 0.305 0.661];
	output_sim = simf_doi_10_1186_1754_1611_5_12_Supplement_page30(sim_x_opt,L_arabinose);
	sim_pair_matrix = [output_sim; GFP_uv];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1186/1754-1611-5-12', 'Supplement-page30', L_arabinose, output_sim, {'pC = 1*{alpha_pC}'; 'araC = {ALPHA_RBS_pC}*pC*{K_L_arabinose}^{n_L_arabinose}/({K_L_arabinose}^{n_L_arabinose} + L_arabinose^{n_L_arabinose})'; 'pBAD = 20*({beta_pBAD} + ({alpha_pBAD} - {beta_pBAD})*{K_pBAD}^{n_pBAD}/({K_pBAD}^{n_pBAD} + araC^{n_pBAD}))'; 'GFP_uv = {scale_GFP_uv_doi_10_1186_1754_1611_5_12}*{ALPHA_RBS_pET_29b}*pBAD'});
	%%%%%%%%%%%%%
	L_arabinose = [0 0.001 0.01 0.1 1 20];
	mRFP1 = [0 0.007 0.05 0.15 0.329 0.464];
	output_sim = simf_doi_10_1186_1754_1611_5_12_Supplement_page31(sim_x_opt,L_arabinose);
	sim_pair_matrix = [output_sim; mRFP1];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1186/1754-1611-5-12', 'Supplement-page31', L_arabinose, output_sim, {'pC = 1*{alpha_pC}'; 'araC = {ALPHA_RBS_pC}*pC*{K_L_arabinose}^{n_L_arabinose}/({K_L_arabinose}^{n_L_arabinose} + L_arabinose^{n_L_arabinose})'; 'pBAD = 15*({beta_pBAD} + ({alpha_pBAD} - {beta_pBAD})*{K_pBAD}^{n_pBAD}/({K_pBAD}^{n_pBAD} + araC^{n_pBAD}))'; 'mRFP1 = {scale_mRFP1_doi_10_1186_1754_1611_5_12}*{ALPHA_RBS_pET_29b}*pBAD'});
	%%%%%%%%%%%%%
	L_arabinose = [0 0.001 0.01 0.1 1 20];
	mRFP1 = [0.004 0.004 0.014 0.032 0.064 0.079];
	output_sim = simf_doi_10_1186_1754_1611_5_12_Supplement_page32(sim_x_opt,L_arabinose);
	sim_pair_matrix = [output_sim; mRFP1];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1186/1754-1611-5-12', 'Supplement-page32', L_arabinose, output_sim, {'pC = 1*{alpha_pC}'; 'araC = {ALPHA_RBS_pC}*pC*{K_L_arabinose}^{n_L_arabinose}/({K_L_arabinose}^{n_L_arabinose} + L_arabinose^{n_L_arabinose})'; 'pBAD = 5*({beta_pBAD} + ({alpha_pBAD} - {beta_pBAD})*{K_pBAD}^{n_pBAD}/({K_pBAD}^{n_pBAD} + araC^{n_pBAD}))'; 'mRFP1 = {scale_mRFP1_doi_10_1186_1754_1611_5_12}*{ALPHA_RBS_pET_29b}*pBAD'});
	%%%%%%%%%%%%%
	IPTG = [0.005 0.01 0.02 0.03 0.05 0.1 0.2 0.3 0.8 1];
	mRFP1 = [0.004 0.038 0.192 0.423 0.75 0.904 0.904 0.908 0.981 1];
	output_sim = simf_doi_10_1371_journal_pone_0039407_figure6B(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; mRFP1];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1371/journal.pone.0039407', 'figure6B', IPTG, output_sim, {'placI = 1*{alpha_placI}'; 'lacI = {ALPHA_RBS_placI}*placI*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'pLAC = 1*({beta_pLAC} + ({alpha_pLAC} - {beta_pLAC})*{K_pLAC}^{n_pLAC}/({K_pLAC}^{n_pLAC} + lacI^{n_pLAC}))'; 'lacZ = {ALPHA_RBS_pLAC}*pLAC'; 'pLlacO_1 = 5*({beta_pLlacO_1} + ({alpha_pLlacO_1} - {beta_pLlacO_1})*{K_pLlacO_1}^{n_pLlacO_1}/({K_pLlacO_1}^{n_pLlacO_1} + lacI^{n_pLlacO_1}))'; 'mRFP1 = {scale_mRFP1_doi_10_1371_journal_pone_0039407}*{ALPHA_BBa_B0034}*pLlacO_1'});
	%%%%%%%%%%%%%
	output_sim = simf_doi_10_1073_pnas_1301301110_sd03_xls_1(sim_x_opt);
	point_pair = [output_sim 	0.422];
	fprintf(sim_pair_fileID, '%f \t %f\n', point_pair);
	print_simulated_data_to_json(json_fileID, 'doi:10.1073/pnas.1301301110', 'sd03.xls-1', [], output_sim, {'BBa_J23117 = 10*{alpha_BBa_J23117}'; 'sfGFP = {scale_sfGFP_doi_10_1073_pnas_1301301110}*{ALPHA_BBa_B0032}*BBa_J23117'});
	%%%%%%%%%%%%%
	output_sim = simf_doi_10_1073_pnas_1301301110_sd03_xls_2(sim_x_opt);
	point_pair = [output_sim 	1];
	fprintf(sim_pair_fileID, '%f \t %f\n', point_pair);
	print_simulated_data_to_json(json_fileID, 'doi:10.1073/pnas.1301301110', 'sd03.xls-2', [], output_sim, {'BBa_J23101 = 10*{alpha_BBa_J23101}'; 'sfGFP = {scale_sfGFP_doi_10_1073_pnas_1301301110}*{ALPHA_BBa_B0032}*BBa_J23101'});
	fclose(sim_pair_fileID);
	fclose(json_fileID);
end

function sim_LOOCV(ensemble_size, ite_num, func_ite_num)
	lb = [0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01 0.03 0.002 0.01 0.01 0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01 0.03 0.01 0.03 0.002 0.01 0.01 0.03 0.002 0.01 0.002 0.002 0.002 0.002 0.002 0.01 0.002 0.002 0.002 0.002 0.01 0.002 0.001 0.002 0 1 0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01 0.03 0.03 0.001 0 1 0.01 0.03 0.001 0.002 0 1 0.01 0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.03 0.001 1 0.01 0.03 0.01 0.002 0.01 0.01 0.03 0.01 0.03 0.01 0.002 0.01 0.03 0.03 0.001 0.002 0 1 0.01 0.01 0.002 0.01 0.01 0.03 0.03 0.001 0.002 0.002 0 1 0.01 0.03 0.03 0.001 0.001 0.002 0.002 0 0 1 1 0.01 0.01 0.01 0.001 0.002 0.002 0 1 0.01 0.001 0.002 0 1 0.001 0.002 0 1 0.03 0.001 0.001 0.002 0 1 1 0.01 0.03 0.001 0.002 0 1 0.01 0.002 0.03 0.001 0.002 0 1 0.01 0.001 0.001 0.001 0.002 0.002 0 0 1 1 1 0.01 0.03 0.01 0.01 0.01 0.03 0.01 0.01 0.03 0.002 0.01 0.001 0.002 0 1 0.01 0.001 0.002 0 1 0.01 0.03 0.01 0.03 0.001 0.002 0.002 0 1 0.03 0.01 0.03 0.01 0.03 0.001 0.002 0 1 0.03 0.001 0.002 0 1 0.01 0.01 0.01 0.03 0.001 0.002 0.002 0 1 0.001 0.002 0 1 0.01 0.01];
	ub = [5 5 50 50 50 50 1 4 4 10 5 50 10 10 5 5 50 50 50 50 1 4 4 10 5 50 50 50 50 1 4 4 10 5 10 5 50 10 10 5 50 10 50 50 50 50 50 10 50 50 50 50 10 50 50 50 1 4 5 5 50 50 50 50 1 4 4 10 5 5 50 1 4 10 5 50 50 1 4 10 5 5 50 50 50 50 1 4 4 5 50 4 10 5 10 50 10 10 5 10 5 10 50 10 5 5 50 50 1 4 10 10 50 10 10 5 5 50 50 50 1 4 10 5 5 50 50 50 50 1 1 4 4 10 10 10 50 50 50 1 4 10 50 50 1 4 50 50 1 4 5 50 50 50 1 4 4 10 5 50 50 1 4 10 50 5 50 50 1 4 10 50 50 50 50 50 1 1 4 4 4 10 5 10 10 10 5 10 10 5 50 10 50 50 1 4 10 50 50 1 4 10 5 10 5 50 50 50 1 4 5 10 5 10 5 50 50 1 4 5 50 50 1 4 10 10 10 5 50 50 50 1 4 50 50 1 4 10 10];
	LOOCV_sim_data_fileID = fopen('Plot_data/LOOCV_sim_data.dat','w');
	LOO_model_name_list = {'doi_10_1186_1754_1611_3_4_figure3_1', 'doi_10_1186_1754_1611_3_4_figure3_2', 'doi_10_1186_1754_1611_3_4_figure3_6', 'doi_10_1186_1754_1611_3_4_figure3_7', 'doi_10_1093_nar_gku593_figure3D_1', 'doi_10_1093_nar_gku593_figure3D_2', 'doi_10_1093_nar_gku593_figure3D_3', 'doi_10_1093_nar_gku593_figure3D_4', 'doi_10_1093_nar_gku593_figure3D_5', 'doi_10_1093_nar_gku593_figure3D_6', 'doi_10_1093_nar_gku593_figureS6', 'doi_10_1093_nar_gku1388_figure2A_1', 'doi_10_1093_nar_gku1388_figure2A_2', 'doi_10_1093_nar_gku1388_figure2A_3', 'doi_10_1093_nar_gku1388_figure2A_4', 'doi_10_1093_nar_gku1388_figure2A_5', 'doi_10_1093_nar_gku1388_figure2A_6', 'doi_10_1093_nar_gku1388_figure2B_1', 'doi_10_1093_nar_gku1388_figure2B_2', 'doi_10_1093_nar_gku964_figure3A_1', 'doi_10_1093_nar_gku964_figure3A_2', 'doi_10_1093_nar_gku964_figure3A_3', 'doi_10_1021_sb300055e_figure4a_1', 'doi_10_1021_sb300055e_figure4a_2', 'doi_10_1021_sb300055e_figure4a_3', 'doi_10_1186_1752_0509_5_111_figure3_a', 'doi_10_1186_1752_0509_5_111_figure3_b', 'doi_10_1038_nmeth_2205_figureS1_2', 'doi_10_1038_nmeth_2205_figureS1_3', 'doi_10_1038_nmeth_2205_figureS1_41', 'doi_10_1038_nmeth_2205_figureS1_42', 'doi_10_1038_nchembio_1411_figure4_1', 'doi_10_1021_sb400131a_figure4c_1', 'doi_10_1021_sb400131a_figure4c_2', 'doi_10_1186_1471_2105_13_S4_S11_figure2_1', 'doi_10_1186_1471_2105_13_S4_S11_figure2_2', 'doi_10_1186_1471_2105_13_S4_S11_figure2_3', 'doi_10_1038_nature03461_figure2b_1', 'doi_10_1038_nature03461_figure2b_2', 'doi_10_1038_nature03461_figure2b_3', 'doi_10_1038_nature09565_figure1C', 'doi_10_1038_nature09565_figureS1C', 'doi_10_1093_nar_gks597_figureS3_left', 'doi_10_1038_nbt_2401_figureS11', 'doi_10_1038_nbt_2401_figureS13', 'doi_10_1038_nature11516_figureS8A', 'doi_10_1038_nature11516_figureS8B', 'doi_10_1038_nature11516_figureS8C', 'doi_10_1186_1754_1611_5_12_Supplement_page1', 'doi_10_1186_1754_1611_5_12_Supplement_page2', 'doi_10_1186_1754_1611_5_12_Supplement_page3', 'doi_10_1186_1754_1611_5_12_Supplement_page4', 'doi_10_1186_1754_1611_5_12_Supplement_page5', 'doi_10_1186_1754_1611_5_12_Supplement_page6', 'doi_10_1186_1754_1611_5_12_Supplement_page7', 'doi_10_1186_1754_1611_5_12_Supplement_page8', 'doi_10_1186_1754_1611_5_12_Supplement_page17', 'doi_10_1186_1754_1611_5_12_Supplement_page18', 'doi_10_1186_1754_1611_5_12_Supplement_page19', 'doi_10_1186_1754_1611_5_12_Supplement_page20', 'doi_10_1186_1754_1611_5_12_Supplement_page21', 'doi_10_1186_1754_1611_5_12_Supplement_page22', 'doi_10_1186_1754_1611_5_12_Supplement_page23', 'doi_10_1186_1754_1611_5_12_Supplement_page24', 'doi_10_1186_1754_1611_5_12_Supplement_page25', 'doi_10_1186_1754_1611_5_12_Supplement_page26', 'doi_10_1186_1754_1611_5_12_Supplement_page27', 'doi_10_1186_1754_1611_5_12_Supplement_page28', 'doi_10_1186_1754_1611_5_12_Supplement_page29', 'doi_10_1186_1754_1611_5_12_Supplement_page30', 'doi_10_1186_1754_1611_5_12_Supplement_page31', 'doi_10_1186_1754_1611_5_12_Supplement_page32', 'doi_10_1073_pnas_1301301110_sd03_xls_1', 'doi_10_1073_pnas_1301301110_sd03_xls_2'};
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
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1093_nar_gku593_figure3D_1, lb, ub);
	desired_output = 0.002;
	output = simf_doi_10_1093_nar_gku593_figure3D_1(x_opt);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1093_nar_gku593_figure3D_1\n');
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
	L_arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2];
	desired_output = [0 0.002 0.041 0.479 0.949 0.997 1 1 1];
	output = simf_doi_10_1093_nar_gku593_figureS6(x_opt, L_arabinose);
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
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1093_nar_gku964_figure3A_1, lb, ub);
	AHL = [0 0 0 0 0.001 0.01 0.1];
	desired_output = [0.024 0.081 0.372 0.763 0.874 0.889 0.891];
	output = simf_doi_10_1093_nar_gku964_figure3A_1(x_opt, AHL);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1093_nar_gku964_figure3A_1\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1093_nar_gku964_figure3A_2, lb, ub);
	AHL = [0 0 0 0 0.001 0.01 0.1];
	desired_output = [0.015 0.052 0.295 0.755 0.911 0.931 0.934];
	output = simf_doi_10_1093_nar_gku964_figure3A_2(x_opt, AHL);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1093_nar_gku964_figure3A_2\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1093_nar_gku964_figure3A_3, lb, ub);
	AHL = [0 0 0 0 0.001 0.01 0.1];
	desired_output = [0.006 0.01 0.05 0.319 0.82 0.98 1];
	output = simf_doi_10_1093_nar_gku964_figure3A_3(x_opt, AHL);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1093_nar_gku964_figure3A_3\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1021_sb300055e_figure4a_1, lb, ub);
	desired_output = 1;
	output = simf_doi_10_1021_sb300055e_figure4a_1(x_opt);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1021_sb300055e_figure4a_1\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1021_sb300055e_figure4a_2, lb, ub);
	desired_output = 0.429;
	output = simf_doi_10_1021_sb300055e_figure4a_2(x_opt);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1021_sb300055e_figure4a_2\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1021_sb300055e_figure4a_3, lb, ub);
	desired_output = 0.119;
	output = simf_doi_10_1021_sb300055e_figure4a_3(x_opt);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1021_sb300055e_figure4a_3\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1186_1752_0509_5_111_figure3_a, lb, ub);
	L_arabinose = [0.03 0.08 0.2 0.5 1.2 3 8 20 50 125 300];
	desired_output = [0.014 0.083 0.181 0.361 0.5 0.694 0.778 0.889 1 0.917 0.778];
	output = simf_doi_10_1186_1752_0509_5_111_figure3_a(x_opt, L_arabinose);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1186_1752_0509_5_111_figure3_a\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1186_1752_0509_5_111_figure3_b, lb, ub);
	L_arabinose = [0.03 0.08 0.2 0.5 1.2 3 8 20 50 125 300];
	desired_output = [0.014 0.014 0.014 0.017 0.021 0.028 0.061 0.167 0.389 0.597 0.556];
	output = simf_doi_10_1186_1752_0509_5_111_figure3_b(x_opt, L_arabinose);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1186_1752_0509_5_111_figure3_b\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1038_nmeth_2205_figureS1_2, lb, ub);
	L_arabinose = [0.017 0.027 0.11 0.27 0.67 1.3 6.7];
	desired_output = [0.007 0.009 0.017 0.06 0.121 0.115 0.103];
	output = simf_doi_10_1038_nmeth_2205_figureS1_2(x_opt, L_arabinose);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1038_nmeth_2205_figureS1_2\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1038_nmeth_2205_figureS1_3, lb, ub);
	aTc = [0.25 1 4 16 64 240];
	desired_output = [0.002 0.002 0.003 0.009 0.036 0.043];
	output = simf_doi_10_1038_nmeth_2205_figureS1_3(x_opt, aTc);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1038_nmeth_2205_figureS1_3\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1038_nmeth_2205_figureS1_41, lb, ub);
	IPTG = [0 0 0.001 0.01 0.1 1];
	desired_output = [0.009 0.034 0.069 0.086 0.095 0.121];
	output = simf_doi_10_1038_nmeth_2205_figureS1_41(x_opt, IPTG);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1038_nmeth_2205_figureS1_41\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1038_nmeth_2205_figureS1_42, lb, ub);
	IPTG = [0 0 0.001 0.01 0.1 1];
	desired_output = [0.129 0.129 0.129 0.345 0.905 1];
	output = simf_doi_10_1038_nmeth_2205_figureS1_42(x_opt, IPTG);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1038_nmeth_2205_figureS1_42\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1038_nchembio_1411_figure4_1, lb, ub);
	IPTG = [0 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1 0.15 0.2 1];
	desired_output = [0.001 0.002 0.003 0.005 0.008 0.016 0.024 0.032 0.06 0.1 0.14 0.2];
	output = simf_doi_10_1038_nchembio_1411_figure4_1(x_opt, IPTG);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1038_nchembio_1411_figure4_1\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1021_sb400131a_figure4c_1, lb, ub);
	aTc = [0.46 4.6 46 460 4600];
	desired_output = [0.025 0.025 0.075 0.25 0.8];
	output = simf_doi_10_1021_sb400131a_figure4c_1(x_opt, aTc);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1021_sb400131a_figure4c_1\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1021_sb400131a_figure4c_2, lb, ub);
	IPTG = [0 0 0.001 0.01 0.1 1];
	desired_output = [0.075 0.075 0.075 0.25 0.875 1];
	output = simf_doi_10_1021_sb400131a_figure4c_2(x_opt, IPTG);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1021_sb400131a_figure4c_2\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1186_1471_2105_13_S4_S11_figure2_1, lb, ub);
	AHL = [0 0 0 0 0 0 0 0 0 0 0];
	desired_output = [0.006 0.086 0.126 0.343 0.571 0.657 0.857 0.914 0.943 0.971 1];
	output = simf_doi_10_1186_1471_2105_13_S4_S11_figure2_1(x_opt, AHL);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1186_1471_2105_13_S4_S11_figure2_1\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1186_1471_2105_13_S4_S11_figure2_2, lb, ub);
	AHL = [0 0 0 0 0 0 0 0 0 0 0];
	desired_output = [0.001 0.011 0.017 0.034 0.057 0.069 0.086 0.097 0.103 0.109 0.114];
	output = simf_doi_10_1186_1471_2105_13_S4_S11_figure2_2(x_opt, AHL);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1186_1471_2105_13_S4_S11_figure2_2\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1186_1471_2105_13_S4_S11_figure2_3, lb, ub);
	AHL = [0 0 0 0 0 0 0 0 0 0 0];
	desired_output = [0.001 0.001 0.001 0.002 0.007 0.014 0.026 0.03 0.033 0.034 0.035];
	output = simf_doi_10_1186_1471_2105_13_S4_S11_figure2_3(x_opt, AHL);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1186_1471_2105_13_S4_S11_figure2_3\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1038_nature03461_figure2b_1, lb, ub);
	AHL = [0.1 0.3 0.5 1 2 4 10 100 1000];
	desired_output = [1 0.947 0.716 0.579 0.158 0.158 0.105 0.053 0.011];
	output = simf_doi_10_1038_nature03461_figure2b_1(x_opt, AHL);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1038_nature03461_figure2b_1\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1038_nature03461_figure2b_2, lb, ub);
	AHL = [0.1 0.3 0.7 1.2 2 5 10 20 30 50 100 500 1000];
	desired_output = [0.947 0.916 0.937 0.916 0.905 0.842 0.684 0.526 0.368 0.284 0.242 0.137 0.126];
	output = simf_doi_10_1038_nature03461_figure2b_2(x_opt, AHL);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1038_nature03461_figure2b_2\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1038_nature03461_figure2b_3, lb, ub);
	AHL = [0.2 0.5 1 2 5 10 50 100 200 300 400 600 1000 2000 3000];
	desired_output = [0.842 0.832 0.8 0.789 0.737 0.695 0.632 0.611 0.526 0.421 0.316 0.263 0.232 0.179 0.126];
	output = simf_doi_10_1038_nature03461_figure2b_3(x_opt, AHL);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1038_nature03461_figure2b_3\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1038_nature09565_figure1C, lb, ub);
	L_arabinose = [0 0.001 0.015 0.05 0.5 5 10];
	desired_output = [0.002 0.002 0.057 0.628 0.999 1 1];
	output = simf_doi_10_1038_nature09565_figure1C(x_opt, L_arabinose);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1038_nature09565_figure1C\n');
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
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1093_nar_gks597_figureS3_left, lb, ub);
	IPTG = [0 0.002 0.01 0.05 0.1 0.25 1];
	desired_output = [0.054 0.107 0.179 0.571 0.693 0.893 1];
	output = simf_doi_10_1093_nar_gks597_figureS3_left(x_opt, IPTG);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1093_nar_gks597_figureS3_left\n');
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
	L_arabinose = [0.1 1 2 5 7 10 12.5 25 37.5 50 62.5];
	desired_output = [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425];
	output = simf_doi_10_1038_nbt_2401_figureS13(x_opt, L_arabinose);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1038_nbt_2401_figureS13\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1038_nature11516_figureS8A, lb, ub);
	L_arabinose = [0 0 0.002 0.01 0.04 0.16 1 4 25];
	desired_output = [0.002 0.002 0.002 0.01 0.137 0.487 0.557 0.558 0.558];
	output = simf_doi_10_1038_nature11516_figureS8A(x_opt, L_arabinose);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1038_nature11516_figureS8A\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1038_nature11516_figureS8B, lb, ub);
	IPTG = [0 0 0 0.001 0.004 0.025 0.1 0.5 2.5];
	desired_output = [0.016 0.016 0.016 0.018 0.04 0.273 0.525 0.588 0.593];
	output = simf_doi_10_1038_nature11516_figureS8B(x_opt, IPTG);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1038_nature11516_figureS8B\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1038_nature11516_figureS8C, lb, ub);
	AHL = [0 0 0 0 0 0 0 0 0];
	desired_output = [0.102 0.102 0.128 0.329 0.89 0.988 0.999 1 1];
	output = simf_doi_10_1038_nature11516_figureS8C(x_opt, AHL);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1038_nature11516_figureS8C\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page1, lb, ub);
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	desired_output = [0.214 0.679 0.764 0.821 0.821 0.836];
	output = simf_doi_10_1186_1754_1611_5_12_Supplement_page1(x_opt, IPTG);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1186_1754_1611_5_12_Supplement_page1\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page2, lb, ub);
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	desired_output = [0.051 0.492 0.61 0.678 0.746 0.847];
	output = simf_doi_10_1186_1754_1611_5_12_Supplement_page2(x_opt, IPTG);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1186_1754_1611_5_12_Supplement_page2\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page3, lb, ub);
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	desired_output = [0.296 0.914 0.964 1 1 1];
	output = simf_doi_10_1186_1754_1611_5_12_Supplement_page3(x_opt, IPTG);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1186_1754_1611_5_12_Supplement_page3\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page4, lb, ub);
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	desired_output = [0.129 0.407 0.443 0.5 0.5 0.5];
	output = simf_doi_10_1186_1754_1611_5_12_Supplement_page4(x_opt, IPTG);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1186_1754_1611_5_12_Supplement_page4\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page5, lb, ub);
	aTc = [0 0.013 0.025 0.05 0.1 0.2];
	desired_output = [0.007 0.429 0.714 0.929 0.986 0.964];
	output = simf_doi_10_1186_1754_1611_5_12_Supplement_page5(x_opt, aTc);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1186_1754_1611_5_12_Supplement_page5\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page6, lb, ub);
	aTc = [0 0.013 0.025 0.05 0.1 0.2];
	desired_output = [0.002 0.288 0.542 0.797 0.932 0.966];
	output = simf_doi_10_1186_1754_1611_5_12_Supplement_page6(x_opt, aTc);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1186_1754_1611_5_12_Supplement_page6\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page7, lb, ub);
	aTc = [0 0.013 0.025 0.05 0.1 0.2];
	desired_output = [0.007 0.271 0.643 0.807 0.929 0.893];
	output = simf_doi_10_1186_1754_1611_5_12_Supplement_page7(x_opt, aTc);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1186_1754_1611_5_12_Supplement_page7\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page8, lb, ub);
	aTc = [0 0.013 0.025 0.05 0.1 0.2];
	desired_output = [0.007 0.429 0.557 0.571 0.557 0.543];
	output = simf_doi_10_1186_1754_1611_5_12_Supplement_page8(x_opt, aTc);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1186_1754_1611_5_12_Supplement_page8\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page17, lb, ub);
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	desired_output = [0.021 0.157 0.2 0.286 0.3 0.279];
	output = simf_doi_10_1186_1754_1611_5_12_Supplement_page17(x_opt, IPTG);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1186_1754_1611_5_12_Supplement_page17\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page18, lb, ub);
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	desired_output = [0.017 0.09 0.115 0.153 0.161 0.183];
	output = simf_doi_10_1186_1754_1611_5_12_Supplement_page18(x_opt, IPTG);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1186_1754_1611_5_12_Supplement_page18\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page19, lb, ub);
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	desired_output = [0.036 0.25 0.343 0.486 0.5 0.536];
	output = simf_doi_10_1186_1754_1611_5_12_Supplement_page19(x_opt, IPTG);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1186_1754_1611_5_12_Supplement_page19\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page20, lb, ub);
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	desired_output = [0.014 0.086 0.107 0.143 0.143 0.143];
	output = simf_doi_10_1186_1754_1611_5_12_Supplement_page20(x_opt, IPTG);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1186_1754_1611_5_12_Supplement_page20\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page21, lb, ub);
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	desired_output = [0.014 0.286 0.421 0.557 0.614 0.571];
	output = simf_doi_10_1186_1754_1611_5_12_Supplement_page21(x_opt, IPTG);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1186_1754_1611_5_12_Supplement_page21\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page22, lb, ub);
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	desired_output = [0.008 0.064 0.105 0.203 0.242 0.339];
	output = simf_doi_10_1186_1754_1611_5_12_Supplement_page22(x_opt, IPTG);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1186_1754_1611_5_12_Supplement_page22\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page23, lb, ub);
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	desired_output = [0.014 0.25 0.4 0.593 0.65 0.664];
	output = simf_doi_10_1186_1754_1611_5_12_Supplement_page23(x_opt, IPTG);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1186_1754_1611_5_12_Supplement_page23\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page24, lb, ub);
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	desired_output = [0.014 0.229 0.307 0.393 0.393 0.414];
	output = simf_doi_10_1186_1754_1611_5_12_Supplement_page24(x_opt, IPTG);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1186_1754_1611_5_12_Supplement_page24\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page25, lb, ub);
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	desired_output = [0.129 0.357 0.729 0.843 0.821 0.714];
	output = simf_doi_10_1186_1754_1611_5_12_Supplement_page25(x_opt, IPTG);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1186_1754_1611_5_12_Supplement_page25\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page26, lb, ub);
	IPTG = [0 0.006 0.013 0.05 0.1 0.5];
	desired_output = [0.034 0.619 0.805 1 0.983 0.881];
	output = simf_doi_10_1186_1754_1611_5_12_Supplement_page26(x_opt, IPTG);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1186_1754_1611_5_12_Supplement_page26\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page27, lb, ub);
	IPTG = [0 0.001 0.004 0.02 0.1 0.5];
	desired_output = [0.093 0.321 0.714 0.943 0.943 0.957];
	output = simf_doi_10_1186_1754_1611_5_12_Supplement_page27(x_opt, IPTG);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1186_1754_1611_5_12_Supplement_page27\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page28, lb, ub);
	IPTG = [0 0.001 0.004 0.02 0.1 0.5];
	desired_output = [0.093 0.321 0.629 0.679 0.586 0.486];
	output = simf_doi_10_1186_1754_1611_5_12_Supplement_page28(x_opt, IPTG);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1186_1754_1611_5_12_Supplement_page28\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page29, lb, ub);
	L_arabinose = [0 0.001 0.01 0.1 1 20];
	desired_output = [0 0.004 0.029 0.086 0.143 0.179];
	output = simf_doi_10_1186_1754_1611_5_12_Supplement_page29(x_opt, L_arabinose);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1186_1754_1611_5_12_Supplement_page29\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page30, lb, ub);
	L_arabinose = [0 0.006 0.013 0.05 0.1 0.5];
	desired_output = [0.008 0.017 0.068 0.153 0.305 0.661];
	output = simf_doi_10_1186_1754_1611_5_12_Supplement_page30(x_opt, L_arabinose);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1186_1754_1611_5_12_Supplement_page30\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page31, lb, ub);
	L_arabinose = [0 0.001 0.01 0.1 1 20];
	desired_output = [0 0.007 0.05 0.15 0.329 0.464];
	output = simf_doi_10_1186_1754_1611_5_12_Supplement_page31(x_opt, L_arabinose);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1186_1754_1611_5_12_Supplement_page31\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page32, lb, ub);
	L_arabinose = [0 0.001 0.01 0.1 1 20];
	desired_output = [0.004 0.004 0.014 0.032 0.064 0.079];
	output = simf_doi_10_1186_1754_1611_5_12_Supplement_page32(x_opt, L_arabinose);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1186_1754_1611_5_12_Supplement_page32\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1073_pnas_1301301110_sd03_xls_1, lb, ub);
	desired_output = 0.422;
	output = simf_doi_10_1073_pnas_1301301110_sd03_xls_1(x_opt);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1073_pnas_1301301110_sd03_xls_1\n');
	for i = 1:length(output)
		fprintf(LOOCV_sim_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_doi_10_1073_pnas_1301301110_sd03_xls_2, lb, ub);
	desired_output = 1;
	output = simf_doi_10_1073_pnas_1301301110_sd03_xls_2(x_opt);
	fprintf(LOOCV_sim_data_fileID, '#doi_10_1073_pnas_1301301110_sd03_xls_2\n');
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
function error = LOOCV_sim_error_doi_10_1093_nar_gku593_figure3D_1(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 4, 4);
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku593_figure3D_2(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 5, 5);
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku593_figure3D_3(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 6, 6);
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku593_figure3D_4(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 7, 7);
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku593_figure3D_5(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 8, 8);
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku593_figure3D_6(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 9, 9);
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku593_figureS6(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 10, 18);
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2A_1(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 19, 19);
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2A_2(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 20, 20);
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2A_3(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 21, 21);
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2A_4(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 22, 22);
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2A_5(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 23, 23);
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2A_6(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 24, 24);
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2B_1(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 25, 31);
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2B_2(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 32, 36);
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku964_figure3A_1(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 37, 43);
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku964_figure3A_2(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 44, 50);
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku964_figure3A_3(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 51, 57);
end
function error = LOOCV_sim_error_doi_10_1021_sb300055e_figure4a_1(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 58, 58);
end
function error = LOOCV_sim_error_doi_10_1021_sb300055e_figure4a_2(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 59, 59);
end
function error = LOOCV_sim_error_doi_10_1021_sb300055e_figure4a_3(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 60, 60);
end
function error = LOOCV_sim_error_doi_10_1186_1752_0509_5_111_figure3_a(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 61, 71);
end
function error = LOOCV_sim_error_doi_10_1186_1752_0509_5_111_figure3_b(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 72, 82);
end
function error = LOOCV_sim_error_doi_10_1038_nmeth_2205_figureS1_2(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 83, 89);
end
function error = LOOCV_sim_error_doi_10_1038_nmeth_2205_figureS1_3(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 90, 95);
end
function error = LOOCV_sim_error_doi_10_1038_nmeth_2205_figureS1_41(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 96, 101);
end
function error = LOOCV_sim_error_doi_10_1038_nmeth_2205_figureS1_42(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 102, 107);
end
function error = LOOCV_sim_error_doi_10_1038_nchembio_1411_figure4_1(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 108, 119);
end
function error = LOOCV_sim_error_doi_10_1021_sb400131a_figure4c_1(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 120, 124);
end
function error = LOOCV_sim_error_doi_10_1021_sb400131a_figure4c_2(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 125, 130);
end
function error = LOOCV_sim_error_doi_10_1186_1471_2105_13_S4_S11_figure2_1(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 131, 141);
end
function error = LOOCV_sim_error_doi_10_1186_1471_2105_13_S4_S11_figure2_2(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 142, 152);
end
function error = LOOCV_sim_error_doi_10_1186_1471_2105_13_S4_S11_figure2_3(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 153, 163);
end
function error = LOOCV_sim_error_doi_10_1038_nature03461_figure2b_1(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 164, 172);
end
function error = LOOCV_sim_error_doi_10_1038_nature03461_figure2b_2(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 173, 185);
end
function error = LOOCV_sim_error_doi_10_1038_nature03461_figure2b_3(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 186, 200);
end
function error = LOOCV_sim_error_doi_10_1038_nature09565_figure1C(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 201, 207);
end
function error = LOOCV_sim_error_doi_10_1038_nature09565_figureS1C(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 208, 214);
end
function error = LOOCV_sim_error_doi_10_1093_nar_gks597_figureS3_left(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 215, 221);
end
function error = LOOCV_sim_error_doi_10_1038_nbt_2401_figureS11(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 222, 231);
end
function error = LOOCV_sim_error_doi_10_1038_nbt_2401_figureS13(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 232, 242);
end
function error = LOOCV_sim_error_doi_10_1038_nature11516_figureS8A(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 243, 251);
end
function error = LOOCV_sim_error_doi_10_1038_nature11516_figureS8B(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 252, 260);
end
function error = LOOCV_sim_error_doi_10_1038_nature11516_figureS8C(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 261, 269);
end
function error = LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page1(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 270, 275);
end
function error = LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page2(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 276, 281);
end
function error = LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page3(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 282, 287);
end
function error = LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page4(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 288, 293);
end
function error = LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page5(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 294, 299);
end
function error = LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page6(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 300, 305);
end
function error = LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page7(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 306, 311);
end
function error = LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page8(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 312, 317);
end
function error = LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page17(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 318, 323);
end
function error = LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page18(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 324, 329);
end
function error = LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page19(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 330, 335);
end
function error = LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page20(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 336, 341);
end
function error = LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page21(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 342, 347);
end
function error = LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page22(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 348, 353);
end
function error = LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page23(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 354, 359);
end
function error = LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page24(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 360, 365);
end
function error = LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page25(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 366, 371);
end
function error = LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page26(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 372, 377);
end
function error = LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page27(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 378, 383);
end
function error = LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page28(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 384, 389);
end
function error = LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page29(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 390, 395);
end
function error = LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page30(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 396, 401);
end
function error = LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page31(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 402, 407);
end
function error = LOOCV_sim_error_doi_10_1186_1754_1611_5_12_Supplement_page32(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 408, 413);
end
function error = LOOCV_sim_error_doi_10_1073_pnas_1301301110_sd03_xls_1(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 414, 414);
end
function error = LOOCV_sim_error_doi_10_1073_pnas_1301301110_sd03_xls_2(x)
	error_all = simultaneous_fitting_error(x) ;
	error = remove_a_sub_array(error_all, 415, 415);
end

