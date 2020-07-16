function real_benchmark_fitting
	%sim_x_opt = simultaneous_fitting(1,1000,50000);
	%print_plotting_sim_data(sim_x_opt);
	%sim_LOOCV(10,100,5000);
	seq_x_opt = sequential_fitting(10, 10);
	seq_LOOCV(10,10);
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
function output = simf_doi_10_1128_AEM_00791_07_figure3a_error(x)
	output = repmat(simf_doi_10_1128_AEM_00791_07_figure3a(x,[0.003 0.005 0.05 0.5 1]) - [0.044 0.089 0.556 1 1],1,5);
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
function output = simf_doi_10_1038_nature12148_figure16a_error(x)
	output = repmat(simf_doi_10_1038_nature12148_figure16a(x,[0 0.001 0.003 0.008 0.025 0.075 0.2 0.7 2 6 20 50]) - [0.015 0.015 0.015 0.02 0.06 0.48 0.84 0.93 0.96 0.98 0.98 1],1,5);
end

function output = simf_doi__10_1093_nar_25_6_1203_figure4a(x, input)
	output = zeros(1,length(input));
	ALPHA_RBSII = x(14);
	ALPHA_RBS_pN25 = x(15);
	K_aTc = x(16);
	K_pLtetO_1 = x(17);
	alpha_pLtetO_1 = x(18);
	alpha_pN25 = x(19);
	beta_pLtetO_1 = x(20);
	n_aTc = x(21);
	n_pLtetO_1 = x(22);
	scale_luc_doi__10_1093_nar_25_6_1203 = x(23);
	for i = 1:length(input)
		pN25 = 1*alpha_pN25;
		aTc = input(i);
		tetR = ALPHA_RBS_pN25*pN25*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pLtetO_1 = 15*(beta_pLtetO_1 + (alpha_pLtetO_1 - beta_pLtetO_1)*K_pLtetO_1^n_pLtetO_1/(K_pLtetO_1^n_pLtetO_1 + tetR^n_pLtetO_1));
		luc = scale_luc_doi__10_1093_nar_25_6_1203*ALPHA_RBSII*pLtetO_1;
		output(i) = luc;
	end
end
function output = simf_doi__10_1093_nar_25_6_1203_figure4a_error(x)
	output = repmat(simf_doi__10_1093_nar_25_6_1203_figure4a(x,[0 1 2 3 4 5 6 7 8 9 10 12 13 17 18 20 25 35 60]) - [0 0.001 0.001 0.001 0.002 0.004 0.007 0.011 0.014 0.032 0.036 0.054 0.089 0.643 0.893 1 1 1 1],1,5);
end

function output = simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6(x, input)
	output = zeros(1,length(input));
	ALPHA_RBSII = x(14);
	ALPHA_RBS_placI = x(24);
	K_IPTG = x(25);
	K_pLAC = x(26);
	alpha_pLAC = x(27);
	alpha_placIq = x(28);
	beta_pLAC = x(29);
	n_IPTG = x(30);
	n_pLAC = x(31);
	scale_eYFP_http___dspace_mit_edu_handle_1721_1_8228 = x(32);
	for i = 1:length(input)
		placIq = 10*alpha_placIq;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placIq*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pLAC = 10*(beta_pLAC + (alpha_pLAC - beta_pLAC)*K_pLAC^n_pLAC/(K_pLAC^n_pLAC + lacI^n_pLAC));
		eYFP = scale_eYFP_http___dspace_mit_edu_handle_1721_1_8228*ALPHA_RBSII*pLAC;
		output(i) = eYFP;
	end
end
function output = simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6_error(x)
	output = repmat(simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6(x,[0.001 0.01 0.03 0.05 0.075 0.1 0.15 0.2 0.3 0.5 1]) - [0.015 0.014 0.023 0.062 0.098 0.187 0.239 0.374 0.572 0.869 1],1,5);
end

function output = simf_doi_10_1073pnas_0408507102_figure2A(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_Unknown_Hooshangi = x(33);
	ALPHA_RBS_placI = x(24);
	K_aTc = x(16);
	K_pLtetO_1 = x(17);
	alpha_pLtetO_1 = x(18);
	alpha_placIq = x(28);
	beta_pLtetO_1 = x(20);
	n_aTc = x(21);
	n_pLtetO_1 = x(22);
	scale_eYFP_doi_10_1073pnas_0408507102 = x(34);
	for i = 1:length(input)
		placIq = 10*alpha_placIq;
		aTc = input(i);
		tetR = ALPHA_RBS_placI*placIq*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pLtetO_1 = 10*(beta_pLtetO_1 + (alpha_pLtetO_1 - beta_pLtetO_1)*K_pLtetO_1^n_pLtetO_1/(K_pLtetO_1^n_pLtetO_1 + tetR^n_pLtetO_1));
		eYFP = scale_eYFP_doi_10_1073pnas_0408507102*ALPHA_RBS_Unknown_Hooshangi*pLtetO_1;
		output(i) = eYFP;
	end
end
function output = simf_doi_10_1073pnas_0408507102_figure2A_error(x)
	output = repmat(simf_doi_10_1073pnas_0408507102_figure2A(x,[9.258 13.887 20.831 32.403 46.29 83.322 162.015 231.45 370.32 648.06 925.8 1388.7]) - [0.003 0.003 0.003 0.003 0.005 0.01 0.05 0.167 0.333 0.667 0.833 1],1,5);
end

function output = simf_doi_10_1073_pnas_0606717104_figure4_6(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pLAC = x(35);
	ALPHA_RBS_placI = x(24);
	K_IPTG = x(25);
	K_pLAC = x(26);
	alpha_pLAC = x(27);
	alpha_placI = x(36);
	beta_pLAC = x(29);
	n_IPTG = x(30);
	n_pLAC = x(31);
	scale_lacZ_doi_10_1073_pnas_0606717104 = x(37);
	for i = 1:length(input)
		placI = 1*alpha_placI;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placI*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pLAC = 1*(beta_pLAC + (alpha_pLAC - beta_pLAC)*K_pLAC^n_pLAC/(K_pLAC^n_pLAC + lacI^n_pLAC));
		lacZ = scale_lacZ_doi_10_1073_pnas_0606717104*ALPHA_RBS_pLAC*pLAC;
		output(i) = lacZ;
	end
end
function output = simf_doi_10_1073_pnas_0606717104_figure4_6_error(x)
	output = repmat(simf_doi_10_1073_pnas_0606717104_figure4_6(x,[0 0.005 0.01 0.02 0.05 0.1 0.2 1]) - [0.001 0.003 0.043 0.5 0.984 0.999 1 1],1,5);
end

function output = simf_doi_10_1186_1754_1611_3_4_figure3_1(x)
	ALPHA_BBa_B0032 = x(38);
	alpha_BBa_J23101 = x(39);
	scale_GFPmut3b_doi_10_1186_1754_1611_3_4 = x(40);
	BBa_J23101 = 10*alpha_BBa_J23101;
	GFPmut3b = scale_GFPmut3b_doi_10_1186_1754_1611_3_4*ALPHA_BBa_B0032*BBa_J23101;
	output = GFPmut3b;
end
function output = simf_doi_10_1186_1754_1611_3_4_figure3_1_error(x)
	output = repmat(simf_doi_10_1186_1754_1611_3_4_figure3_1(x) - 0.69,1,5);
end

function output = simf_doi_10_1186_1754_1611_3_4_figure3_2(x)
	ALPHA_BBa_B0032 = x(38);
	alpha_BBa_J23116 = x(41);
	scale_GFPmut3b_doi_10_1186_1754_1611_3_4 = x(40);
	BBa_J23116 = 10*alpha_BBa_J23116;
	GFPmut3b = scale_GFPmut3b_doi_10_1186_1754_1611_3_4*ALPHA_BBa_B0032*BBa_J23116;
	output = GFPmut3b;
end
function output = simf_doi_10_1186_1754_1611_3_4_figure3_2_error(x)
	output = repmat(simf_doi_10_1186_1754_1611_3_4_figure3_2(x) - 0.018,1,5);
end

function output = simf_doi_10_1186_1754_1611_3_4_figure3_3(x)
	ALPHA_BBa_B0032 = x(38);
	alpha_BBa_J23150 = x(42);
	scale_GFPmut3b_doi_10_1186_1754_1611_3_4 = x(40);
	BBa_J23150 = 10*alpha_BBa_J23150;
	GFPmut3b = scale_GFPmut3b_doi_10_1186_1754_1611_3_4*ALPHA_BBa_B0032*BBa_J23150;
	output = GFPmut3b;
end
function output = simf_doi_10_1186_1754_1611_3_4_figure3_3_error(x)
	output = repmat(simf_doi_10_1186_1754_1611_3_4_figure3_3(x) - 0.193,1,5);
end

function output = simf_doi_10_1186_1754_1611_3_4_figure3_4(x)
	ALPHA_BBa_B0032 = x(38);
	alpha_BBa_J23151 = x(43);
	scale_GFPmut3b_doi_10_1186_1754_1611_3_4 = x(40);
	BBa_J23151 = 10*alpha_BBa_J23151;
	GFPmut3b = scale_GFPmut3b_doi_10_1186_1754_1611_3_4*ALPHA_BBa_B0032*BBa_J23151;
	output = GFPmut3b;
end
function output = simf_doi_10_1186_1754_1611_3_4_figure3_4_error(x)
	output = repmat(simf_doi_10_1186_1754_1611_3_4_figure3_4(x) - 0.4,1,5);
end

function output = simf_doi_10_1186_1754_1611_3_4_figure3_5(x)
	ALPHA_BBa_B0032 = x(38);
	alpha_BBa_J23102 = x(44);
	scale_GFPmut3b_doi_10_1186_1754_1611_3_4 = x(40);
	BBa_J23102 = 10*alpha_BBa_J23102;
	GFPmut3b = scale_GFPmut3b_doi_10_1186_1754_1611_3_4*ALPHA_BBa_B0032*BBa_J23102;
	output = GFPmut3b;
end
function output = simf_doi_10_1186_1754_1611_3_4_figure3_5_error(x)
	output = repmat(simf_doi_10_1186_1754_1611_3_4_figure3_5(x) - 0.655,1,5);
end

function output = simf_doi_10_1186_1754_1611_3_4_figure3_6(x)
	ALPHA_BBa_B0032 = x(38);
	alpha_pLtetO_1 = x(18);
	scale_GFPmut3b_doi_10_1186_1754_1611_3_4 = x(40);
	pLtetO_1 = 10*alpha_pLtetO_1;
	GFPmut3b = scale_GFPmut3b_doi_10_1186_1754_1611_3_4*ALPHA_BBa_B0032*pLtetO_1;
	output = GFPmut3b;
end
function output = simf_doi_10_1186_1754_1611_3_4_figure3_6_error(x)
	output = repmat(simf_doi_10_1186_1754_1611_3_4_figure3_6(x) - 0.897,1,5);
end

function output = simf_doi_10_1186_1754_1611_3_4_figure3_7(x)
	ALPHA_BBa_B0032 = x(38);
	alpha_pLlacO_1 = x(12);
	scale_GFPmut3b_doi_10_1186_1754_1611_3_4 = x(40);
	pLlacO_1 = 10*alpha_pLlacO_1;
	GFPmut3b = scale_GFPmut3b_doi_10_1186_1754_1611_3_4*ALPHA_BBa_B0032*pLlacO_1;
	output = GFPmut3b;
end
function output = simf_doi_10_1186_1754_1611_3_4_figure3_7_error(x)
	output = repmat(simf_doi_10_1186_1754_1611_3_4_figure3_7(x) - 1,1,5);
end

function output = simf_doi_10_1093_nar_gku593_figure3D_1(x)
	ALPHA_BBa_B0030 = x(11);
	alpha_BBa_J23109 = x(45);
	scale_GFPmut3b_doi_10_1093_nar_gku593 = x(46);
	BBa_J23109 = 10*alpha_BBa_J23109;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku593*ALPHA_BBa_B0030*BBa_J23109;
	output = GFPmut3b;
end
function output = simf_doi_10_1093_nar_gku593_figure3D_1_error(x)
	output = repmat(simf_doi_10_1093_nar_gku593_figure3D_1(x) - 0.002,1,5);
end

function output = simf_doi_10_1093_nar_gku593_figure3D_2(x)
	ALPHA_BBa_B0030 = x(11);
	alpha_BBa_J23114 = x(47);
	scale_GFPmut3b_doi_10_1093_nar_gku593 = x(46);
	BBa_J23114 = 10*alpha_BBa_J23114;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku593*ALPHA_BBa_B0030*BBa_J23114;
	output = GFPmut3b;
end
function output = simf_doi_10_1093_nar_gku593_figure3D_2_error(x)
	output = repmat(simf_doi_10_1093_nar_gku593_figure3D_2(x) - 0.018,1,5);
end

function output = simf_doi_10_1093_nar_gku593_figure3D_3(x)
	ALPHA_BBa_B0030 = x(11);
	alpha_BBa_J23116 = x(41);
	scale_GFPmut3b_doi_10_1093_nar_gku593 = x(46);
	BBa_J23116 = 10*alpha_BBa_J23116;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku593*ALPHA_BBa_B0030*BBa_J23116;
	output = GFPmut3b;
end
function output = simf_doi_10_1093_nar_gku593_figure3D_3_error(x)
	output = repmat(simf_doi_10_1093_nar_gku593_figure3D_3(x) - 0.036,1,5);
end

function output = simf_doi_10_1093_nar_gku593_figure3D_4(x)
	ALPHA_BBa_B0030 = x(11);
	alpha_BBa_J23115 = x(48);
	scale_GFPmut3b_doi_10_1093_nar_gku593 = x(46);
	BBa_J23115 = 10*alpha_BBa_J23115;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku593*ALPHA_BBa_B0030*BBa_J23115;
	output = GFPmut3b;
end
function output = simf_doi_10_1093_nar_gku593_figure3D_4_error(x)
	output = repmat(simf_doi_10_1093_nar_gku593_figure3D_4(x) - 0.05,1,5);
end

function output = simf_doi_10_1093_nar_gku593_figure3D_5(x)
	ALPHA_BBa_B0030 = x(11);
	alpha_BBa_J23105 = x(49);
	scale_GFPmut3b_doi_10_1093_nar_gku593 = x(46);
	BBa_J23105 = 10*alpha_BBa_J23105;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku593*ALPHA_BBa_B0030*BBa_J23105;
	output = GFPmut3b;
end
function output = simf_doi_10_1093_nar_gku593_figure3D_5_error(x)
	output = repmat(simf_doi_10_1093_nar_gku593_figure3D_5(x) - 0.086,1,5);
end

function output = simf_doi_10_1093_nar_gku593_figure3D_6(x)
	ALPHA_BBa_B0030 = x(11);
	alpha_BBa_J23106 = x(50);
	scale_GFPmut3b_doi_10_1093_nar_gku593 = x(46);
	BBa_J23106 = 10*alpha_BBa_J23106;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku593*ALPHA_BBa_B0030*BBa_J23106;
	output = GFPmut3b;
end
function output = simf_doi_10_1093_nar_gku593_figure3D_6_error(x)
	output = repmat(simf_doi_10_1093_nar_gku593_figure3D_6(x) - 0.214,1,5);
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
function output = simf_doi_10_1093_nar_gku593_figureS6_error(x)
	output = repmat(simf_doi_10_1093_nar_gku593_figureS6(x,[0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2]) - [0 0.002 0.041 0.479 0.949 0.997 1 1 1],1,5);
end

function output = simf_doi_10_1093_nar_gku1388_figure2A_1(x)
	ALPHA_BBa_B0030 = x(11);
	alpha_BBa_J23101 = x(39);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(51);
	BBa_J23101 = 10*alpha_BBa_J23101;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*BBa_J23101;
	output = GFPmut3b;
end
function output = simf_doi_10_1093_nar_gku1388_figure2A_1_error(x)
	output = repmat(simf_doi_10_1093_nar_gku1388_figure2A_1(x) - 1,1,5);
end

function output = simf_doi_10_1093_nar_gku1388_figure2A_2(x)
	ALPHA_BBa_B0030 = x(11);
	alpha_BBa_J23106 = x(50);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(51);
	BBa_J23106 = 10*alpha_BBa_J23106;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*BBa_J23106;
	output = GFPmut3b;
end
function output = simf_doi_10_1093_nar_gku1388_figure2A_2_error(x)
	output = repmat(simf_doi_10_1093_nar_gku1388_figure2A_2(x) - 0.215,1,5);
end

function output = simf_doi_10_1093_nar_gku1388_figure2A_3(x)
	ALPHA_BBa_B0030 = x(11);
	alpha_BBa_J23105 = x(49);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(51);
	BBa_J23105 = 10*alpha_BBa_J23105;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*BBa_J23105;
	output = GFPmut3b;
end
function output = simf_doi_10_1093_nar_gku1388_figure2A_3_error(x)
	output = repmat(simf_doi_10_1093_nar_gku1388_figure2A_3(x) - 0.095,1,5);
end

function output = simf_doi_10_1093_nar_gku1388_figure2A_4(x)
	ALPHA_BBa_B0030 = x(11);
	alpha_BBa_J23115 = x(48);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(51);
	BBa_J23115 = 10*alpha_BBa_J23115;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*BBa_J23115;
	output = GFPmut3b;
end
function output = simf_doi_10_1093_nar_gku1388_figure2A_4_error(x)
	output = repmat(simf_doi_10_1093_nar_gku1388_figure2A_4(x) - 0.046,1,5);
end

function output = simf_doi_10_1093_nar_gku1388_figure2A_5(x)
	ALPHA_BBa_B0030 = x(11);
	alpha_BBa_J23114 = x(47);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(51);
	BBa_J23114 = 10*alpha_BBa_J23114;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*BBa_J23114;
	output = GFPmut3b;
end
function output = simf_doi_10_1093_nar_gku1388_figure2A_5_error(x)
	output = repmat(simf_doi_10_1093_nar_gku1388_figure2A_5(x) - 0.018,1,5);
end

function output = simf_doi_10_1093_nar_gku1388_figure2A_6(x)
	ALPHA_BBa_B0030 = x(11);
	alpha_BBa_J23117 = x(52);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(51);
	BBa_J23117 = 10*alpha_BBa_J23117;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*BBa_J23117;
	output = GFPmut3b;
end
function output = simf_doi_10_1093_nar_gku1388_figure2A_6_error(x)
	output = repmat(simf_doi_10_1093_nar_gku1388_figure2A_6(x) - 0.003,1,5);
end

function output = simf_doi_10_1093_nar_gku1388_figure2B_1(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0030 = x(11);
	K_aTc = x(16);
	K_pTET_star_ = x(53);
	alpha_BBa_J23117 = x(52);
	alpha_pTET_star_ = x(54);
	beta_pTET_star_ = x(55);
	n_aTc = x(21);
	n_pTET_star_ = x(56);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(51);
	for i = 1:length(input)
		BBa_J23117 = 10*alpha_BBa_J23117;
		aTc = input(i);
		tetR = ALPHA_BBa_B0030*BBa_J23117*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pTET_star_ = 10*(beta_pTET_star_ + (alpha_pTET_star_ - beta_pTET_star_)*K_pTET_star_^n_pTET_star_/(K_pTET_star_^n_pTET_star_ + tetR^n_pTET_star_));
		GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*pTET_star_;
		output(i) = GFPmut3b;
	end
end
function output = simf_doi_10_1093_nar_gku1388_figure2B_1_error(x)
	output = repmat(simf_doi_10_1093_nar_gku1388_figure2B_1(x,[2 4 10 25 50 100 200]) - [0.002 0.014 0.167 0.666 0.817 0.841 0.831],1,5);
end

function output = simf_doi_10_1093_nar_gku1388_figure2B_2(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0030 = x(11);
	K_aTc = x(16);
	K_pTET_star_ = x(53);
	alpha_BBa_J23114 = x(47);
	alpha_pTET_star_ = x(54);
	beta_pTET_star_ = x(55);
	n_aTc = x(21);
	n_pTET_star_ = x(56);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(51);
	for i = 1:length(input)
		BBa_J23114 = 10*alpha_BBa_J23114;
		aTc = input(i);
		tetR = ALPHA_BBa_B0030*BBa_J23114*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pTET_star_ = 10*(beta_pTET_star_ + (alpha_pTET_star_ - beta_pTET_star_)*K_pTET_star_^n_pTET_star_/(K_pTET_star_^n_pTET_star_ + tetR^n_pTET_star_));
		GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*pTET_star_;
		output(i) = GFPmut3b;
	end
end
function output = simf_doi_10_1093_nar_gku1388_figure2B_2_error(x)
	output = repmat(simf_doi_10_1093_nar_gku1388_figure2B_2(x,[10 25 50 100 200]) - [0.001 0.022 0.176 0.411 0.461],1,5);
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
function output = simf_doi_10_1038_nature09565_figure1C_error(x)
	output = repmat(simf_doi_10_1038_nature09565_figure1C(x,[0 0.001 0.015 0.05 0.5 5 10]) - [0.002 0.002 0.057 0.628 0.999 1 1],1,5);
end

function output = simf_doi_10_1038_nature09565_figureS1C(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0032 = x(38);
	ALPHA_BBa_B0033 = x(57);
	K_aTc = x(16);
	K_pLtetO_1 = x(17);
	alpha_BBa_J23117 = x(52);
	alpha_pLtetO_1 = x(18);
	beta_pLtetO_1 = x(20);
	n_aTc = x(21);
	n_pLtetO_1 = x(22);
	scale_eYFP_doi_10_1038_nature09565 = x(59);
	for i = 1:length(input)
		BBa_J23117 = 15*alpha_BBa_J23117;
		aTc = input(i);
		tetR = ALPHA_BBa_B0032*BBa_J23117*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pLtetO_1 = 15*(beta_pLtetO_1 + (alpha_pLtetO_1 - beta_pLtetO_1)*K_pLtetO_1^n_pLtetO_1/(K_pLtetO_1^n_pLtetO_1 + tetR^n_pLtetO_1));
		eYFP = scale_eYFP_doi_10_1038_nature09565*ALPHA_BBa_B0033*pLtetO_1;
		output(i) = eYFP;
	end
end
function output = simf_doi_10_1038_nature09565_figureS1C_error(x)
	output = repmat(simf_doi_10_1038_nature09565_figureS1C(x,[0 0.025 0.25 2.5 25 100 250]) - [0.005 0.006 0.007 0.03 0.297 0.388 0.398],1,5);
end

function output = simf_doi_10_1038_nbt_2401_figureS11(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_A = x(60);
	K_IPTG = x(25);
	K_pTAC = x(61);
	alpha_BBa_J23101 = x(39);
	alpha_pTAC = x(62);
	beta_pTAC = x(63);
	n_IPTG = x(30);
	n_pTAC = x(64);
	scale_sfGFP_doi_10_1038_nbt_2401 = x(65);
	for i = 1:length(input)
		BBa_J23101 = 5*alpha_BBa_J23101;
		IPTG = input(i);
		lacI = ALPHA_RBS_A*BBa_J23101*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pTAC = 5*(beta_pTAC + (alpha_pTAC - beta_pTAC)*K_pTAC^n_pTAC/(K_pTAC^n_pTAC + lacI^n_pTAC));
		sfGFP = scale_sfGFP_doi_10_1038_nbt_2401*ALPHA_RBS_A*pTAC;
		output(i) = sfGFP;
	end
end
function output = simf_doi_10_1038_nbt_2401_figureS11_error(x)
	output = repmat(simf_doi_10_1038_nbt_2401_figureS11(x,[0 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1]) - [0.025 0.025 0.035 0.055 0.115 0.2 0.3 0.4 0.65 1],1,5);
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
function output = simf_doi_10_1038_nbt_2401_figureS13_error(x)
	output = repmat(simf_doi_10_1038_nbt_2401_figureS13(x,[0.1 1 2 5 7 10 12.5 25 37.5 50 62.5]) - [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425],1,5);
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
function output = simf_doi_10_1038_nature11516_figureS8A_error(x)
	output = repmat(simf_doi_10_1038_nature11516_figureS8A(x,[0 0 0.002 0.01 0.04 0.16 1 4 25]) - [0.003 0.003 0.003 0.017 0.23 0.82 0.939 0.941 0.941],1,5);
end

function output = simf_doi_10_1038_nature11516_figureS8B(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_placI = x(24);
	ALPHA_RBS_psicA = x(66);
	K_IPTG = x(25);
	K_pTAC = x(61);
	alpha_pTAC = x(62);
	alpha_placIq = x(28);
	beta_pTAC = x(63);
	n_IPTG = x(30);
	n_pTAC = x(64);
	scale_mRFP1_doi_10_1038_nature11516 = x(67);
	for i = 1:length(input)
		placIq = 15*alpha_placIq;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placIq*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pTAC = 15*(beta_pTAC + (alpha_pTAC - beta_pTAC)*K_pTAC^n_pTAC/(K_pTAC^n_pTAC + lacI^n_pTAC));
		mRFP1 = scale_mRFP1_doi_10_1038_nature11516*ALPHA_RBS_psicA*pTAC;
		output(i) = mRFP1;
	end
end
function output = simf_doi_10_1038_nature11516_figureS8B_error(x)
	output = repmat(simf_doi_10_1038_nature11516_figureS8B(x,[0 0 0 0.001 0.004 0.025 0.1 0.5 2.5]) - [0.027 0.027 0.027 0.03 0.067 0.46 0.884 0.991 1],1,5);
end

function output = simf_doi_10_1038_nbt_2149_Supplement_page17(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pET_29b = x(68);
	ALPHA_RBS_placI = x(24);
	K_IPTG = x(25);
	K_placUV5 = x(69);
	alpha_placI = x(36);
	alpha_placUV5 = x(70);
	beta_placUV5 = x(71);
	n_IPTG = x(30);
	n_placUV5 = x(72);
	scale_mRFP1_doi_10_1038_nbt_2149 = x(73);
	for i = 1:length(input)
		placI = 10*alpha_placI;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placI*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		placUV5 = 10*(beta_placUV5 + (alpha_placUV5 - beta_placUV5)*K_placUV5^n_placUV5/(K_placUV5^n_placUV5 + lacI^n_placUV5));
		mRFP1 = scale_mRFP1_doi_10_1038_nbt_2149*ALPHA_RBS_pET_29b*placUV5;
		output(i) = mRFP1;
	end
end
function output = simf_doi_10_1038_nbt_2149_Supplement_page17_error(x)
	output = repmat(simf_doi_10_1038_nbt_2149_Supplement_page17(x,[0 0.001 0.005 0.01 0.025 0.05 0.1 0.25 0.5]) - [0.02 0.031 0.092 0.306 0.745 0.898 0.969 1 1],1,5);
end

function x_opt = simultaneous_fitting(ensemble_size, ite_num, func_ite_num)
	parameter_name = {'ALPHA_RBS_pBAD'; 'ALPHA_RBS_pC'; 'K_arabinose'; 'K_pBAD'; 'alpha_pBAD'; 'alpha_pC'; 'beta_pBAD'; 'n_arabinose'; 'n_pBAD'; 'scale_GFPuv_doi_10_1128_'; 'ALPHA_BBa_B0030'; 'alpha_pLlacO_1'; 'scale_mCherry_doi_10_103'; 'ALPHA_RBSII'; 'ALPHA_RBS_pN25'; 'K_aTc'; 'K_pLtetO_1'; 'alpha_pLtetO_1'; 'alpha_pN25'; 'beta_pLtetO_1'; 'n_aTc'; 'n_pLtetO_1'; 'scale_luc_doi__10_1093_n'; 'ALPHA_RBS_placI'; 'K_IPTG'; 'K_pLAC'; 'alpha_pLAC'; 'alpha_placIq'; 'beta_pLAC'; 'n_IPTG'; 'n_pLAC'; 'scale_eYFP_http___dspace'; 'ALPHA_RBS_Unknown_Hoosha'; 'scale_eYFP_doi_10_1073pn'; 'ALPHA_RBS_pLAC'; 'alpha_placI'; 'scale_lacZ_doi_10_1073_p'; 'ALPHA_BBa_B0032'; 'alpha_BBa_J23101'; 'scale_GFPmut3b_doi_10_11'; 'alpha_BBa_J23116'; 'alpha_BBa_J23150'; 'alpha_BBa_J23151'; 'alpha_BBa_J23102'; 'alpha_BBa_J23109'; 'scale_GFPmut3b_doi_10_10'; 'alpha_BBa_J23114'; 'alpha_BBa_J23115'; 'alpha_BBa_J23105'; 'alpha_BBa_J23106'; 'scale_GFPmut3b_doi_10_10'; 'alpha_BBa_J23117'; 'K_pTET_star_'; 'alpha_pTET_star_'; 'beta_pTET_star_'; 'n_pTET_star_'; 'ALPHA_BBa_B0033'; 'ALPHA_BBa_B0034'; 'scale_eYFP_doi_10_1038_n'; 'ALPHA_RBS_A'; 'K_pTAC'; 'alpha_pTAC'; 'beta_pTAC'; 'n_pTAC'; 'scale_sfGFP_doi_10_1038_'; 'ALPHA_RBS_psicA'; 'scale_mRFP1_doi_10_1038_'; 'ALPHA_RBS_pET_29b'; 'K_placUV5'; 'alpha_placUV5'; 'beta_placUV5'; 'n_placUV5'; 'scale_mRFP1_doi_10_1038_'};
	lb = [0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01 0.03 0.002 0.01 0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01 0.03 0.01 0.03 0.002 0.01 0.03 0.002 0.01 0.002 0.002 0.002 0.002 0.002 0.01 0.002 0.002 0.002 0.002 0.01 0.002 0.001 0.002 0 1 0.03 0.03 0.01 0.03 0.001 0.002 0 1 0.01 0.03 0.01 0.03 0.001 0.002 0 1 0.01];
	ub = [5 5 50 50 50 50 1 4 4 10 5 50 10 5 5 50 50 50 50 1 4 4 10 5 50 50 50 50 1 4 4 10 5 10 5 50 10 5 50 10 50 50 50 50 50 10 50 50 50 50 10 50 50 50 1 4 5 5 10 5 50 50 1 4 10 5 10 5 50 50 1 4 10];
	[CI_lb, CI_ub, x_opt, solution_ensemble] = fit_a_model(ensemble_size, ite_num, func_ite_num, @simultaneous_fitting_error, lb, ub);
	parameter_fileID = fopen('Plot_data/sim_CI_Hill.dat','w');
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
	SA_fileID = fopen('Plot_data/sensitivity_analysis_Hill.dat','w');
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
	sim_pair_fileID = fopen('Plot_data/sim_pair_Hill.dat','w');
	json_fileID = fopen('JSON/simulation_data_Hill.json','w');
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
	print_simulated_data_to_json(json_fileID, 'doi: 10.1093/nar/25.6.1203', 'figure4a', aTc, output_sim, {'pN25 = 1*{alpha_pN25}'; 'tetR = {ALPHA_RBS_pN25}*pN25*{K_aTc}^{n_aTc}/({K_aTc}^{n_aTc} + aTc^{n_aTc})'; 'pLtetO_1 = 15*({beta_pLtetO_1} + ({alpha_pLtetO_1} - {beta_pLtetO_1})*{K_pLtetO_1}^{n_pLtetO_1}/({K_pLtetO_1}^{n_pLtetO_1} + tetR^{n_pLtetO_1}))'; 'luc = {scale_luc_doi__10_1093_nar_25_6_1203}*{ALPHA_RBSII}*pLtetO_1'});
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
	print_simulated_data_to_json(json_fileID, 'doi:10.1093/nar/gku1388', 'figure2B-1', aTc, output_sim, {'BBa_J23117 = 10*{alpha_BBa_J23117}'; 'tetR = {ALPHA_BBa_B0030}*BBa_J23117*{K_aTc}^{n_aTc}/({K_aTc}^{n_aTc} + aTc^{n_aTc})'; 'pTET_star_ = 10*({beta_pTET_star_} + ({alpha_pTET_star_} - {beta_pTET_star_})*{K_pTET_star_}^{n_pTET_star_}/({K_pTET_star_}^{n_pTET_star_} + tetR^{n_pTET_star_}))'; 'GFPmut3b = {scale_GFPmut3b_doi_10_1093_nar_gku1388}*{ALPHA_BBa_B0030}*pTET_star_'});
	%%%%%%%%%%%%%
	aTc = [10 25 50 100 200];
	GFPmut3b = [0.001 0.022 0.176 0.411 0.461];
	output_sim = simf_doi_10_1093_nar_gku1388_figure2B_2(sim_x_opt,aTc);
	sim_pair_matrix = [output_sim; GFPmut3b];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1093/nar/gku1388', 'figure2B-2', aTc, output_sim, {'BBa_J23114 = 10*{alpha_BBa_J23114}'; 'tetR = {ALPHA_BBa_B0030}*BBa_J23114*{K_aTc}^{n_aTc}/({K_aTc}^{n_aTc} + aTc^{n_aTc})'; 'pTET_star_ = 10*({beta_pTET_star_} + ({alpha_pTET_star_} - {beta_pTET_star_})*{K_pTET_star_}^{n_pTET_star_}/({K_pTET_star_}^{n_pTET_star_} + tetR^{n_pTET_star_}))'; 'GFPmut3b = {scale_GFPmut3b_doi_10_1093_nar_gku1388}*{ALPHA_BBa_B0030}*pTET_star_'});
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
	print_simulated_data_to_json(json_fileID, 'doi:10.1038/nature09565', 'figureS1C', aTc, output_sim, {'BBa_J23117 = 15*{alpha_BBa_J23117}'; 'tetR = {ALPHA_BBa_B0032}*BBa_J23117*{K_aTc}^{n_aTc}/({K_aTc}^{n_aTc} + aTc^{n_aTc})'; 'pLtetO_1 = 15*({beta_pLtetO_1} + ({alpha_pLtetO_1} - {beta_pLtetO_1})*{K_pLtetO_1}^{n_pLtetO_1}/({K_pLtetO_1}^{n_pLtetO_1} + tetR^{n_pLtetO_1}))'; 'eYFP = {scale_eYFP_doi_10_1038_nature09565}*{ALPHA_BBa_B0033}*pLtetO_1'});
	%%%%%%%%%%%%%
	IPTG = [0 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1];
	sfGFP = [0.025 0.025 0.035 0.055 0.115 0.2 0.3 0.4 0.65 1];
	output_sim = simf_doi_10_1038_nbt_2401_figureS11(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; sfGFP];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1038/nbt.2401', 'figureS11', IPTG, output_sim, {'BBa_J23101 = 5*{alpha_BBa_J23101}'; 'lacI = {ALPHA_RBS_A}*BBa_J23101*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'pTAC = 5*({beta_pTAC} + ({alpha_pTAC} - {beta_pTAC})*{K_pTAC}^{n_pTAC}/({K_pTAC}^{n_pTAC} + lacI^{n_pTAC}))'; 'sfGFP = {scale_sfGFP_doi_10_1038_nbt_2401}*{ALPHA_RBS_A}*pTAC'});
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
	print_simulated_data_to_json(json_fileID, 'doi:10.1038/nature11516', 'figureS8B', IPTG, output_sim, {'placIq = 15*{alpha_placIq}'; 'lacI = {ALPHA_RBS_placI}*placIq*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'pTAC = 15*({beta_pTAC} + ({alpha_pTAC} - {beta_pTAC})*{K_pTAC}^{n_pTAC}/({K_pTAC}^{n_pTAC} + lacI^{n_pTAC}))'; 'mRFP1 = {scale_mRFP1_doi_10_1038_nature11516}*{ALPHA_RBS_psicA}*pTAC'});
	%%%%%%%%%%%%%
	IPTG = [0 0.001 0.005 0.01 0.025 0.05 0.1 0.25 0.5];
	mRFP1 = [0.02 0.031 0.092 0.306 0.745 0.898 0.969 1 1];
	output_sim = simf_doi_10_1038_nbt_2149_Supplement_page17(sim_x_opt,IPTG);
	sim_pair_matrix = [output_sim; mRFP1];
	fprintf(sim_pair_fileID, '%f \t %f\n', sim_pair_matrix);
	print_simulated_data_to_json(json_fileID, 'doi:10.1038/nbt.2149', 'Supplement-page17', IPTG, output_sim, {'placI = 10*{alpha_placI}'; 'lacI = {ALPHA_RBS_placI}*placI*{K_IPTG}^{n_IPTG}/({K_IPTG}^{n_IPTG} + IPTG^{n_IPTG})'; 'placUV5 = 10*({beta_placUV5} + ({alpha_placUV5} - {beta_placUV5})*{K_placUV5}^{n_placUV5}/({K_placUV5}^{n_placUV5} + lacI^{n_placUV5}))'; 'mRFP1 = {scale_mRFP1_doi_10_1038_nbt_2149}*{ALPHA_RBS_pET_29b}*placUV5'});
	fclose(sim_pair_fileID);
	fclose(json_fileID);
end

function sim_LOOCV(ensemble_size, ite_num, func_ite_num)
	lb = [0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01 0.03 0.002 0.01 0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01 0.03 0.01 0.03 0.002 0.01 0.03 0.002 0.01 0.002 0.002 0.002 0.002 0.002 0.01 0.002 0.002 0.002 0.002 0.01 0.002 0.001 0.002 0 1 0.03 0.03 0.01 0.03 0.001 0.002 0 1 0.01 0.03 0.01 0.03 0.001 0.002 0 1 0.01];
	ub = [5 5 50 50 50 50 1 4 4 10 5 50 10 5 5 50 50 50 50 1 4 4 10 5 50 50 50 50 1 4 4 10 5 10 5 50 10 5 50 10 50 50 50 50 50 10 50 50 50 50 10 50 50 50 1 4 5 5 10 5 50 50 1 4 10 5 10 5 50 50 1 4 10];
	LOOCV_sim_data_fileID = fopen('Plot_data/LOOCV_sim_data_Hill.dat','w');
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

function parameter_opt = sequential_fitting(N1, N2)
	error_opt = 1e6;
	parameter_opt = zeros(1,73);
	for i = 1:N1
		model_order = randperm(35);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			parameter_opt = parameter;
		end
	end
	parameter_name_list = {'ALPHA_RBS_pBAD'; 'ALPHA_RBS_pC'; 'K_arabinose'; 'K_pBAD'; 'alpha_pBAD'; 'alpha_pC'; 'beta_pBAD'; 'n_arabinose'; 'n_pBAD'; 'scale_GFPuv_doi_10_1128_AEM_00791_07'; 'ALPHA_BBa_B0030'; 'alpha_pLlacO_1'; 'scale_mCherry_doi_10_1038_nature12148'; 'ALPHA_RBSII'; 'ALPHA_RBS_pN25'; 'K_aTc'; 'K_pLtetO_1'; 'alpha_pLtetO_1'; 'alpha_pN25'; 'beta_pLtetO_1'; 'n_aTc'; 'n_pLtetO_1'; 'scale_luc_doi__10_1093_nar_25_6_1203'; 'ALPHA_RBS_placI'; 'K_IPTG'; 'K_pLAC'; 'alpha_pLAC'; 'alpha_placIq'; 'beta_pLAC'; 'n_IPTG'; 'n_pLAC'; 'scale_eYFP_http___dspace_mit_edu_handle_1721_1_8228'; 'ALPHA_RBS_Unknown_Hooshangi'; 'scale_eYFP_doi_10_1073pnas_0408507102'; 'ALPHA_RBS_pLAC'; 'alpha_placI'; 'scale_lacZ_doi_10_1073_pnas_0606717104'; 'ALPHA_BBa_B0032'; 'alpha_BBa_J23101'; 'scale_GFPmut3b_doi_10_1186_1754_1611_3_4'; 'alpha_BBa_J23116'; 'alpha_BBa_J23150'; 'alpha_BBa_J23151'; 'alpha_BBa_J23102'; 'alpha_BBa_J23109'; 'scale_GFPmut3b_doi_10_1093_nar_gku593'; 'alpha_BBa_J23114'; 'alpha_BBa_J23115'; 'alpha_BBa_J23105'; 'alpha_BBa_J23106'; 'scale_GFPmut3b_doi_10_1093_nar_gku1388'; 'alpha_BBa_J23117'; 'K_pTET_star_'; 'alpha_pTET_star_'; 'beta_pTET_star_'; 'n_pTET_star_'; 'ALPHA_BBa_B0033'; 'ALPHA_BBa_B0034'; 'scale_eYFP_doi_10_1038_nature09565'; 'ALPHA_RBS_A'; 'K_pTAC'; 'alpha_pTAC'; 'beta_pTAC'; 'n_pTAC'; 'scale_sfGFP_doi_10_1038_nbt_2401'; 'ALPHA_RBS_psicA'; 'scale_mRFP1_doi_10_1038_nature11516'; 'ALPHA_RBS_pET_29b'; 'K_placUV5'; 'alpha_placUV5'; 'beta_placUV5'; 'n_placUV5'; 'scale_mRFP1_doi_10_1038_nbt_2149'};
	parameter_fileID = fopen('Plot_data/sequential_parameter.dat','w');
	for i = 1:length(parameter_name_list)
		fprintf(parameter_fileID, '%s\t', parameter_name_list{i});
		fprintf(parameter_fileID, '%f\n', parameter_opt(i));
	end
	fclose(parameter_fileID);
end
function best_parameter = sequential_fitting_for_one_order(model_order, N)
	options = optimset('MaxIter', 1000, 'MaxFunEvals', 30000);
	best_parameter = zeros(1,73);
	for model_id = 1:length(model_order)
		switch (model_id)
			case model_order(model_id)
				lb = [0.03 0.002 0.002 0 1 0.001 0.03 0.01 1 0.001];
				ub = [5 50 50 1 4 50 5 10 4 50];
				if (best_parameter(2) > 1e-10)
					lb(1) = best_parameter(2) - 1e-10;
					ub(1) = best_parameter(2) + 1e-10;
				end
				if (best_parameter(6) > 1e-10)
					lb(2) = best_parameter(6) - 1e-10;
					ub(2) = best_parameter(6) + 1e-10;
				end
				if (best_parameter(5) > 1e-10)
					lb(3) = best_parameter(5) - 1e-10;
					ub(3) = best_parameter(5) + 1e-10;
				end
				if (best_parameter(7) > 1e-10)
					lb(4) = best_parameter(7) - 1e-10;
					ub(4) = best_parameter(7) + 1e-10;
				end
				if (best_parameter(9) > 1e-10)
					lb(5) = best_parameter(9) - 1e-10;
					ub(5) = best_parameter(9) + 1e-10;
				end
				if (best_parameter(4) > 1e-10)
					lb(6) = best_parameter(4) - 1e-10;
					ub(6) = best_parameter(4) + 1e-10;
				end
				if (best_parameter(1) > 1e-10)
					lb(7) = best_parameter(1) - 1e-10;
					ub(7) = best_parameter(1) + 1e-10;
				end
				if (best_parameter(10) > 1e-10)
					lb(8) = best_parameter(10) - 1e-10;
					ub(8) = best_parameter(10) + 1e-10;
				end
				if (best_parameter(8) > 1e-10)
					lb(9) = best_parameter(8) - 1e-10;
					ub(9) = best_parameter(8) + 1e-10;
				end
				if (best_parameter(3) > 1e-10)
					lb(10) = best_parameter(3) - 1e-10;
					ub(10) = best_parameter(3) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@simf_doi_10_1128_AEM_00791_07_figure3a_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(2) = x_opt(1);
				best_parameter(6) = x_opt(2);
				best_parameter(5) = x_opt(3);
				best_parameter(7) = x_opt(4);
				best_parameter(9) = x_opt(5);
				best_parameter(4) = x_opt(6);
				best_parameter(1) = x_opt(7);
				best_parameter(10) = x_opt(8);
				best_parameter(8) = x_opt(9);
				best_parameter(3) = x_opt(10);
			case model_order(model_id)
				lb = [0.002 0 1 0.001 0.03 0.01 0.002 1 0.001];
				ub = [50 1 4 50 5 10 50 4 50];
				if (best_parameter(5) > 1e-10)
					lb(1) = best_parameter(5) - 1e-10;
					ub(1) = best_parameter(5) + 1e-10;
				end
				if (best_parameter(7) > 1e-10)
					lb(2) = best_parameter(7) - 1e-10;
					ub(2) = best_parameter(7) + 1e-10;
				end
				if (best_parameter(9) > 1e-10)
					lb(3) = best_parameter(9) - 1e-10;
					ub(3) = best_parameter(9) + 1e-10;
				end
				if (best_parameter(4) > 1e-10)
					lb(4) = best_parameter(4) - 1e-10;
					ub(4) = best_parameter(4) + 1e-10;
				end
				if (best_parameter(11) > 1e-10)
					lb(5) = best_parameter(11) - 1e-10;
					ub(5) = best_parameter(11) + 1e-10;
				end
				if (best_parameter(13) > 1e-10)
					lb(6) = best_parameter(13) - 1e-10;
					ub(6) = best_parameter(13) + 1e-10;
				end
				if (best_parameter(12) > 1e-10)
					lb(7) = best_parameter(12) - 1e-10;
					ub(7) = best_parameter(12) + 1e-10;
				end
				if (best_parameter(8) > 1e-10)
					lb(8) = best_parameter(8) - 1e-10;
					ub(8) = best_parameter(8) + 1e-10;
				end
				if (best_parameter(3) > 1e-10)
					lb(9) = best_parameter(3) - 1e-10;
					ub(9) = best_parameter(3) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@simf_doi_10_1038_nature12148_figure16a_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(5) = x_opt(1);
				best_parameter(7) = x_opt(2);
				best_parameter(9) = x_opt(3);
				best_parameter(4) = x_opt(4);
				best_parameter(11) = x_opt(5);
				best_parameter(13) = x_opt(6);
				best_parameter(12) = x_opt(7);
				best_parameter(8) = x_opt(8);
				best_parameter(3) = x_opt(9);
			case model_order(model_id)
				lb = [0.002 0 1 0.001 0.03 0.01 0.002 0.03 1 0.001];
				ub = [50 1 4 50 5 10 50 5 4 50];
				if (best_parameter(18) > 1e-10)
					lb(1) = best_parameter(18) - 1e-10;
					ub(1) = best_parameter(18) + 1e-10;
				end
				if (best_parameter(20) > 1e-10)
					lb(2) = best_parameter(20) - 1e-10;
					ub(2) = best_parameter(20) + 1e-10;
				end
				if (best_parameter(22) > 1e-10)
					lb(3) = best_parameter(22) - 1e-10;
					ub(3) = best_parameter(22) + 1e-10;
				end
				if (best_parameter(17) > 1e-10)
					lb(4) = best_parameter(17) - 1e-10;
					ub(4) = best_parameter(17) + 1e-10;
				end
				if (best_parameter(14) > 1e-10)
					lb(5) = best_parameter(14) - 1e-10;
					ub(5) = best_parameter(14) + 1e-10;
				end
				if (best_parameter(23) > 1e-10)
					lb(6) = best_parameter(23) - 1e-10;
					ub(6) = best_parameter(23) + 1e-10;
				end
				if (best_parameter(19) > 1e-10)
					lb(7) = best_parameter(19) - 1e-10;
					ub(7) = best_parameter(19) + 1e-10;
				end
				if (best_parameter(15) > 1e-10)
					lb(8) = best_parameter(15) - 1e-10;
					ub(8) = best_parameter(15) + 1e-10;
				end
				if (best_parameter(21) > 1e-10)
					lb(9) = best_parameter(21) - 1e-10;
					ub(9) = best_parameter(21) + 1e-10;
				end
				if (best_parameter(16) > 1e-10)
					lb(10) = best_parameter(16) - 1e-10;
					ub(10) = best_parameter(16) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@simf_doi__10_1093_nar_25_6_1203_figure4a_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(18) = x_opt(1);
				best_parameter(20) = x_opt(2);
				best_parameter(22) = x_opt(3);
				best_parameter(17) = x_opt(4);
				best_parameter(14) = x_opt(5);
				best_parameter(23) = x_opt(6);
				best_parameter(19) = x_opt(7);
				best_parameter(15) = x_opt(8);
				best_parameter(21) = x_opt(9);
				best_parameter(16) = x_opt(10);
			case model_order(model_id)
				lb = [0.03 0.002 0.002 0 1 0.001 0.03 0.01 1 0.001];
				ub = [5 50 50 1 4 50 5 10 4 50];
				if (best_parameter(24) > 1e-10)
					lb(1) = best_parameter(24) - 1e-10;
					ub(1) = best_parameter(24) + 1e-10;
				end
				if (best_parameter(28) > 1e-10)
					lb(2) = best_parameter(28) - 1e-10;
					ub(2) = best_parameter(28) + 1e-10;
				end
				if (best_parameter(27) > 1e-10)
					lb(3) = best_parameter(27) - 1e-10;
					ub(3) = best_parameter(27) + 1e-10;
				end
				if (best_parameter(29) > 1e-10)
					lb(4) = best_parameter(29) - 1e-10;
					ub(4) = best_parameter(29) + 1e-10;
				end
				if (best_parameter(31) > 1e-10)
					lb(5) = best_parameter(31) - 1e-10;
					ub(5) = best_parameter(31) + 1e-10;
				end
				if (best_parameter(26) > 1e-10)
					lb(6) = best_parameter(26) - 1e-10;
					ub(6) = best_parameter(26) + 1e-10;
				end
				if (best_parameter(14) > 1e-10)
					lb(7) = best_parameter(14) - 1e-10;
					ub(7) = best_parameter(14) + 1e-10;
				end
				if (best_parameter(32) > 1e-10)
					lb(8) = best_parameter(32) - 1e-10;
					ub(8) = best_parameter(32) + 1e-10;
				end
				if (best_parameter(30) > 1e-10)
					lb(9) = best_parameter(30) - 1e-10;
					ub(9) = best_parameter(30) + 1e-10;
				end
				if (best_parameter(25) > 1e-10)
					lb(10) = best_parameter(25) - 1e-10;
					ub(10) = best_parameter(25) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(24) = x_opt(1);
				best_parameter(28) = x_opt(2);
				best_parameter(27) = x_opt(3);
				best_parameter(29) = x_opt(4);
				best_parameter(31) = x_opt(5);
				best_parameter(26) = x_opt(6);
				best_parameter(14) = x_opt(7);
				best_parameter(32) = x_opt(8);
				best_parameter(30) = x_opt(9);
				best_parameter(25) = x_opt(10);
			case model_order(model_id)
				lb = [0.002 0.03 0.002 0 1 0.001 0.03 0.01 1 0.001];
				ub = [50 5 50 1 4 50 5 10 4 50];
				if (best_parameter(28) > 1e-10)
					lb(1) = best_parameter(28) - 1e-10;
					ub(1) = best_parameter(28) + 1e-10;
				end
				if (best_parameter(24) > 1e-10)
					lb(2) = best_parameter(24) - 1e-10;
					ub(2) = best_parameter(24) + 1e-10;
				end
				if (best_parameter(18) > 1e-10)
					lb(3) = best_parameter(18) - 1e-10;
					ub(3) = best_parameter(18) + 1e-10;
				end
				if (best_parameter(20) > 1e-10)
					lb(4) = best_parameter(20) - 1e-10;
					ub(4) = best_parameter(20) + 1e-10;
				end
				if (best_parameter(22) > 1e-10)
					lb(5) = best_parameter(22) - 1e-10;
					ub(5) = best_parameter(22) + 1e-10;
				end
				if (best_parameter(17) > 1e-10)
					lb(6) = best_parameter(17) - 1e-10;
					ub(6) = best_parameter(17) + 1e-10;
				end
				if (best_parameter(33) > 1e-10)
					lb(7) = best_parameter(33) - 1e-10;
					ub(7) = best_parameter(33) + 1e-10;
				end
				if (best_parameter(34) > 1e-10)
					lb(8) = best_parameter(34) - 1e-10;
					ub(8) = best_parameter(34) + 1e-10;
				end
				if (best_parameter(21) > 1e-10)
					lb(9) = best_parameter(21) - 1e-10;
					ub(9) = best_parameter(21) + 1e-10;
				end
				if (best_parameter(16) > 1e-10)
					lb(10) = best_parameter(16) - 1e-10;
					ub(10) = best_parameter(16) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@simf_doi_10_1073pnas_0408507102_figure2A_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(28) = x_opt(1);
				best_parameter(24) = x_opt(2);
				best_parameter(18) = x_opt(3);
				best_parameter(20) = x_opt(4);
				best_parameter(22) = x_opt(5);
				best_parameter(17) = x_opt(6);
				best_parameter(33) = x_opt(7);
				best_parameter(34) = x_opt(8);
				best_parameter(21) = x_opt(9);
				best_parameter(16) = x_opt(10);
			case model_order(model_id)
				lb = [0.002 0.03 0.002 0 1 0.001 0.03 0.01 1 0.001];
				ub = [50 5 50 1 4 50 5 10 4 50];
				if (best_parameter(36) > 1e-10)
					lb(1) = best_parameter(36) - 1e-10;
					ub(1) = best_parameter(36) + 1e-10;
				end
				if (best_parameter(24) > 1e-10)
					lb(2) = best_parameter(24) - 1e-10;
					ub(2) = best_parameter(24) + 1e-10;
				end
				if (best_parameter(27) > 1e-10)
					lb(3) = best_parameter(27) - 1e-10;
					ub(3) = best_parameter(27) + 1e-10;
				end
				if (best_parameter(29) > 1e-10)
					lb(4) = best_parameter(29) - 1e-10;
					ub(4) = best_parameter(29) + 1e-10;
				end
				if (best_parameter(31) > 1e-10)
					lb(5) = best_parameter(31) - 1e-10;
					ub(5) = best_parameter(31) + 1e-10;
				end
				if (best_parameter(26) > 1e-10)
					lb(6) = best_parameter(26) - 1e-10;
					ub(6) = best_parameter(26) + 1e-10;
				end
				if (best_parameter(35) > 1e-10)
					lb(7) = best_parameter(35) - 1e-10;
					ub(7) = best_parameter(35) + 1e-10;
				end
				if (best_parameter(37) > 1e-10)
					lb(8) = best_parameter(37) - 1e-10;
					ub(8) = best_parameter(37) + 1e-10;
				end
				if (best_parameter(30) > 1e-10)
					lb(9) = best_parameter(30) - 1e-10;
					ub(9) = best_parameter(30) + 1e-10;
				end
				if (best_parameter(25) > 1e-10)
					lb(10) = best_parameter(25) - 1e-10;
					ub(10) = best_parameter(25) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@simf_doi_10_1073_pnas_0606717104_figure4_6_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(36) = x_opt(1);
				best_parameter(24) = x_opt(2);
				best_parameter(27) = x_opt(3);
				best_parameter(29) = x_opt(4);
				best_parameter(31) = x_opt(5);
				best_parameter(26) = x_opt(6);
				best_parameter(35) = x_opt(7);
				best_parameter(37) = x_opt(8);
				best_parameter(30) = x_opt(9);
				best_parameter(25) = x_opt(10);
			case model_order(model_id)
				lb = [0.002 0.03 0.01];
				ub = [50 5 10];
				if (best_parameter(39) > 1e-10)
					lb(1) = best_parameter(39) - 1e-10;
					ub(1) = best_parameter(39) + 1e-10;
				end
				if (best_parameter(38) > 1e-10)
					lb(2) = best_parameter(38) - 1e-10;
					ub(2) = best_parameter(38) + 1e-10;
				end
				if (best_parameter(40) > 1e-10)
					lb(3) = best_parameter(40) - 1e-10;
					ub(3) = best_parameter(40) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@simf_doi_10_1186_1754_1611_3_4_figure3_1_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(39) = x_opt(1);
				best_parameter(38) = x_opt(2);
				best_parameter(40) = x_opt(3);
			case model_order(model_id)
				lb = [0.002 0.03 0.01];
				ub = [50 5 10];
				if (best_parameter(41) > 1e-10)
					lb(1) = best_parameter(41) - 1e-10;
					ub(1) = best_parameter(41) + 1e-10;
				end
				if (best_parameter(38) > 1e-10)
					lb(2) = best_parameter(38) - 1e-10;
					ub(2) = best_parameter(38) + 1e-10;
				end
				if (best_parameter(40) > 1e-10)
					lb(3) = best_parameter(40) - 1e-10;
					ub(3) = best_parameter(40) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@simf_doi_10_1186_1754_1611_3_4_figure3_2_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(41) = x_opt(1);
				best_parameter(38) = x_opt(2);
				best_parameter(40) = x_opt(3);
			case model_order(model_id)
				lb = [0.002 0.03 0.01];
				ub = [50 5 10];
				if (best_parameter(42) > 1e-10)
					lb(1) = best_parameter(42) - 1e-10;
					ub(1) = best_parameter(42) + 1e-10;
				end
				if (best_parameter(38) > 1e-10)
					lb(2) = best_parameter(38) - 1e-10;
					ub(2) = best_parameter(38) + 1e-10;
				end
				if (best_parameter(40) > 1e-10)
					lb(3) = best_parameter(40) - 1e-10;
					ub(3) = best_parameter(40) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@simf_doi_10_1186_1754_1611_3_4_figure3_3_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(42) = x_opt(1);
				best_parameter(38) = x_opt(2);
				best_parameter(40) = x_opt(3);
			case model_order(model_id)
				lb = [0.002 0.03 0.01];
				ub = [50 5 10];
				if (best_parameter(43) > 1e-10)
					lb(1) = best_parameter(43) - 1e-10;
					ub(1) = best_parameter(43) + 1e-10;
				end
				if (best_parameter(38) > 1e-10)
					lb(2) = best_parameter(38) - 1e-10;
					ub(2) = best_parameter(38) + 1e-10;
				end
				if (best_parameter(40) > 1e-10)
					lb(3) = best_parameter(40) - 1e-10;
					ub(3) = best_parameter(40) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@simf_doi_10_1186_1754_1611_3_4_figure3_4_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(43) = x_opt(1);
				best_parameter(38) = x_opt(2);
				best_parameter(40) = x_opt(3);
			case model_order(model_id)
				lb = [0.002 0.03 0.01];
				ub = [50 5 10];
				if (best_parameter(44) > 1e-10)
					lb(1) = best_parameter(44) - 1e-10;
					ub(1) = best_parameter(44) + 1e-10;
				end
				if (best_parameter(38) > 1e-10)
					lb(2) = best_parameter(38) - 1e-10;
					ub(2) = best_parameter(38) + 1e-10;
				end
				if (best_parameter(40) > 1e-10)
					lb(3) = best_parameter(40) - 1e-10;
					ub(3) = best_parameter(40) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@simf_doi_10_1186_1754_1611_3_4_figure3_5_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(44) = x_opt(1);
				best_parameter(38) = x_opt(2);
				best_parameter(40) = x_opt(3);
			case model_order(model_id)
				lb = [0.002 0.03 0.01];
				ub = [50 5 10];
				if (best_parameter(18) > 1e-10)
					lb(1) = best_parameter(18) - 1e-10;
					ub(1) = best_parameter(18) + 1e-10;
				end
				if (best_parameter(38) > 1e-10)
					lb(2) = best_parameter(38) - 1e-10;
					ub(2) = best_parameter(38) + 1e-10;
				end
				if (best_parameter(40) > 1e-10)
					lb(3) = best_parameter(40) - 1e-10;
					ub(3) = best_parameter(40) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@simf_doi_10_1186_1754_1611_3_4_figure3_6_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(18) = x_opt(1);
				best_parameter(38) = x_opt(2);
				best_parameter(40) = x_opt(3);
			case model_order(model_id)
				lb = [0.002 0.03 0.01];
				ub = [50 5 10];
				if (best_parameter(12) > 1e-10)
					lb(1) = best_parameter(12) - 1e-10;
					ub(1) = best_parameter(12) + 1e-10;
				end
				if (best_parameter(38) > 1e-10)
					lb(2) = best_parameter(38) - 1e-10;
					ub(2) = best_parameter(38) + 1e-10;
				end
				if (best_parameter(40) > 1e-10)
					lb(3) = best_parameter(40) - 1e-10;
					ub(3) = best_parameter(40) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@simf_doi_10_1186_1754_1611_3_4_figure3_7_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(12) = x_opt(1);
				best_parameter(38) = x_opt(2);
				best_parameter(40) = x_opt(3);
			case model_order(model_id)
				lb = [0.002 0.03 0.01];
				ub = [50 5 10];
				if (best_parameter(45) > 1e-10)
					lb(1) = best_parameter(45) - 1e-10;
					ub(1) = best_parameter(45) + 1e-10;
				end
				if (best_parameter(11) > 1e-10)
					lb(2) = best_parameter(11) - 1e-10;
					ub(2) = best_parameter(11) + 1e-10;
				end
				if (best_parameter(46) > 1e-10)
					lb(3) = best_parameter(46) - 1e-10;
					ub(3) = best_parameter(46) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@simf_doi_10_1093_nar_gku593_figure3D_1_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(45) = x_opt(1);
				best_parameter(11) = x_opt(2);
				best_parameter(46) = x_opt(3);
			case model_order(model_id)
				lb = [0.002 0.03 0.01];
				ub = [50 5 10];
				if (best_parameter(47) > 1e-10)
					lb(1) = best_parameter(47) - 1e-10;
					ub(1) = best_parameter(47) + 1e-10;
				end
				if (best_parameter(11) > 1e-10)
					lb(2) = best_parameter(11) - 1e-10;
					ub(2) = best_parameter(11) + 1e-10;
				end
				if (best_parameter(46) > 1e-10)
					lb(3) = best_parameter(46) - 1e-10;
					ub(3) = best_parameter(46) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@simf_doi_10_1093_nar_gku593_figure3D_2_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(47) = x_opt(1);
				best_parameter(11) = x_opt(2);
				best_parameter(46) = x_opt(3);
			case model_order(model_id)
				lb = [0.002 0.03 0.01];
				ub = [50 5 10];
				if (best_parameter(41) > 1e-10)
					lb(1) = best_parameter(41) - 1e-10;
					ub(1) = best_parameter(41) + 1e-10;
				end
				if (best_parameter(11) > 1e-10)
					lb(2) = best_parameter(11) - 1e-10;
					ub(2) = best_parameter(11) + 1e-10;
				end
				if (best_parameter(46) > 1e-10)
					lb(3) = best_parameter(46) - 1e-10;
					ub(3) = best_parameter(46) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@simf_doi_10_1093_nar_gku593_figure3D_3_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(41) = x_opt(1);
				best_parameter(11) = x_opt(2);
				best_parameter(46) = x_opt(3);
			case model_order(model_id)
				lb = [0.002 0.03 0.01];
				ub = [50 5 10];
				if (best_parameter(48) > 1e-10)
					lb(1) = best_parameter(48) - 1e-10;
					ub(1) = best_parameter(48) + 1e-10;
				end
				if (best_parameter(11) > 1e-10)
					lb(2) = best_parameter(11) - 1e-10;
					ub(2) = best_parameter(11) + 1e-10;
				end
				if (best_parameter(46) > 1e-10)
					lb(3) = best_parameter(46) - 1e-10;
					ub(3) = best_parameter(46) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@simf_doi_10_1093_nar_gku593_figure3D_4_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(48) = x_opt(1);
				best_parameter(11) = x_opt(2);
				best_parameter(46) = x_opt(3);
			case model_order(model_id)
				lb = [0.002 0.03 0.01];
				ub = [50 5 10];
				if (best_parameter(49) > 1e-10)
					lb(1) = best_parameter(49) - 1e-10;
					ub(1) = best_parameter(49) + 1e-10;
				end
				if (best_parameter(11) > 1e-10)
					lb(2) = best_parameter(11) - 1e-10;
					ub(2) = best_parameter(11) + 1e-10;
				end
				if (best_parameter(46) > 1e-10)
					lb(3) = best_parameter(46) - 1e-10;
					ub(3) = best_parameter(46) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@simf_doi_10_1093_nar_gku593_figure3D_5_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(49) = x_opt(1);
				best_parameter(11) = x_opt(2);
				best_parameter(46) = x_opt(3);
			case model_order(model_id)
				lb = [0.002 0.03 0.01];
				ub = [50 5 10];
				if (best_parameter(50) > 1e-10)
					lb(1) = best_parameter(50) - 1e-10;
					ub(1) = best_parameter(50) + 1e-10;
				end
				if (best_parameter(11) > 1e-10)
					lb(2) = best_parameter(11) - 1e-10;
					ub(2) = best_parameter(11) + 1e-10;
				end
				if (best_parameter(46) > 1e-10)
					lb(3) = best_parameter(46) - 1e-10;
					ub(3) = best_parameter(46) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@simf_doi_10_1093_nar_gku593_figure3D_6_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(50) = x_opt(1);
				best_parameter(11) = x_opt(2);
				best_parameter(46) = x_opt(3);
			case model_order(model_id)
				lb = [0.03 0.002 0.002 0 1 0.001 0.03 0.01 1 0.001];
				ub = [5 50 50 1 4 50 5 10 4 50];
				if (best_parameter(2) > 1e-10)
					lb(1) = best_parameter(2) - 1e-10;
					ub(1) = best_parameter(2) + 1e-10;
				end
				if (best_parameter(6) > 1e-10)
					lb(2) = best_parameter(6) - 1e-10;
					ub(2) = best_parameter(6) + 1e-10;
				end
				if (best_parameter(5) > 1e-10)
					lb(3) = best_parameter(5) - 1e-10;
					ub(3) = best_parameter(5) + 1e-10;
				end
				if (best_parameter(7) > 1e-10)
					lb(4) = best_parameter(7) - 1e-10;
					ub(4) = best_parameter(7) + 1e-10;
				end
				if (best_parameter(9) > 1e-10)
					lb(5) = best_parameter(9) - 1e-10;
					ub(5) = best_parameter(9) + 1e-10;
				end
				if (best_parameter(4) > 1e-10)
					lb(6) = best_parameter(4) - 1e-10;
					ub(6) = best_parameter(4) + 1e-10;
				end
				if (best_parameter(11) > 1e-10)
					lb(7) = best_parameter(11) - 1e-10;
					ub(7) = best_parameter(11) + 1e-10;
				end
				if (best_parameter(46) > 1e-10)
					lb(8) = best_parameter(46) - 1e-10;
					ub(8) = best_parameter(46) + 1e-10;
				end
				if (best_parameter(8) > 1e-10)
					lb(9) = best_parameter(8) - 1e-10;
					ub(9) = best_parameter(8) + 1e-10;
				end
				if (best_parameter(3) > 1e-10)
					lb(10) = best_parameter(3) - 1e-10;
					ub(10) = best_parameter(3) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@simf_doi_10_1093_nar_gku593_figureS6_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(2) = x_opt(1);
				best_parameter(6) = x_opt(2);
				best_parameter(5) = x_opt(3);
				best_parameter(7) = x_opt(4);
				best_parameter(9) = x_opt(5);
				best_parameter(4) = x_opt(6);
				best_parameter(11) = x_opt(7);
				best_parameter(46) = x_opt(8);
				best_parameter(8) = x_opt(9);
				best_parameter(3) = x_opt(10);
			case model_order(model_id)
				lb = [0.002 0.03 0.01];
				ub = [50 5 10];
				if (best_parameter(39) > 1e-10)
					lb(1) = best_parameter(39) - 1e-10;
					ub(1) = best_parameter(39) + 1e-10;
				end
				if (best_parameter(11) > 1e-10)
					lb(2) = best_parameter(11) - 1e-10;
					ub(2) = best_parameter(11) + 1e-10;
				end
				if (best_parameter(51) > 1e-10)
					lb(3) = best_parameter(51) - 1e-10;
					ub(3) = best_parameter(51) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@simf_doi_10_1093_nar_gku1388_figure2A_1_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(39) = x_opt(1);
				best_parameter(11) = x_opt(2);
				best_parameter(51) = x_opt(3);
			case model_order(model_id)
				lb = [0.002 0.03 0.01];
				ub = [50 5 10];
				if (best_parameter(50) > 1e-10)
					lb(1) = best_parameter(50) - 1e-10;
					ub(1) = best_parameter(50) + 1e-10;
				end
				if (best_parameter(11) > 1e-10)
					lb(2) = best_parameter(11) - 1e-10;
					ub(2) = best_parameter(11) + 1e-10;
				end
				if (best_parameter(51) > 1e-10)
					lb(3) = best_parameter(51) - 1e-10;
					ub(3) = best_parameter(51) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@simf_doi_10_1093_nar_gku1388_figure2A_2_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(50) = x_opt(1);
				best_parameter(11) = x_opt(2);
				best_parameter(51) = x_opt(3);
			case model_order(model_id)
				lb = [0.002 0.03 0.01];
				ub = [50 5 10];
				if (best_parameter(49) > 1e-10)
					lb(1) = best_parameter(49) - 1e-10;
					ub(1) = best_parameter(49) + 1e-10;
				end
				if (best_parameter(11) > 1e-10)
					lb(2) = best_parameter(11) - 1e-10;
					ub(2) = best_parameter(11) + 1e-10;
				end
				if (best_parameter(51) > 1e-10)
					lb(3) = best_parameter(51) - 1e-10;
					ub(3) = best_parameter(51) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@simf_doi_10_1093_nar_gku1388_figure2A_3_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(49) = x_opt(1);
				best_parameter(11) = x_opt(2);
				best_parameter(51) = x_opt(3);
			case model_order(model_id)
				lb = [0.002 0.03 0.01];
				ub = [50 5 10];
				if (best_parameter(48) > 1e-10)
					lb(1) = best_parameter(48) - 1e-10;
					ub(1) = best_parameter(48) + 1e-10;
				end
				if (best_parameter(11) > 1e-10)
					lb(2) = best_parameter(11) - 1e-10;
					ub(2) = best_parameter(11) + 1e-10;
				end
				if (best_parameter(51) > 1e-10)
					lb(3) = best_parameter(51) - 1e-10;
					ub(3) = best_parameter(51) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@simf_doi_10_1093_nar_gku1388_figure2A_4_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(48) = x_opt(1);
				best_parameter(11) = x_opt(2);
				best_parameter(51) = x_opt(3);
			case model_order(model_id)
				lb = [0.002 0.03 0.01];
				ub = [50 5 10];
				if (best_parameter(47) > 1e-10)
					lb(1) = best_parameter(47) - 1e-10;
					ub(1) = best_parameter(47) + 1e-10;
				end
				if (best_parameter(11) > 1e-10)
					lb(2) = best_parameter(11) - 1e-10;
					ub(2) = best_parameter(11) + 1e-10;
				end
				if (best_parameter(51) > 1e-10)
					lb(3) = best_parameter(51) - 1e-10;
					ub(3) = best_parameter(51) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@simf_doi_10_1093_nar_gku1388_figure2A_5_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(47) = x_opt(1);
				best_parameter(11) = x_opt(2);
				best_parameter(51) = x_opt(3);
			case model_order(model_id)
				lb = [0.002 0.03 0.01];
				ub = [50 5 10];
				if (best_parameter(52) > 1e-10)
					lb(1) = best_parameter(52) - 1e-10;
					ub(1) = best_parameter(52) + 1e-10;
				end
				if (best_parameter(11) > 1e-10)
					lb(2) = best_parameter(11) - 1e-10;
					ub(2) = best_parameter(11) + 1e-10;
				end
				if (best_parameter(51) > 1e-10)
					lb(3) = best_parameter(51) - 1e-10;
					ub(3) = best_parameter(51) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@simf_doi_10_1093_nar_gku1388_figure2A_6_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(52) = x_opt(1);
				best_parameter(11) = x_opt(2);
				best_parameter(51) = x_opt(3);
			case model_order(model_id)
				lb = [0.002 0.03 0.002 0 1 0.001 0.01 1 0.001];
				ub = [50 5 50 1 4 50 10 4 50];
				if (best_parameter(52) > 1e-10)
					lb(1) = best_parameter(52) - 1e-10;
					ub(1) = best_parameter(52) + 1e-10;
				end
				if (best_parameter(11) > 1e-10)
					lb(2) = best_parameter(11) - 1e-10;
					ub(2) = best_parameter(11) + 1e-10;
				end
				if (best_parameter(54) > 1e-10)
					lb(3) = best_parameter(54) - 1e-10;
					ub(3) = best_parameter(54) + 1e-10;
				end
				if (best_parameter(55) > 1e-10)
					lb(4) = best_parameter(55) - 1e-10;
					ub(4) = best_parameter(55) + 1e-10;
				end
				if (best_parameter(56) > 1e-10)
					lb(5) = best_parameter(56) - 1e-10;
					ub(5) = best_parameter(56) + 1e-10;
				end
				if (best_parameter(53) > 1e-10)
					lb(6) = best_parameter(53) - 1e-10;
					ub(6) = best_parameter(53) + 1e-10;
				end
				if (best_parameter(51) > 1e-10)
					lb(7) = best_parameter(51) - 1e-10;
					ub(7) = best_parameter(51) + 1e-10;
				end
				if (best_parameter(21) > 1e-10)
					lb(8) = best_parameter(21) - 1e-10;
					ub(8) = best_parameter(21) + 1e-10;
				end
				if (best_parameter(16) > 1e-10)
					lb(9) = best_parameter(16) - 1e-10;
					ub(9) = best_parameter(16) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@simf_doi_10_1093_nar_gku1388_figure2B_1_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(52) = x_opt(1);
				best_parameter(11) = x_opt(2);
				best_parameter(54) = x_opt(3);
				best_parameter(55) = x_opt(4);
				best_parameter(56) = x_opt(5);
				best_parameter(53) = x_opt(6);
				best_parameter(51) = x_opt(7);
				best_parameter(21) = x_opt(8);
				best_parameter(16) = x_opt(9);
			case model_order(model_id)
				lb = [0.002 0.03 0.002 0 1 0.001 0.01 1 0.001];
				ub = [50 5 50 1 4 50 10 4 50];
				if (best_parameter(47) > 1e-10)
					lb(1) = best_parameter(47) - 1e-10;
					ub(1) = best_parameter(47) + 1e-10;
				end
				if (best_parameter(11) > 1e-10)
					lb(2) = best_parameter(11) - 1e-10;
					ub(2) = best_parameter(11) + 1e-10;
				end
				if (best_parameter(54) > 1e-10)
					lb(3) = best_parameter(54) - 1e-10;
					ub(3) = best_parameter(54) + 1e-10;
				end
				if (best_parameter(55) > 1e-10)
					lb(4) = best_parameter(55) - 1e-10;
					ub(4) = best_parameter(55) + 1e-10;
				end
				if (best_parameter(56) > 1e-10)
					lb(5) = best_parameter(56) - 1e-10;
					ub(5) = best_parameter(56) + 1e-10;
				end
				if (best_parameter(53) > 1e-10)
					lb(6) = best_parameter(53) - 1e-10;
					ub(6) = best_parameter(53) + 1e-10;
				end
				if (best_parameter(51) > 1e-10)
					lb(7) = best_parameter(51) - 1e-10;
					ub(7) = best_parameter(51) + 1e-10;
				end
				if (best_parameter(21) > 1e-10)
					lb(8) = best_parameter(21) - 1e-10;
					ub(8) = best_parameter(21) + 1e-10;
				end
				if (best_parameter(16) > 1e-10)
					lb(9) = best_parameter(16) - 1e-10;
					ub(9) = best_parameter(16) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@simf_doi_10_1093_nar_gku1388_figure2B_2_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(47) = x_opt(1);
				best_parameter(11) = x_opt(2);
				best_parameter(54) = x_opt(3);
				best_parameter(55) = x_opt(4);
				best_parameter(56) = x_opt(5);
				best_parameter(53) = x_opt(6);
				best_parameter(51) = x_opt(7);
				best_parameter(21) = x_opt(8);
				best_parameter(16) = x_opt(9);
			case model_order(model_id)
				lb = [0.002 0 1 0.001 0.03 0.01 0.002 0.03 1 0.001];
				ub = [50 1 4 50 5 10 50 5 4 50];
				if (best_parameter(5) > 1e-10)
					lb(1) = best_parameter(5) - 1e-10;
					ub(1) = best_parameter(5) + 1e-10;
				end
				if (best_parameter(7) > 1e-10)
					lb(2) = best_parameter(7) - 1e-10;
					ub(2) = best_parameter(7) + 1e-10;
				end
				if (best_parameter(9) > 1e-10)
					lb(3) = best_parameter(9) - 1e-10;
					ub(3) = best_parameter(9) + 1e-10;
				end
				if (best_parameter(4) > 1e-10)
					lb(4) = best_parameter(4) - 1e-10;
					ub(4) = best_parameter(4) + 1e-10;
				end
				if (best_parameter(57) > 1e-10)
					lb(5) = best_parameter(57) - 1e-10;
					ub(5) = best_parameter(57) + 1e-10;
				end
				if (best_parameter(59) > 1e-10)
					lb(6) = best_parameter(59) - 1e-10;
					ub(6) = best_parameter(59) + 1e-10;
				end
				if (best_parameter(52) > 1e-10)
					lb(7) = best_parameter(52) - 1e-10;
					ub(7) = best_parameter(52) + 1e-10;
				end
				if (best_parameter(58) > 1e-10)
					lb(8) = best_parameter(58) - 1e-10;
					ub(8) = best_parameter(58) + 1e-10;
				end
				if (best_parameter(8) > 1e-10)
					lb(9) = best_parameter(8) - 1e-10;
					ub(9) = best_parameter(8) + 1e-10;
				end
				if (best_parameter(3) > 1e-10)
					lb(10) = best_parameter(3) - 1e-10;
					ub(10) = best_parameter(3) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@simf_doi_10_1038_nature09565_figure1C_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(5) = x_opt(1);
				best_parameter(7) = x_opt(2);
				best_parameter(9) = x_opt(3);
				best_parameter(4) = x_opt(4);
				best_parameter(57) = x_opt(5);
				best_parameter(59) = x_opt(6);
				best_parameter(52) = x_opt(7);
				best_parameter(58) = x_opt(8);
				best_parameter(8) = x_opt(9);
				best_parameter(3) = x_opt(10);
			case model_order(model_id)
				lb = [0.002 0 1 0.001 0.03 0.01 0.002 0.03 1 0.001];
				ub = [50 1 4 50 5 10 50 5 4 50];
				if (best_parameter(18) > 1e-10)
					lb(1) = best_parameter(18) - 1e-10;
					ub(1) = best_parameter(18) + 1e-10;
				end
				if (best_parameter(20) > 1e-10)
					lb(2) = best_parameter(20) - 1e-10;
					ub(2) = best_parameter(20) + 1e-10;
				end
				if (best_parameter(22) > 1e-10)
					lb(3) = best_parameter(22) - 1e-10;
					ub(3) = best_parameter(22) + 1e-10;
				end
				if (best_parameter(17) > 1e-10)
					lb(4) = best_parameter(17) - 1e-10;
					ub(4) = best_parameter(17) + 1e-10;
				end
				if (best_parameter(57) > 1e-10)
					lb(5) = best_parameter(57) - 1e-10;
					ub(5) = best_parameter(57) + 1e-10;
				end
				if (best_parameter(59) > 1e-10)
					lb(6) = best_parameter(59) - 1e-10;
					ub(6) = best_parameter(59) + 1e-10;
				end
				if (best_parameter(52) > 1e-10)
					lb(7) = best_parameter(52) - 1e-10;
					ub(7) = best_parameter(52) + 1e-10;
				end
				if (best_parameter(38) > 1e-10)
					lb(8) = best_parameter(38) - 1e-10;
					ub(8) = best_parameter(38) + 1e-10;
				end
				if (best_parameter(21) > 1e-10)
					lb(9) = best_parameter(21) - 1e-10;
					ub(9) = best_parameter(21) + 1e-10;
				end
				if (best_parameter(16) > 1e-10)
					lb(10) = best_parameter(16) - 1e-10;
					ub(10) = best_parameter(16) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@simf_doi_10_1038_nature09565_figureS1C_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(18) = x_opt(1);
				best_parameter(20) = x_opt(2);
				best_parameter(22) = x_opt(3);
				best_parameter(17) = x_opt(4);
				best_parameter(57) = x_opt(5);
				best_parameter(59) = x_opt(6);
				best_parameter(52) = x_opt(7);
				best_parameter(38) = x_opt(8);
				best_parameter(21) = x_opt(9);
				best_parameter(16) = x_opt(10);
			case model_order(model_id)
				lb = [0.002 0 1 0.001 0.03 0.01 0.002 1 0.001];
				ub = [50 1 4 50 5 10 50 4 50];
				if (best_parameter(62) > 1e-10)
					lb(1) = best_parameter(62) - 1e-10;
					ub(1) = best_parameter(62) + 1e-10;
				end
				if (best_parameter(63) > 1e-10)
					lb(2) = best_parameter(63) - 1e-10;
					ub(2) = best_parameter(63) + 1e-10;
				end
				if (best_parameter(64) > 1e-10)
					lb(3) = best_parameter(64) - 1e-10;
					ub(3) = best_parameter(64) + 1e-10;
				end
				if (best_parameter(61) > 1e-10)
					lb(4) = best_parameter(61) - 1e-10;
					ub(4) = best_parameter(61) + 1e-10;
				end
				if (best_parameter(60) > 1e-10)
					lb(5) = best_parameter(60) - 1e-10;
					ub(5) = best_parameter(60) + 1e-10;
				end
				if (best_parameter(65) > 1e-10)
					lb(6) = best_parameter(65) - 1e-10;
					ub(6) = best_parameter(65) + 1e-10;
				end
				if (best_parameter(39) > 1e-10)
					lb(7) = best_parameter(39) - 1e-10;
					ub(7) = best_parameter(39) + 1e-10;
				end
				if (best_parameter(30) > 1e-10)
					lb(8) = best_parameter(30) - 1e-10;
					ub(8) = best_parameter(30) + 1e-10;
				end
				if (best_parameter(25) > 1e-10)
					lb(9) = best_parameter(25) - 1e-10;
					ub(9) = best_parameter(25) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@simf_doi_10_1038_nbt_2401_figureS11_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(62) = x_opt(1);
				best_parameter(63) = x_opt(2);
				best_parameter(64) = x_opt(3);
				best_parameter(61) = x_opt(4);
				best_parameter(60) = x_opt(5);
				best_parameter(65) = x_opt(6);
				best_parameter(39) = x_opt(7);
				best_parameter(30) = x_opt(8);
				best_parameter(25) = x_opt(9);
			case model_order(model_id)
				lb = [0.002 0 1 0.001 0.03 0.01 0.002 1 0.001];
				ub = [50 1 4 50 5 10 50 4 50];
				if (best_parameter(5) > 1e-10)
					lb(1) = best_parameter(5) - 1e-10;
					ub(1) = best_parameter(5) + 1e-10;
				end
				if (best_parameter(7) > 1e-10)
					lb(2) = best_parameter(7) - 1e-10;
					ub(2) = best_parameter(7) + 1e-10;
				end
				if (best_parameter(9) > 1e-10)
					lb(3) = best_parameter(9) - 1e-10;
					ub(3) = best_parameter(9) + 1e-10;
				end
				if (best_parameter(4) > 1e-10)
					lb(4) = best_parameter(4) - 1e-10;
					ub(4) = best_parameter(4) + 1e-10;
				end
				if (best_parameter(60) > 1e-10)
					lb(5) = best_parameter(60) - 1e-10;
					ub(5) = best_parameter(60) + 1e-10;
				end
				if (best_parameter(65) > 1e-10)
					lb(6) = best_parameter(65) - 1e-10;
					ub(6) = best_parameter(65) + 1e-10;
				end
				if (best_parameter(39) > 1e-10)
					lb(7) = best_parameter(39) - 1e-10;
					ub(7) = best_parameter(39) + 1e-10;
				end
				if (best_parameter(8) > 1e-10)
					lb(8) = best_parameter(8) - 1e-10;
					ub(8) = best_parameter(8) + 1e-10;
				end
				if (best_parameter(3) > 1e-10)
					lb(9) = best_parameter(3) - 1e-10;
					ub(9) = best_parameter(3) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@simf_doi_10_1038_nbt_2401_figureS13_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(5) = x_opt(1);
				best_parameter(7) = x_opt(2);
				best_parameter(9) = x_opt(3);
				best_parameter(4) = x_opt(4);
				best_parameter(60) = x_opt(5);
				best_parameter(65) = x_opt(6);
				best_parameter(39) = x_opt(7);
				best_parameter(8) = x_opt(8);
				best_parameter(3) = x_opt(9);
			case model_order(model_id)
				lb = [0.03 0.002 0.002 0 1 0.001 0.03 0.01 1 0.001];
				ub = [5 50 50 1 4 50 5 10 4 50];
				if (best_parameter(2) > 1e-10)
					lb(1) = best_parameter(2) - 1e-10;
					ub(1) = best_parameter(2) + 1e-10;
				end
				if (best_parameter(6) > 1e-10)
					lb(2) = best_parameter(6) - 1e-10;
					ub(2) = best_parameter(6) + 1e-10;
				end
				if (best_parameter(5) > 1e-10)
					lb(3) = best_parameter(5) - 1e-10;
					ub(3) = best_parameter(5) + 1e-10;
				end
				if (best_parameter(7) > 1e-10)
					lb(4) = best_parameter(7) - 1e-10;
					ub(4) = best_parameter(7) + 1e-10;
				end
				if (best_parameter(9) > 1e-10)
					lb(5) = best_parameter(9) - 1e-10;
					ub(5) = best_parameter(9) + 1e-10;
				end
				if (best_parameter(4) > 1e-10)
					lb(6) = best_parameter(4) - 1e-10;
					ub(6) = best_parameter(4) + 1e-10;
				end
				if (best_parameter(66) > 1e-10)
					lb(7) = best_parameter(66) - 1e-10;
					ub(7) = best_parameter(66) + 1e-10;
				end
				if (best_parameter(67) > 1e-10)
					lb(8) = best_parameter(67) - 1e-10;
					ub(8) = best_parameter(67) + 1e-10;
				end
				if (best_parameter(8) > 1e-10)
					lb(9) = best_parameter(8) - 1e-10;
					ub(9) = best_parameter(8) + 1e-10;
				end
				if (best_parameter(3) > 1e-10)
					lb(10) = best_parameter(3) - 1e-10;
					ub(10) = best_parameter(3) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@simf_doi_10_1038_nature11516_figureS8A_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(2) = x_opt(1);
				best_parameter(6) = x_opt(2);
				best_parameter(5) = x_opt(3);
				best_parameter(7) = x_opt(4);
				best_parameter(9) = x_opt(5);
				best_parameter(4) = x_opt(6);
				best_parameter(66) = x_opt(7);
				best_parameter(67) = x_opt(8);
				best_parameter(8) = x_opt(9);
				best_parameter(3) = x_opt(10);
			case model_order(model_id)
				lb = [0.03 0.002 0.002 0 1 0.001 0.03 0.01 1 0.001];
				ub = [5 50 50 1 4 50 5 10 4 50];
				if (best_parameter(24) > 1e-10)
					lb(1) = best_parameter(24) - 1e-10;
					ub(1) = best_parameter(24) + 1e-10;
				end
				if (best_parameter(28) > 1e-10)
					lb(2) = best_parameter(28) - 1e-10;
					ub(2) = best_parameter(28) + 1e-10;
				end
				if (best_parameter(62) > 1e-10)
					lb(3) = best_parameter(62) - 1e-10;
					ub(3) = best_parameter(62) + 1e-10;
				end
				if (best_parameter(63) > 1e-10)
					lb(4) = best_parameter(63) - 1e-10;
					ub(4) = best_parameter(63) + 1e-10;
				end
				if (best_parameter(64) > 1e-10)
					lb(5) = best_parameter(64) - 1e-10;
					ub(5) = best_parameter(64) + 1e-10;
				end
				if (best_parameter(61) > 1e-10)
					lb(6) = best_parameter(61) - 1e-10;
					ub(6) = best_parameter(61) + 1e-10;
				end
				if (best_parameter(66) > 1e-10)
					lb(7) = best_parameter(66) - 1e-10;
					ub(7) = best_parameter(66) + 1e-10;
				end
				if (best_parameter(67) > 1e-10)
					lb(8) = best_parameter(67) - 1e-10;
					ub(8) = best_parameter(67) + 1e-10;
				end
				if (best_parameter(30) > 1e-10)
					lb(9) = best_parameter(30) - 1e-10;
					ub(9) = best_parameter(30) + 1e-10;
				end
				if (best_parameter(25) > 1e-10)
					lb(10) = best_parameter(25) - 1e-10;
					ub(10) = best_parameter(25) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@simf_doi_10_1038_nature11516_figureS8B_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(24) = x_opt(1);
				best_parameter(28) = x_opt(2);
				best_parameter(62) = x_opt(3);
				best_parameter(63) = x_opt(4);
				best_parameter(64) = x_opt(5);
				best_parameter(61) = x_opt(6);
				best_parameter(66) = x_opt(7);
				best_parameter(67) = x_opt(8);
				best_parameter(30) = x_opt(9);
				best_parameter(25) = x_opt(10);
			case model_order(model_id)
				lb = [0.002 0.03 0.002 0 1 0.001 0.03 0.01 1 0.001];
				ub = [50 5 50 1 4 50 5 10 4 50];
				if (best_parameter(36) > 1e-10)
					lb(1) = best_parameter(36) - 1e-10;
					ub(1) = best_parameter(36) + 1e-10;
				end
				if (best_parameter(24) > 1e-10)
					lb(2) = best_parameter(24) - 1e-10;
					ub(2) = best_parameter(24) + 1e-10;
				end
				if (best_parameter(70) > 1e-10)
					lb(3) = best_parameter(70) - 1e-10;
					ub(3) = best_parameter(70) + 1e-10;
				end
				if (best_parameter(71) > 1e-10)
					lb(4) = best_parameter(71) - 1e-10;
					ub(4) = best_parameter(71) + 1e-10;
				end
				if (best_parameter(72) > 1e-10)
					lb(5) = best_parameter(72) - 1e-10;
					ub(5) = best_parameter(72) + 1e-10;
				end
				if (best_parameter(69) > 1e-10)
					lb(6) = best_parameter(69) - 1e-10;
					ub(6) = best_parameter(69) + 1e-10;
				end
				if (best_parameter(68) > 1e-10)
					lb(7) = best_parameter(68) - 1e-10;
					ub(7) = best_parameter(68) + 1e-10;
				end
				if (best_parameter(73) > 1e-10)
					lb(8) = best_parameter(73) - 1e-10;
					ub(8) = best_parameter(73) + 1e-10;
				end
				if (best_parameter(30) > 1e-10)
					lb(9) = best_parameter(30) - 1e-10;
					ub(9) = best_parameter(30) + 1e-10;
				end
				if (best_parameter(25) > 1e-10)
					lb(10) = best_parameter(25) - 1e-10;
					ub(10) = best_parameter(25) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@simf_doi_10_1038_nbt_2149_Supplement_page17_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(36) = x_opt(1);
				best_parameter(24) = x_opt(2);
				best_parameter(70) = x_opt(3);
				best_parameter(71) = x_opt(4);
				best_parameter(72) = x_opt(5);
				best_parameter(69) = x_opt(6);
				best_parameter(68) = x_opt(7);
				best_parameter(73) = x_opt(8);
				best_parameter(30) = x_opt(9);
				best_parameter(25) = x_opt(10);
		end
	end
end

function seq_LOOCV(N1, N2)
	LOOCV_seq_data_fileID = fopen('Plot_data/LOOCV_sim_data_Hill.dat','w');
	LOO_model_name_list = {'doi_10_1186_1754_1611_3_4_figure3_1', 'doi_10_1186_1754_1611_3_4_figure3_2', 'doi_10_1186_1754_1611_3_4_figure3_6', 'doi_10_1186_1754_1611_3_4_figure3_7', 'doi_10_1093_nar_gku593_figure3D_2', 'doi_10_1093_nar_gku593_figure3D_3', 'doi_10_1093_nar_gku593_figure3D_4', 'doi_10_1093_nar_gku593_figure3D_5', 'doi_10_1093_nar_gku593_figure3D_6', 'doi_10_1093_nar_gku593_figureS6', 'doi_10_1093_nar_gku1388_figure2A_1', 'doi_10_1093_nar_gku1388_figure2A_2', 'doi_10_1093_nar_gku1388_figure2A_3', 'doi_10_1093_nar_gku1388_figure2A_4', 'doi_10_1093_nar_gku1388_figure2A_5', 'doi_10_1093_nar_gku1388_figure2A_6', 'doi_10_1093_nar_gku1388_figure2B_1', 'doi_10_1093_nar_gku1388_figure2B_2', 'doi_10_1038_nature09565_figureS1C', 'doi_10_1038_nbt_2401_figureS11', 'doi_10_1038_nbt_2401_figureS13', 'doi_10_1038_nature11516_figureS8A', 'doi_10_1038_nature11516_figureS8B'};
	LOOCV_error = zeros(23,1);
	error_opt = 1e6;
	x_opt = zeros(1,73);
	for i = 1:N1
		model_order = LOOCV_randperm(35,6);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	desired_output = 0.69;
	output = simf_doi_10_1186_1754_1611_3_4_figure3_1(x_opt);
	fprintf(LOOCV_seq_data_fileID, '#doi_10_1186_1754_1611_3_4_figure3_1\n');
	for i = 1:length(output)
		fprintf(LOOCV_seq_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	error_opt = 1e6;
	x_opt = zeros(1,73);
	for i = 1:N1
		model_order = LOOCV_randperm(35,7);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	desired_output = 0.018;
	output = simf_doi_10_1186_1754_1611_3_4_figure3_2(x_opt);
	fprintf(LOOCV_seq_data_fileID, '#doi_10_1186_1754_1611_3_4_figure3_2\n');
	for i = 1:length(output)
		fprintf(LOOCV_seq_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	error_opt = 1e6;
	x_opt = zeros(1,73);
	for i = 1:N1
		model_order = LOOCV_randperm(35,11);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	desired_output = 0.897;
	output = simf_doi_10_1186_1754_1611_3_4_figure3_6(x_opt);
	fprintf(LOOCV_seq_data_fileID, '#doi_10_1186_1754_1611_3_4_figure3_6\n');
	for i = 1:length(output)
		fprintf(LOOCV_seq_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	error_opt = 1e6;
	x_opt = zeros(1,73);
	for i = 1:N1
		model_order = LOOCV_randperm(35,12);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	desired_output = 1;
	output = simf_doi_10_1186_1754_1611_3_4_figure3_7(x_opt);
	fprintf(LOOCV_seq_data_fileID, '#doi_10_1186_1754_1611_3_4_figure3_7\n');
	for i = 1:length(output)
		fprintf(LOOCV_seq_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	error_opt = 1e6;
	x_opt = zeros(1,73);
	for i = 1:N1
		model_order = LOOCV_randperm(35,14);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	desired_output = 0.018;
	output = simf_doi_10_1093_nar_gku593_figure3D_2(x_opt);
	fprintf(LOOCV_seq_data_fileID, '#doi_10_1093_nar_gku593_figure3D_2\n');
	for i = 1:length(output)
		fprintf(LOOCV_seq_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	error_opt = 1e6;
	x_opt = zeros(1,73);
	for i = 1:N1
		model_order = LOOCV_randperm(35,15);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	desired_output = 0.036;
	output = simf_doi_10_1093_nar_gku593_figure3D_3(x_opt);
	fprintf(LOOCV_seq_data_fileID, '#doi_10_1093_nar_gku593_figure3D_3\n');
	for i = 1:length(output)
		fprintf(LOOCV_seq_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	error_opt = 1e6;
	x_opt = zeros(1,73);
	for i = 1:N1
		model_order = LOOCV_randperm(35,16);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	desired_output = 0.05;
	output = simf_doi_10_1093_nar_gku593_figure3D_4(x_opt);
	fprintf(LOOCV_seq_data_fileID, '#doi_10_1093_nar_gku593_figure3D_4\n');
	for i = 1:length(output)
		fprintf(LOOCV_seq_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	error_opt = 1e6;
	x_opt = zeros(1,73);
	for i = 1:N1
		model_order = LOOCV_randperm(35,17);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	desired_output = 0.086;
	output = simf_doi_10_1093_nar_gku593_figure3D_5(x_opt);
	fprintf(LOOCV_seq_data_fileID, '#doi_10_1093_nar_gku593_figure3D_5\n');
	for i = 1:length(output)
		fprintf(LOOCV_seq_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	error_opt = 1e6;
	x_opt = zeros(1,73);
	for i = 1:N1
		model_order = LOOCV_randperm(35,18);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	desired_output = 0.214;
	output = simf_doi_10_1093_nar_gku593_figure3D_6(x_opt);
	fprintf(LOOCV_seq_data_fileID, '#doi_10_1093_nar_gku593_figure3D_6\n');
	for i = 1:length(output)
		fprintf(LOOCV_seq_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	error_opt = 1e6;
	x_opt = zeros(1,73);
	for i = 1:N1
		model_order = LOOCV_randperm(35,19);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2];
	desired_output = [0 0.002 0.041 0.479 0.949 0.997 1 1 1];
	output = simf_doi_10_1093_nar_gku593_figureS6(x_opt, arabinose);
	fprintf(LOOCV_seq_data_fileID, '#doi_10_1093_nar_gku593_figureS6\n');
	for i = 1:length(output)
		fprintf(LOOCV_seq_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	error_opt = 1e6;
	x_opt = zeros(1,73);
	for i = 1:N1
		model_order = LOOCV_randperm(35,20);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	desired_output = 1;
	output = simf_doi_10_1093_nar_gku1388_figure2A_1(x_opt);
	fprintf(LOOCV_seq_data_fileID, '#doi_10_1093_nar_gku1388_figure2A_1\n');
	for i = 1:length(output)
		fprintf(LOOCV_seq_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	error_opt = 1e6;
	x_opt = zeros(1,73);
	for i = 1:N1
		model_order = LOOCV_randperm(35,21);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	desired_output = 0.215;
	output = simf_doi_10_1093_nar_gku1388_figure2A_2(x_opt);
	fprintf(LOOCV_seq_data_fileID, '#doi_10_1093_nar_gku1388_figure2A_2\n');
	for i = 1:length(output)
		fprintf(LOOCV_seq_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	error_opt = 1e6;
	x_opt = zeros(1,73);
	for i = 1:N1
		model_order = LOOCV_randperm(35,22);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	desired_output = 0.095;
	output = simf_doi_10_1093_nar_gku1388_figure2A_3(x_opt);
	fprintf(LOOCV_seq_data_fileID, '#doi_10_1093_nar_gku1388_figure2A_3\n');
	for i = 1:length(output)
		fprintf(LOOCV_seq_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	error_opt = 1e6;
	x_opt = zeros(1,73);
	for i = 1:N1
		model_order = LOOCV_randperm(35,23);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	desired_output = 0.046;
	output = simf_doi_10_1093_nar_gku1388_figure2A_4(x_opt);
	fprintf(LOOCV_seq_data_fileID, '#doi_10_1093_nar_gku1388_figure2A_4\n');
	for i = 1:length(output)
		fprintf(LOOCV_seq_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	error_opt = 1e6;
	x_opt = zeros(1,73);
	for i = 1:N1
		model_order = LOOCV_randperm(35,24);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	desired_output = 0.018;
	output = simf_doi_10_1093_nar_gku1388_figure2A_5(x_opt);
	fprintf(LOOCV_seq_data_fileID, '#doi_10_1093_nar_gku1388_figure2A_5\n');
	for i = 1:length(output)
		fprintf(LOOCV_seq_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	error_opt = 1e6;
	x_opt = zeros(1,73);
	for i = 1:N1
		model_order = LOOCV_randperm(35,25);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	desired_output = 0.003;
	output = simf_doi_10_1093_nar_gku1388_figure2A_6(x_opt);
	fprintf(LOOCV_seq_data_fileID, '#doi_10_1093_nar_gku1388_figure2A_6\n');
	for i = 1:length(output)
		fprintf(LOOCV_seq_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	error_opt = 1e6;
	x_opt = zeros(1,73);
	for i = 1:N1
		model_order = LOOCV_randperm(35,26);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	aTc = [2 4 10 25 50 100 200];
	desired_output = [0.002 0.014 0.167 0.666 0.817 0.841 0.831];
	output = simf_doi_10_1093_nar_gku1388_figure2B_1(x_opt, aTc);
	fprintf(LOOCV_seq_data_fileID, '#doi_10_1093_nar_gku1388_figure2B_1\n');
	for i = 1:length(output)
		fprintf(LOOCV_seq_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	error_opt = 1e6;
	x_opt = zeros(1,73);
	for i = 1:N1
		model_order = LOOCV_randperm(35,27);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	aTc = [10 25 50 100 200];
	desired_output = [0.001 0.022 0.176 0.411 0.461];
	output = simf_doi_10_1093_nar_gku1388_figure2B_2(x_opt, aTc);
	fprintf(LOOCV_seq_data_fileID, '#doi_10_1093_nar_gku1388_figure2B_2\n');
	for i = 1:length(output)
		fprintf(LOOCV_seq_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	error_opt = 1e6;
	x_opt = zeros(1,73);
	for i = 1:N1
		model_order = LOOCV_randperm(35,29);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	aTc = [0 0.025 0.25 2.5 25 100 250];
	desired_output = [0.005 0.006 0.007 0.03 0.297 0.388 0.398];
	output = simf_doi_10_1038_nature09565_figureS1C(x_opt, aTc);
	fprintf(LOOCV_seq_data_fileID, '#doi_10_1038_nature09565_figureS1C\n');
	for i = 1:length(output)
		fprintf(LOOCV_seq_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	error_opt = 1e6;
	x_opt = zeros(1,73);
	for i = 1:N1
		model_order = LOOCV_randperm(35,30);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	IPTG = [0 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1];
	desired_output = [0.025 0.025 0.035 0.055 0.115 0.2 0.3 0.4 0.65 1];
	output = simf_doi_10_1038_nbt_2401_figureS11(x_opt, IPTG);
	fprintf(LOOCV_seq_data_fileID, '#doi_10_1038_nbt_2401_figureS11\n');
	for i = 1:length(output)
		fprintf(LOOCV_seq_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	error_opt = 1e6;
	x_opt = zeros(1,73);
	for i = 1:N1
		model_order = LOOCV_randperm(35,31);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	arabinose = [0.1 1 2 5 7 10 12.5 25 37.5 50 62.5];
	desired_output = [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425];
	output = simf_doi_10_1038_nbt_2401_figureS13(x_opt, arabinose);
	fprintf(LOOCV_seq_data_fileID, '#doi_10_1038_nbt_2401_figureS13\n');
	for i = 1:length(output)
		fprintf(LOOCV_seq_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	error_opt = 1e6;
	x_opt = zeros(1,73);
	for i = 1:N1
		model_order = LOOCV_randperm(35,32);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	arabinose = [0 0 0.002 0.01 0.04 0.16 1 4 25];
	desired_output = [0.003 0.003 0.003 0.017 0.23 0.82 0.939 0.941 0.941];
	output = simf_doi_10_1038_nature11516_figureS8A(x_opt, arabinose);
	fprintf(LOOCV_seq_data_fileID, '#doi_10_1038_nature11516_figureS8A\n');
	for i = 1:length(output)
		fprintf(LOOCV_seq_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	error_opt = 1e6;
	x_opt = zeros(1,73);
	for i = 1:N1
		model_order = LOOCV_randperm(35,33);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	IPTG = [0 0 0 0.001 0.004 0.025 0.1 0.5 2.5];
	desired_output = [0.027 0.027 0.027 0.03 0.067 0.46 0.884 0.991 1];
	output = simf_doi_10_1038_nature11516_figureS8B(x_opt, IPTG);
	fprintf(LOOCV_seq_data_fileID, '#doi_10_1038_nature11516_figureS8B\n');
	for i = 1:length(output)
		fprintf(LOOCV_seq_data_fileID, '%f\t%f\n', desired_output(i), output(i));
	end
	LOOCV_fileID = fopen('Plot_data/LOOCV_seq.dat','w');
	for i = 1:length(LOO_model_name_list)
		fprintf(LOOCV_fileID, '%s\t', LOO_model_name_list{i});
		fprintf(LOOCV_fileID, '%f\n', LOOCV_error(i));
	end
	fclose(LOOCV_fileID);
end

