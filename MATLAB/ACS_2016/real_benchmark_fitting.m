function real_bechmark_fitting
	tic;
	%separate_fitting(10);
	sep_t = toc;
	tic;
	%seq_x_opt = sequential_fitting(5,10);
	seq_t = toc;
	seq_LOOCV(5,10);
	tic;
	%sim_x_opt = simultaneous_fitting(10);
	sim_t = toc;
	%sim_LOOCV(100);
	%print_a_matrix('Plot_data/running_time.txt', [sep_t seq_t sim_t]);
	%print_plotting_data(seq_x_opt, sim_x_opt);
end
function output = simf_doi_10_1073_pnas_0606717104_figure4_6_1(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pLAC = x(1);
	ALPHA_RBS_placI = x(2);
	K_IPTG = x(3);
	K_pLAC = x(4);
	alpha_pLAC = x(5);
	alpha_placI = x(6);
	beta_pLAC = x(7);
	n_IPTG = x(8);
	n_pLAC = x(9);
	scale_lacZ_doi_10_1073_pnas_0606717104 = x(10);
	for i = 1:length(input);
		placI = 1*alpha_placI;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placI*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pLAC = 1*(beta_pLAC + (alpha_pLAC - beta_pLAC)*K_pLAC^n_pLAC/(K_pLAC^n_pLAC + lacI^n_pLAC));
		lacZ = scale_lacZ_doi_10_1073_pnas_0606717104*ALPHA_RBS_pLAC*pLAC;
		output(i) = lacZ;
	end
end

function output = sepf_doi_10_1073_pnas_0606717104_figure4_6_1(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pLAC = x(1);
	ALPHA_RBS_placI = x(2);
	K_IPTG = x(3);
	K_pLAC = x(4);
	alpha_pLAC = x(5);
	alpha_placI = x(6);
	beta_pLAC = x(7);
	n_IPTG = x(8);
	n_pLAC = x(9);
	scale_lacZ_doi_10_1073_pnas_0606717104 = x(10);
	for i = 1:length(input);
		placI = 1*alpha_placI;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placI*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pLAC = 1*(beta_pLAC + (alpha_pLAC - beta_pLAC)*K_pLAC^n_pLAC/(K_pLAC^n_pLAC + lacI^n_pLAC));
		lacZ = scale_lacZ_doi_10_1073_pnas_0606717104*ALPHA_RBS_pLAC*pLAC;
		output(i) = lacZ;
	end
end

function error = doi_10_1073_pnas_0606717104_figure4_6_1_error(x)
	IPTG = [0 0.005 0.01 0.02 0.05 0.1 0.2 1];
	lacZ = [0.001 0.003 0.043 0.5 0.984 0.999 1 1];
	diff = sepf_doi_10_1073_pnas_0606717104_figure4_6_1(x, IPTG) - lacZ;
	error = repmat(diff,1,5);
end

function output = simf_doi__10_1093_nar_25_6_1203_figure4a_1(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pN25 = x(11);
	ALPHA_RBSII = x(12);
	K_aTc = x(13);
	K_pLtetO_1 = x(14);
	alpha_pLtetO_1 = x(15);
	alpha_pN25 = x(16);
	beta_pLtetO_1 = x(17);
	n_aTc = x(18);
	n_pLtetO_1 = x(19);
	scale_luc_doi__10_1093_nar_25_6_1203 = x(20);
	for i = 1:length(input);
		pN25 = 1*alpha_pN25;
		aTc = input(i);
		tetR = ALPHA_RBS_pN25*pN25*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pLtetO_1 = 15*(beta_pLtetO_1 + (alpha_pLtetO_1 - beta_pLtetO_1)*K_pLtetO_1^n_pLtetO_1/(K_pLtetO_1^n_pLtetO_1 + tetR^n_pLtetO_1));
		luc = scale_luc_doi__10_1093_nar_25_6_1203*ALPHA_RBSII*pLtetO_1;
		output(i) = luc;
	end
end

function output = sepf_doi__10_1093_nar_25_6_1203_figure4a_1(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pN25 = x(1);
	ALPHA_RBSII = x(2);
	K_aTc = x(3);
	K_pLtetO_1 = x(4);
	alpha_pLtetO_1 = x(5);
	alpha_pN25 = x(6);
	beta_pLtetO_1 = x(7);
	n_aTc = x(8);
	n_pLtetO_1 = x(9);
	scale_luc_doi__10_1093_nar_25_6_1203 = x(10);
	for i = 1:length(input);
		pN25 = 1*alpha_pN25;
		aTc = input(i);
		tetR = ALPHA_RBS_pN25*pN25*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pLtetO_1 = 15*(beta_pLtetO_1 + (alpha_pLtetO_1 - beta_pLtetO_1)*K_pLtetO_1^n_pLtetO_1/(K_pLtetO_1^n_pLtetO_1 + tetR^n_pLtetO_1));
		luc = scale_luc_doi__10_1093_nar_25_6_1203*ALPHA_RBSII*pLtetO_1;
		output(i) = luc;
	end
end

function error = doi__10_1093_nar_25_6_1203_figure4a_1_error(x)
	aTc = [0 1 2 3 4 5 6 7 8 9 10 12 13 17 18 20 25 35 60];
	luc = [0 0.001 0.001 0.001 0.002 0.004 0.007 0.011 0.014 0.032 0.036 0.054 0.089 0.643 0.893 1 1 1 1];
	diff = sepf_doi__10_1093_nar_25_6_1203_figure4a_1(x, aTc) - luc;
	error = repmat(diff,1,5);
end

function output = simf_doi_10_1128_AEM_00791_07_figure3a_1(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pBAD = x(21);
	ALPHA_RBS_pC = x(22);
	K_arabinose = x(23);
	K_pBAD = x(24);
	alpha_pBAD = x(25);
	alpha_pC = x(26);
	beta_pBAD = x(27);
	n_arabinose = x(28);
	n_pBAD = x(29);
	scale_GFPuv_doi_10_1128_AEM_00791_07 = x(30);
	for i = 1:length(input);
		pC = 20*alpha_pC;
		arabinose = input(i);
		araC = ALPHA_RBS_pC*pC*K_arabinose^n_arabinose/(K_arabinose^n_arabinose + arabinose^n_arabinose);
		pBAD = 20*(beta_pBAD + (alpha_pBAD - beta_pBAD)*K_pBAD^n_pBAD/(K_pBAD^n_pBAD + araC^n_pBAD));
		GFPuv = scale_GFPuv_doi_10_1128_AEM_00791_07*ALPHA_RBS_pBAD*pBAD;
		output(i) = GFPuv;
	end
end

function output = sepf_doi_10_1128_AEM_00791_07_figure3a_1(x, input)
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
	for i = 1:length(input);
		pC = 20*alpha_pC;
		arabinose = input(i);
		araC = ALPHA_RBS_pC*pC*K_arabinose^n_arabinose/(K_arabinose^n_arabinose + arabinose^n_arabinose);
		pBAD = 20*(beta_pBAD + (alpha_pBAD - beta_pBAD)*K_pBAD^n_pBAD/(K_pBAD^n_pBAD + araC^n_pBAD));
		GFPuv = scale_GFPuv_doi_10_1128_AEM_00791_07*ALPHA_RBS_pBAD*pBAD;
		output(i) = GFPuv;
	end
end

function error = doi_10_1128_AEM_00791_07_figure3a_1_error(x)
	arabinose = [0.003 0.005 0.05 0.5 1];
	GFPuv = [0.044 0.089 0.556 1 1];
	diff = sepf_doi_10_1128_AEM_00791_07_figure3a_1(x, arabinose) - GFPuv;
	error = repmat(diff,1,5);
end

function output = simf_doi_10_1038_nature12148_figure16a_1(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0030 = x(31);
	K_arabinose = x(23);
	K_pBAD = x(24);
	alpha_pBAD = x(25);
	alpha_pLlacO_1 = x(32);
	beta_pBAD = x(27);
	n_arabinose = x(28);
	n_pBAD = x(29);
	scale_mCherry_doi_10_1038_nature12148 = x(33);
	for i = 1:length(input);
		pLlacO_1 = 5*alpha_pLlacO_1;
		arabinose = input(i);
		araC = ALPHA_BBa_B0030*pLlacO_1*K_arabinose^n_arabinose/(K_arabinose^n_arabinose + arabinose^n_arabinose);
		pBAD = 15*(beta_pBAD + (alpha_pBAD - beta_pBAD)*K_pBAD^n_pBAD/(K_pBAD^n_pBAD + araC^n_pBAD));
		mCherry = scale_mCherry_doi_10_1038_nature12148*ALPHA_BBa_B0030*pBAD;
		output(i) = mCherry;
	end
end

function output = sepf_doi_10_1038_nature12148_figure16a_1(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0030 = x(1);
	K_arabinose = x(2);
	K_pBAD = x(3);
	alpha_pBAD = x(4);
	alpha_pLlacO_1 = x(5);
	beta_pBAD = x(6);
	n_arabinose = x(7);
	n_pBAD = x(8);
	scale_mCherry_doi_10_1038_nature12148 = x(9);
	for i = 1:length(input);
		pLlacO_1 = 5*alpha_pLlacO_1;
		arabinose = input(i);
		araC = ALPHA_BBa_B0030*pLlacO_1*K_arabinose^n_arabinose/(K_arabinose^n_arabinose + arabinose^n_arabinose);
		pBAD = 15*(beta_pBAD + (alpha_pBAD - beta_pBAD)*K_pBAD^n_pBAD/(K_pBAD^n_pBAD + araC^n_pBAD));
		mCherry = scale_mCherry_doi_10_1038_nature12148*ALPHA_BBa_B0030*pBAD;
		output(i) = mCherry;
	end
end

function error = doi_10_1038_nature12148_figure16a_1_error(x)
	arabinose = [0 0.001 0.003 0.008 0.025 0.075 0.2 0.7 2 6 20 50];
	mCherry = [0.015 0.015 0.015 0.02 0.06 0.48 0.84 0.93 0.96 0.98 0.98 1];
	diff = sepf_doi_10_1038_nature12148_figure16a_1(x, arabinose) - mCherry;
	error = repmat(diff,1,5);
end

function output = simf_doi_10_1038_nature11516_figureS8A_1(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pC = x(22);
	ALPHA_RBS_psicA = x(34);
	K_arabinose = x(23);
	K_pBAD = x(24);
	alpha_pBAD = x(25);
	alpha_pC = x(26);
	beta_pBAD = x(27);
	n_arabinose = x(28);
	n_pBAD = x(29);
	scale_mRFP1_doi_10_1038_nature11516 = x(35);
	for i = 1:length(input);
		pC = 15*alpha_pC;
		arabinose = input(i);
		araC = ALPHA_RBS_pC*pC*K_arabinose^n_arabinose/(K_arabinose^n_arabinose + arabinose^n_arabinose);
		pBAD = 15*(beta_pBAD + (alpha_pBAD - beta_pBAD)*K_pBAD^n_pBAD/(K_pBAD^n_pBAD + araC^n_pBAD));
		mRFP1 = scale_mRFP1_doi_10_1038_nature11516*ALPHA_RBS_psicA*pBAD;
		output(i) = mRFP1;
	end
end

function output = sepf_doi_10_1038_nature11516_figureS8A_1(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pC = x(1);
	ALPHA_RBS_psicA = x(2);
	K_arabinose = x(3);
	K_pBAD = x(4);
	alpha_pBAD = x(5);
	alpha_pC = x(6);
	beta_pBAD = x(7);
	n_arabinose = x(8);
	n_pBAD = x(9);
	scale_mRFP1_doi_10_1038_nature11516 = x(10);
	for i = 1:length(input);
		pC = 15*alpha_pC;
		arabinose = input(i);
		araC = ALPHA_RBS_pC*pC*K_arabinose^n_arabinose/(K_arabinose^n_arabinose + arabinose^n_arabinose);
		pBAD = 15*(beta_pBAD + (alpha_pBAD - beta_pBAD)*K_pBAD^n_pBAD/(K_pBAD^n_pBAD + araC^n_pBAD));
		mRFP1 = scale_mRFP1_doi_10_1038_nature11516*ALPHA_RBS_psicA*pBAD;
		output(i) = mRFP1;
	end
end

function error = doi_10_1038_nature11516_figureS8A_1_error(x)
	arabinose = [0 0 0.002 0.01 0.04 0.16 1 4 25];
	mRFP1 = [0.001 0.001 0.001 0.003 0.045 0.161 0.184 0.185 0.185];
	diff = sepf_doi_10_1038_nature11516_figureS8A_1(x, arabinose) - mRFP1;
	error = repmat(diff,1,5);
end

function output = simf_doi_10_1038_nature11516_figureS8B_2(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_placI = x(2);
	ALPHA_RBS_psicA = x(34);
	K_IPTG = x(3);
	K_pTAC = x(36);
	alpha_pTAC = x(37);
	alpha_placIq = x(38);
	beta_pTAC = x(39);
	n_IPTG = x(8);
	n_pTAC = x(40);
	scale_mRFP1_doi_10_1038_nature11516 = x(35);
	for i = 1:length(input);
		placIq = 15*alpha_placIq;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placIq*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pTAC = 15*(beta_pTAC + (alpha_pTAC - beta_pTAC)*K_pTAC^n_pTAC/(K_pTAC^n_pTAC + lacI^n_pTAC));
		mRFP1 = scale_mRFP1_doi_10_1038_nature11516*ALPHA_RBS_psicA*pTAC;
		output(i) = mRFP1;
	end
end

function output = sepf_doi_10_1038_nature11516_figureS8B_2(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_placI = x(1);
	ALPHA_RBS_psicA = x(2);
	K_IPTG = x(3);
	K_pTAC = x(4);
	alpha_pTAC = x(5);
	alpha_placIq = x(6);
	beta_pTAC = x(7);
	n_IPTG = x(8);
	n_pTAC = x(9);
	scale_mRFP1_doi_10_1038_nature11516 = x(10);
	for i = 1:length(input);
		placIq = 15*alpha_placIq;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placIq*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pTAC = 15*(beta_pTAC + (alpha_pTAC - beta_pTAC)*K_pTAC^n_pTAC/(K_pTAC^n_pTAC + lacI^n_pTAC));
		mRFP1 = scale_mRFP1_doi_10_1038_nature11516*ALPHA_RBS_psicA*pTAC;
		output(i) = mRFP1;
	end
end

function error = doi_10_1038_nature11516_figureS8B_2_error(x)
	IPTG = [0 0 0 0.001 0.004 0.025 0.1 0.5 2.5];
	mRFP1 = [0.005 0.005 0.005 0.006 0.013 0.09 0.173 0.194 0.196];
	diff = sepf_doi_10_1038_nature11516_figureS8B_2(x, IPTG) - mRFP1;
	error = repmat(diff,1,5);
end

function output = simf_doi_10_1038_nbt_2401_figureS11_1(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_A = x(41);
	K_IPTG = x(3);
	K_pTAC = x(36);
	alpha_BBa_J23101 = x(42);
	alpha_pTAC = x(37);
	beta_pTAC = x(39);
	n_IPTG = x(8);
	n_pTAC = x(40);
	scale_sfGFP_doi_10_1038_nbt_2401 = x(43);
	for i = 1:length(input);
		BBa_J23101 = 5*alpha_BBa_J23101;
		IPTG = input(i);
		lacI = ALPHA_RBS_A*BBa_J23101*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pTAC = 5*(beta_pTAC + (alpha_pTAC - beta_pTAC)*K_pTAC^n_pTAC/(K_pTAC^n_pTAC + lacI^n_pTAC));
		sfGFP = scale_sfGFP_doi_10_1038_nbt_2401*ALPHA_RBS_A*pTAC;
		output(i) = sfGFP;
	end
end

function output = sepf_doi_10_1038_nbt_2401_figureS11_1(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_A = x(1);
	K_IPTG = x(2);
	K_pTAC = x(3);
	alpha_BBa_J23101 = x(4);
	alpha_pTAC = x(5);
	beta_pTAC = x(6);
	n_IPTG = x(7);
	n_pTAC = x(8);
	scale_sfGFP_doi_10_1038_nbt_2401 = x(9);
	for i = 1:length(input);
		BBa_J23101 = 5*alpha_BBa_J23101;
		IPTG = input(i);
		lacI = ALPHA_RBS_A*BBa_J23101*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pTAC = 5*(beta_pTAC + (alpha_pTAC - beta_pTAC)*K_pTAC^n_pTAC/(K_pTAC^n_pTAC + lacI^n_pTAC));
		sfGFP = scale_sfGFP_doi_10_1038_nbt_2401*ALPHA_RBS_A*pTAC;
		output(i) = sfGFP;
	end
end

function error = doi_10_1038_nbt_2401_figureS11_1_error(x)
	IPTG = [0 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1];
	sfGFP = [0.025 0.025 0.035 0.055 0.115 0.2 0.3 0.4 0.65 1];
	diff = sepf_doi_10_1038_nbt_2401_figureS11_1(x, IPTG) - sfGFP;
	error = repmat(diff,1,5);
end

function output = simf_doi_10_1038_nbt_2401_figureS13_2(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_A = x(41);
	K_arabinose = x(23);
	K_pBAD = x(24);
	alpha_BBa_J23101 = x(42);
	alpha_pBAD = x(25);
	beta_pBAD = x(27);
	n_arabinose = x(28);
	n_pBAD = x(29);
	scale_sfGFP_doi_10_1038_nbt_2401 = x(43);
	for i = 1:length(input);
		BBa_J23101 = 5*alpha_BBa_J23101;
		arabinose = input(i);
		araC = ALPHA_RBS_A*BBa_J23101*K_arabinose^n_arabinose/(K_arabinose^n_arabinose + arabinose^n_arabinose);
		pBAD = 5*(beta_pBAD + (alpha_pBAD - beta_pBAD)*K_pBAD^n_pBAD/(K_pBAD^n_pBAD + araC^n_pBAD));
		sfGFP = scale_sfGFP_doi_10_1038_nbt_2401*ALPHA_RBS_A*pBAD;
		output(i) = sfGFP;
	end
end

function output = sepf_doi_10_1038_nbt_2401_figureS13_2(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_A = x(1);
	K_arabinose = x(2);
	K_pBAD = x(3);
	alpha_BBa_J23101 = x(4);
	alpha_pBAD = x(5);
	beta_pBAD = x(6);
	n_arabinose = x(7);
	n_pBAD = x(8);
	scale_sfGFP_doi_10_1038_nbt_2401 = x(9);
	for i = 1:length(input);
		BBa_J23101 = 5*alpha_BBa_J23101;
		arabinose = input(i);
		araC = ALPHA_RBS_A*BBa_J23101*K_arabinose^n_arabinose/(K_arabinose^n_arabinose + arabinose^n_arabinose);
		pBAD = 5*(beta_pBAD + (alpha_pBAD - beta_pBAD)*K_pBAD^n_pBAD/(K_pBAD^n_pBAD + araC^n_pBAD));
		sfGFP = scale_sfGFP_doi_10_1038_nbt_2401*ALPHA_RBS_A*pBAD;
		output(i) = sfGFP;
	end
end

function error = doi_10_1038_nbt_2401_figureS13_2_error(x)
	arabinose = [0.1 1 2 5 7 10 12.5 25 37.5 50 62.5];
	sfGFP = [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425];
	diff = sepf_doi_10_1038_nbt_2401_figureS13_2(x, arabinose) - sfGFP;
	error = repmat(diff,1,5);
end

function output = simf_doi_10_1038_nature09565_figure1C_1(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0033 = x(44);
	ALPHA_BBa_B0034 = x(45);
	K_arabinose = x(23);
	K_pBAD = x(24);
	alpha_BBa_J23117 = x(46);
	alpha_pBAD = x(25);
	beta_pBAD = x(27);
	n_arabinose = x(28);
	n_pBAD = x(29);
	scale_eYFP_doi_10_1038_nature09565 = x(47);
	for i = 1:length(input);
		BBa_J23117 = 15*alpha_BBa_J23117;
		arabinose = input(i);
		araC = ALPHA_BBa_B0034*BBa_J23117*K_arabinose^n_arabinose/(K_arabinose^n_arabinose + arabinose^n_arabinose);
		pBAD = 15*(beta_pBAD + (alpha_pBAD - beta_pBAD)*K_pBAD^n_pBAD/(K_pBAD^n_pBAD + araC^n_pBAD));
		eYFP = scale_eYFP_doi_10_1038_nature09565*ALPHA_BBa_B0033*pBAD;
		output(i) = eYFP;
	end
end

function output = sepf_doi_10_1038_nature09565_figure1C_1(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0033 = x(1);
	ALPHA_BBa_B0034 = x(2);
	K_arabinose = x(3);
	K_pBAD = x(4);
	alpha_BBa_J23117 = x(5);
	alpha_pBAD = x(6);
	beta_pBAD = x(7);
	n_arabinose = x(8);
	n_pBAD = x(9);
	scale_eYFP_doi_10_1038_nature09565 = x(10);
	for i = 1:length(input);
		BBa_J23117 = 15*alpha_BBa_J23117;
		arabinose = input(i);
		araC = ALPHA_BBa_B0034*BBa_J23117*K_arabinose^n_arabinose/(K_arabinose^n_arabinose + arabinose^n_arabinose);
		pBAD = 15*(beta_pBAD + (alpha_pBAD - beta_pBAD)*K_pBAD^n_pBAD/(K_pBAD^n_pBAD + araC^n_pBAD));
		eYFP = scale_eYFP_doi_10_1038_nature09565*ALPHA_BBa_B0033*pBAD;
		output(i) = eYFP;
	end
end

function error = doi_10_1038_nature09565_figure1C_1_error(x)
	arabinose = [0 0.001 0.015 0.05 0.5 5 10];
	eYFP = [0.002 0.002 0.057 0.628 0.999 1 1];
	diff = sepf_doi_10_1038_nature09565_figure1C_1(x, arabinose) - eYFP;
	error = repmat(diff,1,5);
end

function output = simf_doi_10_1038_nature09565_figureS1C_2(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0032 = x(48);
	ALPHA_BBa_B0033 = x(44);
	K_aTc = x(13);
	K_pLtetO_1 = x(14);
	alpha_BBa_J23117 = x(46);
	alpha_pLtetO_1 = x(15);
	beta_pLtetO_1 = x(17);
	n_aTc = x(18);
	n_pLtetO_1 = x(19);
	scale_eYFP_doi_10_1038_nature09565 = x(47);
	for i = 1:length(input);
		BBa_J23117 = 15*alpha_BBa_J23117;
		aTc = input(i);
		tetR = ALPHA_BBa_B0032*BBa_J23117*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pLtetO_1 = 15*(beta_pLtetO_1 + (alpha_pLtetO_1 - beta_pLtetO_1)*K_pLtetO_1^n_pLtetO_1/(K_pLtetO_1^n_pLtetO_1 + tetR^n_pLtetO_1));
		eYFP = scale_eYFP_doi_10_1038_nature09565*ALPHA_BBa_B0033*pLtetO_1;
		output(i) = eYFP;
	end
end

function output = sepf_doi_10_1038_nature09565_figureS1C_2(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0032 = x(1);
	ALPHA_BBa_B0033 = x(2);
	K_aTc = x(3);
	K_pLtetO_1 = x(4);
	alpha_BBa_J23117 = x(5);
	alpha_pLtetO_1 = x(6);
	beta_pLtetO_1 = x(7);
	n_aTc = x(8);
	n_pLtetO_1 = x(9);
	scale_eYFP_doi_10_1038_nature09565 = x(10);
	for i = 1:length(input);
		BBa_J23117 = 15*alpha_BBa_J23117;
		aTc = input(i);
		tetR = ALPHA_BBa_B0032*BBa_J23117*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pLtetO_1 = 15*(beta_pLtetO_1 + (alpha_pLtetO_1 - beta_pLtetO_1)*K_pLtetO_1^n_pLtetO_1/(K_pLtetO_1^n_pLtetO_1 + tetR^n_pLtetO_1));
		eYFP = scale_eYFP_doi_10_1038_nature09565*ALPHA_BBa_B0033*pLtetO_1;
		output(i) = eYFP;
	end
end

function error = doi_10_1038_nature09565_figureS1C_2_error(x)
	aTc = [0 0.025 0.25 2.5 25 100 250];
	eYFP = [0.005 0.006 0.007 0.03 0.297 0.388 0.398];
	diff = sepf_doi_10_1038_nature09565_figureS1C_2(x, aTc) - eYFP;
	error = repmat(diff,1,5);
end

function output = simf_doi_10_1038_nbt_2149_Supplement_page17_1(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pET_29b = x(49);
	ALPHA_RBS_placI = x(2);
	K_IPTG = x(3);
	K_placUV5 = x(50);
	alpha_placI = x(6);
	alpha_placUV5 = x(51);
	beta_placUV5 = x(52);
	n_IPTG = x(8);
	n_placUV5 = x(53);
	scale_mRFP1_doi_10_1038_nbt_2149 = x(54);
	for i = 1:length(input);
		placI = 10*alpha_placI;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placI*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		placUV5 = 10*(beta_placUV5 + (alpha_placUV5 - beta_placUV5)*K_placUV5^n_placUV5/(K_placUV5^n_placUV5 + lacI^n_placUV5));
		mRFP1 = scale_mRFP1_doi_10_1038_nbt_2149*ALPHA_RBS_pET_29b*placUV5;
		output(i) = mRFP1;
	end
end

function output = sepf_doi_10_1038_nbt_2149_Supplement_page17_1(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_pET_29b = x(1);
	ALPHA_RBS_placI = x(2);
	K_IPTG = x(3);
	K_placUV5 = x(4);
	alpha_placI = x(5);
	alpha_placUV5 = x(6);
	beta_placUV5 = x(7);
	n_IPTG = x(8);
	n_placUV5 = x(9);
	scale_mRFP1_doi_10_1038_nbt_2149 = x(10);
	for i = 1:length(input);
		placI = 10*alpha_placI;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placI*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		placUV5 = 10*(beta_placUV5 + (alpha_placUV5 - beta_placUV5)*K_placUV5^n_placUV5/(K_placUV5^n_placUV5 + lacI^n_placUV5));
		mRFP1 = scale_mRFP1_doi_10_1038_nbt_2149*ALPHA_RBS_pET_29b*placUV5;
		output(i) = mRFP1;
	end
end

function error = doi_10_1038_nbt_2149_Supplement_page17_1_error(x)
	IPTG = [0 0.001 0.005 0.01 0.025 0.05 0.1 0.25 0.5];
	mRFP1 = [0.02 0.031 0.092 0.306 0.745 0.898 0.969 1 1];
	diff = sepf_doi_10_1038_nbt_2149_Supplement_page17_1(x, IPTG) - mRFP1;
	error = repmat(diff,1,5);
end

function output = simf_doi_10_1186_1754_1611_3_4_figure3_1(x)
	ALPHA_BBa_B0032 = x(48);
	alpha_BBa_J23101 = x(42);
	scale_GFPmut3b_doi_10_1186_1754_1611_3_4 = x(55);
	BBa_J23101 = 10*alpha_BBa_J23101;
	GFPmut3b = scale_GFPmut3b_doi_10_1186_1754_1611_3_4*ALPHA_BBa_B0032*BBa_J23101;
	output = GFPmut3b;
end

function output = sepf_doi_10_1186_1754_1611_3_4_figure3_1(x)
	ALPHA_BBa_B0032 = x(1);
	alpha_BBa_J23101 = x(2);
	scale_GFPmut3b_doi_10_1186_1754_1611_3_4 = x(3);
	BBa_J23101 = 10*alpha_BBa_J23101;
	GFPmut3b = scale_GFPmut3b_doi_10_1186_1754_1611_3_4*ALPHA_BBa_B0032*BBa_J23101;
	output = GFPmut3b;
end

function error = doi_10_1186_1754_1611_3_4_figure3_1_error(x)
	diff = sepf_doi_10_1186_1754_1611_3_4_figure3_1(x) - 0.69;
	error = repmat(diff,1,5);
end

function output = simf_doi_10_1186_1754_1611_3_4_figure3_2(x)
	ALPHA_BBa_B0032 = x(48);
	alpha_BBa_J23116 = x(56);
	scale_GFPmut3b_doi_10_1186_1754_1611_3_4 = x(55);
	BBa_J23116 = 10*alpha_BBa_J23116;
	GFPmut3b = scale_GFPmut3b_doi_10_1186_1754_1611_3_4*ALPHA_BBa_B0032*BBa_J23116;
	output = GFPmut3b;
end

function output = sepf_doi_10_1186_1754_1611_3_4_figure3_2(x)
	ALPHA_BBa_B0032 = x(1);
	alpha_BBa_J23116 = x(2);
	scale_GFPmut3b_doi_10_1186_1754_1611_3_4 = x(3);
	BBa_J23116 = 10*alpha_BBa_J23116;
	GFPmut3b = scale_GFPmut3b_doi_10_1186_1754_1611_3_4*ALPHA_BBa_B0032*BBa_J23116;
	output = GFPmut3b;
end

function error = doi_10_1186_1754_1611_3_4_figure3_2_error(x)
	diff = sepf_doi_10_1186_1754_1611_3_4_figure3_2(x) - 0.018;
	error = repmat(diff,1,5);
end

function output = simf_doi_10_1186_1754_1611_3_4_figure3_3(x)
	ALPHA_BBa_B0032 = x(48);
	alpha_BBa_J23150 = x(57);
	scale_GFPmut3b_doi_10_1186_1754_1611_3_4 = x(55);
	BBa_J23150 = 10*alpha_BBa_J23150;
	GFPmut3b = scale_GFPmut3b_doi_10_1186_1754_1611_3_4*ALPHA_BBa_B0032*BBa_J23150;
	output = GFPmut3b;
end

function output = sepf_doi_10_1186_1754_1611_3_4_figure3_3(x)
	ALPHA_BBa_B0032 = x(1);
	alpha_BBa_J23150 = x(2);
	scale_GFPmut3b_doi_10_1186_1754_1611_3_4 = x(3);
	BBa_J23150 = 10*alpha_BBa_J23150;
	GFPmut3b = scale_GFPmut3b_doi_10_1186_1754_1611_3_4*ALPHA_BBa_B0032*BBa_J23150;
	output = GFPmut3b;
end

function error = doi_10_1186_1754_1611_3_4_figure3_3_error(x)
	diff = sepf_doi_10_1186_1754_1611_3_4_figure3_3(x) - 0.193;
	error = repmat(diff,1,5);
end

function output = simf_doi_10_1186_1754_1611_3_4_figure3_4(x)
	ALPHA_BBa_B0032 = x(48);
	alpha_BBa_J23151 = x(58);
	scale_GFPmut3b_doi_10_1186_1754_1611_3_4 = x(55);
	BBa_J23151 = 10*alpha_BBa_J23151;
	GFPmut3b = scale_GFPmut3b_doi_10_1186_1754_1611_3_4*ALPHA_BBa_B0032*BBa_J23151;
	output = GFPmut3b;
end

function output = sepf_doi_10_1186_1754_1611_3_4_figure3_4(x)
	ALPHA_BBa_B0032 = x(1);
	alpha_BBa_J23151 = x(2);
	scale_GFPmut3b_doi_10_1186_1754_1611_3_4 = x(3);
	BBa_J23151 = 10*alpha_BBa_J23151;
	GFPmut3b = scale_GFPmut3b_doi_10_1186_1754_1611_3_4*ALPHA_BBa_B0032*BBa_J23151;
	output = GFPmut3b;
end

function error = doi_10_1186_1754_1611_3_4_figure3_4_error(x)
	diff = sepf_doi_10_1186_1754_1611_3_4_figure3_4(x) - 0.4;
	error = repmat(diff,1,5);
end

function output = simf_doi_10_1186_1754_1611_3_4_figure3_5(x)
	ALPHA_BBa_B0032 = x(48);
	alpha_BBa_J23102 = x(59);
	scale_GFPmut3b_doi_10_1186_1754_1611_3_4 = x(55);
	BBa_J23102 = 10*alpha_BBa_J23102;
	GFPmut3b = scale_GFPmut3b_doi_10_1186_1754_1611_3_4*ALPHA_BBa_B0032*BBa_J23102;
	output = GFPmut3b;
end

function output = sepf_doi_10_1186_1754_1611_3_4_figure3_5(x)
	ALPHA_BBa_B0032 = x(1);
	alpha_BBa_J23102 = x(2);
	scale_GFPmut3b_doi_10_1186_1754_1611_3_4 = x(3);
	BBa_J23102 = 10*alpha_BBa_J23102;
	GFPmut3b = scale_GFPmut3b_doi_10_1186_1754_1611_3_4*ALPHA_BBa_B0032*BBa_J23102;
	output = GFPmut3b;
end

function error = doi_10_1186_1754_1611_3_4_figure3_5_error(x)
	diff = sepf_doi_10_1186_1754_1611_3_4_figure3_5(x) - 0.655;
	error = repmat(diff,1,5);
end

function output = simf_doi_10_1186_1754_1611_3_4_figure3_6(x)
	ALPHA_BBa_B0032 = x(48);
	alpha_pLtetO_1 = x(15);
	scale_GFPmut3b_doi_10_1186_1754_1611_3_4 = x(55);
	pLtetO_1 = 10*alpha_pLtetO_1;
	GFPmut3b = scale_GFPmut3b_doi_10_1186_1754_1611_3_4*ALPHA_BBa_B0032*pLtetO_1;
	output = GFPmut3b;
end

function output = sepf_doi_10_1186_1754_1611_3_4_figure3_6(x)
	ALPHA_BBa_B0032 = x(1);
	alpha_pLtetO_1 = x(2);
	scale_GFPmut3b_doi_10_1186_1754_1611_3_4 = x(3);
	pLtetO_1 = 10*alpha_pLtetO_1;
	GFPmut3b = scale_GFPmut3b_doi_10_1186_1754_1611_3_4*ALPHA_BBa_B0032*pLtetO_1;
	output = GFPmut3b;
end

function error = doi_10_1186_1754_1611_3_4_figure3_6_error(x)
	diff = sepf_doi_10_1186_1754_1611_3_4_figure3_6(x) - 0.897;
	error = repmat(diff,1,5);
end

function output = simf_doi_10_1186_1754_1611_3_4_figure3_7(x)
	ALPHA_BBa_B0032 = x(48);
	alpha_pLlacO_1 = x(32);
	scale_GFPmut3b_doi_10_1186_1754_1611_3_4 = x(55);
	pLlacO_1 = 10*alpha_pLlacO_1;
	GFPmut3b = scale_GFPmut3b_doi_10_1186_1754_1611_3_4*ALPHA_BBa_B0032*pLlacO_1;
	output = GFPmut3b;
end

function output = sepf_doi_10_1186_1754_1611_3_4_figure3_7(x)
	ALPHA_BBa_B0032 = x(1);
	alpha_pLlacO_1 = x(2);
	scale_GFPmut3b_doi_10_1186_1754_1611_3_4 = x(3);
	pLlacO_1 = 10*alpha_pLlacO_1;
	GFPmut3b = scale_GFPmut3b_doi_10_1186_1754_1611_3_4*ALPHA_BBa_B0032*pLlacO_1;
	output = GFPmut3b;
end

function error = doi_10_1186_1754_1611_3_4_figure3_7_error(x)
	diff = sepf_doi_10_1186_1754_1611_3_4_figure3_7(x) - 1;
	error = repmat(diff,1,5);
end

function output = simf_doi__10_1073_pnas_1301301110_sd03_xls_1(x)
	ALPHA_BBa_B0032 = x(48);
	alpha_BBa_J23117 = x(46);
	scale_sfGFP_doi__10_1073_pnas_1301301110 = x(60);
	BBa_J23117 = 10*alpha_BBa_J23117;
	sfGFP = scale_sfGFP_doi__10_1073_pnas_1301301110*ALPHA_BBa_B0032*BBa_J23117;
	output = sfGFP;
end

function output = sepf_doi__10_1073_pnas_1301301110_sd03_xls_1(x)
	ALPHA_BBa_B0032 = x(1);
	alpha_BBa_J23117 = x(2);
	scale_sfGFP_doi__10_1073_pnas_1301301110 = x(3);
	BBa_J23117 = 10*alpha_BBa_J23117;
	sfGFP = scale_sfGFP_doi__10_1073_pnas_1301301110*ALPHA_BBa_B0032*BBa_J23117;
	output = sfGFP;
end

function error = doi__10_1073_pnas_1301301110_sd03_xls_1_error(x)
	diff = sepf_doi__10_1073_pnas_1301301110_sd03_xls_1(x) - 0.422;
	error = repmat(diff,1,5);
end

function output = simf_doi__10_1073_pnas_1301301110_sd03_xls_2(x)
	ALPHA_BBa_B0032 = x(48);
	alpha_BBa_J23101 = x(42);
	scale_sfGFP_doi__10_1073_pnas_1301301110 = x(60);
	BBa_J23101 = 10*alpha_BBa_J23101;
	sfGFP = scale_sfGFP_doi__10_1073_pnas_1301301110*ALPHA_BBa_B0032*BBa_J23101;
	output = sfGFP;
end

function output = sepf_doi__10_1073_pnas_1301301110_sd03_xls_2(x)
	ALPHA_BBa_B0032 = x(1);
	alpha_BBa_J23101 = x(2);
	scale_sfGFP_doi__10_1073_pnas_1301301110 = x(3);
	BBa_J23101 = 10*alpha_BBa_J23101;
	sfGFP = scale_sfGFP_doi__10_1073_pnas_1301301110*ALPHA_BBa_B0032*BBa_J23101;
	output = sfGFP;
end

function error = doi__10_1073_pnas_1301301110_sd03_xls_2_error(x)
	diff = sepf_doi__10_1073_pnas_1301301110_sd03_xls_2(x) - 1;
	error = repmat(diff,1,5);
end

function output = simf_doi_10_1093_nar_gku1388_figure2A_1(x)
	ALPHA_BBa_B0030 = x(31);
	alpha_BBa_J23101 = x(42);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(61);
	BBa_J23101 = 10*alpha_BBa_J23101;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*BBa_J23101;
	output = GFPmut3b;
end

function output = sepf_doi_10_1093_nar_gku1388_figure2A_1(x)
	ALPHA_BBa_B0030 = x(1);
	alpha_BBa_J23101 = x(2);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(3);
	BBa_J23101 = 10*alpha_BBa_J23101;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*BBa_J23101;
	output = GFPmut3b;
end

function error = doi_10_1093_nar_gku1388_figure2A_1_error(x)
	diff = sepf_doi_10_1093_nar_gku1388_figure2A_1(x) - 1;
	error = repmat(diff,1,5);
end

function output = simf_doi_10_1093_nar_gku1388_figure2A_2(x)
	ALPHA_BBa_B0030 = x(31);
	alpha_BBa_J23106 = x(62);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(61);
	BBa_J23106 = 10*alpha_BBa_J23106;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*BBa_J23106;
	output = GFPmut3b;
end

function output = sepf_doi_10_1093_nar_gku1388_figure2A_2(x)
	ALPHA_BBa_B0030 = x(1);
	alpha_BBa_J23106 = x(2);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(3);
	BBa_J23106 = 10*alpha_BBa_J23106;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*BBa_J23106;
	output = GFPmut3b;
end

function error = doi_10_1093_nar_gku1388_figure2A_2_error(x)
	diff = sepf_doi_10_1093_nar_gku1388_figure2A_2(x) - 0.215;
	error = repmat(diff,1,5);
end

function output = simf_doi_10_1093_nar_gku1388_figure2A_3(x)
	ALPHA_BBa_B0030 = x(31);
	alpha_BBa_J23105 = x(63);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(61);
	BBa_J23105 = 10*alpha_BBa_J23105;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*BBa_J23105;
	output = GFPmut3b;
end

function output = sepf_doi_10_1093_nar_gku1388_figure2A_3(x)
	ALPHA_BBa_B0030 = x(1);
	alpha_BBa_J23105 = x(2);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(3);
	BBa_J23105 = 10*alpha_BBa_J23105;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*BBa_J23105;
	output = GFPmut3b;
end

function error = doi_10_1093_nar_gku1388_figure2A_3_error(x)
	diff = sepf_doi_10_1093_nar_gku1388_figure2A_3(x) - 0.095;
	error = repmat(diff,1,5);
end

function output = simf_doi_10_1093_nar_gku1388_figure2A_4(x)
	ALPHA_BBa_B0030 = x(31);
	alpha_BBa_J23115 = x(64);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(61);
	BBa_J23115 = 10*alpha_BBa_J23115;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*BBa_J23115;
	output = GFPmut3b;
end

function output = sepf_doi_10_1093_nar_gku1388_figure2A_4(x)
	ALPHA_BBa_B0030 = x(1);
	alpha_BBa_J23115 = x(2);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(3);
	BBa_J23115 = 10*alpha_BBa_J23115;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*BBa_J23115;
	output = GFPmut3b;
end

function error = doi_10_1093_nar_gku1388_figure2A_4_error(x)
	diff = sepf_doi_10_1093_nar_gku1388_figure2A_4(x) - 0.046;
	error = repmat(diff,1,5);
end

function output = simf_doi_10_1093_nar_gku1388_figure2A_5(x)
	ALPHA_BBa_B0030 = x(31);
	alpha_BBa_J23114 = x(65);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(61);
	BBa_J23114 = 10*alpha_BBa_J23114;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*BBa_J23114;
	output = GFPmut3b;
end

function output = sepf_doi_10_1093_nar_gku1388_figure2A_5(x)
	ALPHA_BBa_B0030 = x(1);
	alpha_BBa_J23114 = x(2);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(3);
	BBa_J23114 = 10*alpha_BBa_J23114;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*BBa_J23114;
	output = GFPmut3b;
end

function error = doi_10_1093_nar_gku1388_figure2A_5_error(x)
	diff = sepf_doi_10_1093_nar_gku1388_figure2A_5(x) - 0.018;
	error = repmat(diff,1,5);
end

function output = simf_doi_10_1093_nar_gku1388_figure2A_6(x)
	ALPHA_BBa_B0030 = x(31);
	alpha_BBa_J23117 = x(46);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(61);
	BBa_J23117 = 10*alpha_BBa_J23117;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*BBa_J23117;
	output = GFPmut3b;
end

function output = sepf_doi_10_1093_nar_gku1388_figure2A_6(x)
	ALPHA_BBa_B0030 = x(1);
	alpha_BBa_J23117 = x(2);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(3);
	BBa_J23117 = 10*alpha_BBa_J23117;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*BBa_J23117;
	output = GFPmut3b;
end

function error = doi_10_1093_nar_gku1388_figure2A_6_error(x)
	diff = sepf_doi_10_1093_nar_gku1388_figure2A_6(x) - 0.003;
	error = repmat(diff,1,5);
end

function output = simf_doi_10_1093_nar_gku1388_figure2B_7(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0030 = x(31);
	K_aTc = x(13);
	K_pTET_star_ = x(66);
	alpha_BBa_J23117 = x(46);
	alpha_pTET_star_ = x(67);
	beta_pTET_star_ = x(68);
	n_aTc = x(18);
	n_pTET_star_ = x(69);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(61);
	for i = 1:length(input);
		BBa_J23117 = 10*alpha_BBa_J23117;
		aTc = input(i);
		tetR = ALPHA_BBa_B0030*BBa_J23117*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pTET_star_ = 10*(beta_pTET_star_ + (alpha_pTET_star_ - beta_pTET_star_)*K_pTET_star_^n_pTET_star_/(K_pTET_star_^n_pTET_star_ + tetR^n_pTET_star_));
		GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*pTET_star_;
		output(i) = GFPmut3b;
	end
end

function output = sepf_doi_10_1093_nar_gku1388_figure2B_7(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0030 = x(1);
	K_aTc = x(2);
	K_pTET_star_ = x(3);
	alpha_BBa_J23117 = x(4);
	alpha_pTET_star_ = x(5);
	beta_pTET_star_ = x(6);
	n_aTc = x(7);
	n_pTET_star_ = x(8);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(9);
	for i = 1:length(input);
		BBa_J23117 = 10*alpha_BBa_J23117;
		aTc = input(i);
		tetR = ALPHA_BBa_B0030*BBa_J23117*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pTET_star_ = 10*(beta_pTET_star_ + (alpha_pTET_star_ - beta_pTET_star_)*K_pTET_star_^n_pTET_star_/(K_pTET_star_^n_pTET_star_ + tetR^n_pTET_star_));
		GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*pTET_star_;
		output(i) = GFPmut3b;
	end
end

function error = doi_10_1093_nar_gku1388_figure2B_7_error(x)
	aTc = [0 2 4 10 25 50 100];
	GFPmut3b = [0 0.002 0.014 0.167 0.666 0.817 0.841];
	diff = sepf_doi_10_1093_nar_gku1388_figure2B_7(x, aTc) - GFPmut3b;
	error = repmat(diff,1,5);
end

function output = simf_doi_10_1093_nar_gku1388_figure2B_8(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0030 = x(31);
	K_aTc = x(13);
	K_pTET_star_ = x(66);
	alpha_BBa_J23114 = x(65);
	alpha_pTET_star_ = x(67);
	beta_pTET_star_ = x(68);
	n_aTc = x(18);
	n_pTET_star_ = x(69);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(61);
	for i = 1:length(input);
		BBa_J23114 = 10*alpha_BBa_J23114;
		aTc = input(i);
		tetR = ALPHA_BBa_B0030*BBa_J23114*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pTET_star_ = 10*(beta_pTET_star_ + (alpha_pTET_star_ - beta_pTET_star_)*K_pTET_star_^n_pTET_star_/(K_pTET_star_^n_pTET_star_ + tetR^n_pTET_star_));
		GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*pTET_star_;
		output(i) = GFPmut3b;
	end
end

function output = sepf_doi_10_1093_nar_gku1388_figure2B_8(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0030 = x(1);
	K_aTc = x(2);
	K_pTET_star_ = x(3);
	alpha_BBa_J23114 = x(4);
	alpha_pTET_star_ = x(5);
	beta_pTET_star_ = x(6);
	n_aTc = x(7);
	n_pTET_star_ = x(8);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(9);
	for i = 1:length(input);
		BBa_J23114 = 10*alpha_BBa_J23114;
		aTc = input(i);
		tetR = ALPHA_BBa_B0030*BBa_J23114*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pTET_star_ = 10*(beta_pTET_star_ + (alpha_pTET_star_ - beta_pTET_star_)*K_pTET_star_^n_pTET_star_/(K_pTET_star_^n_pTET_star_ + tetR^n_pTET_star_));
		GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*pTET_star_;
		output(i) = GFPmut3b;
	end
end

function error = doi_10_1093_nar_gku1388_figure2B_8_error(x)
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0 0 0 0.001 0.022 0.176 0.411 0.461];
	diff = sepf_doi_10_1093_nar_gku1388_figure2B_8(x, aTc) - GFPmut3b;
	error = repmat(diff,1,5);
end

function output = simf_doi_10_1093_nar_gku1388_figure2B_9(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0030 = x(31);
	K_aTc = x(13);
	K_pTET_star_ = x(66);
	alpha_BBa_J23115 = x(64);
	alpha_pTET_star_ = x(67);
	beta_pTET_star_ = x(68);
	n_aTc = x(18);
	n_pTET_star_ = x(69);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(61);
	for i = 1:length(input);
		BBa_J23115 = 10*alpha_BBa_J23115;
		aTc = input(i);
		tetR = ALPHA_BBa_B0030*BBa_J23115*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pTET_star_ = 10*(beta_pTET_star_ + (alpha_pTET_star_ - beta_pTET_star_)*K_pTET_star_^n_pTET_star_/(K_pTET_star_^n_pTET_star_ + tetR^n_pTET_star_));
		GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*pTET_star_;
		output(i) = GFPmut3b;
	end
end

function output = sepf_doi_10_1093_nar_gku1388_figure2B_9(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0030 = x(1);
	K_aTc = x(2);
	K_pTET_star_ = x(3);
	alpha_BBa_J23115 = x(4);
	alpha_pTET_star_ = x(5);
	beta_pTET_star_ = x(6);
	n_aTc = x(7);
	n_pTET_star_ = x(8);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(9);
	for i = 1:length(input);
		BBa_J23115 = 10*alpha_BBa_J23115;
		aTc = input(i);
		tetR = ALPHA_BBa_B0030*BBa_J23115*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pTET_star_ = 10*(beta_pTET_star_ + (alpha_pTET_star_ - beta_pTET_star_)*K_pTET_star_^n_pTET_star_/(K_pTET_star_^n_pTET_star_ + tetR^n_pTET_star_));
		GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*pTET_star_;
		output(i) = GFPmut3b;
	end
end

function error = doi_10_1093_nar_gku1388_figure2B_9_error(x)
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.005 0.118 0.277];
	diff = sepf_doi_10_1093_nar_gku1388_figure2B_9(x, aTc) - GFPmut3b;
	error = repmat(diff,1,5);
end

function output = simf_doi_10_1093_nar_gku1388_figure2B_10(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0030 = x(31);
	K_aTc = x(13);
	K_pTET_star_ = x(66);
	alpha_BBa_J23105 = x(63);
	alpha_pTET_star_ = x(67);
	beta_pTET_star_ = x(68);
	n_aTc = x(18);
	n_pTET_star_ = x(69);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(61);
	for i = 1:length(input);
		BBa_J23105 = 10*alpha_BBa_J23105;
		aTc = input(i);
		tetR = ALPHA_BBa_B0030*BBa_J23105*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pTET_star_ = 10*(beta_pTET_star_ + (alpha_pTET_star_ - beta_pTET_star_)*K_pTET_star_^n_pTET_star_/(K_pTET_star_^n_pTET_star_ + tetR^n_pTET_star_));
		GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*pTET_star_;
		output(i) = GFPmut3b;
	end
end

function output = sepf_doi_10_1093_nar_gku1388_figure2B_10(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0030 = x(1);
	K_aTc = x(2);
	K_pTET_star_ = x(3);
	alpha_BBa_J23105 = x(4);
	alpha_pTET_star_ = x(5);
	beta_pTET_star_ = x(6);
	n_aTc = x(7);
	n_pTET_star_ = x(8);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(9);
	for i = 1:length(input);
		BBa_J23105 = 10*alpha_BBa_J23105;
		aTc = input(i);
		tetR = ALPHA_BBa_B0030*BBa_J23105*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pTET_star_ = 10*(beta_pTET_star_ + (alpha_pTET_star_ - beta_pTET_star_)*K_pTET_star_^n_pTET_star_/(K_pTET_star_^n_pTET_star_ + tetR^n_pTET_star_));
		GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*pTET_star_;
		output(i) = GFPmut3b;
	end
end

function error = doi_10_1093_nar_gku1388_figure2B_10_error(x)
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.001 0.001 0.001 0.001 0.001 0.004 0.039 0.147];
	diff = sepf_doi_10_1093_nar_gku1388_figure2B_10(x, aTc) - GFPmut3b;
	error = repmat(diff,1,5);
end

function output = simf_doi_10_1093_nar_gku1388_figure2B_11(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0030 = x(31);
	K_aTc = x(13);
	K_pTET_star_ = x(66);
	alpha_BBa_J23106 = x(62);
	alpha_pTET_star_ = x(67);
	beta_pTET_star_ = x(68);
	n_aTc = x(18);
	n_pTET_star_ = x(69);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(61);
	for i = 1:length(input);
		BBa_J23106 = 10*alpha_BBa_J23106;
		aTc = input(i);
		tetR = ALPHA_BBa_B0030*BBa_J23106*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pTET_star_ = 10*(beta_pTET_star_ + (alpha_pTET_star_ - beta_pTET_star_)*K_pTET_star_^n_pTET_star_/(K_pTET_star_^n_pTET_star_ + tetR^n_pTET_star_));
		GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*pTET_star_;
		output(i) = GFPmut3b;
	end
end

function output = sepf_doi_10_1093_nar_gku1388_figure2B_11(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0030 = x(1);
	K_aTc = x(2);
	K_pTET_star_ = x(3);
	alpha_BBa_J23106 = x(4);
	alpha_pTET_star_ = x(5);
	beta_pTET_star_ = x(6);
	n_aTc = x(7);
	n_pTET_star_ = x(8);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(9);
	for i = 1:length(input);
		BBa_J23106 = 10*alpha_BBa_J23106;
		aTc = input(i);
		tetR = ALPHA_BBa_B0030*BBa_J23106*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pTET_star_ = 10*(beta_pTET_star_ + (alpha_pTET_star_ - beta_pTET_star_)*K_pTET_star_^n_pTET_star_/(K_pTET_star_^n_pTET_star_ + tetR^n_pTET_star_));
		GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*pTET_star_;
		output(i) = GFPmut3b;
	end
end

function error = doi_10_1093_nar_gku1388_figure2B_11_error(x)
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.003 0.004 0.032];
	diff = sepf_doi_10_1093_nar_gku1388_figure2B_11(x, aTc) - GFPmut3b;
	error = repmat(diff,1,5);
end

function output = simf_doi_10_1093_nar_gku1388_figure2B_12(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0030 = x(31);
	K_aTc = x(13);
	K_pTET_star_ = x(66);
	alpha_BBa_J23101 = x(42);
	alpha_pTET_star_ = x(67);
	beta_pTET_star_ = x(68);
	n_aTc = x(18);
	n_pTET_star_ = x(69);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(61);
	for i = 1:length(input);
		BBa_J23101 = 10*alpha_BBa_J23101;
		aTc = input(i);
		tetR = ALPHA_BBa_B0030*BBa_J23101*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pTET_star_ = 10*(beta_pTET_star_ + (alpha_pTET_star_ - beta_pTET_star_)*K_pTET_star_^n_pTET_star_/(K_pTET_star_^n_pTET_star_ + tetR^n_pTET_star_));
		GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*pTET_star_;
		output(i) = GFPmut3b;
	end
end

function output = sepf_doi_10_1093_nar_gku1388_figure2B_12(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0030 = x(1);
	K_aTc = x(2);
	K_pTET_star_ = x(3);
	alpha_BBa_J23101 = x(4);
	alpha_pTET_star_ = x(5);
	beta_pTET_star_ = x(6);
	n_aTc = x(7);
	n_pTET_star_ = x(8);
	scale_GFPmut3b_doi_10_1093_nar_gku1388 = x(9);
	for i = 1:length(input);
		BBa_J23101 = 10*alpha_BBa_J23101;
		aTc = input(i);
		tetR = ALPHA_BBa_B0030*BBa_J23101*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pTET_star_ = 10*(beta_pTET_star_ + (alpha_pTET_star_ - beta_pTET_star_)*K_pTET_star_^n_pTET_star_/(K_pTET_star_^n_pTET_star_ + tetR^n_pTET_star_));
		GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku1388*ALPHA_BBa_B0030*pTET_star_;
		output(i) = GFPmut3b;
	end
end

function error = doi_10_1093_nar_gku1388_figure2B_12_error(x)
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.006 0.006 0.006 0.006 0.006 0.006 0.009 0.021];
	diff = sepf_doi_10_1093_nar_gku1388_figure2B_12(x, aTc) - GFPmut3b;
	error = repmat(diff,1,5);
end

function output = simf_doi_10_1093_nar_gku593_figure3D_1_1(x)
	ALPHA_BBa_B0030 = x(31);
	alpha_BBa_J23109 = x(70);
	scale_GFPmut3b_doi_10_1093_nar_gku593 = x(71);
	BBa_J23109 = 10*alpha_BBa_J23109;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku593*ALPHA_BBa_B0030*BBa_J23109;
	output = GFPmut3b;
end

function output = sepf_doi_10_1093_nar_gku593_figure3D_1_1(x)
	ALPHA_BBa_B0030 = x(1);
	alpha_BBa_J23109 = x(2);
	scale_GFPmut3b_doi_10_1093_nar_gku593 = x(3);
	BBa_J23109 = 10*alpha_BBa_J23109;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku593*ALPHA_BBa_B0030*BBa_J23109;
	output = GFPmut3b;
end

function error = doi_10_1093_nar_gku593_figure3D_1_1_error(x)
	diff = sepf_doi_10_1093_nar_gku593_figure3D_1_1(x) - 0.002;
	error = repmat(diff,1,5);
end

function output = simf_doi_10_1093_nar_gku593_figure3D_2_2(x)
	ALPHA_BBa_B0030 = x(31);
	alpha_BBa_J23114 = x(65);
	scale_GFPmut3b_doi_10_1093_nar_gku593 = x(71);
	BBa_J23114 = 10*alpha_BBa_J23114;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku593*ALPHA_BBa_B0030*BBa_J23114;
	output = GFPmut3b;
end

function output = sepf_doi_10_1093_nar_gku593_figure3D_2_2(x)
	ALPHA_BBa_B0030 = x(1);
	alpha_BBa_J23114 = x(2);
	scale_GFPmut3b_doi_10_1093_nar_gku593 = x(3);
	BBa_J23114 = 10*alpha_BBa_J23114;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku593*ALPHA_BBa_B0030*BBa_J23114;
	output = GFPmut3b;
end

function error = doi_10_1093_nar_gku593_figure3D_2_2_error(x)
	diff = sepf_doi_10_1093_nar_gku593_figure3D_2_2(x) - 0.018;
	error = repmat(diff,1,5);
end

function output = simf_doi_10_1093_nar_gku593_figure3D_3_3(x)
	ALPHA_BBa_B0030 = x(31);
	alpha_BBa_J23116 = x(56);
	scale_GFPmut3b_doi_10_1093_nar_gku593 = x(71);
	BBa_J23116 = 10*alpha_BBa_J23116;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku593*ALPHA_BBa_B0030*BBa_J23116;
	output = GFPmut3b;
end

function output = sepf_doi_10_1093_nar_gku593_figure3D_3_3(x)
	ALPHA_BBa_B0030 = x(1);
	alpha_BBa_J23116 = x(2);
	scale_GFPmut3b_doi_10_1093_nar_gku593 = x(3);
	BBa_J23116 = 10*alpha_BBa_J23116;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku593*ALPHA_BBa_B0030*BBa_J23116;
	output = GFPmut3b;
end

function error = doi_10_1093_nar_gku593_figure3D_3_3_error(x)
	diff = sepf_doi_10_1093_nar_gku593_figure3D_3_3(x) - 0.036;
	error = repmat(diff,1,5);
end

function output = simf_doi_10_1093_nar_gku593_figure3D_4_4(x)
	ALPHA_BBa_B0030 = x(31);
	alpha_BBa_J23115 = x(64);
	scale_GFPmut3b_doi_10_1093_nar_gku593 = x(71);
	BBa_J23115 = 10*alpha_BBa_J23115;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku593*ALPHA_BBa_B0030*BBa_J23115;
	output = GFPmut3b;
end

function output = sepf_doi_10_1093_nar_gku593_figure3D_4_4(x)
	ALPHA_BBa_B0030 = x(1);
	alpha_BBa_J23115 = x(2);
	scale_GFPmut3b_doi_10_1093_nar_gku593 = x(3);
	BBa_J23115 = 10*alpha_BBa_J23115;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku593*ALPHA_BBa_B0030*BBa_J23115;
	output = GFPmut3b;
end

function error = doi_10_1093_nar_gku593_figure3D_4_4_error(x)
	diff = sepf_doi_10_1093_nar_gku593_figure3D_4_4(x) - 0.05;
	error = repmat(diff,1,5);
end

function output = simf_doi_10_1093_nar_gku593_figure3D_5_5(x)
	ALPHA_BBa_B0030 = x(31);
	alpha_BBa_J23105 = x(63);
	scale_GFPmut3b_doi_10_1093_nar_gku593 = x(71);
	BBa_J23105 = 10*alpha_BBa_J23105;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku593*ALPHA_BBa_B0030*BBa_J23105;
	output = GFPmut3b;
end

function output = sepf_doi_10_1093_nar_gku593_figure3D_5_5(x)
	ALPHA_BBa_B0030 = x(1);
	alpha_BBa_J23105 = x(2);
	scale_GFPmut3b_doi_10_1093_nar_gku593 = x(3);
	BBa_J23105 = 10*alpha_BBa_J23105;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku593*ALPHA_BBa_B0030*BBa_J23105;
	output = GFPmut3b;
end

function error = doi_10_1093_nar_gku593_figure3D_5_5_error(x)
	diff = sepf_doi_10_1093_nar_gku593_figure3D_5_5(x) - 0.086;
	error = repmat(diff,1,5);
end

function output = simf_doi_10_1093_nar_gku593_figure3D_6_6(x)
	ALPHA_BBa_B0030 = x(31);
	alpha_BBa_J23106 = x(62);
	scale_GFPmut3b_doi_10_1093_nar_gku593 = x(71);
	BBa_J23106 = 10*alpha_BBa_J23106;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku593*ALPHA_BBa_B0030*BBa_J23106;
	output = GFPmut3b;
end

function output = sepf_doi_10_1093_nar_gku593_figure3D_6_6(x)
	ALPHA_BBa_B0030 = x(1);
	alpha_BBa_J23106 = x(2);
	scale_GFPmut3b_doi_10_1093_nar_gku593 = x(3);
	BBa_J23106 = 10*alpha_BBa_J23106;
	GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku593*ALPHA_BBa_B0030*BBa_J23106;
	output = GFPmut3b;
end

function error = doi_10_1093_nar_gku593_figure3D_6_6_error(x)
	diff = sepf_doi_10_1093_nar_gku593_figure3D_6_6(x) - 0.214;
	error = repmat(diff,1,5);
end

function output = simf_doi_10_1093_nar_gku593_figureS6_7(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0030 = x(31);
	ALPHA_RBS_pC = x(22);
	K_arabinose = x(23);
	K_pBAD = x(24);
	alpha_pBAD = x(25);
	alpha_pC = x(26);
	beta_pBAD = x(27);
	n_arabinose = x(28);
	n_pBAD = x(29);
	scale_GFPmut3b_doi_10_1093_nar_gku593 = x(71);
	for i = 1:length(input);
		pC = 10*alpha_pC;
		arabinose = input(i);
		araC = ALPHA_RBS_pC*pC*K_arabinose^n_arabinose/(K_arabinose^n_arabinose + arabinose^n_arabinose);
		pBAD = 10*(beta_pBAD + (alpha_pBAD - beta_pBAD)*K_pBAD^n_pBAD/(K_pBAD^n_pBAD + araC^n_pBAD));
		GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku593*ALPHA_BBa_B0030*pBAD;
		output(i) = GFPmut3b;
	end
end

function output = sepf_doi_10_1093_nar_gku593_figureS6_7(x, input)
	output = zeros(1,length(input));
	ALPHA_BBa_B0030 = x(1);
	ALPHA_RBS_pC = x(2);
	K_arabinose = x(3);
	K_pBAD = x(4);
	alpha_pBAD = x(5);
	alpha_pC = x(6);
	beta_pBAD = x(7);
	n_arabinose = x(8);
	n_pBAD = x(9);
	scale_GFPmut3b_doi_10_1093_nar_gku593 = x(10);
	for i = 1:length(input);
		pC = 10*alpha_pC;
		arabinose = input(i);
		araC = ALPHA_RBS_pC*pC*K_arabinose^n_arabinose/(K_arabinose^n_arabinose + arabinose^n_arabinose);
		pBAD = 10*(beta_pBAD + (alpha_pBAD - beta_pBAD)*K_pBAD^n_pBAD/(K_pBAD^n_pBAD + araC^n_pBAD));
		GFPmut3b = scale_GFPmut3b_doi_10_1093_nar_gku593*ALPHA_BBa_B0030*pBAD;
		output(i) = GFPmut3b;
	end
end

function error = doi_10_1093_nar_gku593_figureS6_7_error(x)
	arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2];
	GFPmut3b = [0 0.002 0.041 0.479 0.949 0.997 1 1 1];
	diff = sepf_doi_10_1093_nar_gku593_figureS6_7(x, arabinose) - GFPmut3b;
	error = repmat(diff,1,5);
end

function output = simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_placI = x(2);
	ALPHA_RBSII = x(12);
	K_IPTG = x(3);
	K_pLAC = x(4);
	alpha_pLAC = x(5);
	alpha_placIq = x(38);
	beta_pLAC = x(7);
	n_IPTG = x(8);
	n_pLAC = x(9);
	scale_eYFP_http___dspace_mit_edu_handle_1721_1_8228 = x(72);
	for i = 1:length(input);
		placIq = 10*alpha_placIq;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placIq*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pLAC = 10*(beta_pLAC + (alpha_pLAC - beta_pLAC)*K_pLAC^n_pLAC/(K_pLAC^n_pLAC + lacI^n_pLAC));
		eYFP = scale_eYFP_http___dspace_mit_edu_handle_1721_1_8228*ALPHA_RBSII*pLAC;
		output(i) = eYFP;
	end
end

function output = sepf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_placI = x(1);
	ALPHA_RBSII = x(2);
	K_IPTG = x(3);
	K_pLAC = x(4);
	alpha_pLAC = x(5);
	alpha_placIq = x(6);
	beta_pLAC = x(7);
	n_IPTG = x(8);
	n_pLAC = x(9);
	scale_eYFP_http___dspace_mit_edu_handle_1721_1_8228 = x(10);
	for i = 1:length(input);
		placIq = 10*alpha_placIq;
		IPTG = input(i);
		lacI = ALPHA_RBS_placI*placIq*K_IPTG^n_IPTG/(K_IPTG^n_IPTG + IPTG^n_IPTG);
		pLAC = 10*(beta_pLAC + (alpha_pLAC - beta_pLAC)*K_pLAC^n_pLAC/(K_pLAC^n_pLAC + lacI^n_pLAC));
		eYFP = scale_eYFP_http___dspace_mit_edu_handle_1721_1_8228*ALPHA_RBSII*pLAC;
		output(i) = eYFP;
	end
end

function error = http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1_error(x)
	IPTG = [0.001 0.01 0.03 0.05 0.075 0.1 0.15 0.2 0.3 0.5 1];
	eYFP = [0.015 0.014 0.023 0.062 0.098 0.187 0.239 0.374 0.572 0.869 1];
	diff = sepf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1(x, IPTG) - eYFP;
	error = repmat(diff,1,5);
end

function output = simf_doi_10_1073pnas_0408507102_figure2A_1(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_Unknown_Hooshangi = x(73);
	ALPHA_RBS_placI = x(2);
	K_aTc = x(13);
	K_pLtetO_1 = x(14);
	alpha_pLtetO_1 = x(15);
	alpha_placIq = x(38);
	beta_pLtetO_1 = x(17);
	n_aTc = x(18);
	n_pLtetO_1 = x(19);
	scale_eYFP_doi_10_1073pnas_0408507102 = x(74);
	for i = 1:length(input);
		placIq = 10*alpha_placIq;
		aTc = input(i);
		tetR = ALPHA_RBS_placI*placIq*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pLtetO_1 = 10*(beta_pLtetO_1 + (alpha_pLtetO_1 - beta_pLtetO_1)*K_pLtetO_1^n_pLtetO_1/(K_pLtetO_1^n_pLtetO_1 + tetR^n_pLtetO_1));
		eYFP = scale_eYFP_doi_10_1073pnas_0408507102*ALPHA_RBS_Unknown_Hooshangi*pLtetO_1;
		output(i) = eYFP;
	end
end

function output = sepf_doi_10_1073pnas_0408507102_figure2A_1(x, input)
	output = zeros(1,length(input));
	ALPHA_RBS_Unknown_Hooshangi = x(1);
	ALPHA_RBS_placI = x(2);
	K_aTc = x(3);
	K_pLtetO_1 = x(4);
	alpha_pLtetO_1 = x(5);
	alpha_placIq = x(6);
	beta_pLtetO_1 = x(7);
	n_aTc = x(8);
	n_pLtetO_1 = x(9);
	scale_eYFP_doi_10_1073pnas_0408507102 = x(10);
	for i = 1:length(input);
		placIq = 10*alpha_placIq;
		aTc = input(i);
		tetR = ALPHA_RBS_placI*placIq*K_aTc^n_aTc/(K_aTc^n_aTc + aTc^n_aTc);
		pLtetO_1 = 10*(beta_pLtetO_1 + (alpha_pLtetO_1 - beta_pLtetO_1)*K_pLtetO_1^n_pLtetO_1/(K_pLtetO_1^n_pLtetO_1 + tetR^n_pLtetO_1));
		eYFP = scale_eYFP_doi_10_1073pnas_0408507102*ALPHA_RBS_Unknown_Hooshangi*pLtetO_1;
		output(i) = eYFP;
	end
end

function error = doi_10_1073pnas_0408507102_figure2A_1_error(x)
	aTc = [9.258 13.887 20.831 32.403 46.29 83.322 162.015 231.45 370.32 648.06 925.8 1388.7];
	eYFP = [0.003 0.003 0.003 0.003 0.005 0.01 0.05 0.167 0.333 0.667 0.833 1];
	diff = sepf_doi_10_1073pnas_0408507102_figure2A_1(x, aTc) - eYFP;
	error = repmat(diff,1,5);
end

function print_plotting_data(seq_x_opt, sim_x_opt)
	%%%%%%%%%%%%%
	IPTG = [0 0.005 0.01 0.02 0.05 0.1 0.2 1];
	lacZ = [0.001 0.003 0.043 0.5 0.984 0.999 1 1];
	input_log_list = logspace(-4,0.699,1000);
	point_matrix = [IPTG; lacZ];
	point_fileID = fopen('Plot_data/doi_10_1073_pnas_0606717104_figure4_6_1_point.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_matrix);
	fclose(point_fileID);
	output_seq = simf_doi_10_1073_pnas_0606717104_figure4_6_1(seq_x_opt, input_log_list);
	curve_matrix = [input_log_list; output_seq];
	curve_fileID = fopen('Plot_data/doi_10_1073_pnas_0606717104_figure4_6_1_seq_curve.dat','w');
	fprintf(curve_fileID, '%f \t %f\n', curve_matrix);
	fclose(curve_fileID);
	output_sim = simf_doi_10_1073_pnas_0606717104_figure4_6_1(sim_x_opt, input_log_list);
	curve_matrix = [input_log_list; output_sim];
	curve_fileID = fopen('Plot_data/doi_10_1073_pnas_0606717104_figure4_6_1_sim_curve.dat','w');
	fprintf(curve_fileID, '%f \t %f\n', curve_matrix);
	fclose(curve_fileID);
	%%%%%%%%%%%%%
	aTc = [0 1 2 3 4 5 6 7 8 9 10 12 13 17 18 20 25 35 60];
	luc = [0 0.001 0.001 0.001 0.002 0.004 0.007 0.011 0.014 0.032 0.036 0.054 0.089 0.643 0.893 1 1 1 1];
	input_log_list = logspace(-4,3.444,1000);
	point_matrix = [aTc; luc];
	point_fileID = fopen('Plot_data/doi__10_1093_nar_25_6_1203_figure4a_1_point.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_matrix);
	fclose(point_fileID);
	output_seq = simf_doi__10_1093_nar_25_6_1203_figure4a_1(seq_x_opt, input_log_list);
	curve_matrix = [input_log_list; output_seq];
	curve_fileID = fopen('Plot_data/doi__10_1093_nar_25_6_1203_figure4a_1_seq_curve.dat','w');
	fprintf(curve_fileID, '%f \t %f\n', curve_matrix);
	fclose(curve_fileID);
	output_sim = simf_doi__10_1093_nar_25_6_1203_figure4a_1(sim_x_opt, input_log_list);
	curve_matrix = [input_log_list; output_sim];
	curve_fileID = fopen('Plot_data/doi__10_1093_nar_25_6_1203_figure4a_1_sim_curve.dat','w');
	fprintf(curve_fileID, '%f \t %f\n', curve_matrix);
	fclose(curve_fileID);
	%%%%%%%%%%%%%
	arabinose = [0.003 0.005 0.05 0.5 1];
	GFPuv = [0.044 0.089 0.556 1 1];
	input_log_list = logspace(-4,2.097,1000);
	point_matrix = [arabinose; GFPuv];
	point_fileID = fopen('Plot_data/doi_10_1128_AEM_00791_07_figure3a_1_point.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_matrix);
	fclose(point_fileID);
	output_seq = simf_doi_10_1128_AEM_00791_07_figure3a_1(seq_x_opt, input_log_list);
	curve_matrix = [input_log_list; output_seq];
	curve_fileID = fopen('Plot_data/doi_10_1128_AEM_00791_07_figure3a_1_seq_curve.dat','w');
	fprintf(curve_fileID, '%f \t %f\n', curve_matrix);
	fclose(curve_fileID);
	output_sim = simf_doi_10_1128_AEM_00791_07_figure3a_1(sim_x_opt, input_log_list);
	curve_matrix = [input_log_list; output_sim];
	curve_fileID = fopen('Plot_data/doi_10_1128_AEM_00791_07_figure3a_1_sim_curve.dat','w');
	fprintf(curve_fileID, '%f \t %f\n', curve_matrix);
	fclose(curve_fileID);
	%%%%%%%%%%%%%
	arabinose = [0 0.001 0.003 0.008 0.025 0.075 0.2 0.7 2 6 20 50];
	mCherry = [0.015 0.015 0.015 0.02 0.06 0.48 0.84 0.93 0.96 0.98 0.98 1];
	input_log_list = logspace(-4,2.097,1000);
	point_matrix = [arabinose; mCherry];
	point_fileID = fopen('Plot_data/doi_10_1038_nature12148_figure16a_1_point.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_matrix);
	fclose(point_fileID);
	output_seq = simf_doi_10_1038_nature12148_figure16a_1(seq_x_opt, input_log_list);
	curve_matrix = [input_log_list; output_seq];
	curve_fileID = fopen('Plot_data/doi_10_1038_nature12148_figure16a_1_seq_curve.dat','w');
	fprintf(curve_fileID, '%f \t %f\n', curve_matrix);
	fclose(curve_fileID);
	output_sim = simf_doi_10_1038_nature12148_figure16a_1(sim_x_opt, input_log_list);
	curve_matrix = [input_log_list; output_sim];
	curve_fileID = fopen('Plot_data/doi_10_1038_nature12148_figure16a_1_sim_curve.dat','w');
	fprintf(curve_fileID, '%f \t %f\n', curve_matrix);
	fclose(curve_fileID);
	%%%%%%%%%%%%%
	arabinose = [0 0 0.002 0.01 0.04 0.16 1 4 25];
	mRFP1 = [0.001 0.001 0.001 0.003 0.045 0.161 0.184 0.185 0.185];
	input_log_list = logspace(-4,2.097,1000);
	point_matrix = [arabinose; mRFP1];
	point_fileID = fopen('Plot_data/doi_10_1038_nature11516_figureS8A_1_point.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_matrix);
	fclose(point_fileID);
	output_seq = simf_doi_10_1038_nature11516_figureS8A_1(seq_x_opt, input_log_list);
	curve_matrix = [input_log_list; output_seq];
	curve_fileID = fopen('Plot_data/doi_10_1038_nature11516_figureS8A_1_seq_curve.dat','w');
	fprintf(curve_fileID, '%f \t %f\n', curve_matrix);
	fclose(curve_fileID);
	output_sim = simf_doi_10_1038_nature11516_figureS8A_1(sim_x_opt, input_log_list);
	curve_matrix = [input_log_list; output_sim];
	curve_fileID = fopen('Plot_data/doi_10_1038_nature11516_figureS8A_1_sim_curve.dat','w');
	fprintf(curve_fileID, '%f \t %f\n', curve_matrix);
	fclose(curve_fileID);
	%%%%%%%%%%%%%
	IPTG = [0 0 0 0.001 0.004 0.025 0.1 0.5 2.5];
	mRFP1 = [0.005 0.005 0.005 0.006 0.013 0.09 0.173 0.194 0.196];
	input_log_list = logspace(-4,0.699,1000);
	point_matrix = [IPTG; mRFP1];
	point_fileID = fopen('Plot_data/doi_10_1038_nature11516_figureS8B_2_point.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_matrix);
	fclose(point_fileID);
	output_seq = simf_doi_10_1038_nature11516_figureS8B_2(seq_x_opt, input_log_list);
	curve_matrix = [input_log_list; output_seq];
	curve_fileID = fopen('Plot_data/doi_10_1038_nature11516_figureS8B_2_seq_curve.dat','w');
	fprintf(curve_fileID, '%f \t %f\n', curve_matrix);
	fclose(curve_fileID);
	output_sim = simf_doi_10_1038_nature11516_figureS8B_2(sim_x_opt, input_log_list);
	curve_matrix = [input_log_list; output_sim];
	curve_fileID = fopen('Plot_data/doi_10_1038_nature11516_figureS8B_2_sim_curve.dat','w');
	fprintf(curve_fileID, '%f \t %f\n', curve_matrix);
	fclose(curve_fileID);
	%%%%%%%%%%%%%
	IPTG = [0 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1];
	sfGFP = [0.025 0.025 0.035 0.055 0.115 0.2 0.3 0.4 0.65 1];
	input_log_list = logspace(-4,0.699,1000);
	point_matrix = [IPTG; sfGFP];
	point_fileID = fopen('Plot_data/doi_10_1038_nbt_2401_figureS11_1_point.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_matrix);
	fclose(point_fileID);
	output_seq = simf_doi_10_1038_nbt_2401_figureS11_1(seq_x_opt, input_log_list);
	curve_matrix = [input_log_list; output_seq];
	curve_fileID = fopen('Plot_data/doi_10_1038_nbt_2401_figureS11_1_seq_curve.dat','w');
	fprintf(curve_fileID, '%f \t %f\n', curve_matrix);
	fclose(curve_fileID);
	output_sim = simf_doi_10_1038_nbt_2401_figureS11_1(sim_x_opt, input_log_list);
	curve_matrix = [input_log_list; output_sim];
	curve_fileID = fopen('Plot_data/doi_10_1038_nbt_2401_figureS11_1_sim_curve.dat','w');
	fprintf(curve_fileID, '%f \t %f\n', curve_matrix);
	fclose(curve_fileID);
	%%%%%%%%%%%%%
	arabinose = [0.1 1 2 5 7 10 12.5 25 37.5 50 62.5];
	sfGFP = [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425];
	input_log_list = logspace(-4,2.097,1000);
	point_matrix = [arabinose; sfGFP];
	point_fileID = fopen('Plot_data/doi_10_1038_nbt_2401_figureS13_2_point.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_matrix);
	fclose(point_fileID);
	output_seq = simf_doi_10_1038_nbt_2401_figureS13_2(seq_x_opt, input_log_list);
	curve_matrix = [input_log_list; output_seq];
	curve_fileID = fopen('Plot_data/doi_10_1038_nbt_2401_figureS13_2_seq_curve.dat','w');
	fprintf(curve_fileID, '%f \t %f\n', curve_matrix);
	fclose(curve_fileID);
	output_sim = simf_doi_10_1038_nbt_2401_figureS13_2(sim_x_opt, input_log_list);
	curve_matrix = [input_log_list; output_sim];
	curve_fileID = fopen('Plot_data/doi_10_1038_nbt_2401_figureS13_2_sim_curve.dat','w');
	fprintf(curve_fileID, '%f \t %f\n', curve_matrix);
	fclose(curve_fileID);
	%%%%%%%%%%%%%
	arabinose = [0 0.001 0.015 0.05 0.5 5 10];
	eYFP = [0.002 0.002 0.057 0.628 0.999 1 1];
	input_log_list = logspace(-4,2.097,1000);
	point_matrix = [arabinose; eYFP];
	point_fileID = fopen('Plot_data/doi_10_1038_nature09565_figure1C_1_point.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_matrix);
	fclose(point_fileID);
	output_seq = simf_doi_10_1038_nature09565_figure1C_1(seq_x_opt, input_log_list);
	curve_matrix = [input_log_list; output_seq];
	curve_fileID = fopen('Plot_data/doi_10_1038_nature09565_figure1C_1_seq_curve.dat','w');
	fprintf(curve_fileID, '%f \t %f\n', curve_matrix);
	fclose(curve_fileID);
	output_sim = simf_doi_10_1038_nature09565_figure1C_1(sim_x_opt, input_log_list);
	curve_matrix = [input_log_list; output_sim];
	curve_fileID = fopen('Plot_data/doi_10_1038_nature09565_figure1C_1_sim_curve.dat','w');
	fprintf(curve_fileID, '%f \t %f\n', curve_matrix);
	fclose(curve_fileID);
	%%%%%%%%%%%%%
	aTc = [0 0.025 0.25 2.5 25 100 250];
	eYFP = [0.005 0.006 0.007 0.03 0.297 0.388 0.398];
	input_log_list = logspace(-4,3.444,1000);
	point_matrix = [aTc; eYFP];
	point_fileID = fopen('Plot_data/doi_10_1038_nature09565_figureS1C_2_point.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_matrix);
	fclose(point_fileID);
	output_seq = simf_doi_10_1038_nature09565_figureS1C_2(seq_x_opt, input_log_list);
	curve_matrix = [input_log_list; output_seq];
	curve_fileID = fopen('Plot_data/doi_10_1038_nature09565_figureS1C_2_seq_curve.dat','w');
	fprintf(curve_fileID, '%f \t %f\n', curve_matrix);
	fclose(curve_fileID);
	output_sim = simf_doi_10_1038_nature09565_figureS1C_2(sim_x_opt, input_log_list);
	curve_matrix = [input_log_list; output_sim];
	curve_fileID = fopen('Plot_data/doi_10_1038_nature09565_figureS1C_2_sim_curve.dat','w');
	fprintf(curve_fileID, '%f \t %f\n', curve_matrix);
	fclose(curve_fileID);
	%%%%%%%%%%%%%
	IPTG = [0 0.001 0.005 0.01 0.025 0.05 0.1 0.25 0.5];
	mRFP1 = [0.02 0.031 0.092 0.306 0.745 0.898 0.969 1 1];
	input_log_list = logspace(-4,0.699,1000);
	point_matrix = [IPTG; mRFP1];
	point_fileID = fopen('Plot_data/doi_10_1038_nbt_2149_Supplement_page17_1_point.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_matrix);
	fclose(point_fileID);
	output_seq = simf_doi_10_1038_nbt_2149_Supplement_page17_1(seq_x_opt, input_log_list);
	curve_matrix = [input_log_list; output_seq];
	curve_fileID = fopen('Plot_data/doi_10_1038_nbt_2149_Supplement_page17_1_seq_curve.dat','w');
	fprintf(curve_fileID, '%f \t %f\n', curve_matrix);
	fclose(curve_fileID);
	output_sim = simf_doi_10_1038_nbt_2149_Supplement_page17_1(sim_x_opt, input_log_list);
	curve_matrix = [input_log_list; output_sim];
	curve_fileID = fopen('Plot_data/doi_10_1038_nbt_2149_Supplement_page17_1_sim_curve.dat','w');
	fprintf(curve_fileID, '%f \t %f\n', curve_matrix);
	fclose(curve_fileID);
	%%%%%%%%%%%%%
	output_seq = simf_doi_10_1186_1754_1611_3_4_figure3_1(seq_x_opt);
	point_pair = [output_seq 	0.69];
	point_fileID = fopen('Plot_data/doi_10_1186_1754_1611_3_4_figure3_1_seq_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	output_sim = simf_doi_10_1186_1754_1611_3_4_figure3_1(sim_x_opt);
	point_pair = [output_sim 	0.69];
	point_fileID = fopen('Plot_data/doi_10_1186_1754_1611_3_4_figure3_1_sim_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	%%%%%%%%%%%%%
	output_seq = simf_doi_10_1186_1754_1611_3_4_figure3_2(seq_x_opt);
	point_pair = [output_seq 	0.018];
	point_fileID = fopen('Plot_data/doi_10_1186_1754_1611_3_4_figure3_2_seq_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	output_sim = simf_doi_10_1186_1754_1611_3_4_figure3_2(sim_x_opt);
	point_pair = [output_sim 	0.018];
	point_fileID = fopen('Plot_data/doi_10_1186_1754_1611_3_4_figure3_2_sim_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	%%%%%%%%%%%%%
	output_seq = simf_doi_10_1186_1754_1611_3_4_figure3_3(seq_x_opt);
	point_pair = [output_seq 	0.193];
	point_fileID = fopen('Plot_data/doi_10_1186_1754_1611_3_4_figure3_3_seq_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	output_sim = simf_doi_10_1186_1754_1611_3_4_figure3_3(sim_x_opt);
	point_pair = [output_sim 	0.193];
	point_fileID = fopen('Plot_data/doi_10_1186_1754_1611_3_4_figure3_3_sim_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	%%%%%%%%%%%%%
	output_seq = simf_doi_10_1186_1754_1611_3_4_figure3_4(seq_x_opt);
	point_pair = [output_seq 	0.4];
	point_fileID = fopen('Plot_data/doi_10_1186_1754_1611_3_4_figure3_4_seq_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	output_sim = simf_doi_10_1186_1754_1611_3_4_figure3_4(sim_x_opt);
	point_pair = [output_sim 	0.4];
	point_fileID = fopen('Plot_data/doi_10_1186_1754_1611_3_4_figure3_4_sim_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	%%%%%%%%%%%%%
	output_seq = simf_doi_10_1186_1754_1611_3_4_figure3_5(seq_x_opt);
	point_pair = [output_seq 	0.655];
	point_fileID = fopen('Plot_data/doi_10_1186_1754_1611_3_4_figure3_5_seq_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	output_sim = simf_doi_10_1186_1754_1611_3_4_figure3_5(sim_x_opt);
	point_pair = [output_sim 	0.655];
	point_fileID = fopen('Plot_data/doi_10_1186_1754_1611_3_4_figure3_5_sim_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	%%%%%%%%%%%%%
	output_seq = simf_doi_10_1186_1754_1611_3_4_figure3_6(seq_x_opt);
	point_pair = [output_seq 	0.897];
	point_fileID = fopen('Plot_data/doi_10_1186_1754_1611_3_4_figure3_6_seq_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	output_sim = simf_doi_10_1186_1754_1611_3_4_figure3_6(sim_x_opt);
	point_pair = [output_sim 	0.897];
	point_fileID = fopen('Plot_data/doi_10_1186_1754_1611_3_4_figure3_6_sim_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	%%%%%%%%%%%%%
	output_seq = simf_doi_10_1186_1754_1611_3_4_figure3_7(seq_x_opt);
	point_pair = [output_seq 	1];
	point_fileID = fopen('Plot_data/doi_10_1186_1754_1611_3_4_figure3_7_seq_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	output_sim = simf_doi_10_1186_1754_1611_3_4_figure3_7(sim_x_opt);
	point_pair = [output_sim 	1];
	point_fileID = fopen('Plot_data/doi_10_1186_1754_1611_3_4_figure3_7_sim_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	%%%%%%%%%%%%%
	output_seq = simf_doi__10_1073_pnas_1301301110_sd03_xls_1(seq_x_opt);
	point_pair = [output_seq 	0.422];
	point_fileID = fopen('Plot_data/doi__10_1073_pnas_1301301110_sd03_xls_1_seq_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	output_sim = simf_doi__10_1073_pnas_1301301110_sd03_xls_1(sim_x_opt);
	point_pair = [output_sim 	0.422];
	point_fileID = fopen('Plot_data/doi__10_1073_pnas_1301301110_sd03_xls_1_sim_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	%%%%%%%%%%%%%
	output_seq = simf_doi__10_1073_pnas_1301301110_sd03_xls_2(seq_x_opt);
	point_pair = [output_seq 	1];
	point_fileID = fopen('Plot_data/doi__10_1073_pnas_1301301110_sd03_xls_2_seq_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	output_sim = simf_doi__10_1073_pnas_1301301110_sd03_xls_2(sim_x_opt);
	point_pair = [output_sim 	1];
	point_fileID = fopen('Plot_data/doi__10_1073_pnas_1301301110_sd03_xls_2_sim_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	%%%%%%%%%%%%%
	output_seq = simf_doi_10_1093_nar_gku1388_figure2A_1(seq_x_opt);
	point_pair = [output_seq 	1];
	point_fileID = fopen('Plot_data/doi_10_1093_nar_gku1388_figure2A_1_seq_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	output_sim = simf_doi_10_1093_nar_gku1388_figure2A_1(sim_x_opt);
	point_pair = [output_sim 	1];
	point_fileID = fopen('Plot_data/doi_10_1093_nar_gku1388_figure2A_1_sim_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	%%%%%%%%%%%%%
	output_seq = simf_doi_10_1093_nar_gku1388_figure2A_2(seq_x_opt);
	point_pair = [output_seq 	0.215];
	point_fileID = fopen('Plot_data/doi_10_1093_nar_gku1388_figure2A_2_seq_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	output_sim = simf_doi_10_1093_nar_gku1388_figure2A_2(sim_x_opt);
	point_pair = [output_sim 	0.215];
	point_fileID = fopen('Plot_data/doi_10_1093_nar_gku1388_figure2A_2_sim_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	%%%%%%%%%%%%%
	output_seq = simf_doi_10_1093_nar_gku1388_figure2A_3(seq_x_opt);
	point_pair = [output_seq 	0.095];
	point_fileID = fopen('Plot_data/doi_10_1093_nar_gku1388_figure2A_3_seq_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	output_sim = simf_doi_10_1093_nar_gku1388_figure2A_3(sim_x_opt);
	point_pair = [output_sim 	0.095];
	point_fileID = fopen('Plot_data/doi_10_1093_nar_gku1388_figure2A_3_sim_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	%%%%%%%%%%%%%
	output_seq = simf_doi_10_1093_nar_gku1388_figure2A_4(seq_x_opt);
	point_pair = [output_seq 	0.046];
	point_fileID = fopen('Plot_data/doi_10_1093_nar_gku1388_figure2A_4_seq_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	output_sim = simf_doi_10_1093_nar_gku1388_figure2A_4(sim_x_opt);
	point_pair = [output_sim 	0.046];
	point_fileID = fopen('Plot_data/doi_10_1093_nar_gku1388_figure2A_4_sim_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	%%%%%%%%%%%%%
	output_seq = simf_doi_10_1093_nar_gku1388_figure2A_5(seq_x_opt);
	point_pair = [output_seq 	0.018];
	point_fileID = fopen('Plot_data/doi_10_1093_nar_gku1388_figure2A_5_seq_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	output_sim = simf_doi_10_1093_nar_gku1388_figure2A_5(sim_x_opt);
	point_pair = [output_sim 	0.018];
	point_fileID = fopen('Plot_data/doi_10_1093_nar_gku1388_figure2A_5_sim_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	%%%%%%%%%%%%%
	output_seq = simf_doi_10_1093_nar_gku1388_figure2A_6(seq_x_opt);
	point_pair = [output_seq 	0.003];
	point_fileID = fopen('Plot_data/doi_10_1093_nar_gku1388_figure2A_6_seq_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	output_sim = simf_doi_10_1093_nar_gku1388_figure2A_6(sim_x_opt);
	point_pair = [output_sim 	0.003];
	point_fileID = fopen('Plot_data/doi_10_1093_nar_gku1388_figure2A_6_sim_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	%%%%%%%%%%%%%
	aTc = [0 2 4 10 25 50 100];
	GFPmut3b = [0 0.002 0.014 0.167 0.666 0.817 0.841];
	input_log_list = logspace(-4,3.444,1000);
	point_matrix = [aTc; GFPmut3b];
	point_fileID = fopen('Plot_data/doi_10_1093_nar_gku1388_figure2B_7_point.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_matrix);
	fclose(point_fileID);
	output_seq = simf_doi_10_1093_nar_gku1388_figure2B_7(seq_x_opt, input_log_list);
	curve_matrix = [input_log_list; output_seq];
	curve_fileID = fopen('Plot_data/doi_10_1093_nar_gku1388_figure2B_7_seq_curve.dat','w');
	fprintf(curve_fileID, '%f \t %f\n', curve_matrix);
	fclose(curve_fileID);
	output_sim = simf_doi_10_1093_nar_gku1388_figure2B_7(sim_x_opt, input_log_list);
	curve_matrix = [input_log_list; output_sim];
	curve_fileID = fopen('Plot_data/doi_10_1093_nar_gku1388_figure2B_7_sim_curve.dat','w');
	fprintf(curve_fileID, '%f \t %f\n', curve_matrix);
	fclose(curve_fileID);
	%%%%%%%%%%%%%
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0 0 0 0.001 0.022 0.176 0.411 0.461];
	input_log_list = logspace(-4,3.444,1000);
	point_matrix = [aTc; GFPmut3b];
	point_fileID = fopen('Plot_data/doi_10_1093_nar_gku1388_figure2B_8_point.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_matrix);
	fclose(point_fileID);
	output_seq = simf_doi_10_1093_nar_gku1388_figure2B_8(seq_x_opt, input_log_list);
	curve_matrix = [input_log_list; output_seq];
	curve_fileID = fopen('Plot_data/doi_10_1093_nar_gku1388_figure2B_8_seq_curve.dat','w');
	fprintf(curve_fileID, '%f \t %f\n', curve_matrix);
	fclose(curve_fileID);
	output_sim = simf_doi_10_1093_nar_gku1388_figure2B_8(sim_x_opt, input_log_list);
	curve_matrix = [input_log_list; output_sim];
	curve_fileID = fopen('Plot_data/doi_10_1093_nar_gku1388_figure2B_8_sim_curve.dat','w');
	fprintf(curve_fileID, '%f \t %f\n', curve_matrix);
	fclose(curve_fileID);
	%%%%%%%%%%%%%
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.005 0.118 0.277];
	input_log_list = logspace(-4,3.444,1000);
	point_matrix = [aTc; GFPmut3b];
	point_fileID = fopen('Plot_data/doi_10_1093_nar_gku1388_figure2B_9_point.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_matrix);
	fclose(point_fileID);
	output_seq = simf_doi_10_1093_nar_gku1388_figure2B_9(seq_x_opt, input_log_list);
	curve_matrix = [input_log_list; output_seq];
	curve_fileID = fopen('Plot_data/doi_10_1093_nar_gku1388_figure2B_9_seq_curve.dat','w');
	fprintf(curve_fileID, '%f \t %f\n', curve_matrix);
	fclose(curve_fileID);
	output_sim = simf_doi_10_1093_nar_gku1388_figure2B_9(sim_x_opt, input_log_list);
	curve_matrix = [input_log_list; output_sim];
	curve_fileID = fopen('Plot_data/doi_10_1093_nar_gku1388_figure2B_9_sim_curve.dat','w');
	fprintf(curve_fileID, '%f \t %f\n', curve_matrix);
	fclose(curve_fileID);
	%%%%%%%%%%%%%
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.001 0.001 0.001 0.001 0.001 0.004 0.039 0.147];
	input_log_list = logspace(-4,3.444,1000);
	point_matrix = [aTc; GFPmut3b];
	point_fileID = fopen('Plot_data/doi_10_1093_nar_gku1388_figure2B_10_point.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_matrix);
	fclose(point_fileID);
	output_seq = simf_doi_10_1093_nar_gku1388_figure2B_10(seq_x_opt, input_log_list);
	curve_matrix = [input_log_list; output_seq];
	curve_fileID = fopen('Plot_data/doi_10_1093_nar_gku1388_figure2B_10_seq_curve.dat','w');
	fprintf(curve_fileID, '%f \t %f\n', curve_matrix);
	fclose(curve_fileID);
	output_sim = simf_doi_10_1093_nar_gku1388_figure2B_10(sim_x_opt, input_log_list);
	curve_matrix = [input_log_list; output_sim];
	curve_fileID = fopen('Plot_data/doi_10_1093_nar_gku1388_figure2B_10_sim_curve.dat','w');
	fprintf(curve_fileID, '%f \t %f\n', curve_matrix);
	fclose(curve_fileID);
	%%%%%%%%%%%%%
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.003 0.004 0.032];
	input_log_list = logspace(-4,3.444,1000);
	point_matrix = [aTc; GFPmut3b];
	point_fileID = fopen('Plot_data/doi_10_1093_nar_gku1388_figure2B_11_point.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_matrix);
	fclose(point_fileID);
	output_seq = simf_doi_10_1093_nar_gku1388_figure2B_11(seq_x_opt, input_log_list);
	curve_matrix = [input_log_list; output_seq];
	curve_fileID = fopen('Plot_data/doi_10_1093_nar_gku1388_figure2B_11_seq_curve.dat','w');
	fprintf(curve_fileID, '%f \t %f\n', curve_matrix);
	fclose(curve_fileID);
	output_sim = simf_doi_10_1093_nar_gku1388_figure2B_11(sim_x_opt, input_log_list);
	curve_matrix = [input_log_list; output_sim];
	curve_fileID = fopen('Plot_data/doi_10_1093_nar_gku1388_figure2B_11_sim_curve.dat','w');
	fprintf(curve_fileID, '%f \t %f\n', curve_matrix);
	fclose(curve_fileID);
	%%%%%%%%%%%%%
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.006 0.006 0.006 0.006 0.006 0.006 0.009 0.021];
	input_log_list = logspace(-4,3.444,1000);
	point_matrix = [aTc; GFPmut3b];
	point_fileID = fopen('Plot_data/doi_10_1093_nar_gku1388_figure2B_12_point.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_matrix);
	fclose(point_fileID);
	output_seq = simf_doi_10_1093_nar_gku1388_figure2B_12(seq_x_opt, input_log_list);
	curve_matrix = [input_log_list; output_seq];
	curve_fileID = fopen('Plot_data/doi_10_1093_nar_gku1388_figure2B_12_seq_curve.dat','w');
	fprintf(curve_fileID, '%f \t %f\n', curve_matrix);
	fclose(curve_fileID);
	output_sim = simf_doi_10_1093_nar_gku1388_figure2B_12(sim_x_opt, input_log_list);
	curve_matrix = [input_log_list; output_sim];
	curve_fileID = fopen('Plot_data/doi_10_1093_nar_gku1388_figure2B_12_sim_curve.dat','w');
	fprintf(curve_fileID, '%f \t %f\n', curve_matrix);
	fclose(curve_fileID);
	%%%%%%%%%%%%%
	output_seq = simf_doi_10_1093_nar_gku593_figure3D_1_1(seq_x_opt);
	point_pair = [output_seq 	0.002];
	point_fileID = fopen('Plot_data/doi_10_1093_nar_gku593_figure3D_1_1_seq_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	output_sim = simf_doi_10_1093_nar_gku593_figure3D_1_1(sim_x_opt);
	point_pair = [output_sim 	0.002];
	point_fileID = fopen('Plot_data/doi_10_1093_nar_gku593_figure3D_1_1_sim_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	%%%%%%%%%%%%%
	output_seq = simf_doi_10_1093_nar_gku593_figure3D_2_2(seq_x_opt);
	point_pair = [output_seq 	0.018];
	point_fileID = fopen('Plot_data/doi_10_1093_nar_gku593_figure3D_2_2_seq_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	output_sim = simf_doi_10_1093_nar_gku593_figure3D_2_2(sim_x_opt);
	point_pair = [output_sim 	0.018];
	point_fileID = fopen('Plot_data/doi_10_1093_nar_gku593_figure3D_2_2_sim_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	%%%%%%%%%%%%%
	output_seq = simf_doi_10_1093_nar_gku593_figure3D_3_3(seq_x_opt);
	point_pair = [output_seq 	0.036];
	point_fileID = fopen('Plot_data/doi_10_1093_nar_gku593_figure3D_3_3_seq_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	output_sim = simf_doi_10_1093_nar_gku593_figure3D_3_3(sim_x_opt);
	point_pair = [output_sim 	0.036];
	point_fileID = fopen('Plot_data/doi_10_1093_nar_gku593_figure3D_3_3_sim_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	%%%%%%%%%%%%%
	output_seq = simf_doi_10_1093_nar_gku593_figure3D_4_4(seq_x_opt);
	point_pair = [output_seq 	0.05];
	point_fileID = fopen('Plot_data/doi_10_1093_nar_gku593_figure3D_4_4_seq_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	output_sim = simf_doi_10_1093_nar_gku593_figure3D_4_4(sim_x_opt);
	point_pair = [output_sim 	0.05];
	point_fileID = fopen('Plot_data/doi_10_1093_nar_gku593_figure3D_4_4_sim_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	%%%%%%%%%%%%%
	output_seq = simf_doi_10_1093_nar_gku593_figure3D_5_5(seq_x_opt);
	point_pair = [output_seq 	0.086];
	point_fileID = fopen('Plot_data/doi_10_1093_nar_gku593_figure3D_5_5_seq_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	output_sim = simf_doi_10_1093_nar_gku593_figure3D_5_5(sim_x_opt);
	point_pair = [output_sim 	0.086];
	point_fileID = fopen('Plot_data/doi_10_1093_nar_gku593_figure3D_5_5_sim_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	%%%%%%%%%%%%%
	output_seq = simf_doi_10_1093_nar_gku593_figure3D_6_6(seq_x_opt);
	point_pair = [output_seq 	0.214];
	point_fileID = fopen('Plot_data/doi_10_1093_nar_gku593_figure3D_6_6_seq_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	output_sim = simf_doi_10_1093_nar_gku593_figure3D_6_6(sim_x_opt);
	point_pair = [output_sim 	0.214];
	point_fileID = fopen('Plot_data/doi_10_1093_nar_gku593_figure3D_6_6_sim_pair.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_pair);
	fclose(point_fileID);
	%%%%%%%%%%%%%
	arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2];
	GFPmut3b = [0 0.002 0.041 0.479 0.949 0.997 1 1 1];
	input_log_list = logspace(-4,2.097,1000);
	point_matrix = [arabinose; GFPmut3b];
	point_fileID = fopen('Plot_data/doi_10_1093_nar_gku593_figureS6_7_point.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_matrix);
	fclose(point_fileID);
	output_seq = simf_doi_10_1093_nar_gku593_figureS6_7(seq_x_opt, input_log_list);
	curve_matrix = [input_log_list; output_seq];
	curve_fileID = fopen('Plot_data/doi_10_1093_nar_gku593_figureS6_7_seq_curve.dat','w');
	fprintf(curve_fileID, '%f \t %f\n', curve_matrix);
	fclose(curve_fileID);
	output_sim = simf_doi_10_1093_nar_gku593_figureS6_7(sim_x_opt, input_log_list);
	curve_matrix = [input_log_list; output_sim];
	curve_fileID = fopen('Plot_data/doi_10_1093_nar_gku593_figureS6_7_sim_curve.dat','w');
	fprintf(curve_fileID, '%f \t %f\n', curve_matrix);
	fclose(curve_fileID);
	%%%%%%%%%%%%%
	IPTG = [0.001 0.01 0.03 0.05 0.075 0.1 0.15 0.2 0.3 0.5 1];
	eYFP = [0.015 0.014 0.023 0.062 0.098 0.187 0.239 0.374 0.572 0.869 1];
	input_log_list = logspace(-4,0.699,1000);
	point_matrix = [IPTG; eYFP];
	point_fileID = fopen('Plot_data/http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1_point.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_matrix);
	fclose(point_fileID);
	output_seq = simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1(seq_x_opt, input_log_list);
	curve_matrix = [input_log_list; output_seq];
	curve_fileID = fopen('Plot_data/http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1_seq_curve.dat','w');
	fprintf(curve_fileID, '%f \t %f\n', curve_matrix);
	fclose(curve_fileID);
	output_sim = simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1(sim_x_opt, input_log_list);
	curve_matrix = [input_log_list; output_sim];
	curve_fileID = fopen('Plot_data/http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1_sim_curve.dat','w');
	fprintf(curve_fileID, '%f \t %f\n', curve_matrix);
	fclose(curve_fileID);
	%%%%%%%%%%%%%
	aTc = [9.258 13.887 20.831 32.403 46.29 83.322 162.015 231.45 370.32 648.06 925.8 1388.7];
	eYFP = [0.003 0.003 0.003 0.003 0.005 0.01 0.05 0.167 0.333 0.667 0.833 1];
	input_log_list = logspace(-4,3.444,1000);
	point_matrix = [aTc; eYFP];
	point_fileID = fopen('Plot_data/doi_10_1073pnas_0408507102_figure2A_1_point.dat','w');
	fprintf(point_fileID, '%f \t %f\n', point_matrix);
	fclose(point_fileID);
	output_seq = simf_doi_10_1073pnas_0408507102_figure2A_1(seq_x_opt, input_log_list);
	curve_matrix = [input_log_list; output_seq];
	curve_fileID = fopen('Plot_data/doi_10_1073pnas_0408507102_figure2A_1_seq_curve.dat','w');
	fprintf(curve_fileID, '%f \t %f\n', curve_matrix);
	fclose(curve_fileID);
	output_sim = simf_doi_10_1073pnas_0408507102_figure2A_1(sim_x_opt, input_log_list);
	curve_matrix = [input_log_list; output_sim];
	curve_fileID = fopen('Plot_data/doi_10_1073pnas_0408507102_figure2A_1_sim_curve.dat','w');
	fprintf(curve_fileID, '%f \t %f\n', curve_matrix);
	fclose(curve_fileID);
end

function separate_fitting(N)
	options = optimset('MaxIter', 1000, 'MaxFunEvals', 30000);
	combined_CI_table = zeros(41,222);
	%===============
	lb = [0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01];
	ub = [5 5 50 50 50 50 1 4 4 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@doi_10_1073_pnas_0606717104_figure4_6_1_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(1,1) = CI_lb(1);
	combined_CI_table(1,2) = CI_ub(1);
	combined_CI_table(1,3) = x_opt(1);
	combined_CI_table(1,4) = CI_lb(2);
	combined_CI_table(1,5) = CI_ub(2);
	combined_CI_table(1,6) = x_opt(2);
	combined_CI_table(1,7) = CI_lb(3);
	combined_CI_table(1,8) = CI_ub(3);
	combined_CI_table(1,9) = x_opt(3);
	combined_CI_table(1,10) = CI_lb(4);
	combined_CI_table(1,11) = CI_ub(4);
	combined_CI_table(1,12) = x_opt(4);
	combined_CI_table(1,13) = CI_lb(5);
	combined_CI_table(1,14) = CI_ub(5);
	combined_CI_table(1,15) = x_opt(5);
	combined_CI_table(1,16) = CI_lb(6);
	combined_CI_table(1,17) = CI_ub(6);
	combined_CI_table(1,18) = x_opt(6);
	combined_CI_table(1,19) = CI_lb(7);
	combined_CI_table(1,20) = CI_ub(7);
	combined_CI_table(1,21) = x_opt(7);
	combined_CI_table(1,22) = CI_lb(8);
	combined_CI_table(1,23) = CI_ub(8);
	combined_CI_table(1,24) = x_opt(8);
	combined_CI_table(1,25) = CI_lb(9);
	combined_CI_table(1,26) = CI_ub(9);
	combined_CI_table(1,27) = x_opt(9);
	combined_CI_table(1,28) = CI_lb(10);
	combined_CI_table(1,29) = CI_ub(10);
	combined_CI_table(1,30) = x_opt(10);
	%===============
	lb = [0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01];
	ub = [5 5 50 50 50 50 1 4 4 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@doi__10_1093_nar_25_6_1203_figure4a_1_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(2,31) = CI_lb(1);
	combined_CI_table(2,32) = CI_ub(1);
	combined_CI_table(2,33) = x_opt(1);
	combined_CI_table(2,34) = CI_lb(2);
	combined_CI_table(2,35) = CI_ub(2);
	combined_CI_table(2,36) = x_opt(2);
	combined_CI_table(2,37) = CI_lb(3);
	combined_CI_table(2,38) = CI_ub(3);
	combined_CI_table(2,39) = x_opt(3);
	combined_CI_table(2,40) = CI_lb(4);
	combined_CI_table(2,41) = CI_ub(4);
	combined_CI_table(2,42) = x_opt(4);
	combined_CI_table(2,43) = CI_lb(5);
	combined_CI_table(2,44) = CI_ub(5);
	combined_CI_table(2,45) = x_opt(5);
	combined_CI_table(2,46) = CI_lb(6);
	combined_CI_table(2,47) = CI_ub(6);
	combined_CI_table(2,48) = x_opt(6);
	combined_CI_table(2,49) = CI_lb(7);
	combined_CI_table(2,50) = CI_ub(7);
	combined_CI_table(2,51) = x_opt(7);
	combined_CI_table(2,52) = CI_lb(8);
	combined_CI_table(2,53) = CI_ub(8);
	combined_CI_table(2,54) = x_opt(8);
	combined_CI_table(2,55) = CI_lb(9);
	combined_CI_table(2,56) = CI_ub(9);
	combined_CI_table(2,57) = x_opt(9);
	combined_CI_table(2,58) = CI_lb(10);
	combined_CI_table(2,59) = CI_ub(10);
	combined_CI_table(2,60) = x_opt(10);
	%===============
	lb = [0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01];
	ub = [5 5 50 50 50 50 1 4 4 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@doi_10_1128_AEM_00791_07_figure3a_1_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(3,61) = CI_lb(1);
	combined_CI_table(3,62) = CI_ub(1);
	combined_CI_table(3,63) = x_opt(1);
	combined_CI_table(3,64) = CI_lb(2);
	combined_CI_table(3,65) = CI_ub(2);
	combined_CI_table(3,66) = x_opt(2);
	combined_CI_table(3,67) = CI_lb(3);
	combined_CI_table(3,68) = CI_ub(3);
	combined_CI_table(3,69) = x_opt(3);
	combined_CI_table(3,70) = CI_lb(4);
	combined_CI_table(3,71) = CI_ub(4);
	combined_CI_table(3,72) = x_opt(4);
	combined_CI_table(3,73) = CI_lb(5);
	combined_CI_table(3,74) = CI_ub(5);
	combined_CI_table(3,75) = x_opt(5);
	combined_CI_table(3,76) = CI_lb(6);
	combined_CI_table(3,77) = CI_ub(6);
	combined_CI_table(3,78) = x_opt(6);
	combined_CI_table(3,79) = CI_lb(7);
	combined_CI_table(3,80) = CI_ub(7);
	combined_CI_table(3,81) = x_opt(7);
	combined_CI_table(3,82) = CI_lb(8);
	combined_CI_table(3,83) = CI_ub(8);
	combined_CI_table(3,84) = x_opt(8);
	combined_CI_table(3,85) = CI_lb(9);
	combined_CI_table(3,86) = CI_ub(9);
	combined_CI_table(3,87) = x_opt(9);
	combined_CI_table(3,88) = CI_lb(10);
	combined_CI_table(3,89) = CI_ub(10);
	combined_CI_table(3,90) = x_opt(10);
	%===============
	lb = [0.03 0.001 0.001 0.002 0.002 0 1 1 0.01];
	ub = [5 50 50 50 50 1 4 4 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@doi_10_1038_nature12148_figure16a_1_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(4,91) = CI_lb(1);
	combined_CI_table(4,92) = CI_ub(1);
	combined_CI_table(4,93) = x_opt(1);
	combined_CI_table(4,67) = CI_lb(2);
	combined_CI_table(4,68) = CI_ub(2);
	combined_CI_table(4,69) = x_opt(2);
	combined_CI_table(4,70) = CI_lb(3);
	combined_CI_table(4,71) = CI_ub(3);
	combined_CI_table(4,72) = x_opt(3);
	combined_CI_table(4,73) = CI_lb(4);
	combined_CI_table(4,74) = CI_ub(4);
	combined_CI_table(4,75) = x_opt(4);
	combined_CI_table(4,94) = CI_lb(5);
	combined_CI_table(4,95) = CI_ub(5);
	combined_CI_table(4,96) = x_opt(5);
	combined_CI_table(4,79) = CI_lb(6);
	combined_CI_table(4,80) = CI_ub(6);
	combined_CI_table(4,81) = x_opt(6);
	combined_CI_table(4,82) = CI_lb(7);
	combined_CI_table(4,83) = CI_ub(7);
	combined_CI_table(4,84) = x_opt(7);
	combined_CI_table(4,85) = CI_lb(8);
	combined_CI_table(4,86) = CI_ub(8);
	combined_CI_table(4,87) = x_opt(8);
	combined_CI_table(4,97) = CI_lb(9);
	combined_CI_table(4,98) = CI_ub(9);
	combined_CI_table(4,99) = x_opt(9);
	%===============
	lb = [0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01];
	ub = [5 5 50 50 50 50 1 4 4 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@doi_10_1038_nature11516_figureS8A_1_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(5,64) = CI_lb(1);
	combined_CI_table(5,65) = CI_ub(1);
	combined_CI_table(5,66) = x_opt(1);
	combined_CI_table(5,100) = CI_lb(2);
	combined_CI_table(5,101) = CI_ub(2);
	combined_CI_table(5,102) = x_opt(2);
	combined_CI_table(5,67) = CI_lb(3);
	combined_CI_table(5,68) = CI_ub(3);
	combined_CI_table(5,69) = x_opt(3);
	combined_CI_table(5,70) = CI_lb(4);
	combined_CI_table(5,71) = CI_ub(4);
	combined_CI_table(5,72) = x_opt(4);
	combined_CI_table(5,73) = CI_lb(5);
	combined_CI_table(5,74) = CI_ub(5);
	combined_CI_table(5,75) = x_opt(5);
	combined_CI_table(5,76) = CI_lb(6);
	combined_CI_table(5,77) = CI_ub(6);
	combined_CI_table(5,78) = x_opt(6);
	combined_CI_table(5,79) = CI_lb(7);
	combined_CI_table(5,80) = CI_ub(7);
	combined_CI_table(5,81) = x_opt(7);
	combined_CI_table(5,82) = CI_lb(8);
	combined_CI_table(5,83) = CI_ub(8);
	combined_CI_table(5,84) = x_opt(8);
	combined_CI_table(5,85) = CI_lb(9);
	combined_CI_table(5,86) = CI_ub(9);
	combined_CI_table(5,87) = x_opt(9);
	combined_CI_table(5,103) = CI_lb(10);
	combined_CI_table(5,104) = CI_ub(10);
	combined_CI_table(5,105) = x_opt(10);
	%===============
	lb = [0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01];
	ub = [5 5 50 50 50 50 1 4 4 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@doi_10_1038_nature11516_figureS8B_2_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(6,4) = CI_lb(1);
	combined_CI_table(6,5) = CI_ub(1);
	combined_CI_table(6,6) = x_opt(1);
	combined_CI_table(6,100) = CI_lb(2);
	combined_CI_table(6,101) = CI_ub(2);
	combined_CI_table(6,102) = x_opt(2);
	combined_CI_table(6,7) = CI_lb(3);
	combined_CI_table(6,8) = CI_ub(3);
	combined_CI_table(6,9) = x_opt(3);
	combined_CI_table(6,106) = CI_lb(4);
	combined_CI_table(6,107) = CI_ub(4);
	combined_CI_table(6,108) = x_opt(4);
	combined_CI_table(6,109) = CI_lb(5);
	combined_CI_table(6,110) = CI_ub(5);
	combined_CI_table(6,111) = x_opt(5);
	combined_CI_table(6,112) = CI_lb(6);
	combined_CI_table(6,113) = CI_ub(6);
	combined_CI_table(6,114) = x_opt(6);
	combined_CI_table(6,115) = CI_lb(7);
	combined_CI_table(6,116) = CI_ub(7);
	combined_CI_table(6,117) = x_opt(7);
	combined_CI_table(6,22) = CI_lb(8);
	combined_CI_table(6,23) = CI_ub(8);
	combined_CI_table(6,24) = x_opt(8);
	combined_CI_table(6,118) = CI_lb(9);
	combined_CI_table(6,119) = CI_ub(9);
	combined_CI_table(6,120) = x_opt(9);
	combined_CI_table(6,103) = CI_lb(10);
	combined_CI_table(6,104) = CI_ub(10);
	combined_CI_table(6,105) = x_opt(10);
	%===============
	lb = [0.03 0.001 0.001 0.002 0.002 0 1 1 0.01];
	ub = [5 50 50 50 50 1 4 4 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@doi_10_1038_nbt_2401_figureS11_1_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(7,121) = CI_lb(1);
	combined_CI_table(7,122) = CI_ub(1);
	combined_CI_table(7,123) = x_opt(1);
	combined_CI_table(7,7) = CI_lb(2);
	combined_CI_table(7,8) = CI_ub(2);
	combined_CI_table(7,9) = x_opt(2);
	combined_CI_table(7,106) = CI_lb(3);
	combined_CI_table(7,107) = CI_ub(3);
	combined_CI_table(7,108) = x_opt(3);
	combined_CI_table(7,124) = CI_lb(4);
	combined_CI_table(7,125) = CI_ub(4);
	combined_CI_table(7,126) = x_opt(4);
	combined_CI_table(7,109) = CI_lb(5);
	combined_CI_table(7,110) = CI_ub(5);
	combined_CI_table(7,111) = x_opt(5);
	combined_CI_table(7,115) = CI_lb(6);
	combined_CI_table(7,116) = CI_ub(6);
	combined_CI_table(7,117) = x_opt(6);
	combined_CI_table(7,22) = CI_lb(7);
	combined_CI_table(7,23) = CI_ub(7);
	combined_CI_table(7,24) = x_opt(7);
	combined_CI_table(7,118) = CI_lb(8);
	combined_CI_table(7,119) = CI_ub(8);
	combined_CI_table(7,120) = x_opt(8);
	combined_CI_table(7,127) = CI_lb(9);
	combined_CI_table(7,128) = CI_ub(9);
	combined_CI_table(7,129) = x_opt(9);
	%===============
	lb = [0.03 0.001 0.001 0.002 0.002 0 1 1 0.01];
	ub = [5 50 50 50 50 1 4 4 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@doi_10_1038_nbt_2401_figureS13_2_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(8,121) = CI_lb(1);
	combined_CI_table(8,122) = CI_ub(1);
	combined_CI_table(8,123) = x_opt(1);
	combined_CI_table(8,67) = CI_lb(2);
	combined_CI_table(8,68) = CI_ub(2);
	combined_CI_table(8,69) = x_opt(2);
	combined_CI_table(8,70) = CI_lb(3);
	combined_CI_table(8,71) = CI_ub(3);
	combined_CI_table(8,72) = x_opt(3);
	combined_CI_table(8,124) = CI_lb(4);
	combined_CI_table(8,125) = CI_ub(4);
	combined_CI_table(8,126) = x_opt(4);
	combined_CI_table(8,73) = CI_lb(5);
	combined_CI_table(8,74) = CI_ub(5);
	combined_CI_table(8,75) = x_opt(5);
	combined_CI_table(8,79) = CI_lb(6);
	combined_CI_table(8,80) = CI_ub(6);
	combined_CI_table(8,81) = x_opt(6);
	combined_CI_table(8,82) = CI_lb(7);
	combined_CI_table(8,83) = CI_ub(7);
	combined_CI_table(8,84) = x_opt(7);
	combined_CI_table(8,85) = CI_lb(8);
	combined_CI_table(8,86) = CI_ub(8);
	combined_CI_table(8,87) = x_opt(8);
	combined_CI_table(8,127) = CI_lb(9);
	combined_CI_table(8,128) = CI_ub(9);
	combined_CI_table(8,129) = x_opt(9);
	%===============
	lb = [0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01];
	ub = [5 5 50 50 50 50 1 4 4 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@doi_10_1038_nature09565_figure1C_1_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(9,130) = CI_lb(1);
	combined_CI_table(9,131) = CI_ub(1);
	combined_CI_table(9,132) = x_opt(1);
	combined_CI_table(9,133) = CI_lb(2);
	combined_CI_table(9,134) = CI_ub(2);
	combined_CI_table(9,135) = x_opt(2);
	combined_CI_table(9,67) = CI_lb(3);
	combined_CI_table(9,68) = CI_ub(3);
	combined_CI_table(9,69) = x_opt(3);
	combined_CI_table(9,70) = CI_lb(4);
	combined_CI_table(9,71) = CI_ub(4);
	combined_CI_table(9,72) = x_opt(4);
	combined_CI_table(9,136) = CI_lb(5);
	combined_CI_table(9,137) = CI_ub(5);
	combined_CI_table(9,138) = x_opt(5);
	combined_CI_table(9,73) = CI_lb(6);
	combined_CI_table(9,74) = CI_ub(6);
	combined_CI_table(9,75) = x_opt(6);
	combined_CI_table(9,79) = CI_lb(7);
	combined_CI_table(9,80) = CI_ub(7);
	combined_CI_table(9,81) = x_opt(7);
	combined_CI_table(9,82) = CI_lb(8);
	combined_CI_table(9,83) = CI_ub(8);
	combined_CI_table(9,84) = x_opt(8);
	combined_CI_table(9,85) = CI_lb(9);
	combined_CI_table(9,86) = CI_ub(9);
	combined_CI_table(9,87) = x_opt(9);
	combined_CI_table(9,139) = CI_lb(10);
	combined_CI_table(9,140) = CI_ub(10);
	combined_CI_table(9,141) = x_opt(10);
	%===============
	lb = [0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01];
	ub = [5 5 50 50 50 50 1 4 4 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@doi_10_1038_nature09565_figureS1C_2_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(10,142) = CI_lb(1);
	combined_CI_table(10,143) = CI_ub(1);
	combined_CI_table(10,144) = x_opt(1);
	combined_CI_table(10,130) = CI_lb(2);
	combined_CI_table(10,131) = CI_ub(2);
	combined_CI_table(10,132) = x_opt(2);
	combined_CI_table(10,37) = CI_lb(3);
	combined_CI_table(10,38) = CI_ub(3);
	combined_CI_table(10,39) = x_opt(3);
	combined_CI_table(10,40) = CI_lb(4);
	combined_CI_table(10,41) = CI_ub(4);
	combined_CI_table(10,42) = x_opt(4);
	combined_CI_table(10,136) = CI_lb(5);
	combined_CI_table(10,137) = CI_ub(5);
	combined_CI_table(10,138) = x_opt(5);
	combined_CI_table(10,43) = CI_lb(6);
	combined_CI_table(10,44) = CI_ub(6);
	combined_CI_table(10,45) = x_opt(6);
	combined_CI_table(10,49) = CI_lb(7);
	combined_CI_table(10,50) = CI_ub(7);
	combined_CI_table(10,51) = x_opt(7);
	combined_CI_table(10,52) = CI_lb(8);
	combined_CI_table(10,53) = CI_ub(8);
	combined_CI_table(10,54) = x_opt(8);
	combined_CI_table(10,55) = CI_lb(9);
	combined_CI_table(10,56) = CI_ub(9);
	combined_CI_table(10,57) = x_opt(9);
	combined_CI_table(10,139) = CI_lb(10);
	combined_CI_table(10,140) = CI_ub(10);
	combined_CI_table(10,141) = x_opt(10);
	%===============
	lb = [0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01];
	ub = [5 5 50 50 50 50 1 4 4 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@doi_10_1038_nbt_2149_Supplement_page17_1_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(11,145) = CI_lb(1);
	combined_CI_table(11,146) = CI_ub(1);
	combined_CI_table(11,147) = x_opt(1);
	combined_CI_table(11,4) = CI_lb(2);
	combined_CI_table(11,5) = CI_ub(2);
	combined_CI_table(11,6) = x_opt(2);
	combined_CI_table(11,7) = CI_lb(3);
	combined_CI_table(11,8) = CI_ub(3);
	combined_CI_table(11,9) = x_opt(3);
	combined_CI_table(11,148) = CI_lb(4);
	combined_CI_table(11,149) = CI_ub(4);
	combined_CI_table(11,150) = x_opt(4);
	combined_CI_table(11,16) = CI_lb(5);
	combined_CI_table(11,17) = CI_ub(5);
	combined_CI_table(11,18) = x_opt(5);
	combined_CI_table(11,151) = CI_lb(6);
	combined_CI_table(11,152) = CI_ub(6);
	combined_CI_table(11,153) = x_opt(6);
	combined_CI_table(11,154) = CI_lb(7);
	combined_CI_table(11,155) = CI_ub(7);
	combined_CI_table(11,156) = x_opt(7);
	combined_CI_table(11,22) = CI_lb(8);
	combined_CI_table(11,23) = CI_ub(8);
	combined_CI_table(11,24) = x_opt(8);
	combined_CI_table(11,157) = CI_lb(9);
	combined_CI_table(11,158) = CI_ub(9);
	combined_CI_table(11,159) = x_opt(9);
	combined_CI_table(11,160) = CI_lb(10);
	combined_CI_table(11,161) = CI_ub(10);
	combined_CI_table(11,162) = x_opt(10);
	%===============
	lb = [0.03 0.002 0.01];
	ub = [5 50 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@doi_10_1186_1754_1611_3_4_figure3_1_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(12,142) = CI_lb(1);
	combined_CI_table(12,143) = CI_ub(1);
	combined_CI_table(12,144) = x_opt(1);
	combined_CI_table(12,124) = CI_lb(2);
	combined_CI_table(12,125) = CI_ub(2);
	combined_CI_table(12,126) = x_opt(2);
	combined_CI_table(12,163) = CI_lb(3);
	combined_CI_table(12,164) = CI_ub(3);
	combined_CI_table(12,165) = x_opt(3);
	%===============
	lb = [0.03 0.002 0.01];
	ub = [5 50 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@doi_10_1186_1754_1611_3_4_figure3_2_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(13,142) = CI_lb(1);
	combined_CI_table(13,143) = CI_ub(1);
	combined_CI_table(13,144) = x_opt(1);
	combined_CI_table(13,166) = CI_lb(2);
	combined_CI_table(13,167) = CI_ub(2);
	combined_CI_table(13,168) = x_opt(2);
	combined_CI_table(13,163) = CI_lb(3);
	combined_CI_table(13,164) = CI_ub(3);
	combined_CI_table(13,165) = x_opt(3);
	%===============
	lb = [0.03 0.002 0.01];
	ub = [5 50 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@doi_10_1186_1754_1611_3_4_figure3_3_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(14,142) = CI_lb(1);
	combined_CI_table(14,143) = CI_ub(1);
	combined_CI_table(14,144) = x_opt(1);
	combined_CI_table(14,169) = CI_lb(2);
	combined_CI_table(14,170) = CI_ub(2);
	combined_CI_table(14,171) = x_opt(2);
	combined_CI_table(14,163) = CI_lb(3);
	combined_CI_table(14,164) = CI_ub(3);
	combined_CI_table(14,165) = x_opt(3);
	%===============
	lb = [0.03 0.002 0.01];
	ub = [5 50 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@doi_10_1186_1754_1611_3_4_figure3_4_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(15,142) = CI_lb(1);
	combined_CI_table(15,143) = CI_ub(1);
	combined_CI_table(15,144) = x_opt(1);
	combined_CI_table(15,172) = CI_lb(2);
	combined_CI_table(15,173) = CI_ub(2);
	combined_CI_table(15,174) = x_opt(2);
	combined_CI_table(15,163) = CI_lb(3);
	combined_CI_table(15,164) = CI_ub(3);
	combined_CI_table(15,165) = x_opt(3);
	%===============
	lb = [0.03 0.002 0.01];
	ub = [5 50 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@doi_10_1186_1754_1611_3_4_figure3_5_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(16,142) = CI_lb(1);
	combined_CI_table(16,143) = CI_ub(1);
	combined_CI_table(16,144) = x_opt(1);
	combined_CI_table(16,175) = CI_lb(2);
	combined_CI_table(16,176) = CI_ub(2);
	combined_CI_table(16,177) = x_opt(2);
	combined_CI_table(16,163) = CI_lb(3);
	combined_CI_table(16,164) = CI_ub(3);
	combined_CI_table(16,165) = x_opt(3);
	%===============
	lb = [0.03 0.002 0.01];
	ub = [5 50 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@doi_10_1186_1754_1611_3_4_figure3_6_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(17,142) = CI_lb(1);
	combined_CI_table(17,143) = CI_ub(1);
	combined_CI_table(17,144) = x_opt(1);
	combined_CI_table(17,43) = CI_lb(2);
	combined_CI_table(17,44) = CI_ub(2);
	combined_CI_table(17,45) = x_opt(2);
	combined_CI_table(17,163) = CI_lb(3);
	combined_CI_table(17,164) = CI_ub(3);
	combined_CI_table(17,165) = x_opt(3);
	%===============
	lb = [0.03 0.002 0.01];
	ub = [5 50 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@doi_10_1186_1754_1611_3_4_figure3_7_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(18,142) = CI_lb(1);
	combined_CI_table(18,143) = CI_ub(1);
	combined_CI_table(18,144) = x_opt(1);
	combined_CI_table(18,94) = CI_lb(2);
	combined_CI_table(18,95) = CI_ub(2);
	combined_CI_table(18,96) = x_opt(2);
	combined_CI_table(18,163) = CI_lb(3);
	combined_CI_table(18,164) = CI_ub(3);
	combined_CI_table(18,165) = x_opt(3);
	%===============
	lb = [0.03 0.002 0.01];
	ub = [5 50 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@doi__10_1073_pnas_1301301110_sd03_xls_1_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(19,142) = CI_lb(1);
	combined_CI_table(19,143) = CI_ub(1);
	combined_CI_table(19,144) = x_opt(1);
	combined_CI_table(19,136) = CI_lb(2);
	combined_CI_table(19,137) = CI_ub(2);
	combined_CI_table(19,138) = x_opt(2);
	combined_CI_table(19,178) = CI_lb(3);
	combined_CI_table(19,179) = CI_ub(3);
	combined_CI_table(19,180) = x_opt(3);
	%===============
	lb = [0.03 0.002 0.01];
	ub = [5 50 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@doi__10_1073_pnas_1301301110_sd03_xls_2_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(20,142) = CI_lb(1);
	combined_CI_table(20,143) = CI_ub(1);
	combined_CI_table(20,144) = x_opt(1);
	combined_CI_table(20,124) = CI_lb(2);
	combined_CI_table(20,125) = CI_ub(2);
	combined_CI_table(20,126) = x_opt(2);
	combined_CI_table(20,178) = CI_lb(3);
	combined_CI_table(20,179) = CI_ub(3);
	combined_CI_table(20,180) = x_opt(3);
	%===============
	lb = [0.03 0.002 0.01];
	ub = [5 50 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@doi_10_1093_nar_gku1388_figure2A_1_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(21,91) = CI_lb(1);
	combined_CI_table(21,92) = CI_ub(1);
	combined_CI_table(21,93) = x_opt(1);
	combined_CI_table(21,124) = CI_lb(2);
	combined_CI_table(21,125) = CI_ub(2);
	combined_CI_table(21,126) = x_opt(2);
	combined_CI_table(21,181) = CI_lb(3);
	combined_CI_table(21,182) = CI_ub(3);
	combined_CI_table(21,183) = x_opt(3);
	%===============
	lb = [0.03 0.002 0.01];
	ub = [5 50 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@doi_10_1093_nar_gku1388_figure2A_2_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(22,91) = CI_lb(1);
	combined_CI_table(22,92) = CI_ub(1);
	combined_CI_table(22,93) = x_opt(1);
	combined_CI_table(22,184) = CI_lb(2);
	combined_CI_table(22,185) = CI_ub(2);
	combined_CI_table(22,186) = x_opt(2);
	combined_CI_table(22,181) = CI_lb(3);
	combined_CI_table(22,182) = CI_ub(3);
	combined_CI_table(22,183) = x_opt(3);
	%===============
	lb = [0.03 0.002 0.01];
	ub = [5 50 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@doi_10_1093_nar_gku1388_figure2A_3_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(23,91) = CI_lb(1);
	combined_CI_table(23,92) = CI_ub(1);
	combined_CI_table(23,93) = x_opt(1);
	combined_CI_table(23,187) = CI_lb(2);
	combined_CI_table(23,188) = CI_ub(2);
	combined_CI_table(23,189) = x_opt(2);
	combined_CI_table(23,181) = CI_lb(3);
	combined_CI_table(23,182) = CI_ub(3);
	combined_CI_table(23,183) = x_opt(3);
	%===============
	lb = [0.03 0.002 0.01];
	ub = [5 50 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@doi_10_1093_nar_gku1388_figure2A_4_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(24,91) = CI_lb(1);
	combined_CI_table(24,92) = CI_ub(1);
	combined_CI_table(24,93) = x_opt(1);
	combined_CI_table(24,190) = CI_lb(2);
	combined_CI_table(24,191) = CI_ub(2);
	combined_CI_table(24,192) = x_opt(2);
	combined_CI_table(24,181) = CI_lb(3);
	combined_CI_table(24,182) = CI_ub(3);
	combined_CI_table(24,183) = x_opt(3);
	%===============
	lb = [0.03 0.002 0.01];
	ub = [5 50 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@doi_10_1093_nar_gku1388_figure2A_5_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(25,91) = CI_lb(1);
	combined_CI_table(25,92) = CI_ub(1);
	combined_CI_table(25,93) = x_opt(1);
	combined_CI_table(25,193) = CI_lb(2);
	combined_CI_table(25,194) = CI_ub(2);
	combined_CI_table(25,195) = x_opt(2);
	combined_CI_table(25,181) = CI_lb(3);
	combined_CI_table(25,182) = CI_ub(3);
	combined_CI_table(25,183) = x_opt(3);
	%===============
	lb = [0.03 0.002 0.01];
	ub = [5 50 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@doi_10_1093_nar_gku1388_figure2A_6_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(26,91) = CI_lb(1);
	combined_CI_table(26,92) = CI_ub(1);
	combined_CI_table(26,93) = x_opt(1);
	combined_CI_table(26,136) = CI_lb(2);
	combined_CI_table(26,137) = CI_ub(2);
	combined_CI_table(26,138) = x_opt(2);
	combined_CI_table(26,181) = CI_lb(3);
	combined_CI_table(26,182) = CI_ub(3);
	combined_CI_table(26,183) = x_opt(3);
	%===============
	lb = [0.03 0.001 0.001 0.002 0.002 0 1 1 0.01];
	ub = [5 50 50 50 50 1 4 4 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@doi_10_1093_nar_gku1388_figure2B_7_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(27,91) = CI_lb(1);
	combined_CI_table(27,92) = CI_ub(1);
	combined_CI_table(27,93) = x_opt(1);
	combined_CI_table(27,37) = CI_lb(2);
	combined_CI_table(27,38) = CI_ub(2);
	combined_CI_table(27,39) = x_opt(2);
	combined_CI_table(27,196) = CI_lb(3);
	combined_CI_table(27,197) = CI_ub(3);
	combined_CI_table(27,198) = x_opt(3);
	combined_CI_table(27,136) = CI_lb(4);
	combined_CI_table(27,137) = CI_ub(4);
	combined_CI_table(27,138) = x_opt(4);
	combined_CI_table(27,199) = CI_lb(5);
	combined_CI_table(27,200) = CI_ub(5);
	combined_CI_table(27,201) = x_opt(5);
	combined_CI_table(27,202) = CI_lb(6);
	combined_CI_table(27,203) = CI_ub(6);
	combined_CI_table(27,204) = x_opt(6);
	combined_CI_table(27,52) = CI_lb(7);
	combined_CI_table(27,53) = CI_ub(7);
	combined_CI_table(27,54) = x_opt(7);
	combined_CI_table(27,205) = CI_lb(8);
	combined_CI_table(27,206) = CI_ub(8);
	combined_CI_table(27,207) = x_opt(8);
	combined_CI_table(27,181) = CI_lb(9);
	combined_CI_table(27,182) = CI_ub(9);
	combined_CI_table(27,183) = x_opt(9);
	%===============
	lb = [0.03 0.001 0.001 0.002 0.002 0 1 1 0.01];
	ub = [5 50 50 50 50 1 4 4 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@doi_10_1093_nar_gku1388_figure2B_8_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(28,91) = CI_lb(1);
	combined_CI_table(28,92) = CI_ub(1);
	combined_CI_table(28,93) = x_opt(1);
	combined_CI_table(28,37) = CI_lb(2);
	combined_CI_table(28,38) = CI_ub(2);
	combined_CI_table(28,39) = x_opt(2);
	combined_CI_table(28,196) = CI_lb(3);
	combined_CI_table(28,197) = CI_ub(3);
	combined_CI_table(28,198) = x_opt(3);
	combined_CI_table(28,193) = CI_lb(4);
	combined_CI_table(28,194) = CI_ub(4);
	combined_CI_table(28,195) = x_opt(4);
	combined_CI_table(28,199) = CI_lb(5);
	combined_CI_table(28,200) = CI_ub(5);
	combined_CI_table(28,201) = x_opt(5);
	combined_CI_table(28,202) = CI_lb(6);
	combined_CI_table(28,203) = CI_ub(6);
	combined_CI_table(28,204) = x_opt(6);
	combined_CI_table(28,52) = CI_lb(7);
	combined_CI_table(28,53) = CI_ub(7);
	combined_CI_table(28,54) = x_opt(7);
	combined_CI_table(28,205) = CI_lb(8);
	combined_CI_table(28,206) = CI_ub(8);
	combined_CI_table(28,207) = x_opt(8);
	combined_CI_table(28,181) = CI_lb(9);
	combined_CI_table(28,182) = CI_ub(9);
	combined_CI_table(28,183) = x_opt(9);
	%===============
	lb = [0.03 0.001 0.001 0.002 0.002 0 1 1 0.01];
	ub = [5 50 50 50 50 1 4 4 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@doi_10_1093_nar_gku1388_figure2B_9_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(29,91) = CI_lb(1);
	combined_CI_table(29,92) = CI_ub(1);
	combined_CI_table(29,93) = x_opt(1);
	combined_CI_table(29,37) = CI_lb(2);
	combined_CI_table(29,38) = CI_ub(2);
	combined_CI_table(29,39) = x_opt(2);
	combined_CI_table(29,196) = CI_lb(3);
	combined_CI_table(29,197) = CI_ub(3);
	combined_CI_table(29,198) = x_opt(3);
	combined_CI_table(29,190) = CI_lb(4);
	combined_CI_table(29,191) = CI_ub(4);
	combined_CI_table(29,192) = x_opt(4);
	combined_CI_table(29,199) = CI_lb(5);
	combined_CI_table(29,200) = CI_ub(5);
	combined_CI_table(29,201) = x_opt(5);
	combined_CI_table(29,202) = CI_lb(6);
	combined_CI_table(29,203) = CI_ub(6);
	combined_CI_table(29,204) = x_opt(6);
	combined_CI_table(29,52) = CI_lb(7);
	combined_CI_table(29,53) = CI_ub(7);
	combined_CI_table(29,54) = x_opt(7);
	combined_CI_table(29,205) = CI_lb(8);
	combined_CI_table(29,206) = CI_ub(8);
	combined_CI_table(29,207) = x_opt(8);
	combined_CI_table(29,181) = CI_lb(9);
	combined_CI_table(29,182) = CI_ub(9);
	combined_CI_table(29,183) = x_opt(9);
	%===============
	lb = [0.03 0.001 0.001 0.002 0.002 0 1 1 0.01];
	ub = [5 50 50 50 50 1 4 4 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@doi_10_1093_nar_gku1388_figure2B_10_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(30,91) = CI_lb(1);
	combined_CI_table(30,92) = CI_ub(1);
	combined_CI_table(30,93) = x_opt(1);
	combined_CI_table(30,37) = CI_lb(2);
	combined_CI_table(30,38) = CI_ub(2);
	combined_CI_table(30,39) = x_opt(2);
	combined_CI_table(30,196) = CI_lb(3);
	combined_CI_table(30,197) = CI_ub(3);
	combined_CI_table(30,198) = x_opt(3);
	combined_CI_table(30,187) = CI_lb(4);
	combined_CI_table(30,188) = CI_ub(4);
	combined_CI_table(30,189) = x_opt(4);
	combined_CI_table(30,199) = CI_lb(5);
	combined_CI_table(30,200) = CI_ub(5);
	combined_CI_table(30,201) = x_opt(5);
	combined_CI_table(30,202) = CI_lb(6);
	combined_CI_table(30,203) = CI_ub(6);
	combined_CI_table(30,204) = x_opt(6);
	combined_CI_table(30,52) = CI_lb(7);
	combined_CI_table(30,53) = CI_ub(7);
	combined_CI_table(30,54) = x_opt(7);
	combined_CI_table(30,205) = CI_lb(8);
	combined_CI_table(30,206) = CI_ub(8);
	combined_CI_table(30,207) = x_opt(8);
	combined_CI_table(30,181) = CI_lb(9);
	combined_CI_table(30,182) = CI_ub(9);
	combined_CI_table(30,183) = x_opt(9);
	%===============
	lb = [0.03 0.001 0.001 0.002 0.002 0 1 1 0.01];
	ub = [5 50 50 50 50 1 4 4 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@doi_10_1093_nar_gku1388_figure2B_11_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(31,91) = CI_lb(1);
	combined_CI_table(31,92) = CI_ub(1);
	combined_CI_table(31,93) = x_opt(1);
	combined_CI_table(31,37) = CI_lb(2);
	combined_CI_table(31,38) = CI_ub(2);
	combined_CI_table(31,39) = x_opt(2);
	combined_CI_table(31,196) = CI_lb(3);
	combined_CI_table(31,197) = CI_ub(3);
	combined_CI_table(31,198) = x_opt(3);
	combined_CI_table(31,184) = CI_lb(4);
	combined_CI_table(31,185) = CI_ub(4);
	combined_CI_table(31,186) = x_opt(4);
	combined_CI_table(31,199) = CI_lb(5);
	combined_CI_table(31,200) = CI_ub(5);
	combined_CI_table(31,201) = x_opt(5);
	combined_CI_table(31,202) = CI_lb(6);
	combined_CI_table(31,203) = CI_ub(6);
	combined_CI_table(31,204) = x_opt(6);
	combined_CI_table(31,52) = CI_lb(7);
	combined_CI_table(31,53) = CI_ub(7);
	combined_CI_table(31,54) = x_opt(7);
	combined_CI_table(31,205) = CI_lb(8);
	combined_CI_table(31,206) = CI_ub(8);
	combined_CI_table(31,207) = x_opt(8);
	combined_CI_table(31,181) = CI_lb(9);
	combined_CI_table(31,182) = CI_ub(9);
	combined_CI_table(31,183) = x_opt(9);
	%===============
	lb = [0.03 0.001 0.001 0.002 0.002 0 1 1 0.01];
	ub = [5 50 50 50 50 1 4 4 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@doi_10_1093_nar_gku1388_figure2B_12_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(32,91) = CI_lb(1);
	combined_CI_table(32,92) = CI_ub(1);
	combined_CI_table(32,93) = x_opt(1);
	combined_CI_table(32,37) = CI_lb(2);
	combined_CI_table(32,38) = CI_ub(2);
	combined_CI_table(32,39) = x_opt(2);
	combined_CI_table(32,196) = CI_lb(3);
	combined_CI_table(32,197) = CI_ub(3);
	combined_CI_table(32,198) = x_opt(3);
	combined_CI_table(32,124) = CI_lb(4);
	combined_CI_table(32,125) = CI_ub(4);
	combined_CI_table(32,126) = x_opt(4);
	combined_CI_table(32,199) = CI_lb(5);
	combined_CI_table(32,200) = CI_ub(5);
	combined_CI_table(32,201) = x_opt(5);
	combined_CI_table(32,202) = CI_lb(6);
	combined_CI_table(32,203) = CI_ub(6);
	combined_CI_table(32,204) = x_opt(6);
	combined_CI_table(32,52) = CI_lb(7);
	combined_CI_table(32,53) = CI_ub(7);
	combined_CI_table(32,54) = x_opt(7);
	combined_CI_table(32,205) = CI_lb(8);
	combined_CI_table(32,206) = CI_ub(8);
	combined_CI_table(32,207) = x_opt(8);
	combined_CI_table(32,181) = CI_lb(9);
	combined_CI_table(32,182) = CI_ub(9);
	combined_CI_table(32,183) = x_opt(9);
	%===============
	lb = [0.03 0.002 0.01];
	ub = [5 50 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@doi_10_1093_nar_gku593_figure3D_1_1_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(33,91) = CI_lb(1);
	combined_CI_table(33,92) = CI_ub(1);
	combined_CI_table(33,93) = x_opt(1);
	combined_CI_table(33,208) = CI_lb(2);
	combined_CI_table(33,209) = CI_ub(2);
	combined_CI_table(33,210) = x_opt(2);
	combined_CI_table(33,211) = CI_lb(3);
	combined_CI_table(33,212) = CI_ub(3);
	combined_CI_table(33,213) = x_opt(3);
	%===============
	lb = [0.03 0.002 0.01];
	ub = [5 50 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@doi_10_1093_nar_gku593_figure3D_2_2_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(34,91) = CI_lb(1);
	combined_CI_table(34,92) = CI_ub(1);
	combined_CI_table(34,93) = x_opt(1);
	combined_CI_table(34,193) = CI_lb(2);
	combined_CI_table(34,194) = CI_ub(2);
	combined_CI_table(34,195) = x_opt(2);
	combined_CI_table(34,211) = CI_lb(3);
	combined_CI_table(34,212) = CI_ub(3);
	combined_CI_table(34,213) = x_opt(3);
	%===============
	lb = [0.03 0.002 0.01];
	ub = [5 50 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@doi_10_1093_nar_gku593_figure3D_3_3_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(35,91) = CI_lb(1);
	combined_CI_table(35,92) = CI_ub(1);
	combined_CI_table(35,93) = x_opt(1);
	combined_CI_table(35,166) = CI_lb(2);
	combined_CI_table(35,167) = CI_ub(2);
	combined_CI_table(35,168) = x_opt(2);
	combined_CI_table(35,211) = CI_lb(3);
	combined_CI_table(35,212) = CI_ub(3);
	combined_CI_table(35,213) = x_opt(3);
	%===============
	lb = [0.03 0.002 0.01];
	ub = [5 50 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@doi_10_1093_nar_gku593_figure3D_4_4_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(36,91) = CI_lb(1);
	combined_CI_table(36,92) = CI_ub(1);
	combined_CI_table(36,93) = x_opt(1);
	combined_CI_table(36,190) = CI_lb(2);
	combined_CI_table(36,191) = CI_ub(2);
	combined_CI_table(36,192) = x_opt(2);
	combined_CI_table(36,211) = CI_lb(3);
	combined_CI_table(36,212) = CI_ub(3);
	combined_CI_table(36,213) = x_opt(3);
	%===============
	lb = [0.03 0.002 0.01];
	ub = [5 50 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@doi_10_1093_nar_gku593_figure3D_5_5_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(37,91) = CI_lb(1);
	combined_CI_table(37,92) = CI_ub(1);
	combined_CI_table(37,93) = x_opt(1);
	combined_CI_table(37,187) = CI_lb(2);
	combined_CI_table(37,188) = CI_ub(2);
	combined_CI_table(37,189) = x_opt(2);
	combined_CI_table(37,211) = CI_lb(3);
	combined_CI_table(37,212) = CI_ub(3);
	combined_CI_table(37,213) = x_opt(3);
	%===============
	lb = [0.03 0.002 0.01];
	ub = [5 50 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@doi_10_1093_nar_gku593_figure3D_6_6_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(38,91) = CI_lb(1);
	combined_CI_table(38,92) = CI_ub(1);
	combined_CI_table(38,93) = x_opt(1);
	combined_CI_table(38,184) = CI_lb(2);
	combined_CI_table(38,185) = CI_ub(2);
	combined_CI_table(38,186) = x_opt(2);
	combined_CI_table(38,211) = CI_lb(3);
	combined_CI_table(38,212) = CI_ub(3);
	combined_CI_table(38,213) = x_opt(3);
	%===============
	lb = [0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01];
	ub = [5 5 50 50 50 50 1 4 4 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@doi_10_1093_nar_gku593_figureS6_7_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(39,91) = CI_lb(1);
	combined_CI_table(39,92) = CI_ub(1);
	combined_CI_table(39,93) = x_opt(1);
	combined_CI_table(39,64) = CI_lb(2);
	combined_CI_table(39,65) = CI_ub(2);
	combined_CI_table(39,66) = x_opt(2);
	combined_CI_table(39,67) = CI_lb(3);
	combined_CI_table(39,68) = CI_ub(3);
	combined_CI_table(39,69) = x_opt(3);
	combined_CI_table(39,70) = CI_lb(4);
	combined_CI_table(39,71) = CI_ub(4);
	combined_CI_table(39,72) = x_opt(4);
	combined_CI_table(39,73) = CI_lb(5);
	combined_CI_table(39,74) = CI_ub(5);
	combined_CI_table(39,75) = x_opt(5);
	combined_CI_table(39,76) = CI_lb(6);
	combined_CI_table(39,77) = CI_ub(6);
	combined_CI_table(39,78) = x_opt(6);
	combined_CI_table(39,79) = CI_lb(7);
	combined_CI_table(39,80) = CI_ub(7);
	combined_CI_table(39,81) = x_opt(7);
	combined_CI_table(39,82) = CI_lb(8);
	combined_CI_table(39,83) = CI_ub(8);
	combined_CI_table(39,84) = x_opt(8);
	combined_CI_table(39,85) = CI_lb(9);
	combined_CI_table(39,86) = CI_ub(9);
	combined_CI_table(39,87) = x_opt(9);
	combined_CI_table(39,211) = CI_lb(10);
	combined_CI_table(39,212) = CI_ub(10);
	combined_CI_table(39,213) = x_opt(10);
	%===============
	lb = [0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01];
	ub = [5 5 50 50 50 50 1 4 4 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(40,4) = CI_lb(1);
	combined_CI_table(40,5) = CI_ub(1);
	combined_CI_table(40,6) = x_opt(1);
	combined_CI_table(40,34) = CI_lb(2);
	combined_CI_table(40,35) = CI_ub(2);
	combined_CI_table(40,36) = x_opt(2);
	combined_CI_table(40,7) = CI_lb(3);
	combined_CI_table(40,8) = CI_ub(3);
	combined_CI_table(40,9) = x_opt(3);
	combined_CI_table(40,10) = CI_lb(4);
	combined_CI_table(40,11) = CI_ub(4);
	combined_CI_table(40,12) = x_opt(4);
	combined_CI_table(40,13) = CI_lb(5);
	combined_CI_table(40,14) = CI_ub(5);
	combined_CI_table(40,15) = x_opt(5);
	combined_CI_table(40,112) = CI_lb(6);
	combined_CI_table(40,113) = CI_ub(6);
	combined_CI_table(40,114) = x_opt(6);
	combined_CI_table(40,19) = CI_lb(7);
	combined_CI_table(40,20) = CI_ub(7);
	combined_CI_table(40,21) = x_opt(7);
	combined_CI_table(40,22) = CI_lb(8);
	combined_CI_table(40,23) = CI_ub(8);
	combined_CI_table(40,24) = x_opt(8);
	combined_CI_table(40,25) = CI_lb(9);
	combined_CI_table(40,26) = CI_ub(9);
	combined_CI_table(40,27) = x_opt(9);
	combined_CI_table(40,214) = CI_lb(10);
	combined_CI_table(40,215) = CI_ub(10);
	combined_CI_table(40,216) = x_opt(10);
	%===============
	lb = [0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01];
	ub = [5 5 50 50 50 50 1 4 4 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@doi_10_1073pnas_0408507102_figure2A_1_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	combined_CI_table(41,217) = CI_lb(1);
	combined_CI_table(41,218) = CI_ub(1);
	combined_CI_table(41,219) = x_opt(1);
	combined_CI_table(41,4) = CI_lb(2);
	combined_CI_table(41,5) = CI_ub(2);
	combined_CI_table(41,6) = x_opt(2);
	combined_CI_table(41,37) = CI_lb(3);
	combined_CI_table(41,38) = CI_ub(3);
	combined_CI_table(41,39) = x_opt(3);
	combined_CI_table(41,40) = CI_lb(4);
	combined_CI_table(41,41) = CI_ub(4);
	combined_CI_table(41,42) = x_opt(4);
	combined_CI_table(41,43) = CI_lb(5);
	combined_CI_table(41,44) = CI_ub(5);
	combined_CI_table(41,45) = x_opt(5);
	combined_CI_table(41,112) = CI_lb(6);
	combined_CI_table(41,113) = CI_ub(6);
	combined_CI_table(41,114) = x_opt(6);
	combined_CI_table(41,49) = CI_lb(7);
	combined_CI_table(41,50) = CI_ub(7);
	combined_CI_table(41,51) = x_opt(7);
	combined_CI_table(41,52) = CI_lb(8);
	combined_CI_table(41,53) = CI_ub(8);
	combined_CI_table(41,54) = x_opt(8);
	combined_CI_table(41,55) = CI_lb(9);
	combined_CI_table(41,56) = CI_ub(9);
	combined_CI_table(41,57) = x_opt(9);
	combined_CI_table(41,220) = CI_lb(10);
	combined_CI_table(41,221) = CI_ub(10);
	combined_CI_table(41,222) = x_opt(10);
	print_a_matrix('Plot_data/separate_parameter.dat', combined_CI_table);
	parameter_name_list = {'ALPHA_RBS-pLAC'; 'ALPHA_RBS-placI'; 'K_IPTG'; 'K_pLAC'; 'alpha_pLAC'; 'alpha_placI'; 'beta_pLAC'; 'n_IPTG'; 'n_pLAC'; 'scale_lacZ_doi_10_1073_pnas_0606717104'; 'ALPHA_RBS-pN25'; 'ALPHA_RBSII'; 'K_aTc'; 'K_pLtetO-1'; 'alpha_pLtetO-1'; 'alpha_pN25'; 'beta_pLtetO-1'; 'n_aTc'; 'n_pLtetO-1'; 'scale_luc_doi__10_1093_nar_25_6_1203'; 'ALPHA_RBS-pBAD'; 'ALPHA_RBS-pC'; 'K_arabinose'; 'K_pBAD'; 'alpha_pBAD'; 'alpha_pC'; 'beta_pBAD'; 'n_arabinose'; 'n_pBAD'; 'scale_GFPuv_doi_10_1128_AEM_00791_07'; 'ALPHA_BBa_B0030'; 'alpha_pLlacO-1'; 'scale_mCherry_doi_10_1038_nature12148'; 'ALPHA_RBS-psicA'; 'scale_mRFP1_doi_10_1038_nature11516'; 'K_pTAC'; 'alpha_pTAC'; 'alpha_placIq'; 'beta_pTAC'; 'n_pTAC'; 'ALPHA_RBS-A'; 'alpha_BBa_J23101'; 'scale_sfGFP_doi_10_1038_nbt_2401'; 'ALPHA_BBa_B0033'; 'ALPHA_BBa_B0034'; 'alpha_BBa_J23117'; 'scale_eYFP_doi_10_1038_nature09565'; 'ALPHA_BBa_B0032'; 'ALPHA_RBS-pET-29b'; 'K_placUV5'; 'alpha_placUV5'; 'beta_placUV5'; 'n_placUV5'; 'scale_mRFP1_doi_10_1038_nbt_2149'; 'scale_GFPmut3b_doi_10_1186_1754_1611_3_4'; 'alpha_BBa_J23116'; 'alpha_BBa_J23150'; 'alpha_BBa_J23151'; 'alpha_BBa_J23102'; 'scale_sfGFP_doi__10_1073_pnas_1301301110'; 'scale_GFPmut3b_doi_10_1093_nar_gku1388'; 'alpha_BBa_J23106'; 'alpha_BBa_J23105'; 'alpha_BBa_J23115'; 'alpha_BBa_J23114'; 'K_pTET*'; 'alpha_pTET*'; 'beta_pTET*'; 'n_pTET*'; 'alpha_BBa_J23109'; 'scale_GFPmut3b_doi_10_1093_nar_gku593'; 'scale_eYFP_http___dspace_mit_edu_handle_1721_1_8228'; 'ALPHA_RBS-Unknown-Hooshangi'; 'scale_eYFP_doi_10_1073pnas_0408507102'};
	print_union_CI(parameter_name_list, combined_CI_table, 'Plot_data/separate_parameter_CI.dat');
end

function parameter_opt = sequential_fitting(N1, N2)
	error_opt = 1e6;
	parameter_opt = zeros(1,74);
	for i = 1:N1
		model_order = randperm(41);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			parameter_opt = parameter;
		end
	end
	parameter_name_list = {'ALPHA_RBS-pLAC'; 'ALPHA_RBS-placI'; 'K_IPTG'; 'K_pLAC'; 'alpha_pLAC'; 'alpha_placI'; 'beta_pLAC'; 'n_IPTG'; 'n_pLAC'; 'scale_lacZ_doi_10_1073_pnas_0606717104'; 'ALPHA_RBS-pN25'; 'ALPHA_RBSII'; 'K_aTc'; 'K_pLtetO-1'; 'alpha_pLtetO-1'; 'alpha_pN25'; 'beta_pLtetO-1'; 'n_aTc'; 'n_pLtetO-1'; 'scale_luc_doi__10_1093_nar_25_6_1203'; 'ALPHA_RBS-pBAD'; 'ALPHA_RBS-pC'; 'K_arabinose'; 'K_pBAD'; 'alpha_pBAD'; 'alpha_pC'; 'beta_pBAD'; 'n_arabinose'; 'n_pBAD'; 'scale_GFPuv_doi_10_1128_AEM_00791_07'; 'ALPHA_BBa_B0030'; 'alpha_pLlacO-1'; 'scale_mCherry_doi_10_1038_nature12148'; 'ALPHA_RBS-psicA'; 'scale_mRFP1_doi_10_1038_nature11516'; 'K_pTAC'; 'alpha_pTAC'; 'alpha_placIq'; 'beta_pTAC'; 'n_pTAC'; 'ALPHA_RBS-A'; 'alpha_BBa_J23101'; 'scale_sfGFP_doi_10_1038_nbt_2401'; 'ALPHA_BBa_B0033'; 'ALPHA_BBa_B0034'; 'alpha_BBa_J23117'; 'scale_eYFP_doi_10_1038_nature09565'; 'ALPHA_BBa_B0032'; 'ALPHA_RBS-pET-29b'; 'K_placUV5'; 'alpha_placUV5'; 'beta_placUV5'; 'n_placUV5'; 'scale_mRFP1_doi_10_1038_nbt_2149'; 'scale_GFPmut3b_doi_10_1186_1754_1611_3_4'; 'alpha_BBa_J23116'; 'alpha_BBa_J23150'; 'alpha_BBa_J23151'; 'alpha_BBa_J23102'; 'scale_sfGFP_doi__10_1073_pnas_1301301110'; 'scale_GFPmut3b_doi_10_1093_nar_gku1388'; 'alpha_BBa_J23106'; 'alpha_BBa_J23105'; 'alpha_BBa_J23115'; 'alpha_BBa_J23114'; 'K_pTET*'; 'alpha_pTET*'; 'beta_pTET*'; 'n_pTET*'; 'alpha_BBa_J23109'; 'scale_GFPmut3b_doi_10_1093_nar_gku593'; 'scale_eYFP_http___dspace_mit_edu_handle_1721_1_8228'; 'ALPHA_RBS-Unknown-Hooshangi'; 'scale_eYFP_doi_10_1073pnas_0408507102'};
	parameter_fileID = fopen('Plot_data/sequential_parameter.dat','w');
	for i = 1:length(parameter_name_list)
		fprintf(parameter_fileID, '%s\t', parameter_name_list{i});
		fprintf(parameter_fileID, '%f\n', parameter_opt(i));
	end
	fclose(parameter_fileID);
end
function best_parameter = sequential_fitting_for_one_order(model_order, N)
	options = optimset('MaxIter', 1000, 'MaxFunEvals', 30000);
	best_parameter = zeros(1,74);
	for model_id = 1:length(model_order)
		switch (model_id)
			case model_order(model_id)
				lb = [0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01];
				ub = [5 5 50 50 50 50 1 4 4 10];
				if (best_parameter(1) > 1e-10)
					lb(1) = best_parameter(1) - 1e-10;
					ub(1) = best_parameter(1) + 1e-10;
				end
				if (best_parameter(2) > 1e-10)
					lb(2) = best_parameter(2) - 1e-10;
					ub(2) = best_parameter(2) + 1e-10;
				end
				if (best_parameter(3) > 1e-10)
					lb(3) = best_parameter(3) - 1e-10;
					ub(3) = best_parameter(3) + 1e-10;
				end
				if (best_parameter(4) > 1e-10)
					lb(4) = best_parameter(4) - 1e-10;
					ub(4) = best_parameter(4) + 1e-10;
				end
				if (best_parameter(5) > 1e-10)
					lb(5) = best_parameter(5) - 1e-10;
					ub(5) = best_parameter(5) + 1e-10;
				end
				if (best_parameter(6) > 1e-10)
					lb(6) = best_parameter(6) - 1e-10;
					ub(6) = best_parameter(6) + 1e-10;
				end
				if (best_parameter(7) > 1e-10)
					lb(7) = best_parameter(7) - 1e-10;
					ub(7) = best_parameter(7) + 1e-10;
				end
				if (best_parameter(8) > 1e-10)
					lb(8) = best_parameter(8) - 1e-10;
					ub(8) = best_parameter(8) + 1e-10;
				end
				if (best_parameter(9) > 1e-10)
					lb(9) = best_parameter(9) - 1e-10;
					ub(9) = best_parameter(9) + 1e-10;
				end
				if (best_parameter(10) > 1e-10)
					lb(10) = best_parameter(10) - 1e-10;
					ub(10) = best_parameter(10) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@doi_10_1073_pnas_0606717104_figure4_6_1_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(1) = x_opt(1);
				best_parameter(2) = x_opt(2);
				best_parameter(3) = x_opt(3);
				best_parameter(4) = x_opt(4);
				best_parameter(5) = x_opt(5);
				best_parameter(6) = x_opt(6);
				best_parameter(7) = x_opt(7);
				best_parameter(8) = x_opt(8);
				best_parameter(9) = x_opt(9);
				best_parameter(10) = x_opt(10);
			case model_order(model_id)
				lb = [0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01];
				ub = [5 5 50 50 50 50 1 4 4 10];
				if (best_parameter(11) > 1e-10)
					lb(1) = best_parameter(11) - 1e-10;
					ub(1) = best_parameter(11) + 1e-10;
				end
				if (best_parameter(12) > 1e-10)
					lb(2) = best_parameter(12) - 1e-10;
					ub(2) = best_parameter(12) + 1e-10;
				end
				if (best_parameter(13) > 1e-10)
					lb(3) = best_parameter(13) - 1e-10;
					ub(3) = best_parameter(13) + 1e-10;
				end
				if (best_parameter(14) > 1e-10)
					lb(4) = best_parameter(14) - 1e-10;
					ub(4) = best_parameter(14) + 1e-10;
				end
				if (best_parameter(15) > 1e-10)
					lb(5) = best_parameter(15) - 1e-10;
					ub(5) = best_parameter(15) + 1e-10;
				end
				if (best_parameter(16) > 1e-10)
					lb(6) = best_parameter(16) - 1e-10;
					ub(6) = best_parameter(16) + 1e-10;
				end
				if (best_parameter(17) > 1e-10)
					lb(7) = best_parameter(17) - 1e-10;
					ub(7) = best_parameter(17) + 1e-10;
				end
				if (best_parameter(18) > 1e-10)
					lb(8) = best_parameter(18) - 1e-10;
					ub(8) = best_parameter(18) + 1e-10;
				end
				if (best_parameter(19) > 1e-10)
					lb(9) = best_parameter(19) - 1e-10;
					ub(9) = best_parameter(19) + 1e-10;
				end
				if (best_parameter(20) > 1e-10)
					lb(10) = best_parameter(20) - 1e-10;
					ub(10) = best_parameter(20) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@doi__10_1093_nar_25_6_1203_figure4a_1_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(11) = x_opt(1);
				best_parameter(12) = x_opt(2);
				best_parameter(13) = x_opt(3);
				best_parameter(14) = x_opt(4);
				best_parameter(15) = x_opt(5);
				best_parameter(16) = x_opt(6);
				best_parameter(17) = x_opt(7);
				best_parameter(18) = x_opt(8);
				best_parameter(19) = x_opt(9);
				best_parameter(20) = x_opt(10);
			case model_order(model_id)
				lb = [0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01];
				ub = [5 5 50 50 50 50 1 4 4 10];
				if (best_parameter(21) > 1e-10)
					lb(1) = best_parameter(21) - 1e-10;
					ub(1) = best_parameter(21) + 1e-10;
				end
				if (best_parameter(22) > 1e-10)
					lb(2) = best_parameter(22) - 1e-10;
					ub(2) = best_parameter(22) + 1e-10;
				end
				if (best_parameter(23) > 1e-10)
					lb(3) = best_parameter(23) - 1e-10;
					ub(3) = best_parameter(23) + 1e-10;
				end
				if (best_parameter(24) > 1e-10)
					lb(4) = best_parameter(24) - 1e-10;
					ub(4) = best_parameter(24) + 1e-10;
				end
				if (best_parameter(25) > 1e-10)
					lb(5) = best_parameter(25) - 1e-10;
					ub(5) = best_parameter(25) + 1e-10;
				end
				if (best_parameter(26) > 1e-10)
					lb(6) = best_parameter(26) - 1e-10;
					ub(6) = best_parameter(26) + 1e-10;
				end
				if (best_parameter(27) > 1e-10)
					lb(7) = best_parameter(27) - 1e-10;
					ub(7) = best_parameter(27) + 1e-10;
				end
				if (best_parameter(28) > 1e-10)
					lb(8) = best_parameter(28) - 1e-10;
					ub(8) = best_parameter(28) + 1e-10;
				end
				if (best_parameter(29) > 1e-10)
					lb(9) = best_parameter(29) - 1e-10;
					ub(9) = best_parameter(29) + 1e-10;
				end
				if (best_parameter(30) > 1e-10)
					lb(10) = best_parameter(30) - 1e-10;
					ub(10) = best_parameter(30) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@doi_10_1128_AEM_00791_07_figure3a_1_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(21) = x_opt(1);
				best_parameter(22) = x_opt(2);
				best_parameter(23) = x_opt(3);
				best_parameter(24) = x_opt(4);
				best_parameter(25) = x_opt(5);
				best_parameter(26) = x_opt(6);
				best_parameter(27) = x_opt(7);
				best_parameter(28) = x_opt(8);
				best_parameter(29) = x_opt(9);
				best_parameter(30) = x_opt(10);
			case model_order(model_id)
				lb = [0.03 0.001 0.001 0.002 0.002 0 1 1 0.01];
				ub = [5 50 50 50 50 1 4 4 10];
				if (best_parameter(31) > 1e-10)
					lb(1) = best_parameter(31) - 1e-10;
					ub(1) = best_parameter(31) + 1e-10;
				end
				if (best_parameter(23) > 1e-10)
					lb(2) = best_parameter(23) - 1e-10;
					ub(2) = best_parameter(23) + 1e-10;
				end
				if (best_parameter(24) > 1e-10)
					lb(3) = best_parameter(24) - 1e-10;
					ub(3) = best_parameter(24) + 1e-10;
				end
				if (best_parameter(25) > 1e-10)
					lb(4) = best_parameter(25) - 1e-10;
					ub(4) = best_parameter(25) + 1e-10;
				end
				if (best_parameter(32) > 1e-10)
					lb(5) = best_parameter(32) - 1e-10;
					ub(5) = best_parameter(32) + 1e-10;
				end
				if (best_parameter(27) > 1e-10)
					lb(6) = best_parameter(27) - 1e-10;
					ub(6) = best_parameter(27) + 1e-10;
				end
				if (best_parameter(28) > 1e-10)
					lb(7) = best_parameter(28) - 1e-10;
					ub(7) = best_parameter(28) + 1e-10;
				end
				if (best_parameter(29) > 1e-10)
					lb(8) = best_parameter(29) - 1e-10;
					ub(8) = best_parameter(29) + 1e-10;
				end
				if (best_parameter(33) > 1e-10)
					lb(9) = best_parameter(33) - 1e-10;
					ub(9) = best_parameter(33) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@doi_10_1038_nature12148_figure16a_1_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(31) = x_opt(1);
				best_parameter(23) = x_opt(2);
				best_parameter(24) = x_opt(3);
				best_parameter(25) = x_opt(4);
				best_parameter(32) = x_opt(5);
				best_parameter(27) = x_opt(6);
				best_parameter(28) = x_opt(7);
				best_parameter(29) = x_opt(8);
				best_parameter(33) = x_opt(9);
			case model_order(model_id)
				lb = [0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01];
				ub = [5 5 50 50 50 50 1 4 4 10];
				if (best_parameter(22) > 1e-10)
					lb(1) = best_parameter(22) - 1e-10;
					ub(1) = best_parameter(22) + 1e-10;
				end
				if (best_parameter(34) > 1e-10)
					lb(2) = best_parameter(34) - 1e-10;
					ub(2) = best_parameter(34) + 1e-10;
				end
				if (best_parameter(23) > 1e-10)
					lb(3) = best_parameter(23) - 1e-10;
					ub(3) = best_parameter(23) + 1e-10;
				end
				if (best_parameter(24) > 1e-10)
					lb(4) = best_parameter(24) - 1e-10;
					ub(4) = best_parameter(24) + 1e-10;
				end
				if (best_parameter(25) > 1e-10)
					lb(5) = best_parameter(25) - 1e-10;
					ub(5) = best_parameter(25) + 1e-10;
				end
				if (best_parameter(26) > 1e-10)
					lb(6) = best_parameter(26) - 1e-10;
					ub(6) = best_parameter(26) + 1e-10;
				end
				if (best_parameter(27) > 1e-10)
					lb(7) = best_parameter(27) - 1e-10;
					ub(7) = best_parameter(27) + 1e-10;
				end
				if (best_parameter(28) > 1e-10)
					lb(8) = best_parameter(28) - 1e-10;
					ub(8) = best_parameter(28) + 1e-10;
				end
				if (best_parameter(29) > 1e-10)
					lb(9) = best_parameter(29) - 1e-10;
					ub(9) = best_parameter(29) + 1e-10;
				end
				if (best_parameter(35) > 1e-10)
					lb(10) = best_parameter(35) - 1e-10;
					ub(10) = best_parameter(35) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@doi_10_1038_nature11516_figureS8A_1_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(22) = x_opt(1);
				best_parameter(34) = x_opt(2);
				best_parameter(23) = x_opt(3);
				best_parameter(24) = x_opt(4);
				best_parameter(25) = x_opt(5);
				best_parameter(26) = x_opt(6);
				best_parameter(27) = x_opt(7);
				best_parameter(28) = x_opt(8);
				best_parameter(29) = x_opt(9);
				best_parameter(35) = x_opt(10);
			case model_order(model_id)
				lb = [0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01];
				ub = [5 5 50 50 50 50 1 4 4 10];
				if (best_parameter(2) > 1e-10)
					lb(1) = best_parameter(2) - 1e-10;
					ub(1) = best_parameter(2) + 1e-10;
				end
				if (best_parameter(34) > 1e-10)
					lb(2) = best_parameter(34) - 1e-10;
					ub(2) = best_parameter(34) + 1e-10;
				end
				if (best_parameter(3) > 1e-10)
					lb(3) = best_parameter(3) - 1e-10;
					ub(3) = best_parameter(3) + 1e-10;
				end
				if (best_parameter(36) > 1e-10)
					lb(4) = best_parameter(36) - 1e-10;
					ub(4) = best_parameter(36) + 1e-10;
				end
				if (best_parameter(37) > 1e-10)
					lb(5) = best_parameter(37) - 1e-10;
					ub(5) = best_parameter(37) + 1e-10;
				end
				if (best_parameter(38) > 1e-10)
					lb(6) = best_parameter(38) - 1e-10;
					ub(6) = best_parameter(38) + 1e-10;
				end
				if (best_parameter(39) > 1e-10)
					lb(7) = best_parameter(39) - 1e-10;
					ub(7) = best_parameter(39) + 1e-10;
				end
				if (best_parameter(8) > 1e-10)
					lb(8) = best_parameter(8) - 1e-10;
					ub(8) = best_parameter(8) + 1e-10;
				end
				if (best_parameter(40) > 1e-10)
					lb(9) = best_parameter(40) - 1e-10;
					ub(9) = best_parameter(40) + 1e-10;
				end
				if (best_parameter(35) > 1e-10)
					lb(10) = best_parameter(35) - 1e-10;
					ub(10) = best_parameter(35) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@doi_10_1038_nature11516_figureS8B_2_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(2) = x_opt(1);
				best_parameter(34) = x_opt(2);
				best_parameter(3) = x_opt(3);
				best_parameter(36) = x_opt(4);
				best_parameter(37) = x_opt(5);
				best_parameter(38) = x_opt(6);
				best_parameter(39) = x_opt(7);
				best_parameter(8) = x_opt(8);
				best_parameter(40) = x_opt(9);
				best_parameter(35) = x_opt(10);
			case model_order(model_id)
				lb = [0.03 0.001 0.001 0.002 0.002 0 1 1 0.01];
				ub = [5 50 50 50 50 1 4 4 10];
				if (best_parameter(41) > 1e-10)
					lb(1) = best_parameter(41) - 1e-10;
					ub(1) = best_parameter(41) + 1e-10;
				end
				if (best_parameter(3) > 1e-10)
					lb(2) = best_parameter(3) - 1e-10;
					ub(2) = best_parameter(3) + 1e-10;
				end
				if (best_parameter(36) > 1e-10)
					lb(3) = best_parameter(36) - 1e-10;
					ub(3) = best_parameter(36) + 1e-10;
				end
				if (best_parameter(42) > 1e-10)
					lb(4) = best_parameter(42) - 1e-10;
					ub(4) = best_parameter(42) + 1e-10;
				end
				if (best_parameter(37) > 1e-10)
					lb(5) = best_parameter(37) - 1e-10;
					ub(5) = best_parameter(37) + 1e-10;
				end
				if (best_parameter(39) > 1e-10)
					lb(6) = best_parameter(39) - 1e-10;
					ub(6) = best_parameter(39) + 1e-10;
				end
				if (best_parameter(8) > 1e-10)
					lb(7) = best_parameter(8) - 1e-10;
					ub(7) = best_parameter(8) + 1e-10;
				end
				if (best_parameter(40) > 1e-10)
					lb(8) = best_parameter(40) - 1e-10;
					ub(8) = best_parameter(40) + 1e-10;
				end
				if (best_parameter(43) > 1e-10)
					lb(9) = best_parameter(43) - 1e-10;
					ub(9) = best_parameter(43) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@doi_10_1038_nbt_2401_figureS11_1_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(41) = x_opt(1);
				best_parameter(3) = x_opt(2);
				best_parameter(36) = x_opt(3);
				best_parameter(42) = x_opt(4);
				best_parameter(37) = x_opt(5);
				best_parameter(39) = x_opt(6);
				best_parameter(8) = x_opt(7);
				best_parameter(40) = x_opt(8);
				best_parameter(43) = x_opt(9);
			case model_order(model_id)
				lb = [0.03 0.001 0.001 0.002 0.002 0 1 1 0.01];
				ub = [5 50 50 50 50 1 4 4 10];
				if (best_parameter(41) > 1e-10)
					lb(1) = best_parameter(41) - 1e-10;
					ub(1) = best_parameter(41) + 1e-10;
				end
				if (best_parameter(23) > 1e-10)
					lb(2) = best_parameter(23) - 1e-10;
					ub(2) = best_parameter(23) + 1e-10;
				end
				if (best_parameter(24) > 1e-10)
					lb(3) = best_parameter(24) - 1e-10;
					ub(3) = best_parameter(24) + 1e-10;
				end
				if (best_parameter(42) > 1e-10)
					lb(4) = best_parameter(42) - 1e-10;
					ub(4) = best_parameter(42) + 1e-10;
				end
				if (best_parameter(25) > 1e-10)
					lb(5) = best_parameter(25) - 1e-10;
					ub(5) = best_parameter(25) + 1e-10;
				end
				if (best_parameter(27) > 1e-10)
					lb(6) = best_parameter(27) - 1e-10;
					ub(6) = best_parameter(27) + 1e-10;
				end
				if (best_parameter(28) > 1e-10)
					lb(7) = best_parameter(28) - 1e-10;
					ub(7) = best_parameter(28) + 1e-10;
				end
				if (best_parameter(29) > 1e-10)
					lb(8) = best_parameter(29) - 1e-10;
					ub(8) = best_parameter(29) + 1e-10;
				end
				if (best_parameter(43) > 1e-10)
					lb(9) = best_parameter(43) - 1e-10;
					ub(9) = best_parameter(43) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@doi_10_1038_nbt_2401_figureS13_2_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(41) = x_opt(1);
				best_parameter(23) = x_opt(2);
				best_parameter(24) = x_opt(3);
				best_parameter(42) = x_opt(4);
				best_parameter(25) = x_opt(5);
				best_parameter(27) = x_opt(6);
				best_parameter(28) = x_opt(7);
				best_parameter(29) = x_opt(8);
				best_parameter(43) = x_opt(9);
			case model_order(model_id)
				lb = [0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01];
				ub = [5 5 50 50 50 50 1 4 4 10];
				if (best_parameter(44) > 1e-10)
					lb(1) = best_parameter(44) - 1e-10;
					ub(1) = best_parameter(44) + 1e-10;
				end
				if (best_parameter(45) > 1e-10)
					lb(2) = best_parameter(45) - 1e-10;
					ub(2) = best_parameter(45) + 1e-10;
				end
				if (best_parameter(23) > 1e-10)
					lb(3) = best_parameter(23) - 1e-10;
					ub(3) = best_parameter(23) + 1e-10;
				end
				if (best_parameter(24) > 1e-10)
					lb(4) = best_parameter(24) - 1e-10;
					ub(4) = best_parameter(24) + 1e-10;
				end
				if (best_parameter(46) > 1e-10)
					lb(5) = best_parameter(46) - 1e-10;
					ub(5) = best_parameter(46) + 1e-10;
				end
				if (best_parameter(25) > 1e-10)
					lb(6) = best_parameter(25) - 1e-10;
					ub(6) = best_parameter(25) + 1e-10;
				end
				if (best_parameter(27) > 1e-10)
					lb(7) = best_parameter(27) - 1e-10;
					ub(7) = best_parameter(27) + 1e-10;
				end
				if (best_parameter(28) > 1e-10)
					lb(8) = best_parameter(28) - 1e-10;
					ub(8) = best_parameter(28) + 1e-10;
				end
				if (best_parameter(29) > 1e-10)
					lb(9) = best_parameter(29) - 1e-10;
					ub(9) = best_parameter(29) + 1e-10;
				end
				if (best_parameter(47) > 1e-10)
					lb(10) = best_parameter(47) - 1e-10;
					ub(10) = best_parameter(47) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@doi_10_1038_nature09565_figure1C_1_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(44) = x_opt(1);
				best_parameter(45) = x_opt(2);
				best_parameter(23) = x_opt(3);
				best_parameter(24) = x_opt(4);
				best_parameter(46) = x_opt(5);
				best_parameter(25) = x_opt(6);
				best_parameter(27) = x_opt(7);
				best_parameter(28) = x_opt(8);
				best_parameter(29) = x_opt(9);
				best_parameter(47) = x_opt(10);
			case model_order(model_id)
				lb = [0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01];
				ub = [5 5 50 50 50 50 1 4 4 10];
				if (best_parameter(48) > 1e-10)
					lb(1) = best_parameter(48) - 1e-10;
					ub(1) = best_parameter(48) + 1e-10;
				end
				if (best_parameter(44) > 1e-10)
					lb(2) = best_parameter(44) - 1e-10;
					ub(2) = best_parameter(44) + 1e-10;
				end
				if (best_parameter(13) > 1e-10)
					lb(3) = best_parameter(13) - 1e-10;
					ub(3) = best_parameter(13) + 1e-10;
				end
				if (best_parameter(14) > 1e-10)
					lb(4) = best_parameter(14) - 1e-10;
					ub(4) = best_parameter(14) + 1e-10;
				end
				if (best_parameter(46) > 1e-10)
					lb(5) = best_parameter(46) - 1e-10;
					ub(5) = best_parameter(46) + 1e-10;
				end
				if (best_parameter(15) > 1e-10)
					lb(6) = best_parameter(15) - 1e-10;
					ub(6) = best_parameter(15) + 1e-10;
				end
				if (best_parameter(17) > 1e-10)
					lb(7) = best_parameter(17) - 1e-10;
					ub(7) = best_parameter(17) + 1e-10;
				end
				if (best_parameter(18) > 1e-10)
					lb(8) = best_parameter(18) - 1e-10;
					ub(8) = best_parameter(18) + 1e-10;
				end
				if (best_parameter(19) > 1e-10)
					lb(9) = best_parameter(19) - 1e-10;
					ub(9) = best_parameter(19) + 1e-10;
				end
				if (best_parameter(47) > 1e-10)
					lb(10) = best_parameter(47) - 1e-10;
					ub(10) = best_parameter(47) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@doi_10_1038_nature09565_figureS1C_2_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(48) = x_opt(1);
				best_parameter(44) = x_opt(2);
				best_parameter(13) = x_opt(3);
				best_parameter(14) = x_opt(4);
				best_parameter(46) = x_opt(5);
				best_parameter(15) = x_opt(6);
				best_parameter(17) = x_opt(7);
				best_parameter(18) = x_opt(8);
				best_parameter(19) = x_opt(9);
				best_parameter(47) = x_opt(10);
			case model_order(model_id)
				lb = [0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01];
				ub = [5 5 50 50 50 50 1 4 4 10];
				if (best_parameter(49) > 1e-10)
					lb(1) = best_parameter(49) - 1e-10;
					ub(1) = best_parameter(49) + 1e-10;
				end
				if (best_parameter(2) > 1e-10)
					lb(2) = best_parameter(2) - 1e-10;
					ub(2) = best_parameter(2) + 1e-10;
				end
				if (best_parameter(3) > 1e-10)
					lb(3) = best_parameter(3) - 1e-10;
					ub(3) = best_parameter(3) + 1e-10;
				end
				if (best_parameter(50) > 1e-10)
					lb(4) = best_parameter(50) - 1e-10;
					ub(4) = best_parameter(50) + 1e-10;
				end
				if (best_parameter(6) > 1e-10)
					lb(5) = best_parameter(6) - 1e-10;
					ub(5) = best_parameter(6) + 1e-10;
				end
				if (best_parameter(51) > 1e-10)
					lb(6) = best_parameter(51) - 1e-10;
					ub(6) = best_parameter(51) + 1e-10;
				end
				if (best_parameter(52) > 1e-10)
					lb(7) = best_parameter(52) - 1e-10;
					ub(7) = best_parameter(52) + 1e-10;
				end
				if (best_parameter(8) > 1e-10)
					lb(8) = best_parameter(8) - 1e-10;
					ub(8) = best_parameter(8) + 1e-10;
				end
				if (best_parameter(53) > 1e-10)
					lb(9) = best_parameter(53) - 1e-10;
					ub(9) = best_parameter(53) + 1e-10;
				end
				if (best_parameter(54) > 1e-10)
					lb(10) = best_parameter(54) - 1e-10;
					ub(10) = best_parameter(54) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@doi_10_1038_nbt_2149_Supplement_page17_1_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(49) = x_opt(1);
				best_parameter(2) = x_opt(2);
				best_parameter(3) = x_opt(3);
				best_parameter(50) = x_opt(4);
				best_parameter(6) = x_opt(5);
				best_parameter(51) = x_opt(6);
				best_parameter(52) = x_opt(7);
				best_parameter(8) = x_opt(8);
				best_parameter(53) = x_opt(9);
				best_parameter(54) = x_opt(10);
			case model_order(model_id)
				lb = [0.03 0.002 0.01];
				ub = [5 50 10];
				if (best_parameter(48) > 1e-10)
					lb(1) = best_parameter(48) - 1e-10;
					ub(1) = best_parameter(48) + 1e-10;
				end
				if (best_parameter(42) > 1e-10)
					lb(2) = best_parameter(42) - 1e-10;
					ub(2) = best_parameter(42) + 1e-10;
				end
				if (best_parameter(55) > 1e-10)
					lb(3) = best_parameter(55) - 1e-10;
					ub(3) = best_parameter(55) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@doi_10_1186_1754_1611_3_4_figure3_1_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(48) = x_opt(1);
				best_parameter(42) = x_opt(2);
				best_parameter(55) = x_opt(3);
			case model_order(model_id)
				lb = [0.03 0.002 0.01];
				ub = [5 50 10];
				if (best_parameter(48) > 1e-10)
					lb(1) = best_parameter(48) - 1e-10;
					ub(1) = best_parameter(48) + 1e-10;
				end
				if (best_parameter(56) > 1e-10)
					lb(2) = best_parameter(56) - 1e-10;
					ub(2) = best_parameter(56) + 1e-10;
				end
				if (best_parameter(55) > 1e-10)
					lb(3) = best_parameter(55) - 1e-10;
					ub(3) = best_parameter(55) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@doi_10_1186_1754_1611_3_4_figure3_2_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(48) = x_opt(1);
				best_parameter(56) = x_opt(2);
				best_parameter(55) = x_opt(3);
			case model_order(model_id)
				lb = [0.03 0.002 0.01];
				ub = [5 50 10];
				if (best_parameter(48) > 1e-10)
					lb(1) = best_parameter(48) - 1e-10;
					ub(1) = best_parameter(48) + 1e-10;
				end
				if (best_parameter(57) > 1e-10)
					lb(2) = best_parameter(57) - 1e-10;
					ub(2) = best_parameter(57) + 1e-10;
				end
				if (best_parameter(55) > 1e-10)
					lb(3) = best_parameter(55) - 1e-10;
					ub(3) = best_parameter(55) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@doi_10_1186_1754_1611_3_4_figure3_3_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(48) = x_opt(1);
				best_parameter(57) = x_opt(2);
				best_parameter(55) = x_opt(3);
			case model_order(model_id)
				lb = [0.03 0.002 0.01];
				ub = [5 50 10];
				if (best_parameter(48) > 1e-10)
					lb(1) = best_parameter(48) - 1e-10;
					ub(1) = best_parameter(48) + 1e-10;
				end
				if (best_parameter(58) > 1e-10)
					lb(2) = best_parameter(58) - 1e-10;
					ub(2) = best_parameter(58) + 1e-10;
				end
				if (best_parameter(55) > 1e-10)
					lb(3) = best_parameter(55) - 1e-10;
					ub(3) = best_parameter(55) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@doi_10_1186_1754_1611_3_4_figure3_4_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(48) = x_opt(1);
				best_parameter(58) = x_opt(2);
				best_parameter(55) = x_opt(3);
			case model_order(model_id)
				lb = [0.03 0.002 0.01];
				ub = [5 50 10];
				if (best_parameter(48) > 1e-10)
					lb(1) = best_parameter(48) - 1e-10;
					ub(1) = best_parameter(48) + 1e-10;
				end
				if (best_parameter(59) > 1e-10)
					lb(2) = best_parameter(59) - 1e-10;
					ub(2) = best_parameter(59) + 1e-10;
				end
				if (best_parameter(55) > 1e-10)
					lb(3) = best_parameter(55) - 1e-10;
					ub(3) = best_parameter(55) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@doi_10_1186_1754_1611_3_4_figure3_5_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(48) = x_opt(1);
				best_parameter(59) = x_opt(2);
				best_parameter(55) = x_opt(3);
			case model_order(model_id)
				lb = [0.03 0.002 0.01];
				ub = [5 50 10];
				if (best_parameter(48) > 1e-10)
					lb(1) = best_parameter(48) - 1e-10;
					ub(1) = best_parameter(48) + 1e-10;
				end
				if (best_parameter(15) > 1e-10)
					lb(2) = best_parameter(15) - 1e-10;
					ub(2) = best_parameter(15) + 1e-10;
				end
				if (best_parameter(55) > 1e-10)
					lb(3) = best_parameter(55) - 1e-10;
					ub(3) = best_parameter(55) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@doi_10_1186_1754_1611_3_4_figure3_6_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(48) = x_opt(1);
				best_parameter(15) = x_opt(2);
				best_parameter(55) = x_opt(3);
			case model_order(model_id)
				lb = [0.03 0.002 0.01];
				ub = [5 50 10];
				if (best_parameter(48) > 1e-10)
					lb(1) = best_parameter(48) - 1e-10;
					ub(1) = best_parameter(48) + 1e-10;
				end
				if (best_parameter(32) > 1e-10)
					lb(2) = best_parameter(32) - 1e-10;
					ub(2) = best_parameter(32) + 1e-10;
				end
				if (best_parameter(55) > 1e-10)
					lb(3) = best_parameter(55) - 1e-10;
					ub(3) = best_parameter(55) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@doi_10_1186_1754_1611_3_4_figure3_7_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(48) = x_opt(1);
				best_parameter(32) = x_opt(2);
				best_parameter(55) = x_opt(3);
			case model_order(model_id)
				lb = [0.03 0.002 0.01];
				ub = [5 50 10];
				if (best_parameter(48) > 1e-10)
					lb(1) = best_parameter(48) - 1e-10;
					ub(1) = best_parameter(48) + 1e-10;
				end
				if (best_parameter(46) > 1e-10)
					lb(2) = best_parameter(46) - 1e-10;
					ub(2) = best_parameter(46) + 1e-10;
				end
				if (best_parameter(60) > 1e-10)
					lb(3) = best_parameter(60) - 1e-10;
					ub(3) = best_parameter(60) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@doi__10_1073_pnas_1301301110_sd03_xls_1_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(48) = x_opt(1);
				best_parameter(46) = x_opt(2);
				best_parameter(60) = x_opt(3);
			case model_order(model_id)
				lb = [0.03 0.002 0.01];
				ub = [5 50 10];
				if (best_parameter(48) > 1e-10)
					lb(1) = best_parameter(48) - 1e-10;
					ub(1) = best_parameter(48) + 1e-10;
				end
				if (best_parameter(42) > 1e-10)
					lb(2) = best_parameter(42) - 1e-10;
					ub(2) = best_parameter(42) + 1e-10;
				end
				if (best_parameter(60) > 1e-10)
					lb(3) = best_parameter(60) - 1e-10;
					ub(3) = best_parameter(60) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@doi__10_1073_pnas_1301301110_sd03_xls_2_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(48) = x_opt(1);
				best_parameter(42) = x_opt(2);
				best_parameter(60) = x_opt(3);
			case model_order(model_id)
				lb = [0.03 0.002 0.01];
				ub = [5 50 10];
				if (best_parameter(31) > 1e-10)
					lb(1) = best_parameter(31) - 1e-10;
					ub(1) = best_parameter(31) + 1e-10;
				end
				if (best_parameter(42) > 1e-10)
					lb(2) = best_parameter(42) - 1e-10;
					ub(2) = best_parameter(42) + 1e-10;
				end
				if (best_parameter(61) > 1e-10)
					lb(3) = best_parameter(61) - 1e-10;
					ub(3) = best_parameter(61) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@doi_10_1093_nar_gku1388_figure2A_1_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(31) = x_opt(1);
				best_parameter(42) = x_opt(2);
				best_parameter(61) = x_opt(3);
			case model_order(model_id)
				lb = [0.03 0.002 0.01];
				ub = [5 50 10];
				if (best_parameter(31) > 1e-10)
					lb(1) = best_parameter(31) - 1e-10;
					ub(1) = best_parameter(31) + 1e-10;
				end
				if (best_parameter(62) > 1e-10)
					lb(2) = best_parameter(62) - 1e-10;
					ub(2) = best_parameter(62) + 1e-10;
				end
				if (best_parameter(61) > 1e-10)
					lb(3) = best_parameter(61) - 1e-10;
					ub(3) = best_parameter(61) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@doi_10_1093_nar_gku1388_figure2A_2_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(31) = x_opt(1);
				best_parameter(62) = x_opt(2);
				best_parameter(61) = x_opt(3);
			case model_order(model_id)
				lb = [0.03 0.002 0.01];
				ub = [5 50 10];
				if (best_parameter(31) > 1e-10)
					lb(1) = best_parameter(31) - 1e-10;
					ub(1) = best_parameter(31) + 1e-10;
				end
				if (best_parameter(63) > 1e-10)
					lb(2) = best_parameter(63) - 1e-10;
					ub(2) = best_parameter(63) + 1e-10;
				end
				if (best_parameter(61) > 1e-10)
					lb(3) = best_parameter(61) - 1e-10;
					ub(3) = best_parameter(61) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@doi_10_1093_nar_gku1388_figure2A_3_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(31) = x_opt(1);
				best_parameter(63) = x_opt(2);
				best_parameter(61) = x_opt(3);
			case model_order(model_id)
				lb = [0.03 0.002 0.01];
				ub = [5 50 10];
				if (best_parameter(31) > 1e-10)
					lb(1) = best_parameter(31) - 1e-10;
					ub(1) = best_parameter(31) + 1e-10;
				end
				if (best_parameter(64) > 1e-10)
					lb(2) = best_parameter(64) - 1e-10;
					ub(2) = best_parameter(64) + 1e-10;
				end
				if (best_parameter(61) > 1e-10)
					lb(3) = best_parameter(61) - 1e-10;
					ub(3) = best_parameter(61) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@doi_10_1093_nar_gku1388_figure2A_4_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(31) = x_opt(1);
				best_parameter(64) = x_opt(2);
				best_parameter(61) = x_opt(3);
			case model_order(model_id)
				lb = [0.03 0.002 0.01];
				ub = [5 50 10];
				if (best_parameter(31) > 1e-10)
					lb(1) = best_parameter(31) - 1e-10;
					ub(1) = best_parameter(31) + 1e-10;
				end
				if (best_parameter(65) > 1e-10)
					lb(2) = best_parameter(65) - 1e-10;
					ub(2) = best_parameter(65) + 1e-10;
				end
				if (best_parameter(61) > 1e-10)
					lb(3) = best_parameter(61) - 1e-10;
					ub(3) = best_parameter(61) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@doi_10_1093_nar_gku1388_figure2A_5_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(31) = x_opt(1);
				best_parameter(65) = x_opt(2);
				best_parameter(61) = x_opt(3);
			case model_order(model_id)
				lb = [0.03 0.002 0.01];
				ub = [5 50 10];
				if (best_parameter(31) > 1e-10)
					lb(1) = best_parameter(31) - 1e-10;
					ub(1) = best_parameter(31) + 1e-10;
				end
				if (best_parameter(46) > 1e-10)
					lb(2) = best_parameter(46) - 1e-10;
					ub(2) = best_parameter(46) + 1e-10;
				end
				if (best_parameter(61) > 1e-10)
					lb(3) = best_parameter(61) - 1e-10;
					ub(3) = best_parameter(61) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@doi_10_1093_nar_gku1388_figure2A_6_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(31) = x_opt(1);
				best_parameter(46) = x_opt(2);
				best_parameter(61) = x_opt(3);
			case model_order(model_id)
				lb = [0.03 0.001 0.001 0.002 0.002 0 1 1 0.01];
				ub = [5 50 50 50 50 1 4 4 10];
				if (best_parameter(31) > 1e-10)
					lb(1) = best_parameter(31) - 1e-10;
					ub(1) = best_parameter(31) + 1e-10;
				end
				if (best_parameter(13) > 1e-10)
					lb(2) = best_parameter(13) - 1e-10;
					ub(2) = best_parameter(13) + 1e-10;
				end
				if (best_parameter(66) > 1e-10)
					lb(3) = best_parameter(66) - 1e-10;
					ub(3) = best_parameter(66) + 1e-10;
				end
				if (best_parameter(46) > 1e-10)
					lb(4) = best_parameter(46) - 1e-10;
					ub(4) = best_parameter(46) + 1e-10;
				end
				if (best_parameter(67) > 1e-10)
					lb(5) = best_parameter(67) - 1e-10;
					ub(5) = best_parameter(67) + 1e-10;
				end
				if (best_parameter(68) > 1e-10)
					lb(6) = best_parameter(68) - 1e-10;
					ub(6) = best_parameter(68) + 1e-10;
				end
				if (best_parameter(18) > 1e-10)
					lb(7) = best_parameter(18) - 1e-10;
					ub(7) = best_parameter(18) + 1e-10;
				end
				if (best_parameter(69) > 1e-10)
					lb(8) = best_parameter(69) - 1e-10;
					ub(8) = best_parameter(69) + 1e-10;
				end
				if (best_parameter(61) > 1e-10)
					lb(9) = best_parameter(61) - 1e-10;
					ub(9) = best_parameter(61) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@doi_10_1093_nar_gku1388_figure2B_7_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(31) = x_opt(1);
				best_parameter(13) = x_opt(2);
				best_parameter(66) = x_opt(3);
				best_parameter(46) = x_opt(4);
				best_parameter(67) = x_opt(5);
				best_parameter(68) = x_opt(6);
				best_parameter(18) = x_opt(7);
				best_parameter(69) = x_opt(8);
				best_parameter(61) = x_opt(9);
			case model_order(model_id)
				lb = [0.03 0.001 0.001 0.002 0.002 0 1 1 0.01];
				ub = [5 50 50 50 50 1 4 4 10];
				if (best_parameter(31) > 1e-10)
					lb(1) = best_parameter(31) - 1e-10;
					ub(1) = best_parameter(31) + 1e-10;
				end
				if (best_parameter(13) > 1e-10)
					lb(2) = best_parameter(13) - 1e-10;
					ub(2) = best_parameter(13) + 1e-10;
				end
				if (best_parameter(66) > 1e-10)
					lb(3) = best_parameter(66) - 1e-10;
					ub(3) = best_parameter(66) + 1e-10;
				end
				if (best_parameter(65) > 1e-10)
					lb(4) = best_parameter(65) - 1e-10;
					ub(4) = best_parameter(65) + 1e-10;
				end
				if (best_parameter(67) > 1e-10)
					lb(5) = best_parameter(67) - 1e-10;
					ub(5) = best_parameter(67) + 1e-10;
				end
				if (best_parameter(68) > 1e-10)
					lb(6) = best_parameter(68) - 1e-10;
					ub(6) = best_parameter(68) + 1e-10;
				end
				if (best_parameter(18) > 1e-10)
					lb(7) = best_parameter(18) - 1e-10;
					ub(7) = best_parameter(18) + 1e-10;
				end
				if (best_parameter(69) > 1e-10)
					lb(8) = best_parameter(69) - 1e-10;
					ub(8) = best_parameter(69) + 1e-10;
				end
				if (best_parameter(61) > 1e-10)
					lb(9) = best_parameter(61) - 1e-10;
					ub(9) = best_parameter(61) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@doi_10_1093_nar_gku1388_figure2B_8_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(31) = x_opt(1);
				best_parameter(13) = x_opt(2);
				best_parameter(66) = x_opt(3);
				best_parameter(65) = x_opt(4);
				best_parameter(67) = x_opt(5);
				best_parameter(68) = x_opt(6);
				best_parameter(18) = x_opt(7);
				best_parameter(69) = x_opt(8);
				best_parameter(61) = x_opt(9);
			case model_order(model_id)
				lb = [0.03 0.001 0.001 0.002 0.002 0 1 1 0.01];
				ub = [5 50 50 50 50 1 4 4 10];
				if (best_parameter(31) > 1e-10)
					lb(1) = best_parameter(31) - 1e-10;
					ub(1) = best_parameter(31) + 1e-10;
				end
				if (best_parameter(13) > 1e-10)
					lb(2) = best_parameter(13) - 1e-10;
					ub(2) = best_parameter(13) + 1e-10;
				end
				if (best_parameter(66) > 1e-10)
					lb(3) = best_parameter(66) - 1e-10;
					ub(3) = best_parameter(66) + 1e-10;
				end
				if (best_parameter(64) > 1e-10)
					lb(4) = best_parameter(64) - 1e-10;
					ub(4) = best_parameter(64) + 1e-10;
				end
				if (best_parameter(67) > 1e-10)
					lb(5) = best_parameter(67) - 1e-10;
					ub(5) = best_parameter(67) + 1e-10;
				end
				if (best_parameter(68) > 1e-10)
					lb(6) = best_parameter(68) - 1e-10;
					ub(6) = best_parameter(68) + 1e-10;
				end
				if (best_parameter(18) > 1e-10)
					lb(7) = best_parameter(18) - 1e-10;
					ub(7) = best_parameter(18) + 1e-10;
				end
				if (best_parameter(69) > 1e-10)
					lb(8) = best_parameter(69) - 1e-10;
					ub(8) = best_parameter(69) + 1e-10;
				end
				if (best_parameter(61) > 1e-10)
					lb(9) = best_parameter(61) - 1e-10;
					ub(9) = best_parameter(61) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@doi_10_1093_nar_gku1388_figure2B_9_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(31) = x_opt(1);
				best_parameter(13) = x_opt(2);
				best_parameter(66) = x_opt(3);
				best_parameter(64) = x_opt(4);
				best_parameter(67) = x_opt(5);
				best_parameter(68) = x_opt(6);
				best_parameter(18) = x_opt(7);
				best_parameter(69) = x_opt(8);
				best_parameter(61) = x_opt(9);
			case model_order(model_id)
				lb = [0.03 0.001 0.001 0.002 0.002 0 1 1 0.01];
				ub = [5 50 50 50 50 1 4 4 10];
				if (best_parameter(31) > 1e-10)
					lb(1) = best_parameter(31) - 1e-10;
					ub(1) = best_parameter(31) + 1e-10;
				end
				if (best_parameter(13) > 1e-10)
					lb(2) = best_parameter(13) - 1e-10;
					ub(2) = best_parameter(13) + 1e-10;
				end
				if (best_parameter(66) > 1e-10)
					lb(3) = best_parameter(66) - 1e-10;
					ub(3) = best_parameter(66) + 1e-10;
				end
				if (best_parameter(63) > 1e-10)
					lb(4) = best_parameter(63) - 1e-10;
					ub(4) = best_parameter(63) + 1e-10;
				end
				if (best_parameter(67) > 1e-10)
					lb(5) = best_parameter(67) - 1e-10;
					ub(5) = best_parameter(67) + 1e-10;
				end
				if (best_parameter(68) > 1e-10)
					lb(6) = best_parameter(68) - 1e-10;
					ub(6) = best_parameter(68) + 1e-10;
				end
				if (best_parameter(18) > 1e-10)
					lb(7) = best_parameter(18) - 1e-10;
					ub(7) = best_parameter(18) + 1e-10;
				end
				if (best_parameter(69) > 1e-10)
					lb(8) = best_parameter(69) - 1e-10;
					ub(8) = best_parameter(69) + 1e-10;
				end
				if (best_parameter(61) > 1e-10)
					lb(9) = best_parameter(61) - 1e-10;
					ub(9) = best_parameter(61) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@doi_10_1093_nar_gku1388_figure2B_10_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(31) = x_opt(1);
				best_parameter(13) = x_opt(2);
				best_parameter(66) = x_opt(3);
				best_parameter(63) = x_opt(4);
				best_parameter(67) = x_opt(5);
				best_parameter(68) = x_opt(6);
				best_parameter(18) = x_opt(7);
				best_parameter(69) = x_opt(8);
				best_parameter(61) = x_opt(9);
			case model_order(model_id)
				lb = [0.03 0.001 0.001 0.002 0.002 0 1 1 0.01];
				ub = [5 50 50 50 50 1 4 4 10];
				if (best_parameter(31) > 1e-10)
					lb(1) = best_parameter(31) - 1e-10;
					ub(1) = best_parameter(31) + 1e-10;
				end
				if (best_parameter(13) > 1e-10)
					lb(2) = best_parameter(13) - 1e-10;
					ub(2) = best_parameter(13) + 1e-10;
				end
				if (best_parameter(66) > 1e-10)
					lb(3) = best_parameter(66) - 1e-10;
					ub(3) = best_parameter(66) + 1e-10;
				end
				if (best_parameter(62) > 1e-10)
					lb(4) = best_parameter(62) - 1e-10;
					ub(4) = best_parameter(62) + 1e-10;
				end
				if (best_parameter(67) > 1e-10)
					lb(5) = best_parameter(67) - 1e-10;
					ub(5) = best_parameter(67) + 1e-10;
				end
				if (best_parameter(68) > 1e-10)
					lb(6) = best_parameter(68) - 1e-10;
					ub(6) = best_parameter(68) + 1e-10;
				end
				if (best_parameter(18) > 1e-10)
					lb(7) = best_parameter(18) - 1e-10;
					ub(7) = best_parameter(18) + 1e-10;
				end
				if (best_parameter(69) > 1e-10)
					lb(8) = best_parameter(69) - 1e-10;
					ub(8) = best_parameter(69) + 1e-10;
				end
				if (best_parameter(61) > 1e-10)
					lb(9) = best_parameter(61) - 1e-10;
					ub(9) = best_parameter(61) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@doi_10_1093_nar_gku1388_figure2B_11_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(31) = x_opt(1);
				best_parameter(13) = x_opt(2);
				best_parameter(66) = x_opt(3);
				best_parameter(62) = x_opt(4);
				best_parameter(67) = x_opt(5);
				best_parameter(68) = x_opt(6);
				best_parameter(18) = x_opt(7);
				best_parameter(69) = x_opt(8);
				best_parameter(61) = x_opt(9);
			case model_order(model_id)
				lb = [0.03 0.001 0.001 0.002 0.002 0 1 1 0.01];
				ub = [5 50 50 50 50 1 4 4 10];
				if (best_parameter(31) > 1e-10)
					lb(1) = best_parameter(31) - 1e-10;
					ub(1) = best_parameter(31) + 1e-10;
				end
				if (best_parameter(13) > 1e-10)
					lb(2) = best_parameter(13) - 1e-10;
					ub(2) = best_parameter(13) + 1e-10;
				end
				if (best_parameter(66) > 1e-10)
					lb(3) = best_parameter(66) - 1e-10;
					ub(3) = best_parameter(66) + 1e-10;
				end
				if (best_parameter(42) > 1e-10)
					lb(4) = best_parameter(42) - 1e-10;
					ub(4) = best_parameter(42) + 1e-10;
				end
				if (best_parameter(67) > 1e-10)
					lb(5) = best_parameter(67) - 1e-10;
					ub(5) = best_parameter(67) + 1e-10;
				end
				if (best_parameter(68) > 1e-10)
					lb(6) = best_parameter(68) - 1e-10;
					ub(6) = best_parameter(68) + 1e-10;
				end
				if (best_parameter(18) > 1e-10)
					lb(7) = best_parameter(18) - 1e-10;
					ub(7) = best_parameter(18) + 1e-10;
				end
				if (best_parameter(69) > 1e-10)
					lb(8) = best_parameter(69) - 1e-10;
					ub(8) = best_parameter(69) + 1e-10;
				end
				if (best_parameter(61) > 1e-10)
					lb(9) = best_parameter(61) - 1e-10;
					ub(9) = best_parameter(61) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@doi_10_1093_nar_gku1388_figure2B_12_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(31) = x_opt(1);
				best_parameter(13) = x_opt(2);
				best_parameter(66) = x_opt(3);
				best_parameter(42) = x_opt(4);
				best_parameter(67) = x_opt(5);
				best_parameter(68) = x_opt(6);
				best_parameter(18) = x_opt(7);
				best_parameter(69) = x_opt(8);
				best_parameter(61) = x_opt(9);
			case model_order(model_id)
				lb = [0.03 0.002 0.01];
				ub = [5 50 10];
				if (best_parameter(31) > 1e-10)
					lb(1) = best_parameter(31) - 1e-10;
					ub(1) = best_parameter(31) + 1e-10;
				end
				if (best_parameter(70) > 1e-10)
					lb(2) = best_parameter(70) - 1e-10;
					ub(2) = best_parameter(70) + 1e-10;
				end
				if (best_parameter(71) > 1e-10)
					lb(3) = best_parameter(71) - 1e-10;
					ub(3) = best_parameter(71) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@doi_10_1093_nar_gku593_figure3D_1_1_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(31) = x_opt(1);
				best_parameter(70) = x_opt(2);
				best_parameter(71) = x_opt(3);
			case model_order(model_id)
				lb = [0.03 0.002 0.01];
				ub = [5 50 10];
				if (best_parameter(31) > 1e-10)
					lb(1) = best_parameter(31) - 1e-10;
					ub(1) = best_parameter(31) + 1e-10;
				end
				if (best_parameter(65) > 1e-10)
					lb(2) = best_parameter(65) - 1e-10;
					ub(2) = best_parameter(65) + 1e-10;
				end
				if (best_parameter(71) > 1e-10)
					lb(3) = best_parameter(71) - 1e-10;
					ub(3) = best_parameter(71) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@doi_10_1093_nar_gku593_figure3D_2_2_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(31) = x_opt(1);
				best_parameter(65) = x_opt(2);
				best_parameter(71) = x_opt(3);
			case model_order(model_id)
				lb = [0.03 0.002 0.01];
				ub = [5 50 10];
				if (best_parameter(31) > 1e-10)
					lb(1) = best_parameter(31) - 1e-10;
					ub(1) = best_parameter(31) + 1e-10;
				end
				if (best_parameter(56) > 1e-10)
					lb(2) = best_parameter(56) - 1e-10;
					ub(2) = best_parameter(56) + 1e-10;
				end
				if (best_parameter(71) > 1e-10)
					lb(3) = best_parameter(71) - 1e-10;
					ub(3) = best_parameter(71) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@doi_10_1093_nar_gku593_figure3D_3_3_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(31) = x_opt(1);
				best_parameter(56) = x_opt(2);
				best_parameter(71) = x_opt(3);
			case model_order(model_id)
				lb = [0.03 0.002 0.01];
				ub = [5 50 10];
				if (best_parameter(31) > 1e-10)
					lb(1) = best_parameter(31) - 1e-10;
					ub(1) = best_parameter(31) + 1e-10;
				end
				if (best_parameter(64) > 1e-10)
					lb(2) = best_parameter(64) - 1e-10;
					ub(2) = best_parameter(64) + 1e-10;
				end
				if (best_parameter(71) > 1e-10)
					lb(3) = best_parameter(71) - 1e-10;
					ub(3) = best_parameter(71) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@doi_10_1093_nar_gku593_figure3D_4_4_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(31) = x_opt(1);
				best_parameter(64) = x_opt(2);
				best_parameter(71) = x_opt(3);
			case model_order(model_id)
				lb = [0.03 0.002 0.01];
				ub = [5 50 10];
				if (best_parameter(31) > 1e-10)
					lb(1) = best_parameter(31) - 1e-10;
					ub(1) = best_parameter(31) + 1e-10;
				end
				if (best_parameter(63) > 1e-10)
					lb(2) = best_parameter(63) - 1e-10;
					ub(2) = best_parameter(63) + 1e-10;
				end
				if (best_parameter(71) > 1e-10)
					lb(3) = best_parameter(71) - 1e-10;
					ub(3) = best_parameter(71) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@doi_10_1093_nar_gku593_figure3D_5_5_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(31) = x_opt(1);
				best_parameter(63) = x_opt(2);
				best_parameter(71) = x_opt(3);
			case model_order(model_id)
				lb = [0.03 0.002 0.01];
				ub = [5 50 10];
				if (best_parameter(31) > 1e-10)
					lb(1) = best_parameter(31) - 1e-10;
					ub(1) = best_parameter(31) + 1e-10;
				end
				if (best_parameter(62) > 1e-10)
					lb(2) = best_parameter(62) - 1e-10;
					ub(2) = best_parameter(62) + 1e-10;
				end
				if (best_parameter(71) > 1e-10)
					lb(3) = best_parameter(71) - 1e-10;
					ub(3) = best_parameter(71) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@doi_10_1093_nar_gku593_figure3D_6_6_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(31) = x_opt(1);
				best_parameter(62) = x_opt(2);
				best_parameter(71) = x_opt(3);
			case model_order(model_id)
				lb = [0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01];
				ub = [5 5 50 50 50 50 1 4 4 10];
				if (best_parameter(31) > 1e-10)
					lb(1) = best_parameter(31) - 1e-10;
					ub(1) = best_parameter(31) + 1e-10;
				end
				if (best_parameter(22) > 1e-10)
					lb(2) = best_parameter(22) - 1e-10;
					ub(2) = best_parameter(22) + 1e-10;
				end
				if (best_parameter(23) > 1e-10)
					lb(3) = best_parameter(23) - 1e-10;
					ub(3) = best_parameter(23) + 1e-10;
				end
				if (best_parameter(24) > 1e-10)
					lb(4) = best_parameter(24) - 1e-10;
					ub(4) = best_parameter(24) + 1e-10;
				end
				if (best_parameter(25) > 1e-10)
					lb(5) = best_parameter(25) - 1e-10;
					ub(5) = best_parameter(25) + 1e-10;
				end
				if (best_parameter(26) > 1e-10)
					lb(6) = best_parameter(26) - 1e-10;
					ub(6) = best_parameter(26) + 1e-10;
				end
				if (best_parameter(27) > 1e-10)
					lb(7) = best_parameter(27) - 1e-10;
					ub(7) = best_parameter(27) + 1e-10;
				end
				if (best_parameter(28) > 1e-10)
					lb(8) = best_parameter(28) - 1e-10;
					ub(8) = best_parameter(28) + 1e-10;
				end
				if (best_parameter(29) > 1e-10)
					lb(9) = best_parameter(29) - 1e-10;
					ub(9) = best_parameter(29) + 1e-10;
				end
				if (best_parameter(71) > 1e-10)
					lb(10) = best_parameter(71) - 1e-10;
					ub(10) = best_parameter(71) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@doi_10_1093_nar_gku593_figureS6_7_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(31) = x_opt(1);
				best_parameter(22) = x_opt(2);
				best_parameter(23) = x_opt(3);
				best_parameter(24) = x_opt(4);
				best_parameter(25) = x_opt(5);
				best_parameter(26) = x_opt(6);
				best_parameter(27) = x_opt(7);
				best_parameter(28) = x_opt(8);
				best_parameter(29) = x_opt(9);
				best_parameter(71) = x_opt(10);
			case model_order(model_id)
				lb = [0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01];
				ub = [5 5 50 50 50 50 1 4 4 10];
				if (best_parameter(2) > 1e-10)
					lb(1) = best_parameter(2) - 1e-10;
					ub(1) = best_parameter(2) + 1e-10;
				end
				if (best_parameter(12) > 1e-10)
					lb(2) = best_parameter(12) - 1e-10;
					ub(2) = best_parameter(12) + 1e-10;
				end
				if (best_parameter(3) > 1e-10)
					lb(3) = best_parameter(3) - 1e-10;
					ub(3) = best_parameter(3) + 1e-10;
				end
				if (best_parameter(4) > 1e-10)
					lb(4) = best_parameter(4) - 1e-10;
					ub(4) = best_parameter(4) + 1e-10;
				end
				if (best_parameter(5) > 1e-10)
					lb(5) = best_parameter(5) - 1e-10;
					ub(5) = best_parameter(5) + 1e-10;
				end
				if (best_parameter(38) > 1e-10)
					lb(6) = best_parameter(38) - 1e-10;
					ub(6) = best_parameter(38) + 1e-10;
				end
				if (best_parameter(7) > 1e-10)
					lb(7) = best_parameter(7) - 1e-10;
					ub(7) = best_parameter(7) + 1e-10;
				end
				if (best_parameter(8) > 1e-10)
					lb(8) = best_parameter(8) - 1e-10;
					ub(8) = best_parameter(8) + 1e-10;
				end
				if (best_parameter(9) > 1e-10)
					lb(9) = best_parameter(9) - 1e-10;
					ub(9) = best_parameter(9) + 1e-10;
				end
				if (best_parameter(72) > 1e-10)
					lb(10) = best_parameter(72) - 1e-10;
					ub(10) = best_parameter(72) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(2) = x_opt(1);
				best_parameter(12) = x_opt(2);
				best_parameter(3) = x_opt(3);
				best_parameter(4) = x_opt(4);
				best_parameter(5) = x_opt(5);
				best_parameter(38) = x_opt(6);
				best_parameter(7) = x_opt(7);
				best_parameter(8) = x_opt(8);
				best_parameter(9) = x_opt(9);
				best_parameter(72) = x_opt(10);
			case model_order(model_id)
				lb = [0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01];
				ub = [5 5 50 50 50 50 1 4 4 10];
				if (best_parameter(73) > 1e-10)
					lb(1) = best_parameter(73) - 1e-10;
					ub(1) = best_parameter(73) + 1e-10;
				end
				if (best_parameter(2) > 1e-10)
					lb(2) = best_parameter(2) - 1e-10;
					ub(2) = best_parameter(2) + 1e-10;
				end
				if (best_parameter(13) > 1e-10)
					lb(3) = best_parameter(13) - 1e-10;
					ub(3) = best_parameter(13) + 1e-10;
				end
				if (best_parameter(14) > 1e-10)
					lb(4) = best_parameter(14) - 1e-10;
					ub(4) = best_parameter(14) + 1e-10;
				end
				if (best_parameter(15) > 1e-10)
					lb(5) = best_parameter(15) - 1e-10;
					ub(5) = best_parameter(15) + 1e-10;
				end
				if (best_parameter(38) > 1e-10)
					lb(6) = best_parameter(38) - 1e-10;
					ub(6) = best_parameter(38) + 1e-10;
				end
				if (best_parameter(17) > 1e-10)
					lb(7) = best_parameter(17) - 1e-10;
					ub(7) = best_parameter(17) + 1e-10;
				end
				if (best_parameter(18) > 1e-10)
					lb(8) = best_parameter(18) - 1e-10;
					ub(8) = best_parameter(18) + 1e-10;
				end
				if (best_parameter(19) > 1e-10)
					lb(9) = best_parameter(19) - 1e-10;
					ub(9) = best_parameter(19) + 1e-10;
				end
				if (best_parameter(74) > 1e-10)
					lb(10) = best_parameter(74) - 1e-10;
					ub(10) = best_parameter(74) + 1e-10;
				end
				parameter_num = length(lb);
				x0 = zeros(1, parameter_num);
				x_opt = zeros(1, parameter_num);
				error_opt = 1e6;
				for i = 1:N
					for j = 1:parameter_num
						x0(j) = lb(j) + (ub(j) - lb(j))*rand;
					end
					[x_min, error_min] = lsqnonlin(@doi_10_1073pnas_0408507102_figure2A_1_error, x0, lb, ub, options);
					if (error_min < error_opt)
						error_opt = error_min;
						x_opt = x_min;
					end
				end
				best_parameter(73) = x_opt(1);
				best_parameter(2) = x_opt(2);
				best_parameter(13) = x_opt(3);
				best_parameter(14) = x_opt(4);
				best_parameter(15) = x_opt(5);
				best_parameter(38) = x_opt(6);
				best_parameter(17) = x_opt(7);
				best_parameter(18) = x_opt(8);
				best_parameter(19) = x_opt(9);
				best_parameter(74) = x_opt(10);
		end
	end
end

function seq_LOOCV(N1, N2)
	LOO_model_name_list = {'doi_10_1038_nature11516_figureS8A_1', 'doi_10_1038_nature11516_figureS8B_2', 'doi_10_1038_nbt_2401_figureS11_1', 'doi_10_1038_nbt_2401_figureS13_2', 'doi_10_1038_nature09565_figureS1C_2', 'doi_10_1186_1754_1611_3_4_figure3_1', 'doi_10_1186_1754_1611_3_4_figure3_2', 'doi_10_1186_1754_1611_3_4_figure3_6', 'doi_10_1186_1754_1611_3_4_figure3_7', 'doi__10_1073_pnas_1301301110_sd03_xls_1', 'doi__10_1073_pnas_1301301110_sd03_xls_2', 'doi_10_1093_nar_gku1388_figure2A_1', 'doi_10_1093_nar_gku1388_figure2A_2', 'doi_10_1093_nar_gku1388_figure2A_3', 'doi_10_1093_nar_gku1388_figure2A_4', 'doi_10_1093_nar_gku1388_figure2A_5', 'doi_10_1093_nar_gku1388_figure2A_6', 'doi_10_1093_nar_gku1388_figure2B_7', 'doi_10_1093_nar_gku1388_figure2B_8', 'doi_10_1093_nar_gku1388_figure2B_9', 'doi_10_1093_nar_gku1388_figure2B_10', 'doi_10_1093_nar_gku1388_figure2B_11', 'doi_10_1093_nar_gku1388_figure2B_12', 'doi_10_1093_nar_gku593_figure3D_2_2', 'doi_10_1093_nar_gku593_figure3D_3_3', 'doi_10_1093_nar_gku593_figure3D_4_4', 'doi_10_1093_nar_gku593_figure3D_5_5', 'doi_10_1093_nar_gku593_figure3D_6_6', 'doi_10_1093_nar_gku593_figureS6_7'};
	LOOCV_error = zeros(29,1);
	error_opt = 1e6;
	x_opt = zeros(1,74);
	for i = 1:N1
		model_order = LOOCV_randperm(41,4);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	arabinose = [0 0 0.002 0.01 0.04 0.16 1 4 25];
	mRFP1 = [0.001 0.001 0.001 0.003 0.045 0.161 0.184 0.185 0.185];
	output = simf_doi_10_1038_nature11516_figureS8A_1(x_opt, arabinose) - mRFP1;
	LOOCV_error(1) = norm(output);
	error_opt = 1e6;
	x_opt = zeros(1,74);
	for i = 1:N1
		model_order = LOOCV_randperm(41,5);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	IPTG = [0 0 0 0.001 0.004 0.025 0.1 0.5 2.5];
	mRFP1 = [0.005 0.005 0.005 0.006 0.013 0.09 0.173 0.194 0.196];
	output = simf_doi_10_1038_nature11516_figureS8B_2(x_opt, IPTG) - mRFP1;
	LOOCV_error(2) = norm(output);
	error_opt = 1e6;
	x_opt = zeros(1,74);
	for i = 1:N1
		model_order = LOOCV_randperm(41,6);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	IPTG = [0 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1];
	sfGFP = [0.025 0.025 0.035 0.055 0.115 0.2 0.3 0.4 0.65 1];
	output = simf_doi_10_1038_nbt_2401_figureS11_1(x_opt, IPTG) - sfGFP;
	LOOCV_error(3) = norm(output);
	error_opt = 1e6;
	x_opt = zeros(1,74);
	for i = 1:N1
		model_order = LOOCV_randperm(41,7);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	arabinose = [0.1 1 2 5 7 10 12.5 25 37.5 50 62.5];
	sfGFP = [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425];
	output = simf_doi_10_1038_nbt_2401_figureS13_2(x_opt, arabinose) - sfGFP;
	LOOCV_error(4) = norm(output);
	error_opt = 1e6;
	x_opt = zeros(1,74);
	for i = 1:N1
		model_order = LOOCV_randperm(41,9);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	aTc = [0 0.025 0.25 2.5 25 100 250];
	eYFP = [0.005 0.006 0.007 0.03 0.297 0.388 0.398];
	output = simf_doi_10_1038_nature09565_figureS1C_2(x_opt, aTc) - eYFP;
	LOOCV_error(5) = norm(output);
	error_opt = 1e6;
	x_opt = zeros(1,74);
	for i = 1:N1
		model_order = LOOCV_randperm(41,11);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	output = simf_doi_10_1186_1754_1611_3_4_figure3_1(x_opt) - 0.69;
	LOOCV_error(6) = norm(output);
	error_opt = 1e6;
	x_opt = zeros(1,74);
	for i = 1:N1
		model_order = LOOCV_randperm(41,12);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	output = simf_doi_10_1186_1754_1611_3_4_figure3_2(x_opt) - 0.018;
	LOOCV_error(7) = norm(output);
	error_opt = 1e6;
	x_opt = zeros(1,74);
	for i = 1:N1
		model_order = LOOCV_randperm(41,16);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	output = simf_doi_10_1186_1754_1611_3_4_figure3_6(x_opt) - 0.897;
	LOOCV_error(8) = norm(output);
	error_opt = 1e6;
	x_opt = zeros(1,74);
	for i = 1:N1
		model_order = LOOCV_randperm(41,17);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	output = simf_doi_10_1186_1754_1611_3_4_figure3_7(x_opt) - 1;
	LOOCV_error(9) = norm(output);
	error_opt = 1e6;
	x_opt = zeros(1,74);
	for i = 1:N1
		model_order = LOOCV_randperm(41,18);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	output = simf_doi__10_1073_pnas_1301301110_sd03_xls_1(x_opt) - 0.422;
	LOOCV_error(10) = norm(output);
	error_opt = 1e6;
	x_opt = zeros(1,74);
	for i = 1:N1
		model_order = LOOCV_randperm(41,19);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	output = simf_doi__10_1073_pnas_1301301110_sd03_xls_2(x_opt) - 1;
	LOOCV_error(11) = norm(output);
	error_opt = 1e6;
	x_opt = zeros(1,74);
	for i = 1:N1
		model_order = LOOCV_randperm(41,20);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	output = simf_doi_10_1093_nar_gku1388_figure2A_1(x_opt) - 1;
	LOOCV_error(12) = norm(output);
	error_opt = 1e6;
	x_opt = zeros(1,74);
	for i = 1:N1
		model_order = LOOCV_randperm(41,21);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	output = simf_doi_10_1093_nar_gku1388_figure2A_2(x_opt) - 0.215;
	LOOCV_error(13) = norm(output);
	error_opt = 1e6;
	x_opt = zeros(1,74);
	for i = 1:N1
		model_order = LOOCV_randperm(41,22);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	output = simf_doi_10_1093_nar_gku1388_figure2A_3(x_opt) - 0.095;
	LOOCV_error(14) = norm(output);
	error_opt = 1e6;
	x_opt = zeros(1,74);
	for i = 1:N1
		model_order = LOOCV_randperm(41,23);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	output = simf_doi_10_1093_nar_gku1388_figure2A_4(x_opt) - 0.046;
	LOOCV_error(15) = norm(output);
	error_opt = 1e6;
	x_opt = zeros(1,74);
	for i = 1:N1
		model_order = LOOCV_randperm(41,24);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	output = simf_doi_10_1093_nar_gku1388_figure2A_5(x_opt) - 0.018;
	LOOCV_error(16) = norm(output);
	error_opt = 1e6;
	x_opt = zeros(1,74);
	for i = 1:N1
		model_order = LOOCV_randperm(41,25);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	output = simf_doi_10_1093_nar_gku1388_figure2A_6(x_opt) - 0.003;
	LOOCV_error(17) = norm(output);
	error_opt = 1e6;
	x_opt = zeros(1,74);
	for i = 1:N1
		model_order = LOOCV_randperm(41,26);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	aTc = [0 2 4 10 25 50 100];
	GFPmut3b = [0 0.002 0.014 0.167 0.666 0.817 0.841];
	output = simf_doi_10_1093_nar_gku1388_figure2B_7(x_opt, aTc) - GFPmut3b;
	LOOCV_error(18) = norm(output);
	error_opt = 1e6;
	x_opt = zeros(1,74);
	for i = 1:N1
		model_order = LOOCV_randperm(41,27);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0 0 0 0.001 0.022 0.176 0.411 0.461];
	output = simf_doi_10_1093_nar_gku1388_figure2B_8(x_opt, aTc) - GFPmut3b;
	LOOCV_error(19) = norm(output);
	error_opt = 1e6;
	x_opt = zeros(1,74);
	for i = 1:N1
		model_order = LOOCV_randperm(41,28);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.005 0.118 0.277];
	output = simf_doi_10_1093_nar_gku1388_figure2B_9(x_opt, aTc) - GFPmut3b;
	LOOCV_error(20) = norm(output);
	error_opt = 1e6;
	x_opt = zeros(1,74);
	for i = 1:N1
		model_order = LOOCV_randperm(41,29);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.001 0.001 0.001 0.001 0.001 0.004 0.039 0.147];
	output = simf_doi_10_1093_nar_gku1388_figure2B_10(x_opt, aTc) - GFPmut3b;
	LOOCV_error(21) = norm(output);
	error_opt = 1e6;
	x_opt = zeros(1,74);
	for i = 1:N1
		model_order = LOOCV_randperm(41,30);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.003 0.004 0.032];
	output = simf_doi_10_1093_nar_gku1388_figure2B_11(x_opt, aTc) - GFPmut3b;
	LOOCV_error(22) = norm(output);
	error_opt = 1e6;
	x_opt = zeros(1,74);
	for i = 1:N1
		model_order = LOOCV_randperm(41,31);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.006 0.006 0.006 0.006 0.006 0.006 0.009 0.021];
	output = simf_doi_10_1093_nar_gku1388_figure2B_12(x_opt, aTc) - GFPmut3b;
	LOOCV_error(23) = norm(output);
	error_opt = 1e6;
	x_opt = zeros(1,74);
	for i = 1:N1
		model_order = LOOCV_randperm(41,33);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	output = simf_doi_10_1093_nar_gku593_figure3D_2_2(x_opt) - 0.018;
	LOOCV_error(24) = norm(output);
	error_opt = 1e6;
	x_opt = zeros(1,74);
	for i = 1:N1
		model_order = LOOCV_randperm(41,34);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	output = simf_doi_10_1093_nar_gku593_figure3D_3_3(x_opt) - 0.036;
	LOOCV_error(25) = norm(output);
	error_opt = 1e6;
	x_opt = zeros(1,74);
	for i = 1:N1
		model_order = LOOCV_randperm(41,35);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	output = simf_doi_10_1093_nar_gku593_figure3D_4_4(x_opt) - 0.05;
	LOOCV_error(26) = norm(output);
	error_opt = 1e6;
	x_opt = zeros(1,74);
	for i = 1:N1
		model_order = LOOCV_randperm(41,36);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	output = simf_doi_10_1093_nar_gku593_figure3D_5_5(x_opt) - 0.086;
	LOOCV_error(27) = norm(output);
	error_opt = 1e6;
	x_opt = zeros(1,74);
	for i = 1:N1
		model_order = LOOCV_randperm(41,37);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	output = simf_doi_10_1093_nar_gku593_figure3D_6_6(x_opt) - 0.214;
	LOOCV_error(28) = norm(output);
	error_opt = 1e6;
	x_opt = zeros(1,74);
	for i = 1:N1
		model_order = LOOCV_randperm(41,38);
		parameter = sequential_fitting_for_one_order(model_order, N2);
		error = norm(simultaneous_fitting_error(parameter));
		if (error < error_opt)
			error_opt = error;
			x_opt = parameter;
		end
	end
	arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2];
	GFPmut3b = [0 0.002 0.041 0.479 0.949 0.997 1 1 1];
	output = simf_doi_10_1093_nar_gku593_figureS6_7(x_opt, arabinose) - GFPmut3b;
	LOOCV_error(29) = norm(output);
	LOOCV_fileID = fopen('Plot_data/LOOCV_seq.dat','w');
	for i = 1:length(LOO_model_name_list)
		fprintf(LOOCV_fileID, '%s\t', LOO_model_name_list{i});
		fprintf(LOOCV_fileID, '%f\n', LOOCV_error(i));
	end
	fclose(LOOCV_fileID);
end

function x_opt = simultaneous_fitting(N)
	options = optimset('MaxIter', 1000, 'MaxFunEvals', 30000);
	parameter_name = {'ALPHA-RBS-pLAC'; 'ALPHA-RBS-placI'; 'K-IPTG'; 'K-pLAC'; 'alpha-pLAC'; 'alpha-placI'; 'beta-pLAC'; 'n-IPTG'; 'n-pLAC'; 'scale-lacZ-doi-10-1073-p'; 'ALPHA-RBS-pN25'; 'ALPHA-RBSII'; 'K-aTc'; 'K-pLtetO-1'; 'alpha-pLtetO-1'; 'alpha-pN25'; 'beta-pLtetO-1'; 'n-aTc'; 'n-pLtetO-1'; 'scale-luc-doi--10-1093-n'; 'ALPHA-RBS-pBAD'; 'ALPHA-RBS-pC'; 'K-arabinose'; 'K-pBAD'; 'alpha-pBAD'; 'alpha-pC'; 'beta-pBAD'; 'n-arabinose'; 'n-pBAD'; 'scale-GFPuv-doi-10-1128-'; 'ALPHA-BBa-B0030'; 'alpha-pLlacO-1'; 'scale-mCherry-doi-10-103'; 'ALPHA-RBS-psicA'; 'scale-mRFP1-doi-10-1038-'; 'K-pTAC'; 'alpha-pTAC'; 'alpha-placIq'; 'beta-pTAC'; 'n-pTAC'; 'ALPHA-RBS-A'; 'alpha-BBa-J23101'; 'scale-sfGFP-doi-10-1038-'; 'ALPHA-BBa-B0033'; 'ALPHA-BBa-B0034'; 'alpha-BBa-J23117'; 'scale-eYFP-doi-10-1038-n'; 'ALPHA-BBa-B0032'; 'ALPHA-RBS-pET-29b'; 'K-placUV5'; 'alpha-placUV5'; 'beta-placUV5'; 'n-placUV5'; 'scale-mRFP1-doi-10-1038-'; 'scale-GFPmut3b-doi-10-11'; 'alpha-BBa-J23116'; 'alpha-BBa-J23150'; 'alpha-BBa-J23151'; 'alpha-BBa-J23102'; 'scale-sfGFP-doi--10-1073'; 'scale-GFPmut3b-doi-10-10'; 'alpha-BBa-J23106'; 'alpha-BBa-J23105'; 'alpha-BBa-J23115'; 'alpha-BBa-J23114'; 'K-pTET*'; 'alpha-pTET*'; 'beta-pTET*'; 'n-pTET*'; 'alpha-BBa-J23109'; 'scale-GFPmut3b-doi-10-10'; 'scale-eYFP-http---dspace'; 'ALPHA-RBS-Unknown-Hoosha'; 'scale-eYFP-doi-10-1073pn'};
	lb = [0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01 0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01 0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01 0.03 0.002 0.01 0.03 0.01 0.001 0.002 0.002 0 1 0.03 0.002 0.01 0.03 0.03 0.002 0.01 0.03 0.03 0.001 0.002 0 1 0.01 0.01 0.002 0.002 0.002 0.002 0.01 0.01 0.002 0.002 0.002 0.002 0.001 0.002 0 1 0.002 0.01 0.01 0.03 0.01];
	ub = [5 5 50 50 50 50 1 4 4 10 5 5 50 50 50 50 1 4 4 10 5 5 50 50 50 50 1 4 4 10 5 50 10 5 10 50 50 50 1 4 5 50 10 5 5 50 10 5 5 50 50 1 4 10 10 50 50 50 50 10 10 50 50 50 50 50 50 1 4 50 10 10 5 10];
	parameter_num = length(lb);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@simultaneous_fitting_error, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	parameter_fileID = fopen('Plot_data/simultaneous_parameter_CI.dat','w');
	for i = 1:length(parameter_name)
		fprintf(parameter_fileID, '%s\t', parameter_name{i});
		fprintf(parameter_fileID, '%f\t', lb(i));
		fprintf(parameter_fileID, '%f\t', ub(i));
		fprintf(parameter_fileID, '%f\t', CI_lb(i));
		fprintf(parameter_fileID, '%f\t', CI_ub(i));
		fprintf(parameter_fileID, '%f\n', x_opt(i));
	end
	fclose(parameter_fileID);
end
function error = simultaneous_fitting_error(x)
	error = zeros(206,1);
	current_base_index = 0;
	%==============
	IPTG = [0 0.005 0.01 0.02 0.05 0.1 0.2 1];
	lacZ = [0.001 0.003 0.043 0.5 0.984 0.999 1 1];
	output = simf_doi_10_1073_pnas_0606717104_figure4_6_1(x, IPTG) - lacZ;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 1 2 3 4 5 6 7 8 9 10 12 13 17 18 20 25 35 60];
	luc = [0 0.001 0.001 0.001 0.002 0.004 0.007 0.011 0.014 0.032 0.036 0.054 0.089 0.643 0.893 1 1 1 1];
	output = simf_doi__10_1093_nar_25_6_1203_figure4a_1(x, aTc) - luc;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.003 0.005 0.05 0.5 1];
	GFPuv = [0.044 0.089 0.556 1 1];
	output = simf_doi_10_1128_AEM_00791_07_figure3a_1(x, arabinose) - GFPuv;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.003 0.008 0.025 0.075 0.2 0.7 2 6 20 50];
	mCherry = [0.015 0.015 0.015 0.02 0.06 0.48 0.84 0.93 0.96 0.98 0.98 1];
	output = simf_doi_10_1038_nature12148_figure16a_1(x, arabinose) - mCherry;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0 0.002 0.01 0.04 0.16 1 4 25];
	mRFP1 = [0.001 0.001 0.001 0.003 0.045 0.161 0.184 0.185 0.185];
	output = simf_doi_10_1038_nature11516_figureS8A_1(x, arabinose) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0 0 0.001 0.004 0.025 0.1 0.5 2.5];
	mRFP1 = [0.005 0.005 0.005 0.006 0.013 0.09 0.173 0.194 0.196];
	output = simf_doi_10_1038_nature11516_figureS8B_2(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1];
	sfGFP = [0.025 0.025 0.035 0.055 0.115 0.2 0.3 0.4 0.65 1];
	output = simf_doi_10_1038_nbt_2401_figureS11_1(x, IPTG) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.1 1 2 5 7 10 12.5 25 37.5 50 62.5];
	sfGFP = [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425];
	output = simf_doi_10_1038_nbt_2401_figureS13_2(x, arabinose) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.015 0.05 0.5 5 10];
	eYFP = [0.002 0.002 0.057 0.628 0.999 1 1];
	output = simf_doi_10_1038_nature09565_figure1C_1(x, arabinose) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 0.025 0.25 2.5 25 100 250];
	eYFP = [0.005 0.006 0.007 0.03 0.297 0.388 0.398];
	output = simf_doi_10_1038_nature09565_figureS1C_2(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.025 0.05 0.1 0.25 0.5];
	mRFP1 = [0.02 0.031 0.092 0.306 0.745 0.898 0.969 1 1];
	output = simf_doi_10_1038_nbt_2149_Supplement_page17_1(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_1(x) - 0.69;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_3(x) - 0.193;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_4(x) - 0.4;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_5(x) - 0.655;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_6(x) - 0.897;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_7(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_1(x) - 0.422;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_2(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_1(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_2(x) - 0.215;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_3(x) - 0.095;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_4(x) - 0.046;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_5(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_6(x) - 0.003;
	current_base_index = current_base_index + 1;
	%==============
	aTc = [0 2 4 10 25 50 100];
	GFPmut3b = [0 0.002 0.014 0.167 0.666 0.817 0.841];
	output = simf_doi_10_1093_nar_gku1388_figure2B_7(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0 0 0 0.001 0.022 0.176 0.411 0.461];
	output = simf_doi_10_1093_nar_gku1388_figure2B_8(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.005 0.118 0.277];
	output = simf_doi_10_1093_nar_gku1388_figure2B_9(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.001 0.001 0.001 0.001 0.001 0.004 0.039 0.147];
	output = simf_doi_10_1093_nar_gku1388_figure2B_10(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.003 0.004 0.032];
	output = simf_doi_10_1093_nar_gku1388_figure2B_11(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.006 0.006 0.006 0.006 0.006 0.006 0.009 0.021];
	output = simf_doi_10_1093_nar_gku1388_figure2B_12(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_1_1(x) - 0.002;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_2_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_3_3(x) - 0.036;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_4_4(x) - 0.05;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_5_5(x) - 0.086;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_6_6(x) - 0.214;
	current_base_index = current_base_index + 1;
	%==============
	arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2];
	GFPmut3b = [0 0.002 0.041 0.479 0.949 0.997 1 1 1];
	output = simf_doi_10_1093_nar_gku593_figureS6_7(x, arabinose) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0.001 0.01 0.03 0.05 0.075 0.1 0.15 0.2 0.3 0.5 1];
	eYFP = [0.015 0.014 0.023 0.062 0.098 0.187 0.239 0.374 0.572 0.869 1];
	output = simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1(x, IPTG) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [9.258 13.887 20.831 32.403 46.29 83.322 162.015 231.45 370.32 648.06 925.8 1388.7];
	eYFP = [0.003 0.003 0.003 0.003 0.005 0.01 0.05 0.167 0.333 0.667 0.833 1];
	output = simf_doi_10_1073pnas_0408507102_figure2A_1(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
end

function sim_LOOCV(N)
	options = optimset('MaxIter', 1000, 'MaxFunEvals', 30000);
	lb = [0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01 0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01 0.03 0.03 0.001 0.001 0.002 0.002 0 1 1 0.01 0.03 0.002 0.01 0.03 0.01 0.001 0.002 0.002 0 1 0.03 0.002 0.01 0.03 0.03 0.002 0.01 0.03 0.03 0.001 0.002 0 1 0.01 0.01 0.002 0.002 0.002 0.002 0.01 0.01 0.002 0.002 0.002 0.002 0.001 0.002 0 1 0.002 0.01 0.01 0.03 0.01];
	ub = [5 5 50 50 50 50 1 4 4 10 5 5 50 50 50 50 1 4 4 10 5 5 50 50 50 50 1 4 4 10 5 50 10 5 10 50 50 50 1 4 5 50 10 5 5 50 10 5 5 50 50 1 4 10 10 50 50 50 50 10 10 50 50 50 50 50 50 1 4 50 10 10 5 10];
	parameter_num = length(lb);
	LOO_model_name_list = {'doi_10_1038_nature11516_figureS8A_1', 'doi_10_1038_nature11516_figureS8B_2', 'doi_10_1038_nbt_2401_figureS11_1', 'doi_10_1038_nbt_2401_figureS13_2', 'doi_10_1038_nature09565_figureS1C_2', 'doi_10_1186_1754_1611_3_4_figure3_1', 'doi_10_1186_1754_1611_3_4_figure3_2', 'doi_10_1186_1754_1611_3_4_figure3_6', 'doi_10_1186_1754_1611_3_4_figure3_7', 'doi__10_1073_pnas_1301301110_sd03_xls_1', 'doi__10_1073_pnas_1301301110_sd03_xls_2', 'doi_10_1093_nar_gku1388_figure2A_1', 'doi_10_1093_nar_gku1388_figure2A_2', 'doi_10_1093_nar_gku1388_figure2A_3', 'doi_10_1093_nar_gku1388_figure2A_4', 'doi_10_1093_nar_gku1388_figure2A_5', 'doi_10_1093_nar_gku1388_figure2A_6', 'doi_10_1093_nar_gku1388_figure2B_7', 'doi_10_1093_nar_gku1388_figure2B_8', 'doi_10_1093_nar_gku1388_figure2B_9', 'doi_10_1093_nar_gku1388_figure2B_10', 'doi_10_1093_nar_gku1388_figure2B_11', 'doi_10_1093_nar_gku1388_figure2B_12', 'doi_10_1093_nar_gku593_figure3D_2_2', 'doi_10_1093_nar_gku593_figure3D_3_3', 'doi_10_1093_nar_gku593_figure3D_4_4', 'doi_10_1093_nar_gku593_figure3D_5_5', 'doi_10_1093_nar_gku593_figure3D_6_6', 'doi_10_1093_nar_gku593_figureS6_7'};
	LOOCV_error = zeros(29,1);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@LOOCV_sim_error_doi_10_1038_nature11516_figureS8A_1, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	arabinose = [0 0 0.002 0.01 0.04 0.16 1 4 25];
	mRFP1 = [0.001 0.001 0.001 0.003 0.045 0.161 0.184 0.185 0.185];
	output = simf_doi_10_1038_nature11516_figureS8A_1(x_opt, arabinose) - mRFP1;
	LOOCV_error(1) = norm(output);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@LOOCV_sim_error_doi_10_1038_nature11516_figureS8B_2, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	IPTG = [0 0 0 0.001 0.004 0.025 0.1 0.5 2.5];
	mRFP1 = [0.005 0.005 0.005 0.006 0.013 0.09 0.173 0.194 0.196];
	output = simf_doi_10_1038_nature11516_figureS8B_2(x_opt, IPTG) - mRFP1;
	LOOCV_error(2) = norm(output);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@LOOCV_sim_error_doi_10_1038_nbt_2401_figureS11_1, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	IPTG = [0 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1];
	sfGFP = [0.025 0.025 0.035 0.055 0.115 0.2 0.3 0.4 0.65 1];
	output = simf_doi_10_1038_nbt_2401_figureS11_1(x_opt, IPTG) - sfGFP;
	LOOCV_error(3) = norm(output);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@LOOCV_sim_error_doi_10_1038_nbt_2401_figureS13_2, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	arabinose = [0.1 1 2 5 7 10 12.5 25 37.5 50 62.5];
	sfGFP = [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425];
	output = simf_doi_10_1038_nbt_2401_figureS13_2(x_opt, arabinose) - sfGFP;
	LOOCV_error(4) = norm(output);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@LOOCV_sim_error_doi_10_1038_nature09565_figureS1C_2, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	aTc = [0 0.025 0.25 2.5 25 100 250];
	eYFP = [0.005 0.006 0.007 0.03 0.297 0.388 0.398];
	output = simf_doi_10_1038_nature09565_figureS1C_2(x_opt, aTc) - eYFP;
	LOOCV_error(5) = norm(output);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@LOOCV_sim_error_doi_10_1186_1754_1611_3_4_figure3_1, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	output = simf_doi_10_1186_1754_1611_3_4_figure3_1(x_opt) - 0.69;
	LOOCV_error(6) = norm(output);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@LOOCV_sim_error_doi_10_1186_1754_1611_3_4_figure3_2, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	output = simf_doi_10_1186_1754_1611_3_4_figure3_2(x_opt) - 0.018;
	LOOCV_error(7) = norm(output);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@LOOCV_sim_error_doi_10_1186_1754_1611_3_4_figure3_6, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	output = simf_doi_10_1186_1754_1611_3_4_figure3_6(x_opt) - 0.897;
	LOOCV_error(8) = norm(output);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@LOOCV_sim_error_doi_10_1186_1754_1611_3_4_figure3_7, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	output = simf_doi_10_1186_1754_1611_3_4_figure3_7(x_opt) - 1;
	LOOCV_error(9) = norm(output);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@LOOCV_sim_error_doi__10_1073_pnas_1301301110_sd03_xls_1, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	output = simf_doi__10_1073_pnas_1301301110_sd03_xls_1(x_opt) - 0.422;
	LOOCV_error(10) = norm(output);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@LOOCV_sim_error_doi__10_1073_pnas_1301301110_sd03_xls_2, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	output = simf_doi__10_1073_pnas_1301301110_sd03_xls_2(x_opt) - 1;
	LOOCV_error(11) = norm(output);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2A_1, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	output = simf_doi_10_1093_nar_gku1388_figure2A_1(x_opt) - 1;
	LOOCV_error(12) = norm(output);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2A_2, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	output = simf_doi_10_1093_nar_gku1388_figure2A_2(x_opt) - 0.215;
	LOOCV_error(13) = norm(output);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2A_3, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	output = simf_doi_10_1093_nar_gku1388_figure2A_3(x_opt) - 0.095;
	LOOCV_error(14) = norm(output);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2A_4, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	output = simf_doi_10_1093_nar_gku1388_figure2A_4(x_opt) - 0.046;
	LOOCV_error(15) = norm(output);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2A_5, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	output = simf_doi_10_1093_nar_gku1388_figure2A_5(x_opt) - 0.018;
	LOOCV_error(16) = norm(output);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2A_6, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	output = simf_doi_10_1093_nar_gku1388_figure2A_6(x_opt) - 0.003;
	LOOCV_error(17) = norm(output);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2B_7, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	aTc = [0 2 4 10 25 50 100];
	GFPmut3b = [0 0.002 0.014 0.167 0.666 0.817 0.841];
	output = simf_doi_10_1093_nar_gku1388_figure2B_7(x_opt, aTc) - GFPmut3b;
	LOOCV_error(18) = norm(output);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2B_8, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0 0 0 0.001 0.022 0.176 0.411 0.461];
	output = simf_doi_10_1093_nar_gku1388_figure2B_8(x_opt, aTc) - GFPmut3b;
	LOOCV_error(19) = norm(output);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2B_9, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.005 0.118 0.277];
	output = simf_doi_10_1093_nar_gku1388_figure2B_9(x_opt, aTc) - GFPmut3b;
	LOOCV_error(20) = norm(output);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2B_10, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.001 0.001 0.001 0.001 0.001 0.004 0.039 0.147];
	output = simf_doi_10_1093_nar_gku1388_figure2B_10(x_opt, aTc) - GFPmut3b;
	LOOCV_error(21) = norm(output);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2B_11, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.003 0.004 0.032];
	output = simf_doi_10_1093_nar_gku1388_figure2B_11(x_opt, aTc) - GFPmut3b;
	LOOCV_error(22) = norm(output);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2B_12, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.006 0.006 0.006 0.006 0.006 0.006 0.009 0.021];
	output = simf_doi_10_1093_nar_gku1388_figure2B_12(x_opt, aTc) - GFPmut3b;
	LOOCV_error(23) = norm(output);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@LOOCV_sim_error_doi_10_1093_nar_gku593_figure3D_2_2, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	output = simf_doi_10_1093_nar_gku593_figure3D_2_2(x_opt) - 0.018;
	LOOCV_error(24) = norm(output);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@LOOCV_sim_error_doi_10_1093_nar_gku593_figure3D_3_3, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	output = simf_doi_10_1093_nar_gku593_figure3D_3_3(x_opt) - 0.036;
	LOOCV_error(25) = norm(output);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@LOOCV_sim_error_doi_10_1093_nar_gku593_figure3D_4_4, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	output = simf_doi_10_1093_nar_gku593_figure3D_4_4(x_opt) - 0.05;
	LOOCV_error(26) = norm(output);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@LOOCV_sim_error_doi_10_1093_nar_gku593_figure3D_5_5, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	output = simf_doi_10_1093_nar_gku593_figure3D_5_5(x_opt) - 0.086;
	LOOCV_error(27) = norm(output);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@LOOCV_sim_error_doi_10_1093_nar_gku593_figure3D_6_6, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	output = simf_doi_10_1093_nar_gku593_figure3D_6_6(x_opt) - 0.214;
	LOOCV_error(28) = norm(output);
	solution_ensemble = zeros(N, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:N
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(@LOOCV_sim_error_doi_10_1093_nar_gku593_figureS6_7, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	[CI_lb, CI_ub, x_opt] = extract_CI(solution_ensemble);
	arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2];
	GFPmut3b = [0 0.002 0.041 0.479 0.949 0.997 1 1 1];
	output = simf_doi_10_1093_nar_gku593_figureS6_7(x_opt, arabinose) - GFPmut3b;
	LOOCV_error(29) = norm(output);
	LOOCV_fileID = fopen('Plot_data/LOOCV_sim.dat','w');
	for i = 1:length(LOO_model_name_list)
		fprintf(LOOCV_fileID, '%s\t', LOO_model_name_list{i});
		fprintf(LOOCV_fileID, '%f\n', LOOCV_error(i));
	end
	fclose(LOOCV_fileID);
end
function error = LOOCV_sim_error_doi_10_1038_nature11516_figureS8A_1(x)
	error = zeros(197,1);
	current_base_index = 0;
	%==============
	IPTG = [0 0.005 0.01 0.02 0.05 0.1 0.2 1];
	lacZ = [0.001 0.003 0.043 0.5 0.984 0.999 1 1];
	output = simf_doi_10_1073_pnas_0606717104_figure4_6_1(x, IPTG) - lacZ;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 1 2 3 4 5 6 7 8 9 10 12 13 17 18 20 25 35 60];
	luc = [0 0.001 0.001 0.001 0.002 0.004 0.007 0.011 0.014 0.032 0.036 0.054 0.089 0.643 0.893 1 1 1 1];
	output = simf_doi__10_1093_nar_25_6_1203_figure4a_1(x, aTc) - luc;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.003 0.005 0.05 0.5 1];
	GFPuv = [0.044 0.089 0.556 1 1];
	output = simf_doi_10_1128_AEM_00791_07_figure3a_1(x, arabinose) - GFPuv;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.003 0.008 0.025 0.075 0.2 0.7 2 6 20 50];
	mCherry = [0.015 0.015 0.015 0.02 0.06 0.48 0.84 0.93 0.96 0.98 0.98 1];
	output = simf_doi_10_1038_nature12148_figure16a_1(x, arabinose) - mCherry;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0 0 0.001 0.004 0.025 0.1 0.5 2.5];
	mRFP1 = [0.005 0.005 0.005 0.006 0.013 0.09 0.173 0.194 0.196];
	output = simf_doi_10_1038_nature11516_figureS8B_2(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1];
	sfGFP = [0.025 0.025 0.035 0.055 0.115 0.2 0.3 0.4 0.65 1];
	output = simf_doi_10_1038_nbt_2401_figureS11_1(x, IPTG) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.1 1 2 5 7 10 12.5 25 37.5 50 62.5];
	sfGFP = [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425];
	output = simf_doi_10_1038_nbt_2401_figureS13_2(x, arabinose) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.015 0.05 0.5 5 10];
	eYFP = [0.002 0.002 0.057 0.628 0.999 1 1];
	output = simf_doi_10_1038_nature09565_figure1C_1(x, arabinose) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 0.025 0.25 2.5 25 100 250];
	eYFP = [0.005 0.006 0.007 0.03 0.297 0.388 0.398];
	output = simf_doi_10_1038_nature09565_figureS1C_2(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.025 0.05 0.1 0.25 0.5];
	mRFP1 = [0.02 0.031 0.092 0.306 0.745 0.898 0.969 1 1];
	output = simf_doi_10_1038_nbt_2149_Supplement_page17_1(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_1(x) - 0.69;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_3(x) - 0.193;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_4(x) - 0.4;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_5(x) - 0.655;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_6(x) - 0.897;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_7(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_1(x) - 0.422;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_2(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_1(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_2(x) - 0.215;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_3(x) - 0.095;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_4(x) - 0.046;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_5(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_6(x) - 0.003;
	current_base_index = current_base_index + 1;
	%==============
	aTc = [0 2 4 10 25 50 100];
	GFPmut3b = [0 0.002 0.014 0.167 0.666 0.817 0.841];
	output = simf_doi_10_1093_nar_gku1388_figure2B_7(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0 0 0 0.001 0.022 0.176 0.411 0.461];
	output = simf_doi_10_1093_nar_gku1388_figure2B_8(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.005 0.118 0.277];
	output = simf_doi_10_1093_nar_gku1388_figure2B_9(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.001 0.001 0.001 0.001 0.001 0.004 0.039 0.147];
	output = simf_doi_10_1093_nar_gku1388_figure2B_10(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.003 0.004 0.032];
	output = simf_doi_10_1093_nar_gku1388_figure2B_11(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.006 0.006 0.006 0.006 0.006 0.006 0.009 0.021];
	output = simf_doi_10_1093_nar_gku1388_figure2B_12(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_1_1(x) - 0.002;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_2_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_3_3(x) - 0.036;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_4_4(x) - 0.05;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_5_5(x) - 0.086;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_6_6(x) - 0.214;
	current_base_index = current_base_index + 1;
	%==============
	arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2];
	GFPmut3b = [0 0.002 0.041 0.479 0.949 0.997 1 1 1];
	output = simf_doi_10_1093_nar_gku593_figureS6_7(x, arabinose) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0.001 0.01 0.03 0.05 0.075 0.1 0.15 0.2 0.3 0.5 1];
	eYFP = [0.015 0.014 0.023 0.062 0.098 0.187 0.239 0.374 0.572 0.869 1];
	output = simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1(x, IPTG) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [9.258 13.887 20.831 32.403 46.29 83.322 162.015 231.45 370.32 648.06 925.8 1388.7];
	eYFP = [0.003 0.003 0.003 0.003 0.005 0.01 0.05 0.167 0.333 0.667 0.833 1];
	output = simf_doi_10_1073pnas_0408507102_figure2A_1(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
end
function error = LOOCV_sim_error_doi_10_1038_nature11516_figureS8B_2(x)
	error = zeros(197,1);
	current_base_index = 0;
	%==============
	IPTG = [0 0.005 0.01 0.02 0.05 0.1 0.2 1];
	lacZ = [0.001 0.003 0.043 0.5 0.984 0.999 1 1];
	output = simf_doi_10_1073_pnas_0606717104_figure4_6_1(x, IPTG) - lacZ;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 1 2 3 4 5 6 7 8 9 10 12 13 17 18 20 25 35 60];
	luc = [0 0.001 0.001 0.001 0.002 0.004 0.007 0.011 0.014 0.032 0.036 0.054 0.089 0.643 0.893 1 1 1 1];
	output = simf_doi__10_1093_nar_25_6_1203_figure4a_1(x, aTc) - luc;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.003 0.005 0.05 0.5 1];
	GFPuv = [0.044 0.089 0.556 1 1];
	output = simf_doi_10_1128_AEM_00791_07_figure3a_1(x, arabinose) - GFPuv;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.003 0.008 0.025 0.075 0.2 0.7 2 6 20 50];
	mCherry = [0.015 0.015 0.015 0.02 0.06 0.48 0.84 0.93 0.96 0.98 0.98 1];
	output = simf_doi_10_1038_nature12148_figure16a_1(x, arabinose) - mCherry;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0 0.002 0.01 0.04 0.16 1 4 25];
	mRFP1 = [0.001 0.001 0.001 0.003 0.045 0.161 0.184 0.185 0.185];
	output = simf_doi_10_1038_nature11516_figureS8A_1(x, arabinose) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1];
	sfGFP = [0.025 0.025 0.035 0.055 0.115 0.2 0.3 0.4 0.65 1];
	output = simf_doi_10_1038_nbt_2401_figureS11_1(x, IPTG) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.1 1 2 5 7 10 12.5 25 37.5 50 62.5];
	sfGFP = [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425];
	output = simf_doi_10_1038_nbt_2401_figureS13_2(x, arabinose) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.015 0.05 0.5 5 10];
	eYFP = [0.002 0.002 0.057 0.628 0.999 1 1];
	output = simf_doi_10_1038_nature09565_figure1C_1(x, arabinose) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 0.025 0.25 2.5 25 100 250];
	eYFP = [0.005 0.006 0.007 0.03 0.297 0.388 0.398];
	output = simf_doi_10_1038_nature09565_figureS1C_2(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.025 0.05 0.1 0.25 0.5];
	mRFP1 = [0.02 0.031 0.092 0.306 0.745 0.898 0.969 1 1];
	output = simf_doi_10_1038_nbt_2149_Supplement_page17_1(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_1(x) - 0.69;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_3(x) - 0.193;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_4(x) - 0.4;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_5(x) - 0.655;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_6(x) - 0.897;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_7(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_1(x) - 0.422;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_2(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_1(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_2(x) - 0.215;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_3(x) - 0.095;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_4(x) - 0.046;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_5(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_6(x) - 0.003;
	current_base_index = current_base_index + 1;
	%==============
	aTc = [0 2 4 10 25 50 100];
	GFPmut3b = [0 0.002 0.014 0.167 0.666 0.817 0.841];
	output = simf_doi_10_1093_nar_gku1388_figure2B_7(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0 0 0 0.001 0.022 0.176 0.411 0.461];
	output = simf_doi_10_1093_nar_gku1388_figure2B_8(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.005 0.118 0.277];
	output = simf_doi_10_1093_nar_gku1388_figure2B_9(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.001 0.001 0.001 0.001 0.001 0.004 0.039 0.147];
	output = simf_doi_10_1093_nar_gku1388_figure2B_10(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.003 0.004 0.032];
	output = simf_doi_10_1093_nar_gku1388_figure2B_11(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.006 0.006 0.006 0.006 0.006 0.006 0.009 0.021];
	output = simf_doi_10_1093_nar_gku1388_figure2B_12(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_1_1(x) - 0.002;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_2_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_3_3(x) - 0.036;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_4_4(x) - 0.05;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_5_5(x) - 0.086;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_6_6(x) - 0.214;
	current_base_index = current_base_index + 1;
	%==============
	arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2];
	GFPmut3b = [0 0.002 0.041 0.479 0.949 0.997 1 1 1];
	output = simf_doi_10_1093_nar_gku593_figureS6_7(x, arabinose) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0.001 0.01 0.03 0.05 0.075 0.1 0.15 0.2 0.3 0.5 1];
	eYFP = [0.015 0.014 0.023 0.062 0.098 0.187 0.239 0.374 0.572 0.869 1];
	output = simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1(x, IPTG) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [9.258 13.887 20.831 32.403 46.29 83.322 162.015 231.45 370.32 648.06 925.8 1388.7];
	eYFP = [0.003 0.003 0.003 0.003 0.005 0.01 0.05 0.167 0.333 0.667 0.833 1];
	output = simf_doi_10_1073pnas_0408507102_figure2A_1(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
end
function error = LOOCV_sim_error_doi_10_1038_nbt_2401_figureS11_1(x)
	error = zeros(196,1);
	current_base_index = 0;
	%==============
	IPTG = [0 0.005 0.01 0.02 0.05 0.1 0.2 1];
	lacZ = [0.001 0.003 0.043 0.5 0.984 0.999 1 1];
	output = simf_doi_10_1073_pnas_0606717104_figure4_6_1(x, IPTG) - lacZ;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 1 2 3 4 5 6 7 8 9 10 12 13 17 18 20 25 35 60];
	luc = [0 0.001 0.001 0.001 0.002 0.004 0.007 0.011 0.014 0.032 0.036 0.054 0.089 0.643 0.893 1 1 1 1];
	output = simf_doi__10_1093_nar_25_6_1203_figure4a_1(x, aTc) - luc;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.003 0.005 0.05 0.5 1];
	GFPuv = [0.044 0.089 0.556 1 1];
	output = simf_doi_10_1128_AEM_00791_07_figure3a_1(x, arabinose) - GFPuv;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.003 0.008 0.025 0.075 0.2 0.7 2 6 20 50];
	mCherry = [0.015 0.015 0.015 0.02 0.06 0.48 0.84 0.93 0.96 0.98 0.98 1];
	output = simf_doi_10_1038_nature12148_figure16a_1(x, arabinose) - mCherry;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0 0.002 0.01 0.04 0.16 1 4 25];
	mRFP1 = [0.001 0.001 0.001 0.003 0.045 0.161 0.184 0.185 0.185];
	output = simf_doi_10_1038_nature11516_figureS8A_1(x, arabinose) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0 0 0.001 0.004 0.025 0.1 0.5 2.5];
	mRFP1 = [0.005 0.005 0.005 0.006 0.013 0.09 0.173 0.194 0.196];
	output = simf_doi_10_1038_nature11516_figureS8B_2(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.1 1 2 5 7 10 12.5 25 37.5 50 62.5];
	sfGFP = [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425];
	output = simf_doi_10_1038_nbt_2401_figureS13_2(x, arabinose) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.015 0.05 0.5 5 10];
	eYFP = [0.002 0.002 0.057 0.628 0.999 1 1];
	output = simf_doi_10_1038_nature09565_figure1C_1(x, arabinose) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 0.025 0.25 2.5 25 100 250];
	eYFP = [0.005 0.006 0.007 0.03 0.297 0.388 0.398];
	output = simf_doi_10_1038_nature09565_figureS1C_2(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.025 0.05 0.1 0.25 0.5];
	mRFP1 = [0.02 0.031 0.092 0.306 0.745 0.898 0.969 1 1];
	output = simf_doi_10_1038_nbt_2149_Supplement_page17_1(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_1(x) - 0.69;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_3(x) - 0.193;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_4(x) - 0.4;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_5(x) - 0.655;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_6(x) - 0.897;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_7(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_1(x) - 0.422;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_2(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_1(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_2(x) - 0.215;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_3(x) - 0.095;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_4(x) - 0.046;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_5(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_6(x) - 0.003;
	current_base_index = current_base_index + 1;
	%==============
	aTc = [0 2 4 10 25 50 100];
	GFPmut3b = [0 0.002 0.014 0.167 0.666 0.817 0.841];
	output = simf_doi_10_1093_nar_gku1388_figure2B_7(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0 0 0 0.001 0.022 0.176 0.411 0.461];
	output = simf_doi_10_1093_nar_gku1388_figure2B_8(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.005 0.118 0.277];
	output = simf_doi_10_1093_nar_gku1388_figure2B_9(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.001 0.001 0.001 0.001 0.001 0.004 0.039 0.147];
	output = simf_doi_10_1093_nar_gku1388_figure2B_10(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.003 0.004 0.032];
	output = simf_doi_10_1093_nar_gku1388_figure2B_11(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.006 0.006 0.006 0.006 0.006 0.006 0.009 0.021];
	output = simf_doi_10_1093_nar_gku1388_figure2B_12(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_1_1(x) - 0.002;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_2_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_3_3(x) - 0.036;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_4_4(x) - 0.05;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_5_5(x) - 0.086;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_6_6(x) - 0.214;
	current_base_index = current_base_index + 1;
	%==============
	arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2];
	GFPmut3b = [0 0.002 0.041 0.479 0.949 0.997 1 1 1];
	output = simf_doi_10_1093_nar_gku593_figureS6_7(x, arabinose) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0.001 0.01 0.03 0.05 0.075 0.1 0.15 0.2 0.3 0.5 1];
	eYFP = [0.015 0.014 0.023 0.062 0.098 0.187 0.239 0.374 0.572 0.869 1];
	output = simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1(x, IPTG) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [9.258 13.887 20.831 32.403 46.29 83.322 162.015 231.45 370.32 648.06 925.8 1388.7];
	eYFP = [0.003 0.003 0.003 0.003 0.005 0.01 0.05 0.167 0.333 0.667 0.833 1];
	output = simf_doi_10_1073pnas_0408507102_figure2A_1(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
end
function error = LOOCV_sim_error_doi_10_1038_nbt_2401_figureS13_2(x)
	error = zeros(195,1);
	current_base_index = 0;
	%==============
	IPTG = [0 0.005 0.01 0.02 0.05 0.1 0.2 1];
	lacZ = [0.001 0.003 0.043 0.5 0.984 0.999 1 1];
	output = simf_doi_10_1073_pnas_0606717104_figure4_6_1(x, IPTG) - lacZ;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 1 2 3 4 5 6 7 8 9 10 12 13 17 18 20 25 35 60];
	luc = [0 0.001 0.001 0.001 0.002 0.004 0.007 0.011 0.014 0.032 0.036 0.054 0.089 0.643 0.893 1 1 1 1];
	output = simf_doi__10_1093_nar_25_6_1203_figure4a_1(x, aTc) - luc;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.003 0.005 0.05 0.5 1];
	GFPuv = [0.044 0.089 0.556 1 1];
	output = simf_doi_10_1128_AEM_00791_07_figure3a_1(x, arabinose) - GFPuv;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.003 0.008 0.025 0.075 0.2 0.7 2 6 20 50];
	mCherry = [0.015 0.015 0.015 0.02 0.06 0.48 0.84 0.93 0.96 0.98 0.98 1];
	output = simf_doi_10_1038_nature12148_figure16a_1(x, arabinose) - mCherry;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0 0.002 0.01 0.04 0.16 1 4 25];
	mRFP1 = [0.001 0.001 0.001 0.003 0.045 0.161 0.184 0.185 0.185];
	output = simf_doi_10_1038_nature11516_figureS8A_1(x, arabinose) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0 0 0.001 0.004 0.025 0.1 0.5 2.5];
	mRFP1 = [0.005 0.005 0.005 0.006 0.013 0.09 0.173 0.194 0.196];
	output = simf_doi_10_1038_nature11516_figureS8B_2(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1];
	sfGFP = [0.025 0.025 0.035 0.055 0.115 0.2 0.3 0.4 0.65 1];
	output = simf_doi_10_1038_nbt_2401_figureS11_1(x, IPTG) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.015 0.05 0.5 5 10];
	eYFP = [0.002 0.002 0.057 0.628 0.999 1 1];
	output = simf_doi_10_1038_nature09565_figure1C_1(x, arabinose) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 0.025 0.25 2.5 25 100 250];
	eYFP = [0.005 0.006 0.007 0.03 0.297 0.388 0.398];
	output = simf_doi_10_1038_nature09565_figureS1C_2(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.025 0.05 0.1 0.25 0.5];
	mRFP1 = [0.02 0.031 0.092 0.306 0.745 0.898 0.969 1 1];
	output = simf_doi_10_1038_nbt_2149_Supplement_page17_1(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_1(x) - 0.69;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_3(x) - 0.193;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_4(x) - 0.4;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_5(x) - 0.655;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_6(x) - 0.897;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_7(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_1(x) - 0.422;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_2(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_1(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_2(x) - 0.215;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_3(x) - 0.095;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_4(x) - 0.046;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_5(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_6(x) - 0.003;
	current_base_index = current_base_index + 1;
	%==============
	aTc = [0 2 4 10 25 50 100];
	GFPmut3b = [0 0.002 0.014 0.167 0.666 0.817 0.841];
	output = simf_doi_10_1093_nar_gku1388_figure2B_7(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0 0 0 0.001 0.022 0.176 0.411 0.461];
	output = simf_doi_10_1093_nar_gku1388_figure2B_8(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.005 0.118 0.277];
	output = simf_doi_10_1093_nar_gku1388_figure2B_9(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.001 0.001 0.001 0.001 0.001 0.004 0.039 0.147];
	output = simf_doi_10_1093_nar_gku1388_figure2B_10(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.003 0.004 0.032];
	output = simf_doi_10_1093_nar_gku1388_figure2B_11(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.006 0.006 0.006 0.006 0.006 0.006 0.009 0.021];
	output = simf_doi_10_1093_nar_gku1388_figure2B_12(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_1_1(x) - 0.002;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_2_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_3_3(x) - 0.036;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_4_4(x) - 0.05;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_5_5(x) - 0.086;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_6_6(x) - 0.214;
	current_base_index = current_base_index + 1;
	%==============
	arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2];
	GFPmut3b = [0 0.002 0.041 0.479 0.949 0.997 1 1 1];
	output = simf_doi_10_1093_nar_gku593_figureS6_7(x, arabinose) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0.001 0.01 0.03 0.05 0.075 0.1 0.15 0.2 0.3 0.5 1];
	eYFP = [0.015 0.014 0.023 0.062 0.098 0.187 0.239 0.374 0.572 0.869 1];
	output = simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1(x, IPTG) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [9.258 13.887 20.831 32.403 46.29 83.322 162.015 231.45 370.32 648.06 925.8 1388.7];
	eYFP = [0.003 0.003 0.003 0.003 0.005 0.01 0.05 0.167 0.333 0.667 0.833 1];
	output = simf_doi_10_1073pnas_0408507102_figure2A_1(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
end
function error = LOOCV_sim_error_doi_10_1038_nature09565_figureS1C_2(x)
	error = zeros(199,1);
	current_base_index = 0;
	%==============
	IPTG = [0 0.005 0.01 0.02 0.05 0.1 0.2 1];
	lacZ = [0.001 0.003 0.043 0.5 0.984 0.999 1 1];
	output = simf_doi_10_1073_pnas_0606717104_figure4_6_1(x, IPTG) - lacZ;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 1 2 3 4 5 6 7 8 9 10 12 13 17 18 20 25 35 60];
	luc = [0 0.001 0.001 0.001 0.002 0.004 0.007 0.011 0.014 0.032 0.036 0.054 0.089 0.643 0.893 1 1 1 1];
	output = simf_doi__10_1093_nar_25_6_1203_figure4a_1(x, aTc) - luc;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.003 0.005 0.05 0.5 1];
	GFPuv = [0.044 0.089 0.556 1 1];
	output = simf_doi_10_1128_AEM_00791_07_figure3a_1(x, arabinose) - GFPuv;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.003 0.008 0.025 0.075 0.2 0.7 2 6 20 50];
	mCherry = [0.015 0.015 0.015 0.02 0.06 0.48 0.84 0.93 0.96 0.98 0.98 1];
	output = simf_doi_10_1038_nature12148_figure16a_1(x, arabinose) - mCherry;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0 0.002 0.01 0.04 0.16 1 4 25];
	mRFP1 = [0.001 0.001 0.001 0.003 0.045 0.161 0.184 0.185 0.185];
	output = simf_doi_10_1038_nature11516_figureS8A_1(x, arabinose) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0 0 0.001 0.004 0.025 0.1 0.5 2.5];
	mRFP1 = [0.005 0.005 0.005 0.006 0.013 0.09 0.173 0.194 0.196];
	output = simf_doi_10_1038_nature11516_figureS8B_2(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1];
	sfGFP = [0.025 0.025 0.035 0.055 0.115 0.2 0.3 0.4 0.65 1];
	output = simf_doi_10_1038_nbt_2401_figureS11_1(x, IPTG) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.1 1 2 5 7 10 12.5 25 37.5 50 62.5];
	sfGFP = [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425];
	output = simf_doi_10_1038_nbt_2401_figureS13_2(x, arabinose) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.015 0.05 0.5 5 10];
	eYFP = [0.002 0.002 0.057 0.628 0.999 1 1];
	output = simf_doi_10_1038_nature09565_figure1C_1(x, arabinose) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.025 0.05 0.1 0.25 0.5];
	mRFP1 = [0.02 0.031 0.092 0.306 0.745 0.898 0.969 1 1];
	output = simf_doi_10_1038_nbt_2149_Supplement_page17_1(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_1(x) - 0.69;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_3(x) - 0.193;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_4(x) - 0.4;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_5(x) - 0.655;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_6(x) - 0.897;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_7(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_1(x) - 0.422;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_2(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_1(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_2(x) - 0.215;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_3(x) - 0.095;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_4(x) - 0.046;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_5(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_6(x) - 0.003;
	current_base_index = current_base_index + 1;
	%==============
	aTc = [0 2 4 10 25 50 100];
	GFPmut3b = [0 0.002 0.014 0.167 0.666 0.817 0.841];
	output = simf_doi_10_1093_nar_gku1388_figure2B_7(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0 0 0 0.001 0.022 0.176 0.411 0.461];
	output = simf_doi_10_1093_nar_gku1388_figure2B_8(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.005 0.118 0.277];
	output = simf_doi_10_1093_nar_gku1388_figure2B_9(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.001 0.001 0.001 0.001 0.001 0.004 0.039 0.147];
	output = simf_doi_10_1093_nar_gku1388_figure2B_10(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.003 0.004 0.032];
	output = simf_doi_10_1093_nar_gku1388_figure2B_11(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.006 0.006 0.006 0.006 0.006 0.006 0.009 0.021];
	output = simf_doi_10_1093_nar_gku1388_figure2B_12(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_1_1(x) - 0.002;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_2_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_3_3(x) - 0.036;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_4_4(x) - 0.05;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_5_5(x) - 0.086;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_6_6(x) - 0.214;
	current_base_index = current_base_index + 1;
	%==============
	arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2];
	GFPmut3b = [0 0.002 0.041 0.479 0.949 0.997 1 1 1];
	output = simf_doi_10_1093_nar_gku593_figureS6_7(x, arabinose) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0.001 0.01 0.03 0.05 0.075 0.1 0.15 0.2 0.3 0.5 1];
	eYFP = [0.015 0.014 0.023 0.062 0.098 0.187 0.239 0.374 0.572 0.869 1];
	output = simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1(x, IPTG) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [9.258 13.887 20.831 32.403 46.29 83.322 162.015 231.45 370.32 648.06 925.8 1388.7];
	eYFP = [0.003 0.003 0.003 0.003 0.005 0.01 0.05 0.167 0.333 0.667 0.833 1];
	output = simf_doi_10_1073pnas_0408507102_figure2A_1(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
end
function error = LOOCV_sim_error_doi_10_1186_1754_1611_3_4_figure3_1(x)
	error = zeros(205,1);
	current_base_index = 0;
	%==============
	IPTG = [0 0.005 0.01 0.02 0.05 0.1 0.2 1];
	lacZ = [0.001 0.003 0.043 0.5 0.984 0.999 1 1];
	output = simf_doi_10_1073_pnas_0606717104_figure4_6_1(x, IPTG) - lacZ;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 1 2 3 4 5 6 7 8 9 10 12 13 17 18 20 25 35 60];
	luc = [0 0.001 0.001 0.001 0.002 0.004 0.007 0.011 0.014 0.032 0.036 0.054 0.089 0.643 0.893 1 1 1 1];
	output = simf_doi__10_1093_nar_25_6_1203_figure4a_1(x, aTc) - luc;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.003 0.005 0.05 0.5 1];
	GFPuv = [0.044 0.089 0.556 1 1];
	output = simf_doi_10_1128_AEM_00791_07_figure3a_1(x, arabinose) - GFPuv;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.003 0.008 0.025 0.075 0.2 0.7 2 6 20 50];
	mCherry = [0.015 0.015 0.015 0.02 0.06 0.48 0.84 0.93 0.96 0.98 0.98 1];
	output = simf_doi_10_1038_nature12148_figure16a_1(x, arabinose) - mCherry;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0 0.002 0.01 0.04 0.16 1 4 25];
	mRFP1 = [0.001 0.001 0.001 0.003 0.045 0.161 0.184 0.185 0.185];
	output = simf_doi_10_1038_nature11516_figureS8A_1(x, arabinose) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0 0 0.001 0.004 0.025 0.1 0.5 2.5];
	mRFP1 = [0.005 0.005 0.005 0.006 0.013 0.09 0.173 0.194 0.196];
	output = simf_doi_10_1038_nature11516_figureS8B_2(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1];
	sfGFP = [0.025 0.025 0.035 0.055 0.115 0.2 0.3 0.4 0.65 1];
	output = simf_doi_10_1038_nbt_2401_figureS11_1(x, IPTG) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.1 1 2 5 7 10 12.5 25 37.5 50 62.5];
	sfGFP = [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425];
	output = simf_doi_10_1038_nbt_2401_figureS13_2(x, arabinose) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.015 0.05 0.5 5 10];
	eYFP = [0.002 0.002 0.057 0.628 0.999 1 1];
	output = simf_doi_10_1038_nature09565_figure1C_1(x, arabinose) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 0.025 0.25 2.5 25 100 250];
	eYFP = [0.005 0.006 0.007 0.03 0.297 0.388 0.398];
	output = simf_doi_10_1038_nature09565_figureS1C_2(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.025 0.05 0.1 0.25 0.5];
	mRFP1 = [0.02 0.031 0.092 0.306 0.745 0.898 0.969 1 1];
	output = simf_doi_10_1038_nbt_2149_Supplement_page17_1(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_3(x) - 0.193;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_4(x) - 0.4;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_5(x) - 0.655;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_6(x) - 0.897;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_7(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_1(x) - 0.422;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_2(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_1(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_2(x) - 0.215;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_3(x) - 0.095;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_4(x) - 0.046;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_5(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_6(x) - 0.003;
	current_base_index = current_base_index + 1;
	%==============
	aTc = [0 2 4 10 25 50 100];
	GFPmut3b = [0 0.002 0.014 0.167 0.666 0.817 0.841];
	output = simf_doi_10_1093_nar_gku1388_figure2B_7(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0 0 0 0.001 0.022 0.176 0.411 0.461];
	output = simf_doi_10_1093_nar_gku1388_figure2B_8(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.005 0.118 0.277];
	output = simf_doi_10_1093_nar_gku1388_figure2B_9(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.001 0.001 0.001 0.001 0.001 0.004 0.039 0.147];
	output = simf_doi_10_1093_nar_gku1388_figure2B_10(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.003 0.004 0.032];
	output = simf_doi_10_1093_nar_gku1388_figure2B_11(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.006 0.006 0.006 0.006 0.006 0.006 0.009 0.021];
	output = simf_doi_10_1093_nar_gku1388_figure2B_12(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_1_1(x) - 0.002;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_2_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_3_3(x) - 0.036;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_4_4(x) - 0.05;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_5_5(x) - 0.086;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_6_6(x) - 0.214;
	current_base_index = current_base_index + 1;
	%==============
	arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2];
	GFPmut3b = [0 0.002 0.041 0.479 0.949 0.997 1 1 1];
	output = simf_doi_10_1093_nar_gku593_figureS6_7(x, arabinose) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0.001 0.01 0.03 0.05 0.075 0.1 0.15 0.2 0.3 0.5 1];
	eYFP = [0.015 0.014 0.023 0.062 0.098 0.187 0.239 0.374 0.572 0.869 1];
	output = simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1(x, IPTG) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [9.258 13.887 20.831 32.403 46.29 83.322 162.015 231.45 370.32 648.06 925.8 1388.7];
	eYFP = [0.003 0.003 0.003 0.003 0.005 0.01 0.05 0.167 0.333 0.667 0.833 1];
	output = simf_doi_10_1073pnas_0408507102_figure2A_1(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
end
function error = LOOCV_sim_error_doi_10_1186_1754_1611_3_4_figure3_2(x)
	error = zeros(205,1);
	current_base_index = 0;
	%==============
	IPTG = [0 0.005 0.01 0.02 0.05 0.1 0.2 1];
	lacZ = [0.001 0.003 0.043 0.5 0.984 0.999 1 1];
	output = simf_doi_10_1073_pnas_0606717104_figure4_6_1(x, IPTG) - lacZ;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 1 2 3 4 5 6 7 8 9 10 12 13 17 18 20 25 35 60];
	luc = [0 0.001 0.001 0.001 0.002 0.004 0.007 0.011 0.014 0.032 0.036 0.054 0.089 0.643 0.893 1 1 1 1];
	output = simf_doi__10_1093_nar_25_6_1203_figure4a_1(x, aTc) - luc;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.003 0.005 0.05 0.5 1];
	GFPuv = [0.044 0.089 0.556 1 1];
	output = simf_doi_10_1128_AEM_00791_07_figure3a_1(x, arabinose) - GFPuv;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.003 0.008 0.025 0.075 0.2 0.7 2 6 20 50];
	mCherry = [0.015 0.015 0.015 0.02 0.06 0.48 0.84 0.93 0.96 0.98 0.98 1];
	output = simf_doi_10_1038_nature12148_figure16a_1(x, arabinose) - mCherry;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0 0.002 0.01 0.04 0.16 1 4 25];
	mRFP1 = [0.001 0.001 0.001 0.003 0.045 0.161 0.184 0.185 0.185];
	output = simf_doi_10_1038_nature11516_figureS8A_1(x, arabinose) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0 0 0.001 0.004 0.025 0.1 0.5 2.5];
	mRFP1 = [0.005 0.005 0.005 0.006 0.013 0.09 0.173 0.194 0.196];
	output = simf_doi_10_1038_nature11516_figureS8B_2(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1];
	sfGFP = [0.025 0.025 0.035 0.055 0.115 0.2 0.3 0.4 0.65 1];
	output = simf_doi_10_1038_nbt_2401_figureS11_1(x, IPTG) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.1 1 2 5 7 10 12.5 25 37.5 50 62.5];
	sfGFP = [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425];
	output = simf_doi_10_1038_nbt_2401_figureS13_2(x, arabinose) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.015 0.05 0.5 5 10];
	eYFP = [0.002 0.002 0.057 0.628 0.999 1 1];
	output = simf_doi_10_1038_nature09565_figure1C_1(x, arabinose) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 0.025 0.25 2.5 25 100 250];
	eYFP = [0.005 0.006 0.007 0.03 0.297 0.388 0.398];
	output = simf_doi_10_1038_nature09565_figureS1C_2(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.025 0.05 0.1 0.25 0.5];
	mRFP1 = [0.02 0.031 0.092 0.306 0.745 0.898 0.969 1 1];
	output = simf_doi_10_1038_nbt_2149_Supplement_page17_1(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_1(x) - 0.69;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_3(x) - 0.193;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_4(x) - 0.4;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_5(x) - 0.655;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_6(x) - 0.897;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_7(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_1(x) - 0.422;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_2(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_1(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_2(x) - 0.215;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_3(x) - 0.095;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_4(x) - 0.046;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_5(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_6(x) - 0.003;
	current_base_index = current_base_index + 1;
	%==============
	aTc = [0 2 4 10 25 50 100];
	GFPmut3b = [0 0.002 0.014 0.167 0.666 0.817 0.841];
	output = simf_doi_10_1093_nar_gku1388_figure2B_7(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0 0 0 0.001 0.022 0.176 0.411 0.461];
	output = simf_doi_10_1093_nar_gku1388_figure2B_8(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.005 0.118 0.277];
	output = simf_doi_10_1093_nar_gku1388_figure2B_9(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.001 0.001 0.001 0.001 0.001 0.004 0.039 0.147];
	output = simf_doi_10_1093_nar_gku1388_figure2B_10(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.003 0.004 0.032];
	output = simf_doi_10_1093_nar_gku1388_figure2B_11(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.006 0.006 0.006 0.006 0.006 0.006 0.009 0.021];
	output = simf_doi_10_1093_nar_gku1388_figure2B_12(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_1_1(x) - 0.002;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_2_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_3_3(x) - 0.036;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_4_4(x) - 0.05;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_5_5(x) - 0.086;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_6_6(x) - 0.214;
	current_base_index = current_base_index + 1;
	%==============
	arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2];
	GFPmut3b = [0 0.002 0.041 0.479 0.949 0.997 1 1 1];
	output = simf_doi_10_1093_nar_gku593_figureS6_7(x, arabinose) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0.001 0.01 0.03 0.05 0.075 0.1 0.15 0.2 0.3 0.5 1];
	eYFP = [0.015 0.014 0.023 0.062 0.098 0.187 0.239 0.374 0.572 0.869 1];
	output = simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1(x, IPTG) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [9.258 13.887 20.831 32.403 46.29 83.322 162.015 231.45 370.32 648.06 925.8 1388.7];
	eYFP = [0.003 0.003 0.003 0.003 0.005 0.01 0.05 0.167 0.333 0.667 0.833 1];
	output = simf_doi_10_1073pnas_0408507102_figure2A_1(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
end
function error = LOOCV_sim_error_doi_10_1186_1754_1611_3_4_figure3_6(x)
	error = zeros(205,1);
	current_base_index = 0;
	%==============
	IPTG = [0 0.005 0.01 0.02 0.05 0.1 0.2 1];
	lacZ = [0.001 0.003 0.043 0.5 0.984 0.999 1 1];
	output = simf_doi_10_1073_pnas_0606717104_figure4_6_1(x, IPTG) - lacZ;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 1 2 3 4 5 6 7 8 9 10 12 13 17 18 20 25 35 60];
	luc = [0 0.001 0.001 0.001 0.002 0.004 0.007 0.011 0.014 0.032 0.036 0.054 0.089 0.643 0.893 1 1 1 1];
	output = simf_doi__10_1093_nar_25_6_1203_figure4a_1(x, aTc) - luc;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.003 0.005 0.05 0.5 1];
	GFPuv = [0.044 0.089 0.556 1 1];
	output = simf_doi_10_1128_AEM_00791_07_figure3a_1(x, arabinose) - GFPuv;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.003 0.008 0.025 0.075 0.2 0.7 2 6 20 50];
	mCherry = [0.015 0.015 0.015 0.02 0.06 0.48 0.84 0.93 0.96 0.98 0.98 1];
	output = simf_doi_10_1038_nature12148_figure16a_1(x, arabinose) - mCherry;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0 0.002 0.01 0.04 0.16 1 4 25];
	mRFP1 = [0.001 0.001 0.001 0.003 0.045 0.161 0.184 0.185 0.185];
	output = simf_doi_10_1038_nature11516_figureS8A_1(x, arabinose) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0 0 0.001 0.004 0.025 0.1 0.5 2.5];
	mRFP1 = [0.005 0.005 0.005 0.006 0.013 0.09 0.173 0.194 0.196];
	output = simf_doi_10_1038_nature11516_figureS8B_2(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1];
	sfGFP = [0.025 0.025 0.035 0.055 0.115 0.2 0.3 0.4 0.65 1];
	output = simf_doi_10_1038_nbt_2401_figureS11_1(x, IPTG) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.1 1 2 5 7 10 12.5 25 37.5 50 62.5];
	sfGFP = [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425];
	output = simf_doi_10_1038_nbt_2401_figureS13_2(x, arabinose) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.015 0.05 0.5 5 10];
	eYFP = [0.002 0.002 0.057 0.628 0.999 1 1];
	output = simf_doi_10_1038_nature09565_figure1C_1(x, arabinose) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 0.025 0.25 2.5 25 100 250];
	eYFP = [0.005 0.006 0.007 0.03 0.297 0.388 0.398];
	output = simf_doi_10_1038_nature09565_figureS1C_2(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.025 0.05 0.1 0.25 0.5];
	mRFP1 = [0.02 0.031 0.092 0.306 0.745 0.898 0.969 1 1];
	output = simf_doi_10_1038_nbt_2149_Supplement_page17_1(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_1(x) - 0.69;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_3(x) - 0.193;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_4(x) - 0.4;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_5(x) - 0.655;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_7(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_1(x) - 0.422;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_2(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_1(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_2(x) - 0.215;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_3(x) - 0.095;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_4(x) - 0.046;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_5(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_6(x) - 0.003;
	current_base_index = current_base_index + 1;
	%==============
	aTc = [0 2 4 10 25 50 100];
	GFPmut3b = [0 0.002 0.014 0.167 0.666 0.817 0.841];
	output = simf_doi_10_1093_nar_gku1388_figure2B_7(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0 0 0 0.001 0.022 0.176 0.411 0.461];
	output = simf_doi_10_1093_nar_gku1388_figure2B_8(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.005 0.118 0.277];
	output = simf_doi_10_1093_nar_gku1388_figure2B_9(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.001 0.001 0.001 0.001 0.001 0.004 0.039 0.147];
	output = simf_doi_10_1093_nar_gku1388_figure2B_10(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.003 0.004 0.032];
	output = simf_doi_10_1093_nar_gku1388_figure2B_11(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.006 0.006 0.006 0.006 0.006 0.006 0.009 0.021];
	output = simf_doi_10_1093_nar_gku1388_figure2B_12(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_1_1(x) - 0.002;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_2_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_3_3(x) - 0.036;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_4_4(x) - 0.05;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_5_5(x) - 0.086;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_6_6(x) - 0.214;
	current_base_index = current_base_index + 1;
	%==============
	arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2];
	GFPmut3b = [0 0.002 0.041 0.479 0.949 0.997 1 1 1];
	output = simf_doi_10_1093_nar_gku593_figureS6_7(x, arabinose) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0.001 0.01 0.03 0.05 0.075 0.1 0.15 0.2 0.3 0.5 1];
	eYFP = [0.015 0.014 0.023 0.062 0.098 0.187 0.239 0.374 0.572 0.869 1];
	output = simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1(x, IPTG) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [9.258 13.887 20.831 32.403 46.29 83.322 162.015 231.45 370.32 648.06 925.8 1388.7];
	eYFP = [0.003 0.003 0.003 0.003 0.005 0.01 0.05 0.167 0.333 0.667 0.833 1];
	output = simf_doi_10_1073pnas_0408507102_figure2A_1(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
end
function error = LOOCV_sim_error_doi_10_1186_1754_1611_3_4_figure3_7(x)
	error = zeros(205,1);
	current_base_index = 0;
	%==============
	IPTG = [0 0.005 0.01 0.02 0.05 0.1 0.2 1];
	lacZ = [0.001 0.003 0.043 0.5 0.984 0.999 1 1];
	output = simf_doi_10_1073_pnas_0606717104_figure4_6_1(x, IPTG) - lacZ;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 1 2 3 4 5 6 7 8 9 10 12 13 17 18 20 25 35 60];
	luc = [0 0.001 0.001 0.001 0.002 0.004 0.007 0.011 0.014 0.032 0.036 0.054 0.089 0.643 0.893 1 1 1 1];
	output = simf_doi__10_1093_nar_25_6_1203_figure4a_1(x, aTc) - luc;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.003 0.005 0.05 0.5 1];
	GFPuv = [0.044 0.089 0.556 1 1];
	output = simf_doi_10_1128_AEM_00791_07_figure3a_1(x, arabinose) - GFPuv;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.003 0.008 0.025 0.075 0.2 0.7 2 6 20 50];
	mCherry = [0.015 0.015 0.015 0.02 0.06 0.48 0.84 0.93 0.96 0.98 0.98 1];
	output = simf_doi_10_1038_nature12148_figure16a_1(x, arabinose) - mCherry;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0 0.002 0.01 0.04 0.16 1 4 25];
	mRFP1 = [0.001 0.001 0.001 0.003 0.045 0.161 0.184 0.185 0.185];
	output = simf_doi_10_1038_nature11516_figureS8A_1(x, arabinose) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0 0 0.001 0.004 0.025 0.1 0.5 2.5];
	mRFP1 = [0.005 0.005 0.005 0.006 0.013 0.09 0.173 0.194 0.196];
	output = simf_doi_10_1038_nature11516_figureS8B_2(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1];
	sfGFP = [0.025 0.025 0.035 0.055 0.115 0.2 0.3 0.4 0.65 1];
	output = simf_doi_10_1038_nbt_2401_figureS11_1(x, IPTG) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.1 1 2 5 7 10 12.5 25 37.5 50 62.5];
	sfGFP = [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425];
	output = simf_doi_10_1038_nbt_2401_figureS13_2(x, arabinose) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.015 0.05 0.5 5 10];
	eYFP = [0.002 0.002 0.057 0.628 0.999 1 1];
	output = simf_doi_10_1038_nature09565_figure1C_1(x, arabinose) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 0.025 0.25 2.5 25 100 250];
	eYFP = [0.005 0.006 0.007 0.03 0.297 0.388 0.398];
	output = simf_doi_10_1038_nature09565_figureS1C_2(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.025 0.05 0.1 0.25 0.5];
	mRFP1 = [0.02 0.031 0.092 0.306 0.745 0.898 0.969 1 1];
	output = simf_doi_10_1038_nbt_2149_Supplement_page17_1(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_1(x) - 0.69;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_3(x) - 0.193;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_4(x) - 0.4;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_5(x) - 0.655;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_6(x) - 0.897;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_1(x) - 0.422;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_2(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_1(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_2(x) - 0.215;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_3(x) - 0.095;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_4(x) - 0.046;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_5(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_6(x) - 0.003;
	current_base_index = current_base_index + 1;
	%==============
	aTc = [0 2 4 10 25 50 100];
	GFPmut3b = [0 0.002 0.014 0.167 0.666 0.817 0.841];
	output = simf_doi_10_1093_nar_gku1388_figure2B_7(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0 0 0 0.001 0.022 0.176 0.411 0.461];
	output = simf_doi_10_1093_nar_gku1388_figure2B_8(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.005 0.118 0.277];
	output = simf_doi_10_1093_nar_gku1388_figure2B_9(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.001 0.001 0.001 0.001 0.001 0.004 0.039 0.147];
	output = simf_doi_10_1093_nar_gku1388_figure2B_10(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.003 0.004 0.032];
	output = simf_doi_10_1093_nar_gku1388_figure2B_11(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.006 0.006 0.006 0.006 0.006 0.006 0.009 0.021];
	output = simf_doi_10_1093_nar_gku1388_figure2B_12(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_1_1(x) - 0.002;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_2_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_3_3(x) - 0.036;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_4_4(x) - 0.05;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_5_5(x) - 0.086;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_6_6(x) - 0.214;
	current_base_index = current_base_index + 1;
	%==============
	arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2];
	GFPmut3b = [0 0.002 0.041 0.479 0.949 0.997 1 1 1];
	output = simf_doi_10_1093_nar_gku593_figureS6_7(x, arabinose) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0.001 0.01 0.03 0.05 0.075 0.1 0.15 0.2 0.3 0.5 1];
	eYFP = [0.015 0.014 0.023 0.062 0.098 0.187 0.239 0.374 0.572 0.869 1];
	output = simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1(x, IPTG) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [9.258 13.887 20.831 32.403 46.29 83.322 162.015 231.45 370.32 648.06 925.8 1388.7];
	eYFP = [0.003 0.003 0.003 0.003 0.005 0.01 0.05 0.167 0.333 0.667 0.833 1];
	output = simf_doi_10_1073pnas_0408507102_figure2A_1(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
end
function error = LOOCV_sim_error_doi__10_1073_pnas_1301301110_sd03_xls_1(x)
	error = zeros(205,1);
	current_base_index = 0;
	%==============
	IPTG = [0 0.005 0.01 0.02 0.05 0.1 0.2 1];
	lacZ = [0.001 0.003 0.043 0.5 0.984 0.999 1 1];
	output = simf_doi_10_1073_pnas_0606717104_figure4_6_1(x, IPTG) - lacZ;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 1 2 3 4 5 6 7 8 9 10 12 13 17 18 20 25 35 60];
	luc = [0 0.001 0.001 0.001 0.002 0.004 0.007 0.011 0.014 0.032 0.036 0.054 0.089 0.643 0.893 1 1 1 1];
	output = simf_doi__10_1093_nar_25_6_1203_figure4a_1(x, aTc) - luc;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.003 0.005 0.05 0.5 1];
	GFPuv = [0.044 0.089 0.556 1 1];
	output = simf_doi_10_1128_AEM_00791_07_figure3a_1(x, arabinose) - GFPuv;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.003 0.008 0.025 0.075 0.2 0.7 2 6 20 50];
	mCherry = [0.015 0.015 0.015 0.02 0.06 0.48 0.84 0.93 0.96 0.98 0.98 1];
	output = simf_doi_10_1038_nature12148_figure16a_1(x, arabinose) - mCherry;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0 0.002 0.01 0.04 0.16 1 4 25];
	mRFP1 = [0.001 0.001 0.001 0.003 0.045 0.161 0.184 0.185 0.185];
	output = simf_doi_10_1038_nature11516_figureS8A_1(x, arabinose) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0 0 0.001 0.004 0.025 0.1 0.5 2.5];
	mRFP1 = [0.005 0.005 0.005 0.006 0.013 0.09 0.173 0.194 0.196];
	output = simf_doi_10_1038_nature11516_figureS8B_2(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1];
	sfGFP = [0.025 0.025 0.035 0.055 0.115 0.2 0.3 0.4 0.65 1];
	output = simf_doi_10_1038_nbt_2401_figureS11_1(x, IPTG) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.1 1 2 5 7 10 12.5 25 37.5 50 62.5];
	sfGFP = [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425];
	output = simf_doi_10_1038_nbt_2401_figureS13_2(x, arabinose) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.015 0.05 0.5 5 10];
	eYFP = [0.002 0.002 0.057 0.628 0.999 1 1];
	output = simf_doi_10_1038_nature09565_figure1C_1(x, arabinose) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 0.025 0.25 2.5 25 100 250];
	eYFP = [0.005 0.006 0.007 0.03 0.297 0.388 0.398];
	output = simf_doi_10_1038_nature09565_figureS1C_2(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.025 0.05 0.1 0.25 0.5];
	mRFP1 = [0.02 0.031 0.092 0.306 0.745 0.898 0.969 1 1];
	output = simf_doi_10_1038_nbt_2149_Supplement_page17_1(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_1(x) - 0.69;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_3(x) - 0.193;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_4(x) - 0.4;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_5(x) - 0.655;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_6(x) - 0.897;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_7(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_2(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_1(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_2(x) - 0.215;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_3(x) - 0.095;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_4(x) - 0.046;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_5(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_6(x) - 0.003;
	current_base_index = current_base_index + 1;
	%==============
	aTc = [0 2 4 10 25 50 100];
	GFPmut3b = [0 0.002 0.014 0.167 0.666 0.817 0.841];
	output = simf_doi_10_1093_nar_gku1388_figure2B_7(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0 0 0 0.001 0.022 0.176 0.411 0.461];
	output = simf_doi_10_1093_nar_gku1388_figure2B_8(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.005 0.118 0.277];
	output = simf_doi_10_1093_nar_gku1388_figure2B_9(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.001 0.001 0.001 0.001 0.001 0.004 0.039 0.147];
	output = simf_doi_10_1093_nar_gku1388_figure2B_10(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.003 0.004 0.032];
	output = simf_doi_10_1093_nar_gku1388_figure2B_11(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.006 0.006 0.006 0.006 0.006 0.006 0.009 0.021];
	output = simf_doi_10_1093_nar_gku1388_figure2B_12(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_1_1(x) - 0.002;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_2_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_3_3(x) - 0.036;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_4_4(x) - 0.05;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_5_5(x) - 0.086;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_6_6(x) - 0.214;
	current_base_index = current_base_index + 1;
	%==============
	arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2];
	GFPmut3b = [0 0.002 0.041 0.479 0.949 0.997 1 1 1];
	output = simf_doi_10_1093_nar_gku593_figureS6_7(x, arabinose) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0.001 0.01 0.03 0.05 0.075 0.1 0.15 0.2 0.3 0.5 1];
	eYFP = [0.015 0.014 0.023 0.062 0.098 0.187 0.239 0.374 0.572 0.869 1];
	output = simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1(x, IPTG) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [9.258 13.887 20.831 32.403 46.29 83.322 162.015 231.45 370.32 648.06 925.8 1388.7];
	eYFP = [0.003 0.003 0.003 0.003 0.005 0.01 0.05 0.167 0.333 0.667 0.833 1];
	output = simf_doi_10_1073pnas_0408507102_figure2A_1(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
end
function error = LOOCV_sim_error_doi__10_1073_pnas_1301301110_sd03_xls_2(x)
	error = zeros(205,1);
	current_base_index = 0;
	%==============
	IPTG = [0 0.005 0.01 0.02 0.05 0.1 0.2 1];
	lacZ = [0.001 0.003 0.043 0.5 0.984 0.999 1 1];
	output = simf_doi_10_1073_pnas_0606717104_figure4_6_1(x, IPTG) - lacZ;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 1 2 3 4 5 6 7 8 9 10 12 13 17 18 20 25 35 60];
	luc = [0 0.001 0.001 0.001 0.002 0.004 0.007 0.011 0.014 0.032 0.036 0.054 0.089 0.643 0.893 1 1 1 1];
	output = simf_doi__10_1093_nar_25_6_1203_figure4a_1(x, aTc) - luc;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.003 0.005 0.05 0.5 1];
	GFPuv = [0.044 0.089 0.556 1 1];
	output = simf_doi_10_1128_AEM_00791_07_figure3a_1(x, arabinose) - GFPuv;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.003 0.008 0.025 0.075 0.2 0.7 2 6 20 50];
	mCherry = [0.015 0.015 0.015 0.02 0.06 0.48 0.84 0.93 0.96 0.98 0.98 1];
	output = simf_doi_10_1038_nature12148_figure16a_1(x, arabinose) - mCherry;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0 0.002 0.01 0.04 0.16 1 4 25];
	mRFP1 = [0.001 0.001 0.001 0.003 0.045 0.161 0.184 0.185 0.185];
	output = simf_doi_10_1038_nature11516_figureS8A_1(x, arabinose) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0 0 0.001 0.004 0.025 0.1 0.5 2.5];
	mRFP1 = [0.005 0.005 0.005 0.006 0.013 0.09 0.173 0.194 0.196];
	output = simf_doi_10_1038_nature11516_figureS8B_2(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1];
	sfGFP = [0.025 0.025 0.035 0.055 0.115 0.2 0.3 0.4 0.65 1];
	output = simf_doi_10_1038_nbt_2401_figureS11_1(x, IPTG) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.1 1 2 5 7 10 12.5 25 37.5 50 62.5];
	sfGFP = [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425];
	output = simf_doi_10_1038_nbt_2401_figureS13_2(x, arabinose) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.015 0.05 0.5 5 10];
	eYFP = [0.002 0.002 0.057 0.628 0.999 1 1];
	output = simf_doi_10_1038_nature09565_figure1C_1(x, arabinose) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 0.025 0.25 2.5 25 100 250];
	eYFP = [0.005 0.006 0.007 0.03 0.297 0.388 0.398];
	output = simf_doi_10_1038_nature09565_figureS1C_2(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.025 0.05 0.1 0.25 0.5];
	mRFP1 = [0.02 0.031 0.092 0.306 0.745 0.898 0.969 1 1];
	output = simf_doi_10_1038_nbt_2149_Supplement_page17_1(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_1(x) - 0.69;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_3(x) - 0.193;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_4(x) - 0.4;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_5(x) - 0.655;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_6(x) - 0.897;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_7(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_1(x) - 0.422;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_1(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_2(x) - 0.215;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_3(x) - 0.095;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_4(x) - 0.046;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_5(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_6(x) - 0.003;
	current_base_index = current_base_index + 1;
	%==============
	aTc = [0 2 4 10 25 50 100];
	GFPmut3b = [0 0.002 0.014 0.167 0.666 0.817 0.841];
	output = simf_doi_10_1093_nar_gku1388_figure2B_7(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0 0 0 0.001 0.022 0.176 0.411 0.461];
	output = simf_doi_10_1093_nar_gku1388_figure2B_8(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.005 0.118 0.277];
	output = simf_doi_10_1093_nar_gku1388_figure2B_9(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.001 0.001 0.001 0.001 0.001 0.004 0.039 0.147];
	output = simf_doi_10_1093_nar_gku1388_figure2B_10(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.003 0.004 0.032];
	output = simf_doi_10_1093_nar_gku1388_figure2B_11(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.006 0.006 0.006 0.006 0.006 0.006 0.009 0.021];
	output = simf_doi_10_1093_nar_gku1388_figure2B_12(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_1_1(x) - 0.002;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_2_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_3_3(x) - 0.036;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_4_4(x) - 0.05;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_5_5(x) - 0.086;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_6_6(x) - 0.214;
	current_base_index = current_base_index + 1;
	%==============
	arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2];
	GFPmut3b = [0 0.002 0.041 0.479 0.949 0.997 1 1 1];
	output = simf_doi_10_1093_nar_gku593_figureS6_7(x, arabinose) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0.001 0.01 0.03 0.05 0.075 0.1 0.15 0.2 0.3 0.5 1];
	eYFP = [0.015 0.014 0.023 0.062 0.098 0.187 0.239 0.374 0.572 0.869 1];
	output = simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1(x, IPTG) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [9.258 13.887 20.831 32.403 46.29 83.322 162.015 231.45 370.32 648.06 925.8 1388.7];
	eYFP = [0.003 0.003 0.003 0.003 0.005 0.01 0.05 0.167 0.333 0.667 0.833 1];
	output = simf_doi_10_1073pnas_0408507102_figure2A_1(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2A_1(x)
	error = zeros(205,1);
	current_base_index = 0;
	%==============
	IPTG = [0 0.005 0.01 0.02 0.05 0.1 0.2 1];
	lacZ = [0.001 0.003 0.043 0.5 0.984 0.999 1 1];
	output = simf_doi_10_1073_pnas_0606717104_figure4_6_1(x, IPTG) - lacZ;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 1 2 3 4 5 6 7 8 9 10 12 13 17 18 20 25 35 60];
	luc = [0 0.001 0.001 0.001 0.002 0.004 0.007 0.011 0.014 0.032 0.036 0.054 0.089 0.643 0.893 1 1 1 1];
	output = simf_doi__10_1093_nar_25_6_1203_figure4a_1(x, aTc) - luc;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.003 0.005 0.05 0.5 1];
	GFPuv = [0.044 0.089 0.556 1 1];
	output = simf_doi_10_1128_AEM_00791_07_figure3a_1(x, arabinose) - GFPuv;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.003 0.008 0.025 0.075 0.2 0.7 2 6 20 50];
	mCherry = [0.015 0.015 0.015 0.02 0.06 0.48 0.84 0.93 0.96 0.98 0.98 1];
	output = simf_doi_10_1038_nature12148_figure16a_1(x, arabinose) - mCherry;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0 0.002 0.01 0.04 0.16 1 4 25];
	mRFP1 = [0.001 0.001 0.001 0.003 0.045 0.161 0.184 0.185 0.185];
	output = simf_doi_10_1038_nature11516_figureS8A_1(x, arabinose) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0 0 0.001 0.004 0.025 0.1 0.5 2.5];
	mRFP1 = [0.005 0.005 0.005 0.006 0.013 0.09 0.173 0.194 0.196];
	output = simf_doi_10_1038_nature11516_figureS8B_2(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1];
	sfGFP = [0.025 0.025 0.035 0.055 0.115 0.2 0.3 0.4 0.65 1];
	output = simf_doi_10_1038_nbt_2401_figureS11_1(x, IPTG) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.1 1 2 5 7 10 12.5 25 37.5 50 62.5];
	sfGFP = [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425];
	output = simf_doi_10_1038_nbt_2401_figureS13_2(x, arabinose) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.015 0.05 0.5 5 10];
	eYFP = [0.002 0.002 0.057 0.628 0.999 1 1];
	output = simf_doi_10_1038_nature09565_figure1C_1(x, arabinose) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 0.025 0.25 2.5 25 100 250];
	eYFP = [0.005 0.006 0.007 0.03 0.297 0.388 0.398];
	output = simf_doi_10_1038_nature09565_figureS1C_2(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.025 0.05 0.1 0.25 0.5];
	mRFP1 = [0.02 0.031 0.092 0.306 0.745 0.898 0.969 1 1];
	output = simf_doi_10_1038_nbt_2149_Supplement_page17_1(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_1(x) - 0.69;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_3(x) - 0.193;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_4(x) - 0.4;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_5(x) - 0.655;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_6(x) - 0.897;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_7(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_1(x) - 0.422;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_2(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_2(x) - 0.215;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_3(x) - 0.095;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_4(x) - 0.046;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_5(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_6(x) - 0.003;
	current_base_index = current_base_index + 1;
	%==============
	aTc = [0 2 4 10 25 50 100];
	GFPmut3b = [0 0.002 0.014 0.167 0.666 0.817 0.841];
	output = simf_doi_10_1093_nar_gku1388_figure2B_7(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0 0 0 0.001 0.022 0.176 0.411 0.461];
	output = simf_doi_10_1093_nar_gku1388_figure2B_8(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.005 0.118 0.277];
	output = simf_doi_10_1093_nar_gku1388_figure2B_9(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.001 0.001 0.001 0.001 0.001 0.004 0.039 0.147];
	output = simf_doi_10_1093_nar_gku1388_figure2B_10(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.003 0.004 0.032];
	output = simf_doi_10_1093_nar_gku1388_figure2B_11(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.006 0.006 0.006 0.006 0.006 0.006 0.009 0.021];
	output = simf_doi_10_1093_nar_gku1388_figure2B_12(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_1_1(x) - 0.002;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_2_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_3_3(x) - 0.036;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_4_4(x) - 0.05;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_5_5(x) - 0.086;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_6_6(x) - 0.214;
	current_base_index = current_base_index + 1;
	%==============
	arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2];
	GFPmut3b = [0 0.002 0.041 0.479 0.949 0.997 1 1 1];
	output = simf_doi_10_1093_nar_gku593_figureS6_7(x, arabinose) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0.001 0.01 0.03 0.05 0.075 0.1 0.15 0.2 0.3 0.5 1];
	eYFP = [0.015 0.014 0.023 0.062 0.098 0.187 0.239 0.374 0.572 0.869 1];
	output = simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1(x, IPTG) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [9.258 13.887 20.831 32.403 46.29 83.322 162.015 231.45 370.32 648.06 925.8 1388.7];
	eYFP = [0.003 0.003 0.003 0.003 0.005 0.01 0.05 0.167 0.333 0.667 0.833 1];
	output = simf_doi_10_1073pnas_0408507102_figure2A_1(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2A_2(x)
	error = zeros(205,1);
	current_base_index = 0;
	%==============
	IPTG = [0 0.005 0.01 0.02 0.05 0.1 0.2 1];
	lacZ = [0.001 0.003 0.043 0.5 0.984 0.999 1 1];
	output = simf_doi_10_1073_pnas_0606717104_figure4_6_1(x, IPTG) - lacZ;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 1 2 3 4 5 6 7 8 9 10 12 13 17 18 20 25 35 60];
	luc = [0 0.001 0.001 0.001 0.002 0.004 0.007 0.011 0.014 0.032 0.036 0.054 0.089 0.643 0.893 1 1 1 1];
	output = simf_doi__10_1093_nar_25_6_1203_figure4a_1(x, aTc) - luc;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.003 0.005 0.05 0.5 1];
	GFPuv = [0.044 0.089 0.556 1 1];
	output = simf_doi_10_1128_AEM_00791_07_figure3a_1(x, arabinose) - GFPuv;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.003 0.008 0.025 0.075 0.2 0.7 2 6 20 50];
	mCherry = [0.015 0.015 0.015 0.02 0.06 0.48 0.84 0.93 0.96 0.98 0.98 1];
	output = simf_doi_10_1038_nature12148_figure16a_1(x, arabinose) - mCherry;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0 0.002 0.01 0.04 0.16 1 4 25];
	mRFP1 = [0.001 0.001 0.001 0.003 0.045 0.161 0.184 0.185 0.185];
	output = simf_doi_10_1038_nature11516_figureS8A_1(x, arabinose) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0 0 0.001 0.004 0.025 0.1 0.5 2.5];
	mRFP1 = [0.005 0.005 0.005 0.006 0.013 0.09 0.173 0.194 0.196];
	output = simf_doi_10_1038_nature11516_figureS8B_2(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1];
	sfGFP = [0.025 0.025 0.035 0.055 0.115 0.2 0.3 0.4 0.65 1];
	output = simf_doi_10_1038_nbt_2401_figureS11_1(x, IPTG) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.1 1 2 5 7 10 12.5 25 37.5 50 62.5];
	sfGFP = [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425];
	output = simf_doi_10_1038_nbt_2401_figureS13_2(x, arabinose) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.015 0.05 0.5 5 10];
	eYFP = [0.002 0.002 0.057 0.628 0.999 1 1];
	output = simf_doi_10_1038_nature09565_figure1C_1(x, arabinose) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 0.025 0.25 2.5 25 100 250];
	eYFP = [0.005 0.006 0.007 0.03 0.297 0.388 0.398];
	output = simf_doi_10_1038_nature09565_figureS1C_2(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.025 0.05 0.1 0.25 0.5];
	mRFP1 = [0.02 0.031 0.092 0.306 0.745 0.898 0.969 1 1];
	output = simf_doi_10_1038_nbt_2149_Supplement_page17_1(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_1(x) - 0.69;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_3(x) - 0.193;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_4(x) - 0.4;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_5(x) - 0.655;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_6(x) - 0.897;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_7(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_1(x) - 0.422;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_2(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_1(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_3(x) - 0.095;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_4(x) - 0.046;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_5(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_6(x) - 0.003;
	current_base_index = current_base_index + 1;
	%==============
	aTc = [0 2 4 10 25 50 100];
	GFPmut3b = [0 0.002 0.014 0.167 0.666 0.817 0.841];
	output = simf_doi_10_1093_nar_gku1388_figure2B_7(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0 0 0 0.001 0.022 0.176 0.411 0.461];
	output = simf_doi_10_1093_nar_gku1388_figure2B_8(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.005 0.118 0.277];
	output = simf_doi_10_1093_nar_gku1388_figure2B_9(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.001 0.001 0.001 0.001 0.001 0.004 0.039 0.147];
	output = simf_doi_10_1093_nar_gku1388_figure2B_10(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.003 0.004 0.032];
	output = simf_doi_10_1093_nar_gku1388_figure2B_11(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.006 0.006 0.006 0.006 0.006 0.006 0.009 0.021];
	output = simf_doi_10_1093_nar_gku1388_figure2B_12(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_1_1(x) - 0.002;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_2_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_3_3(x) - 0.036;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_4_4(x) - 0.05;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_5_5(x) - 0.086;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_6_6(x) - 0.214;
	current_base_index = current_base_index + 1;
	%==============
	arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2];
	GFPmut3b = [0 0.002 0.041 0.479 0.949 0.997 1 1 1];
	output = simf_doi_10_1093_nar_gku593_figureS6_7(x, arabinose) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0.001 0.01 0.03 0.05 0.075 0.1 0.15 0.2 0.3 0.5 1];
	eYFP = [0.015 0.014 0.023 0.062 0.098 0.187 0.239 0.374 0.572 0.869 1];
	output = simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1(x, IPTG) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [9.258 13.887 20.831 32.403 46.29 83.322 162.015 231.45 370.32 648.06 925.8 1388.7];
	eYFP = [0.003 0.003 0.003 0.003 0.005 0.01 0.05 0.167 0.333 0.667 0.833 1];
	output = simf_doi_10_1073pnas_0408507102_figure2A_1(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2A_3(x)
	error = zeros(205,1);
	current_base_index = 0;
	%==============
	IPTG = [0 0.005 0.01 0.02 0.05 0.1 0.2 1];
	lacZ = [0.001 0.003 0.043 0.5 0.984 0.999 1 1];
	output = simf_doi_10_1073_pnas_0606717104_figure4_6_1(x, IPTG) - lacZ;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 1 2 3 4 5 6 7 8 9 10 12 13 17 18 20 25 35 60];
	luc = [0 0.001 0.001 0.001 0.002 0.004 0.007 0.011 0.014 0.032 0.036 0.054 0.089 0.643 0.893 1 1 1 1];
	output = simf_doi__10_1093_nar_25_6_1203_figure4a_1(x, aTc) - luc;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.003 0.005 0.05 0.5 1];
	GFPuv = [0.044 0.089 0.556 1 1];
	output = simf_doi_10_1128_AEM_00791_07_figure3a_1(x, arabinose) - GFPuv;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.003 0.008 0.025 0.075 0.2 0.7 2 6 20 50];
	mCherry = [0.015 0.015 0.015 0.02 0.06 0.48 0.84 0.93 0.96 0.98 0.98 1];
	output = simf_doi_10_1038_nature12148_figure16a_1(x, arabinose) - mCherry;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0 0.002 0.01 0.04 0.16 1 4 25];
	mRFP1 = [0.001 0.001 0.001 0.003 0.045 0.161 0.184 0.185 0.185];
	output = simf_doi_10_1038_nature11516_figureS8A_1(x, arabinose) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0 0 0.001 0.004 0.025 0.1 0.5 2.5];
	mRFP1 = [0.005 0.005 0.005 0.006 0.013 0.09 0.173 0.194 0.196];
	output = simf_doi_10_1038_nature11516_figureS8B_2(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1];
	sfGFP = [0.025 0.025 0.035 0.055 0.115 0.2 0.3 0.4 0.65 1];
	output = simf_doi_10_1038_nbt_2401_figureS11_1(x, IPTG) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.1 1 2 5 7 10 12.5 25 37.5 50 62.5];
	sfGFP = [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425];
	output = simf_doi_10_1038_nbt_2401_figureS13_2(x, arabinose) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.015 0.05 0.5 5 10];
	eYFP = [0.002 0.002 0.057 0.628 0.999 1 1];
	output = simf_doi_10_1038_nature09565_figure1C_1(x, arabinose) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 0.025 0.25 2.5 25 100 250];
	eYFP = [0.005 0.006 0.007 0.03 0.297 0.388 0.398];
	output = simf_doi_10_1038_nature09565_figureS1C_2(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.025 0.05 0.1 0.25 0.5];
	mRFP1 = [0.02 0.031 0.092 0.306 0.745 0.898 0.969 1 1];
	output = simf_doi_10_1038_nbt_2149_Supplement_page17_1(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_1(x) - 0.69;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_3(x) - 0.193;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_4(x) - 0.4;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_5(x) - 0.655;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_6(x) - 0.897;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_7(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_1(x) - 0.422;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_2(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_1(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_2(x) - 0.215;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_4(x) - 0.046;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_5(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_6(x) - 0.003;
	current_base_index = current_base_index + 1;
	%==============
	aTc = [0 2 4 10 25 50 100];
	GFPmut3b = [0 0.002 0.014 0.167 0.666 0.817 0.841];
	output = simf_doi_10_1093_nar_gku1388_figure2B_7(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0 0 0 0.001 0.022 0.176 0.411 0.461];
	output = simf_doi_10_1093_nar_gku1388_figure2B_8(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.005 0.118 0.277];
	output = simf_doi_10_1093_nar_gku1388_figure2B_9(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.001 0.001 0.001 0.001 0.001 0.004 0.039 0.147];
	output = simf_doi_10_1093_nar_gku1388_figure2B_10(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.003 0.004 0.032];
	output = simf_doi_10_1093_nar_gku1388_figure2B_11(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.006 0.006 0.006 0.006 0.006 0.006 0.009 0.021];
	output = simf_doi_10_1093_nar_gku1388_figure2B_12(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_1_1(x) - 0.002;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_2_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_3_3(x) - 0.036;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_4_4(x) - 0.05;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_5_5(x) - 0.086;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_6_6(x) - 0.214;
	current_base_index = current_base_index + 1;
	%==============
	arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2];
	GFPmut3b = [0 0.002 0.041 0.479 0.949 0.997 1 1 1];
	output = simf_doi_10_1093_nar_gku593_figureS6_7(x, arabinose) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0.001 0.01 0.03 0.05 0.075 0.1 0.15 0.2 0.3 0.5 1];
	eYFP = [0.015 0.014 0.023 0.062 0.098 0.187 0.239 0.374 0.572 0.869 1];
	output = simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1(x, IPTG) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [9.258 13.887 20.831 32.403 46.29 83.322 162.015 231.45 370.32 648.06 925.8 1388.7];
	eYFP = [0.003 0.003 0.003 0.003 0.005 0.01 0.05 0.167 0.333 0.667 0.833 1];
	output = simf_doi_10_1073pnas_0408507102_figure2A_1(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2A_4(x)
	error = zeros(205,1);
	current_base_index = 0;
	%==============
	IPTG = [0 0.005 0.01 0.02 0.05 0.1 0.2 1];
	lacZ = [0.001 0.003 0.043 0.5 0.984 0.999 1 1];
	output = simf_doi_10_1073_pnas_0606717104_figure4_6_1(x, IPTG) - lacZ;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 1 2 3 4 5 6 7 8 9 10 12 13 17 18 20 25 35 60];
	luc = [0 0.001 0.001 0.001 0.002 0.004 0.007 0.011 0.014 0.032 0.036 0.054 0.089 0.643 0.893 1 1 1 1];
	output = simf_doi__10_1093_nar_25_6_1203_figure4a_1(x, aTc) - luc;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.003 0.005 0.05 0.5 1];
	GFPuv = [0.044 0.089 0.556 1 1];
	output = simf_doi_10_1128_AEM_00791_07_figure3a_1(x, arabinose) - GFPuv;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.003 0.008 0.025 0.075 0.2 0.7 2 6 20 50];
	mCherry = [0.015 0.015 0.015 0.02 0.06 0.48 0.84 0.93 0.96 0.98 0.98 1];
	output = simf_doi_10_1038_nature12148_figure16a_1(x, arabinose) - mCherry;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0 0.002 0.01 0.04 0.16 1 4 25];
	mRFP1 = [0.001 0.001 0.001 0.003 0.045 0.161 0.184 0.185 0.185];
	output = simf_doi_10_1038_nature11516_figureS8A_1(x, arabinose) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0 0 0.001 0.004 0.025 0.1 0.5 2.5];
	mRFP1 = [0.005 0.005 0.005 0.006 0.013 0.09 0.173 0.194 0.196];
	output = simf_doi_10_1038_nature11516_figureS8B_2(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1];
	sfGFP = [0.025 0.025 0.035 0.055 0.115 0.2 0.3 0.4 0.65 1];
	output = simf_doi_10_1038_nbt_2401_figureS11_1(x, IPTG) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.1 1 2 5 7 10 12.5 25 37.5 50 62.5];
	sfGFP = [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425];
	output = simf_doi_10_1038_nbt_2401_figureS13_2(x, arabinose) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.015 0.05 0.5 5 10];
	eYFP = [0.002 0.002 0.057 0.628 0.999 1 1];
	output = simf_doi_10_1038_nature09565_figure1C_1(x, arabinose) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 0.025 0.25 2.5 25 100 250];
	eYFP = [0.005 0.006 0.007 0.03 0.297 0.388 0.398];
	output = simf_doi_10_1038_nature09565_figureS1C_2(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.025 0.05 0.1 0.25 0.5];
	mRFP1 = [0.02 0.031 0.092 0.306 0.745 0.898 0.969 1 1];
	output = simf_doi_10_1038_nbt_2149_Supplement_page17_1(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_1(x) - 0.69;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_3(x) - 0.193;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_4(x) - 0.4;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_5(x) - 0.655;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_6(x) - 0.897;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_7(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_1(x) - 0.422;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_2(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_1(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_2(x) - 0.215;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_3(x) - 0.095;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_5(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_6(x) - 0.003;
	current_base_index = current_base_index + 1;
	%==============
	aTc = [0 2 4 10 25 50 100];
	GFPmut3b = [0 0.002 0.014 0.167 0.666 0.817 0.841];
	output = simf_doi_10_1093_nar_gku1388_figure2B_7(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0 0 0 0.001 0.022 0.176 0.411 0.461];
	output = simf_doi_10_1093_nar_gku1388_figure2B_8(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.005 0.118 0.277];
	output = simf_doi_10_1093_nar_gku1388_figure2B_9(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.001 0.001 0.001 0.001 0.001 0.004 0.039 0.147];
	output = simf_doi_10_1093_nar_gku1388_figure2B_10(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.003 0.004 0.032];
	output = simf_doi_10_1093_nar_gku1388_figure2B_11(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.006 0.006 0.006 0.006 0.006 0.006 0.009 0.021];
	output = simf_doi_10_1093_nar_gku1388_figure2B_12(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_1_1(x) - 0.002;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_2_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_3_3(x) - 0.036;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_4_4(x) - 0.05;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_5_5(x) - 0.086;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_6_6(x) - 0.214;
	current_base_index = current_base_index + 1;
	%==============
	arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2];
	GFPmut3b = [0 0.002 0.041 0.479 0.949 0.997 1 1 1];
	output = simf_doi_10_1093_nar_gku593_figureS6_7(x, arabinose) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0.001 0.01 0.03 0.05 0.075 0.1 0.15 0.2 0.3 0.5 1];
	eYFP = [0.015 0.014 0.023 0.062 0.098 0.187 0.239 0.374 0.572 0.869 1];
	output = simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1(x, IPTG) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [9.258 13.887 20.831 32.403 46.29 83.322 162.015 231.45 370.32 648.06 925.8 1388.7];
	eYFP = [0.003 0.003 0.003 0.003 0.005 0.01 0.05 0.167 0.333 0.667 0.833 1];
	output = simf_doi_10_1073pnas_0408507102_figure2A_1(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2A_5(x)
	error = zeros(205,1);
	current_base_index = 0;
	%==============
	IPTG = [0 0.005 0.01 0.02 0.05 0.1 0.2 1];
	lacZ = [0.001 0.003 0.043 0.5 0.984 0.999 1 1];
	output = simf_doi_10_1073_pnas_0606717104_figure4_6_1(x, IPTG) - lacZ;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 1 2 3 4 5 6 7 8 9 10 12 13 17 18 20 25 35 60];
	luc = [0 0.001 0.001 0.001 0.002 0.004 0.007 0.011 0.014 0.032 0.036 0.054 0.089 0.643 0.893 1 1 1 1];
	output = simf_doi__10_1093_nar_25_6_1203_figure4a_1(x, aTc) - luc;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.003 0.005 0.05 0.5 1];
	GFPuv = [0.044 0.089 0.556 1 1];
	output = simf_doi_10_1128_AEM_00791_07_figure3a_1(x, arabinose) - GFPuv;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.003 0.008 0.025 0.075 0.2 0.7 2 6 20 50];
	mCherry = [0.015 0.015 0.015 0.02 0.06 0.48 0.84 0.93 0.96 0.98 0.98 1];
	output = simf_doi_10_1038_nature12148_figure16a_1(x, arabinose) - mCherry;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0 0.002 0.01 0.04 0.16 1 4 25];
	mRFP1 = [0.001 0.001 0.001 0.003 0.045 0.161 0.184 0.185 0.185];
	output = simf_doi_10_1038_nature11516_figureS8A_1(x, arabinose) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0 0 0.001 0.004 0.025 0.1 0.5 2.5];
	mRFP1 = [0.005 0.005 0.005 0.006 0.013 0.09 0.173 0.194 0.196];
	output = simf_doi_10_1038_nature11516_figureS8B_2(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1];
	sfGFP = [0.025 0.025 0.035 0.055 0.115 0.2 0.3 0.4 0.65 1];
	output = simf_doi_10_1038_nbt_2401_figureS11_1(x, IPTG) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.1 1 2 5 7 10 12.5 25 37.5 50 62.5];
	sfGFP = [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425];
	output = simf_doi_10_1038_nbt_2401_figureS13_2(x, arabinose) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.015 0.05 0.5 5 10];
	eYFP = [0.002 0.002 0.057 0.628 0.999 1 1];
	output = simf_doi_10_1038_nature09565_figure1C_1(x, arabinose) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 0.025 0.25 2.5 25 100 250];
	eYFP = [0.005 0.006 0.007 0.03 0.297 0.388 0.398];
	output = simf_doi_10_1038_nature09565_figureS1C_2(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.025 0.05 0.1 0.25 0.5];
	mRFP1 = [0.02 0.031 0.092 0.306 0.745 0.898 0.969 1 1];
	output = simf_doi_10_1038_nbt_2149_Supplement_page17_1(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_1(x) - 0.69;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_3(x) - 0.193;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_4(x) - 0.4;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_5(x) - 0.655;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_6(x) - 0.897;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_7(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_1(x) - 0.422;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_2(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_1(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_2(x) - 0.215;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_3(x) - 0.095;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_4(x) - 0.046;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_6(x) - 0.003;
	current_base_index = current_base_index + 1;
	%==============
	aTc = [0 2 4 10 25 50 100];
	GFPmut3b = [0 0.002 0.014 0.167 0.666 0.817 0.841];
	output = simf_doi_10_1093_nar_gku1388_figure2B_7(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0 0 0 0.001 0.022 0.176 0.411 0.461];
	output = simf_doi_10_1093_nar_gku1388_figure2B_8(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.005 0.118 0.277];
	output = simf_doi_10_1093_nar_gku1388_figure2B_9(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.001 0.001 0.001 0.001 0.001 0.004 0.039 0.147];
	output = simf_doi_10_1093_nar_gku1388_figure2B_10(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.003 0.004 0.032];
	output = simf_doi_10_1093_nar_gku1388_figure2B_11(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.006 0.006 0.006 0.006 0.006 0.006 0.009 0.021];
	output = simf_doi_10_1093_nar_gku1388_figure2B_12(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_1_1(x) - 0.002;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_2_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_3_3(x) - 0.036;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_4_4(x) - 0.05;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_5_5(x) - 0.086;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_6_6(x) - 0.214;
	current_base_index = current_base_index + 1;
	%==============
	arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2];
	GFPmut3b = [0 0.002 0.041 0.479 0.949 0.997 1 1 1];
	output = simf_doi_10_1093_nar_gku593_figureS6_7(x, arabinose) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0.001 0.01 0.03 0.05 0.075 0.1 0.15 0.2 0.3 0.5 1];
	eYFP = [0.015 0.014 0.023 0.062 0.098 0.187 0.239 0.374 0.572 0.869 1];
	output = simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1(x, IPTG) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [9.258 13.887 20.831 32.403 46.29 83.322 162.015 231.45 370.32 648.06 925.8 1388.7];
	eYFP = [0.003 0.003 0.003 0.003 0.005 0.01 0.05 0.167 0.333 0.667 0.833 1];
	output = simf_doi_10_1073pnas_0408507102_figure2A_1(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2A_6(x)
	error = zeros(205,1);
	current_base_index = 0;
	%==============
	IPTG = [0 0.005 0.01 0.02 0.05 0.1 0.2 1];
	lacZ = [0.001 0.003 0.043 0.5 0.984 0.999 1 1];
	output = simf_doi_10_1073_pnas_0606717104_figure4_6_1(x, IPTG) - lacZ;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 1 2 3 4 5 6 7 8 9 10 12 13 17 18 20 25 35 60];
	luc = [0 0.001 0.001 0.001 0.002 0.004 0.007 0.011 0.014 0.032 0.036 0.054 0.089 0.643 0.893 1 1 1 1];
	output = simf_doi__10_1093_nar_25_6_1203_figure4a_1(x, aTc) - luc;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.003 0.005 0.05 0.5 1];
	GFPuv = [0.044 0.089 0.556 1 1];
	output = simf_doi_10_1128_AEM_00791_07_figure3a_1(x, arabinose) - GFPuv;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.003 0.008 0.025 0.075 0.2 0.7 2 6 20 50];
	mCherry = [0.015 0.015 0.015 0.02 0.06 0.48 0.84 0.93 0.96 0.98 0.98 1];
	output = simf_doi_10_1038_nature12148_figure16a_1(x, arabinose) - mCherry;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0 0.002 0.01 0.04 0.16 1 4 25];
	mRFP1 = [0.001 0.001 0.001 0.003 0.045 0.161 0.184 0.185 0.185];
	output = simf_doi_10_1038_nature11516_figureS8A_1(x, arabinose) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0 0 0.001 0.004 0.025 0.1 0.5 2.5];
	mRFP1 = [0.005 0.005 0.005 0.006 0.013 0.09 0.173 0.194 0.196];
	output = simf_doi_10_1038_nature11516_figureS8B_2(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1];
	sfGFP = [0.025 0.025 0.035 0.055 0.115 0.2 0.3 0.4 0.65 1];
	output = simf_doi_10_1038_nbt_2401_figureS11_1(x, IPTG) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.1 1 2 5 7 10 12.5 25 37.5 50 62.5];
	sfGFP = [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425];
	output = simf_doi_10_1038_nbt_2401_figureS13_2(x, arabinose) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.015 0.05 0.5 5 10];
	eYFP = [0.002 0.002 0.057 0.628 0.999 1 1];
	output = simf_doi_10_1038_nature09565_figure1C_1(x, arabinose) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 0.025 0.25 2.5 25 100 250];
	eYFP = [0.005 0.006 0.007 0.03 0.297 0.388 0.398];
	output = simf_doi_10_1038_nature09565_figureS1C_2(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.025 0.05 0.1 0.25 0.5];
	mRFP1 = [0.02 0.031 0.092 0.306 0.745 0.898 0.969 1 1];
	output = simf_doi_10_1038_nbt_2149_Supplement_page17_1(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_1(x) - 0.69;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_3(x) - 0.193;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_4(x) - 0.4;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_5(x) - 0.655;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_6(x) - 0.897;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_7(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_1(x) - 0.422;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_2(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_1(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_2(x) - 0.215;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_3(x) - 0.095;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_4(x) - 0.046;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_5(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	aTc = [0 2 4 10 25 50 100];
	GFPmut3b = [0 0.002 0.014 0.167 0.666 0.817 0.841];
	output = simf_doi_10_1093_nar_gku1388_figure2B_7(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0 0 0 0.001 0.022 0.176 0.411 0.461];
	output = simf_doi_10_1093_nar_gku1388_figure2B_8(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.005 0.118 0.277];
	output = simf_doi_10_1093_nar_gku1388_figure2B_9(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.001 0.001 0.001 0.001 0.001 0.004 0.039 0.147];
	output = simf_doi_10_1093_nar_gku1388_figure2B_10(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.003 0.004 0.032];
	output = simf_doi_10_1093_nar_gku1388_figure2B_11(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.006 0.006 0.006 0.006 0.006 0.006 0.009 0.021];
	output = simf_doi_10_1093_nar_gku1388_figure2B_12(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_1_1(x) - 0.002;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_2_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_3_3(x) - 0.036;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_4_4(x) - 0.05;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_5_5(x) - 0.086;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_6_6(x) - 0.214;
	current_base_index = current_base_index + 1;
	%==============
	arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2];
	GFPmut3b = [0 0.002 0.041 0.479 0.949 0.997 1 1 1];
	output = simf_doi_10_1093_nar_gku593_figureS6_7(x, arabinose) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0.001 0.01 0.03 0.05 0.075 0.1 0.15 0.2 0.3 0.5 1];
	eYFP = [0.015 0.014 0.023 0.062 0.098 0.187 0.239 0.374 0.572 0.869 1];
	output = simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1(x, IPTG) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [9.258 13.887 20.831 32.403 46.29 83.322 162.015 231.45 370.32 648.06 925.8 1388.7];
	eYFP = [0.003 0.003 0.003 0.003 0.005 0.01 0.05 0.167 0.333 0.667 0.833 1];
	output = simf_doi_10_1073pnas_0408507102_figure2A_1(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2B_7(x)
	error = zeros(199,1);
	current_base_index = 0;
	%==============
	IPTG = [0 0.005 0.01 0.02 0.05 0.1 0.2 1];
	lacZ = [0.001 0.003 0.043 0.5 0.984 0.999 1 1];
	output = simf_doi_10_1073_pnas_0606717104_figure4_6_1(x, IPTG) - lacZ;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 1 2 3 4 5 6 7 8 9 10 12 13 17 18 20 25 35 60];
	luc = [0 0.001 0.001 0.001 0.002 0.004 0.007 0.011 0.014 0.032 0.036 0.054 0.089 0.643 0.893 1 1 1 1];
	output = simf_doi__10_1093_nar_25_6_1203_figure4a_1(x, aTc) - luc;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.003 0.005 0.05 0.5 1];
	GFPuv = [0.044 0.089 0.556 1 1];
	output = simf_doi_10_1128_AEM_00791_07_figure3a_1(x, arabinose) - GFPuv;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.003 0.008 0.025 0.075 0.2 0.7 2 6 20 50];
	mCherry = [0.015 0.015 0.015 0.02 0.06 0.48 0.84 0.93 0.96 0.98 0.98 1];
	output = simf_doi_10_1038_nature12148_figure16a_1(x, arabinose) - mCherry;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0 0.002 0.01 0.04 0.16 1 4 25];
	mRFP1 = [0.001 0.001 0.001 0.003 0.045 0.161 0.184 0.185 0.185];
	output = simf_doi_10_1038_nature11516_figureS8A_1(x, arabinose) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0 0 0.001 0.004 0.025 0.1 0.5 2.5];
	mRFP1 = [0.005 0.005 0.005 0.006 0.013 0.09 0.173 0.194 0.196];
	output = simf_doi_10_1038_nature11516_figureS8B_2(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1];
	sfGFP = [0.025 0.025 0.035 0.055 0.115 0.2 0.3 0.4 0.65 1];
	output = simf_doi_10_1038_nbt_2401_figureS11_1(x, IPTG) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.1 1 2 5 7 10 12.5 25 37.5 50 62.5];
	sfGFP = [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425];
	output = simf_doi_10_1038_nbt_2401_figureS13_2(x, arabinose) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.015 0.05 0.5 5 10];
	eYFP = [0.002 0.002 0.057 0.628 0.999 1 1];
	output = simf_doi_10_1038_nature09565_figure1C_1(x, arabinose) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 0.025 0.25 2.5 25 100 250];
	eYFP = [0.005 0.006 0.007 0.03 0.297 0.388 0.398];
	output = simf_doi_10_1038_nature09565_figureS1C_2(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.025 0.05 0.1 0.25 0.5];
	mRFP1 = [0.02 0.031 0.092 0.306 0.745 0.898 0.969 1 1];
	output = simf_doi_10_1038_nbt_2149_Supplement_page17_1(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_1(x) - 0.69;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_3(x) - 0.193;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_4(x) - 0.4;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_5(x) - 0.655;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_6(x) - 0.897;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_7(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_1(x) - 0.422;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_2(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_1(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_2(x) - 0.215;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_3(x) - 0.095;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_4(x) - 0.046;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_5(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_6(x) - 0.003;
	current_base_index = current_base_index + 1;
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0 0 0 0.001 0.022 0.176 0.411 0.461];
	output = simf_doi_10_1093_nar_gku1388_figure2B_8(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.005 0.118 0.277];
	output = simf_doi_10_1093_nar_gku1388_figure2B_9(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.001 0.001 0.001 0.001 0.001 0.004 0.039 0.147];
	output = simf_doi_10_1093_nar_gku1388_figure2B_10(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.003 0.004 0.032];
	output = simf_doi_10_1093_nar_gku1388_figure2B_11(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.006 0.006 0.006 0.006 0.006 0.006 0.009 0.021];
	output = simf_doi_10_1093_nar_gku1388_figure2B_12(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_1_1(x) - 0.002;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_2_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_3_3(x) - 0.036;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_4_4(x) - 0.05;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_5_5(x) - 0.086;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_6_6(x) - 0.214;
	current_base_index = current_base_index + 1;
	%==============
	arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2];
	GFPmut3b = [0 0.002 0.041 0.479 0.949 0.997 1 1 1];
	output = simf_doi_10_1093_nar_gku593_figureS6_7(x, arabinose) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0.001 0.01 0.03 0.05 0.075 0.1 0.15 0.2 0.3 0.5 1];
	eYFP = [0.015 0.014 0.023 0.062 0.098 0.187 0.239 0.374 0.572 0.869 1];
	output = simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1(x, IPTG) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [9.258 13.887 20.831 32.403 46.29 83.322 162.015 231.45 370.32 648.06 925.8 1388.7];
	eYFP = [0.003 0.003 0.003 0.003 0.005 0.01 0.05 0.167 0.333 0.667 0.833 1];
	output = simf_doi_10_1073pnas_0408507102_figure2A_1(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2B_8(x)
	error = zeros(198,1);
	current_base_index = 0;
	%==============
	IPTG = [0 0.005 0.01 0.02 0.05 0.1 0.2 1];
	lacZ = [0.001 0.003 0.043 0.5 0.984 0.999 1 1];
	output = simf_doi_10_1073_pnas_0606717104_figure4_6_1(x, IPTG) - lacZ;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 1 2 3 4 5 6 7 8 9 10 12 13 17 18 20 25 35 60];
	luc = [0 0.001 0.001 0.001 0.002 0.004 0.007 0.011 0.014 0.032 0.036 0.054 0.089 0.643 0.893 1 1 1 1];
	output = simf_doi__10_1093_nar_25_6_1203_figure4a_1(x, aTc) - luc;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.003 0.005 0.05 0.5 1];
	GFPuv = [0.044 0.089 0.556 1 1];
	output = simf_doi_10_1128_AEM_00791_07_figure3a_1(x, arabinose) - GFPuv;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.003 0.008 0.025 0.075 0.2 0.7 2 6 20 50];
	mCherry = [0.015 0.015 0.015 0.02 0.06 0.48 0.84 0.93 0.96 0.98 0.98 1];
	output = simf_doi_10_1038_nature12148_figure16a_1(x, arabinose) - mCherry;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0 0.002 0.01 0.04 0.16 1 4 25];
	mRFP1 = [0.001 0.001 0.001 0.003 0.045 0.161 0.184 0.185 0.185];
	output = simf_doi_10_1038_nature11516_figureS8A_1(x, arabinose) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0 0 0.001 0.004 0.025 0.1 0.5 2.5];
	mRFP1 = [0.005 0.005 0.005 0.006 0.013 0.09 0.173 0.194 0.196];
	output = simf_doi_10_1038_nature11516_figureS8B_2(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1];
	sfGFP = [0.025 0.025 0.035 0.055 0.115 0.2 0.3 0.4 0.65 1];
	output = simf_doi_10_1038_nbt_2401_figureS11_1(x, IPTG) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.1 1 2 5 7 10 12.5 25 37.5 50 62.5];
	sfGFP = [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425];
	output = simf_doi_10_1038_nbt_2401_figureS13_2(x, arabinose) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.015 0.05 0.5 5 10];
	eYFP = [0.002 0.002 0.057 0.628 0.999 1 1];
	output = simf_doi_10_1038_nature09565_figure1C_1(x, arabinose) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 0.025 0.25 2.5 25 100 250];
	eYFP = [0.005 0.006 0.007 0.03 0.297 0.388 0.398];
	output = simf_doi_10_1038_nature09565_figureS1C_2(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.025 0.05 0.1 0.25 0.5];
	mRFP1 = [0.02 0.031 0.092 0.306 0.745 0.898 0.969 1 1];
	output = simf_doi_10_1038_nbt_2149_Supplement_page17_1(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_1(x) - 0.69;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_3(x) - 0.193;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_4(x) - 0.4;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_5(x) - 0.655;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_6(x) - 0.897;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_7(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_1(x) - 0.422;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_2(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_1(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_2(x) - 0.215;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_3(x) - 0.095;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_4(x) - 0.046;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_5(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_6(x) - 0.003;
	current_base_index = current_base_index + 1;
	%==============
	aTc = [0 2 4 10 25 50 100];
	GFPmut3b = [0 0.002 0.014 0.167 0.666 0.817 0.841];
	output = simf_doi_10_1093_nar_gku1388_figure2B_7(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.005 0.118 0.277];
	output = simf_doi_10_1093_nar_gku1388_figure2B_9(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.001 0.001 0.001 0.001 0.001 0.004 0.039 0.147];
	output = simf_doi_10_1093_nar_gku1388_figure2B_10(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.003 0.004 0.032];
	output = simf_doi_10_1093_nar_gku1388_figure2B_11(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.006 0.006 0.006 0.006 0.006 0.006 0.009 0.021];
	output = simf_doi_10_1093_nar_gku1388_figure2B_12(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_1_1(x) - 0.002;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_2_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_3_3(x) - 0.036;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_4_4(x) - 0.05;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_5_5(x) - 0.086;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_6_6(x) - 0.214;
	current_base_index = current_base_index + 1;
	%==============
	arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2];
	GFPmut3b = [0 0.002 0.041 0.479 0.949 0.997 1 1 1];
	output = simf_doi_10_1093_nar_gku593_figureS6_7(x, arabinose) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0.001 0.01 0.03 0.05 0.075 0.1 0.15 0.2 0.3 0.5 1];
	eYFP = [0.015 0.014 0.023 0.062 0.098 0.187 0.239 0.374 0.572 0.869 1];
	output = simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1(x, IPTG) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [9.258 13.887 20.831 32.403 46.29 83.322 162.015 231.45 370.32 648.06 925.8 1388.7];
	eYFP = [0.003 0.003 0.003 0.003 0.005 0.01 0.05 0.167 0.333 0.667 0.833 1];
	output = simf_doi_10_1073pnas_0408507102_figure2A_1(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2B_9(x)
	error = zeros(198,1);
	current_base_index = 0;
	%==============
	IPTG = [0 0.005 0.01 0.02 0.05 0.1 0.2 1];
	lacZ = [0.001 0.003 0.043 0.5 0.984 0.999 1 1];
	output = simf_doi_10_1073_pnas_0606717104_figure4_6_1(x, IPTG) - lacZ;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 1 2 3 4 5 6 7 8 9 10 12 13 17 18 20 25 35 60];
	luc = [0 0.001 0.001 0.001 0.002 0.004 0.007 0.011 0.014 0.032 0.036 0.054 0.089 0.643 0.893 1 1 1 1];
	output = simf_doi__10_1093_nar_25_6_1203_figure4a_1(x, aTc) - luc;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.003 0.005 0.05 0.5 1];
	GFPuv = [0.044 0.089 0.556 1 1];
	output = simf_doi_10_1128_AEM_00791_07_figure3a_1(x, arabinose) - GFPuv;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.003 0.008 0.025 0.075 0.2 0.7 2 6 20 50];
	mCherry = [0.015 0.015 0.015 0.02 0.06 0.48 0.84 0.93 0.96 0.98 0.98 1];
	output = simf_doi_10_1038_nature12148_figure16a_1(x, arabinose) - mCherry;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0 0.002 0.01 0.04 0.16 1 4 25];
	mRFP1 = [0.001 0.001 0.001 0.003 0.045 0.161 0.184 0.185 0.185];
	output = simf_doi_10_1038_nature11516_figureS8A_1(x, arabinose) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0 0 0.001 0.004 0.025 0.1 0.5 2.5];
	mRFP1 = [0.005 0.005 0.005 0.006 0.013 0.09 0.173 0.194 0.196];
	output = simf_doi_10_1038_nature11516_figureS8B_2(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1];
	sfGFP = [0.025 0.025 0.035 0.055 0.115 0.2 0.3 0.4 0.65 1];
	output = simf_doi_10_1038_nbt_2401_figureS11_1(x, IPTG) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.1 1 2 5 7 10 12.5 25 37.5 50 62.5];
	sfGFP = [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425];
	output = simf_doi_10_1038_nbt_2401_figureS13_2(x, arabinose) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.015 0.05 0.5 5 10];
	eYFP = [0.002 0.002 0.057 0.628 0.999 1 1];
	output = simf_doi_10_1038_nature09565_figure1C_1(x, arabinose) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 0.025 0.25 2.5 25 100 250];
	eYFP = [0.005 0.006 0.007 0.03 0.297 0.388 0.398];
	output = simf_doi_10_1038_nature09565_figureS1C_2(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.025 0.05 0.1 0.25 0.5];
	mRFP1 = [0.02 0.031 0.092 0.306 0.745 0.898 0.969 1 1];
	output = simf_doi_10_1038_nbt_2149_Supplement_page17_1(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_1(x) - 0.69;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_3(x) - 0.193;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_4(x) - 0.4;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_5(x) - 0.655;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_6(x) - 0.897;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_7(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_1(x) - 0.422;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_2(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_1(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_2(x) - 0.215;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_3(x) - 0.095;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_4(x) - 0.046;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_5(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_6(x) - 0.003;
	current_base_index = current_base_index + 1;
	%==============
	aTc = [0 2 4 10 25 50 100];
	GFPmut3b = [0 0.002 0.014 0.167 0.666 0.817 0.841];
	output = simf_doi_10_1093_nar_gku1388_figure2B_7(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0 0 0 0.001 0.022 0.176 0.411 0.461];
	output = simf_doi_10_1093_nar_gku1388_figure2B_8(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.001 0.001 0.001 0.001 0.001 0.004 0.039 0.147];
	output = simf_doi_10_1093_nar_gku1388_figure2B_10(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.003 0.004 0.032];
	output = simf_doi_10_1093_nar_gku1388_figure2B_11(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.006 0.006 0.006 0.006 0.006 0.006 0.009 0.021];
	output = simf_doi_10_1093_nar_gku1388_figure2B_12(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_1_1(x) - 0.002;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_2_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_3_3(x) - 0.036;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_4_4(x) - 0.05;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_5_5(x) - 0.086;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_6_6(x) - 0.214;
	current_base_index = current_base_index + 1;
	%==============
	arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2];
	GFPmut3b = [0 0.002 0.041 0.479 0.949 0.997 1 1 1];
	output = simf_doi_10_1093_nar_gku593_figureS6_7(x, arabinose) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0.001 0.01 0.03 0.05 0.075 0.1 0.15 0.2 0.3 0.5 1];
	eYFP = [0.015 0.014 0.023 0.062 0.098 0.187 0.239 0.374 0.572 0.869 1];
	output = simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1(x, IPTG) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [9.258 13.887 20.831 32.403 46.29 83.322 162.015 231.45 370.32 648.06 925.8 1388.7];
	eYFP = [0.003 0.003 0.003 0.003 0.005 0.01 0.05 0.167 0.333 0.667 0.833 1];
	output = simf_doi_10_1073pnas_0408507102_figure2A_1(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2B_10(x)
	error = zeros(198,1);
	current_base_index = 0;
	%==============
	IPTG = [0 0.005 0.01 0.02 0.05 0.1 0.2 1];
	lacZ = [0.001 0.003 0.043 0.5 0.984 0.999 1 1];
	output = simf_doi_10_1073_pnas_0606717104_figure4_6_1(x, IPTG) - lacZ;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 1 2 3 4 5 6 7 8 9 10 12 13 17 18 20 25 35 60];
	luc = [0 0.001 0.001 0.001 0.002 0.004 0.007 0.011 0.014 0.032 0.036 0.054 0.089 0.643 0.893 1 1 1 1];
	output = simf_doi__10_1093_nar_25_6_1203_figure4a_1(x, aTc) - luc;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.003 0.005 0.05 0.5 1];
	GFPuv = [0.044 0.089 0.556 1 1];
	output = simf_doi_10_1128_AEM_00791_07_figure3a_1(x, arabinose) - GFPuv;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.003 0.008 0.025 0.075 0.2 0.7 2 6 20 50];
	mCherry = [0.015 0.015 0.015 0.02 0.06 0.48 0.84 0.93 0.96 0.98 0.98 1];
	output = simf_doi_10_1038_nature12148_figure16a_1(x, arabinose) - mCherry;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0 0.002 0.01 0.04 0.16 1 4 25];
	mRFP1 = [0.001 0.001 0.001 0.003 0.045 0.161 0.184 0.185 0.185];
	output = simf_doi_10_1038_nature11516_figureS8A_1(x, arabinose) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0 0 0.001 0.004 0.025 0.1 0.5 2.5];
	mRFP1 = [0.005 0.005 0.005 0.006 0.013 0.09 0.173 0.194 0.196];
	output = simf_doi_10_1038_nature11516_figureS8B_2(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1];
	sfGFP = [0.025 0.025 0.035 0.055 0.115 0.2 0.3 0.4 0.65 1];
	output = simf_doi_10_1038_nbt_2401_figureS11_1(x, IPTG) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.1 1 2 5 7 10 12.5 25 37.5 50 62.5];
	sfGFP = [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425];
	output = simf_doi_10_1038_nbt_2401_figureS13_2(x, arabinose) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.015 0.05 0.5 5 10];
	eYFP = [0.002 0.002 0.057 0.628 0.999 1 1];
	output = simf_doi_10_1038_nature09565_figure1C_1(x, arabinose) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 0.025 0.25 2.5 25 100 250];
	eYFP = [0.005 0.006 0.007 0.03 0.297 0.388 0.398];
	output = simf_doi_10_1038_nature09565_figureS1C_2(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.025 0.05 0.1 0.25 0.5];
	mRFP1 = [0.02 0.031 0.092 0.306 0.745 0.898 0.969 1 1];
	output = simf_doi_10_1038_nbt_2149_Supplement_page17_1(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_1(x) - 0.69;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_3(x) - 0.193;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_4(x) - 0.4;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_5(x) - 0.655;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_6(x) - 0.897;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_7(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_1(x) - 0.422;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_2(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_1(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_2(x) - 0.215;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_3(x) - 0.095;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_4(x) - 0.046;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_5(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_6(x) - 0.003;
	current_base_index = current_base_index + 1;
	%==============
	aTc = [0 2 4 10 25 50 100];
	GFPmut3b = [0 0.002 0.014 0.167 0.666 0.817 0.841];
	output = simf_doi_10_1093_nar_gku1388_figure2B_7(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0 0 0 0.001 0.022 0.176 0.411 0.461];
	output = simf_doi_10_1093_nar_gku1388_figure2B_8(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.005 0.118 0.277];
	output = simf_doi_10_1093_nar_gku1388_figure2B_9(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.003 0.004 0.032];
	output = simf_doi_10_1093_nar_gku1388_figure2B_11(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.006 0.006 0.006 0.006 0.006 0.006 0.009 0.021];
	output = simf_doi_10_1093_nar_gku1388_figure2B_12(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_1_1(x) - 0.002;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_2_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_3_3(x) - 0.036;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_4_4(x) - 0.05;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_5_5(x) - 0.086;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_6_6(x) - 0.214;
	current_base_index = current_base_index + 1;
	%==============
	arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2];
	GFPmut3b = [0 0.002 0.041 0.479 0.949 0.997 1 1 1];
	output = simf_doi_10_1093_nar_gku593_figureS6_7(x, arabinose) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0.001 0.01 0.03 0.05 0.075 0.1 0.15 0.2 0.3 0.5 1];
	eYFP = [0.015 0.014 0.023 0.062 0.098 0.187 0.239 0.374 0.572 0.869 1];
	output = simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1(x, IPTG) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [9.258 13.887 20.831 32.403 46.29 83.322 162.015 231.45 370.32 648.06 925.8 1388.7];
	eYFP = [0.003 0.003 0.003 0.003 0.005 0.01 0.05 0.167 0.333 0.667 0.833 1];
	output = simf_doi_10_1073pnas_0408507102_figure2A_1(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2B_11(x)
	error = zeros(198,1);
	current_base_index = 0;
	%==============
	IPTG = [0 0.005 0.01 0.02 0.05 0.1 0.2 1];
	lacZ = [0.001 0.003 0.043 0.5 0.984 0.999 1 1];
	output = simf_doi_10_1073_pnas_0606717104_figure4_6_1(x, IPTG) - lacZ;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 1 2 3 4 5 6 7 8 9 10 12 13 17 18 20 25 35 60];
	luc = [0 0.001 0.001 0.001 0.002 0.004 0.007 0.011 0.014 0.032 0.036 0.054 0.089 0.643 0.893 1 1 1 1];
	output = simf_doi__10_1093_nar_25_6_1203_figure4a_1(x, aTc) - luc;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.003 0.005 0.05 0.5 1];
	GFPuv = [0.044 0.089 0.556 1 1];
	output = simf_doi_10_1128_AEM_00791_07_figure3a_1(x, arabinose) - GFPuv;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.003 0.008 0.025 0.075 0.2 0.7 2 6 20 50];
	mCherry = [0.015 0.015 0.015 0.02 0.06 0.48 0.84 0.93 0.96 0.98 0.98 1];
	output = simf_doi_10_1038_nature12148_figure16a_1(x, arabinose) - mCherry;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0 0.002 0.01 0.04 0.16 1 4 25];
	mRFP1 = [0.001 0.001 0.001 0.003 0.045 0.161 0.184 0.185 0.185];
	output = simf_doi_10_1038_nature11516_figureS8A_1(x, arabinose) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0 0 0.001 0.004 0.025 0.1 0.5 2.5];
	mRFP1 = [0.005 0.005 0.005 0.006 0.013 0.09 0.173 0.194 0.196];
	output = simf_doi_10_1038_nature11516_figureS8B_2(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1];
	sfGFP = [0.025 0.025 0.035 0.055 0.115 0.2 0.3 0.4 0.65 1];
	output = simf_doi_10_1038_nbt_2401_figureS11_1(x, IPTG) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.1 1 2 5 7 10 12.5 25 37.5 50 62.5];
	sfGFP = [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425];
	output = simf_doi_10_1038_nbt_2401_figureS13_2(x, arabinose) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.015 0.05 0.5 5 10];
	eYFP = [0.002 0.002 0.057 0.628 0.999 1 1];
	output = simf_doi_10_1038_nature09565_figure1C_1(x, arabinose) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 0.025 0.25 2.5 25 100 250];
	eYFP = [0.005 0.006 0.007 0.03 0.297 0.388 0.398];
	output = simf_doi_10_1038_nature09565_figureS1C_2(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.025 0.05 0.1 0.25 0.5];
	mRFP1 = [0.02 0.031 0.092 0.306 0.745 0.898 0.969 1 1];
	output = simf_doi_10_1038_nbt_2149_Supplement_page17_1(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_1(x) - 0.69;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_3(x) - 0.193;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_4(x) - 0.4;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_5(x) - 0.655;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_6(x) - 0.897;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_7(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_1(x) - 0.422;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_2(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_1(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_2(x) - 0.215;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_3(x) - 0.095;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_4(x) - 0.046;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_5(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_6(x) - 0.003;
	current_base_index = current_base_index + 1;
	%==============
	aTc = [0 2 4 10 25 50 100];
	GFPmut3b = [0 0.002 0.014 0.167 0.666 0.817 0.841];
	output = simf_doi_10_1093_nar_gku1388_figure2B_7(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0 0 0 0.001 0.022 0.176 0.411 0.461];
	output = simf_doi_10_1093_nar_gku1388_figure2B_8(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.005 0.118 0.277];
	output = simf_doi_10_1093_nar_gku1388_figure2B_9(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.001 0.001 0.001 0.001 0.001 0.004 0.039 0.147];
	output = simf_doi_10_1093_nar_gku1388_figure2B_10(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.006 0.006 0.006 0.006 0.006 0.006 0.009 0.021];
	output = simf_doi_10_1093_nar_gku1388_figure2B_12(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_1_1(x) - 0.002;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_2_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_3_3(x) - 0.036;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_4_4(x) - 0.05;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_5_5(x) - 0.086;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_6_6(x) - 0.214;
	current_base_index = current_base_index + 1;
	%==============
	arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2];
	GFPmut3b = [0 0.002 0.041 0.479 0.949 0.997 1 1 1];
	output = simf_doi_10_1093_nar_gku593_figureS6_7(x, arabinose) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0.001 0.01 0.03 0.05 0.075 0.1 0.15 0.2 0.3 0.5 1];
	eYFP = [0.015 0.014 0.023 0.062 0.098 0.187 0.239 0.374 0.572 0.869 1];
	output = simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1(x, IPTG) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [9.258 13.887 20.831 32.403 46.29 83.322 162.015 231.45 370.32 648.06 925.8 1388.7];
	eYFP = [0.003 0.003 0.003 0.003 0.005 0.01 0.05 0.167 0.333 0.667 0.833 1];
	output = simf_doi_10_1073pnas_0408507102_figure2A_1(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku1388_figure2B_12(x)
	error = zeros(198,1);
	current_base_index = 0;
	%==============
	IPTG = [0 0.005 0.01 0.02 0.05 0.1 0.2 1];
	lacZ = [0.001 0.003 0.043 0.5 0.984 0.999 1 1];
	output = simf_doi_10_1073_pnas_0606717104_figure4_6_1(x, IPTG) - lacZ;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 1 2 3 4 5 6 7 8 9 10 12 13 17 18 20 25 35 60];
	luc = [0 0.001 0.001 0.001 0.002 0.004 0.007 0.011 0.014 0.032 0.036 0.054 0.089 0.643 0.893 1 1 1 1];
	output = simf_doi__10_1093_nar_25_6_1203_figure4a_1(x, aTc) - luc;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.003 0.005 0.05 0.5 1];
	GFPuv = [0.044 0.089 0.556 1 1];
	output = simf_doi_10_1128_AEM_00791_07_figure3a_1(x, arabinose) - GFPuv;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.003 0.008 0.025 0.075 0.2 0.7 2 6 20 50];
	mCherry = [0.015 0.015 0.015 0.02 0.06 0.48 0.84 0.93 0.96 0.98 0.98 1];
	output = simf_doi_10_1038_nature12148_figure16a_1(x, arabinose) - mCherry;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0 0.002 0.01 0.04 0.16 1 4 25];
	mRFP1 = [0.001 0.001 0.001 0.003 0.045 0.161 0.184 0.185 0.185];
	output = simf_doi_10_1038_nature11516_figureS8A_1(x, arabinose) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0 0 0.001 0.004 0.025 0.1 0.5 2.5];
	mRFP1 = [0.005 0.005 0.005 0.006 0.013 0.09 0.173 0.194 0.196];
	output = simf_doi_10_1038_nature11516_figureS8B_2(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1];
	sfGFP = [0.025 0.025 0.035 0.055 0.115 0.2 0.3 0.4 0.65 1];
	output = simf_doi_10_1038_nbt_2401_figureS11_1(x, IPTG) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.1 1 2 5 7 10 12.5 25 37.5 50 62.5];
	sfGFP = [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425];
	output = simf_doi_10_1038_nbt_2401_figureS13_2(x, arabinose) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.015 0.05 0.5 5 10];
	eYFP = [0.002 0.002 0.057 0.628 0.999 1 1];
	output = simf_doi_10_1038_nature09565_figure1C_1(x, arabinose) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 0.025 0.25 2.5 25 100 250];
	eYFP = [0.005 0.006 0.007 0.03 0.297 0.388 0.398];
	output = simf_doi_10_1038_nature09565_figureS1C_2(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.025 0.05 0.1 0.25 0.5];
	mRFP1 = [0.02 0.031 0.092 0.306 0.745 0.898 0.969 1 1];
	output = simf_doi_10_1038_nbt_2149_Supplement_page17_1(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_1(x) - 0.69;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_3(x) - 0.193;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_4(x) - 0.4;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_5(x) - 0.655;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_6(x) - 0.897;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_7(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_1(x) - 0.422;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_2(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_1(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_2(x) - 0.215;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_3(x) - 0.095;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_4(x) - 0.046;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_5(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_6(x) - 0.003;
	current_base_index = current_base_index + 1;
	%==============
	aTc = [0 2 4 10 25 50 100];
	GFPmut3b = [0 0.002 0.014 0.167 0.666 0.817 0.841];
	output = simf_doi_10_1093_nar_gku1388_figure2B_7(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0 0 0 0.001 0.022 0.176 0.411 0.461];
	output = simf_doi_10_1093_nar_gku1388_figure2B_8(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.005 0.118 0.277];
	output = simf_doi_10_1093_nar_gku1388_figure2B_9(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.001 0.001 0.001 0.001 0.001 0.004 0.039 0.147];
	output = simf_doi_10_1093_nar_gku1388_figure2B_10(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.003 0.004 0.032];
	output = simf_doi_10_1093_nar_gku1388_figure2B_11(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_1_1(x) - 0.002;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_2_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_3_3(x) - 0.036;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_4_4(x) - 0.05;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_5_5(x) - 0.086;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_6_6(x) - 0.214;
	current_base_index = current_base_index + 1;
	%==============
	arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2];
	GFPmut3b = [0 0.002 0.041 0.479 0.949 0.997 1 1 1];
	output = simf_doi_10_1093_nar_gku593_figureS6_7(x, arabinose) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0.001 0.01 0.03 0.05 0.075 0.1 0.15 0.2 0.3 0.5 1];
	eYFP = [0.015 0.014 0.023 0.062 0.098 0.187 0.239 0.374 0.572 0.869 1];
	output = simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1(x, IPTG) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [9.258 13.887 20.831 32.403 46.29 83.322 162.015 231.45 370.32 648.06 925.8 1388.7];
	eYFP = [0.003 0.003 0.003 0.003 0.005 0.01 0.05 0.167 0.333 0.667 0.833 1];
	output = simf_doi_10_1073pnas_0408507102_figure2A_1(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku593_figure3D_2_2(x)
	error = zeros(205,1);
	current_base_index = 0;
	%==============
	IPTG = [0 0.005 0.01 0.02 0.05 0.1 0.2 1];
	lacZ = [0.001 0.003 0.043 0.5 0.984 0.999 1 1];
	output = simf_doi_10_1073_pnas_0606717104_figure4_6_1(x, IPTG) - lacZ;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 1 2 3 4 5 6 7 8 9 10 12 13 17 18 20 25 35 60];
	luc = [0 0.001 0.001 0.001 0.002 0.004 0.007 0.011 0.014 0.032 0.036 0.054 0.089 0.643 0.893 1 1 1 1];
	output = simf_doi__10_1093_nar_25_6_1203_figure4a_1(x, aTc) - luc;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.003 0.005 0.05 0.5 1];
	GFPuv = [0.044 0.089 0.556 1 1];
	output = simf_doi_10_1128_AEM_00791_07_figure3a_1(x, arabinose) - GFPuv;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.003 0.008 0.025 0.075 0.2 0.7 2 6 20 50];
	mCherry = [0.015 0.015 0.015 0.02 0.06 0.48 0.84 0.93 0.96 0.98 0.98 1];
	output = simf_doi_10_1038_nature12148_figure16a_1(x, arabinose) - mCherry;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0 0.002 0.01 0.04 0.16 1 4 25];
	mRFP1 = [0.001 0.001 0.001 0.003 0.045 0.161 0.184 0.185 0.185];
	output = simf_doi_10_1038_nature11516_figureS8A_1(x, arabinose) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0 0 0.001 0.004 0.025 0.1 0.5 2.5];
	mRFP1 = [0.005 0.005 0.005 0.006 0.013 0.09 0.173 0.194 0.196];
	output = simf_doi_10_1038_nature11516_figureS8B_2(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1];
	sfGFP = [0.025 0.025 0.035 0.055 0.115 0.2 0.3 0.4 0.65 1];
	output = simf_doi_10_1038_nbt_2401_figureS11_1(x, IPTG) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.1 1 2 5 7 10 12.5 25 37.5 50 62.5];
	sfGFP = [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425];
	output = simf_doi_10_1038_nbt_2401_figureS13_2(x, arabinose) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.015 0.05 0.5 5 10];
	eYFP = [0.002 0.002 0.057 0.628 0.999 1 1];
	output = simf_doi_10_1038_nature09565_figure1C_1(x, arabinose) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 0.025 0.25 2.5 25 100 250];
	eYFP = [0.005 0.006 0.007 0.03 0.297 0.388 0.398];
	output = simf_doi_10_1038_nature09565_figureS1C_2(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.025 0.05 0.1 0.25 0.5];
	mRFP1 = [0.02 0.031 0.092 0.306 0.745 0.898 0.969 1 1];
	output = simf_doi_10_1038_nbt_2149_Supplement_page17_1(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_1(x) - 0.69;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_3(x) - 0.193;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_4(x) - 0.4;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_5(x) - 0.655;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_6(x) - 0.897;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_7(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_1(x) - 0.422;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_2(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_1(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_2(x) - 0.215;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_3(x) - 0.095;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_4(x) - 0.046;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_5(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_6(x) - 0.003;
	current_base_index = current_base_index + 1;
	%==============
	aTc = [0 2 4 10 25 50 100];
	GFPmut3b = [0 0.002 0.014 0.167 0.666 0.817 0.841];
	output = simf_doi_10_1093_nar_gku1388_figure2B_7(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0 0 0 0.001 0.022 0.176 0.411 0.461];
	output = simf_doi_10_1093_nar_gku1388_figure2B_8(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.005 0.118 0.277];
	output = simf_doi_10_1093_nar_gku1388_figure2B_9(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.001 0.001 0.001 0.001 0.001 0.004 0.039 0.147];
	output = simf_doi_10_1093_nar_gku1388_figure2B_10(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.003 0.004 0.032];
	output = simf_doi_10_1093_nar_gku1388_figure2B_11(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.006 0.006 0.006 0.006 0.006 0.006 0.009 0.021];
	output = simf_doi_10_1093_nar_gku1388_figure2B_12(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_1_1(x) - 0.002;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_3_3(x) - 0.036;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_4_4(x) - 0.05;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_5_5(x) - 0.086;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_6_6(x) - 0.214;
	current_base_index = current_base_index + 1;
	%==============
	arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2];
	GFPmut3b = [0 0.002 0.041 0.479 0.949 0.997 1 1 1];
	output = simf_doi_10_1093_nar_gku593_figureS6_7(x, arabinose) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0.001 0.01 0.03 0.05 0.075 0.1 0.15 0.2 0.3 0.5 1];
	eYFP = [0.015 0.014 0.023 0.062 0.098 0.187 0.239 0.374 0.572 0.869 1];
	output = simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1(x, IPTG) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [9.258 13.887 20.831 32.403 46.29 83.322 162.015 231.45 370.32 648.06 925.8 1388.7];
	eYFP = [0.003 0.003 0.003 0.003 0.005 0.01 0.05 0.167 0.333 0.667 0.833 1];
	output = simf_doi_10_1073pnas_0408507102_figure2A_1(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku593_figure3D_3_3(x)
	error = zeros(205,1);
	current_base_index = 0;
	%==============
	IPTG = [0 0.005 0.01 0.02 0.05 0.1 0.2 1];
	lacZ = [0.001 0.003 0.043 0.5 0.984 0.999 1 1];
	output = simf_doi_10_1073_pnas_0606717104_figure4_6_1(x, IPTG) - lacZ;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 1 2 3 4 5 6 7 8 9 10 12 13 17 18 20 25 35 60];
	luc = [0 0.001 0.001 0.001 0.002 0.004 0.007 0.011 0.014 0.032 0.036 0.054 0.089 0.643 0.893 1 1 1 1];
	output = simf_doi__10_1093_nar_25_6_1203_figure4a_1(x, aTc) - luc;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.003 0.005 0.05 0.5 1];
	GFPuv = [0.044 0.089 0.556 1 1];
	output = simf_doi_10_1128_AEM_00791_07_figure3a_1(x, arabinose) - GFPuv;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.003 0.008 0.025 0.075 0.2 0.7 2 6 20 50];
	mCherry = [0.015 0.015 0.015 0.02 0.06 0.48 0.84 0.93 0.96 0.98 0.98 1];
	output = simf_doi_10_1038_nature12148_figure16a_1(x, arabinose) - mCherry;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0 0.002 0.01 0.04 0.16 1 4 25];
	mRFP1 = [0.001 0.001 0.001 0.003 0.045 0.161 0.184 0.185 0.185];
	output = simf_doi_10_1038_nature11516_figureS8A_1(x, arabinose) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0 0 0.001 0.004 0.025 0.1 0.5 2.5];
	mRFP1 = [0.005 0.005 0.005 0.006 0.013 0.09 0.173 0.194 0.196];
	output = simf_doi_10_1038_nature11516_figureS8B_2(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1];
	sfGFP = [0.025 0.025 0.035 0.055 0.115 0.2 0.3 0.4 0.65 1];
	output = simf_doi_10_1038_nbt_2401_figureS11_1(x, IPTG) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.1 1 2 5 7 10 12.5 25 37.5 50 62.5];
	sfGFP = [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425];
	output = simf_doi_10_1038_nbt_2401_figureS13_2(x, arabinose) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.015 0.05 0.5 5 10];
	eYFP = [0.002 0.002 0.057 0.628 0.999 1 1];
	output = simf_doi_10_1038_nature09565_figure1C_1(x, arabinose) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 0.025 0.25 2.5 25 100 250];
	eYFP = [0.005 0.006 0.007 0.03 0.297 0.388 0.398];
	output = simf_doi_10_1038_nature09565_figureS1C_2(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.025 0.05 0.1 0.25 0.5];
	mRFP1 = [0.02 0.031 0.092 0.306 0.745 0.898 0.969 1 1];
	output = simf_doi_10_1038_nbt_2149_Supplement_page17_1(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_1(x) - 0.69;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_3(x) - 0.193;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_4(x) - 0.4;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_5(x) - 0.655;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_6(x) - 0.897;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_7(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_1(x) - 0.422;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_2(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_1(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_2(x) - 0.215;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_3(x) - 0.095;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_4(x) - 0.046;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_5(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_6(x) - 0.003;
	current_base_index = current_base_index + 1;
	%==============
	aTc = [0 2 4 10 25 50 100];
	GFPmut3b = [0 0.002 0.014 0.167 0.666 0.817 0.841];
	output = simf_doi_10_1093_nar_gku1388_figure2B_7(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0 0 0 0.001 0.022 0.176 0.411 0.461];
	output = simf_doi_10_1093_nar_gku1388_figure2B_8(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.005 0.118 0.277];
	output = simf_doi_10_1093_nar_gku1388_figure2B_9(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.001 0.001 0.001 0.001 0.001 0.004 0.039 0.147];
	output = simf_doi_10_1093_nar_gku1388_figure2B_10(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.003 0.004 0.032];
	output = simf_doi_10_1093_nar_gku1388_figure2B_11(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.006 0.006 0.006 0.006 0.006 0.006 0.009 0.021];
	output = simf_doi_10_1093_nar_gku1388_figure2B_12(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_1_1(x) - 0.002;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_2_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_4_4(x) - 0.05;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_5_5(x) - 0.086;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_6_6(x) - 0.214;
	current_base_index = current_base_index + 1;
	%==============
	arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2];
	GFPmut3b = [0 0.002 0.041 0.479 0.949 0.997 1 1 1];
	output = simf_doi_10_1093_nar_gku593_figureS6_7(x, arabinose) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0.001 0.01 0.03 0.05 0.075 0.1 0.15 0.2 0.3 0.5 1];
	eYFP = [0.015 0.014 0.023 0.062 0.098 0.187 0.239 0.374 0.572 0.869 1];
	output = simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1(x, IPTG) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [9.258 13.887 20.831 32.403 46.29 83.322 162.015 231.45 370.32 648.06 925.8 1388.7];
	eYFP = [0.003 0.003 0.003 0.003 0.005 0.01 0.05 0.167 0.333 0.667 0.833 1];
	output = simf_doi_10_1073pnas_0408507102_figure2A_1(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku593_figure3D_4_4(x)
	error = zeros(205,1);
	current_base_index = 0;
	%==============
	IPTG = [0 0.005 0.01 0.02 0.05 0.1 0.2 1];
	lacZ = [0.001 0.003 0.043 0.5 0.984 0.999 1 1];
	output = simf_doi_10_1073_pnas_0606717104_figure4_6_1(x, IPTG) - lacZ;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 1 2 3 4 5 6 7 8 9 10 12 13 17 18 20 25 35 60];
	luc = [0 0.001 0.001 0.001 0.002 0.004 0.007 0.011 0.014 0.032 0.036 0.054 0.089 0.643 0.893 1 1 1 1];
	output = simf_doi__10_1093_nar_25_6_1203_figure4a_1(x, aTc) - luc;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.003 0.005 0.05 0.5 1];
	GFPuv = [0.044 0.089 0.556 1 1];
	output = simf_doi_10_1128_AEM_00791_07_figure3a_1(x, arabinose) - GFPuv;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.003 0.008 0.025 0.075 0.2 0.7 2 6 20 50];
	mCherry = [0.015 0.015 0.015 0.02 0.06 0.48 0.84 0.93 0.96 0.98 0.98 1];
	output = simf_doi_10_1038_nature12148_figure16a_1(x, arabinose) - mCherry;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0 0.002 0.01 0.04 0.16 1 4 25];
	mRFP1 = [0.001 0.001 0.001 0.003 0.045 0.161 0.184 0.185 0.185];
	output = simf_doi_10_1038_nature11516_figureS8A_1(x, arabinose) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0 0 0.001 0.004 0.025 0.1 0.5 2.5];
	mRFP1 = [0.005 0.005 0.005 0.006 0.013 0.09 0.173 0.194 0.196];
	output = simf_doi_10_1038_nature11516_figureS8B_2(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1];
	sfGFP = [0.025 0.025 0.035 0.055 0.115 0.2 0.3 0.4 0.65 1];
	output = simf_doi_10_1038_nbt_2401_figureS11_1(x, IPTG) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.1 1 2 5 7 10 12.5 25 37.5 50 62.5];
	sfGFP = [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425];
	output = simf_doi_10_1038_nbt_2401_figureS13_2(x, arabinose) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.015 0.05 0.5 5 10];
	eYFP = [0.002 0.002 0.057 0.628 0.999 1 1];
	output = simf_doi_10_1038_nature09565_figure1C_1(x, arabinose) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 0.025 0.25 2.5 25 100 250];
	eYFP = [0.005 0.006 0.007 0.03 0.297 0.388 0.398];
	output = simf_doi_10_1038_nature09565_figureS1C_2(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.025 0.05 0.1 0.25 0.5];
	mRFP1 = [0.02 0.031 0.092 0.306 0.745 0.898 0.969 1 1];
	output = simf_doi_10_1038_nbt_2149_Supplement_page17_1(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_1(x) - 0.69;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_3(x) - 0.193;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_4(x) - 0.4;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_5(x) - 0.655;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_6(x) - 0.897;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_7(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_1(x) - 0.422;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_2(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_1(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_2(x) - 0.215;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_3(x) - 0.095;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_4(x) - 0.046;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_5(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_6(x) - 0.003;
	current_base_index = current_base_index + 1;
	%==============
	aTc = [0 2 4 10 25 50 100];
	GFPmut3b = [0 0.002 0.014 0.167 0.666 0.817 0.841];
	output = simf_doi_10_1093_nar_gku1388_figure2B_7(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0 0 0 0.001 0.022 0.176 0.411 0.461];
	output = simf_doi_10_1093_nar_gku1388_figure2B_8(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.005 0.118 0.277];
	output = simf_doi_10_1093_nar_gku1388_figure2B_9(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.001 0.001 0.001 0.001 0.001 0.004 0.039 0.147];
	output = simf_doi_10_1093_nar_gku1388_figure2B_10(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.003 0.004 0.032];
	output = simf_doi_10_1093_nar_gku1388_figure2B_11(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.006 0.006 0.006 0.006 0.006 0.006 0.009 0.021];
	output = simf_doi_10_1093_nar_gku1388_figure2B_12(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_1_1(x) - 0.002;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_2_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_3_3(x) - 0.036;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_5_5(x) - 0.086;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_6_6(x) - 0.214;
	current_base_index = current_base_index + 1;
	%==============
	arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2];
	GFPmut3b = [0 0.002 0.041 0.479 0.949 0.997 1 1 1];
	output = simf_doi_10_1093_nar_gku593_figureS6_7(x, arabinose) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0.001 0.01 0.03 0.05 0.075 0.1 0.15 0.2 0.3 0.5 1];
	eYFP = [0.015 0.014 0.023 0.062 0.098 0.187 0.239 0.374 0.572 0.869 1];
	output = simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1(x, IPTG) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [9.258 13.887 20.831 32.403 46.29 83.322 162.015 231.45 370.32 648.06 925.8 1388.7];
	eYFP = [0.003 0.003 0.003 0.003 0.005 0.01 0.05 0.167 0.333 0.667 0.833 1];
	output = simf_doi_10_1073pnas_0408507102_figure2A_1(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku593_figure3D_5_5(x)
	error = zeros(205,1);
	current_base_index = 0;
	%==============
	IPTG = [0 0.005 0.01 0.02 0.05 0.1 0.2 1];
	lacZ = [0.001 0.003 0.043 0.5 0.984 0.999 1 1];
	output = simf_doi_10_1073_pnas_0606717104_figure4_6_1(x, IPTG) - lacZ;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 1 2 3 4 5 6 7 8 9 10 12 13 17 18 20 25 35 60];
	luc = [0 0.001 0.001 0.001 0.002 0.004 0.007 0.011 0.014 0.032 0.036 0.054 0.089 0.643 0.893 1 1 1 1];
	output = simf_doi__10_1093_nar_25_6_1203_figure4a_1(x, aTc) - luc;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.003 0.005 0.05 0.5 1];
	GFPuv = [0.044 0.089 0.556 1 1];
	output = simf_doi_10_1128_AEM_00791_07_figure3a_1(x, arabinose) - GFPuv;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.003 0.008 0.025 0.075 0.2 0.7 2 6 20 50];
	mCherry = [0.015 0.015 0.015 0.02 0.06 0.48 0.84 0.93 0.96 0.98 0.98 1];
	output = simf_doi_10_1038_nature12148_figure16a_1(x, arabinose) - mCherry;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0 0.002 0.01 0.04 0.16 1 4 25];
	mRFP1 = [0.001 0.001 0.001 0.003 0.045 0.161 0.184 0.185 0.185];
	output = simf_doi_10_1038_nature11516_figureS8A_1(x, arabinose) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0 0 0.001 0.004 0.025 0.1 0.5 2.5];
	mRFP1 = [0.005 0.005 0.005 0.006 0.013 0.09 0.173 0.194 0.196];
	output = simf_doi_10_1038_nature11516_figureS8B_2(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1];
	sfGFP = [0.025 0.025 0.035 0.055 0.115 0.2 0.3 0.4 0.65 1];
	output = simf_doi_10_1038_nbt_2401_figureS11_1(x, IPTG) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.1 1 2 5 7 10 12.5 25 37.5 50 62.5];
	sfGFP = [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425];
	output = simf_doi_10_1038_nbt_2401_figureS13_2(x, arabinose) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.015 0.05 0.5 5 10];
	eYFP = [0.002 0.002 0.057 0.628 0.999 1 1];
	output = simf_doi_10_1038_nature09565_figure1C_1(x, arabinose) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 0.025 0.25 2.5 25 100 250];
	eYFP = [0.005 0.006 0.007 0.03 0.297 0.388 0.398];
	output = simf_doi_10_1038_nature09565_figureS1C_2(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.025 0.05 0.1 0.25 0.5];
	mRFP1 = [0.02 0.031 0.092 0.306 0.745 0.898 0.969 1 1];
	output = simf_doi_10_1038_nbt_2149_Supplement_page17_1(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_1(x) - 0.69;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_3(x) - 0.193;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_4(x) - 0.4;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_5(x) - 0.655;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_6(x) - 0.897;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_7(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_1(x) - 0.422;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_2(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_1(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_2(x) - 0.215;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_3(x) - 0.095;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_4(x) - 0.046;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_5(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_6(x) - 0.003;
	current_base_index = current_base_index + 1;
	%==============
	aTc = [0 2 4 10 25 50 100];
	GFPmut3b = [0 0.002 0.014 0.167 0.666 0.817 0.841];
	output = simf_doi_10_1093_nar_gku1388_figure2B_7(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0 0 0 0.001 0.022 0.176 0.411 0.461];
	output = simf_doi_10_1093_nar_gku1388_figure2B_8(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.005 0.118 0.277];
	output = simf_doi_10_1093_nar_gku1388_figure2B_9(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.001 0.001 0.001 0.001 0.001 0.004 0.039 0.147];
	output = simf_doi_10_1093_nar_gku1388_figure2B_10(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.003 0.004 0.032];
	output = simf_doi_10_1093_nar_gku1388_figure2B_11(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.006 0.006 0.006 0.006 0.006 0.006 0.009 0.021];
	output = simf_doi_10_1093_nar_gku1388_figure2B_12(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_1_1(x) - 0.002;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_2_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_3_3(x) - 0.036;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_4_4(x) - 0.05;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_6_6(x) - 0.214;
	current_base_index = current_base_index + 1;
	%==============
	arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2];
	GFPmut3b = [0 0.002 0.041 0.479 0.949 0.997 1 1 1];
	output = simf_doi_10_1093_nar_gku593_figureS6_7(x, arabinose) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0.001 0.01 0.03 0.05 0.075 0.1 0.15 0.2 0.3 0.5 1];
	eYFP = [0.015 0.014 0.023 0.062 0.098 0.187 0.239 0.374 0.572 0.869 1];
	output = simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1(x, IPTG) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [9.258 13.887 20.831 32.403 46.29 83.322 162.015 231.45 370.32 648.06 925.8 1388.7];
	eYFP = [0.003 0.003 0.003 0.003 0.005 0.01 0.05 0.167 0.333 0.667 0.833 1];
	output = simf_doi_10_1073pnas_0408507102_figure2A_1(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku593_figure3D_6_6(x)
	error = zeros(205,1);
	current_base_index = 0;
	%==============
	IPTG = [0 0.005 0.01 0.02 0.05 0.1 0.2 1];
	lacZ = [0.001 0.003 0.043 0.5 0.984 0.999 1 1];
	output = simf_doi_10_1073_pnas_0606717104_figure4_6_1(x, IPTG) - lacZ;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 1 2 3 4 5 6 7 8 9 10 12 13 17 18 20 25 35 60];
	luc = [0 0.001 0.001 0.001 0.002 0.004 0.007 0.011 0.014 0.032 0.036 0.054 0.089 0.643 0.893 1 1 1 1];
	output = simf_doi__10_1093_nar_25_6_1203_figure4a_1(x, aTc) - luc;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.003 0.005 0.05 0.5 1];
	GFPuv = [0.044 0.089 0.556 1 1];
	output = simf_doi_10_1128_AEM_00791_07_figure3a_1(x, arabinose) - GFPuv;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.003 0.008 0.025 0.075 0.2 0.7 2 6 20 50];
	mCherry = [0.015 0.015 0.015 0.02 0.06 0.48 0.84 0.93 0.96 0.98 0.98 1];
	output = simf_doi_10_1038_nature12148_figure16a_1(x, arabinose) - mCherry;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0 0.002 0.01 0.04 0.16 1 4 25];
	mRFP1 = [0.001 0.001 0.001 0.003 0.045 0.161 0.184 0.185 0.185];
	output = simf_doi_10_1038_nature11516_figureS8A_1(x, arabinose) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0 0 0.001 0.004 0.025 0.1 0.5 2.5];
	mRFP1 = [0.005 0.005 0.005 0.006 0.013 0.09 0.173 0.194 0.196];
	output = simf_doi_10_1038_nature11516_figureS8B_2(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1];
	sfGFP = [0.025 0.025 0.035 0.055 0.115 0.2 0.3 0.4 0.65 1];
	output = simf_doi_10_1038_nbt_2401_figureS11_1(x, IPTG) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.1 1 2 5 7 10 12.5 25 37.5 50 62.5];
	sfGFP = [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425];
	output = simf_doi_10_1038_nbt_2401_figureS13_2(x, arabinose) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.015 0.05 0.5 5 10];
	eYFP = [0.002 0.002 0.057 0.628 0.999 1 1];
	output = simf_doi_10_1038_nature09565_figure1C_1(x, arabinose) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 0.025 0.25 2.5 25 100 250];
	eYFP = [0.005 0.006 0.007 0.03 0.297 0.388 0.398];
	output = simf_doi_10_1038_nature09565_figureS1C_2(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.025 0.05 0.1 0.25 0.5];
	mRFP1 = [0.02 0.031 0.092 0.306 0.745 0.898 0.969 1 1];
	output = simf_doi_10_1038_nbt_2149_Supplement_page17_1(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_1(x) - 0.69;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_3(x) - 0.193;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_4(x) - 0.4;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_5(x) - 0.655;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_6(x) - 0.897;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_7(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_1(x) - 0.422;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_2(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_1(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_2(x) - 0.215;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_3(x) - 0.095;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_4(x) - 0.046;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_5(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_6(x) - 0.003;
	current_base_index = current_base_index + 1;
	%==============
	aTc = [0 2 4 10 25 50 100];
	GFPmut3b = [0 0.002 0.014 0.167 0.666 0.817 0.841];
	output = simf_doi_10_1093_nar_gku1388_figure2B_7(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0 0 0 0.001 0.022 0.176 0.411 0.461];
	output = simf_doi_10_1093_nar_gku1388_figure2B_8(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.005 0.118 0.277];
	output = simf_doi_10_1093_nar_gku1388_figure2B_9(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.001 0.001 0.001 0.001 0.001 0.004 0.039 0.147];
	output = simf_doi_10_1093_nar_gku1388_figure2B_10(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.003 0.004 0.032];
	output = simf_doi_10_1093_nar_gku1388_figure2B_11(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.006 0.006 0.006 0.006 0.006 0.006 0.009 0.021];
	output = simf_doi_10_1093_nar_gku1388_figure2B_12(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_1_1(x) - 0.002;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_2_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_3_3(x) - 0.036;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_4_4(x) - 0.05;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_5_5(x) - 0.086;
	current_base_index = current_base_index + 1;
	%==============
	arabinose = [0 0.001 0.005 0.021 0.083 0.33 1.33 5.32 21.2];
	GFPmut3b = [0 0.002 0.041 0.479 0.949 0.997 1 1 1];
	output = simf_doi_10_1093_nar_gku593_figureS6_7(x, arabinose) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0.001 0.01 0.03 0.05 0.075 0.1 0.15 0.2 0.3 0.5 1];
	eYFP = [0.015 0.014 0.023 0.062 0.098 0.187 0.239 0.374 0.572 0.869 1];
	output = simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1(x, IPTG) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [9.258 13.887 20.831 32.403 46.29 83.322 162.015 231.45 370.32 648.06 925.8 1388.7];
	eYFP = [0.003 0.003 0.003 0.003 0.005 0.01 0.05 0.167 0.333 0.667 0.833 1];
	output = simf_doi_10_1073pnas_0408507102_figure2A_1(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
end
function error = LOOCV_sim_error_doi_10_1093_nar_gku593_figureS6_7(x)
	error = zeros(197,1);
	current_base_index = 0;
	%==============
	IPTG = [0 0.005 0.01 0.02 0.05 0.1 0.2 1];
	lacZ = [0.001 0.003 0.043 0.5 0.984 0.999 1 1];
	output = simf_doi_10_1073_pnas_0606717104_figure4_6_1(x, IPTG) - lacZ;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 1 2 3 4 5 6 7 8 9 10 12 13 17 18 20 25 35 60];
	luc = [0 0.001 0.001 0.001 0.002 0.004 0.007 0.011 0.014 0.032 0.036 0.054 0.089 0.643 0.893 1 1 1 1];
	output = simf_doi__10_1093_nar_25_6_1203_figure4a_1(x, aTc) - luc;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.003 0.005 0.05 0.5 1];
	GFPuv = [0.044 0.089 0.556 1 1];
	output = simf_doi_10_1128_AEM_00791_07_figure3a_1(x, arabinose) - GFPuv;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.003 0.008 0.025 0.075 0.2 0.7 2 6 20 50];
	mCherry = [0.015 0.015 0.015 0.02 0.06 0.48 0.84 0.93 0.96 0.98 0.98 1];
	output = simf_doi_10_1038_nature12148_figure16a_1(x, arabinose) - mCherry;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0 0.002 0.01 0.04 0.16 1 4 25];
	mRFP1 = [0.001 0.001 0.001 0.003 0.045 0.161 0.184 0.185 0.185];
	output = simf_doi_10_1038_nature11516_figureS8A_1(x, arabinose) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0 0 0.001 0.004 0.025 0.1 0.5 2.5];
	mRFP1 = [0.005 0.005 0.005 0.006 0.013 0.09 0.173 0.194 0.196];
	output = simf_doi_10_1038_nature11516_figureS8B_2(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.07 0.1];
	sfGFP = [0.025 0.025 0.035 0.055 0.115 0.2 0.3 0.4 0.65 1];
	output = simf_doi_10_1038_nbt_2401_figureS11_1(x, IPTG) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0.1 1 2 5 7 10 12.5 25 37.5 50 62.5];
	sfGFP = [0.015 0.02 0.025 0.05 0.075 0.11 0.15 0.25 0.325 0.375 0.425];
	output = simf_doi_10_1038_nbt_2401_figureS13_2(x, arabinose) - sfGFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	arabinose = [0 0.001 0.015 0.05 0.5 5 10];
	eYFP = [0.002 0.002 0.057 0.628 0.999 1 1];
	output = simf_doi_10_1038_nature09565_figure1C_1(x, arabinose) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 0.025 0.25 2.5 25 100 250];
	eYFP = [0.005 0.006 0.007 0.03 0.297 0.388 0.398];
	output = simf_doi_10_1038_nature09565_figureS1C_2(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	IPTG = [0 0.001 0.005 0.01 0.025 0.05 0.1 0.25 0.5];
	mRFP1 = [0.02 0.031 0.092 0.306 0.745 0.898 0.969 1 1];
	output = simf_doi_10_1038_nbt_2149_Supplement_page17_1(x, IPTG) - mRFP1;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_1(x) - 0.69;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_3(x) - 0.193;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_4(x) - 0.4;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_5(x) - 0.655;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_6(x) - 0.897;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1186_1754_1611_3_4_figure3_7(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_1(x) - 0.422;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi__10_1073_pnas_1301301110_sd03_xls_2(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_1(x) - 1;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_2(x) - 0.215;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_3(x) - 0.095;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_4(x) - 0.046;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_5(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku1388_figure2A_6(x) - 0.003;
	current_base_index = current_base_index + 1;
	%==============
	aTc = [0 2 4 10 25 50 100];
	GFPmut3b = [0 0.002 0.014 0.167 0.666 0.817 0.841];
	output = simf_doi_10_1093_nar_gku1388_figure2B_7(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0 0 0 0.001 0.022 0.176 0.411 0.461];
	output = simf_doi_10_1093_nar_gku1388_figure2B_8(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.005 0.118 0.277];
	output = simf_doi_10_1093_nar_gku1388_figure2B_9(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.001 0.001 0.001 0.001 0.001 0.004 0.039 0.147];
	output = simf_doi_10_1093_nar_gku1388_figure2B_10(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.003 0.003 0.003 0.003 0.003 0.003 0.004 0.032];
	output = simf_doi_10_1093_nar_gku1388_figure2B_11(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [0 2 4 10 25 50 100 200];
	GFPmut3b = [0.006 0.006 0.006 0.006 0.006 0.006 0.009 0.021];
	output = simf_doi_10_1093_nar_gku1388_figure2B_12(x, aTc) - GFPmut3b;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_1_1(x) - 0.002;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_2_2(x) - 0.018;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_3_3(x) - 0.036;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_4_4(x) - 0.05;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_5_5(x) - 0.086;
	current_base_index = current_base_index + 1;
	%==============
	error(current_base_index + 1) =  simf_doi_10_1093_nar_gku593_figure3D_6_6(x) - 0.214;
	current_base_index = current_base_index + 1;
	%==============
	IPTG = [0.001 0.01 0.03 0.05 0.075 0.1 0.15 0.2 0.3 0.5 1];
	eYFP = [0.015 0.014 0.023 0.062 0.098 0.187 0.239 0.374 0.572 0.869 1];
	output = simf_http___dspace_mit_edu_handle_1721_1_8228_figure4_6_1(x, IPTG) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
	current_base_index = current_base_index + length(output);
	%==============
	aTc = [9.258 13.887 20.831 32.403 46.29 83.322 162.015 231.45 370.32 648.06 925.8 1388.7];
	eYFP = [0.003 0.003 0.003 0.003 0.005 0.01 0.05 0.167 0.333 0.667 0.833 1];
	output = simf_doi_10_1073pnas_0408507102_figure2A_1(x, aTc) - eYFP;
	for i = 1:length(output)
		error(current_base_index + i) = output(i);
	end
end

