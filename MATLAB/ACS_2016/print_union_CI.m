function print_union_CI(parameter_name_list, CI_combination, file_name)
	parameter_number = length(parameter_name_list);
	final_CI = zeros(parameter_number,2);
	for i = 1:size(CI_combination,1)
		for j = 1:parameter_number
			if (CI_combination(i,j*3-2) > 0 && CI_combination(i,j*3-1) > 0)
				if (final_CI(j,1) > CI_combination(i,j*3-2) || final_CI(j,1) == 0)
					final_CI(j,1) = CI_combination(i,j*3-2);
				end
				if (final_CI(j,2) < CI_combination(i,j*3-1) || final_CI(j,2) == 0)
					final_CI(j,2) = CI_combination(i,j*3-1);
				end
			end
		end
	end
	CI_fileID = fopen(file_name,'w');
	for i = 1:parameter_number
		fprintf(CI_fileID, '%s \t', parameter_name_list{i});
		fprintf(CI_fileID, '%f \t', final_CI(i,1));
		fprintf(CI_fileID, '%f \n', final_CI(i,2));		
	end
end
