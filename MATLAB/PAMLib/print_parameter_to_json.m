function print_parameter_to_json(filename, lb, ub, CI_lb, CI_ub, x_opt, solution_ensemble, parameter_name_list)
	json_fileID = fopen(filename,'w');
	n_row = size(solution_ensemble,1);
	for k = 1:length(parameter_name_list)
		if (k > 1)
			fprintf(json_fileID, '\n');
		end
		fprintf(json_fileID, '{"name": "%s", ', parameter_name_list{k});
		fprintf(json_fileID, '"lb": %f, ', lb(k));
		fprintf(json_fileID, '"ub": %f, ', ub(k));
		fprintf(json_fileID, '"CI_lb": %f, ', CI_lb(k));
		fprintf(json_fileID, '"CI_ub": %f, ', CI_ub(k));
		fprintf(json_fileID, '"opt": %f, ', x_opt(k));
		fprintf(json_fileID, '"ensemble": [');
		for i = 1:n_row - 1
			fprintf(json_fileID, '%f, ', solution_ensemble(i,k));
		end
		fprintf(json_fileID, '%f]}', solution_ensemble(n_row,k));
	end
	fclose(json_fileID);	
end
