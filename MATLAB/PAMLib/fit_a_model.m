function [CI_lb, CI_ub, x_opt, solution_ensemble] = fit_a_model(ensemble_size, ite_num, func_num, error_function, lb, ub)
	options = optimset('MaxIter', ite_num, 'MaxFunEvals', func_num);
	parameter_num = length(lb);
	solution_ensemble = zeros(ensemble_size, parameter_num + 1);
	x0 = zeros(1, parameter_num);
	for i = 1:ensemble_size
		for j = 1:parameter_num
			x0(j) = lb(j) + (ub(j) - lb(j))*rand;
		end
		[x_min, error_min] = lsqnonlin(error_function, x0, lb, ub, options);
		for j = 1:parameter_num
			solution_ensemble(i,j) = x_min(1,j);
		end
		solution_ensemble(i, parameter_num + 1) = error_min;
	end
	n_row = size(solution_ensemble,1);
	n_col = size(solution_ensemble,2);
	min_error = 1e6;
	CI_lb = zeros(1,n_col);
	CI_ub = zeros(1,n_col);
	x_opt = zeros(1,n_col);
	for i = 1:n_row
		if (solution_ensemble(i,n_col) < min_error)
			min_error = solution_ensemble(i,n_col);
			for j = 1:(n_col - 1)
				CI_lb(j) = solution_ensemble(i,j);
				CI_ub(j) = solution_ensemble(i,j);
				x_opt(j) = solution_ensemble(i,j);
			end
		end
	end
	for i = 1:n_row
		if (solution_ensemble(i,n_col) < 1.5*min_error)
			for j = 1:(n_col - 1)
				if (solution_ensemble(i,j) < CI_lb(j))
					CI_lb(j) = solution_ensemble(i,j);
				end
				if (solution_ensemble(i,j) > CI_ub(j))
					CI_ub(j) = solution_ensemble(i,j);
				end
			end			
		end
	end
end
