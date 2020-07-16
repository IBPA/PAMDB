function s = LOOCV_randperm(N, x)
	tmp = randperm(N);
	s = zeros(N - 1,1);
	offset = 0;
	for i = 1:N
		if (tmp(i) == x)
			offset = 1;
		else
			s(i-offset) = tmp(i);
		end
	end
end
