function y = remove_a_sub_array(x, p1, p2)
	y = zeros(1, length(x) - (p2 - p1 + 1));
	index = 1;
	for i = 1: length(x)
		if (i < p1 || i > p2)
		y(index) = x(i);
		index = index + 1;
	end
end
