function print_a_matrix(file_name, M)
	point_fileID = fopen(file_name,'w');
	for i = 1:size(M,1)
		for j = 1:size(M,2)
			fprintf(point_fileID, '%f \t', M(i,j));
		end
		fprintf(point_fileID, '\n');
	end
	fclose(point_fileID);
end
