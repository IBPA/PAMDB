function print_simulated_data_to_json(json_file, uri, exp_ref, input, output, model_code)
	fprintf(json_file, '{"uri": "%s", "experiment": "%s", "model": [', uri, exp_ref);
	for i = 1:length(model_code)
		if (i < length(model_code))
			fprintf(json_file, '"%s", ', model_code{i});
		else
			fprintf(json_file, '"%s"], "input_output": [', model_code{i});
		end
	end
	if (length(input) > 0) 
		for i = 1:length(input)
			fprintf(json_file, '[%f, %f]', input(i), output(i));
			if (i < length(input))
				fprintf(json_file, ', ');
			else
				fprintf(json_file, ']}\n');
			end
		end
	else
		fprintf(json_file, '[%f]]}\n', output(1));
	end
end
