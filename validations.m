function rms_error = validations(Y_valid, Ezy_valid, Vzy_valid)

dzy_valid = sqrt(Vzy_valid);
n = length(Y_valid);
for alpha = linspace(0.05, 0.5, 10)
  lambda = norminv(1 - alpha / 2);
  dalpha = lambda * dzy_valid;
  lim_min = Ezy_valid - dalpha;
  lim_max = Ezy_valid + dalpha;
  in = lim_min < Y_valid & Y_valid < lim_max;
  p_in = sum(in) / n;
  fprintf(1, 'p_in: %11.4e,  in: %11.4e\n', p_in, 1 - alpha);
end

rms_error = sqrt(mean(Vzy_valid .* (Ezy_valid - Y_valid).^2));

end
