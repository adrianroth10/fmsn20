load swissRainfall.mat

normplot(swissRain(:, 1))
print('proj1/normplot_data', '-deps')
close
normplot(sqrt(swissRain(:, 1)))
xlabel('sqrt(Data)')
print('proj1/normplot_sqrt_data', '-deps')
close
normplot(log(swissRain(:, 1) + 1))
xlabel('log(Data + 1)')
print('proj1/normplot_log_data', '-deps')
close
