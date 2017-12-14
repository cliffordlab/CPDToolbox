function googleColors = loadGoogleColors()

% Load Google Material Design colors into a structure.
% These are much nicer than Matlab default colors.
% RGB values should be in % or fraction out of 100
%
% Uses hex2rgb.m
%
% By Erik Reinertsen <er@gatech.edu>

googleColors.blue = hex2rgb('2196F3');
googleColors.red = hex2rgb('F44336');
googleColors.green = hex2rgb('4CAF50');
googleColors.purple = hex2rgb('9C27B0');
googleColors.orange = hex2rgb('FF9800');
googleColors.gray = hex2rgb('9E9E9E');
googleColors.grey = hex2rgb('9E9E9E');
googleColors.yellow = hex2rgb('FFEB3B');
googleColors.amber = hex2rgb('FFC107');
googleColors.brown = hex2rgb('795548');

end