% test point spread function
axial_size = 7;
lateral_size = 13;
elevational_size = 7;
transducer_freq = 4.5; %MHz
resolution_axial = 0.145; %mm
var_axial = 0.145;
var_elevational = 0.3;
var_lateral = 0.3;
half_axial = axial_size * resolution_axial/2;
half_lateral = lateral_size * resolution_axial/2;
half_elevational = elevational_size * resolution_axial/2
axial_kernel = [];
lateral_kernel = [];
elevational_kernel = [];
for i = 1:axial_size
    x = (i-1)*resolution_axial - half_axial;
    axial_kernel = [axial_kernel exp(-0.5 * x*x/var_axial)* cos(2*pi*transducer_freq*x)];
end

for i = 1:lateral_size
    y = (i-1)*resolution_axial - half_lateral;
    lateral_kernel = [lateral_kernel exp(-0.5 * y * y/var_lateral)];
end

for i = 1: elevational_size
    z = (i-1)*resolution_axial - half_elevational;
    elevational_kernel = [elevational_kernel exp(-0.5 * z * z/var_elevational)];
end

convolution1 = axial_kernel'*lateral_kernel;
mesh(convolution1)
% [X,Y] = meshgrid(-3:1:3);
% mesh(X,Y,convolution1(1:7, 1:7))