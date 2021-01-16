% Test transducer
angle_radian = -pi*30/180;
radius = 30; %cm
tranducer_elements = 512
tranducer_element_separation = (60*pi/180)*radius/tranducer_elements
[radius*sin(angle_radian), radius*cos(angle_radian), 0]

p = [1;1;1];
position = [-135;0;0];
rotz(0)
radius* rotz(-90*pi/180)*p

radius*rotx(0*pi/180)*p

radius*roty(0*pi/180)*p

%theta = (-30/512)*pi/180;
theta = (-30 + 30*2/512)*pi/180;
pos = position + radius*roty(0*pi/180)*rotx(0*pi/180)*rotz(-90*pi/180)*[sin(theta); cos(theta); 0]

amplitude = tranducer_element_separation/radius

angle = -(amplitude * tranducer_elements / 2) + amplitude/2
%%





