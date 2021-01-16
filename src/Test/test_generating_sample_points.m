% test generating sample points
radius = 30; % mm
transducer_element = 512;
transducer_amplitude = 70*pi/180;
transducer_element_separation = transducer_amplitude * radius / transducer_element;
amp = transducer_element_separation / radius;
angle_center_of_element = amp/2;
angle = - (amp * transducer_element /2) + angle_center_of_element;
position1 = [];
position2 = [];
position3 = [];
center_point = [0, 30, 0];

for i = 1 : 512
    delta_x = random('Normal',0,0.005);
    delta_z = random('Normal',0,0.005);
    center_point_noise = [delta_x, 30, delta_z];
%     position1 = [position1; radius * sin(angle), radius * cos(angle), 0];
    temp = rotz(-angle)*center_point';
    position2 = [position2; temp'];
    temp_1 = rotz(-angle)*center_point_noise';
    position3 = [position3; temp_1'];
    angle = angle + amp;
end
figure
scatter3(position2(:,1), position2(:,2), position2(:,3))
hold on
scatter3(position3(:,1), position3(:,2), position3(:,3))





