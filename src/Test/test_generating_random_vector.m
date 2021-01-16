v = [0,0,1];
figure
quiver3(0,0,0, v(1), v(2), v(3),  'b', 'LineWidth', 2, 'ShowArrowHead', 'on','AutoScaleFactor',1 )
xlim([-1,1])
ylim([-1,1])
zlim([0,1])
random_vector1 = [];
roughness = 100000
title ('Roughness(Cos) = 100000')
for i = 1: 100
    cos_theta = power_cosine_variate(roughness);
    random_vector = random_unit_vector(v, cos_theta);
    random_vector1 = [random_vector1; random_vector];
    hold on
    quiver3(0,0,0, random_vector(1), random_vector(2), random_vector(3), 'r','LineWidth', 0.1, 'AutoScaleFactor',1)
end


function cos_theta = power_cosine_variate(x)
    number = random('Uniform', 0, 1);
    indice = x + 1;
    exponente = 1.0 / indice;
    cos_theta = number.^exponente;
end

function random_vector = random_unit_vector(v, cos_theta)
    flag = false;
    p = 0.5;
    px = 0;
    py = 0;
    while (p > 0.25)
        a = random('Uniform', 0, 1)*2*pi;
        r = 0.5 * sqrt(random('Uniform', 0, 1));
        px = r * cos(a);
        py = r * sin(a);
        p = px*px + py*py;
    end
    if abs(v(1)) > abs(v(2))
        temp = v(1);
        v(1) = v(2);
        v(2) = temp;
        flag = true;
    end
    
    b = 1 - v(1)*v(1);
    radicando = 1 - cos_theta*cos_theta;
    radicando = radicando / (p * b);
    c = sqrt(radicando);
    px = px *c;
    py = py *c;
    d = cos_theta - v(1) * px;
    wx = v(1)*cos_theta - b * px;
    wy = v(2)*d + v(3) * py;
    wz = v(3)*d - v(2) * py;
    if flag
       temp1 = wy;
       wy = wx;
       wx = temp1;
    end
    random_vector = [wx, wy, wz]; 
end


