s = 1:1000;
e = ones(1, 1000)*0.3;

y = power(e, 1./s);

figure
plot(s, y)

%%
R = 0.5;
figure(4)
for i = 1 : 1000

a = random('Uniform', 0, 1)* 2* pi;
r = R * sqrt(random('Uniform', 0, 1));

x = r * cos(a);
y = r * sin(a);
if (x*x+ y*y > 0.25)
    scatter(x, y)
    hold on
end
end
xlim([-2 2])
ylim([-2 2])
