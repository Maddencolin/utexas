# Colin Madden Applied Orbital Mechanics Semester Project
clc
clear
format longE

j2 = 0.0010826267;
j3 = -2.5327 * 10^-6; 
radius_earth = 6378.1363; %km
radius_sun = 700000; %km
mu_earth = 398600.4415; %km^3/s^2
mu_sun = 1.327 * 10^11; %km^3/s^2
mu_moon = 3903; %km^3/s^2
area_mass = 1 * 10^-8; %km^2/kg
c_r = 1.5;
c_d = 2.0;
earth_rot_rate = 2 * pi / 86164; %rad/s
JD_UTC = 367 * 2020 - floor((7 * (2020 + floor((3 + 9) / 12)) / 4)) + floor((275 * 3) / 9) + 1 + 1721013.5 + (1 / 24) * (12 + (1 / 60) * (0 + (0 / 60)));
JD_UT1 = JD_UTC;

%% Task 1 LEO
[r0_LEO, v0_LEO] = pos_vel(6763, 0.001, 50, 0, 0, 0, mu_earth);
odeoptions = odeset('RelTol', 1 * 10^-12, 'AbsTol', 1 * 10^-20);
y_vec = [r0_LEO, v0_LEO];
t0 = (0:10 * 60:24 * 60 * 60 * 5)';

[T_2b, r_2b] = ode45(@(T_2b, r_2b)odeProp(T_2b, r_2b, j2, j3, radius_earth, radius_sun, mu_earth, mu_sun, mu_moon, c_r, area_mass, earth_rot_rate, c_d, JD_UT1, 0, 0, 0, 0, 0, 0), t0, y_vec, odeoptions);

[T_2b_j2, r_2b_j2] = ode45(@(T_2b_j2, r_2b_j2)odeProp(T_2b_j2, r_2b_j2, j2, j3, radius_earth, radius_sun, mu_earth, mu_sun, mu_moon, c_r, area_mass, earth_rot_rate, c_d, JD_UT1, 1, 0, 0, 0, 0, 0), t0, y_vec, odeoptions);
delta_2b_j2 = r_2b_j2 - r_2b;
delta_2b_j2 = sqrt(delta_2b_j2(:, 1).^2 + delta_2b_j2(:, 2).^2 + delta_2b_j2(:, 3).^2);

[T_2b_j3, r_2b_j3] = ode45(@(T_2b_j3, r_2b_j3)odeProp(T_2b_j3, r_2b_j3, j2, j3, radius_earth, radius_sun, mu_earth, mu_sun, mu_moon, c_r, area_mass, earth_rot_rate, c_d, JD_UT1, 0, 1, 0, 0, 0, 0), t0, y_vec, odeoptions);
delta_2b_j3 = r_2b_j3 - r_2b;
delta_2b_j3 = sqrt(delta_2b_j3(:, 1).^2 + delta_2b_j3(:, 2).^2 + delta_2b_j3(:, 3).^2);

[T_2b_d, r_2b_d] = ode45(@(T_2b_d, r_2b_d)odeProp(T_2b_d, r_2b_d, j2, j3, radius_earth, radius_sun, mu_earth, mu_sun, mu_moon, c_r, area_mass, earth_rot_rate, c_d, JD_UT1, 0, 0, 1, 0, 0, 0), t0, y_vec, odeoptions);
delta_2b_d = r_2b_d - r_2b;
delta_2b_d = sqrt(delta_2b_d(:, 1).^2 + delta_2b_d(:, 2).^2 + delta_2b_d(:, 3).^2);

[T_2b_sun, r_2b_sun] = ode45(@(T_2b_sun, r_2b_sun)odeProp(T_2b_sun, r_2b_sun, j2, j3, radius_earth, radius_sun, mu_earth, mu_sun, mu_moon, c_r, area_mass, earth_rot_rate, c_d, JD_UT1, 0, 0, 0, 1, 0, 0), t0, y_vec, odeoptions);
delta_2b_sun = r_2b_sun - r_2b;
delta_2b_sun = sqrt(delta_2b_sun(:, 1).^2 + delta_2b_sun(:, 2).^2 + delta_2b_sun(:, 3).^2);

[T_2b_moon, r_2b_moon] = ode45(@(T_2b_moon, r_2b_moon)odeProp(T_2b_moon, r_2b_moon, j2, j3, radius_earth, radius_sun, mu_earth, mu_sun, mu_moon, c_r, area_mass, earth_rot_rate, c_d, JD_UT1, 0, 0, 0, 0, 1, 0), t0, y_vec, odeoptions);
delta_2b_moon = r_2b_moon - r_2b;
delta_2b_moon = sqrt(delta_2b_moon(:, 1).^2 + delta_2b_moon(:, 2).^2 + delta_2b_moon(:, 3).^2);

[T_2b_srp, r_2b_srp] = ode45(@(T_2b_srp, r_2b_srp)odeProp(T_2b_srp, r_2b_srp, j2, j3, radius_earth, radius_sun, mu_earth, mu_sun, mu_moon, c_r, area_mass, earth_rot_rate, c_d, JD_UT1, 0, 0, 0, 0, 0, 1), t0, y_vec, odeoptions);
delta_2b_srp = r_2b_srp - r_2b;
delta_2b_srp = sqrt(delta_2b_srp(:, 1).^2 + delta_2b_srp(:, 2).^2 + delta_2b_srp(:, 3).^2);

% [T_2b_all, r_2b_all] = ode45(@(T_2b_all, r_2b_all)odeProp(T_2b_all, r_2b_all, j2, j3, radius_earth, radius_sun, mu_earth, mu_sun, mu_moon, c_r, area_mass, earth_rot_rate, c_d, JD_UT1, 1, 1, 1, 1, 1, 1), t0, y_vec, odeoptions);
% delta_2b_all = r_2b_all - r_2b;
% delta_2b_all = sqrt(delta_2b_all(:, 1).^2 + delta_2b_all(:, 2).^2 + delta_2b_all(:, 3).^2);

t0 = (0:10:24 * 60 * 5)';
figure(1)
semilogy(t0, delta_2b_j2)
hold on
semilogy(t0, delta_2b_j3)
hold on
semilogy(t0, delta_2b_d)
hold on
semilogy(t0, delta_2b_sun)
hold on
semilogy(t0, delta_2b_moon)
hold on
semilogy(t0, delta_2b_srp)
hold on
% semilogy(t0, delta_2b_all)
ylabel('RSS Error (km)')
xlabel('Time (Minutes)')
legend('2 Body + J2', '2 Body + J3', '2 Body + Drag', '2 Body + Sun', '2 Body + Moon', '2 Body + SRP') % , 'Full Dynamics')
title('\Deltar(t) for LEO')

%% Task 1 MEO
[r0_MEO, v0_MEO] = pos_vel(26560, 0.001, 55, 0, 0, 0, mu_earth);
odeoptions = odeset('RelTol', 1 * 10^-12, 'AbsTol', 1 * 10^-20);
y_vec = [r0_MEO, v0_MEO];
t0 = (0:10 * 60:24 * 60 * 60 * 5)';

[T_2b, r_2b] = ode45(@(T_2b, r_2b)odeProp(T_2b, r_2b, j2, j3, radius_earth, radius_sun, mu_earth, mu_sun, mu_moon, c_r, area_mass, earth_rot_rate, c_d, JD_UT1, 0, 0, 0, 0, 0, 0), t0, y_vec, odeoptions);

[T_2b_j2, r_2b_j2] = ode45(@(T_2b_j2, r_2b_j2)odeProp(T_2b_j2, r_2b_j2, j2, j3, radius_earth, radius_sun, mu_earth, mu_sun, mu_moon, c_r, area_mass, earth_rot_rate, c_d, JD_UT1, 1, 0, 0, 0, 0, 0), t0, y_vec, odeoptions);
delta_2b_j2 = r_2b_j2 - r_2b;
delta_2b_j2 = sqrt(delta_2b_j2(:, 1).^2 + delta_2b_j2(:, 2).^2 + delta_2b_j2(:, 3).^2);

[T_2b_j3, r_2b_j3] = ode45(@(T_2b_j3, r_2b_j3)odeProp(T_2b_j3, r_2b_j3, j2, j3, radius_earth, radius_sun, mu_earth, mu_sun, mu_moon, c_r, area_mass, earth_rot_rate, c_d, JD_UT1, 0, 1, 0, 0, 0, 0), t0, y_vec, odeoptions);
delta_2b_j3 = r_2b_j3 - r_2b;
delta_2b_j3 = sqrt(delta_2b_j3(:, 1).^2 + delta_2b_j3(:, 2).^2 + delta_2b_j3(:, 3).^2);

[T_2b_d, r_2b_d] = ode45(@(T_2b_d, r_2b_d)odeProp(T_2b_d, r_2b_d, j2, j3, radius_earth, radius_sun, mu_earth, mu_sun, mu_moon, c_r, area_mass, earth_rot_rate, c_d, JD_UT1, 0, 0, 1, 0, 0, 0), t0, y_vec, odeoptions);
delta_2b_d = r_2b_d - r_2b;
delta_2b_d = sqrt(delta_2b_d(:, 1).^2 + delta_2b_d(:, 2).^2 + delta_2b_d(:, 3).^2);

[T_2b_sun, r_2b_sun] = ode45(@(T_2b_sun, r_2b_sun)odeProp(T_2b_sun, r_2b_sun, j2, j3, radius_earth, radius_sun, mu_earth, mu_sun, mu_moon, c_r, area_mass, earth_rot_rate, c_d, JD_UT1, 0, 0, 0, 1, 0, 0), t0, y_vec, odeoptions);
delta_2b_sun = r_2b_sun - r_2b;
delta_2b_sun = sqrt(delta_2b_sun(:, 1).^2 + delta_2b_sun(:, 2).^2 + delta_2b_sun(:, 3).^2);

[T_2b_moon, r_2b_moon] = ode45(@(T_2b_moon, r_2b_moon)odeProp(T_2b_moon, r_2b_moon, j2, j3, radius_earth, radius_sun, mu_earth, mu_sun, mu_moon, c_r, area_mass, earth_rot_rate, c_d, JD_UT1, 0, 0, 0, 0, 1, 0), t0, y_vec, odeoptions);
delta_2b_moon = r_2b_moon - r_2b;
delta_2b_moon = sqrt(delta_2b_moon(:, 1).^2 + delta_2b_moon(:, 2).^2 + delta_2b_moon(:, 3).^2);

[T_2b_srp, r_2b_srp] = ode45(@(T_2b_srp, r_2b_srp)odeProp(T_2b_srp, r_2b_srp, j2, j3, radius_earth, radius_sun, mu_earth, mu_sun, mu_moon, c_r, area_mass, earth_rot_rate, c_d, JD_UT1, 0, 0, 0, 0, 0, 1), t0, y_vec, odeoptions);
delta_2b_srp = r_2b_srp - r_2b;
delta_2b_srp = sqrt(delta_2b_srp(:, 1).^2 + delta_2b_srp(:, 2).^2 + delta_2b_srp(:, 3).^2);

% [T_2b_all, r_2b_all] = ode45(@(T_2b_all, r_2b_all)odeProp(T_2b_all, r_2b_all, j2, j3, radius_earth, radius_sun, mu_earth, mu_sun, mu_moon, c_r, area_mass, earth_rot_rate, c_d, JD_UT1, 1, 1, 1, 1, 1, 1), t0, y_vec, odeoptions);
% delta_2b_all = r_2b_all - r_2b;
% delta_2b_all = sqrt(delta_2b_all(:, 1).^2 + delta_2b_all(:, 2).^2 + delta_2b_all(:, 3).^2);

t0 = (0:10:24 * 60 * 5)';
figure(2)
semilogy(t0, delta_2b_j2)
hold on
semilogy(t0, delta_2b_j3)
hold on
semilogy(t0, delta_2b_d)
hold on
semilogy(t0, delta_2b_sun)
hold on
semilogy(t0, delta_2b_moon)
hold on
semilogy(t0, delta_2b_srp)
hold on
% semilogy(t0, delta_2b_all)
ylabel('RSS Error (km)')
xlabel('Time (minutes)')
legend('2 Body + J2', '2 Body + J3', '2 Body + Drag', '2 Body + Sun', '2 Body + Moon', '2 Body + SRP') % , 'Full Dynamics')
title('\Deltar(t) for MEO')

%% Task 1 GEO
[r0_GEO, v0_GEO] = pos_vel(42164, 0.01, 0.5, -120, 0, 0, mu_earth);
odeoptions = odeset('RelTol', 1 * 10^-12, 'AbsTol', 1 * 10^-20);
y_vec = [r0_GEO, v0_GEO];
t0 = (0:10 * 60:24 * 60 * 60 * 5)';

[T_2b, r_2b] = ode45(@(T_2b, r_2b)odeProp(T_2b, r_2b, j2, j3, radius_earth, radius_sun, mu_earth, mu_sun, mu_moon, c_r, area_mass, earth_rot_rate, c_d, JD_UT1, 0, 0, 0, 0, 0, 0), t0, y_vec, odeoptions);

[T_2b_j2, r_2b_j2] = ode45(@(T_2b_j2, r_2b_j2)odeProp(T_2b_j2, r_2b_j2, j2, j3, radius_earth, radius_sun, mu_earth, mu_sun, mu_moon, c_r, area_mass, earth_rot_rate, c_d, JD_UT1, 1, 0, 0, 0, 0, 0), t0, y_vec, odeoptions);
delta_2b_j2 = r_2b_j2 - r_2b;
delta_2b_j2 = sqrt(delta_2b_j2(:, 1).^2 + delta_2b_j2(:, 2).^2 + delta_2b_j2(:, 3).^2);

[T_2b_j3, r_2b_j3] = ode45(@(T_2b_j3, r_2b_j3)odeProp(T_2b_j3, r_2b_j3, j2, j3, radius_earth, radius_sun, mu_earth, mu_sun, mu_moon, c_r, area_mass, earth_rot_rate, c_d, JD_UT1, 0, 1, 0, 0, 0, 0), t0, y_vec, odeoptions);
delta_2b_j3 = r_2b_j3 - r_2b;
delta_2b_j3 = sqrt(delta_2b_j3(:, 1).^2 + delta_2b_j3(:, 2).^2 + delta_2b_j3(:, 3).^2);

[T_2b_d, r_2b_d] = ode45(@(T_2b_d, r_2b_d)odeProp(T_2b_d, r_2b_d, j2, j3, radius_earth, radius_sun, mu_earth, mu_sun, mu_moon, c_r, area_mass, earth_rot_rate, c_d, JD_UT1, 0, 0, 1, 0, 0, 0), t0, y_vec, odeoptions);
delta_2b_d = r_2b_d - r_2b;
delta_2b_d = sqrt(delta_2b_d(:, 1).^2 + delta_2b_d(:, 2).^2 + delta_2b_d(:, 3).^2);

[T_2b_sun, r_2b_sun] = ode45(@(T_2b_sun, r_2b_sun)odeProp(T_2b_sun, r_2b_sun, j2, j3, radius_earth, radius_sun, mu_earth, mu_sun, mu_moon, c_r, area_mass, earth_rot_rate, c_d, JD_UT1, 0, 0, 0, 1, 0, 0), t0, y_vec, odeoptions);
delta_2b_sun = r_2b_sun - r_2b;
delta_2b_sun = sqrt(delta_2b_sun(:, 1).^2 + delta_2b_sun(:, 2).^2 + delta_2b_sun(:, 3).^2);

[T_2b_moon, r_2b_moon] = ode45(@(T_2b_moon, r_2b_moon)odeProp(T_2b_moon, r_2b_moon, j2, j3, radius_earth, radius_sun, mu_earth, mu_sun, mu_moon, c_r, area_mass, earth_rot_rate, c_d, JD_UT1, 0, 0, 0, 0, 1, 0), t0, y_vec, odeoptions);
delta_2b_moon = r_2b_moon - r_2b;
delta_2b_moon = sqrt(delta_2b_moon(:, 1).^2 + delta_2b_moon(:, 2).^2 + delta_2b_moon(:, 3).^2);

[T_2b_srp, r_2b_srp] = ode45(@(T_2b_srp, r_2b_srp)odeProp(T_2b_srp, r_2b_srp, j2, j3, radius_earth, radius_sun, mu_earth, mu_sun, mu_moon, c_r, area_mass, earth_rot_rate, c_d, JD_UT1, 0, 0, 0, 0, 0, 1), t0, y_vec, odeoptions);
delta_2b_srp = r_2b_srp - r_2b;
delta_2b_srp = sqrt(delta_2b_srp(:, 1).^2 + delta_2b_srp(:, 2).^2 + delta_2b_srp(:, 3).^2);

% [T_2b_all, r_2b_all] = ode45(@(T_2b_all, r_2b_all)odeProp(T_2b_all, r_2b_all, j2, j3, radius_earth, radius_sun, mu_earth, mu_sun, mu_moon, c_r, area_mass, earth_rot_rate, c_d, JD_UT1, 1, 1, 1, 1, 1, 1), t0, y_vec, odeoptions);
% delta_2b_all = r_2b_all - r_2b;
% delta_2b_all = sqrt(delta_2b_all(:, 1).^2 + delta_2b_all(:, 2).^2 + delta_2b_all(:, 3).^2);

t0 = (0:10:24 * 60 * 5)';
figure(3)
semilogy(t0, delta_2b_j2)
hold on
semilogy(t0, delta_2b_j3)
hold on
semilogy(t0, delta_2b_d)
hold on
semilogy(t0, delta_2b_sun)
hold on
semilogy(t0, delta_2b_moon)
hold on
semilogy(t0, delta_2b_srp)
hold on
% semilogy(t0, delta_2b_all)
ylabel('RSS Error (km)')
xlabel('Time (minutes)')
legend('2 Body + J2', '2 Body + J3', '2 Body + Drag', '2 Body + Sun', '2 Body + Moon', '2 Body + SRP') % , 'Full Dynamics')
title('\Deltar(t) for GEO')

%% Task 1 Molniya
[r0_Molniya, v0_Molniya] = pos_vel(26000, 0.72, 75, 90, -90, 0, mu_earth);
odeoptions = odeset('RelTol', 1 * 10^-12, 'AbsTol', 1 * 10^-20);
y_vec = [r0_Molniya, v0_Molniya];
t0 = (0:10 * 60:24 * 60 * 60 * 5)';

[T_2b, r_2b] = ode45(@(T_2b, r_2b)odeProp(T_2b, r_2b, j2, j3, radius_earth, radius_sun, mu_earth, mu_sun, mu_moon, c_r, area_mass, earth_rot_rate, c_d, JD_UT1, 0, 0, 0, 0, 0, 0), t0, y_vec, odeoptions);

[T_2b_j2, r_2b_j2] = ode45(@(T_2b_j2, r_2b_j2)odeProp(T_2b_j2, r_2b_j2, j2, j3, radius_earth, radius_sun, mu_earth, mu_sun, mu_moon, c_r, area_mass, earth_rot_rate, c_d, JD_UT1, 1, 0, 0, 0, 0, 0), t0, y_vec, odeoptions);
delta_2b_j2 = r_2b_j2 - r_2b;
delta_2b_j2 = sqrt(delta_2b_j2(:, 1).^2 + delta_2b_j2(:, 2).^2 + delta_2b_j2(:, 3).^2);

[T_2b_j3, r_2b_j3] = ode45(@(T_2b_j3, r_2b_j3)odeProp(T_2b_j3, r_2b_j3, j2, j3, radius_earth, radius_sun, mu_earth, mu_sun, mu_moon, c_r, area_mass, earth_rot_rate, c_d, JD_UT1, 0, 1, 0, 0, 0, 0), t0, y_vec, odeoptions);
delta_2b_j3 = r_2b_j3 - r_2b;
delta_2b_j3 = sqrt(delta_2b_j3(:, 1).^2 + delta_2b_j3(:, 2).^2 + delta_2b_j3(:, 3).^2);

[T_2b_d, r_2b_d] = ode45(@(T_2b_d, r_2b_d)odeProp(T_2b_d, r_2b_d, j2, j3, radius_earth, radius_sun, mu_earth, mu_sun, mu_moon, c_r, area_mass, earth_rot_rate, c_d, JD_UT1, 0, 0, 1, 0, 0, 0), t0, y_vec, odeoptions);
delta_2b_d = r_2b_d - r_2b;
delta_2b_d = sqrt(delta_2b_d(:, 1).^2 + delta_2b_d(:, 2).^2 + delta_2b_d(:, 3).^2);

[T_2b_sun, r_2b_sun] = ode45(@(T_2b_sun, r_2b_sun)odeProp(T_2b_sun, r_2b_sun, j2, j3, radius_earth, radius_sun, mu_earth, mu_sun, mu_moon, c_r, area_mass, earth_rot_rate, c_d, JD_UT1, 0, 0, 0, 1, 0, 0), t0, y_vec, odeoptions);
delta_2b_sun = r_2b_sun - r_2b;
delta_2b_sun = sqrt(delta_2b_sun(:, 1).^2 + delta_2b_sun(:, 2).^2 + delta_2b_sun(:, 3).^2);

[T_2b_moon, r_2b_moon] = ode45(@(T_2b_moon, r_2b_moon)odeProp(T_2b_moon, r_2b_moon, j2, j3, radius_earth, radius_sun, mu_earth, mu_sun, mu_moon, c_r, area_mass, earth_rot_rate, c_d, JD_UT1, 0, 0, 0, 0, 1, 0), t0, y_vec, odeoptions);
delta_2b_moon = r_2b_moon - r_2b;
delta_2b_moon = sqrt(delta_2b_moon(:, 1).^2 + delta_2b_moon(:, 2).^2 + delta_2b_moon(:, 3).^2);

[T_2b_srp, r_2b_srp] = ode45(@(T_2b_srp, r_2b_srp)odeProp(T_2b_srp, r_2b_srp, j2, j3, radius_earth, radius_sun, mu_earth, mu_sun, mu_moon, c_r, area_mass, earth_rot_rate, c_d, JD_UT1, 0, 0, 0, 0, 0, 1), t0, y_vec, odeoptions);
delta_2b_srp = r_2b_srp - r_2b;
delta_2b_srp = sqrt(delta_2b_srp(:, 1).^2 + delta_2b_srp(:, 2).^2 + delta_2b_srp(:, 3).^2);

% [T_2b_all, r_2b_all] = ode45(@(T_2b_all, r_2b_all)odeProp(T_2b_all, r_2b_all, j2, j3, radius_earth, radius_sun, mu_earth, mu_sun, mu_moon, c_r, area_mass, earth_rot_rate, c_d, JD_UT1, 1, 1, 1, 1, 1, 1), t0, y_vec, odeoptions);
% delta_2b_all = r_2b_all - r_2b;
% delta_2b_all = sqrt(delta_2b_all(:, 1).^2 + delta_2b_all(:, 2).^2 + delta_2b_all(:, 3).^2);

t0 = (0:10:24 * 60 * 5)';
figure(4)
semilogy(t0, delta_2b_j2)
hold on
semilogy(t0, delta_2b_j3)
hold on
semilogy(t0, delta_2b_d)
hold on
semilogy(t0, delta_2b_sun)
hold on
semilogy(t0, delta_2b_moon)
hold on
semilogy(t0, delta_2b_srp)
hold on
% semilogy(t0, delta_2b_all)
ylabel('RSS Error (km)')
xlabel('Time (minutes)')
legend('2 Body + J2', '2 Body + J3', '2 Body + Drag', '2 Body + Sun', '2 Body + Moon', '2 Body + SRP') % , 'Full Dynamics')
title('\Deltar(t) for Molniya Orbit')

%% Task 2 NASA data
clc
a = 22081775.58 * 0.3048 / 1000; %km
e = 0.0219118;
nu = 23.969; % deg
phi1 = 28.5; % deg N
lambda1 = 279.45; % deg E
phi2 = 34; % deg N
lambda2 = 241; % deg E
n = 3;
numIterations = 3;
[a, e, i, arg_p, RAAN, nu] = burnout(a, e, nu, phi1, lambda1, phi2, lambda2, n, numIterations, mu_earth, radius_earth, earth_rot_rate, j2);
% T = 2 * pi * sqrt(a^3 / mu_earth); % period (s)
% E_a = 2 * atan(sqrt((1 - e) / (1 + e)) * tand(nu1 / 2)); % Eccentric anomaly at burnout
% t_nu = T / (2 * pi) * (E_a - e * sin(E_a)); % Time since periapsis at burnout
% r1 = a * (1 - e * cos(E_a));
% v1 = sqrt(mu_earth * (2 / r1 - 1 / a));
% p = a * (1 - e^2);
% lambda_2e = lambda2 + n * earth_rot_rate * 180 / pi * T;
% change_time1 = T / 360 * (lambda_2e - lambda1);
% 
% % 1st Iteration
% delta_lambda_12e = (lambda2 + n * earth_rot_rate * 180 / pi * T) + (earth_rot_rate * 180 / pi) * change_time1 - lambda1;
% nu_2e = acosd(sind(phi2) * sind(phi1) + cosd(phi2) * cosd(phi1) * cosd(delta_lambda_12e)) + nu1;
% E_a_2e = 2 * atan(sqrt((1 - e) / (1 + e)) * tand(nu_2e / 2));
% t_nu_2e = T / (2 * pi) * (E_a_2e - e * sin(E_a_2e));
% psi = asind(sind(delta_lambda_12e) * cosd(phi2) / sind(nu_2e - nu1));
% i = acosd(sind(psi) * cosd(phi1));
% arg_p = asind(sind(phi1) / sind(i)) - nu1;
% 
% delta_argp = 3 * sqrt(mu_earth / a^3) * radius_earth^2 * j2 / (4 * p^2) * (4 - 5 * sind(i)^2) * 180 / pi * T * n;
% delta_RAAN = -3 * sqrt(mu_earth / a^3) * radius_earth^2 * j2 / (2 * p^2) * cosd(i) * 180 / pi * T * n;
% delta_phi2 = sind(i) * cosd(arg_p + nu_2e) / sind(nu_2e - nu1) * delta_argp;
% delta_lambda2 = cosd(i) * secd(arg_p + nu_2e)^2 / (1 + cosd(i)^2 *tand(arg_p + nu_2e)^2) * delta_argp + delta_RAAN;
% 
% phi2 = phi2 - delta_phi2;
% lambda2 = lambda2 - delta_lambda2;
% change_time1 = t_nu_2e - t_nu;
% 
% % 2nd Iteration
% delta_lambda_12e = (lambda2 + n * earth_rot_rate * 180 / pi * T) + (earth_rot_rate * 180 / pi) * change_time1 - lambda1;
% nu_2e = acosd(sind(phi2) * sind(phi1) + cosd(phi2) * cosd(phi1) * cosd(delta_lambda_12e)) + nu1;
% E_a_2e = 2 * atan(sqrt((1 - e) / (1 + e)) * tand(nu_2e / 2));
% t_nu_2e = T / (2 * pi) * (E_a_2e - e * sin(E_a_2e));
% psi = asind(sind(delta_lambda_12e) * cosd(phi2) / sind(nu_2e - nu1));
% i = acosd(sind(psi) * cosd(phi1));
% arg_p = asind(sind(phi1) / sind(i)) - nu1;
% 
% delta_argp = 3 * sqrt(mu_earth / a^3) * radius_earth^2 * j2 / (4 * p^2) * (4 - 5 * sind(i)^2) * 180 / pi * T * n;
% delta_RAAN = -3 * sqrt(mu_earth / a^3) * radius_earth^2 * j2 / (2 * p^2) * cosd(i) * 180 / pi * T * n;
% delta_phi2 = sind(i) * cosd(arg_p + nu_2e) / sind(nu_2e - nu1) * delta_argp;
% delta_lambda2 = cosd(i) * secd(arg_p + nu_2e)^2 / (1 + cosd(i)^2 *tand(arg_p + nu_2e)^2) * delta_argp + delta_RAAN;
% 
% phi2 = phi2 - delta_phi2;
% lambda2 = lambda2 - delta_lambda2;
% change_time1 = t_nu_2e - t_nu;
% 
% % 3rd Iteration
% delta_lambda_12e = (lambda2 + n * earth_rot_rate * 180 / pi * T) + (earth_rot_rate * 180 / pi) * change_time1 - lambda1;
% nu_2e = acosd(sind(phi2) * sind(phi1) + cosd(phi2) * cosd(phi1) * cosd(delta_lambda_12e)) + nu1;
% E_a_2e = 2 * atan(sqrt((1 - e) / (1 + e)) * tand(nu_2e / 2));
% t_nu_2e = T / (2 * pi) * (E_a_2e - e * sin(E_a_2e));
% psi = asind(sind(delta_lambda_12e) * cosd(phi2) / sind(nu_2e - nu1));
% i = acosd(sind(psi) * cosd(phi1));
% arg_p = asind(sind(phi1) / sind(i)) - nu1;
% 
% delta_argp = 3 * sqrt(mu_earth / a^3) * radius_earth^2 * j2 / (4 * p^2) * (4 - 5 * sind(i)^2) * 180 / pi * T * n;
% delta_RAAN = -3 * sqrt(mu_earth / a^3) * radius_earth^2 * j2 / (2 * p^2) * cosd(i) * 180 / pi * T * n;
% delta_phi2 = sind(i) * cosd(arg_p + nu_2e) / sind(nu_2e - nu1) * delta_argp;
% delta_lambda2 = cosd(i) * secd(arg_p + nu_2e)^2 / (1 + cosd(i)^2 *tand(arg_p + nu_2e)^2) * delta_argp + delta_RAAN;
% 
% phi2 = phi2 - delta_phi2;
% lambda2 = lambda2 - delta_lambda2;
% change_time1 = t_nu_2e - t_nu;

%% Task 2 Project data
clc
a = 6500;
e = 0.001;
nu = 20; % deg
n = 10.25;
phi1 = 30.3; % deg N
lambda1 = 360 - 120.6; % deg E
phi2 = 30.3; % deg N
lambda2 = 360 - 97.7; % deg E
numIterations = 4;
[a, e, i, arg_p, RAAN, nu] = burnout(a, e, nu, phi1, lambda1, phi2, lambda2, n, numIterations, mu_earth, radius_earth, earth_rot_rate, j2);
[r0_burnout, v0_burnout] = pos_vel(a, e, i, RAAN, arg_p, nu, mu_earth);
period = 2 * pi * sqrt(a^3 / mu_earth);
t0 = 0:n * period;
JD_UTC = 2451545.0;
JD_UT1 = JD_UTC;
a
e
i
arg_p
RAAN
nu
r0_burnout
v0_burnout
%%
y0 = [r0_burnout, v0_burnout];
odeoptions = odeset('RelTol', 1 * 10^-12, 'AbsTol', 1 * 10^-20);
[T_2b_j2_burnout, r_2b_j2_burnout] = ode45(@(T_2b_j2_burnout, r_2b_j2_burnout)odeProp(T_2b_j2_burnout, r_2b_j2_burnout, j2, j3, radius_earth, radius_sun, mu_earth, mu_sun, mu_moon, c_r, area_mass, earth_rot_rate, c_d, JD_UT1, 1, 0, 0, 0, 0, 0), t0, y0, odeoptions);
r_2b_j2_burnout = [r_2b_j2_burnout(:, 1), r_2b_j2_burnout(:, 2), r_2b_j2_burnout(:, 3)];
j = 1;
lambda_burnout = zeros(length(r_2b_j2_burnout), 1);
phi_burnout = zeros(length(r_2b_j2_burnout), 1);
while j < length(r_2b_j2_burnout) + 1
    nu_g = earth_rot_rate * (j - 1);
    Q = [cos(nu_g), sin(nu_g), 0; -sin(nu_g), cos(nu_g), 0; 0, 0, 1];
    r_2b_j2_burnout(j, :) = Q * r_2b_j2_burnout(j, :)'; 
    j = j + 1;
end

lambda_burnout(:, 1) = atan2d(r_2b_j2_burnout(:, 2), r_2b_j2_burnout(:, 1));
phi_burnout(:, 1) = asind(r_2b_j2_burnout(:, 3) ./ sqrt(r_2b_j2_burnout(:, 1).^2 + r_2b_j2_burnout(:, 2).^2 + r_2b_j2_burnout(:, 3).^2));

figure(5)
load earth_coastline.mat
plot(earth_coastline(:,1),earth_coastline(:,2),'k')
hold on
axis equal
xlim([-180,180])
ylim([-90,90])
xlabel('Longitude (degrees)')
ylabel('Latitude (degrees)')
title('Burnout Orbit Calculated Using j2 + 2 Body Dynamics')
plot(lambda_burnout, phi_burnout, '.');

n = 1;
[T_2b_j2_burnout, r_2b_j2_burnout] = ode45(@(T_2b_j2_burnout, r_2b_j2_burnout)odeProp(T_2b_j2_burnout, r_2b_j2_burnout, j2, j3, radius_earth, radius_sun, mu_earth, mu_sun, mu_moon, c_r, area_mass, earth_rot_rate, c_d, JD_UT1, 1, 1, 1, 1, 1, 1), t0, y0, odeoptions);
v_2b_j2_burnout = [r_2b_j2_burnout(:, 4), r_2b_j2_burnout(:, 5), r_2b_j2_burnout(:, 6)];
r_2b_j2_burnout = [r_2b_j2_burnout(:, 1), r_2b_j2_burnout(:, 2), r_2b_j2_burnout(:, 3)];
j = 1;
lambda_burnout = zeros(length(r_2b_j2_burnout), 1);
phi_burnout = zeros(length(r_2b_j2_burnout), 1);
while j < length(r_2b_j2_burnout) + 1
    nu_g = earth_rot_rate * (j - 1);
    Q = [cos(nu_g), sin(nu_g), 0; -sin(nu_g), cos(nu_g), 0; 0, 0, 1];
    r_2b_j2_burnout(j, :) = Q * r_2b_j2_burnout(j, :)'; 
    j = j + 1;
end
lambda_burnout(:, 1) = atan2d(real(r_2b_j2_burnout(:, 2)), real(r_2b_j2_burnout(:, 1)));
phi_burnout(:, 1) = asind(real(r_2b_j2_burnout(:, 3) ./ sqrt(r_2b_j2_burnout(:, 1).^2 + r_2b_j2_burnout(:, 2).^2 + r_2b_j2_burnout(:, 3).^2)));

figure(6)
load earth_coastline.mat
plot(earth_coastline(:,1),earth_coastline(:,2),'k')
hold on
axis equal
xlim([-180,180])
ylim([-90,90])
xlabel('Longitude (degrees)')
ylabel('Latitude (degrees)')
title('Burnout Orbit Calculated Using a Full-Fidelity Model')
plot(lambda_burnout, phi_burnout, '.');

%% Functions

function [a, e, i, arg_p, RAAN, nu] = burnout(a, e, nu, phi1, lambda1, phi2, lambda2, n, numIterations, mu_earth, radius_earth, earth_rot_rate, j2)
    T = 2 * pi * sqrt(a^3 / mu_earth); % period (s)
    E_a = 2 * atan(sqrt((1 - e) / (1 + e)) * tand(nu / 2)); % Eccentric anomaly at burnout
    p = a * (1 - e^2);
    lambda_2e = lambda2 + n * earth_rot_rate * 180 / pi * T;
    change_time1 = T / 360 * (lambda_2e - lambda1);
    t_nu = T / (2 * pi) * (E_a - e * sin(E_a));
    j = 0;
    while j < numIterations
        delta_lambda_12e = (lambda2 + n * earth_rot_rate * 180 / pi * T) + (earth_rot_rate * 180 / pi) * change_time1 - lambda1;
        nu_2e = acosd(sind(phi2) * sind(phi1) + cosd(phi2) * cosd(phi1) * cosd(delta_lambda_12e)) + nu;
        E_a_2e = 2 * atan(sqrt((1 - e) / (1 + e)) * tand(nu_2e / 2));
        t_nu_2e = T / (2 * pi) * (E_a_2e - e * sin(E_a_2e));
        psi = asind(sind(delta_lambda_12e) * cosd(phi2) / sind(nu_2e - nu))
        i = acosd(sind(psi) * cosd(phi1))
        arg_p = asind(sind(phi1) / sind(i)) - nu
        if j == 0
            delta_argp = 3 * sqrt(mu_earth / a^3) * radius_earth^2 * j2 / (4 * p^2) * (4 - 5 * sind(i)^2) * 180 / pi * (n* T + change_time1);
            delta_RAAN = -3 * sqrt(mu_earth / a^3) * radius_earth^2 * j2 / (2 * p^2) * cosd(i) * 180 / pi * (n * T + change_time1);
            delta_phi2 = sind(i) * cosd(arg_p + nu_2e) / cosd(phi2) * delta_argp;
            delta_lambda2 = cosd(i) * secd(arg_p + nu_2e)^2 / (1 + cosd(i)^2 *tand(arg_p + nu_2e)^2) * delta_argp + delta_RAAN;
            phi2 = phi2 - delta_phi2;
            lambda2 = lambda2 - delta_lambda2;
        end
        change_time1 = t_nu_2e - t_nu;
        j = j + 1;
    end
    RAAN = lambda1 - atand(sind(phi1) * tand(psi));
end

function a_j2 = acc_j2(r_vec, j2, radius_earth, mu)
    r = norm(r_vec);
    r_i = r_vec(1);
    r_j = r_vec(2);
    r_k = r_vec(3);
    a_j2 = ((3 * mu * j2 * radius_earth^2) / (2 * r^5)) * [r_i * ((5 * (r_k / r)^2) - 1), r_j * ((5 * (r_k / r)^2) - 1), r_k * ((5 * (r_k / r)^2) - 3)];
end

function a_j3 = acc_j3(r_vec, j3, radius_earth, mu)
    r = norm(r_vec);
    r_i = r_vec(1);
    r_j = r_vec(2);
    r_k = r_vec(3);
    a_j3 = ((-5 * j3 * mu * radius_earth^3) / (2 * r^7)) * [r_i * ((3 * r_k) - (7 * r_k^3 / r^2)), r_j * ((3 * r_k) - (7 * r_k^3 / r^2)), (6 * r_k^2) - (7 * r_k^4 / r^2) - (3 * r^2 / 5)];
end

function a_d = acc_drag(c_d, area_mass, r, v_rel_vec, radius_earth)
    v_rel = norm(v_rel_vec);
    [rho0, h0, H] = getDensityParams(r - radius_earth);
    density = rho0 * 10^9 * exp(-((r - radius_earth) - h0) / H);
    a_d = -0.5 * c_d * area_mass * density * v_rel_vec * v_rel;
end

function a_srp = acc_srp(r_vec, sunPos_vec, c_r, area_mass, radius_earth, radius_sun)
    r_sc_sun_vec = sunPos_vec - r_vec;
    r_sc_sun = norm(r_sc_sun_vec);
    SF = 1367;
    c = 3 * 10^8;
    p_srp = SF / c;
    a = asin(radius_sun / r_sc_sun);
    b = asin(radius_earth / norm(r_vec));
    c = acos(dot(-r_vec, r_sc_sun_vec) / (norm(r_vec) * r_sc_sun));
    if c < abs(a - b)
        gamma = 0;
    elseif a + b <= c
        gamma = 1;
    else
        if abs(a - b) < c && c < a + b
            x = (c^2 + a^2 - b^2) / (2 * c);
            y = sqrt(a^2 - x^2);
            area = b^2 * acos((c - x) / b) + a^2 * acos(x / a) - c * y;
        else
            area = pi * b^2;
        end
        gamma = 1 - area / (pi * a^2);
    end
    a_srp = -p_srp * gamma * c_r * area_mass * r_sc_sun_vec / r_sc_sun * 1000;
end

function [rho0,h0,H] = getDensityParams(altitude)
    %
    % function [rho0,h0,H] = getDensityParams( altitude )
    % ---------------------------------------------------------------------
    % 
    % Description:
    %
    %  Get the density model parameters for the 1976 Standard Atmosphere
    %  exponential model
    % 
    % Inputs:
    %
    %  altitude - Altitude above the Earth's surface in kilometers
    % 
    % Outputs:
    % 
    %  rho0 - Nominal density in kg/m^3
    %  h0   - Base altitude (not radius!) in km
    %  H    - Scale height in km
    %
    % Assumptions/References:
    %
    %  Vallado and McClain, Third Edition, P. 564
    %
    % Dependencies:
    %
    %  None
    %
    % Modification History:
    % 
    %   18jan17     Brandon A. Jones      original version
    %
    % ---------------------------------------------------------------------
    % Copyright University of Texas at Austin, 2017
    %
    if altitude > 1000
        rho0 = 3.019e-15;
        h0 = 1000;
        H = 268;
    elseif altitude > 900
        rho0 = 5.245e-15;
        h0 = 900;
        H = 181.05;
    elseif altitude > 800
        rho0 = 1.170e-14;
        h0 = 800;
        H = 124.64;
    elseif altitude > 700
        rho0 = 3.614e-14;
        h0 = 700;
        H = 88.667;
    elseif altitude > 600
        rho0 = 1.454e-13;
        h0 = 600;
        H = 71.835;
    elseif altitude > 500
        rho0 = 6.967e-13;
        h0 = 500;
        H = 63.822;
    elseif altitude > 450
        rho0 = 1.585e-12;
        h0 = 450;
        H = 60.828;
    elseif altitude > 400
        rho0 = 3.725e-12;
        h0 = 400;
        H = 58.515;
    elseif altitude > 350
        rho0 = 9.518e-12;
        h0 = 350;
        H = 53.298;
    elseif altitude > 300
        rho0 = 2.418e-11;
        h0 = 300;
        H = 53.628;
    elseif altitude > 250
        rho0 = 7.248e-11;
        h0 = 250;
        H = 45.546;
    elseif altitude > 200
        rho0 = 2.789e-10;
        h0 = 200;
        H = 37.105;
    elseif altitude > 180
        rho0 = 5.464e-10;
        h0 = 180;
        H = 29.740;
    elseif altitude > 150
        rho0 = 2.070e-9;
        h0 = 150;
        H = 22.523;
    elseif altitude > 140
        rho0 = 3.845e-9;
        h0 = 140;
        H = 16.149;
    elseif altitude > 130
        rho0 = 8.484e-9;
        h0 = 130;
        H = 12.636;
    elseif altitude > 120
        rho0 = 2.438e-8;
        h0 = 120;
        H = 9.473;
    elseif altitude > 110
        rho0 = 9.661e-8;
        h0 = 110;
        H = 7.263;
    elseif altitude > 100
        rho0 = 5.297e-7;
        h0 = 100;
        H = 5.877;
    elseif altitude > 90
        rho0 = 3.396e-6;
        h0 = 90;
        H = 5.382;
    elseif altitude > 80
        rho0 = 1.905e-5;
        h0 = 80;
        H = 5.799;
    elseif altitude > 70
        rho0 = 8.770e-5;
        h0 = 70;
        H = 6.549;
    elseif altitude > 60
        rho0 = 3.206e-4;
        h0 = 60;
        H = 7.714;
    elseif altitude > 50
        rho0 = 1.057e-3;
        h0 = 50;
        H = 8.382;
    elseif altitude > 40
        rho0 = 3.972e-3;
        h0 = 40;
        H = 7.554;
    elseif altitude > 30
        rho0 = 1.774e-2;
        h0 = 30;
        H = 6.682;
    elseif altitude > 25
        rho0 = 3.899e-2;
        h0 = 25;
        H = 6.349;
    else
        rho0 = 1.225;
        h0 = 0;
        H = 7.249;
    end
end

function sunPosition = Sun(JD_UT1)
    T_UT1 = (JD_UT1 - 2451545) / 36525;
    T_TDB = T_UT1;
    T_TT = T_UT1;
    lambda_M_sun = 280.460 + 36000.771 * T_UT1; %degrees
    M_sun = 357.52772333 + 35999.05034 * T_TDB; %degrees
    M_sun = M_sun * pi / 180; %rad
    lambda_ecl = lambda_M_sun + 1.914666471 * sin(M_sun) + 0.019994643 * sin(2 * M_sun); %degrees
    lambda_M_sun = lambda_M_sun * pi / 180; %rad
    lambda_ecl = lambda_ecl * pi / 180; %rad
    r_sun = 1.000140612 - 0.016708617 * cos(M_sun) - 0.000139589 * cos(2 * M_sun); %AU
    epsilon = 23.439291 - 0.0130042 * T_TDB; %degrees
    epsilon = epsilon * pi / 180; %rad 
    r_sun_vec_MOD = [r_sun .* cos(lambda_ecl), r_sun .* cos(epsilon) .* sin(lambda_ecl), r_sun .* sin(epsilon) .* sin(lambda_ecl)];
    zeta = 2306.2181 * T_TT + 0.30188 * T_TT^2 + 0.017998 * T_TT^3; %arcseconds
    zeta = zeta * (pi / (180 * 3600)); %rad
    theta = 2004.3109 * T_TT - 0.42665 * T_TT^2 - 0.041833 * T_TT^3; %arcseconds
    theta = theta * (pi / (180 * 3600)); %rad
    z = 2306.2181 * T_TT + 1.09468 * T_TT^2 + 0.018203 * T_TT^3; %arcseconds
    z = z * (pi / (180 * 3600)); %rad
    Q_GCRF_MOD = [cos(zeta), sin(zeta), 0; -sin(zeta), cos(zeta), 0; 0, 0, 1] * [cos(-theta), 0, -sin(-theta); 0, 1, 0; sin(-theta), 0, cos(-theta)] * [cos(z), sin(z), 0; -sin(z), cos(z), 0; 0, 0, 1];
    r_sun_vec_GCRF = Q_GCRF_MOD * r_sun_vec_MOD'; % AU
    sunPosition = r_sun_vec_GCRF * 149597870.7; %km
end

function a_2b = acc_2b(y0, mu_earth)
    y0 = y0(1:3);
    r = norm(y0);
    a_2b = -mu_earth / (r^3) * y0(1:3);
end

function a_sun = acc_sun(r_sc_sun_vec, r_sun, mu_sun)
    r_sc_sun = norm(r_sc_sun_vec);
    a_sun = mu_sun * ((r_sc_sun_vec / r_sc_sun^3) - (r_sun' / norm(r_sun)^3));
end

function a_moon = acc_moon(r_sc_moon_vec, r_moon, mu_moon)
    r_sc_moon = norm(r_sc_moon_vec);
    a_moon = mu_moon * ((r_sc_moon_vec / r_sc_moon^3) - (r_moon / norm(r_moon)^3));
end

function moonPosition = Moon(JD_UT1, radius_earth)
    T_UT1 = (JD_UT1 - 2451545) / 36525;
    T_TDB = T_UT1;
    lambda_ecl = 218.32 + 481267.8813 * T_TDB + 6.29 * sind(134.9 + 477198.85 * T_TDB) - 1.27 * sind(259.2 - 413335.38 * T_TDB) + 0.66 * sind(235.7 + 890534.23 * T_TDB) + 0.21 * sind(269.9 + 954397.70 * T_TDB) - 0.19 * sind(357.5 + 35999.05 * T_TDB) - 0.11 * sind(186.6 + 966404.05 * T_TDB);
    %lambda_ecl = lambda_ecl * pi / 180; %rad
    phi_ecl = 5.13 * sind(93.3 + 483202.03 * T_TDB) + 0.28 * sind(228.2 + 960400.87 * T_TDB) - 0.28 * sind(318.3 + 6003.18 * T_TDB) - 0.17 * sind(217.6 - 407332.20 * T_TDB);
    %phi_ecl = phi_ecl * pi / 1809; %rad
    p = 0.9508 + 0.0518 * cosd(134.9 + 477198.85 * T_TDB) + 0.0095 * cosd(259.2 - 413335.38 * T_TDB) + 0.0078 * cosd(235.7 + 890534.23 * T_TDB) + 0.0028 * cosd(269.9 + 954397.70 * T_TDB);
    %p = p * pi / 180; %rad
    epsilon_bar = 23.439291 - 0.0130042 * T_TDB - (1.64 * 10^-7) * T_TDB^2 + (5.04 * 10 ^-7) * T_TDB^3;
    %epsilon_bar = epsilon_bar * pi / 180; %rad
    r_moon = radius_earth / sind(p);
    moonPosition = r_moon * [cosd(phi_ecl) * cosd(lambda_ecl), cosd(epsilon_bar) * cosd(phi_ecl) * sind(lambda_ecl) - sind(epsilon_bar) * sind(phi_ecl), sind(epsilon_bar) * cosd(phi_ecl) * sind(lambda_ecl) + cosd(epsilon_bar) * sind(phi_ecl)];
    delta_alpha_0 = 0.0146 / 3600; %deg
    zeta_0 = -0.16617 / 3600; %deg
    eta_0 = -0.0068192 / 3600; %deg
    Q_GCRF_J200 = [cosd(-delta_alpha_0), sind(-delta_alpha_0), 0; -sind(-delta_alpha_0), cosd(-delta_alpha_0), 0; 0, 0, 1] * [cosd(-zeta_0), 0, -sind(-zeta_0); 0, 1, 0; sind(-zeta_0), 0, cosd(-zeta_0)] * [1, 0, 0; 0, cosd(eta_0), sind(eta_0); 0, -sind(eta_0), cosd(eta_0)];
    moonPosition = (Q_GCRF_J200 * moonPosition')';
end

function [r, v] = pos_vel(a, e, i, RAAN, arg_p, true_anomaly, mu)
    i = i * pi / 180; % rad
    RAAN = RAAN * pi / 180; % rad
    arg_p = arg_p * pi / 180; % rad
    true_anomaly = true_anomaly * pi / 180; % rad
    p = a * (1 - e^2);
    r_pqw = p / (1 + e * cos(true_anomaly)) * [cos(true_anomaly), sin(true_anomaly), 0];
    v_pqw = sqrt(mu / p) * [-sin(true_anomaly), e + cos(true_anomaly), 0];
    Q_ijk_pqw = [cos(-RAAN), sin(-RAAN), 0; -sin(-RAAN), cos(-RAAN), 0; 0, 0, 1] * [1, 0, 0; 0, cos(-i), sin(-i); 0, -sin(-i), cos(-i)] * [cos(-arg_p), sin(-arg_p), 0; -sin(-arg_p), cos(-arg_p), 0; 0, 0, 1];
    r = Q_ijk_pqw * r_pqw';
    v = Q_ijk_pqw * v_pqw';
end

function y = odeProp(t0, y0, j2, j3, radius_earth, radius_sun, mu_earth, mu_sun, mu_moon, c_r, area_mass, earth_rot_rate, c_d, JD_UT1, x_j2, x_j3, x_drag, x_sun, x_moon, x_srp)
    t0 = JD_UT1 + (t0 / 24 / 60 / 60);
    r = sqrt(y0(1)^2 + y0(2)^2 + y0(3)^2);
    y = zeros(6, 1);
    r_sun = Sun(t0)';
    r_moon = Moon(t0, radius_earth)';
    r_sc_sun_vec = r_sun' - y0(1:3);
    r_sc_moon_vec = r_moon - y0(1:3);
    if x_j2 > 0
        a_j2 = acc_j2([y0(1), y0(2), y0(3)], j2, radius_earth, mu_earth);
    else 
        a_j2 = 0;
    end
    if x_j3 > 0
        a_j3 = acc_j3([y0(1), y0(2), y0(3)], j3, radius_earth, mu_earth);
    else
        a_j3 = 0;
    end
    if x_srp > 0
        a_srp = acc_srp([y0(1), y0(2), y0(3)], r_sun, c_r, area_mass, radius_earth, radius_sun);
    else 
        a_srp = 0;
    end
    if x_drag > 0
        v_rel_vec = [y0(4), y0(5), y0(6)] - cross([0, 0, earth_rot_rate], [y0(1), y0(2), y0(3)]);
        a_d = acc_drag(c_d, area_mass, r, v_rel_vec, radius_earth);
    else
        a_d = 0;
    end
    a_2b = acc_2b(y0, mu_earth);
    if x_sun > 0
        a_sun = acc_sun(r_sc_sun_vec, r_sun, mu_sun);
    else 
        a_sun = 0;
    end
    if x_moon > 0
        a_moon = acc_moon(r_sc_moon_vec, r_moon, mu_moon);
    else
        a_moon = 0;
    end
    y(1) = y0(4);
    y(2) = y0(5);
    y(3) = y0(6);
    y(4:6) = a_2b + a_sun + a_moon + a_j2' + a_j3' + a_srp' + a_d';
end

