clear 
clc
parameters
omega = 314;

syms Ipoc I2 Im Ib Vb Ppoc Sb Pb

Cc_prime = n^2*Cc;
Rc_prime = n^2*Rc;
Lc_prime = n^2*Lc;
Z1 = R1 + 1j*omega*L1;
Z2_prime = (R2_prime + Rc_prime) + 1j*omega*(L2_prime + Lc_prime);
Zm =  1j*omega*Lm;
Zc = 1/(1j*omega*Cc_prime);

Vpoc = 33000/sqrt(3)*n;
Ppoc = Ipoc*Vpoc;
Ic = Vpoc/Zc;
I2 = Ipoc + Ic;
Vab = Vpoc + I2*Z2_prime;
Im = Vab/Zm;
Ib = Im + I2;
Vb = Vab + Z1*Ib;
Sb = Vb*conj(Ib);
Pb = real(Sb);
current = [0:10:1000];
pb_vect = eval(subs(simplify(Pb), Ipoc, current));
ppoc_vect = eval(subs(simplify(Ppoc), Ipoc, current));

figure()
plot(current, pb_vect)
hold on
plot(current, ppoc_vect)
hold off
%vpasolve(Pb == Ppoc, Ipoc)

% Ipoc_min = -500;
% Ipoc_max = 500;
% for i = 1:1000
% Ipoc = (Ipoc_max + Ipoc_min)*0.5;
% Vpoc = 33000/sqrt(3)*n;
% Ppoc = Ipoc*Vpoc;
% Ic = Vpoc/Zc;
% I2 = Ipoc + Ic;
% Vab = Vpoc + I2*Z2_prime;
% Im = Vab/Zm;
% Ib = Im + I2;
% Vb = Vab + Z1*Ib;
% Sb = Vb*conj(Ib);
% Pb = real(Sb);
% 
% if abs(Pb - Ppoc) < 0.001 && (Ipoc_max - Ipoc_min)<1e-10 
%   disp('convergence')
%   break
% else 
%   if (Pb - Ppoc) > 0
%     Ipoc_max = Ipoc;
%   else
%     Ipoc_min = Ipoc;
%   end
% end
% % if abs(Pb - Ppoc) < 0.1
% %   disp('Ok convergence')
% % else
% %   if Ppoc > Pb
% %     Ipoc = Ipoc - (Ipoc - Ipoc_old)/2;
% %   else
% %     Ipoc = Ipoc + (Ipoc - Ipoc_old)/2;
% %   end
% % end
% % Ipoc_old = Ipoc;
% disp(strcat('Ppoc = ', num2str(Ppoc), 'Pb=', num2str(Pb) , 'diff', num2str(Pb-Ppoc)))
% end
