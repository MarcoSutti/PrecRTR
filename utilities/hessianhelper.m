function Hess = hessianhelper(hessian, udot, norm_xdot)

% 2*hessian(hessian(udot))/norm_xdot;

H = hessian(hessian(udot));

Hess.M = 2/norm_xdot * H.M;
Hess.Up = 2/norm_xdot * H.Up;
Hess.Vp = 2/norm_xdot * H.Vp;

end