restart;
read cat(kernelopts(homedir),"/opt/FD/FD.mpl"):  # change to your FD path

Clean_FD();
Make_FD();

compactificationX   := d*exp(c*x) - d*exp(-c*x);
compactificationX_i := solve(X = compactificationX, x)[1];
compactificationY   := d*exp(c*y) - d*exp(-c*y);
compactificationY_i := solve(Y = compactificationY, y)[1];
compactificationZ   := d*exp(c*z) - d*exp(-c*z);
compactificationZ_i := solve(Z = compactificationZ, z)[1];

grid_functions := {phi1,pi1,phi2,pi2,At,Bt,Ax,Bx,Ay,By,Az,Bz,chi,xi,modphi,qden,eden,jdenx,lorenz,elecx,elecy,elecz,magx,magy,magz,eem,gaussE,gaussB,dxdX,dxdX2,dydY,dydY2,dzdZ,dzdZ2,Xx,Yy,Zz,jac};

eqn_list := [eqn_phi1, eqn_pi1_poly, eqn_pi1_log, eqn_phi2, eqn_pi2_poly, eqn_pi2_log, eqn_At, eqn_Bt, eqn_Ax, eqn_Bx, eqn_Ay, eqn_By, eqn_Az, eqn_Bz, eqn_chi, eqn_xi, eqn_chi_init, eqn_xi_init, eqn_modphi, eqn_qden, eqn_eden_poly, eqn_eden_log, eqn_jdenx, eqn_lorenz, eqn_elecx, eqn_elecy, eqn_elecz, eqn_magx, eqn_magy, eqn_magz, eqn_eem, eqn_gaussE, eqn_gaussB, eqn_phi1_ire_poly, eqn_phi1_ire_log, eqn_phi2_ire_poly, eqn_phi2_ire_log, eqn_At_ire, eqn_Ax_ire, eqn_Ay_ire, eqn_Az_ire];

Vpoly := m^2*(phi1(t,X,Y,Z)^2 + phi2(t,X,Y,Z)^2) - (g/2)*(phi1(t,X,Y,Z)^2 + phi2(t,X,Y,Z)^2)^2 + (h/3)*(phi1(t,X,Y,Z)^2 + phi2(t,X,Y,Z)^2)^3;
Vlog  := -mu^2*(phi1(t,X,Y,Z)^2 + phi2(t,X,Y,Z)^2)*log(beta^2*(phi1(t,X,Y,Z)^2 + phi2(t,X,Y,Z)^2));
Vpoly_diff1 := Physics :- diff(Vpoly,phi1(t,X,Y,Z));
Vpoly_diff2 := Physics :- diff(Vpoly,phi2(t,X,Y,Z));
Vlog_diff1  := Physics :- diff(Vlog,phi1(t,X,Y,Z));
Vlog_diff2  := Physics :- diff(Vlog,phi2(t,X,Y,Z));

# redefine these quantities to deal with underflow issues in log argument
eps := 10^(-100);
Vlog  := -mu^2*(phi1(t,X,Y,Z)^2 + phi2(t,X,Y,Z)^2)*log(beta^2*(phi1(t,X,Y,Z)^2 + phi2(t,X,Y,Z)^2) + eps);
Vlog_diff1  := -2*mu^2*phi1(t, X, Y, Z)*ln(beta^2*(phi1(t, X, Y, Z)^2 + phi2(t, X, Y, Z)^2) + eps) - 2*mu^2*phi1(t, X, Y, Z);
Vlog_diff2  := -2*mu^2*phi2(t, X, Y, Z)*ln(beta^2*(phi1(t, X, Y, Z)^2 + phi2(t, X, Y, Z)^2) + eps) - 2*mu^2*phi2(t, X, Y, Z);

eqn_phi1 := diff(phi1(t,X,Y,Z),t) = pi1(t,X,Y,Z);
eqn_phi2 := diff(phi2(t,X,Y,Z),t) = pi2(t,X,Y,Z);
eqn_At   := diff(At(t,X,Y,Z),t)   = Bt(t,X,Y,Z);
eqn_Ax   := diff(Ax(t,X,Y,Z),t)   = Bx(t,X,Y,Z);
eqn_Ay   := diff(Ay(t,X,Y,Z),t)   = By(t,X,Y,Z);
eqn_Az   := diff(Az(t,X,Y,Z),t)   = Bz(t,X,Y,Z);

eqn_pi1_poly := diff(pi1(t, X, Y, Z), t) = (At(t, X, Y, Z)^2*phi1(t, X, Y, Z) - phi1(t, X, Y, Z)*Ax(t, X, Y, Z)^2 - phi1(t, X, Y, Z)*Az(t, X, Y, Z)^2 - phi1(t, X, Y, Z)*Ay(t, X, Y, Z)^2)*e^2 + (-2*At(t, X, Y, Z)*pi2(t, X, Y, Z) + 2*Ax(t, X, Y, Z)*diff(phi2(t, X, Y, Z), X) + 2*Az(t, X, Y, Z)*diff(phi2(t, X, Y, Z), Z) + 2*Ay(t, X, Y, Z)*diff(phi2(t, X, Y, Z), Y))*e + diff(phi1(t, X, Y, Z), X, X) + diff(phi1(t, X, Y, Z), Y, Y) + diff(phi1(t, X, Y, Z), Z, Z) - (1/2)*Vpoly_diff1 - cc*chi(t,X,Y,Z)^2*phi1(t,X,Y,Z);
eqn_pi2_poly := diff(pi2(t, X, Y, Z), t) = (phi2(t, X, Y, Z)*At(t, X, Y, Z)^2 - phi2(t, X, Y, Z)*Ax(t, X, Y, Z)^2 - phi2(t, X, Y, Z)*Az(t, X, Y, Z)^2 - phi2(t, X, Y, Z)*Ay(t, X, Y, Z)^2)*e^2 + (2*pi1(t, X, Y, Z)*At(t, X, Y, Z) - 2*diff(phi1(t, X, Y, Z), Z)*Az(t, X, Y, Z) - 2*diff(phi1(t, X, Y, Z), X)*Ax(t, X, Y, Z) - 2*diff(phi1(t, X, Y, Z), Y)*Ay(t, X, Y, Z))*e + diff(phi2(t, X, Y, Z), Y, Y) + diff(phi2(t, X, Y, Z), Z, Z) + diff(phi2(t, X, Y, Z), X, X) - (1/2)*Vpoly_diff2 - cc*chi(t,X,Y,Z)^2*phi2(t,X,Y,Z);
eqn_pi1_log := diff(pi1(t, X, Y, Z), t) = (At(t, X, Y, Z)^2*phi1(t, X, Y, Z) - phi1(t, X, Y, Z)*Ax(t, X, Y, Z)^2 - phi1(t, X, Y, Z)*Az(t, X, Y, Z)^2 - phi1(t, X, Y, Z)*Ay(t, X, Y, Z)^2)*e^2 + (-2*At(t, X, Y, Z)*pi2(t, X, Y, Z) + 2*Ax(t, X, Y, Z)*diff(phi2(t, X, Y, Z), X) + 2*Az(t, X, Y, Z)*diff(phi2(t, X, Y, Z), Z) + 2*Ay(t, X, Y, Z)*diff(phi2(t, X, Y, Z), Y))*e + diff(phi1(t, X, Y, Z), X, X) + diff(phi1(t, X, Y, Z), Y, Y) + diff(phi1(t, X, Y, Z), Z, Z) - (1/2)*Vlog_diff1 - cc*chi(t,X,Y,Z)^2*phi1(t,X,Y,Z);
eqn_pi2_log := diff(pi2(t, X, Y, Z), t) = (phi2(t, X, Y, Z)*At(t, X, Y, Z)^2 - phi2(t, X, Y, Z)*Ax(t, X, Y, Z)^2 - phi2(t, X, Y, Z)*Az(t, X, Y, Z)^2 - phi2(t, X, Y, Z)*Ay(t, X, Y, Z)^2)*e^2 + (2*pi1(t, X, Y, Z)*At(t, X, Y, Z) - 2*diff(phi1(t, X, Y, Z), Z)*Az(t, X, Y, Z) - 2*diff(phi1(t, X, Y, Z), X)*Ax(t, X, Y, Z) - 2*diff(phi1(t, X, Y, Z), Y)*Ay(t, X, Y, Z))*e + diff(phi2(t, X, Y, Z), Y, Y) + diff(phi2(t, X, Y, Z), Z, Z) + diff(phi2(t, X, Y, Z), X, X) - (1/2)*Vlog_diff2 - cc*chi(t,X,Y,Z)^2*phi2(t,X,Y,Z);

eqn_Bt := diff(Bt(t, X, Y, Z), t) = -(2*phi1(t, X, Y, Z)^2*At(t, X, Y, Z) + 2*phi2(t, X, Y, Z)^2*At(t, X, Y, Z))*e^2 - (-2*phi1(t, X, Y, Z)*pi2(t, X, Y, Z) + 2*phi2(t, X, Y, Z)*pi1(t, X, Y, Z))*e + diff(At(t, X, Y, Z), X, X) + diff(At(t, X, Y, Z), Y, Y) + diff(At(t, X, Y, Z), Z, Z);
eqn_Bx := diff(Bx(t, X, Y, Z), t) = -(2*Ax(t, X, Y, Z)*phi1(t, X, Y, Z)^2 + 2*Ax(t, X, Y, Z)*phi2(t, X, Y, Z)^2)*e^2 - (-2*phi1(t, X, Y, Z)*diff(phi2(t, X, Y, Z), X) + 2*phi2(t, X, Y, Z)*diff(phi1(t, X, Y, Z), X))*e + diff(Ax(t, X, Y, Z), X, X) + diff(Ax(t, X, Y, Z), Y, Y) + diff(Ax(t, X, Y, Z), Z, Z);
eqn_By := diff(By(t, X, Y, Z), t) = -(2*phi1(t, X, Y, Z)^2*Ay(t, X, Y, Z) + 2*phi2(t, X, Y, Z)^2*Ay(t, X, Y, Z))*e^2 - (-2*phi1(t, X, Y, Z)*diff(phi2(t, X, Y, Z), Y) + 2*phi2(t, X, Y, Z)*diff(phi1(t, X, Y, Z), Y))*e + diff(Ay(t, X, Y, Z), Y, Y) + diff(Ay(t, X, Y, Z), X, X) + diff(Ay(t, X, Y, Z), Z, Z);
eqn_Bz := diff(Bz(t, X, Y, Z), t) = -(2*Az(t, X, Y, Z)*phi1(t, X, Y, Z)^2 + 2*Az(t, X, Y, Z)*phi2(t, X, Y, Z)^2)*e^2 - (-2*phi1(t, X, Y, Z)*diff(phi2(t, X, Y, Z), Z) + 2*phi2(t, X, Y, Z)*diff(phi1(t, X, Y, Z), Z))*e + diff(Az(t, X, Y, Z), Z, Z) + diff(Az(t, X, Y, Z), X, X) + diff(Az(t, X, Y, Z), Y, Y);

eqn_chi := diff(chi(t,X,Y,Z),t) = xi(t,X,Y,Z);
eqn_xi  := diff(xi(t,X,Y,Z),t)  = diff(chi(t,X,Y,Z),X,X) + diff(chi(t,X,Y,Z),Y,Y) + diff(chi(t,X,Y,Z),Z,Z) - cc*(phi1(t,X,Y,Z)^2 + phi2(t,X,Y,Z)^2)*chi(t,X,Y,Z);
eqn_chi_init := amp*exp(-((sqrt(((X - chiX0)/wwX)^2 + ((Y - chiY0)/wwY)^2 + ((Z - chiZ0)/wwZ)^2) - r0)/delta)^2);
eqn_xi_init  := (diff(chi(t,X,Y,Z),X)*X + diff(chi(t,X,Y,Z),Y)*Y + diff(chi(t,X,Y,Z),Z)*Z + chi(t,X,Y,Z))/(sqrt(X^2+Y^2+Z^2 + 10^(-20)));  # avoid div by 0 @ X=Y=Z=0

eqn_modphi := modphi(t,X,Y,Z) = sqrt(phi1(t,X,Y,Z)^2 + phi2(t,X,Y,Z)^2);
eqn_qden   := qden(t,X,Y,Z) = 2*(phi2(t,X,Y,Z)*pi1(t,X,Y,Z) - phi1(t,X,Y,Z)*pi2(t,X,Y,Z)) + 2*e*At(t,X,Y,Z)*modphi(t,X,Y,Z)^2;
eqn_eden_poly := eden(t,X,Y,Z) = -2*phi1(t, X, Y, Z)*pi2(t, X, Y, Z)*At(t, X, Y, Z)*e + 2*phi2(t, X, Y, Z)*pi1(t, X, Y, Z)*At(t, X, Y, Z)*e - 2*phi1(t, X, Y, Z)*diff(phi2(t, X, Y, Z), Y)*Ay(t, X, Y, Z)*e + 2*phi2(t, X, Y, Z)*diff(phi1(t, X, Y, Z), Y)*Ay(t, X, Y, Z)*e - 2*phi1(t, X, Y, Z)*diff(phi2(t, X, Y, Z), X)*Ax(t, X, Y, Z)*e + 2*phi2(t, X, Y, Z)*diff(phi1(t, X, Y, Z), X)*Ax(t, X, Y, Z)*e - 2*phi1(t, X, Y, Z)*diff(phi2(t, X, Y, Z), Z)*Az(t, X, Y, Z)*e + 2*phi2(t, X, Y, Z)*diff(phi1(t, X, Y, Z), Z)*Az(t, X, Y, Z)*e + pi2(t, X, Y, Z)^2 + pi1(t, X, Y, Z)^2 + diff(Ay(t, X, Y, Z), X)^2/2 + diff(Ax(t, X, Y, Z), Y)^2/2 + diff(At(t, X, Y, Z), Y)^2/2 + diff(At(t, X, Y, Z), Z)^2/2 + Bz(t, X, Y, Z)^2/2 + diff(At(t, X, Y, Z), X)^2/2 - Bx(t, X, Y, Z)*diff(At(t, X, Y, Z), X) - By(t, X, Y, Z)*diff(At(t, X, Y, Z), Y) - Bz(t, X, Y, Z)*diff(At(t, X, Y, Z), Z) - diff(Ay(t, X, Y, Z), X)*diff(Ax(t, X, Y, Z), Y) - diff(Az(t, X, Y, Z), X)*diff(Ax(t, X, Y, Z), Z) - diff(Az(t, X, Y, Z), Y)*diff(Ay(t, X, Y, Z), Z) + diff(Az(t, X, Y, Z), Y)^2/2 + diff(Ay(t, X, Y, Z), Z)^2/2 + diff(Az(t, X, Y, Z), X)^2/2 + diff(Ax(t, X, Y, Z), Z)^2/2 + Bx(t, X, Y, Z)^2/2 + By(t, X, Y, Z)^2/2 + diff(phi2(t, X, Y, Z), Z)^2 + diff(phi2(t, X, Y, Z), Y)^2 + diff(phi2(t, X, Y, Z), X)^2 + diff(phi1(t, X, Y, Z), Z)^2 + diff(phi1(t, X, Y, Z), Y)^2 + diff(phi1(t, X, Y, Z), X)^2 + e^2*(Az(t, X, Y, Z)^2 + At(t, X, Y, Z)^2 + Ax(t, X, Y, Z)^2 + Ay(t, X, Y, Z)^2)*phi2(t, X, Y, Z)^2 + e^2*(Az(t, X, Y, Z)^2 + At(t, X, Y, Z)^2 + Ax(t, X, Y, Z)^2 + Ay(t, X, Y, Z)^2)*phi1(t, X, Y, Z)^2 + Vpoly;
eqn_eden_log := eden(t,X,Y,Z) = -2*phi1(t, X, Y, Z)*pi2(t, X, Y, Z)*At(t, X, Y, Z)*e + 2*phi2(t, X, Y, Z)*pi1(t, X, Y, Z)*At(t, X, Y, Z)*e - 2*phi1(t, X, Y, Z)*diff(phi2(t, X, Y, Z), Y)*Ay(t, X, Y, Z)*e + 2*phi2(t, X, Y, Z)*diff(phi1(t, X, Y, Z), Y)*Ay(t, X, Y, Z)*e - 2*phi1(t, X, Y, Z)*diff(phi2(t, X, Y, Z), X)*Ax(t, X, Y, Z)*e + 2*phi2(t, X, Y, Z)*diff(phi1(t, X, Y, Z), X)*Ax(t, X, Y, Z)*e - 2*phi1(t, X, Y, Z)*diff(phi2(t, X, Y, Z), Z)*Az(t, X, Y, Z)*e + 2*phi2(t, X, Y, Z)*diff(phi1(t, X, Y, Z), Z)*Az(t, X, Y, Z)*e + pi2(t, X, Y, Z)^2 + pi1(t, X, Y, Z)^2 + diff(Ay(t, X, Y, Z), X)^2/2 + diff(Ax(t, X, Y, Z), Y)^2/2 + diff(At(t, X, Y, Z), Y)^2/2 + diff(At(t, X, Y, Z), Z)^2/2 + Bz(t, X, Y, Z)^2/2 + diff(At(t, X, Y, Z), X)^2/2 - Bx(t, X, Y, Z)*diff(At(t, X, Y, Z), X) - By(t, X, Y, Z)*diff(At(t, X, Y, Z), Y) - Bz(t, X, Y, Z)*diff(At(t, X, Y, Z), Z) - diff(Ay(t, X, Y, Z), X)*diff(Ax(t, X, Y, Z), Y) - diff(Az(t, X, Y, Z), X)*diff(Ax(t, X, Y, Z), Z) - diff(Az(t, X, Y, Z), Y)*diff(Ay(t, X, Y, Z), Z) + diff(Az(t, X, Y, Z), Y)^2/2 + diff(Ay(t, X, Y, Z), Z)^2/2 + diff(Az(t, X, Y, Z), X)^2/2 + diff(Ax(t, X, Y, Z), Z)^2/2 + Bx(t, X, Y, Z)^2/2 + By(t, X, Y, Z)^2/2 + diff(phi2(t, X, Y, Z), Z)^2 + diff(phi2(t, X, Y, Z), Y)^2 + diff(phi2(t, X, Y, Z), X)^2 + diff(phi1(t, X, Y, Z), Z)^2 + diff(phi1(t, X, Y, Z), Y)^2 + diff(phi1(t, X, Y, Z), X)^2 + e^2*(Az(t, X, Y, Z)^2 + At(t, X, Y, Z)^2 + Ax(t, X, Y, Z)^2 + Ay(t, X, Y, Z)^2)*phi2(t, X, Y, Z)^2 + e^2*(Az(t, X, Y, Z)^2 + At(t, X, Y, Z)^2 + Ax(t, X, Y, Z)^2 + Ay(t, X, Y, Z)^2)*phi1(t, X, Y, Z)^2 + Vlog;
eqn_jdenx := jdenx(t,X,Y,Z) = (Y*diff(Ax(t, X, Y, Z), Z) - Y*diff(Az(t, X, Y, Z), X) - Z*diff(Ax(t, X, Y, Z), Y) + Z*diff(Ay(t, X, Y, Z), X))*diff(At(t, X, Y, Z), X) + (-By(t, X, Y, Z)*Y - Bz(t, X, Y, Z)*Z + Y*diff(At(t, X, Y, Z), Y) + Z*diff(At(t, X, Y, Z), Z))*diff(Ay(t, X, Y, Z), Z) + (By(t, X, Y, Z)*Y + Bz(t, X, Y, Z)*Z - Y*diff(At(t, X, Y, Z), Y) - Z*diff(At(t, X, Y, Z), Z))*diff(Az(t, X, Y, Z), Y) + 2*Z*(phi2(t, X, Y, Z)*e*At(t, X, Y, Z) + pi1(t, X, Y, Z))*diff(phi1(t, X, Y, Z), Y) - 2*Y*(phi2(t, X, Y, Z)*e*At(t, X, Y, Z) + pi1(t, X, Y, Z))*diff(phi1(t, X, Y, Z), Z) + 2*Z*(-phi1(t, X, Y, Z)*e*At(t, X, Y, Z) + pi2(t, X, Y, Z))*diff(phi2(t, X, Y, Z), Y) - 2*Y*(-phi1(t, X, Y, Z)*e*At(t, X, Y, Z) + pi2(t, X, Y, Z))*diff(phi2(t, X, Y, Z), Z) + Bx(t, X, Y, Z)*diff(Ax(t, X, Y, Z), Y)*Z - Bx(t, X, Y, Z)*diff(Ax(t, X, Y, Z), Z)*Y - Bx(t, X, Y, Z)*diff(Ay(t, X, Y, Z), X)*Z + diff(Az(t, X, Y, Z), X)*Bx(t, X, Y, Z)*Y - 2*(Y*Az(t, X, Y, Z) - Z*Ay(t, X, Y, Z))*e*(e*(phi1(t, X, Y, Z)^2 + phi2(t, X, Y, Z)^2)*At(t, X, Y, Z) + pi1(t, X, Y, Z)*phi2(t, X, Y, Z) - pi2(t, X, Y, Z)*phi1(t, X, Y, Z));
eqn_lorenz := lorenz(t,X,Y,Z) = -Bt(t,X,Y,Z) + diff(Ax(t,X,Y,Z), X) + diff(Ay(t,X,Y,Z), Y) + diff(Az(t,X,Y,Z), Z);
eqn_elecx  := elecx(t,X,Y,Z) = diff(At(t,X,Y,Z), X) - Bx(t,X,Y,Z);
eqn_elecy  := elecy(t,X,Y,Z) = diff(At(t,X,Y,Z), Y) - By(t,X,Y,Z);
eqn_elecz  := elecz(t,X,Y,Z) = diff(At(t,X,Y,Z), Z) - Bz(t,X,Y,Z);
eqn_magx   := magx(t,X,Y,Z) = diff(Az(t,X,Y,Z), Y) - diff(Ay(t,X,Y,Z), Z);
eqn_magy   := magy(t,X,Y,Z) = diff(Ax(t,X,Y,Z), Z) - diff(Az(t,X,Y,Z), X);
eqn_magz   := magz(t,X,Y,Z) = diff(Ay(t,X,Y,Z), X) - diff(Ax(t,X,Y,Z), Y);
eqn_eem    := eem(t,X,Y,Z) = 0.5*(elecx(t,X,Y,Z)^2+elecy(t,X,Y,Z)^2+elecz(t,X,Y,Z)^2+magx(t,X,Y,Z)^2+magy(t,X,Y,Z)^2+magz(t,X,Y,Z)^2);
eqn_gaussE := gaussE(t,X,Y,Z) = diff(At(t,X,Y,Z), X, X) - diff(Bx(t,X,Y,Z), X) + diff(At(t,X,Y,Z), Y, Y) - diff(By(t,X,Y,Z), Y) + diff(At(t,X,Y,Z), Z, Z) - diff(Bz(t,X,Y,Z), Z) - e*qden(t,X,Y,Z);
eqn_gaussB := gaussB(t,X,Y,Z) = diff(magx(t,X,Y,Z), X) + diff(magy(t,X,Y,Z), Y) + diff(magz(t,X,Y,Z), Z);

eqn_phi1_ire_poly := diff(phi1(t,X,Y,Z),t,t) = rhs(eqn_pi1_poly);
eqn_phi1_ire_log  := diff(phi1(t,X,Y,Z),t,t) = rhs(eqn_pi1_log);
eqn_phi2_ire_poly := diff(phi2(t,X,Y,Z),t,t) = rhs(eqn_pi2_poly);
eqn_phi2_ire_log  := diff(phi2(t,X,Y,Z),t,t) = rhs(eqn_pi2_log);
eqn_At_ire := diff(At(t,X,Y,Z),t,t) = rhs(eqn_Bt);
eqn_Ax_ire := diff(Ax(t,X,Y,Z),t,t) = rhs(eqn_Bx);
eqn_Ay_ire := diff(Ay(t,X,Y,Z),t,t) = rhs(eqn_By);
eqn_Az_ire := diff(Az(t,X,Y,Z),t,t) = rhs(eqn_Bz);

chain:={
  diff(f(t,X,Y,Z),X)   = diff(f(t,x(X),y(Y),z(Z)),X),
  diff(f(t,X,Y,Z),Y)   = diff(f(t,x(X),y(Y),z(Z)),Y),
  diff(f(t,X,Y,Z),Z)   = diff(f(t,x(X),y(Y),z(Z)),Z),
  diff(f(t,X,Y,Z),X,X) = diff(f(t,x(X),y(Y),z(Z)),X,X),
  diff(f(t,X,Y,Z),Y,Y) = diff(f(t,x(X),y(Y),z(Z)),Y,Y),
  diff(f(t,X,Y,Z),Z,Z) = diff(f(t,x(X),y(Y),z(Z)),Z,Z),
  diff(f(t,X,Y,Z),t)   = diff(f(t,x,y,z),t)
};

subvars1 := subs(
  D[2](f)(t,x(X),y(Y),z(Z))   = diff(f(t,x,y,z),x),
  D[3](f)(t,x(X),y(Y),z(Z))   = diff(f(t,x,y,z),y),
  D[4](f)(t,x(X),y(Y),z(Z))   = diff(f(t,x,y,z),z),
  D[2,2](f)(t,x(X),y(Y),z(Z)) = diff(f(t,x,y,z),x,x),
  D[3,3](f)(t,x(X),y(Y),z(Z)) = diff(f(t,x,y,z),y,y),
  D[4,4](f)(t,x(X),y(Y),z(Z)) = diff(f(t,x,y,z),z,z),
  diff(x(X),X,X) = dxdX2(t,x,y,z),
  diff(x(X),X)   = dxdX(t,x,y,z),
  diff(y(Y),Y,Y) = dydY2(t,x,y,z),
  diff(y(Y),Y)   = dydY(t,x,y,z),
  diff(z(Z),Z,Z) = dzdZ2(t,x,y,z),
  diff(z(Z),Z)   = dzdZ(t,x,y,z),
chain);

subvars2 := f(t,X,Y,Z) = f(t,x,y,z);

subvars3 := X = compactificationX, Y = compactificationY, Z = compactificationZ;

count := 1:
for f in grid_functions do
  out1[count] := subvars1[];
  out2[count] := subvars2;
  out3[count] := subvars3;
  count := count+1:
end do:
fullvars1 := convert(out1,set);
fullvars2 := convert(out2,set);
fullvars3 := convert(out3,set);
f := 'f':

count := 1:
for f in eqn_list do
  eqn[count] :=
    subs(
      fullvars1,
      fullvars2,
      fullvars3,
    f);
  print(%):
  count := count+1;
end do:
f := 'f':

eqn_phi1     := eqn[1];
eqn_pi1_poly := eqn[2];
eqn_pi1_log  := eqn[3];
eqn_phi2     := eqn[4];
eqn_pi2_poly := eqn[5];
eqn_pi2_log  := eqn[6];
eqn_At   := eqn[7];
eqn_Bt   := eqn[8];
eqn_Ax   := eqn[9];
eqn_Bx   := eqn[10];
eqn_Ay   := eqn[11];
eqn_By   := eqn[12];
eqn_Az   := eqn[13];
eqn_Bz   := eqn[14];
eqn_chi  := eqn[15];
eqn_xi   := eqn[16];
eqn_chi_init  := eqn[17];
eqn_xi_init   := eqn[18];
eqn_modphi    := eqn[19];
eqn_qden      := eqn[20];
eqn_eden_poly := eqn[21];
eqn_eden_log  := eqn[22];
eqn_jdenx    := eqn[23];
eqn_lorenz  := eqn[24];
eqn_elecx := eqn[25];
eqn_elecy := eqn[26];
eqn_elecz := eqn[27];
eqn_magx := eqn[28];
eqn_magy := eqn[29];
eqn_magz := eqn[30];
eqn_eem  := eqn[31];
eqn_gaussE := eqn[32];
eqn_gaussB := eqn[33];
eqn_phi1_ire_poly := eqn[34];
eqn_phi1_ire_log  := eqn[35];
eqn_phi2_ire_poly := eqn[36];
eqn_phi2_ire_log  := eqn[37];
eqn_At_ire := eqn[38];
eqn_Ax_ire := eqn[39];
eqn_Ay_ire := eqn[40];
eqn_Az_ire := eqn[41];

# independent residual equations will use default second-order stencils
Gen_Res_Code(lhs(eqn_phi1_ire_poly)-rhs(eqn_phi1_ire_poly), input="c", proc_name="ire_phi1_poly");
Gen_Res_Code(lhs(eqn_phi1_ire_log)-rhs(eqn_phi1_ire_log),   input="c", proc_name="ire_phi1_log");
Gen_Res_Code(lhs(eqn_phi2_ire_poly)-rhs(eqn_phi2_ire_poly), input="c", proc_name="ire_phi2_poly");
Gen_Res_Code(lhs(eqn_phi2_ire_log)-rhs(eqn_phi2_ire_log),   input="c", proc_name="ire_phi2_log");
Gen_Res_Code(lhs(eqn_At_ire)-rhs(eqn_At_ire), input="c", proc_name="ire_at");
Gen_Res_Code(lhs(eqn_Ax_ire)-rhs(eqn_Ax_ire), input="c", proc_name="ire_ax");
Gen_Res_Code(lhs(eqn_Ay_ire)-rhs(eqn_Ay_ire), input="c", proc_name="ire_ay");
Gen_Res_Code(lhs(eqn_Az_ire)-rhs(eqn_Az_ire), input="c", proc_name="ire_az");

# switch to fourth-order stencils from now on
fds:=table([ t=[-1,-1],x=[-1,-1],y=[-1,-1],z=[-1,-1] ]):
Update_FD_Table(4,fds):
SFDT();

FD_table[t] := [[0],[0,1]];  # use only n and np1 for time integration (RK4)
NT := f -> FD(f,[[1],[0,0,0]]);

eqn_phi1_R_D     := Gen_Sten(rhs(eqn_phi1));
eqn_pi1_poly_R_D := Gen_Sten(rhs(eqn_pi1_poly));
eqn_pi1_log_R_D  := Gen_Sten(rhs(eqn_pi1_log));
eqn_phi2_R_D     := Gen_Sten(rhs(eqn_phi2));
eqn_pi2_poly_R_D := Gen_Sten(rhs(eqn_pi2_poly));
eqn_pi2_log_R_D  := Gen_Sten(rhs(eqn_pi2_log));

eqn_phi1_L_D     := Gen_Sten(lhs(eqn_phi1));
eqn_pi1_poly_L_D := Gen_Sten(lhs(eqn_pi1_poly));
eqn_pi1_log_L_D  := Gen_Sten(lhs(eqn_pi1_log));
eqn_phi2_L_D     := Gen_Sten(lhs(eqn_phi2));
eqn_pi2_poly_L_D := Gen_Sten(lhs(eqn_pi2_poly));
eqn_pi2_log_L_D  := Gen_Sten(lhs(eqn_pi2_log));

eqn_At_R_D := Gen_Sten(rhs(eqn_At));
eqn_Bt_R_D := Gen_Sten(rhs(eqn_Bt));
eqn_Ax_R_D := Gen_Sten(rhs(eqn_Ax));
eqn_Bx_R_D := Gen_Sten(rhs(eqn_Bx));
eqn_Ay_R_D := Gen_Sten(rhs(eqn_Ay));
eqn_By_R_D := Gen_Sten(rhs(eqn_By));
eqn_Az_R_D := Gen_Sten(rhs(eqn_Az));
eqn_Bz_R_D := Gen_Sten(rhs(eqn_Bz));

eqn_At_L_D := Gen_Sten(lhs(eqn_At));
eqn_Bt_L_D := Gen_Sten(lhs(eqn_Bt));
eqn_Ax_L_D := Gen_Sten(lhs(eqn_Ax));
eqn_Bx_L_D := Gen_Sten(lhs(eqn_Bx));
eqn_Ay_L_D := Gen_Sten(lhs(eqn_Ay));
eqn_By_L_D := Gen_Sten(lhs(eqn_By));
eqn_Az_L_D := Gen_Sten(lhs(eqn_Az));
eqn_Bz_L_D := Gen_Sten(lhs(eqn_Bz));

eqn_chi_R_D  := Gen_Sten(rhs(eqn_chi));
eqn_xi_R_D   := Gen_Sten(rhs(eqn_xi));

eqn_chi_L_D  := Gen_Sten(lhs(eqn_chi));
eqn_xi_L_D   := Gen_Sten(lhs(eqn_xi));

spec_evol_phi1 := [
  { i=[3,Nx-2,1] , j=[3,Ny-2,1] , k=[3,Nz-2,1] }        = eqn_phi1_R_D,
  { i=[1,1,1]    , j=[1,Ny,1]   , k=[1,Nz,1],  b=xmin } = phi1(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[Nx,Nx,1]  , j=[1,Ny,1]   , k=[1,Nz,1],  b=xmax } = phi1(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,1,1]    , k=[1,Nz,1],  b=ymin } = phi1(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[Ny,Ny,1]  , k=[1,Nz,1],  b=ymax } = phi1(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,Ny,1]   , k=[1,1,1],   b=zmin } = phi1(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,Ny,1]   , k=[Nz,Nz,1], b=zmax } = phi1(n+1,i,j,k) - myzero*x(i)*y(j)*z(k)
];

spec_evol_pi1_poly := [
  { i=[3,Nx-2,1] , j=[3,Ny-2,1] , k=[3,Nz-2,1] }        = eqn_pi1_poly_R_D,
  { i=[1,1,1]    , j=[1,Ny,1]   , k=[1,Nz,1],  b=xmin } = pi1(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[Nx,Nx,1]  , j=[1,Ny,1]   , k=[1,Nz,1],  b=xmax } = pi1(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,1,1]    , k=[1,Nz,1],  b=ymin } = pi1(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[Ny,Ny,1]  , k=[1,Nz,1],  b=ymax } = pi1(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,Ny,1]   , k=[1,1,1],   b=zmin } = pi1(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,Ny,1]   , k=[Nz,Nz,1], b=zmax } = pi1(n+1,i,j,k) - myzero*x(i)*y(j)*z(k)
];

spec_evol_pi1_log := [
  { i=[3,Nx-2,1] , j=[3,Ny-2,1] , k=[3,Nz-2,1] }        = eqn_pi1_log_R_D,
  { i=[1,1,1]    , j=[1,Ny,1]   , k=[1,Nz,1],  b=xmin } = pi1(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[Nx,Nx,1]  , j=[1,Ny,1]   , k=[1,Nz,1],  b=xmax } = pi1(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,1,1]    , k=[1,Nz,1],  b=ymin } = pi1(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[Ny,Ny,1]  , k=[1,Nz,1],  b=ymax } = pi1(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,Ny,1]   , k=[1,1,1],   b=zmin } = pi1(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,Ny,1]   , k=[Nz,Nz,1], b=zmax } = pi1(n+1,i,j,k) - myzero*x(i)*y(j)*z(k)
];

spec_evol_phi2 := [
  { i=[3,Nx-2,1] , j=[3,Ny-2,1] , k=[3,Nz-2,1] }        = eqn_phi2_R_D,
  { i=[1,1,1]    , j=[1,Ny,1]   , k=[1,Nz,1],  b=xmin } = phi2(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[Nx,Nx,1]  , j=[1,Ny,1]   , k=[1,Nz,1],  b=xmax } = phi2(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,1,1]    , k=[1,Nz,1],  b=ymin } = phi2(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[Ny,Ny,1]  , k=[1,Nz,1],  b=ymax } = phi2(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,Ny,1]   , k=[1,1,1],   b=zmin } = phi2(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,Ny,1]   , k=[Nz,Nz,1], b=zmax } = phi2(n+1,i,j,k) - myzero*x(i)*y(j)*z(k)
];

spec_evol_pi2_poly := [
  { i=[3,Nx-2,1] , j=[3,Ny-2,1] , k=[3,Nz-2,1] }        = eqn_pi2_poly_R_D,
  { i=[1,1,1]    , j=[1,Ny,1]   , k=[1,Nz,1],  b=xmin } = pi2(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[Nx,Nx,1]  , j=[1,Ny,1]   , k=[1,Nz,1],  b=xmax } = pi2(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,1,1]    , k=[1,Nz,1],  b=ymin } = pi2(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[Ny,Ny,1]  , k=[1,Nz,1],  b=ymax } = pi2(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,Ny,1]   , k=[1,1,1],   b=zmin } = pi2(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,Ny,1]   , k=[Nz,Nz,1], b=zmax } = pi2(n+1,i,j,k) - myzero*x(i)*y(j)*z(k)
];

spec_evol_pi2_log := [
  { i=[3,Nx-2,1] , j=[3,Ny-2,1] , k=[3,Nz-2,1] }        = eqn_pi2_log_R_D,
  { i=[1,1,1]    , j=[1,Ny,1]   , k=[1,Nz,1],  b=xmin } = pi2(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[Nx,Nx,1]  , j=[1,Ny,1]   , k=[1,Nz,1],  b=xmax } = pi2(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,1,1]    , k=[1,Nz,1],  b=ymin } = pi2(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[Ny,Ny,1]  , k=[1,Nz,1],  b=ymax } = pi2(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,Ny,1]   , k=[1,1,1],   b=zmin } = pi2(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,Ny,1]   , k=[Nz,Nz,1], b=zmax } = pi2(n+1,i,j,k) - myzero*x(i)*y(j)*z(k)
];

spec_evol_At := [
  { i=[3,Nx-2,1] , j=[3,Ny-2,1] , k=[3,Nz-2,1] }        = eqn_At_R_D,
  { i=[1,1,1]    , j=[1,Ny,1]   , k=[1,Nz,1],  b=xmin } = At(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[Nx,Nx,1]  , j=[1,Ny,1]   , k=[1,Nz,1],  b=xmax } = At(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,1,1]    , k=[1,Nz,1],  b=ymin } = At(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[Ny,Ny,1]  , k=[1,Nz,1],  b=ymax } = At(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,Ny,1]   , k=[1,1,1],   b=zmin } = At(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,Ny,1]   , k=[Nz,Nz,1], b=zmax } = At(n+1,i,j,k) - myzero*x(i)*y(j)*z(k)
];

spec_evol_Bt := [
  { i=[3,Nx-2,1] , j=[3,Ny-2,1] , k=[3,Nz-2,1] }        = eqn_Bt_R_D,
  { i=[1,1,1]    , j=[1,Ny,1]   , k=[1,Nz,1],  b=xmin } = Bt(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[Nx,Nx,1]  , j=[1,Ny,1]   , k=[1,Nz,1],  b=xmax } = Bt(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,1,1]    , k=[1,Nz,1],  b=ymin } = Bt(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[Ny,Ny,1]  , k=[1,Nz,1],  b=ymax } = Bt(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,Ny,1]   , k=[1,1,1],   b=zmin } = Bt(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,Ny,1]   , k=[Nz,Nz,1], b=zmax } = Bt(n+1,i,j,k) - myzero*x(i)*y(j)*z(k)
];

spec_evol_Ax := [
  { i=[3,Nx-2,1] , j=[3,Ny-2,1] , k=[3,Nz-2,1] }        = eqn_Ax_R_D,
  { i=[1,1,1]    , j=[1,Ny,1]   , k=[1,Nz,1],  b=xmin } = Ax(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[Nx,Nx,1]  , j=[1,Ny,1]   , k=[1,Nz,1],  b=xmax } = Ax(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,1,1]    , k=[1,Nz,1],  b=ymin } = Ax(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[Ny,Ny,1]  , k=[1,Nz,1],  b=ymax } = Ax(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,Ny,1]   , k=[1,1,1],   b=zmin } = Ax(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,Ny,1]   , k=[Nz,Nz,1], b=zmax } = Ax(n+1,i,j,k) - myzero*x(i)*y(j)*z(k)
];

spec_evol_Bx := [
  { i=[3,Nx-2,1] , j=[3,Ny-2,1] , k=[3,Nz-2,1] }        = eqn_Bx_R_D,
  { i=[1,1,1]    , j=[1,Ny,1]   , k=[1,Nz,1],  b=xmin } = Bx(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[Nx,Nx,1]  , j=[1,Ny,1]   , k=[1,Nz,1],  b=xmax } = Bx(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,1,1]    , k=[1,Nz,1],  b=ymin } = Bx(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[Ny,Ny,1]  , k=[1,Nz,1],  b=ymax } = Bx(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,Ny,1]   , k=[1,1,1],   b=zmin } = Bx(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,Ny,1]   , k=[Nz,Nz,1], b=zmax } = Bx(n+1,i,j,k) - myzero*x(i)*y(j)*z(k)
];

spec_evol_Ay := [
  { i=[3,Nx-2,1] , j=[3,Ny-2,1] , k=[3,Nz-2,1] }        = eqn_Ay_R_D,
  { i=[1,1,1]    , j=[1,Ny,1]   , k=[1,Nz,1],  b=xmin } = Ay(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[Nx,Nx,1]  , j=[1,Ny,1]   , k=[1,Nz,1],  b=xmax } = Ay(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,1,1]    , k=[1,Nz,1],  b=ymin } = Ay(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[Ny,Ny,1]  , k=[1,Nz,1],  b=ymax } = Ay(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,Ny,1]   , k=[1,1,1],   b=zmin } = Ay(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,Ny,1]   , k=[Nz,Nz,1], b=zmax } = Ay(n+1,i,j,k) - myzero*x(i)*y(j)*z(k)
];

spec_evol_By := [
  { i=[3,Nx-2,1] , j=[3,Ny-2,1] , k=[3,Nz-2,1] }        = eqn_By_R_D,
  { i=[1,1,1]    , j=[1,Ny,1]   , k=[1,Nz,1],  b=xmin } = By(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[Nx,Nx,1]  , j=[1,Ny,1]   , k=[1,Nz,1],  b=xmax } = By(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,1,1]    , k=[1,Nz,1],  b=ymin } = By(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[Ny,Ny,1]  , k=[1,Nz,1],  b=ymax } = By(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,Ny,1]   , k=[1,1,1],   b=zmin } = By(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,Ny,1]   , k=[Nz,Nz,1], b=zmax } = By(n+1,i,j,k) - myzero*x(i)*y(j)*z(k)
];

spec_evol_Az := [
  { i=[3,Nx-2,1] , j=[3,Ny-2,1] , k=[3,Nz-2,1] }        = eqn_Az_R_D,
  { i=[1,1,1]    , j=[1,Ny,1]   , k=[1,Nz,1],  b=xmin } = Az(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[Nx,Nx,1]  , j=[1,Ny,1]   , k=[1,Nz,1],  b=xmax } = Az(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,1,1]    , k=[1,Nz,1],  b=ymin } = Az(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[Ny,Ny,1]  , k=[1,Nz,1],  b=ymax } = Az(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,Ny,1]   , k=[1,1,1],   b=zmin } = Az(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,Ny,1]   , k=[Nz,Nz,1], b=zmax } = Az(n+1,i,j,k) - myzero*x(i)*y(j)*z(k)
];

spec_evol_Bz := [
  { i=[3,Nx-2,1] , j=[3,Ny-2,1] , k=[3,Nz-2,1] }        = eqn_Bz_R_D,
  { i=[1,1,1]    , j=[1,Ny,1]   , k=[1,Nz,1],  b=xmin } = Bz(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[Nx,Nx,1]  , j=[1,Ny,1]   , k=[1,Nz,1],  b=xmax } = Bz(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,1,1]    , k=[1,Nz,1],  b=ymin } = Bz(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[Ny,Ny,1]  , k=[1,Nz,1],  b=ymax } = Bz(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,Ny,1]   , k=[1,1,1],   b=zmin } = Bz(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,Ny,1]   , k=[Nz,Nz,1], b=zmax } = Bz(n+1,i,j,k) - myzero*x(i)*y(j)*z(k)
];

spec_evol_chi := [
  { i=[3,Nx-2,1] , j=[3,Ny-2,1] , k=[3,Nz-2,1] }        = eqn_chi_R_D,
  { i=[1,1,1]    , j=[1,Ny,1]   , k=[1,Nz,1],  b=xmin } = chi(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[Nx,Nx,1]  , j=[1,Ny,1]   , k=[1,Nz,1],  b=xmax } = chi(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,1,1]    , k=[1,Nz,1],  b=ymin } = chi(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[Ny,Ny,1]  , k=[1,Nz,1],  b=ymax } = chi(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,Ny,1]   , k=[1,1,1],   b=zmin } = chi(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,Ny,1]   , k=[Nz,Nz,1], b=zmax } = chi(n+1,i,j,k) - myzero*x(i)*y(j)*z(k)
];

spec_evol_xi := [
  { i=[3,Nx-2,1] , j=[3,Ny-2,1] , k=[3,Nz-2,1] }        = eqn_xi_R_D,
  { i=[1,1,1]    , j=[1,Ny,1]   , k=[1,Nz,1],  b=xmin } = xi(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[Nx,Nx,1]  , j=[1,Ny,1]   , k=[1,Nz,1],  b=xmax } = xi(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,1,1]    , k=[1,Nz,1],  b=ymin } = xi(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[Ny,Ny,1]  , k=[1,Nz,1],  b=ymax } = xi(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,Ny,1]   , k=[1,1,1],   b=zmin } = xi(n+1,i,j,k) - myzero*x(i)*y(j)*z(k),
  { i=[1,Nx,1]   , j=[1,Ny,1]   , k=[Nz,Nz,1], b=zmax } = xi(n+1,i,j,k) - myzero*x(i)*y(j)*z(k)
];

A_Gen_Eval_Code(spec_evol_phi1,     input="d", proc_name="u_phi1_rk");
A_Gen_Eval_Code(spec_evol_pi1_poly, input="d", proc_name="u_pi1_poly_rk");
A_Gen_Eval_Code(spec_evol_pi1_log,  input="d", proc_name="u_pi1_log_rk");
A_Gen_Eval_Code(spec_evol_phi2,     input="d", proc_name="u_phi2_rk");
A_Gen_Eval_Code(spec_evol_pi2_poly, input="d", proc_name="u_pi2_poly_rk");
A_Gen_Eval_Code(spec_evol_pi2_log,  input="d", proc_name="u_pi2_log_rk");

A_Gen_Eval_Code(spec_evol_At, input="d", proc_name="u_at_rk");
A_Gen_Eval_Code(spec_evol_Bt, input="d", proc_name="u_bt_rk");
A_Gen_Eval_Code(spec_evol_Ax, input="d", proc_name="u_ax_rk");
A_Gen_Eval_Code(spec_evol_Bx, input="d", proc_name="u_bx_rk");
A_Gen_Eval_Code(spec_evol_Ay, input="d", proc_name="u_ay_rk");
A_Gen_Eval_Code(spec_evol_By, input="d", proc_name="u_by_rk");
A_Gen_Eval_Code(spec_evol_Az, input="d", proc_name="u_az_rk");
A_Gen_Eval_Code(spec_evol_Bz, input="d", proc_name="u_bz_rk");

A_Gen_Eval_Code(spec_evol_chi,  input="d", proc_name="u_chi_rk");
A_Gen_Eval_Code(spec_evol_xi,   input="d", proc_name="u_xi_rk");

eqn_modphi_R_D := Gen_Sten(rhs(eqn_modphi));
eqn_qden_R_D   := Gen_Sten(rhs(eqn_qden));
eqn_eden_poly_R_D := Gen_Sten(rhs(eqn_eden_poly));
eqn_eden_log_R_D  := Gen_Sten(rhs(eqn_eden_log));
eqn_jdenx_R_D   := Gen_Sten(rhs(eqn_jdenx));
eqn_lorenz_R_D := Gen_Sten(rhs(eqn_lorenz));
eqn_elecx_R_D  := Gen_Sten(rhs(eqn_elecx));
eqn_elecy_R_D  := Gen_Sten(rhs(eqn_elecy));
eqn_elecz_R_D  := Gen_Sten(rhs(eqn_elecz));
eqn_magx_R_D   := Gen_Sten(rhs(eqn_magx));
eqn_magy_R_D   := Gen_Sten(rhs(eqn_magy));
eqn_magz_R_D   := Gen_Sten(rhs(eqn_magz));
eqn_eem_R_D    := Gen_Sten(rhs(eqn_eem));
eqn_gaussE_R_D := Gen_Sten(rhs(eqn_gaussE));
eqn_gaussB_R_D := Gen_Sten(rhs(eqn_gaussB));

eqn_modphi_L_D := Gen_Sten(lhs(eqn_modphi));
eqn_qden_L_D   := Gen_Sten(lhs(eqn_qden));
eqn_eden_poly_L_D := Gen_Sten(lhs(eqn_eden_poly));
eqn_eden_log_L_D  := Gen_Sten(lhs(eqn_eden_log));
eqn_jdenx_L_D   := Gen_Sten(lhs(eqn_jdenx));
eqn_lorenz_L_D := Gen_Sten(lhs(eqn_lorenz));
eqn_elecx_L_D  := Gen_Sten(lhs(eqn_elecx));
eqn_elecy_L_D  := Gen_Sten(lhs(eqn_elecy));
eqn_elecz_L_D  := Gen_Sten(lhs(eqn_elecz));
eqn_magx_L_D   := Gen_Sten(lhs(eqn_magx));
eqn_magy_L_D   := Gen_Sten(lhs(eqn_magy));
eqn_magz_L_D   := Gen_Sten(lhs(eqn_magz));
eqn_eem_L_D    := Gen_Sten(lhs(eqn_eem));
eqn_gaussE_L_D := Gen_Sten(lhs(eqn_gaussE));
eqn_gaussB_L_D := Gen_Sten(lhs(eqn_gaussB));

eqn_modphi_D := NT(eqn_modphi_L_D - eqn_modphi_R_D);
eqn_qden_D   := NT(eqn_qden_L_D - eqn_qden_R_D);
eqn_eden_poly_D := NT(eqn_eden_poly_L_D - eqn_eden_poly_R_D);
eqn_eden_log_D  := NT(eqn_eden_log_L_D - eqn_eden_log_R_D);
eqn_jdenx_D   := NT(eqn_jdenx_L_D - eqn_jdenx_R_D);
eqn_lorenz_D := NT(eqn_lorenz_L_D - eqn_lorenz_R_D);
eqn_elecx_D  := NT(eqn_elecx_L_D - eqn_elecx_R_D);
eqn_elecy_D  := NT(eqn_elecy_L_D - eqn_elecy_R_D);
eqn_elecz_D  := NT(eqn_elecz_L_D - eqn_elecz_R_D);
eqn_magx_D   := NT(eqn_magx_L_D - eqn_magx_R_D);
eqn_magy_D   := NT(eqn_magy_L_D - eqn_magy_R_D);
eqn_magz_D   := NT(eqn_magz_L_D - eqn_magz_R_D);
eqn_eem_D    := NT(eqn_eem_L_D - eqn_eem_R_D);
eqn_gaussE_D := NT(eqn_gaussE_L_D - eqn_gaussE_R_D);
eqn_gaussB_D := NT(eqn_gaussB_L_D - eqn_gaussB_R_D);

spec_evol_modphi := [
  { i=[1,Nx,1],   j=[1,Ny,1],   k=[1,Nz,1] }   = eqn_modphi_D
];

spec_evol_qden   := [
  { i=[1,Nx,1],   j=[1,Ny,1],   k=[1,Nz,1] }   = eqn_qden_D
];

spec_evol_eden_poly := [
  { i=[3,Nx-2,1], j=[3,Ny-2,1], k=[3,Nz-2,1] } = eqn_eden_poly_D
];

spec_evol_eden_log  := [
  { i=[3,Nx-2,1], j=[3,Ny-2,1], k=[3,Nz-2,1] } = eqn_eden_log_D
];

spec_evol_jdenx  := [
  { i=[3,Nx-2,1], j=[3,Ny-2,1], k=[3,Nz-2,1] } = eqn_jdenx_D
];

spec_evol_lorenz  := [
  { i=[3,Nx-2,1], j=[3,Ny-2,1], k=[3,Nz-2,1] } = eqn_lorenz_D
];

spec_evol_elecx  := [
  { i=[3,Nx-2,1], j=[3,Ny-2,1], k=[3,Nz-2,1] } = eqn_elecx_D
];

spec_evol_elecy  := [
  { i=[3,Nx-2,1], j=[3,Ny-2,1], k=[3,Nz-2,1] } = eqn_elecy_D
];

spec_evol_elecz  := [
  { i=[3,Nx-2,1], j=[3,Ny-2,1], k=[3,Nz-2,1] } = eqn_elecz_D
];

spec_evol_magx  := [
  { i=[3,Nx-2,1], j=[3,Ny-2,1], k=[3,Nz-2,1] } = eqn_magx_D
];

spec_evol_magy  := [
  { i=[3,Nx-2,1], j=[3,Ny-2,1], k=[3,Nz-2,1] } = eqn_magy_D
];

spec_evol_magz  := [
  { i=[3,Nx-2,1], j=[3,Ny-2,1], k=[3,Nz-2,1] } = eqn_magz_D
];

spec_evol_eem  := [
  { i=[3,Nx-2,1], j=[3,Ny-2,1], k=[3,Nz-2,1] } = eqn_eem_D
];

spec_evol_gaussE  := [
  { i=[3,Nx-2,1], j=[3,Ny-2,1], k=[3,Nz-2,1] } = eqn_gaussE_D
];

spec_evol_gaussB  := [
  { i=[3,Nx-2,1], j=[3,Ny-2,1], k=[3,Nz-2,1] } = eqn_gaussB_D
];

A_Gen_Solve_Code(spec_evol_modphi, {modphi(n+1,i,j,k)}, input="d", proc_name="u_modphi");
A_Gen_Solve_Code(spec_evol_qden,   {qden(n+1,i,j,k)},   input="d", proc_name="u_qden");
A_Gen_Solve_Code(spec_evol_eden_poly, {eden(n+1,i,j,k)}, input="d", proc_name="u_eden_poly");
A_Gen_Solve_Code(spec_evol_eden_log,  {eden(n+1,i,j,k)}, input="d", proc_name="u_eden_log");
A_Gen_Solve_Code(spec_evol_jdenx,  {jdenx(n+1,i,j,k)},   input="d", proc_name="u_jdenx");
A_Gen_Solve_Code(spec_evol_lorenz, {lorenz(n+1,i,j,k)}, input="d", proc_name="u_lorenz");
A_Gen_Solve_Code(spec_evol_elecx,  {elecx(n+1,i,j,k)},  input="d", proc_name="u_elecx");
A_Gen_Solve_Code(spec_evol_elecy,  {elecy(n+1,i,j,k)},  input="d", proc_name="u_elecy");
A_Gen_Solve_Code(spec_evol_elecz,  {elecz(n+1,i,j,k)},  input="d", proc_name="u_elecz");
A_Gen_Solve_Code(spec_evol_magx,   {magx(n+1,i,j,k)},   input="d", proc_name="u_magx");
A_Gen_Solve_Code(spec_evol_magy,   {magy(n+1,i,j,k)},   input="d", proc_name="u_magy");
A_Gen_Solve_Code(spec_evol_magz,   {magz(n+1,i,j,k)},   input="d", proc_name="u_magz");
A_Gen_Solve_Code(spec_evol_eem,    {eem(n+1,i,j,k)},    input="d", proc_name="u_eem");
A_Gen_Solve_Code(spec_evol_gaussE, {gaussE(n+1,i,j,k)}, input="d", proc_name="u_gausse");
A_Gen_Solve_Code(spec_evol_gaussB, {gaussB(n+1,i,j,k)}, input="d", proc_name="u_gaussb");

# initial data
init_modphi := rhs(eqn_modphi);
init_qden   := rhs(eqn_qden);
init_eden_poly := rhs(eqn_eden_poly);
init_eden_log  := rhs(eqn_eden_log);
init_jdenx   := rhs(eqn_jdenx);
init_lorenz := rhs(eqn_lorenz);
init_elecx  := rhs(eqn_elecx);
init_elecy  := rhs(eqn_elecy);
init_elecz  := rhs(eqn_elecz);
init_magx   := rhs(eqn_magx);
init_magy   := rhs(eqn_magy);
init_magz   := rhs(eqn_magz);
init_eem    := rhs(eqn_eem);
init_gaussE := rhs(eqn_gaussE);
init_gaussB := rhs(eqn_gaussB);
Gen_Eval_Code(init_modphi, input="c", proc_name="init_modphi");
Gen_Eval_Code(init_qden,   input="c", proc_name="init_qden");
Gen_Eval_Code(init_eden_poly, input="c", proc_name="init_eden_poly");
Gen_Eval_Code(init_eden_log,  input="c", proc_name="init_eden_log");
Gen_Eval_Code(init_jdenx,  input="c", proc_name="init_jdenx");
Gen_Eval_Code(init_lorenz, input="c", proc_name="init_lorenz");
Gen_Eval_Code(init_elecx,  input="c", proc_name="init_elecx");
Gen_Eval_Code(init_elecy,  input="c", proc_name="init_elecy");
Gen_Eval_Code(init_elecz,  input="c", proc_name="init_elecz");
Gen_Eval_Code(init_magx,   input="c", proc_name="init_magx");
Gen_Eval_Code(init_magy,   input="c", proc_name="init_magy");
Gen_Eval_Code(init_magz,   input="c", proc_name="init_magz");
Gen_Eval_Code(init_eem,    input="c", proc_name="init_eem");
Gen_Eval_Code(init_gaussE, input="c", proc_name="init_gausse");
Gen_Eval_Code(init_gaussB, input="c", proc_name="init_gaussb");

init_chi := eqn_chi_init;
init_xi  := eqn_xi_init;
Gen_Eval_Code(init_chi, input="c", proc_name="init_chi");
Gen_Eval_Code(init_xi,  input="c", proc_name="init_xi");

init_Xx    := compactificationX + myzero*y*z;
init_Yy    := compactificationY + myzero*x*z;
init_Zz    := compactificationZ + myzero*x*y;
init_dxdX  := simplify(subs(X = compactificationX, diff(compactificationX_i, X))) + myzero*y*z;
init_dxdX2 := simplify(subs(X = compactificationX, diff(compactificationX_i, X, X))) + myzero*y*z;
init_dydY  := simplify(subs(Y = compactificationY, diff(compactificationY_i, Y))) + myzero*x*z;
init_dydY2 := simplify(subs(Y = compactificationY, diff(compactificationY_i, Y, Y))) + myzero*x*z;
init_dzdZ  := simplify(subs(Z = compactificationZ, diff(compactificationZ_i, Z))) + myzero*x*y;
init_dzdZ2 := simplify(subs(Z = compactificationZ, diff(compactificationZ_i, Z, Z))) + myzero*x*y;
init_jac   := simplify(diff(compactificationX,x)*diff(compactificationY,y)*diff(compactificationZ,z));

Gen_Eval_Code(init_Xx,    input="c", proc_name="init_xx");
Gen_Eval_Code(init_dxdX,  input="c", proc_name="init_dxdx");
Gen_Eval_Code(init_dxdX2, input="c", proc_name="init_dxdx2");
Gen_Eval_Code(init_Yy,    input="c", proc_name="init_yy");
Gen_Eval_Code(init_dydY,  input="c", proc_name="init_dydy");
Gen_Eval_Code(init_dydY2, input="c", proc_name="init_dydy2");
Gen_Eval_Code(init_Zz,    input="c", proc_name="init_zz");
Gen_Eval_Code(init_dzdZ,  input="c", proc_name="init_dzdz");
Gen_Eval_Code(init_dzdZ2, input="c", proc_name="init_dzdz2");
Gen_Eval_Code(init_jac,   input="c", proc_name="init_jac");
