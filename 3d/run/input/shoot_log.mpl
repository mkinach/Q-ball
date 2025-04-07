with(ArrayTools):
with(CurveFitting):

interface(displayprecision=-1):
CWD := currentdir();
Digits := 50;  # software precision

# parameters to be used for shooting
etilde  := 1.1;
gtilde0 := 1.0;
low0    := 0.6;
high0   := 0.7;
ftilde0 := (low0 + high0) / 2;
mu      := 1.0;
beta    := 1.0;
V       := -mu^2*modphi^2*log(beta^2*modphi^2);
numiters   := 75:
xmax       := 25:
xmax_f     := 10000;
mystepsize := 0.01:
mystepsize_f := mystepsize;

numelements := ceil(xmax/mystepsize);
xlist := [seq(i*mystepsize, i=0..numelements)]:

print(0, "Trying initial ftilde0:", ftilde0);

low  := low0:
high := high0:

for j from 1 to numiters do
  exitvar := 0:

  ics := {ftilde(0) = ftilde0,
          gtilde(0) = gtilde0,
          utilde(0) = 0,
          ntilde(0) = 0}:

  sys_ode := {
    diff(ftilde(x), x) = utilde(x),
    diff(gtilde(x), x) = ntilde(x),
    diff(utilde(x), x) = piecewise(x = 0, (1/3)*(-gtilde(x)^2*ftilde(x) - ftilde(x)*ln(ftilde(x)^2) - ftilde(x)),
                                         -(2/x)*utilde(x) - gtilde(x)^2*ftilde(x) - ftilde(x)*ln(ftilde(x)^2) - ftilde(x)),
    diff(ntilde(x), x) = piecewise(x = 0, (2/3)*etilde^2*gtilde(x)*ftilde(x)^2,
                                         -(2/x)*ntilde(x) + 2*etilde^2*ftilde(x)^2*gtilde(x))
  }:

  sol := dsolve( [ sys_ode[], ics[] ],
                   numeric,
                   method   = classical[rk4],
                   output   = Array(xlist),
                   stepsize = mystepsize,
                   maxfun   = 0
               ):

  sol(1):
  x_vals := sol(2)[..,1]:
  ftilde_vals := sol(2)[..,2]:
  gtilde_vals := sol(2)[..,3]:
  ntilde_vals := sol(2)[..,4]:
  utilde_vals := sol(2)[..,5]:

  # works on some parts of the branch
  if utilde_vals(2) > 0 then  # stepsize must be sufficiently small!
    count := 0:
    for z from 1 to numelements do
      if utilde_vals(z)*utilde_vals(z+1) < 0 then
        count := count + 1:
        if count = 2 then
          if ftilde_vals(z) > 0 then
            print(j, "Last val too high. Trying new ftilde0:", ftilde0);
            high := (low + high) / 2:
            ftilde0 := (low + high)/2:
            break:
          else
            print(j, "Last val too high. Trying new ftilde0:", ftilde0);
            low := (low + high) / 2:
            ftilde0 := (low+high) / 2:
            break:
          end if:
        end if:
      end if:
    end do:
    next:
  end if:

  # works on other parts of the branch
  for z from 1 to numelements + 1 do
    if utilde_vals(z) > 0 then
      if ftilde_vals(z) > 0 then
        low := (low + high) / 2:
        print(j, "Last val too low. Trying new ftilde0:", (low + high) / 2);
      else
        high := (low + high) / 2:
        print(j, "Last val too high. Trying new ftilde0:", (low + high) / 2);
      end if;
      ftilde0 := (low + high) / 2;
      break;
    elif z = numelements + 1 then
      print("INTEGRATION REGION NOT LARGE ENOUGH TO EVER GET utilde_vals > 0.");
      exitvar := 1;
    end if;
  end do:

  if exitvar = 1 then
    break;
  end if:

end do:

# apply a linear fit near the point where logdiff changes sign or becomes imaginary
ftilde_max := 0:
for j from 1 to numelems(ftilde_vals) do
  if ftilde_vals(j) > ftilde_max then
    ftilde_max := ftilde_vals(j);
  end if;
end do:
logpt_max := -1:
logpt_min := -1:
logdiff_last := log(ftilde_vals(2)) - log(ftilde_vals(1)):
for i from 2 to numelems(ftilde_vals) - 1 do
  logdiff := log(ftilde_vals(i+1)) - log(ftilde_vals(i));
  if (logdiff_last < 0 and logdiff > 0) or (Im(log(ftilde_vals(i+1))) <> 0) then
    if ftilde_vals(i)/ftilde_max < 0.01 then
      logpt_max := x_vals(i)*0.975;
      logpt_min := x_vals(i)*0.95;
      break;
    end if;
  end if;
  logdiff_last := logdiff;
end do:

# check that the relative slope is small near the fit point
relslope := abs(utilde_vals(i)/ftilde_max);
print("relsope:", relslope);
if relslope > 1e-9 then
  print("NO Q-BALL SOLUTION DETECTED!");
end if;

# round to the nearest factor of mystepsize
logpt_min := round(logpt_min/mystepsize)*mystepsize;
logpt_max := round(logpt_max/mystepsize)*mystepsize;

# perform a linear fit
fit := LeastSquares(x_vals[ceil(logpt_min/mystepsize)..ceil(logpt_max/mystepsize)], map(log, ftilde_vals[ceil(logpt_min/mystepsize)..ceil(logpt_max/mystepsize)]), x);
fit_slope := diff(%, x);
fit_lower := LeastSquares(x_vals[ceil(0.75*logpt_min/mystepsize)..ceil(0.75*logpt_max/mystepsize)], map(log, x_vals[ceil(0.75*logpt_min/mystepsize)..ceil(0.75*logpt_max/mystepsize)].~ftilde_vals[ceil(0.75*logpt_min/mystepsize)..ceil(0.75*logpt_max/mystepsize)]), x);
fit_lower_slope := diff(%, x);
linearity_test := abs(fit_lower_slope/fit_slope);
if (linearity_test < 0.75) or (linearity_test > 1.25) then
  print("NO Q-BALL SOLUTION DETECTED!");
end if;

# check scalar field value at the fit point
ftilde_vals(ceil(logpt_min/mystepsize));
expfit := exp(fit)/x;
d_expfit := diff(expfit, x):
expfit_f := unapply(expfit, x):
d_expfit_f := unapply(d_expfit, x):

# re-solve the equations up to some large radius
fit_start := round(logpt_max/mystepsize_f)*mystepsize_f;
fit_start_pt := round(fit_start/mystepsize_f);
numelements_f := ceil((xmax_f - fit_start)/mystepsize_f);
xlist_f := [seq(fit_start+i*mystepsize_f, i=0..numelements_f)]:

ftilde_f := proc(x);
  if not type(evalf(x), 'numeric') then
    'procname'(x);
  else
    expfit_f(x);
  end if;
end proc;

ics := {
  gtilde_f(fit_start) = gtilde_vals(fit_start_pt+1),
  utilde_f(fit_start) = utilde_vals(fit_start_pt+1),
  ntilde_f(fit_start) = ntilde_vals(fit_start_pt+1)
}:

sys_ode := {
  diff(gtilde_f(x), x) = ntilde_f(x),
  diff(utilde_f(x), x) = -(2/x)*utilde_f(x) - gtilde_f(x)^2*ftilde_f(x) - ftilde_f(x)*ln(ftilde_f(x)^2) - ftilde_f(x),
  diff(ntilde_f(x), x) = -(2/x)*ntilde_f(x) + 2*etilde^2*ftilde_f(x)^2*gtilde_f(x)
}:

sol := dsolve( [ sys_ode[], ics[] ],
                 numeric,
                 method = classical[rk4],
                 output = Array(xlist_f),
                 stepsize = mystepsize_f,
                 maxfun = 0,
                 known = ftilde_f
             ):

sol(1):
x_vals_f := sol(2)[..,1]:
gtilde_vals_f := sol(2)[..,2]:
ntilde_vals_f := sol(2)[..,3]:
utilde_vals_f := sol(2)[..,4]:
xvals_all       := Concatenate(1, x_vals[1..fit_start_pt], x_vals_f):
ftilde_vals_all := Concatenate(1, ftilde_vals[1..fit_start_pt], map(ftilde_f, x_vals_f)):
gtilde_vals_all := Concatenate(1, gtilde_vals[1..fit_start_pt], gtilde_vals_f):

# extract the value of omega from the asymptotic value of gtilde
numpts := Size(ftilde_vals_all)[1]:
x_lastfew := xvals_all[numpts-100..numpts]:
gtilde_lastfew := gtilde_vals_all[numpts-100..numpts]:
gtilde_fit := LeastSquares(x_lastfew, gtilde_lastfew, v, curve=a+b/v):
wtilde := gtilde_vals_all[Size(gtilde_vals_all)[1]];  # for comparison purposes
wtilde := limit(gtilde_fit, v=infinity);

# convert gtilde into the physically-relevant quantity A0
if etilde = 0 then
  A0tilde := 0.0*~gtilde_vals_all:
else
  A0tilde := (-wtilde +~ gtilde_vals_all):
end if:
A0tilde_lastfew := A0tilde[numpts-100..numpts]:
A0tilde_fit := LeastSquares(x_lastfew, A0tilde_lastfew, v, curve=a+b/v);
A0tilde_constant := limit(A0tilde_fit, v=infinity);
A0tilde := A0tilde -~ A0tilde_constant:

# convert the rest of the quantities from tilde to non-tilde
r := xvals_all/~mu:
e := etilde*mu*beta;
w := mu*wtilde;
f := ftilde_vals_all/~beta:
A0 := (-mu/etilde)*~A0tilde:

# compute the total energy of the Q-ball
A0_diff := 0*~A0:
f_diff  := 0*~f:
A0_diff(1) := (A0(2) - A0(1))/mystepsize:
f_diff(1)  := (f(2) - f(1))/mystepsize:
for ii from 2 to Size(A0)[1] - 1 do
  A0_diff(ii) := (A0(ii+1) - A0(ii-1))/(2*mystepsize);
  f_diff(ii)  := (f(ii+1) - f(ii-1))/(2*mystepsize);
end do:
maxpoints := Size(A0)[1]:

V_vals_all := -mu^2*~ftilde_vals_all^~2*~map(log, beta^2*~ftilde_vals_all^~2):

Eden_E := (1/2)*~A0_diff^~2:
Eden_G := f_diff^~2:
Eden_P := V_vals_all:
Eden_T := gtilde_vals_all^~2*~ftilde_vals_all^~2:

Etot_E := 0:
Etot_G := 0:
Etot_P := 0:
Etot_T := 0:

for z from 1 to maxpoints - 1 do
  Etot_E := Etot_E + (xvals_all(z+1) - xvals_all(z))*(Eden_E(z+1)*xvals_all(z+1)^2 + Eden_E(z)*xvals_all(z)^2)/2;
  Etot_G := Etot_G + (xvals_all(z+1) - xvals_all(z))*(Eden_G(z+1)*xvals_all(z+1)^2 + Eden_G(z)*xvals_all(z)^2)/2;
  Etot_P := Etot_P + (xvals_all(z+1) - xvals_all(z))*(Eden_P(z+1)*xvals_all(z+1)^2 + Eden_P(z)*xvals_all(z)^2)/2;
  Etot_T := Etot_T + (xvals_all(z+1) - xvals_all(z))*(Eden_T(z+1)*xvals_all(z+1)^2 + Eden_T(z)*xvals_all(z)^2)/2;
end do:

Etot_E := 4*Pi*Etot_E;
Etot_G := 4*Pi*Etot_G;
Etot_P := 4*Pi*Etot_P;
Etot_T := 4*Pi*Etot_T;

Etot := Etot_E + Etot_G + Etot_P + Etot_T;
Eden_vals := Eden_E +~ Eden_G +~ Eden_P +~ Eden_T:

# compute the total Noether charge of the Q-ball
Qden_vals := -2*~ftilde_vals^~2*~gtilde_vals:

Qtot := 0:
intmax := round(logpt_min/mystepsize);

for z from 1 to intmax do
  Qtot := Qtot + (xvals_all(z+1) - xvals_all(z))*(Qden_vals(z+1)*xvals_all(z+1)^2 + Qden_vals(z)*xvals_all(z)^2)/2;
end do:

Qtot := 4*Pi*Qtot;

# write data to file
exportvar := 1:
if exportvar = 1 then
  if mystepsize - mystepsize_f = 0 then
    dataarray := < r | f | A0 >:
    outputname := cat(CWD,"/initdata");
    ExportMatrix(outputname, dataarray, delimiter="\t"):
    paramfilename := cat(CWD, "/initdata.fparam"):
    Digits := 15:
    params:=[ ["e :=", e], ["w :=", w], ["mu :=", mu], ["beta :=", beta] ]:
    writedata(paramfilename, params, [string, float]):
    print("Data exported");
  else
    print("Error: mystepsize != mystepsize_f (needed for interpolation)")
  end if:
end if:
