package JouleThompsonCoefficient
  function JTCoeff
    extends Modelica.Icons.Function;
    import Simulator.Files.Thermodynamic_Functions.*;
    input Real P, T, a;
    input Real Pc, Tc, Cp, rho;
    Real Tr, gamma;
    output Real mu;
  algorithm
    Tr := T / Tc;
    gamma := Cp / (Cp - 8.314);
    if a == 1 then
      mu := 0.0048823 * Tc * (18 / (Tr ^ 2 - 1)) / (Pc * Cp * gamma);
    elseif a == 2 then
      mu := -1 / (1000 * rho * Cp);
    else
      mu := 0;
    end if;
//compound is a gas
  end JTCoeff;

      model Test
              import Simulator.*;
              import data = Simulator.Files.ChemsepDatabase;
              parameter Integer Nc = 2;
              parameter data.Water eth;
              parameter data.Pyridine benz;
              parameter data.GeneralProperties C[Nc] = {eth, benz};
              Real Pc[Nc], Tc[Nc], Tb[Nc], Tm[Nc], rho[Nc], Cp[Nc];
              Integer a[Nc];
              Real P, T;
              Real mu[Nc], x[Nc] = {0.4, 0.6};
              Real JTC;
      equation
              P = 101325;
              T = 298.15;
              for i in 1:Nc loop
                      Pc[i] = C[i].Pc;
                      Tc[i] = C[i].Tc;
                      Tb[i] = C[i].Tb;
                      Tm[i] = C[i].Tm;
                      if T >= Tb[i] then
                              a[i] = 1;
                              rho[i] = 0;
                              Cp[i] = Files.ThermodynamicFunctions.VapCpId(C[i].VapCp, T);
                      elseif T < Tb[i] and T >= Tm[i] then
                              a[i] = 2;
                              Cp[i] = Files.ThermodynamicFunctions.LiqCpId(C[i].LiqCp, T);
                              rho[i] = Files.ThermodynamicFunctions.Dens(C[i].LiqDen, Tc[i], T, P);
                      else
                              a[i] = 3;
                              Cp[i] = 0;
                              rho[i] = 0;
                      end if;
                      mu[i] = JouleThompsonCoefficient.JTCoeff(P, T, a[i], Pc[i], Tc[i], Cp[i], rho[i]);
              end for;
              JTC = sum(x[:] .* mu[:]);
      end Test;
  
  
end JouleThompsonCoefficient;
