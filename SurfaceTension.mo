package SurfaceTension
  function SurfTens
    extends Modelica.Icons.Function;
    import Simulator.Files.Thermodynamic_Functions.*;
    input Real T, coeff[6];
    output Real sigma;
  algorithm
    sigma := coeff[2] + exp(coeff[3] / T + coeff[4] + coeff[5] * T + coeff[6] * T * T);
  end SurfTens;

      model Test
              import Simulator.*;
              import data = Simulator.Files.ChemsepDatabase;
              parameter data.Benzene ace;
              parameter data.Phenol meth;
              parameter data.Aniline benz;
              parameter Integer Nc = 3;
              parameter data.GeneralProperties C[Nc] = {ace, meth, benz};
              Real ST;
              parameter Real x[Nc] = {0.33, 0.33, 0.34}, P = 101325, T = 300;
              Real R[Nc];
      equation
              for i in 1:Nc loop
                      if T > C[i].SigmaT[1] and T < C[i].SigmaT[2] then
                              R[i] = SurfaceTension.SurfTens(T, C[i].Sigma);
                      else
                              R[i] = SurfaceTension.SurfTensEmp(P, T, C[i].Pc / 100000, C[i].Tc, C[i].Tb) / 1000;
                      end if;
              end for;
              ST = sum(x[:] .* R[:]);
      end Test;
  
  

  function SurfTensEmp
    extends Modelica.Icons.Function;
    import Simulator.Files.Thermodynamic_Functions.*;
    input Real P, T;
    input Real Pc, Tc, Tb;
    output Real sigma;
  protected
    Real Tr, Tbr;
    Real alpha;
  algorithm
    Tr := T / Tc;
    Tbr := Tb / Tc;
    alpha := 0.9076 * (1 + Tbr * log(Pc / 1.01325) / (1 - Tbr));
    sigma := Pc ^ (2 / 3) * Tc ^ (1 / 3) * (0.132 * alpha - 0.279) * (1 - Tr) ^ (11 / 9);
  end SurfTensEmp;
end SurfaceTension;
