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
                  
                      R[i] = SurfaceTension.SurfTens(T, C[i].Sigma);
              end for;
              ST = sum(x[:] .* R[:]);
      end Test;
  
  
end SurfaceTension;
