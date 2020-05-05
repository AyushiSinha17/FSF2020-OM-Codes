package SurfaceTension
  function SurfTens
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
  
  Real x[Nc]={0.33,0.33,0.34};
  Real R[Nc];
                                  Real Pc[Nc]; 
                                  Real Tc[Nc] ;
                                  Real Tb[Nc]; 
  equation
  
                             for i in 1:Nc loop
                                   Pc[i] = C[i].Pc / 100000;
                                   Tc[i] = C[i].Tc;
                                   Tb[i]= C[i].Tb;    
                                                                                 
                                 R[i] = SurfaceTension.SurfTens(101325, 280, Pc[i], Tc[i], Tb[i]) /1000;
                                                                   end for;
                                                                  
                                                   ST = sum(x[:].*R[:]);
                                                  
                                                                 
  
  
  end Test;
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
end SurfaceTension;
