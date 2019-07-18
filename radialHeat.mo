within ;
model radialHeat
  import SI = Modelica.SIunits;
  import Modelica.Constants.pi;
  import Modelica.Math.log;

  parameter SI.Radius r1 = 0.025;
  parameter SI.Radius r2 = 0.035;
  parameter SI.Radius r3 = 0.06;
  constant SI.ThermalConductivity k1 = 45;
  constant SI.ThermalConductivity k2 = 0.087;

  constant SI.CoefficientOfHeatTransfer h1 = 11350;
  SI.CoefficientOfHeatTransfer h2;

  parameter SI.Temperature T1 = 418;
  SI.Temperature T2;
  SI.Temperature T3;
  SI.Temperature T4;
  SI.Temperature T5;
  SI.HeatFlowRate q;

  parameter SI.Time valTime[:]=  0:3600:86400;
  parameter SI.Temperature valTemp[:] = {291.5,293.0,295.5,300.0,302.5,305.5,
  306.0,306.5,308.5,310.0,312.5,313.5,313.5,312.5,310.0,309.0,
  305.5,300.5,298.5,293.5,290.5,286.0,280.5,275.5,273.0};

  function LinearInterpolation
    input Real x;
    input Real tableX[:];
    input Real tableY[:];
    output Real y;

  protected
      Integer n;
      Real slope;

  algorithm
    n := size(tableX, 1);
    assert(size(tableX, 1) == size(tableY, 1), "Error: tableX and tableY with different size");
    assert(x >= tableX[1] and x <= tableX[n], "Error: independent variable is out of range");
    for i in 1:n - 1 loop
      if x >= tableX[i] and x <= tableX[i + 1] then
        slope := (tableY[i+1] - tableY[i])/(tableX[i+1] - tableX[i]);
        y := tableY[i] + slope*(x - tableX[i]);
      end if;
    end for;
  end LinearInterpolation;

equation
  h2 = 1.32*(abs(T4 - T5)/(2*r3))^0.25;
  T1 - T2 = q/(h1*2*pi*r1);
  T2 - T3 = q*log(r2/r1)/(2*pi*k1);
  T3 - T4 = q*log(r3/r2)/(2*pi*k2);
  T4 - T5 = q/(h2*2*pi*r3);
  // T5 = 291.5 - 8.5*sin(2*pi*time/86400);
  T5 = LinearInterpolation(x = time, tableX = valTime, tableY = valTemp);
end radialHeat;
