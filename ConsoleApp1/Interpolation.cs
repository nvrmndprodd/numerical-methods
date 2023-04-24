using System.Drawing;

namespace ConsoleApp1;

public class Interpolation
{
    public Interpolation()
    {
        Execute();
    }

    private void Execute()
    {
        var a = -1;
        var b = 4;
        var n = 8;
        
        var x = ParticularNode(a, b, 100);
        var yOrigin = GetFunctionValues(x, F);

        var (randomX, randomY) = RandomNodes(a, b, n);
        var yRandomLagrange = LagrangeValues(randomX, randomY, x);

        var (particularX, particularY) = ParticularNodes(a, b, n);
        var yParticularLagrange = LagrangeValues(particularX, particularY, x);

        var deltaOriginAndRandomValues = new double[100];
        var deltaOriginAndParticularValues = new double[100];
        for (var i = 0; i < x.Length; ++i)
        {
            deltaOriginAndRandomValues[i] = Double.Abs(yOrigin[i] - yRandomLagrange[i]);
            deltaOriginAndParticularValues[i] = Double.Abs(yOrigin[i] - yParticularLagrange[i]);
        }

        var plt = new ScottPlot.Plot(1000, 1000);
        plt.AddScatter(x, yOrigin, color: Color.Cyan, label: "origin");
        plt.AddScatter(x, yRandomLagrange, color: Color.Chartreuse, label: "random");
        plt.AddScatter(x, yParticularLagrange, color: Color.Black, label: "particular");
        plt.Grid(true);
        plt.Legend();
        plt.SaveFig("Interpolation.png");
        
        var plt1 = new ScottPlot.Plot(600, 200);
        plt1.AddScatter(x, deltaOriginAndRandomValues, Color.Black, label: "F - random");
        plt1.AddScatter(x, deltaOriginAndParticularValues, Color.HotPink, label: "F - particular");
        plt1.Grid(true);
        plt1.Legend();
        plt1.SaveFig("InterpolationDelta.png");
    }

    private double[] GetFunctionValues(double[] x, Func<double, double> func)
    {
        var values = new double[x.Length];
        
        for (int i = 0; i < x.Length; i++) 
            values[i] = func(x[i]);
        
        return values;
    }

    private double[] LagrangeValues(double[] x, double[] y, double[] xValues)
    {
        var yValues = new double[xValues.Length];

        for (var i = 0; i < xValues.Length; i++)
        {
            yValues[i] = LagrangeValue(x, y, xValues[i]);
        }

        return yValues;
    }

    private double LagrangeValue(double[] x, double[] y, double xValue)
    {
        var yValue = 0d;
        var products = y[0];

        for (var i = 0; i < x.Length; ++i)
        {
            products = y[i];
            
            for (int j = 0; j < x.Length; j++)
            {
                if (i != j) 
                    products *= (xValue - x[j]) / (x[i] - x[j]);
            }

            yValue += products;
        }

        return yValue;
    }

    private Tuple<double[], double[]> ParticularNodes(double a, double b, int n)
    {
        var step = (b - a) / (n - 1);
        var x = new double[n];
        var y = new double[n];

        for (var i = 0; i < n; ++i)
        {
            x[i] = a + i * step;
            y[i] = F(x[i]);
        }

        return new Tuple<double[], double[]>(x, y);
    }

    private double[] ParticularNode(double a, double b, int n)
    {
        var step = (b - a) / n;
        var x = new double[n];

        for (var i = 0; i < n; ++i)
        {
            x[i] = a + i * step;
        }

        return x;
    }

    private Tuple<double[], double[]> RandomNodes(double a, double b, int n)
    {
        var x = new double[n];
        var y = new double[n];

        for (var i = 0; i < n; ++i)
        {
            x[i] = 0.5d * ((b - a) * Math.Cos((double)(1 + 2 * i) / (2 * n) * Math.PI) + b + a );
            y[i] = F(x[i]);
        }

        return new Tuple<double[], double[]>(x, y);
    }

    private double F2(double x) => 
        Math.Abs(x) * F(x);
    private double F(double x) => 
        Math.Pow(x, 3) - Math.Exp(x) + 1;
}