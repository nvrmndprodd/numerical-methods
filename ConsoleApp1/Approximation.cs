using System.Drawing;
using MathNet.Numerics.Integration;
using MathNet.Numerics.LinearAlgebra;

namespace ConsoleApp1;

public class Approximation
{
    private double[] _x;
    private double[] _y;

    public Approximation()
    {
        //_x = new[] { -1, -0.5, 0, 0.5, 1 };
        _x = SplitSegment(-1, 1, 100);
        _y = new double[101];
        
        for (var i = 0; i < _x.Length; i++) 
            _y[i] = F(_x[i]);

        Execute();
    }

    private void Execute()
    {
        var coefficients = Mnk(_x, _y, 3);
        var x = SplitSegment(-1, 1, 100);
        
        var yOrigin = GetFunctionValues(x, F);
        var yApproximation = PolynomialValues(x, coefficients);
        
        var ck = new double[4];
        for (var i = 0; i < 4; ++i) 
            ck[i] = ScalarProductOne(i) / ScalarProductTwo(i);

        var yLegendre = LegendreValues(x, ck);
        
        var plt = new ScottPlot.Plot(600, 600);
        plt.AddScatter(x, yOrigin, color: Color.Cyan, label: "origin");
        plt.AddScatter(x, yApproximation, color: Color.Chartreuse, label: "approximation");
        plt.AddScatter(x, yLegendre, color: Color.Chocolate, label: "legendre");
        plt.Grid(true);
        plt.Legend();
        plt.SaveFig("Approximation.png");

        var deltaOriginAndApproximation = new double[101];
        var deltaOriginAndLegendre = new double[101];
        for (var i = 0; i < 101; i++)
        {
            deltaOriginAndApproximation[i] = double.Abs(yOrigin[i] - yApproximation[i]);
            deltaOriginAndLegendre[i] = double.Abs(yOrigin[i] - yLegendre[i]);
        }
        var plt1 = new ScottPlot.Plot(600, 200);
        plt1.AddScatter(x, deltaOriginAndApproximation, Color.Black, label: "F - approximation");
        plt1.AddScatter(x, deltaOriginAndLegendre, Color.HotPink, label: "F - legendre");
        plt1.Grid(true);
        plt1.Legend();
        plt1.SaveFig("ApproximationDelta.png");
    }

    // вернет коэффициенты для наименьшеквадратичного полинома
    private double[] Mnk(double[] x, double[] y, int n)
    {
        var b = Vector<double>.Build.Dense(n + 1);
        for (int i = 0; i < n + 1; i++)
        {
            var sum = 0d;

            for (var j = 0; j < x.Length; ++j) 
                sum += double.Pow(x[j], i) * y[j];

            b[i] = sum;
        }

        var xPows = new double[2 * n + 1];

        for (int i = 0; i < 2 * n + 1; i++) 
            xPows[i] = x.Sum(xi => double.Pow(xi, i));

        var A = Matrix<double>.Build.Dense(n + 1, n + 1);
        for (var i = 0; i < A.RowCount; i++)
            for (var j = 0; j < A.ColumnCount; j++)
                A[i, j] = xPows[i + j];

        return A.Solve(b).ToArray();
    }

    private double[] PolynomialValues(double[] x, double[] coefficients)
    {
        var y = new double[x.Length];
        var polynomAsString = "";
        var flag = true;

        for (var i = 0; i < x.Length; ++i)
        {
            var r = 0d; // степени помноженные на коэффициенты
            for (var j = 0; j < coefficients.Length; j++)
            {
                if (flag) 
                    polynomAsString += $"x ^ {j} * {coefficients[j]} + ";
                r += double.Pow(x[i], j) * coefficients[j];
            }

            if (flag)
                Console.WriteLine(polynomAsString);
            flag = false;

            y[i] = r;
        }
        return y;
    }

    private double[] LegendreValues(double[] x, double[] coefficients)
    {
        var y = new double[x.Length];
        var polynomAsString = "";
        var flag = true;

        for (var i = 0; i < x.Length; ++i)
        {
            var r = coefficients[0];
            for (var j = 1; j < coefficients.Length; j++)
            {
                if (flag) 
                    polynomAsString += $"x ^ {j} * {coefficients[j]} + ";
                r += LegendrePolynomial(x[i], j) * coefficients[j];
            }

            if (flag)
                Console.WriteLine(polynomAsString);
            flag = false;

            y[i] = r;
        }
        return y;
    }

    private double ScalarProductOne(int n) => 
        SimpsonRule.IntegrateComposite(x => LegendrePolynomial(x, n) * F(x), -1, 1, 20);
    
    private double ScalarProductTwo(int n) => 
        SimpsonRule.IntegrateComposite(x => LegendrePolynomial(x, n) * LegendrePolynomial(x, n), -1, 1, 20);
    
    private double LegendrePolynomial(double x, int n) =>
        n switch
        {
            0 => 1d,
            1 => x,
            2 => (3 * (double.Pow(x, 2)) - 1) / 2,
            3 => (5 * (double.Pow(x, 3)) - 3 * x) / 2,
            _ => 0
        };

    private double[] SplitSegment(double a, double b, int n)
    {
        var step = (b - a) / n;
        var x = new double[n + 1];
        
        for (var i = 0; i < n + 1; ++i) 
            x[i] = a + i * step;

        return x;
    }
    
    private double[] GetFunctionValues(double[] x, Func<double, double> func)
    {
        var values = new double[x.Length];
        
        for (int i = 0; i < x.Length; i++) 
            values[i] = func(x[i]);
        
        return values;
    }

    // private double F(double x) => 
    //     Double.Sqrt(Double.Pow(x, 2)) + Double.Log(Double.Pow(x, 2));

    private double F(double x) => 
         double.Pow(x, 2) * double.Tan(x);
    
    // private double F(double x) 
    //     => double.Pow(x, 3) + double.Exp(x);
}