using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace ConsoleApp1;

public class MinimizationOfQuadraticFunction
{
    private const int N = 16;
    private const double epsilon = 10e-6;

    private Matrix<double> A;
    private Matrix<double> B;

    public MinimizationOfQuadraticFunction()
    {
        A = DenseMatrix.OfArray(new[,]
            {
                { 4, 1, 1 },
                { 1, 6 + 0.2 * N, -1 },
                { 1, -1, 8 + 0.2 * N }
            });

        B = DenseMatrix.OfArray(new double[,]
            {
                { 1 },
                { -2 },
                { 3 }
            });
        
        GradientDescent(A, B);
        Console.WriteLine("---------------------------------------------------------------------------------------------");
        Console.WriteLine("Coordinate descent");
        CoordinateDescent(A, B);
    }

    private void GradientDescent(Matrix<double> a, Matrix<double> b)
    {
        Matrix<double> x0 = DenseMatrix.OfArray(new double[,]
        {
            { 1 },
            { 0 },
            { 0 }
        });
    
        Matrix<double> x1 = x0;
        var i = 0;
        
        while (true)
        {
            ++i;
            var q = A * x0 + B;

            var m = -(q.Transpose() * q) * (q.Transpose() * (A * q)).Inverse();
            x1 = x0 + m[0,0] * q;
            
            if (double.Abs(F(x1) - F(x0)) <= epsilon)
                break;
            
            x0 = x1;
        }

        Console.WriteLine($"Iterations: {i}");
        Console.WriteLine($"Extremum Point: {x1}");
        Console.WriteLine($"Value at the extremum point: {F(x1)}");
    }

    private void CoordinateDescent(Matrix<double> a, Matrix<double> b)
    {
        Matrix<double> x0 = DenseMatrix.OfArray(new double[,]
        {
            { 0 },
            { 0 },
            { 0 }
        });
        
        Matrix<double> e = DenseMatrix.OfArray(new double[,]
        {
            { 1 },
            { 0 },
            { 0 }
        });

        Matrix<double> x1 = x0;
        var i = 0;

        while (true)
        {
            var m = -(e.Transpose() * (A * x0 + B)) * (e.Transpose() * (A * e)).Inverse();
            x1 = x0 + m[0, 0] * e;
            
            if ((x1 - x0).L1Norm() <= epsilon)
                break;

            x0 = x1;
            e[i % 3, 0] = 0;
            ++i;
            e[i % 3, 0] = 1;
        }
        
        Console.WriteLine($"Iterations: {i}");
        Console.WriteLine($"Extremum Point: {x1}");
        Console.WriteLine($"Value at the extremum point: {F(x1)}");
    }

    private double F(Matrix<double> x) =>
        2 * Math.Pow(x[0, 0], 2) + (3 + 0.1 * N) * Math.Pow(x[1, 0], 2) + (4 + 0.1 * N) * Math.Pow(x[2, 0], 2) +
        x[0, 0] * x[1, 0] - x[1, 0] * x[2, 0] + x[0, 0] * x[1, 0] + x[0, 0] - 2 * x[1, 0] + 3 * x[2, 0] + N;
}