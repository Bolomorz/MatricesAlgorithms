using MatricesAlgorithms.MatrixObjects;

namespace MatricesAlgorithms.Algorithms;

internal class BairstowAlgorithm
{
    private class QuadraticRoots
    {
        internal required Complex Root1 { get; set; }
        internal required Complex Root2 { get; set; }
    }
    private class QuadFactor
    {
        internal required Complex U { get; set; }
        internal required Complex V { get; set; }
        internal required Complex[] Result { get; set; }
    }
    private const int ITERMAX = 1000;
    protected Complex[] EigenValues { get; set; }

    internal BairstowAlgorithm(CharacteristicPolynomial polynomial)
    {
        EigenValues = Bairstow(polynomial.GetResult()).ToArray();
    }

    private List<Complex> Bairstow(Complex[] A)
    {
        List<Complex> roots = new();

        while(A.Length > 2)
        {
            int n = A.Length;
            Complex u = A[n - 2] / A[n - 1];
            Complex v = A[n - 3] / A[n - 1];

            QuadFactor qf = SolveQuadFactor(A, u, v);
            QuadraticRoots qr = SolveQuadraticEquation(new(1), qf.U, qf.V);

            roots.Add(qr.Root1);
            roots.Add(qr.Root2);

            A = qf.Result;
        }

        if(A.Length == 2) roots.Add(-(A[0] / A[1]));

        return roots;
    }
    private QuadFactor SolveQuadFactor(Complex[] A, Complex U, Complex V)
    {
        int n = A.Length - 1;
        Complex u = new(U.Re, U.Im);
        Complex v = new(V.Re, V.Im);
        Complex c = new(1);
        Complex d = new(1);
        int iteration = 0;

        Complex[] B = new Complex[n + 1];
        while( c * c + d * d >= Complex.Tol && iteration < ITERMAX )
        {
            Complex[] F = new Complex[n + 1];
            B = new Complex[n + 1];
            for(int i = 0; i < n + 1; i++)
            {
                F[i] = new();
                B[i] = new();
            }

            for(int i = n - 2; i >= 0; i--)
            {
                B[i] = A[i + 2] - u * B[i + 1] - v * B[i + 2];
                F[i] = B[i + 2] - u * F[i + 1] - v * F[i + 2];
            }

            c = A[1] - u * B[0] - v * B[1];
            d = A[0] - v * B[0];
            Complex g = B[1] - u * F[0] - v * F[1];
            Complex h = B[0] - v * F[0];
            Complex det = v * g * g + h * (h - u * g);
            u -= (-h * c + g * d) / det;
            v -= (-g * v * c + (g * u - h) * d) / det;
            iteration++;
        }

        Complex[] result = new Complex[n - 1];
        for(int i = 0; i < n; i++) result[i] = B[i];

        return new(){ Result = result, U = u, V = v};
    }
    private QuadraticRoots SolveQuadraticEquation(Complex A, Complex B, Complex C)
    {
        Complex det = B * B - 4 * A * C;
        return new(){ Root1 = (-B + det.SquareRoot()) / (2.0 * A), Root2 = -(B / A + (-B + det.SquareRoot()) / (2.0 * A))};
    }

    internal Complex[] GetResult() => EigenValues;
}