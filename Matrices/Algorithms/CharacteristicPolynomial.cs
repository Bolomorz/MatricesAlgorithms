using MatricesAlgorithms.MatrixObjects;

namespace MatricesAlgorithms.Algorithms;

internal class CharacteristicPolynomial
{
    protected Complex[] Polynomial { get; set; }

    public CharacteristicPolynomial(Matrix matrix)
    {
        if(!matrix.IsQuadratic()) throw new Exception("cannot calculate CharacteristicPolynomial of non quadratic matrix.");
        Polynomial = FaddeevLeVerrier(matrix);
    }

    private Complex[] FaddeevLeVerrier(Matrix H)
    {
        int n = H.GetRows();
        Complex[] C = new Complex[n];
        Matrix[] M = new Matrix[n];
        int k = 2;

        C[n - 1] = new(1);
        C[n - 2] = -H.Trace();
        M[0] = new(SpecialQuadratic.Identity, n);

        while(k <= n)
        {
            M[k - 1] = H * M[k - 2] + C[n - k] * new Matrix(SpecialQuadratic.Identity, n); 
            C[n - k - 1] = -(1/k) * (H * M[k - 1]).Trace();
            k++;
        }

        return C;
    }

    internal Complex[] GetResult() => Polynomial;
}