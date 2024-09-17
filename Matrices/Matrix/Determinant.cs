namespace MatricesAlgorithms.MatrixObjects;

internal class Determinant
{
    private class Decomposition
    {
        internal required bool Success { get; set; }
        internal required Complex[,] Decompose { get; set; }
        internal required int[] P { get; set; }
    }
    internal Complex Value { get; }

    internal Determinant(Matrix matrix)
    {
        if(!matrix.IsQuadratic()) throw new Exception("cannot calculate determinant of non quadratric matrix.");
        var decomposition = LUPDecompose(matrix.GetValues(), matrix.GetRows(), Complex.Tol);
        Value = decomposition.Success ? LUPDeterminant(decomposition, matrix.GetRows()) : CalculateDeterminant(matrix.GetValues(), matrix.GetRows());
    }

    #region LUP Algorithm - Faster
    private Complex LUPDeterminant(Decomposition decomposition, int n)
    {
        Complex d = decomposition.Decompose[0, 0];
        for(int i = 1; i < n; i++) d *= decomposition.Decompose[i, i];
        return (decomposition.P[n] - n) % 2 == 0 ? d : -d;
    }
    private Decomposition LUPDecompose(Complex[,] A, int n, Complex Tol)
    {
        int i, j, k, imax;
        Complex maxA, absA;
        Complex[] ptr = new Complex[n];

        int[] P = new int[n + 1];
        Complex[,] decompose = new Complex[n, n];

        for(i = 0; i < n; i++) for(j = 0; j < n; j++) decompose[i, j] = A[i, j];

        for(i = 0; i <= n; i++) P[i] = i;

        for(i = 0; i < n; i++)
        {
            maxA = new Complex();
            imax = i;

            for(k = i; k < n; k++)
            {
                absA = new Complex(decompose[k, i].Absolute());
                if(absA > maxA)
                {
                    maxA = absA;
                    imax = k;
                }
            }

            if(maxA < Tol) return new(){ Decompose = decompose, P = P, Success = false };

            if(imax != i)
            {
                j = P[i];
                P[i] = P[imax];
                P[imax] = j;

                for(k = 0; k < n; k++) ptr[k] = decompose[i, k];
                for(k = 0; k < n; k++) decompose[i, k] = decompose[imax, k];
                for(k = 0; k < n; k++) decompose[imax, k] = ptr[k];

                P[n]++;
            }

            for(j = i+1; j < n; j++)
            {
                decompose [j, i] /= decompose[i, i];
                for(k = i+1; k < n; k++) decompose[j, k] -= decompose[j, i] * decompose[i, k];
            }
        }
        return new(){ Decompose = decompose, P = P, Success = true };
    }
    #endregion

    #region SubDeterminant Algorithm - Very Slow
    private Complex[,] SubMatrix(Complex[,] m, int i, int k, int n)
    {
        Complex[,] sub = new Complex[n - 1, n - 1];
        int subrow = 0;
        int subcol;
        i--; k--;
        for(int row = 0; row < n; row++)
        {
            if(row != i)
            {
                subcol = 0;
                for(int col = 0; col < n; col++) if(col != k) sub[subrow, subcol++] = m[row, col];
                subrow++;
            }
        }
        return sub;
    }
    private Complex CalculateDeterminant(Complex[,] m, int n)
    {
        if(n == 2) return m[0,0] * m[1,1] - m[0,1] * m[1, 0];
        Complex sum = new();
        List<Complex> dets = new();
        for(int i = 1; i <= n; i++) dets.Add(CalculateDeterminant(SubMatrix(m, 1, i, n), n - 1));
        for(int col = 1; col <= n; col++)
        {
            Complex sub = dets[col - 1];
            sum += (col + 1) % 2 == 0 ? m[0, col - 1] * sub : -1 * m[0, col - 1] * sub;
        }
        return sum;
    }
    #endregion
}