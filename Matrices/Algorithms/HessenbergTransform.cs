using MatricesAlgorithms.MatrixObjects;

namespace MatricesAlgorithms.Algorithms;

internal class HessenbergTransform
{
    protected Matrix Matrix { get; set; }

    internal HessenbergTransform(Matrix matrix)
    {
        if(!matrix.IsQuadratic()) throw new Exception("cannot calculate HessenbergTransform of non quadratic matrix.");
        Matrix = CalculateHessenbergTransform(matrix);
    }

    private Matrix CalculateHessenbergTransform(Matrix m)
    {
        int n = m.GetRows();
        Complex[,] hbtransform = new Complex[n, n];

        for(int i = 1; i <= n; i++) for(int j = 1; j <= n; j++) hbtransform[i - 1, j - 1] = m.GetValue(i, j);

        for(int j = 1; j <= n - 2; j++)
        {
            for(int i = j + 2; i <= n; i++)
            {
                Complex aij = hbtransform[i - 1, j - 1];
                Complex aj1j = hbtransform[j, j - 1];
                if(aij != Complex.Zero)
                {
                    Complex w, c, s;
                    if(aj1j.Absolute() < double.Epsilon * aij.Absolute())
                    {
                        w = -m.GetValue(i, j);
                        c = new();
                        s = new(1);
                    }
                    else
                    {
                        Complex det = aj1j * aj1j + aij * aij;
                        w = aj1j.Sign() * det.SquareRoot();
                        c = aj1j / w;
                        s = -(aij / w);
                    }

                    hbtransform[j, j - 1] = w;
                    hbtransform[i - 1, j - 1] = new();

                    for(int k = j + 1; k <= n; k++)
                    {
                        Complex h = c * hbtransform[j, k - 1] - s * hbtransform[i - 1, k - 1];
                        hbtransform[i - 1, k - 1] = s * hbtransform[j, k - 1] + c * hbtransform[i - 1, k - 1];
                        hbtransform[j, k - 1] = h;
                    }
                    for(int k = 1; k <= n; k++)
                    {
                        Complex h = c * hbtransform[k - 1, j] - s * hbtransform[k - 1, i - 1];
                        hbtransform[k - 1, i - 1] = s * hbtransform[k - 1, j] + c * hbtransform[k - 1, i - 1];
                        hbtransform[k - 1, j] = h;
                    }
                }
            }
        }
        return new(hbtransform);
    }

    internal Matrix GetResult() => Matrix;
}