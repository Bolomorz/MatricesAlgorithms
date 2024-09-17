using MatricesAlgorithms.MatrixObjects;

namespace MatricesAlgorithms.Algorithms;

internal class QRTransform
{
    protected Matrix Matrix { get; set; }

    internal QRTransform(HessenbergTransform hessenberg)
    {
        Matrix hbt = hessenberg.GetResult();
        throw new NotImplementedException();
    }

    internal Matrix GetResult() => Matrix; 
}