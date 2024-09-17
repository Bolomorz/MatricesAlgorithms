namespace MatricesAlgorithms.MatrixObjects;

internal class Complex
{
    internal double Re { get; set; }
    internal double Im { get; set; }

    internal static Complex Tol = new(double.Epsilon);
    internal static Complex Zero = new();

    internal Complex(double re, double im)
    {
        Re = re;
        Im = im;
    }
    internal Complex(double re)
    {
        Re = re;
        Im = 0;
    }
    internal Complex()
    {
        Re = 0;
        Im = 0;
    }

    internal Complex SquareRoot()
    {
        if(Im == 0)
        {
            return Re < 0 ? new(0, Math.Sqrt(-Re)) : new(Math.Sqrt(Re));
        }
        else
        {
            double re = Math.Sqrt((Re + Math.Sqrt(Square())) / 2);
            double im = (Im / Math.Abs(Im)) * Math.Sqrt((-Re + Math.Sqrt(Square())) / 2);
            return new(re, im);
        }
    }
    internal Complex Sign() => new(Re/Absolute(), Im/Absolute());
    internal double Absolute() => Math.Sqrt(Square());
    internal double Square() => Re * Re + Im * Im;
    internal Complex Conjugate() => new(Re, -Im);

    public override string ToString() => Im == 0 ? string.Format("{0}", Math.Round(Re, 5)) : string.Format("{0} + i * {1}", Math.Round(Re, 5), Math.Round(Im, 5));
    public override int GetHashCode() => Re.GetHashCode() + Im.GetHashCode();
    public override bool Equals(object? obj) => (obj == null || !(obj is Complex)) ? false : this == (Complex)obj;

    public static Complex operator +(Complex A) => A;
    public static Complex operator -(Complex A) => new(-A.Re, -A.Im);

    public static Complex operator +(Complex A, Complex B) => new(A.Re + B.Re, A.Im + B.Im);
    public static Complex operator +(double A, Complex B) => new(A + B.Re, B.Im);
    public static Complex operator +(Complex A, double B) => new(A.Re + B, A.Im);

    public static Complex operator -(Complex A, Complex B) => new(A.Re - B.Re, A.Im - B.Im);
    public static Complex operator -(double A, Complex B) => new(A - B.Re, B.Im);
    public static Complex operator -(Complex A, double B) => new(A.Re - B, A.Im);

    public static Complex operator *(Complex A, Complex B) => new(A.Re * B.Re - A.Im * B.Im, A.Re * B.Im + B.Re * A.Im);
    public static Complex operator *(double A, Complex B) => new(A * B.Re, A * B.Im);
    public static Complex operator *(Complex A, double B) => new(A.Re * B, B * A.Im);

    public static Complex operator /(Complex A, Complex B)
    {
        if(B.Re == 0 && B.Im == 0) throw new DivideByZeroException();
        return new((A.Re * B.Re + A.Im * B.Im) / B.Square(), (B.Re * A.Im - A.Re * B.Im) / B.Square());
    }
    public static Complex operator /(double A, Complex B)
    {
        if(B.Re == 0 && B.Im == 0) throw new DivideByZeroException();
        return new(A * B.Re / B.Square(), A * B.Im / B.Square());
    }
    public static Complex operator /(Complex A, double B)
    {
        if(B == 0) throw new DivideByZeroException();
        return new(A.Re / B, A.Im / B);
    }

    public static bool operator ==(Complex A, Complex B) => A.Re == B.Re && A.Im == B.Im;
    public static bool operator !=(Complex A, Complex B) => !(A == B);

    public static bool operator >=(Complex A, Complex B) => A.Re >= B.Re;
    public static bool operator <=(Complex A, Complex B) => A.Re <= B.Re;

    public static bool operator >(Complex A, Complex B) => A.Re > B.Re;
    public static bool operator <(Complex A, Complex B) => A.Re < B.Re;
}