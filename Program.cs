using AALab_5;
using System.Diagnostics;

uint swaps = 0;
Random r = new Random();
Matrix[] Matrixs = new Matrix[50];
for (int i = 0; i < 50; i++)
{
    Matrixs[i] = GenMatrix();
}

Matrix[] Copy(Matrix[] A)
{
    Matrix[] res = new Matrix[A.Length];
    for (int i = 0; i < A.Length; i++)
    {
        res[i] = new Matrix(A[i].M.GetLength(0), A[i].M.GetLength(1));
        for (int j = 0; j < A[i].M.GetLength(0); j++)
        {
            for (int k = 0; k < A[i].M.GetLength(1); k++)
            {
                res[i].M[j, k] = A[i].M[j, k];
            }
        }
    }
    return res;
}
void Test(double d, StreamWriter datawriter)
{
    Stopwatch sw = new Stopwatch();
    double sumT = 0;
    double sumL = 0;
    double sumK = 0;
    double l = 0;
    double k = 0;
    Matrix[] B = Copy(Matrixs);
    for (int i = 0; i < 50; i++)
    {
        for (int j = 0; j < B[i].M.GetLength(0); j++)
        {
            for (int j2 = 0; j2 < B[i].M.GetLength(1); j2++)
            {
                datawriter.Write(B[i].M[j, j2] + " ");
            }
            datawriter.WriteLine();
        }
        sw.Restart();
        B[i] = LLL(B[i], d);
        sw.Stop();
        l = Length(B[i]);
        k = CoefAdamara(B[i]);
        datawriter.WriteLine(sw.Elapsed.Milliseconds + " " + l + " " + k + " " + swaps + " " + d);
        sumT += sw.Elapsed.Milliseconds;
        sumL += l;
        sumK += k;
        Console.Write("|");
    }
    Console.WriteLine();
    Console.WriteLine(d + ": " + sumT / 50 + " " + sumL / 50 + " " + sumK / 50 + " " + swaps / 50);
}

Parallel.Invoke(
    () => Test(0.50, new StreamWriter("data1.txt")),
    () => Test(0.75, new StreamWriter("data2.txt")),
    () => Test(0.90, new StreamWriter("data3.txt")),
    () => Test(0.95, new StreamWriter("data4.txt")),
    () => Test(0.99, new StreamWriter("data5.txt"))
);


double Length(Matrix B)
{
    double res = 0;
    for (int i = 0; i < B.M.GetLength(0); i++)
    {
        res += B.M[i, 0] * B.M[i, 0];
    }
    return Math.Sqrt(res);
}

double CoefAdamara(Matrix B)
{
    double Logres = 0;
    for (int i = 0; i < B.M.GetLength(1); i++)
    {
        Logres = Logres + 0.5 * Math.Log(Scal(B, i, B, i));
    }
    return Math.Exp((det(B) - Logres) / B.M.GetLength(1));
}

double det(Matrix B)
{
    Matrix L;
    Matrix U;
    Matrix P;
    (L, U, P) = Matrix.decomposition(B);
    double det = 1;
    for (int i = 0; i < U.M.GetLength(0); i++)
    {
        det = det + Math.Log(Math.Abs(U.M[i, i]));
    }
    return det;
}

Matrix LLL(Matrix B, double d)
{
    Matrix res = B;
    int m = res.M.GetLength(1);
    int k = 1;
    Matrix Bstar;
    double[,] mu;
    double[] norms;
    (Bstar, mu, norms) = GrammaShmidta(res);
    while (k < m)
    {
        for (int j = k - 1; j >= 0; j--)
        {
            if (Math.Abs(mu[k, j]) > 0.5)
            {
                for (int i = 0; i < res.M.GetLength(0); i++)
                {
                    res.M[i, k] = res.M[i, k] - Math.Round(mu[k, j]) * res.M[i, j];
                }
                (Bstar, mu, norms) = GrammaShmidta(res);
            }
        }
        if (norms[k] < (d - mu[k, k - 1] * mu[k, k - 1]) * norms[k - 1])
        {
            Swap(res, k - 1, k);
            (Bstar, mu, norms) = GrammaShmidta(res);
            k = Math.Max(k - 1, 1);
        }   
        else
        {
            k++;
        }
    }
    return res;
}
(Matrix, double[,], double[]) GrammaShmidta(Matrix B)
{
    int n = B.M.GetLength(0);
    int m = B.M.GetLength(1);
    Matrix Bstar = new Matrix(n, m);
    double[,] mu = new double[m, m];
    double[] norms = new double[m];
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            Bstar.M[j, i] = B.M[j, i];
        }
        for (int j = 0; j < i; j++)
        {
            mu[i, j] = Scal(B, i, Bstar, j) / norms[j];
            for (int k = 0; k < n; k++)
            {
                Bstar.M[k, i] = Bstar.M[k, i] - mu[i, j] * Bstar.M[k, j];
            }
        }
        norms[i] = Scal(Bstar, i, Bstar, i);
    }
    return (Bstar, mu, norms);
}

double Scal(Matrix A, int x, Matrix B, int y)
{
    double res = 0;
    for (int i = 0; i < A.M.GetLength(0); i++)
    {
        res = res + A.M[i, x] * B.M[i, y];
    }
    return res;
}

void Swap(Matrix A, int x, int y)
{
    swaps++;
    double temp;
    for (int i = 0; i < A.M.GetLength(0); i++)
    {
        temp = A.M[i, x];
        A.M[i, x] = A.M[i, y];
        A.M[i, y] = temp;
    }
}

Matrix GenMatrix(int n = 30)
{ 
    Matrix res = new Matrix(n, n);
    for (int i = 0; i < n; i++)
    {
        res.M[i, i] = 1;
    }
    int x, y, c;
    for (int i = 0; i < 160; i++)
    {
        x = r.Next(0, n);
        y = r.Next(0, n);
        c = r.Next(-10, 10);
        if (x != y && c != 0)
        { 
            for (int j = 0; j < n; j++)
            {
                res.M[j, y] = res.M[j, y] + c * res.M[j, x];
            }
        }
    }
    for (int i = 0; i < 160; i++)
    {
        x = r.Next(0, n);
        y = r.Next(0, n);
        c = r.Next(-10, 10);
        if (x != y && c != 0)
        {
            for (int j = 0; j < n; j++)
            {
                res.M[j, y] = res.M[j, y] + c * res.M[j, x];
            }
        }
    }
    return res;
}