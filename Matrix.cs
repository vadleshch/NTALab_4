namespace AALab_5
{
    internal class Matrix
    {
        public double[,] M;
        public Matrix (int n, int m = 1) 
        { 
            M = new double[n, m];
        }
        public static Matrix add(Matrix A, Matrix B)
        {
            Matrix C = new Matrix(A.M.GetLength(0), A.M.GetLength(1));
            for (int i = 0; i < A.M.GetLength(0); i++)
            {
                for (int j = 0; j < A.M.GetLength(1); j++)
                {
                    C.M[i, j] = A.M[i, j] + B.M[i, j];
                }
            }
            return C;
        }
        public static Matrix sub(Matrix A, Matrix B)
        {
            Matrix C = new Matrix(A.M.GetLength(0), A.M.GetLength(1));
            for (int i = 0; i < A.M.GetLength(0); i++)
            {
                for (int j = 0; j < A.M.GetLength(1); j++)
                {
                    C.M[i, j] = A.M[i, j] - B.M[i, j];
                }
            }
            return C;
        }
        public static Matrix mul (Matrix A, Matrix B)
        {
            Matrix C = new Matrix(A.M.GetLength(0), B.M.GetLength(1));
            int n = A.M.GetLength(1);
            
            for (int i = 0; i < A.M.GetLength(0); i++)
            {
                for (int j = 0; j < B.M.GetLength(1); j++)
                {
                    for (int k = 0; k < n; k++)
                    {
                        C.M[i, j] += A.M[i, k] * B.M[k, j];
                    }
                }
            }
            return C;
        }
        public static (Matrix L, Matrix U, Matrix P) decomposition(Matrix M)
        {
            //A*P=L*U
            int n = M.M.GetLength(0);
            Matrix A = new Matrix(n, n);
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    A.M[i, j] = M.M[i, j];
                }
            }
            int[] p = new int[n];
            Matrix L = new Matrix(n, n);
            Matrix U = new Matrix(n, n);
            Matrix P = new Matrix(n, n);
            for (int i = 0; i < n; i++)
            {
                p[i] = i;
            }
            double v = 0;
            int k1 = 0;
            for (int k = 0; k < n; k++)
            {
                v = 0;
                k1 = k;
                for (int i = k; i < n; i++)
                {
                    if (v < Math.Abs(A.M[i, k]))
                    {
                        v = Math.Abs(A.M[i, k]);
                        k1 = i;
                    }
                }
                if (v == 0)
                { 
                    return (L, U, P);
                }
                if (k1 != k)
                {
                    for (int j = 0; j < n; j++)
                    {
                        A.M[k, j] = A.M[k, j] + A.M[k1, j];
                        A.M[k1, j] = A.M[k, j] - A.M[k1, j];
                        A.M[k, j] = A.M[k, j] - A.M[k1, j];
                    }
                    p[k] = p[k] + p[k1];
                    p[k1] = p[k] - p[k1];
                    p[k] = p[k] - p[k1];
                }
                for (int i = k + 1; i < n; i++)
                {
                    A.M[i, k] = A.M[i, k] / A.M[k, k];
                    for (int j = k + 1; j < n; j++)
                    {
                        A.M[i, j] = A.M[i, j] - A.M[i, k] * A.M[k, j];
                    }
                }
            }
            for (int i = 0; i < n; i++)
            {
                L.M[i, i] = 1;
                for (int j = 0; j < i; j++)
                {
                    L.M[i, j] = A.M[i, j];
                }
            }
            for (int i = 0; i < n; i++)
            {
                for (int j = i; j < n; j++)
                {
                    U.M[i, j] = A.M[i, j];
                }
            }
            for (int i = 0; i < n; i++)
            {
                P.M[i, p[i]] = 1;
            }
            return (L, U, P);
        }
    }
}
