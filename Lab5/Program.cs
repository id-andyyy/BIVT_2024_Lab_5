using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Numerics;
using System.Reflection;
using System.Reflection.Metadata;
using System.Text;
using static System.Runtime.InteropServices.JavaScript.JSType;

public class Program
{
    public void PrintMatrix(int[,] matrix)
    {
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                Console.Write($"{matrix[i, j]}\t");
            }
            Console.WriteLine();
        }
        Console.WriteLine();
    }
    
    public static void Main()
    {
        Program program = new Program();
        double[,] one = null, two = null;
        program.Task_3_1(ref one, ref two);
    }
    
    #region Level 1

    public long Combinations(int n, int k)
    {
        return Factorial(n) / (Factorial(k) * Factorial(n - k));
    }

    public long Factorial(int n)
    {
        long answer = 1;

        for (int i = 2; i <= n; i++)
        {
            answer *= i;
        }

        return answer;
    }
    
    public long Task_1_1(int n, int k)
    {
        if (n < k || n < 0 || k < 0) return 0;
        
        long answer = Combinations(n, k);

        return answer;
    }

    public bool CheckTriangle(double a, double b, double c)
    {
        if (a + b > c && a + c > b && b + c > a) return true;
        return false;
    }

    public double GeronArea(double a, double b, double c)
    {
        double p = (a + b + c) / 2;

        return Math.Sqrt(p * (p - a) * (p - b) * (p - c));
    }

    public int Task_1_2(double[] first, double[] second)
    {
        double a1 = first[0], a2 = first[1], a3 = first[2], b1 = second[0], b2 = second[1], b3 = second[2];

        if (!CheckTriangle(a1, a2, a3) || !CheckTriangle(b1, b2, b3)) return -1;
        
        double S1 = GeronArea(a1, a2, a3), S2 = GeronArea(b1, b2, b3);

        if (S1 > S2) return 1;
        if (S1 < S2) return 2;
        return 0;
    }

    public double GetDistance(double v, double a, int t)
    {
        return v * t + a * t * t / 2;
    }
    
    public int Task_1_3a(double v1, double a1, double v2, double a2, int time)
    {
        double s1 = GetDistance(v1, a1, time), s2 = GetDistance(v2, a2, time);

        if (s1 > s2) return 1;
        if (s1 < s2) return 2;
        return 0;
    }

    public int Task_1_3b(double v1, double a1, double v2, double a2)
    {
        int time = 1;
        while (GetDistance(v1, a1, time) > GetDistance(v2, a2, time)) time++;

        return time;
    }
    #endregion

    #region Level 2

    public void FindMaxIndex(int[,] matrix, out int maxI, out int maxJ)
    {
        int max = matrix[0, 0];
        maxI = 0;
        maxJ = 0;

        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                if (matrix[i, j] > max)
                {
                    max = matrix[i, j];
                    maxI = i;
                    maxJ = j;
                }
            }
        }
    }
    
    public void Task_2_1(int[,] A, int[,] B)
    {
        int maxAI, maxAJ, maxBI, maxBJ;
        FindMaxIndex(A, out maxAI, out maxAJ);
        FindMaxIndex(B, out maxBI, out maxBJ);

        (A[maxAI, maxAJ], B[maxBI, maxBJ]) = (B[maxBI, maxBJ], A[maxAI, maxAJ]);
    }

    public void Task_2_2(double[] A, double[] B)
    {
        // code here

        // create and use FindMaxIndex(array);
        // only 1 array has to be changed!

        // end
    }

    public int FindDiagonalMaxIndex(int[,] matrix)
    {
        int max = matrix[0, 0], index = 0;

        for (int i = 1; i < matrix.GetLength(0); i++)
        {
            if (matrix[i, i] > max)
            {
                max = matrix[i, i];
                index = i;
            }
        }

        return index;
    }

    public void DeleteMatrixColumn(ref int[,] matrix, int index)
    {
        int[,] temp = new int[matrix.GetLength(0) - 1, matrix.GetLength(1)];

        for (int i = 0; i < temp.GetLength(0); i++)
        {
            for (int j = 0; j < temp.GetLength(1); j++)
            {
                if (i < index) temp[i, j] = matrix[i, j];
                else if (i > index - 1) temp[i, j] = matrix[i + 1, j];
            }
        }
        
        matrix = temp;
    }

    public void Task_2_3(ref int[,] B, ref int[,] C)
    {
        DeleteMatrixColumn(ref B, FindDiagonalMaxIndex(B));
        DeleteMatrixColumn(ref C, FindDiagonalMaxIndex(C));
    }

    public void Task_2_4(int[,] A, int[,] B)
    {
        // code here

        //  create and use method FindDiagonalMaxIndex(matrix); like in Task_2_3

        // end
    }

    public void FindMaxColumn(int[,] matrix, int columnIndex, out int rowIndex)
    {
        int max = matrix[0, 1];
        rowIndex = 0;

        for (int i = 1; i < matrix.GetLength(0); i++)
        {
            if (matrix[i, 1] > max)
            {
                max = matrix[i, 1];
                rowIndex = i;
            }
        }
    }

    public void Task_2_5(int[,] A, int[,] B)
    {
        int rowA, rowB;
        FindMaxColumn(A, 0, out rowA);
        FindMaxColumn(B, 0, out rowB);

        for (int j = 0; j < A.GetLength(1); j++)
        {
            (A[rowA, j], B[rowB, j]) = (B[rowB, j], A[rowA, j]);
        }
    }

    public void Task_2_6(ref int[] A, int[] B)
    {
        // code here

        // create and use FindMax(matrix, out row, out column); like in Task_2_1
        // create and use DeleteElement(array, index);

        // end
    }

    public int CountRowPositive(int[,] matrix, int rowIndex)
    {
        int cnt = 0;

        for (int j = 0; j < matrix.GetLength(1); j++)
        {
            if (matrix[rowIndex, j] > 0) cnt++;
        }

        return cnt;
    }

    public int CountColumnPositive(int[,] matrix, int colIndex)
    {
        int cnt = 0;

        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            if (matrix[i, colIndex] > 0) cnt++;
        }

        return cnt;
    }

    public void Task_2_7(ref int[,] B, int[,] C)
    {
        int maxB = int.MinValue, maxIndexB = -1;

        for (int i = 0; i < B.GetLength(0); i++)
        {
            int cnt = CountRowPositive(B, i);
            if (cnt > maxB)
            {
                maxB = cnt;
                maxIndexB = i;
            }
        }

        int maxC = int.MinValue, maxIndexC = -1;

        for (int j = 0; j < C.GetLength(1); j++)
        {
            int cnt = CountColumnPositive(C, j);
            if (cnt > maxC)
            {
                maxC = cnt;
                maxIndexC = j;
            }
        }

        int[,] temp = new int[B.GetLength(0) + 1, B.GetLength(1)];
        
        for (int i = 0; i < temp.GetLength(0); i++)
        {
            for (int j = 0; j < temp.GetLength(1); j++)
            {
                if (i <= maxIndexB) temp[i, j] = B[i, j];
                else if (i == maxIndexB + 1) temp[i, j] = C[j, maxIndexC];
                else temp[i, j] = B[i - 1, j];
            }
        }

        B = temp;
    }
    
    public void Task_2_8(int[] A, int[] B)
    {
        // code here

        // create and use SortArrayPart(array, startIndex);

        // end
    }

    public int[] SumPositiveElementsInColumns(int[,] matrix)
    {
        int[] array = new int[matrix.GetLength(1)];

        for (int j = 0; j < matrix.GetLength(1); j++)
        {
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                if (matrix[i, j] > 0) array[j] += matrix[i, j];
            }
        }

        return array;
    }

    public int[] Task_2_9(int[,] A, int[,] C)
    {
        int[] answer = new int[A.GetLength(1) + C.GetLength(1)];

        int[] arrayA = SumPositiveElementsInColumns(A), arrayC = SumPositiveElementsInColumns(C);

        for (int i = 0; i < arrayA.Length; i++)
        {
            answer[i] = arrayA[i];
        }

        for (int i = arrayA.Length; i < arrayA.Length + arrayC.Length; i++)
        {
            answer[i] = arrayC[i - arrayA.Length];
        }
        
        return answer;
    }

    public void Task_2_10(ref int[,] matrix)
    {
        // code here

        // create and use RemoveColumn(matrix, columnIndex);

        // end
    }

    public void Task_2_11(int[,] A, int[,] B)
    {
        int maxAI, maxAJ, maxBI, maxBJ;
        
        FindMaxIndex(A, out maxAI, out maxAJ);
        FindMaxIndex(B, out maxBI, out maxBJ);

        (A[maxAI, maxAJ], B[maxBI, maxBJ]) = (B[maxBI, maxBJ], A[maxAI, maxAJ]);
    }
    public void Task_2_12(int[,] A, int[,] B)
    {
        // code here

        // create and use FindMaxColumnIndex(matrix);

        // end
    }

    public void RemoveRow(ref int[,] matrix, int rowIndex)
    {
        int[,] temp = new int[matrix.GetLength(0) - 1, matrix.GetLength(1)];

        for (int i = 0; i < temp.GetLength(0); i++)
        {
            for (int j = 0; j < temp.GetLength(1); j++)
            {
                if (i < rowIndex) temp[i, j] = matrix[i, j];
                else if (i >= rowIndex) temp[i, j] = matrix[i + 1, j];
            }
        }

        matrix = temp;
    }

    public void Task_2_13(ref int[,] matrix)
    {
        int max = matrix[0, 0], maxRow = 0, min = matrix[0, 0], minRow = 0;

        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                if (matrix[i, j] > max)
                {
                    max = matrix[i, j];
                    maxRow = i;
                }

                if (matrix[i, j] < min)
                {
                    min = matrix[i, j];
                    minRow = i;
                }
            }
        }
        if (minRow < maxRow)
        {
            RemoveRow(ref matrix, maxRow);
            RemoveRow(ref matrix, minRow);
        }
        else if (minRow > maxRow)
        {
            RemoveRow(ref matrix, minRow);
            RemoveRow(ref matrix, maxRow);
        }
        else
        {
            RemoveRow(ref matrix, minRow);
        }
    }

    public void Task_2_14(int[,] matrix)
    {
        // code here

        // create and use SortRow(matrix, rowIndex);

        // end
    }

    public double GetAverageWithoutMinMax(int[,] matrix)
    {
        int max = matrix[0, 0], indexMaxI = 0, indexMaxJ = 0, min = matrix[0, 0], indexMinI = 0, indexMinJ = 0;
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                if (matrix[i, j] > max)
                {
                    max = matrix[i, j];
                    indexMaxI = i;
                    indexMaxJ = j;
                }

                if (matrix[i, j] < min)
                {
                    min = matrix[i, j];
                    indexMinI = i;
                    indexMinJ = j;
                }
            }
        }

        double sum = 0;
        int cnt = 0;
        
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                if ((i == indexMaxI || i == indexMinI) && (j == indexMaxJ || j == indexMinJ)) continue;
                sum += matrix[i, j];
                cnt++;
            }
        }

        return sum / cnt;
    }

    public int Task_2_15(int[,] A, int[,] B, int[,] C)
    {
        double[] array = { GetAverageWithoutMinMax(A), GetAverageWithoutMinMax(B), GetAverageWithoutMinMax(C) };

        if (array[0] < array[1] && array[1] < array[2]) return 1;
        if (array[0] > array[1] && array[1] > array[2]) return -1;
        return 0;
    }

    public void Task_2_16(int[] A, int[] B)
    {
        // code here

        // create and use SortNegative(array);

        // end
    }

    public void SortRowsByMaxElement(int[,] matrix)
    {
        int[] maxs = new int[matrix.GetLength(0)];

        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            maxs[i] = matrix[i, 0];
            for (int j = 1; j < matrix.GetLength(1); j++)
            {
                maxs[i] = (matrix[i, j] > maxs[i]) ? matrix[i, j] : maxs[i];
            }
        }

        for (int i = 1; i < maxs.Length; i++)
        {
            int k = i - 1;

            while (k >= 0 && maxs[k] < maxs[k + 1])
            {
                (maxs[k], maxs[k + 1]) = (maxs[k + 1], maxs[k]);
                
                for (int j = 0; j < matrix.GetLength(1); j++)
                {
                    (matrix[k, j], matrix[k + 1, j]) = (matrix[k + 1, j], matrix[k, j]);
                }

                k--;
            }
        }
    }

    public void Task_2_17(int[,] A, int[,] B)
    {
        SortRowsByMaxElement(A);
        SortRowsByMaxElement(B);
    }

    public void Task_2_18(int[,] A, int[,] B)
    {
        // code here

        // create and use SortDiagonal(matrix);

        // end
    }

    public void Task_2_19(ref int[,] matrix)
    {
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            bool flag = false;
            
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                if (matrix[i, j] == 0)
                {
                    flag = true;
                    break;
                }
            }
            
            if (!flag) continue;
            
            RemoveRow(ref matrix, i);
            i--;
        }
    }
    public void Task_2_20(ref int[,] A, ref int[,] B)
    {
        // code here

        // use RemoveColumn(matrix, columnIndex); from 2_10

        // end
    }

    public int[] CreateArrayFromMins(int[,] matrix)
    {
        int[] array = new int[matrix.GetLength(0)];
        
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            array[i] = matrix[i, i];
            for (int j = i + 1; j < matrix.GetLength(1); j++)
            {
                array[i] = (matrix[i, j] < array[i]) ? matrix[i, j] : array[i];
            }
        }

        return array;
    }

    public void Task_2_21(int[,] A, int[,] B, out int[] answerA, out int[] answerB)
    {
        answerA = CreateArrayFromMins(A);
        answerB = CreateArrayFromMins(B);
    }

    public void Task_2_22(int[,] matrix, out int[] rows, out int[] cols)
    {
        rows = null;
        cols = null;

        // code here

        // create and use CountNegativeInRow(matrix, rowIndex);
        // create and use FindMaxNegativePerColumn(matrix);

        // end
    }

    public void MatrixValuesChange(double[,] matrix)
    {
        double[] array = new double[matrix.GetLength(0) * matrix.GetLength(1)];

        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                array[i * matrix.GetLength(0) + j] = matrix[i, j];
            }
        }

        for (int i = 1; i < array.Length; i++)
        {
            int j = i - 1;

            while (j >= 0 && array[j] < array[j + 1])
            {
                (array[j], array[j + 1]) = (array[j + 1], array[j]);
                j--;
            }
        }
        
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                bool flag = false;
                
                for (int k = 0; k < 5; k++)
                {
                    if (matrix[i, j] == array[k])
                    {
                        flag = true;
                        matrix[i, j] = (matrix[i, j] > 0) ? matrix[i, j] * 2 : matrix[i, j] / 2;
                        break;
                    }
                }
                
                if (flag) continue;

                matrix[i, j] = (matrix[i, j] > 0) ? matrix[i, j] / 2 : matrix[i, j] * 2;
            }
        }
    }
    
    public void Task_2_23(double[,] A, double[,] B)
    {
        MatrixValuesChange(A);
        MatrixValuesChange(B);
    }

    public void Task_2_24(int[,] A, int[,] B)
    {
        // code here

        // use FindMaxIndex(matrix, out row, out column); like in 2_1
        // create and use SwapColumnDiagonal(matrix, columnIndex);

        // end
    }

    public int FindRowWithMaxNegativeCount(int[,] matrix)
    {
        int max = int.MinValue, index = -1;

        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            int cnt = CountNegativeInRow(matrix, i);

            if (cnt > max)
            {
                max = cnt;
                index = i;
            }
        }

        return index;
    }

    public int CountNegativeInRow(int[,] matrix, int rowIndex)
    {
        int cnt = 0;
        
        for (int j = 0; j < matrix.GetLength(1); j++)
        {
            if (matrix[rowIndex, j] < 0) cnt++;
        }

        return cnt;
    }

    public void Task_2_25(int[,] A, int[,] B, out int maxA, out int maxB)
    {
        maxA = FindRowWithMaxNegativeCount(A);
        maxB = FindRowWithMaxNegativeCount(B);
    }

    public void Task_2_26(int[,] A, int[,] B)
    {
        // code here

        // create and use FindRowWithMaxNegativeCount(matrix); like in 2_25
        // in FindRowWithMaxNegativeCount use CountNegativeInRow(matrix, rowIndex); from 2_22

        // end
    }

    public void FindRowMaxIndex(int[,] matrix, int rowIndex, out int columnIndex)
    {
        int max = matrix[rowIndex, 0];
        columnIndex = 0;

        for (int j = 0; j < matrix.GetLength(1); j++)
        {
            if (matrix[rowIndex, j] > max)
            {
                max = matrix[rowIndex, j];
                columnIndex = j;
            }
        }
    }

    public void ReplaceMaxElementOdd(int[,] matrix, int row, int column)
    {
        matrix[row, column] *= column + 1;
    }

    public void ReplaceMaxElementEven(int[,] matrix, int row, int column)
    {
        matrix[row, column] = 0;
    }

    public void ReplaceMatrixElements(int[,] matrix)
    {
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            int index;
            if ((i + 1) % 2 != 0)
            {
                FindRowMaxIndex(matrix, i, out index);
                ReplaceMaxElementOdd(matrix, i, index);
            }
            else
            {
                FindRowMaxIndex(matrix, i, out index);
                ReplaceMaxElementEven(matrix, i, index);
            }
        }
    }

    public void Task_2_27(int[,] A, int[,] B)
    {
        ReplaceMatrixElements(A);
        ReplaceMatrixElements(B);
    }

    public void Task_2_28a(int[] first, int[] second, ref int answerFirst, ref int answerSecond)
    {
        // code here

        // create and use FindSequence(array, A, B); // 1 - increasing, 0 - no sequence,  -1 - decreasing
        // A and B - start and end indexes of elements from array for search

        // end
    }

    public void Task_2_28b(int[] first, int[] second, ref int[,] answerFirst, ref int[,] answerSecond)
    {
        // code here

        // use FindSequence(array, A, B); from Task_2_28a or entirely Task_2_28a
        // A and B - start and end indexes of elements from array for search

        // end
    }

    public void Task_2_28c(int[] first, int[] second, ref int[] answerFirst, ref int[] answerSecond)
    {
        // code here

        // use FindSequence(array, A, B); from Task_2_28a or entirely Task_2_28a or Task_2_28b
        // A and B - start and end indexes of elements from array for search

        // end
    }
    #endregion

    #region Level 3

    public delegate double SumFunction(double x, int i, ref int iCalculated);
    public delegate double YFunction(double x);

    public double SumFunction1(double x, int i, ref int iFactorial)
    {
        if (i != 0) iFactorial *= i;
        
        return Math.Cos(i * x) / iFactorial;
    }

    public double SumFunction2(double x, int i, ref int iSign)
    {
        iSign *= -1;

        return iSign * Math.Cos(i * x) / (i * i);
    }

    public double YFunction1(double x)
    {
        return Math.Exp(Math.Cos(x)) * Math.Cos(Math.Sin(x));
    }

    public double YFunction2(double x)
    {
        return (x * x - Math.PI * Math.PI / 3) / 4;
    }

    public double CalculateSum(SumFunction sumFunction, double x, int i)
    {
        int iCalculated = 1;
        double sum = 0, summand = sumFunction(x, i++, ref iCalculated);
        
        while (Math.Abs(summand) > 1e-4)
        {
            sum += summand;
            summand = sumFunction(x, i++, ref iCalculated);
        }

        return sum;
    }
    
    public void GetSumAndY(double[,] sumAndY, SumFunction sumFunction, YFunction yFunction, double a, double h, int firstI)
    {
        for (int i = 0; i < sumAndY.GetLength(0); i++)
        {
            double x = a + h * i;

            sumAndY[i, 0] = CalculateSum(sumFunction, x, firstI);
            sumAndY[i, 1] = yFunction(x);
        }
    }
    
    public void Task_3_1(ref double[,] firstSumAndY, ref double[,] secondSumAndY)
    {
        double a1 = 0.1, b1 = 1, h1 = 0.1;
        firstSumAndY = new double[(int) ((b1 - a1) / h1) + 1, 2];
        GetSumAndY(firstSumAndY, SumFunction1, YFunction1, a1, h1, 0);
        
        double a2 = Math.PI / 5, b2 = Math.PI, h2 = Math.PI / 25;
        secondSumAndY = new double[(int)((b2 - a2) / h2) + 1, 2];
        GetSumAndY(secondSumAndY, SumFunction2, YFunction2, a2, h2, 1);
    }

    public void Task_3_2(int[,] matrix)
    {
        // SortRowStyle sortStyle = default(SortRowStyle); - uncomment me

        // code here

        // create and use public delegate SortRowStyle(matrix, rowIndex);
        // create and use methods SortAscending(matrix, rowIndex) and SortDescending(matrix, rowIndex)
        // change method in variable sortStyle in the loop here and use it for row sorting

        // end
    }
    
    public delegate void SwapDirector(double[] array);

    public void SwapRight(double[] array)  // start swapping from the last element
    {
        for (int i = array.Length - 1; i > 0; i -= 2)
        {
            (array[i], array[i - 1]) = (array[i - 1], array[i]);
        }
    }

    public void SwapLeft(double[] array)  // start swapping from the first element
    {
        for (int i = 0; i < array.Length - 1; i += 2)
        {
            (array[i], array[i + 1]) = (array[i + 1], array[i]);
        }
    }

    public double GetSum(double[] array, int start, int h)
    {
        double sum = 0;
        
        for (int i = start; i < array.Length; i += h)
        {
            sum += array[i];
        }

        return sum;
    }

    public double Task_3_3(double[] array)
    {
        SwapDirector swapper;

        double avg = 0;

        foreach (double el in array)
        {
            avg += el;
        }

        avg /= array.Length;

        if (array[0] > avg) swapper = SwapLeft;
        else swapper = SwapRight;

        swapper(array);

        return GetSum(array, 1, 2);
    }

    public int Task_3_4(int[,] matrix, bool isUpperTriangle)
    {
        int answer = 0;

        // code here

        // create and use public delegate GetTriangle(matrix);
        // create and use methods GetUpperTriange(array) and GetLowerTriange(array)
        // create and use GetSum(GetTriangle, matrix)

        // end

        return answer;
    }

    public int CountSignFlips(YFunction yFunction, double a, double b, double h)
    {
        double lastY = yFunction(a);
        int cnt = 0;
        for (double i = a + h; i <= b; i += h)
        {
            double curY = yFunction(i);
            if (curY > 0 && lastY < 0 || curY < 0 && lastY > 0) cnt++;
            lastY = curY;
        }

        return cnt + 1;
    }

    public double FunctionY1(double x)
    {
        return x * x - Math.Sin(x);
    }

    public double FunctionY2(double x)
    {
        return Math.Exp(x) - 1;
    }
    
    public void Task_3_5(out int func1, out int func2)
    {
        func1 = CountSignFlips(FunctionY1, 0, 2, 0.1);
        func2 = CountSignFlips(FunctionY2, -1, 1, 0.2);
    }

    public void Task_3_6(int[,] matrix)
    {
        // code here

        // create and use public delegate FindElementDelegate(matrix);
        // use method FindDiagonalMaxIndex(matrix) like in Task_2_3;
        // create and use method FindFirstRowMaxIndex(matrix);
        // create and use method SwapColumns(matrix, FindDiagonalMaxIndex, FindFirstRowMaxIndex);

        // end
    }

    public delegate int CountPositive(int[,] matrix, int index);

    public void InsertColumn(ref int[,] matrixB, int CountRow, int[,] matrixC, int CountColumn)
    {
        int[,] temp = new int[matrixB.GetLength(0) + 1, matrixC.GetLength(0)];

        for (int i = 0; i < temp.GetLength(0); i++)
        {
            for (int j = 0; j < temp.GetLength(1); j++)
            {
                if (i < CountRow) temp[i, j] = matrixB[i, j];
                else if (i == CountRow) temp[i, j] = matrixC[j, CountColumn];
                else temp[i, j] = matrixB[i - 1, j];
            }
        }

        matrixB = temp;
    }

    public void Task_3_7(ref int[,] B, int[,] C)
    {
        CountPositive countPositive;

        int max = int.MinValue, indexColumn = -1;
        countPositive = CountRowPositive;
        for (int i = 0; i < B.GetLength(0); i++)
        {
            int cnt = countPositive(B, i);
            if (cnt > max)
            {
                max = cnt;
                indexColumn = i;
            }
        }

        max = int.MinValue;
        int indexRow = -1;
        countPositive = CountColumnPositive;
        for (int j = 0; j < C.GetLength(1); j++)
        {
            int cnt = countPositive(C, j);
            if (cnt > max)
            {
                max = cnt;
                indexRow = j;
            }
        }
        
        InsertColumn(ref B, indexRow + 1, C, indexColumn);
    }

    public void Task_3_10(ref int[,] matrix)
    {
        // FindIndex searchArea = default(FindIndex); - uncomment me

        // code here

        // create and use public delegate FindIndex(matrix);
        // create and use method FindMaxBelowDiagonalIndex(matrix);
        // create and use method FindMinAboveDiagonalIndex(matrix);
        // use RemoveColumn(matrix, columnIndex) from Task_2_10
        // create and use method RemoveColumns(matrix, findMaxBelowDiagonalIndex, findMinAboveDiagonalIndex)

        // end
    }

    public void FindMinIndex(int[,] matrix, out int indexI, out int indexJ)
    {
        int min = matrix[0, 0];
        indexI = 0;
        indexJ = 0;

        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                if (matrix[i, j] < min)
                {
                    min = matrix[i, j];
                    indexI = i;
                    indexJ = j;
                }
            }
        }
    }

    public delegate void FindElementDelegate(int[,] matrix, out int indexI, out int indexJ);

    public void RemoveRows(ref int[,] matrix, FindElementDelegate findMax, FindElementDelegate findMin)
    {
        int maxI, maxJ, minI, minJ;
        findMax(matrix, out maxI, out maxJ);
        findMin(matrix, out minI, out minJ);
        
        if (minI < maxI)
        {
            RemoveRow(ref matrix, maxI);
            RemoveRow(ref matrix, minI);
        }
        else if (minI > maxI)
        {
            RemoveRow(ref matrix, minI);
            RemoveRow(ref matrix, maxI);
        }
        else
        {
            RemoveRow(ref matrix, minI);
        }
    }
    
    public void Task_3_13(ref int[,] matrix)
    {
        RemoveRows(ref matrix, FindMaxIndex, FindMinIndex);
    }

    public void Task_3_22(int[,] matrix, out int[] rows, out int[] cols)
    {

        rows = null;
        cols = null;

        // code here

        // create and use public delegate GetNegativeArray(matrix);
        // use GetNegativeCountPerRow(matrix) from Task_2_22
        // use GetMaxNegativePerColumn(matrix) from Task_2_22
        // create and use method FindNegatives(matrix, searcherRows, searcherCols, out rows, out cols);

        // end
    }

    public delegate void ReplaceMaxElement(int[,] matrix, int rowIndex, int max);

    public void EvenOddRowsTransform(int[,] matrix, ReplaceMaxElement replaceMaxElementOdd,
        ReplaceMaxElement replaceMaxElementEven)
    {
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            int index;
            if ((i + 1) % 2 != 0)
            {
                FindRowMaxIndex(matrix, i, out index);
                replaceMaxElementOdd(matrix, i, index);
            }
            else
            {
                FindRowMaxIndex(matrix, i, out index);
                replaceMaxElementEven(matrix, i, index);
            }
        }
    }

    public void Task_3_27(int[,] A, int[,] B)
    {
        EvenOddRowsTransform(A, ReplaceMaxElementOdd, ReplaceMaxElementEven);
        EvenOddRowsTransform(B, ReplaceMaxElementOdd, ReplaceMaxElementEven);
    }

    public void Task_3_28a(int[] first, int[] second, ref int answerFirst, ref int answerSecond)
    {
        // code here

        // create public delegate IsSequence(array, left, right);
        // create and use method FindIncreasingSequence(array, A, B); similar to FindSequence(array, A, B) in Task_2_28a
        // create and use method FindDecreasingSequence(array, A, B); similar to FindSequence(array, A, B) in Task_2_28a
        // create and use method DefineSequence(array, findIncreasingSequence, findDecreasingSequence);

        // end
    }

    public void Task_3_28c(int[] first, int[] second, ref int[] answerFirstIncrease, ref int[] answerFirstDecrease, ref int[] answerSecondIncrease, ref int[] answerSecondDecrease)
    {
        // code here

        // create public delegate IsSequence(array, left, right);
        // use method FindIncreasingSequence(array, A, B); from Task_3_28a
        // use method FindDecreasingSequence(array, A, B); from Task_3_28a
        // create and use method FindLongestSequence(array, sequence);

        // end
    }
    #endregion
    #region bonus part

    public delegate void MatrixConverter(double[,] matrix);

    public void ToUpperTriangular(double[,] matrix)
    {
        for (int j = 0; j < matrix.GetLength(1); j++)
        {
            for (int i = j + 1; i < matrix.GetLength(0); i++)
            {
                double k = matrix[i, j] / matrix[j, j];

                for (int jj = j; jj < matrix.GetLength(0); jj++)
                {
                    matrix[i, jj] -= k * matrix[j, jj];
                }
            }
        }
    }

    public void ToLowerTriangular(double[,] matrix)
    {
        for (int j = matrix.GetLength(1) - 1; j >= 0; j--)
        {
            for (int i = j - 1; i >= 0; i--)
            {
                double k = matrix[i, j] / matrix[j, j];

                for (int jj = j; jj >= 0; jj--)
                {
                    matrix[i, jj] -= k * matrix[j, jj];
                }
            }
        }
    }

    public void ToLeftDiagonal(double[,] matrix)
    {
        ToUpperTriangular(matrix);
        ToLowerTriangular(matrix);
    }

    public void ToRightDiagonal(double[,] matrix)
    {
        ToLowerTriangular(matrix);
        ToUpperTriangular(matrix);
    }
    
    public double[,] Task_4(double[,] matrix, int index)
    {
        MatrixConverter[] mc = new MatrixConverter[] { ToUpperTriangular, ToLowerTriangular, ToLeftDiagonal, ToRightDiagonal };
        
        mc[index](matrix);

        return matrix;
    }
    
    #endregion
}
