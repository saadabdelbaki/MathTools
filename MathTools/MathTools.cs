using System;
using System.Collections;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Numerics;
using System.Security.Cryptography;

namespace MathTools
{
    public static class MathTools
    {
        //tolerance for checking powers and roots being equal when using the double type for math
        private const double Tolerance = 0.000001;

        /// <summary>
        ///     Find all unique factorizations of a number N
        /// </summary>
        /// <param name="n"></param>
        /// <returns></returns>
        public static IEnumerable<IEnumerable<int>> FactorizationsOf(int n)
        {
            List<int> primeFactors = PrimeFactorization(n);
            var factorizations = GetAllPartitions(primeFactors.ToArray()).
                Select(partitions => partitions.
                    Select(partition => partition.Product()));
            var uniqueFactorizations = new HashSet<IEnumerable<int>>(factorizations, new MultiSetComparer<int>());
            uniqueFactorizations.RemoveWhere(x => x.Count() == 1 && x.First() == n);
            foreach (IEnumerable<int> factorization in uniqueFactorizations)
            {
                if (factorization.ToArray().Product() != n)
                {
                    throw new InvalidOperationException("factorization is not equal to number");
                }

            }
            return uniqueFactorizations;
        }


        public static BigInteger BigSqrt(BigInteger x)
        {
            int b = 15; // this is the next bit we try 
            BigInteger r = 0; // r will contain the result
            BigInteger r2 = 0; // here we maintain r squared
            while (b >= 0)
            {
                BigInteger sr2 = r2;
                BigInteger sr = r;
                // compute (r+(1<<b))**2, we have r**2 already.
                r2 += (uint)((r << (1 + b)) + (1 << (b + b)));
                r += (uint)(1 << b);
                if (r2 > x)
                {
                    r = sr;
                    r2 = sr2;
                }
                b--;
            }
            return r;
        }

        public static string[] ReadCsvFile(string filename)
        {
            var r = new StreamReader(filename);
            string line = r.ReadToEnd();
            r.Close();

            string[] values = line.Split(',');

            for (int i = 0; i < values.Length; i++)
            {
                values[i] = values[i].Trim('"');
            }

            return values;
        }

        public static int AlphaValueSum(string s)
        {
            return s.Sum(letter => Convert.ToInt32(letter) - 64);
        }

        public static BigInteger DigitSum(BigInteger number)
        {
            BigInteger sum = 0;
            BigInteger temp = number;
            while (temp > 0)
            {
                sum += temp % 10;
                temp = temp / 10;
            }
            return sum;
        }

        public static int DigitFactorialSum(int number)
        {
            int[] f = { 1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880 };
            int sum = 0;
            int temp = number;
            while (temp > 0)
            {
                sum += f[temp % 10];
                temp = temp / 10;
            }
            return sum;
        }

        public static BigInteger PowBigInteger(BigInteger b, BigInteger p)
        {
            BigInteger result = b;
            for (int i = 1; i < p; i++)
            {
                result *= b;
            }
            return result;
        }

        public static int PowInt(int b, int p)
        {
            int result = b;
            for (int i = 1; i < p; i++)
            {
                result *= b;
            }
            if (result < 0)
            {
                throw new OverflowException();
            }
            return result;
        }


        public static List<int> SpiralDiagonals(int size)
        {
            var diagonals = new List<int>();
            int skip = 1;
            int n = 1;
            int count = 0;
            while (n < size * size)
            {
                diagonals.Add(n);
                count++;
                n += skip + 1;
                if (count % 4 == 0)
                {
                    skip += 2;
                }
            }
            diagonals.Add(n);


            return diagonals;
        }

        public static BigInteger XToPowerN(int x, int n)
        {
            BigInteger result = x;
            for (int i = 1; i < n; i++)
            {
                result *= x;
            }
            return result;
        }

        /// <summary>
        ///     Calculates the first N digits of Sqrt(s)
        /// </summary>
        /// <param name="s">square</param>
        /// <param name="l">digits</param>
        /// <returns></returns>
        public static int[] FirstNDigitsOfSqrt(int s, int l)
        {
            var sqrtDigits = new List<int>();
            int digitPair;

            BigInteger remainder = 0;

            BigInteger p;
            int x;
            BigInteger y = 0;

            List<int> square = IntToDigitArray(s).ToList();

            if (square.Count % 2 == 1)
            {
                square.Insert(0, 0);
            }

            while (sqrtDigits.Count < l)
            {
                if (square.Any())
                {
                    digitPair = square[0] * 10 + square[1];
                    square.RemoveRange(0, 2);
                }
                else
                {
                    digitPair = 0;
                }

                BigInteger c = remainder * 100 + digitPair;

                if (!sqrtDigits.Any())
                {
                    p = 0;
                }
                else
                {
                    p = DigitArrayToBigInt(sqrtDigits.ToArray());
                }

                for (x = 0; x < 10; x++)
                {
                    BigInteger yNext = (x + 1) * (20 * p + (x + 1)); //next value of y
                    if (yNext > c)
                    {
                        y = x * (20 * p + x);
                        break;
                    }
                }
                sqrtDigits.Add(x);
                remainder = c - y;
            }
            return sqrtDigits.ToArray();
        }

        public static void PrintList<T>(IEnumerable<T> list)
        {
            foreach (T item in list)
            {
                Console.Write("{0}, ", item);
            }
            Console.WriteLine();
        }

        public static void PrintListOfLists<T>(IEnumerable<IEnumerable<T>> lists)
        {
            foreach (var list in lists)
            {
                PrintList(list);
            }
            Console.WriteLine();
        }

        public static int Mod(int x, int m)
        {
            return (x % m + m) % m;
        }

        #region Extension Methods

        public static int Product(this int[] numbers)
        {
            return numbers.Aggregate(1, (current, number) => current * number);
        }

        #endregion

        #region Primes

        //prime stuff
        public static bool IsPrime(int n)
        {
            if (n <= 1) return false;
            if (n == 2) return true;
            if (n % 2 == 0) return false;
            if (n < 9) return true;
            if (n % 3 == 0) return false;

            int counter = 5;
            while ((counter * counter) <= n)
            {
                if (n % counter == 0) return false;
                if (n % (counter + 2) == 0) return false;
                counter += 6;
            }
            return true;
        }

        public static bool BigIsPrime(BigInteger n)
        {
            if (n <= 1) return false;
            if (n == 2) return true;
            if (n % 2 == 0) return false;
            if (n < 9) return true;
            if (n % 3 == 0) return false;

            BigInteger counter = 5;
            while ((counter * counter) <= n)
            {
                if (n % counter == 0) return false;
                if (n % (counter + 2) == 0) return false;
                counter += 6;
            }
            return true;
        }

        public static int NthPrime(int n)
        {
            int p = 1;
            int pthPrime = 2; //first prime is 2
            int i = pthPrime + 1;

            while (p < n)
            {
                if (IsPrime(i))
                {
                    pthPrime = i;
                    p++;
                }

                i += 2;
            }
            return pthPrime;
        }

        public static List<int> PrimesBelowN(int n) //using Sieve of Eratosthenes
        {
            var primes = new List<int>();

            var sieve = new bool[n]; //make an empy array wih n elements
            for (int i = 0; i < sieve.Length; i++) //initialise all elements to true
            {
                sieve[i] = true;
            }

            for (int i = 2; i < Math.Sqrt(n); i++) //starting at 2
            {
                if (!sieve[i]) continue;
                for (int j = i * i; j < n; j += i)
                {
                    sieve[j] = false;
                }
            }
            // sieve is now only true for sieve[i] where i is prime

            //put primes in list
            for (int i = 2; i < sieve.Length; i++)
            {
                if (sieve[i])
                {
                    primes.Add(i);
                }
            }

            return primes;
        }

        public static int[] ESieve(int upperLimit)
        {
            if (upperLimit == 2)
            {
                return new[] { 2 };
            }
            int sieveBound = (upperLimit - 1) / 2;
            int upperSqrt = ((int)Math.Sqrt(upperLimit) - 1) / 2;

            var primeBits = new BitArray(sieveBound + 1, true);

            for (int i = 1; i <= upperSqrt; i++)
            {
                if (!primeBits.Get(i)) continue;
                for (int j = i * 2 * (i + 1); j <= sieveBound; j += 2 * i + 1)
                {
                    primeBits.Set(j, false);
                }
            }

            var numbers = new List<int>((int)(upperLimit / (Math.Log(upperLimit) - 1.08366))) { 2 };

            for (int i = 1; i <= sieveBound; i++)
            {
                if (primeBits.Get(i))
                {
                    numbers.Add(2 * i + 1);
                }
            }

            return numbers.ToArray();
        } //optimized Sieve of Eratosthenes

        public static List<int> PrimeFactorization(int n)
        {
            var primeFactors = new List<int>();

            //add all the 2s that divide n to the list
            while (n % 2 == 0)
            {
                primeFactors.Add(2);
                n = n / 2;
            }
            // n must be odd at this point.  So we can skip one element (Note i = i +2)
            for (int i = 3; i <= Math.Sqrt(n); i = i + 2)
            {
                // While i divides n, add i to the list and divide n
                while (n % i == 0)
                {
                    primeFactors.Add(i);
                    n = n / i;
                }
            }

            // This condition is to handle the case whien n is a prime number
            // greater than 2
            if (n > 2)
            {
                primeFactors.Add(n);
            }

            return primeFactors;
        }

        #endregion

        #region Divisors

        //divisor stuff
        public static int NumberOfDivisors(long n)
        {
            int numDivisors = 2; //all numbers n have at least factors 1 and n
            var sqrt = (int)Math.Sqrt(n);

            for (int i = 1; i < sqrt; i++)
            {
                if (n % i == 0)
                {
                    numDivisors += 2; //add 2 because the divisors are i and n/i
                }
            }
            //Correction if the number is a perfect square
            if (sqrt * sqrt == n)
            {
                numDivisors--;
                //added 2 earlier, but the square root would be counted twice, so subtract 1 from the total
            }
            return numDivisors;
        }

        public static int SumOfDivisors(int n)
        {
            return ListDivisors(n).Sum();
        }

        public static List<int> ListDivisors(int n)
        {
            var divisors = new List<int> { 1 };

            var sqrt = (int)Math.Sqrt(n);

            for (int i = 2; i <= sqrt; i++)
            {
                if (n % i == 0)
                {
                    //divisors are i and n/i
                    divisors.Add(i);
                    divisors.Add(n / i);
                }
            }
            return divisors.Distinct().ToList();
        }

        public static int GreatestCommonFactor(int a, int b)
        {
            int y;
            int x;

            if (a > b)
            {
                x = a;
                y = b;
            }
            else
            {
                x = b;
                y = a;
            }

            while (x % y != 0)
            {
                int temp = x;
                x = y;
                y = temp % x;
            }
            return y;
        }

        public static long GreatestCommonFactor(long a, long b)
        {
            long y;
            long x;

            if (a > b)
            {
                x = a;
                y = b;
            }
            else
            {
                x = b;
                y = a;
            }

            while (x % y != 0)
            {
                long temp = x;
                x = y;
                y = temp % x;
            }
            return y;
        }

        public static bool Coprime(int a, int b)
        {
            return GreatestCommonFactor(a, b) == 1;
        }

        public static bool Coprime(long a, long b)
        {
            return GreatestCommonFactor(a, b) == 1;
        }

        #endregion

        #region Combinatorics and Permutations

        //combinatorics
        public static long BinomialCoefficient(int n, int k)
        {
            //formulaa is ( Factorial(n) ) / ( Factorial(n - k) * Factorial(k) );
            //but the multiplicitive formula is more efficient for individual coefficients
            if (k < 0 || k > n)
            {
                return 0;
            }
            if (k == 0 || k == n)
            {
                return 1;
            }
            k = Math.Min(k, n - k); //take advantage of symmetry
            long c = 1;
            for (long i = 0; i < k; i++)
            {
                c = c * (n - i) / (i + 1);
            }
            return c;
        }

        public static BigInteger BigBinomialCoefficient(int n, int k)
        {
            //formulaa is ( Factorial(n) ) / ( Factorial(n - k) * Factorial(k) );
            //but the multiplicitive formula is more efficient for individual coefficients
            if (k < 0 || k > n)
            {
                return 0;
            }
            if (k == 0 || k == n)
            {
                return 1;
            }
            k = Math.Min(k, n - k); //take advantage of symmetry
            BigInteger c = 1;
            for (long i = 0; i < k; i++)
            {
                c = c * (n - i) / (i + 1);
            }
            return c;
        }

        public static BigInteger Factorial(long n)
        {
            BigInteger product = 1;

            for (long i = 1; i <= n; i++)
            {
                product = product * i;
            }
            return product;
        }

        public static List<int> NthPermuatation(int n, List<int> elements)
        {
            elements.Sort();
            var maxPermutations = (int)Factorial(elements.Count);
            //returns BigInt, cast to normal - make sure not to use this with big numbers
            if (n > maxPermutations)
            {
                throw new InvalidOperationException("n must be less than possible number of permutations");
            }
            if (n <= 0)
            {
                throw new InvalidOperationException("n must be greater than or equal to 1");
            }

            //start actual algorithm
            var permutedList = new List<int>();
            var toBePermuted = new List<int>(elements);
            int remainder = n - 1;

            for (int i = 1; i < elements.Count; i++)
            {
                var factorial = (int)Factorial(elements.Count - i);
                int factorialMultiplier = remainder / factorial;
                remainder = remainder % factorial;

                int selectedElement = toBePermuted[factorialMultiplier];
                toBePermuted.Remove(selectedElement);
                permutedList.Add(selectedElement);

                if (remainder == 0)
                {
                    break;
                }
            }

            permutedList.AddRange(toBePermuted);

            return permutedList;
        }

        public static bool IsPermutationOf(int a, int b)
        {
            var arr = new int[10];

            int temp = a;
            while (temp > 0)
            {
                arr[temp % 10]++;
                temp /= 10;
            }

            temp = b;
            while (temp > 0)
            {
                arr[temp % 10]--;
                temp /= 10;
            }

            for (int i = 0; i < 10; i++)
            {
                if (arr[i] != 0)
                {
                    return false;
                }
            }
            return true;
        }

        public static bool IsPermutationOf(long a, long b)
        {
            var arr = new int[10];

            long temp = a;
            while (temp > 0)
            {
                arr[temp % 10]++;
                temp /= 10;
            }

            temp = b;
            while (temp > 0)
            {
                arr[temp % 10]--;
                temp /= 10;
            }

            for (int i = 0; i < 10; i++)
            {
                if (arr[i] != 0)
                {
                    return false;
                }
            }
            return true;
        }

        public static List<List<int>> PermutationsOf(List<int> elements)
        {
            var maxPermutations = (int)Factorial(elements.Count);
            var permutations = new List<List<int>>();

            for (int i = 1; i <= maxPermutations; i++)
            {
                permutations.Add(NthPermuatation(i, elements));
            }
            return permutations;
        }

        public static List<List<int>> UniquePermutationsOf(List<int> elements)
        {
            var permutations = new List<List<int>>(PermutationsOf(elements));
            List<List<int>> unique = UniqueListsIn(permutations);
            return unique;
        }

        public static List<List<int>> UniqueListsIn(List<List<int>> lists)
        {
            var unique = new List<List<int>>(lists);
            for (int i = 0; i < unique.Count; i++)
            {
                List<int> perm = unique[i];
                for (int j = i + 1; j < unique.Count; j++)
                {
                    if (unique[j].SequenceEqual(perm))
                    {
                        unique.RemoveAt(j);
                    }
                }
            }
            return unique;
        }

        #endregion

        #region Arrays

        //array stuff
        public static int[] IntToDigitArray(long n)
        {
            int numDigits = n.ToString(CultureInfo.InvariantCulture).Length;

            var digitArray = new int[numDigits];
            long temp = n;
            for (int i = 0; i < numDigits; i++)
            {
                digitArray[numDigits - 1 - i] = (int)(temp % 10);
                temp = temp / 10;
            }
            return digitArray;
        }

        public static long DigitArrayToInt(int[] array)
        {
            long number = 0;
            long place = 1;
            int[] temp = array.Reverse().ToArray();
            foreach (int digit in temp)
            {
                number += digit * place;
                place *= 10;
            }
            return number;
        }

        public static BigInteger DigitArrayToBigInt(int[] array)
        {
            BigInteger number = 0;
            BigInteger place = 1;
            int[] temp = array.Reverse().ToArray();
            foreach (int digit in temp)
            {
                number += digit * place;
                place *= 10;
            }
            return number;
        }

        public static void AddDigitsToList(int n, List<int> list)
        {
            while (n > 0)
            {
                int digit = n % 10;
                n = n / 10;
                list.Add(digit);
            }
        }

        public static int[][] CopyJaggedArray(int[][] source)
        {
            return source.Select(s => s.ToArray()).ToArray();
        }

        #endregion

        #region Number Properties

        //number properties
        public static bool IsPandigital(int n)
        {
            int digits = 0;
            int count = 0;

            for (; n > 0; n /= 10, ++count)
            {
                if ((digits) == (digits |= 1 << (n - ((n / 10) * 10) - 1)))
                    return false;
            }

            return digits == (1 << count) - 1;
        }

        public static bool IsAmicable(int n)
        {
            //Let d(n) be defined as the sum of proper divisors of n (numbers less than n which divide evenly into n).
            //If d(a) = b and d(b) = a, where a ≠ b, then a and b are an amicable pair and each of a and b are called amicable numbers.     
            int temp = SumOfDivisors(n);
            if (SumOfDivisors(temp) == n && temp != n) //If d(a) = b and d(b) = a, where a ≠ b
            {
                return true;
            }
            return false;
        }

        public static bool IsCube(int n)
        {
            double cubeTest = Math.Pow(n, 1 / 3.0);
            return Math.Abs(cubeTest - ((int)cubeTest)) < Tolerance;
        }

        #endregion

        #region Euler's Totient / Phi

        //Euler's Totient / Phi
        public static int Totient(int n)
        {
            double totient = n;
            IEnumerable<int> primes = PrimeFactorization(n).Distinct();
            foreach (int prime in primes)
            {
                totient *= 1 - (1.0 / prime);
            }
            return (int)totient;
        }

        public static double PhiRatio(int n)
        {
            double phi = Phi(n);
            return n / phi;
        }

        public static long Phi(int n)
        {
            return PhiSieve(n)[n];
        }

        public static int[] PhiSieve(int limit)
        {
            int[] phi = Enumerable.Range(0, limit + 1).ToArray();
            for (int i = 2; i <= limit; i++)
            {
                if (phi[i] == i)
                {
                    for (int j = i; j <= limit; j += i)
                    {
                        phi[j] = phi[j] / i * (i - 1);
                    }
                }
            }
            return phi;
        }

        #endregion

        #region Geometric Numbers

        //geometric number stuff
        public static int TriangleNumber(int n)
        {
            //Triangle 	  	Tn=n(n+1)/2 	  	1, 3, 6, 10, 15, ...
            if (n < 1)
            {
                throw new InvalidOperationException("n must be 1 or greater");
            }
            return n * (n + 1) / 2;
        }

        public static long TriangleNumber(long n)
        {
            //Triangle 	  	Tn=n(n+1)/2 	  	1, 3, 6, 10, 15, ...
            if (n < 1)
            {
                throw new InvalidOperationException("n must be 1 or greater");
            }
            return n * (n + 1) / 2;
        }

        public static bool IsTriangular(int n)
        {
            double triTest = (Math.Sqrt(1 + 8 * n) - 1.0) / 2.0;
            return Math.Abs(triTest - ((int)triTest)) < Tolerance;
        }

        public static int SquareNumber(int n)
        {
            //Square 	  	Sn=n^2 	  	1, 4, 9, 16, 25, ...
            if (n < 1)
            {
                throw new InvalidOperationException("n must be 1 or greater");
            }
            return n * n;
        }

        public static long SquareNumber(long n)
        {
            //Square 	  	Sn=n^2 	  	1, 4, 9, 16, 25, ...
            if (n < 1)
            {
                throw new InvalidOperationException("n must be 1 or greater");
            }
            return n * n;
        }

        public static bool IsSquare(int n)
        {
            double squareTest = Math.Sqrt(n);
            return Math.Abs(squareTest - ((int)squareTest)) < Tolerance;
        }

        public static int PentagonalNumber(int n)
        {
            //Pentagonal 	  	Pn=n(3n−1)/2 	  	1, 5, 12, 22, 35, ...
            if (n < 1)
            {
                throw new InvalidOperationException("n must be 1 or greater");
            }
            return n * (3 * n - 1) / 2;
        }

        public static bool IsPentagonal(int n)
        {
            double penTest = (Math.Sqrt(1 + 24 * n) + 1.0) / 6.0;
            return Math.Abs(penTest - ((int)penTest)) < Tolerance;
        }

        public static bool IsPentagonal(long n)
        {
            double penTest = (Math.Sqrt(1 + 24 * n) + 1.0) / 6.0;
            return Math.Abs(penTest - ((long)penTest)) < Tolerance;
        }

        public static int HexagonalNumber(int n)
        {
            //Hexagonal 	  	Hn=n(2n−1) 	  	1, 6, 15, 28, 45, ...
            if (n < 1)
            {
                throw new InvalidOperationException("n must be 1 or greater");
            }
            return n * (2 * n - 1);
        }

        public static bool IsHexagonal(int n)
        {
            double hexTest = (Math.Sqrt(1 + 8 * n) + 1.0) / 4.0;
            return Math.Abs(hexTest - ((int)hexTest)) < Tolerance;
        }

        public static bool IsHexagonal(long n)
        {
            double hexTest = (Math.Sqrt(1 + 8 * n) + 1.0) / 4.0;
            return Math.Abs(hexTest - ((long)hexTest)) < Tolerance;
        }

        public static int HeptagonalNumber(int n)
        {
            //Heptagonal 	  	P7,n=n(5n−3)/2 	  	1, 7, 18, 34, 55, ...
            if (n < 1)
            {
                throw new InvalidOperationException("n must be 1 or greater");
            }
            return n * (5 * n - 3) / 2;
        }

        public static bool IsHeptagonal(int n)
        {
            //Heptagonal 	  	P7,n=n(5n−3)/2 	  	1, 7, 18, 34, 55, ...
            double heptTest = (Math.Sqrt(9 + 40 * n) + 3.0) / 10.0;
            return Math.Abs(heptTest - ((int)heptTest)) < Tolerance;
        }

        public static bool IsHeptagonal(long n)
        {
            //Heptagonal 	  	P7,n=n(5n−3)/2 	  	1, 7, 18, 34, 55, ...
            double heptTest = (Math.Sqrt(9 + 40 * n) + 3.0) / 10.0;
            return Math.Abs(heptTest - ((long)heptTest)) < Tolerance;
        }

        public static int OctagonalNumber(int n)
        {
            //Octagonal 	  	P8,n=n(3n−2) 	  	1, 8, 21, 40, 65, ...
            if (n < 1)
            {
                throw new InvalidOperationException("n must be 1 or greater");
            }
            return n * (3 * n - 2);
        }

        public static bool IsOctagonal(int n)
        {
            //Octagonal 	  	P8,n=n(3n−2) 	  	1, 8, 21, 40, 65, ...
            double octTest = (Math.Sqrt(1 + 3 * n) + 1.0) / 3.0;
            return Math.Abs(octTest - ((int)octTest)) < Tolerance;
        }

        public static bool IsOctagonal(long n)
        {
            //Octagonal 	  	P8,n=n(3n−2) 	  	1, 8, 21, 40, 65, ...
            double octTest = (Math.Sqrt(1 + 3 * n) + 1.0) / 3.0;
            return Math.Abs(octTest - ((long)octTest)) < Tolerance;
        }

        #endregion

        #region Palindromes

        //palindromes
        public static bool IsPalindrome(int n)
        {
            const int b = 10;
            int reversed = 0;
            int k = n;

            while (k > 0)
            {
                reversed = b * reversed + k % b;
                k /= b;
            }
            return n == reversed;
        }

        public static bool IsPalindrome(int n, int b)
        {
            int reversed = 0;
            int k = n;

            while (k > 0)
            {
                reversed = b * reversed + k % b;
                k /= b;
            }
            return n == reversed;
        }

        public static int ReverseNumber(int n)
        {
            const int b = 10;
            int reversed = 0;
            int k = n;

            while (k > 0)
            {
                reversed = b * reversed + k % b;
                k /= b;
            }
            return reversed;
        }

        public static BigInteger ReverseNumber(BigInteger number)
        {
            char[] k = number.ToString().ToCharArray();
            Array.Reverse(k);
            return BigInteger.Parse(new string(k));
        }

        public static bool IsPalindrome(BigInteger number)
        {
            return number == ReverseNumber(number);
        }

        #endregion

        #region Partitions

        public static IEnumerable<T[][]> GetAllPartitions<T>(T[] elements)
        {
            return GetAllPartitions(new T[][] { }, elements);
        }

        private static IEnumerable<T[][]> GetAllPartitions<T>(
            T[][] fixedParts, T[] suffixElements)
        {
            // A trivial partition consists of the fixed parts
            // followed by all suffix elements as one block
            yield return fixedParts.Concat(new[] { suffixElements }).ToArray();

            // Get all two-group-partitions of the suffix elements
            // and sub-divide them recursively
            IEnumerable<Tuple<T[], T[]>> suffixPartitions = GetTuplePartitions(suffixElements);
            foreach (var suffixPartition in suffixPartitions)
            {
                IEnumerable<T[][]> subPartitions = GetAllPartitions(
                    fixedParts.Concat(new[] { suffixPartition.Item1 }).ToArray(),
                    suffixPartition.Item2);
                foreach (var subPartition in subPartitions)
                {
                    yield return subPartition;
                }
            }
        }

        private static IEnumerable<Tuple<T[], T[]>> GetTuplePartitions<T>(
            T[] elements)
        {
            // No result if less than 2 elements
            if (elements.Length < 2) yield break;

            // Generate all 2-part partitions
            for (int pattern = 1; pattern < 1 << (elements.Length - 1); pattern++)
            {
                // Create the two result sets and
                // assign the first element to the first set
                List<T>[] resultSets =
                {
                    new List<T> {elements[0]}, new List<T>()
                };
                // Distribute the remaining elements
                for (int index = 1; index < elements.Length; index++)
                {
                    resultSets[(pattern >> (index - 1)) & 1].Add(elements[index]);
                }

                yield return Tuple.Create(
                    resultSets[0].ToArray(), resultSets[1].ToArray());
            }
        }

        #endregion

        /// <summary>
        ///     A comparer for lists of lists
        /// </summary>
        /// <typeparam name="T"></typeparam>
        public class MultiSetComparer<T> : IEqualityComparer<IEnumerable<T>>
        {
            public bool Equals(IEnumerable<T> first, IEnumerable<T> second)
            {
                if (first == null)
                    return second == null;

                if (second == null)
                    return false;

                if (ReferenceEquals(first, second))
                    return true;

                var firstCollection = first as ICollection<T>;
                var secondCollection = second as ICollection<T>;
                if (firstCollection != null && secondCollection != null)
                {
                    if (firstCollection.Count != secondCollection.Count)
                        return false;

                    if (firstCollection.Count == 0)
                        return true;
                }

                return !HaveMismatchedElement(first, second);
            }

            private static bool HaveMismatchedElement(IEnumerable<T> first, IEnumerable<T> second)
            {
                int firstCount;
                int secondCount;

                var firstElementCounts = GetElementCounts(first, out firstCount);
                var secondElementCounts = GetElementCounts(second, out secondCount);

                if (firstCount != secondCount)
                    return true;

                foreach (var kvp in firstElementCounts)
                {
                    firstCount = kvp.Value;
                    secondElementCounts.TryGetValue(kvp.Key, out secondCount);

                    if (firstCount != secondCount)
                        return true;
                }

                return false;
            }

            private static Dictionary<T, int> GetElementCounts(IEnumerable<T> enumerable, out int nullCount)
            {
                var dictionary = new Dictionary<T, int>();
                nullCount = 0;

                foreach (T element in enumerable)
                {
                    if (element == null)
                    {
                        nullCount++;
                    }
                    else
                    {
                        int num;
                        dictionary.TryGetValue(element, out num);
                        num++;
                        dictionary[element] = num;
                    }
                }

                return dictionary;
            }

            public int GetHashCode(IEnumerable<T> enumerable)
            {
                int hash = 17;

                foreach (T val in enumerable.OrderBy(x => x))
                    hash = hash * 23 + val.GetHashCode();

                return hash;
            }
        }

        //shuffle a list of objects
        public static void Shuffle<T>(this IList<T> list)
        {
            RNGCryptoServiceProvider provider = new RNGCryptoServiceProvider();
            int n = list.Count;
            while (n > 1)
            {
                byte[] box = new byte[1];
                do provider.GetBytes(box);
                while (!(box[0] < n * (Byte.MaxValue / n)));
                int k = (box[0] % n);
                n--;
                T value = list[k];
                list[k] = list[n];
                list[n] = value;
            }
        }

    }
}