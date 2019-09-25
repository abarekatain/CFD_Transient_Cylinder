using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CFDP3
{
    class Program
    {
        //Reference Type for each Node
        struct Node
        {
            public float r, Tet, T, aE, aP, aN, aW, aS, Sc, Sp, b,ap0,Tp0;

        }
        static float t = 0;
        static float dr, dtet , dt = 0.1f;
        static float r1 = 0.02f, r2 = 0.03f, r3 = 0.05f;
        static float GridScale = 0.5f;
        static int Ni = 90, Nj = 30;
        static float K = 1;
        static float h = 5,RoC = 1,q = 3000;
        static float T1 = 400, Tinf = 300, InitTemp = 350;
        static Node[,] Nodes;
        static float[,] OldTemps;
        static float Res = 100, Epsilon = 0.0003f;
        static float W = 1f;
        static float Q , QR,QU,QB,QS;
        static int Counter = 0,n=0;
        static float S = 500000;

        static void Main(string[] args)
        {
            Grid();
            AssignTemps();
            AssignSource();
            AssignCoeffs();

            Solve();
            CalculateQ();
            Console.ReadLine();

            //Export();
        }
        static void Grid()
        {
            Ni = (int)(Ni / GridScale);
            Nj = (int)(Nj / GridScale);

            dtet = (float)Math.PI/2/Ni;
            dr = (r3-r1) / Nj;
            Nodes = new Node[Ni + 2, Nj + 2];
            OldTemps = new float[Ni + 2, Nj + 2];
            //OldTwall = new float[Ni + 2];


            for (int i = 0; i < Ni + 2; i++)
            {

                for (int j = 0; j < Nj + 2; j++)
                {
                    if (i == 0) Nodes[i, j].Tet = 0;
                    else if (i == 1) Nodes[i, j].Tet = dtet / 2;
                    else if (i == Ni + 1) Nodes[i, j].Tet = Nodes[i - 1, j].Tet + dtet / 2;
                    else Nodes[i, j].Tet = Nodes[i - 1, j].Tet + dtet;

                    if (j == 0) Nodes[i, j].r = r1;
                    else if (j == 1) Nodes[i, j].r = r1 + dr / 2;
                    else if (j == Nj + 1) Nodes[i, j].r = Nodes[i, j - 1].r + dr / 2;
                    else Nodes[i, j].r = Nodes[i, j - 1].r + dr;
                }
            }

        }

        static void AssignTemps()
        {
            for (int i = 0; i < Ni + 2; i++)
            {
                for (int j = 0; j < Nj + 2; j++)
                {
                    if (j == 0) Nodes[i, 0].T = T1;
                    else
                        Nodes[i, j].T = InitTemp;
                }

            }
        }

        static void AssignSource()
        {
            for (int i = 0; i < Ni + 2; i++)
            {
                for (int j = 0; j < Nj + 2; j++)
                {
                    if (Nodes[i, j].r <= r2)
                    {
                        Nodes[i, j].Sc = S;
                        Nodes[i, j].Sp = 0;
                    }
                    else
                    {
                        Nodes[i, j].Sc = 0;
                        Nodes[i, j].Sp = 0;
                    }
                }
            }
        }

        static void AssignCoeffs()
        {
            for (int i = 1; i < Ni + 1; i++)
            {
                for (int j = 1; j < Nj + 1; j++)
                {
                    Nodes[i, j].aN = K * (Nodes[i,j].r + dr/2) * dtet / (Nodes[i, j + 1].r - Nodes[i, j].r);
                    Nodes[i, j].aS = K * (Nodes[i, j].r - dr / 2) * dtet / (Nodes[i, j].r - Nodes[i, j - 1].r);

                    if (i != 1)
                        Nodes[i, j].aW = K * dr / (Nodes[i, j].Tet - Nodes[i - 1, j].Tet) / Nodes[i,j].r;
                    else
                        Nodes[i, j].aW = 0;

                    if (i != Ni)
                        Nodes[i, j].aE = K * dr / (Nodes[i+1, j].Tet - Nodes[i, j].Tet) / Nodes[i, j].r;
                    else
                        Nodes[i, j].aE = 0;

                    Nodes[i, j].ap0 = RoC * (Nodes[i, j].r) * dr * dtet / dt ;
                    Nodes[i, j].aP = Nodes[i, j].aN + Nodes[i, j].aS + Nodes[i, j].aW + Nodes[i, j].aE + Nodes[i, j].ap0 - Nodes[i, j].Sp * (Nodes[i, j].r) * dr * dtet;
                    if (i != Ni)
                    {
                        Nodes[i, j].b = Nodes[i, j].Sc * (Nodes[i, j].r) * dr * dtet + Nodes[i, j].ap0 * InitTemp;
                    }
                    else
                    Nodes[i, j].b = Nodes[i, j].Sc * (Nodes[i, j].r) * dr * dtet + Nodes[i, j].ap0 * InitTemp + q * dr ;
                }
            }
        }


        static void Assignb()
        {
            for (int i = 1; i < Ni + 1; i++)
            {
                for (int j = 1; j < Nj + 1; j++)
                {
                    if (i != Ni)
                    {
                        Nodes[i, j].b = Nodes[i, j].Sc * (Nodes[i, j].r) * dr * dtet + Nodes[i, j].ap0 * Nodes[i, j].T;
                    }
                    else
                        Nodes[i, j].b = Nodes[i, j].Sc * (Nodes[i, j].r) * dr * dtet + Nodes[i, j].ap0 * Nodes[i,j].T + q * dr;
                }
            }
        }

        static void Solve()
        {
            while (t<=100)
            {
                t += dt;
                GaussSeidel();
                for (int i = 1; i < Ni+1; i++)
                {
                    Nodes[i, Nj + 1].T = (2 * K / dr * Nodes[i, Nj].T + h * Tinf) / (2 * K / dr + h);
                }
                Assignb();
                Export();
                Console.WriteLine("t= {0} Node[45,15] : T = {1} , b = {2}", t, Nodes[45, 15].T, Nodes[45, 15].b);
            }
        }

        static void GaussSeidel()
        {
            Res = 100;
            while (Res > Epsilon)
            {
                Counter++;
                for (int i = 1; i < Ni + 1; i++)
                {
                    for (int j = 1; j < Nj + 1; j++)
                    {

                        OldTemps[i, j] = Nodes[i, j].T;

                    }
                }

                for (int j = 1; j < Nj + 1; j++)
                {
                    for (int i = 1; i < Ni + 1; i++)
                    {
                        Nodes[i, j].T = (1 / Nodes[i, j].aP) * (Nodes[i, j].aE * Nodes[i + 1, j].T + Nodes[i, j].aW * Nodes[i - 1, j].T + Nodes[i, j].aN * Nodes[i, j + 1].T + Nodes[i, j].aS * Nodes[i, j - 1].T + Nodes[i, j].b);
                    }
                }
                CalculateResidual();
                for (int j = 1; j < Nj + 1; j++)
                {
                    for (int i = 1; i < Ni + 1; i++)
                    {
                        Nodes[i, j].T += (W - 1) * (Nodes[i, j].T - OldTemps[i, j]);
                    }
                }
                
            }
        }

        static void CalculateQ()
        {
            Q = 0;
            for (int i = 0; i < Ni + 2; i++)
            {
                for (int j = 0; j < Nj + 2; j++)
                {
                    if (j == 0)
                        QB += K * (Nodes[i, 1].T - Nodes[i, 0].T) / (dr / 2) * r1 * dtet;
                    else if(j==Nj+1)
                        QU += K * (Nodes[i, Nj].T - Nodes[i, Nj+1].T) / (dr / 2) * r3 * dtet;

                }

            }
            QR = q * (r3 - r1);
            QS = S * (r1+r2) * (r2-r1) * (float)Math.PI/4;
            Q = QB + QR + QS + QU;

            Console.WriteLine(QB + "\n" + QR + "\n" + QS + "\n" + QU + "\n" + Q);
        }

        static void CalculateResidual()
        {
            Res = 0;
            for (int i = 1; i < Ni + 1; i++)
            {
                for (int j = 1; j < Nj + 1; j++)
                {
                    Res += Math.Abs(Nodes[i, j].T - OldTemps[i, j]);
                }
            }
            Res /= Ni * Nj;
        }



        static void Export()
        {
            string[] Data = new string[(Ni + 2) * (Nj + 2) + 2];
            Data[0] = "VARIABLES = X, Y, T";
            Data[1] = string.Format("ZONE I = {0} , J = {1} , SOLUTIONTIME = {2}", Nj + 2, Ni + 2,t);
            for (int i = 0; i < Ni + 2; i++)
            {
                for (int j = 0; j < Nj + 2; j++)
                {
                    Data[i * (Nj + 2) + j + 2] = (Nodes[i, j].r * Math.Sin(Nodes[i, j].Tet)) + "\t" + (Nodes[i, j].r * Math.Cos(Nodes[i, j].Tet)) + "\t" + Nodes[i, j].T;
                }
            }
            File.WriteAllLines(string.Format(@"F:\Plots\{0}.plt",n++), Data);

        }


    }
}
