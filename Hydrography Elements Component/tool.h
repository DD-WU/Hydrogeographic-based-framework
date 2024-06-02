#pragma once
#include <math.h>
#define G 9.8
#define C 0.5
#define PI 3.1415926
class Tool
{
public:
	//---------------------------------------------------------------------------
// GENERAL UTILITY BOBBINS
	static double norm(double* x, int n)
	{
		double res = 0.;
		for (int i = 0; i < n; i++) res += x[i] * x[i];
		return sqrt(res / n);
	}
	static double getmax(double a, double b)
	{
		return(a > b ? a : b);
	}
	static double getmin(double a, double b)
	{
		return(a < b ? a : b);
	}
	static double maxmag(double a, double b)
	{
		double res;
		if (fabs(a) > fabs(b)) res = a; else res = b;
		return res;
	}
	static double minmag(double a, double b)
	{
		double res;
		if (fabs(a) < fabs(b)) res = a; else res = b;
		return res;
	}
	static double InterpBC(double* varlist, double t) {
		int i = 0;
		double dt, a1, a2;
		double res = 0;

		// for values less than the start of the array - use 1st value
		if (t < varlist[1]) return(varlist[0]);

		while (varlist[(i + 1) * 2 + 1] > -0.9)
		{
			if (varlist[i * 2 + 1] <= t && varlist[i * 2 + 3] > t)
			{
				dt = varlist[i * 2 + 3] - varlist[i * 2 + 1];
				a1 = (t - varlist[i * 2 + 1]) / dt;
				a2 = 1.0 - a1;

				res = a1 * varlist[i * 2 + 2] + a2 * varlist[i * 2 + 0];
			}

			i++;
		}

		// for values greater than the end of the array - use last value
		if (t >= varlist[i * 2 + 1]) res = varlist[i * 2];

		return(res);
	};
	static double CalcA(double n, double s, double w, double Q)
	{
		double A;
		Q = fabs(Q);
		A = pow(n * pow(w, (2. / 3.)) * Q / fabs(s), (3. / 5.));

		return(A);
	}
	static int MaskTest(int m1, int m2)
	{
		if (m1 == -1 && m2 == -1) return(1);
		if (m1 == -1 && m2 > 0) return(1);
		if (m1 > 0 && m2 == -1) return(1);
		return(0);
	}
	static double waterDepthToInflow(double waterDepth, double A, double C1) {
		double r = sqrt(A / PI);
		if (r / waterDepth > C1) {
			return C * 2 * PI * r * sqrt(2 * G * waterDepth) * waterDepth; //m³/s
		}
		else
		{
			return C * A * sqrt(2 * G * waterDepth); //m³/s
		}
	}
private:

};


