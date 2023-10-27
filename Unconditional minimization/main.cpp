#include <iostream>
#include <math.h>

#include <process.h>
#include <conio.h>

#include <limits>

using namespace std;

double epsilon{0.001};

double tau{0.1 * sqrt(epsilon)};


double f(double x, double y)
{
	return x*x+2*y*y;
}

double df_dx_1(double x, double y)
{
	return (f(x + tau, y) - f(x - tau, y)) / (2 * tau);
}

double df_dx_2(double x, double y)
{
	return (f(x, y + tau) - f(x, y - tau)) / (2 * tau);
}

double d2f_dx1_dx1(double _x, double _y)
{
	return (f(_x+tau,_y)-2*f(_x,_y)+f(_x-tau,_y))/(tau*tau);
}

double d2f_dx2_dx2(double _x, double _y)
{
	return (f(_x,_y-tau)-2*f(_x,_y)+f(_x,_y+tau))/(tau*tau);
}

double d2f_dx1_dx2(double _x, double _y)
{
	return (f(_x + tau, _y + tau) - f(_x + tau, _y - tau) - f(_x - tau, _y + tau) + f(_x - tau, _y - tau)) / (4 * tau * tau);
}


double Euclid_norm(double A_1, double A_2, double B_1, double B_2)
{
	double vector_1 = B_1 - A_1;
	double vector_2 = B_2 - A_2;

	return sqrt(pow(vector_1, 2) + pow(vector_2, 2));

}

double norm_of_vector(double x, double y)
{
	return sqrt(pow(x, 2) + pow(y, 2));
}



double phi(double t, double _current_point__1, double _current_point__2, double H__1, double H__2)
{
	return f(_current_point__1 + t * H__1, _current_point__2 + t * H__2);
}

double find_alpha_k(double current_point__1, double current_point__2, double H__1, double H__2)
{
	double a = 0; 
	double b = 3;

	double n = 1;

	double a_n_1 = a;
	double b_n_1 = b;


	while ((b - a - epsilon) / (pow(2, n + 1)) > (epsilon / 2))
	{
		
		double x_1_n_1 = ((a_n_1 + b_n_1) / 2) - epsilon / 2;
		double x_2_n_1 = ((a_n_1 + b_n_1) / 2) + epsilon / 2;

		

		if(phi(x_1_n_1, current_point__1, current_point__2, H__1, H__2) < phi(x_2_n_1, current_point__1, current_point__2, H__1, H__2))
		{
			b_n_1 = x_2_n_1;
			n++;
			continue;
		}

		if (phi(x_1_n_1, current_point__1, current_point__2, H__1, H__2) > phi(x_2_n_1, current_point__1, current_point__2, H__1, H__2))
		{
			a_n_1 = x_1_n_1;
			n++;
			continue;
		}

		n++;

	}

	return (b_n_1 + a_n_1) / 2;
}


double find_alpha_k__(double current_point__1, double current_point__2, double H__1, double H__2)
{

	double beta = 0.9;
	double lambda = 0.4;
	double mu = 2;

	double alpha = beta;

link2:
	if (f(current_point__1 + alpha * H__1, current_point__2 + alpha * H__2) < f(current_point__1, current_point__2))
	{
		if (alpha==beta)
		{

			alpha = mu * beta;

link:
			if (f(current_point__1+alpha*H__1, current_point__2+alpha*H__2)<f(current_point__1+beta*H__1, current_point__2+beta*H__2))
			{
				alpha = alpha * mu;
				goto link;
			}
			else
			{
				return alpha;
			}

		}
		else
		{
			return alpha;
		}

	}
	else
	{
		alpha = lambda * beta;
		goto link2;
	}
	
	abort();
}

int main()
{
	
	std::cout << "Input the start point:\n";

	double initial_point__1;
	double initial_point__2;

	cin >> initial_point__1;
	cin >> initial_point__2;

	double current_point__1 = initial_point__1;
	double current_point__2 = initial_point__2;
	
	int k = 1;

	for (;;k++)
	{
		double H__1 = df_dx_1(current_point__1,current_point__2);
		double H__2 = df_dx_2(current_point__1, current_point__2);

		H__1 = -1.0 * H__1;
		H__2 = -1.0 * H__2;

		double alpha = find_alpha_k(current_point__1,current_point__2,H__1,H__2);

		double next_point__1 = current_point__1 + alpha * H__1;
		double next_point__2 = current_point__2 + alpha * H__2;

		if(Euclid_norm(next_point__1,next_point__2,current_point__1,current_point__2)<sqrt(epsilon) && norm_of_vector(H__1,H__2)<sqrt(epsilon) && abs(f(current_point__1,current_point__2)-f(next_point__1,next_point__2))<sqrt(epsilon))
		{
			current_point__1 = next_point__1;
			current_point__2 = next_point__2;
			break;
		}

		current_point__1 = next_point__1;
		current_point__2 = next_point__2;
	}

	// START OF 2 METHOD

	cout << '\n' << '\n' << '\n' << '\n';

	double current_point__1__ = current_point__1;
	double current_point__2__ = current_point__2;

	int l = 1;

	double g_11;
	double g_12;
	double g_21;
	double g_22;

	for (;; l++)
	{
		g_11 = d2f_dx1_dx1(current_point__1__, current_point__2__);
		g_12 = d2f_dx1_dx2(current_point__1__, current_point__2__);
		g_21 = d2f_dx1_dx2(current_point__1__, current_point__2__);
		g_22 = d2f_dx2_dx2(current_point__1__, current_point__2__);

		
		double H__1__ = (g_12 * df_dx_2(current_point__1__,current_point__2__) - g_22 * df_dx_1(current_point__1__, current_point__2__)) / (g_11 * g_22 - g_12 * g_21);
		double H__2__ = (g_11 * df_dx_2(current_point__1__, current_point__2__) - g_21 * df_dx_1(current_point__1__, current_point__2__)) / (g_21 * g_12 - g_11 * g_22);
		
		double alpha__ = find_alpha_k__(current_point__1__, current_point__2__, H__1__, H__2__);

		double next_point__1__ = current_point__1__ + alpha__ * H__1__;
		double next_point__2__ = current_point__2__ + alpha__ * H__2__;

		if(Euclid_norm(next_point__1__, next_point__2__, current_point__1__, current_point__2__) < epsilon && norm_of_vector(H__1__, H__2__) < epsilon && abs(f(current_point__1__, current_point__2__) - f(next_point__1__, next_point__2__)) < epsilon)
		{
			current_point__1__ = next_point__1__;
			current_point__2__ = next_point__2__;
			break;
		}

		current_point__1__ = next_point__1__;
		current_point__2__ = next_point__2__;
		
	}

	std::cout << "----------------------------------\n";
	std::cout << "1 method\n----------------------------------\nPoint: ";
	std::cout << "(" << current_point__1 << "," << current_point__2 << ")" << "\nNumber of iterations: ";
	std::cout << k;
	std::cout << "\n\n";

	std::cout << "----------------------------------"<<"\n";
	std::cout << "2 method\n----------------------------------\nPoint: ";
	std::cout << "(" << current_point__1__ << "," << current_point__2__ << ")" << "\nNumber of iterations: ";
	std::cout << l;
	std::cout << "\n\n";
	
	return 1;
}

//--------------------------------------------
// x*x+2*y*y
// (1,1)
// a = 0
// b = 3
// beta = 0.4
// lambda = 0.3
// mu = 2
// -------------------------------------------
// exp((x-1)*(x-1)+0.1*(y-2)*(y-2))
// (1,1), (0,0)
// a = 0
// b = 3
// beta = 0.9
// lambda = 0.4
// mu = 2
// -------------------------------------------
// sin((x-1)*(x-1)+0.1*(y-2)*(y-2))
// (3,1), (0,0)
// a = 0
// b = 3
// beta = 0.9
// lambda = 0.4
// mu = 2
// -------------------------------------------
// (x-1)*(x-1)+0.1*(y-2)*(y-2)
// (1,1), (4,4), (-1,4), (0,0)
// a = 0
// b = 3
// beta = 0.9
// lambda = 0.4
// mu = 2
// -------------------------------------------
//
//
//
//
//
//
//
//
//
// 
//





