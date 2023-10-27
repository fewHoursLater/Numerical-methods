import math
import numpy

a = 0.0;
b = 1.0;

epsilon = 0.000001;

def K(x: float, t: float):
    return x*t*t;

def f(x):
    return 1;

def exact_solution(x: float):
    return 1+4*x/9;


def U_n(x: float, n: int):
    
    n = n*10;

    h = (b-a)/(n-1);

    C_0 = h/2;
    C_i = h;
    C_last = h/2;

    value = f(x);

    var = Gauss_solve(n-1);

    for i in range(n):

        t_i=a+i*h;

        if i==0:
            value += C_0*K(x, t_i)*var[0];
            continue;

        if i==(n-1):
            value += C_last*K(x, t_i)*var[n-1];
            continue;

        value += C_i*K(x, t_i)*var[i];

    return value;

def ksi(x: float, i: int):
    
    return pow(U_n(x, i+1)-U_n(x, i),2);

def diff(x: float, k: int):
    return pow(U_n(x,k)-exact_solution(x),2);

def Gauss_solve(number_of_segmetns: int): 
    h = (b-a)/number_of_segmetns;         
    right_parts = [];

    n=0;

    while a + n*h<=b:
        x = a+n*h;
        right_parts.append(f(x));
        n+=1;

    Row = number_of_segmetns +1;
    Column = number_of_segmetns+1;

    matrix = [];

    for j in range(Row):   
        
        temp_ = [];
        
        t_first_argument=a+j*h;

        for i in range(Column):
            
            t_second_argument = a+i*h;

            if i==0:

                if i==j:
                    temp_.append(1-(h/2)*K(t_first_argument,t_second_argument));
                    continue;

                else:
                    temp_.append((-h/2)*K(t_first_argument,t_second_argument));
                    continue;
            

            if i==number_of_segmetns:

                if i==j:
                    temp_.append(1-(h/2)*K(t_first_argument,t_second_argument));
                    continue;

                else:
                    temp_.append((-h/2)*K(t_first_argument,t_second_argument));
                    continue;


            if i==j:
                temp_.append(1-h*K(t_first_argument,t_second_argument));
                continue;

            else:
                temp_.append(-h*K(t_first_argument,t_second_argument));
                continue;
                
        matrix.append(temp_)
    
    G = numpy.array(matrix);
    v = numpy.array(right_parts);

    return numpy.linalg.solve(G, v);

def integrate(k: int, left_bound: float, right_bound: float):     
    C_1 = 5/9;
    C_2 = 8/9;
    C_3 = 5/9;
    
    t_1 = -1*math.sqrt(15)/5;
    t_3 = math.sqrt(15)/5;

    x_1 = ((right_bound+left_bound)/2)+((right_bound-left_bound)/2)*t_1;
    x_2 = ((right_bound+left_bound)/2);
    x_3 = ((right_bound+left_bound)/2)+((right_bound-left_bound)/2)*t_3;

    return ((right_bound-left_bound)/2)*(C_1*ksi(x_1,k)+C_2*ksi(x_2,k)+C_3*ksi(x_3,k));


def integrate_diff(k: int, left_bound: float, right_bound: float):     
    C_1 = 5/9;
    C_2 = 8/9;
    C_3 = 5/9;
    
    t_1 = -1*math.sqrt(15)/5;
    t_3 = math.sqrt(15)/5;

    x_1 = ((right_bound+left_bound)/2)+((right_bound-left_bound)/2)*t_1;
    x_2 = ((right_bound+left_bound)/2);
    x_3 = ((right_bound+left_bound)/2)+((right_bound-left_bound)/2)*t_3;

    return ((right_bound-left_bound)/2)*(C_1*diff(x_1,k)+C_2*diff(x_2,k)+C_3*diff(x_3,k));



def integrate_for_exact_solution(k: int):     
    
    first_step = (b-a)/10;

    current_step = first_step;

    left_bound = a;

    I=0;

    while left_bound+current_step<=b:  

        I_h = integrate_diff(k,left_bound,left_bound+current_step); 

        I_h_2 = integrate_diff(k,left_bound,left_bound+current_step/2) + integrate_diff(k,left_bound+current_step/2,left_bound+current_step);

        delta = (I_h_2-I_h)/(pow(2,5)-1);        

        if delta > max(0.0001, 0.0001*abs(I_h_2)): 
            current_step = current_step/2;
            continue;

        I=I+I_h;

        if current_step<1.0e+8:
            break;

        if delta <= (1/pow(2,5))*max(0.000001, 0.000001*abs(I_h_2)):
            left_bound = left_bound+current_step;
            current_step = current_step*2;

            if left_bound+current_step>b:
                current_step = b - left_bound;

            continue;
        
        else:
            left_bound = left_bound+current_step;
        
            if left_bound+current_step>b:
                current_step=b-left_bound;

            continue;

    return I;


first_step = (b-a)/10;

current_step = first_step;

left_bound = a;

k=1;

while True:
    
    I=0;

    while left_bound+current_step<=b:  

        I_h = integrate(k,left_bound,left_bound+current_step); 

        I_h_2 = integrate(k,left_bound,left_bound+current_step/2) + integrate(k,left_bound+current_step/2,left_bound+current_step);

        delta = (I_h_2-I_h)/(pow(2,5)-1);        

        if delta > max(0.0001, 0.0001*abs(I_h_2)): 
            current_step = current_step/2;
            continue;

        I=I+I_h;

        if current_step<1.0e+8:
            break;

        if delta <= (1/pow(2,5))*max(0.000001, 0.000001*abs(I_h_2)):
            left_bound = left_bound+current_step;
            current_step = current_step*2;

            if left_bound+current_step>b:
                current_step = b - left_bound;

            continue;
        
        else:
            left_bound = left_bound+current_step;
        
            if left_bound+current_step>b:
                current_step=b-left_bound;

            continue;

        

    if I>pow(epsilon,2):
        print('||u_',10*k,'-','u_',10*(k+1),'|| = %.16f' % (math.sqrt(I)));
        k = k+1;
        continue;
    else:
        print('||u_',10*k,'-','u_',10*(k+1),'|| = %.16f' % (math.sqrt(I)));
        print('\n');
        k = k+1;
        print('||u_%(order)d-u(x)|| = %(value).16f' % {"order": k*10, "value": math.sqrt(integrate_for_exact_solution(k))});

        break;

print('----------------------------------\n');









    
        


