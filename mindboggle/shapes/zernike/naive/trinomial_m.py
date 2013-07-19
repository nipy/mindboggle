'''VALID'''

def trinomial(i,j,k) :                             #function t=trinomial(i,j,k)
                                                   #% Computes the trinomial of
                                                   #% the three input arguments

    assert isinstance(i,int) and isinstance(j,int) and isinstance(k,int)
    assert i >= 0 and j >= 0 and k >= 0
    aux_1 = factorial(i+j+k)                       #aux_1=factorial(i+j+k);
    aux_2 = factorial(i)*factorial(j)*factorial(k) #aux_2=factorial(i)*factorial(j)*factorial(k);
    return float(aux_1)/float(aux_2)               #t = aux_1/aux_2;

