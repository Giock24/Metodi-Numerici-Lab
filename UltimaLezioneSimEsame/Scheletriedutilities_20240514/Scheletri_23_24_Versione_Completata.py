import numpy as np
import math
import matplotlib.pyplot as plt
import sympy as sym
from sympy.utilities.lambdify import lambdify

def sign(x):
  """
  Funzione segno che restituisce 1 se x è positivo, 0 se x è zero e -1 se x è negativo.
  """
  return math.copysign(1, x)

def metodo_bisezione(fname, a, b, tolx,tolf):
    """
    Implementa il metodo di bisezione per il calcolo degli zeri di un'equazione non lineare.
    Parametri:
    f: La funzione da cui si vuole calcolare lo zero.
    a: L'estremo sinistro dell'intervallo di ricerca.
    b: L'estremo destro dell'intervallo di ricerca.
    tol: La tolleranza di errore.
    Restituisce:
    Lo zero approssimato della funzione, il numero di iterazioni e la lista di valori intermedi.
    """
    fa=fname(a);
    fb=fname(b);
    if  fa * fb > 0:
        print("Non è possibile applicare il metodo di bisezione \n")
        return None, None,None

    it = 0
    v_xk = []

    maxit = math.ceil(math.log((b - a) / tolx) / math.log(2))-1

 
    while it < maxit and np.abs(b - a) > tolx:  
        xk =  (b+a)/2
        v_xk.append(xk)
        it += 1
        fxk=fname(xk)
        if fxk==0:
            return xk, it, v_xk

     
        if sign(fa)*sign(fxk)>0:   
            a = xk
            fa = fxk
        elif sign(fxk)*sign(fb)>0:    
            b = xk
            fb = fxk

 
    return xk, it, v_xk




def falsi(fname, a, b, maxit, tolx,tolf):
    """
    Implementa il metodo di falsa posizione per il calcolo degli zeri di un'equazione non lineare
    Parametri:
        f: La funzione da cui si vuole calcolare lo zero.
        a: L'estremo sinistro dell'intervallo di ricerca.
        b: L'estremo destro dell'intervallo di ricerca.
        tol: La tolleranza di errore
    Restituisce:
        Lo zero approssimato della funzione, il numero di iterazioni e la lista di valori intermedi.
    """
    fa=fname(a);
    fb=fname(b);
 
    if  fa * fb > 0: 
        print("Non è possibile applicare il metodo di falsa posizione \n")
        return None, None,None

    it = 0
    v_xk = []
 
    fxk=10

 
    while  it < maxit and np.abs(b - a) > tolx and np.abs(fxk) > tolf: 
        xk = a - fa * (b - a)/(fb - fa)
        v_xk.append(xk)
        it += 1
        fxk=fname(xk)
        if fxk==0:
            return xk, it, v_xk
    
        if sign(fa)*sign(fxk)>0:   
            a = xk
            fa = fxk
        elif sign(fxk)*sign(fb)>0:    
            b = xk
            fb = fxk

 
    return xk, it, v_xk


def corde(fname,m,x0,tolx,tolf,nmax):
    """
    Implementa il metodo delle corde per il calcolo degli zeri di un'equazione non linear
    Parametri:
        fname: La funzione da cui si vuole calcolare lo zero.
        m: coefficiente angolare della retta che rimane fisso per tutte le iterazioni
        tolx: La tolleranza di errore tra due iterati successivi
        tolf: tolleranza sul valore della funzione
        nmax: numero massimo di iterazio
    Restituisce:
        Lo zero approssimato della funzione, il numero di iterazioni e la lista degli iterati intermedi.
    """
    xk=[]
    fx0=fname(x0)
    d=fx0/m
    x1=x0-d
    fx1=fname(x1)
    xk.append(x1)
    it=1
    
    # errore relativo x1 = np.abs(x1 - x0)/np.abs(x1)
    while it < nmax and np.abs(x1 - x0)/np.abs(x1) > tolx and np.abs(fx1) > tolf:
        x0=x1
        fx0=fname(x0)
        d=fx0/m
        """
        #x1= ascissa del punto di intersezione tra  la retta che passa per il punto
        (xi,f(xi)) e ha pendenza uguale a m  e l'asse x
        """
        x1=x0-d
        fx1=fname(x1)
        it=it+1
        
        xk.append(x1)
          
        if it==nmax:
            print('raggiunto massimo numero di iterazioni \n')
            
        
    return x1,it,xk

def newton(fname,fpname,x0,tolx,tolf,nmax):
    """
    Implementa il metodo di Newton per il calcolo degli zeri di un'equazione non lineare.
    Parametri:
     fname: La funzione di cui si vuole calcolare lo zero.
     fpname: La derivata prima della funzione di  cui si vuole calcolare lo zero.
     x0: iterato iniziale
     tolx: La tolleranza di errore tra due iterati successivi
     tolf: tolleranza sul valore della funzione
     nmax: numero massimo di iterazione
    Restituisce:
     Lo zero approssimato della funzione, il numero di iterazioni e la lista degli iterati intermedi.
    """ 
    xk=[]
    fx0=fname(x0)
    if np.abs(fpname(x0)) <= np.spacing(1): #Se la derivata prima e' pià piccola della precisione di macchina stop
        print(" derivata prima nulla in x0")
        return None, None,None
    
    d=fx0/fpname(x0)
    x1=x0-d
    
    fx1=fname(x1)
    xk.append(x1)
    it=1
    
    while it < nmax and np.abs(x1 - x0)/np.abs(x1) > tolx and np.abs(fx1) > tolf:
        x0=x1
        fx0=fname(x0)
        if np.abs(fpname(x0)) <= np.spacing(1): #Se la derivata prima e' pià piccola della precisione di macchina stop
             print(" derivata prima nulla in x0")
             return None, None,None
        d=fx0/fpname(x0)
         
        x1=x0-d
        fx1=fname(x1)
        it=it+1
        
        xk.append(x1)
          
        if it==nmax:
            print('raggiunto massimo numero di iterazioni \n')
            
        
    return x1,it,xk

def secanti(fname,xm1,x0,tolx,tolf,nmax):
    """
    Implementa il metodo delle secanti per il calcolo degli zeri di un'equazione non lineare.
    Parametri:
    fname: La funzione di cui si vuole calcolare lo zero.
    xm1, x0: primi due iterati
    tolx: La tolleranza di errore tra due iterati successivi
    tolf: tolleranza sul valore della funzione
    nmax: numero massimo di iterazione
    Restituisce:
    Lo zero approssimato della funzione, il numero di iterazioni e la lista degli iterati intermedi.
    """
    xk=[]
    fxm1=fname(xm1)
    fx0=fname(x0)
    d=fx0*(x0-xm1)/(fx0-fxm1) # ATTENZIONE: questo perchè xm1 è considerato come primo iterato mentre x0 il secondo
    # guarda nei commenti come ha scritto i parametri
    x1=x0-d
    xk.append(x1)
    fx1=fname(x1)
    it=1
       
    while it<nmax and abs(fx1)>=tolf and abs(d)>=tolx*abs(x1):
        xm1=x0
        x0=x1
        fxm1=fname(xm1)
        fx0=fname(x0)
        d=fx0*(x0-xm1)/(fx0-fxm1)
        x1=x0-d
        fx1=fname(x1)
        xk.append(x1);
        it=it+1;
           
       
        if it==nmax:
            print('Secanti: raggiunto massimo numero di iterazioni \n')
        
    return x1,it,xk

def newton_mod(fname,fpname,m,x0,tolx,tolf,nmax):
    """
    Implementa il metodo di Newton modificato da utilizzato per il calcolo degli zeri di un'equazione non lineare
    nel caso di zeri multipli.
    Parametri:
    fname: La funzione di cui si vuole calcolare lo zero.
    fpname: La derivata prima della funzione di  cui si vuole calcolare lo zero.
    m: molteplicità della radice
    x0: iterato iniziale
    tolx: La tolleranza di errore tra due iterati successivi
    tolf: tolleranza sul valore della funzione
    nmax: numero massimo di iterazione
    Restituisce:
    Lo zero approssimato della funzione, il numero di iterazioni e la lista degli iterati intermedi.
    """ 
 
    xk=[]
    fx0=fname(x0)
    
    if np.abs(fpname(x0)) <= np.spacing(1) :
        print(" derivata prima nulla in x0")
        return None, None,None
    
    d=fx0/fpname(x0)
    x1= x0 - m * d
    
    fx1=fname(x1)
    xk.append(x1)
    it=1
        
    while it <= nmax and np.abs(x1 - x0)/np.abs(x1) > tolx and np.abs(fx1) > tolf :
        x0=x1
        fx0=fname(x0)
        
        if np.abs(fpname(x0)) <= np.spacing(1): #Se la derivata prima e' pià piccola della precisione di macchina stop
            print(" derivata prima nulla in x0")
            return None, None,None
        
        d=fx0/fpname(x0) 
        x1= x0 - m * d
        fx1=fname(x1)
        it=it+1
     
        xk.append(x1)
      
    if it==nmax:
        print('raggiunto massimo numero di iterazioni \n')
        
    
    return x1,it,xk
    
def stima_ordine(xk,iterazioni):
    #Vedi dispensa allegata per la spiegazione
    k=iterazioni-4
    p=np.log(abs(xk[k+2]-xk[k+3])/abs(xk[k+1]-xk[k+2]))/np.log(abs(xk[k+1]-xk[k+2])/abs(xk[k]-xk[k+1]));

    ordine=p
    return ordine

def my_newtonSys(fun, jac, x0, tolx, tolf, nmax):

    """
    Funzione per la risoluzione del sistema F(x)=0
    mediante il metodo di Newton.
    Parametri
    ----------
    fun : funzione vettoriale contenente ciascuna equazione non lineare del sistema.
    jac : funzione che calcola la matrice Jacobiana della funzione vettoriale.
    x0 : array
      Vettore contenente l'approssimazione iniziale della soluzione.
    tolx : float
      Parametro di tolleranza per l'errore assoluto.
    tolf : float
      Parametro di tolleranza per l'errore relativo.
    nmax : int
      Numero massimo di iterazioni.
    Restituisce
    -------
    x : array
      Vettore soluzione del sistema (o equazione) non lineare.
    it : int
      Numero di iterazioni fatte per ottenere l'approssimazione desiderata.
    Xm : array
      Vettore contenente la norma dell'errore relativo tra due iterati successivi.
    """

    matjac = jac(x0)
    if np.linalg.det(matjac) == 0:
        print("La matrice dello Jacobiano calcolata nell'iterato precedente non è a rango massimo")
        return None, None,None

    s = np.linalg.solve(matjac, -fun(x0))
    # Aggiornamento della soluzione
    it = 1
    x1 = x0 + s
    fx1 = fun(x1)
    Xm = [np.linalg.norm(s, 1)/np.linalg.norm(x1,1)]

    while it <= nmax and np.linalg.norm(x1 - x0, 1) > tolx and np.linalg.norm(fx1, 1) > tolf:
        x0 = x1
        it += 1
        matjac = jac(x0)
        if np.linalg.det(matjac) == 0:
            print("La matrice dello Jacobiano calcolata nell'iterato precedente non è a rango massimo")
            return None, None,None

   
        s = np.linalg.solve(matjac, -fun(x0))

        # Aggiornamento della soluzione
        x1 = x0 + s
        fx1 = fun(x1)
        Xm.append(np.linalg.norm(s, 1)/np.linalg.norm(x1,1))

    return x1, it, Xm


def my_newtonSys_corde(fun, jac, x0, tolx, tolf, nmax):

    """
    Funzione per la risoluzione del sistema f(x)=0
    mediante il metodo di Newton, con variante delle corde, in cui lo Jacobiano non viene calcolato
    ad ogni iterazione, ma rimane fisso, calcolato nell'iterato iniziale x0.
  
    Parametri
    ----------
    fun : funzione vettoriale contenente ciascuna equazione non lineare del sistema.
    jac : funzione che calcola la matrice Jacobiana della funzione vettoriale.
    x0 : array
        Vettore contenente l'approssimazione iniziale della soluzione.
    tolx : float
        Parametro di tolleranza per l'errore tra due soluzioni successive.
    tolf : float
        Parametro di tolleranza sul valore della funzione.
    nmax : int
        Numero massimo di iterazioni.
    
    Restituisce
    -------
    x : array
        Vettore soluzione del sistema (o equazione) non lineare.
    it : int
        Numero di iterazioni fatte per ottenere l'approssimazione desiderata.
    Xm : array
        Vettore contenente la norma dell'errore relativo tra due iterati successivi.
    """

    matjac = jac(x0)   
    if np.linalg.det(matjac) == 0:
        print("La matrice dello Jacobiano calcolata nell'iterato precedente non è a rango massimo")
        return None, None,None
    s = np.linalg.solve(matjac, -fun(x0))
    # Aggiornamento della soluzione
    it = 1
    x1 = x0 + s
    fx1 = fun(x1)
    Xm = [np.linalg.norm(s, 1)/np.linalg.norm(x1,1)]

    while it <= nmax and np.linalg.norm(x1 - x0, 1) > tolx and np.linalg.norm(fx1, 1) > tolf:
        x0 = x1
        it += 1
        
        
        if np.linalg.det(matjac) == 0:
            print("La matrice dello Jacobiano calcolata nell'iterato precedente non è a rango massimo")
            return None, None,None
        
         
        
        s = np.linalg.solve(matjac, -fun(x0))
    
        # Aggiornamento della soluzione
        x1 =  x0 + s
        fx1 = fun(x1)
        Xm.append(np.linalg.norm(s, 1)/np.linalg.norm(x1,1))

    return x1, it, Xm

def my_newtonSys_sham(fun, jac, x0, tolx, tolf, nmax):

    """
  Funzione per la risoluzione del sistema f(x)=0
  mediante il metodo di Newton, con variante delle shamanski, in cui lo Jacobiano viene
  aggiornato ogni un tot di iterazioni, deciso dall'utente.

  Parametri
  ----------
  fun : funzione vettoriale contenente ciascuna equazione non lineare del sistema.
  jac : funzione che calcola la matrice Jacobiana della funzione vettoriale.
  x0 : array
    Vettore contenente l'approssimazione iniziale della soluzione.
  tolx : float
    Parametro di tolleranza per l'errore tra due soluzioni successive.
  tolf : float
    Parametro di tolleranza sul valore della funzione.
  nmax : int
    Numero massimo di iterazioni.

  Restituisce
  -------
  x : array
    Vettore soluzione del sistema (o equazione) non lineare.
  it : int
    Numero di iterazioni fatte per ottenere l'approssimazione desiderata.
  Xm : array
      Vettore contenente la norma dell'errore relativo tra due iterati successivi.
    """

    matjac = jac(x0)
    if np.linalg.det(matjac) == 0:
        print("La matrice dello Jacobiano calcolata nell'iterato precedente non è a rango massimo")
        return None,None,None

    s = np.linalg.solve(matjac, -fun(x0))
    # Aggiornamento della soluzione
    it = 1
    x1 = x0 + s
    fx1 = fun(x1)

    Xm = [np.linalg.norm(s, 1)/np.linalg.norm(x1,1)]
    update=10  #Numero di iterazioni durante le quali non si aggiorna la valutazione dello Jacobiano nell'iterato attuale
    while it <= nmax and np.linalg.norm(x1 - x0, 1) > tolx and np.linalg.norm(fx1, 1) > tolf:
        x0 = x1
        it += 1
        if it%update==0:   #Valuto la matrice di iterazione nel nuovo iterato ogni "update" iterazioni
            matjac = jac(x0)
       
            if np.linalg.det(matjac) == 0:
                print("La matrice dello Jacobiano calcolata nell'iterato precedente non è a rango massimo")
                return None,None,None
            else:
             
               s = np.linalg.solve(matjac, -fun(x0))
        else:
              
            s = np.linalg.solve(matjac, -fun(x0))

        # Aggiornamento della soluzione
        x1 = x0 + s
        fx1 = fun(x1)
        Xm.append(np.linalg.norm(s, 1)/np.linalg.norm(x1,1))

    return x1, it, Xm


def my_newton_minimo(gradiente, Hess, x0, tolx, tolf, nmax):

    """
    DA UTILIZZARE NEL CASO IN CUI CALCOLATE DRIVATE PARZIALI PER GRADIENTE ED HESSIANO SENZA UTILIZZO DI SYMPY
    
    Funzione di newton-raphson per calcolare il minimo di una funzione in più variabili
    Parametri
    ----------
    gradiente : 
      Nome della funzione che calcola il gradiente della funzione non lineare.
    Hess :  
      Nome della funzione che calcola la matrice Hessiana della funzione non lineare.
    x0 : array
      Vettore contenente l'approssimazione iniziale della soluzione.
    tolx : float
      Parametro di tolleranza per l'errore assoluto.
    tolf : float
      Parametro di tolleranza per l'errore relativo.
    nmax : int
      Numero massimo di iterazioni.
    Restituisce
    -------
    x : array
        Vettore soluzione del sistema (o equazione) non lineare.
    it : int
        Numero di iterazioni fatte per ottenere l'approssimazione desiderata.
    Xm : array
        Vettore contenente la norma del passo ad ogni iterazione.
    """

    matHess = Hess(x0)
    if np.linalg.det(matHess) == 0:
        print("La matrice Hessiana calcolata nell'iterato precedente non è a rango massimo")
        return None, None, None
    grad_fx0= gradiente(x0)    
    s = np.linalg.solve(matHess, -grad_fx0)
    # Aggiornamento della soluzione
    it = 1
    x1 = x0 + s
    grad_fx1 = gradiente(x1)
    Xm = [np.linalg.norm(s, 1)]
  
    while it <= nmax and np.linalg.norm(x1 - x0, 1)/np.linalg.norm(x0, 1) > tolx and np.linalg.norm(grad_fx1) > tolf:
     
        x0 = x1
        it += 1
        matHess = Hess(x0)
        grad_fx0=grad_fx1
         
        if np.linalg.det(matHess) == 0:
            print("La matrice Hessiana calcolata nell'iterato precedente non è a rango massimo")
            return None, None, None
          
     
        s = np.linalg.solve(matHess, -grad_fx0)
         
        # Aggiornamento della soluzione
        x1 = x0 + s
    
        #Calcolo del gradiente nel nuovo iterato
        grad_fx1  = gradiente(x1)
        print(np.linalg.norm(s, 1))
        Xm.append(np.linalg.norm(s, 1))

    return x1, it, Xm

def my_newton_minimo_MOD(gradiente, Hess, x0, tolx, tolf, nmax):

    """
    Funzione di newton-raphson per calcolare il minimo di una funzione in più variabili, 
    modificato nel caso in cui si utilizzando sympy per calcolare Gradiente ed Hessiano. 
    Rispetto alla precedente versione cambia esclusivamente il modo di valutare 
    il vettore gradiente e la matrice Hessiana in un punto 
    Parametri
    ----------
    fun : 
    Nome della funzione che calcola il gradiente della funzione non lineare.
    Hess :  
    Nome della funzione che calcola la matrice Hessiana della funzione non lineare.
    x0 : array
    Vettore contenente l'approssimazione iniziale della soluzione.
    tolx : float
    Parametro di tolleranza per l'errore assoluto.
    tolf : float
    Parametro di tolleranza per l'errore relativo.
    nmax : int
    Numero massimo di iterazioni.
    Restituisce
    -------
    x : array
    Vettore soluzione del sistema (o equazione) non lineare.
    it : int
    Numero di iterazioni fatte per ottenere l'approssimazione desiderata.
    Xm : array
    Vettore contenente la norma del passo ad ogni iterazione.
    """

    
    matHess = np.array([[Hess[0, 0](x0[0], x0[1]), Hess[0, 1](x0[0], x0[1])],
                      [Hess[1, 0](x0[0], x0[1]), Hess[1, 1](x0[0], x0[1])]])
 

    gradiente_x0=np.array([gradiente[0](x0[0], x0[1]),gradiente[1](x0[0], x0[1])])
   
    if np.linalg.det(matHess) == 0:
        print("La matrice Hessiana calcolata nell'iterato precedente non è a rango massimo")
        return None, None, None
      
    s = np.linalg.solve(matHess, -gradiente_x0)
  
    # Aggiornamento della soluzione
    it = 1
    x1 = x0 + s
    grad_fx1=np.array([gradiente[0](x1[0],x1[1]),gradiente[1](x1[0],x1[1])])
    Xm = [np.linalg.norm(s, 1)]
  
    while it < nmax and np.linalg.norm(x1 - x0, 1)/np.linalg.norm(x1, 1) > tolx and np.linalg.norm(grad_fx0) > tolf:
     
        x0 = x1
        it += 1
        matHess = np.array([[Hess[0, 0](x0[0], x0[1]), Hess[0, 1](x0[0], x0[1])],
                      [Hess[1, 0](x0[0], x0[1]), Hess[1, 1](x0[0], x0[1])]])
        grad_fx0=grad_fx1
          
        if np.linalg.det(matHess) == 0:
            print("La matrice Hessiana calcolata nell'iterato precedente non è a rango massimo")
            return None, None, None
          
        
        s = np.linalg.solve(matHess, -grad_fx0)
         
        # Aggiornamento della soluzione
        x1 = x0 + s
        #Aggiorno il gradiente per la prossima iterazione 
        grad_fx1=np.array([gradiente[0](x1[0],x1[1]),gradiente[1](x1[0],x1[1])])
        print(np.linalg.norm(s, 1))
        Xm.append(np.linalg.norm(s, 1))

    return x1, it, Xm
#---------------------------------------------
# Metodi diretti aggiunti tra gli schemi
# originariamente non ci sono
import SolveTriangular
from scipy.linalg import lu

def my_LU(A,b):
    PT, L, U = lu(A)
    P = PT.T
    x = np.zeros_like(b)
    
    y, flag = SolveTriangular.Lsolve(L, P@b)
    if (flag == 1):
        return x
    else:
        x, flag = SolveTriangular.Usolve(U,y)
    
    return x


from scipy.linalg import qr

def my_QRsolve(A, b):
    Q, R = qr(A)
    y = Q.T@b # puoi scriverlo così perchè Q è una matrice ortogonale
    x, flag = SolveTriangular.Usolve(R, y)
    if (flag == 1):
        return np.zeros_like(b)
    
    return x

from scipy.linalg import cholesky

def my_cholesky(A,b):
    L = cholesky(A)
    LT = L.T
    
    y, flag = SolveTriangular.Lsolve(L, b)
    if (flag == 1):
        return x
    else:
        x, flag = SolveTriangular.Usolve(LT,y)
    
    return x
#---------------------------------------------
def jacobi(A,b,x0,toll,it_max):
    errore=1000
    d=np.diag(A)
    n=A.shape[0]
    invM=np.diag(1/d)
    E=np.tril(A, -1)
    F=np.triu(A,+1)
    N=-(E+F)
    T=invM@N
    autovalori=np.linalg.eigvals(T)
    raggiospettrale=np.sqrt(np.max(autovalori))
    print("raggio spettrale jacobi", raggiospettrale)
    it=0
    
    er_vet=[]
    while it <= it_max and errore > toll:
        x=T@x0+invM@b
        errore=np.linalg.norm(x-x0)/np.linalg.norm(x)
        er_vet.append(errore)
        x0=x.copy()
        it=it+1
    return x,it,er_vet

def gauss_seidel(A,b,x0,toll,it_max):
    errore=1000
    d=np.diag(A)
    D=np.diag(d)
    E=np.tril(A, -1)
    F=np.triu(A, +1)
    M=(E+D)
    N=-F
    T=np.linalg.inv(M)@N
    autovalori=np.linalg.eigvals(T)
    raggiospettrale=np.max(np.abs(autovalori))
    print("raggio spettrale Gauss-Seidel ",raggiospettrale)
    it=0
    er_vet=[]
    while it < it_max and errore > toll:
        #temp=#to do
        x=T@x0 + np.linalg.inv(M)@b
        errore=np.linalg.norm(x-x0)/np.linalg.norm(x)
        er_vet.append(errore)
        x0=x.copy()
        it=it+1
    return x,it,er_vet

def gauss_seidel_sor(A,b,x0,toll,it_max,omega):
    errore=1000
    d=np.diag(A)
    D=np.diag(d)
    Dinv=np.linalg.inv(D)
    E=np.tril(A,-1)
    F=np.triu(A,+1)
    Momega=D+omega*E
    Nomega=(1-omega)*D-omega*F
    T=np.linalg.inv(Momega)@Nomega
    autovalori=np.linalg.eigvals(T)
    raggiospettrale=np.max(np.abs(autovalori))
    print("raggio spettrale Gauss-Seidel SOR ", raggiospettrale)
    
    M=D+E
    N=-F
    it=0
    xold=x0.copy()
    xnew=x0.copy()
    er_vet=[]
    while it < it_max and errore > toll:
        #temp=#to do
        xtilde = (np.linalg.inv(M)@N@xold + np.linalg.inv(M)@b) - xold
        xnew= xold + xtilde*omega
        errore=np.linalg.norm(xnew-xold)/np.linalg.norm(xnew)
        er_vet.append(errore)
        xold=xnew.copy()
        it=it+1
    return xnew,it,er_vet

def steepestdescent(A,b,x0,itmax,tol):
 
    n,m=A.shape
    if n!=m:
        print("Matrice non quadrata")
        return [],[]
    
    
   # inizializzare le variabili necessarie
    x = x0

     
    r = A@x-b
    p = -r
    it = 0
    nb=np.linalg.norm(b)
    errore=np.linalg.norm(r)/nb
    vec_sol=[]
    vec_sol.append(x)
    vet_r=[]
    vet_r.append(errore)
     
    # utilizzare il metodo del gradiente per trovare la soluzione
    while it < itmax and errore > tol:
        it=it+1
        Ap=A@p
       
        alpha = (r.T@r)/((A@r).T@r)
        x = x+alpha*p
         
        vec_sol.append(x)
        r=r+alpha*Ap
        errore=np.linalg.norm(r)/nb
        vet_r.append(errore)
        p = -A@x + b
        
     
    return x,vet_r,vec_sol,it


def conjugate_gradient(A,b,x0,itmax,tol):
    n,m=A.shape
    if n!=m:
        print("Matrice non quadrata")
        return [],[]
    
    
   # inizializzare le variabili necessarie
    x = x0
    
    r = A@x-b
    p = -r
    it = 0
    nb=np.linalg.norm(b)
    errore=np.linalg.norm(r)/nb
    vec_sol=[]
    vec_sol.append(x0)
    vet_r=[]
    vet_r.append(errore)
# utilizzare il metodo del gradiente coniugato per calcolare la soluzione
    while errore >= tol and it< itmax:
        it=it+1
        Ap=#to do
        alpha = -#to do
        x =#to do
        vec_sol.append(x)
        rtr_old=r.T@r
        r=r+alpha*Ap
        gamma= 
        errore=np.linalg.norm(r)/nb
        vet_r.append(errore)
        p =  #to do
   
    
    return x,vet_r,vec_sol,it

def eqnorm(A,b):
#Risolve un sistema sovradeterminato con il metodo delle equazioni normali
    G= A.T@A
     
    f= A.T@b
    
    L= cholesky(G,lower=True)
    U=L.T
        
   
    z= SolveTriangular.Lsolve(L,f)
    x= SolveTriangular.Usolve(U,z)
    
    return x
    
def qrLS(A,b):
#Risolve un sistema sovradeterminato con il metodo QR-LS
    n=A.shape[1]  # numero di colonne di A
    Q,R=spLin.qr(A)
    h=Q.T@b
    x,flag=SolveTriangular.Usolve(R[:n,:],h[:n])
    residuo=np.linalg.norm(h[n:])**2
    return x,residuo

def SVDLS(A,b):
    #Risolve un sistema sovradeterminato con il metodo SVD-LS
    m,n=A.shape  #numero di righe e  numero di colonne di A
    U,s,VT=spLin.svd(A)  #Attenzione : Restituisce U, il numpy-array 1d che contiene la diagonale della matrice Sigma e VT=VTrasposta)
    #Quindi 
    V=VT.T
    thresh=np.spacing(1)*m*s[0] ##Calcolo del rango della matrice, numero dei valori singolari maggiori di una soglia
    k=np.count_nonzero(s>thresh)
    print("rango=",k)
    d=U.T@b
    d1=d[:k].reshape(k,1)
    s1=s[:k].reshape(k,1)
    #Risolve il sistema diagonale di dimensione kxk avene come matrice dei coefficienti la matrice Sigma
    c=d1/s1
    x=V[:,:k]@c 
    residuo=np.linalg.norm(d[k:])**2
    return x,residuo

def plagr(xnodi,j):
    """
    Restituisce i coefficienti del j-esimo pol di
    Lagrange associato ai punti del vettore xnodi
    """
    xzeri=np.zeros_like(xnodi)
    n=xnodi.size
    if j==0:
        xzeri=xnodi[1:n]
    else:
        xzeri=np.append(xnodi[0:j],xnodi[j+1:n])
    
    num=np.poly(xzeri)
    den=np.polyval(num, xnodi[j])
    
    p=num/den
    
    return p

def InterpL(x, y, xx):
     """"
        %funzione che determina in un insieme di punti il valore del polinomio
        %interpolante ottenuto dalla formula di Lagrange.
        % DATI INPUT
        %  x  vettore con i nodi dell'interpolazione
        %  f  vettore con i valori dei nodi 
        %  xx vettore con i punti in cui si vuole calcolare il polinomio
        % DATI OUTPUT
        %  y vettore contenente i valori assunti dal polinomio interpolante
        %
     """
    n=x.size
    m=xx.size
    L=np.zeros((m,n))
    for j in range(n):
        p=plagr(x,j)
        L[:,j]=np.polyval(p,xx)
    return L@y
'''
