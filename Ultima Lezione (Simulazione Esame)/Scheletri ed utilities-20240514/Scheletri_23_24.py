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
 if   #to do
     print("Non è possibile applicare il metodo di bisezione \n")
     return None, None,None

 it = 0
 v_xk = []

 maxit = math.ceil(math.log((b - a) / tolx) / math.log(2))-1

 
 while :  #to do
    xk =  #to do
    v_xk.append(xk)
    it += 1
    fxk=fname(xk)
    if fxk==0:
      return xk, it, v_xk

     
    if sign(fa)*sign(fxk)>0:   
       # to do
    elif sign(fxk)*sign(fb)>0:    
       # to do

 
 return xk, it, v_xk


def falsi(fname, a, b, maxit, tolx,tolf):
 """
 Implementa il metodo di falsa posizione per il calcolo degli zeri di un'equazione non lineare.

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
 
if  : #to do
     print("Non è possibile applicare il metodo di falsa posizione \n")
     return None, None,None

 it = 0
 v_xk = []
 
 fxk=10

 
 while  : # to do
    xk = #to do
    v_xk.append(xk)
    it += 1
    fxk=fname(xk)
    if fxk==0:
      return xk, it, v_xk

    # 
    if sign(fa)*sign(fxk)>0:   
       #to do
        #to do
    elif sign(fxk)*sign(fb)>0:    
      #to do
      #to do

 
 return xk, it, v_xk


def corde(fname,m,x0,tolx,tolf,nmax):
 """
 Implementa il metodo delle corde per il calcolo degli zeri di un'equazione non lineare.

 Parametri:
  fname: La funzione da cui si vuole calcolare lo zero.
  m: coefficiente angolare della retta che rimane fisso per tutte le iterazioni
  tolx: La tolleranza di errore tra due iterati successivi
  tolf: tolleranza sul valore della funzione
  nmax: numero massimo di iterazione

 Restituisce:
  Lo zero approssimato della funzione, il numero di iterazioni e la lista degli iterati intermedi.
 """
        xk=[]
        fx0=#to do
        d=#to do
        x1=#to do
        fx1=fname(x1)
        xk.append(x1)
        it=1
        
        while  :
           x0= # to do
           fx0= #to do
           d= #to do
           '''
           #x1= ascissa del punto di intersezione tra  la retta che passa per il punto
           (xi,f(xi)) e ha pendenza uguale a m  e l'asse x
           '''
           x1=#to do  
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
        if : #to do
            print(" derivata prima nulla in x0")
            return None, None,None
        
        d=#to do
        x1=#to do
        
        fx1=fname(x1)
        xk.append(x1)
        it=1
        
        while #to do :
           x0= #to do
           fx0= #to do
           if #to do: #Se la derivata prima e' pià piccola della precisione di macchina stop
                print(" derivata prima nulla in x0")
                return None, None,None
           d=#to do
            
           x1=#to do
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
        fxm1=#to do
        fx0=#to do
        d=#to do
        x1=#to do;
        xk.append(x1)
        fx1=fname(x1);
        it=1
       
        while it<nmax and abs(fx1)>=tolf and abs(d)>=tolx*abs(x1):
            xm1=#to do
            x0=#to do
            fxm1=#to do)
            fx0=#to do 
            d=#to do
            x1=#to do
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
        fx0=#to do
        if #to do :
            print(" derivata prima nulla in x0")
            return None, None,None

        d=#to do
        x1=#to do
        
        fx1=#to do
        xk.append(x1)
        it=1
        
        while #to do :
           x0=#to do
           fx0=#to do
           if #to do: #Se la derivata prima e' pià piccola della precisione di macchina stop
                print(" derivata prima nulla in x0")
                return None, None,None
           d=#to do
            
           x1=#to do 
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
  if #to do:
    print("La matrice dello Jacobiano calcolata nell'iterato precedente non è a rango massimo")
    return None, None,None

  s = #to do
  # Aggiornamento della soluzione
  it = 1
  x1 =#to do
  fx1 = fun(x1)

  Xm = [np.linalg.norm(s, 1)/np.linalg.norm(x1,1)]

  while#to do:
    x0 =#to do
    it += 1
    matjac = jac(x0)
    if #to do:
            print("La matrice dello Jacobiano calcolata nell'iterato precedente non è a rango massimo")
            return None, None,None

   
    s =#to do

    # Aggiornamento della soluzione
    x1 = #to do
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
  if #to do:
    print("La matrice dello Jacobiano calcolata nell'iterato precedente non è a rango massimo")
    return None, None,None
  s = #to do
  # Aggiornamento della soluzione
  it = 1
  x1 = #to do
  fx1 = fun(x1)

  Xm = [np.linalg.norm(s, 1)/np.linalg.norm(x1,1)]

  while #to do:
    x0 = #to do
    it += 1
   
   
    if #to do:
        print("La matrice dello Jacobiano calcolata nell'iterato precedente non è a rango massimo")
        return None, None,None
    
     
    
    s = #to do

    # Aggiornamento della soluzione
    x1 =  #to do
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
  if #to do:
    print("La matrice dello Jacobiano calcolata nell'iterato precedente non è a rango massimo")
    return None,None,None

  s = #to do
  # Aggiornamento della soluzione
  it = 1
  x1 = x0 + s#to do
  fx1 = fun(x1)

  Xm = [np.linalg.norm(s, 1)/np.linalg.norm(x1,1)]
  update=10  #Numero di iterazioni durante le quali non si aggiorna la valutazione dello Jacobiano nell'iterato attuale
  while #to do:
    x0 =  #to do
    it += 1
    if it%update==0:   #Valuto la matrice di iterazione nel nuovo iterato ogni "update" iterazioni
        #to do
   
        if #to do == 0:
           print("La matrice dello Jacobiano calcolata nell'iterato precedente non è a rango massimo")
           return None,None,None
        else:
         
           s = #to do
    else:
          
           s = -#to do

    # Aggiornamento della soluzione
    x1 = #to do
    fx1 = fun(x1)
    Xm.append(np.linalg.norm(s, 1)/np.linalg.norm(x1,1))

  return x1, it, Xm


def my_newton_minimo(gradiente, Hess, x0, tolx, tolf, nmax):

  """
  DA UTILIZZARE NEL CASO IN CUI CALCOLATE DRIVATE PARZIALI PER GRADIENTE ED HESSIANO SENZA UTILIZZO DI SYMPY
  
  Funzione di newton-raphson per calcolare il minimo di una funzione in più variabili

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

  matHess = #to do
  if #to do:
    print("La matrice Hessiana calcolata nell'iterato precedente non è a rango massimo")
    return None, None, None
  grad_fx0= gradiente(x0)    
  s = #to do
  # Aggiornamento della soluzione
  it = 1
  x1 =#to do
  grad_fx1 = gradiente(x1)
  Xm = [np.linalg.norm(s, 1)]
  
  while#to do:
     
    x0 = #to do
    it += 1
    matHess = #to do
    grad_fx0=grad_fx1
     
    if #to do:
       
      print("La matrice Hessiana calcolata nell'iterato precedente non è a rango massimo")
      return None, None, None
      
 
    s = #to do
     
    # Aggiornamento della soluzione
    x1 = #to do

    #Calcolo del gradiente nel nuovo iterato
    grad_fx1  = gradiente(x1)
    print(np.linalg.norm(s, 1))
    Xm.append(np.linalg.norm(s, 1))

  return x1, it, Xm

def my_newton_minimo_MOD(gradiente, Hess, x0, tolx, tolf, nmax):

  """
  Funzione di newton-raphson per calcolare il minimo di una funzione in più variabili, modificato nel caso in cui si utilizzando sympy 
  per calcolare Gradiente ed Hessiano. Rispetto alla precedente versione cambia esclusivamente il modo di valutare il vettore gradiente e la matrice Hessiana in un punto 
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
   
  if #to do
    print("La matrice Hessiana calcolata nell'iterato precedente non è a rango massimo")
    return None, None, None
      
  s = #to do
  
  # Aggiornamento della soluzione
  it = 1
  x1 = #to do
  grad_fx1=np.array([gradiente[0](x1[0],x1[1]),gradiente[1](x1[0],x1[1])])
  Xm = [np.linalg.norm(s, 1)]
  
  while #to do:
     
    x0 = x1
    it += 1
    matHess = #to do
    grad_fx0=grad_fx1
      
    if np.linalg.det(matHess) == 0:
       
      print("La matrice Hessiana calcolata nell'iterato precedente non è a rango massimo")
      return None, None, None
      

    
    s = #to do
     
    # Aggiornamento della soluzione
    x1 =#to do
    #Aggiorno il gradiente per la prossima iterazione 
    grad_fx1=np.array([gradiente[0](x1[0],x1[1]),gradiente[1](x1[0],x1[1])])
    print(np.linalg.norm(s, 1))
    Xm.append(np.linalg.norm(s, 1))

  return x1, it, Xm

def jacobi(A,b,x0,toll,it_max):
    errore=1000
    d=np.diag(A)
    n=A.shape[0]
    invM=np.diag(1/d)
    E=#to do
    F=#to do
    N=#to do
    T=#to do
    autovalori=np.linalg.eigvals(T)
    raggiospettrale=#to do
    print("raggio spettrale jacobi", raggiospettrale)
    it=0
    
    er_vet=[]
    while #to do:
        x=#to do
        errore=np.linalg.norm(x-x0)/np.linalg.norm(x)
        er_vet.append(errore)
        x0=x.copy()
        it=it+1
    return x,it,er_vet

def gauss_seidel(A,b,x0,toll,it_max):
    errore=1000
    d=np.diag(A)
    D=#to do
    E=#to do
    F=#to do
    M=#to do
    N=#to do
    T=#to do
    autovalori=np.linalg.eigvals(T)
    raggiospettrale=#to do
    print("raggio spettrale Gauss-Seidel ",raggiospettrale)
    it=0
    er_vet=[]
    while #to do:
        temp=#to do
        x= #to do 
        errore=np.linalg.norm(x-x0)/np.linalg.norm(x)
        er_vet.append(errore)
        x0=x.copy()
        it=it+1
    return x,it,er_vet

def gauss_seidel_sor(A,b,x0,toll,it_max,omega):
    errore=1000
    d=#to donp.diag(A)
    D=#to do
    Dinv=#to do
    E=#to do
    F=#to do
    Momega=D+omega*E
    Nomega=(1-omega)*D-omega*F
    T=#to do
    autovalori=np.linalg.eigvals(T)
    raggiospettrale=np.max(np.abs(autovalori))
    print("raggio spettrale Gauss-Seidel SOR ", raggiospettrale)
    
    M=D+E
    N=-F
    it=0
    xold=x0.copy()
    xnew=x0.copy()
    er_vet=[]
    while #to do
        temp=#to do
        xtilde#to do
        xnew=#to do
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
    p =  
    it = 0
    nb=np.linalg.norm(b)
    errore=np.linalg.norm(r)/nb
    vec_sol=[]
    vec_sol.append(x)
    vet_r=[]
    vet_r.append(errore)
     
# utilizzare il metodo del gradiente per trovare la soluzione
    while #to do:
        it=it+1
        Ap= #to do
       
        alpha = #to do
                
        x =  
        
         
        vec_sol.append(x)
        r=r+alpha*Ap
        errore=np.linalg.norm(r)/nb
        vet_r.append(errore)
        p =#to do
        
     
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
    G= 
     
    f= 
    
    L= 
    U=L.T
        
   
    z=
    x=
    
    return x
    
def qrLS(A,b):
#Risolve un sistema sovradeterminato con il metodo QR-LS
    n=A.shape[1]  # numero di colonne di A
    Q,R=spLin.qr(A)
    h=#to do
    x,flag=SolveTriangular.Usolve( #to do)
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
    d=#to do
    d1=#to do
    s1=#to do
    #Risolve il sistema diagonale di dimensione kxk avene come matrice dei coefficienti la matrice Sigma
    c=#to do
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
       xzeri=np.append(#to do)
    
    num=#to do
    den=#to do
    
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
        p=#to do
        L[:,j]=#to do
    
    
     return L@y