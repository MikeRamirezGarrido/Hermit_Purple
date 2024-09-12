!Understanding Quantum Mechanics: From analytical to numerical analysis
!Hydrogen atom's electron probability densisty and radiuses
!Laguerre RungeKutta
     Function Lafirst(MLa)
             double precision, intent(IN) :: MLa
             double precision :: Lafirst
             Lafirst = MLa
     End Function
         

     Function Lasec(n,l,x,La,MLa)
             double precision , intent(IN) :: n,l
             double precision, intent(IN) :: x, La, MLa
             double precision :: Lasec
             Lasec = -((2*l+2-x)*MLa+(n-l-1)*La)/(x)
     End Function

!Legendre RungeKutta    
     Function Lefirst(QLe)
             double precision, intent(IN) :: QLe
             double precision :: Lefirst
             Lefirst = QLe
     End Function
      
     Function Lesec(m,l,x,Le,QLe)
             double precision , intent(IN) :: m, l
             double precision, intent(IN) :: x, Le, QLe
             double precision :: Lesec, r
             r = ((m**2)/(1-x**2))-l*(l+1)
             Lesec = ((2*x)*QLe+r*Le)/(1-x**2)
     End Function
      
!Auxiliar coefficient function      
     Function why(m,l)
             double precision , intent(IN) :: m, l
             double precision :: why, pi
             pi = 3.14159265358979323D0
             why = SQRT(((2*l+1)/(4*pi))*(gamma(l-m+1)/gamma(l+m+1)))
     End Function
      
!Auxiliar function     
     Function casio(n,l,x)
             double precision, intent(IN) :: n, l, x
             double precision :: a0, e, casio
             e = 2.7182818D0
             a0 =  0.0005291772109D0
             casio = (2/(n**2))*(SQRT((gamma(n-l)/(gamma(n+l+1)*a0**3))))*((x)**l)*e**(-0.5*x)
     End Function

     Program HydrogenAtom
             double precision :: n, m, l
             double precision :: Lasec, Lafirst
             double precision :: ak1, ak2, ak3, ak4, x0, y0, w0
             double precision :: h, i, cl, ya, xa, wa
             double precision :: am1, am2, am3, am4, ET
             double precision ::  ek1, ek2, ek3, ek4, ex0, ey0, ew0
             double precision :: em1, em2, em3, em4, eh, ei
             double precision :: eya, exa, ewa, ecl, a0
             double precision :: Lesec, Lefirst, pi, T
             double precision :: why, casio, rho, rhomax, xmax
             double precision :: xmaxError, riemannError, Low, Tlow, r1, r2
!Data output instruntions in line 100
!Divided into two files for size management
!            Open(4,file='data1.txt')
!            Open(5,file='data2.txt')

             Write(*,*) "Understanding Quantum Mechanics: From analytical to numerical analysis"
             Write(*,*) "Hydrogen atom's electron probability densisty and radiuses"
             pi = 3.14159265D0
!Bohr radius x10^7
             a0 = 0.0005291772109D0
             x0 = 0.001D0
             ex0 = -0.99999D0
!Write here
!Write quantic numbers
             n = 4.0D0
             l = 2.0D0
             m = 1.0D0
             
             rhomax = 0
             riemannError = 0
             xmaxError = 0
             r1 = 0
             r2 = 0
             ET = 0
             Low = 0
             Tlow = 0
             T = 0
             
!Write intial value L (Laguerre_{n-l-1}^{2l+1}) at x = 0.001       
             y0 = 5.999D0
!Write intial value L' (Derivative Laguerre_{n-l-1}^{2l+1}) at x = 0.001
             w0 = -1.0D0
         
             xa = x0
             ya = y0
             wa = w0
             cl = 0

!If you want the data output change stepsize h to 0.01 and the number of repetitions 
!to around 2000 for efficiency. Uncomment lines 64, 65, 181 and 182.           
2            IF(cl .LE. 15155) THEN
             h = 0.001D0
             
                     If(cl .EQ. 0) THEN
                     h = 0.0D0
                     End If

             i = h/2
      
             ak1 = Lafirst(wa)
             am1 = Lasec(n,l,xa,ya,wa)
             ak2 = Lafirst(wa+am1*i)
             am2 = Lasec(n,l,xa+i, ya+i*ak1, wa+i*am1)
             ak3 = Lafirst(wa+am2*i)
             am3 = Lasec(n,l,xa+i,ya+i*ak2, wa+i*am2)
             ak4 = Lafirst(wa+am3*h)
             am4 = Lasec(n,l,xa+h, ya+h*ak3, wa+h*am3)
             
             xa = xa + h
             ya = ya + h*(ak1+2*ak2+2*ak3+ak4)/6
             wa = wa + h*(am1+2*am2+2*am3+am4)/6
         
             cl = cl + 1


!Riemann sum       
                rho = ((0.5*n*xa)**2.0)*((casio(n,l,xa)*ya)**2.0)
                 
             If (rhomax .LT. rho) Then
                r2 = r1
                r1 = xmax
                rhomax  = rho
                
                xmax = 0.5*n*xa
             End If

             If(Low .GT. rho*0.5*n*h) Then
                Low = rho*0.5*n*h        
             End IF     

             T = rho*0.5*n*h + T
             Tlow = Low + Tlow
             Low = rho*0.5*n*h


!Write intial value P (Legendre_{l}^{m}) at x = -0.99999        
             ey0 = 0.0134162401601939D0
!Write intial value P' (Derivative Legendre_{l}^{m}) at x = -0.99999
             ew0 = 670.7952376D0
         
             exa = ex0
             eya = ey0
             ewa = ew0
             ecl = 0
             
             
3            IF(ecl .LE. 499) THEN       
                eh = -exa+COS(pi-ecl*pi/500)
                If(ecl .EQ. 0) THEN
                     eh = 0
                End IF 
              

                ei = eh/2.0D0
      
                ek1 = Lefirst(ewa)
                em1 = Lesec(m,l,exa,eya,ewa)
                ek2 = Lefirst(ewa+em1*ei)
                em2 = Lesec(m,l,exa+ei, eya+ei*ek1, ewa+ei*em1)
                ek3 = Lefirst(ewa+em2*ei)
                em3 = Lesec(m,l,exa+ei,eya+ei*ek2, ewa+ei*em2)
                ek4 = Lefirst(ewa+em3*eh)
                em4 = Lesec(m,l,exa+eh, eya+eh*ek3, ewa+eh*em3)
             
                exa = exa + eh
                eya = eya + eh*(ek1+2*ek2+2*ek3+ek4)/6
                ewa = ewa + eh*(em1+2*em2+2*em3+em4)/6        
             

!             Write(4,*) 0.5*n*xa*exa," , ",0.5*n*xa*SQRT(1-exa**2)," , ", (casio(n,l,xa)*ya*why(m,l)*eya)**2
!             Write(5,*) 0.5*n*xa*exa," , ",-0.5*n*xa*SQRT(1-exa**2)," , ", (casio(n,l,xa)*ya*why(m,l)*eya)**2

              
                ecl = ecl + 1
             
                GO TO 3     
                End If
             GO TO 2
             End If

             riemannError=2*(T-Tlow)*a0**3
             xmaxError = xmax-r2

             Write(*,*) "Max radius and error, r, p(r) and p(r) error:"
             Write(*,*) xmax, xmaxError, xa*0.5*n, T*a0**3, riemannError
!Range of the riemann sum             
!             Write(*,*) T*a0**3-riemannError, T*a0**3+riemannError
         
     End Program
