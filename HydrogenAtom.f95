       Function factorial(x)
              double precision , intent(IN) :: x
              double precision  :: factorial, s,i
              s = 1
              i = 1
1              If(i.LE.x) Then
                     s = s*i 
                     i = i + 1
                     Go To 1
              End If
              factorial  = s
       End Function

!Laguerre
       Function Lafirst(MLa)
              double precision, intent(IN) :: MLa
              double precision :: Lafirst
              Lafirst = MLa
       End Function
              
!Runge Kutta second function w' = 2xw-2ny   
       Function Lasec(n,l,x,La,MLa)
              double precision , intent(IN) :: n,l 
              double precision, intent(IN) :: x, La, MLa
              double precision :: Lasec
              Lasec = -((2*l+2-x)*MLa+(n+l)*La)/(x)
       End Function
       
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
       
       Function why(m,l)
              double precision , intent(IN) :: m, l
              double precision :: why, pi
              double precision  :: factorial
              pi = 3.141592653589793238
              why = ((-1)**m)*SQRT(((2*l+1)/(4*pi))*(factorial(l-m)/factorial(l+m)))
       End Function
       
       Function casio(n,l,x)
              double precision, intent(IN) :: n, l, x
              double precision :: a0, e, factorial, casio
              e = 2.14159265358979
              a0 = 0.000529462966
              casio = -(2/n**2)*(SQRT((factorial(n-l-1)/(factorial(n+l)*a0)**3)))*((x)**l)*e**(-0.5*x)
       End Function 

       Program HydrogenAtom
              double precision  :: n, m, l, factorial
              double precision :: x, y, w, Lasec, Lafirst
              double precision :: ak1, ak2, ak3, ak4, x0, y0, w0
              double precision :: h, i, cl, ya, xa, wa
              double precision :: am1, am2, am3, am4
              double precision ::  ek1, ek2, ek3, ek4, ex0, ey0, ew0
              double precision :: em1, em2, em3, em4, eh, ei
              double precision :: eya, exa, ewa, ecl
              double precision :: Lesec, Lefirst
              double precision :: why, casio, rho
              Open(4,file='data.txt')
              Open(5,file='data2.txt')
              Write(*,*) "Hermit Purple by Ramirez Garrido"
              Write(*,*) "Runge_Kutta aproximation for Hermit polynomial"

              x0 = 0.01
              ex0 = -0.999999
!Write here
!Write degree
              n = 7
              l = 1
              m = 1

!Write intial value H_n at x = 0           
              y0 = 161.723023137
!Write intial value H_n' at x = 0,   H_{n}'(x) = 2*n*H_{n-1} (x) if that helps. 
              w0 = -325.4031
!thank you              
              
              xa = x0
              ya = y0
              wa = w0
              cl = 0
  
2             IF(cl .LT. 300) THEN 
              h = 0.08
          
              
              
              ak1 = 0
              ak2 = 0
              ak3 = 0
              ak4 = 0
              am1 = 0
              am2 = 0
              am3 = 0
              am4 = 0
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
             
              
              
              !Write intial value H_n at x = 0           
              ey0 = -0.00141421320881966026
!Write intial value H_n' at x = 0,   H_{n}'(x) = 2*n*H_{n-1} (x) if that helps. 
              ew0 = -707.1062509
!thank you              
              
              exa = ex0
              eya = ey0
              ewa = ew0
              ecl = 0
              
3             IF(ecl .LT. 19999) THEN           

              eh = 0.0001
              
              ek1 = 0
              ek2 = 0
              ek3 = 0
              ek4 = 0
              em1 = 0
              em2 = 0
              em3 = 0
              em4 = 0
              ei = eh/2
       
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
              

              If(xa .LT. 4 .AND. xa .GE. 2 ) Then
              Write(4,*) 0.5*n*xa*exa," , ",0.5*n*xa*SQRT(1-exa**2)," , ", (casio(n,l,xa)*ya*why(m,l)*eya)**2
              Write(5,*) 0.5*n*xa*exa," , ",-0.5*n*xa*SQRT(1-exa**2)," , ", (casio(n,l,xa)*ya*why(m,l)*eya)**2
              End If
              ecl = ecl + 1
              GO TO 3
              End If
              GO TO 2
              End If
       End Program
