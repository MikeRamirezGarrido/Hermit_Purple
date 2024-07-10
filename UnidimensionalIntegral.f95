!Ramirez Garrido
!Hydrogen atom's electron probability densisty and radiuses
!Laguerre
  !Runge Kutta first function y' = w
       Function first(w)
              double precision, intent(IN) :: w
              double precision :: first
              first = w
       End Function
              
!Runge Kutta second function w' = 2xw-2ny   
       Function sec(n,x,y,w)
              double precision, intent(IN) :: n, x, y, w     
              double precision :: sec
              sec = 2*x*w-2*n*y
       End Function
       
  	 Program HydrogenAtom
     		 double precision  :: n, m
     		 double precision :: sec, first
     		 double precision :: ak1, ak2, ak3, ak4, x0, y0, w0
     		 double precision :: h, i, cl, ya, xa, wa
     		 double precision :: am1, am2, am3, am4, ET
     		 double precision ::  ek1, ek2, ek3, ek4, ex0, ey0, ew0
     		 double precision :: em1, em2, em3, em4, eh, ei
     		 double precision :: eya, exa, ewa, ecl, a0
     		 double precision :: pi, masse, omega, e, xs, psi, axs, apsi, T
!     		 Open(4,file='data.txt')
!     		 Open(5,file='data2.txt')

     		 Write(*,*) "Ramirez Garrido"
     		 Write(*,*) "Hydrogen atom's electron probability densisty and radiuses"
     		 pi = 3.14159265D0
               masse = 1.0D0
               omega = 0.00031D0
               hbar = 0.000115767D0
               e = 2.718281828D0
               T = 0.0D0
     		 x0 = -6.0D0
     		 ex0 = -6.0D0
!Write here
!Write quantic numbers
     		 n = 2.0D0
     		
     		 m = 8.0D0
    		 
     		 
    		 
!Write intial value L at x = 0.001		 
     		 y0 = 142.0D0
!Write intial value L' at x = 0.001
     		 w0 = -48.0D0
   		 
     		 xa = x0
     		 ya = y0
     		 wa = w0
     		 cl = 0
     		 
 !Write intial value P at x = -0.99999		 
     		 ey0 = 279702672.0D0
!Write intial value P' at x = -0.99999
     		 ew0 = -417544704.0D0
   		 
     		 exa = ex0
     		 eya = ey0
     		 ewa = ew0
 
2    		 IF(cl .LE. 1200) THEN
     		 h = 0.01D0
    		 
            		 

     		 i = h/2
      
     		 ak1 = first(wa)
     		 am1 = sec(n,xa,ya,wa)
     		 ak2 = first(wa+am1*i)
     		 am2 = sec(n,xa+i, ya+i*ak1, wa+i*am1)
     		 ak3 = first(wa+am2*i)
     		 am3 = sec(n,xa+i,ya+i*ak2, wa+i*am2)
     		 ak4 = first(wa+am3*h)
     		 am4 = sec(n,xa+h, ya+h*ak3, wa+h*am3)
     		 
     		 xa = xa + h
     		 ya = ya + h*(ak1+2*ak2+2*ak3+ak4)/6
     		 wa = wa + h*(am1+2*am2+2*am3+am4)/6
   		 
     		 cl = cl + 1

                     
              xs = SQRT(hbar/(masse*omega))*xa
           	psi = ( 1/((2**n)*gamma(n+1))**0.5)*((masse*omega/(pi*hbar))**0.25)*(e**(-masse*omega*xs**2/(2*hbar)))	 
     	


     		 
     		 
    		 

     		 ek1 = first(ewa)
     		 em1 = sec(m,exa,eya,ewa)
     		 ek2 = first(ewa+em1*i)
     		 em2 = sec(m,exa+i, eya+i*ek1, ewa+i*em1)
     		 ek3 = first(ewa+em2*i)
     		 em3 = sec(m,exa+i,eya+i*ek2, ewa+i*em2)
     		 ek4 = first(ewa+em3*h)
     		 em4 = sec(m,exa+h, eya+h*ek3, ewa+h*em3)
     	
     		 exa = exa + h
     		 eya = eya + h*(ek1+2*ek2+2*ek3+ek4)/6
     		 ewa = ewa + h*(em1+2*em2+2*em3+em4)/6
     	
    		axs = SQRT(hbar/(masse*omega))*xa
           	apsi = ( 1/((2**n)*gamma(n+1))**0.5)*((masse*omega/(pi*hbar))**0.25)*(e**(-masse*omega*axs**2/(2*hbar)))
              T = T + ya*eya*psi*apsi*xs*h/xa
     	
    		 
     		 GO TO 2
    		 
     		 End If
   		 
     		 
    	
     		 
    		 
     	
   		 
!Last elements of the Laguerre's and Legendre's Runge-Kutta results. To compare  		 
     		 Write(*,*) xa, ya
     		 Write(*,*) exa, eya
     		 Write(*,*) T
     		 
  	 End Program
