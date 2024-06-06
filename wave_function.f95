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
       
       Program RungeKutta
              double precision :: n, x, y, w, sec, first
              double precision :: k1, k2, k3, k4, x0, y0, w0
              double precision :: m, h, i, l, ya, xa, wa
              double precision :: m1, m2, m3, m4
              double precision :: hbar, omega, masse, e, pi, xs, psi, T
              Write(*,*) "Hermit Purple by Ramirez Garrido"
              Write(*,*) "Write initial values in source code"
              x0 = -7
              pi = 3.141592654
              e = 2.718281828
              

              masse = 1
              omega = 0.000115767
              hbar = 0.000115767
              
              
!Write here
!Write degree
              n = 5
             
!Write intial value H_n at x = 0     
              y0 = -483784
!Write intial value H_n' at x = 0,   H_{n}'(x) = 2*n*H_{n-1} (x) if that helps. 
              w0 = 2*5*(36076)
              T = 0
              Open(3, file="data.txt")
              xa = x0
              ya = y0
              wa = w0
              l = 1.0
2             IF(l .LT. 14000.0) THEN           
              h = 0.001
              
              k1 = 0
              k2 = 0
              k3 = 0
              k4 = 0
              m1 = 0
              m2 = 0
              m3 = 0
              m4 = 0
              i = h/2
          
       
	        k1 = first(wa)
		m1 = sec(n,xa,ya,wa)
       	      	k2 = first(wa+m1*i)
       	        m2 = sec(n,xa+i, ya+i*k1, wa+i*m1)
       	        k3 = first(wa+m2*i)
       	        m3 = sec(n,xa+i,ya+i*k2, wa+i*m2)
       	        k4 = first(wa+m3*h)
       	        m4 = sec(n,xa+h, ya+h*k3, wa+h*m3)
       	       
                     xa = xa + h
                     ya = ya + h*(k1+2.0*k2+2.0*k3+k4)/6.0
                     wa = wa + h*(m1+2.0*m2+2.0*m3+m4)/6.0
              
              	
              xs = SQRT(hbar/(masse*omega))*xa
              psi = ( 1/((2**n)*gamma(n+1))**0.5)*((masse*omega/(pi*hbar))**0.25)*(e**(-masse*omega*xs**2/(2*hbar)))
	      
              Write(3,*) xs,",", ya,",", ya*psi, ",", (ya*psi)**2
              T = T + ((ya*psi)**2)*xs*h/xa
              l = l + 1.0
              GO TO 2
              End If
              Write(*,*) T
       End Program
