!Uni-dimesional quantum harmonic oscillator
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
              double precision :: hbar, omega, masse, e, pi, xs, psi
              double precision :: T, Error
              Write(*,*) "Hermit Purple by Ramirez Garrido"
              Write(*,*) "Uni-dimesional quantum harmonic oscillator"
              Write(*,*) "Write initial values in source code"
              pi = 3.141592654D0
              e = 2.718281828D0
              Error = 0.0D0
!hbar using electron mass as mass unit               
              hbar = 0.000115767D0
!Write here              
              x0 = -7.0D0
              masse = 1.0D0
              omega = 0.00031D0
              n = 8.0D0
!Write intial value H_n at x = x0     
              y0 = 1085747600.0D0
!Write intial value H_n' at x = x0,   H_{n}'(x) = 2*n*H_{n-1} (x) if that helps. 
              w0 = 2*8*(-83965616.0D0)
              T = 0.0D0
!              Open(2, file="data.txt")
              xa = x0
              ya = y0
              wa = w0
              l = 1.0
1             IF(l .LT. 14000) THEN           
              h = 0.001D0
              
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
	      
!              Write(2,*) xs, ",", ya*psi, ",", (ya*psi)**2
              T = T + ((ya*psi)**2)*xs*h/xa
              Error = Error + ((ya*psi)**2)*xs*h*(0.0005)/xa
              l = l + 1.0
              GO TO 1
              End If
              Write(*,*) T, Error
       End Program
