      PROGRAM MODELISATION_NUMERIQUE

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                            c
c MODELISATION NUMERIQUE D'UN DISQUE D'ACCRETION AUTOUR D'UN TROU NOIR :     c
c                                                                            c
c APPARENCE & SPECTRE                                                        c
c                                                                            c
c Programme realise par : Herve AUSSEL                                       c
c                         Emmanuel QUEMENER                                  c
c                                                                            c
c Dans le cadre de l'unite de modelisation numerique du DEA                  c
c d'Astrophysique & Techniques Spatiales (Paris 7 & 11)                      c
c                                                                            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      IMPLICIT NONE
      REAL*8 pi
      PARAMETER (pi=3.141592654)

c
c PARAMETRES PHYSIQUES DU SYSTEME DISQUE, TROU-NOIR, OBSERVATEUR
c

      REAL*8 m                  ! masse reduite du trou noir
      REAL*8 rs,ri,re           ! rayons de Schwarzschild,interieur et exterieur
      REAL*8 tho                ! angle normale-disque & direction-observateur
      REAL*8 ro                 ! distance entre l'observateur et le trou noir
      REAL*8 q               ! parametre de spectre
      REAL*4 bss                ! borne superieure du spectre

c
c PARAMETRES DE LA PLAQUE PHOTOGRAPHIQUE
c

      INTEGER*4 dim             ! dimension de la plaque photographique (en pixels)
      PARAMETER (dim=1024)
      REAL*4 stp                ! pas d'iteration en pixels sur la plaque
      INTEGER*4 nmx             ! nombre total d'iterations
      REAL*8 d                  ! distance exploree du centre de la plaque
      REAL*8 bmx                ! parametre d'impact maximal

c
c PARAMETRES INITIAUX D' INTEGRATION
c

      REAL*8 db                 ! avancement du parametre d'impact
      REAL*8 b                  ! parametre d'impact du photon
      REAL*8 h                  ! avancement angulaire initial

c
c PARAMETRES D'INTEGRATION DU SYSTEME DIFFERENTIEL
c

      REAL*8 up,vp,pp           ! parametres du SED , precedents
      REAL*8 us,vs,ps           ! parametres du SED , suivants

c
c RESULTATS DE L'INTEGRATION DE L'EQUATION DIFFERENTIELLE
c

      REAL*8 rp(2048)    ! vecteur polaire
      INTEGER*4 nh              ! nombre d'avancements angulaires

c
c RESULTATS DE LA SIMULATION
c

      REAL*4 zp(dim,dim) 	! image en decalage spectral
      REAL*4 fp(dim,dim) 	! image en flux

c
c AUTRES VARIABLES...
c

      INTEGER*4 pc              ! pourcentage du balayage realise
      REAL*8 phi                ! angle de la plaque photographique
      REAL*8 thi                ! angle trou noir-observateur & direction emission
      REAL*8 thx                ! angle thi maximal
      REAL*8 phd                ! angle donnant l'intersection des plans
      REAL*8 php                ! angle d'intersection intermediaire
      REAL*8 nr                 ! indice d'intersection reel du vecteur
      INTEGER*4 ni              ! indice d'intersection entier du vecteur
      INTEGER*4 ii              ! indice de 1' image (primaire, secondaire, ...)
      REAL*8 r                  ! rayon moyen d'impact sur le disque
      INTEGER*4 nbr             ! nombre de bandes de rapport de frequence
      INTEGER*4 i,imx,j,n
      LOGICAL tst
      REAL*8 atanp              ! fonction arctangente : resultat entre 0 et 2*pi
      LOGICAL raie            ! variable interactive de type de spectre
      LOGICAL dist            ! variable interactive de distance

      REAL*4 t_start,t_stop
c
c TYPE DE SPECTRE
c

c      raie=.TRUE.               ! TRUE -> raie infiniment fine
      raie=.FALSE.               ! FALSE -> corps noir

c
c TYPE DE DISTANCE
c

c      dist=.TRUE.               ! TRUE -> distance infinie
      dist=.FALSE.               ! FALSE -> distance finie

c
c INITIALISATION DES TABLEAUX
c

      DO i=1,dim
         DO j=1,dim
            zp(i,j)=0.
            fp(i,j)=0.
         END DO
      END DO

c
c DEFINITION DES PARAMETRES PHYSIQUES
c

      m=1.
      rs=2.*m
      ri=3.*rs
      re=12.
      ro=100.
      tho=pi/180.*80.
      IF (raie .EQV. .TRUE.) THEN
         q=-2
      ELSE
         q=-3/4
      ENDIF

c
c DEFINITION DES PARAMETRES NUMERIQUES
c

      nmx=dim
      stp=dim/(2.*nmx)
      bmx=1.25*re
      b=0.
      thx=asin(bmx/ro)
      pc=0
      nbr=256
      IF (raie .EQV. .TRUE.) THEN
         bss=2.
      ELSE
         bss=1.e19
      ENDIF

c
c BOUCLE PRINCIPALE DE LA SIMULATION
c

      CALL TIMER(t_start)
      
      DO n=1,nmx

c
c Affectation des conditions initiales de la resolution Runge-Kutta
c ( Premier etage d'imbrication -> n )
c
         h=4*pi/2048.
         d=stp*n
         IF (dist .EQV. .TRUE.) THEN
            db=bmx/nmx
            b=db*n
            up=0.
            vp=1.
         ELSE
            thi=thx/nmx*n
            db=ro*DSIN(thi)-b
            b=ro*DSIN(thi)
            up=DSIN(thi)
            vp=DCOS(thi)
         ENDIF
         pp=0.
         nh=1
         CALL rungekutta(ps,us,vs,pp,up,vp,h,m,b)
         rp(nh)=ABS(b/us)

c
c Resolution de l'equation differentielle pour un parametre b donne
c Les resultats sont stockes dans le vecteur polaire rp(n)
c ( Deuxieme etage d'imbrication -> n & nh )
c

         DO WHILE ((rp(nh).GE.rs).AND.(rp(nh).LE.rp(1)))
            nh=nh+1
            pp=ps
            up=us
            vp=vs
            CALL rungekutta(ps,us,vs,pp,up,vp,h,m,b)
            rp(nh)=b/us
         END DO

c
c Mise a zero de tous les termes non redefinis du vecteur polaire
c ( Deuxieme etage d'imbrication -> n & 1 )
c

         DO i=nh+1,2048
            rp(i)=0. 
         END DO

c
c Determination des points de la trajectoire coupant le plan du disque
c par echantillonnage de 1' angle phi sur un cercle, a pas radial constant
c ( Deuxieme etage d'imbrication -> n & 1 )
c

         imx=INT(8*d)
         DO i=0,imx
            phi=2*pi/imx*i

c Calcul de l'angle "brut" d'intersection dans la base polaire

            phd=atanp(DCOS(phi)*DSIN(tho),DCOS(tho))
            phd=phd-INT(phd/pi)*pi
            ii=0
            tst=.FALSE.

c
c Boucle tant que - que l'image tertiaire n'a pas ete teste
c - le rayon n'est pas dans le domaine du disque
c ( Troisieme etage d'imbrication -> n , 1 , (ii, tst) , r )
c

            DO WHILE((ii.LE.2).AND.(tst.EQV..FALSE.))

c Calcul de l'angle d'intersection pour chacune des images

               php=phd+ii*pi
               nr=php/h
               ni=INT(nr)

c Interpolation lineaire de la distance photon-trou noir

               IF (ni.LT.nh) THEN
                  r=(rp(ni+1)-rp(ni))*(nr-ni*1.)+rp(ni)
               ELSE
                  r=rp(ni)
               ENDIF


c
c Test d'impact sur le disque
c ( Quatrieme etage d'imbrication -> n , i , (ii,tst) , r )
c

               IF ((r.LE.re).AND.(r.GE.ri)) then
                  tst=.TRUE.

c S'il y a impact calcul - du rapport de frequence
c - du spectre de photons et de flux

                  CALL impact(raie,d,phi,dim,r,b,tho,m,zp,fp,
     & q,db,h,nbr,bss)

               END IF

               ii=ii+1
            END DO
         END DO
      END DO


      CALL TIMER(t_stop)

      WRITE (*,*) 'Elapsed time = ', t_stop - t_start
      
c Appel de la routine d'affichage

      CALL affichage(zp,fp,dim,nbr,bss,bmx,m,ri,re,tho,pi)

      WRITE(*,*)
      WRITE(*,*) 'C''est TERMINE ! '
      WRITE(*,*) 

      STOP
      END PROGRAM

c
c Cette fonction permet de trouver un arctangente entre 0 et 2*pi
c

      REAL*8 function atanp(x,y)
      REAL*8 x,y,pi,eps
      pi=3.141592654
      eps=1.e-20
      IF (ABS(x).LE.eps) then
         IF (ABS(y).EQ.Y) then
            atanp=pi/2.
         ELSE
            atanp=3.*pi/2.
         END IF
      ELSE
         IF (ABS(x).EQ.x) then
            IF (ABS(y).EQ.Y) then
               atanp=DATAN(y/x)
            ELSE
               atanp=DATAN(y/x)+2.*pi
            END IF
         ELSE
            atanp=DATAN(y/x)+pi
         END IF
      END IF
      RETURN
      END

c
c
c
      SUBROUTINE timer(ytime)
      REAL*4 ytime
      INTEGER*4 clk_max,clk_rate,clk_read

      CALL SYSTEM_CLOCK(clk_read,clk_rate,clk_max)

      ytime=FLOAT(clk_read)/FLOAT(clk_rate)

      RETURN
      END
c
c Premiere fonction du systeme -> ( d(u)/d(phi)=f(phi,u,v) )
c

      REAL*8 function f(v)
      REAL*8 v
      f=v
      RETURN
      END

c
c Deuxieme fonction du systeme -> ( d(v)/d(phi)=g(phi,u,v) )
c

      REAL*8 function g(u,m,b)
      REAL*8 u, m, b
      g=3.*m/b*u**2-u
      RETURN
      END

c
c Routine de d'intergration par la methode de Runge-Kutta
c

      SUBROUTINE rungekutta(ps,us,vs,pp,up,vp,h,m,b)
      REAL*8 ps,us,vs,pp,up,vp,h,m,b
      CALL calcul(us,vs,up,vp,h,m,b)
      ps=pp+h
      RETURN
      END

      SUBROUTINE calcul(us,vs,up,vp,h,m,b)
      REAL*8 us,vs,up,vp,h,m,b
      REAL*8 c(4),d(4),f,g
      c(1)=h*f(vp)
      c(2)=h* f(vp+c(1)/2.)
      c(3)=h*f(vp+c(2)/2.)
      c(4)=h*f(vp+c(3))
      d(1)=h*g(up,m,b)
      d(2)=h*g(up+d(1)/2.,m,b)
      d(3)=h*g(up+d(2)/2.,m,b)
      d(4)=h*g(up+d(3),m,b)
      us=up+(c(1)+2.*c(2)+2.*c(3)+c(4))/6.
      vs=vp+(d(1)+2.*d(2)+2.*d(3)+d(4))/6.
      RETURN
      END

c
c En cas d'impact, cette procedure :
c - etablit la position du pixel consisdere sur les images
c - appelle la routine de calcul de decalage spectral
c - appelle la routine de calcul du flux et du spectre
c - affecte a chacun des pixels sa valeur calculee
c

      SUBROUTINE impact(raie,d,phi,dim,r,b,tho,m,zp,fp, 
     & q,db,h,nbr,bss)
      REAL*8 d,phi,r,b,tho,m,db,h,q
      INTEGER*4 dim
      REAL*4 zp(dim,dim),fp(dim,dim),bss
      REAL*4 xe,ye
      INTEGER*4 xi,yi
      REAL*8 flx,rf
      LOGICAL raie
      xe=d*DSIN(phi)
      ye=-d*DCOS(phi)

c Calcul des coordonnees (x,y) du pixel sur chacune des images

      xi=INT(xe)+INT(dim/2.)
      yi=INT(ye)+INT(dim/2.)

c Appel de la routine de decalage spectral

      CALL decalage_spectral(rf,r,b,phi,tho,m)

c Appel de la routine de calcul de flux et de modification de spectre
c pour le cas d'une raie infiniment fine ou un corps noir

      IF (raie .EQV. .TRUE.) THEN
         CALL spectre(rf,q,b,db,h,r,m,bss,flx)
      ELSE
         CALL spectre_bb(rf,b,db,h,r,m,nbr,bss,flx)
      END IF

c Affectation sur chacune des images du decalage spectral et du flux

      IF(zp(xi,yi).EQ.0.) zp(xi,yi)=1./rf
      IF(fp(xi,yi).EQ.0.) fp(xi,yi)=flx

      RETURN
      END

c Calcul du rapport entre la frequence recue et la frequence emise

      SUBROUTINE decalage_spectral(rf,r,b,phi,tho,m)
      REAL*8 rf,r,b,phi,tho,m
      rf=sqrt(1-3*m/r)/(1+sqrt(m/r**3)*b*sin(tho)*sin(phi))
      RETURN
      END

c Calcul du flux puis du spectre pour une raie infiniment fine

      SUBROUTINE spectre(rf,q,b,db,h,r,m,bss,flx)
      REAL*8 rf,b,db,h,r,m,flx,q
      REAL*4 bss
      flx=(r/m)**q*rf**4*b*db*h
      RETURN
      END

c Calcul du flux puis du spectre pour un corps noir

      SUBROUTINE spectre_bb(rf,b,db,h,r,m,nbr,bss,flx)
      REAL*8 nu_rec,nu_em,qu,v
      REAL*8 rf,b,db,h,r,m,flx,temp,temp_em,flux_int
      REAL*8 m_point
      REAL*4 bss
      INTEGER*4 nbr
      INTEGER*4 fi,posfreq
      REAL*8, PARAMETER :: planck=6.62e-34 ! definition des constantes
      REAL*8, PARAMETER :: c=3.e8
      REAL*8, PARAMETER :: k=1.38e-23

c
c Definition des parametres pour le calcul de la temperature d'emission
c puis calcul
c

      temp=3.e7
      m_point=1.

c v -> C du rapport

      v=1.-3./r

c qu -> Q du rapport

      qu=((1-3./r)**-0.5)*(r**-0.5)*(r**0.5-6.**0.5+
     & (3.**0.5/2.)*log((r**0.5+3.**0.5)/(r**0.5-3**0.5)*0.1715728752))
      qu=abs(qu)
      temp_em=temp*m**0.5
      temp_em=temp_em*m_point**0.25
      temp_em=temp_em*r**(-0.75)
      temp_em=temp_em*v**(-0.125)
      temp_em=temp_em*qu**0.25

c 
c initialisation des compteurs de flux
c flux int : flux recu pour le pixel dans la tranque posfreq
c flx : flux "bolometrique", integre de 0. a bss
c

      flux_int=0.
      flx=0.

c
c Balayage en frequence du spectre
c

      DO fi=1,nbr
         nu_em=bss*fi/nbr 	!frequence d' emission
         nu_rec=nu_em*rf 	!frequence de reception
         posfreq=INT(nu_rec*FLOAT(nbr)/bss) !case ou se trouve nu rec

c
c test pour savoir si la frequence de reception est bien dans le
c domaine du spectre calcule
c

         if ((posfreq.GE.1).AND.(posfreq.LE.nbr)) THEN

c Loi du corps noir

            flux_int=2.*nu_em**3/c**2*planck/
     & (exp(planck*nu_em/(k*temp_em))-1.)

c Integration sur le pixel

            flux_int=flux_int*rf**3*b*db*h

c Integration bolometrique

            flx=flx+flux_int
         ENDIF
      END DO

      RETURN
      END

c
c AFFICHAGE DES RESULTATS : 6 fenetres de sortie
c
c - Image en decalage spectral
c - Image en flux
c - Dessin de parametres physiques interessants
c - Dessin des limites du disque a l'infini
c - Dessin du spectre en flux
c - Dessin du spectre en nombre de photons
c

      SUBROUTINE affichage(zp,fp,dim,nbr,bss,bmx,m,ri,re,tho,pi)
      INTEGER*4 dim
      REAL*4 zp(dim,dim),fp(dim,dim),bss
      REAL*8 bmx,m,ri,re,tho,pi
      INTEGER*4 np,i,j
      REAL*8 fm,nm
      REAL*4 fmx,zmx,x(1024),y(1024),tr(6),c(20),bmp
      LOGICAL tg
      LOGICAL ps

      tg=.TRUE.                 ! TRUE : affichage en tons de gris
c                               ! FALSE : affichage en contours

      ps=.TRUE.                ! TRUE : affichage postscript 
c                               ! FALSE : affichage fenetre Xwindow


c
c TRANSFORMATION DE L'IMAGE DES FLUX : FLUX -> FLUX RELATIFS
c

c
c DETERMINATION DU DECALAGE SPECTRAL MAXIMUM ET DU FLUX MAXIMUM
c

      fmx=0.
      zmx=0.
      DO i=1,dim
         DO j=1,dim
            IF (fmx.LT.fp(i,j)) fmx=fp(i,j)
            IF (zmx.LT.zp(i,j)) zmx=zp(i,j)
         END DO
      END DO

      write(*,*) 'Flux max :',fmx
      write(*,*) 'Z max :',zmx
      
      
      DATA tr/0.,1.,0.,0.,0.,1./

c
c AFFICHAGE DE L'IMAGE DES DECALAGES SPECTRAUX
c

      if (ps .EQV. .TRUE.) THEN
         CALL PGBEGIN (0,'z.ps/ps',1,1)
      ELSE
         CALL PGBEGIN (0,'/xwindow',1,1)
      END IF
      CALL PGENV (1.,dim/1.,1.,dim/1.,1,-1)
      CALL PGLABEL('Image des decalages spectraux','',
     &'DISQUE D''ACCRETION AUTOUR D''UN TROU NOIR')
      IF (tg .EQV. .TRUE.) THEN
         DO i=0,25
            DO j=0,25
               zp(i+25,j+25)=1.
            END DO
         END DO
         CALL PGGRAY(zp,dim,dim,2,dim-1,2,dim-1,zmx,0.,tr)
      ELSE
         DO i=1,20
            c(i)=0.2*(i-1)
         END DO
         CALL PGCONS(zp,dim,dim,1,dim,1,dim,c,20,tr)
      END IF
      CALL PGEND

c
c AFFICHAGE DE L'IMAGE DES FLUX
c

      if (ps .EQV. .TRUE.) THEN
         CALL PGBEGIN (0,'flux.ps/ps',1,1)
      ELSE
         CALL PGBEGIN (0,'/xwindow',1,1)
      END IF
      CALL PGENV(1.,dim/1.,1.,dim/1.,1,-1)
      CALL PGLABEL(' Image des flux','',
     &'DISQUE D''ACCRETION AUTOUR D''UN TROU NOIR')
      if (tg) THEN
         DO i=0,25
            DO j=0,25
               fp(i+25,j+25)=1.
            END DO
         END DO
         CALL PGGRAY(fp,dim,dim,2,dim-1,2,dim-1,fmx,0.,tr)
      ELSE
         DO i=1,8
            c(i)=0.0625*2.**(i-1)
         END DO
         CALL PGCONS(fp,dim,dim,1,dim,1,dim,c,8,tr)
      END IF
      CALL PGEND

      END
c compilation : 
c gfortran --std=gnu -o trou_noir_Fortran trou_noir.f -lpgplot
c SORTIE GRAPHIQUE : -> Envoi sur l'ecran : /xwindow
c ( dans PGBEGIN ) -> Envoi en format postscript : nom image.ps/ps





