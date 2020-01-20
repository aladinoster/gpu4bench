/*
	Programme original realise en Fortran 77 en mars 1994
	pour les Travaux Pratiques de Modelisation Numerique
	DEA d'astrophysique et techniques spatiales de Meudon

		par Herve Aussel et Emmanuel Quemener

	Conversion en C par Emmanuel Quemener en aout 1997
        Modification par Emmanuel Quemener en aout 2019

	Remerciements a :

	- Herve Aussel pour sa procedure sur le spectre de corps noir
	- Didier Pelat pour l'aide lors de ce travail
	- Jean-Pierre Luminet pour son article de 1979
	- Le Numerical Recipies pour ses recettes de calcul
	- Luc Blanchet pour sa disponibilite lors de mes interrogations en RG

	Compilation sous gcc ( Compilateur GNU sous Linux ) :

		gcc -O2 -o trou_noir trou_noir.c -lm
*/ 

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#define nbr 256 /* Nombre de colonnes du spectre */

#define PI 3.14159265359

#define TRACKPOINTS 2048

#if TYPE == FP32
#define MYFLOAT float
#else
#define MYFLOAT double
#endif

MYFLOAT atanp(MYFLOAT x,MYFLOAT y)
{
  MYFLOAT angle;

  angle=atan2(y,x);

  if (angle<0)
    {
      angle+=2*PI;
    }

  return angle;
}


MYFLOAT f(MYFLOAT v)
{
  return v;
}

MYFLOAT g(MYFLOAT u,MYFLOAT m,MYFLOAT b)
{
  return (3.*m/b*pow(u,2)-u);
}


void calcul(MYFLOAT *us,MYFLOAT *vs,MYFLOAT up,MYFLOAT vp,
	    MYFLOAT h,MYFLOAT m,MYFLOAT b)
{
  MYFLOAT c[4],d[4];

  c[0]=h*f(vp);
  c[1]=h*f(vp+c[0]/2.);
  c[2]=h*f(vp+c[1]/2.);
  c[3]=h*f(vp+c[2]);
  d[0]=h*g(up,m,b);
  d[1]=h*g(up+d[0]/2.,m,b);
  d[2]=h*g(up+d[1]/2.,m,b);
  d[3]=h*g(up+d[2],m,b);

  *us=up+(c[0]+2.*c[1]+2.*c[2]+c[3])/6.;
  *vs=vp+(d[0]+2.*d[1]+2.*d[2]+d[3])/6.;
}

void rungekutta(MYFLOAT *ps,MYFLOAT *us,MYFLOAT *vs,
		MYFLOAT pp,MYFLOAT up,MYFLOAT vp,
		MYFLOAT h,MYFLOAT m,MYFLOAT b)
{
  calcul(us,vs,up,vp,h,m,b);
  *ps=pp+h;
}


MYFLOAT decalage_spectral(MYFLOAT r,MYFLOAT b,MYFLOAT phi,
			 MYFLOAT tho,MYFLOAT m)
{
  return (sqrt(1-3*m/r)/(1+sqrt(m/pow(r,3))*b*sin(tho)*sin(phi)));
}

MYFLOAT spectre(MYFLOAT rf,MYFLOAT q,MYFLOAT b,MYFLOAT db,
	     MYFLOAT h,MYFLOAT r,MYFLOAT m,MYFLOAT bss)
{
  MYFLOAT flx;

  flx=exp(q*log(r/m))*pow(rf,4)*b*db*h;
  return(flx);
}

MYFLOAT spectre_cn(MYFLOAT rf,MYFLOAT b,MYFLOAT db,
		MYFLOAT h,MYFLOAT r,MYFLOAT m,MYFLOAT bss)
{
  
  MYFLOAT flx;
  MYFLOAT nu_rec,nu_em,qu,temp_em,flux_int;
  int fi,posfreq;

#define planck 6.62e-34
#define k 1.38e-23
#define c2 9.e16
#define temp 3.e7
#define m_point 1.

#define lplanck (log(6.62)-34.*log(10.))
#define lk (log(1.38)-23.*log(10.))
#define lc2 (log(9.)+16.*log(10.))
  
  MYFLOAT v=1.-3./r;

  qu=1./sqrt((1.-3./r)*r)*(sqrt(r)-sqrt(6.)+sqrt(3.)/2.*log((sqrt(r)+sqrt(3.))/(sqrt(r)-sqrt(3.))* 0.17157287525380988 ));

  temp_em=temp*sqrt(m)*exp(0.25*log(m_point)-0.75*log(r)-0.125*log(v)+0.25*log(fabs(qu)));

  flux_int=0.;
  flx=0.;

  for (fi=0;fi<nbr;fi++)
    {
      nu_em=bss*(MYFLOAT)fi/(MYFLOAT)nbr;
      nu_rec=nu_em*rf;
      posfreq=(int)(nu_rec*(MYFLOAT)nbr/bss);
      if ((posfreq>0)&&(posfreq<nbr))
  	{
	  flux_int=2.*planck/c2*pow(nu_em,3)/(exp(planck*nu_em/(k*temp_em))-1.)*exp(3.*log(rf))*b*db*h;
  	  flx+=flux_int;
  	}
    }

  return((MYFLOAT)flx);
}

void impact(MYFLOAT d,MYFLOAT phi,int dim,MYFLOAT r,MYFLOAT b,MYFLOAT tho,MYFLOAT m,
	    MYFLOAT **zp,MYFLOAT **fp,
	    MYFLOAT q,MYFLOAT db,
	    MYFLOAT h,MYFLOAT bss,int raie)
{
  MYFLOAT xe,ye;
  int xi,yi;
  MYFLOAT flx,rf;
  xe=d*sin(phi);
  ye=-d*cos(phi);

  xi=(int)xe+dim/2;
  yi=(int)ye+dim/2;

  rf=decalage_spectral(r,b,phi,tho,m);

  if (raie==0)
    {
      bss=1.e19;
      flx=spectre_cn(rf,b,db,h,r,m,bss);
    }
  else
    {
      bss=2.;
      flx=spectre(rf,q,b,db,h,r,m,bss);
    }
  
  if (zp[xi][yi]==0.)
    {
      zp[xi][yi]=1./rf;
    }
  
  if (fp[xi][yi]==0.)
    {
      fp[xi][yi]=flx;
    }

}

void sauvegarde_pgm(char nom[24],unsigned int **image,int dim)
{
  FILE            *sortie;
  unsigned long   i,j;
  
  sortie=fopen(nom,"w");
  
  fprintf(sortie,"P5\n");
  fprintf(sortie,"%i %i\n",dim,dim);
  fprintf(sortie,"255\n");

  for (j=0;j<dim;j++) for (i=0;i<dim;i++)
    {
      fputc(image[i][j],sortie);
    }

  fclose(sortie);
}

int main(int argc,char *argv[])
{

  MYFLOAT m,rs,ri,re,tho;
  int q;

  MYFLOAT bss,stp;
  int nmx,dim;
  MYFLOAT d,bmx,db,b,h;
  MYFLOAT up,vp,pp;
  MYFLOAT us,vs,ps;
  MYFLOAT rp[TRACKPOINTS];
  MYFLOAT **zp,**fp;
  unsigned int **izp,**ifp;
  MYFLOAT zmx,fmx,zen,fen;
  MYFLOAT flux_tot,impc_tot;
  MYFLOAT phi,thi,thx,phd,php,nr,r;
  int ni,ii,i,imx,j,n,tst,dist,raie,pc,fcl,zcl;
  MYFLOAT nh;

  if (argc==2)
    {
      if (strcmp(argv[1],"-aide")==0)
	{
	  printf("\nSimulation d'un disque d'accretion autour d'un trou noir\n");
	  printf("\nParametres a definir :\n\n");
	  printf("  1) Dimension de l'Image\n");
	  printf("  2) Masse relative du trou noir\n");
	  printf("  3) Dimension du disque exterieur\n");
	  printf("  4) Inclinaison par rapport au disque (en degres)\n");
	  printf("  5) Rayonnement de disque MONOCHROMATIQUE ou CORPS_NOIR\n");
	  printf("  6) Impression des images NEGATIVE ou POSITIVE\n"); 
	  printf("  7) Nom de l'image des Flux\n");
	  printf("  8) Nom de l'image des decalages spectraux\n");
	  printf("\nSi aucun parametre defini, parametres par defaut :\n\n");
	  printf("  1) Dimension de l'image : 1024 pixels de cote\n");
	  printf("  2) Masse relative du trou noir : 1\n");
	  printf("  3) Dimension du disque exterieur : 12 \n");
	  printf("  4) Inclinaison par rapport au disque (en degres) : 10\n");
	  printf("  5) Rayonnement de disque CORPS_NOIR\n");
	  printf("  6) Impression des images NEGATIVE ou POSITIVE\n"); 
       	  printf("  7) Nom de l'image des flux : flux.pgm\n");
	  printf("  8) Nom de l'image des z : z.pgm\n");
	}
    }
  
  if (argc==9)
    {
      printf("# Utilisation les valeurs definies par l'utilisateur\n");
      
      dim=atoi(argv[1]);
      m=atof(argv[2]);
      re=atof(argv[3]);
      tho=PI/180.*(90-atof(argv[4]));
      
      rs=2.*m;
      ri=3.*rs;

      if (strcmp(argv[5],"CORPS_NOIR")==0)
	{
	  raie=0;
	}
      else
	{
	  raie=1;
	}

    }
  else
    {
      printf("# Utilisation les valeurs par defaut\n");
      
      dim=1024;
      m=1.;
      rs=2.*m;
      ri=3.*rs;
      re=12.;
      tho=PI/180.*80;
      // Corps noir
      raie=0;
    }

  if (raie==1)
    {
      bss=2.;
      q=-2;
    }
  else
    {
      bss=1.e19;
      q=-0.75;
    }

      printf("# Dimension de l'image : %i\n",dim);
      printf("# Masse : %f\n",m);
      printf("# Rayon singularite : %f\n",rs);
      printf("# Rayon interne : %f\n",ri);
      printf("# Rayon externe : %f\n",re);
      printf("# Inclinaison a la normale en radian : %f\n",tho);
  
  zp=(MYFLOAT**)calloc(dim,sizeof(MYFLOAT*));
  zp[0]=(MYFLOAT*)calloc(dim*dim,sizeof(MYFLOAT));
  
  fp=(MYFLOAT**)calloc(dim,sizeof(MYFLOAT*));
  fp[0]=(MYFLOAT*)calloc(dim*dim,sizeof(MYFLOAT));

  izp=(unsigned int**)calloc(dim,sizeof(unsigned int*));
  izp[0]=(unsigned int*)calloc(dim*dim,sizeof(unsigned int));
  
  ifp=(unsigned int**)calloc(dim,sizeof(unsigned int*));
  ifp[0]=(unsigned int*)calloc(dim*dim,sizeof(unsigned int));

  for (i=1;i<dim;i++)
    {
      zp[i]=zp[i-1]+dim;
      fp[i]=fp[i-1]+dim;
      izp[i]=izp[i-1]+dim;
      ifp[i]=ifp[i-1]+dim;
    }

  nmx=dim;
  stp=dim/(2.*nmx);
  bmx=1.25*re;
  b=0.;
  pc=0;

  struct timeval tv1,tv2;
  struct timezone tz;
      
  // Set start timer
  gettimeofday(&tv1, &tz);

  for (n=1;n<=nmx;n++)
    {     
      h=4.*PI/(MYFLOAT)TRACKPOINTS;
      d=stp*n;

      db=bmx/(MYFLOAT)nmx;
      b=db*(MYFLOAT)n;
      up=0.;
      vp=1.;
      
      pp=0.;
      nh=1;

      rungekutta(&ps,&us,&vs,pp,up,vp,h,m,b);
    
      rp[(int)nh]=fabs(b/us);
      
      do
	{
	  nh++;
	  pp=ps;
	  up=us;
	  vp=vs;
	  rungekutta(&ps,&us,&vs,pp,up,vp,h,m,b);
	  
	  rp[(int)nh]=b/us;
	  
	} while ((rp[(int)nh]>=rs)&&(rp[(int)nh]<=rp[1]));
      
      for (i=nh+1;i<TRACKPOINTS;i++)
	{
	  rp[i]=0.; 
	}
      
      imx=(int)(8*d);
      
      for (i=0;i<=imx;i++)
	{
	  phi=2.*PI/(MYFLOAT)imx*(MYFLOAT)i;
	  phd=atanp(cos(phi)*sin(tho),cos(tho));
	  phd=fmod(phd,PI);
	  ii=0;
	  tst=0;
	  
	  do
	    {
	      php=phd+(MYFLOAT)ii*PI;
	      nr=php/h;
	      ni=(int)nr;

	      if ((MYFLOAT)ni<nh)
		{
		  r=(rp[ni+1]-rp[ni])*(nr-ni*1.)+rp[ni];
		}
	      else
		{
		  r=rp[ni];
		}
	   
	      if ((r<=re)&&(r>=ri))
		{
		  tst=1;
		  impact(d,phi,dim,r,b,tho,m,zp,fp,q,db,h,bss,raie);
		}
	      
	      ii++;
	    } while ((ii<=2)&&(tst==0));
	}
    }

  // Set stop timer
  gettimeofday(&tv2, &tz);
  
  double elapsed=(double)((tv2.tv_sec-tv1.tv_sec) * 1000000L +
			  (tv2.tv_usec-tv1.tv_usec))/1000000;  
  
  fmx=fp[0][0];
  zmx=zp[0][0];
  
  int zimx=0,zjmx=0,fimx=0,fjmx=0;

  for (i=0;i<dim;i++) for (j=0;j<dim;j++)
    {
      if (fmx<fp[i][j])
	{
	  fimx=i;
	  fjmx=j;
	  fmx=fp[i][j];
	}
      
      if (zmx<zp[i][j])
	{
	  zimx=i;
	  zjmx=j;
	  zmx=zp[i][j];
	}
    }

  printf("\nElapsed Time : %lf",(double)elapsed);
  printf("\nZ max @(%i,%i) : %f",zimx,zjmx,zmx);
  printf("\nFlux max @(%i,%i) : %f\n\n",fimx,fjmx,fmx);

  for (i=0;i<dim;i++) for (j=0;j<dim;j++)
    {
      zcl=(int)(255/zmx*zp[i][dim-1-j]);
      fcl=(int)(255/fmx*fp[i][dim-1-j]);

      if (strcmp(argv[6],"NEGATIVE")==0)
	{
	  if (zcl>255)
	    {
	      izp[i][j]=0;
	    }
	  else
	    {
	      izp[i][j]=255-zcl;
	    }
	  
	  if (fcl>255)
	    {
	      ifp[i][j]=0;
	    }
	  else
	    {
	      ifp[i][j]=255-fcl;
	    } 
	  
	}
      else
	{
	  if (zcl>255)
	    {
	      izp[i][j]=255;
	    }
	  else
	    {
	      izp[i][j]=zcl;
	    }
	  
	  if (fcl>255)
	    {
	      ifp[i][j]=255;
	    }
	  else
	    {
	      ifp[i][j]=fcl;
	    } 
	  
	}
	
    }

  if (argc==9)
   {
     sauvegarde_pgm(argv[7],ifp,dim);
     sauvegarde_pgm(argv[8],izp,dim);
   }
  else
    {
      sauvegarde_pgm("z.pgm",izp,dim);
      sauvegarde_pgm("flux.pgm",ifp,dim);
    }

  free(zp[0]);
  free(zp);
  free(fp[0]);
  free(fp);

  free(izp[0]);
  free(izp);
  free(ifp[0]);
  free(ifp);

}


