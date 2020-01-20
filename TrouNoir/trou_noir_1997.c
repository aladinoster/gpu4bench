/*
	Programme original realise en Fortran 77 en mars 1994
	pour les Travaux Pratiques de Modelisation Numerique
	DEA d'astrophysique et techniques spatiales de Meudon

		par Herve Aussel et Emmanuel Quemener

	Conversion en C par Emmanuel Quemener en aout 1997

	Remerciements a :

	- Herve Aussel pour sa procedure sur le spectre de corps noir
	- Didier Pelat pour l'aide lors de ce travail
	- Jean-Pierre Luminet pour son article de 1979
	- Le Numerical Recipies pour ses recettes de calcul
	- Luc Blanchet pour sa disponibilite lors de mes interrogations en RG

	Mes Coordonnees :	Emmanuel Quemener
				Departement Optique
				ENST de Bretagne
				BP 832
				29285 BREST Cedex

				Emmanuel.Quemener@enst-bretagne.fr

	Compilation sous gcc ( Compilateur GNU sous Linux ) :

		gcc -O6 -m486 -o trou_noir trou_noir.c -lm
*/ 

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define nbr 200 /* Nombre de colonnes du spectre */

#define PI 3.14159265359

double atanp(double x,double y)
{
  double angle;

  angle=atan2(y,x);

  if (angle<0)
    {
      angle+=2*PI;
    }

  return angle;
}


double f(double v)
{
  return v;
}

double g(double u,double m,double b)
{
// return (3.*m/b*pow(u,2)-u);
  return (3.*m/b*pow(u,2)-u);
}


void calcul(double *us,double *vs,double up,double vp,
	    double h,double m,double b)
{
  double c[4],d[4];

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

void rungekutta(double *ps,double *us,double *vs,
		double pp,double up,double vp,
		double h,double m,double b)
{
  calcul(us,vs,up,vp,h,m,b);
  *ps=pp+h;
}


double decalage_spectral(double r,double b,double phi,
			 double tho,double m)
{
  return (sqrt(1-3*m/r)/(1+sqrt(m/pow(r,3))*b*sin(tho)*sin(phi)));
}

void spectre(int nt[nbr],double fx[nbr],double rf,int q,double b,double db,
	     double h,double r,double m,double bss,double *flx)
{
  int fi;

  fi=(int)(rf*nbr/bss);
  nt[fi]+=1;
  *flx=pow(r/m,q)*pow(rf,4)*b*db*h;
  fx[fi]=fx[fi]+*flx;
}

void spectre_cn(int nt[nbr],double fx[nbr],double rf,double b,double db,
		double h,double r,double m,double bss,double *flx)
{
  double nu_rec,nu_em,qu,v,temp,temp_em,flux_int,m_point,planck,c,k;
  int fi,posfreq;

  planck=6.62e-34;
  k=1.38e-23;
  temp=3.e7;
  // m_point=1.e14;
  m_point=10.;
  v=1.-3./r;

  qu=1/sqrt(1-3./r)/sqrt(r)*(sqrt(r)-sqrt(6)+sqrt(3)/2*log((sqrt(r)+sqrt(3))/(sqrt(r)-sqrt(3))*(sqrt(6)-sqrt(3))/(sqrt(6)+sqrt(3))));

  temp_em=temp*sqrt(m)*exp(0.25*log(m_point))*exp(-0.75*log(r))*
    exp(-0.125*log(v))*exp(0.25*log(qu));

  flux_int=0;
  *flx=0;

  for (fi=1;fi<nbr;fi++)
    {
      nu_em=bss*fi/nbr;
      nu_rec=nu_em*rf; 
      posfreq=1./bss*nu_rec*nbr;
      if ((posfreq>0)&&(posfreq<nbr))
	{
	  flux_int=2*planck/9e16*pow(nu_em,3)/(exp(planck*nu_em/k/temp_em)-1.)*pow(rf,3)*b*db*h;
	  fx[posfreq]+=flux_int;
	  *flx+=flux_int;
	  nt[posfreq]+=1;
	}
    }

  printf("%f %f %f %f\n",b,db,r,*flx);
}

void impact(double d,double phi,int dim,double r,double b,double tho,double m,
	    double **zp,double **fp,
	    int nt[200],double fx[200],int q,double db,
	    double h,double bss,int raie)
{
  double xe,ye;
  int xi,yi;
  double flx,rf;
  xe=d*sin(phi);
  ye=-d*cos(phi);

  xi=(int)xe+dim/2;
  yi=(int)ye+dim/2;

  rf=decalage_spectral(r,b,phi,tho,m);

  if (raie==0)
    {
      spectre(nt,fx,rf,q,b,db,h,r,m,bss,&flx);
    }
  else
    {
      spectre_cn(nt,fx,rf,b,db,h,r,m,bss,&flx);
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

void sauvegarde_dat(char nom[24],double tableau[3][nbr],int raie)
{
  FILE            *sortie;
  unsigned long   i;
  
  sortie=fopen(nom,"w");

  fprintf(sortie,"# Trou Noir entoure d'un Disque d'Accretion\n");
 
  if (raie==0)
    {
      fprintf(sortie,"# Colonne 1 : Frequence_Recue/Frequence_Emise\n");
    }
  else
    {
      fprintf(sortie,"# Colonne 1 : Fréquence d'Emission en Hertz\n");
    }

  fprintf(sortie,"# Colonne 2 : Intensite Normalisee\n");
  fprintf(sortie,"# Colonne 3 : Nombre d'Impacts Normalise\n");
  
  for (i=1;i<nbr;i++)
    {
      fprintf(sortie,"%f\t%f\t%f\n",tableau[0][i],tableau[1][i],tableau[2][i]);
    }

  fclose(sortie);
}

int main(int argc,char *argv[])
{

  double m,rs,ri,re,tho,ro;
  int q;

  double bss,stp;
  int nmx,dim;
  double d,bmx,db,b,h;
  double up,vp,pp;
  double us,vs,ps;
  double rp[2000];
  double **zp,**fp;
  unsigned int **izp,**ifp;
  double zmx,fmx,zen,fen;
  double flux_tot,impc_tot;
  double fx[nbr];
  int nt[nbr];
  double tableau[3][nbr];
  double phi,thi,thx,phd,php,nr,r;
  int ni,ii,i,imx,j,n,tst,dist,raie,pc,fcl,zcl;
  double nh;

  if (argc==2)
    {
      if (strcmp(argv[1],"-aide")==0)
	{
	  printf("\nSimulation d'un disque d'accretion autour d'un trou noir\n");
	  printf("\nParametres a definir :\n\n");
	  printf("   1) Dimension de l'Image\n");
	  printf("   2) Masse relative du trou noir\n");
	  printf("   3) Dimension du disque exterieur\n");
	  printf("   4) Distance de l'observateur\n");
	  printf("   5) Inclinaison par rapport au disque (en degres)\n");
	  printf("   6) Observation a distance FINIE ou INFINIE\n");
	  printf("   7) Rayonnement de disque MONOCHROMATIQUE ou CORPS_NOIR\n");
	  printf("   8) Normalisation des flux INTERNE ou EXTERNE\n");
	  printf("   9) Normalisation de z INTERNE ou EXTERNE\n"); 
	  printf("  10) Impression des images NEGATIVE ou POSITIVE\n"); 
	  printf("  11) Nom de l'image des Flux\n");
	  printf("  12) Nom de l'image des decalages spectraux\n");
	  printf("  13) Nom du fichier contenant le spectre\n");
	  printf("  14) Valeur de normalisation des flux\n");
	  printf("  15) Valeur de normalisation des decalages spectraux\n");
	  printf("\nSi aucun parametre defini, parametres par defaut :\n\n");
	  printf("   1) Dimension de l'image : 256 pixels de cote\n");
	  printf("   2) Masse relative du trou noir : 1\n");
	  printf("   3) Dimension du disque exterieur : 12 \n");
	  printf("   4) Distance de l'observateur : 100 \n");
	  printf("   5) Inclinaison par rapport au disque (en degres) : 10\n");
	  printf("   6) Observation a distance FINIE\n");
	  printf("   7) Rayonnement de disque MONOCHROMATIQUE\n");
	  printf("   8) Normalisation des flux INTERNE\n");
	  printf("   9) Normalisation des z INTERNE\n");
	  printf("  10) Impression des images NEGATIVE ou POSITIVE\n"); 
       	  printf("  11) Nom de l'image des flux : flux.pgm\n");
	  printf("  12) Nom de l'image des z : z.pgm\n");
	  printf("  13) Nom du fichier contenant le spectre : spectre.dat\n");
	  printf("  14) <non definie>\n");
	  printf("  15) <non definie>\n");
	}
    }
  
  if ((argc==14)||(argc==16))
    {
      printf("# Utilisation les valeurs definies par l'utilisateur\n");
      
      dim=atoi(argv[1]);
      m=atof(argv[2]);
      re=atof(argv[3]);
      ro=atof(argv[4]);
      tho=PI/180.*(90-atof(argv[5]));
      
      rs=2.*m;
      ri=3.*rs;
      q=-2;

      if (strcmp(argv[6],"FINIE")==0)
	{
	  dist=0;
	}
      else
	{
	  dist=1;
	}

      if (strcmp(argv[7],"MONOCHROMATIQUE")==0)
	{
	  raie=0;
	}
      else
	{
	  raie=1;
	}

      if (strcmp(argv[8],"EXTERNE")==0)
	{
	  fen=atof(argv[14]);
	}
      
      if (strcmp(argv[9],"EXTERNE")==0)
	{
	  zen=atof(argv[15]);
	}
      
    }
  else
    {
      printf("# Utilisation les valeurs par defaut\n");
      
      dim=256;
      m=1.;
      rs=2.*m;
      ri=3.*rs;
      re=12.;
      ro=100.;
      tho=PI/180.*80;
      q=-2;
      dist=0;
      raie=0;
    }
      
      printf("# Dimension de l'image : %i\n",dim);
      printf("# Masse : %f\n",m);
      printf("# Rayon singularite : %f\n",rs);
      printf("# Rayon interne : %f\n",ri);
      printf("# Rayon externe : %f\n",re);
      printf("# Distance de l'observateur : %f\n",ro);
      printf("# Inclinaison a la normale en radian : %f\n",tho);
  
  for (i=0;i<nbr;i++)
    {
      fx[i]=0.;
      nt[i]=0;
    }  

  zp=(double**)calloc(dim,sizeof(double*));
  zp[0]=(double*)calloc(dim*dim,sizeof(double));
  
  fp=(double**)calloc(dim,sizeof(double*));
  fp[0]=(double*)calloc(dim*dim,sizeof(double));

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
  thx=asin(bmx/ro);
  pc=0;

  if (raie==0)
    {
      bss=2;
    }
  else
    {
      bss=3e21;
    }

  for (n=1;n<=nmx;n++)
    {     
      h=PI/500.;
      d=stp*n;

      if (dist==1)
	{
	  db=bmx/(double)nmx;
	  b=db*(double)n;
	  up=0.;
	  vp=1.;
	}
      else
	{
	  thi=thx/(double)nmx*(double)n;
	  db=ro*sin(thi)-b;
	  b=ro*sin(thi);
	  up=sin(thi);
	  vp=cos(thi);
	}
      
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
      
      for (i=nh+1;i<2000;i++)
	{
	  rp[i]=0.; 
	}
      
      imx=(int)(8*d);
      
      for (i=0;i<=imx;i++)
	{
	  phi=2.*PI/(double)imx*(double)i;
	  phd=atanp(cos(phi)*sin(tho),cos(tho));
	  phd=fmod(phd,PI);
	  ii=0;
	  tst=0;
	  
	  do
	    {
	      php=phd+(double)ii*PI;
	      nr=php/h;
	      ni=(int)nr;

	      if ((double)ni<nh)
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
		  impact(d,phi,dim,r,b,tho,m,zp,fp,nt,fx,q,db,h,bss,raie);
		}
	      
	      ii++;
	    } while ((ii<=2)&&(tst==0));
	}
    }

  fmx=fp[0][0];
  zmx=zp[0][0];
  
  for (i=0;i<dim;i++) for (j=0;j<dim;j++)
    {
      if (fmx<fp[i][j])
	{
	  fmx=fp[i][j];
	}
      
      if (zmx<zp[i][j])
	{
	  zmx=zp[i][j];
	}
    }

  printf("\nLe flux maximal detecte est de %f",fmx);
  printf("\nLe decalage spectral maximal detecte est de %f\n\n",zmx);

  if (strcmp(argv[8],"EXTERNE")==0)
    {
      fmx=fen;
    }

  if (strcmp(argv[9],"EXTERNE")==0)
    {  
      zmx=zen;
    }

  for (i=0;i<dim;i++) for (j=0;j<dim;j++)
    {
      zcl=(int)(255/zmx*zp[i][dim-1-j]);
      fcl=(int)(255/fmx*fp[i][dim-1-j]);

      if (strcmp(argv[8],"NEGATIVE")==0)
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

  flux_tot=0;
  impc_tot=0;

  for (i=1;i<nbr;i++)
    {
      flux_tot+=fx[i];
      impc_tot+=nt[i];
    }

  for (i=1;i<nbr;i++)
    {
      tableau[0][i]=bss*i/nbr;
      tableau[1][i]=fx[i]/flux_tot;
      tableau[2][i]=(double)nt[i]/(double)impc_tot;
    }

  if ((argc==14)||(argc==16))
   {
     sauvegarde_pgm(argv[11],ifp,dim);
     sauvegarde_pgm(argv[12],izp,dim);
     sauvegarde_dat(argv[13],tableau,raie);
   }
  else
    {
      sauvegarde_pgm("z.pgm",izp,dim);
      sauvegarde_pgm("flux.pgm",ifp,dim);
      sauvegarde_dat("spectre.dat",tableau,raie);
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


