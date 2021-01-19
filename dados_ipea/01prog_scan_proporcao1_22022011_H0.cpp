
// Programa exemplo para o DEV
// 

// Inclusão de bibliotecas
#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>
#include <time.h>

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define pi (3.1415926535897932)

// Definindo Constantes
#define max 1000
#define nxx 100000
#define nsim 10000

// Declarando as funcoes
double logpow( double x1, double x2 );
double  power( double x1, double x2 );
double rnorm(double mu, double sigma);
double sqr( double x );
int multinomial(long n1);
float ran1(long idum);
void ordena0(long esq, long dir);
void particao0(long esq, long dir, long *i, long *j);
void ordena2(long esq, long dir);
void particao2(long esq, long dir, long *i, long *j);

clock_t t_inicio, t_fim;

// Criando Arquivos
FILE *dados;
FILE *saidas;
FILE *saidas_H0;


// declarando Variáveis
int i,j, id[max], pop[max], cas[max], soma, m, pop0, cas0, id0;
double  cx[max], cy[max], x0, y0_;
double  dist[max][max];
double limite;
double popt, cast;
int print1=1;// imprime
int print0=0;// não imprime
long idum;
double xxx[nxx];
long rapid[max];
long cassim[max];// vetor dos casos simulados
double distord[max][max];
int inddistord[max][max];
int popz, casz, w, popz2, casz2;
double miz, miz2;
double llr, llrz, llrz2;
int centro, raio;
int cast2;
double multin[max];    // vetor que auxilia a geração multinomial
long p, imc, nrep;  

double popac, vetllrmax[nsim]; 

int vetcentro[nsim];
int vetraio[nsim];
int vetcasz[nsim];
double vetmiz[nsim];
int vetpopz[nsim];
double pvalor;

// programa Principal
main()
{  
      t_inicio= clock();
      m=853; //número de regiões do mapa
      nrep=1001;// número de réplicas de Monte Carlo
      limite=0.30;//tamanho do raio de cada centróide
      dados=fopen("basemg.txt", "r");//abre o arquivo para leitura
      //dados=fopen("basemg_13regioes.txt", "r");//abre o arquivo para leitura
      //saidas=fopen("resultados_H0.txt","w");//abre o arquivo para escrita
      saidas_H0=fopen("result_H0.txt","w");//abre o arquvio para escrita das verossimilhancas sob H0
      popt=0.0;//população total
      cast=0.0;//casos totais
      
    if(print0)  fprintf(saidas,"id  cx      cy      pop      cas \n");
      
      for (i=1;i<=m;i++)
      {
         fscanf(dados,"%d %lg %lg %d %d",&id0, &x0, &y0_, &pop0, &cas0);//escanear os dados do arquivo aberto e associá-los as respectivas variáveis
         id[i]=id0; //identificador da regiao
         cx[i]=x0;
         cy[i]=y0_;
         pop[i]=pop0;
         cas[i]=cas0;
         popt+=pop[i];//popt= popt+pop[i]
         cast+=cas[i];//cast= cast+cas[i]
        if(print0) fprintf(saidas,"%2d %6.2f %6.2f %6d %6d \n" ,  id[i], cx[i], cy[i], pop[i], cas[i]);
       } 
       cast2=int(cast);
       if(print0) fprintf(saidas,"\n popt=%f cast=%f\n\n",popt,cast);//não imprime saídas quando print for 0.
       // Calculo de distancias euclidianas
       for (i=1; i<m;i++)
       {
           dist[i][i]=0.0;
           for (j=i+1;j<=m;j++)
           {
             dist[i][j]=sqrt((cx[i]-cx[j])*(cx[i]-cx[j])+(cy[i]-cy[j])*(cy[i]-cy[j])); 
             dist[j][i]=dist[i][j]; 
            }        
       }      
 // imprimindo a matriz e distancia 
       for (i=1; i<=m;i++)
       {
             for (j=1;j<=m;j++)
           {        
              if(print0) fprintf(saidas,"%5.2f ", dist[i][j]); 
           }
         if(print0) fprintf(saidas,"\n");  
       }  
       if(print0) fprintf(saidas,"\n"); 
       
   //distâncias e índices ordenados entre as m regiões
    for (i=1; i<=m;i++)
       {
           for (j=1;j<=m;j++)
           {
              xxx[j]=dist[i][j];
              rapid[j]=j; 
            }  
            ordena0(1,m);
            for (j=1;j<=m;j++)
            {
                distord[i][j]=xxx[j];
                inddistord[i][j]=rapid[j];
            }         
       } 
       //imprime a matriz de dist. ordenada
        for (i=1; i<=m;i++)
       {
             for (j=1;j<=m;j++)
           {        
              if(print0) fprintf(saidas,"%5.2f ", distord[i][j]); 
           }
         if(print0) fprintf(saidas,"\n");  
       }  
       if(print0) fprintf(saidas,"\n");  
       //imprime a matriz dos indices que ordena as dist.
        for (i=1; i<=m;i++)
       {
             for (j=1;j<=m;j++)
           {        
              if(print0) fprintf(saidas,"%d ", inddistord[i][j]); 
           }
         if(print0) fprintf(saidas,"\n");  
       }  
       if(print0) fprintf(saidas,"\n"); 

  popac=0.0;
   for(p=1;p<=m;p++) // auxilia à funcao multinomial
   { 
       popac+=pop[p];
       multin[p]=popac/popt;
   }
       
  pvalor=0.0;

 //Grande loop para simulação Monte Carlo  
    
for (imc=0; imc<nrep; imc++)
{
    llrz2=0.0;
    vetllrmax[imc]=0.0;
   
    if(imc > 0)  /* simulacao via multinomial */
    {
//       cast2=12676;
       multinomial(cast2);
       for(p=1;p<=m;p++) 
       {
          cas[p]=cassim[p];
          if (print0) fprintf(saidas,"imc=%d cas[%d]=%d \n", imc, p, cas[p]);        
       }
     }
       
       //calcula as LLR(logaritmo da razao de verossimilhança)de cada zona z
       for(i=1;i<=m;i++)
       {
                        popz=0;
                        casz=0;
                        w=1;
                        do{
                           popz+=pop[inddistord[i][w]];
                           casz+=cas[inddistord[i][w]];
                           miz=cast*(popz/(float)popt);
                           if (casz>miz)
                           {
                              llrz= logpow(casz/miz,casz)+ logpow((cast-casz)/(cast-miz),cast-casz); // 1
                              //llrz=pow((casz/(float)popz)+1.0,sqrt(popz));  // 2
                              //llrz=pow((casz/(float)popz)+1.0,casz);        // 3
                              //llrz=pow((casz/(float)popz)+1.0,sqrt(casz));  // 4
                           }
                           else
                           {
                               llrz=0.0;
                           }
                           if (print0) fprintf(saidas,"%12d %6d %6.2f %8.4f \n", popz, casz, miz, llrz);// 8 representa os espaços, .4 representa o número de casas decimais  
                           // armazenando o llr maximo
                           if (print0) fprintf(saidas,"====================================== \n");
                           if (llrz>llrz2)
                           {
                              llrz2=llrz;
                              centro=i;
                              raio=w;
                              casz2=casz;
                              popz2=popz;
                              miz2=miz;                     
                              //imprime todas as zonas cujo llr  e diferente de zero
                              if (print0) fprintf(saidas,"centro=%d raio=%d casz=%6d miz=%6.2f popz=%12d llrz=%f regioes: ", centro, raio,casz, miz, popz, llrz);            
                              for (j=1;j<=raio;j++)
                              {
                                  if (print0) fprintf(saidas,"%d-",inddistord[centro][j]);
                              }
                              if (print0) fprintf(saidas,"\n");
                            } 
                                                    
                           w++;        
                           }while(popz<limite*popt);
                                                     
        }
 
 // imprime as informaçoes apenas da zona que maximiza a llr
         if (print0) fprintf(saidas,"====================================== \n");
         if (print0) fprintf(saidas,"centro=%d raio=%d llrz=%f regioes: ", centro, raio, llrz2);            
         for (j=1;j<=raio;j++)
         {
            if (print0) fprintf(saidas,"%d-",inddistord[centro][j]);
         }
         if (print0) fprintf(saidas,"\n"); 
      
      vetllrmax[imc]=llrz2;
      vetcentro[imc]=centro;
      vetraio[imc]=raio;
      vetcasz[imc]=casz2;
      vetmiz[imc]=miz2;
      vetpopz[imc]=popz2;
      if (print1) fprintf(saidas,"imc=%d centro=%d raio=%d casz=%d miz=%f popz=%d llrz=%f \n", imc, vetcentro[imc], vetraio[imc], vetcasz[imc], vetmiz[imc], vetpopz[imc], vetllrmax[imc]);      
      if (print0) fprintf(saidas," %f \n", vetllrmax[imc]);      
         
         for (j=1;j<=vetraio[imc];j++)
         {
            if (print1) fprintf(saidas,"%d-",inddistord[vetcentro[imc]][j]);
         }
      if (print1) fprintf(saidas,"\n \n"); 
      if (print0) printf(" simulacao %d \n", imc); 
      
       if(imc>0)
   {
      if (print1) fprintf(saidas_H0,"%f \n", vetllrmax[imc]);
      if(vetllrmax[imc]>vetllrmax[0]) 
        pvalor++;
   }

if (print1) printf("%d \n", imc);
}// Fecha o loop de Monte Carlo 

pvalor=pvalor/nrep;
if (print1) fprintf(saidas, "pvalor=%f", pvalor);

t_fim= clock();

printf("Tempo de execucao (em segundos): %lf \n", ((double) (t_fim-t_inicio)/CLOCKS_PER_SEC));                            
                              
   fclose(dados); 
   fclose(saidas);  
   fclose(saidas_H0);
   
   system("PAUSE");
}

// ##############################################################
//  FUNÇÕES DIVERSAS A PARTIR DAQUI
// ##############################################################

double logpow(double x1, double x2)//calcula o logaritmo de uma potência
{
  if(x2>1e-20) return( (log(x1)*x2) );
  else return(0);
}

double power(double x1, double x2)// calcula uma potencia x1^x2
{
 /* power(2,3)=8 */
  return( exp(log(x1)*x2) );
}

double sqr( double x )
{
  return (x*x);
}

float ran1(long idum0)// gera número aleatório entre (0,1)
{
   int j;
   long k;
   float temp;
   static long iy=0;
   static long iv[NTAB];
   if( idum0<=0 || !iy ){
     if( -(idum0)<1 ) idum0=1;
     else idum0 = -(idum0);
     for(j=NTAB+7;j>=0;j--){
       k= (idum0)/IQ;
       idum0= IA * (idum0 - k * IQ) - IR * k;
       if(idum0 < 0) idum0 +=IM;
       if(j<NTAB) iv[j]=idum0;
     }
     iy=iv[0];
    }
    k= (idum0)/IQ;
    idum0 = IA * (idum0 - k * IQ) - IR * k;
    if(idum0 < 0) idum0 += IM;
    j=iy/NDIV;
    iy=iv[j];
    iv[j] = idum0;
    temp=AM*iy;
    idum=idum0;
    if(temp>RNMX)return RNMX;
    else return temp;
}


void particao0(long esq, long dir, long *i, long *j){
     double temp,pivot;
     long tempi;
     *i=esq;
     *j=dir;
     pivot=xxx[((*i)+(*j))/2]; //encontra valor do pivot
     do{
       while (xxx[*i]<pivot)
         (*i)++;
       while (xxx[*j]>pivot)
         (*j)--;
       if ((*i)<=(*j)){
         //coloca distxy em ordem crescente
         temp=xxx[*i];
         xxx[*i]=xxx[*j];
         xxx[*j]=temp;
         tempi=rapid[*i];
         rapid[*i]=rapid[*j];
         rapid[*j]=tempi;
         (*i)++;
         (*j)--;
       }
     }while ((*i)<=(*j));
}

void ordena0(long esq, long dir)
    {
      long i,j;
      particao0(esq,dir,&i,&j);
      if (esq<j)
        ordena0(esq,j);
      if (i<dir)
        ordena0(i,dir);
    }

void particao2(long esq, long dir, long *i, long *j){
     double temp,pivot;
     long tempi;
     *i=esq;
     *j=dir;
     pivot=xxx[((*i)+(*j))/2]; //encontra valor do pivot
     do{
       while (xxx[*i]<pivot)
         (*i)++;
       while (xxx[*j]>pivot)
         (*j)--;
       if ((*i)<=(*j)){
         //coloca distxy em ordem crescente
         temp=xxx[*i];
         xxx[*i]=xxx[*j];
         xxx[*j]=temp;       
         (*i)++;
         (*j)--;
       }
     }while ((*i)<=(*j));
}

void ordena2(long esq, long dir)
    {
      long i,j;
      particao2(esq,dir,&i,&j);
      if (esq<j)
        ordena2(esq,j);
      if (i<dir)
        ordena2(i,dir);
    }

int multinomial(long n)
{
  long i,j;
  for(j=1;j<=m;j++)
    cassim[j]=0;
  for(i=1;i<=n;i++)
    xxx[i]=ran1(idum);
  ordena2(1,n);
  j=1;
  for(i=1;i<=n;i++)
  {
    while(xxx[i]>multin[j]) j++;
    cassim[j]++;
  }
  return(0);
}


      
