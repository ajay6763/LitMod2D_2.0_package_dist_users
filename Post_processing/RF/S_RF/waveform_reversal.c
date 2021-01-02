//header files
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>
main()
{
  int len;
 int i,index;
 float data[5000],d[5000],time[5000],t[5000],t_1,t_2;   
len=0;
  char file[50];
  printf("enter file to work on: \n");
  scanf("%s \n",file);
  printf("%s\n",file );
  FILE *fr,*fw;
  fr=fopen(file,"r");
  fw=fopen("out.txt","w+");
  while ( !feof(fr))
	{    
		//printf("\n%d",len);
		fscanf(fr,"%f%f \n",&t_1,&t_2); 
		//printf("\n%f %f",t_1,t_2);
                data[len]=t_2;
                time[len]=t_1;
		len=len++;
	}
  

printf("\n %d\n",len);
       //printf("\nworking good\n");

for (i=0;i <len; i++)
{       //printf("working good\n");
	index=(len-1)-i;
	d[i]=data[index];
	d[i]=-1*d[i];        
	t[i]=time[i]; 
        //printf("%f\n",d[i]);
        fprintf(fw,"%f\t%f\n",t[i],d[i]);
}

printf("\nworking good\n");
       
fclose(fr);

//fprintf(fw,"%d\n",t);		
	        
fclose(fw);
}
