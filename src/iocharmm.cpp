#include <iostream>
#include <cstring>

#include "iocharmm.hpp"

#define _FORMAT103_ "%d %s %d %s %s %d %lf %lf %d %lf %lf %lf %lf"
#define _FORMAT102_ "%d %s %d %s %s %d %lf %lf %d %lf %lf"
#define _FORMAT101_ "%d %s %d %s %s %d %lf %lf %d"

using namespace std;

void Iocharmm::readInput(FILE *input,vector<string>& listDcd,
			 vector<Atom>& atom,int& nAtom,int& firstFrame,
			 int& lastFrame,int& freqFrame,int& fAtom1,
			 int& lAtom1,int& fAtom2,int& lAtom2,int& nDcd,
			 double& rMin,double& rMax,vector<double>& box)
{
  char buff1[1024]="", *buff2=NULL,*buff3=NULL;
  char buffName[1024];
  
  while(fgets(buff1,1024,input)!=NULL)
  {
    buff2=strtok(buff1," \n\t");
    
    if(buff2==NULL)
      continue;
    
    if(!strcmp(buff2,"psf"))
    {
      fscanf(input,"%s",buffName);
      string namePsf(buffName);
      
      readPsf(namePsf,nAtom,atom);
    }
    else if(!strcmp(buff2,"frame"))
    {
      fscanf(input,"%d %d %d",&firstFrame,&lastFrame,&freqFrame);
    }
    else if(!strcmp(buff2,"dcd"))
    {
      buff3=strtok(NULL," \n\t");
      nDcd=atoi(buff3);
      
      for(int i=0;i<nDcd;i++)
      {
	fscanf(input,"%s",buffName);
	listDcd.push_back(buffName);
      }
    }
    else if(!strcmp(buff2,"selec"))
    {
      fscanf(input,"%d %d",&fAtom1,&lAtom1);
      fscanf(input,"%d %d",&fAtom2,&lAtom2);
      fscanf(input,"%lf %lf",&rMin,&rMax);
    }
    else if(!strcmp(buff2,"pbc"))
    {
      fscanf(input,"%lf %lf %lf",&(box[0]),&(box[1]),&(box[2]));
    }
  }
}

void Iocharmm::readPsf(string& namePsf,int& nAtom,vector<Atom>& atom)
{
  
  FILE *psfFile=NULL;
  
  char buff1[1024]="", *buff2=NULL;
  
  char resn[5]="",segn[5]="",label[5]="";
  
  int idx,ires,type,imove;
  
  double q,m,ech,eha,alphadp,tholei;
  
  bool lExt(false),lDrude(false),lCmap(false),lCheq(false);

  psfFile=fopen(namePsf.c_str(),"r");

  if(psfFile==NULL)
  {
    cout<< "Cannot find PSF file: "<<namePsf<<endl;
    exit(10);
  }

  while(fgets(buff1,1024,psfFile)!=NULL)
  {
    
    if(strstr(buff1,"EXT")!=NULL)
    {
      lExt=true;
    }
    
    if(strstr(buff1,"CMAP")!=NULL)
    {
      lCmap=true;
    }
    
    if(strstr(buff1,"CHEQ")!=NULL)
    {
      lCheq=true;
    }
    
    if(strstr(buff1,"DRUDE")!=NULL)
    {
      lDrude=true;
    }
    
    if(strstr(buff1,"NATOM")!=NULL)
    {
	buff2=strtok(buff1," \n\t");
	nAtom=atoi(buff2);
	break;
    }
  }

  for(int i=0; i<nAtom;i++)
  {
    if(lCheq)
    {
      if(lDrude)
      {
	fscanf(psfFile,_FORMAT103_,&idx,segn,&ires,resn,label,&type,&q,&m,&imove,&ech,&eha,&alphadp,&tholei);
      }
      else
      {
	fscanf(psfFile,_FORMAT102_,&idx,segn,&ires,resn,label,&type,&q,&m,&imove,&ech,&eha);
      }
    }
    else
    {
      fscanf(psfFile,_FORMAT101_,&idx,segn,&ires,resn,label,&type,&q,&m,&imove);
    }
    
    atom.push_back(Atom(segn,(ires-1),resn,label,type,q,m));
    
  }

  fclose(psfFile);
}