/*
 *  read_dcd : c++ class + main file example for reading a CHARMM dcd file
 *  Copyright (C) 2013  Florent Hedin
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <string>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <vector>

#include "dcd.hpp"
#include "atom.hpp"
#include "iocharmm.hpp"

using namespace std;

int main(int argc,char* argv[])
{
  FILE *input=NULL;
  
  ofstream output;
  
  bool ***u;
  double ***v,***w;
  
  int nAtom,nAtomC,nAtomS,nDcd;
  int firstFrame,lastFrame,freqFrame;
  int fAtom1,lAtom1,fAtom2,lAtom2;
  
  vector<double> box(3,0.);
  
  double rMin(0.),rMax(0.);
  
  vector<string> listDcd;
  vector<Atom> atom;
  
  if(argc>2)
  {
    input=fopen(argv[1],"r");
    output.open(argv[2]);
  }
  else
    cout<<"Usage: "<<argv[0]<< " input_file out_file"<<endl;
  
  Iocharmm::readInput(input,listDcd,atom,nAtom,firstFrame,lastFrame,
		      freqFrame,fAtom1,lAtom1,fAtom2,lAtom2,nDcd,
		      rMin,rMax,box);
  
  firstFrame--;
  lastFrame--;
  fAtom1--;
  lAtom1--;
  fAtom2--;
  lAtom2--;
  
  nAtomC=lAtom1-fAtom1+1;
  nAtomS=lAtom2-fAtom2+1;
  
  rMin=rMin*rMin;
  rMax=rMax*rMax;
  
  vector<DCD*> dcds;
  
  for(int i=0;i<nDcd;i++)
  {
    dcds.push_back(new DCD(listDcd.at(i).c_str(),false));
  }
  
  cout<<"DCDs are opened."<<endl;
  
  int nResidC(atom.at(lAtom1).getIres()+1);

  int nResidS(atom.at(lAtom2).getIres()+1);
  
  int nFrame((lastFrame-firstFrame)/freqFrame+1);
  
  u=new bool** [nFrame];
  v=new double** [nFrame];
  w=new double** [nFrame];
  
  for(int i=0;i<nFrame;i++)
  {
    u[i]=new bool* [nResidC];
    v[i]=new double* [nAtomC];
    w[i]=new double* [nAtomS];
    
    for(int j=0;j<nResidC;j++)
    {
      u[i][j]=new bool [nResidS];
      
      for(int k=0;k<nResidS;k++)
	u[i][j][k]=false;
    }
    
    for(int j=0;j<nAtomC;j++)
    {
      v[i][j]=new double [3];
      
      for(int k=0;k<3;k++)
	v[i][j][k]=0.;
    }
    
    for(int j=0;j<nAtomS;j++)
    {
      w[i][j]=new double [3];
      
      for(int k=0;k<3;k++)
	w[i][j][k]=0.;
    }
  }
  
  vector<int>life(nFrame,0);
  vector<int>life2(nFrame,0);
  
  vector< vector<int> > time(nResidC,vector<int>(nResidS,-1));
  vector< vector<int> > time2(nResidC,vector<int>(nResidS,-1));
  
  vector< vector<bool> > remove(nResidC,vector<bool>(nResidS,false));
  vector< vector<bool> > remove2(nResidC,vector<bool>(nResidS,false));
  
  vector<double>corr(nFrame,0.);
  vector<double>nCorr(nFrame,0.);
  
  float *x=NULL,*y=NULL,*z=NULL;
  
  for(int i=0;i<nDcd;i++)
  {
    
    cout<<"Frame: "<<i<<endl;
    
    for(int j=0;j<firstFrame;j++)
    {
      dcds.at(i)->read_oneFrame();
    }
    
    for(int j=firstFrame;j<=lastFrame;j++)
    {
      int jj((j-firstFrame)/freqFrame);
      
      dcds.at(i)->read_oneFrame();
      
      if( (j-firstFrame)%freqFrame == 0 )
      {
	x=dcds.at(i)->getX();
	y=dcds.at(i)->getY();
	z=dcds.at(i)->getZ();
	
	Atom::atomSelec(atom,u,v,w,box,x,y,z,rMin,rMax,
			fAtom1,lAtom1,fAtom2,lAtom2,jj);
	
      } // if( (j-firstFrame)%freqFrame == 0 )
      
    } // for(int j=firstFrame;j<=lastFrame;j++)
    
    Atom::lifeComp(u,life,life2,time,time2,remove,remove2,
		   nResidC,nResidS,nFrame);
    
    Atom::corrComp(atom,u,v,w,corr,nCorr,box,fAtom1,lAtom1,fAtom2,lAtom2,nFrame);
    
  } // for(int i=0;i<nDcd;i++)
  
  for(int i=0;i<nFrame;i++)
  {
    corr[i]/=nCorr[i];
    
    output << i << '\t' << corr[i] << '\t' << (double)life[i]/((double)life[0]) << '\t' << (double)life2[i]/((double)life2[0]) <<endl;
  }
  
  // free mem  //
  
  for(int i=0;i<nDcd;i++)
  {
    delete dcds.at(i);
  }
  
  for(int i=0;i<nFrame;i++)
  {
    
    for(int j=0;j<nResidC;j++)
      delete [] u[i][j];
    
    for(int j=0;j<nAtomC;j++)
      delete [] v[i][j];
    
    for(int j=0;j<nAtomS;j++)
      delete [] w[i][j];
    
    delete [] u[i];
    delete [] v[i];
    delete [] w[i];
  }
  
  /*delete [] u;
  delete [] v;
  delete [] w;*/
  
  output.close();
  fclose(input);
  
  return EXIT_SUCCESS;
}
