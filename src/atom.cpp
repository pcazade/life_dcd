/*
 * <one line to give the program's name and a brief idea of what it does.>
 * Copyright (C) 2013  Pierre-Andre Cazade <email>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <cstdlib>
#include <iostream>

#include "atom.hpp"

using namespace std;

Atom::Atom()
{
    m_ires=0;
    m_resi=0;
    m_type=0;

    m_mass=0.;
    m_charge=0.;

    m_segmentName="";
    m_residueName="";
    m_label="";
}


Atom::Atom(string& segn,int ires,string& resn,string& label,int type,double q,double m)
{
    m_ires=ires;
    m_resi=0;
    m_type=type;

    m_mass=m;
    m_charge=q;

    m_segmentName=segn;
    m_residueName=resn;
    m_label=label;
}

Atom::Atom(char segn[],int ires,char resn[],char label[],int type,double q,double m)
{
    m_ires=ires;
    m_resi=0;
    m_type=type;

    m_mass=m;
    m_charge=q;

    m_segmentName=segn;
    m_residueName=resn;
    m_label=label;
}

Atom::Atom(string& segn,int ires,int resi,string& resn,string& label,int type,double q,double m)
{
    m_ires=ires;
    m_resi=resi;
    m_type=type;

    m_mass=m;
    m_charge=q;

    m_segmentName=segn;
    m_residueName=resn;
    m_label=label;
}

string Atom::getSegn() const
{
    return(m_segmentName);
}

string Atom::getResn() const
{
    return(m_residueName);
}

string Atom::getLabel() const
{
    return(m_label);
}

int Atom::getResi() const
{
    return(m_resi);
}

int Atom::getIres() const
{
    return(m_ires);
}

int Atom::getType() const
{
    return(m_type);
}

double Atom::getMass() const
{
    return(m_mass);
}

double Atom::getCharge() const
{
    return(m_charge);
}

void Atom::setSegn(std::string& segn)
{
    m_segmentName=segn;
}

void Atom::setResn(std::string& resn)
{
    m_residueName=resn;
}

void Atom::setLabel(std::string& label)
{
    m_label=label;
}

void Atom::setResi(int resi)
{
    m_resi=resi;
}

void Atom::setIres(int ires)
{
    m_ires=ires;
}

void Atom::setType(int type)
{
    m_type=type;
}

void Atom::setMass(double m)
{
    m_mass=m;
}

void Atom::setCharge(double q)
{
    m_charge=q;
}

void Atom::atomSelec(const vector<Atom>& atom,bool** u[],double** v[],double** w[],
                     const vector<double>& box,const float x[],const float y[],
                     const float z[],const double rMin,const double rMax,
                     const int fAtom1,const int lAtom1,const int fAtom2,
                     const int lAtom2,const int frame)
{
    int ii,jj;
    double r(0.);
    double ux(0.),uy(0.),uz(0.);
    double comx(0.),comy(0.),comz(0.);
    double comox(0.),comoy(0.),comoz(0.);
    double mass(0.);

    int nAtomC=lAtom1-fAtom1+1;
    int nAtomS=lAtom2-fAtom2+1;

    int nResidC(atom.at(lAtom1).getIres()+1);
    int nResidS(atom.at(lAtom2).getIres()+1);

    int nAtomCr=nAtomC/nResidC;
    int nAtomSr=nAtomS/nResidS;

    for(int i=0; i<nResidC; i++)
    {
        for(int j=0; j<nResidS; j++)
        {
            u[frame][i][j]=false;
        }
    }

    for(int i=0; i<nResidC; i++)
    {
        comx=0.;
        comy=0.;
        comz=0.;

        comox=0.;
        comoy=0.;
        comoz=0.;

        mass=0.;

        if(frame>0)
        {
	    
	    ii=i*nAtomCr;
            for(int j=0; j<nAtomCr; j++)
            {
                jj=fAtom1+ii;

                comx+=atom.at(jj).getMass()*x[jj];
                comy+=atom.at(jj).getMass()*y[jj];
                comz+=atom.at(jj).getMass()*z[jj];

                comox+=atom.at(jj).getMass()*v[frame-1][ii][0];
                comoy+=atom.at(jj).getMass()*v[frame-1][ii][0];
                comoz+=atom.at(jj).getMass()*v[frame-1][ii][0];

                mass+=atom.at(jj).getMass();
		
		ii++;
            }

            comx/=mass;
            comy/=mass;
            comz/=mass;

            comox/=mass;
            comoy/=mass;
            comoz/=mass;

            ux=comx-comox;
            uy=comy-comoy;
            uz=comz-comoz;

            ux-=box.at(0)*nint(ux/box.at(0));
            uy-=box.at(1)*nint(uy/box.at(1));
            uz-=box.at(2)*nint(uz/box.at(2));

            comx=ux+comox;
            comy=uy+comoy;
            comz=uz+comoz;
	    
	    ii=i*nAtomCr;
            for(int j=0; j<nAtomCr; j++)
            {
                jj=fAtom1+ii;

                ux=x[jj]-comx;
                uy=y[jj]-comy;
                uz=z[jj]-comz;

                ux-=box.at(0)*nint(ux/box.at(0));
                uy-=box.at(1)*nint(uy/box.at(1));
                uz-=box.at(2)*nint(uz/box.at(2));

                v[frame][ii][0]=ux+comx;
                v[frame][ii][1]=uy+comy;
                v[frame][ii][2]=uz+comz;
		
		ii++;
            }

        } // if(frame>0)
        else
        {
	    ii=i*nAtomCr;
            for(int j=0; j<nAtomCr; j++)
            {
                jj=fAtom1+ii;

                comx+=atom.at(jj).getMass()*x[jj];
                comy+=atom.at(jj).getMass()*y[jj];
                comz+=atom.at(jj).getMass()*z[jj];

                mass+=atom.at(jj).getMass();
		
		ii++;
            }

            comx/=mass;
            comy/=mass;
            comz/=mass;

            comx-=box.at(0)*nint(comx/box.at(0));
            comy-=box.at(1)*nint(comy/box.at(1));
            comz-=box.at(2)*nint(comz/box.at(2));
	    
	    ii=i*nAtomCr;
            for(int j=0; j<nAtomCr; j++)
            {
                jj=fAtom1+ii;

                ux=x[jj]-comx;
                uy=y[jj]-comy;
                uz=z[jj]-comz;

                ux-=box.at(0)*nint(ux/box.at(0));
                uy-=box.at(1)*nint(uy/box.at(1));
                uz-=box.at(2)*nint(uz/box.at(2));

                v[frame][ii][0]=ux+comx;
                v[frame][ii][1]=uy+comy;
                v[frame][ii][2]=uz+comz;
		
		ii++;
            }

        }  // else if(frame>0)

    } // for(int i=0;i<nResidC;i++)

    for(int i=0; i<nResidS; i++)
    {
        comx=0.;
        comy=0.;
        comz=0.;

        comox=0.;
        comoy=0.;
        comoz=0.;

        mass=0.;

        if(frame>0)
        {
	    
	    ii=i*nAtomSr;
            for(int j=0; j<nAtomSr; j++)
            {
                jj=fAtom2+ii;

                comx+=atom.at(jj).getMass()*x[jj];
                comy+=atom.at(jj).getMass()*y[jj];
                comz+=atom.at(jj).getMass()*z[jj];

                comox+=atom.at(jj).getMass()*w[frame-1][ii][0];
                comoy+=atom.at(jj).getMass()*w[frame-1][ii][1];
                comoz+=atom.at(jj).getMass()*w[frame-1][ii][2];

                mass+=atom.at(jj).getMass();
		
		ii++;
            }

            comx/=mass;
            comy/=mass;
            comz/=mass;

            comox/=mass;
            comoy/=mass;
            comoz/=mass;

            ux=comx-comox;
            uy=comy-comoy;
            uz=comz-comoz;

            ux-=box.at(0)*nint(ux/box.at(0));
            uy-=box.at(1)*nint(uy/box.at(1));
            uz-=box.at(2)*nint(uz/box.at(2));

            comx=ux+comox;
            comy=uy+comoy;
            comz=uz+comoz;
	    
	    ii=i*nAtomSr;
            for(int j=0; j<nAtomSr; j++)
            {
                jj=fAtom2+ii;

                ux=x[jj]-comx;
                uy=y[jj]-comy;
                uz=z[jj]-comz;

                ux-=box.at(0)*nint(ux/box.at(0));
                uy-=box.at(1)*nint(uy/box.at(1));
                uz-=box.at(2)*nint(uz/box.at(2));

                w[frame][ii][0]=ux+comx;
                w[frame][ii][1]=uy+comy;
                w[frame][ii][2]=uz+comz;
		
		ii++;
            }

        } // if(frame>0)
        else
        {
	    ii=i*nAtomSr;
            for(int j=0; j<nAtomSr; j++)
            {
                jj=fAtom2+ii;

                comx+=atom.at(jj).getMass()*x[jj];
                comy+=atom.at(jj).getMass()*y[jj];
                comz+=atom.at(jj).getMass()*z[jj];

                mass+=atom.at(jj).getMass();
		
		ii++;
            }

            comx/=mass;
            comy/=mass;
            comz/=mass;

            comx-=box.at(0)*nint(comx/box.at(0));
            comy-=box.at(1)*nint(comy/box.at(1));
            comz-=box.at(2)*nint(comz/box.at(2));
	    
	    ii=i*nAtomSr;
            for(int j=0; j<nAtomSr; j++)
            {
                jj=fAtom2+ii;

                ux=x[jj]-comx;
                uy=y[jj]-comy;
                uz=z[jj]-comz;

                ux-=box.at(0)*nint(ux/box.at(0));
                uy-=box.at(1)*nint(uy/box.at(1));
                uz-=box.at(2)*nint(uz/box.at(2));

                w[frame][ii][0]=ux+comx;
                w[frame][ii][1]=uy+comy;
                w[frame][ii][2]=uz+comz;
		
		ii++;
            }

        }  // else if(frame>0)

    } // for(int i=0;i<nResidS;i++)

    for(int i=fAtom1; i<=lAtom1; i++)
    {
        for(int j=fAtom2; j<=lAtom2; j++)
        {
            ux=x[j]-x[i];
            uy=y[j]-y[i];
            uz=z[j]-z[i];

            ux-=box.at(0)*nint(ux/box.at(0));
            uy-=box.at(1)*nint(uy/box.at(1));
            uz-=box.at(2)*nint(uz/box.at(2));

            r=(ux*ux)+(uy*uy)+(uz*uz);

            if( (r>rMin) && (r<=rMax) )
            {
                u[frame][atom.at(i).getIres()][atom.at(j).getIres()]=true;
            }
        }
    }

}

void Atom::lifeComp(bool** u[], vector<int>&life,vector<int>&life2,
                    vector< vector<int> >& time,vector< vector<int> >& time2,
                    vector< vector<bool> >& remove,vector< vector<bool> >& remove2,
                    const int nResidC,const int nResidS,const int nFrame)
{

    for(int i=0; i<nFrame; i++)
    {

        for(int j=0; j<nResidC; j++)
        {

            for(int k=0; k<nResidS; k++)
            {

                if(u[i][j][k])
                {
                    time.at(j).at(k)++;
                    life.at(time.at(j).at(k))++;

                    if(remove2.at(j).at(k))
                    {
                        for(int l=0; l<=time2.at(j).at(k); l++)
                        {
                            life2.at(l)--;
                        }

                        remove2.at(j).at(k)=false;

                    } // if(remove2.at(j).at(k))

                    time2.at(j).at(k)=-1;

                    if(i==0)
                        remove.at(j).at(k)=true;

                } // if(u[i][j][k])
                else
                {
                    time2.at(j).at(k)++;
                    life2.at(time2.at(j).at(k))++;

                    if(remove.at(j).at(k))
                    {
                        for(int l=0; l<=time.at(j).at(k); l++)
                        {
                            life.at(l)--;
                        }

                        remove.at(j).at(k)=false;

                    } // if(remove2.at(j).at(k))

                    time.at(j).at(k)=-1;

                    if(i==0)
                        remove2.at(j).at(k)=true;

                } // if(u.at(i).at(j).at(k)) else

            } // for(int k=0;k<nResidS;k++)

        } // for(int j=0;j<nResidC;j++)

    } // for(int i=0;i<nFrame;i++)


    for(int j=0; j<nResidC; j++)
    {

        for(int k=0; k<nResidS; k++)
        {

            for(int i=0; i<=time.at(j).at(k); i++)
            {
                life.at(i)--;
            } // for(int i=0;i<=time.at(j).at(k);i++)

            for(int i=0; i<=time2.at(j).at(k); i++)
            {
                life2.at(i)--;
            } // for(int i=0;i<=time2.at(j).at(k);i++)

        } // for(int k=0;k<nResidS;k++)

    } // for(int j=0;j<nResidC;j++)

}

void Atom::corrComp(const vector<Atom>& atom,bool** u[],double** v[],
                    double** w[],vector<double>&corr,vector<double>&nCorr,
                    const vector<double>& box,const int fAtom1,const int lAtom1,
                    const int fAtom2,const int lAtom2,const int nFrame)
{

    bool test(false);

    int tt(0);

    int kk,ll,m,n,mm,nn;

    double mass(0.);
    double ux(0.),uy(0.),uz(0.);
    double ax(0.),ay(0.),az(0.);
    double bx(0.),by(0.),bz(0.);
    double cx(0.),cy(0.),cz(0.);
    double dx(0.),dy(0.),dz(0.);

    double acomx(0.),acomy(0.),acomz(0.);
    double bcomx(0.),bcomy(0.),bcomz(0.);

    int nAtomC=lAtom1-fAtom1+1;
    int nAtomS=lAtom2-fAtom2+1;

    int nResidC(atom.at(lAtom1).getIres()+1);
    int nResidS(atom.at(lAtom2).getIres()+1);

    int nAtomCr=nAtomC/nResidC;
    int nAtomSr=nAtomS/nResidS;

    for(int k=0; k<nResidC; k++)
    {

        for(int l=0; l<nResidS; l++)
        {

            for(int i=0; i<nFrame; i++)
            {

                if(u[i][k][l])
                {

                    acomx=0.;
                    acomy=0.;
                    acomz=0.;

                    mass=0.;
		    
		    kk=k*nAtomCr;
                    for(int m=0; m<nAtomCr; m++)
                    {
                        mm=fAtom1+kk;
                        acomx+=atom.at(mm).getMass()*v[i][kk][0];
                        acomy+=atom.at(mm).getMass()*v[i][kk][1];
                        acomz+=atom.at(mm).getMass()*v[i][kk][2];

                        mass+=atom.at(mm).getMass();
			
			kk++;
                    }

                    acomx/=mass;
                    acomy/=mass;
                    acomz/=mass;

                    cx=box.at(0)*nint(acomx/box.at(0));
                    cy=box.at(1)*nint(acomy/box.at(1));
                    cz=box.at(2)*nint(acomz/box.at(2));

                    ax=acomx-cx;
                    ay=acomy-cy;
                    az=acomz-cz;

                    bcomx=0.;
                    bcomy=0.;
                    bcomz=0.;

                    mass=0.;
		    
		    ll=l*nAtomSr;
                    for(int n=0; n<nAtomSr; n++)
                    {
                        nn=fAtom2+ll;
                        bcomx+=atom.at(nn).getMass()*w[i][ll][0];
                        bcomy+=atom.at(nn).getMass()*w[i][ll][1];
                        bcomz+=atom.at(nn).getMass()*w[i][ll][2];

                        mass+=atom.at(nn).getMass();
			
			ll++;
                    }

                    bcomx/=mass;
                    bcomy/=mass;
                    bcomz/=mass;

                    ux=bcomx-ax;
                    uy=bcomy-ay;
                    uz=bcomz-az;

                    bx=box.at(0)*nint(ux/box.at(0));
                    by=box.at(1)*nint(uy/box.at(1));
                    bz=box.at(2)*nint(uz/box.at(2));

                    dx=bx-cx;
                    dy=by-cy;
                    dz=bz-cz;

                    for(int j=i; j<nFrame; j++)
                    {

                        test=false;
			
			kk=k*nAtomCr;
                        for(int m=0; m<nAtomCr; m++)
                        {
			    
			    ll=l*nAtomSr;
                            for(int n=0; n<nAtomSr; n++)
                            {

                                ux=w[j][ll][0]-v[j][kk][0]-dx;
                                uy=w[j][ll][1]-v[j][kk][1]-dy;
                                uz=w[j][ll][2]-v[j][kk][2]-dz;

                                tt=abs(inint(ux/box.at(0)))+abs(inint(uy/box.at(1)))+abs(inint(uz/box.at(2)));

                                if(tt==0)
                                {
                                    test=true;
                                    break;
                                }
                                
                                ll++;

                            } // for(int n=0;n<nAtomSr;n++)

                            if(test)
                                break;
			    
			    kk++;

                        } // for(int m=0;m<nAtomCr;m++)

                        if(test)
                        {
                            corr[j-i]+=(double)(u[i][k][l]*u[j][k][l]);
                            nCorr[j-i]+=1.;
                        }
                        else
                        {
			    /*if((j-i)==0)
			      cout << i << ' ' << j << ' ' << k << ' ' << nAtomCr << ' ' << l << ' ' << nAtomSr << ' ' << tt<< endl;*/
                            corr[j-i]+=0.;
                            nCorr[j-i]+=1.;
                        }

                    } // for(int j=i;j<nFrame;j++)

                } // if(u[i][k][l])

            } // for(int i=0;i<nFrame;i++)

        } // for(int l=0;l<nResidS;l++)

    } // for(int k=0;k<nResidC;k++)
}

double Atom::nint(const double x)
{
    int m;
    m=( (x)>=0.?(int)((x)+0.5):(int)((x)-0.5) );
    return ( (double)m );
}

int Atom::inint(const double x)
{
    int m;
    m=( (x)>=0.?(int)((x)+0.5):(int)((x)-0.5) );
    return ( m );
}

Atom::~Atom()
{

}

