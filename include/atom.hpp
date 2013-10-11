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

#ifndef ATOM_H
#define ATOM_H

#include <string>
#include <vector>

class Atom
{
public:
  
  Atom();
  
  Atom(std::string& segn,int ires,std::string& resn,std::string& label,int type,
       double q,double m);
  
  Atom(char segn[],int ires,char resn[],char label[],int type,double q,double m);
  
  Atom(std::string& segn,int ires,int resi,std::string& resn,std::string& label,int type,
       double q,double m);
  
  std::string getSegn() const;
  std::string getResn() const;
  std::string getLabel() const;
  
  int getResi() const;
  int getIres() const;
  int getType() const;
  
  double getMass() const;
  double getCharge() const;
  
  void setSegn(std::string& segn);
  void setResn(std::string& resn);
  void setLabel(std::string& label);
  
  void setResi(int resi);
  void setIres(int ires);
  void setType(int type);
  
  void setMass(double m);
  void setCharge(double q);
  
  static void atomSelec(const std::vector<Atom>& atom,bool** u[],double** v[],double** w[],
			const std::vector<double>& box,const float x[],const float y[],
			const float z[],const double rMin,const double rMax,
			const int fAtom1,const int lAtom1,const int fAtom2,
			const int lAtom2,const int frame);
  
  static void lifeComp(bool** u[], std::vector<int>&life,std::vector<int>&life2,
		      std::vector< std::vector<int> >& time,std::vector< std::vector<int> >& time2,
		      std::vector< std::vector<bool> >& remove,std::vector< std::vector<bool> >& remove2,
		      const int nResidC,const int nResidS,const int nFrame);
  
  static void corrComp(const std::vector<Atom>& atom,bool** u[],double** v[],
		       double** w[],std::vector<double>&corr,std::vector<double>&nCorr,
		       const std::vector<double>& box,const int fAtom1,const int lAtom1,
		       const int fAtom2,const int lAtom2,const int nFrame);
  
  static double nint(const double x);
  
  static int inint(const double x);
  
  ~Atom();
  
private:
  
  int m_ires;
  int m_resi;
  int m_type;
  
  double m_mass;
  double m_charge;
  
  std::string m_segmentName;
  std::string m_residueName;
  std::string m_label;

};

#endif // ATOM_H
