/* 
 * File:   iocharmm.hpp
 * Author: cazade
 *
 * Created on June 26, 2013, 3:20 PM
 */

#ifndef IOCHARMM_HPP
#define	IOCHARMM_HPP

#include <cstdlib>
#include <cstdio>
#include <string>
#include <vector>

#include "atom.hpp"

class Iocharmm
{
public:
  
  static void readInput(FILE *input,std::vector<std::string>& listDcd,
			std::vector<Atom>& atom,int& nAtom,int& firstFrame,
			int& lastFrame,int& freqFrame,int& fAtom1,
			int& lAtom1,int& fAtom2,int& lAtom2,int& nDcd,
			double& rMin,double& rMax,std::vector<double>& box);
  
  static void readPsf(std::string& namePsf,int& nAtom,std::vector<Atom>& atom);
  
};

#endif