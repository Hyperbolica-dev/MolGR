#pragma once

/**********************************************************************
forcefielduff.h - UFF force field.

Copyright (C) 2007 by Geoffrey Hutchison
Some portions Copyright (C) 2006-2007 by Tim Vandermeersch

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <tuple>

#include <openbabel/forcefield.h>
#include <openbabel/base.h>
#include <openbabel/mol.h>

namespace OpenBabel
{
  class MolgrFFBondCalculationUFF : public OBFFCalculation2
  {
    public:
      double bt = 0.0; // bond order (e.g., 1.41 for amide)
      double kb = 0.0;
      double r0 = 0.0;
      double rab = 0.0;
      double delta = 0.0;

      template<bool> void Compute();
  };

  class MolgrFFAngleCalculationUFF : public OBFFCalculation3
  {
    public:
      int at = 0; //angletype (ATIJK)
      bool linear = false;
      double ka = 0.0;
      double theta0 = 0.0;
      double theta = 0.0;
      double delta = 0.0;
      double c0 = 0.0;
      double c1 = 0.0;
      double c2 = 0.0;
      double zi = 0.0;
      double zk = 0.0;
      double rij = 0.0;
      double rjk = 0.0;
      double rik = 0.0;
      double cosT0 = 0.0; // cos theta0
      int coord = 0;
      int n = 0;

      template<bool> void Compute();
  };

  class MolgrFFTorsionCalculationUFF : public OBFFCalculation4
  {
    public:
      int n = 0;
      double tt = 0.0; //torsiontype (i.e. b-c bond order)
      double V = 0.0;
      double tor = 0.0;
      double cosNPhi0 = 0.0;

      template<bool> void Compute();

  };

  class MolgrFFOOPCalculationUFF : public OBFFCalculation4
  {
    public:
      double koop = 0.0;
      double angle = 0.0;
      double c0 = 0.0;
      double c1 = 0.0;
      double c2 = 0.0;

      template<bool> void Compute();
  };

  class MolgrFFVDWCalculationUFF : public OBFFCalculation2
  {
    public:
      bool is14 = false;
      bool samering = false;
      double ka = 0.0;
      double kaSquared = 0.0;
      double Ra = 0.0;
      double kb = 0.0;
      double Rb = 0.0;
      double kab = 0.0;
      double rab = 0.0;

      template<bool> void Compute();
  };

  class MolgrFFElectrostaticCalculationUFF : public OBFFCalculation2
  {
    public:
      double qq = 0.0;
      double rab = 0.0;

      template<bool> void Compute();
  };

  // Class MolgrForceFieldUFF
  // class introduction in forcefieldUFF.cpp
  class MolgrForceFieldUFF: public OBForceField
  {
  protected:
    //!  Parses the parameter file
    bool ParseParamFile();
    //!  Sets atomtypes to UFF types in _mol
    bool SetTypes();
    //!  Fill OBFFXXXCalculation vectors
    bool SetupCalculations();
    //! Setup pointers in OBFFXXXCalculation vectors
    bool SetupPointers();
    bool SetupVDWCalculation(OBAtom *a, OBAtom *b, MolgrFFVDWCalculationUFF &vdwcalc);
    //!  By default, electrostatic terms are disabled
    //!  This is discouraged, since the parameterization is not designed for it
    //!  But if you want, we give you the option.
    bool SetupElectrostatics();
    //! Same as OBForceField::GetParameter, but simpler
    OBFFParameter* GetParameterUFF(std::string a, std::vector<OBFFParameter> &parameter);

    // OBFFParameter vectors to contain the parameters
    std::vector<OBFFParameter> _ffparams;
    std::unordered_map<std::string, std::size_t> _ffparam_index;

    // OBFFXXXCalculationYYY vectors to contain the calculations
    std::vector<MolgrFFBondCalculationUFF>          _bondcalculations;
    std::vector<MolgrFFAngleCalculationUFF>         _anglecalculations;
    std::vector<MolgrFFTorsionCalculationUFF>       _torsioncalculations;
    std::vector<MolgrFFOOPCalculationUFF>           _oopcalculations;
    std::vector<MolgrFFVDWCalculationUFF>           _vdwcalculations;
    std::vector<MolgrFFElectrostaticCalculationUFF> _electrostaticcalculations;

    OBFFConstraints constraints_;
    unsigned int fix_atom_ = 0;
    unsigned int ignore_atom_ = 0;
    bool use_atom_typing_cache_ = false;
    std::string atom_typing_cache_key_;

    void ActivateThreadLocalInstance();
    static MolgrForceFieldUFF* ActiveInstance();
    static thread_local MolgrForceFieldUFF* active_instance_;

  public:
    //! Constructor
    MolgrForceFieldUFF() : OBForceField("", false)
    {
      _validSetup = false;
      _init = false;
      _rvdw = 7.0;
      _rele = 15.0;
      _epsilon = 1.0; // electrostatics not used
      _pairfreq = 10;
      _cutoff = false;
      _linesearch = LineSearchType::Newton2Num;
    }

    //! Destructor
    ~MolgrForceFieldUFF() override;

     //!Clone the current instance. May be desirable in multithreaded environments
    MolgrForceFieldUFF* MakeNewInstance() override
    {
       return new MolgrForceFieldUFF();
    }

    bool Setup(OBMol &mol);
    bool Setup(OBMol &mol, OBFFConstraints &constraints);

    //! Assignment
    MolgrForceFieldUFF &operator = (MolgrForceFieldUFF &);

    OBFFConstraints& GetConstraints();
    void SetConstraints(OBFFConstraints& constraints);
    void SetFixAtom(int index);
    void UnsetFixAtom();
    void SetIgnoreAtom(int index);
    void UnsetIgnoreAtom();
    static bool IgnoreCalculation(int a, int b);
    static bool IgnoreCalculation(int a, int b, int c);
    static bool IgnoreCalculation(int a, int b, int c, int d);
    void ConfigureAtomTypingCache(bool enabled, std::string cache_key);
    void CopyPerceivedHybridizationTo(OBMol &mol) const;
    std::vector<std::string> DebugAtomTypes() const;

    //! Get the description for this force field
    const char* Description()
    {
      return "Universal Force Field.";
    }

    //! Get the unit in which the energy is expressed
    std::string GetUnit()
      {
        return std::string("kJ/mol");  // Note that we convert from kcal/mol internally
      }

    //! \return that analytical gradients are implemented for UFF
    bool HasAnalyticalGradients() { return true; }

    //! \return total energy
    double Energy(bool gradients = true);
    //! \return the bond stretching energy
    template<bool> double E_Bond();
    double E_Bond(bool gradients = true)
    {
      return gradients ? E_Bond<true>() : E_Bond<false>();
    }
    //! Returns the angle bending energy
    template<bool> double E_Angle();
    double E_Angle(bool gradients = true)
    {
      return gradients ? E_Angle<true>() : E_Angle<false>();
    }
    //! Returns the torsional energy
    template<bool> double E_Torsion();
    double E_Torsion(bool gradients = true)
    {
      return gradients ? E_Torsion<true>() : E_Torsion<false>();
    }
    //! Returns the out-of-plane bending energy
    template<bool> double E_OOP();
    double E_OOP(bool gradients = true)
    {
      return gradients ? E_OOP<true>() : E_OOP<false>();
    }
    //! Returns the Van der Waals energy (Buckingham potential)
    template<bool> double E_VDW();
    double E_VDW(bool gradients = true)
    {
      return gradients ? E_VDW<true>() : E_VDW<false>();
    }
    //! Returns the dipole-dipole interaction energy
    template<bool> double E_Electrostatic();
    double E_Electrostatic(bool gradients = true)
    {
      return gradients ? E_Electrostatic<true>() : E_Electrostatic<false>();
    }

    //! Compare and print the numerical and analytical gradients
    bool ValidateGradients();

  }; // class MolgrForceFieldUFF

  void ClearMolgrUffAtomTypeAssignmentCache();
  std::tuple<std::size_t, std::size_t, std::size_t> MolgrUffAtomTypeAssignmentCacheInfo();

}// namespace OpenBabel

//! \file forcefieldUFF.h
//! \brief UFF force field
