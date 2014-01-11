#include "pot_2d.hpp"

class C2DFFTPotential : public C2DPotential
{
public:
  C2DFFTPotential(const ConfigReaderPtr &configReader);
  virtual void DisplayParams();
  virtual void makeSlices(int nlayer, char *fileName, atom *center);
  virtual void AddAtomToSlices(std::vector<atom>::iterator &atom, 
                               float_tt atomX, float_tt atomY, float_tt atomZ);
protected:
  virtual void AddAtomPeriodic(std::vector<atom>::iterator &atom, 
                         float_tt atomBoxX, unsigned int ix, 
                         float_tt atomBoxY, unsigned int iy, 
                         float_tt atomZ);
  virtual void AddAtomNonPeriodic(std::vector<atom>::iterator &atom, 
                         float_tt atomBoxX, unsigned int ix, 
                         float_tt atomBoxY, unsigned int iy, 
                         float_tt atomZ);
  complex_tt *GetAtomPotential2D(int Znum, double B);
private:
  unsigned m_nyAtBox, m_nxyAtBox, m_nyAtBox2, m_nxyAtBox2; //Size of atom box in pixels
};
