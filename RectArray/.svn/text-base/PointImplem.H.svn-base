#include <assert.h>

inline PPoint getOnes()
{
  int ones[DIM];
  for (int k = 0 ; k < DIM; k++)
    {
      ones[k] = 1;
    }
  return PPoint(ones);
};
inline PPoint getZeros()
{
  int zeros[DIM];
  for (int k = 0 ; k < DIM; k++)
    {
      zeros[k] = 0;
    }
  return PPoint(zeros);
};
// PPoint getUnitv(int a_dir)
// {
//   assert((a_dir < DIM) && (a_dir >= 0));
//   int unitv[DIM];
//   for (int k = 0 ; k < DIM; k++)
//     {
//       unitv[k] = 0;
//     }
//  unitv[a_dir] = 1;
//  return PPoint(unitv);
// };
PPoint::PPoint(){};
PPoint::PPoint(const int a_tuple[DIM])
{
  for (int j = 0; j < DIM; j++)
    {
      m_tuple[j]=a_tuple[j];
    }
};
PPoint::PPoint(const PPoint& a_pt)
{
  *this = a_pt;
};
PPoint PPoint::operator+(const PPoint& a_rhsPoint) const
{
  int tuple[DIM];
  for (int j = 0; j < DIM; j++)
    {
      tuple[j] = m_tuple[j]+a_rhsPoint.m_tuple[j]; //a_tuple
    }
  PPoint ret(tuple);
  return ret;
};
PPoint PPoint::operator-(const PPoint& a_rhsPoint) const
{
  int tuple[DIM];
  for (int j = 0; j < DIM; j++)
    {
      tuple[j] = m_tuple[j]-a_rhsPoint.m_tuple[j]; 
    }
  PPoint ret(tuple);
  return ret;
};
PPoint PPoint::operator*(const PPoint& a_rhsPoint) const
{
  int tuple[DIM];
  for (int j = 0; j < DIM; j++)
    {
      tuple[j] = m_tuple[j]*a_rhsPoint.m_tuple[j]; //a_tuple
    }
  PPoint ret(tuple);
  return ret;
}; 
PPoint PPoint::operator/(const PPoint& a_rhsPoint) const
{
  int tuple[DIM];
  for (int j = 0; j < DIM; j++)
    {
      assert(a_rhsPoint[j]!=0);
      if (m_tuple[j]*a_rhsPoint[j] >= 0)
        {
          tuple[j] = m_tuple[j]/a_rhsPoint[j];
        }
      else
        {
          tuple[j] = -(-m_tuple[j]/a_rhsPoint[j]);
        }
    }
  PPoint ret(tuple);
  return ret;

};
PPoint PPoint::operator*(int a_nref) const
{
  int tuple[DIM];
  for (int j = 0; j < DIM; j++)
    {
      tuple[j] = m_tuple[j]*a_nref;
    }
  PPoint ret(tuple);
  return ret;
};    
PPoint PPoint::operator/(int a_nref) const
{
  int tuple[DIM];
  assert(a_nref != 0);
  for (int j = 0; j < DIM; j++)
    {
      
      if (m_tuple[j]*a_nref >= 0)
        {
          tuple[j] = m_tuple[j]/a_nref;
        }
      else
        {
          tuple[j] = -(-m_tuple[j]/a_nref);
        }
    }
  PPoint ret(tuple);
  return ret;
};
void PPoint::operator+=(const PPoint& a_rhsPoint)
{
  for (int j = 0; j < DIM; j++)
    {
      m_tuple[j] = m_tuple[j]+a_rhsPoint.m_tuple[j]; 
    }
};
void PPoint::operator-=(const PPoint& a_rhsPoint)
{
  for (int j = 0; j < DIM; j++)
    {
      m_tuple[j] = m_tuple[j]-a_rhsPoint.m_tuple[j]; 
    }
};
void PPoint::operator*=(const PPoint& a_rhsPoint)
{
  for (int j = 0; j < DIM; j++)
    {
      m_tuple[j] = m_tuple[j]*a_rhsPoint.m_tuple[j]; 
    }
};
void PPoint::operator/=(const PPoint& a_rhsPoint)
{
  for (int j = 0; j < DIM; j++)
    {
      assert(a_rhsPoint[j]!=0);
      if (m_tuple[j]*a_rhsPoint[j] >= 0)
        {
          m_tuple[j] = m_tuple[j]/a_rhsPoint[j];
        }
      else
        {
          m_tuple[j] = -(-m_tuple[j]/a_rhsPoint[j]);
        }
    }
};
void PPoint::operator+=(int a_nref)
{
  for (int j = 0; j < DIM; j++)
    {
      m_tuple[j] = m_tuple[j]+a_nref; 
    }
};
void PPoint::operator-=(int a_nref)
{
  for (int j = 0; j < DIM; j++)
    {
      m_tuple[j] = m_tuple[j]-a_nref; 
    }
};
void PPoint::operator*=(int a_nref)
{
  for (int j = 0; j < DIM; j++)
    {
      m_tuple[j] = m_tuple[j]*a_nref; 
    }
};
void PPoint::operator/=(int a_nref)
{
  assert(a_nref!=0);
  for (int j = 0; j < DIM; j++)
    {
      
      if (m_tuple[j]*a_nref >= 0)
        {
          m_tuple[j] = m_tuple[j]/a_nref;
        }
      else
        {
          m_tuple[j] = -(-m_tuple[j]/a_nref);
        }
    }
};
bool PPoint::operator==(const PPoint& a_rhsPoint) const 
{
  bool ret = true;
  for (int j = 0; j < DIM; j++)
    {
      if (m_tuple[j] != a_rhsPoint[j]) ret=false;
    }
  return ret;
};

// const int& PPoint::operator[](const int& a_index) const
// {
//   assert(a_index < DIM);
//   assert(a_index >= 0);
//   return m_tuple[a_index];
// };
