# include "Box.H"
#include <cassert>
Box::Box()
{
  m_lowCorner = getZeros();
  m_highCorner= getOnes();
  m_highCorner *= -1;
  m_size = -1;
};
Box::Box(const Box& a_box)
{
  m_lowCorner = a_box.m_lowCorner;
  m_highCorner = a_box.m_highCorner;
  m_size = a_box.m_size;
};
Box::Box(const PPoint& a_lowCorner,const PPoint& a_highCorner)
{
 m_lowCorner = a_lowCorner;
 m_highCorner = a_highCorner;
 m_size = 1;
  for (int dir = 0;dir < DIM;dir++)
    {
      m_size=(m_highCorner[dir] - m_lowCorner[dir] + 1)*m_size;
    }
};
//Box returns a new box, should have put const to make clear.
Box Box::operator*(const Box& a_rightBox) const
{
  int newLo[DIM];
  int newHi[DIM];
  for (int i = 0 ; i < DIM; i++)
    {
      newLo[i] = m_lowCorner[i];
      if (m_lowCorner[i] < a_rightBox.m_lowCorner[i])
        {
          newLo[i] = a_rightBox.m_lowCorner[i];
        }
      newHi[i] = m_highCorner[i];
      if (m_highCorner[i] > a_rightBox.m_highCorner[i])
        {
          newHi[i] = a_rightBox.m_highCorner[i];
        }
    }
  PPoint newLowCorner(newLo);
  PPoint newHighCorner(newHi); 
  Box ret(newLowCorner,newHighCorner);
  return ret;
};
Box Box::shift(int a_direction, int a_offset)
{
  Box returnBox = Box(*this);
  returnBox.m_lowCorner += getUnitv(a_direction)*a_offset;
  returnBox.m_highCorner += getUnitv(a_direction)*a_offset;
  return returnBox;
};
Box Box::shift(const PPoint& a_pt)
{
  Box returnBox = Box(*this);
  returnBox.m_lowCorner += a_pt;
  returnBox.m_highCorner += a_pt;
  return returnBox;
};
Box Box::grow(int a_offset)
{
  Box returnBox(*this);
  returnBox.m_lowCorner -= getOnes()*a_offset;
  returnBox.m_highCorner += getOnes()*a_offset;
  returnBox.m_size = 1;
  for (int dir = 0;dir < DIM;dir++)
    {
      returnBox.m_size=
        (returnBox.m_highCorner[dir] - 
         returnBox.m_lowCorner[dir] + 1)*returnBox.m_size;
    }
  return returnBox;
};
Box Box::grow(const PPoint& a_offset)
{
  Box returnBox(*this);
  returnBox.m_lowCorner -= a_offset;
  returnBox.m_highCorner += a_offset;
   returnBox.m_size = 1;
  for (int dir = 0;dir < DIM;dir++)
    {
      returnBox.m_size=
        (returnBox.m_highCorner[dir] - 
         returnBox.m_lowCorner[dir] + 1)*returnBox.m_size;
    }
  return returnBox;
};
Box Box::coarsen(int a_nref)
{
  Box returnBox(*this);
  returnBox.m_lowCorner /=a_nref;
  returnBox.m_highCorner /=a_nref;
  returnBox.m_size = 1;
  for (int dir = 0;dir < DIM;dir++)
    {
      returnBox.m_size=
        (returnBox.m_highCorner[dir] - 
         returnBox.m_lowCorner[dir] + 1)*returnBox.m_size;
    }
  return returnBox;
};
Box Box::coarsen(const PPoint& a_pt)
{
  Box returnBox(*this);
  returnBox.m_lowCorner /=a_pt;
  returnBox.m_highCorner /= a_pt;
   returnBox.m_size = 1;
  for (int dir = 0;dir < DIM;dir++)
    {
      returnBox.m_size=
        (returnBox.m_highCorner[dir] - 
         returnBox.m_lowCorner[dir] + 1)*returnBox.m_size;
    }
  return returnBox;
};
Box Box::refine(int a_nref)
{
  Box returnBox(*this);
  returnBox.m_lowCorner *=a_nref;
  returnBox.m_highCorner *= a_nref;
  returnBox.m_size = 1;
  for (int dir = 0;dir < DIM;dir++)
    {
      returnBox.m_size=
        (returnBox.m_highCorner[dir] - 
         returnBox.m_lowCorner[dir] + 1)*returnBox.m_size;
    }
  return returnBox;
};
Box Box::refine(const PPoint& a_pt)
{
  Box returnBox(*this);
  returnBox.m_lowCorner *=a_pt;
  returnBox.m_highCorner *= a_pt;
  returnBox.m_size = 1;
  for (int dir = 0;dir < DIM;dir++)
    {
      returnBox.m_size=
        (returnBox.m_highCorner[dir] - 
         returnBox.m_lowCorner[dir] + 1)*returnBox.m_size;
    }  
  return returnBox;
};
// int Box::getIndex(const PPoint& a_pt) const
// {
//   int factor = 1;
//   int linIndex = a_pt[0] - m_lowCorner[0];
//   for (int i = 1;i < DIM;i++)
//     { 
//       factor = factor*(m_highCorner[i-1] - m_lowCorner[i-1]+1);
//       linIndex = linIndex + (a_pt[i] - m_lowCorner[i])*factor;
//     }
//   return linIndex;
// };
void Box::increment(PPoint& a_pt) const
{
  int k = getIndex(a_pt);
  assert(k < sizeOf());
  assert(k >= 0);
  if (k < sizeOf()-1) 
    {
      a_pt = getPoint(k+1);
    }
  else
    {
      a_pt+=getOnes();
    }
};
bool Box::notDone(const PPoint& a_pt) const
{
  int k = getIndex(a_pt);
  return (k < sizeOf());
};
PPoint Box::getPoint(int k) const
{
  assert(k >= 0);
  int tuple[DIM];
  for (int i = 0;i < DIM; i++)
    {      
      int factor = (m_highCorner[i] - m_lowCorner[i] + 1);
      int kred = k%factor;
      tuple[i] = kred + m_lowCorner[i];
      k = (k - kred)/factor;
    }
  PPoint pt(tuple);
  return pt;
};
bool Box::contains(const PPoint& a_pt) const
{
  int k=getIndex(a_pt);
  return (k < sizeOf());
};
