/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkTypedArray.h.in

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkTypeUInt16Array - dynamic, self-adjusting array of vtkTypeUInt16
// .SECTION Description
// vtkTypeUInt16Array is an array of values of type vtkTypeUInt16.  It
// provides methods for insertion and retrieval of values and will
// automatically resize itself to hold new data.

#ifndef __vtkTypeUInt16Array_h
#define __vtkTypeUInt16Array_h

#include "vtkUnsignedShortArray.h"

class VTK_COMMON_EXPORT vtkTypeUInt16Array : public vtkUnsignedShortArray
{
public:
  static vtkTypeUInt16Array* New();
  vtkTypeMacro(vtkTypeUInt16Array,vtkUnsignedShortArray);
  void PrintSelf(ostream& os, vtkIndent indent);

protected:
  vtkTypeUInt16Array(vtkIdType numComp=1);
  ~vtkTypeUInt16Array();

private:
  vtkTypeUInt16Array(const vtkTypeUInt16Array&);  // Not implemented.
  void operator=(const vtkTypeUInt16Array&);  // Not implemented.
};

#endif
