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
// .NAME vtkTypeUInt64Array - dynamic, self-adjusting array of vtkTypeUInt64
// .SECTION Description
// vtkTypeUInt64Array is an array of values of type vtkTypeUInt64.  It
// provides methods for insertion and retrieval of values and will
// automatically resize itself to hold new data.

#ifndef __vtkTypeUInt64Array_h
#define __vtkTypeUInt64Array_h

#include "vtkUnsignedLongLongArray.h"

class VTK_COMMON_EXPORT vtkTypeUInt64Array : public vtkUnsignedLongLongArray
{
public:
  static vtkTypeUInt64Array* New();
  vtkTypeMacro(vtkTypeUInt64Array,vtkUnsignedLongLongArray);
  void PrintSelf(ostream& os, vtkIndent indent);

protected:
  vtkTypeUInt64Array(vtkIdType numComp=1);
  ~vtkTypeUInt64Array();

private:
  vtkTypeUInt64Array(const vtkTypeUInt64Array&);  // Not implemented.
  void operator=(const vtkTypeUInt64Array&);  // Not implemented.
};

#endif
