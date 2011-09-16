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
// .NAME vtkTypeUInt32Array - dynamic, self-adjusting array of vtkTypeUInt32
// .SECTION Description
// vtkTypeUInt32Array is an array of values of type vtkTypeUInt32.  It
// provides methods for insertion and retrieval of values and will
// automatically resize itself to hold new data.

#ifndef __vtkTypeUInt32Array_h
#define __vtkTypeUInt32Array_h

#include "vtkUnsignedIntArray.h"

class VTK_COMMON_EXPORT vtkTypeUInt32Array : public vtkUnsignedIntArray
{
public:
  static vtkTypeUInt32Array* New();
  vtkTypeMacro(vtkTypeUInt32Array,vtkUnsignedIntArray);
  void PrintSelf(ostream& os, vtkIndent indent);

protected:
  vtkTypeUInt32Array(vtkIdType numComp=1);
  ~vtkTypeUInt32Array();

private:
  vtkTypeUInt32Array(const vtkTypeUInt32Array&);  // Not implemented.
  void operator=(const vtkTypeUInt32Array&);  // Not implemented.
};

#endif
