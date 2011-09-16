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
// .NAME vtkTypeFloat32Array - dynamic, self-adjusting array of vtkTypeFloat32
// .SECTION Description
// vtkTypeFloat32Array is an array of values of type vtkTypeFloat32.  It
// provides methods for insertion and retrieval of values and will
// automatically resize itself to hold new data.

#ifndef __vtkTypeFloat32Array_h
#define __vtkTypeFloat32Array_h

#include "vtkFloatArray.h"

class VTK_COMMON_EXPORT vtkTypeFloat32Array : public vtkFloatArray
{
public:
  static vtkTypeFloat32Array* New();
  vtkTypeMacro(vtkTypeFloat32Array,vtkFloatArray);
  void PrintSelf(ostream& os, vtkIndent indent);

protected:
  vtkTypeFloat32Array(vtkIdType numComp=1);
  ~vtkTypeFloat32Array();

private:
  vtkTypeFloat32Array(const vtkTypeFloat32Array&);  // Not implemented.
  void operator=(const vtkTypeFloat32Array&);  // Not implemented.
};

#endif
