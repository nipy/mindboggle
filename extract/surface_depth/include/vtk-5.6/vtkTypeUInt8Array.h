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
// .NAME vtkTypeUInt8Array - dynamic, self-adjusting array of vtkTypeUInt8
// .SECTION Description
// vtkTypeUInt8Array is an array of values of type vtkTypeUInt8.  It
// provides methods for insertion and retrieval of values and will
// automatically resize itself to hold new data.

#ifndef __vtkTypeUInt8Array_h
#define __vtkTypeUInt8Array_h

#include "vtkUnsignedCharArray.h"

class VTK_COMMON_EXPORT vtkTypeUInt8Array : public vtkUnsignedCharArray
{
public:
  static vtkTypeUInt8Array* New();
  vtkTypeMacro(vtkTypeUInt8Array,vtkUnsignedCharArray);
  void PrintSelf(ostream& os, vtkIndent indent);

protected:
  vtkTypeUInt8Array(vtkIdType numComp=1);
  ~vtkTypeUInt8Array();

private:
  vtkTypeUInt8Array(const vtkTypeUInt8Array&);  // Not implemented.
  void operator=(const vtkTypeUInt8Array&);  // Not implemented.
};

#endif
