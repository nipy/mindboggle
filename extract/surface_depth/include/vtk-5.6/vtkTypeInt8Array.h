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
// .NAME vtkTypeInt8Array - dynamic, self-adjusting array of vtkTypeInt8
// .SECTION Description
// vtkTypeInt8Array is an array of values of type vtkTypeInt8.  It
// provides methods for insertion and retrieval of values and will
// automatically resize itself to hold new data.

#ifndef __vtkTypeInt8Array_h
#define __vtkTypeInt8Array_h

#include "vtkCharArray.h"

class VTK_COMMON_EXPORT vtkTypeInt8Array : public vtkCharArray
{
public:
  static vtkTypeInt8Array* New();
  vtkTypeMacro(vtkTypeInt8Array,vtkCharArray);
  void PrintSelf(ostream& os, vtkIndent indent);

protected:
  vtkTypeInt8Array(vtkIdType numComp=1);
  ~vtkTypeInt8Array();

private:
  vtkTypeInt8Array(const vtkTypeInt8Array&);  // Not implemented.
  void operator=(const vtkTypeInt8Array&);  // Not implemented.
};

#endif
