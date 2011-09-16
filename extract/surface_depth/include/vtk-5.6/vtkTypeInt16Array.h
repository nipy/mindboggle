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
// .NAME vtkTypeInt16Array - dynamic, self-adjusting array of vtkTypeInt16
// .SECTION Description
// vtkTypeInt16Array is an array of values of type vtkTypeInt16.  It
// provides methods for insertion and retrieval of values and will
// automatically resize itself to hold new data.

#ifndef __vtkTypeInt16Array_h
#define __vtkTypeInt16Array_h

#include "vtkShortArray.h"

class VTK_COMMON_EXPORT vtkTypeInt16Array : public vtkShortArray
{
public:
  static vtkTypeInt16Array* New();
  vtkTypeMacro(vtkTypeInt16Array,vtkShortArray);
  void PrintSelf(ostream& os, vtkIndent indent);

protected:
  vtkTypeInt16Array(vtkIdType numComp=1);
  ~vtkTypeInt16Array();

private:
  vtkTypeInt16Array(const vtkTypeInt16Array&);  // Not implemented.
  void operator=(const vtkTypeInt16Array&);  // Not implemented.
};

#endif
