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
// .NAME vtkTypeFloat64Array - dynamic, self-adjusting array of vtkTypeFloat64
// .SECTION Description
// vtkTypeFloat64Array is an array of values of type vtkTypeFloat64.  It
// provides methods for insertion and retrieval of values and will
// automatically resize itself to hold new data.

#ifndef __vtkTypeFloat64Array_h
#define __vtkTypeFloat64Array_h

#include "vtkDoubleArray.h"

class VTK_COMMON_EXPORT vtkTypeFloat64Array : public vtkDoubleArray
{
public:
  static vtkTypeFloat64Array* New();
  vtkTypeMacro(vtkTypeFloat64Array,vtkDoubleArray);
  void PrintSelf(ostream& os, vtkIndent indent);

protected:
  vtkTypeFloat64Array(vtkIdType numComp=1);
  ~vtkTypeFloat64Array();

private:
  vtkTypeFloat64Array(const vtkTypeFloat64Array&);  // Not implemented.
  void operator=(const vtkTypeFloat64Array&);  // Not implemented.
};

#endif
